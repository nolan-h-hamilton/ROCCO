#include "wls_backend.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static void rocco_swap_double(double *left, double *right)
{
    double temp = *left;
    *left = *right;
    *right = temp;
}

static size_t rocco_partition_f64(
    double *values,
    size_t left,
    size_t right,
    size_t pivot_index)
{
    double pivot_value = values[pivot_index];
    size_t store_index = left;
    size_t idx = 0;

    rocco_swap_double(&values[pivot_index], &values[right]);
    for (idx = left; idx < right; ++idx)
    {
        if (values[idx] < pivot_value)
        {
            rocco_swap_double(&values[store_index], &values[idx]);
            ++store_index;
        }
    }
    rocco_swap_double(&values[right], &values[store_index]);
    return store_index;
}

static void rocco_partition3_f64(
    double *values,
    size_t left,
    size_t right,
    size_t pivot_index,
    size_t *equal_left_out,
    size_t *equal_right_out)
{
    double pivot_value = values[pivot_index];
    size_t less = left;
    size_t idx = left;
    size_t greater = right;

    while (idx <= greater)
    {
        if (values[idx] < pivot_value)
        {
            rocco_swap_double(&values[less], &values[idx]);
            ++less;
            ++idx;
        }
        else if (values[idx] > pivot_value)
        {
            rocco_swap_double(&values[idx], &values[greater]);
            --greater;
        }
        else
        {
            ++idx;
        }
    }

    *equal_left_out = less;
    *equal_right_out = greater;
}

static double rocco_select_kth_f64(
    double *values,
    size_t left,
    size_t right,
    size_t kth_index)
{
    size_t pivot_index = 0;
    size_t equal_left = 0;
    size_t equal_right = 0;

    while (left < right)
    {
        pivot_index = left + ((right - left) / 2U);
        rocco_partition3_f64(
            values,
            left,
            right,
            pivot_index,
            &equal_left,
            &equal_right);
        if (kth_index < equal_left)
        {
            right = equal_left - 1U;
        }
        else if (kth_index > equal_right)
        {
            left = equal_right + 1U;
        }
        else
        {
            return values[kth_index];
        }
    }
    return values[left];
}

static double rocco_median_f64(double *values, size_t value_count)
{
    size_t upper_index = 0;
    double upper = 0.0;
    double lower = 0.0;

    if (value_count == 0U)
    {
        return 0.0;
    }
    if (value_count == 1U)
    {
        return values[0];
    }

    upper_index = value_count / 2U;
    if ((value_count & 1U) == 1U)
    {
        return rocco_select_kth_f64(values, 0U, value_count - 1U, upper_index);
    }

    upper = rocco_select_kth_f64(values, 0U, value_count - 1U, upper_index);
    lower = rocco_select_kth_f64(values, 0U, upper_index - 1U, upper_index - 1U);
    return 0.5 * (lower + upper);
}

static double rocco_small_median_f64(double *values, size_t value_count)
{
    size_t idx = 0;
    size_t inner = 0;
    double key = 0.0;

    if (value_count == 0U)
    {
        return 0.0;
    }
    for (idx = 1U; idx < value_count; ++idx)
    {
        key = values[idx];
        inner = idx;
        while (inner > 0U && values[inner - 1U] > key)
        {
            values[inner] = values[inner - 1U];
            --inner;
        }
        values[inner] = key;
    }

    if ((value_count & 1U) == 1U)
    {
        return values[value_count / 2U];
    }
    return 0.5 * (values[(value_count / 2U) - 1U] + values[value_count / 2U]);
}

static double rocco_robust_scale_f64(double *work_buffer, size_t value_count)
{
    size_t idx = 0;
    double median = 0.0;
    double mad = 0.0;

    if (value_count == 0U)
    {
        return 1.0e-6;
    }
    median = rocco_median_f64(work_buffer, value_count);
    for (idx = 0U; idx < value_count; ++idx)
    {
        work_buffer[idx] = fabs(work_buffer[idx] - median);
    }
    mad = rocco_median_f64(work_buffer, value_count);
    mad *= 1.4826;
    if (!(mad > 1.0e-6))
    {
        return 1.0e-6;
    }
    return mad;
}

int rocco_score_centered_wls_f64(
    const double *centered_matrix,
    size_t sample_count,
    size_t locus_count,
    double lower_bound_z,
    double prior_df,
    double *sample_weights_out,
    double *mean_out,
    double *standard_error_out,
    double *z_scores_out,
    double *scores_out)
{
    double *consensus = NULL;
    double *sample_scales = NULL;
    double *weighted_var = NULL;
    double *work_buffer = NULL;
    double *column_buffer = NULL;
    size_t sample_idx = 0U;
    size_t locus_idx = 0U;
    size_t positive_count = 0U;
    double scale_floor = 1.0e-6;
    double sum_w = 0.0;
    double sum_w2 = 0.0;
    double effective_n = 1.0;
    double variance_denom = 1.0;
    double global_var = 1.0e-6;

    if (centered_matrix == NULL || sample_weights_out == NULL || mean_out == NULL ||
        standard_error_out == NULL || z_scores_out == NULL || scores_out == NULL)
    {
        return -2;
    }
    if (sample_count == 0U || locus_count == 0U)
    {
        return -2;
    }

    consensus = (double *)malloc(locus_count * sizeof(double));
    sample_scales = (double *)malloc(sample_count * sizeof(double));
    weighted_var = (double *)malloc(locus_count * sizeof(double));
    work_buffer = (double *)malloc(locus_count * sizeof(double));
    column_buffer = (double *)malloc(sample_count * sizeof(double));
    if (consensus == NULL || sample_scales == NULL || weighted_var == NULL ||
        work_buffer == NULL || column_buffer == NULL)
    {
        free(consensus);
        free(sample_scales);
        free(weighted_var);
        free(work_buffer);
        free(column_buffer);
        return -1;
    }

    if (sample_count == 1U)
    {
        for (locus_idx = 0U; locus_idx < locus_count; ++locus_idx)
        {
            mean_out[locus_idx] = centered_matrix[locus_idx];
            work_buffer[locus_idx] = centered_matrix[locus_idx];
        }
        scale_floor = rocco_robust_scale_f64(work_buffer, locus_count);
        sample_weights_out[0] = 1.0;
        for (locus_idx = 0U; locus_idx < locus_count; ++locus_idx)
        {
            standard_error_out[locus_idx] = scale_floor;
            z_scores_out[locus_idx] = mean_out[locus_idx] / fmax(scale_floor, 1.0e-8);
            scores_out[locus_idx] = mean_out[locus_idx] - (lower_bound_z * scale_floor);
        }
        free(consensus);
        free(sample_scales);
        free(weighted_var);
        free(work_buffer);
        free(column_buffer);
        return 0;
    }

    for (locus_idx = 0U; locus_idx < locus_count; ++locus_idx)
    {
        for (sample_idx = 0U; sample_idx < sample_count; ++sample_idx)
        {
            column_buffer[sample_idx] =
                centered_matrix[(sample_idx * locus_count) + locus_idx];
        }
        consensus[locus_idx] = rocco_small_median_f64(column_buffer, sample_count);
    }

    for (sample_idx = 0U; sample_idx < sample_count; ++sample_idx)
    {
        for (locus_idx = 0U; locus_idx < locus_count; ++locus_idx)
        {
            work_buffer[locus_idx] = fabs(
                centered_matrix[(sample_idx * locus_count) + locus_idx] -
                consensus[locus_idx]);
        }
        sample_scales[sample_idx] = rocco_median_f64(work_buffer, locus_count);
    }

    positive_count = 0U;
    for (sample_idx = 0U; sample_idx < sample_count; ++sample_idx)
    {
        if (sample_scales[sample_idx] > 0.0)
        {
            column_buffer[positive_count] = sample_scales[sample_idx];
            ++positive_count;
        }
    }
    if (positive_count > 0U)
    {
        scale_floor = rocco_small_median_f64(column_buffer, positive_count);
    }
    if (!(scale_floor > 1.0e-6))
    {
        scale_floor = 1.0e-6;
    }

    for (sample_idx = 0U; sample_idx < sample_count; ++sample_idx)
    {
        double sample_var = fmax(sample_scales[sample_idx], scale_floor);
        sample_var *= sample_var;
        sample_weights_out[sample_idx] = 1.0 / sample_var;
        sum_w += sample_weights_out[sample_idx];
        sum_w2 += sample_weights_out[sample_idx] * sample_weights_out[sample_idx];
    }

    effective_n = (sum_w * sum_w) / fmax(sum_w2, 1.0e-12);
    variance_denom = fmax(sum_w - (sum_w2 / sum_w), 1.0e-8);

    for (locus_idx = 0U; locus_idx < locus_count; ++locus_idx)
    {
        double weighted_sum = 0.0;
        double var_sum = 0.0;
        for (sample_idx = 0U; sample_idx < sample_count; ++sample_idx)
        {
            double value = centered_matrix[(sample_idx * locus_count) + locus_idx];
            weighted_sum += sample_weights_out[sample_idx] * value;
        }
        mean_out[locus_idx] = weighted_sum / sum_w;
        for (sample_idx = 0U; sample_idx < sample_count; ++sample_idx)
        {
            double resid =
                centered_matrix[(sample_idx * locus_count) + locus_idx] -
                mean_out[locus_idx];
            var_sum += sample_weights_out[sample_idx] * resid * resid;
        }
        weighted_var[locus_idx] = var_sum / variance_denom;
    }

    positive_count = 0U;
    for (locus_idx = 0U; locus_idx < locus_count; ++locus_idx)
    {
        if (weighted_var[locus_idx] > 0.0)
        {
            work_buffer[positive_count] = weighted_var[locus_idx];
            ++positive_count;
        }
    }
    if (positive_count > 0U)
    {
        global_var = rocco_median_f64(work_buffer, positive_count);
    }
    else
    {
        memcpy(work_buffer, weighted_var, locus_count * sizeof(double));
        global_var = rocco_median_f64(work_buffer, locus_count);
    }
    if (!(global_var > 1.0e-6))
    {
        global_var = 1.0e-6;
    }

    for (locus_idx = 0U; locus_idx < locus_count; ++locus_idx)
    {
        double posterior_var =
            ((fmax(effective_n - 1.0, 0.0) * weighted_var[locus_idx]) +
             (prior_df * global_var)) /
            fmax((effective_n - 1.0) + prior_df, 1.0);
        standard_error_out[locus_idx] = sqrt(posterior_var / sum_w);
        z_scores_out[locus_idx] =
            mean_out[locus_idx] / fmax(standard_error_out[locus_idx], 1.0e-8);
        scores_out[locus_idx] =
            mean_out[locus_idx] - (lower_bound_z * standard_error_out[locus_idx]);
    }

    free(consensus);
    free(sample_scales);
    free(weighted_var);
    free(work_buffer);
    free(column_buffer);
    return 0;
}
