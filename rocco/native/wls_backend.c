#include "wls_backend.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* swap two doubles (comparators needed for qsort, etc.) */
static void rocco_swap_double(double *left, double *right)
{
    double temp = *left;
    *left = *right;
    *right = temp;
}

/* partition an array of doubles around a pivot (for quick-select/sort) */
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
            /* move to the left partition */
            rocco_swap_double(&values[store_index], &values[idx]);
            ++store_index;
        }
    }
    rocco_swap_double(&values[right], &values[store_index]);
    return store_index;
}

/* three-way partition an array of doubles around a pivot (for quick-select/sort, robustness to equal values) */
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

    /** Cases: (i) value < pivot_value, (ii) value > pivot_value, (iii) value == pivot_value */

    while (idx <= greater)
    {
        if (values[idx] < pivot_value)
        {
            /* move to the left partition */
            rocco_swap_double(&values[less], &values[idx]);
            ++less;
            ++idx;
        }
        else if (values[idx] > pivot_value)
        {
            /* move to the right partition */
            rocco_swap_double(&values[idx], &values[greater]);
            --greater;
        }
        else
        {
            /* equal value, just move forward */
            ++idx;
        }
    }
    /* point at boundaries of the equal-val partition */
    *equal_left_out = less;
    *equal_right_out = greater;
}

/* select the k-th smallest element in an array of doubles */
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
    // insertion for small array medians (stable/fast relative to qselect)
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

typedef struct rocco_xy_pair
{
    double x;
    double y;
} rocco_xy_pair;

static int rocco_compare_xy_pair_by_x(
    const void *left,
    const void *right)
{
    const rocco_xy_pair *left_pair = (const rocco_xy_pair *)left;
    const rocco_xy_pair *right_pair = (const rocco_xy_pair *)right;
    if (left_pair->x < right_pair->x)
    {
        return -1;
    }
    if (left_pair->x > right_pair->x)
    {
        return 1;
    }
    if (left_pair->y < right_pair->y)
    {
        return -1;
    }
    if (left_pair->y > right_pair->y)
    {
        return 1;
    }
    return 0;
}

static size_t rocco_resolve_spatial_window(
    size_t locus_count,
    int requested_window)
{
    size_t window = 0U;
    /* odd window or nothing */
    if (locus_count < 5U)
    {
        return 0U;
    }
    window = requested_window > 0 ? (size_t)requested_window : 31U;
    if (window < 5U)
    {
        window = 5U;
    }
    if (window > locus_count)
    {
        window = locus_count;
    }
    if ((window & 1U) == 0U)
    {
        window = (window == locus_count) ? (window - 1U) : (window + 1U);
    }
    if (window < 5U)
    {
        return 0U;
    }
    return window;
}

static int rocco_pava_f64(
    const double *values,
    const double *weights,
    size_t value_count,
    double *out)
{
    /* pool bins until the trend stops dropping */
    double *block_values = NULL;
    double *block_weights = NULL;
    size_t *block_lengths = NULL;
    size_t block_count = 0U;
    size_t idx = 0U;
    size_t cursor = 0U;

    if (values == NULL || weights == NULL || out == NULL)
    {
        return -2;
    }
    if (value_count == 0U)
    {
        return 0;
    }

    block_values = (double *)malloc(value_count * sizeof(double));
    block_weights = (double *)malloc(value_count * sizeof(double));
    block_lengths = (size_t *)malloc(value_count * sizeof(size_t));
    if (block_values == NULL || block_weights == NULL || block_lengths == NULL)
    {
        free(block_values);
        free(block_weights);
        free(block_lengths);
        return -1;
    }

    for (idx = 0U; idx < value_count; ++idx)
    {
        block_values[block_count] = values[idx];
        block_weights[block_count] = fmax(weights[idx], 1.0e-8);
        block_lengths[block_count] = 1U;
        ++block_count;

        /* merge until monotone */
        while (block_count >= 2U &&
               block_values[block_count - 2U] > block_values[block_count - 1U])
        {
            /* merge the last two blocks */
            double total_weight =
                block_weights[block_count - 2U] + block_weights[block_count - 1U];
            double merged_value =
                ((block_values[block_count - 2U] * block_weights[block_count - 2U]) +
                 (block_values[block_count - 1U] * block_weights[block_count - 1U])) /
                total_weight;
            /* write back to the second-to-last block and drop the last */
            block_values[block_count - 2U] = merged_value;
            block_weights[block_count - 2U] = total_weight;
            block_lengths[block_count - 2U] += block_lengths[block_count - 1U];
            --block_count;
        }
    }

    /* now we expand blocks back into the output array */
    cursor = 0U;
    for (idx = 0U; idx < block_count; ++idx)
    {
        size_t repeat_idx = 0U;
        for (repeat_idx = 0U; repeat_idx < block_lengths[idx]; ++repeat_idx)
        {
            out[cursor] = block_values[idx];
            ++cursor;
        }
    }

    free(block_values);
    free(block_weights);
    free(block_lengths);
    return 0;
}

/** linear interpolation between grid pts for monotone fit */
static double rocco_linear_interp_f64(
    const double *x_values,
    const double *y_values,
    size_t value_count,
    double target_x)
{
    size_t left = 0U;
    size_t right = 0U;
    size_t mid = 0U;
    double x_left = 0.0;
    double x_right = 0.0;
    double weight = 0.0;

    if (value_count == 0U)
    {
        return 1.0e-8;
    }
    if (value_count == 1U || target_x <= x_values[0])
    {
        return y_values[0];
    }
    if (target_x >= x_values[value_count - 1U])
    {
        return y_values[value_count - 1U];
    }

    left = 0U;
    right = value_count - 1U;
    /* bisect */
    while ((right - left) > 1U)
    {
        mid = left + ((right - left) / 2U);
        if (x_values[mid] <= target_x)
        {
            left = mid;
        }
        else
        {
            right = mid;
        }
    }

    x_left = x_values[left];
    x_right = x_values[right];
    if (x_right <= x_left)
    {
        return fmax(y_values[right], y_values[left]);
    }
    weight = (target_x - x_left) / (x_right - x_left);
    return y_values[left] + (weight * (y_values[right] - y_values[left]));
}

/* fit |mean| ~ variance trend */
static int rocco_fit_monotone_variance_trend_f64(
    const double *covariate,
    const double *raw_variance,
    size_t locus_count,
    double *trend_out)
{
    rocco_xy_pair *pairs = NULL;
    double *bin_cov = NULL;
    double *bin_var = NULL;
    double *bin_weight = NULL;
    double *fit_var = NULL;
    double *work_buffer = NULL;
    double *knot_cov = NULL;
    double *knot_var = NULL;
    size_t finite_count = 0U;
    size_t idx = 0U;
    size_t bin_count = 0U;
    size_t used_bins = 0U;
    size_t knot_count = 0U;
    double fallback = 1.0e-6;

    if (covariate == NULL || raw_variance == NULL || trend_out == NULL)
    {
        return -2;
    }

    pairs = (rocco_xy_pair *)malloc(locus_count * sizeof(rocco_xy_pair));
    work_buffer = (double *)malloc(locus_count * sizeof(double));
    if (pairs == NULL || work_buffer == NULL)
    {
        free(pairs);
        free(work_buffer);
        return -1;
    }

    for (idx = 0U; idx < locus_count; ++idx)
    {
        if (isfinite(covariate[idx]) && isfinite(raw_variance[idx]))
        {
            pairs[finite_count].x = fabs(covariate[idx]);
            pairs[finite_count].y = fmax(raw_variance[idx], 1.0e-8);
            work_buffer[finite_count] = pairs[finite_count].y;
            ++finite_count;
        }
    }
    if (finite_count > 0U)
    {
        fallback = fmax(rocco_median_f64(work_buffer, finite_count), 1.0e-8);
    }
    if (finite_count < 4U)
    {
        for (idx = 0U; idx < locus_count; ++idx)
        {
            trend_out[idx] = fallback;
        }
        free(pairs);
        free(work_buffer);
        return 0;
    }

    qsort(pairs, finite_count, sizeof(rocco_xy_pair), rocco_compare_xy_pair_by_x);
    /* start with coarse binning */
    bin_count = (size_t)fmax(4.0, floor(1.0 + (log((double)finite_count + 1.0) / log(2.0))));
    bin_cov = (double *)malloc(bin_count * sizeof(double));
    bin_var = (double *)malloc(bin_count * sizeof(double));
    bin_weight = (double *)malloc(bin_count * sizeof(double));
    fit_var = (double *)malloc(bin_count * sizeof(double));
    knot_cov = (double *)malloc(bin_count * sizeof(double));
    knot_var = (double *)malloc(bin_count * sizeof(double));
    if (bin_cov == NULL || bin_var == NULL || bin_weight == NULL || fit_var == NULL ||
        knot_cov == NULL || knot_var == NULL)
    {
        free(pairs);
        free(work_buffer);
        free(bin_cov);
        free(bin_var);
        free(bin_weight);
        free(fit_var);
        free(knot_cov);
        free(knot_var);
        return -1;
    }

    used_bins = 0U;
    for (idx = 0U; idx < bin_count; ++idx)
    {
        /* pick out representative pair from from each bin */
        size_t left = (idx * finite_count) / bin_count;
        size_t right = ((idx + 1U) * finite_count) / bin_count;
        size_t width = 0U;
        if (right <= left)
        {
            continue;
        }
        width = right - left;
        if ((width & 1U) == 1U)
        {
            bin_cov[used_bins] = pairs[left + (width / 2U)].x;
        }
        else
        {
            bin_cov[used_bins] =
                0.5 * (pairs[left + (width / 2U) - 1U].x + pairs[left + (width / 2U)].x);
        }
        for (size_t inner = 0U; inner < width; ++inner)
        {
            work_buffer[inner] = pairs[left + inner].y;
        }
        bin_var[used_bins] = rocco_median_f64(work_buffer, width);
        bin_weight[used_bins] = (double)width;
        ++used_bins;
    }

    if (used_bins == 0U)
    {
        for (idx = 0U; idx < locus_count; ++idx)
        {
            trend_out[idx] = fallback;
        }
        free(pairs);
        free(work_buffer);
        free(bin_cov);
        free(bin_var);
        free(bin_weight);
        free(fit_var);
        free(knot_cov);
        free(knot_var);
        return 0;
    }
    if (used_bins == 1U)
    {
        double constant = fmax(bin_var[0], 1.0e-8);
        for (idx = 0U; idx < locus_count; ++idx)
        {
            trend_out[idx] = constant;
        }
        free(pairs);
        free(work_buffer);
        free(bin_cov);
        free(bin_var);
        free(bin_weight);
        free(fit_var);
        free(knot_cov);
        free(knot_var);
        return 0;
    }

    if (rocco_pava_f64(bin_var, bin_weight, used_bins, fit_var) != 0)
    {
        free(pairs);
        free(work_buffer);
        free(bin_cov);
        free(bin_var);
        free(bin_weight);
        free(fit_var);
        free(knot_cov);
        free(knot_var);
        return -1;
    }

    knot_count = 0U;
    for (idx = 0U; idx < used_bins; ++idx)
    {
        double cov_value = bin_cov[idx];
        double var_value = fmax(fit_var[idx], 1.0e-8);
        if (knot_count > 0U && cov_value <= knot_cov[knot_count - 1U])
        {
            knot_var[knot_count - 1U] = fmax(knot_var[knot_count - 1U], var_value);
            continue;
        }
        knot_cov[knot_count] = cov_value;
        knot_var[knot_count] = var_value;
        ++knot_count;
    }

    if (knot_count == 0U)
    {
        for (idx = 0U; idx < locus_count; ++idx)
        {
            trend_out[idx] = fallback;
        }
    }
    else if (knot_count == 1U)
    {
        double constant = fmax(knot_var[0], 1.0e-8);
        for (idx = 0U; idx < locus_count; ++idx)
        {
            trend_out[idx] = constant;
        }
    }
    else
    {
        for (idx = 0U; idx < locus_count; ++idx)
        {
            if (!isfinite(covariate[idx]))
            {
                trend_out[idx] = fallback;
                continue;
            }
            trend_out[idx] = fmax(
                rocco_linear_interp_f64(knot_cov, knot_var, knot_count, fabs(covariate[idx])),
                1.0e-8);
        }
    }

    free(pairs);
    free(work_buffer);
    free(bin_cov);
    free(bin_var);
    free(bin_weight);
    free(fit_var);
    free(knot_cov);
    free(knot_var);
    return 0;
}

static int rocco_rolling_ar1_innovation_variance_f64(
    const double *values,
    size_t value_count,
    size_t window,
    double *variance_out)
{
    /* local munc track */
    double *var_at_start = NULL;
    size_t half_window = 0U;
    size_t max_start = 0U;
    size_t start_idx = 0U;
    size_t region_idx = 0U;
    double sum_y = 0.0;
    double sum_sq_y = 0.0;
    double sum_lag_prod = 0.0;
    double window_double = 0.0;
    double pair_count = 0.0;

    if (values == NULL || variance_out == NULL)
    {
        return -2;
    }
    if (value_count == 0U)
    {
        return -2;
    }

    window = rocco_resolve_spatial_window(value_count, (int)window);
    if (window == 0U || value_count < 4U)
    {
        memset(variance_out, 0, value_count * sizeof(double));
        return 0;
    }

    half_window = window / 2U;
    max_start = value_count - window;
    var_at_start = (double *)malloc((max_start + 1U) * sizeof(double));
    if (var_at_start == NULL)
    {
        return -1;
    }

    for (region_idx = 0U; region_idx < window; ++region_idx)
    {
        double current_value = values[region_idx];
        sum_y += current_value;
        sum_sq_y += current_value * current_value;
        if (region_idx < (window - 1U))
        {
            sum_lag_prod += current_value * values[region_idx + 1U];
        }
    }

    window_double = (double)window;
    pair_count = (double)(window - 1U);
    for (start_idx = 0U; start_idx <= max_start; ++start_idx)
    {
        double leaving_value = values[start_idx];
        double entering_value = values[start_idx + window - 1U];
        double sum_x_seq = sum_y - entering_value;
        double sum_y_seq = sum_y - leaving_value;
        double mean_all = sum_y / window_double;
        double gamma0_num = sum_sq_y - (window_double * mean_all * mean_all);
        double gamma1_num = 0.0;
        double lambda_eff = 1.0 / (window_double + 1.0);
        double scale_floor = 0.0;
        double denom = 0.0;
        double eps = 0.0;
        double beta1 = 0.0;
        double gamma0 = 0.0;
        double one_minus_beta_sq = 0.0;
        if (gamma0_num < 0.0)
        {
            gamma0_num = 0.0;
        }

        gamma1_num =
            sum_lag_prod - (mean_all * sum_x_seq) - (mean_all * sum_y_seq) +
            (pair_count * mean_all * mean_all);
        scale_floor = 1.0e-4 * (gamma0_num + 1.0);
        denom = (gamma0_num * (1.0 + lambda_eff)) + scale_floor;
        eps = 1.0e-12 * (gamma0_num + 1.0);
        if (denom > eps)
        {
            beta1 = gamma1_num / denom;
        }
        if (beta1 > 0.99)
        {
            beta1 = 0.99;
        }
        else if (beta1 < 0.0)
        {
            beta1 = 0.0;
        }
        gamma0 = gamma0_num / window_double;
        one_minus_beta_sq = 1.0 - (beta1 * beta1);
        if (one_minus_beta_sq < 0.0)
        {
            one_minus_beta_sq = 0.0;
        }
        var_at_start[start_idx] = fmax(gamma0 * one_minus_beta_sq, 0.0);

        if (start_idx < max_start)
        {
            double next_value = values[start_idx + window];
            double next_lag_left = values[start_idx + window - 1U];
            double next_lag_right = values[start_idx + 1U];
            sum_y = (sum_y - leaving_value) + next_value;
            sum_sq_y =
                sum_sq_y - (leaving_value * leaving_value) + (next_value * next_value);
            sum_lag_prod =
                sum_lag_prod - (leaving_value * next_lag_right) +
                (next_lag_left * next_value);
        }
    }

    for (region_idx = 0U; region_idx < value_count; ++region_idx)
    {
        ptrdiff_t candidate = (ptrdiff_t)region_idx - (ptrdiff_t)half_window;
        if (candidate < 0)
        {
            candidate = 0;
        }
        else if ((size_t)candidate > max_start)
        {
            candidate = (ptrdiff_t)max_start;
        }
        variance_out[region_idx] = var_at_start[candidate];
    }

    free(var_at_start);
    return 0;
}

int rocco_score_centered_wls_f64(
    const double *centered_matrix,
    size_t sample_count,
    size_t locus_count,
    double lower_bound_z,
    double prior_df,
    double min_effect,
    int use_min_effect,
    int spatial_window,
    double precision_floor_ratio,
    double *mean_out,
    double *raw_variance_out,
    double *prior_variance_out,
    double *moderated_variance_out,
    double *standard_error_out,
    double *scores_out,
    double *degrees_of_freedom_out,
    int *resolved_window_out)
{
    double *obs_variance = NULL;
    double *prior_track = NULL;
    double *weighted_sum = NULL;
    double *precision_sum = NULL;
    double *raw_precision_sum = NULL;
    double *prior_precision_sum = NULL;
    double *work_buffer = NULL;
    size_t sample_idx = 0U;
    size_t locus_idx = 0U;
    size_t resolved_window = 0U;
    double local_df = 0.0;
    double prior_df_pos = 0.0;
    double total_df = 0.0;
    double precision_floor_ratio_pos = 0.0;
    int status = 0;

    if (centered_matrix == NULL || mean_out == NULL || raw_variance_out == NULL ||
        prior_variance_out == NULL || moderated_variance_out == NULL ||
        standard_error_out == NULL || scores_out == NULL)
    {
        return -2;
    }
    if (sample_count == 0U || locus_count == 0U)
    {
        return -2;
    }

    prior_df_pos = fmax(prior_df, 0.0);
    precision_floor_ratio_pos = fmax(precision_floor_ratio, 0.0);
    resolved_window = rocco_resolve_spatial_window(locus_count, spatial_window);
    local_df = resolved_window > 0U ? fmax(4.0, (double)resolved_window - 3.0) : 1.0;
    total_df = local_df + prior_df_pos;
    if (degrees_of_freedom_out != NULL)
    {
        *degrees_of_freedom_out = total_df;
    }
    if (resolved_window_out != NULL)
    {
        *resolved_window_out = (int)resolved_window;
    }

    obs_variance = (double *)malloc(locus_count * sizeof(double));
    prior_track = (double *)malloc(locus_count * sizeof(double));
    weighted_sum = (double *)malloc(locus_count * sizeof(double));
    precision_sum = (double *)malloc(locus_count * sizeof(double));
    raw_precision_sum = (double *)malloc(locus_count * sizeof(double));
    prior_precision_sum = (double *)malloc(locus_count * sizeof(double));
    work_buffer = (double *)malloc(locus_count * sizeof(double));
    if (obs_variance == NULL || prior_track == NULL || weighted_sum == NULL ||
        precision_sum == NULL || raw_precision_sum == NULL ||
        prior_precision_sum == NULL || work_buffer == NULL)
    {
        free(obs_variance);
        free(prior_track);
        free(weighted_sum);
        free(precision_sum);
        free(raw_precision_sum);
        free(prior_precision_sum);
        free(work_buffer);
        return -1;
    }

    memset(weighted_sum, 0, locus_count * sizeof(double));
    memset(precision_sum, 0, locus_count * sizeof(double));
    memset(raw_precision_sum, 0, locus_count * sizeof(double));
    memset(prior_precision_sum, 0, locus_count * sizeof(double));

    /* one track at a time in scratch, stays O(n) */
    for (sample_idx = 0U; sample_idx < sample_count; ++sample_idx)
    {
        const double *row = centered_matrix + (sample_idx * locus_count);
        if (resolved_window == 0U || locus_count < 4U)
        {
            double scale_floor = 1.0e-6;
            for (locus_idx = 0U; locus_idx < locus_count; ++locus_idx)
            {
                work_buffer[locus_idx] = row[locus_idx];
            }
            scale_floor = rocco_robust_scale_f64(work_buffer, locus_count);
            scale_floor = fmax(scale_floor * scale_floor, 1.0e-8);
            for (locus_idx = 0U; locus_idx < locus_count; ++locus_idx)
            {
                obs_variance[locus_idx] = scale_floor;
                prior_track[locus_idx] = scale_floor;
            }
        }
        else
        {
            status = rocco_rolling_ar1_innovation_variance_f64(
                row,
                locus_count,
                resolved_window,
                obs_variance);
            if (status != 0)
            {
                free(obs_variance);
                free(prior_track);
                free(weighted_sum);
                free(precision_sum);
                free(raw_precision_sum);
                free(prior_precision_sum);
                free(work_buffer);
                return status == -1 ? -1 : -2;
            }
            for (locus_idx = 0U; locus_idx < locus_count; ++locus_idx)
            {
                obs_variance[locus_idx] = fmax(obs_variance[locus_idx], 1.0e-8);
            }
            status = rocco_fit_monotone_variance_trend_f64(
                row,
                obs_variance,
                locus_count,
                prior_track);
            if (status != 0)
            {
                free(obs_variance);
                free(prior_track);
                free(weighted_sum);
                free(precision_sum);
                free(raw_precision_sum);
                free(prior_precision_sum);
                free(work_buffer);
                return status == -1 ? -1 : -2;
            }
        }

        for (locus_idx = 0U; locus_idx < locus_count; ++locus_idx)
        {
            double posterior_variance = 0.0;
            double posterior_precision = 0.0;
            double variance_floor = 0.0;
            double obs_value = fmax(obs_variance[locus_idx], 1.0e-8);
            double prior_value = fmax(prior_track[locus_idx], 1.0e-8);

            posterior_variance =
                ((local_df * obs_value) + (prior_df_pos * prior_value)) /
                fmax(total_df, 1.0);
            variance_floor = precision_floor_ratio_pos * prior_value;
            if (posterior_variance < variance_floor)
            {
                posterior_variance = variance_floor;
            }
            posterior_variance = fmax(posterior_variance, 1.0e-8);
            posterior_precision = 1.0 / posterior_variance;
            raw_precision_sum[locus_idx] += 1.0 / obs_value;
            prior_precision_sum[locus_idx] += 1.0 / prior_value;
            precision_sum[locus_idx] += posterior_precision;
            weighted_sum[locus_idx] += posterior_precision * row[locus_idx];
        }
    }

    /* final m x n precision combine */
    for (locus_idx = 0U; locus_idx < locus_count; ++locus_idx)
    {
        double z_score = 0.0;
        double locus_precision = fmax(precision_sum[locus_idx], 1.0e-8);
        mean_out[locus_idx] = weighted_sum[locus_idx] / locus_precision;
        raw_variance_out[locus_idx] =
            (double)sample_count / fmax(raw_precision_sum[locus_idx], 1.0e-8);
        prior_variance_out[locus_idx] =
            (double)sample_count / fmax(prior_precision_sum[locus_idx], 1.0e-8);
        moderated_variance_out[locus_idx] = (double)sample_count / locus_precision;
        standard_error_out[locus_idx] = sqrt(1.0 / locus_precision);
        z_score = mean_out[locus_idx] / fmax(standard_error_out[locus_idx], 1.0e-8);
        if (use_min_effect != 0)
        {
            scores_out[locus_idx] =
                (mean_out[locus_idx] - fmax(min_effect, 0.0)) /
                fmax(standard_error_out[locus_idx], 1.0e-8);
        }
        else
        {
            scores_out[locus_idx] = z_score - lower_bound_z;
        }
    }

    free(obs_variance);
    free(prior_track);
    free(weighted_sum);
    free(precision_sum);
    free(raw_precision_sum);
    free(prior_precision_sum);
    free(work_buffer);
    return 0;
}
