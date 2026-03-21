#include "baseline_backend.h"

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/*
 * Native cross-fit Whittaker baseline backend for ROCCO.
 *
 * This implementation is adapted from the Consenrich baseline routines:
 *   https://github.com/nolan-h-hamilton/Consenrich
 *   https://github.com/nolan-h-hamilton/Consenrich/blob/main/src/consenrich/cconsenrich.pyx
 */

typedef struct rocco_baseline_state
{
    double *a0;
    double *a1;
    double *a2;
    double *d;
    double *l1;
    double *l2;
    double *forward_rhs;
    double *backward_rhs;
    double *other_baseline;
    double *buffer;
} rocco_baseline_state;

static int rocco_init_baseline_state(
    size_t value_count,
    rocco_baseline_state *state)
{
    size_t off1_count = value_count > 1 ? (value_count - 1) : 1;
    size_t off2_count = value_count > 2 ? (value_count - 2) : 1;
    size_t total_count =
        (5 * value_count) + (2 * off1_count) + (2 * off2_count);
    double *buffer = NULL;

    if (state == NULL)
    {
        return -1;
    }

    memset(state, 0, sizeof(*state));
    buffer = (double *)malloc(total_count * sizeof(double));
    if (buffer == NULL)
    {
        return -1;
    }

    state->buffer = buffer;
    state->a0 = buffer;
    state->a1 = state->a0 + value_count;
    state->a2 = state->a1 + off1_count;
    state->d = state->a2 + off2_count;
    state->l1 = state->d + value_count;
    state->l2 = state->l1 + off1_count;
    state->forward_rhs = state->l2 + off2_count;
    state->backward_rhs = state->forward_rhs + value_count;
    state->other_baseline = state->backward_rhs + value_count;
    return 0;
}

static void rocco_free_baseline_state(
    rocco_baseline_state *state)
{
    if (state == NULL)
    {
        return;
    }
    if (state->buffer != NULL)
    {
        free(state->buffer);
    }
    memset(state, 0, sizeof(*state));
}

static void rocco_solve_sod_ldl_f64(
    const double *a0,
    const double *a1,
    const double *a2,
    const double *rhs,
    double *solution,
    double *d,
    double *l1,
    double *l2,
    double *forward_rhs,
    double *backward_rhs,
    size_t value_count)
{
    ptrdiff_t idx = 0;
    double t1 = 0.0;
    double t2 = 0.0;

    if (value_count == 0)
    {
        return;
    }
    if (value_count == 1)
    {
        solution[0] = rhs[0] / a0[0];
        return;
    }

    d[0] = a0[0];
    l1[0] = a1[0] / d[0];
    if (value_count > 2)
    {
        l2[0] = a2[0] / d[0];
    }

    d[1] = a0[1] - ((l1[0] * l1[0]) * d[0]);
    if (value_count > 2)
    {
        t1 = ((l2[0] * d[0]) * l1[0]);
        l1[1] = (a1[1] - t1) / d[1];
    }
    if (value_count > 3)
    {
        l2[1] = a2[1] / d[1];
    }

    for (idx = 2; idx < (ptrdiff_t)value_count; ++idx)
    {
        t1 = ((l1[idx - 1] * l1[idx - 1]) * d[idx - 1]);
        t2 = ((l2[idx - 2] * l2[idx - 2]) * d[idx - 2]);
        d[idx] = a0[idx] - t1 - t2;

        if (idx <= ((ptrdiff_t)value_count - 2))
        {
            t1 = ((l2[idx - 1] * d[idx - 1]) * l1[idx - 1]);
            l1[idx] = (a1[idx] - t1) / d[idx];
        }

        if (idx <= ((ptrdiff_t)value_count - 3))
        {
            l2[idx] = a2[idx] / d[idx];
        }
    }

    /* Forward solve :math:`Ly = b`. */
    forward_rhs[0] = rhs[0];
    forward_rhs[1] = rhs[1] - (l1[0] * forward_rhs[0]);
    for (idx = 2; idx < (ptrdiff_t)value_count; ++idx)
    {
        t1 = l1[idx - 1] * forward_rhs[idx - 1];
        t2 = l2[idx - 2] * forward_rhs[idx - 2];
        forward_rhs[idx] = rhs[idx] - t1 - t2;
    }

    /* Diagonal solve :math:`Dz = y`. */
    for (idx = 0; idx < (ptrdiff_t)value_count; ++idx)
    {
        backward_rhs[idx] = forward_rhs[idx] / d[idx];
    }

    /* Backward solve :math:`L^\top x = z`. */
    solution[value_count - 1] = backward_rhs[value_count - 1];
    solution[value_count - 2] =
        backward_rhs[value_count - 2] -
        (l1[value_count - 2] * solution[value_count - 1]);
    if (value_count == 2)
    {
        return;
    }
    for (idx = (ptrdiff_t)value_count - 3; idx >= 0; --idx)
    {
        t1 = l1[idx] * solution[idx + 1];
        t2 = l2[idx] * solution[idx + 2];
        solution[idx] = backward_rhs[idx] - t1 - t2;
    }
}

static void rocco_loc_baseline_masked_f64(
    const double *y_values,
    size_t value_count,
    int parity,
    double penalty_lambda,
    double *baseline_out,
    rocco_baseline_state *state)
{
    size_t idx = 0;
    double weight = 0.0;

    if (value_count < 3)
    {
        if (value_count > 0)
        {
            memcpy(
                baseline_out,
                y_values,
                value_count * sizeof(double));
        }
        return;
    }

    /*
     * Solve :math:`(W + \lambda D^\top D)b = Wy` with a parity mask in
     * :math:`W` so one parity is fit from the other.
     */
    state->a0[0] = ((parity == 0) ? 1.0 : 0.0) + penalty_lambda;
    state->a0[1] = ((parity == 1) ? 1.0 : 0.0) + (5.0 * penalty_lambda);
    baseline_out[0] = (parity == 0) ? y_values[0] : 0.0;
    baseline_out[1] = (parity == 1) ? y_values[1] : 0.0;

    for (idx = 2; idx < value_count - 2; ++idx)
    {
        weight = ((idx & 1U) == (size_t)parity) ? 1.0 : 0.0;
        state->a0[idx] = weight + (6.0 * penalty_lambda);
        baseline_out[idx] = weight * y_values[idx];
    }

    state->a0[value_count - 2] =
        ((((value_count - 2) & 1U) == (size_t)parity) ? 1.0 : 0.0) +
        (5.0 * penalty_lambda);
    state->a0[value_count - 1] =
        ((((value_count - 1) & 1U) == (size_t)parity) ? 1.0 : 0.0) +
        penalty_lambda;
    baseline_out[value_count - 2] =
        ((((value_count - 2) & 1U) == (size_t)parity) ? y_values[value_count - 2] : 0.0);
    baseline_out[value_count - 1] =
        ((((value_count - 1) & 1U) == (size_t)parity) ? y_values[value_count - 1] : 0.0);

    /* These are the three bands of :math:`\lambda D^\top D`. */
    state->a1[0] = -2.0 * penalty_lambda;
    for (idx = 1; idx < value_count - 2; ++idx)
    {
        state->a1[idx] = -4.0 * penalty_lambda;
    }
    state->a1[value_count - 2] = -2.0 * penalty_lambda;

    for (idx = 0; idx < value_count - 2; ++idx)
    {
        state->a2[idx] = penalty_lambda;
    }

    rocco_solve_sod_ldl_f64(
        state->a0,
        state->a1,
        state->a2,
        baseline_out,
        baseline_out,
        state->d,
        state->l1,
        state->l2,
        state->forward_rhs,
        state->backward_rhs,
        value_count);
}

int rocco_crossfit_whittaker_baseline_f64(
    const double *y_values,
    size_t value_count,
    double penalty_lambda,
    double *baseline_out)
{
    rocco_baseline_state state;
    size_t idx = 0;

    if (y_values == NULL || baseline_out == NULL)
    {
        return -1;
    }

    if (value_count < 25)
    {
        for (idx = 0; idx < value_count; ++idx)
        {
            baseline_out[idx] = 0.0;
        }
        return 0;
    }

    if (rocco_init_baseline_state(value_count, &state) != 0)
    {
        return -1;
    }

    /* Cross-fit average :math:`\hat b = (\hat b_{\mathrm{even}} + \hat b_{\mathrm{odd}})/2`. */
    rocco_loc_baseline_masked_f64(
        y_values,
        value_count,
        0,
        penalty_lambda,
        baseline_out,
        &state);
    rocco_loc_baseline_masked_f64(
        y_values,
        value_count,
        1,
        penalty_lambda,
        state.other_baseline,
        &state);

    for (idx = 0; idx < value_count; ++idx)
    {
        baseline_out[idx] = 0.5 * (baseline_out[idx] + state.other_baseline[idx]);
    }

    rocco_free_baseline_state(&state);
    return 0;
}

int rocco_crossfit_whittaker_baseline_matrix_f64(
    const double *matrix_values,
    size_t row_count,
    size_t column_count,
    double penalty_lambda,
    double *baseline_out)
{
    size_t row_idx = 0;

    if (matrix_values == NULL || baseline_out == NULL)
    {
        return -1;
    }

    for (row_idx = 0; row_idx < row_count; ++row_idx)
    {
        const double *row_values = matrix_values + (row_idx * column_count);
        double *row_out = baseline_out + (row_idx * column_count);
        if (rocco_crossfit_whittaker_baseline_f64(
                row_values,
                column_count,
                penalty_lambda,
                row_out) != 0)
        {
            return -1;
        }
    }

    return 0;
}
