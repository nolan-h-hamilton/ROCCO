#ifndef ROCCO_BASELINE_BACKEND_H
#define ROCCO_BASELINE_BACKEND_H

#include <stddef.h>

#ifdef __cplusplus
extern "C"
{
#endif

    int rocco_crossfit_whittaker_baseline_f64(
        const double *y_values,
        size_t value_count,
        double penalty_lambda,
        double *baseline_out);

    int rocco_crossfit_whittaker_baseline_matrix_f64(
        const double *matrix_values,
        size_t row_count,
        size_t column_count,
        double penalty_lambda,
        double *baseline_out);

#ifdef __cplusplus
}
#endif

#endif
