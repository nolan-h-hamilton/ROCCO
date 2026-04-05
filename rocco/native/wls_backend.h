#ifndef ROCCO_WLS_BACKEND_H
#define ROCCO_WLS_BACKEND_H

#include <stddef.h>

#ifdef __cplusplus
extern "C"
{
#endif

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
        int *resolved_window_out);

#ifdef __cplusplus
}
#endif

#endif
