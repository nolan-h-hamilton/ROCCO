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
        double *sample_weights_out,
        double *mean_out,
        double *standard_error_out,
        double *z_scores_out,
        double *scores_out);

#ifdef __cplusplus
}
#endif

#endif
