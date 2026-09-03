import importlib as _importlib

from ._version import __version__
from .constants import GENOME_DICT
from .dp import (
    build_switch_costs,
    objective_value,
    solve_chrom_exact,
    solve_penalized_chain,
)
from .dependence import (
    choose_dependence_span,
    estimate_dependence_radius_for_window,
)
from .inference import (
    candidate_mask_from_wls,
    estimateStationaryBootstrapBudget,
    estimate_empirical_bayes_budgets,
    fit_beta_prior_mle,
    score_loci_wls,
)
from .readtracks import (
    check_type_bam_files,
    generate_chrom_matrix,
    get_bam_chrom_reads,
    get_bigwig_chrom_scores,
    get_chroms_and_sizes,
    get_track_type,
)
from .scores import (
    EmpiricalNull,
    get_ecdf,
    get_read_length,
    multi_ecdf,
    raw_count_matrix,
    score_peaks,
)

_ROCCO_EXPORTS = {
    "chrom_solution_to_bed",
    "combine_chrom_results",
    "cscores_quantiles",
    "json_config",
    "main",
    "resolve_config",
    "score_central_tendency_chrom",
}
_SUBMODULES = {
    "constants",
    "dependence",
    "dp",
    "inference",
    "readtracks",
    "rocco",
    "scores",
}

__all__ = [
    "EmpiricalNull",
    "GENOME_DICT",
    "__version__",
    "build_switch_costs",
    "candidate_mask_from_wls",
    "check_type_bam_files",
    "chrom_solution_to_bed",
    "choose_dependence_span",
    "combine_chrom_results",
    "cscores_quantiles",
    "estimateStationaryBootstrapBudget",
    "estimate_dependence_radius_for_window",
    "estimate_empirical_bayes_budgets",
    "fit_beta_prior_mle",
    "generate_chrom_matrix",
    "get_bam_chrom_reads",
    "get_bigwig_chrom_scores",
    "get_chroms_and_sizes",
    "get_ecdf",
    "get_read_length",
    "get_track_type",
    "json_config",
    "main",
    "multi_ecdf",
    "objective_value",
    "raw_count_matrix",
    "resolve_config",
    "score_central_tendency_chrom",
    "score_loci_wls",
    "score_peaks",
    "solve_chrom_exact",
    "solve_penalized_chain",
]


def __getattr__(name):
    if name in _ROCCO_EXPORTS:
        rocco_module = _importlib.import_module(".rocco", __name__)
        value = getattr(rocco_module, name)
        globals()[name] = value
        for module_name in (
            "constants",
            "dependence",
            "dp",
            "inference",
            "readtracks",
            "rocco",
            "scores",
        ):
            globals().pop(module_name, None)
        return value
    if name in _SUBMODULES:
        return _importlib.import_module(f".{name}", __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(__all__)


for _module_name in (
    "constants",
    "dependence",
    "dp",
    "inference",
    "readtracks",
    "scores",
):
    globals().pop(_module_name, None)
del _module_name
