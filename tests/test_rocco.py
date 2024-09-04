import os
import numpy as np
import pybedtools as pbt
import pytest
    
from rocco import *

@pytest.fixture
def test_setup():
    chromosomes = ['chr19', 'chr21', 'chrX']
    chrom_budget_dict = {'chr19': 0.045, 'chr21': 0.02, 'chrX': 0.015}
    chrom_ref_results = {'chr19': 'ref_chr19.bed', 'chr21': 'ref_chr21.bed', 'chrX': 'ref_chrX.bed'}
    combined_ref_file = 'combined_ref.bed'
    matrices = np.load('test_data.npz')
    intervals = np.load('test_intervals.npz')
    ID_ = 'test'
    C1 = 1.0
    C2 = -1.0
    C3 = 1.0

    min_gap = 0.95
    min_jaccard = 0.99
    
    return {
        "chromosomes": chromosomes,
        "chrom_budget_dict": chrom_budget_dict,
        "chrom_ref_results": chrom_ref_results,
        "combined_ref_file": combined_ref_file,
        "matrices": matrices,
        "intervals": intervals,
        "ID_": ID_,
        "C1": C1,
        "C2": C2,
        "C3": C3,
        "min_gap": min_gap,
        "min_jaccard": min_jaccard
    }


@pytest.mark.consistency
def test_consistency_multiple_chromosomes(test_setup, seed=42):
    np.random.seed(seed)  # Set the random seed for consistency
    
    # Unpack the setup data
    chromosomes = test_setup["chromosomes"]
    chrom_budget_dict = test_setup["chrom_budget_dict"]
    chrom_ref_results = test_setup["chrom_ref_results"]
    matrices = test_setup["matrices"]
    intervals = test_setup["intervals"]
    ID_ = test_setup["ID_"]
    C1 = test_setup["C1"]
    C2 = test_setup["C2"]
    C3 = test_setup["C3"]
    min_gap = test_setup["min_gap"]
    min_jaccard = test_setup["min_jaccard"]

    for chrom_ in chromosomes:
        # use default scoring scheme and parameters for this test
        ct_scores = rocco.score_central_tendency_chrom(matrices[chrom_], method='quantile', quantile=0.50, power=1.0)
        disp_scores = rocco.score_dispersion_chrom(matrices[chrom_], method='mad', power=1.0)
        boundary_scores = rocco.score_boundary_chrom(ct_scores, denom=1.0, power=1.0)
        chrom_scores = C1*ct_scores + C2*disp_scores + C3*boundary_scores

        # solve/refine with default parameters
        chrom_lp_sol, lp_bound = rocco.solve_relaxation_chrom_pdlp(chrom_scores, budget=chrom_budget_dict[chrom_], gamma=1.0, verbose=True)
        rround_sol, rround_objval = rocco.get_rround_sol(chrom_lp_sol, chrom_scores, chrom_budget_dict[chrom_], gamma=1.0)
        chrom_outfile = chrom_solution_to_bed(chrom_, intervals[chrom_], rround_sol, ID=ID_)
        assert os.path.exists(chrom_outfile), f'{chrom_}: Solution file not found'
        chrom_pbt = pbt.BedTool(chrom_outfile)
        chrom_ref_pbt = pbt.BedTool(chrom_ref_results[chrom_])
        chrom_jaccard = round(chrom_pbt.jaccard(chrom_ref_pbt)['jaccard'],5)
        gap_ = round(rround_objval / lp_bound,5)
        assert gap_ >= min_gap, f'{chrom_}: Performance gap with LP-ideal is below threshold {gap_}: < {min_gap}'
        assert chrom_jaccard > min_jaccard, f'{chrom_}: Jaccard index for results with reference is below threshold {chrom_jaccard} < {min_jaccard}'


@pytest.mark.correctness
def test_combine_chrom_results_no_names(test_setup):
    # combine the chromosome-specific reference bed files in the repo and combine to create a new file
    combined_outfile = combine_chrom_results([str(x) for x in test_setup["chrom_ref_results"].values()], output_file='test_combined.bed')
    assert os.path.exists(combined_outfile), f'Combined solution file {combined_outfile} not found'
    # load the combined file created by combine_chrom_results into a BedTool object
    combined_pbt = pbt.BedTool(combined_outfile)

    # load the already-combined reference combined file that is in the repo to a BedTool object
    combined_ref_pbt = pbt.BedTool(test_setup['combined_ref_file'])
    
    # check the jaccard index between the combined solution and the reference
    combined_jaccard = round(combined_pbt.jaccard(combined_ref_pbt)['jaccard'],5)

    assert combined_jaccard > test_setup["min_jaccard"], f'Jaccard index for combined results with reference is below threshold {combined_jaccard} < {test_setup["min_jaccard"]}'
    
    try:
        os.remove(combined_outfile)
    except:
        pass


@pytest.mark.correctness
def test_combine_chrom_results_add_names(test_setup):
    # combine the chromosome-specific reference bed files in the repo and combine to create a new file
    combined_outfile = combine_chrom_results([str(x) for x in test_setup["chrom_ref_results"].values()], output_file='test_combined.bed', name_features=True)
    assert os.path.exists(combined_outfile), f'Combined solution file {combined_outfile} not found'
    # load the combined file created by combine_chrom_results into a BedTool object
    combined_pbt = pbt.BedTool(combined_outfile)

    # load the already-combined reference combined file that is in the repo to a BedTool object
    combined_ref_pbt = pbt.BedTool(test_setup['combined_ref_file'])
    
    # check the jaccard index between the combined solution and the reference
    combined_jaccard = round(combined_pbt.jaccard(combined_ref_pbt)['jaccard'],5)
    assert combined_jaccard > test_setup["min_jaccard"], f'Jaccard index for combined results with reference is below threshold {combined_jaccard} < {test_setup["min_jaccard"]}'
    
    # now check the naming
    assert combined_pbt.field_count() > 3, f'Combined solution file does not have the expected number of fields: {combined_pbt.field_count()}'
    np.random.seed(42)
    for i in np.random.choice(range(len(combined_pbt)-1), size=50):
        idx_ = int(i)
        assert combined_pbt[idx_].name.split('_')[0] == combined_pbt[idx_].chrom, f'Incorrect chromosome name assigned to feature {combined_pbt[idx_].name}'
        assert combined_pbt[idx_].name.split('_')[1] == str(combined_pbt[idx_].start), f'Incorrect start position assigned to feature {combined_pbt[idx_].name}'
        assert combined_pbt[idx_].name.split('_')[2] == str(combined_pbt[idx_].end), f'Incorrect end position assigned to feature {combined_pbt[idx_].name}'
    try:
        os.remove(combined_outfile)
    except:
        pass

@pytest.mark.correctness
def test_minmax_scale():

    # case 1
    x = [1, 2, 3, 4, 5]
    minval_ = 0
    maxval_ = 10
    x_minmax = minmax_scale(x, min_val=minval_, max_val=maxval_)
    assert min(x_minmax) == 0, f'Min-max scaling failed : {x, x_minmax, minval_, maxval_}'
    assert max(x_minmax) == 10, f'Min-max scaling failed : {x, x_minmax, minval_, maxval_}'
    
    # case 2
    x = [-1, 0, 1, 2, 3]
    minval_ = 0
    maxval_ = 4
    x_minmax = minmax_scale(x, min_val=minval_, max_val=maxval_)
    assert min(x_minmax) == 0, f'Min-max scaling failed : {x, x_minmax, minval_, maxval_}'
    assert max(x_minmax) == 4, f'Min-max scaling failed : {x, x_minmax, minval_, maxval_}'
