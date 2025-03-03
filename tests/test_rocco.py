import os
import numpy as np
import pybedtools as pbt
import pytest
import subprocess
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
    r"""Test the consistency of the solution across multiple chromosomes
    
    Track performance against original reference 

    """
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
    r"""Test :func:`rocco.combine_chrom_results` with no names"""
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
def test_score_central_tendency_chrom():
    r"""Test :func:`rocco.score_central_tendency_chrom`"""
    # Case I
    X = np.array([[1, 2, 3, 4, 5], [2, 3, 4, 5, 6], [3, 4, 5, 6, 10]])
    method = 'quantile'
    scores = score_central_tendency_chrom(X, method=method)
    for i,x in enumerate([2,3,4,5,6]):
        assert scores[i] == x, f'Central tendency scoring failed : {X, scores, method}'


    # Case II
    method = 'mean'
    scores = score_central_tendency_chrom(X, method=method)
    for i,x in enumerate([2,3,4,5,7]):
        assert scores[i] == x, f'Central tendency scoring failed : {X, scores, method}'
    
    # Case III
    X = np.array([[0, 0, 0, 0, 0], [2, 3, 4, 5, 6],
                  [1, 2, 3, 4, 5], [2, 3, 4, 5, 6], 
                  [1, 2, 3, 4, 5], [2, 3, 4, 5, 6],
                  [1, 2, 3, 4, 5], [5, 3, 4, 5, 6],
                  [1, 2, 3, 4, 5], [1000, 1000, 1000, 1000, 1000]])
    method = 'tmean'
    scores = score_central_tendency_chrom(X, method=method, tprop=.11)
    # output should be a list of 5 values -- the mean of the center 8 values in each column
    for i,x in enumerate([1.875, 2.5, 3.5, 4.5, 5.5]):
        assert scores[i] == x, f'Central tendency scoring failed : {X, scores, method}'


@pytest.mark.correctness
def test_score_dispersion_chrom():
    r"""Test :func:`rocco.score_dispersion_chrom`"""
    
    # Case I
    X = np.zeros(shape=(11,3))
    X[:,0] = [0, 1, 1, 1, 1, 5, 25, 25, 25, 25, 100]
    X[:,1] = np.ones(11)
    X[:,2] = [1,1,1,0,1,0,1,0,1,0,1]


    method = 'mad'
    scores = score_dispersion_chrom(X, method=method)
    for i,x in enumerate([5,0,0]):
        assert scores[i] == x, f'Dispersion method failed: {X, scores, method}'

    # Case II
    method = 'std'
    scores = score_dispersion_chrom(X, method=method)
    for i,x in enumerate([27.89265, 0, 0.481045]):
        assert np.isclose(x,scores[i], rtol=1e-4, atol=1e-4), f'Dispersion method failed : {X, scores, method}'

    # Case III
    method = 'iqr'
    scores = score_dispersion_chrom(X, method=method)
    for i,x in enumerate([24, 0, 1]):
        assert scores[i] == x, f'Dispersion method failed : {X, scores, method}'


@pytest.mark.correctness
def test_score_boundary_chrom():
    r"""Test :func:`rocco.score_boundary_chrom` with default parameters"""
    vec = np.array([-10,-5,0,3,10])
    denom = 1
    scores = score_boundary_chrom(vec, denom=denom)
    for i,x in enumerate([5/11.0, 5/6.0, 5.0, 7/4.0, 7/11.0]):
        assert scores[i] == x, f'Boundary scoring failed : {vec, denom}'


@pytest.mark.correctness
def test_parsig_default_parameters():
    r"""Test :func:`rocco.parsig`"""
    scores = np.array([-50, 1, 2, 3, 4, 5, 50]) 
    B_ = 0.80
    M_ = 10
    R_ = 5
    # remove_min is True by default
    transformed_scores = parsig(scores, gamma=1.0, parsig_B=B_, parsig_M=M_, parsig_R=R_)
    assert np.min(transformed_scores) >= 0, f'Negative values found in transformed scores: {transformed_scores}'
    assert np.max(transformed_scores) <= 10, f'Maximum value of transformed scores is unexpected: {np.max(transformed_scores)}'


@pytest.mark.correctness
def test_get_floor_eps_sol_basic():
    r"""Test :func:`rocco.get_floor_eps_sol` with basic input"""
    chrom_lp_sol = np.array([1,0,.9,.9,.96,.97,0,0,0,.95])
    budget = 0.5
    int_tol = 1e-6
    eps_mult = 1.01
    result = get_floor_eps_sol(chrom_lp_sol, budget, int_tol, eps_mult)
    assert np.all(result >= 0) and np.all(result <= 1), f"feasible region for each dvar should be [0,1]: {result}"
    assert np.sum(result) <= np.floor(len(chrom_lp_sol) * budget), f"Solution should respect the budget: {result}, {np.sum(result)}"
    assert all([np.isclose(abs(int(x) - x),0) for x in result]), f"All solutions should be binary (0 or 1), {result}"
    assert np.allclose(np.array([1, 0, 0, 0, 1, 1, 0, 0, 0, 1,]), result), f"Expected solution: {result}"


@pytest.mark.correctness
def test_get_floor_eps_sol_allint():
    r"""Test :func:`rocco.get_floor_eps_sol` with all integer values"""
    chrom_lp_sol = np.array([1,1,1,1,1,0,0,0,0,0])
    budget = 0.5
    int_tol = 1e-6
    eps_mult = 1.01
    result = get_floor_eps_sol(chrom_lp_sol, budget, int_tol, eps_mult)
    assert np.all(result >= 0) and np.all(result <= 1), f"feasible region for each dvar should be [0,1]: {result}"
    assert np.sum(result) <= np.floor(len(chrom_lp_sol) * budget), f"Solution should respect the budget: {result}, {np.sum(result)}"
    assert all([np.isclose(x,1) for x in result[0:5]]), f"First five dvars should be 1, {result}"


@pytest.mark.correctness
def test_get_floor_eps_sol_fullbudget():
    r"""Test :func:`rocco.get_floor_eps_sol` with a budget of 1.0"""
    chrom_lp_sol = np.ones(10)*.50
    budget = 1.0
    int_tol = 1e-6
    eps_mult = 1.01
    result = get_floor_eps_sol(chrom_lp_sol, budget, int_tol, eps_mult)
    assert np.all(result >= 0) and np.all(result <= 1), f"feasible region for each dvar should be [0,1]: {result}"
    assert np.sum(result) <= np.floor(len(chrom_lp_sol) * budget), f"Solution should respect the budget: {result}, {np.sum(result)}"
    assert all([np.isclose(x,1) for x in result]), f"Solution should be ones vec, {result}"


@pytest.mark.correctness
def test_no_input_no_args():
    result = subprocess.run(['rocco'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert result.returncode == 0, f'Expected return code 0, got {result.returncode}'
    assert 'usage:' in result.stdout.decode(), f'Expected help message, got {result.stdout.decode()}'


@pytest.mark.correctness
def test_no_input_listed():
    result = subprocess.run(['rocco', '--input_files'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert result.returncode != 0, f'Expected non-zero return code, got {result.returncode}'
    assert 'usage:' in result.stderr.decode(), f'Expected help message, got {result.stderr.decode()}'


@pytest.mark.correctness
def test_unrecognized_arg():
    result = subprocess.run(['rocco', '--unrecognized_arg'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert result.returncode != 0, f'Expected non-zero return code, got {result.returncode}'
    assert 'unrecognized' in result.stderr.decode(), f'Expected unrecognized argument message, got {result.stderr.decode()}'


@pytest.mark.correctness
def test_bedgraph_input():
    result = subprocess.run(['rocco', '--input_files', 'test.bedgraph'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Should raise a value error with message containing "Please convert"
    assert result.returncode != 0, f'Expected non-zero return code, got {result.returncode}'
    assert 'convert' in result.stderr.decode().lower(), f'Expected error message requesting conversion to bigwig or use of BAM input, got {result.stderr.decode()}'


@pytest.mark.correctness
def test_resolve_transformations():
    r"""Test :func:`rocco.resolve_transformations`"""
    
    # Case 1: Conflicting transformations explicitly specified
    args = {
        'transform_savgol': True,
        'transform_savgol_window_bp': 100,
        'transform_savgol_order': 3,
        'transform_medfilt': True,
        'transform_medfilt_kernel': 5
    }
    resolved_args = resolve_transformations(args)
    assert resolved_args['transform_savgol'] == False, f'Expected False, got {resolved_args["transform_savgol"]}'
    assert resolved_args['transform_savgol_window_bp'] is None, f'Expected None, got {resolved_args["transform_savgol_window_bp"]}'
    assert resolved_args['transform_savgol_order'] is None, f'Expected None, got {resolved_args["transform_savgol_order"]}'
    assert resolved_args['transform_medfilt'] == True, f'Expected True, got {resolved_args["transform_medfilt"]}'
    assert resolved_args['transform_medfilt_kernel'] == 5, f'Expected 5, got {resolved_args["transform_medfilt_kernel"]}'
    
    # Case II: Conflicting transformations implicated only by their respective flags
    args2 = {
        'transform_savgol': False,
        'transform_savgol_window_bp': 100,
        'transform_savgol_order': 3,
        'transform_medfilt': False,
        'transform_medfilt_kernel': 5
    }
    resolved_args2 = resolve_transformations(args2)
    assert resolved_args2['transform_savgol'] == False, f'Expected False, got {resolved_args2["transform_savgol"]}'
    assert resolved_args2['transform_savgol_window_bp'] is None, f'Expected None, got {resolved_args2["transform_savgol_window_bp"]}'
    assert resolved_args2['transform_savgol_order'] is None, f'Expected None, got {resolved_args2["transform_savgol_order"]}'
    assert resolved_args2['transform_medfilt'] == True, f'Expected True, got {resolved_args2["transform_medfilt"]}'
    assert resolved_args2['transform_medfilt_kernel'] == 5, f'Expected 5, got {resolved_args2["transform_medfilt_kernel"]}'
