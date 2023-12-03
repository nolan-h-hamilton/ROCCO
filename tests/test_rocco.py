import os
import subprocess
import shutil
import pytest

solver = 'CLARABEL'
peakfile_minsize = 100
jaccard_min = .95
rr_iter = '100'
chr19_start = '20000000'
chr19_end = '50000000'
chr19_budget= '0.035'
chr20_start = '20000000'
chr20_end = '50000000'
chr20_budget = '0.035'

group_A_bedfile = 'group_A.bed'
group_B_bedfile = 'group_B.bed'

tracks_chrom19_wig = os.path.relpath('data/tracks_chr19/chr19_test.bam.bw.wig.orig')
tracks_chrom19_wig_fs = os.path.relpath('data/tracks_chr19_fixedStep/chr19_test_bam.bw.wig.fs.orig')
tracks_chrom20_wig = os.path.relpath('data/tracks_chr20/chr20_test.bam.bw.wig.orig')

chr19_outfile = f'ROCCO_out_chr19_{chr19_budget}_1.0_0.0_1.0_1.0_1.0.bed'
chr19_ref_file = os.path.relpath('data/ref/ref_chr19.bed')

chr20_outfile = f'ROCCO_out_chr20_{chr20_budget}_1.0_0.0_1.0_1.0_1.0.bed'
chr20_ref_file = os.path.relpath('data/ref/ref_chr20.bed')

gwide_outfile = 'gwide_out.bed'
gwide_outfile_fs = 'gwide_out_fs.bed'
gwide_outdir = os.path.relpath('output/combined')
gwide_outdir_fs = os.path.relpath('output/combined_fs')
gwide_ref_file = os.path.relpath('data/ref/ref_gwide.bed')

metadata_file = 'test_metadata.tsv'
identifiers_file = 'test_identifiers.txt'
identifiers_outdir = 'test_identifiers_outdir'
identifiers_outdir_fs = 'test_identifiers_outdir_fs'

@pytest.mark.regular
def test_sim(var=1, samples=10):
    """
    Add noise to existing wiggle file to generate varying tracks for testing purposes.
    `samples` signal tracks are generated.
    """
    proc1 = subprocess.run(['python', 'sim.py', tracks_chrom19_wig, str(var), str(samples)],
                           capture_output=True, encoding='utf-8')
    proc2 = subprocess.run(['python', 'sim.py', tracks_chrom20_wig, str(var), str(samples)],
                           capture_output=True, encoding='utf-8')
    assert proc1.returncode == 0, f'sim.py failed for {tracks_chrom19_wig}'
    assert proc2.returncode == 0, f'sim.py failed for {tracks_chrom20_wig}'

@pytest.mark.regular
def test_rocco_exists():
    """
    Ensure command line utility 'rocco' is available
    """
    assert shutil.which('rocco') is not None

@pytest.mark.regular
def test_chrom_help():
    """
    Ensure 'rocco chrom -h' returns help message
    """
    # check if return code is 0
    proc = subprocess.run(['rocco', 'chrom', '-h'],capture_output=True, encoding='utf-8')
    assert proc.returncode == 0
    assert 'usage: rocco chrom' in proc.stdout

@pytest.mark.regular
def test_gwide_help():
    """
    Ensure 'rocco gwide -h' returns help message
    """
    # check if return code is 0
    proc = subprocess.run(['rocco', 'gwide', '-h'],capture_output=True, encoding='utf-8')
    assert proc.returncode == 0
    assert 'usage: rocco gwide' in proc.stdout

@pytest.mark.regular
def test_budgets_help():
    """
    Ensure 'rocco budgets -h' returns help message
    """
    # check if return code is 0
    proc = subprocess.run(['rocco', 'budgets', '-h'],capture_output=True, encoding='utf-8')
    assert proc.returncode == 0
    assert 'usage: rocco budgets' in proc.stdout

@pytest.mark.regular
def test_prep_help():
    """
    Ensure 'rocco prep -h' returns help message
    """
    # check if return code is 0
    proc = subprocess.run(['rocco', 'prep', '-h'],capture_output=True, encoding='utf-8')
    assert proc.returncode == 0
    assert 'usage: rocco prep' in proc.stdout

@pytest.mark.regular
def test_run_chrom_chr19():
    """
    Ensure 'rocco chrom' runs successfully for chr19
    """
    # check if command returns exit code '0'
    assert subprocess.run(['rocco', 'chrom','--chrom', 'chr19', '--wig_path',  os.path.dirname(tracks_chrom19_wig),
                           '--budget', chr19_budget, '-N', rr_iter, '--start', chr19_start, '--end', chr19_end,
                           '--solver', solver],
                          capture_output=True, encoding='utf-8').returncode == 0, 'rocco chrom returncode != 0'
    # check if output file was created properly
    assert chr19_outfile in os.listdir()
    # ensure output file is nonempty (at least 'x' bytes)
    assert os.stat(chr19_outfile).st_size > peakfile_minsize

@pytest.mark.regular
def test_comp_ref_chr19():
    """
    Ensure the generated BED files are reasonably consistent with constant references.
    Wont match perfectly due the random simulation, but should exhibit broad similarity
    """
    cmd = ["bedtools", "jaccard", "-a", chr19_outfile , "-b", chr19_ref_file ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    data = proc.stdout.splitlines()[1]
    jaccard_chr19 = float(data.split("\t")[2])
    assert jaccard_chr19 >= jaccard_min

@pytest.mark.regular
def test_run_chrom_chr19_identifiers():
    """
    Ensure 'rocco chrom --identifiers' runs successfully and uses correct subset
    of samples in `test_identifiers.txt`
    """
    # check if command returns exit code '0'
    assert subprocess.run(['rocco', 'chrom','--chrom', 'chr19', '--wig_path',  os.path.dirname(tracks_chrom19_wig),
                           '--budget', chr19_budget, '-N', rr_iter, '--start', '10000000', '--end', '20000000',
                           '--solver', solver, '--identifiers', identifiers_file, '--outdir', identifiers_outdir],
                          capture_output=True, encoding='utf-8').returncode == 0, 'rocco chrom returncode != 0'
    # check if output file was created in the output directory `identifiers_outdir`
    assert chr19_outfile in os.listdir(identifiers_outdir)
    # ensure output file is nonempty (at least 'x' bytes)
    assert os.stat(os.path.join(identifiers_outdir,chr19_outfile)).st_size > peakfile_minsize

    identifier_tags = []
    id_wig_match_dict = {}
    for line in open(identifiers_file, mode='r', encoding='utf-8'):
        identifier_tags.append(line.strip())
    for line in open(os.path.join(identifiers_outdir, 'chr19_sample_wig_match.log'), mode='r', encoding='utf-8'):
        parsed_line = line.strip()
        parsed_line = parsed_line.split('\t')
        id_wig_match_dict.update({parsed_line[0]: parsed_line[1]})

    for samp_tag in identifier_tags:
        assert samp_tag in id_wig_match_dict.keys(), f'no wig file found for sample {samp_tag}'
    assert len(id_wig_match_dict) == len(identifier_tags), f'a sample was missed or included erroneously:\n{str(id_wig_match_dict)}'
    assert len(list(set(id_wig_match_dict.values()))) == len(list(id_wig_match_dict.values())), f'multiple distinct samples erroneously matching to same wig file:\n{str(id_wig_match_dict)}'

@pytest.mark.regular
def test_run_chrom_chr20():
    """
    Ensure 'rocco chrom' runs successfully for chr20
    """
    # check if command returns exit code '0'
    assert subprocess.run(['rocco', 'chrom', '--chrom', 'chr20', '--wig_path', os.path.dirname(tracks_chrom20_wig),
                            '--budget', chr20_budget, '-N', rr_iter, '--start', chr20_start, '--end', chr20_end,  '--solver', solver],
                          capture_output=True, encoding='utf-8').returncode == 0,'rocco chrom returncode != 0'
    # check if output file was created properly
    assert chr20_outfile in os.listdir()
    # ensure output file is nonempty (at least 'x' bytes)
    assert os.stat(chr20_outfile).st_size > peakfile_minsize

@pytest.mark.regular
def test_comp_ref_chr20():
    """
    Ensure the generated BED files are reasonably consistent with constant references.
    Wont match perfectly due the random simulation, but should exhibit broad similarity
    """
    cmd = ["bedtools", "jaccard", "-a", chr20_outfile , "-b", chr20_ref_file ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    data = proc.stdout.splitlines()[1]
    jaccard_chr20 = float(data.split("\t")[2])
    assert jaccard_chr20 >= jaccard_min

@pytest.mark.extra
def test_run_gwide():
    """
    Ensure `rocco gwide -p test_params.csv` runs successfully
    """
    # check if command returns exit code '0'
    assert subprocess.run(["rocco", "gwide", "-p", "test_params.csv","--outdir", gwide_outdir, "--combine", gwide_outfile, "-N", rr_iter,  '--solver', solver],
                          capture_output=True, encoding='utf-8').returncode == 0, 'rocco gwide returncode != 0'
    # check if output file was created properly
    assert gwide_outfile in os.listdir()
    # ensure output file is nonempty (at least 'x' bytes)
    assert os.stat(gwide_outfile).st_size > peakfile_minsize
    chrom_files = [x for x in os.listdir(gwide_outdir) if x.split('.')[-1] == 'bed']
    # ensure outdir contains chromosome-specific BED files
    assert len(chrom_files) == 2

@pytest.mark.extra
def test_comp_ref_gwide():
    """
    Ensure the generated BED files are reasonably consistent with constant references.
    Wont match perfectly due the random simulation, but should exhibit broad similarity
    """
    cmd = ["bedtools", "jaccard", "-a", gwide_outfile, "-b", gwide_ref_file]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    data = proc.stdout.splitlines()[1]
    jaccard_gwide = float(data.split("\t")[2])
    assert jaccard_gwide >= jaccard_min

@pytest.mark.extra
def test_gwide_metadata_and_group_comp():
    """
    Ensure correct splitting of samples by group when calling `--coldata`
    """
    # check if command returns exit code '0'
    assert subprocess.run(["rocco", "gwide", "-p", "test_params.csv",  '--coldata', metadata_file, '--group_col', 'status'],
                          capture_output=True, encoding='utf-8').returncode == 0, 'rocco gwide --coldata != 0'
    assert os.stat(group_A_bedfile).st_size > peakfile_minsize
    assert os.stat(group_B_bedfile).st_size > peakfile_minsize
    """
    Since the division of samples in the metadata file is arbitrary,
    we can expect consistency between `A` and `B` groups'
    results.
    """
    cmd = ["bedtools", "jaccard", "-a", group_A_bedfile, '-b', group_B_bedfile]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    data = proc.stdout.splitlines()[1]
    jaccard_group = float(data.split("\t")[2])
    assert jaccard_group >= jaccard_min

@pytest.mark.regular
def test_gwide_exclude_chroms():
    """
    Running `rocco gwide -p test_params.csv --exclude_chroms chr19,chr20` should not create
    any new BED files since `test_params.csv` only contains entries for these chromosomes.
    """
    num_beds_orig = len([x for x in os.listdir() if x.split('.')[-1] == 'bed'])
    # check if command returns exit code '0'
    assert subprocess.run(["rocco", "gwide", "-p", "test_params.csv",  '--exclude_chroms', 'chr19,chr20'],
                          capture_output=True, encoding='utf-8').returncode == 0, 'rocco gwide --exclude_chroms != 0'
    num_beds_curr = len([x for x in os.listdir() if x.split('.')[-1] == 'bed'])
    assert num_beds_orig == num_beds_curr

@pytest.mark.fixedStep
def test_sim_fs_chr19(var=1, samples=10):
    """
    simulate data in fixedStep format for chrom19
    """
    proc1 = subprocess.run(['python', 'sim.py', tracks_chrom19_wig_fs, str(var), str(samples), 'fs'],
                           capture_output=True, encoding='utf-8')
    assert proc1.returncode == 0, f'sim.py failed for {tracks_chrom19_wig_fs}'

@pytest.mark.fixedStep
def test_run_chrom_chr19_fixedStep():
    """
    Ensure 'rocco chrom' runs successfully for chr19 with fixedStep input
    """
    # remove output from previous run if present
    if os.path.exists(chr19_outfile):
        os.remove(chr19_outfile)
    # check if command returns exit code '0'
    proc = subprocess.run(['rocco', 'chrom','--chrom', 'chr19', '--wig_path',  os.path.dirname(tracks_chrom19_wig_fs),
                           '--budget', chr19_budget, '-N', rr_iter, '--start', chr19_start, '--end', chr19_end,
                           '--solver', solver, '--fixedStep', '--verbose'],
                          capture_output=True, encoding='utf-8')
    assert proc.returncode == 0, 'rocco chrom failed'
    # check if output file was created properly
    assert chr19_outfile in os.listdir()
    # ensure output file is nonempty (at least 'x' bytes)
    assert os.stat(chr19_outfile).st_size > peakfile_minsize

@pytest.mark.fixedStep
def test_comp_ref_chr19_fs():
    """
    same test as in `comp_ref_chr19` but assuming fixedStep format

    note `chr19_outfile` should exist after running `chrom --fixedStep`
    """
    cmd = ["bedtools", "jaccard", "-a", chr19_outfile , "-b", chr19_ref_file ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    data = proc.stdout.splitlines()[1]
    jaccard_chr19 = float(data.split("\t")[2])
    assert jaccard_chr19 >= jaccard_min


@pytest.mark.fixedStep
def test_run_chrom_chr19_identifiers_fs():
    """
    Ensure 'rocco chrom --identifiers' runs successfully and uses correct subset
    of samples in `test_identifiers.txt`
    """
    # check if command returns exit code '0'
    assert subprocess.run(['rocco', 'chrom','--chrom', 'chr19', '--wig_path',  os.path.dirname(tracks_chrom19_wig_fs),
                           '--budget', chr19_budget, '-N', rr_iter, '--start', '10000000', '--end', '20000000',
                           '--solver', solver, '--identifiers', identifiers_file, '--outdir', identifiers_outdir_fs, '--fixedStep'],
                          capture_output=True, encoding='utf-8').returncode == 0, 'rocco chrom returncode != 0'
    # check if output file was created in the output directory `identifiers_outdir`
    assert chr19_outfile in os.listdir(identifiers_outdir)
    # ensure output file is nonempty (at least 'x' bytes)
    assert os.stat(os.path.join(identifiers_outdir,chr19_outfile)).st_size > peakfile_minsize

    identifier_tags = []
    id_wig_match_dict = {}
    for line in open(identifiers_file, mode='r', encoding='utf-8'):
        identifier_tags.append(line.strip())
    for line in open(os.path.join(identifiers_outdir, 'chr19_sample_wig_match.log'), mode='r', encoding='utf-8'):
        parsed_line = line.strip()
        parsed_line = parsed_line.split('\t')
        id_wig_match_dict.update({parsed_line[0]: parsed_line[1]})

    for samp_tag in identifier_tags:
        assert samp_tag in id_wig_match_dict.keys(), f'no wig file found for sample {samp_tag}'
    assert len(id_wig_match_dict) == len(identifier_tags), f'a sample was missed or included erroneously:\n{str(id_wig_match_dict)}'
    assert len(list(set(id_wig_match_dict.values()))) == len(list(id_wig_match_dict.values())), f'multiple distinct samples erroneously matching to same wig file:\n{str(id_wig_match_dict)}'

@pytest.mark.fixedStep
def test_run_gwide_fs():
    """
    Ensure `rocco gwide -p test_params.csv` runs successfully
    """
    # check if command returns exit code '0'
    assert subprocess.run(["rocco", "gwide", "-p", "test_params.csv","--outdir", gwide_outdir_fs, "--combine", gwide_outfile_fs, "-N", rr_iter,  '--solver', solver],
                          capture_output=True, encoding='utf-8').returncode == 0, 'rocco gwide returncode != 0'
    # check if output file was created properly
    assert gwide_outfile_fs in os.listdir()
    # ensure output file is nonempty (at least 'x' bytes)
    assert os.stat(gwide_outfile_fs).st_size > peakfile_minsize
    chrom_files = [x for x in os.listdir(gwide_outdir_fs) if x.split('.')[-1] == 'bed']
    # ensure outdir contains chromosome-specific BED files
    assert len(chrom_files) == 2

@pytest.mark.fixedStep
def test_comp_ref_gwide_fs():
    """
    Can use same ref file as in comp_ref_gwide
    """
    cmd = ["bedtools", "jaccard", "-a", gwide_outfile_fs, "-b", gwide_ref_file]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    data = proc.stdout.splitlines()[1]
    jaccard_gwide = float(data.split("\t")[2])
    assert jaccard_gwide >= jaccard_min