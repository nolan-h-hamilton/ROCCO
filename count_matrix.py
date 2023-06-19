import os
import argparse
import pandas as pd
from pybedtools import BedTool
from pysam import head, SamtoolsError

def is_alignment(filepath):
    if 'bai' == filepath.split('.')[-1]:
        return False
    try:
        head("-n 1", filepath)
        return True
    except SamtoolsError:
        return False

def proc_bed(bedfile, header=False):
    i = 0
    bed_string = ""
    with open(bedfile, 'r') as bf:
        for line in bf:
            if header and i == 0:
                i += 1
                continue

            line = line.strip()
            line = line.split('\t')
            if len(line) == 6:
                return BedTool(bedfile)

            if len(line) > 6:
                new_line = "\t".join(line[0:5]) + '\n'
                bed_string += new_line
                continue

            new_line = "\t".join(line[0:3])
            if len(line) == 3:
                # create name for peak: chr_start_end
                name = '_'.join([line[0],line[1],line[2]])
                # assign 1000 for score and . for strand
                new_line += "\t{}\t{}\t{}".format(name,'1000','.')

            elif len(line) == 4:
                newline += "\t{}".format(line[3])
                new_line += "\t{}\t{}".format('1000','.')

            elif len(line) == 5:
                newline += "\t{}".format(line[3])
                newline += "\t{}".format(line[4])
                new_line += "\t{}".format('.')

            new_line += '\n'
            bed_string += new_line

    return BedTool(bed_string, from_string=True)


def peak_count_matrix(bedfile, bams=None, outfile=None,  bdir='.'):
    """Runs `multicov` for a specified `bedfile` and BAMs"""
    bamfiles = []
    if bdir.split('/')[-1] == '/':
        bdir = bdir[0:-1]

    if bams is not None:
        bam_df = pd.read_csv(bams,
                             header=None,
                             index_col=False,
                             skiprows=[0],
                             sep="\t")

        IDs = [x for x in list(bam_df.iloc[:,0])]
        init_IDs= [x for x in IDs]
        for i,ID in enumerate(IDs):
            if not os.path.exists(ID):
                for f_ in os.listdir(bdir):
                    f_ = f_.strip()
                    f_ = bdir + '/' + f_
                    if ID in f_ and  is_alignment(f_):
                        IDs[i] = f_

        bamfiles = IDs
    else:
        print('no sample metadata was provided...exiting')
        return None
    bed = proc_bed(bedfile)
    bed.multi_bam_coverage(bams=bamfiles,
                               output=outfile)

    cols = ['chrom','start','end','name','score','strand']
    for bf in init_IDs:
        cols.append(bf)
    count_df = pd.read_csv(outfile,sep='\t',names=cols)
    count_df = count_df.drop(
        columns=['chrom','start','end','score','strand'])
    count_df.to_csv(outfile,sep='\t',index=False)
    return count_df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--peakfile', help="BED file containing peaks")
    parser.add_argument('-d', '--bamfiles', default=None,
                        help="sample metadata file. leftmost column should list sample\
                          identifiers for which there exists an alignment file named \
                          accordingly in the directory `--")
    parser.add_argument('--bamdir', type=str, default='.')
    parser.add_argument('-o','--outfile', default='ROCCO_peak_counts.tsv')
    args = vars(parser.parse_args())
    peak_count_matrix(args['peakfile'], args['bamfiles'],
                      args['outfile'],bdir=args['bamdir'])

if __name__ == "__main__":
    main()
