import os
import argparse
import subprocess


def sort_combine_bed(outfile,dir_='.'):
    """sorts and combines chromosome-specific bedfiles"""
    if dir_[-1] == '/':
        dir_ = dir_[0:-1]
    filenames = [f_ for f_ in os.listdir(dir_)
                 if f_.split('.')[-1] == 'bed' and 'ROCCO_out' in f_]
    filenames = sorted(filenames,key=lambda x: int(val) if (val:=x.split('_')[2][3:]).isnumeric() else ord(val))
    
     
    filenames = [dir_ + '/' + fname for fname in filenames]
    with open(outfile, 'w') as outfile_:
        for fname in filenames:
            p = subprocess.Popen(('cat', fname),
                             stdout=outfile_.fileno())
            p.wait()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--dir',
                        default='.',
                        help="directory containing bed files to be combined")
    parser.add_argument('-o', '--outfile',
                        default='ROCCO_combined.bed',
                        help="file containing sorted/merged results")
    
    args = vars(parser.parse_args())
    sort_combine_bed(args['outfile'], args['dir'])
    
if __name__ == "__main__":
    main()
