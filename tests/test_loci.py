import sys
sys.path.append('../')
from locus import Locus
from loci import Loci
import numpy as np
import pandas as pd

def log(txt, script):
    print(f'{script}: {txt}')

signal_matrix = np.array([[0,0,0,0,0,2,10,2,0,0], [0,0,0,0,0,1,9,1,0,0],[0,0,0,0,0,3,11,3,0,0]])
signal_matrix = pd.DataFrame(signal_matrix)
signal_matrix.columns = [i*50 for i in range(0,10)]
loci_object = Loci()
size_ = 50
for i, loc in enumerate(signal_matrix.columns):
    new_loc = Locus(position=0 + i * size_,
                    size=size_,
                    sig_data=signal_matrix[loc])
    if i in [5,6,7]:
        new_loc.accessible = 1
    loci_object.append_locus(new_loc)
assert loci_object.is_empty() == False, "loci_object empty"
assert loci_object.head.position == 0, "incorrect position for first locus/first locus d.n.e"
assert loci_object.head.right.position == 50, "first locus points to wrong locus/first locus `right` d.n.e"
assert loci_object.head.right.left.position == 0, "locus left/right attributes not set correctly"


list_rep = loci_object.as_list()
assert len(list_rep) == 10, "list representation missing a locus"
for loc_obj in list_rep:
    assert len(loc_obj.sig_data) == 3, "missing locus signal data"


loci_object.combine_selected()
list_rep = loci_object.as_list()
i = 0
for loc_obj in loci_object:
    if i == 5:
        assert loc_obj.size == 150, "combined locus has incorrect size"
        i += 1
    else:
        assert loc_obj.size == 50, "combined locus has incorrect size"
        i += 1
assert i == 8, "combined loci object should have length 8"

log('all cases passed', sys.argv[0])