"""
Locus object

Args:
    position: nucleotide position of the locus.
    left (Locus): points to previous adjacent locus
    right (Locus): points to subsequent adjacent locus
    size (int): size (nt) of locus.
    parent (Locus): parent locus of the current locus\
        defined for initial loci after Loci.combine()
    subloci (list): list of subloci contained within the current locus\
        defined for merged loci after Loci.combine()
    sig_data (list): list of replicate's signal values at the locus
    accessible (float): A value between [0,1] measuring accessibility\
        restricted to 0,1 in a complete (integral) solution.

"""
class Locus:
    def __init__(self, position=None,left=None,right=None,
                 size=None, parent=None, subloci=None,
                 sig_data=None, accessible=None):
        self.position = position
        self.left = left
        self.right = right
        self.size = size
        self.parent = parent
        self.subloci = subloci
        self.sig_data = sig_data
        if accessible is None:
            self.accessible = -1
        else:
            self.accessible = accessible


    def __str__(self):
        """
        Returns a string representation of the locus

        Returns:
            a string repr. of a dictionary with a key for each\
                field of the Locus object
        """
        return str(self.__dict__)

