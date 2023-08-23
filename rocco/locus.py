class Locus:
    r"""
    `Locus` objects are the $n \approx \frac{|\mathscr{L}|}{L}$ nodes
    dividing genomic region $\mathscr{L}$ into segments of $L$ nucleotides.

    Attributes:
        position (int): Nucleotide position of Locus `self`.
        left (Locus): Points to the previous adjacent locus.
        right (Locus): Points to the subsequent adjacent locus.
        size (int): Size (in nucleotides) of the Locus `self`.
            Corresponds to $L$ before merging.
        parent (Locus): Parent locus of the Locus `self`.
            Only defined for Locus objects merged by `Loci.combine()`.
        subloci (list): List of subloci comprising Locus objet `self`.
            Only defined for Locus objects that were created from merged
            loci with `Loci.combine()`.
        sig_data (np.ndarray): $K \times 1$ vector of replicates' signal
            values at Locus `self`. Corresponds to $\mathbf{s}_i$.
        accessible (float): a value between 0 and 1. Corresponds to
            $\ell_i$.

    Note:
        After optimization, ROCCO merges contiguous, accessible Locus objects.
            Meaning `Locus.size` may be larger than the initial locus size, and
            is in this way inconsistent with the definition of $L$ in the paper.
    """
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


