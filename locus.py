import numpy as np

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
        return str(self.__dict__)

