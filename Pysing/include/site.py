# -*- coding: utf-8 -*-
# Site class used in the ising project.

from mathtools import *

# class that contains the matrices A(s=up) and
# A(s=down) for a given site s
class Site():

    # constructor
    def __init__(self, D=None, random=False, zeros=False, up=None, down=None):
        # set up and down to real symmetric matrices
        # with equal eigenvalues (spin symmetry)
        if random:
            self.up, self.down = enforce_equal_eigenvalues(
                rand_sym_mat(D),rand_sym_mat(D), symmetric=True)
        elif zeros:
            self.up, self.down = enforce_equal_eigenvalues(
                np.zeros((D,D)),np.zeros((D,D)), symmetric=True)
        elif up is not None and down is not None:
            self.up = np.copy(up)
            self.down = np.copy(down)

    # predefined constructor cases
    @classmethod
    def zeros(cls, D):
        return cls(D=D, zeros=True)
    @classmethod
    def random(cls, D):
        return cls(D=D, random=True)
    @classmethod
    def compose(cls,up,down):
        return cls(up=up, down=down)


    # get the up (i=0) or down (i=1) matrix
    def get(self,i):
        return [self.up, self.down][i]
    # set the up (i=0) or down (i=1) matrix
    def set(self,i,value):
        if i == 0:
            self.up = value
        elif i == 1:
            self.down = value
    # get a random matrix
    def rand(self):
        return self.get(random.randint(2))


    # '+=' operator
    def __iadd__(self,other):
        self.up += other.up
        self.down += other.down
        return self

    # '*' operator (scalar multiplication)
    def __mul__(self, scalar):
        up = self.up * scalar
        down = self.down * scalar
        return Site.compose(up,down)
    __rmul__ = __mul__

    # '+' operator
    def __add__(self,other):
        up = self.up + other.up
        down = self.down + other.down
        return Site.compose(up,down)

    # operator [] for getting
    # key: 0 == up, 1 == down
    def __getitem__(self, key):
        return [self.up, self.down][key]
    # operator [] for assignment
    # key: 0 == up, 1 == down
    def __setitem__(self, key, value):
        if key == 0:
            self.up = value
        elif key == 1:
            self.down = value

    # representation (for printing)
    def __repr__(self):
        return "Up: %r\nDown:%r\n" \
                % (self.up, self.down)
