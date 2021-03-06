# -*- coding: utf-8 -*-
# Site class used in the ising project.

from mathtools import *

# class that contains the matrices A(s=up) and
# A(s=down) for a given site s
class Site():

    # constructor
    def __init__(self, D=None, random=False, spinsym=True, \
                normalize=False, zeros=False, up=None, down=None):
        # set up and down to real symmetric matrices
        # with equal eigenvalues (spin symmetry)
        self.D = D
        if random:
            if spinsym:
                self.up, self.down = enforce_equal_eigenvalues(
                    rand_sym_mat(D),rand_sym_mat(D), symmetric=True)
            else:
                self.up, self.down = rand_sym_mat(D),rand_sym_mat(D)
        elif zeros:
            self.up, self.down = np.zeros((D,D)),np.zeros((D,D))
        elif up is not None and down is not None:
            self.up = np.copy(up)  # WAS np.copy(up)
            self.down = np.copy(down)  # WAS np.copy(down)

    # predefined constructor cases:
    # - fill with zeros
    @classmethod
    def zeros(cls, D):
        return cls(D=D, zeros=True)
    # - fill with random r in [1,0)
    @classmethod
    def random(cls, D, spinsym=True, normalize=False):
        return cls(D=D, random=True, spinsym=spinsym, \
                normalize=normalize)
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

    # '*' operator (scalar multiplication)
    # @profile
    def __mul__(self, scalar):
        up = self.up * scalar
        down = self.down * scalar
        return Site.compose(up,down)
    __rmul__ = __mul__

    # '+=' operator
    # @profile
    def __iadd__(self,other):
        self.up += other.up
        self.down += other.down
        return self
    # '-=' operator
    def __isub__(self,other):
        self.up -= other.up
        self.down -= other.down
        return self

    # '+' operator
    # @profile
    def __add__(self,other):
        up = self.up + other.up
        down = self.down + other.down
        return Site.compose(up,down)
    # '-' operator
    def __sub__(self, other):
        up = self.up - other.up
        down = self.down - other.down
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

# generalization of the numpy.sign function to arrays and Site
def sign(input, iterable=False, isSite=False):
    if iterable:
        return np.array([sign(element, iterable=False, \
                isSite=isSite) for element in input])
    elif isSite:
        return Site.compose(np.sign(input.up), np.sign(input.down))
    else:
        return np.sign(object)

# generalization of the numpy.multiply function to arrays and Site
def multiply(a,b, iterable=False, isSite=False):
    if iterable:
        if len(a) != len(b):
            raise ValueError('Input arrays should have equel length.')
        else:
            return np.array([multiply(a[i],b[i], iterable=False, \
                    isSite=isSite) for i in range(len(a))])
    elif isSite:
        return Site.compose(np.multiply(a.up, b.up), \
                    np.multiply(a.down, b.down))
    else:
        return np.multiply(a,b)

def main():
    x = Site.random(2)
    y = Site.random(2)
    z = Site.random(2)
    print('x:\n', x, '\ny:\n', y, '\nz:\n', z)

    z += 3 * x
    print('z:\n', z)

    y += x
    print('y:\n', y)
    x = Site.zeros(2)
    print('y:\n', y)
    print('z:\n', z)

    # # some debugging for multiply
    # x = Site.random(2)
    # y = Site.random(2)
    # print('x:\n', x, '\ny:\n', y)
    # a = np.multiply(x.up, y.up)
    # b = np.multiply(x.down, y.down)
    # print('a:\n', a, '\nb:\n', b)
    # z1 = Site.compose(a,b)
    # z2 = multiply(x,y, isSite=True)
    # print('z1:\n', z1)
    # print('z2:\n', z2)

# don't run main if included 
if __name__ == "__main__":
    main()
