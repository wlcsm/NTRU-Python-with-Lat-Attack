from sympy.abc import x
from sympy import ZZ, Poly
import numpy as np


class NTRU:
    """ Stores information about the generic NTRU parameters """
    def __init__(self, N, p, q):
        self.N = N
        self.p = p
        self.q = q
        self.R = Poly(x**N - 1, x).set_domain(ZZ)


def load_NTRU(fileName, isPriv):
    """ Loads the NTRU system for the Private and Public domains """
    param_file = np.load(fileName, allow_pickle=True)
    ntru = NTRU(*[int(param_file[i]) for i in ['N', 'p', 'q']])
    if isPriv: # Private domain, used for decryption
        f = asPoly(param_file['f'].astype(np.int))
        f_p = asPoly(param_file['f_p'].astype(np.int))
        return ntru, f, f_p
    else: # Public domain, used for encryption
        h = asPoly(param_file['h'].astype(np.int))
        return ntru, h


def asPoly(arr):
    """ Convert a NumPy array to a polynomial """
    return Poly(arr[::-1],x).set_domain(ZZ)


def asArr(poly, N):
    """ Convert polynomial to a NumPy array """
    tmp = poly.all_coeffs()[::-1]
    return np.pad(tmp, (0, N - len(tmp)))
