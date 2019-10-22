from sympy.abc import x
from sympy import ZZ, Poly


class NTRUPoly:
    """ Polynomial Class for NTRU system """
    coeff = []

    def __init__(self, generator, N, k):
        """ Initializes the polynomial type """
        self.R = Poly(x ** N - 1, x).set_domain(ZZ)
        if type(generator) == sympy.polys.polytools.Poly:
            self.poly = (generator % self.R).trunc(k)
        else
            self.poly = Poly(generator[::-1], x).set_domain(ZZ)

    def toArr(self):
        """ Export the polynomial as a numpy array """
        poly_as_arr = self.poly.all_coeffs()[::-1]
        return np.array(poly_as_arr, (0, self.N - len(poly_as_arr)))

    def trunc(self, k):
        self.poly = self.poly.trunc(k)


