class Potential(object):
    pass

class ChaoticPotential(Potential):
    def __init__(self, m=1):
        self.m = m

    def __call__(self, phi):
        return self.m**2*phi**2/2

    def d(self, phi):
        return self.m**2*phi

    def dd(self, phi):
        return self.m**2
