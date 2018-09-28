import numpy
from scipy.constants import c, hbar, G

lp = numpy.sqrt(hbar*G/c**3)
mp = numpy.sqrt(hbar*c/G) 
tp = numpy.sqrt(hbar*G/c**5) 
Mpc = 3.086e+22
H0_h = 100 *1000 / Mpc * tp
DA = 13.91726e3
kstar = 0.05 * Mpc**-1 /( lp**2 / tp)
