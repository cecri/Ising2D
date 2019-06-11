import scipy.special as sp
import numpy as np
def energy(beta):
    k = 1.0/(np.sinh(2*beta)**2)
    return -1.0/np.tanh(2.0*beta)*(1.0 + 2.0/np.pi *(2.0*np.tanh(2*beta)**2-1.0)*sp.ellipk(4*k/(1+k)**2))

if __name__ == '__main__':
    for i in range(1,101):
        beta = i*0.01
        print("{:.2f}\t{:.5f}".format(beta, energy(beta)))

