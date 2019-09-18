import numpy as np
import numpy.linalg as LA

def to_array(sigma):
    n = []
    for i in range(4):
        n.append(sigma%2)
        sigma = sigma//2
    return 1-2*np.array(n)

def Ham(s):
    return -s[0]*s[1]-s[1]*s[3]-s[0]*s[2]-s[2]*s[3]

def oloc(s, beta):
    w = np.arccosh(np.exp(beta))/2
    W = np.array([[w,w,0,0],[0,w,0,w],[w,0,w,0],[0,0,w,w]])
    theta = np.matmul(W,s)
    return np.concatenate([s, np.tanh(theta), np.kron(s,np.tanh(theta))])

def ent_ee(evals):
    evals = np.array([e for e in evals if e > 1e-10])
    return -sum(evals*np.log2(evals))

if __name__ == '__main__':
    beta = 1.4
    w = np.arccosh(np.exp(beta))/2
    W = np.array([[w,w,0,0],[0,w,0,w],[w,0,w,0],[0,0,w,w]])

    n = 4
    m = 4
    
    s = np.array([to_array(sigma) for sigma in range(16)])

    Z = sum(np.exp(-beta*Ham(si)) for si in s)

    rbm_prob = np.zeros(16)
    for sidx, si in enumerate(s):
        rbm_prob[sidx] = np.prod(np.cosh(np.matmul(W,si)))

    rbm_prob = rbm_prob*rbm_prob
    rbm_prob /= sum(rbm_prob)
    
    s2 = np.zeros((n+m+n*m, n+m+n*m))
    omean = np.zeros((1,n+m+n*m))
    for sidx, si in enumerate(s):
        p = np.exp(-beta*Ham(si))/Z
        o = oloc(si, beta)
        o = o[np.newaxis,:]
        s2 += p*np.matmul(np.transpose(o),o)
        omean += p*o

    smat = s2 - np.matmul(np.transpose(omean),omean)

    w,v = LA.eigh(smat)
    for i in range(5):
        ev = v[-i-1,n+m:]
        ev /= LA.norm(ev)
        psi = ev.reshape(4,4)
        rho = np.matmul(np.transpose(psi),psi)
        evals = LA.eigvalsh(rho)
        print(w[-i-1], evals)
