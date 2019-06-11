import numpy as np
import matplotlib.pyplot as plt


if __name__=='__main__':
    fig = plt.figure(figsize=(5,4))

    ax = plt.gca()

    betas = [0.1, 0.45, 0.58, 0.65, 0.9]
    plots = []

    ax.set_xlim(0,5000)
    ax.set_ylim(1e-8, 1e+3)
    ax.set_yscale('log')

    for beta in betas:
        print(int(beta*100))
        filename = "SMAT_{:03d}.dat".format(int(beta*100))
        d = np.loadtxt(filename)
        ys = d[1:]
        beta = d[0]
        p, = ax.plot(range(ys.size), ys[::-1], linewidth=2, label=r"$\beta={:.2f}$".format(beta))

    ax.tick_params(axis='both', which='major', labelsize=20)
    plt.legend(fontsize=18)
    #plt.show()
    plt.savefig('plots.eps',bbox_inches="tight")
