import numpy as np
import matplotlib.pyplot as plt

from matplotlib import animation, rc



if __name__=='__main__':
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate = 1800)

    fig = plt.figure()

    ax = plt.gca()

    plots = []

    ax.set_xlim(0,7168)
    ax.set_ylim(1e-22, 1e+3)
    ax.set_yscale('log')

    for i in range(10,101):
        filename = "SMAT_{:03d}.dat".format(i)
        d = np.loadtxt(filename)
        ys = d[1:]
        beta = d[0]
        p, = ax.plot(range(ys.size), ys[::-1], color='b')
        beta_t = r"$\beta={:.2f}$".format(beta)
        t = ax.text(0.6, 0.8, beta_t, transform=ax.transAxes, fontsize=18)

        plots.append([p,t])

    ani = animation.ArtistAnimation(fig, plots, interval=50, blit=True)
    ani.save('smat.mp4', writer = writer)
    #plt.show()
