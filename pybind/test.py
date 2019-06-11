import Ising2D
import numpy as np

model = Ising2D.Ising2D(10,10)

for beta in np.linspace(0.05, 1.0, 20):
    sampler = Ising2D.WolffSampler(model, beta)
    sampler.set_seed(2)
    sampler.randomize_conf()

    #thermalizing
    for i in range(100):
        sampler.sweep()

    #sampling
    energies = []
    for i in range(5000):
        sampler.sweep()
        energies.append(model.energy(sampler.conf))

    print("{}\t{}".format(beta, np.mean(energies)/model.size))


