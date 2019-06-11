import Ising2D

model = Ising2D.Ising2D(4,5)
print(model.neighbors(0))
sampler = Ising2D.WolffSampler(model, 100.0)
sampler.set_seed(2)
sampler.randomize_conf()
for i in range(100):
    sampler.sweep()
    print(sampler.conf)

