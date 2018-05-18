import numpy as np

file = open("parameters", "w")

file.truncate()

file.write("size discretization iterations initial_population birth_rate death_rate birth_variance death_variance seed\n")

seed = 1
for i in np.linspace(1e-4, 1e-2, 5):
    for j in np.linspace(1e-4, 1e-2, 5):
        file.write("10 50000 50000 1000 5e-2 1e-3 " + str(i) + " " + str(j) + " " + str(seed) + "\n")
        seed += 1

file.close()
