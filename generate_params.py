import numpy as np

file = open("parameters", "w")

file.truncate()

file.write("size discretization iterations initial_population birth_rate death_rate birth_variance death_variance seed\n")

seed = 1
for i in np.linspace(1e-5, 1e-3, 30):
    for j in np.linspace(1e-5, 1e-3, 30):
        file.write("1 5000 5000 300 5e-3 1e-4 " + str(i) + " " + str(j) + " " + str(seed) + "\n")
        seed += 1

file.close()
