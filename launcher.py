#!/usr/bin/env python
from multiprocessing.pool import ThreadPool
import os
import time
import subprocess
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def launch(filename, params):
    print('Started')
    out = subprocess.check_output(['../' + filename] + list(params), universal_newlines=True)
    out = np.array(out.split()).reshape(-1, 3)

    f = plt.figure()
    plt.title('br: {}, dr: {}, bv: {}, dv: {}'.format(
        params['birth_rate'],
        params['death_rate'],
        params['birth_variance'],
        params['death_variance'],
    ))
    plt.xlabel('iteration')
    plt.ylabel('population')
    plt.plot(np.arange(int(params['iterations'])), out[:,1].astype(np.int))
    f.savefig('{}-{}-{}-{}.png'.format(
        params['birth_rate'],
        params['death_rate'],
        params['birth_variance'],
        params['death_variance'],
    ))

program_file = sys.argv[1]
params_file = sys.argv[2]
try:
    cores = int(sys.argv[3])
except (IndexError) as e:
    cores = 1
params_list = pd.read_csv(params_file, sep='\s+').astype(str)

directory = 'Sim_' + time.strftime('%d.%m.%y-%H.%M.%S')
if not os.path.exists(directory):
    os.makedirs(directory)
os.chdir(directory)

tp = ThreadPool(cores)
for i, params in params_list.iterrows():
    # launch(program_file, params)
    tp.apply_async(launch, [program_file, params])

tp.close()
tp.join()