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
from string import ascii_lowercase, ascii_uppercase, digits

dictionary = digits + ascii_lowercase + ascii_uppercase
of = False

def b36encode(arr):
    result = ''
    for num in arr:
        try:
            result += dictionary[int(num)]
        except IndexError as e:
            result += '+'
            if not of:
                print('Overflow: {}'.format(num))
                of = True
    return result

def launch(filename, params, repeats):
    path = '{}-{}-{}-{}'.format(
        params['birth_rate'],
        params['death_rate'],
        params['birth_variance'],
        params['death_variance'],
    )
    if not os.path.exists(path):
        os.makedirs(path)
    # os.chdir(path)

    print('Started')
    for i in range(repeats):
        out = subprocess.check_output(['../' + filename] + list(params), universal_newlines=True)
        out = np.array(out.split()).reshape(-1, 3 + int(params['discretization']))
        with open('{}/{}_{}_traj'.format(path, path, i), 'w+') as f:
            for row in out:
                f.write(b36encode(row[3:]) + '\n')
        with open('{}/{}_{}_pops'.format(path, path, i), 'w+') as f:
            for row in out:
                f.write(row[1] + '\n')


    # f = plt.figure()
    # plt.title('{}, br: {}, dr: {}, bv: {}, dv: {}'.format(
    #     params['birth_rate'],
    #     params['death_rate'],
    #     params['birth_variance'],
    #     params['death_variance'],
    # ))
    # plt.xlabel('iteration')
    # plt.ylabel('population')
    # plt.plot(np.arange(int(params['iterations'])), out[:,1].astype(np.int))
    # f.savefig(path + '.png')


if __name__ == '__main__':
    usage_string = 'Usage: ' + sys.argv[0] + ' program_file parameters_file repeats [num_of_cores]'
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        raise ValueError('Invalid number of command line arguments. ' + usage_string)


    program_file = sys.argv[1]
    params_file = sys.argv[2]
    repeats = int(sys.argv[3])
    try:
        cores = int(sys.argv[4])
    except (IndexError) as e:
        cores = 1
    params_list = pd.read_csv(params_file, sep='\s+').astype(str)

    directory = 'Sim_' + time.strftime('%d.%m.%y-%H.%M.%S')
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.chdir(directory)

    tp = ThreadPool(cores)
    for i, params in params_list.iterrows():
        # launch(program_file, params, repeats)
        tp.apply_async(launch, [program_file, params, repeats])

    tp.close()
    tp.join()