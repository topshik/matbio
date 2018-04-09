#!/usr/bin/env python3
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
import time
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
    start_time = time.time()
    print('Started')
    for i in range(repeats):
        params['seed'] = str(int(params['seed']) + 1)
        out = subprocess.Popen(['../' + filename] + list(params), universal_newlines=True, stdout=subprocess.PIPE)
        with open('{}_{}'.format(
                params['birth_variance'],
                params['death_variance'],
        ), 'w+') as f:
            f.write('{} {}\n'.format(
                params['birth_variance'],
                params['death_variance'],
            ))
            while True:
                row = out.stdout.readline()
                if row == '':
                    break
                f.write('{}\n'.format(row.split()[1]))
                f.flush()
    print('Ended: {}'.format(time.time() - start_time))

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
        launch(program_file, params, repeats)
        # tp.apply_async(launch, [program_file, params, repeats])

    tp.close()
    tp.join()
