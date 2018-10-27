#!/usr/bin/env python3
import sys
import subprocess
import pandas as pd
import numpy as np

def write_dap_full(x,y,prefix,r):
    names = np.array([('geno', i+1, f'group{r}') for i in range(x.shape[1])])
    with open(f'{prefix}.data', 'w') as f:
        print(*(['pheno', 'pheno', f'group{r}'] + list(np.array(y).ravel())), file=f)
        np.savetxt(f, np.hstack((names, x.T)), fmt = '%s', delimiter = ' ')
        
def run_dap_full(prefix, args):
    cmd = ['dap-g', '-d', f'{prefix}.data', '-o', f'{prefix}.result', '--output_all'] + ' '.join(args).split()
    subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()    
           
def write_dap_ss(z,prefix):
    '''z-score vesion of dap input is the same as FINEMAP'''
    ids = np.array([str(i+1) for i in range(z.shape[0])])
    with open(f'{prefix}.z', 'w') as f:
        np.savetxt(f,  np.vstack((ids, z)).T, fmt = '%s', delimiter = ' ')

def run_dap_z(ld, prefix, args):
    cmd = ['dap-g', '-d_z', f'{prefix}.z', '-d_ld', ld, '-o', f'{prefix}.result', '--all'] + ' '.join(args).split()
    subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()    
    
def extract_dap_output(prefix):
    out = [x.strip().split() for x in open(f'{prefix}.result').readlines()]
    pips = []
    clusters = []
    still_pip = True
    for line in out:
        if len(line) == 0:
            continue
        if len(line) > 2 and line[2] == 'cluster_pip':
            still_pip = False
            continue
        if still_pip and (not line[0].startswith('((')):
            continue
        if still_pip:
            pips.append([line[1], float(line[2]), float(line[3]), int(line[4])])
        else:
            clusters.append([len(clusters) + 1, float(line[2]), float(line[3])])
    pips = pd.DataFrame(pips, columns = ['snp', 'snp_prob', 'snp_log10bf', 'cluster'])
    clusters = pd.DataFrame(clusters, columns = ['cluster', 'cluster_prob', 'cluster_avg_r2'])
    clusters = pd.merge(clusters, pips.groupby(['cluster'])['snp'].apply(','.join).reset_index(), on = 'cluster')
    return {'snp': pips, 'set': clusters}

def dap_single(x, y, prefix, r, args):
    write_dap_full(x,y,prefix,r)
    run_dap_full(prefix,args)
    return extract_dap_output(prefix)

def dap_single_z(z, ld, prefix, args):
    write_dap_ss(z,prefix)
    run_dap_z(ld,prefix,args)
    return extract_dap_output(prefix)

def dap_batch(X, Y, prefix, *args):
    return dict([(f'V{r+1}', dap_single(X, Y[:,r], f'{prefix}_condition_{r+1}', r+1, args)) for r in range(Y.shape[1])])

def dap_batch_z(z, ld, prefix, *args):
    return dict([(f'V{r+1}', dap_single_z(z[:,r], ld, f'{prefix}_condition_{r+1}', args)) for r in range(z.shape[1])])

import os
from dsc.dsc_io import load_rds, save_rds
import tempfile
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")

input_file = os.path.expanduser(sys.argv[1])
output_file = os.path.expanduser(sys.argv[2])
args = sys.argv[3:]
data = load_rds(input_file)['data']
cache = tempfile.NamedTemporaryFile(suffix = '.dap')
posterior = dap_batch(data['X'], data['Y'], cache.name, ' '.join(args))
save_rds(posterior, output_file + '.rds')
