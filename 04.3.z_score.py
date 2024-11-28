#! /usr/bin/env python
import re
import numpy as np
from scipy import stats

ppi = {}

all_ratio = 'all_lung.ratio'
shark_ratio = 'all_dogshark.ratio'
out_file = 'output.z_score'

with open(all_ratio,'r') as F:
    F.readline()
    for line in F:
        values = re.split('\s+', line.strip())
        ppi_name = values[0]
        hvg_tag = values[1]
        percents = [float(x) for x in values[2:] if x != 'NA']

        mean = np.mean(percents)
        std = np.std(percents)

        if hvg_tag not in ppi:
            ppi[hvg_tag] = {}
        ppi[hvg_tag][ppi_name] = {
            'mean': mean,
            'std': std,
            'percents': percents
        }

with open(shark_ratio,'r') as F, open(out_file,'w') as OH:
    F.readline()
    for line in F:
        values = re.split('\s+', line.strip())
        ppi_name = values[0]
        percents = [float(x) for x in values[1:] if x != 'NA']
        if len(percents) == 0: continue
        max_percent = max(percents)

        for hvg_tag in ppi:
            if ppi_name not in ppi[hvg_tag]:
                continue
            mean = ppi[hvg_tag][ppi_name]['mean']
            std = ppi[hvg_tag][ppi_name]['std']

            # z_score = (max_percent - mean)/std
            # p_value = stats.norm.cdf(z_score)

            if std != 0:
                z_score = (max_percent - mean) / std
                p_value = stats.norm.cdf(z_score)
            else:
                z_score = float('inf') if max_percent > mean else float('-inf')
                p_value = 1 if max_percent > mean else 0

            fold_change = 999 if max_percent == 0 else mean / max_percent

            OH.write(f'{ppi_name}\t{hvg_tag}\t{max_percent}\t{mean}\t{std}\t{z_score}\t{p_value}\t{fold_change}\n')
