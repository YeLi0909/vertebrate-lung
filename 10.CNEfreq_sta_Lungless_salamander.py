#!/usr/bin/env python
import os
import re
from collections import defaultdict

blast_dir = 'Blastout'
out_file = 'cne_freq.txt'
evalue_threshold = 1e-10

score = defaultdict(int)
species_num = 0

for filename in os.listdir(blast_dir):
    if filename.endswith('.out'):
        species_num += 1
        filepath = os.path.join(blast_dir, filename)
        file_cnes = set()  # 用于跟踪每个文件中的唯一CNE
        with open(filepath, 'r') as f:
            for line in f:
                values = re.split(r'\s+', line.strip())
                if len(values) >= 11:  # 确保有足够的列
                    query_id = values[0]
                    evalue = float(values[10])  # E-value 通常在第11列 (索引10)
                    if evalue < evalue_threshold:
                        match = re.match(r".*:(cne\d+)", query_id)
                        if match:
                            cne_id = match.group(1)
                            file_cnes.add(cne_id)
        
        # 对文件中的每个唯一CNE增加计数
        for cne_id in file_cnes:
            score[cne_id] += 1

with open(out_file, 'w') as oh:
    oh.write('cne_id\tfrequency\trate\n')
    for cne_id, frequency in score.items():
        rate = frequency / species_num
        oh.write(f'{cne_id}\t{frequency}\t{rate:.4f}\n')