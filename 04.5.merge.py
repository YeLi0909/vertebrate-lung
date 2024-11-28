#! /usr/bin/python
import csv
from collections import defaultdict

# 定义输入文件和输出文件
file1 = 'output.z_score'
file2 = 'all_lungfish.ratio.zscore'
file3 = 'all_bichir.ratio.zscore'
file4 = 'all_mouse.ratio.zscore'
file5 = '../18.ppi_significant_multi/all_lung.ratio.global_fdr.validate'
file6 = '../22.abca3_sftpb/min_p_value_combinations.tsv'
output_file = 'merge.out'

# 使用字典来存储数据
data = defaultdict(lambda: [''] * 10)

# 读取文件并处理数据
def process_file(filename, column_indices, data_index):
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        if filename != file1:
            next(reader)  # 跳过标题行
        for row in reader:
            key = (row[0], row[1].replace('LUNG_', ''))
            for i, col_index in enumerate(column_indices):
                data[key][data_index + i] = row[col_index]

# 处理五个文件
process_file(file1, [0, 1, 2, 3, 5, 7], 0)
process_file(file2, [-1], 6)
process_file(file3, [-1], 7)
process_file(file4, [-1], 8)
process_file(file5, [2], 9)

# 处理min_p_value_combinations.tsv文件
with open(file6, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    next(reader)  # 跳过标题行
    min_p_value_data = {row[0]: row[1:] for row in reader}  # 使用 combination 作为唯一的 key，存储剩余数据

# 写入输出文件
with open(output_file, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t')

    # 写入标题行
    writer.writerow(['combination', 'hvg_tag', 'shark_ratio', 'mean_ratio', 'z_score', 'fold_change', 'Lungfish_ZScore', 'Bichir_ZScore', 'Mouse_Zscore', 'validate_tissue_count', 'cell_type', 'ratio', 'p_value', 'stats'])

    # 写入所有行，将空值替换为NA
    for key, row in data.items():
        combination_key = key[0]  # 获取 combination 作为唯一的 key
        row = ['NA' if cell == '' else cell for cell in row]  # 替换空值为NA
        if combination_key in min_p_value_data:
            min_p_values = ['NA' if cell == '' else cell for cell in min_p_value_data[combination_key]]  # 替换空值为NA
            writer.writerow(row + min_p_values)
        else:
            writer.writerow(row)

print(f"合并完成，输出文件为: {output_file}")