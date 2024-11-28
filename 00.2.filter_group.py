import os
from ete3 import Tree
from collections import defaultdict

# 创建输出目录
output_dir = "all_genes_processed"
os.makedirs(output_dir, exist_ok=True)

# 初始化数据结构
data = defaultdict(set)

# 打开并读取文件
orthogroup_id = ''
with open("orthogroup_summary.txt", "r") as file, open("putative_orthogroups.txt","w") as output:
    for line in file:
        if line.startswith("Orthogroup:"):
            parts = line.strip().split(" ")
            orthogroup_id = parts[1]
            continue
        if line.startswith("Node:"):
            # 使用\t分隔行
            parts = line.strip().split("\t")
            node_id = parts[0].split(";")[0].split(" ")[1]
            if len(parts) > 1:
                tree = Tree(parts[1], format=1)

                all_leafs = tree.get_leaves()
               
                # 获取所有叶子节点
                all_species = set()
                for leaf in all_leafs:
                    if "|" in leaf.name:
                        species, gene = leaf.name.split("|")
                        all_species.add(species) if species != 'dogshark' else None

                if len(all_species) < 6: continue
                if not {"bichir", "lungfish"} & set(all_species): continue
                output.write(orthogroup_id + "_" + node_id + "\t" + line.strip() + "\n")
                for leaf in all_leafs:
                    if "|" in leaf.name:
                        species, gene = leaf.name.split("|")
                        data[species].add(gene)
                # break

# 输出每个物种的基因到单独的文件
for species, genes in data.items():
    output_file = os.path.join(output_dir, f"{species}.genes")
    with open(output_file, "w") as out_file:
        for gene in sorted(genes):
            out_file.write(f"{gene}\n")

# print("处理完成。输出文件保存在 'all_genes_processed' 目录中。")
