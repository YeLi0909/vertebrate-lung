#! /usr/bin/python
import os, re
from collections import defaultdict
from itertools import combinations

def read_gene_list(gene_lst):
    gene_group_mapping = defaultdict(str)
    group_gene_info = defaultdict(lambda: {"gene_name": None, "hvg_tag": set(), "species": defaultdict(set)})
    species_set = set()

    with open(gene_lst,'r') as file:
        file.readline()

        for line in file:
            parts = re.split("\t",line.strip())
            group_id, tag = parts[0], parts[1].upper()
            gene_name_counter = defaultdict(int)

            for i in range(2,len(parts)):
                if not parts[i]: continue

                for basic_element in parts[i].split(','):
                    # basic_element = basic_element.rstrip('*').lower()
                    basic_element = basic_element.rstrip('*')
                    species_id, gene_id = basic_element.split('|')
                    if species_id == "dogshark": continue
                    species_set.add(species_id)
                    raw_name = gene_id.upper()
                    gene_name_counter[raw_name] += 1
                    gene_group_mapping[basic_element] = group_id
                    group_gene_info[group_id]["species"][species_id].add(gene_id)
                    # if gene_id.upper() == "SFTPB":
                    #     print(basic_element)

            real_name = max(gene_name_counter, key = gene_name_counter.get)
            group_gene_info[group_id]["gene_name"] = real_name
            group_gene_info[group_id]["hvg_tag"].add(tag)

    return gene_group_mapping, group_gene_info, species_set

def read_cell_details(data_dir, species_set):
    cell_detail = defaultdict(lambda: defaultdict(str))
    cell_type_count = defaultdict(lambda: defaultdict(int))

    for species in species_set:
        file_name = f"{species}_meta.tsv"
        file_path = os.path.join(data_dir, file_name)

        if not os.path.exists(file_path):
            print(f"Warning: File {file_path} does not exist. Skipping.")
            continue

        with open(file_path, 'r') as file:
            next(file, None)

            for line in file:
                if re.search(r'Ciliated cell', line):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    print(f"Warning: Invalid line in {file_name}: {line.strip()}")
                    continue

                cell_id = parts[0]
                cell_type = parts[-1].upper()
                cell_id = re.sub(r'-', '.', cell_id) # for very odd and inconsistent format
                cell_detail[species][cell_id] = cell_type
                cell_type_count[cell_type][species] += 1

    return cell_detail, cell_type_count

def read_exp_data(data_dir, species_set, gene_group_mapping, cell_detail):
    gene_exp = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

    for species in species_set:
        file_name = f"{species}_exp.tsv"
        file_path = os.path.join(data_dir, file_name)

        if not os.path.exists(file_path):
            print(f"Warning: File {file_path} does not exist. Skipping.")
            continue

        with open(file_path, 'r') as file:
            header = file.readline().strip().split('\t')
            cell_ids = header[1:]

            for line in file:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    print(f"Warning: Invalid line in {file_name}: {line.strip()}")
                    continue

                # gene_id = parts[1].lower()
                gene_id = parts[1]
                # if gene_id == "bircir021097":
                #     print(species, gene_id)
                expression_levels = parts[2:]

                raw_id = species + "|" + gene_id
                if raw_id not in gene_group_mapping:
                    continue
                for cell_id, exp_level in zip(cell_ids, expression_levels):
                    try:
                        exp_value = float(exp_level)
                        if exp_value > 0:
                            cell_id = re.sub(r'-', '.', cell_id)
                            if cell_id not in cell_detail[species]: continue
                            cell_type = cell_detail[species][cell_id]
                            if len(cell_type) == 0:
                                print(species, cell_id)
                                break
                            gene_exp[species][gene_id][cell_type].add(cell_id)
                    except ValueError:
                        print(f"Warning: Invalid expression value in {file_name} for gene {gene_id}, cell {cell_id}: {exp_level}")

    return gene_exp

def process_gene_pair(group1, group2, group_gene_info, species_set, gene_exp):
    gene_name1 = group_gene_info[group1]["gene_name"]
    gene_name2 = group_gene_info[group2]["gene_name"]
    comb_name = f"{gene_name1}_{gene_name2}|{group1}_{group2}"

    all_tags = group_gene_info[group1]["hvg_tag"] | group_gene_info[group2]["hvg_tag"]
    result = defaultdict(lambda: defaultdict(lambda: "NA"))

    for hvg_tag in all_tags:
        for species in species_set:
            genes_for_group1 = group_gene_info[group1]["species"].get(species, set())
            genes_for_group2 = group_gene_info[group2]["species"].get(species, set())

            if not genes_for_group1 or not genes_for_group2:
                result[hvg_tag][species] = "NA"
            else:
                max_count = 0
                for gene1 in genes_for_group1:
                    for gene2 in genes_for_group2:
                        if gene1 in gene_exp[species] and gene2 in gene_exp[species]:
                            coexp_cells = gene_exp[species][gene1][hvg_tag] & gene_exp[species][gene2][hvg_tag]
                            count = len(coexp_cells)
                            max_count = max(max_count, count)
                
                result[hvg_tag][species] = max_count

    return comb_name, result

def analysis_data(group_gene_info, species_set, gene_exp, cell_type_count, out_file):
    exp_cell_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: "NA")))

    for group1, group2 in combinations(sorted(group_gene_info), 2):
        comb_name, result = process_gene_pair(group1, group2, group_gene_info, species_set, gene_exp)
        for hvg_tag in result:
            for species in result[hvg_tag]:
                exp_cell_counts[hvg_tag][comb_name][species] = result[hvg_tag][species]

    with open(out_file, 'w') as output:
        header = ["combination", "hvg_tag"] + sorted(species_set)
        output.write("\t".join(header) + "\n")

        for hvg_tag in sorted(exp_cell_counts.keys()):
            for comb_name in sorted(exp_cell_counts[hvg_tag].keys()):
                content = [comb_name, hvg_tag]

                for species in sorted(species_set):
                    exp_count = exp_cell_counts[hvg_tag][comb_name][species]
                    if exp_count == "NA":
                        content.append("NA")
                    else:
                        all_count = cell_type_count[hvg_tag][species]
                        ratio = exp_count / all_count if all_count != 0 else 0
                        content.append(f"{ratio:.4f}")
                
                output.write("\t".join(map(str, content)) + "\n")

def output_gene_exp(gene_exp, output_file):
    with open(output_file, 'w') as f:
        f.write("species\tgene_id\tcell_type\tcell_count\n")
        for species in gene_exp:
            for gene_id in gene_exp[species]:
                for cell_type in gene_exp[species][gene_id]:
                    cell_count = len(gene_exp[species][gene_id][cell_type])
                    f.write(f"{species}\t{gene_id}\t{cell_type}\t{cell_count}\n")

def output_cell_type_count(cell_type_count, output_file):
    with open(output_file, 'w') as f:
        f.write("cell_type\tspecies\tcount\n")
        for cell_type in cell_type_count:
            for species in cell_type_count[cell_type]:
                count = cell_type_count[cell_type][species]
                f.write(f"{cell_type}\t{species}\t{count}\n")

def main():
    data_dir = 'data_lung'
    gene_lst = 'bony_conserved_lung_genes.txt'
    out_file = 'all_lung.ratio'

    gene_group_mapping, group_gene_info, species_set = read_gene_list(gene_lst)
    print("gene list loaded")
    cell_detail, cell_type_count = read_cell_details(data_dir, species_set)
    print("meta data loaded")
    gene_exp = read_exp_data(data_dir, species_set, gene_group_mapping, cell_detail)
    print("expression data loaded")

    gene_exp_out_file = 'gene_exp_output.tsv'
    output_gene_exp(gene_exp, gene_exp_out_file)
    print(f"gene_exp data written to {gene_exp_out_file}")

    cell_type_count_out_file = 'cell_type_count_output.tsv'
    output_cell_type_count(cell_type_count, cell_type_count_out_file)
    print(f"cell_type_count data written to {cell_type_count_out_file}")

    analysis_data(group_gene_info, species_set, gene_exp, cell_type_count, out_file)

if __name__ == '__main__':
    main()