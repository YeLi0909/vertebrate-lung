#! /usr/bin/python
import os, re
from collections import defaultdict
from itertools import combinations

def read_gene_list(gene_lst, target_species):
    gene_group_mapping = defaultdict(str)
    group_gene_info = defaultdict(lambda: {"gene_name": None, "hvg_tag": set(), "all_gene": set()})

    with open(gene_lst,'r') as file:
        file.readline()

        for line in file:
            parts = re.split("\t",line.strip())
            group_id, tag = parts[0], parts[1].upper()
            gene_name_counter = defaultdict(int)

            for i in range(2,len(parts)):
                if not parts[i]: continue

                for basic_element in parts[i].split(','):
                    basic_element = basic_element.rstrip('*').lower()
                    species_id, gene_id = basic_element.split('|')
                    raw_name = gene_id.upper()
                    gene_name_counter[raw_name] += 1
                    if species_id != target_species: continue
                    gene_group_mapping[basic_element] = group_id
                    group_gene_info[group_id]["all_gene"].add(gene_id)
                    # if gene_id.upper() == "SFTPB":
                    #     print(basic_element)

            real_name = max(gene_name_counter, key = gene_name_counter.get)
            group_gene_info[group_id]["gene_name"] = real_name
            group_gene_info[group_id]["hvg_tag"].add(tag)

    return gene_group_mapping, group_gene_info

def read_cell_details(data_dir, target_species):
    cell_detail = defaultdict(lambda: defaultdict(str))
    cell_type_count = defaultdict(int)

    file_name = f"{target_species}_meta.tsv"
    file_path = os.path.join(data_dir, file_name)

    if not os.path.exists(file_path):
        print(f"Warning: File {file_path} does not exist. Skipping.")
        exit()

    with open(file_path, 'r') as file:
        next(file)

        for line in file:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                print(f"Warning: Invalid line in {file_name}: {line.strip()}")
                continue

            cell_id = parts[0]
            cell_type = parts[-1].upper()
            # cell_id = re.sub(r'-', '.', cell_id) # for very odd and inconsistent format
            cell_detail[cell_id] = cell_type
            cell_type_count[cell_type] += 1

    return cell_detail, cell_type_count

def read_exp_data(data_dir, target_species, gene_group_mapping, cell_detail):
    gene_exp = defaultdict(lambda: defaultdict(set))

    file_name = f"{target_species}_exp.tsv"
    file_path = os.path.join(data_dir, file_name)

    if not os.path.exists(file_path):
        print(f"Warning: File {file_path} does not exist. Skipping.")
        exit()

    with open(file_path, 'r') as file:
        header = file.readline().strip().split('\t')
        cell_ids = header[1:]

        for line in file:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                print(f"Warning: Invalid line in {file_name}: {line.strip()}")
                continue

            gene_id = parts[1].lower()
            expression_levels = parts[2:]

            raw_id = target_species + "|" + gene_id
            if raw_id not in gene_group_mapping:
                continue
            for cell_id, exp_level in zip(cell_ids, expression_levels):
                try:
                    exp_value = float(exp_level)
                    if exp_value > 0:
                        cell_type = cell_detail[cell_id]
                        if len(cell_type) == 0:
                            print(target_species, cell_id)
                            break
                        gene_exp[gene_id][cell_type].add(cell_id)
                except ValueError:
                    print(f"Warning: Invalid expression value in {file_name} for gene {gene_id}, cell {cell_id}: {exp_level}")

    return gene_exp

def process_gene_pair(group1, group2, group_gene_info, gene_exp, all_cell_type):
    gene_name1 = group_gene_info[group1]["gene_name"]
    gene_name2 = group_gene_info[group2]["gene_name"]
    comb_name = f"{gene_name1}_{gene_name2}|{group1}_{group2}"

    result = defaultdict(lambda: "NA")

    for cell_type in all_cell_type:
        genes_for_group1 = group_gene_info[group1].get("all_gene", set())
        genes_for_group2 = group_gene_info[group2].get("all_gene", set())

        if not genes_for_group1 or not genes_for_group2:
            result[cell_type] = "NA"
        else:
            max_count = 0
            for gene1 in genes_for_group1:
                for gene2 in genes_for_group2:
                    if gene1 in gene_exp and gene2 in gene_exp:
                        coexp_cells = gene_exp[gene1][cell_type] & gene_exp[gene2][cell_type]
                        count = len(coexp_cells)
                        max_count = max(max_count, count)
            
            result[cell_type] = max_count

    return comb_name, result

def analysis_data(group_gene_info, gene_exp, cell_type_count, out_file):
    exp_cell_counts = defaultdict(lambda: defaultdict(lambda: "NA"))

    all_cell_type = cell_type_count.keys()
    for group1, group2 in combinations(sorted(group_gene_info), 2):
        comb_name, result = process_gene_pair(group1, group2, group_gene_info, gene_exp, all_cell_type)
        for cell_type in result:
            exp_cell_counts[comb_name][cell_type] = result[cell_type]

    with open(out_file, 'w') as output:
        header = ["combination"] + sorted(all_cell_type)
        output.write("\t".join(header) + "\n")

        for comb_name in sorted(exp_cell_counts.keys()):
            content = [comb_name]

            for cell_type in sorted(all_cell_type):
                exp_count = exp_cell_counts[comb_name][cell_type]
                if exp_count == "NA":
                    content.append("NA")
                else:
                    all_count = cell_type_count[cell_type]
                    # print(comb_name, cell_type, exp_count, all_count)
                    ratio = exp_count / all_count if all_count != 0 else 0
                    content.append(f"{ratio:.4f}")
                
            output.write("\t".join(map(str, content)) + "\n")

def output_gene_exp(gene_exp, output_file):
    with open(output_file, 'w') as f:
        f.write("gene_id\tcell_type\tcell_count\n")
        for gene_id in gene_exp:
            for cell_type in gene_exp[gene_id]:
                cell_count = len(gene_exp[gene_id][cell_type])
                f.write(f"{gene_id}\t{cell_type}\t{cell_count}\n")

def output_cell_type_count(cell_type_count, output_file):
    with open(output_file, 'w') as f:
        f.write("cell_type\tcount\n")
        for cell_type in cell_type_count:
            count = cell_type_count[cell_type]
            f.write(f"{cell_type}\t{count}\n")

def main():
    data_dir = 'data_lung'
    gene_lst = 'bony_conserved_lung_genes.txt'
    target_species = "dogshark"
    out_file = f'all_{target_species}.ratio'

    gene_group_mapping, group_gene_info = read_gene_list(gene_lst, target_species)
    print("gene list loaded")
    cell_detail, cell_type_count = read_cell_details(data_dir, target_species)
    print("meta data loaded")
    gene_exp = read_exp_data(data_dir, target_species, gene_group_mapping, cell_detail)
    print("expression data loaded")

    gene_exp_out_file = f'{target_species}_gene_exp_output.tsv'
    output_gene_exp(gene_exp, gene_exp_out_file)
    print(f"gene_exp data written to {gene_exp_out_file}")

    cell_type_count_out_file = f'{target_species}_cell_type_count_output.tsv'
    output_cell_type_count(cell_type_count, cell_type_count_out_file)
    print(f"cell_type_count data written to {cell_type_count_out_file}")

    analysis_data(group_gene_info, gene_exp, cell_type_count, out_file)

if __name__ == '__main__':
    main()