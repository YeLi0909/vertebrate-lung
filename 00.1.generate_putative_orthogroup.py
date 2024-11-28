import os
import re
from ete3 import Tree
from collections import defaultdict

def parse_newick(file_path):
    with open(file_path, 'r') as file:
        newick_str = file.read().strip()
        newick_str = re.sub(r'(\w+)_(\w+\|)', r'\2', newick_str)
    return newick_str

def is_outgroup(parent, grandparent):
    # print(parent, grandparent)
    parent_node = species_tree.get_common_ancestor(*parent)
    
    for species in grandparent:
        if species not in parent:
            species_node = species_tree.search_nodes(name=species)[0]
            if parent_node not in species_node.get_ancestors():
                return True
    
    return False

def process_intermediate_nodes(tree):
    top_nodes = set()
    visited_nodes = set()
    trace_info = []
    initial_nodes = []

    def initialize_node(node):
        if not hasattr(node, 'genes'):
            if node.is_leaf():
                species = node.name.split('|')[0]
                node.add_features(
                    genes=set([node.name]),
                    species=set([species])
                )
            else:
                for child in node.children:
                    initialize_node(child)
                node.add_features(
                    genes=set.union(*[child.genes for child in node.children]),
                    species=set.union(*[child.species for child in node.children])
                )

    def check_condition(node, parent):
        species_increase = len(parent.species) - len(node.species)
        gene_increase = len(parent.genes) - len(node.genes)
        return species_increase > 0 and gene_increase - species_increase <= 3

    for node in tree.traverse("postorder"):
        initialize_node(node)

        if not node.is_leaf() and len(node.species) >= 4 and len(node.genes) / len(node.species) <= 2:
            initial_nodes.append(node)

    # Refine initial nodes
    refined_initial_nodes = []
    for node in initial_nodes:
        is_redundant = False
        for other_node in initial_nodes:
            if node != other_node and node in other_node.get_ancestors() and node.genes == other_node.genes:
                is_redundant = True
                break
        if not is_redundant:
            refined_initial_nodes.append(node)

    # Process refined initial nodes
    for seed_node in refined_initial_nodes:
        trace = [seed_node.name]
        while seed_node.up and seed_node.up.name not in visited_nodes:
            parent = seed_node.up
            initialize_node(parent)
            
            if check_condition(seed_node, parent):
                trace.append(parent.name)
                seed_node = parent
            else:
                if parent.up and parent.up.name not in visited_nodes:
                    grandparent = parent.up
                    initialize_node(grandparent)
                    if check_condition(parent, grandparent):
                        if is_outgroup(parent.species, grandparent.species):
                            trace.append(grandparent.name)
                            seed_node = grandparent
                        else:
                            break
                    else:
                        break
                else:
                    break
        
        if seed_node.name not in visited_nodes:
            top_nodes.add(seed_node)
            visited_nodes.add(seed_node.name)
            for descendant in seed_node.iter_descendants():
                visited_nodes.add(descendant.name)
            trace_info.append(f"From {trace[0]} to {trace[-1]}: {' -> '.join(trace)}")

    highest_level_nodes = set()
    for node in top_nodes:
        is_highest = True
        for other_node in top_nodes:
            if node != other_node and node in other_node.get_ancestors():
                is_highest = False
                break
        if is_highest:
            highest_level_nodes.add(node)

    return list(highest_level_nodes)

def format_node_info(node):
    subtree = node.copy()
    for n in subtree.traverse():
        n.dist = 0
        if not n.is_leaf() and not n.name:
            n.name = "Internal"
    node_name = node.name if node.name else "Internal"
    return f"Node: {node_name}; Genes: {len(node.genes)}; Species: {len(node.species)};\t{subtree.write(format=8)}"

newick = "(dogshark,(bichir,(lungfish,(frog,((lizard,chicken),((pig,human),(rat,mouse)))))));"
species_tree = Tree(newick, format=1)

dst_dir = "/data02/liye/project/analys_sr/28.orthfider/result/Results_Aug13/Resolved_Gene_Trees"
output_file = "orthogroup_summary.txt"
table_output_file = "species_genes_table.txt"

all_species = set()
table_data = []

with open(output_file, 'w') as output:
    for filename in os.listdir(dst_dir):
        if filename.endswith('_tree.txt'):
            ortho_id, _ = filename.split('_tree')
            tree_path = os.path.join(dst_dir, filename)
            
            try:
                newick = parse_newick(tree_path)
                tree = Tree(newick, format=1)
                
                output.write(f"\nOrthogroup: {ortho_id}\n")
                top_nodes = process_intermediate_nodes(tree)

                valid_nodes = [node for node in top_nodes if len(node.species) >= 3]
                if valid_nodes:
                    for i, node in enumerate(valid_nodes, 1):
                        output.write(format_node_info(node) + "\n")
                        
                        node_name = node.name if node.name else "Internal"
                        node_id = f"{ortho_id}_{node_name}"

                        species_genes = defaultdict(list)
                        for gene in node.genes:
                            species, gene_name = gene.split('|')
                            # species = species.split('_')[-1]
                            # species_genes[species].append(gene_name)
                            species_genes[species].append(gene)
                            all_species.add(species)
                        
                        table_data.append((node_id, species_genes))
            except Exception as e:
                output.write(f"Error processing {filename}: {str(e)}\n")
                raise

print(f"Analysis complete. Output file in {output_file}")

with open(table_output_file, 'w') as table_output:
    header = ['Orthogroup_node_id'] + sorted(all_species)
    table_output.write('\t'.join(header) + '\n')
    
    for node_id, species_genes in table_data:
        row = [node_id]
        for species in sorted(all_species):
            genes = ','.join(species_genes[species]) if species in species_genes else ''
            row.append(genes)
        table_output.write('\t'.join(row) + '\n')

print(f"Table output file created: {table_output_file}")