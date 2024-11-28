from collections import Counter, defaultdict
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
import networkx as nx
import matplotlib.pyplot as plt

def read_file(filename):
    gene_counts = Counter()
    gene_cne_map = defaultdict(set)
    all_cnes = set()
    with open(filename, 'r') as file:
        next(file)  # Skip header
        for line in file:
            gene, cne = line.split()[:2]
            gene_counts[gene] += 1
            gene_cne_map[gene].add(cne)
            all_cnes.add(cne)
    return gene_counts, gene_cne_map, all_cnes

def fisher_test(foreground_count, background_count, total_foreground, total_background):
    table = [
        [foreground_count, background_count],
        [total_foreground - foreground_count, total_background - background_count]
    ]
    odds_ratio, p_value = stats.fisher_exact(table, alternative='greater')
    return p_value, odds_ratio, table

# Read input files
background, bg_gene_cne_map, bg_all_cnes = read_file('shared.cnes.2m.filter.distance')
foreground, fg_gene_cne_map, fg_all_cnes = read_file('active_cne.distance')

all_genes = set(background.keys()) | set(foreground.keys())

total_background = len(bg_all_cnes)
total_foreground = len(fg_all_cnes)

# Perform Fisher's exact test for each gene
results = []
for gene in all_genes:
    bg_count = len(bg_gene_cne_map.get(gene, set()))
    fg_count = len(fg_gene_cne_map.get(gene, set()))
    p_value, odds_ratio, table = fisher_test(fg_count, bg_count, total_foreground, total_background)
    results.append((gene, bg_count, fg_count, p_value, odds_ratio, table))

# Sort results by p-value
results.sort(key=lambda x: x[3])

# Calculate FDR
_, fdr = fdrcorrection([r[3] for r in results])

# Add FDR to results
results_with_fdr = [r + (f,) for r, f in zip(results, fdr)]

# Save results to file
with open('enrichment_results.tsv', 'w') as outfile:
    outfile.write("Gene\tBackground_Count\tForeground_Count\tP_value\tOdds_Ratio\tFDR\tForeground_In\tBackground_In\tForeground_Out\tBackground_Out\n")
    for result in results_with_fdr:
        gene, bg_count, fg_count, p_value, odds_ratio, table, fdr = result
        fg_in, bg_in = table[0]
        fg_out, bg_out = table[1]
        outfile.write(f"{gene}\t{bg_count}\t{fg_count}\t{p_value}\t{odds_ratio}\t{fdr}\t{fg_in}\t{bg_in}\t{fg_out}\t{bg_out}\n")

print("Enrichment analysis complete. Results saved to 'enrichment_results.tsv'.")

# Get significant genes (FDR < 0.05)
significant_genes = [result[0] for result in results_with_fdr if result[6] < 0.05]
with open("significant_genes.txt", "w") as oh:
    for gene in significant_genes:
        oh.write(gene + "\n")

# Create a graph of significant genes sharing CNEs
G = nx.Graph()
for gene1 in significant_genes:
    for gene2 in significant_genes:
        if gene1 < gene2:  # Avoid duplicate edges
            shared_cnes = fg_gene_cne_map[gene1].intersection(fg_gene_cne_map[gene2])
            if shared_cnes:
                G.add_edge(gene1, gene2, weight=len(shared_cnes))

# Remove isolated nodes
G.remove_nodes_from(list(nx.isolates(G)))

# Draw the graph
plt.figure(figsize=(20, 20))
pos = nx.spring_layout(G, k=0.5, iterations=50)
nx.draw(G, pos, node_color='lightblue', node_size=500, font_size=8, with_labels=True)
edge_weights = nx.get_edge_attributes(G, 'weight')
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_weights, font_size=6)

# Save the graph as PDF
plt.savefig('gene_cne_network.pdf', format='pdf', dpi=300, bbox_inches='tight')
print("Gene-CNE network graph saved as 'gene_cne_network.pdf'.")

# Save edge information to a text file
with open('gene_cne_network_edges.txt', 'w') as f:
    f.write("Gene1\tGene2\tShared_CNEs\n")
    for edge in G.edges(data=True):
        f.write(f"{edge[0]}\t{edge[1]}\t{edge[2]['weight']}\n")

print("Edge information saved to 'gene_cne_network_edges.txt'.")