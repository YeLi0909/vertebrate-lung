import re
from collections import defaultdict
from scipy import stats
import numpy as np

def read_file(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 4:
                data.append(parts)
    return data

def count_entries(data, group, condition, value_condition=None):
    count = 0
    for entry in data:
        if entry[1] == group and (condition(entry[2]) if callable(condition) else entry[2] in condition):
            if value_condition is None or value_condition(float(entry[3])):
                count += 1
    return count

def calculate_cramers_v(confusion_matrix):
    chi2 = stats.chi2_contingency(confusion_matrix)[0]
    n = np.sum(confusion_matrix)
    min_dim = min(confusion_matrix.shape) - 1
    return np.sqrt(chi2 / (n * min_dim))

# Read the data
data = read_file('ultra_cne.lst.uniq')

# Define conditions
# active_conditions = ['strong', 'enhancer']
active_conditions = ['strong', 'enhancer']
is_active = lambda x: x in active_conditions
is_passive = lambda x: x not in active_conditions
value_condition = lambda x: x < 3

# Count entries for each category
jawed_active = count_entries(data, 'jawed', is_active)
jawed_active_cond = count_entries(data, 'jawed', is_active, value_condition)
jawed_passive = count_entries(data, 'jawed', is_passive)
jawed_passive_cond = count_entries(data, 'jawed', is_passive, value_condition)

bony_active = count_entries(data, 'bony', is_active)
bony_active_cond = count_entries(data, 'bony', is_active, value_condition)
bony_passive = count_entries(data, 'bony', is_passive)
bony_passive_cond = count_entries(data, 'bony', is_passive, value_condition)

# Create confusion matrices
jawed_matrix = np.array([[jawed_active_cond, jawed_active - jawed_active_cond], 
                         [jawed_passive_cond, jawed_passive - jawed_passive_cond]])

bony_matrix = np.array([[bony_active_cond, bony_active - bony_active_cond], 
                        [bony_passive_cond, bony_passive - bony_passive_cond]])

# Perform Fisher's exact test
jawed_oddsratio, jawed_pvalue = stats.fisher_exact(jawed_matrix)
bony_oddsratio, bony_pvalue = stats.fisher_exact(bony_matrix)

# Calculate Cramer's V
jawed_cramers_v = calculate_cramers_v(jawed_matrix)
bony_cramers_v = calculate_cramers_v(bony_matrix)

# Print results
print("Jawed Fish Results:")
print(f"Confusion Matrix:\n{jawed_matrix}")
print(f"Fisher's Exact Test p-value: {jawed_pvalue}")
print(f"Cramer's V: {jawed_cramers_v}")

print("\nBony Fish Results:")
print(f"Confusion Matrix:\n{bony_matrix}")
print(f"Fisher's Exact Test p-value: {bony_pvalue}")
print(f"Cramer's V: {bony_cramers_v}")
