"""
Author: Dr. techn. Sebastian Raubitzek MSc. BSc.; SBA Research, Complexity and Resilience Research Group
This Python script demonstrates the use of the AssemblyAnalysis class to generate and analyze multiple autocatalytic sets, focusing on understanding how assembly and autocatalytic properties manifest across different reaction networks.

### Problem Description:
We are simulating and analyzing the formation of multiple autocatalytic sets, which represent networks of reactions where products catalyze other reactions within the network. The objective is to observe the variability in the assembly index and to assess the presence of autocatalytic subnetworks.

### Model Overview:
The AssemblyAnalysis class is used to generate a specified number of reaction networks, each leading to a final product. The final product ('AB') is assembled through a series of reactions that may or may not be autocatalytic, depending on the environmental parameters, such as the probability of assigning catalysts to reactions.

### Parameters:
- **Final Product**: 'AB' - A sequence representing the target object to be assembled through autocatalytic reactions.
- **Number of Sets (n)**: 10 - The number of autocatalytic sets (reaction networks) to generate and analyze.
- **Catalyst Probability**: 0.8 - The probability that a reaction in the network will have a catalyst assigned. A higher probability indicates a more catalyst-rich environment, which may affect the structure of the assembly.
- **Catalyst Strategy (switch_cat)**: 2 (Weighted random assignment with reuse) - Catalysts are assigned based on their weighted probability and can be reused across multiple reactions. This reflects more realistic biochemical or chemical systems where catalysts are often reused.
- **Autocatalytic Subnetwork Switch (switch_aut_prob)**: True - If set to True, the assembly index for each set is counted multiple times based on the number of autocatalytic subnetworks within that set. This provides insight into how the presence of autocatalytic subsets affects the overall assembly process.

### Expected Outcome:
The script generates multiple reaction networks and evaluates the assembly list, which quantifies the assembly indices of the networks. The script also checks and prints the statistics of the autocatalytic subnetworks, helping to understand the distribution and prevalence of autocatalysis in the generated networks.

### Implementation:
- The script initializes the AssemblyAnalysis class with the given parameters and generates multiple reaction networks.
- It calculates the assembly list (a measure of the complexity and depth of the assembly process) and prints basic statistics, such as the minimum, mean, maximum, and standard deviation of the assembly indices.
- Additionally, it checks for autocatalytic subnetworks, providing counts and mean values, which are critical in understanding how autocatalysis impacts the overall network structure.

This implementation allows for an exploration of how different catalyst assignment strategies and autocatalytic tendencies influence the complexity and organization of reaction networks.
"""

from class_AssemblyAnalysis import *
import numpy as np
import os

# PARAMETERS ###########################################################################################################
# The final product of the assembly analysis, represented as a string of characters.
# Each character can be thought of as a different element or molecule involved in the reactions.
#final_product = 'ABCDEFGH'

final_product = 'AB'


# The number of different autocatalytic sets to generate.
n = 10

# The probability of assigning a catalyst to each reaction, represented as a float between 0 and 1.
# For example, a value of 0.8 means there is an 80% chance that any given reaction will have a catalyst.
catalyst_probability = 0.8

# The strategy used for assigning catalysts to reactions, represented as an integer:
# 0 - Random assignment without reuse: Catalysts are assigned randomly to reactions without reusing the same catalyst.
# 1 - Weighted random assignment without reuse: Catalysts are assigned based on their probability, weighted by the length
#     of the catalyst string, without allowing reuse.
# 2 - Weighted random assignment with reuse: Similar to option 1, but catalysts can be reused across different reactions.
switch_cat = 2

# Boolean switch to determine if the assembly index should be counted multiple times based on the number of autocatalytic
# subnetworks. If True, the assembly index is counted multiple times; otherwise, it is counted once.
switch_aut_prob = True
########################################################################################################################

# Initialize the AssemblyAnalysis with the final product and the specified parameters
assembly_analysis = AssemblyAnalysis(final_product, n=n, catalyst_probability=catalyst_probability, switch_cat=switch_cat, switch_aut_prob=switch_aut_prob, print_out=True)

# Find the assembly list for the generated reaction networks
assembly_list = assembly_analysis.find_assembly_list()

# Check if the assembly list is not empty and print statistics
if assembly_list:
    print(
        f'Assembly List, Min: {np.min(np.array(assembly_list))} Mean: {np.mean(np.array(assembly_list))}, Max: {np.max(np.array(assembly_list))}, Std: {np.std(np.array(assembly_list))}, Length: {len(assembly_list)}')
else:
    print('Assembly List is empty.')

# Print the counts of autocatalytic subnetworks
print(f'Autocatalytic subnetwork counts: {assembly_analysis.autocatalytic_subnetwork_counts}')
# Print the mean of the autocatalytic subnetwork counts, handling the case where the list might be empty
print(
    f'Mean autocatalytic subnetwork count: {np.mean(assembly_analysis.autocatalytic_subnetwork_counts) if assembly_analysis.autocatalytic_subnetwork_counts else float("nan")}')
