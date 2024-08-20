"""
Author: Dr. techn. Sebastian Raubitzek MSc. BSc.; SBA Research, Complexity and Resilience Research Group
The AssemblyAnalysis class is designed to perform extensive analysis on autocatalytic sets, focusing on how different configurations of catalyst assignment, probability, and network structures influence the assembly process of complex products. This class builds upon the foundation provided by the AutocatalyticSet class and is a key component of the study presented in "Autocatalytic Sets and Assembly Theory: A Toy Model Perspective" by Sebastian Raubitzek et al.

### Purpose:
The AssemblyAnalysis class systematically generates multiple instances of autocatalytic sets and analyzes their assembly indices. It allows for the exploration of how different factors, such as catalyst probability and autocatalytic subnetwork presence, affect the depth and complexity of reaction networks.

### Key Features:
- **Multiple Instance Generation**: The class generates numerous autocatalytic sets based on the specified final product, allowing for robust statistical analysis across different configurations.
- **Catalyst Assignment and Autocatalytic Subnetwork Analysis**: By leveraging different catalyst assignment strategies and considering autocatalytic subnetworks, the class provides insights into how these factors impact the assembly process.
- **Assembly Index Calculation**: The class calculates the assembly index for each generated network, which represents the maximum depth of the network, adjusted based on the presence of autocatalytic subnetworks and the overall autocatalytic nature of the set.
- **Statistical Analysis**: It aggregates data across multiple generated networks to produce meaningful statistics, which can be used to identify trends and patterns in the assembly process.

### Attributes:
- **final_product (str)**: The final product that the reaction networks aim to synthesize.
- **n (int)**: The number of different autocatalytic sets to generate and analyze.
- **catalyst_probability (float)**: The probability of assigning a catalyst to each reaction.
- **switch_cat (int)**: Determines the strategy used for catalyst assignment (0: random, 1: weighted random without reuse, 2: weighted random with reuse).
- **switch_aut_prob (bool/str)**: Determines how the presence of autocatalytic subnetworks influences the assembly index (True: linear, "square": exponential).
- **autocatalytic_subnetwork_counts (list)**: Stores the count of autocatalytic subnetworks for each generated set.
- **subnetwork_counts (list)**: Stores the total count of subnetworks for each generated set.
- **max_depths (list)**: Stores the maximum depth (complexity) of each generated network.

### Methods:
- **find_assembly_index()**: Generates multiple autocatalytic sets, processes them to calculate the assembly index for each set, and returns the minimum assembly index across all sets.
- **process_autocatalytic_set(autocatalytic_set)**: Processes a single AutocatalyticSet instance, extracting and analyzing its subnetworks, and determining the set's maximum depth and subnetwork counts.
- **get_max_layer(subnetwork_df)**: Determines the maximum layer (depth) within a given subnetwork.
- **find_assembly_list(return_all_counts=False)**: Generates multiple autocatalytic sets, calculates the assembly indices, and returns a list of these indices. Optionally returns detailed counts of subnetworks and depths.

### Connection to AutocatalyticSet Class:
The AssemblyAnalysis class is closely connected to the AutocatalyticSet class. Each instance of an autocatalytic set, generated and analyzed by AssemblyAnalysis, is an instance of the AutocatalyticSet class. The AssemblyAnalysis class relies on AutocatalyticSet to:
- **Build Reaction Networks**: AutocatalyticSet constructs the reaction networks that AssemblyAnalysis evaluates.
- **Extract and Analyze Subnetworks**: AutocatalyticSet provides the tools for identifying and analyzing subnetworks within the broader reaction network.
- **Determine Autocatalytic Properties**: AssemblyAnalysis uses the methods in AutocatalyticSet to assess whether the entire network or specific subnetworks are autocatalytic.

The role of AssemblyAnalysis is to scale up the analysis by systematically generating multiple instances of these autocatalytic sets, enabling a comprehensive statistical evaluation of how different parameters influence the assembly process. This connection allows for a deeper exploration of the emergent properties of autocatalytic systems, as discussed in the paper.

### Context and Usage:
The AssemblyAnalysis class is intended for use in large-scale analyses where multiple instances of autocatalytic sets need to be generated and compared. It is particularly valuable in exploring how variations in catalyst assignment strategies, probabilities, and autocatalytic subnetwork handling affect the complexity and depth of the resulting networks. By aggregating data across thousands of generated sets, the class enables researchers to draw significant conclusions about the underlying mechanisms of autocatalytic assembly.

This class plays a crucial role in the broader project, as it builds upon the fundamental principles of autocatalysis and assembly theory discussed in the paper, providing a practical tool for conducting extensive computational experiments.
"""


import numpy as np
from class_AutocatalyticSet import AutocatalyticSet
import pandas as pd

class AssemblyAnalysis:
    def __init__(self, final_product, n=10000, catalyst_probability=1.0, switch_cat=2, switch_aut_prob=False, print_out=False):
        self.final_product = final_product
        self.n = n
        self.catalyst_probability = catalyst_probability
        self.switch_cat = switch_cat
        self.switch_aut_prob = switch_aut_prob
        self.autocatalytic_subnetwork_counts = []
        self.subnetwork_counts = []
        self.max_depths = []

        self.print_out = print_out

    def find_assembly_index(self):
        assembly_indices = []

        for _ in range(self.n):
            autocatalytic_set = AutocatalyticSet(self.final_product, self.catalyst_probability,
                                                 switch_cat=self.switch_cat, print_out=self.print_out)
            autocatalytic_set.build_reactions()
            subnetwork_count, max_depth = self.process_autocatalytic_set(autocatalytic_set)

            self.autocatalytic_subnetwork_counts.append(subnetwork_count)

            if self.switch_aut_prob:
                assembly_indices.extend([max_depth] * subnetwork_count)
            else:
                assembly_indices.append(max_depth)

        if not assembly_indices:
            return float('inf')

        return min(assembly_indices)

    def process_autocatalytic_set(self, autocatalytic_set):
        subnetwork_count = 0 #counts all subnetworks
        autocatalytic_subnetwork_count = 0 #counts all autocatalytic subnetworks
        max_depth = 0
        autocatalytic_set.extract_subnetworks()

        #get max layer
        for reaction in autocatalytic_set.reactions:
            layer = reaction['Layer']  # This is already an integer
            if layer > max_depth:
                max_depth = layer

        #layer only counts the layer on the product level, but we count layer on the Constituent level, which is just +1
        max_depth = max_depth + 1

        for product, subnetwork_df in autocatalytic_set.subnetwork_dict.items(): # going thorugh all the subnetworks and assess if they are autocatalytic
            #print(f'Product: {product}')
            #print(f'subnetwork_df: {subnetwork_df}')
            #print(subnetwork_df)
            #print(f'Is this subnetwork autocatalytic? {autocatalytic_set.is_autocatalytic_subnetwork(subnetwork_df)}')
            if autocatalytic_set.is_autocatalytic_subnetwork(subnetwork_df):
                autocatalytic_subnetwork_count += 1
            subnetwork_count += 1

        return subnetwork_count, autocatalytic_subnetwork_count, max_depth

    def get_max_layer(self, subnetwork_df):
        max_layer = 0
        for _, reaction in subnetwork_df.iterrows():
            layer = reaction['Layer']  # This is already an integer
            if layer > max_layer:
                max_layer = layer
        return max_layer + 1 # because the layer is always the layer of the final product, it's just a matter a counting

    def find_assembly_list(self, return_all_counts=False):
        assembly_indices = []

        for _ in range(self.n):
            autocatalytic_set = AutocatalyticSet(self.final_product, self.catalyst_probability,
                                                 switch_cat=self.switch_cat, print_out=self.print_out)
            autocatalytic_set.build_reactions()
            subnetwork_count, autocatalytic_subnetwork_count, max_depth = self.process_autocatalytic_set(autocatalytic_set)

            whole_set_autocatalytic = 0

            if autocatalytic_set.is_autocatalytic_subnetwork(autocatalytic_set.reactions_df):
                whole_set_autocatalytic = 1

            if self.print_out:
                print(
                    f'Is the whole set is autocatalytic? {autocatalytic_set.is_autocatalytic_subnetwork(autocatalytic_set.reactions_df)}')
                print(autocatalytic_set.reactions_df)
                print(
                    f'There are {subnetwork_count} subnetworks, and {autocatalytic_subnetwork_count} of them are autocatalytic.')
                print(f'The max_depth of the current network is {max_depth} ')
                print(f'The Multiplier, i.e. 1 + (is autocatlayitc) + (#autocatalytic subsets) is {(1 + whole_set_autocatalytic + autocatalytic_subnetwork_count)}')
                print('##########################################################################')
                print('##########################################################################')
                print('##########################################################################')
                print('##########################################################################')

            self.autocatalytic_subnetwork_counts.append(autocatalytic_subnetwork_count)
            self.subnetwork_counts.append(subnetwork_count)
            self.max_depths.append(max_depth)

            if self.switch_aut_prob == True:
                assembly_indices.extend([max_depth] * (1 + whole_set_autocatalytic + autocatalytic_subnetwork_count)) #the maximal depth is multiplied with 1+ the number of autocatalytic subsets + 1 if the whole set is autocatalytic. This means that for each network we take into account the max depth at least once; however if there are sveral autocatalytic subsets amogn the network this is increased; i.e. one additional copy of the assembly index for each autocatalytic subset
            elif self.switch_aut_prob == "square":
                assembly_indices.extend([max_depth] * pow((1 + whole_set_autocatalytic + autocatalytic_subnetwork_count),2)) #the maximal depth is multiplied with 1+ the number of autocatalytic subsets + 1 if the whole set is autocatalytic. This means that for each network we take into account the max depth at least once; however if there are sveral autocatalytic subsets amogn the network this is increased; i.e. one additional copy of the assembly index for each autocatalytic subset
            else:
                assembly_indices.append(max_depth)
        if return_all_counts:
            return assembly_indices, self.subnetwork_counts, self.autocatalytic_subnetwork_counts, self.max_depths

        return assembly_indices
