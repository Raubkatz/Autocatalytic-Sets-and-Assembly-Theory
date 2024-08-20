"""
Author: Dr. techn. Sebastian Raubitzek MSc. BSc.; SBA Research, Complexity and Resilience Research Group
This Python script demonstrates the use of the AutocatalyticSet class to simulate and analyze autocatalytic reaction networks, as part of a simplified model discussed in the paper "Autocatalytic Sets and Assembly Theory: A Toy Model Perspective" by Sebastian Raubitzek et al.

### Problem Description:
We are modeling the formation of a complex final product (e.g., 'ABCDEFGH') through a network of reactions that are potentially autocatalytic. The goal is to understand how autocatalysis influences the assembly of such products, particularly in environments with varying levels of catalytic activity.

### Mathematical Model:
The core concept revolves around the idea of autocatalytic sets, where the reactions in the network are catalyzed by products of other reactions within the same network. This creates a self-sustaining system that can potentially explain the emergence of complex biochemical structures.

### Numerical Method:
The script builds a reaction network for a predefined final product using a probabilistic approach to assign catalysts:
- **Catalyst Assignment**: Catalysts are assigned to reactions based on a specified probability, simulating different environmental conditions.
- **Catalyst Strategies**: Three strategies for catalyst assignment are implemented:
  - Random assignment without reuse.
  - Weighted random assignment without reuse.
  - Weighted random assignment with reuse (most realistic for complex systems).

### Parameters Used:
- **Final Product**: 'ABCDEFGH' - A sequence representing the complex object to be assembled.
- **Catalyst Probability**: 0.75 - The likelihood that a given reaction will have an associated catalyst.
- **Catalyst Strategy**: 2 (Weighted random assignment with reuse).

### Expected Outcome:
The script will generate a network of reactions leading to the final product, assess whether the network and its subnetworks are autocatalytic, and visualize the structure of the network. The outcome will provide insights into the role of autocatalysis in complex system formation.

### Implementation:
The script initializes the AutocatalyticSet with the specified parameters, builds the reaction network, checks for autocatalytic properties, and visualizes the network and its subnetworks.
"""

import os
from class_AutocatalyticSet import *

# PARAMETERS ###########################################################################################################
# The final product of the autocatalytic set, represented as a string of characters.
# Each character can be thought of as a different element or molecule involved in the reactions.
#final_product = 'ABCDEFGαβΓγΔδ'

final_product = 'ABCDEFGH'


# The probability of assigning a catalyst to each reaction, represented as a float between 0 and 1.
# For example, a value of 0.5 means there is a 50% chance that any given reaction will have a catalyst.
catalyst_probability = 0.75

# The strategy used for assigning catalysts to reactions, represented as an integer:
# 0 - Random assignment without reuse: Catalysts are assigned randomly to reactions without reusing the same catalyst.
# 1 - Weighted random assignment without reuse: Catalysts are assigned based on their probability, weighted by the length
#     of the catalyst string, without allowing reuse.
# 2 - Weighted random assignment with reuse: Similar to option 1, but catalysts can be reused across different reactions.
switch_cat = 2

# The directory path where the results (like the reactions file) will be saved.
# Ensure this directory exists or the script has permissions to create it.
save_path = "autocatalytic_demonstration"
########################################################################################################################

# Check if the save_path directory exists, and create it if it does not
if not os.path.exists(save_path):
    os.makedirs(save_path)

# Initialize the AutocatalyticSet with the final product and the specified parameters
autocatalytic_set = AutocatalyticSet(final_product, catalyst_probability=catalyst_probability, switch_cat=switch_cat)

# Build the reactions for the autocatalytic set
autocatalytic_set.build_reactions()

# Save the reactions to a file
autocatalytic_set.save_reactions_to_file(f'./{save_path}/autocatalytic_reactions.txt')

# Display the reactions in a readable format
print(f"The whole network for {final_product}:")
autocatalytic_set.display_reactions()

# Display the reactions in a readable format
print(f"Is the whole set autocatalytic? {autocatalytic_set.is_autocatalytic_subnetwork(autocatalytic_set.reactions_df)}")


# Extract and analyze subnetworks within the autocatalytic set
autocatalytic_set.extract_subnetworks()

# Draw the whole network
autocatalytic_set.draw_network()

print("All non-empty subnetworks:")
# Iterate through the extracted subnetworks
for product, df in autocatalytic_set.subnetwork_dict.items():
    print(f"Subnetwork for {product}:")
    print(df)

    # Check if the subnetwork is autocatalytic
    is_autocatalytic = autocatalytic_set.is_autocatalytic_subnetwork(df)
    print(f"Is subnetwork for {product} autocatalytic? {is_autocatalytic}")

    # Draw the subnetwork
    autocatalytic_set.draw_network_sub(reactions_df=df, key=product)

# Check if each subnetwork is autocatalytic and display the result
for product, subnetwork_df in autocatalytic_set.subnetwork_dict.items():
    is_autocatalytic = autocatalytic_set.is_autocatalytic_subnetwork(subnetwork_df)
    print(f"Is subnetwork for {product} autocatalytic? {is_autocatalytic}")
