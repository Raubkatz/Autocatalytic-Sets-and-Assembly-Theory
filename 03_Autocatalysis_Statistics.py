"""
Author: Dr. techn. Sebastian Raubitzek MSc. BSc.
This Python script utilizes the AssemblyAnalysis class to generate and analyze multiple autocatalytic sets with varying parameters. The script is designed to explore how different combinations of factors influence the assembly process and the emergence of autocatalytic subnetworks.

### Problem Description:
We are conducting a series of experiments to analyze the formation of complex final products (e.g., sequences of up to 26 letters) through autocatalytic reaction networks. The goal is to understand how varying the catalyst probability, catalyst assignment strategy, and the consideration of autocatalytic subnetworks affect the assembly process.

### Parameters Used:
- **Run Name**: 'run_10000_26_letters_squared_missed' - A unique identifier for this experiment.
- **Final Products**: ['ABCDEFGHIJKLMNOPQRSTUVWXY', 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'] - A list of sequences representing the complex products to be assembled.
- **Catalyst Probabilities**: [0.0, 0.1, 0.2, ..., 1.0] - A range of probabilities for assigning catalysts to reactions.
- **Catalyst Assignment Strategies**: [0, 1, 2] - Different methods for assigning catalysts:
  - 0: Random assignment without reuse.
  - 1: Weighted random assignment without reuse.
  - 2: Weighted random assignment with reuse.
- **Autocatalytic Probability Switches**: [True, False, "square"] - Determines how the assembly index considers autocatalytic subnetworks.
- **Number of Sets**: 10,000 - The number of different autocatalytic sets to generate and analyze.

### Expected Outcome:
The script will iterate through all combinations of the provided parameters, generate the corresponding reaction networks, and analyze their assembly indices and autocatalytic subnetworks. Results, including various statistical measures (e.g., mean, median, variance), will be saved for each configuration, providing comprehensive insights into the impact of the different parameters.

### Implementation:
The script initializes the AssemblyAnalysis class for each set of parameters, generates the reaction networks, computes the assembly indices, and analyzes the autocatalytic subnetworks. The results are stored in a structured format (CSV and pickle files) for further analysis.
"""

from class_AssemblyAnalysis import *
import numpy as np
import os
import pandas as pd
import pickle
import string
import itertools

# PARAMETERS ###########################################################################################################
# Lists of parameters to iterate through for assembly analysis

# Name of the experiment
#run = "run_1000_26letters"
run = "run_10000_26_letters_squared_missed"

#final_products = ['ABCDEF', 'ABCDEFG', 'ABCDEFGH', 'ABCDEFGHI', 'ABCDEFGHIJ', 'ABCDEFGHIJK', 'ABCDEFGHIJKL', 'ABCDEFGHIJKLM', 'ABCDEFGHIJKLMN', 'ABCDEFGHIJKLMNO', 'ABCDEFGHIJKLMNOP', 'ABCDEFGHIJKLMNOPQ', 'ABCDEFGHIJKLMNOPQR', 'ABCDEFGHIJKLMNOPQRS', 'ABCDEFGHIJKLMNOPQRST', 'ABCDEFGHIJKLMNOPQRSTU', 'ABCDEFGHIJKLMNOPQRSTUV', 'ABCDEFGHIJKLMNOPQRSTUVW', 'ABCDEFGHIJKLMNOPQRSTUVWX', 'ABCDEFGHIJKLMNOPQRSTUVWXY', 'ABCDEFGHIJKLMNOPQRSTUVWXYZ']
final_products = ['ABCDEFGHIJKLMNOPQRSTUVWXY', 'ABCDEFGHIJKLMNOPQRSTUVWXYZ']



# List of final products, each represented as a string of characters.
#final_products = ['AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSs', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsT', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTt', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtU', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUu', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuV', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVv', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvW', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWw', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwX', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXx', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxY', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYy', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZ', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz']
#final_products = ['Aa', 'AaB', 'AaBb', 'AaBbC', 'AaBbCc', 'AaBbCcD', 'AaBbCcDd', 'AaBbCcDdE', 'AaBbCcDdEe', 'AaBbCcDdEeF', 'AaBbCcDdEeFf', 'AaBbCcDdEeFfG', 'AaBbCcDdEeFfGg', 'AaBbCcDdEeFfGgH', 'AaBbCcDdEeFfGgHh', 'AaBbCcDdEeFfGgHhI', 'AaBbCcDdEeFfGgHhIi', 'AaBbCcDdEeFfGgHhIiJ', 'AaBbCcDdEeFfGgHhIiJj', 'AaBbCcDdEeFfGgHhIiJjK', 'AaBbCcDdEeFfGgHhIiJjKk', 'AaBbCcDdEeFfGgHhIiJjKkL', 'AaBbCcDdEeFfGgHhIiJjKkLl', 'AaBbCcDdEeFfGgHhIiJjKkLlM', 'AaBbCcDdEeFfGgHhIiJjKkLlMm', 'AaBbCcDdEeFfGgHhIiJjKkLlMmN', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNn', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnO', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOo', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoP', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPp', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQ', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQq', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqR', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRr', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrS', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSs', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsT', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTt', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtU', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUu', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuV', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVv', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvW', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWw', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwX', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXx', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxY', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYy', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZ', 'AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz']
#final_products = ['AB', 'ABC', 'ABCD', 'ABCDE', 'ABCDEF', 'ABCDEFG', 'ABCDEFGH', 'ABCDEFGHI', 'ABCDEFGHIJ', 'ABCDEFGHIJK', 'ABCDEFGHIJKL', 'ABCDEFGHIJKLM', 'ABCDEFGHIJKLMN', 'ABCDEFGHIJKLMNO', 'ABCDEFGHIJKLMNOP', 'ABCDEFGHIJKLMNOPQ', 'ABCDEFGHIJKLMNOPQR', 'ABCDEFGHIJKLMNOPQRS', 'ABCDEFGHIJKLMNOPQRST', 'ABCDEFGHIJKLMNOPQRSTU', 'ABCDEFGHIJKLMNOPQRSTUV', 'ABCDEFGHIJKLMNOPQRSTUVW', 'ABCDEFGHIJKLMNOPQRSTUVWX', 'ABCDEFGHIJKLMNOPQRSTUVWXY', 'ABCDEFGHIJKLMNOPQRSTUVWXYZ']


# List of catalyst probabilities, each represented as a float between 0 and 1.
catalyst_probabilities = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

# List of strategies for assigning catalysts to reactions.
# 0 - Random assignment without reuse
# 1 - Weighted random assignment without reuse
# 2 - Weighted random assignment with reuse
switch_cats = [0, 1, 2]

# List of boolean switches to determine if the assembly index should be counted multiple times based on the number of autocatalytic subnetworks.
switch_aut_probs = [True, False, "square"]

# The number of different autocatalytic sets to generate.
n = 10000

print(f'Run: {run}')
# The directory path where the results will be saved.
save_path = "assembly_analysis_results"
########################################################################################################################

add = "_" + str(run)

# Check if the save_path directory exists, and create it if it does not
if not os.path.exists(save_path):
    os.makedirs(save_path)

# Save initial parameters to a text file
with open(f'{save_path}/parameters{add}.txt', 'w') as param_file:
    param_file.write(f'Run: {run}\n')
    param_file.write(f'Final Products: {final_products}\n')
    param_file.write(f'Catalyst Probabilities: {catalyst_probabilities}\n')
    param_file.write(f'Switch Catalyst Assignment: {switch_cats}\n')
    param_file.write(f'Switch Autocatalytic Probabilities: {switch_aut_probs}\n')
    param_file.write(f'Number of Sets to Generate: {n}\n')

# Initialize a list to store the results
results = []

# Iterate through all combinations of parameters
for final_product in final_products:
    for catalyst_probability in catalyst_probabilities:
        print(f'Run: {run}')
        for switch_cat in switch_cats:
            for switch_aut_prob in switch_aut_probs:

                # Print current configuration details
                print(f"\nRunning analysis with the following configuration:")
                print(f"Final product: {final_product}, Length: {len(final_product)}")
                print(f"Catalyst probability: {catalyst_probability}")
                print(f"Switch for catalyst assignment (switch_cat): {switch_cat}")
                print(f"Switch for autocatalytic probability consideration (switch_aut_prob): {switch_aut_prob}")

                # Explanation for the current switch_aut_prob
                if switch_aut_prob == True:
                    print(
                        "The current run uses switch_aut_prob = True, which means the autocatalytic subsets are considered as additional elements of the statistics.")
                elif switch_aut_prob == "square":
                    print(
                        "The current run uses switch_aut_prob = 'square', which means the autocatalytic subsets are considered sqaured as additional elements of the statistics.")
                else:
                    print(
                        "The current run uses switch_aut_prob = False, which means the autocatalytic subsets are not considered as additional elements of the statistics.")

                # Initialize the AssemblyAnalysis with the current set of parameters
                assembly_analysis = AssemblyAnalysis(
                    final_product,
                    n=n,
                    catalyst_probability=catalyst_probability,
                    switch_cat=switch_cat,
                    switch_aut_prob=switch_aut_prob,
                    print_out=False
                )

                # Find the assembly list for the generated reaction networks
                assembly_list = assembly_analysis.find_assembly_list()

                print(assembly_list)

                # Collect the statistics for the current configuration
                if assembly_list:
                    min_val = np.min(np.array(assembly_list))
                    mean_val = np.mean(np.array(assembly_list))
                    median_val = np.median(np.array(assembly_list))
                    max_val = np.max(np.array(assembly_list))
                    std_val = np.std(np.array(assembly_list))
                    var_val = np.var(np.array(assembly_list))
                    range_val = np.ptp(np.array(assembly_list))
                    q1_val = np.percentile(np.array(assembly_list), 25)
                    q3_val = np.percentile(np.array(assembly_list), 75)
                    iqr_val = q3_val - q1_val
                    length_val = len(assembly_list)
                else:
                    min_val = mean_val = median_val = max_val = std_val = var_val = range_val = q1_val = q3_val = iqr_val = length_val = float('nan')

                # Check if the assembly list is not empty and print statistics
                if assembly_list:
                    print(
                        f'Assembly List, Min: {min_val}, Mean: {mean_val}, Median: {median_val}, Max: {max_val}, Std: {std_val}, Var: {var_val}, Range: {range_val}, Q1: {q1_val}, Q3: {q3_val}, IQR: {iqr_val}, Length: {length_val}')
                else:
                    print('Assembly List is empty.')

                # Collect the counts of autocatalytic subnetworks
                subnetwork_counts = assembly_analysis.autocatalytic_subnetwork_counts
                mean_subnetwork_count = np.mean(subnetwork_counts) if subnetwork_counts else float('nan')

                print(f'Final product: {final_product}, Length: {len(final_product)}')

                # Append the results to the list
                results.append({
                    'final_product': final_product,
                    'len_final_product': len(final_product),
                    'catalyst_probability': catalyst_probability,
                    'switch_cat': switch_cat,
                    'switch_aut_prob': switch_aut_prob,
                    'min_assembly_index': min_val,
                    'mean_assembly_index': mean_val,
                    'median_assembly_index': median_val,
                    'max_assembly_index': max_val,
                    'std_assembly_index': std_val,
                    'var_assembly_index': var_val,
                    'range_assembly_index': range_val,
                    'q1_assembly_index': q1_val,
                    'q3_assembly_index': q3_val,
                    'iqr_assembly_index': iqr_val,
                    'length_assembly_list': length_val,
                    'mean_subnetwork_count': mean_subnetwork_count
                })

                # Save the assembly list as a pickle file
                pickle_filename = f'{save_path}/assembly_list_{run}_{final_product}_{catalyst_probability}_{switch_cat}_{switch_aut_prob}.pkl'
                with open(pickle_filename, 'wb') as pickle_file:
                    pickle.dump(assembly_list, pickle_file)

                # Convert the results to a DataFrame and save to a CSV file
                df = pd.DataFrame(results)
                df.to_csv(f'{save_path}/assembly_analysis_results{add}.csv', index=False)

print(df)
print("Assembly analysis completed and results saved.")
