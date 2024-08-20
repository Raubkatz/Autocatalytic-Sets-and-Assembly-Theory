"""
Author: Dr. techn. Sebastian Raubitzek MSc. BSc.; SBA Research, Complexity and Resilience Research Group
This Python script generates visualizations to analyze the impact of different catalyst probabilities, catalyst assignment strategies, and autocatalytic amplification methods on the assembly of autocatalytic sets. The script produces line plots that compare the mean, standard deviation, and minimum assembly indices across various experimental conditions.

### Problem Description:
The goal is to visually examine how different experimental configurations affect the depth and variability of reaction networks generated in an assembly analysis. By plotting these metrics against the length of the final product, the script provides insights into the underlying patterns and differences induced by varying catalyst probabilities and other parameters.

### Parameters Used:
- **Run Name**: 'run_10000_26_letters_squared' - The specific experiment being analyzed.
- **Discard Filters**:
  - `discard_1`: True - Configurations with `switch_cat = 1` are excluded.
  - `discard_linear`: True - Configurations with `switch_aut_prob = "True"` are excluded.
- **Catalyst Probabilities**: Extracted from the dataset - Each unique catalyst probability value is used to create separate plots.

### Visualization:
The script generates the following plots for each catalyst probability:
1. **Mean Assembly Index**: Plots the average reaction/event depth against the length of the final product.
2. **Standard Deviation of Assembly Index**: Plots the variability in reaction/event depth against the length of the final product.
3. **Minimum Assembly Index**: Plots the minimum depth required for assembly against the length of the final product.

### Titles and Labels:
- **Switch Catalyst Titles**: Descriptions for different catalyst assignment strategies are provided in the plots.
  - 0: "Random catalyst assignment without reuse"
  - 1: "Weighted random catalyst assignment without reuse"
  - 2: "Weighted random catalyst assignment with reuse"
- **Switch Autocatalytic Probability Titles**: Descriptions for how autocatalytic amplification is handled:
  - "True": "Linear autocatalysis amplification"
  - "False": "No autocatalysis amplification"
  - "square": "Exponential autocatalysis amplification"

### Output:
- **PNG and EPS Files**: Each plot is saved as both PNG and EPS files in the `assembly_analysis_results` directory.

### Implementation:
The script reads the processed assembly analysis results, filters the data based on the specified conditions, and iteratively generates plots for each combination of catalyst probability, catalyst assignment strategy, and autocatalytic amplification method.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set the color palette to match the earlier surface plots
sns.set_palette('viridis')

run = "run_10000_26_letters_squared"

discard_1 = True
discard_linear = True
add = ""

if discard_1:
    add = add + "_discard1"
if discard_linear:
    add = add + "_discardlin"

df = pd.read_csv(f'./assembly_analysis_results/assembly_analysis_results_{run}_{add}_updated.csv')

# Title mappings
switch_cat_titles = {
    0: "Random catalyst assignment without reuse",
    1: "Weighted random catalyst assignment without reuse",
    2: "Weighted random catalyst assignment with reuse"
}

switch_aut_prob_titles = {
    "True": "Linear autocatalysis amplification",
    "False": "No autocatalysis amplification",
    "square": "Exponential autocatalysis amplification"
}

catalyst_probabilities = df['catalyst_probability'].unique()

for catalyst_probability in catalyst_probabilities:
    subset_cp = df[df['catalyst_probability'] == catalyst_probability]

    # Plotting the mean assembly index
    plt.figure(figsize=(12, 8))
    for switch_cat in subset_cp['switch_cat'].unique():
        for switch_aut_prob in subset_cp['switch_aut_prob'].unique():
            subset = subset_cp[(subset_cp['switch_cat'] == switch_cat) & (subset_cp['switch_aut_prob'] == switch_aut_prob)]
            label = f"{switch_cat_titles[switch_cat]}\n{switch_aut_prob_titles[switch_aut_prob]}"
            plt.plot(subset['len_final_product'], subset['mean_assembly_index'], marker='o', label=label)

    plt.xlabel('Length of Final Product')
    plt.ylabel('Average Reaction/Event Depth')
    plt.title(f'Average Reaction/Event Depth vs Length of Final Product\n(Catalyst Probability = {catalyst_probability})')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'./assembly_analysis_results/line_plot_mean_{str(catalyst_probability).replace(".","")}.png')
    plt.savefig(f'./assembly_analysis_results/line_plot_mean_{str(catalyst_probability).replace(".","")}.eps')
    plt.show()

    # Plotting the std assembly index
    plt.figure(figsize=(12, 8))
    for switch_cat in subset_cp['switch_cat'].unique():
        for switch_aut_prob in subset_cp['switch_aut_prob'].unique():
            subset = subset_cp[(subset_cp['switch_cat'] == switch_cat) & (subset_cp['switch_aut_prob'] == switch_aut_prob)]
            label = f"{switch_cat_titles[switch_cat]}\n{switch_aut_prob_titles[switch_aut_prob]}"
            plt.plot(subset['len_final_product'], subset['std_assembly_index'], marker='o', label=label)

    plt.xlabel('Length of Final Product')
    plt.ylabel('StdDev Reaction/Event')
    plt.title(f'Standard Deviation of Reaction/Event Depth vs Length of Final Product\n(Catalyst Probability = {catalyst_probability})')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'./assembly_analysis_results/line_plot_std_{str(catalyst_probability).replace(".","")}.png')
    plt.savefig(f'./assembly_analysis_results/line_plot_std_{str(catalyst_probability).replace(".","")}.eps')
    plt.show()

    # Plotting the minimum assembly index
    plt.figure(figsize=(12, 8))
    for switch_cat in subset_cp['switch_cat'].unique():
        for switch_aut_prob in subset_cp['switch_aut_prob'].unique():
            subset = subset_cp[(subset_cp['switch_cat'] == switch_cat) & (subset_cp['switch_aut_prob'] == switch_aut_prob)]
            label = f"{switch_cat_titles[switch_cat]}\n{switch_aut_prob_titles[switch_aut_prob]}"
            plt.plot(subset['len_final_product'], subset['min_assembly_index'], marker='o', label=label)

    plt.xlabel('Length of Final Product')
    plt.ylabel('Minimum Assembly Index')
    plt.title(f'Minimum Assembly Index vs Length of Final Product\n(Catalyst Probability = {catalyst_probability})')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'./assembly_analysis_results/line_plot_min_{str(catalyst_probability).replace(".","")}.png')
    plt.savefig(f'./assembly_analysis_results/line_plot_min_{str(catalyst_probability).replace(".","")}.eps')
    plt.show()
