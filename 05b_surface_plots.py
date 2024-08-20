"""
Author: Dr. techn. Sebastian Raubitzek MSc. BSc.
This Python script generates 3D surface plots to visualize the impact of different catalyst probabilities, catalyst assignment strategies, and autocatalytic amplification methods on the assembly of autocatalytic sets. The script focuses on analyzing the differences in mean and standard deviation of assembly indices across various configurations.

### Problem Description:
The script visualizes how different configurations in the assembly process affect the mean and standard deviation of reaction/event depths. By creating 3D surface plots, the script provides a clear representation of how these factors interact, offering insights into the complexity and variability of the resulting networks.

### Parameters Used:
- **Run Name**: 'run_10000_26_letters_squared' - The identifier for the specific experiment being analyzed.
- **Discard Filters**:
  - `discard_1`: True - Configurations with `switch_cat = 1` are excluded.
  - `discard_linear`: True - Configurations with `switch_aut_prob = "True"` are excluded.
- **Title Mappings**:
  - **Catalyst Assignment Strategies (switch_cat)**:
    - 0: "Random catalyst assignment without reuse"
    - 1: "Weighted random catalyst assignment without reuse"
    - 2: "Weighted random catalyst assignment with reuse"
  - **Autocatalytic Probability Handling (switch_aut_prob)**:
    - "True": "Linear autocatalysis amplification"
    - "False": "No autocatalysis amplification"
    - "square": "Exponential autocatalysis amplification"

### Visualization:
The script produces the following 3D surface plots:
1. **Mean Assembly Index Difference (mean_assembly_diff)**: Visualizes the difference in mean assembly indices compared to a baseline.
2. **Standard Deviation of Assembly Index Difference (std_assembly_diff)**: Visualizes the difference in standard deviation of assembly indices compared to a baseline.

### Output:
- **PNG and EPS Files**: Each 3D surface plot is saved as both PNG and EPS files in the `assembly_analysis_results` directory.

### Implementation:
1. **Data Preparation**: The script reads the processed assembly analysis results and filters the data based on the specified conditions.
2. **3D Surface Plot Generation**:
   - The script creates a 3D surface plot for the mean assembly index differences (`mean_assembly_diff`).
   - It also creates a 3D surface plot for the standard deviation differences (`std_assembly_diff`).
3. **Loop Through Configurations**: The script iterates over each unique combination of `switch_cat` and `switch_aut_prob`, generating the corresponding plots.

### Execution:
To generate the plots, the script pivots the data to form a 3D grid, then uses Matplotlib to plot the surfaces, and finally saves the plots for further analysis.
"""


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

run = "run_10000_26_letters_squared"

discard_1 = True
discard_linear = True
add = ""

if discard_1:
    add = add + "_discard1"
if discard_linear:
    add = add + "_discardlin"


df = pd.read_csv(f'./assembly_analysis_results/assembly_analysis_results_{run}_{add}_updated.csv')


# Define the title mappings for switch_cat and switch_aut_prob
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


# Define a function to create a 3D surface plot for mean_assembly_diff
def create_3d_surface_plot_mean(df, switch_cat, switch_aut_prob):
    filtered_df = df[(df['switch_cat'] == switch_cat) & (df['switch_aut_prob'] == switch_aut_prob)]

    if filtered_df.empty:
        print(f"No data available for switch_cat={switch_cat} and switch_aut_prob={switch_aut_prob}")
        return

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Pivot the data to create a 3D grid
    pivot_table = filtered_df.pivot_table(index='len_final_product', columns='catalyst_probability',
                                          values='mean_assembly_diff')
    X, Y = np.meshgrid(pivot_table.columns, pivot_table.index)
    Z = pivot_table.values

    # Plot the surface
    ax.plot_surface(X, Y, Z, cmap='viridis')

    # Set labels
    ax.set_ylabel('Length of Final Product')
    ax.set_xlabel('Catalyst Probability')
    ax.set_zlabel('Average Reaction/Event Depth')

    # Set title
    ax.set_title(f"Average Reaction/Event Depth\n{switch_cat_titles[switch_cat]}\n{switch_aut_prob_titles[switch_aut_prob]}", pad=20)

    plt.tight_layout()
    plt.savefig(f'./assembly_analysis_results/surface_plot_mean_{str(switch_cat)}_{str(switch_aut_prob)}.png')
    plt.savefig(f'./assembly_analysis_results/surface_plot_mean_{str(switch_cat)}_{str(switch_aut_prob)}.eps')
    plt.show()


# Define a function to create a 3D surface plot for std_assembly_diff
def create_3d_surface_plot_std(df, switch_cat, switch_aut_prob):
    filtered_df = df[(df['switch_cat'] == switch_cat) & (df['switch_aut_prob'] == switch_aut_prob)]

    if filtered_df.empty:
        print(f"No data available for switch_cat={switch_cat} and switch_aut_prob={switch_aut_prob}")
        return

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Pivot the data to create a 3D grid
    pivot_table = filtered_df.pivot_table(index='len_final_product', columns='catalyst_probability',
                                          values='std_assembly_diff')
    X, Y = np.meshgrid(pivot_table.columns, pivot_table.index)
    Z = pivot_table.values

    # Plot the surface
    ax.plot_surface(X, Y, Z, cmap='viridis')

    # Set labels
    ax.set_ylabel('Length of Final Product')
    ax.set_xlabel('Catalyst Probability')
    ax.set_zlabel('StdDev Reaction/Event Depth')

    # Set title
    ax.set_title(f"StdDev Reaction/Event Depth\n{switch_cat_titles[switch_cat]}\n{switch_aut_prob_titles[switch_aut_prob]}", pad=20)

    plt.tight_layout()
    plt.savefig(f'./assembly_analysis_results/surface_plot_std_{str(switch_cat)}_{str(switch_aut_prob)}.png')
    plt.savefig(f'./assembly_analysis_results/surface_plot_std_{str(switch_cat)}_{str(switch_aut_prob)}.eps')
    plt.show()


# Loop through each combination of switch_cat and switch_aut_prob and create plots
for switch_cat in df['switch_cat'].unique():
    for switch_aut_prob in df['switch_aut_prob'].unique():
        create_3d_surface_plot_mean(df, switch_cat, switch_aut_prob)
        create_3d_surface_plot_std(df, switch_cat, switch_aut_prob)
