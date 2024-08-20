"""
Author: Dr. techn. Sebastian Raubitzek MSc. BSc.; SBA Research, Complexity and Resilience Research Group
This Python script processes the results from an assembly analysis experiment, focusing on how different parameter configurations affect the assembly process of autocatalytic sets. The script filters, adjusts, and compares the assembly indices against a baseline, ultimately saving the processed data for further analysis.

### Problem Description:
We aim to analyze the results of an experiment that generated autocatalytic sets with varying parameters. Specifically, the goal is to compare each configuration's assembly indices against a defined baseline to identify significant deviations or trends.

### Parameters Used:
- **Run Name**: 'run_10000_26_letters_squared' - The identifier for the specific experiment being analyzed.
- **Discard Filters**:
  - `discard_1`: True - If set, configurations with `switch_cat = 1` are excluded from the analysis.
  - `discard_linear`: True - If set, configurations with `switch_aut_prob = "True"` are excluded from the analysis.
- **Baseline Definition**: The baseline is defined as configurations where `switch_cat = 0` and `switch_aut_prob = "False"`, representing a standard comparison point.

### Data Processing:
1. **Filtering**: The script loads the experiment data, filters it based on the specified `len_final_product` range (6 to 20 characters), and applies the discard filters as needed.
2. **Baseline Extraction**: A subset of the data matching the baseline criteria is extracted for use in comparisons.
3. **Baseline Assignment**: For each unique combination of `catalyst_probability`, `len_final_product`, `switch_cat`, and `switch_aut_prob`, the script calculates the differences between the assembly indices and their respective baseline values.
4. **Difference Calculation**: The script computes the difference between each configuration's mean and standard deviation of assembly indices and the corresponding baseline values.

### Output:
- **CSV File**: The processed DataFrame, including the baseline-adjusted assembly indices and their differences, is saved to a new CSV file for further analysis.

### Implementation:
The script filters the input data, extracts baseline values, assigns these baselines to the relevant configurations, and calculates the differences. The final results are saved and displayed for verification.
"""



import pandas as pd

run = "run_10000_26_letters_squared"

discard_1 = True
discard_linear = True
add = ""

if discard_1:
    add = add + "_discard1"
if discard_linear:
    add = add + "_discardlin"

# Load the CSV file into a DataFrame
df = pd.read_csv(f'./assembly_analysis_results/assembly_analysis_results_{run}.csv')

# Filter the DataFrame to the necessary range
df = df[df['len_final_product'].between(6, 20)]


# Discard rows where switch_cat = 1 if discard_1 is set to True
if discard_1:
    df = df[df['switch_cat'] != 1]
if discard_linear:
    df = df[df['switch_aut_prob'] != "True"]

# Extract the baseline DataFrame
baseline_df = df[(df['switch_cat'] == 0) & (df['switch_aut_prob'] == "False")]

print(baseline_df)

# Define a function to assign baselines and calculate differences
def assign_baseline(df, baseline_df):
    # Initialize baseline columns
    df['baseline_mean_assembly_index'] = 0
    df['baseline_std_assembly_index'] = 0
    df['mean_assembly_diff'] = 0
    df['std_assembly_diff'] = 0

    # Iterate over each unique combination of catalyst_probability and len_final_product
    for catalyst_probability in df['catalyst_probability'].unique():
        for len_final_product in df['len_final_product'].unique():
            for switch_cat in df['switch_cat'].unique():
                for switch_aut_prob in df['switch_aut_prob'].unique():

                    # Find the baseline values for the current grid point
                    baseline_values = baseline_df[(baseline_df['catalyst_probability'] == catalyst_probability) & (baseline_df['len_final_product'] == len_final_product)]

                    baseline_mean = baseline_values['mean_assembly_index'].values
                    baseline_std = baseline_values['std_assembly_index'].values

                    print(baseline_values)

                    current_assembly_values = df[(df['catalyst_probability'] == catalyst_probability) & (df['len_final_product'] == len_final_product) & (df['switch_cat'] == switch_cat) & (df['switch_aut_prob'] == switch_aut_prob)]

                    # Assign baseline values
                    df.loc[current_assembly_values.index, 'baseline_mean_assembly_index'] = baseline_mean
                    df.loc[current_assembly_values.index, 'baseline_std_assembly_index'] = baseline_std

                    # Calculate differences
                    df.loc[current_assembly_values.index, 'mean_assembly_diff'] = current_assembly_values['mean_assembly_index'] - baseline_mean
                    df.loc[current_assembly_values.index, 'std_assembly_diff'] = current_assembly_values['std_assembly_index'] - baseline_std

    return df

# Apply the function to the DataFrame
df = assign_baseline(df, baseline_df)

# Save the updated DataFrame to a CSV file

output_file_path = f'./assembly_analysis_results/assembly_analysis_results_{run}_{add}_updated.csv'

df.to_csv(output_file_path, index=False)

# Display the first few rows of the updated DataFrame
print(df.head())
