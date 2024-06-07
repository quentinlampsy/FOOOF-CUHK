import os
# ensure working directory
cwd = os.path.dirname(os.path.realpath(__file__))
os.chdir(cwd)
import pandas as pd
import numpy as np
import shutil
import re

# set up input folder
dir_path = f'5 alpha topo/'
output_dir = f'5 topo cleaned/'

# Specify the frequency range for alpha range
cut_off_component = 15 # the last 15 Hz is not used
max_Hz = 40 - cut_off_component
freqs = np.array(range(2, 41-cut_off_component)) # 2 to 25 Hz


# Check if output_dir is not empty and it exists
if output_dir and os.path.exists(output_dir):
    # Remove the directory and its contents
    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir)
    else:
        os.remove(output_dir)
# Recreate the directory
os.makedirs(output_dir)

# Get a list of all files in the directory
all_files = os.listdir(dir_path)

# Filter the list to only include files with 'alpha' in the name and end with '.csv'
linear_result_files = [file for file in all_files if 'alpha' in file and file.endswith('.csv')]
linear_result_files = sorted(linear_result_files, key=lambda file: int(re.search(r'\d+', file).group()))

# Initialize an empty DataFrame to store all data
stat_data = pd.DataFrame()

# Loop over the linear_result files
for file in linear_result_files:
    file_name, _ = os.path.splitext(file)
    file_name = file_name.replace('_alpha_result', '')
    file_name = file_name.replace('_all_raw', '')
    # Read the data from the file
    linear_result = pd.read_csv(f'{dir_path}/{file}')


    no_peak = linear_result[linear_result['CF'].isna()]
    peak_else = linear_result[(linear_result['CF'] < 8.2) | (linear_result['CF'] > 13.5)]
    peak_else = pd.concat([peak_else, no_peak], axis=0)

    # Replace 'CF' values with NaN and 'PW' with 0
    linear_result.loc[peak_else.index, 'CF'] = np.nan
    linear_result.loc[peak_else.index, 'PW'] = 0
    linear_result.loc[peak_else.index, 'BW'] = 0
    # linear_result = linear_result.drop(no_peak.index)

    # keep the one with largest PW
    linear_result = linear_result.loc[linear_result.groupby(['PID', 'exp', 'load', 'age'])['PW'].idxmax()]
    linear_result.reset_index(drop=True, inplace=True)

    linear_result['file_name'] = file_name
    linear_result['band'] = 'alpha'
    linear_result = linear_result.sort_values(by=['PID', 'exp', 'load', 'age'])
    stat_data = pd.concat([stat_data, linear_result])

# Specify the frequency range for alpha range
cut_off_component = 24 # the last 24 Hz is not used
max_Hz = 40 - cut_off_component
freqs = np.array(range(2, 41-cut_off_component)) # 2 to 16 Hz

# do the analysis on the specfic folder
dir_path = f'5 theta topo/'


# Get a list of all files in the directory
all_files = os.listdir(dir_path)

# Filter the list to only include files with 'linear_result' in the name and end with '.csv'
linear_result_files = [file for file in all_files if 'theta' in file and file.endswith('.csv')]
linear_result_files = sorted(linear_result_files, key=lambda file: int(re.search(r'\d+', file).group()))
# Loop over the linear_result files
for file in linear_result_files:
    file_name, _ = os.path.splitext(file)
    file_name = file_name.replace('_theta_result', '')
    file_name = file_name.replace('_all_raw', '')
    # Read the data from the file
    linear_result = pd.read_csv(f'{dir_path}/{file}')

    no_peak = linear_result[linear_result['CF'].isna()]
    peak_else = linear_result[(linear_result['CF'] < 3.2) | (linear_result['CF'] > 7.8)]
    peak_else = pd.concat([peak_else, no_peak], axis=0)

    # Replace 'CF' values with NaN and 'PW' with 0
    linear_result.loc[peak_else.index, 'CF'] = np.nan
    linear_result.loc[peak_else.index, 'PW'] = 0
    linear_result.loc[peak_else.index, 'BW'] = 0

    # keep the one with largest PW
    linear_result = linear_result.loc[linear_result.groupby(['PID', 'exp', 'load', 'age'])['PW'].idxmax()]
    linear_result.reset_index(drop=True, inplace=True)

    linear_result['file_name'] = file_name
    linear_result['band'] = 'theta'
    linear_result = linear_result.sort_values(by=['PID', 'exp', 'load', 'age'])
    stat_data = pd.concat([stat_data, linear_result])

# Group the data by 'age' and 'exp'
stat_data['PW'] = stat_data['PW'].fillna(0)
grouped = stat_data.groupby(['age', 'exp', 'load'])

print(stat_data.head())

stat_data.to_csv(f'{output_dir}/permu_all_sub_data.csv')
# Loop over the groups
for (age, exp, load), group in grouped:
    # Pivot the data
    for para in ['CF', 'PW', 'offset', 'exponent']:
        pivot_df = group.pivot_table(index=['PID', 'band'], columns='file_name', values=para)
        
        # Reset the index
        pivot_df = pivot_df.reset_index()
        # Extract the columns starting with 'EEG'
        eeg_cols = [col for col in pivot_df.columns if col.startswith('EEG')]
        # Sort the columns starting with 'EEG'
        eeg_cols = sorted(eeg_cols, key=lambda col: int(re.search(r'\d+', col).group()))
        
        # Rearrange the columns
        pivot_df = pivot_df[['PID', 'band'] + eeg_cols]
        
        # Flatten the MultiIndex columns
        if age == 0:
            age_name = 'child'
        elif age == 1:
            age_name = 'adult'
        file_name = f'{age_name}_{exp}_{load}_{para}'
        
        # Save the pivoted data to a CSV file
        pivot_df.to_csv(f'{output_dir}/{file_name}.csv')

print(pivot_df.head())