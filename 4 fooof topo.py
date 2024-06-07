# Import required code for visualizing example models
import os, io, shutil
# ensure working directory
cwd = os.path.dirname(os.path.realpath(__file__))
os.chdir(cwd)
import pandas as pd
import numpy as np
# the preprocess package
from fooof import FOOOF
# from fooof.sim.utils import set_random_seed
from fooof.plts.spectra import plot_spectra
from fooof.plts.annotate import plot_annotated_model
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from statsmodels.stats.anova import AnovaRM
from statsmodels.regression.mixed_linear_model import MixedLM
from scipy import stats
import re

# Specify the directory
cut_off_component = 15 # the last 15 Hz is not used
max_Hz = 40 - cut_off_component
freqs = np.array(range(2, 41-cut_off_component)) # 2 to 40 Hz

# Specify the output directory
dir_path = f'5 alpha topo/'
# Specify the input directory
input_dir = 'topography data'

# Check if dir_path is not empty and it exists
if dir_path and os.path.exists(dir_path):
    # Remove the directory and its contents
    if os.path.isdir(dir_path):
        shutil.rmtree(dir_path)
    else:
        os.remove(dir_path)
# Recreate the directory
os.makedirs(dir_path)


# create the linear result dataframe to save the results
linear_result = pd.DataFrame()


def get_peak_params(fm):
    # in case the FOOOF did not detect any peaks in a channel
    if fm.peak_params_.size != 0:
        peak = pd.DataFrame(fm.peak_params_, columns=['CF', 'PW', 'BW'])
    else:
        peak = pd.DataFrame({'CF': [np.nan], 'PW': [0], 'BW': [0]})
    peak[['offset', 'exponent']] = pd.DataFrame(fm.aperiodic_params_.reshape(1, -1), columns=['offset', 'exponent'])
    peak['PID'] = np.full(peak.shape[0], pid)
    peak['exp'] = np.full(peak.shape[0], exp)
    peak['load'] = np.full(peak.shape[0], load)
    peak['age'] = np.full(peak.shape[0], group)
    peak['r_squared'] = np.full(peak.shape[0], fm.r_squared_)
    peak['error'] = np.full(peak.shape[0], fm.error_)
    # move PID to the first column
    cols = list(peak.columns)
    cols.insert(0, cols.pop(cols.index('PID')))
    cols.insert(1, cols.pop(cols.index('exp')))
    cols.insert(2, cols.pop(cols.index('load')))
    cols.insert(3, cols.pop(cols.index('age')))
    peak = peak.loc[:, cols]
    return peak


# Get a list of all files in the data directory
all_files = os.listdir(input_dir)

# Filter the list to only include files with 'EEG' in the name
alpha_files = [file for file in all_files if 'EEG' in file]
alpha_files = sorted(alpha_files, key=lambda file: int(re.search(r'\d+', file).group()))

# Loop over the alpha files
for file in alpha_files:
    file_name, _ = os.path.splitext(file)
    # Read the data from the file
    data = pd.read_csv(f'{input_dir}/{file}')
    data = data.iloc[:, 1:] 
    data = data.sort_values(['id', 'exp', 'load']).reset_index(drop=True)

    linear_result = pd.DataFrame()
    for index, row in data.iterrows():
        sub_powers = np.array(row[1:-(4+cut_off_component)]) + 1e-10 # add a small number to avoid zero
        pid = row['id']
        exp = row['exp']
        load = row['load']
        group = row['age']
        # use fixed model, the peak_width_limits is set to 3-7 Hz, min_peak_height is set to 0.15, peak_threshold is set to 2, max_n_peaks is set to 4
        fm2 = FOOOF(peak_width_limits=[3, 7], min_peak_height=0.15, peak_threshold=2, max_n_peaks=4, aperiodic_mode='fixed', verbose=False)
        fm2.fit(freqs, sub_powers)
        peak_l = get_peak_params(fm2)
        peak_l['offset'] = peak_l['offset'].fillna(method='ffill')
        peak_l['exponent'] = peak_l['exponent'].fillna(method='ffill')
        # Add a 'file_name' column to peak_l
        peak_l['file_name'] = file_name
        linear_result = pd.concat([linear_result, peak_l], axis=0)

    linear_result.to_csv(os.path.join(dir_path, f'{file_name}_alpha_result.csv'), index=False)
    print(f'{file_name} alpha finished.')

# for the theta
# Specify the directory
cut_off_component = 24 # the last 15 Hz is not used
max_Hz = 40 - cut_off_component
freqs = np.array(range(2, 41-cut_off_component)) # 2 to 16 Hz

dir_path = f'5 theta topo/'

# Check if dir_path is not empty and it exists
if dir_path and os.path.exists(dir_path):
    # Remove the directory and its contents
    if os.path.isdir(dir_path):
        shutil.rmtree(dir_path)
    else:
        os.remove(dir_path)
# Recreate the directory
os.makedirs(dir_path)


# create the linear result dataframe to save the results
linear_result = pd.DataFrame()

# Filter the list to only include files with 'EEG' in the name
theta_files = [file for file in all_files if 'EEG' in file]
theta_files = sorted(theta_files, key=lambda file: int(re.search(r'\d+', file).group()))
# Loop over the alpha files
for file in theta_files:
    file_name, _ = os.path.splitext(file)
    # Read the data from the file
    data = pd.read_csv(f'{input_dir}/{file}')
    data = data.iloc[:, 1:] 
    data = data.sort_values(['id', 'exp', 'load']).reset_index(drop=True)

    linear_result = pd.DataFrame()
    for index, row in data.iterrows():
        sub_powers = np.array(row[1:-(4+cut_off_component)]) + 1e-10 # add a small number to avoid zero
        pid = row['id']
        exp = row['exp']
        load = row['load']
        group = row['age']
        # use fixed model, the peak width limit is set to 1.5-7 Hz, min_peak_height is set to 0, peak_threshold is set to 1, max_n_peaks is set to 4
        fm2 = FOOOF(peak_width_limits=[1.5, 7], min_peak_height=0, peak_threshold=1, max_n_peaks=4, aperiodic_mode='fixed', verbose=False)
        fm2.fit(freqs, sub_powers)
        peak_l = get_peak_params(fm2)
        peak_l['offset'] = peak_l['offset'].fillna(method='ffill')
        peak_l['exponent'] = peak_l['exponent'].fillna(method='ffill')
        # Add a 'file_name' column to peak_l
        peak_l['file_name'] = file_name
        linear_result = pd.concat([linear_result, peak_l], axis=0)

    linear_result.to_csv(os.path.join(dir_path, f'{file_name}_theta_result.csv'), index=False)
    print(f'{file_name} theta finished.')