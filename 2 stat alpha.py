import os
# ensure working directory
cwd = os.path.dirname(os.path.realpath(__file__))
os.chdir(cwd)
import pandas as pd
import numpy as np
import shutil
from scipy import stats
import matplotlib.pyplot as plt
from PIL import Image
import itertools

# do the analysis on the specfic folder
dir_path = f'2 alpha/'
output_dir = f'3 alpha_stat/'
output_pic_dir = os.path.join(cwd, '4 alpha pic')

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
os.makedirs(output_pic_dir, exist_ok=True)

# Get a list of all files in the directory
all_files = os.listdir(dir_path)

# Filter the list to only include files with 'linear_result' in the name and end with '.csv'
linear_result_files = [file for file in all_files if 'linear_result' in file and file.endswith('.csv')]

# Initialize an empty DataFrame to store all data
stat_data = pd.DataFrame()

# Loop over the linear_result files
for file in linear_result_files:
    file_name, _ = os.path.splitext(file)
    file_name = file_name.replace('_linear_result', '')
    # Read the data from the file
    linear_result = pd.read_csv(f'{dir_path}/{file}')
    linear_result['PID'] = linear_result['PID'].astype(str)
    linear_result['offset'] = linear_result['offset'].fillna(method='ffill')
    linear_result['exponent'] = linear_result['exponent'].fillna(method='ffill')
    # print(linear_result.head())

    no_peak = linear_result[linear_result['CF'].isna()]

    peak_else = linear_result[(linear_result['CF'] < 8.2) | (linear_result['CF'] > 13.5)]
    peak_else = pd.concat([peak_else, no_peak], axis=0)
    peak_else.to_csv(f'{output_dir}{file_name}_peak_else_sub.csv')

    # Replace 'CF' values with NaN and 'PW' with 0
    linear_result.loc[peak_else.index, 'CF'] = np.nan
    linear_result.loc[peak_else.index, 'PW'] = 0
    linear_result.loc[peak_else.index, 'BW'] = 0
    linear_result = linear_result.drop(no_peak.index)
    # filter out the participants with r_squared < 0.7
    linear_result = linear_result[linear_result['r_squared'] > 0.7]
    # keep the one with largest PW, if there are multiple peaks
    linear_result = linear_result.loc[linear_result.groupby(['PID', 'exp', 'load', 'age'])['PW'].idxmax()]
    linear_result.reset_index(drop=True, inplace=True)
    linear_result.to_csv(f'{output_dir}{file_name}_stat_data.csv')
    # Append the data to linear_result
    stat_data = pd.concat([stat_data, linear_result])

# print(stat_data.head())
# Group the data by 'PID', 'exp', 'load', and 'age', and calculate the mean
mean_values = stat_data.groupby(['PID', 'exp', 'load', 'age'])[['CF', 'PW', 'BW', 'offset', 'exponent']].mean()

# Reset the index to make 'PID', 'exp', 'load', and 'age' regular columns again
mean_values.reset_index(inplace=True)

# Count the number of occurrences of each 'PID' value
pid_counts = mean_values['PID'].value_counts()
if any(pid_counts != 4):
    print("Alert: Some PID values do not appear exactly 4 times! Which indicates data missing.")
else:
    print("All PID values appear exactly 4 times.")

mean_values.to_csv(f'{output_dir}mean_alpha_params.csv', index=False)
print(mean_values.head())


# plot the figures
# one sample t test
group_1t = mean_values.groupby(['PID', 'age', 'exp', 'load'])
# Create a new DataFrame to store the results
result = pd.DataFrame(columns=['PID', 'age', 'exp', 'CF_diff', 
                               'PW_diff', 'BW_diff', 'offset_diff', 'exponent_diff'])

# Loop over the groups
for name, group in group_1t:
    # If the load is '1b', store the values
    if name[3] == '1b':
        cf_1b = group['CF'].values
        pw_1b = group['PW'].values
        bw_1b = group['BW'].values
        offset_1b = group['offset'].values
        exponent_1b = group['exponent'].values
    # If the load is '2b', calculate the difference and store the result
    elif name[3] == '2b':
        cf_diff = group['CF'].values - cf_1b
        pw_diff = group['PW'].values - pw_1b
        bw_diff = group['BW'].values - bw_1b
        offset_diff = group['offset'].values - offset_1b
        exponent_diff = group['exponent'].values - exponent_1b
        # Create a temporary DataFrame from the dictionary
        temp_df = pd.DataFrame({'PID': name[0], 'age': name[1], 'exp': name[2], 
                                'CF_diff': cf_diff, 'PW_diff': pw_diff, 'BW_diff': bw_diff, 
                                'offset_diff': offset_diff, 'exponent_diff': exponent_diff})
       # Append the temporary DataFrame to the result DataFrame
        result = pd.concat([result, temp_df], ignore_index=True)

# print(result.head())
result.to_csv(f'{output_dir}mean_alpha_delta_nback_params.csv', index=False)
# We want to test if the population mean is 0.0
popmean = 0.0
# Group the data by 'age' and 'exp'
groups_1test = result.groupby(['age', 'exp'])
# Initialize an empty DataFrame to store the results
results_df = pd.DataFrame(columns=['age', 'exp', 'column', 'T statistic', 'P value'])
for (age, exp), group in groups_1test:
    # print(f"Age: {age}, Exp: {exp}")
    # Loop over the columns in the result DataFrame
    for column in result.columns:
        # Skip the 'PID', 'age', and 'exp' columns
        if column in ['PID', 'age', 'exp']:
            continue
        # Remove NaN values by row through each column
        group = group.dropna(subset=[column])

        # Perform the one-sample t-test
        t_statistic, p_value = stats.ttest_1samp(group[column], popmean)

        # Store the results in the DataFrame
        temp_df = pd.DataFrame({
            'age': age,
            'exp': exp,
            'column': column,
            'T statistic': round(t_statistic, 3),
            'P value': round(p_value, 3)
        }, index=[0])
        # Concatenate the new data to the results DataFrame
        results_df = pd.concat([results_df, temp_df], ignore_index=True)
print(results_df)
results_df.to_csv(f'{output_dir}one_sample_t_test.csv', index=False)



# the fundtion to calculate the aperiodic power values
def calculate_y_ap(group):
    offset = np.repeat(group['offset'].values, len(freqs))
    exp = np.repeat(group['exponent'].values, len(freqs))
    y_ap = offset - np.log10(freqs ** exp)
    return y_ap

# the fundtion to calculate the periodic power values
def calculate_y_p(group):
    cf = np.repeat(group['CF'].values, len(freqs))
    pw = np.repeat(group['PW'].values, len(freqs))
    bw = np.repeat(group['BW'].values, len(freqs))

    ys = np.zeros_like(freqs)
    y_p = ys + pw * np.exp(-(freqs-cf)**2 / (2*bw**2))
    return y_p


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


groups = mean_values.groupby(['PID', 'age', 'load', 'exp'])
# Apply the function to each group
ap_power = groups.apply(calculate_y_ap).reset_index()
p_power = groups.apply(calculate_y_p).reset_index()

sub_info = ap_power[['PID', 'age', 'load', 'exp']]
# Create a new DataFrame with columns named from 'V2' to 'V{max_Hz}'
columns = ['V' + str(i) for i in range(2, max_Hz + 1)]
ap_power = pd.DataFrame(ap_power[0].to_list(), columns=columns)
p_power = pd.DataFrame(p_power[0].to_list(), columns=columns)
ap_power_all = pd.concat([sub_info, ap_power], axis=1)
p_power_all = pd.concat([sub_info, p_power], axis=1)
# print(ap_power_all.head())

###########################################
# Alpha FOOOFed oscilation plot in Fig. 4 #
###########################################

# Merge raw_data_eeg and ap_power_all on the group columns
merged_df = pd.merge(p_power_all, ap_power_all, on=['PID', 'age', 'load', 'exp'], suffixes=('_p', '_ap'))
# Subtract the ap_power_all columns from the raw_data_eeg columns
for i in range(2, max_Hz+1):
    merged_df[f'V{i}'] = merged_df[f'V{i}_p'] + merged_df[f'V{i}_ap'] #- merged_df[f'V{i}_ap']

# Drop the original columns
merged_df.drop(columns=[f'V{i}_p' for i in range(2, max_Hz+1)] + [f'V{i}_ap' for i in range(2, max_Hz+1)], inplace=True)
merged_df.to_csv(f'{output_dir}/p+ap_alpha_load.csv', index=False)

# calcualte the mean value for each column startswith V by groups_mean
groups_mean = merged_df.groupby(['age', 'load', 'exp'])

# Select columns that start with 'V'
v_columns = [col for col in merged_df.columns if col.startswith('V')]
# Calculate the mean value for each column that starts with 'V'
mean_df = groups_mean[v_columns].mean()
mean_df.to_csv(f'{output_dir}/mean_alpha_load.csv')
# print(mean_df.head())

# Get unique combinations of 'exp' and 'load'
conditions = mean_df.reset_index()[['exp', 'age']].drop_duplicates().values.tolist()

# Create a 2x2 grid of subplots
# Determine the global minimum and maximum values across all your data
global_min = mean_df.min().min()
global_max = mean_df.max().max()
fig, axs = plt.subplots(2, 2, figsize=(18, 18))
axs = axs.flatten()  # Flatten the array of axes to make it easier to iterate over

# Read the head map pic and change it to gray scale
img = Image.open(os.path.join(cwd, 'montage pic', 'alpha.png')).convert('LA')

# Loop over the conditions and plot each condition in a separate subplot
for i, (ax, (exp, age)) in enumerate(zip(axs, conditions)):
    # Select the data for this condition
    data = mean_df.xs((exp, age), level=('exp', 'age'))
    # Determine the age group
    if age == 0:
        age_group = 'Children'
    elif age == 1:
        age_group = 'Adults'
    else:
        age_group = 'Unknown age group'

    # Determine the experiment
    if exp == 'wmc':
        experiment = 'Chinese character experiment'
    elif exp == 'wmv':
        experiment = 'visual pattern experiment'
    else:
        experiment = 'Unknown experiment'
    

    # Loop over the age groups and plot each one
    for load in data.index.unique():
        # Determine the load
        if load == '1b':
            load_n = '1-back'
        elif load == '2b':
            load_n = '2-back'
        else:
            load_n = 'Unknown condition'
        ax.plot(data.loc[load], label=f'{load_n}')
    
    

    # Set the title
    ax.set_title(f'{age_group} {experiment}', fontsize=24)
    ax.set_xlabel('Frequency (Hz)', fontsize=20)
    if i not in [1, 3]:  # Dont plot y title for the right side of the subplot
        ax.set_ylabel(r'Peak Power ($\mu V^{2})$', fontsize=20)
    new_ticks = np.arange(0, max_Hz-1, 2)
    # Generate new labels for these tick marks
    new_labels = new_ticks + 2 # starts from 2 Hz
    # Set the new x-axis tick marks
    ax.set_xticks(new_ticks)
    # Set the tick labels with larger font size
    ax.tick_params(axis='both', which='major', labelsize=14)
    # Set the x-axis labels
    ax.set_xticklabels(new_labels)
    ax.legend(fontsize=18, loc='upper right')
    # Set the same y-axis limits for all subplots
    ax.set_ylim(global_min, global_max)

    # Add an overall y-axis title
    fig.text(0.04, 0.5, 'Occipital-parietal Area', va='center', rotation='vertical', fontsize=28)

    # Get the t, p values for the plot
    off_t = results_df.loc[(results_df['age'] == age) & (results_df['exp'] == exp) & (results_df['column'] == 'offset_diff'), 'T statistic'].values[0]
    off_p = results_df.loc[(results_df['age'] == age) & (results_df['exp'] == exp) & (results_df['column'] == 'offset_diff'), 'P value'].values[0]
    exp_t = results_df.loc[(results_df['age'] == age) & (results_df['exp'] == exp) & (results_df['column'] == 'exponent_diff'), 'T statistic'].values[0]
    exp_p = results_df.loc[(results_df['age'] == age) & (results_df['exp'] == exp) & (results_df['column'] == 'exponent_diff'), 'P value'].values[0]
    cf_t = results_df.loc[(results_df['age'] == age) & (results_df['exp'] == exp) & (results_df['column'] == 'CF_diff'), 'T statistic'].values[0]
    cf_p = results_df.loc[(results_df['age'] == age) & (results_df['exp'] == exp) & (results_df['column'] == 'CF_diff'), 'P value'].values[0]
    pw_t = results_df.loc[(results_df['age'] == age) & (results_df['exp'] == exp) & (results_df['column'] == 'PW_diff'), 'T statistic'].values[0]
    pw_p = results_df.loc[(results_df['age'] == age) & (results_df['exp'] == exp) & (results_df['column'] == 'PW_diff'), 'P value'].values[0]
    bw_t = results_df.loc[(results_df['age'] == age) & (results_df['exp'] == exp) & (results_df['column'] == 'BW_diff'), 'T statistic'].values[0]
    bw_p = results_df.loc[(results_df['age'] == age) & (results_df['exp'] == exp) & (results_df['column'] == 'BW_diff'), 'P value'].values[0]
    
    # add significance mark(*) manually
    if i < 3:
        cell_text = [[r'$\Delta ap_{off}$', round(result.loc[(result['age'] == age) & (result['exp'] == exp), 'offset_diff'].mean(), 2)],
                [r'$\Delta ap_{exp}$', round(result.loc[(result['age'] == age) & (result['exp'] == exp), 'exponent_diff'].mean(), 2)],
                [r'$\Delta \alpha_{cf}$', round(result.loc[(result['age'] == age) & (result['exp'] == exp), 'CF_diff'].mean(), 2)],
                [r'$\Delta \alpha_{pw}$', round(result.loc[(result['age'] == age) & (result['exp'] == exp), 'PW_diff'].mean(), 2)]]
            #  [r'$\Delta \alpha_{bw}$', round(result.loc[(result['age'] == age) & (result['exp'] == exp), 'BW_diff'].values[0], 2)]]
    else:
        cell_text = [[r'$\Delta ap_{off}$', round(result.loc[(result['age'] == age) & (result['exp'] == exp), 'offset_diff'].mean(), 2)],
                [r'$\Delta ap_{exp}$', round(result.loc[(result['age'] == age) & (result['exp'] == exp), 'exponent_diff'].mean(), 2)],
                [r'$\Delta \alpha_{cf}$', round(result.loc[(result['age'] == age) & (result['exp'] == exp), 'CF_diff'].mean(), 2)],
                [r'$\Delta \alpha_{pw}$', str(round(result.loc[(result['age'] == age) & (result['exp'] == exp), 'PW_diff'].mean(), 2)) + '*']]
        
    table = ax.table(cellText=cell_text, loc='bottom left', bbox=[0, 0, 0.2, 0.4])
    # Set the font size of the table
    table.set_fontsize(20)
    
    cells = table.get_celld()
    for cell in cells.values():
        cell.set_text_props(horizontalalignment='center', weight='bold')

    # Center get bg color for the significance
    for key, cell in table.get_celld().items():
        t_list = [off_t, exp_t, cf_t, pw_t, bw_t]
        p_list = [off_p, exp_p, cf_p, pw_p, bw_p]
        for i, (t, p) in enumerate(zip(t_list, p_list)):
            if key[1] == 1 and key[0] == i:  # key[1] is the column index, key[0] is the row index
                if p <= 0.05 and t > 0:
                    cell.set_facecolor('darkorange')
                elif p <= 0.05 and t < 0:
                    cell.set_facecolor('darkturquoise')

    # Add the head map image to each subplot
    first_subplot = True
    for ax in axs.flat:
        if first_subplot: 
            # Create a new axes for the image in the current subplot
            img_ax = ax.inset_axes([0.62, 0.48, 0.35, 0.35])
            img_ax.imshow(img)
            img_ax.axis('off')  # Hide the axis
            first_subplot = False

plt.savefig(os.path.join(output_pic_dir, f'Fig4 alpha.png'), dpi=500, bbox_inches='tight')
# plt.show()

###################################################
# plot the aperiodic power load in Fig. 3 #
###################################################

merged_df = pd.merge(p_power_all, ap_power_all, on=['PID', 'age', 'load', 'exp'], suffixes=('_p', '_ap'))
# Subtract the ap_power_all columns from the raw_data_eeg columns
for i in range(2, max_Hz+1):
    merged_df[f'V{i}'] = merged_df[f'V{i}_p'] + merged_df[f'V{i}_ap'] - merged_df[f'V{i}_p']

# Drop the original columns
merged_df.drop(columns=[f'V{i}_p' for i in range(2, max_Hz+1)] + [f'V{i}_ap' for i in range(2, max_Hz+1)], inplace=True)
# print(merged_df.head())   
merged_df.to_csv(f'{output_dir}/ap_alpha_load.csv', index=False)
# calcualte the mean value for each column startswith V by groups_mean
groups_mean = merged_df.groupby(['age', 'load', 'exp'])

# Select columns that start with 'V'
v_columns = [col for col in merged_df.columns if col.startswith('V')]
# Calculate the mean value for each column that starts with 'V'
mean_df = groups_mean[v_columns].mean()
# Get unique combinations of 'age'
conditions = mean_df.reset_index()['age'].drop_duplicates().values.tolist()
# Loop over the unique PIDs in merged_df
PID_df = merged_df.set_index(['age', 'load', 'exp', 'PID'])
global_min = PID_df.min().min()
global_max = PID_df.max().max()
fig, axs = plt.subplots(1, 2, figsize=(18, 9))
axs = axs.flatten()  # Flatten the array of axes to make it easier to iterate over
# Loop over the conditions and plot each condition in a separate subplot
for i, (ax, (age)) in enumerate(zip(axs, conditions)):
    # Select the data for this age
    data = mean_df.xs(age, level='age')
    # Define the colors
    colors = ['#797BB7', '#E79397', '#80BA8A', '#51B1B7']
    # Create a cycle of colors
    color_cycle = itertools.cycle(colors)
    # Lighten the colors by 70%
    light_colors = [lighten_color(color, 0.3) for color in colors]
    # Create a cycle of light colors
    light_color_cycle = itertools.cycle(light_colors)    
    
    # Each subject data
    PID_data = PID_df.xs(age, level='age')

    # Group by the first two levels of the index
    for (level0_value, level1_value), group in PID_data.groupby(level=[0, 1]):
        # 'group' is a subset of PID_data where the first level of the index equals 'level0_value'
        # and the second level of the index equals 'level1_value'
        # Get the next color in the cycle
        color = next(light_color_cycle)
        for pid in group.index:
            # Plot the data for this group
            ax.plot(group.loc[(pid, slice(None))], color=color, alpha=0.7)

    # Determine the age group
    if age == 0:
        age_group = 'Children'
    elif age == 1:
        age_group = 'Adults'
    else:
        age_group = 'Unknown age'

    for loadnexp in data.index:
        # Determine the load
        if loadnexp[0] == '1b':
            load_n = '1-back'
        elif loadnexp[0] == '2b':
            load_n = '2-back'
        else:
            load_n = 'Unknown condition'

        # Determine the experiment
        if loadnexp[1] == 'wmc':
            experiment = 'Chinese character'
        elif loadnexp[1] == 'wmv':
            experiment = 'visual pattern'
        else:
            experiment = 'Unknown'

        # Get the next color in the cycle
        color = next(color_cycle)
        # Plot the data for this 'exp' and 'load'
        ax.plot(data.loc[(loadnexp, slice(None))], color=color, label=f'{load_n}, {experiment}')

    # Set the title
    ax.set_title(f'{age_group}', fontsize=24)
    ax.set_xlabel('Frequency (Hz)', fontsize=20)
    if i not in [1]:  # Dont plot y title for the right side of the subplot
        ax.set_ylabel(r'Peak Power ($\mu V^{2})$', fontsize=20)

    # Generate the new locations for the x-axis tick marks
    new_ticks = np.arange(0, max_Hz-1, 2)
    # Generate new labels for these tick marks
    new_labels = new_ticks + 2
    # Set the new x-axis tick marks
    ax.set_xticks(new_ticks)
    # Set the tick labels with larger font size
    ax.tick_params(axis='both', which='major', labelsize=14)
    # Set the x-axis labels
    ax.set_xticklabels(new_labels)
    ax.legend(fontsize=14, loc='upper right')


plt.savefig(os.path.join(output_pic_dir, f'Fig3 a.png'), dpi=500, bbox_inches='tight')

###################################################
# plot the periodic power load in Fig. 3 #
###################################################

merged_df = pd.merge(p_power_all, ap_power_all, on=['PID', 'age', 'load', 'exp'], suffixes=('_p', '_ap'))
# Subtract the ap_power_all columns from the raw_data_eeg columns
for i in range(2, max_Hz+1):
    merged_df[f'V{i}'] = merged_df[f'V{i}_p'] + merged_df[f'V{i}_ap'] - merged_df[f'V{i}_ap']

# Drop the original columns
merged_df.drop(columns=[f'V{i}_p' for i in range(2, max_Hz+1)] + [f'V{i}_ap' for i in range(2, max_Hz+1)], inplace=True)
# print(merged_df.head())   
merged_df.to_csv(f'{output_dir}/p_alpha_load.csv', index=False)
# calcualte the mean value for each column startswith V by groups_mean
groups_mean = merged_df.groupby(['age', 'load', 'exp'])

# Select columns that start with 'V'
v_columns = [col for col in merged_df.columns if col.startswith('V')]
# Calculate the mean value for each column that starts with 'V'
mean_df = groups_mean[v_columns].mean()
# Get unique combinations of 'age'
conditions = mean_df.reset_index()['age'].drop_duplicates().values.tolist()
# Loop over the unique PIDs in merged_df
PID_df = merged_df.set_index(['age', 'load', 'exp', 'PID'])
global_min = PID_df.min().min()
global_max = PID_df.max().max()
fig, axs = plt.subplots(1, 2, figsize=(18, 9))
axs = axs.flatten()  # Flatten the array of axes to make it easier to iterate over
# Loop over the conditions and plot each condition in a separate subplot
for i, (ax, (age)) in enumerate(zip(axs, conditions)):
    # Select the data for this age
    data = mean_df.xs(age, level='age')
    # Define the colors
    colors = ['#797BB7', '#E79397', '#80BA8A', '#51B1B7']
    # Create a cycle of colors
    color_cycle = itertools.cycle(colors)
    # Lighten the colors by 70%
    light_colors = [lighten_color(color, 0.3) for color in colors]
    # Create a cycle of light colors
    light_color_cycle = itertools.cycle(light_colors)    
    
    # Each subject data
    PID_data = PID_df.xs(age, level='age')

    # Group by the first two levels of the index
    for (level0_value, level1_value), group in PID_data.groupby(level=[0, 1]):
        # 'group' is a subset of PID_data where the first level of the index equals 'level0_value'
        # and the second level of the index equals 'level1_value'
        # Get the next color in the cycle
        color = next(light_color_cycle)
        for pid in group.index:
            # Plot the data for this group
            ax.plot(group.loc[(pid, slice(None))], color=color, alpha=0.7)

    # Determine the age group
    if age == 0:
        age_group = 'Children'
    elif age == 1:
        age_group = 'Adults'
    else:
        age_group = 'Unknown age'

    for loadnexp in data.index:
        # Determine the load
        if loadnexp[0] == '1b':
            load_n = '1-back'
        elif loadnexp[0] == '2b':
            load_n = '2-back'
        else:
            load_n = 'Unknown condition'

        # Determine the experiment
        if loadnexp[1] == 'wmc':
            experiment = 'Chinese character'
        elif loadnexp[1] == 'wmv':
            experiment = 'visual pattern'
        else:
            experiment = 'Unknown'

        # Get the next color in the cycle
        color = next(color_cycle)
        # Plot the data for this 'exp' and 'load'
        ax.plot(data.loc[(loadnexp, slice(None))], color=color, label=f'{load_n}, {experiment}')

    # Set the title
    ax.set_title(f'{age_group}', fontsize=24)
    ax.set_xlabel('Frequency (Hz)', fontsize=20)
    if i not in [1]:  # Dont plot y title for the right side of the subplot
        ax.set_ylabel(r'Peak Power ($\mu V^{2})$', fontsize=20)

    # Generate the new locations for the x-axis tick marks
    new_ticks = np.arange(0, max_Hz-1, 2)
    # Generate new labels for these tick marks
    new_labels = new_ticks + 2
    # Set the new x-axis tick marks
    ax.set_xticks(new_ticks)
    # Set the tick labels with larger font size
    ax.tick_params(axis='both', which='major', labelsize=14)
    # Set the x-axis labels
    ax.set_xticklabels(new_labels)
    ax.legend(fontsize=14, loc='upper right')
    # Set the tick labels with larger font size
    ax.tick_params(axis='both', which='major', labelsize=14)
    # Set the x-axis labels
    ax.set_xticklabels(new_labels)
    ax.legend(fontsize=14, loc='upper right')


plt.savefig(os.path.join(output_pic_dir, f'Fig3 g.png'), dpi=500, bbox_inches='tight')
