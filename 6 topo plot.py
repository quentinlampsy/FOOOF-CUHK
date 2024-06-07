import os, shutil
import numpy as np
import pandas as pd
import mne
import matplotlib.pyplot as plt
from PIL import Image
from sklearn.preprocessing import StandardScaler

cwd = os.path.dirname(os.path.realpath(__file__))
os.chdir(cwd)

data_dir = os.path.join(cwd, '5 topo cleaned')

output_dir_alpha = os.path.join(cwd, '6 alpha topo')
# Check if output_dir_alpha is not empty and it exists
if output_dir_alpha and os.path.exists(output_dir_alpha):
    # Remove the directory and its contents
    if os.path.isdir(output_dir_alpha):
        shutil.rmtree(output_dir_alpha)
    else:
        os.remove(output_dir_alpha)
# Recreate the directory
os.makedirs(output_dir_alpha, exist_ok=True)

output_dir_theta = os.path.join(cwd, '6 theta topo')
# Check if output_dir_theta is not empty and it exists
if output_dir_theta and os.path.exists(output_dir_theta):
    # Remove the directory and its contents
    if os.path.isdir(output_dir_theta):
        shutil.rmtree(output_dir_theta)
    else:
        os.remove(output_dir_theta)
# Recreate the directory
os.makedirs(output_dir_theta, exist_ok=True)

# Specify the path to the montage file
montage_file_path = os.path.join(cwd, os.path.join('montage files', 'GSN-HydroCel-128.sfp'))

# The head of our montage file looks like this:
# E1	5.787677636	5.520863216	-2.577468644
# E2	5.291804727	6.709097557	0.307434896
# E3	3.864122447	7.63424051	3.067770143
# E4	2.868837559	7.145708546	4.989564557
# E5	1.479340453	5.68662139	6.812878187

# Read the file into a DataFrame
montage_df = pd.read_csv(montage_file_path, sep='\t', header=None, names=['ch_name', 'x', 'y', 'z'])

# Extract the channel names and coordinates
ch_names = montage_df['ch_name'].tolist()
x = montage_df['x'].values
y = montage_df['y'].values
z = montage_df['z'].values

# # Create a DigMontage
montage = mne.channels.make_dig_montage(dict(zip(ch_names, zip(x, y, z))), coord_frame='head')


# function to fit eeg data to montage
def montage_fit(data):
    # Convert the DataFrame to a 3D numpy array
    data = data.loc[:, 'EEG1': 'EEG128'].values[:, np.newaxis]

    # Initialize a new StandardScaler instance
    scaler = StandardScaler()

    # # Fit the scaler to the data and transform the data for normalizing data by rows
    # data = data.squeeze()  # remove singleton dimensions
    # data_zscore = scaler.fit_transform(data.T)

    # # Add an extra dimension to match the original data shape
    # data_zscore = data_zscore.T[:, np.newaxis]

    # Create an Info object for the EEG data
    info = mne.create_info(montage.ch_names, sfreq=1.0, ch_types='eeg')

    # # Average the data across trials (the first dimension)
    # avg_data = np.mean(data_zscore, axis=0)
    avg_data = np.mean(data, axis=0)

    # Reshape avg_data to (128, 1)
    avg_data = np.reshape(avg_data, (128, 1))

    # Create an Evoked object for the EEG data
    fits = mne.EvokedArray(avg_data, info)

    # Apply the montage to the Raw object
    fits.set_montage(montage)

    return fits


def plot_topo(fits, group, exp, load, i):
    # reshape to an 1D array
    fits_1D = np.squeeze(fits.data)
    # set up the figure
    fig, ax = plt.subplots()
    im, _ = mne.viz.plot_topomap(fits_1D, fits.info, cmap='RdBu_r', axes=ax, show=False)
    if group == 'adult':
        group = 'Adults'
    elif group == 'child':
        group = 'Children'
    if exp == 'wmc':
        exp = 'CC'
    elif exp == 'wmv':
        exp = 'VP'
    
    title = f'{group} {exp} {load}'
    ax.set_title(title, fontsize=24)

    print(group, load, exp, i)
    # Add a colorbar to the right of the plot
    if i in [4, 10, 11]:
        fig.colorbar(im, ax=ax)

    return fig, title


# List comprehension to get all CSV files with '_PW' in the filename
csv_files = [f for f in os.listdir(data_dir) if f.endswith('.csv') and '_PW' in f]

i = 0
for f in csv_files:
    if i >= 1:
        group_temp = group
        exp_temp = exp
        load_temp = load
        # fits_alpha_temp = fits_alpha
        # fits_theta_temp = fits_theta
        data_alpha_temp = data_alpha
        data_theta_temp = data_theta


    group = f.split('_')[0]
    exp = f.split('_')[1]
    load = f.split('_')[2]
    # print(group, exp, load) # adult wmc 1b

    data = pd.read_csv(os.path.join(data_dir, f))
    data_alpha = data.loc[data['band'] == 'alpha', :]
    data_theta = data.loc[data['band'] == 'theta', :]

    fits_alpha = montage_fit(data_alpha)
    fits_theta = montage_fit(data_theta)

    fig_alpha, title_alpha = plot_topo(fits_alpha, group, exp, load, i)
    fig_theta, title_theta = plot_topo(fits_theta, group, exp, load, i)

    # Save the figure
    fig_alpha.savefig(os.path.join(output_dir_alpha, f'{title_alpha}.png'))
    fig_theta.savefig(os.path.join(output_dir_theta, f'{title_theta}.png'))

    plt.close(fig_alpha)
    plt.close(fig_theta)
    
    i += 1
    if i >= 2:
        # first minus then average, this is for the difference plot
        if group == group_temp and exp == exp_temp:
            data_delta_load_alpha = data_alpha.loc[:, 'EEG1':].values - data_alpha_temp.loc[:, 'EEG1':].values
            data_delta_load_theta = data_theta.loc[:, 'EEG1':].values - data_theta_temp.loc[:, 'EEG1':].values

            fits_delta_load_alpha = montage_fit(pd.DataFrame(data_delta_load_alpha, columns=data_alpha.columns[3:]))
            fits_delta_load_theta = montage_fit(pd.DataFrame(data_delta_load_theta, columns=data_theta.columns[3:]))

            fig_delta_load_alpha, title_delta_load_alpha = plot_topo(fits_delta_load_alpha, group, exp, f'{load} - {load_temp}', i)
            fig_delta_load_theta, title_delta_load_theta = plot_topo(fits_delta_load_theta, group, exp, f'{load} - {load_temp}', i)

            fig_delta_load_alpha.savefig(os.path.join(output_dir_alpha, f'{title_delta_load_alpha}.png'))
            fig_delta_load_theta.savefig(os.path.join(output_dir_theta, f'{title_delta_load_theta}.png'))

            plt.close(fig_alpha)
            plt.close(fig_theta)

            i += 1
        # # first average then minus, this is for the difference plot
        # if group == group_temp and exp == exp_temp:
        #     fits_delta_load_alpha_data = fits_alpha.data - fits_alpha_temp.data
        #     fits_delta_load_theta_data = fits_theta.data - fits_theta_temp.data

        #     # Create new EvokedArray objects with the subtracted data
        #     fits_delta_load_alpha = mne.EvokedArray(fits_delta_load_alpha_data, fits_alpha.info)
        #     fits_delta_load_theta = mne.EvokedArray(fits_delta_load_theta_data, fits_theta.info)

        #     fig_delta_load_alpha, title_delta_load_alpha = plot_topo(fits_delta_load_alpha, group, exp, f'{load} - {load_temp}')
        #     fig_delta_load_theta, title_delta_load_theta = plot_topo(fits_delta_load_theta, group, exp, f'{load} - {load_temp}')

        #     fig_delta_load_alpha.savefig(os.path.join(output_dir_alpha, f'{title_delta_load_alpha}.png'))
        #     fig_delta_load_theta.savefig(os.path.join(output_dir_theta, f'{title_delta_load_theta}.png'))


# Get all PNG files in the directories and sort them
alpha_pics = sorted([f for f in os.listdir(output_dir_alpha) if f.endswith('.png')], key=lambda x: (len(x), x))
theta_pics = sorted([f for f in os.listdir(output_dir_theta) if f.endswith('.png')], key=lambda x: (len(x), x))

# Create a 3-row, 4-column subplot structure for alpha images
fig, axs = plt.subplots(3, 4, figsize=(15, 15))

# Plot all alpha images
for i, f in enumerate(alpha_pics):
    img = Image.open(os.path.join(output_dir_alpha, f))
    row = i // 4  # Determine row index
    col = i % 4  # Determine column index
    axs[row, col].imshow(img)
    axs[row, col].axis('off')  # Hide axes

# Adjust layout to be tighter
plt.tight_layout()

# Add a title to the joint plot
plt.suptitle('Topographies for alpha band in each condition each group', fontsize=24)

# Save the alpha joint plot
plt.savefig(os.path.join(output_dir_alpha, 'Alpha joint.png'))
plt.close(fig)
# Clear the current figure to start a new one for theta
plt.clf()

# Create a new 3-row, 4-column subplot structure for theta images
fig, axs = plt.subplots(3, 4, figsize=(15, 15))

# Plot all theta images
for i, f in enumerate(theta_pics):
    img = Image.open(os.path.join(output_dir_theta, f))
    row = i // 4  # Determine row index
    col = i % 4  # Determine column index
    axs[row, col].imshow(img)
    axs[row, col].axis('off')  # Hide axes

# Adjust layout to be tighter
plt.tight_layout()

# Add a title to the joint plot
plt.suptitle('Topographies for theta band in each condition each group', fontsize=24)

# Save the theta joint plot
plt.savefig(os.path.join(output_dir_theta, 'Theta joint.png'))
plt.close(fig)
# Clear the current figure
plt.clf()

# Plot the sensor locations with channel names
fig, ax = plt.subplots()
mne.viz.plot_sensors(fits_alpha.info, show_names=True, axes=ax)
ax.set_title('Sensor Locations with Channel Names')
fig.savefig(os.path.join('montage files', 'channel location.png'))
plt.close(fig)
# plt.show()

