##### First plot #####

#############################################################################
#This script hold whole dataset including the reads with or without tag-jump#
#############################################################################

#load the required libraries
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

# Create dataframes from the provided data
uchime_denovo_data = {
    "FILE": ["ERR6454461.fasta", "ERR6454462.fasta", "ERR6454463.fasta", "ERR6454464.fasta", "ERR6454465.fasta",
             "ERR6454466.fasta", "ERR6454467.fasta", "ERR6454468.fasta", "ERR6454469.fasta", "ERR6454470.fasta",
             "ERR6454471.fasta", "ERR6454472.fasta", "ERR6454473.fasta", "ERR6454474.fasta", "ERR6454475.fasta",
             "ERR6454476.fasta", "ERR6454477.fasta", "ERR6454478.fasta"],
    "Non-chimeric reads without tag-jump": [176535, 87889, 87646, 261847, 349319, 111274, 146681, 361163, 360439, 998035, 160357, 95272, 273924, 245273, 372754, 313935, 334612, 94765],
    "Non-chimeric reads with tag-jump": [176255, 63595, 65824, 235319, 289378, 104501, 110319, 324860, 252619, 996721, 152601, 68418, 211543, 210658, 192483, 193727, 216007, 54275],
    "Chimeric reads without tag-jump": [1286, 339, 319, 2518, 2354, 957, 951, 3400, 2939, 10220, 3227, 1321, 25431, 1642, 2779, 5239, 4714, 526],
    "Chimeric reads with tag-jump": [1285, 316, 367, 3732, 1099, 884, 1454, 3122, 3004, 10218, 3077, 1482, 25588, 1483, 3012, 4807, 4768, 619],
    "BLASTn recovery without tag-jump": [422, 139, 19, 1262, 1277, 538, 639, 2679, 2568, 7046, 1104, 317, 701, 133, 183, 107, 622, 29],
    "BLASTn recovery with tag-jump": [421, 111, 110, 2495, 714, 489, 694, 2372, 2611, 7045, 1124, 320, 840, 611, 455, 146, 356, 45]
}

chimeras_denovo_data = {
    "FILE": ["ERR6454461.fasta", "ERR6454462.fasta", "ERR6454463.fasta", "ERR6454464.fasta", "ERR6454465.fasta",
             "ERR6454466.fasta", "ERR6454467.fasta", "ERR6454468.fasta", "ERR6454469.fasta", "ERR6454470.fasta",
             "ERR6454471.fasta", "ERR6454472.fasta", "ERR6454473.fasta", "ERR6454474.fasta", "ERR6454475.fasta",
             "ERR6454476.fasta", "ERR6454477.fasta", "ERR6454478.fasta"],
    "Non-chimeric reads without tag-jump": [171203, 84490, 82366, 255591, 331088, 108861, 138419, 347544, 341398, 897334, 159342, 89266, 272363, 236590, 348272, 295857, 321177, 86558],
    "Non-chimeric reads with tag-jump": [170939, 60830, 61212, 231143, 275937, 102035, 104817, 311684, 237071, 896021, 151542, 64185, 213087, 202583, 177376, 187675, 206184, 48076],
    "Chimeric reads without tag-jump": [6618, 3738, 5599, 8774, 20585, 3370, 9213, 17019, 21980, 110921, 4242, 7327, 26992, 10325, 27261, 23317, 18149, 8733],
    "Chimeric reads with tag-jump": [6601, 3081, 4979, 7908, 14540, 3350, 6956, 16298, 18552, 110918, 4136, 5715, 24044, 9558, 18119, 10859, 14591, 6818],
    "BLASTn recovery without tag-jump": [3870, 3415, 3081, 6544, 13470, 2129, 6334, 14657, 21798, 43038, 805, 176, 6212, 1648, 472, 16594, 1172, 7619],
    "BLASTn recovery with tag-jump": [3856, 2776, 2973, 5946, 9876, 2127, 4369, 14001, 18376, 43035, 803, 156, 5787, 1504, 186, 6995, 186, 5940]
}

dada2_data = {
    "FILE": ["ERR6454461.fasta", "ERR6454462.fasta", "ERR6454463.fasta", "ERR6454464.fasta", "ERR6454465.fasta",
             "ERR6454466.fasta", "ERR6454467.fasta", "ERR6454468.fasta", "ERR6454469.fasta", "ERR6454470.fasta",
             "ERR6454471.fasta", "ERR6454472.fasta", "ERR6454473.fasta", "ERR6454474.fasta", "ERR6454475.fasta",
             "ERR6454476.fasta", "ERR6454477.fasta", "ERR6454478.fasta"],
    "Non-chimeric reads without tag-jump": [161791, 80051, 77648, 242208, 305617, 104112, 130522, 332200, 321809, 790551, 152390, 81329, 243909, 224077, 316278, 275414, 296732, 79453],
    "Non-chimeric reads with tag-jump": [252601, 90955, 94351, 365009, 439636, 164875, 166830, 508833, 374654, 966078, 229647, 87814, 306617, 259678, 247668, 298874, 311622, 71099],
    "Chimeric reads without tag-jump": [16030, 8177, 10317, 22157, 44339, 8119, 17110, 32363, 41569, 217704, 11194, 15264, 55446, 22838, 59255, 43760, 42594, 15838],
    "Chimeric reads with tag-jump": [20593, 10245, 13502, 26106, 45094, 9758, 19784, 45443, 51898, 303374, 15053, 18272, 63023, 29093, 57004, 33134, 43393, 18486],
    "BLASTn recovery without tag-jump": [8776, 7114, 5675, 16003, 30668, 5237, 11927, 27933, 41160, 102566, 2215, 556, 9864, 4496, 1666, 28578, 2833, 13035],
    "BLASTn recovery with tag-jump": [11333, 9061, 8372, 19118, 31403, 6806, 12555, 38963, 51350, 137185, 2676, 599, 14183, 5003, 741, 19305, 787, 15641]
}

# Convert to pandas dataframes
df_uchime_denovo = pd.DataFrame(uchime_denovo_data)
df_chimeras_denovo = pd.DataFrame(chimeras_denovo_data)
df_dada2 = pd.DataFrame(dada2_data)

# Combine dataframes for easier plotting
df_uchime_denovo['Method'] = 'Uchime Denovo'
df_chimeras_denovo['Method'] = 'Chimeras Denovo'
df_dada2['Method'] = 'DADA2'

combined_df = pd.concat([df_uchime_denovo, df_chimeras_denovo, df_dada2])

# Define plotting function
def plot_multi_panel_line_plot(data):
    methods = data['Method'].unique()
    categories = [
        'Non-chimeric reads without tag-jump',
        'Non-chimeric reads with tag-jump',
        'Chimeric reads without tag-jump',
        'Chimeric reads with tag-jump',
        'BLASTn recovery without tag-jump',
        'BLASTn recovery with tag-jump'
    ]
    colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown']
    linestyles = ['-', '--', '-', '--', '-', '--']

    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(14, 18), sharex=True)
    axes = axes.flatten()

    for i, method in enumerate(methods):
        ax = axes[i]
        method_data = data[data['Method'] == method]
        for j, category in enumerate(categories):
            linestyle = '--' if 'without tag-jump' in category else '-'
            ax.plot(method_data['FILE'], method_data[category], label=category, color=colors[j], linestyle=linestyle)
        
        ax.set_title(method)
        ax.set_ylabel('Number of Reads')
        ax.legend(loc='upper left')
        ax.grid(True)
        ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: f'{int(x):,}'))

    plt.xlabel('File')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

# Plot the multi-panel line plot
plot_multi_panel_line_plot(combined_df)


######################## Second plot #################################

########################################################################
#This script designed for the dataset including the reads with tag-jump#
########################################################################

# Load the required libraries
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

# Create dataframes from the provided data
uchime_denovo_data = {
    "FILE": ["ERR6454461.fasta", "ERR6454462.fasta", "ERR6454463.fasta", "ERR6454464.fasta", "ERR6454465.fasta",
             "ERR6454466.fasta", "ERR6454467.fasta", "ERR6454468.fasta", "ERR6454469.fasta", "ERR6454470.fasta",
             "ERR6454471.fasta", "ERR6454472.fasta", "ERR6454473.fasta", "ERR6454474.fasta", "ERR6454475.fasta",
             "ERR6454476.fasta", "ERR6454477.fasta", "ERR6454478.fasta"],
    "Non-chimeric reads": [176255, 63595, 65824, 235319, 289378, 104501, 110319, 324860, 252619, 996721, 152601, 68418, 211543, 210658, 192483, 193727, 216007, 54275],
    "Chimeric reads": [1285, 316, 367, 3732, 1099, 884, 1454, 3122, 3004, 10218, 3077, 1482, 25588, 1483, 3012, 4807, 4768, 619],
    "BLASTn recovery": [421, 111, 110, 2495, 714, 489, 694, 2372, 2611, 7045, 1124, 320, 840, 611, 455, 146, 356, 45]
}

chimeras_denovo_data = {
    "FILE": ["ERR6454461.fasta", "ERR6454462.fasta", "ERR6454463.fasta", "ERR6454464.fasta", "ERR6454465.fasta",
             "ERR6454466.fasta", "ERR6454467.fasta", "ERR6454468.fasta", "ERR6454469.fasta", "ERR6454470.fasta",
             "ERR6454471.fasta", "ERR6454472.fasta", "ERR6454473.fasta", "ERR6454474.fasta", "ERR6454475.fasta",
             "ERR6454476.fasta", "ERR6454477.fasta", "ERR6454478.fasta"],
    "Non-chimeric reads": [170939, 60830, 61212, 231143, 275937, 102035, 104817, 311684, 237071, 896021, 151542, 64185, 213087, 202583, 177376, 187675, 206184, 48076],
    "Chimeric reads": [6601, 3081, 4979, 7908, 14540, 3350, 6956, 16298, 18552, 110918, 4136, 5715, 24044, 9558, 18119, 10859, 14591, 6818],
    "BLASTn recovery": [3856, 2776, 2973, 5946, 9876, 2127, 4369, 14001, 18376, 43035, 803, 156, 5787, 1504, 186, 6995, 186, 5940]
}

dada2_data = {
    "FILE": ["ERR6454461.fasta", "ERR6454462.fasta", "ERR6454463.fasta", "ERR6454464.fasta", "ERR6454465.fasta",
             "ERR6454466.fasta", "ERR6454467.fasta", "ERR6454468.fasta", "ERR6454469.fasta", "ERR6454470.fasta",
             "ERR6454471.fasta", "ERR6454472.fasta", "ERR6454473.fasta", "ERR6454474.fasta", "ERR6454475.fasta",
             "ERR6454476.fasta", "ERR6454477.fasta", "ERR6454478.fasta"],
    "Non-chimeric reads": [252601, 90955, 94351, 365009, 439636, 164875, 166830, 508833, 374654, 966078, 229647, 87814, 306617, 259678, 247668, 298874, 311622, 71099],
    "Chimeric reads": [20593, 10245, 13502, 26106, 45094, 9758, 19784, 45443, 51898, 303374, 15053, 18272, 63023, 29093, 57004, 33134, 43393, 18486],
    "BLASTn recovery": [11333, 9061, 8372, 19118, 31403, 6806, 12555, 38963, 51350, 137185, 2676, 599, 14183, 5003, 741, 19305, 787, 15641]
}

# Convert to pandas dataframes
df_uchime_denovo = pd.DataFrame(uchime_denovo_data)
df_chimeras_denovo = pd.DataFrame(chimeras_denovo_data)
df_dada2 = pd.DataFrame(dada2_data)

# Combine dataframes for easier plotting
df_uchime_denovo['Method'] = 'Uchime Denovo'
df_chimeras_denovo['Method'] = 'Chimeras Denovo'
df_dada2['Method'] = 'DADA2'

combined_df = pd.concat([df_uchime_denovo, df_chimeras_denovo, df_dada2])

# Define plotting function
def plot_multi_panel_line_plot(data):
    methods = data['Method'].unique()
    categories = [
        'Non-chimeric reads',
        'Chimeric reads',
        'BLASTn recovery'
    ]
    colors = ['blue', 'green', 'purple']
    linestyles = ['-', '--', '-']

    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(14, 18), sharex=True)
    axes = axes.flatten()

    for i, method in enumerate(methods):
        ax = axes[i]
        method_data = data[data['Method'] == method]
        for j, category in enumerate(categories):
            linestyle = linestyles[j]
            ax.plot(method_data['FILE'], method_data[category], label=category, color=colors[j], linestyle=linestyle)
        
        ax.set_title(method)
        ax.set_ylabel('Number of Reads')
        ax.legend(loc='upper left')
        ax.grid(True)
        ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: f'{int(x):,}'))

    plt.xlabel('File')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

# Plot the multi-panel line plot
plot_multi_panel_line_plot(combined_df)


############################## THIRD PLOT #########################

#############################################################################################################################
#This script is related to the queries which tag-jump filtering applied on them and to show logaritmic value embeded in them#
#############################################################################################################################

#load the required libraries
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

# Create dataframes from the provided data
uchime_denovo_data = {
    "FILE": ["ERR6454461.fasta", "ERR6454462.fasta", "ERR6454463.fasta", "ERR6454464.fasta", "ERR6454465.fasta",
             "ERR6454466.fasta", "ERR6454467.fasta", "ERR6454468.fasta", "ERR6454469.fasta", "ERR6454470.fasta",
             "ERR6454471.fasta", "ERR6454472.fasta", "ERR6454473.fasta", "ERR6454474.fasta", "ERR6454475.fasta",
             "ERR6454476.fasta", "ERR6454477.fasta", "ERR6454478.fasta"],
    "Non-chimeric reads": [176255, 63595, 65824, 235319, 289378, 104501, 110319, 324860, 252619, 996721, 152601, 68418, 211543, 210658, 192483, 193727, 216007, 54275],
    "Chimeric reads": [1285, 316, 367, 3732, 1099, 884, 1454, 3122, 3004, 10218, 3077, 1482, 25588, 1483, 3012, 4807, 4768, 619],
    "BLASTn recovered reads": [421, 111, 110, 2495, 714, 489, 694, 2372, 2611, 7045, 1124, 320, 840, 611, 455, 146, 356, 45]
}

chimeras_denovo_data = {
    "FILE": ["ERR6454461.fasta", "ERR6454462.fasta", "ERR6454463.fasta", "ERR6454464.fasta", "ERR6454465.fasta",
             "ERR6454466.fasta", "ERR6454467.fasta", "ERR6454468.fasta", "ERR6454469.fasta", "ERR6454470.fasta",
             "ERR6454471.fasta", "ERR6454472.fasta", "ERR6454473.fasta", "ERR6454474.fasta", "ERR6454475.fasta",
             "ERR6454476.fasta", "ERR6454477.fasta", "ERR6454478.fasta"],
    "Non-chimeric reads": [170939, 60830, 61212, 231143, 275937, 102035, 104817, 311684, 237071, 896021, 151542, 64185, 213087, 202583, 177376, 187675, 206184, 48076],
    "Chimeric reads": [6601, 3081, 4979, 7908, 14540, 3350, 6956, 16298, 18552, 110918, 4136, 5715, 24044, 9558, 18119, 10859, 14591, 6818],
    "BLASTn recovered reads": [3856, 2776, 2973, 5946, 9876, 2127, 4369, 14001, 18376, 43035, 803, 156, 5787, 1504, 186, 6995, 186, 5940]
}

dada2_data = {
    "FILE": ["ERR6454461.fasta", "ERR6454462.fasta", "ERR6454463.fasta", "ERR6454464.fasta", "ERR6454465.fasta",
             "ERR6454466.fasta", "ERR6454467.fasta", "ERR6454468.fasta", "ERR6454469.fasta", "ERR6454470.fasta",
             "ERR6454471.fasta", "ERR6454472.fasta", "ERR6454473.fasta", "ERR6454474.fasta", "ERR6454475.fasta",
             "ERR6454476.fasta", "ERR6454477.fasta", "ERR6454478.fasta"],
    "Non-chimeric reads": [252601, 90955, 94351, 365009, 439636, 164875, 166830, 508833, 374654, 966078, 229647, 87814, 306617, 259678, 247668, 298874, 311622, 71099],
    "Chimeric reads": [20593, 10245, 13502, 26106, 45094, 9758, 19784, 45443, 51898, 303374, 15053, 18272, 63023, 29093, 57004, 33134, 43393, 18486],
    "BLASTn recovered reads": [11333, 9061, 8372, 19118, 31403, 6806, 12555, 38963, 51350, 137185, 2676, 599, 14183, 5003, 741, 19305, 787, 15641]
}

# Convert to pandas dataframes
df_uchime_denovo = pd.DataFrame(uchime_denovo_data)
df_chimeras_denovo = pd.DataFrame(chimeras_denovo_data)
df_dada2 = pd.DataFrame(dada2_data)

# Combine dataframes for easier plotting
df_uchime_denovo['Method'] = 'Uchime Denovo'
df_chimeras_denovo['Method'] = 'Chimeras Denovo'
df_dada2['Method'] = 'DADA2'

combined_df = pd.concat([df_uchime_denovo, df_chimeras_denovo, df_dada2])

# Convert abundances to logarithmic scale
for col in combined_df.columns:
    if col not in ['FILE', 'Method']:
        combined_df[col] = combined_df[col].apply(lambda x: np.log10(x + 1))

# Define plotting function
def plot_multi_panel_line_plot(data):
    methods = data['Method'].unique()
    categories = [
        'Non-chimeric reads',
        'Chimeric reads',
        'BLASTn recovered reads'
    ]
    colors = ['blue', 'orange', 'green']
    linestyles = ['-', '--', '-']

    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(14, 18), sharex=True)
    axes = axes.flatten()

    for i, method in enumerate(methods):
        ax = axes[i]
        method_data = data[data['Method'] == method]
        for j, category in enumerate(categories):
            linestyle = linestyles[j]
            ax.plot(method_data['FILE'], method_data[category], label=category, color=colors[j], linestyle=linestyle)
        
        ax.set_title(method)
        ax.set_ylabel('Log10(Number of Reads + 1)')
        ax.legend(loc='upper left')
        ax.grid(True)

    plt.xlabel('File')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

# Plot the multi-panel line plot
plot_multi_panel_line_plot(combined_df)
