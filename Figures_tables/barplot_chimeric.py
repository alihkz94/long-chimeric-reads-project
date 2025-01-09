"""
This script generates stacked bar plots for the number of chimeric 
and false positive chimeric sequences detected by different chimera 
detection modules.
"""

#importlibraries
import pandas as pd
import matplotlib.pyplot as plt
import os

# Optional: Enhance plot aesthetics with seaborn
import seaborn as sns
sns.set(style="whitegrid")

# Create output directory for the images
output_dir = 'chimera_detection_plots_separate'
os.makedirs(output_dir, exist_ok=True)

# Function to create dataframes for each module
def create_module_dataframe(module_name):
    if module_name == "Chimeras_denovo":
        data = {
            'Sample ID': [
                'ERR6454461', 'ERR6454462', 'ERR6454463', 'ERR6454464',
                'ERR6454465', 'ERR6454466', 'ERR6454467', 'ERR6454468',
                'ERR6454469', 'ERR6454470', 'ERR6454471', 'ERR6454472',
                'ERR6454473', 'ERR6454474', 'ERR6454475', 'ERR6454476',
                'ERR6454477', 'ERR6454478'
            ],
            'Chimeric': [
                1161, 422, 1207, 2043, 4558, 1015, 1190, 2960,
                2448, 36708, 1669, 3335, 10027, 2416, 6873, 4458,
                5046, 2176
            ],
            'False positive chimeras': [
                4214, 2314, 2935, 4435, 8435, 1663, 5203, 11993,
                14463, 30170, 1823, 1362, 5470, 5253, 8639, 4213,
                6099, 3279
            ]
        }
    elif module_name == "UCHIME_denovo":
        data = {
            'Sample ID': [
                'ERR6454461', 'ERR6454462', 'ERR6454463', 'ERR6454464',
                'ERR6454465', 'ERR6454466', 'ERR6454467', 'ERR6454468',
                'ERR6454469', 'ERR6454470', 'ERR6454471', 'ERR6454472',
                'ERR6454473', 'ERR6454474', 'ERR6454475', 'ERR6454476',
                'ERR6454477', 'ERR6454478'
            ],
            'Chimeric': [
                855, 185, 216, 1278, 428, 430, 688, 1539,
                1479, 5231, 1638, 686, 14253, 328, 1873, 3092,
                1915, 452
            ],
            'False positive chimeras': [
                45, 107, 93, 316, 381, 86, 542, 240,
                597, 333, 276, 162, 186, 186, 267, 330,
                359, 15
            ]
        }
    elif module_name == "DADA2":
        data = {
            'Sample ID': [
                'ERR6454461', 'ERR6454462', 'ERR6454463', 'ERR6454464',
                'ERR6454465', 'ERR6454466', 'ERR6454467', 'ERR6454468',
                'ERR6454469', 'ERR6454470', 'ERR6454471', 'ERR6454472',
                'ERR6454473', 'ERR6454474', 'ERR6454475', 'ERR6454476',
                'ERR6454477', 'ERR6454478'
            ],
            'Chimeric': [
                1985, 775, 1483, 4158, 7444, 1507, 1935, 5152,
                3299, 59982, 3054, 5072, 15867, 4118, 10288, 6762,
                8379, 2771
            ],
            'False positive chimeras': [
                6810, 3660, 4356, 7265, 12579, 2766, 7153, 15863,
                20093, 40793, 3339, 2384, 7574, 7715, 12758, 7062,
                8947, 4884
            ]
        }
    else:
        raise ValueError("Invalid module name.")
    
    df = pd.DataFrame(data)
    return df

# Function to create and save a stacked bar chart for a given module
def create_and_save_plot(module_name, df):
    # Define colors
    colors = {
        'Chimeric': '#D55E00',              # Vermillion
        'False positive chimeras': '#009E73' # Bluish green
    }
    
    # Set figure size close to A4 dimensions (8.27 x 11.69 inches)
    fig, ax = plt.subplots(figsize=(11.69, 8.27))  # Landscape orientation for better readability
    
    # Create stacked bar chart
    chimeric_bars = ax.bar(df['Sample ID'], df['Chimeric'], label='Chimeric', color=colors['Chimeric'])
    false_pos_bars = ax.bar(df['Sample ID'], df['False positive chimeras'], bottom=df['Chimeric'],
                            label='False positive chimeras', color=colors['False positive chimeras'])
    
    # Customize the plot
    ax.set_title(f'Chimera Detection Results - {module_name}', fontsize=20, fontweight='bold')
    ax.set_xlabel('Sample ID', fontsize=16, fontweight='bold')
    ax.set_ylabel('Number of Sequences', fontsize=16, fontweight='bold')
    
    # Customize tick labels to be bold and rotated
    ax.set_xticklabels(df['Sample ID'], rotation=45, ha='right', fontsize=12, fontweight='bold')
    ax.tick_params(axis='y', labelsize=14)
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')
    
    # Customize legend labels to be bold
    legend = ax.legend(fontsize=14, loc='upper right')
    for text in legend.get_texts():
        text.set_fontweight('bold')
    
    # Adjust layout to prevent clipping of tick-labels
    plt.tight_layout()
    
    # Save the figure with high resolution
    output_path = os.path.join(output_dir, f'{module_name}_Chimera_Detection.png')
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close()
    
    print(f"Plot for '{module_name}' has been saved as '{output_path}' with 600 DPI resolution.")

# Define modules
modules = ["UCHIME_denovo", "Chimeras_denovo", "DADA2"]

# Generate and save plots for each module
for module in modules:
    df = create_module_dataframe(module)
    create_and_save_plot(module, df)