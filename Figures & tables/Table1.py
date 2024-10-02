import pandas as pd
import imgkit

# Data
data = {
    "Module": ["UCHIME_denovo", "UCHIME_denovo", "Chimeras_denovo", "Chimeras_denovo", "DADA2"],
    "Settings": ["Custom", "Default", "Custom", "Default", "Default"],
    "OTUs Number Before ReChime": [46552, 48046, 49466, 49126, 49290],
    "OTUs Number After ReChime": [49112, 48011, 49466, 49112, 49275],
    "Change in OTUs": [2560, -35, 0, -14, -15],
    "Percentage Change in OTUs": [5.5, -0.073, 0, -0.029, -0.03],
    "Num of Seqs Before ReChime": [4831720, 4878887, 4567719, 4188648, 4216091],
    "Num of Seqs After ReChime": [4832431, 4881201, 4591864, 4258779, 4288884],
    "Change in Seqs": [711, 2314, 24145, 70131, 72793],
    "Percentage Change in Seqs": [0.015, 0.047, 0.529, 1.673, 1.727]
}

df = pd.DataFrame(data)

# Format specific columns
df['Percentage Change in OTUs'] = df['Percentage Change in OTUs'].map(lambda x: '{:.3g}'.format(x))
df['Percentage Change in Seqs'] = df['Percentage Change in Seqs'].map(lambda x: '{:.3g}'.format(x))

# Create the HTML table with styles
html = df.style.set_table_styles(
    [
        {'selector': 'thead th', 'props': [('background-color', '#4CAF50'), ('color', 'white'), ('text-align', 'center')]},
        {'selector': 'tbody td', 'props': [('text-align', 'center')]},
        {'selector': 'tbody td:nth-child(4), tbody td:nth-child(6), tbody td:nth-child(8), tbody td:nth-child(10)', 'props': [('color', 'red')]}
    ]
).set_properties(**{'border': '1px solid black', 'padding': '5px'}).render()

# Save the HTML to a file
html_file = 'styled_ReChime_OTU_Richness_Table.html'
with open(html_file, 'w') as f:
    f.write(html)

# Define imgkit configuration to use xvfb
config = imgkit.config(xvfb='/usr/bin/xvfb')

# Convert the HTML file to a high-resolution PNG image
imgkit.from_file(html_file, 'styled_ReChime_OTU_Richness_Table_high_res.png', options={'quality': '100', 'width': '1200', 'zoom': '2'}, config=config)

# Alternatively, convert to a JPG image
# imgkit.from_file(html_file, 'styled_ReChime_OTU_Richness_Table_high_res.jpg', options={'quality': '100', 'width': '1200', 'zoom': '2'}, config=config)
