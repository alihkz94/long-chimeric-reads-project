# FILTER THE TABLE FOR GETTING AN OUTPUT WITH SPECIES

import pandas as pd

# Read the CSV file into a pandas DataFrame without column names
df = pd.read_csv('chimera.csv', header=None)

# Drop rows with any empty cell in the DataFrame
df.dropna(inplace=True)

# Filter rows where all 10 columns have non-null values
filtered_df = df[df.notnull().all(axis=1)]

# Save the filtered DataFrame to a new CSV file
filtered_df.to_csv('filtered_file.csv', index=False, header=False)


