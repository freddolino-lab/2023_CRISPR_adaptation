#This script is used to handle technical artifacts like PCR jackpots from native adaptation samples. 

#usage python3 phageAD-nativeAD-outlier-IQR.py strain1_spacer_count.txt

#Input file is a text file with counts for each new spacer acquired for a given sample. One number per line. 

#Output file:  First line: the input filename; Second line: sum of spacer counts excluding outlier; Third line: mean of spacer counts excluding outlier; Fourth line: number of excluded values; Fifth line: product of excluded count and mean; Sixth line: cleaned spacer counts.

import numpy as np
import argparse
import os

# Set up argument parsing to get the filename from the command line
parser = argparse.ArgumentParser(description="Detect outliers using the IQR method.")
parser.add_argument('filename', type=str, help='The path to the input txt file containing the numbers')
args = parser.parse_args()

# Read numbers from the specified txt file
filename = args.filename

# Load the numbers from the file, assuming one number per line
data = np.loadtxt(filename)

# Calculate Q1 (25th percentile) and Q3 (75th percentile)
Q1 = np.percentile(data, 25)
Q3 = np.percentile(data, 75)

# Calculate the Interquartile Range (IQR)
IQR = Q3 - Q1

# Determine the lower and upper bounds for outliers
lower_bound = Q1 - 1.5 * IQR
upper_bound = Q3 + 1.5 * IQR

# Identify outliers
outliers = [x for x in data if x < lower_bound or x > upper_bound]

# Count the number of outliers excluded
number_of_excluded = len(outliers)

# Filter the data to remove outliers
cleaned_data = [x for x in data if lower_bound <= x <= upper_bound]

# Calculate the mean and sum of the cleaned data
if cleaned_data:
    mean_cleaned_data = np.mean(cleaned_data)
    sum_cleaned_data = np.sum(cleaned_data)
else:
    mean_cleaned_data = 0  # Avoid division by zero if no cleaned data
    sum_cleaned_data = 0   # If no cleaned data, sum should be zero

# Calculate the product of the number of excluded values and the mean of the cleaned data
product_excluded_mean = number_of_excluded * mean_cleaned_data
clean_total = product_excluded_mean + sum_cleaned_data

# Generate the output file name
output_filename = f"{os.path.splitext(filename)[0]}_output.txt"

# Write the results to the output file
with open(output_filename, 'w') as output_file:
    output_file.write(f"{filename}\n")  # First line: the input filename
    output_file.write(f"{sum_cleaned_data}\n")  # Second line: sum of spacer counts excluding outlier
    output_file.write(f"{mean_cleaned_data}\n")  # Third line: mean of spacer counts excluding outlier
    output_file.write(f"{number_of_excluded}\n")  # Fourth line: number of excluded values
    output_file.write(f"{product_excluded_mean}\n")  # Fifth line: product of excluded count and mean
    output_file.write(f"{clean_total}\n") # Sixth line: cleaned spacer counts

# Print the file path for confirmation
print(f"Results saved to: {output_filename}")
