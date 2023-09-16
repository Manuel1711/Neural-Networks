import pandas as pd
import numpy as np
import sys

input_file = sys.argv[1]

# Read the dataset into a Pandas DataFrame
df = pd.read_csv(input_file, names=["Time", "Value"])

# Group values by time
grouped_data = df.groupby("Time")["Value"]
print(grouped_data)
# Define the number of bootstraps
num_bootstraps = 100

# Perform bootstrapping
bootstrap_means = []

for _ in range(num_bootstraps):
    bootstrap_sample = []

    # Resample values while keeping the same time
    for time, values in grouped_data:
        resampled_values = np.random.choice(values, len(values), replace=True)
        bootstrap_sample.extend(resampled_values)

    # Calculate the mean of the bootstrap sample
    bootstrap_mean = np.mean(bootstrap_sample)
    bootstrap_means.append(bootstrap_mean)


if len(sys.argv) > 2:
    output_file = sys.argv[2]
else:
    output_file = "bootstrap_results.dat"


with open(output_file, "w") as file:
    file.write(f"Time \t \t Mean_Bootstrap\n")
    for i, mean in enumerate(bootstrap_means):
        file.write(f"{i+1} \t \t {mean:.9f}\n")

print(f"Bootstrap results have been written to {output_file}")

