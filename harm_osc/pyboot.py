import pandas as pd
import numpy as np
import sys
from matplotlib import pyplot as plt

sys.path.append(r'/home/ftarantelli/Desktop/projects/NNetwork/harm_osc/statanalysis-master')
from bootstrap import bootstrap_for_primary as btsp

input_file = sys.argv[1]

# Read the dataset into a Pandas DataFrame
#df = pd.read_csv(input_file, names=["Time", "Value"])
time, corr = np.loadtxt(input_file, unpack='True')

time_max = np.max(time)
length = len(time)

num_conf = int(length / (time_max+1))
bin_time = int(time_max / 2 + 1)

new_corr = np.ndarray( shape=(int(time_max / 2 + 1), num_conf) )

for i in range(0, num_conf):
    for k in range(0, bin_time):
        aux = int(i*(time_max+1) + k)
        aux1 = int(i*(time_max+1) + time_max - k)
        new_corr[k,i] = 0.5*corr[aux] + 0.5*corr[aux1]

mean_corr = np.array([0.]*bin_time)
err_corr = np.array([0.]*bin_time)
for ii in range(0, bin_time):
    mean_corr[ii], err_corr[ii] = btsp(lambda x : x, new_corr[ii], 100, 1000)

#print(bin_time)
#time = np.linspace(0, 99, bin_time)
#plt.errorbar(time, mean_corr, yerr = err_corr, fmt = '.')
#plt.show()

if len(sys.argv) > 2:
    output_file = sys.argv[2]
else:
    output_file = "bootstrap_results.dat"


with open(output_file, "w") as file:
    file.write(f"# BootstrapNum \t \t Time \t \t Mean_Corr\n")
    for i, mean in enumerate(mean_corr):
        file.write(f"1 \t \t {i} \t \t {mean:.9f} \n")
    for i, mean in enumerate(mean_corr):
        file.write(f"1 \t \t {i+100} \t \t {mean_corr[99-i]:.9f} \n")

print(f"Bootstrap results have been written to {output_file}")

'''
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

'''
