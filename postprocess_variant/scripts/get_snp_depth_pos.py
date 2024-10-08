import numpy as np
import matplotlib.pyplot as plt

# Read the output file generated by the bcftools query command
with open(snakemake.input["depth"], 'r') as f:
#with open('41-Callipepla_depth.txt', 'r') as f: 
    lines = f.readlines()

# Parse the positions, chromosomes, and depths into separate lists
chromosomes = []
positions = []
depths = []
for line in lines:
    fields = line.strip().split('\t')
    chromosomes.append(fields[0])
    positions.append(int(fields[1]))
    depths.append(int(fields[2]))

# Calculate the average depth and standard deviation
avg_depth = np.mean(depths)
std_depth = np.std(depths)

# Exclude positions that exceed 2 standard deviations above and below the mean
threshold_high = avg_depth + (2 * std_depth)
threshold_low = avg_depth - (2 * std_depth)
exceed_positions = [positions[i] for i in range(len(positions)) if depths[i] > threshold_high or depths[i] < threshold_low]
exceed_chromosomes = [chromosomes[i] for i in range(len(chromosomes)) if depths[i] > threshold_high or depths[i] < threshold_low]


# Calculate the percentage of sites that are removed by the filter
num_sites = len(lines)
percent_removed = len(exceed_positions) / num_sites * 100

# Print the percentage of sites that are removed by the filter
print(f'{percent_removed:.2f}% of sites are removed by the filter.')

print(avg_depth)
print(std_depth)
print(threshold_high)
print(threshold_low)

# Write the excluded positions to a new file
with open(snakemake.output["pos"], 'w') as f:
#with open('41-Callipepla_test.txt', 'w') as f:
    for i in range(len(exceed_positions)):
        f.write(f'{exceed_chromosomes[i]}\t{exceed_positions[i]}\n')

# Plot a histogram of the depth distribution with the cutoffs marked as a vertical line
plt.hist(depths, bins=np.logspace(np.log10(1), np.log10(max(depths)), 50))
plt.xscale('log')
plt.axvline(x=threshold_high, color='r', linestyle='dashed', linewidth=2)
plt.axvline(x=threshold_low, color='r', linestyle='dashed', linewidth=2)
plt.xlabel('Depth')
plt.ylabel('Count')
plt.title('Depth Distribution')
plt.savefig('depth_histogram.png', dpi=300, bbox_inches='tight')
