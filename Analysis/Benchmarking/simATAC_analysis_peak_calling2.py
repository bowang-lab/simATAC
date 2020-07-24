import re
from collections import Counter
import numpy as np


# This function compares the coordinates of peak regions from the list of bulk peaks,
# list of peak bins from real bin by cell matrix, and the list of peak bins from simATAC's 
# generated bin by cell matrix, and prints the percentage of regions that have overlap.
# Inputs:
# bulk_peak_file: file address of bulk peak regions.
# real_peak_file: file address of peak bins from real bin by cell matrix.
# sim_peak_file: file address of peak bins from simulated bin by cell matrix.
# Output: -
def peakAnalysis(bulk_peak_file, real_peak_file, sim_peak_file, data):

	bulk_peak = list()
	real = list()
	sim = list()

	# Read the chromosome, starting position, and ending position of peak regions.
	with open (bulk_peak_file, "r") as fileHandler:
	    for line in fileHandler:
	    	part = line.rstrip().split("\t")
	    	if(data == "B"):
	    		bulk_peak.append(part[0:3])
	    	else:
	    		bulk_peak.append(part)

	with open (real_peak_file, "r") as fileHandler:
	    for line in fileHandler:
	    	part = re.split(':|-',line.rstrip()) 
	    	real.append(part)


	with open (sim_peak_file, "r") as fileHandler:
	    for line in fileHandler:
	    	part = re.split(':|-',line.rstrip()) 
	    	sim.append(part)


	# Counts the percentage of common bins between real and simulated peak bins.
	counter = 0
	for r in range(len(real)):
		for s in range(len(sim)):
			if real[r][0] == sim[s][0] and real[r][1] == sim[s][1] and real[r][2] == sim[s][2]:
				counter = counter + 1

	print("Number of common bins in real and simulated peak bins:")
	print(counter)
	print("Percentage of common bins in real and simulated peak bins:")
	print(round(counter/len(real), 2))
	print(round(counter/len(sim), 2))


	# Analyze pair of (bulk peaks, simualted peak bins)
	signr = np.zeros(len(real))
	signs = np.zeros(len(sim))
	signp = np.zeros(len(bulk_peak))
	counter = 0
	for p in range(len(bulk_peak)):
		for s in range(len(sim)):
			if bulk_peak[p][0] == sim[s][0] and (int(bulk_peak[p][1]) >= int(sim[s][1]) or int(bulk_peak[p][2]) <= int(sim[s][2])):
				signp[p] = 1
				signs[s] = 1

	print("Percentage of regions (in bulk peaks) having intersection with simualted peak bins:")
	print(round(np.count_nonzero(signp)/len(bulk_peak), 2))
	print("Percentage of bins (in simulated peak bins) having intersection with bulk peaks:")
	print(round(np.count_nonzero(signs)/len(sim), 2))


	# Analyze pair of (bulk peaks, real peak bins)
	signr = np.zeros(len(real))
	signs = np.zeros(len(sim))
	signp = np.zeros(len(bulk_peak))
	counter = 0
	for p in range(len(bulk_peak)):
		for r in range(len(real)):
			if bulk_peak[p][0] == real[r][0] and (int(bulk_peak[p][1]) >= int(real[r][1]) or int(bulk_peak[p][2]) <= int(real[r][2])):
				signp[p] = 1
				signr[r] = 1

	print("Percentage of regions (in bulk peaks) having intersection with real peak bins:")
	print(round(np.count_nonzero(signp)/len(bulk_peak), 2))
	print("Percentage of bins (in real peak bins) having intersection with bulk peaks:")
	print(round(np.count_nonzero(signr)/len(real), 2))


# Perform analysis on Buenrostro2018
print("Buenrostro2018")
peakAnalysis("../Benchmarking_Data/Buenrostro2018_bulk_peak.bed",
	"../Benchmarking_Data/Buenrostro2018_peaks_real.txt",
	"../Benchmarking_Data/Buenrostro2018_peaks_simulated.txt", "B")

# Perform analysis on Cusanovich2018
print("Cusanovich2018")
peakAnalysis("../Benchmarking_Data/Cusanovich2018_bulk_peak.bed",
	"../Benchmarking_Data/Cusanovich2018_peaks_real.txt",
	"../Benchmarking_Data/Cusanovich2018_peaks_simulated.txt", "C")

# Perform analysis on PBMCs
print("PBMCs")
peakAnalysis("../Benchmarking_Data/BMCs_bulk_peak.bed",
	"../Benchmarking_Data/PBMCs_peaks_real.txt",
	"../Benchmarking_Data/PBMCs_peaks_simulated.txt", "P")
