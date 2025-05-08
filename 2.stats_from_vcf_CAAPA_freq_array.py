import calc_stats_for_swifr_revised_Khomanionly_mod
import os
import sys

folder = sys.argv[1]
file = sys.argv[2]
start = sys.argv[3]
end = sys.argv[4]
samplesize = sys.argv[5]
freq1 = sys.argv[6]
out = sys.argv[7]
array_num = sys.argv[8]
last = sys.argv[9]


#if it's the first or the last array, the function scan_targets() will remove 25kb from start and from end (we loose part of the variants, but we ensure the statistics are all calculated in the same window size), but we don't want that for the rest of array (create overlapping windows)
N = 120000 
if int(array_num) > 1 and int(array_num) < int(last): 
	start_ov = int(start) - N//2
	end_ov = int(end) + N//2
elif int(array_num) == 1:
	start_ov = int(start)
	end_ov = int(end) + N//2
elif int(array_num) == int(last):
	start_ov = int(start) - N//2
	end_ov = int(end)


for f in os.listdir(folder):
    file_name = os.path.join(folder, f)
    if os.path.isfile(file_name) and (file in file_name):
        print(file_name)
        vcf = calc_stats_for_swifr_revised.VCF(file_name, n=samplesize)
        stat_matrix, snps, pos = vcf.create_stat(N=int(N), target_freq=(float(freq1),0.99), target_list=None, target_range=(start_ov,end_ov), scale=False, pca=False, output_name=out)


