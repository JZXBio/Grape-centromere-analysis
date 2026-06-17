#!/usr/bin/python
# -*- coding: UTF-8 -*-
import sys
argvs=sys.argv
argv_len=len(sys.argv)
if argv_len < 2 or argvs[1] in ("h", "-h", "help", "-help") or not argvs[1].startswith("step") :
    print ("samples_satellite.py-----help:")
    print ("")
    print ("Usage：")

    print ("step18          Segment the genome, calculate kmers, and analyze homogenization")          
    
    print ("-thread \t\tNumber of threads (default 50). Some steps use multiprocessing.")
    print ("-i      \t\tRequired for step0. Input fasta file.")
    print (" ")
    sys.exit()
argv1=argvs[1]

import subprocess
import csv
import os
import math
import time
import timeit 
import copy
import re #处理正则
from multiprocessing import Pool, cpu_count

###分离参数
args_dict = {}
i = 1
while i < argv_len:
    if argvs[i].startswith('-'):
        arg_name = argvs[i][1:]
        if i + 1 < argv_len and not argvs[i + 1].startswith('-'):
            args_dict[arg_name] = argvs[i + 1]
            i += 1
        else:
            args_dict[arg_name] = True
    i += 1
#print(args_dict)
##
time_start=timeit.default_timer()
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  

if "thread" in args_dict:  thread=int(args_dict["thread"])  
else:thread=70



if  os.path.exists('./samples_satellite')==False:
    subprocess.run(["mkdir samples_satellite"], shell=True)

        
  
### Exploring homogenization
if argv1=="stepall" or "step18" in argv1:  
    kmer_list=[5,7,9,11,14,17, 21, 31, 41, 51, 61, 71, 81, 91, 101] 
    #kmer_list=[31] 
    #kmer_list=[10,11,12,13,14] 
    #kmer_list=[11] 
    #kmer_list=[11] 
    #window=1000;jianju=10000   # Lowest std is at kmer=9
    #block_min=50000;window=5000;jianju=2000    # Lowest std is at 12  
    #window=10000;jianju=1000  # Lowest std is at kmer=11       ## This spacing only affects plotting; with smoothing curve, taking 1000 points per display doesn't cause lag
    #block_min=100000;window=10000;jianju=2000    # Lowest std is at 13
    block_min=20000;window=10000;jianju=5000    # Lowest std is at 11
    #block_min=20000;window=10000;jianju=2000    # Lowest std is at 11
    if  os.path.exists('./samples_satellite/18_Homogenization')==False:
        subprocess.run(["mkdir ./samples_satellite/18_Homogenization"], shell=True)
    if argv1=="stepall" or argv1=="step18_readme":
        print('Printing instructions')
        with open('./samples_satellite/18_Homogenization/readme','w') as f:
            txt=r'''
            18 — Only identify homogenization exceeding 10000bp
            18.0 — Extract sequences
            18.1 — Segment sequences
            18.2 — Calculate various HI indices for different kmers
            18.3 — Overall analysis, determine an appropriate kmer
            18.4 — Overall chromosome differences
            18.4b — Differences per chromosome

                
            18.10 — All chromosomes, and any chromosome with any sample
            18.11 — Analyze the relationship between HOR and homogenization, count HORs within each region
            
            '''
            f.write(txt)     
    if argv1=="stepall" or argv1=="step18" or argv1=="step18.0":  
        print('step18.0 — Extract sequences')
        if  os.path.exists('./samples_satellite/18_Homogenization/0_seqs')==False:
            subprocess.run(["mkdir ./samples_satellite/18_Homogenization/0_seqs"], shell=True)            
        # Generate reverse complement sequence
        def reverse_complement(sequence):
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
            reverse_complement_seq = ''.join(complement[base] for base in reversed(sequence))
            return reverse_complement_seq    
        #    
        region_name_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if eachline_arr[0]=='':continue
                one_region_name=f"{eachline_arr[0]}:{eachline_arr[1]}{eachline_arr[6]}:{eachline_arr[3]}-{eachline_arr[4]}"
                #if one_region_name!="PN40024:region_3-:13509958-14079218" and one_region_name!="PN40024:region_1+:15445074-18499067" :continue
                region_name_list.append(one_region_name)
        print(len(region_name_list))  
        
        def run_step(one_region):
            sample=one_region.split(':')[0]
            chromosome=one_region.split(':')[1][:-1].replace('region_','Chr')
            strand=one_region.split(':')[1][-1]
            pos=one_region.split(':')[-1]
            pos_arr=pos.split('-')
            region_len=int(pos_arr[1])-int(pos_arr[0])+1
            #print(sample,chromosome,strand,pos)
            subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample}/{sample}.fasta {chromosome}:{pos} > ./samples_satellite/18_Homogenization/0_seqs/{one_region}_tmp"], shell=True)
            with open(f"./samples_satellite/18_Homogenization/0_seqs/{one_region}_tmp",'r') as f:
                next(f)
                seq=f.read().replace('\n','')
            if strand=='-':
                seq=reverse_complement(seq)
            with open(f"./samples_satellite/18_Homogenization/0_seqs/{one_region}",'w') as f:  
                f.write(seq)
            seq_len=len(seq)
            if region_len!=seq_len:print(f'error\t{region_len}\t{seq_len}')
            
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, region_name_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(region_name_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()       
    if argv1=="stepall" or argv1=="step18"  or argv1=="step18.1":  
        print('step18.1 — Segment sequences')
        if  os.path.exists('./samples_satellite/18_Homogenization/1_seqs')==False:
            subprocess.run(["mkdir ./samples_satellite/18_Homogenization/1_seqs"], shell=True)            
        region_name_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if eachline_arr[0]=='':continue
                one_region_name=f"{eachline_arr[0]}:{eachline_arr[1]}{eachline_arr[6]}:{eachline_arr[3]}-{eachline_arr[4]}"
                #if one_region_name!="PN40024:region_3-:13509958-14079218" and one_region_name!="PN40024:region_1+:15445074-18499067" :continue
                region_name_list.append(one_region_name)
        print(len(region_name_list))  
        
        # Create summary file
        summary_file = "./samples_satellite/18_Homogenization/1_infosum"
        os.makedirs(os.path.dirname(summary_file), exist_ok=True)  # Ensure directory exists
        with open(summary_file, 'w') as f:
            f.write("Original_Region\tregion_length\tStart\tEnd\tLength\tSegment_File\n")  # Write header
        # Create filtered summary file
        summary_file_filter = "./samples_satellite/18_Homogenization/1_infosum_filter"
        os.makedirs(os.path.dirname(summary_file), exist_ok=True)  # Ensure directory exists
        with open(summary_file_filter, 'w') as f:
            f.write("Original_Region\tregion_length\tStart\tEnd\tLength\tSegment_File\n")  # Write header
            
        region_name_list = []
        with open('./samples_satellite/2_good_regions', 'r') as f:
            for line in f:
                eachline_arr = line.strip().split('\t')
                if eachline_arr[0] == 'sample': continue
                if eachline_arr[0] == '': continue
                one_region_name = f"{eachline_arr[0]}:{eachline_arr[1]}{eachline_arr[6]}:{eachline_arr[3]}-{eachline_arr[4]}"
                region_name_list.append(one_region_name)
        print(len(region_name_list))  
         
        def run_step(one_region):
            region_strand = one_region.split(':')[-2][-1]
            region_pos = one_region.split(':')[-1].split('-')
            pos_start = int(region_pos[0])
            pos_end = int(region_pos[1])
            region_length = pos_end - pos_start + 1
            
            with open(f"./samples_satellite/18_Homogenization/0_seqs/{one_region}", 'r') as f:
                seq = f.read()
            
            segment_info = []  # Store all segment info for the current region
            
            if region_strand == '+':
                tmp_delta = 0
                while tmp_delta < region_length:
                    part_start_index = tmp_delta
                    part_end_index = tmp_delta + window
                    if part_end_index > region_length:
                        part_end_index = region_length
                    part_seq = seq[part_start_index:part_end_index]
                    part_start_pos = pos_start + part_start_index
                    part_end_pos = pos_start + part_end_index - 1
                    segment_file = f"./samples_satellite/18_Homogenization/1_seqs/{one_region}___{part_start_pos}-{part_end_pos}"
                    with open(segment_file, 'w') as f:
                        f.write(f'>{one_region}___{part_start_pos}-{part_end_pos}\n{part_seq}')
                    # Record segment info
                    segment_info.append([one_region,region_length, part_start_pos, part_end_pos, len(part_seq), segment_file])
                    tmp_delta += window
                
                # Handle the last segment (if not fully covered)
                if (region_length % window) != 0:
                    part_start_index = max(0, region_length - window)
                    part_end_index = region_length
                    part_seq = seq[part_start_index:part_end_index]
                    part_start_pos = pos_start + part_start_index
                    part_end_pos = pos_start + part_end_index - 1
                    segment_file = f"./samples_satellite/18_Homogenization/1_seqs/{one_region}___{part_start_pos}-{part_end_pos}"
                    with open(segment_file, 'w') as f:
                        f.write(f'>{one_region}___{part_start_pos}-{part_end_pos}\n{part_seq}')
                    segment_info.append([one_region,region_length, part_start_pos, part_end_pos, len(part_seq), segment_file])
            
            elif region_strand == '-':
                tmp_delta = region_length
                while tmp_delta > 0:
                    part_start_index = max(0, tmp_delta - window)
                    part_end_index = tmp_delta
                    part_seq = seq[part_start_index:part_end_index]
                    part_start_pos = pos_end - part_start_index
                    part_end_pos = pos_end - part_end_index + 1
                    segment_file = f"./samples_satellite/18_Homogenization/1_seqs/{one_region}___{part_start_pos}-{part_end_pos}"
                    with open(segment_file, 'w') as f:
                        f.write(f'>{one_region}___{part_start_pos}-{part_end_pos}\n{part_seq}')
                    segment_info.append([one_region,region_length, part_start_pos, part_end_pos, len(part_seq), segment_file])
                    tmp_delta -= window
                
                # Handle the last segment (if not fully covered)
                if (region_length % window) != 0:
                    part_start_index = 0
                    part_end_index = min(window, region_length)
                    part_seq = seq[part_start_index:part_end_index]
                    part_start_pos = pos_end - part_start_index
                    part_end_pos = pos_end - part_end_index + 1
                    segment_file = f"./samples_satellite/18_Homogenization/1_seqs/{one_region}___{part_start_pos}-{part_end_pos}"
                    with open(segment_file, 'w') as f:
                        f.write(f'>{one_region}___{part_start_pos}-{part_end_pos}\n{part_seq}')
                    segment_info.append([one_region,region_length, part_start_pos, part_end_pos, len(part_seq), segment_file])
                    
            # Convert to tuples for deduplication
            unique_segments = [list(seg) for seg in set(tuple(sublist) for sublist in segment_info)]
             
            # Sort by a field (e.g., by part_start_pos)
            unique_segments.sort(key=lambda x: x[2])
            unique_segments.sort(key=lambda x: x[0])
            # Write the current region's segment info to the summary file
            with open(summary_file_filter, 'a') as f2:
                with open(summary_file, 'a') as f:
                    for info in unique_segments:
                        f.write("\t".join(map(str, info)) + "\n")
                        
                        if int(info[1])>=block_min and int(info[4])==window:
                            f2.write("\t".join(map(str, info)) + "\n")
            
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, region_name_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(region_name_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()               
            
    if argv1=="stepall" or argv1=="step18"  or argv1=="step18.2":
        print('Take kmers, calculate HI indices, including main kmer percentage (0-100), kmer diversity (0-100), HI index (standardized hi entropy, 0-100)')
        if  os.path.exists('./samples_satellite/18_Homogenization/2_kmer')==False:
            subprocess.run(["mkdir ./samples_satellite/18_Homogenization/2_kmer"], shell=True)      

                
        import pandas as pd
        import numpy as np
        from collections import Counter
         
        def calculate_hi_indices(seq, k):
            """
            Calculate centromere homogenization indices (HI_dominant, HI_diversity, HI_entropy)
            
            Parameters:
                seq (str): DNA sequence string
                k (int): k-mer length, default is 21
            
            Returns:
                dict: Dictionary containing three indices
            """
            # 1. Extract all k-mers and count frequencies
            kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
            kmer_counts = Counter(kmers)
            
            # Convert to DataFrame for analysis
            kmers_df = pd.DataFrame(list(kmer_counts.items()), columns=["kmer", "count"])
            
            # 2. Calculate total k-mer count and unique k-mer count
            total_count = kmers_df["count"].sum()
            unique_count = len(kmers_df)
            
            if total_count == 0:
                return {"HI_dominant": 0, "HI_diversity": 0, "HI_entropy": 0}
            
            # 3. Calculate HI_dominant (dominant k-mer percentage)
            dominant_count = kmers_df["count"].max()
            hi_dominant = round((dominant_count / total_count),4)
            
            # 4. Calculate HI_diversity (unique k-mer percentage)
            hi_diversity = round((unique_count / total_count) ,4)
            
            # 5. Calculate HI_entropy (normalized entropy)
            p = kmers_df["count"] / total_count
            entropy = -np.sum(p * np.log2(p))
            hi_entropy = 1 - (entropy / np.log2(unique_count)) if unique_count > 1 else 1.0
            
        
            # 6. Calculate Coverage (corrected formula: Coverage = Total_kmers / Unique_kmers)
            Coverage = round(total_count / unique_count, 2)
            
            # 7. Calculate rare k-mer ratio (k-mers with count ≤ 2)
            rare_kmer_ratio = round((kmers_df[kmers_df["count"] <= 2].shape[0] / unique_count), 4)
            
            return [total_count,dominant_count,unique_count,hi_dominant,hi_diversity,round(entropy,3),round(hi_entropy,5),Coverage,rare_kmer_ratio]
  
                    
        summary_file_filter = "./samples_satellite/18_Homogenization/1_infosum_filter" 
        input_list=[]
        with open(summary_file_filter,'r') as f:
            next(f)
            for line in f:
                Original_Region,region_length,Start,End,Length,Segment_File=line.strip().split('\t')
                input_list.append(Segment_File)
        input_list_num=len(input_list)        
        print(f'Input count: {input_list_num}')
        
            
        #####
        def run_step(info)    :    
            one_segment_file,kmer=info
            region=one_segment_file.split('/')[-1].split('___')[0]
            chromosome=one_segment_file.split(':')[-2][:-1].replace('region_','Chr')
            with open(one_segment_file,'r')  as f:
                next(f)
                seq=f.read()
                result=calculate_hi_indices(seq,kmer)
                #print(result)
                result_str='\t'.join([str(x) for x in result])
                #print(result)
                with open(f'./samples_satellite/18_Homogenization/2_kmer/kmer_{kmer}','a') as f:
                    f.write(f"{chromosome}\t{region}\t{one_segment_file}\t{result_str}\n")
                    
        kmer_list_len=len(kmer_list)
        i=0
        for kmer in kmer_list:
            i+=1
            print(f'\n\nkmer processing progress: {i}/{kmer_list_len}\t\tkmer{kmer}')
            if os.path.exists(f'./samples_satellite/18_Homogenization/2_kmer/kmer_{kmer}')==True:
                continue
            with open(f'./samples_satellite/18_Homogenization/2_kmer/kmer_{kmer}','w') as f:
                f.write(f"chromosome\tregion\tsegment_file\ttotal_count\tdominant_count\tunique_count\thi_dominant\thi_diversity\tentropy\thi_entropy\tCoverage\trare_kmer_ratio\n")           
            input_list2=[]
            for one in input_list:
                input_list2.append([one,kmer])
            with Pool(processes=thread) as pool:
                # Use imap to get results one by one
                for i, result in enumerate(pool.imap(run_step, input_list2), start=1):
                    # Results can be processed here, e.g., stored or printed
                    progress = (i / len(input_list2)) * 100
                    sys.stdout.write(f"\rProgress: {progress:.2f}%")
                    sys.stdout.flush()                      
    if argv1=="stepall" or argv1=="step18"  or argv1=="step18.3":   
        import numpy as np        
        if  os.path.exists('./samples_satellite/18_Homogenization/3_stat')==False:
            subprocess.run(["mkdir ./samples_satellite/18_Homogenization/3_stat"], shell=True) 
        print('Determining the appropriate kmer value')
        kmer_list_len=len(kmer_list)
        dict_kmer_info={}
        for kmer in kmer_list:
            print(kmer)
            if kmer not in dict_kmer_info:
                dict_kmer_info[kmer]={}
                dict_kmer_info[kmer]['hi_entropy']=[]
                dict_kmer_info[kmer]['Coverage']=[]
                dict_kmer_info[kmer]['rare_kmer_ratio']=[]
            with open(f"./samples_satellite/18_Homogenization/2_kmer/kmer_{kmer}",'r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    hi_entropy=         float(eachline_arr[9])
                    Coverage=           float(eachline_arr[10])
                    rare_kmer_ratio=    float(eachline_arr[11])
                    dict_kmer_info[kmer]['hi_entropy'].append(hi_entropy)
                    dict_kmer_info[kmer]['Coverage'].append(Coverage)
                    dict_kmer_info[kmer]['rare_kmer_ratio'].append(rare_kmer_ratio)
                  
        min_length = min(len(info['Coverage']) for info in dict_kmer_info.values())
        print("Minimum list length:", min_length)
        #print(len(dict_kmer_info[9]['rare_kmer_ratio'])            );sys.exit()  
        
        
        with open('./samples_satellite/18_Homogenization/3_stat/hi_entropy_sd','w') as f:            
            f.write('kmer\tkmer_num\thi_entropy_sd\n')
        

        maxstd_kmer='';
        for kmer in kmer_list:
            info=dict_kmer_info[kmer]
            hi_entropy_list=info['hi_entropy']
            Coverage_list=info['Coverage']
            rare_kmer_ratio_list=info['rare_kmer_ratio']
            dict_kmer_info[kmer]['hi_entropy'].sort()
            dict_kmer_info[kmer]['Coverage'].sort()
            dict_kmer_info[kmer]['rare_kmer_ratio'].sort() 
            hi_entropy_sd = round(np.std(hi_entropy_list),6)
            ### Record the kmer with the highest std
            if maxstd_kmer=='':
                maxstd_kmer=kmer
                maxstd_kmer_value=hi_entropy_sd
            if hi_entropy_sd>maxstd_kmer_value:
                maxstd_kmer=kmer
                maxstd_kmer_value=hi_entropy_sd                
            ###
            with open('./samples_satellite/18_Homogenization/3_stat/hi_entropy_sd','a') as f:
                f.write(f'kmer_{kmer}\t{kmer}\t{hi_entropy_sd}\n')
            '''with open(f'./samples_satellite/18_Homogenization/3_stat/Coverage_list___kmer{kmer}','w') as f:    
                for one in Coverage_list:    
                    f.write(one+'\n')                
            with open(f'./samples_satellite/18_Homogenization/3_stat/rare_kmer_ratio_list___kmer{kmer}','w') as f:    
                for one in rare_kmer_ratio_list:    
                    f.write(one+'\n')                    '''
        
        ## Determine the top and bottom 15% boundary points
        hi_entropy_list=[]
        with open(f"./samples_satellite/18_Homogenization/2_kmer/kmer_{maxstd_kmer}",'r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                hi_entropy=         float(eachline_arr[9])
                hi_entropy_list.append(hi_entropy)
        # Calculate top and bottom boundaries (15% quantiles)
        bound_14 = np.percentile(hi_entropy_list, 14)  # Lower boundary (15%)
        bound_43 = np.percentile(hi_entropy_list, 43)  # Lower boundary (15%)
        bound_57 = np.percentile(hi_entropy_list, 57)  # Lower boundary (15%)
        bound_86 = np.percentile(hi_entropy_list, 86)  # Upper boundary (85%)                
        with open('./samples_satellite/18_Homogenization/3_stat/maxstd_kmer','w') as f:
            f.write(f"maxstd_kmer\thi_entropy_sd\tbound_14\tbound_43\tbound_57\tbound_86\n")
            f.write(f"{maxstd_kmer}\t{hi_entropy_sd}\t{bound_14}\t{bound_43}\t{bound_57}\t{bound_86}\n")
            
        with open('./samples_satellite/18_Homogenization/3_stat/infos','w') as f:
            f.write(f"kmer\tkmer_num\tindex\thi_entropy\tCoverage\trare_kmer_ratio\n")
            kk=0;
            while kk<min_length   : 
                for kmer in kmer_list:
                    hi_entropy=dict_kmer_info[kmer]['hi_entropy'][kk]
                    Coverage=dict_kmer_info[kmer]['Coverage'][kk]
                    rare_kmer_ratio=dict_kmer_info[kmer]['rare_kmer_ratio'][kk]
                    if kk%jianju==0:
                        f.write(f"kmer_{kmer}\t{kmer}\t{kk}\t{hi_entropy}\t{Coverage}\t{rare_kmer_ratio}\n")
                ##
                kk+=1
                    
    if argv1=="stepall" or argv1=="step18"  or argv1=="step18.3" or argv1=="step18.3p":
        print('Plotting')            
        R_txt=f'''library(ggplot2)
library(dplyr)
#install.packages("stringr") 
library("stringr")

# Set working directory
setwd('./')

# Read data
input_file1 <- read.table('infos', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_line(data = input_file1, aes(x = index, y = hi_entropy, color = kmer_num, group = kmer_num))
p <- p + theme_classic()
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("hi_entropy.pdf", width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()

# Create plot
p <- ggplot()
p <- p + geom_line(data = input_file1, aes(x = index, y = Coverage, color = kmer_num, group = kmer_num))
p <- p + theme_classic()
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("Coverage.pdf", width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()

# Create plot
p <- ggplot()
p <- p + geom_line(data = input_file1, aes(x = index, y = rare_kmer_ratio, color = kmer_num, group = kmer_num))
p <- p + theme_classic()
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("rare_kmer_ratio.pdf", width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()


# Read data
input_file2 <- read.table('hi_entropy_sd', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')
# Create plot
p <- ggplot()
p <- p + geom_col(data = input_file2, aes(x = kmer_num, y = hi_entropy_sd))
p <- p + theme_classic()
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("hi_entropy_sd.pdf", width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()

p <- ggplot()
p <- p + geom_line(data = input_file2, aes(x = kmer_num, y = hi_entropy_sd))
p <- p + geom_point(data = input_file2, aes(x = kmer_num, y = hi_entropy_sd,color=kmer_num), shape = 16,stroke = 0,size=3)
p <- p + theme_classic()
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("hi_entropy_sd_line.pdf", width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()
'''
        with open(f'./samples_satellite/18_Homogenization/3_stat/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = f"./samples_satellite/18_Homogenization/3_stat/"
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../../')            
        
    if argv1=="stepall" or argv1=="step18"  or argv1=="step18.4":
        print('18.4 — Overall chromosome differences')
        import numpy as np
        if  os.path.exists('./samples_satellite/18_Homogenization/4_chr')==False:
            subprocess.run(["mkdir ./samples_satellite/18_Homogenization/4_chr"], shell=True)          
        #    
        with open('./samples_satellite/18_Homogenization/3_stat/maxstd_kmer','r') as f:   
            next(f)
            for line in f:
                maxstd_kmer,hi_entropy_sd,bound_14,bound_43,bound_57,bound_86=line.strip().split('\t')
                bound_14,bound_43,bound_57,bound_86=float(bound_14),float(bound_43),float(bound_57),float(bound_86)
        print(f"maxstd_kmer：{maxstd_kmer}")   
        #    
        kmer_list_len=len(kmer_list)
        dict_kmer_info={}
        for kmer in kmer_list:
            print(kmer)
            if str(kmer) != maxstd_kmer:continue
            print(f'Only calculate this kmer for kmer{kmer}')
            #########
            if kmer not in dict_kmer_info:
                dict_kmer_info[kmer]={}
            with open(f"./samples_satellite/18_Homogenization/2_kmer/kmer_{kmer}",'r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    chromosome=         eachline_arr[0]
                    hi_entropy=         float(eachline_arr[9])
                    
                    if chromosome not in dict_kmer_info[kmer]:dict_kmer_info[kmer][chromosome]=[]
                    dict_kmer_info[kmer][chromosome].append(hi_entropy)
                  
            with open(f"./samples_satellite/18_Homogenization/4_chr/kmer_{kmer}",'w') as f:    
                f.write(f"chromosome\tregion_number\tmean\tstd_dev\tlow_percent\tmiddle_percent\thigh_percent\n")
                dict_chr_info=dict_kmer_info[kmer]
                chromosome_list= dict_chr_info.keys()
                sorted_chromosomes = sorted(chromosome_list, key=lambda x: int(x[3:]) if x[3:].isdigit() else float('inf'))
                
                for chromosome in sorted_chromosomes:
                    value_list=dict_chr_info[chromosome]
                    #num
                    number=len(value_list)
                    # Calculate mean
                    from scipy import stats

                    # Calculate trimmed mean, e.g., removing 10% from both ends
                    mean = stats.trim_mean(value_list, proportiontocut=0.15)
                    #mean = np.mean(value_list)
                     
                    # Calculate standard deviation (population standard deviation, denominator is n)
                    std_dev = np.std(value_list)
                    ##
                    low_num = len([x for x in value_list if x < bound_14])
                    middle_num= len([x for x in value_list if x>bound_43 and x < bound_57])
                    high_num= len([x for x in value_list if x > bound_86])
                    low_percent=round(low_num/number,3)
                    middle_percent=round(middle_num/number,3)
                    high_percent=round(high_num/number,3)
                    ##
                    f.write(f"{chromosome}\t{number}\t{round(mean,3)}\t{round(std_dev,3)}\t{low_percent}\t{middle_percent}\t{high_percent}\n")
    if argv1=="stepall" or argv1=="step18"  or argv1=="step18.4" or argv1=="step18.4b":
        print('18.4b — Differences per chromosome')
        import numpy as np
        if  os.path.exists('./samples_satellite/18_Homogenization/4_chr')==False:
            subprocess.run(["mkdir ./samples_satellite/18_Homogenization/4_chr"], shell=True)          
        #
        print('Loading a sample information table ./samples_satellite/sample_info')
        dict_sample_type={}
        with open('./samples_satellite/sample_info') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=6:continue
                sample_name=eachline_arr[0]
                sample_type=eachline_arr[3]
                sample_type2=eachline_arr[5]
                if sample_type=="Eurasian" and sample_type2=="Table":           sample_type3="Table"
                elif sample_type=="Eurasian" and sample_type2=="Wine":           sample_type3="Wine"
                elif sample_type=="East_Asia":                                  sample_type3="East_Asia"
                elif sample_type=="America":                                  sample_type3="America"
                else:sample_type=="other"

                dict_sample_type[sample_name]=sample_type3
        #
        with open('./samples_satellite/18_Homogenization/3_stat/maxstd_kmer','r') as f:   
            next(f)
            for line in f:
                maxstd_kmer,hi_entropy_sd,bound_14,bound_43,bound_57,bound_86=line.strip().split('\t')
                bound_14,bound_43,bound_57,bound_86=float(bound_14),float(bound_43),float(bound_57),float(bound_86)
        print(f"maxstd_kmer：{maxstd_kmer}")  

        #    
        kmer_list_len=len(kmer_list)
        dict_kmer_info={}
        for kmer in kmer_list:
            print(kmer)
            if str(kmer) != maxstd_kmer:continue
            print(f'Only calculate this kmer for kmer{kmer}')
            #########
            if kmer not in dict_kmer_info:
                dict_kmer_info[kmer]={}
            with open(f"./samples_satellite/18_Homogenization/2_kmer/kmer_{kmer}",'r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    sample=             eachline_arr[1].split(':')[0]
                    chromosome=         eachline_arr[0]
                    samplechr=          sample+'___'+chromosome
                    hi_entropy=         float(eachline_arr[9])
                    
                    if samplechr not in dict_kmer_info[kmer]:dict_kmer_info[kmer][samplechr]=[]
                    dict_kmer_info[kmer][samplechr].append(hi_entropy)
                  
            with open(f"./samples_satellite/18_Homogenization/4_chr/kmer_{kmer}_b",'w') as f:    
                f.write(f"samplechr\tsample\tchromosome\tregion_number\tmean\tstd_dev\tlow_percent\tmiddle_percent\thigh_percent\tsample_type\n")
                dict_chr_info=dict_kmer_info[kmer]
                samplechr_list= dict_chr_info.keys()
                sorted_samplechr = sorted(samplechr_list, key=lambda x: int(x[3:]) if x[3:].isdigit() else float('inf'))
                
                for samplechr in sorted_samplechr:
                    value_list=dict_chr_info[samplechr]
                    #num
                    number=len(value_list)
                    # Calculate mean
                    #mean = np.mean(value_list)
                    mean = stats.trim_mean(value_list, proportiontocut=0.15) 
                    # Calculate standard deviation (population standard deviation, denominator is n)
                    std_dev = np.std(value_list)
                    ##
                    low_num = len([x for x in value_list if x < bound_14])
                    middle_num= len([x for x in value_list if x>bound_43 and x < bound_57])
                    high_num= len([x for x in value_list if x > bound_86])
                    low_percent=round(low_num/number,3)
                    middle_percent=round(middle_num/number,3)
                    high_percent=round(high_num/number,3)
                    ##
                    sample=         samplechr.split('___')[0]
                    sample_type=        dict_sample_type[sample]
                    if sample_type=='other':continue
                    chromosome=     samplechr.split('___')[1]
                    f.write(f"{samplechr}\t{sample}\t{chromosome}\t{number}\t{round(mean,3)}\t{round(std_dev,3)}\t{low_percent}\t{middle_percent}\t{high_percent}\t{sample_type}\n")
    if argv1=="stepall" or argv1=="step18"  or argv1=="step18.4" or argv1=="step18.4b" or argv1=="step18.4bp":
        with open('./samples_satellite/18_Homogenization/3_stat/maxstd_kmer','r') as f:   
            next(f)
            for line in f:
                maxstd_kmer,hi_entropy_sd,bound_14,bound_43,bound_57,bound_86=line.strip().split('\t')
                bound_14,bound_43,bound_57,bound_86=float(bound_14),float(bound_43),float(bound_57),float(bound_86)
        print(f"maxstd_kmer：{maxstd_kmer}")         
        print('Plotting')            
        R_txt=f'''library(ggplot2)
library(dplyr)
#install.packages("stringr") 
library("stringr")

# Set working directory
setwd('./')

# Read data
input_file1 <- read.table('kmer_{maxstd_kmer}_b', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')
###
input_file2 <- input_file1 %>%arrange(as.numeric(sub("Chr", "", chromosome)))
chr_order <- unique(input_file2$chromosome)
input_file2$chromosome <- factor(input_file2$chromosome, levels = chr_order)


# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file2, aes(x = chromosome, y = mean,color=sample_type,size=region_number*0.01),position = position_jitter(width = 0.3),alpha=0.95,  shape = 16,stroke = 0)
p <- p + geom_boxplot(data = input_file2, aes(x = chromosome, y = mean),size=1,outlier.shape = NA)

p <- p +scale_color_manual(name = "Variant Mode", 
                           values = c("Table" = "#E5E54C",   #8fc31f
                                      "Wine" = "#7e318e",  
                                      "East_Asia" = "#0066ff",  
                                      "America" = "#29a329"))  #33cc33
p <- p + theme_classic()
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("chromosome.pdf", width = 40 / 2.54, height = 15 / 2.54)
print(p)
dev.off()


p <- ggplot()
p <- p + geom_point(data = input_file2, aes(x = region_number, y = mean,color=sample_type,size=region_number*0.01),position = position_jitter(width = 0.3),alpha=0.95,  shape = 16,stroke = 0)
p <- p +scale_color_manual(name = "Variant Mode", 
                           values = c("Table" = "#E5E54C",   #8fc31f
                                      "Wine" = "#7e318e",  
                                      "East_Asia" = "#0066ff",  
                                      "America" = "#29a329"))  #33cc33
p <- p + theme_classic()
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("chromosome_size2mean.pdf", width = 20 / 2.54, height = 16 / 2.54)
print(p)
dev.off()
'''
        with open(f'./samples_satellite/18_Homogenization/4_chr/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = f"./samples_satellite/18_Homogenization/4_chr/"
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../../')                            
    if argv1=="stepall" or argv1=="step18"  or argv1=="step18.4" or argv1=="step18.4c":     
        print('Calculating differences between different groups')
        import pingouin as pg
        import numpy as np        
        with open('./samples_satellite/18_Homogenization/3_stat/maxstd_kmer','r') as f:   
            next(f)
            for line in f:
                maxstd_kmer,hi_entropy_sd,bound_14,bound_43,bound_57,bound_86=line.strip().split('\t')
                bound_14,bound_43,bound_57,bound_86=float(bound_14),float(bound_43),float(bound_57),float(bound_86)
        print(f"maxstd_kmer：{maxstd_kmer}") 
        
        dict_chr_type_info={}
        with open(f"./samples_satellite/18_Homogenization/4_chr/kmer_{maxstd_kmer}_b",'r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                #sample=         eachline_arr[1].replace('_hap','.hap')
                chromosome=     eachline_arr[2]
                sample_type=    eachline_arr[9]
                mean=           eachline_arr[4]
                #if sample[-4:] not in ['.hap1','.hap2']:print('error');sys.exit()
                #sample_pure=sample[-4:]
                #if sample_pure not in dict_samplepure_chr_type_info: dict_samplepure_chr_type_info[sample_pure]={}
                if chromosome not in dict_chr_type_info:                    dict_chr_type_info[chromosome]={}
                if sample_type not in dict_chr_type_info[chromosome]:       dict_chr_type_info[chromosome][sample_type]=[]
                dict_chr_type_info[chromosome][sample_type].append(mean)
        #############        
        with open(f"./samples_satellite/18_Homogenization/4_chr/kmer_{maxstd_kmer}_c2",'w') as f2:     
            with open(f"./samples_satellite/18_Homogenization/4_chr/kmer_{maxstd_kmer}_c",'w') as f:     
                f2.write('chromosome\tgroup\tU\tP_value\tRBC\tCliffdelta\n')
                f.write('chromosome\tTable\tWine\tEast_Asia\tAmerica\n')
                for chromosome,dict_type_info  in dict_chr_type_info.items():
                    if chromosome=='Chr18':continue
                    if "Table" not in dict_type_info:       Table_str=''; Table_group=[]
                    else:                                   Table_group=dict_type_info["Table"];                Table_str=",".join(Table_group)
                    if "Wine" not in dict_type_info:        Wine_str=''; Wine_group=[]
                    else:                                   Wine_group=dict_type_info["Wine"];                  Wine_str=",".join(Wine_group)
                    if "East_Asia" not in dict_type_info:   East_Asia_str=''; East_Asia_group=[]
                    else:                                   East_Asia_group=dict_type_info["East_Asia"];        East_Asia_str=",".join(East_Asia_group)
                    if "America" not in dict_type_info:     America_str=''; America_group=[]
                    else:                                   America_group=dict_type_info["America"];            America_str=",".join(America_group)      
                    ##
                    Table_group=[float(x) for x in Table_group]
                    Wine_group=[float(x) for x in Wine_group]
                    East_Asia_group=[float(x) for x in East_Asia_group]
                    America_group=[float(x) for x in America_group]
                    #
                    cultivar_group=     Table_group+Wine_group
                    wild_group=         East_Asia_group+America_group
    
                    # Define function: calculate Mann-Whitney U test, RBC, and Cliff's delta
                    def compute_mwu_effects(group1, group2):
                        # Calculate Mann-Whitney U test
                        result = pg.mwu(group1, group2, alternative='two-sided')
                        U = result['U-val'].values[0]
                        p_val = result['p-val'].values[0]
                        rbc = result['RBC'].values[0]
                        
                        # Calculate Cliff's delta
                        n1, n2 = len(group1), len(group2)
                        cliff_delta = (2 * U / (n1 * n2)) - 1  # Normalized to [-1, 1]
                        
                        return U, p_val, rbc, cliff_delta
                    
                    # Calculate effect sizes for each comparison group
                    Table2Wine_U, Table2Wine_p, Table2Wine_rbc, Table2Wine_cliffdelta = compute_mwu_effects(Table_group,Wine_group)
                    EA2Am_U, EA2Am_p, EA2Am_rbc, EA2Am_cliffdelta = compute_mwu_effects(East_Asia_group,America_group)
                    Cult2Wild_U, Cult2Wild_p, Cult2Wild_rbc, Cult2Wild_cliffdelta = compute_mwu_effects(cultivar_group,wild_group)
                    
                    # Write to file f (summary data)
                    f.write(f"{chromosome}\t"
                            f"{Table_str}\t{Wine_str}\t{East_Asia_str}\t{America_str}\n")
                    
                    # Write to file f2 (sub-item data)
                    f2.write(f"{chromosome}\tTable2Wine\t{Table2Wine_U}\t{Table2Wine_p}\t{Table2Wine_rbc}\t{Table2Wine_cliffdelta}\n")
                    f2.write(f"{chromosome}\tEast_Asia2America\t{EA2Am_U}\t{EA2Am_p}\t{EA2Am_rbc}\t{EA2Am_cliffdelta}\n")
                    f2.write(f"{chromosome}\tCult2Wildivar\t{Cult2Wild_U}\t{Cult2Wild_p}\t{Cult2Wild_rbc}\t{Cult2Wild_cliffdelta}\n")
    if argv1=="stepall" or argv1=="step18"  or argv1=="step18.4" or argv1=="step18.4c" or argv1=="step18.4cp":
        with open('./samples_satellite/18_Homogenization/3_stat/maxstd_kmer','r') as f:   
            next(f)
            for line in f:
                maxstd_kmer,hi_entropy_sd,bound_14,bound_43,bound_57,bound_86=line.strip().split('\t')
                bound_14,bound_43,bound_57,bound_86=float(bound_14),float(bound_43),float(bound_57),float(bound_86)
        print(f"maxstd_kmer：{maxstd_kmer}")         
        print('Plotting')            
        R_txt=f'''library(ggplot2)
library(dplyr)
#install.packages("stringr") 
library("stringr")

# Set working directory
setwd('./')

# Read data
input_file1 <- read.table('kmer_{maxstd_kmer}_c2', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')
###
input_file2 <- input_file1 %>%arrange(as.numeric(sub("Chr", "", chromosome)))
chr_order <- unique(input_file2$chromosome)
input_file2$chromosome <- factor(input_file2$chromosome, levels = chr_order)

input_file2$group <- factor(
  input_file2$group,
  levels = c("Table2Wine", "East_Asia2America", "Cult2Wildivar")
) 

# Create plot
p <- ggplot()
p <- p + geom_col(data = input_file2,aes(x = chromosome, y = RBC, fill = group),position = position_dodge(width = 0.7),  width = 0.6 ) 
p <- p + theme_classic()+scale_fill_manual(
                           values = c("Table2Wine"="#e1e15d","East_Asia2America"="#107693","Cult2Wildivar"="#b25062")) 
    
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("chromosome_RBC.pdf", width = 20 / 2.54, height = 5 / 2.54)
print(p)
dev.off()

# Create plot
p <- ggplot()
p <- p + geom_col(data = input_file2,aes(x = chromosome, y = Cliffdelta, fill = group),position = position_dodge(width = 0.7),  width = 0.6 ) 
p <- p + theme_classic()
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("chromosome_cliffdelta.pdf", width = 20 / 2.54, height = 5 / 2.54)
print(p)
dev.off()

# Create plot
p <- ggplot()
p <- p + geom_col(data = input_file2,aes(x = chromosome, y = U, fill = group),position = position_dodge(width = 0.7),  width = 0.6 ) 
p <- p + theme_classic()
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("chromosome_U.pdf", width = 20 / 2.54, height = 5 / 2.54)
print(p)
dev.off()

# Create plot
p <- ggplot()
p <- p + geom_col(data = input_file2,aes(x = chromosome, y = -log10(P_value), fill = group),position = position_dodge(width = 0.7),  width = 0.6 ) 
p <- p + theme_classic()
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("chromosome_P_value.pdf", width = 20 / 2.54, height = 5 / 2.54)
print(p)
dev.off()
'''
        with open(f'./samples_satellite/18_Homogenization/4_chr/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = f"./samples_satellite/18_Homogenization/4_chr/"
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../../')                                                                              
                    
    ### The plus versions are for convenient editing 
    if argv1=="stepall" or argv1=="step18"  or argv1=="step18.10": 
        print('Visualizing all chromosomes')
        if  os.path.exists('./samples_satellite/18_Homogenization/10_plot')==False:
            subprocess.run(["mkdir ./samples_satellite/18_Homogenization/10_plot"], shell=True)
        with open(f'./samples_satellite/18_Homogenization/10_plot/kmer_all_info','w') as ff:
            ff.write(f"kmer\tsample\tchromosome\tstrand\tstart\tend\thi_entropy\tymin\tymax\n")
        def run_step(one_kmer):
            with open(f'./samples_satellite/18_Homogenization/10_plot/kmer_{one_kmer}_info','w') as f2:
                f2.write(f"sample\tchromosome\tstrand\tstart\tend\thi_entropy\tymin\tymax\n")   ## Note: Hi_entropy
                with open(f'./samples_satellite/18_Homogenization/2_kmer/kmer_{one_kmer}','r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        chromosome,region,segment_file,total_count,dominant_count,unique_count,hi_dominant,hi_diversity,entropy,hi_entropy,Coverage,rare_kmer_ratio=eachline_arr
                        sample=region.split(':')[0]
                        if sample not in ["Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap1","V059.hap2","V063.hap1","V063.hap2","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap1","Chardonnay_hap2","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap1","V074.hap2","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap1","V062.hap2","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap1","PinotNoir_hap2","PinotNoir2_hap1","PinotNoir2_hap2","V069.hap1","V069.hap2","PN40024","PN40024_hap1","PN40024_hap2","V087.hap1","V087.hap2","V061.hap1","V061.hap2","V065.hap1","V065.hap2","V064.hap1","V064.hap2","ThompsonSeedless_hap1","ThompsonSeedless_hap2","ThompsonSeedless2_hap1","ThompsonSeedless2_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap1","V067.hap2","V058.hap1","V058.hap2","V107.hap1","V107.hap2","V108.hap1","V108.hap2","V081.hap1","V081.hap2","V015.hap1","V015.hap2","V105.hap1","V105.hap2","V126.hap1","V126.hap2","V112.hap1","V112.hap2","DavidiiGrape_hap1","DavidiiGrape_hap2","V007.hap1","V007.hap2","V008.hap1","V008.hap2","V012.hap1","V012.hap2","V100.hap1","V100.hap2","V099.hap1","V099.hap2","V098.hap1","V098.hap2","V117.hap1","V117.hap2","V120.hap1","V120.hap2","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap1","V124.hap2","V125.hap1","V125.hap2","WoollyGrape_hap1","WoollyGrape_hap2","V106.hap1","V106.hap2","V102.hap1","V102.hap2","V079.hap1","V079.hap2","V030.hap1","V030.hap2","V031.hap1","V031.hap2","V005.hap1","V005.hap2","V048.hap1","V048.hap2","V050.hap1","V050.hap2","V051.hap1","V051.hap2","V049.hap1","V049.hap2","V038.hap1","V038.hap2","V043.hap1","V043.hap2","V052.hap1","V052.hap2","V034.hap1","V034.hap2","V023.hap1","V023.hap2","V032.hap1","V032.hap2","V033.hap1","V033.hap2","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap1","V022.hap2","V040.hap1","V040.hap2","V019.hap1","V019.hap2","V041.hap1","V041.hap2","V077.hap1","V077.hap2","V037.hap1","V037.hap2","V096.hap1","V096.hap2","ShineMuscat_hap1","ShineMuscat_hap2","V036.hap1","V036.hap2","V076.hap1","V076.hap2","V055.hap1","V055.hap2","V018.hap1","V018.hap2","V072.hap1","V072.hap2"]:continue
                        strand=region.split(":")[1][-1]
                        pos_arr=segment_file.split('___')[1].split('-')
                        start=pos_arr[0]
                        end=pos_arr[1]
                        ymin='' ;ymax=''
                        if strand=='+': ymin=0 ;ymax=1
                        elif strand=='-': ymin=-1 ;ymax=0
                        f2.write(f"{sample}\t{chromosome}\t{strand}\t{start}\t{end}\t{hi_entropy}\t{ymin}\t{ymax}\n")
                        #####
                        with open(f'./samples_satellite/18_Homogenization/10_plot/kmer_all_info','a') as ff:
                            ff.write(f"{one_kmer}\t{sample}\t{chromosome}\t{strand}\t{start}\t{end}\t{hi_entropy}\t{ymin}\t{ymax}\n")
            Plot_txt=f"""
            library(ggplot2)
            library(dplyr)
            # Get all command line arguments
            
            ### Monomer
            print("")

              input_file=read.table('kmer_{one_kmer}_info', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
    
              input_file$chromosome <- factor(
                  input_file$chromosome,
                  levels = c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10','Chr11','Chr12','Chr13','Chr14','Chr15','Chr16','Chr17','Chr18','Chr19')
                )  
    input_file$sample <- factor(input_file$sample, levels = c(
            "Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap1","V059.hap2","V063.hap1","V063.hap2","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap1","Chardonnay_hap2","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap1","V074.hap2","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap1","V062.hap2","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap1","PinotNoir_hap2","PinotNoir2_hap1","PinotNoir2_hap2","V069.hap1","V069.hap2","PN40024","PN40024_hap1","PN40024_hap2","V087.hap1","V087.hap2","V061.hap1","V061.hap2","V065.hap1","V065.hap2","V064.hap1","V064.hap2","ThompsonSeedless_hap1","ThompsonSeedless_hap2","ThompsonSeedless2_hap1","ThompsonSeedless2_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap1","V067.hap2","V058.hap1","V058.hap2","V107.hap1","V107.hap2","V108.hap1","V108.hap2","V081.hap1","V081.hap2","V015.hap1","V015.hap2","V105.hap1","V105.hap2","V126.hap1","V126.hap2","V112.hap1","V112.hap2","DavidiiGrape_hap1","DavidiiGrape_hap2","V007.hap1","V007.hap2","V008.hap1","V008.hap2","V012.hap1","V012.hap2","V100.hap1","V100.hap2","V099.hap1","V099.hap2","V098.hap1","V098.hap2","V117.hap1","V117.hap2","V120.hap1","V120.hap2","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap1","V124.hap2","V125.hap1","V125.hap2","WoollyGrape_hap1","WoollyGrape_hap2","V106.hap1","V106.hap2","V102.hap1","V102.hap2","V079.hap1","V079.hap2","V030.hap1","V030.hap2","V031.hap1","V031.hap2","V005.hap1","V005.hap2","V048.hap1","V048.hap2","V050.hap1","V050.hap2","V051.hap1","V051.hap2","V049.hap1","V049.hap2","V038.hap1","V038.hap2","V043.hap1","V043.hap2","V052.hap1","V052.hap2","V034.hap1","V034.hap2","V023.hap1","V023.hap2","V032.hap1","V032.hap2","V033.hap1","V033.hap2","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap1","V022.hap2","V040.hap1","V040.hap2","V019.hap1","V019.hap2","V041.hap1","V041.hap2","V077.hap1","V077.hap2","V037.hap1","V037.hap2","V096.hap1","V096.hap2","ShineMuscat_hap1","ShineMuscat_hap2","V036.hap1","V036.hap2","V076.hap1","V076.hap2","V055.hap1","V055.hap2","V018.hap1","V018.hap2","V072.hap1","V072.hap2"))
             # Check data and handle missing values
            input_file <- input_file %>%
              filter(!is.na(start), !is.na(end), !is.na(hi_entropy)) %>%  # Filter invalid rows
              group_by(sample, chromosome) %>%
              mutate(
                ymin = as.numeric(factor(sample)) - 0.4,  # Assign y-axis position for each sample
                ymax = as.numeric(factor(sample)) + 0.4
              )
            
            # Create plot object
            p <- ggplot() +
              geom_rect(
                data = input_file,
                aes(
                  xmin = start / 1000000,  # Convert to Mb
                  xmax = end / 1000000,
                  ymin = ymin,
                  ymax = ymax,
                  fill = hi_entropy  # Fill with hi_entropy
                )
              ) +
              scale_fill_gradientn(
                name = "Hi-C Entropy",
                colors = c("#440154", "#404387", "#29788e", "#23a883", "#7ad150", "#fde724"),
                limits = range(input_file$hi_entropy, na.rm = TRUE),
                values = scales::rescale(quantile(input_file$hi_entropy, probs = seq(0, 1, length.out = 6), na.rm = TRUE))
              ) +
              guides(fill = guide_colorbar(title.position = "top", barwidth = 10)) +
              facet_grid(sample ~ chromosome,  drop = FALSE) +
              theme_classic() +
              theme(
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                strip.text.y = element_text(angle = 0)  # Display sample names horizontally
              )
            
            # Save as PDF
            pdf(file = "{one_kmer}_sum_for_plot.pdf", width = 400 / 2.54, height = 150 / 2.54)
            print(p)
            dev.off()
            
            
            input_file <- input_file[input_file$sample == "PN40024", ]
          # Create plot object
            p <- ggplot() +
              geom_rect(
                data = input_file,
                aes(
                  xmin = start / 1000000,  # Convert to Mb
                  xmax = end / 1000000,
                  ymin = ymin,
                  ymax = ymax,
                  fill = hi_entropy  # Fill with hi_entropy
                )
              ) +
              scale_fill_gradientn(
                name = "Hi-C Entropy",
                colors = c("#440154", "#404387", "#29788e", "#23a883", "#7ad150", "#fde724"),
                limits = range(input_file$hi_entropy, na.rm = TRUE),
                values = scales::rescale(quantile(input_file$hi_entropy, probs = seq(0, 1, length.out = 6), na.rm = TRUE))
              ) +
              guides(fill = guide_colorbar(title.position = "top", barwidth = 10)) +
              facet_grid(sample ~ chromosome,  drop = FALSE) +
              theme_classic() +
              theme(
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                strip.text.y = element_text(angle = 0)  # Display sample names horizontally
              )
            
            # Save as PDF
            pdf(file = "{one_kmer}_ref_plot.pdf", width = 20 / 2.54, height = 2 / 2.54)
            print(p)
            dev.off()
            """
            with open(f'./samples_satellite/18_Homogenization/10_plot/plot_{one_kmer}.R','w') as f:
                f.write(Plot_txt)
            os.chdir(f'./samples_satellite/18_Homogenization/10_plot/')
            subprocess.run([f'Rscript plot_{one_kmer}.R  '], shell=True)  #>null 2>&1 
            os.chdir('../../../')                 
                
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, kmer_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(kmer_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()         
    if argv1=="stepall" or argv1=="step18"  or argv1=="step18.10_ref":             
        Plot_txt=f"""
            library(ggplot2)
            library(dplyr)
            # Get all command line arguments
            
            ### Monomer
            print("")

              input_file=read.table('kmer_all_info', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
    
              input_file$chromosome <- factor(
                  input_file$chromosome,
                  levels = c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10','Chr11','Chr12','Chr13','Chr14','Chr15','Chr16','Chr17','Chr18','Chr19')
                )  
    input_file$sample <- factor(input_file$sample, levels = c(
            "Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap1","V059.hap2","V063.hap1","V063.hap2","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap1","Chardonnay_hap2","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap1","V074.hap2","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap1","V062.hap2","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap1","PinotNoir_hap2","PinotNoir2_hap1","PinotNoir2_hap2","V069.hap1","V069.hap2","PN40024","PN40024_hap1","PN40024_hap2","V087.hap1","V087.hap2","V061.hap1","V061.hap2","V065.hap1","V065.hap2","V064.hap1","V064.hap2","ThompsonSeedless_hap1","ThompsonSeedless_hap2","ThompsonSeedless2_hap1","ThompsonSeedless2_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap1","V067.hap2","V058.hap1","V058.hap2","V107.hap1","V107.hap2","V108.hap1","V108.hap2","V081.hap1","V081.hap2","V015.hap1","V015.hap2","V105.hap1","V105.hap2","V126.hap1","V126.hap2","V112.hap1","V112.hap2","DavidiiGrape_hap1","DavidiiGrape_hap2","V007.hap1","V007.hap2","V008.hap1","V008.hap2","V012.hap1","V012.hap2","V100.hap1","V100.hap2","V099.hap1","V099.hap2","V098.hap1","V098.hap2","V117.hap1","V117.hap2","V120.hap1","V120.hap2","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap1","V124.hap2","V125.hap1","V125.hap2","WoollyGrape_hap1","WoollyGrape_hap2","V106.hap1","V106.hap2","V102.hap1","V102.hap2","V079.hap1","V079.hap2","V030.hap1","V030.hap2","V031.hap1","V031.hap2","V005.hap1","V005.hap2","V048.hap1","V048.hap2","V050.hap1","V050.hap2","V051.hap1","V051.hap2","V049.hap1","V049.hap2","V038.hap1","V038.hap2","V043.hap1","V043.hap2","V052.hap1","V052.hap2","V034.hap1","V034.hap2","V023.hap1","V023.hap2","V032.hap1","V032.hap2","V033.hap1","V033.hap2","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap1","V022.hap2","V040.hap1","V040.hap2","V019.hap1","V019.hap2","V041.hap1","V041.hap2","V077.hap1","V077.hap2","V037.hap1","V037.hap2","V096.hap1","V096.hap2","ShineMuscat_hap1","ShineMuscat_hap2","V036.hap1","V036.hap2","V076.hap1","V076.hap2","V055.hap1","V055.hap2","V018.hap1","V018.hap2","V072.hap1","V072.hap2"))
             # Check data and handle missing values
            input_file <- input_file %>%
              filter(!is.na(start), !is.na(end), !is.na(hi_entropy)) %>%  # Filter invalid rows
              group_by(sample, chromosome) %>%
              mutate(
                ymin = as.numeric(factor(sample)) - 0.4,  # Assign y-axis position for each sample
                ymax = as.numeric(factor(sample)) + 0.4
              )
  
            input_file <- input_file[input_file$sample == "PN40024", ]
          # Create plot object
            p <- ggplot() +
              geom_rect(
                data = input_file,
                aes(
                  xmin = start / 1000000,  # Convert to Mb
                  xmax = end / 1000000,
                  ymin = ymin,
                  ymax = ymax,
                  fill = hi_entropy  # Fill with hi_entropy
                )
              ) +
              scale_fill_gradientn(
                name = "Hi-C Entropy",
                colors = c("#440154", "#404387", "#29788e", "#23a883", "#7ad150", "#fde724"),
                limits = range(input_file$hi_entropy, na.rm = TRUE),
                values = scales::rescale(quantile(input_file$hi_entropy, probs = seq(0, 1, length.out = 6), na.rm = TRUE))
              ) +
              guides(fill = guide_colorbar(title.position = "top", barwidth = 10)) +
              facet_grid(kmer ~ chromosome) +
              theme_classic() +
              theme(
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                strip.text.y = element_text(angle = 0)  # Display sample names horizontally
              )
            
            # Save as PDF
            pdf(file = "kmerall_ref_plot.pdf", width = 80 / 2.54, height = 20 / 2.54)
            print(p)
            dev.off()
            """
        with open(f'./samples_satellite/18_Homogenization/10_plot/plot_kmerall.R','w') as f:
            f.write(Plot_txt)
        os.chdir(f'./samples_satellite/18_Homogenization/10_plot/')
        subprocess.run([f'Rscript plot_kmerall.R  '], shell=True)  #>null 2>&1 
        os.chdir('../../../')                   
    
    #The code for other Chr is similar. 
    if argv1=="stepall" or 'step18.11' in argv1:
        print('Calculate HOR length and HI index')
        subprocess.run(["mkdir ./samples_satellite/18_Homogenization/11_HOR2homo"], shell=True)
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.11"  or argv1=="step18.11s00":
            print('step18.11s00-04: All HOR length obtained is 2G, including HOR overlaps, double-counted')
            print('step18.11s10-11: HOR + noHOR length obtained is 0.8G, because overlaps were removed')
            print('Refer to step18.10, ./samples_satellite/18_Homogenization/10_plot/kmer_all_info is similar to this file')
    
            one_kmer=31
            with open(f'./samples_satellite/18_Homogenization/11_HOR2homo/kmer_info','w') as f2:
                f2.write(f"sample\tchromosome\tstrand\tstart\tend\thi_entropy\tymin\tymax\n")   ## Note: Hi_entropy
                with open(f'./samples_satellite/18_Homogenization/2_kmer/kmer_{one_kmer}','r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        chromosome,region,segment_file,total_count,dominant_count,unique_count,hi_dominant,hi_diversity,entropy,hi_entropy,Coverage,rare_kmer_ratio=eachline_arr
                        sample=region.split(':')[0]
                        if sample not in ["Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap1","V059.hap2","V063.hap1","V063.hap2","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap1","Chardonnay_hap2","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap1","V074.hap2","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap1","V062.hap2","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap1","PinotNoir_hap2","PinotNoir2_hap1","PinotNoir2_hap2","V069.hap1","V069.hap2","PN40024","PN40024_hap1","PN40024_hap2","V087.hap1","V087.hap2","V061.hap1","V061.hap2","V065.hap1","V065.hap2","V064.hap1","V064.hap2","ThompsonSeedless_hap1","ThompsonSeedless_hap2","ThompsonSeedless2_hap1","ThompsonSeedless2_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap1","V067.hap2","V058.hap1","V058.hap2","V107.hap1","V107.hap2","V108.hap1","V108.hap2","V081.hap1","V081.hap2","V015.hap1","V015.hap2","V105.hap1","V105.hap2","V126.hap1","V126.hap2","V112.hap1","V112.hap2","DavidiiGrape_hap1","DavidiiGrape_hap2","V007.hap1","V007.hap2","V008.hap1","V008.hap2","V012.hap1","V012.hap2","V100.hap1","V100.hap2","V099.hap1","V099.hap2","V098.hap1","V098.hap2","V117.hap1","V117.hap2","V120.hap1","V120.hap2","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap1","V124.hap2","V125.hap1","V125.hap2","WoollyGrape_hap1","WoollyGrape_hap2","V106.hap1","V106.hap2","V102.hap1","V102.hap2","V079.hap1","V079.hap2","V030.hap1","V030.hap2","V031.hap1","V031.hap2","V005.hap1","V005.hap2","V048.hap1","V048.hap2","V050.hap1","V050.hap2","V051.hap1","V051.hap2","V049.hap1","V049.hap2","V038.hap1","V038.hap2","V043.hap1","V043.hap2","V052.hap1","V052.hap2","V034.hap1","V034.hap2","V023.hap1","V023.hap2","V032.hap1","V032.hap2","V033.hap1","V033.hap2","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap1","V022.hap2","V040.hap1","V040.hap2","V019.hap1","V019.hap2","V041.hap1","V041.hap2","V077.hap1","V077.hap2","V037.hap1","V037.hap2","V096.hap1","V096.hap2","ShineMuscat_hap1","ShineMuscat_hap2","V036.hap1","V036.hap2","V076.hap1","V076.hap2","V055.hap1","V055.hap2","V018.hap1","V018.hap2","V072.hap1","V072.hap2"]:continue
                        strand=region.split(":")[1][-1]
                        pos_arr=segment_file.split('___')[1].split('-')
                        start1=pos_arr[0]
                        end1=pos_arr[1]
                        start=min(start1,end1)
                        end=max(start1,end1)
                        ymin='' ;ymax=''
                        if strand=='+': ymin=0 ;ymax=1
                        elif strand=='-': ymin=-1 ;ymax=0
                        f2.write(f"{sample}\t{chromosome}\t{strand}\t{start}\t{end}\t{hi_entropy}\t{ymin}\t{ymax}\n")
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.11"  or argv1=="step18.11s0":    
            
            input_list=[]
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/kmer_info','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if len(eachline_arr)!=8:print(line)
                    sample,chromosome,strand,start,end,hi_entropy,ymin,ymax=eachline_arr
                    #if sample!="MuscatHamburg_hap1" or chromosome!=	"Chr19":continue
    
                    input_list.append([sample,chromosome,strand,start,end,hi_entropy])
            #print(input_list)        
            #sys.exit()
            region_name_list=[]
            dict_region_samplechr={}
            with open('./samples_satellite/2_good_regions','r') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if eachline_arr[0]=='sample':continue
                    if eachline_arr[0]=='':continue
                    chromosome=eachline_arr[1].replace('region_','Chr')
                    sample=eachline_arr[0]
                    if sample not in ["Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap1","V059.hap2","V063.hap1","V063.hap2","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap1","Chardonnay_hap2","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap1","V074.hap2","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap1","V062.hap2","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap1","PinotNoir_hap2","PinotNoir2_hap1","PinotNoir2_hap2","V069.hap1","V069.hap2","PN40024","PN40024_hap1","PN40024_hap2","V087.hap1","V087.hap2","V061.hap1","V061.hap2","V065.hap1","V065.hap2","V064.hap1","V064.hap2","ThompsonSeedless_hap1","ThompsonSeedless_hap2","ThompsonSeedless2_hap1","ThompsonSeedless2_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap1","V067.hap2","V058.hap1","V058.hap2","V107.hap1","V107.hap2","V108.hap1","V108.hap2","V081.hap1","V081.hap2","V015.hap1","V015.hap2","V105.hap1","V105.hap2","V126.hap1","V126.hap2","V112.hap1","V112.hap2","DavidiiGrape_hap1","DavidiiGrape_hap2","V007.hap1","V007.hap2","V008.hap1","V008.hap2","V012.hap1","V012.hap2","V100.hap1","V100.hap2","V099.hap1","V099.hap2","V098.hap1","V098.hap2","V117.hap1","V117.hap2","V120.hap1","V120.hap2","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap1","V124.hap2","V125.hap1","V125.hap2","WoollyGrape_hap1","WoollyGrape_hap2","V106.hap1","V106.hap2","V102.hap1","V102.hap2","V079.hap1","V079.hap2","V030.hap1","V030.hap2","V031.hap1","V031.hap2","V005.hap1","V005.hap2","V048.hap1","V048.hap2","V050.hap1","V050.hap2","V051.hap1","V051.hap2","V049.hap1","V049.hap2","V038.hap1","V038.hap2","V043.hap1","V043.hap2","V052.hap1","V052.hap2","V034.hap1","V034.hap2","V023.hap1","V023.hap2","V032.hap1","V032.hap2","V033.hap1","V033.hap2","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap1","V022.hap2","V040.hap1","V040.hap2","V019.hap1","V019.hap2","V041.hap1","V041.hap2","V077.hap1","V077.hap2","V037.hap1","V037.hap2","V096.hap1","V096.hap2","ShineMuscat_hap1","ShineMuscat_hap2","V036.hap1","V036.hap2","V076.hap1","V076.hap2","V055.hap1","V055.hap2","V018.hap1","V018.hap2","V072.hap1","V072.hap2"]:continue
                    one_region_name=f"{eachline_arr[0]}:{eachline_arr[1]}{eachline_arr[6]}:{eachline_arr[3]}-{eachline_arr[4]}"
                    #if one_region_name!='MuscatHamburg_hap1:region_19-:16767395-18114955':continue
                    region_name_list.append(one_region_name)
                    dict_region_samplechr[one_region_name]=[eachline_arr[0],chromosome]
                    
            ALL_LENGTH=0
            dict_samplechr_info={}     
            ALL_MARK_LEN=0
            for one_region in  region_name_list:    
                sample,chromosome=dict_region_samplechr[one_region]
                #if sample!="MuscatHamburg_hap1" or chromosome!=	"Chr19":continue
                samplechr=f"{sample}|{chromosome}"
                if samplechr not in dict_samplechr_info:
                    dict_samplechr_info[samplechr]=[]
                with open(f'./samples_satellite/12_Hor_stat/1_blocks_stat//{one_region}','r')   as f:
                    next(f)
                    for line in f:
                        #print(line)
                        eachline_arr=line.strip().split('\t')
                        mer,segment,segment_type,segment_type_num,segment_length,HOR_type,circ_start,circ_end,pos_start1,pos_end1,markserial_start,markserial_end,layer,father_layer,pos_length,markserial_num,HOR_repeat_num=eachline_arr
                        pos_start1=int(pos_start1)
                        pos_end1=int(pos_end1)
                        pos_start=min(pos_start1,pos_end1)
                        pos_end=max(pos_start1,pos_end1)
                        pos_length=int(pos_length)
                        #if float(HOR_repeat_num)==150:print(eachline_arr)
                        ALL_LENGTH+=pos_length
                        dict_samplechr_info[samplechr].append([mer,segment,segment_type,segment_length,HOR_type,pos_start,pos_end,pos_length,markserial_num,HOR_repeat_num])
                        #if segment=='28-23|28-28-23|28-28-28-23|28-28-23':
                        #    ALL_MARK_LEN+=pos_length
                        #    if pos_length>1000:print(samplechr,eachline_arr)
            print(f'ALL_LENGTH:{ALL_LENGTH}')
            #print(f'ALL_MARK_LEN:{ALL_MARK_LEN}')
            
            def run_step(one):    
                sample,chromosome,strand,start,end,hi_entropy=one
                #print(one)
                start=int(start)
                end=int(end)
                region_len=end-start+1
                samplechr=f"{sample}|{chromosome}"
                HOR_list=dict_samplechr_info[samplechr]
                good_HOR_list=[]
                for one_HOR in HOR_list:
                    #if segment=='28-23|28-28-23|28-28-28-23|28-28-23':print(one_HOR)
                    mer,segment,segment_type,segment_length,HOR_type,pos_start1,pos_end1,pos_length,markserial_num,HOR_repeat_num=one_HOR
                    #if float(HOR_repeat_num)==150:print(one_HOR)
                    pos_start=int(pos_start1)
                    pos_end=int(pos_end1)
                    #if pos_end-pos_start<0:print(11111111111)
                    one_HOR_len=pos_end-pos_start+1
                    if pos_start>=end or pos_end<=start:continue
                    elif pos_start>=start and pos_end<=end:
                        good_HOR_list.append([one,one_HOR,1]);#print(1)
                        
                    elif pos_start<=start and pos_end>end:  
                        percent=region_len/one_HOR_len
                        good_HOR_list.append([one,one_HOR,percent]);#print(1)
                    elif pos_start<=start and pos_end>=start and pos_end<=end :      
                        percent=(pos_end-start+1)/one_HOR_len
                        good_HOR_list.append([one,one_HOR,percent]);#print(1)
                    elif pos_start>=start and pos_start<=end and pos_end>=end :      
                        percent=(end-pos_start+1)/one_HOR_len     
                        good_HOR_list.append([one,one_HOR,percent]);#print(1)
                    else:
                        print('error, this situation should not occur')
                #print(HOR_list)        
                return good_HOR_list 
                    
            
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/0_sum','w') as f:
                f.write(f"sample\tchromosome\tstrand\tstart\tend\thi_entropy\tmer\tsegment\tsegment_type\tsegment_length\tHOR_type\tpos_start\tpos_end\tpos_length\tmarkserial_num\tHOR_repeat_num\tpercent\n")
                with Pool(processes=thread) as pool:
                    # Use imap to get results one by one
                    for i, result in enumerate(pool.imap(run_step, input_list), start=1):
                        # Results can be processed here, e.g., stored or printed
                        for one_result in result:
                            HIregion,one_HOR,percent=one_result
                            ####
                            sample,chromosome,strand,Start,end,hi_entropy=HIregion
                            mer,segment,segment_type,segment_length,HOR_type,pos_start,pos_end,pos_length,markserial_num,HOR_repeat_num=one_HOR
                            #if segment=='28-23|28-28-23|28-28-28-23|28-28-23':print(one_HOR)
                            #####
                            f.write(f"{sample}\t{chromosome}\t{strand}\t{Start}\t{end}\t{hi_entropy}\t{mer}\t{segment}\t{segment_type}\t{segment_length}\t{HOR_type}\t{pos_start}\t{pos_end}\t{pos_length}\t{markserial_num}\t{HOR_repeat_num}\t{percent}\n")
                            
                        ###########    
                        progress = (i / len(input_list)) * 100
                        sys.stdout.write(f"\rProgress: {progress:.2f}%")
                        sys.stdout.flush()   
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.11"  or argv1=="step18.11s1":  
            def normalize_hor_pattern_rotation_only(sequence, separator='|'):
                """
                Normalize HOR monomer sequence, only handling rotation, not considering reverse
                
                Example: "28-28-28-23|28-28-23|28-28-28-28-23"
                   -> Normalize to unique string
                
                Strategy: Find the smallest lexicographic order among all rotations
                """
                # 1. Split by separator
                parts = sequence.split(separator)
                n = len(parts)
                
                # 2. Generate all rotations
                candidates = []
                doubled = parts + parts  # Double sequence for rotation
                
                for i in range(n):
                    rotated = doubled[i:i+n]
                    candidate = separator.join(rotated)
                    candidates.append(candidate)
                
                # 3. Return the smallest lexicographic string
                return min(candidates)    
                
            dict_HORtype_length={}
            dict_segmentlen_length={}
            dict_mer_length={}
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/1_sum_filter','w') as f2:
                f2.write(f"horname\tsample\tchromosome\tstrand\tstart\tend\thi_entropy\tmer\tsegment\tsegment_type\tsegment_length\tHOR_type\tpos_start\tpos_end\tpos_length\tmarkserial_num\tHOR_repeat_num\tpercent\n")
                with open('./samples_satellite/18_Homogenization/11_HOR2homo/0_sum','r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t')            
                        sample,chromosome,strand,start,end,hi_entropy,mer,segment,segment_type,segment_length,HOR_type,pos_start,pos_end,pos_length,markserial_num,HOR_repeat_num,percent=eachline_arr
    
                        ##############
    
                        #### Count mer sum
                        dict_change = {
                                "28-23|28-28-28-23": "MON51-MON107",
                                "28-28-23|28-28-28-23": "MON79-MON107",
                                "28-28-28-23|28-28-28-23": "(MON107)2",
                                "28-28-28-23|28-28-28-28-23": "MON107-MON135",
                                "28-28-28-23|28-28-28-28-28-23": "MON107-MON163",
                                "28-28-28-28-23|28-28-28-28-23": "(MON135)2",
                                "28-23|28-28-28-23|28-28-28-23": "MON51-(MON107)2",           # After normalization
                                "28-28-23|28-28-23|28-28-28-23": "(MON79)2-MON107",      # After normalization
                                "28-28-23|28-28-28-23|28-28-28-23": "MON79-(MON107)2",      # After normalization (note difference from above)
                                "28-28-23|28-28-28-23|28-28-28-28-23": "MON79-MON107-MON135",
                                "28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)3",
                                "28-28-28-28-23|28-28-28-28-23|28-28-28-28-23": "(MON135)3",
                                "28-23|28-28-23|28-28-28-23|28-28-23": "MON51-(MON79)2-MON107",  # After normalization
                                "28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)4",
                                "28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)5",
                                "28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)6",
                                "28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)10"
                            }
                        if segment not in dict_change: name_new='other'
                        else: name_new=   dict_change[segment] 
                        #if segment=='28-23|28-28-23|28-28-28-23|28-28-23':print(eachline_arr)
                        #if name_new=='MON51-(MON79)2-MON107':print(eachline_arr)
                        f2.write(f"{name_new}\t{sample}\t{chromosome}\t{strand}\t{start}\t{end}\t{hi_entropy}\t{mer}\t{segment}\t{segment_type}\t{segment_length}\t{HOR_type}\t{pos_start}\t{pos_end}\t{pos_length}\t{markserial_num}\t{HOR_repeat_num}\t{percent}\n")
                        ###
                        length_good=int(pos_length)*float(percent)
                        if mer not in dict_mer_length: dict_mer_length[mer]=[]
                        dict_mer_length[mer].append([length_good,hi_entropy])
                        if segment_length not in dict_segmentlen_length: dict_segmentlen_length[segment_length]=[]
                        dict_segmentlen_length[segment_length].append([length_good,hi_entropy])
                        if name_new not in dict_HORtype_length:dict_HORtype_length[name_new]=[]
                        dict_HORtype_length[name_new].append([length_good,hi_entropy])
            def get_entropy_range(hi_entropy):
                """
                Generate interval range based on hi_entropy value
                Example: 0.125 -> "0.12-0.13"
                     0.111 -> "0.11-0.12"
                     0.119 -> "0.11-0.12"
                     0.120 -> "0.12-0.13"
                """
                hi_entropy = float(hi_entropy)
                # Round down to 0.01
                lower = math.floor(hi_entropy * 100) / 100
                upper = lower + 0.01
                return f"{lower:.2f}-{upper:.2f}"       
         
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/1_plot_mer','w') as f:
                f.write('mer\tlength\tHIindex\tHIentropy\tHIentropy_class\n')
                for  mer,length_HIentropy in  dict_mer_length.items():     
                    for one in length_HIentropy:
                        length,HIentropy=one
                        HIentropy_class=get_entropy_range(HIentropy)
                        f.write(f"{mer}\t{length}\t{HIentropy}\t{HIentropy_class}\n")
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/1_plot_segmentlen','w') as f:
                f.write('segmentlen\tlength\tHIentropy\tHIentropy_class\n')
                for  segmentlen,length_HIentropy in  dict_segmentlen_length.items():          
                    for one in length_HIentropy:
                        length,HIentropy=one
                        HIentropy_class=get_entropy_range(HIentropy)
                        f.write(f"{segmentlen}\t{length}\t{HIentropy}\t{HIentropy_class}\n")
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/1_plot_HORtype','w') as f:
                f.write('HORtype\tlength\tHIentropy\tHIentropy_class\n')
                for  HORtype,length_HIentropy in  dict_HORtype_length.items():          
                    for one in length_HIentropy:
                        length,HIentropy=one
                        HIentropy_class=get_entropy_range(HIentropy)
                        f.write(f"{HORtype}\t{length}\t{HIentropy}\t{HIentropy_class}\n")              
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.11"  or argv1=="step18.11s2":     
            # Read and output by HIentropy_class
            dict_class_len={}
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/1_plot_mer', 'r') as f:
                header = f.readline()
                for line in f:
                    parts = line.strip().split('\t')
                    mer, length,  HIentropy, HIentropy_class = parts
                    newclass=f"{HIentropy_class}|||||{mer}"
                    if newclass not in dict_class_len:
                        dict_class_len[newclass]=0
                    dict_class_len[newclass]   +=float(length)
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/2_plot_mer', 'w') as f:     
                f.write('mer\tHIentropy_class\tlength\n')
                for one,length in dict_class_len.items():
                    HIentropy_class,mer=one.split('|||||')
                    f.write(f"{mer}\t{HIentropy_class}\t{length:.0f}\n")
            ###############
            dict_class_len={}
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/1_plot_segmentlen', 'r') as f:
                header = f.readline()
                for line in f:
                    parts = line.strip().split('\t')
                    segmentlen, length,  HIentropy, HIentropy_class = parts
                    newclass=f"{HIentropy_class}|||||{segmentlen}"
                    if newclass not in dict_class_len:
                        dict_class_len[newclass]=0
                    dict_class_len[newclass]   +=float(length)
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/2_plot_segmentlen', 'w') as f:     
                f.write('segmentlen\tHIentropy_class\tlength\n')
                for one,length in dict_class_len.items():
                    HIentropy_class,segmentlen=one.split('|||||')
                    f.write(f"{segmentlen}\t{HIentropy_class}\t{length:.0f}\n")        
            ###############
            dict_class_len={}
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/1_plot_HORtype', 'r') as f:
                header = f.readline()
                for line in f:
                    parts = line.strip().split('\t')
                    HORtype, length,  HIentropy, HIentropy_class = parts
                    newclass=f"{HIentropy_class}|||||{HORtype}"
                    if newclass not in dict_class_len:
                        dict_class_len[newclass]=0
                    dict_class_len[newclass]   +=float(length)
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/2_plot_HORtype', 'w') as f:     
                f.write('HORtype\tHIentropy_class\tlength\n')
                for one,length in dict_class_len.items():
                    HIentropy_class,HORtype=one.split('|||||')
                    f.write(f"{HORtype}\t{HIentropy_class}\t{length:.0f}\n")        
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.11"  or argv1=="step18.11s3":            
            print('Plotting')            
            R_txt=f'''library(ggplot2)
    
    # Set working directory
    setwd('./')
    
    # Read data
    input_file1 <- read.table('2_plot_mer', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')
    
    # Create plot
    p <- ggplot()
    p <- p + geom_point(data = input_file1, aes(x = mer, y = length,color=HIentropy_class))
    p <- p + theme_classic()
    # Add title and axis labels
    p <- p + theme(
      #axis.title = element_blank(),  # Hide axis titles
      #axis.text = element_blank(),   # Hide axis tick text
      #legend.position = "none",
      #axis.ticks = element_blank(),  # Hide axis ticks
      #axis.line = element_blank()    # Hide axis lines
    )+ coord_flip()
    pdf("3_plot_mer.pdf", width = 20 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()
    '''
            with open(f'./samples_satellite/18_Homogenization/11_HOR2homo/plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)
            new_directory = f"./samples_satellite/18_Homogenization/11_HOR2homo/"
            os.chdir(new_directory)
            subprocess.run(['Rscript plot.R'], shell=True)   
            
            print('Plotting')            
            R_txt=f'''library(ggplot2)
    
    # Set working directory
    setwd('./')
    
    # Read data
    input_file1 <- read.table('2_plot_HORtype', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')
    
    # Create plot
    p <- ggplot()
    p <- p + geom_point(data = input_file1, aes(x = HORtype, y = length,color=HIentropy_class))
    p <- p + theme_classic()
    # Add title and axis labels
    p <- p + theme(
      #axis.title = element_blank(),  # Hide axis titles
      #axis.text = element_blank(),   # Hide axis tick text
      #legend.position = "none",
      #axis.ticks = element_blank(),  # Hide axis ticks
      #axis.line = element_blank()    # Hide axis lines
    )+ coord_flip()
    pdf("3_plot_HORtype.pdf", width = 20 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()
    '''
            with open(f'plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)
            subprocess.run(['Rscript plot.R'], shell=True)           
            
            
            print('Plotting')            
            R_txt=f'''library(ggplot2)
    
    # Set working directory
    setwd('./')
    
    # Read data
    input_file1 <- read.table('2_plot_segmentlen', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')
    
    # Create plot
    p <- ggplot()
    p <- p + geom_point(data = input_file1, aes(x = segmentlen, y = length,color=HIentropy_class))
    p <- p + theme_classic()
    # Add title and axis labels
    p <- p + theme(
      #axis.title = element_blank(),  # Hide axis titles
      #axis.text = element_blank(),   # Hide axis tick text
      #legend.position = "none",
      #axis.ticks = element_blank(),  # Hide axis ticks
      #axis.line = element_blank()    # Hide axis lines
    )+ coord_flip()
    pdf("3_plot_segmentlen.pdf", width = 20 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()
    '''
            with open(f'plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)
            subprocess.run(['Rscript plot.R'], shell=True)           
                    
            
            
            
            
            
            
            
            os.chdir('../../../')     
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.11"  or argv1=="step18.11s4": 
            print('Convert to wide format')
            dict_HORtype_HIClass={}
            HORtype_set=set()
            HIentropy_class_set=set()
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/2_plot_HORtype','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    HORtype,HIentropy_class,length=eachline_arr
                    dict_HORtype_HIClass[f"{HORtype}||{HIentropy_class}"]=length
                    HORtype_set.add(HORtype)
                    HIentropy_class_set.add(HIentropy_class)
            HORtype_list=list(HORtype_set);  HORtype_list.sort()
            HIentropy_class_list=list(HIentropy_class_set);  HIentropy_class_list.sort()
            print(HORtype_list,HIentropy_class_list)
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/4_plot_HORtype','w') as f:
                head='one_type\t'+'\t'.join(HIentropy_class_list)+'\n'
                f.write(head)
                for one_HORtype in HORtype_list:
                    line_arr=[]
                    for one_HIclass in HIentropy_class_list:
                        name = f"{one_HORtype}||{one_HIclass}"
                        if name not in dict_HORtype_HIClass: length=0
                        else :length=dict_HORtype_HIClass[name]
                        line_arr.append(length)
                    line_str=    one_HORtype+'\t'+'\t'.join(str(x) for x in line_arr)+'\n'
                    f.write(line_str)
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.11"  or argv1=="step18.11s10":
            print('Re-analyze, HOR and non-HOR')
            def get_entropy_range(hi_entropy):
                """
                Generate interval range based on hi_entropy value
                Example: 0.125 -> "0.12-0.13"
                     0.111 -> "0.11-0.12"
                     0.119 -> "0.11-0.12"
                     0.120 -> "0.12-0.13"
                """
                hi_entropy = float(hi_entropy)
                # Round down to 0.01
                lower = math.floor(hi_entropy * 100) / 100
                upper = lower + 0.01
                return f"{lower:.2f}-{upper:.2f}"  
            input_list=[]
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/kmer_info','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    sample,chromosome,strand,start,end,hi_entropy,ymin,ymax=eachline_arr
                    hi_entropy_class=get_entropy_range(hi_entropy)
                    input_list.append([sample,chromosome,strand,start,end,hi_entropy,hi_entropy_class])
         
            region_name_list=[]
            dict_region_samplechr={}
            with open('./samples_satellite/2_good_regions','r') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if eachline_arr[0]=='sample':continue
                    if eachline_arr[0]=='':continue
                    chromosome=eachline_arr[1].replace('region_','Chr')
                    one_region_name=f"{eachline_arr[0]}:{eachline_arr[1]}{eachline_arr[6]}:{eachline_arr[3]}-{eachline_arr[4]}"
                    region_name_list.append(one_region_name)
                    dict_region_samplechr[one_region_name]=[eachline_arr[0],chromosome,int(eachline_arr[3]),int(eachline_arr[4])]
                    
            def get_noHOR_regions(region_start, region_end, HOR_list):
                """
                Get regions within the interval that are not covered by any HOR
                
                Parameters:
                    region_start: Region start position
                    region_end: Region end position
                    HOR_list: List, each element is [pos_start, pos_end]
                
                Returns:
                    noHOR_list: List, each element is [start, end], representing areas with no HOR coverage
                """
                if not HOR_list:
                    return [[region_start, region_end, 'noHOR']]
                
                # 1. Sort by start position
                sorted_hor = sorted(HOR_list, key=lambda x: x[0])
                
                # 2. Merge overlapping or adjacent HORs (adjacent also merged because there's no gap)
                merged_hor = []
                for start, end,_ in sorted_hor:
                    if not merged_hor:
                        merged_hor.append([start, end])
                    else:
                        last_start, last_end = merged_hor[-1]
                        # If overlapping or adjacent (end+1 >= last_start indicates overlap or adjacency)
                        if start <= last_end + 1:
                            merged_hor[-1][1] = max(last_end, end)
                        else:
                            merged_hor.append([start, end])
                
                # 3. Find non-HOR regions
                noHOR_list = []
                current = region_start
                
                for hor_start, hor_end in merged_hor:
                    # HOR completely to the right of region, exit
                    if hor_start > region_end:
                        break
                    # HOR completely to the left of region, skip
                    if hor_end < region_start:
                        continue
                    
                    # Non-HOR region before current HOR
                    if current < hor_start:
                        no_start = max(current, region_start)
                        no_end = min(hor_start - 1, region_end)
                        if no_start <= no_end:
                            noHOR_list.append([no_start, no_end,'noHOR'])
                    
                    # Update current position to after HOR end
                    current = max(current, hor_end + 1)
                
                # 4. Region after the last HOR
                if current <= region_end:
                    noHOR_list.append([current, region_end,'noHOR'])
                
                return noHOR_list               
            def merge_HOR_list(HOR_list):
                """
                Merge overlapping or adjacent regions in HOR_list
                
                Parameters:
                    HOR_list: List, each element is [pos_start, pos_end, 'HOR'] or [pos_start, pos_end]
                
                Returns:
                    merged_list: Merged list, each element is [pos_start, pos_end, 'HOR']
                """
                if not HOR_list:
                    return []
                
                # 1. Extract coordinates (compatible with both formats)
                coords = []
                for item in HOR_list:
                    if len(item) == 3:
                        coords.append([item[0], item[1]])
                    else:
                        coords.append([item[0], item[1]])
                
                # 2. Sort by start position
                sorted_hor = sorted(coords, key=lambda x: x[0])
                
                # 3. Merge overlapping or adjacent regions
                merged = []
                for start, end in sorted_hor:
                    if not merged:
                        merged.append([start, end])
                    else:
                        last_start, last_end = merged[-1]
                        # If overlapping or adjacent (end+1 >= last_start indicates overlap or adjacency)
                        if start <= last_end + 1:
                            merged[-1][1] = max(last_end, end)
                        else:
                            merged.append([start, end])
                
                # 4. Add 'HOR' label
                merged_list = [[start, end, 'HOR'] for start, end in merged]
                
                return merged_list
                
            dict_samplechr_info={}       
            for one_region in  region_name_list:    
                sample,chromosome,region_start,region_end=dict_region_samplechr[one_region]
                samplechr=f"{sample}|{chromosome}"
                #if samplechr not in dict_samplechr_info:
                 #   dict_samplechr_info[samplechr]=[]
                    
                HOR_list=[]    
                with open(f'./samples_satellite/12_Hor_stat/1_blocks_stat//{one_region}','r')   as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        mer,segment,segment_type,segment_type_num,segment_length,HOR_type,circ_start,circ_end,pos_start1,pos_end1,markserial_start,markserial_end,layer,father_layer,pos_length,markserial_num,HOR_repeat_num=eachline_arr
    
                        pos_start1=int(pos_start1)
                        pos_end1=int(pos_end1)
                        pos_start=min(pos_start1,pos_end1)
                        pos_end=max(pos_start1,pos_end1)
                        #if float(HOR_repeat_num)<3:continue 
                        #if int(mer)<2:continue
                        HOR_list.append([pos_start,pos_end,'HOR'])        
                HOR_list= merge_HOR_list(HOR_list)
                
                noHOR_list=get_noHOR_regions(region_start, region_end, HOR_list)        
                all_HOR_noHOR_list=HOR_list+noHOR_list
                dict_samplechr_info[samplechr]=all_HOR_noHOR_list
                
                
            def run_step(one):    
                sample, chromosome, strand, start, end, hi_entropy, hi_entropy_class = one
                start = int(start)
                end = int(end)
                samplechr = f"{sample}|{chromosome}"
                all_HOR_noHOR_list = dict_samplechr_info.get(samplechr, [])
                good_list = []
                
                for one_HOR_noHOR in all_HOR_noHOR_list:
                    pos_start, pos_end, one_type = one_HOR_noHOR
                    pos_start = int(pos_start)
                    pos_end = int(pos_end)
                    
                    # Calculate overlap length
                    overlap_start = max(start, pos_start)
                    overlap_end = min(end, pos_end)
                    
                    if overlap_start > overlap_end:
                        continue  # No overlap
                    
                    overlap_len = overlap_end - overlap_start + 1
                    hor_len = pos_end - pos_start + 1
                    percent = overlap_len / hor_len
                    
                    good_list.append([one, one_HOR_noHOR, percent])
                
                return good_list
                    
            
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/10_sum','w') as f:
                f.write(f"sample\tchromosome\tstrand\tstart\tend\thi_entropy\thi_entropy_class\tpos_start\tpos_end\tpos_length\tone_type\tpercent\n")
                with Pool(processes=thread) as pool:
                    # Use imap to get results one by one
                    for i, result in enumerate(pool.imap(run_step, input_list), start=1):
                        # Results can be processed here, e.g., stored or printed
                        for one_result in result:
                            HIregion,one_HOR_noHOR,percent=one_result
                            ####
                            sample,chromosome,strand,start,end,hi_entropy,hi_entropy_class=HIregion
                            pos_start,pos_end,one_type=one_HOR_noHOR
                            pos_length=pos_end-pos_start+1
                            #####
                            f.write(f"{sample}\t{chromosome}\t{strand}\t{start}\t{end}\t{hi_entropy}\t{hi_entropy_class}\t{pos_start}\t{pos_end}\t{pos_length}\t{one_type}\t{percent}\n")
                            
                        ###########    
                        progress = (i / len(input_list)) * 100
                        sys.stdout.write(f"\rProgress: {progress:.2f}%")
                        sys.stdout.flush()           
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.11"  or argv1=="step18.11s11":        
            dict_name_length={}
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/10_sum','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    sample,chromosome,strand,start,end,hi_entropy,hi_entropy_class,pos_start,pos_end,pos_length,one_type,percent=eachline_arr
                    length=int(pos_length)*float(percent)
                    name=f"{hi_entropy_class}|||{one_type}"
                    if name not in  dict_name_length:dict_name_length[name]= 0
                    dict_name_length[name]+=length
            with open('./samples_satellite/18_Homogenization/11_HOR2homo/11_sum','w') as f:  
                f.write(f"hi_entropy_class\tone_type\tlength\n")
                for name ,length in dict_name_length.items():
                    hi_entropy_class,one_type=name.split('|||')
                    f.write(f"{hi_entropy_class}\t{one_type}\t{length}\n")
            import pandas as pd
            # Read data
            df = pd.read_csv('./samples_satellite/18_Homogenization/11_HOR2homo/11_sum', sep='\t')
            
            # Convert to wide format: rows as one_type, columns as hi_entropy_class
            df_wide = df.pivot(index='one_type', columns='hi_entropy_class', values='length').fillna(0)
            
            # Reset index to make one_type a column
            df_wide = df_wide.reset_index()
            
            # Save
            df_wide.to_csv('./samples_satellite/18_Homogenization/11_HOR2homo/11_sum_format', sep='\t', index=False)
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.11"  or argv1=="step18.11s20":         
            import pandas as pd
            import numpy as np
            
            # Read file, set first column as index
            df1 = pd.read_csv('./samples_satellite/18_Homogenization/11_HOR2homo/4_plot_HORtype', sep='\t', index_col=0)
            df2 = pd.read_csv('./samples_satellite/18_Homogenization/11_HOR2homo/11_sum_format', sep='\t', index_col=0)
            
            # Ensure all numeric columns are numeric type, convert non-numeric to NaN, then fill with 0
            df1 = df1.apply(pd.to_numeric, errors='coerce').fillna(0)
            df2 = df2.apply(pd.to_numeric, errors='coerce').fillna(0)
            
            # Add, fill missing values with 0 (ensure index alignment)
            df_merged = df1.add(df2, fill_value=0)
            
            # Ensure no NaN values again
            df_merged = df_merged.fillna(0)
            
            # Save original merged data (keep)
            df_merged_reset = df_merged.reset_index()
            # Ensure no NaN when saving
            df_merged_reset = df_merged_reset.fillna(0)
            df_merged_reset.to_csv('./samples_satellite/18_Homogenization/11_HOR2homo/20_plot_HORtype_merged', sep='\t', index=False, float_format='%.0f')
            
            # Calculate row sums
            row_sums = df_merged.sum(axis=1)
            
            # Avoid division by zero; rows that are all zeros remain as 0
            df_percentage = df_merged.div(row_sums, axis=0).where(row_sums != 0, 0)
            
            # Ensure no NaN in percentage data again
            df_percentage = df_percentage.fillna(0)
            
            # Reset index
            df_percentage = df_percentage.reset_index()
            
            # Save percentage data
            df_percentage.to_csv('./samples_satellite/18_Homogenization/11_HOR2homo/20_plot_HORtype_merged_percentage', sep='\t', index=False, float_format='%.6f')
            
            # Print statistics
            print(f"Saved: 20_plot_HORtype_merged (original data)")
            print(f"Saved: 20_plot_HORtype_merged_percentage (percentage data, each row sums to 1)")
            print(f"\nData statistics:")
            print(f"df1 shape: {df1.shape}, missing value count: {df1.isna().sum().sum()}")
            print(f"df2 shape: {df2.shape}, missing value count: {df2.isna().sum().sum()}")
            print(f"df_merged shape: {df_merged.shape}, missing value count: {df_merged.isna().sum().sum()}")
            print(f"df_percentage shape: {df_percentage.shape}, missing value count: {df_percentage.isna().sum().sum()}")
            print(f"\nRow sum check (first 5 rows):")
            print(df_percentage.iloc[:, 1:].sum(axis=1).head())
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.11"  or argv1=="step18.11s20p":   
            print('Plotting')            
            R_txt=f'''library(pheatmap)

library(pheatmap)

# Read data
input_file1 <- read.table('20_plot_HORtype_merged_percentage', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL, sep = '\t')

# Set first column as row names
rownames(input_file1) <- input_file1[,1]
input_file1 <- input_file1[,-1]

# Specify row order
desired_order <- c("noHOR", "HOR", "(MON107)2", "(MON107)3", "(MON107)4", 
                   "(MON107)5", "(MON107)6", "(MON107)10", "(MON135)2", 
                   "(MON135)3", "MON51-MON107", "MON79-MON107", "MON107-MON135", 
                   "MON107-MON163", "MON51-(MON107)2", "(MON79)2-MON107", 
                   "MON79-(MON107)2", "MON79-MON107-MON135", 
                   "MON51-(MON79)2-MON107", "other")

# Rearrange in specified order
input_file1 <- input_file1[desired_order, ]

# Convert to matrix
data_matrix <- as.matrix(input_file1)

# Calculate log10 transformation
log_matrix <- log10(data_matrix + 0.01)

# Find positions where original data is exactly 0
zero_positions <- data_matrix == 0

# Set these positions to NA (so pheatmap will display them in gray or white)
log_matrix[zero_positions] <- NA

# Define colors (excluding white because NA will automatically display as white)
log_colors <- colorRampPalette(c("navy", "blue", "cyan", "yellow", "red", "darkred"))(100)

# Draw heatmap
pheatmap(log_matrix,
         main = "HOR Type Distribution (log10 scale)\n(White = 0%)",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = log_colors,
         display_numbers = FALSE,
         fontsize_row = 8,
         fontsize_col = 8,
         angle_col = 45,
         cellwidth = 15,
         cellheight = 10,
         na_col = "white",  # NA values displayed as white
         filename = "20_plot_HORtype_merged_percentage_log.pdf")
    
    '''
            with open(f'./samples_satellite/18_Homogenization/11_HOR2homo/plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)
            new_directory = f"./samples_satellite/18_Homogenization/11_HOR2homo/"
            os.chdir(new_directory)
            subprocess.run(['Rscript plot.R'], shell=True)    
            os.chdir('../../../')           
    
    if argv1=="stepall" or 'step18.12' in argv1:
        print('Calculate HOR repeat count and HI index')
        subprocess.run(["mkdir ./samples_satellite/18_Homogenization/12_HOR2homo"], shell=True)
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.12"  or argv1=="step18.12s00":
            print('Refer to step18.10, ./samples_satellite/18_Homogenization/10_plot/kmer_all_info is similar to this file')
    
            one_kmer=31
            with open(f'./samples_satellite/18_Homogenization/12_HOR2homo/kmer_info','w') as f2:
                f2.write(f"sample\tchromosome\tstrand\tstart\tend\thi_entropy\tymin\tymax\n")   ## Note: Hi_entropy
                with open(f'./samples_satellite/18_Homogenization/2_kmer/kmer_{one_kmer}','r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        chromosome,region,segment_file,total_count,dominant_count,unique_count,hi_dominant,hi_diversity,entropy,hi_entropy,Coverage,rare_kmer_ratio=eachline_arr
                        sample=region.split(':')[0]
                        if sample not in ["Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap1","V059.hap2","V063.hap1","V063.hap2","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap1","Chardonnay_hap2","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap1","V074.hap2","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap1","V062.hap2","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap1","PinotNoir_hap2","PinotNoir2_hap1","PinotNoir2_hap2","V069.hap1","V069.hap2","PN40024","PN40024_hap1","PN40024_hap2","V087.hap1","V087.hap2","V061.hap1","V061.hap2","V065.hap1","V065.hap2","V064.hap1","V064.hap2","ThompsonSeedless_hap1","ThompsonSeedless_hap2","ThompsonSeedless2_hap1","ThompsonSeedless2_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap1","V067.hap2","V058.hap1","V058.hap2","V107.hap1","V107.hap2","V108.hap1","V108.hap2","V081.hap1","V081.hap2","V015.hap1","V015.hap2","V105.hap1","V105.hap2","V126.hap1","V126.hap2","V112.hap1","V112.hap2","DavidiiGrape_hap1","DavidiiGrape_hap2","V007.hap1","V007.hap2","V008.hap1","V008.hap2","V012.hap1","V012.hap2","V100.hap1","V100.hap2","V099.hap1","V099.hap2","V098.hap1","V098.hap2","V117.hap1","V117.hap2","V120.hap1","V120.hap2","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap1","V124.hap2","V125.hap1","V125.hap2","WoollyGrape_hap1","WoollyGrape_hap2","V106.hap1","V106.hap2","V102.hap1","V102.hap2","V079.hap1","V079.hap2","V030.hap1","V030.hap2","V031.hap1","V031.hap2","V005.hap1","V005.hap2","V048.hap1","V048.hap2","V050.hap1","V050.hap2","V051.hap1","V051.hap2","V049.hap1","V049.hap2","V038.hap1","V038.hap2","V043.hap1","V043.hap2","V052.hap1","V052.hap2","V034.hap1","V034.hap2","V023.hap1","V023.hap2","V032.hap1","V032.hap2","V033.hap1","V033.hap2","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap1","V022.hap2","V040.hap1","V040.hap2","V019.hap1","V019.hap2","V041.hap1","V041.hap2","V077.hap1","V077.hap2","V037.hap1","V037.hap2","V096.hap1","V096.hap2","ShineMuscat_hap1","ShineMuscat_hap2","V036.hap1","V036.hap2","V076.hap1","V076.hap2","V055.hap1","V055.hap2","V018.hap1","V018.hap2","V072.hap1","V072.hap2"]:continue
                        strand=region.split(":")[1][-1]
                        pos_arr=segment_file.split('___')[1].split('-')
                        start1=pos_arr[0]
                        end1=pos_arr[1]
                        start=min(start1,end1)
                        end=max(start1,end1)
                        ymin='' ;ymax=''
                        if strand=='+': ymin=0 ;ymax=1
                        elif strand=='-': ymin=-1 ;ymax=0
                        f2.write(f"{sample}\t{chromosome}\t{strand}\t{start}\t{end}\t{hi_entropy}\t{ymin}\t{ymax}\n")
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.12"  or argv1=="step18.12s0":    
            
            input_list=[]
            with open('./samples_satellite/18_Homogenization/12_HOR2homo/kmer_info','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if len(eachline_arr)!=8:print(line)
                    sample,chromosome,strand,start,end,hi_entropy,ymin,ymax=eachline_arr
                    #if sample!="MuscatHamburg_hap1" or chromosome!=	"Chr19":continue
    
                    input_list.append([sample,chromosome,strand,start,end,hi_entropy])
            #print(input_list)        
            #sys.exit()
            region_name_list=[]
            dict_region_samplechr={}
            with open('./samples_satellite/2_good_regions','r') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if eachline_arr[0]=='sample':continue
                    if eachline_arr[0]=='':continue
                    chromosome=eachline_arr[1].replace('region_','Chr')
                    sample=eachline_arr[0]
                    if sample not in ["Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap1","V059.hap2","V063.hap1","V063.hap2","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap1","Chardonnay_hap2","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap1","V074.hap2","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap1","V062.hap2","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap1","PinotNoir_hap2","PinotNoir2_hap1","PinotNoir2_hap2","V069.hap1","V069.hap2","PN40024","PN40024_hap1","PN40024_hap2","V087.hap1","V087.hap2","V061.hap1","V061.hap2","V065.hap1","V065.hap2","V064.hap1","V064.hap2","ThompsonSeedless_hap1","ThompsonSeedless_hap2","ThompsonSeedless2_hap1","ThompsonSeedless2_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap1","V067.hap2","V058.hap1","V058.hap2","V107.hap1","V107.hap2","V108.hap1","V108.hap2","V081.hap1","V081.hap2","V015.hap1","V015.hap2","V105.hap1","V105.hap2","V126.hap1","V126.hap2","V112.hap1","V112.hap2","DavidiiGrape_hap1","DavidiiGrape_hap2","V007.hap1","V007.hap2","V008.hap1","V008.hap2","V012.hap1","V012.hap2","V100.hap1","V100.hap2","V099.hap1","V099.hap2","V098.hap1","V098.hap2","V117.hap1","V117.hap2","V120.hap1","V120.hap2","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap1","V124.hap2","V125.hap1","V125.hap2","WoollyGrape_hap1","WoollyGrape_hap2","V106.hap1","V106.hap2","V102.hap1","V102.hap2","V079.hap1","V079.hap2","V030.hap1","V030.hap2","V031.hap1","V031.hap2","V005.hap1","V005.hap2","V048.hap1","V048.hap2","V050.hap1","V050.hap2","V051.hap1","V051.hap2","V049.hap1","V049.hap2","V038.hap1","V038.hap2","V043.hap1","V043.hap2","V052.hap1","V052.hap2","V034.hap1","V034.hap2","V023.hap1","V023.hap2","V032.hap1","V032.hap2","V033.hap1","V033.hap2","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap1","V022.hap2","V040.hap1","V040.hap2","V019.hap1","V019.hap2","V041.hap1","V041.hap2","V077.hap1","V077.hap2","V037.hap1","V037.hap2","V096.hap1","V096.hap2","ShineMuscat_hap1","ShineMuscat_hap2","V036.hap1","V036.hap2","V076.hap1","V076.hap2","V055.hap1","V055.hap2","V018.hap1","V018.hap2","V072.hap1","V072.hap2"]:continue
                    one_region_name=f"{eachline_arr[0]}:{eachline_arr[1]}{eachline_arr[6]}:{eachline_arr[3]}-{eachline_arr[4]}"
                    #if one_region_name!='MuscatHamburg_hap1:region_19-:16767395-18114955':continue
                    region_name_list.append(one_region_name)
                    dict_region_samplechr[one_region_name]=[eachline_arr[0],chromosome]
                    
            ALL_LENGTH=0
            dict_samplechr_info={}     
            ALL_MARK_LEN=0
            for one_region in  region_name_list:    
                sample,chromosome=dict_region_samplechr[one_region]
                #if sample!="MuscatHamburg_hap1" or chromosome!=	"Chr19":continue
                samplechr=f"{sample}|{chromosome}"
                if samplechr not in dict_samplechr_info:
                    dict_samplechr_info[samplechr]=[]
                with open(f'./samples_satellite/12_Hor_stat/1_blocks_stat//{one_region}','r')   as f:
                    next(f)
                    for line in f:
                        #print(line)
                        eachline_arr=line.strip().split('\t')
                        mer,segment,segment_type,segment_type_num,segment_length,HOR_type,circ_start,circ_end,pos_start1,pos_end1,markserial_start,markserial_end,layer,father_layer,pos_length,markserial_num,HOR_repeat_num=eachline_arr
                        pos_start1=int(pos_start1)
                        pos_end1=int(pos_end1)
                        pos_start=min(pos_start1,pos_end1)
                        pos_end=max(pos_start1,pos_end1)
                        pos_length=int(pos_length)
                        #if float(HOR_repeat_num)==150:print(eachline_arr)
                        ALL_LENGTH+=pos_length
                        dict_samplechr_info[samplechr].append([mer,segment,segment_type,segment_length,HOR_type,pos_start,pos_end,pos_length,markserial_num,HOR_repeat_num])
                        #if segment=='28-23|28-28-23|28-28-28-23|28-28-23':
                        #    ALL_MARK_LEN+=pos_length
                        #    if pos_length>1000:print(samplechr,eachline_arr)
            print(f'ALL_LENGTH:{ALL_LENGTH}')
            #print(f'ALL_MARK_LEN:{ALL_MARK_LEN}')
            
            def run_step(one):    
                sample,chromosome,strand,start,end,hi_entropy=one
                #print(one)
                start=int(start)
                end=int(end)
                region_len=end-start+1
                samplechr=f"{sample}|{chromosome}"
                HOR_list=dict_samplechr_info[samplechr]
                good_HOR_list=[]
                for one_HOR in HOR_list:
                    #if segment=='28-23|28-28-23|28-28-28-23|28-28-23':print(one_HOR)
                    mer,segment,segment_type,segment_length,HOR_type,pos_start1,pos_end1,pos_length,markserial_num,HOR_repeat_num=one_HOR
                    #if float(HOR_repeat_num)==150:print(one_HOR)
                    pos_start=int(pos_start1)
                    pos_end=int(pos_end1)
                    #if pos_end-pos_start<0:print(11111111111)
                    one_HOR_len=pos_end-pos_start+1
                    if pos_start>=end or pos_end<=start:continue
                    elif pos_start>=start and pos_end<=end:
                        good_HOR_list.append([one,one_HOR,1]);#print(1)
                        
                    elif pos_start<=start and pos_end>end:  
                        percent=region_len/one_HOR_len
                        good_HOR_list.append([one,one_HOR,percent]);#print(1)
                    elif pos_start<=start and pos_end>=start and pos_end<=end :      
                        percent=(pos_end-start+1)/one_HOR_len
                        good_HOR_list.append([one,one_HOR,percent]);#print(1)
                    elif pos_start>=start and pos_start<=end and pos_end>=end :      
                        percent=(end-pos_start+1)/one_HOR_len     
                        good_HOR_list.append([one,one_HOR,percent]);#print(1)
                    else:
                        print('error, this situation should not occur')
                #print(HOR_list)        
                return good_HOR_list 
                    
            
            with open('./samples_satellite/18_Homogenization/12_HOR2homo/0_sum','w') as f:
                f.write(f"sample\tchromosome\tstrand\tstart\tend\thi_entropy\tmer\tsegment\tsegment_type\tsegment_length\tHOR_type\tpos_start\tpos_end\tpos_length\tmarkserial_num\tHOR_repeat_num\tpercent\n")
                with Pool(processes=thread) as pool:
                    # Use imap to get results one by one
                    for i, result in enumerate(pool.imap(run_step, input_list), start=1):
                        # Results can be processed here, e.g., stored or printed
                        for one_result in result:
                            HIregion,one_HOR,percent=one_result
                            ####
                            sample,chromosome,strand,Start,end,hi_entropy=HIregion
                            mer,segment,segment_type,segment_length,HOR_type,pos_start,pos_end,pos_length,markserial_num,HOR_repeat_num=one_HOR
                            #if segment=='28-23|28-28-23|28-28-28-23|28-28-23':print(one_HOR)
                            #####
                            f.write(f"{sample}\t{chromosome}\t{strand}\t{Start}\t{end}\t{hi_entropy}\t{mer}\t{segment}\t{segment_type}\t{segment_length}\t{HOR_type}\t{pos_start}\t{pos_end}\t{pos_length}\t{markserial_num}\t{HOR_repeat_num}\t{percent}\n")
                            
                        ###########    
                        progress = (i / len(input_list)) * 100
                        sys.stdout.write(f"\rProgress: {progress:.2f}%")
                        sys.stdout.flush()   
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.12"  or argv1=="step18.12s1":  
            def normalize_hor_pattern_rotation_only(sequence, separator='|'):
                """
                Normalize HOR monomer sequence, only handling rotation, not considering reverse
                
                Example: "28-28-28-23|28-28-23|28-28-28-28-23"
                   -> Normalize to unique string
                
                Strategy: Find the smallest lexicographic order among all rotations
                """
                # 1. Split by separator
                parts = sequence.split(separator)
                n = len(parts)
                
                # 2. Generate all rotations
                candidates = []
                doubled = parts + parts  # Double sequence for rotation
                
                for i in range(n):
                    rotated = doubled[i:i+n]
                    candidate = separator.join(rotated)
                    candidates.append(candidate)
                
                # 3. Return the smallest lexicographic string
                return min(candidates)    
                
            dict_HORtype_length={}
            dict_segmentlen_length={}
            dict_mer_length={}
            with open('./samples_satellite/18_Homogenization/12_HOR2homo/1_sum_filter','w') as f2:
                f2.write(f"horname\tsample\tchromosome\tstrand\tstart\tend\thi_entropy\tmer\tsegment\tsegment_type\tsegment_length\tHOR_type\tpos_start\tpos_end\tpos_length\tmarkserial_num\tHOR_repeat_num\tpercent\n")
                with open('./samples_satellite/18_Homogenization/12_HOR2homo/0_sum','r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t')            
                        sample,chromosome,strand,start,end,hi_entropy,mer,segment,segment_type,segment_length,HOR_type,pos_start,pos_end,pos_length,markserial_num,HOR_repeat_num,percent=eachline_arr
    
                        ##############
    
                        #### Count mer sum
                        dict_change = {
                                "28-23|28-28-28-23": "MON51-MON107",
                                "28-28-23|28-28-28-23": "MON79-MON107",
                                "28-28-28-23|28-28-28-23": "(MON107)2",
                                "28-28-28-23|28-28-28-28-23": "MON107-MON135",
                                "28-28-28-23|28-28-28-28-28-23": "MON107-MON163",
                                "28-28-28-28-23|28-28-28-28-23": "(MON135)2",
                                "28-23|28-28-28-23|28-28-28-23": "MON51-(MON107)2",           # After normalization
                                "28-28-23|28-28-23|28-28-28-23": "(MON79)2-MON107",      # After normalization
                                "28-28-23|28-28-28-23|28-28-28-23": "MON79-(MON107)2",      # After normalization (note difference from above)
                                "28-28-23|28-28-28-23|28-28-28-28-23": "MON79-MON107-MON135",
                                "28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)3",
                                "28-28-28-28-23|28-28-28-28-23|28-28-28-28-23": "(MON135)3",
                                "28-23|28-28-23|28-28-28-23|28-28-23": "MON51-(MON79)2-MON107",  # After normalization
                                "28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)4",
                                "28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)5",
                                "28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)6",
                                "28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)10"
                            }
                        if segment not in dict_change: name_new='other'
                        else: name_new=   dict_change[segment] 
                        #if segment=='28-23|28-28-23|28-28-28-23|28-28-23':print(eachline_arr)
                        #if name_new=='MON51-(MON79)2-MON107':print(eachline_arr)
                        f2.write(f"{name_new}\t{sample}\t{chromosome}\t{strand}\t{start}\t{end}\t{hi_entropy}\t{mer}\t{segment}\t{segment_type}\t{segment_length}\t{HOR_type}\t{pos_start}\t{pos_end}\t{pos_length}\t{markserial_num}\t{HOR_repeat_num}\t{percent}\n")
                        ###
                        repeat_num_good=float(HOR_repeat_num)*float(percent)
                        if mer not in dict_mer_length: dict_mer_length[mer]=[]
                        dict_mer_length[mer].append([repeat_num_good,hi_entropy])
                        if segment_length not in dict_segmentlen_length: dict_segmentlen_length[segment_length]=[]
                        dict_segmentlen_length[segment_length].append([repeat_num_good,hi_entropy])
                        if name_new not in dict_HORtype_length:dict_HORtype_length[name_new]=[]
                        dict_HORtype_length[name_new].append([repeat_num_good,hi_entropy])
            def get_entropy_range(hi_entropy):
                """
                Generate interval range based on hi_entropy value
                Example: 0.125 -> "0.12-0.13"
                     0.111 -> "0.11-0.12"
                     0.119 -> "0.11-0.12"
                     0.120 -> "0.12-0.13"
                """
                hi_entropy = float(hi_entropy)
                # Round down to 0.01
                lower = math.floor(hi_entropy * 100) / 100
                upper = lower + 0.01
                return f"{lower:.2f}-{upper:.2f}"       
         
            with open('./samples_satellite/18_Homogenization/12_HOR2homo/1_plot_mer','w') as f:
                f.write('mer\tlength\tHIindex\tHIentropy\tHIentropy_class\n')
                for  mer,length_HIentropy in  dict_mer_length.items():     
                    for one in length_HIentropy:
                        repeat_num_good,HIentropy=one
                        HIentropy_class=get_entropy_range(HIentropy)
                        f.write(f"{mer}\t{repeat_num_good}\t{HIentropy}\t{HIentropy_class}\n")
            with open('./samples_satellite/18_Homogenization/12_HOR2homo/1_plot_segmentlen','w') as f:
                f.write('segmentlen\tlength\tHIentropy\tHIentropy_class\n')
                for  segmentlen,length_HIentropy in  dict_segmentlen_length.items():          
                    for one in length_HIentropy:
                        repeat_num_good,HIentropy=one
                        HIentropy_class=get_entropy_range(HIentropy)
                        f.write(f"{segmentlen}\t{repeat_num_good}\t{HIentropy}\t{HIentropy_class}\n")
            with open('./samples_satellite/18_Homogenization/12_HOR2homo/1_plot_HORtype','w') as f:
                f.write('HORtype\tlength\tHIentropy\tHIentropy_class\n')
                for  HORtype,length_HIentropy in  dict_HORtype_length.items():          
                    for one in length_HIentropy:
                        repeat_num_good,HIentropy=one
                        HIentropy_class=get_entropy_range(HIentropy)
                        f.write(f"{HORtype}\t{repeat_num_good}\t{HIentropy}\t{HIentropy_class}\n")              
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.12"  or argv1=="step18.12s2":     
            # Read and output by HIentropy_class
            dict_class_len={}
            with open('./samples_satellite/18_Homogenization/12_HOR2homo/1_plot_mer', 'r') as f:
                header = f.readline()
                for line in f:
                    parts = line.strip().split('\t')
                    mer, repeat_num_good,  HIentropy, HIentropy_class = parts
                    newclass=f"{HIentropy_class}|||||{mer}"
                    if newclass not in dict_class_len:
                        dict_class_len[newclass]=0
                    dict_class_len[newclass]   +=float(repeat_num_good)
            with open('./samples_satellite/18_Homogenization/12_HOR2homo/2_plot_mer', 'w') as f:     
                f.write('mer\tHIentropy_class\trepeat_num_good\n')
                for one,repeat_num_good in dict_class_len.items():
                    HIentropy_class,mer=one.split('|||||')
                    f.write(f"{mer}\t{HIentropy_class}\t{repeat_num_good:.1f}\n")
            ###############
            dict_class_len={}
            with open('./samples_satellite/18_Homogenization/12_HOR2homo/1_plot_segmentlen', 'r') as f:
                header = f.readline()
                for line in f:
                    parts = line.strip().split('\t')
                    segmentlen, repeat_num_good,  HIentropy, HIentropy_class = parts
                    newclass=f"{HIentropy_class}|||||{segmentlen}"
                    if newclass not in dict_class_len:
                        dict_class_len[newclass]=0
                    dict_class_len[newclass]   +=float(repeat_num_good)
            with open('./samples_satellite/18_Homogenization/12_HOR2homo/2_plot_segmentlen', 'w') as f:     
                f.write('segmentlen\tHIentropy_class\trepeat_num_good\n')
                for one,repeat_num_good in dict_class_len.items():
                    HIentropy_class,segmentlen=one.split('|||||')
                    f.write(f"{segmentlen}\t{HIentropy_class}\t{repeat_num_good:.1f}\n")        
            ###############
            dict_class_len={}
            with open('./samples_satellite/18_Homogenization/12_HOR2homo/1_plot_HORtype', 'r') as f:
                header = f.readline()
                for line in f:
                    parts = line.strip().split('\t')
                    HORtype, repeat_num_good,  HIentropy, HIentropy_class = parts
                    newclass=f"{HIentropy_class}|||||{HORtype}"
                    if newclass not in dict_class_len:
                        dict_class_len[newclass]=0
                    dict_class_len[newclass]   +=float(repeat_num_good)
            with open('./samples_satellite/18_Homogenization/12_HOR2homo/2_plot_HORtype', 'w') as f:     
                f.write('HORtype\tHIentropy_class\trepeat_num_good\n')
                for one,repeat_num_good in dict_class_len.items():
                    HIentropy_class,HORtype=one.split('|||||')
                    f.write(f"{HORtype}\t{HIentropy_class}\t{repeat_num_good:.1f}\n")        
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.12"  or argv1=="step18.12s3":            
            print('Plotting')            
            R_txt=f'''library(ggplot2)
    
    # Set working directory
    setwd('./')
    
    # Read data
    input_file1 <- read.table('2_plot_mer', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')
    
    # Create plot
    p <- ggplot()
    p <- p + geom_point(data = input_file1, aes(x = mer, y = repeat_num_good,color=HIentropy_class))
    p <- p + theme_classic()
    # Add title and axis labels
    p <- p + theme(
      #axis.title = element_blank(),  # Hide axis titles
      #axis.text = element_blank(),   # Hide axis tick text
      #legend.position = "none",
      #axis.ticks = element_blank(),  # Hide axis ticks
      #axis.line = element_blank()    # Hide axis lines
    )+ coord_flip()
    pdf("3_plot_mer.pdf", width = 20 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()
    '''
            with open(f'./samples_satellite/18_Homogenization/12_HOR2homo/plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)
            new_directory = f"./samples_satellite/18_Homogenization/12_HOR2homo/"
            os.chdir(new_directory)
            subprocess.run(['Rscript plot.R'], shell=True)   
            
            print('Plotting')            
            R_txt=f'''library(ggplot2)
    
    # Set working directory
    setwd('./')
    
    # Read data
    input_file1 <- read.table('2_plot_HORtype', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')
    
    # Create plot
    p <- ggplot()
    p <- p + geom_point(data = input_file1, aes(x = HORtype, y = repeat_num_good,color=HIentropy_class))
    p <- p + theme_classic()
    # Add title and axis labels
    p <- p + theme(
      #axis.title = element_blank(),  # Hide axis titles
      #axis.text = element_blank(),   # Hide axis tick text
      #legend.position = "none",
      #axis.ticks = element_blank(),  # Hide axis ticks
      #axis.line = element_blank()    # Hide axis lines
    )+ coord_flip()
    pdf("3_plot_HORtype.pdf", width = 20 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()
    '''
            with open(f'plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)
            subprocess.run(['Rscript plot.R'], shell=True)           
            
            
            print('Plotting')            
            R_txt=f'''library(ggplot2)
    
    # Set working directory
    setwd('./')
    
    # Read data
    input_file1 <- read.table('2_plot_segmentlen', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')
    
    # Create plot
    p <- ggplot()
    p <- p + geom_point(data = input_file1, aes(x = segmentlen, y = repeat_num_good,color=HIentropy_class))
    p <- p + theme_classic()
    # Add title and axis labels
    p <- p + theme(
      #axis.title = element_blank(),  # Hide axis titles
      #axis.text = element_blank(),   # Hide axis tick text
      #legend.position = "none",
      #axis.ticks = element_blank(),  # Hide axis ticks
      #axis.line = element_blank()    # Hide axis lines
    )+ coord_flip()
    pdf("3_plot_segmentlen.pdf", width = 20 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()
    '''
            with open(f'plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)
            subprocess.run(['Rscript plot.R'], shell=True)           
                    
            
            
            
            
            
            
            
            os.chdir('../../../')     
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.12"  or argv1=="step18.12s4": 
            print('Convert to wide format')
            dict_HORtype_HIClass={}
            HORtype_set=set()
            HIentropy_class_set=set()
            with open('./samples_satellite/18_Homogenization/12_HOR2homo/2_plot_HORtype','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    HORtype,HIentropy_class,repeat_num_good=eachline_arr
                    dict_HORtype_HIClass[f"{HORtype}||{HIentropy_class}"]=repeat_num_good
                    HORtype_set.add(HORtype)
                    HIentropy_class_set.add(HIentropy_class)
            HORtype_list=list(HORtype_set);  HORtype_list.sort()
            HIentropy_class_list=list(HIentropy_class_set);  HIentropy_class_list.sort()
            print(HORtype_list,HIentropy_class_list)
            with open('./samples_satellite/18_Homogenization/12_HOR2homo/4_plot_HORtype','w') as f:
                head='one_type\t'+'\t'.join(HIentropy_class_list)+'\n'
                f.write(head)
                for one_HORtype in HORtype_list:
                    line_arr=[]
                    for one_HIclass in HIentropy_class_list:
                        name = f"{one_HORtype}||{one_HIclass}"
                        if name not in dict_HORtype_HIClass: repeat_num_good=0
                        else :repeat_num_good=dict_HORtype_HIClass[name]
                        line_arr.append(repeat_num_good)
                    line_str=    one_HORtype+'\t'+'\t'.join(str(x) for x in line_arr)+'\n'
                    f.write(line_str)
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.12" or argv1=="step18.12s4b":         
            import pandas as pd
            import numpy as np
            
            # Read original data
            df = pd.read_csv('./samples_satellite/18_Homogenization/12_HOR2homo/4_plot_HORtype', sep='\t', index_col=0)
            
            # Ensure all numeric columns are numeric type
            df = df.apply(pd.to_numeric, errors='coerce').fillna(0)

            
            # Calculate row sums, convert to percentage (each row sums to 100)
            row_sums = df.sum(axis=1)
            df_percentage = df.div(row_sums, axis=0).where(row_sums != 0, 0) * 100
            
            # Save percentage data
            df_percentage.reset_index().to_csv('./samples_satellite/18_Homogenization/12_HOR2homo/4_plot_HORtype_percentage', sep='\t', index=False, float_format='%.2f')
            
            print(f"Saved: 20_plot_HORtype_merged (original data)")
            print(f"Saved: 20_plot_HORtype_merged_percentage (percentage data, each row sums to 100%)")                    
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.12"  or argv1=="step18.12s10":
            print('Re-analyze, HOR and non-HOR, by HOR density')
            def get_entropy_range(hi_entropy):
                """
                Generate interval range based on hi_entropy value
                Example: 0.125 -> "0.12-0.13"
                     0.111 -> "0.11-0.12"
                     0.119 -> "0.11-0.12"
                     0.120 -> "0.12-0.13"
                """
                hi_entropy = float(hi_entropy)
                # Round down to 0.01
                lower = math.floor(hi_entropy * 100) / 100
                upper = lower + 0.01
                return f"{lower:.2f}-{upper:.2f}"
                
            dict_pos_info={}
            with open('./samples_satellite/18_Homogenization/12_HOR2homo/0_sum','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    sample,chromosome,strand,start,end,hi_entropy,mer,segment,segment_type,segment_length,HOR_type,pos_start,pos_end,pos_length,markserial_num,HOR_repeat_num,percent=eachline_arr
                    length=int(end)-int(start)+1
                    density=float(HOR_repeat_num)*float(percent)/(length/10000)
                    #hi_entropy_class=get_entropy_range(hi_entropy)
                    pos=f"{sample}||{chromosome}||{strand}||{start}||{end}"
                    if pos not in dict_pos_info:
                        dict_pos_info[pos]={}
                        dict_pos_info[pos]['HOR_repeat_num']=0
                        dict_pos_info[pos]['hi_entropy']=float(hi_entropy)
                    dict_pos_info[pos]['HOR_repeat_num']+=density
            with open('./samples_satellite/18_Homogenization/12_HOR2homo/10_sum','w') as f2:   
                f2.write(f"pos\tHOR_repeat_num\thi_entropy\n")
                for pos,info in dict_pos_info.items():
                    HOR_repeat_num=info['HOR_repeat_num']
                    hi_entropy=info['hi_entropy']    
      
                    f2.write(f"{pos}\t{HOR_repeat_num}\t{hi_entropy}\n")
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.12"  or argv1=="step18.12s10p1":   
            print('Plotting')            
            R_txt=f'''library(ggplot2)
    
    
    # Set working directory
    setwd('./')
    
    library(pheatmap)
    
    # Read data
    input_file1 <- read.table('10_sum', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL, sep = '\t')
    
   # Create plot	
    p <- ggplot()
    p <- p + geom_point(data = input_file1, aes(x = HOR_repeat_num, y = hi_entropy,color=hi_entropy))
    p <- p + theme_classic()
    # Add title and axis labels
    p <- p + theme(
      #axis.title = element_blank(),  # Hide axis titles
      #axis.text = element_blank(),   # Hide axis tick text
      #legend.position = "none",
      #axis.ticks = element_blank(),  # Hide axis ticks
      #axis.line = element_blank()    # Hide axis lines
    )
    pdf("10_sum.pdf", width = 20 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()
    
    '''
            with open(f'./samples_satellite/18_Homogenization/12_HOR2homo/plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)
            new_directory = f"./samples_satellite/18_Homogenization/12_HOR2homo/"
            os.chdir(new_directory)
            subprocess.run(['Rscript plot.R'], shell=True)    
            os.chdir('../../../')                       
        if argv1=="stepall" or argv1=="step18" or argv1=="step18.12"  or argv1=="step18.12s10p2":   
            print('Plotting')            
            R_txt=f'''
    library(ggplot2)
library(hexbin)

# Set working directory
setwd('./')

# Read data
input_file1 <- read.table('10_sum', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL, sep = '\t')

# Hexagonal binning
p <- ggplot(input_file1, aes(x = HOR_repeat_num, y = hi_entropy)) +
  geom_hex(bins = 30) +  # bins controls the number of hexagons
  scale_fill_gradient(low = "white", high = "#006699", trans = "log") +  # log transformation for better appearance
  theme_classic() +
  labs(x = "HOR Repeat Number", y = "HI Entropy", fill = "Count")

pdf("10_sum_hexbin.pdf", width = 8, height = 6)
print(p)
dev.off()
    
    '''
            with open(f'./samples_satellite/18_Homogenization/12_HOR2homo/plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)
            new_directory = f"./samples_satellite/18_Homogenization/12_HOR2homo/"
            os.chdir(new_directory)
            subprocess.run(['Rscript plot.R'], shell=True)    
            os.chdir('../../../')
                    
        
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))      
time_end=timeit.default_timer()
print('All the running time: %.0f Seconds'%(time_end-time_start))






















