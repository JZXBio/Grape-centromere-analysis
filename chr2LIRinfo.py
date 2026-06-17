#!/usr/bin/python
# -*- coding: UTF-8 -*-
import sys
argvs = sys.argv
argv_len = len(sys.argv)
if argv_len < 2 or argvs[1] in ("h", "-h", "help", "-help"):
    print("chr2LIRinfo.py-----help:")
    print("")
    print("Usage：")
    print("chr2LIRinfo.py step1   Count various types of VSat1-inter LIR, some LIR seem to be missing")
    print("chr2LIRinfo.py step1b   Directly count chromosome ratios")
    print("chr2LIRinfo.py step1s   Summarize step_LIR_stat results, convert to wide format")
    
    print("chr2LIRinfo.py step_LIR_stat_time   Calculate similarity between LIR arms, can infer formation time, correct based on multiple files, first look at connecting lines, if no connecting lines take the center point (e.g., Chr18 basically has no connecting lines), then manually check low similarity ones, found that above 93% are mostly correct, below 93% drop sharply to around 85, all wrong... The code uses two layers of with Pool(processes=thread) as pool")
    print("chr2LIRinfo.py step_LIR_stat_time_sum   Summarize, need to manually filter results, delete unreasonable rows, filename: LIR_stat1s_revised_manual")
    print("chr2LIRinfo.py step_LIR_stat_time_plot   Re-plot")
    
    print("chr2LIRinfo.py step3_LIRVSat6_TE_diff_revise   Compare VSat6 and LIRSat1")
    
    #print("chr2blast.py step_LIR_stat3    Reclassify, containing LISat3 is ClassIII, only LISat1 is ClassI, only LISat2 is ClassII, with LISat2/1 is ClassII")
    
    print('The final match_percent may be greater than 1 because the alignment length can be well controlled, but there are overlapping alignments, causing some matches to be counted twice, so this value is high, considered to have limited impact, for reference only')
    print("-thread \t\tThreads (default 70), some steps use multiprocessing")
    print("-i      \t\tStep0 requires, input fasta file")
    print(" ")
    sys.exit()

import subprocess
import csv
import os
import math
import time
import timeit
import copy
import re  # For regular expressions
from multiprocessing import Pool, cpu_count
time_start = timeit.default_timer()
print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))

argv1 = argvs[1]

if os.path.exists('./chr2LIRinfo') == False:
    subprocess.run(["mkdir chr2LIRinfo"], shell=True)

if argv1 == 'step1':
    subprocess.run(["mkdir chr2LIRinfo/1_stat"], shell=True)
    dict_samplechr_infolist = {}
    with open('./chr2blast/1_blastn/sum', 'r') as f:
        next(f)
        for line in f:
            eachline_arr = line.strip().split('\t')
            sample = eachline_arr[0]
            cenid = eachline_arr[1]
            if cenid not in ['LIRSat1', 'LIRSat2', 'LIRSat3']:
                continue
            chromosome = eachline_arr[2]
            strand = eachline_arr[3]
            start = int(eachline_arr[4])
            end = int(eachline_arr[5])
            length = int(eachline_arr[6])
            one_id = chromosome
            
            #
            one_id_lower = one_id.lower()
            if one_id.isdigit():
                one_id = "Chr" + one_id
            elif one_id_lower.endswith('chr01') or one_id_lower.endswith('chr1') or one_id_lower == 'region_1':
                one_id = "Chr1"
            elif one_id_lower.endswith('chr02') or one_id_lower.endswith('chr2') or one_id_lower == 'region_2':
                one_id = "Chr2"
            elif one_id_lower.endswith('chr03') or one_id_lower.endswith('chr3') or one_id_lower == 'region_3':
                one_id = "Chr3"
            elif one_id_lower.endswith('chr04') or one_id_lower.endswith('chr4') or one_id_lower == 'region_4':
                one_id = "Chr4"
            elif one_id_lower.endswith('chr05') or one_id_lower.endswith('chr5') or one_id_lower == 'region_5':
                one_id = "Chr5"
            elif one_id_lower.endswith('chr06') or one_id_lower.endswith('chr6') or one_id_lower == 'region_6':
                one_id = "Chr6"
            elif one_id_lower.endswith('chr07') or one_id_lower.endswith('chr7') or one_id_lower == 'region_7':
                one_id = "Chr7"
            elif one_id_lower.endswith('chr08') or one_id_lower.endswith('chr8') or one_id_lower == 'region_8':
                one_id = "Chr8"
            elif one_id_lower.endswith('chr09') or one_id_lower.endswith('chr9') or one_id_lower == 'region_9':
                one_id = "Chr9"
            elif one_id_lower.endswith('chr10') or one_id_lower == 'region_10':
                one_id = "Chr10"
            elif one_id_lower.endswith('chr11') or one_id_lower == 'region_11':
                one_id = "Chr11"
            elif one_id_lower.endswith('chr12') or one_id_lower == 'region_12':
                one_id = "Chr12"
            elif one_id_lower.endswith('chr13') or one_id_lower == 'region_13':
                one_id = "Chr13"
            elif one_id_lower.endswith('chr14') or one_id_lower == 'region_14':
                one_id = "Chr14"
            elif one_id_lower.endswith('chr15') or one_id_lower == 'region_15':
                one_id = "Chr15"
            elif one_id_lower.endswith('chr16') or one_id_lower == 'region_16':
                one_id = "Chr16"
            elif one_id_lower.endswith('chr17') or one_id_lower == 'region_17':
                one_id = "Chr17"
            elif one_id_lower.endswith('chr18') or one_id_lower == 'region_18':
                one_id = "Chr18"
            elif one_id_lower.endswith('chr19') or one_id_lower == 'region_19':
                one_id = "Chr19"
            elif one_id_lower.endswith('chr20') or one_id_lower == 'region_20':
                one_id = "Chr20"
            else:
                continue
            chromosome = one_id
            samplechr = sample + '___' + chromosome
            if samplechr not in dict_samplechr_infolist:
                dict_samplechr_infolist[samplechr] = []
            dict_samplechr_infolist[samplechr].append([cenid, start, end, length])
    #print(dict_samplechr_infolist)
    print('Loading sample info table ../samples_satellite/sample_info')
    dict_sample_type = {}
    with open('../samples_satellite/sample_info') as f:
        for line in f:
            eachline_arr = line.strip().split('\t')
            if len(eachline_arr) != 6:
                continue
            sample_name = eachline_arr[0]
            sample_type = eachline_arr[3]
            dict_sample_type[sample_name] = sample_type
    manual_list = ['Chr3', 'Chr13']
    manual_list = []
    with open('./chr2LIRinfo/1_stat/LIR_stat1', 'w') as f2:
        f2.write(f"sample\tsample_type\tchromosome\tborder1\tborder2\tonetype\tnum\n")
        with open('../samples_satellite/2_good_regions_interarray_LIR', 'r') as f:
            next(f)
            for line in f:
                eachline_arr = line.strip().split('\t')
                sample = eachline_arr[0]
                #
                sample_type = dict_sample_type[sample]
                if sample_type not in ["Eurasian", "East_Asia", "America"]:
                    continue
                if sample in ['074.hap1', '074.hap2']:
                    continue
                #
                chromosome = eachline_arr[1]
                if chromosome in manual_list:
                    continue
                if chromosome == 'Chr18':
                    continue
                LIR_interarrays = eachline_arr[13]
                if LIR_interarrays == 'NA':
                    continue
                LIR_interarrays_list = LIR_interarrays.split('|')
                for one in LIR_interarrays_list:
                    arr = one.split('-')
                    if len(arr) != 2:
                        print(arr)
                        print('error,238')
                    start, end = arr
                    #start,end=int(start),int(end)
                    border1 = int(start) - 50000
                    border2 = int(end) + 50000
                    dict_type_num = {}
                    ############
                    samplechr = sample + '___' + chromosome
                    if samplechr not in dict_samplechr_infolist:
                        continue
                    info_list = dict_samplechr_infolist[samplechr]
                    for one_info in info_list:
                        one_info_type, one_info_start, one_info_end, one_info_length = one_info
                        if one_info_start > border1 and one_info_end < border2:
                            if one_info_type not in dict_type_num:
                                dict_type_num[one_info_type] = 0
                            dict_type_num[one_info_type] += one_info_length
                    #print(dict_type_num)
                    for onetype, num in dict_type_num.items():
                        f2.write(f"{sample}\t{sample_type}\t{chromosome}\t{border1}\t{border2}\t{onetype}\t{num}\n")
    with open('./chr2LIRinfo/1_stat/LIR_stat1', 'a') as f2:
        with open('../samples_satellite/chr18_LIR', 'r') as f:
            next(f)
            for line in f:
                eachline_arr = line.strip().split('\t')
                if len(eachline_arr) != 4:
                    continue
                sample, start, end, one_type = eachline_arr
                #
                sample_type = dict_sample_type[sample]
                if sample_type not in ["Eurasian", "East_Asia", "America"]:
                    continue
                if sample in ['074.hap1', '074.hap2']:
                    continue
                #
                if one_type == '.':
                    continue
                start, end = float(start), float(end)
                start = 10000000 + start * 1000000
                end = 10000000 + end * 1000000
                chromosome = 'Chr18'
                border1 = int(start) - 50000
                border2 = int(end) + 50000
                dict_type_num = {}
                ############
                samplechr = sample + '___' + chromosome
                if samplechr not in dict_samplechr_infolist:
                    continue
                info_list = dict_samplechr_infolist[samplechr]
                for one_info in info_list:
                    one_info_type, one_info_start, one_info_end, one_info_length = one_info
                    if one_info_start > border1 and one_info_end < border2:
                        if one_info_type not in dict_type_num:
                            dict_type_num[one_info_type] = 0
                        dict_type_num[one_info_type] += one_info_length
                #print(dict_type_num)
                for onetype, num in dict_type_num.items():
                    f2.write(f"{sample}\t{sample_type}\t{chromosome}\t{border1}\t{border2}\t{onetype}\t{num}\n")
    for one_manual_chr in manual_list:
        with open('./chr2LIRinfo/1_stat/LIR_stat1', 'a') as f2:
            with open(f'../samples_satellite/{one_manual_chr}_LIR', 'r') as f:
                next(f)
                for line in f:
                    eachline_arr = line.strip().split('\t')
                    if len(eachline_arr) != 4:
                        continue
                    sample, start, end, one_type = eachline_arr
                    #
                    sample_type = dict_sample_type[sample]
                    if sample_type not in ["Eurasian", "East_Asia", "America"]:
                        continue
                    if sample in ['074.hap1', '074.hap2']:
                        continue
                    #
                    start, end = float(start), float(end)
                    #start=10000000+start*1000000
                    #end=10000000+end*1000000
                    chromosome = one_manual_chr
                    #if one_type=='.':continue
                    border1 = int(start) - 50000
                    border2 = int(end) + 50000
                    dict_type_num = {}
                    ############
                    samplechr = sample + '___' + chromosome
                    if samplechr not in dict_samplechr_infolist:
                        continue
                    info_list = dict_samplechr_infolist[samplechr]
                    for one_info in info_list:
                        one_info_type, one_info_start, one_info_end, one_info_length = one_info
                        print(sample, one_info, border1, border2)
                        if one_info_start > border1 and one_info_end < border2:
                            if one_info_type not in dict_type_num:
                                dict_type_num[one_info_type] = 0
                            dict_type_num[one_info_type] += one_info_length
                    #print(dict_type_num)
                    for onetype, num in dict_type_num.items():
                        f2.write(f"{sample}\t{sample_type}\t{chromosome}\t{border1}\t{border2}\t{onetype}\t{num}\n")
    dict_LIR_type_num = {}
    with open('./chr2LIRinfo/1_stat/LIR_stat1', 'r') as f:
        next(f)
        for line in f:
            eachline = line.strip()
            eachline_arr = eachline.split('\t')
            sample, sample_type, chromosome, border1, border2, onetype, num = eachline_arr
            LIR = sample + '___' + chromosome + '___' + border1 + '___' + border2
            if LIR not in dict_LIR_type_num:
                dict_LIR_type_num[LIR] = {}
            dict_LIR_type_num[LIR][onetype] = num

    dict_sampletype_chromosome_mark_num = {}
    with open('./chr2LIRinfo/1_stat/LIR_stat1b', 'w') as f:
        f.write(f"sample\tsample_type\tchromosome\tborder1\tborder2\tLIRSat1\tLIRSat2\tLIRSat3\tmark\n")
        for LIR, dict_type_num in dict_LIR_type_num.items():
            sample, chromosome, border1, border2 = LIR.split('___')
            if 'LIRSat1' in dict_type_num:
                LIRSat1 = dict_type_num['LIRSat1']
            else:
                LIRSat1 = 0
            if 'LIRSat2' in dict_type_num:
                LIRSat2 = dict_type_num['LIRSat2']
            else:
                LIRSat2 = 0
            if 'LIRSat3' in dict_type_num:
                LIRSat3 = dict_type_num['LIRSat3']
            else:
                LIRSat3 = 0
            LIRSat1 = int(LIRSat1)
            LIRSat2 = int(LIRSat2)
            LIRSat3 = int(LIRSat3)

            if LIRSat1 > 5000 and LIRSat2 < 5000 and LIRSat3 < 5000:
                mark = 'LIRSat1-main'
            elif LIRSat1 < 5000 and LIRSat2 > 5000 and LIRSat3 < 5000:
                mark = 'LIRSat2-main'
            elif LIRSat1 < 5000 and LIRSat2 < 5000 and LIRSat3 > 5000:
                mark = 'LIRSat3-main'
            elif LIRSat1 > 5000 and LIRSat2 > 5000 and LIRSat3 < 5000:
                mark = 'LIRSat_12'
            elif LIRSat1 > 5000 and LIRSat2 < 5000 and LIRSat3 > 5000:
                mark = 'LIRSat_13'
            elif LIRSat1 < 5000 and LIRSat2 > 5000 and LIRSat3 > 5000:
                mark = 'LIRSat_23'
            elif LIRSat1 > 5000 and LIRSat2 > 5000 and LIRSat3 > 5000:
                mark = 'LIRSat_123'
            else:
                mark = 'little_LIRSat'
            """
            if LIRSat1>5000 and LIRSat2<5000 and LIRSat3<5000: mark='LIRSat1-main'
            elif LIRSat1<5000 and LIRSat2>5000 and LIRSat3<5000: mark='LIRSat2-main'
            elif LIRSat3>5000: mark='LIRSat3-main'
            elif LIRSat1>5000 and LIRSat2>5000 and LIRSat3<5000: mark='LIRSat_12'
            else:mark='little_LIRSat'
            """
            sample_type = dict_sample_type[sample]
            f.write(f"{sample}\t{sample_type}\t{chromosome}\t{border1}\t{border2}\t{LIRSat1}\t{LIRSat2}\t{LIRSat3}\t{mark}\n")
            ###
            if sample_type not in dict_sampletype_chromosome_mark_num:
                dict_sampletype_chromosome_mark_num[sample_type] = {}
            if chromosome not in dict_sampletype_chromosome_mark_num[sample_type]:
                dict_sampletype_chromosome_mark_num[sample_type][chromosome] = {}
            if mark not in dict_sampletype_chromosome_mark_num[sample_type][chromosome]:
                dict_sampletype_chromosome_mark_num[sample_type][chromosome][mark] = 0
            dict_sampletype_chromosome_mark_num[sample_type][chromosome][mark] += 1

    with open('./chr2LIRinfo/1_stat/LIR_stat1c', 'w') as f:
        f.write(f"sampletype\tchromosome\tsampletype_mark\tmark\tnum\n")
        for sampletype, dict_chromosome_mark_num in dict_sampletype_chromosome_mark_num.items():
            for chromosome, dict_mark_num in dict_chromosome_mark_num.items():
                for mark, num in dict_mark_num.items():
                    if mark == 'little_LIRSat':
                        continue
                    if sampletype not in ['America', 'East_Asia', 'Eurasian']:
                        continue
                    sampletype_mark = f"{sampletype}___{mark}"
                    f.write(f"{sampletype}\t{chromosome}\t{sampletype_mark}\t{mark}\t{num}\n")

    R_txt = f'''
            library(ggplot2)
            library(dplyr)
            
            print("")
           
              input_file=read.table('LIR_stat1c', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
            input_file$chromosome <- factor(
              input_file$chromosome,
              levels = c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10','Chr11','Chr12','Chr13','Chr14','Chr15','Chr16','Chr17','Chr18','Chr19')
            )
  
    color_values=c(
        "LIRSat1-main"="#94caec",
        "LIRSat2-main"="#5ea899",
        "LIRSat3-main"="#dccd7c",
        "LIRSat_12"="#2151b2",
        "LIRSat_13"="#a53c3b",
        "LIRSat_23"="#adad42",
        "LIRSat_123"="#ff0000",
        "no_LIRSat"="#cecece"
    )        
            # Create plot object
            p <- ggplot() 
         #p = p +  geom_bar(data = input_file, aes(x = chromosome, y = num, fill = mark),stat = "identity", position = "stack", width = 0.8)   # Stacked bar chart
         # Calculate percentage (group by chromosome)
input_file_percent <- input_file %>%
  group_by(chromosome) %>%
  mutate(percent = num / sum(num) * 100) %>%
  ungroup()

# Draw percentage stacked bar chart
p <- ggplot() +
  geom_bar(
    data = input_file_percent,
    aes(x = chromosome, y = percent, fill =mark),  #,group=sampletype
    stat = "identity",
    position = "stack",  # Default is "stack", can be omitted
    width = 0.8
  )+scale_fill_manual(values = color_values, drop = FALSE)
                #+coord_fixed()+
                  #p = p + facet_grid(sampletype ~ ., scales = "free", space = "free")
            p =p+theme_classic() +         
                theme(
                  #axis.ticks.y = element_blank(),
                  #axis.text.y = element_blank(),
                  #legend.position = "none",
                  #axis.text.x = element_blank()
                ) 
             
              
              # Save as PDF
              pdf(file = paste0('LIR_stat1a', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
              print(p)
              dev.off()
   input_file$sampletype <- factor(
      input_file$sampletype,
      levels = c("Eurasian","East_Asia","America")
    )    
# Draw grouped bar chart (facet by chromosome)
p <- ggplot(input_file) +
  geom_bar(
    aes(x = sampletype, y = num, fill = mark),
    stat = "identity",  # Key fix: explicitly specify stat = "identity"
    width = 0.8,
    position = "stack"  # Default stack, change to "dodge" for side-by-side
  ) +scale_fill_manual(values = color_values, drop = FALSE)+
  facet_grid(. ~ chromosome) +  # Facet by chromosome
  theme_classic() +
  theme(
    # Optional: adjust facet title and label styles
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0),
    # Optional: adjust legend position
    legend.position = "right"
  ) +
  labs(x = "Sample Type", y = "Count", fill = "Mark")  # Add axis labels and legend title

# Save as PDF
pdf(file = "LIR_stat1b.pdf", width = 20 / 2.54, height = 20 / 2.54)
print(p)
dev.off()
    '''
    with open(f'./chr2LIRinfo/1_stat/plot.R', 'w', encoding='utf-8') as f:
        f.write(R_txt)
    new_directory = f'./chr2LIRinfo/1_stat/'
    os.chdir(new_directory)
    subprocess.run(['Rscript plot.R'], shell=True)
    os.chdir('../../../')

if argv1 == 'step1b':
    subprocess.run(["mkdir chr2LIRinfo/1_stat_b_chr"], shell=True)
    dict_samplechr_type_num = {}
    with open('./chr2blast/1_blastn/sum', 'r') as f:
        next(f)
        for line in f:
            eachline_arr = line.strip().split('\t')
            sample = eachline_arr[0]
            cenid = eachline_arr[1]
            #if cenid not in ['LIRSat1','LIRSat2','LIRSat3']:continue
            chromosome = eachline_arr[2]
            strand = eachline_arr[3]
            start = int(eachline_arr[4])
            end = int(eachline_arr[5])
            length = int(eachline_arr[6])
            one_id = chromosome
            
            #
            one_id_lower = one_id.lower()
            if one_id.isdigit():
                one_id = "Chr" + one_id
            elif one_id_lower.endswith('chr01') or one_id_lower.endswith('chr1') or one_id_lower == 'region_1':
                one_id = "Chr1"
            elif one_id_lower.endswith('chr02') or one_id_lower.endswith('chr2') or one_id_lower == 'region_2':
                one_id = "Chr2"
            elif one_id_lower.endswith('chr03') or one_id_lower.endswith('chr3') or one_id_lower == 'region_3':
                one_id = "Chr3"
            elif one_id_lower.endswith('chr04') or one_id_lower.endswith('chr4') or one_id_lower == 'region_4':
                one_id = "Chr4"
            elif one_id_lower.endswith('chr05') or one_id_lower.endswith('chr5') or one_id_lower == 'region_5':
                one_id = "Chr5"
            elif one_id_lower.endswith('chr06') or one_id_lower.endswith('chr6') or one_id_lower == 'region_6':
                one_id = "Chr6"
            elif one_id_lower.endswith('chr07') or one_id_lower.endswith('chr7') or one_id_lower == 'region_7':
                one_id = "Chr7"
            elif one_id_lower.endswith('chr08') or one_id_lower.endswith('chr8') or one_id_lower == 'region_8':
                one_id = "Chr8"
            elif one_id_lower.endswith('chr09') or one_id_lower.endswith('chr9') or one_id_lower == 'region_9':
                one_id = "Chr9"
            elif one_id_lower.endswith('chr10') or one_id_lower == 'region_10':
                one_id = "Chr10"
            elif one_id_lower.endswith('chr11') or one_id_lower == 'region_11':
                one_id = "Chr11"
            elif one_id_lower.endswith('chr12') or one_id_lower == 'region_12':
                one_id = "Chr12"
            elif one_id_lower.endswith('chr13') or one_id_lower == 'region_13':
                one_id = "Chr13"
            elif one_id_lower.endswith('chr14') or one_id_lower == 'region_14':
                one_id = "Chr14"
            elif one_id_lower.endswith('chr15') or one_id_lower == 'region_15':
                one_id = "Chr15"
            elif one_id_lower.endswith('chr16') or one_id_lower == 'region_16':
                one_id = "Chr16"
            elif one_id_lower.endswith('chr17') or one_id_lower == 'region_17':
                one_id = "Chr17"
            elif one_id_lower.endswith('chr18') or one_id_lower == 'region_18':
                one_id = "Chr18"
            elif one_id_lower.endswith('chr19') or one_id_lower == 'region_19':
                one_id = "Chr19"
            elif one_id_lower.endswith('chr20') or one_id_lower == 'region_20':
                one_id = "Chr20"
            else:
                continue
            chromosome = one_id
            samplechr = sample + '___' + chromosome
            if samplechr not in dict_samplechr_type_num:
                dict_samplechr_type_num[samplechr] = {}
            if cenid not in dict_samplechr_type_num[samplechr]:
                dict_samplechr_type_num[samplechr][cenid] = 0
            dict_samplechr_type_num[samplechr][cenid] += int(length)

    print('Loading sample info table ../samples_satellite/sample_info')
    dict_sample_type = {}
    with open('../samples_satellite/sample_info') as f:
        for line in f:
            eachline_arr = line.strip().split('\t')
            if len(eachline_arr) != 6:
                continue
            sample_name = eachline_arr[0]
            sample_type = eachline_arr[3]
            dict_sample_type[sample_name] = sample_type

    with open('./chr2LIRinfo/1_stat_b_chr/LIR_stat2b', 'w') as f2:
        f2.write(f"sample\tsample_type\tchromosome\tmark\n")
        dict_chromosome_mark_type_info = {}
        for one_samplechr, dict_type_num in dict_samplechr_type_num.items():
            sample, chromosome = one_samplechr.split('___')
            sample_type = dict_sample_type[sample]
            if sample_type not in ["Eurasian", "East_Asia", "America"]:
                continue
            if sample in ['074.hap1', '074.hap2']:
                continue
            if 'LIRSat1' not in dict_type_num:
                LIRSat1 = 0
            else:
                LIRSat1 = dict_type_num['LIRSat1']
            if 'LIRSat2' not in dict_type_num:
                LIRSat2 = 0
            else:
                LIRSat2 = dict_type_num['LIRSat2']
            if 'LIRSat3' not in dict_type_num:
                LIRSat3 = 0
            else:
                LIRSat3 = dict_type_num['LIRSat3']

            if LIRSat1 > 5000 and LIRSat2 < 5000 and LIRSat3 < 5000:
                mark = 'LIRSat1-main'
            elif LIRSat1 < 5000 and LIRSat2 > 5000 and LIRSat3 < 5000:
                mark = 'LIRSat2-main'
            elif LIRSat1 < 5000 and LIRSat2 < 5000 and LIRSat3 > 5000:
                mark = 'LIRSat3-main'
            elif LIRSat1 > 5000 and LIRSat2 > 5000 and LIRSat3 < 5000:
                mark = 'LIRSat_12'
            elif LIRSat1 > 5000 and LIRSat2 < 5000 and LIRSat3 > 5000:
                mark = 'LIRSat_13'
            elif LIRSat1 < 5000 and LIRSat2 > 5000 and LIRSat3 > 5000:
                mark = 'LIRSat_23'
            elif LIRSat1 > 5000 and LIRSat2 > 5000 and LIRSat3 > 5000:
                mark = 'LIRSat_123'
            elif LIRSat1 < 5000 and LIRSat2 < 5000 and LIRSat3 < 5000:
                mark = 'no_LIRSat'
            else:
                print('error')
                sys.exit()
            """
            if LIRSat1>5000 and LIRSat2<5000 and LIRSat3<5000: mark='LIRSat1-main'
            elif LIRSat1<5000 and LIRSat2>5000 and LIRSat3<5000: mark='LIRSat2-main'
            elif LIRSat3>5000: mark='LIRSat3-main'
            elif LIRSat1>5000 and LIRSat2>5000 and LIRSat3<5000: mark='LIRSat_12'
            elif LIRSat1<5000 and LIRSat2<5000 and LIRSat3<5000:mark='no_LIRSat'
            else:print('error');sys.exit()"""
            if chromosome not in dict_chromosome_mark_type_info:
                dict_chromosome_mark_type_info[chromosome] = {}
            if mark not in dict_chromosome_mark_type_info[chromosome]:
                dict_chromosome_mark_type_info[chromosome][mark] = {}
            if sample_type not in dict_chromosome_mark_type_info[chromosome][mark]:
                dict_chromosome_mark_type_info[chromosome][mark][sample_type] = 0
            dict_chromosome_mark_type_info[chromosome][mark][sample_type] += 1
            f2.write(f"{sample}\t{sample_type}\t{chromosome}\t{mark}\n")

    with open('./chr2LIRinfo/1_stat_b_chr/LIR_stat2', 'w') as f:
        f.write(f"sampletype\tchromosome\tsampletype_mark\tmark\tnum\n")
        for chromosome, dict_mark_type_num in dict_chromosome_mark_type_info.items():
            for mark, dict_type_num in dict_mark_type_num.items():
                for sampletype, num in dict_type_num.items():
                    #if mark =='little_LIRSat':continue
                    #if sampletype not in ['America','East_Asia','Eurasian']:continue
                    sampletype_mark = f"{sampletype}___{mark}"
                    f.write(f"{sampletype}\t{chromosome}\t{sampletype_mark}\t{mark}\t{num}\n")
    R_txt = f'''
            library(ggplot2)
            library(dplyr)
            
            print("")
           
              input_file=read.table('LIR_stat2', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
            input_file$chromosome <- factor(
              input_file$chromosome,
              levels = c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10','Chr11','Chr12','Chr13','Chr14','Chr15','Chr16','Chr17','Chr18','Chr19')
            )
            input_file$sampletype <- factor(
              input_file$sampletype,
              levels = c("Eurasian","East_Asia","America")
            )            
            # Create plot object
            p <- ggplot() 
         #p = p +  geom_bar(data = input_file, aes(x = chromosome, y = num, fill = mark),stat = "identity", position = "stack", width = 0.8)   # Stacked bar chart
         # Calculate percentage (group by chromosome)
input_file_percent <- input_file %>%
  group_by(chromosome) %>%
  mutate(percent = num / sum(num) * 100) %>%
  ungroup()
color_values=c(
    "LIRSat1-main"="#94caec",
    "LIRSat2-main"="#5ea899",
    "LIRSat3-main"="#dccd7c",
    "LIRSat_12"="#2151b2",
    "LIRSat_13"="#a53c3b",
    "LIRSat_23"="#adad42",
    "LIRSat_123"="#ff0000",
    "no_LIRSat"="#cecece"
)

# Draw percentage stacked bar chart
p <- ggplot() +
  geom_bar(
    data = input_file_percent,
    aes(x = chromosome, y = percent, fill =mark),  #,group=sampletype
    stat = "identity",
    position = "stack",  # Default is "stack", can be omitted
    width = 0.8
  )+scale_fill_manual(values = color_values, drop = FALSE)#+coord_fixed()+
                  #p = p + facet_grid(sampletype ~ ., scales = "free", space = "free")
            p =p+theme_classic() +         
                theme(
                  #axis.ticks.y = element_blank(),
                  #axis.text.y = element_blank(),
                  #legend.position = "none",
                  #axis.text.x = element_blank()
                ) 
             
              
              # Save as PDF
              pdf(file = paste0('LIR_stat2a', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
              print(p)
              dev.off()
              
# Draw grouped bar chart (facet by chromosome)
p <- ggplot(input_file) +
  geom_bar(
    aes(x = sampletype, y = num, fill = mark),
    stat = "identity",  # Key fix: explicitly specify stat = "identity"
    width = 0.8,
    position = "stack"  # Default stack, change to "dodge" for side-by-side
  )+ scale_fill_manual(values = color_values, drop = FALSE)+
  facet_grid(. ~ chromosome) +  # Facet by chromosome
  theme_classic() +
  theme(
    # Optional: adjust facet title and label styles
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0),
    # Optional: adjust legend position
    legend.position = "right"
  ) +
  labs(x = "Sample Type", y = "Count", fill = "Mark")  # Add axis labels and legend title

# Save as PDF
pdf(file = "LIR_stat2b.pdf", width = 20 / 2.54, height = 20 / 2.54)
print(p)
dev.off()
    '''
    with open(f'./chr2LIRinfo/1_stat_b_chr/plot.R', 'w', encoding='utf-8') as f:
        f.write(R_txt)
    new_directory = f'./chr2LIRinfo/1_stat_b_chr//'
    os.chdir(new_directory)
    subprocess.run(['Rscript plot.R'], shell=True)
    os.chdir('../../../')

if argv1 == 'step1s':
    dict_type_num = {}
    with open('./chr2LIRinfo/1_stat/LIR_stat1', 'r') as f:
        next(f)
        for line in f:
            eachline_arr = line.strip().split('\t')
            sample, sample_type, chromosome, border1, border2, onetype, num = eachline_arr
            name = sample + '|||' + sample_type + '|||' + chromosome + '|||' + border1 + '|||' + border2
            if name not in dict_type_num:
                dict_type_num[name] = {}
            if onetype not in dict_type_num[name]:
                dict_type_num[name][onetype] = 0
            dict_type_num[name][onetype] += int(num)
    kk = 0
    with open('./chr2LIRinfo/1_stat/LIR_stat1s', 'w') as f:
        f.write(f"ID\tsample\tsample_type\tchromosome\tborder1\tborder2\tLIRSat1\tLIRSat2\tLIRSat3\tmark\tLIRclass\n")
        for name, dict_num in dict_type_num.items():
            sample, sample_type, chromosome, border1, border2 = name.split('|||')
            LIRSat1 = dict_num.get('LIRSat1', 0)
            LIRSat2 = dict_num.get('LIRSat2', 0)
            LIRSat3 = dict_num.get('LIRSat3', 0)
            if LIRSat1 > 5000 and LIRSat2 < 5000 and LIRSat3 < 5000:
                mark = 'LIRSat1-main'
                LIRclass = 'ClassI'
            elif LIRSat1 < 5000 and LIRSat2 > 5000 and LIRSat3 < 5000:
                mark = 'LIRSat2-main'
                LIRclass = 'ClassII'
            elif LIRSat1 < 5000 and LIRSat2 < 5000 and LIRSat3 > 5000:
                mark = 'LIRSat3-main'
                LIRclass = 'ClassIII'
            elif LIRSat1 > 5000 and LIRSat2 > 5000 and LIRSat3 < 5000:
                mark = 'LIRSat_12'
                LIRclass = 'ClassIV'
            elif LIRSat1 > 5000 and LIRSat2 < 5000 and LIRSat3 > 5000:
                mark = 'LIRSat_13'
                LIRclass = 'ClassV'
            elif LIRSat1 < 5000 and LIRSat2 > 5000 and LIRSat3 > 5000:
                mark = 'LIRSat_23'
                LIRclass = 'OTHER'
            elif LIRSat1 > 5000 and LIRSat2 > 5000 and LIRSat3 > 5000:
                mark = 'LIRSat_123'
                LIRclass = 'OTHER'
            else:
                mark = 'little_LIRSat'
                LIRclass = 'OTHER'
            kk += 1
            f.write(f"{kk}\t{sample}\t{sample_type}\t{chromosome}\t{border1}\t{border2}\t{LIRSat1}\t{LIRSat2}\t{LIRSat3}\t{mark}\t{LIRclass}\n")

if argv1 == 'step2_time':
    print('LIR from bioinformatics analysis, find connecting lines results, then analyze center point positions, then chr2aligntime analyze similarity')
    input_list = []
    with open('./chr2LIRinfo/1_stat/LIR_stat1s', 'r') as f:
        next(f)
        for line in f:
            eachline_arr = line.strip().split('\t')
            ID, sample, sample_type, chromosome, border1, border2, LIRSat1, LIRSat2, LIRSat3, mark, LIRclass = eachline_arr
            input_list.append(eachline_arr)
    input_list_num = len(input_list)
    
    ##
    plotlines = []
    with open("./chr2moddotplot3/2_slashout_updown", "r") as f:
        next(f)
        for line in f:
            eachline_arr = line.strip().split('\t')
            sample, chromosome, onetype, start_x, start_y, end_x, end_y, seq1_len, seq2_len = eachline_arr[:9]
            if onetype == 'down':
                if abs(int(end_x) - int(end_y)) > 500000:
                    continue
                plotlines.append([sample, chromosome, int(start_x), int(start_y), int(end_x), int(end_y)])
    ##
    #kk=0
    #for one in input_list:
    def run_step1(one):
        #kk+=1
        #print(f"{one}\t\t\t{kk}/{input_list_num}")
        print(one)
        ID, sample, sample_type, chromosome, border1, border2, LIRSat1, LIRSat2, LIRSat3, mark, LIRclass = one
        border1 = int(border1)
        border2 = int(border2)
        #center1=round(int(border1)+(int(border2)-int(border1))/2-5000,0)
        #center2=round(int(border1)+(int(border2)-int(border1))/2+5000,0)
        ##### Calculate center point
        tmp_lenth = 0
        tmp_line = []
        for oneline in plotlines:
            #if oneline[0]=='Baimunage_hap1' :print(oneline)
            if oneline[0] != sample or oneline[1] != chromosome:
                continue
            start_x, start_y, end_x, end_y = oneline[-4:]
            #print(oneline)
            if end_x > border1 and end_x < border2 and end_y > border1 and end_y < border2:
                if end_x - start_x > tmp_lenth:
                    tmp_lenth = end_x - start_x
                    tmp_line = [start_x, start_y, end_x, end_y]
        if tmp_line != []:
            start_x, start_y, end_x, end_y = tmp_line
            center1 = end_x
            center2 = end_y
            border1_revise = max(border1, start_x)
            border2_revise = min(border2, start_y)
            ###############################
            #print(center1,center2)
            para = f"{sample}:{chromosome}:{border1_revise}-{center1:.0f} {sample}:{chromosome}:{center2:.0f}-{border2_revise} "
            name = f"{ID}___{sample}:{chromosome}:{border1_revise}-{border2_revise}"
            if os.path.exists(f'./chr2LIRinfo/alignmtime/{name}/ALL_SUM') == True:
                return False
            subprocess.run([f'mkdir -p ./chr2LIRinfo/alignmtime/{name}'], shell=True)
            with open(f'./chr2LIRinfo/alignmtime/{name}/region_info', 'w') as ff:
                ff.write(f"ID\tsample\tchromosome\tborder1_revise\tcenter1_revise\tcenter2_revise\tborder2_revise\n")
                ff.write(f"{ID}\t{sample}\t{chromosome}\t{border1_revise}\t{center1:.0f}\t{center2:.0f}\t{border2_revise}\n")
            subprocess.run([f'python /home/lain/aaa_data/run0/new_work_dir/chr2aligntime.py {para} ./chr2LIRinfo/alignmtime/{name}'], shell=True)
        else:
            print('Missing connecting lines');
            center1 = round(int(border1) + (int(border2) - int(border1)) / 2 - 5000, 0)
            center2 = round(int(border1) + (int(border2) - int(border1)) / 2 + 5000, 0)
            #return False
            para = f"{sample}:{chromosome}:{border1}-{center1:.0f} {sample}:{chromosome}:{center2:.0f}-{border2} "
            name = f"{ID}___{sample}:{chromosome}:{border1}-{border2}"
            if os.path.exists(f'./chr2LIRinfo/alignmtime/{name}/ALL_SUM') == True:
                return False
            subprocess.run([f'mkdir -p ./chr2LIRinfo/alignmtime/{name}'], shell=True)
            with open(f'./chr2LIRinfo/alignmtime/{name}/region_info', 'w') as ff:
                ff.write(f"ID\tsample\tchromosome\tborder1_revise\tcenter1_revise\tcenter2_revise\tborder2_revise\n")
                ff.write(f"{ID}\t{sample}\t{chromosome}\t.\t{center1:.0f}\t{center2:.0f}\t.\n")
            subprocess.run([f'python /home/lain/aaa_data/run0/new_work_dir/chr2aligntime.py {para} ./chr2LIRinfo/alignmtime/{name}'], shell=True)

    with Pool(processes=5) as pool:
        # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step1, input_list), start=1):
            # Process results here, e.g., store or print
            progress = (i / len(input_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush()

if argv1 == 'step2_time_sum':
    files = os.listdir('./chr2LIRinfo/alignmtime/')
    files.sort()
    subprocess.run([f'mkdir ./chr2LIRinfo/2_FRtime'], shell=True)
    dict_id_info = {}
    for one in files:
        ID = one.split('___')[0]
        #print(ID)
        #dict_id_info[ID]={}
        if os.path.exists(f'./chr2LIRinfo/alignmtime/{one}/ALL_SUM') == False:
            print(one)
            continue
        if os.path.exists(f'./chr2LIRinfo/alignmtime/{one}/region_info') == False:
            print(one)
            continue
        with open(f'./chr2LIRinfo/alignmtime/{one}/ALL_SUM', 'r') as f:
            length, var, percent = f.read().strip().split('\t')
        with open(f'./chr2LIRinfo/alignmtime/{one}/region_info', 'r') as f:
            next(f)
            for line in f:
                eachline_arr = line.strip().split('\t')
                ID, sample, chromosome, border1_revise, center1_revise, center2_revise, border2_revise = eachline_arr
        dict_id_info[ID] = [sample, chromosome, border1_revise, center1_revise, center2_revise, border2_revise, length, var, percent]
    #print(dict_id_info)

    with open('./chr2LIRinfo/2_FRtime/LIR_stat1s', 'w') as f2:
        f2.write(f"ID\tsample\tsample_type\tchromosome\tborder1\tborder2\tLIRSat1\tLIRSat2\tLIRSat3\tmark\tLIRclass\tborder1_revise\tcenter1_revise\tcenter2_revise\tborder2_revise\tlength\tvar\tpercent\n")
    with open('./chr2LIRinfo/1_stat/LIR_stat1s', 'r') as f:
        next(f)
        for line in f:
            eachline = line.strip()
            eachline_arr = eachline.split('\t')
            ID, sample, sample_type, chromosome, border1, border2, LIRSat1, LIRSat2, LIRSat3, mark, LIRclass = eachline_arr
            #print(ID)
            if ID in dict_id_info:
                _, _, border1_revise, center1_revise, center2_revise, border2_revise, length, var, percent = dict_id_info[ID]
                with open('./chr2LIRinfo/2_FRtime/LIR_stat1s', 'a') as f2:
                    f2.write(f"{eachline}\t{border1_revise}\t{center1_revise}\t{center2_revise}\t{border2_revise}\t{length}\t{var}\t{percent}\n")

    R_txt = f'''
                library(ggplot2)
                library(dplyr)
                
                print("")
               
      input_file=read.table('LIR_stat1s', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file, aes(x = LIRclass, y = percent),position = position_jitter(width = 0.3),alpha=0.95,  shape = 16,stroke = 0)
p <- p + geom_boxplot(data = input_file, aes(x = LIRclass, y = percent),size=1,outlier.shape = NA)

    # Save as PDF
    pdf(file = "1_LIRClass_stat.pdf", width = 20 / 2.54, height = 20 / 2.54)
    print(p)
    dev.off()
        '''
    with open(f'./chr2LIRinfo/2_FRtime//plot.R', 'w', encoding='utf-8') as f:
        f.write(R_txt)
    new_directory = f'./chr2LIRinfo/2_FRtime///'
    os.chdir(new_directory)
    subprocess.run(['Rscript plot.R'], shell=True)
    os.chdir('../../../')
if argv1 == 'step2_time_plot':
    R_txt = f'''
                library(ggplot2)
                library(dplyr)
                
                print("")
               
      input_file=read.table('LIR_stat1s_revised_manual', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
      input_file$chromosome <- factor(
          input_file$chromosome,
          levels = c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10','Chr11','Chr12','Chr13','Chr14','Chr15','Chr16','Chr17','Chr18','Chr19')
        )  
# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file, aes(x = LIRclass, y = percent),position = position_jitter(width = 0.3),alpha=0.95,  shape = 16,stroke = 0)
p <- p + geom_boxplot(data = input_file, aes(x = LIRclass, y = percent),size=1,outlier.shape = NA)

    # Save as PDF
    pdf(file = "2_LIR_Class_stat.pdf", width = 20 / 2.54, height = 20 / 2.54)
    print(p)
    dev.off()
    
# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file, aes(x = chromosome, y = percent),position = position_jitter(width = 0.3),alpha=0.95,  shape = 16,stroke = 0)
p <- p + geom_boxplot(data = input_file, aes(x = chromosome, y = percent),size=1,outlier.shape = NA)

    # Save as PDF
    pdf(file = "2_LIR_Chr_stat.pdf", width = 20 / 2.54, height = 20 / 2.54)
    print(p)
    dev.off() 
    
    
# Create plot
p <- ggplot()
#p <- p + geom_boxplot(data = input_file, aes(x = chromosome, y = percent),size=1,outlier.shape = NA)
p <- p + geom_point(data = input_file, aes(x = chromosome, y = percent,color=LIRclass),position = position_jitter(width = 0.3),alpha=0.95,  shape = 16,stroke = 0)
p=p+theme_classic()
    # Save as PDF
    pdf(file = "2_LIR_ChrClass_stat.pdf", width = 10 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()  
    
# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file, 
                    aes(x = percent, y = chromosome, color = LIRclass),
                    position = position_jitter(height = 0.3),  # Note: height instead of width
                    alpha = 0.95, shape = 16, stroke = 0, size = 1)
p <- p + theme_classic()+
scale_color_manual(values = c(
    "ClassI"="#65beea",
    "ClassII"="#5ea899",
    "ClassIII"="#db762a","ClassV"="#db762a",
   "ClassIV"="#2151b2"), drop = FALSE)
p <- p + scale_y_discrete(limits = rev(levels(input_file$chromosome)))  # Keep Chr1 at the top

    # Save as PDF
    pdf(file = "2_LIR_ChrClass_stat2.pdf", width = 10 / 2.54, height = 8 / 2.54)
    print(p)
    dev.off()    
        '''
    with open(f'./chr2LIRinfo/2_FRtime//plot.R', 'w', encoding='utf-8') as f:
        f.write(R_txt)
    new_directory = f'./chr2LIRinfo/2_FRtime///'
    os.chdir(new_directory)
    subprocess.run(['Rscript plot.R'], shell=True)
    os.chdir('../../../')

if argv1 == 'step2_time_plot_tekayinsert':
    print('Find TE insertion time inside LIR, only analyze tekay.')
    subprocess.run(["mkdir -p chr2LIRinfo/2_FRtime_TE/tekay"], shell=True)
    TE_list = []
    with open('/home/lain/aaa_data/run0/new_work_dir/chr2EDTA/5_tree/Tekay/0_info_LTR', 'r') as f:
        next(f)
        for line in f:
            eachline_arr = line.strip().split('\t')
            sample, sample_type, chromosome, start, end, strand, TE_name, Class, Class1, Class2, Class3, Complete, Domains, dom_names, seqs, seqs_used, mafftinput, seq_index, LTR1_start, LTR1_end, LTR1_len, seq1, LTR1_TSD, LTR2_start, LTR2_end, LTR2_len, seq2, LTR2_TSD, similarity, time_years_Ma = eachline_arr
            TE_list.append([sample, chromosome, int(start), int(end), similarity])
    with open('./chr2LIRinfo/2_FRtime_TE/tekay/5_LIR_TE', 'w') as f2:
        f2.write(f"sample\tchromosome\tLIRclass\tLIR_left\tLIR_right\tLIR_similarity\tTE_start\tTE_end\tTE_similarity\n")
        with open('./chr2LIRinfo/2_FRtime/LIR_stat1s_revised_manual', 'r') as f:
            next(f)
            for line in f:
                eachline = line.strip()
                eachline_arr = eachline.split('\t')

                if len(eachline_arr) != 18:
                    continue

                ID, sample, sample_type, chromosome, border1, border2, LIRSat1, LIRSat2, LIRSat3, mark, LIRclass, border1_revise, center1_revise, center2_revise, border2_revise, length, var, percent = eachline_arr
                if border1_revise != '.':
                    border_left = int(border1_revise) - 50000
                else:
                    border_left = int(border1)
                if border2_revise != '.':
                    border_right = int(border2_revise) + 50000
                else:
                    border_right = int(border2)
                ##
                for one_TE in TE_list:
                    TE_sample, TE_chromosome, TE_start, TE_end, TE_similarity = one_TE
                    if TE_sample == sample and TE_chromosome == chromosome:
                        if TE_start > border_left and TE_end < border_right:
                            percent_revise = float(percent) * 100
                            #if float(TE_similarity)<90:continue
                            f2.write(f"{sample}\t{chromosome}\t{LIRclass}\t{border_left}\t{border_right}\t{percent_revise}\t{TE_start}\t{TE_end}\t{TE_similarity}\n")

    R_txt = f'''
                library(ggplot2)
                library(dplyr)
                
                print("")
               
      input_file=read.table('5_LIR_TE', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
 
# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file, aes(x = LIR_similarity, y = TE_similarity),position = position_jitter(width = 0.3),alpha=0.95,  shape = 16,stroke = 0)

    # Save as PDF
    pdf(file = "5_LIRtime_TEtime.pdf", width = 10 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()
    
# Create plot
p <- ggplot()
p <- p + geom_boxplot(data = input_file, aes(x = LIRclass, y = TE_similarity),size=1,outlier.shape = NA)
p <- p + geom_point(data = input_file, aes(x = LIRclass, y = TE_similarity),position = position_jitter(width = 0.3),alpha=0.95,  shape = 16,stroke = 0)


    # Save as PDF
    pdf(file = "5_LIRClass_TEtime.pdf", width = 10 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()   
        '''
    with open(f'./chr2LIRinfo/2_FRtime_TE/tekay/plot.R', 'w', encoding='utf-8') as f:
        f.write(R_txt)
    new_directory = f'./chr2LIRinfo/2_FRtime_TE/tekay/'
    os.chdir(new_directory)
    subprocess.run(['Rscript plot.R'], shell=True)
    os.chdir('../../../')
    print('5_LIRtime_TEtime.pdf shows no clear trend, new TEs are concentrated in new LIRs which is fine, but some TE LTR similarities are too low')
    print('5_LIRClass_TEtime.pdf also has mediocre results, overall indicating TE insertion time estimates are not very accurate')

if argv1 == 'step3_LIRVSat6_TE_diff_revise':
    subprocess.run(["mkdir -p chr2LIRinfo/3_LIRVSat6_TE_diff"], shell=True)
    print('Calculate TE density differences between VSat6 and LIR, plus chromosome arms')
    
    def merge_intervals(intervals):
        """Merge overlapping intervals"""
        if not intervals:
            return []
        intervals.sort(key=lambda x: x[0])
        merged = [list(intervals[0])]
        for current_start, current_end in intervals[1:]:
            last_start, last_end = merged[-1]
            if current_start <= last_end:
                merged[-1][1] = max(last_end, current_end)
            else:
                merged.append([current_start, current_end])
        return merged
    
    def subtract_intervals(main_interval, exclude_intervals):
        """Subtract exclude intervals from main interval, return remaining intervals"""
        if not exclude_intervals:
            return [main_interval]
        
        start, end = main_interval
        # Only keep exclude intervals that overlap with main interval, and sort
        valid_excludes = []
        for ex_start, ex_end in exclude_intervals:
            if ex_start < end and ex_end > start:  # Overlap exists
                valid_excludes.append([max(start, ex_start), min(end, ex_end)])
        
        if not valid_excludes:
            return [[start, end]]
        
        # Sort by start position
        valid_excludes.sort(key=lambda x: x[0])
        
        # Merge overlapping exclude intervals
        merged_excludes = []
        cur_start, cur_end = valid_excludes[0]
        for ex_start, ex_end in valid_excludes[1:]:
            if ex_start <= cur_end:
                cur_end = max(cur_end, ex_end)
            else:
                merged_excludes.append([cur_start, cur_end])
                cur_start, cur_end = ex_start, ex_end
        merged_excludes.append([cur_start, cur_end])
        
        # Calculate remaining intervals
        remaining = []
        current_pos = start
        for ex_start, ex_end in merged_excludes:
            if current_pos < ex_start:
                remaining.append([current_pos, ex_start])
            current_pos = max(current_pos, ex_end)
        
        if current_pos < end:
            remaining.append([current_pos, end])
        
        return remaining
    
    def calculate_TE_coverage(TE_list, region_start, region_end, exclude_intervals=None):
        """Calculate actual TE coverage length within region (after merging overlaps), can exclude specific intervals"""
        overlapping_intervals = []
        for te_type, te_start, te_end in TE_list:
            # Calculate overlap between TE and target region
            overlap_start = max(te_start, region_start)
            overlap_end = min(te_end, region_end)
            if overlap_start < overlap_end:
                # If there are intervals to exclude
                if exclude_intervals:
                    # Subtract exclude intervals from TE overlap portion
                    remaining_intervals = subtract_intervals([overlap_start, overlap_end], exclude_intervals)
                    for start, end in remaining_intervals:
                        overlapping_intervals.append([start, end])
                else:
                    overlapping_intervals.append([overlap_start, overlap_end])
        
        # Merge overlapping intervals
        merged = merge_intervals(overlapping_intervals)
        
        # Calculate total coverage length and TE count
        total_covered_length = sum(end - start for start, end in merged)
        te_count = len(merged)
        
        return te_count, total_covered_length

    print('First record allowed range, 1M around VSat1')
    dict_samplechr_region = {}
    with open('/home/lain/aaa_data/run0/samples_satellite/2_good_regions_main', 'r') as f:
        next(f)
        for line in f:
            sample, chromosome, start, end, length = line.strip().split('\t')
            samplechr = f"{sample}|||{chromosome}"
            dict_samplechr_region[samplechr] = [int(start) - 1000000, int(end) + 1000000]
    
    print('VSat6 within this range will be used')
    target_VSat6_list = []
    USED_SAMPLECHR = set()
    with open('/home/lain/aaa_data/run0/stat_plot/0-region2info', 'r') as f:
        next(f)
        for line in f:
            sample, chromosome, chromosome_new, centype, chr_start, chr_end, length, strand, match_percent = line.strip().split('\t')
            samplechr = f"{sample}|||{chromosome_new}"
            if samplechr not in dict_samplechr_region:
                continue
            border1, border2 = dict_samplechr_region[samplechr]
            if centype == 'cen_383':
                if int(chr_start) > border1 and int(chr_end) < border2:
                    if samplechr in USED_SAMPLECHR:
                        continue
                    USED_SAMPLECHR.add(samplechr)
                    target_VSat6_list.append(['VSat6', sample, chromosome_new, int(chr_start), int(chr_end)])
    
    print('Store positions of VSat1/6/5')
    dict_samplechr_VSat1poslist = {}
    with open('/home/lain/aaa_data/run0/stat_plot/0-region2info', 'r') as f:
        next(f)
        for line in f:
            sample, chromosome, chromosome_new, centype, chr_start, chr_end, length, strand, match_percent = line.strip().split('\t')
            samplechr = f"{sample}|||{chromosome_new}"
            if samplechr not in dict_samplechr_region:
                continue
            if samplechr not in dict_samplechr_VSat1poslist:
                dict_samplechr_VSat1poslist[samplechr] = []
            #if centype not in ['cen_107','cen_355','cen_383']:continue
            if centype not in ['cen_107']:
                continue
            dict_samplechr_VSat1poslist[samplechr].append([int(chr_start), int(chr_end)])
    
    print('Get LIR positions, directly use the symmetric ones I obtained from LIRSat, those that already had F and R arm differences')
    target_LIR_list = []
    with open('/home/lain/aaa_data/run0/new_work_dir//chr2LIRinfo/2_FRtime/LIR_stat1s_revised_manual', 'r') as f:
        next(f)
        for line in f:
            eachline_arr = line.strip().split('\t')
            if len(eachline_arr) < 11:
                continue
            ID, sample, sample_type, chromosome, border1, border2, LIRSat1, LIRSat2, LIRSat3, mark, LIRclass = eachline_arr[:11]
            target_LIR_list.append(['LIR', sample, chromosome, int(border1), int(border2)])
    
    print('Example')
    print(target_VSat6_list[0])
    print(target_LIR_list[0])
    
    print('Iterate through each LIR/VSat6')
    
    def run_step(target_one):
        target_type, sample, chromosome, start, end = target_one
        samplechr = f"{sample}|||{chromosome}"
        
        ## Read TE file
        TE_list = []
        try:
            with open(f"./chr2EDTA/1_EDTA/{sample}/{sample}.fasta.mod.EDTA.TEanno.gff3", 'r') as f:
                next(f)
                for line in f:
                    eachline_arr = line.strip().split('\t')
                    if len(eachline_arr) != 9:
                        continue
                    one_chromosome = eachline_arr[0]
                    one_type = eachline_arr[2]
                    one_start = int(eachline_arr[3])
                    one_end = int(eachline_arr[4])
                    if one_end - one_start > 30000:
                        continue
                    if one_type not in ["Gypsy_LTR_retrotransposon", "LTR_retrotransposon", "PIF_Harbinger_TIR_transposon",
                                       "Copia_LTR_retrotransposon", "CACTA_TIR_transposon", "hAT_TIR_transposon",
                                       "helitron", "Mutator_TIR_transposon", "Tc1_Mariner_TIR_transposon"]:
                        continue
                    if one_chromosome == chromosome:
                        TE_list.append([one_type, one_start, one_end])
        except FileNotFoundError:
            print(f"Warning: TE file not found for {sample}")
            return False
        
        satellite_list = dict_samplechr_VSat1poslist.get(samplechr, [])
        
        if target_type == "LIR":
            center = start + (end - start) / 2
            start = center - 50000
            end = center + 50000
            region_length = 100000
            te_count, te_total_length = calculate_TE_coverage(TE_list, start, end)
            return sample, chromosome, target_type, region_length, te_count, te_total_length
            
        elif target_type == "VSat6":
            region_length = 100000
            
            # Left region (start-50000 to start)
            left_start = start - 50000
            left_end = start
            # Calculate overlap with VSat1 (exclude intervals)
            left_exclude_intervals = []
            for vsat1_start, vsat1_end in satellite_list:
                overlap_start = max(left_start, vsat1_start)
                overlap_end = min(left_end, vsat1_end)
                if overlap_start < overlap_end:
                    left_exclude_intervals.append([overlap_start, overlap_end])
            # Calculate TE coverage in left region (excluding VSat1)
            te_count1, te_total_length1 = calculate_TE_coverage(TE_list, left_start, left_end, left_exclude_intervals)
            
            # Right region (end to end+50000)
            right_start = end
            right_end = end + 50000
            # Calculate overlap with VSat1 (exclude intervals)
            right_exclude_intervals = []
            for vsat1_start, vsat1_end in satellite_list:
                overlap_start = max(right_start, vsat1_start)
                overlap_end = min(right_end, vsat1_end)
                if overlap_start < overlap_end:
                    right_exclude_intervals.append([overlap_start, overlap_end])
            # Calculate TE coverage in right region (excluding VSat1)
            te_count2, te_total_length2 = calculate_TE_coverage(TE_list, right_start, right_end, right_exclude_intervals)
            
            te_count = te_count1 + te_count2
            te_total_length = te_total_length1 + te_total_length2
            return sample, chromosome, target_type, region_length, te_count, te_total_length

    with open('./chr2LIRinfo/3_LIRVSat6_TE_diff/1_TEdensity', 'w') as f:
        f.write('sample\tchromosome\ttype\tregion_length\tTEnum\tTElen\n')
    
    input_list = target_VSat6_list + target_LIR_list
    from multiprocessing import Pool
    import sys
    
    with Pool(processes=50) as pool:
        for i, result in enumerate(pool.imap(run_step, input_list), start=1):
            if result != False:
                sample, chromosome, target_type, result_length, result_TEnum, result_TElen = result
                with open('./chr2LIRinfo/3_LIRVSat6_TE_diff/1_TEdensity', 'a') as f:
                    f.write(f'{sample}\t{chromosome}\t{target_type}\t{result_length}\t{result_TEnum}\t{result_TElen}\n')
            progress = (i / len(input_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush()
    
    with open('./chr2LIRinfo/3_LIRVSat6_TE_diff/1_TEdensity', 'a') as f2:
        with open('./chr2EDTA/allaverage', 'r') as f:
            next(f)
            for line in f:
                eachline_arr = line.strip().split('\t')
                one_sample, chrom, chr_total_len, te_total_count, te_total_len, region_type = eachline_arr
                f2.write(f"{one_sample}\t{chrom}\tchr_average\t{chr_total_len}\t{te_total_count}\t{te_total_len}\n")

    with open('./chr2LIRinfo/3_LIRVSat6_TE_diff/1_TEdensity', 'a') as f2:
        with open('./chr2EDTA/allaverage_arm', 'r') as f:
            next(f)
            for line in f:
                eachline_arr = line.strip().split('\t')
                one_sample, chrom, chr_total_len, te_total_count, te_total_len, region_type = eachline_arr
                f2.write(f"{one_sample}\t{chrom}\tchr_arm_average\t{chr_total_len}\t{te_total_count}\t{te_total_len}\n")

    with open('./chr2LIRinfo/3_LIRVSat6_TE_diff/1_TEdensity', 'a') as f2:
        with open('./chr2EDTA/allaverage_cen', 'r') as f:
            next(f)
            for line in f:
                eachline_arr = line.strip().split('\t')
                one_sample, chrom, chr_total_len, te_total_count, te_total_len, region_type = eachline_arr
                f2.write(f"{one_sample}\t{chrom}\tchr_cen_average\t{chr_total_len}\t{te_total_count}\t{te_total_len}\n")

if argv1 == 'step3_LIRVSat6_TE_diff_norm':
    dict_samplechr_info = {}
    with open('./chr2LIRinfo/3_LIRVSat6_TE_diff/1_TEdensity', 'r') as f:
        next(f)
        for line in f:
            eachline_arr = line.strip().split('\t')
            sample, chromosome, regiontype, region_length, TEnum, TElen = eachline_arr
            if regiontype == 'chr_average':
                dict_samplechr_info[f"{sample}|||{chromosome}"] = [int(TEnum) / int(region_length), float(TElen) / int(region_length)]
    with open('./chr2LIRinfo/3_LIRVSat6_TE_diff/1_TEdensity_norm', 'w') as f2:
        f2.write(f"sample\tchromosome\ttype\tregion_length\tTEnum\tTElen\tTEnum_norm\tTElen_norm\n")
        with open('./chr2LIRinfo/3_LIRVSat6_TE_diff/1_TEdensity', 'r') as f:
            next(f)
            for line in f:
                eachline_arr = line.strip().split('\t')
                sample, chromosome, regiontype, region_length, TEnum, TElen = eachline_arr
                if regiontype == 'chr_average':
                    continue
                chraverage_num, chraverage_len = dict_samplechr_info[f"{sample}|||{chromosome}"]
                density_num_norm = round(int(TEnum) / int(region_length) / chraverage_num, 4)
                #print(TElen)
                density_len_norm = round(float(TElen) / int(region_length) / chraverage_len, 4)
                f2.write(f"{sample}\t{chromosome}\t{regiontype}\t{region_length}\t{TEnum}\t{TElen}\t{density_num_norm}\t{density_len_norm}\n")

if argv1 == 'step3_LIRVSat6_TE_diff' or argv1 == 'step3_LIRVSat6_TE_diff_plot':
    R_txt = f'''
            library(ggplot2)
            library(dplyr)
            
            print("")
           
              input_file=read.table('1_TEdensity', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
#input_file <- input_file %>% filter(chromosome %in% c("Chr10", "Chr11", "Chr15", "Chr16") | type %in% c("LIR"))
    color_values=c(

    )        

# Draw percentage stacked bar chart
p <- ggplot()
p <- p + geom_point(data = input_file, aes(x = TEnum/region_length*1000000, y = TElen/region_length*1000000,color=type), shape = 16,stroke = 0)

    p =p+theme_classic() +         
        theme(
          #axis.ticks.y = element_blank(),
          #axis.text.y = element_blank(),
          #legend.position = "none",
          #axis.text.x = element_blank()
        ) 
      
      # Save as PDF
      pdf(file = paste0('1_TEdensity', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
      print(p)
      dev.off()

# Correct box plot
p <- ggplot(input_file, aes(x = type, y = TEnum/region_length*1000000, fill = type)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_discrete(name = "TE Type")  # If fill needs to be set, use scale_fill_*
  # No need for scale_x_continuous()
  facet_wrap(~metric, scales = "free_y")
  # Save as PDF
  pdf(file = paste0('1_TEdensity_boxplot_TEnum', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
  print(p)
  dev.off()  
  
# Correct box plot
p <- ggplot(input_file, aes(x = type, y = TElen/region_length*1000000, fill = type)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_discrete(name = "TE Type")  # If fill needs to be set, use scale_fill_*
  # No need for scale_x_continuous()
  facet_wrap(~metric, scales = "free_y")
  # Save as PDF
  pdf(file = paste0('1_TEdensity_boxplot_TElen', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
  print(p)
  dev.off()   
  
  
  library(ggplot2)
library(dplyr)

# Read data (assuming data contains columns: type, TEnum, TElen, region_length, TEnum_se, TElen_se)
input_file <- read.table('1_TEdensity', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')

# Calculate normalized values (TE count/length per million base pairs)
input_file <- input_file %>%
  mutate(
    TEnum_norm = TEnum / region_length * 1000000,
    TElen_norm = TElen / region_length * 1000000
  )

# Double error bar scatter plot
p_error <- ggplot(input_file, aes(x = TEnum_norm, y = TElen_norm, color = type)) +
  geom_point(shape = 16, stroke = 0, size = 3) +  # Scatter points
  geom_errorbar(aes(ymin = TElen_norm - TElen_se, ymax = TElen_norm + TElen_se), 
                width = 0.2, alpha = 0.7) +       # Vertical error bars (y direction)
  geom_errorbarh(aes(xmin = TEnum_norm - TEnum_se, xmax = TEnum_norm + TEnum_se), 
                 height = 0.2, alpha = 0.7) +     # Horizontal error bars (x direction)
  theme_classic() +
  labs(
    x = "TE Number per Mb",
    y = "TE Length per Mb",
    color = "TE Type"
  ) +
  theme(
    legend.position = "right",
    text = element_text(size = 12)
  )

# Save as PDF
pdf(file = paste0('1_TEdensity_errorbar_plot', ".pdf"), width = 20 / 2.54, height = 20 / 2.54)
print(p_error)
dev.off()



    '''
    with open(f'./chr2LIRinfo/3_LIRVSat6_TE_diff/plot.R', 'w', encoding='utf-8') as f:
        f.write(R_txt)
    new_directory = f'./chr2LIRinfo/3_LIRVSat6_TE_diff/'
    os.chdir(new_directory)
    subprocess.run(['Rscript plot.R'], shell=True)
    os.chdir('../../../')

if argv1 == 'step3_LIRVSat6_TE_diff_plot2D':
    R_txt = f'''
          library(ggplot2)
library(dplyr)

# Read data
input_file <- read.table(
  "1_TEdensity",  # Replace with your filename
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  sep = "\t"
)

# Calculate normalized values (TE count/length per million base pairs)
input_file <- input_file %>%
  mutate(
    TEnum_norm = TEnum / region_length * 1000000,
    TElen_norm = TElen / region_length * 1000000
  )

# Calculate mean and standard deviation for each type
summary_data <- input_file %>%
  group_by(type) %>%
  summarise(
    mean_TEnum = mean(TEnum_norm, na.rm = TRUE),
    mean_TElen = mean(TElen_norm, na.rm = TRUE),
    sd_TEnum = sd(TEnum_norm, na.rm = TRUE),
    sd_TElen = sd(TElen_norm, na.rm = TRUE)
  )
p_summary <- ggplot(summary_data, aes(x = mean_TEnum, y = mean_TElen, color = type)) +
  geom_point(shape = 16, stroke = 0, size = 4) +  # Thicker mean points
  geom_errorbar(
    aes(ymin = mean_TElen - sd_TElen, ymax = mean_TElen + sd_TElen),
    width = 0.2, alpha = 0.7
  ) +  # Vertical error bars (y direction)
  geom_errorbarh(
    aes(xmin = mean_TEnum - sd_TEnum, xmax = mean_TEnum + sd_TEnum),
    height = 0.2, alpha = 0.7
  ) +  # Horizontal error bars (x direction)
  theme_classic() +
  labs(
    x = "Mean TE Number per Mb",
    y = "Mean TE Length per Mb",
    color = "TE Type"
  ) +
  theme(
    legend.position = "right",
    text = element_text(size = 12)
  )

# Save as PDF
pdf(file = "1_TEdensity_summary_errorbar_plot.pdf", width = 20 / 2.54, height = 20 / 2.54)
print(p_summary)
dev.off()



    '''
    with open(f'./chr2LIRinfo/3_LIRVSat6_TE_diff/plot.R', 'w', encoding='utf-8') as f:
        f.write(R_txt)
    new_directory = f'./chr2LIRinfo/3_LIRVSat6_TE_diff/'
    os.chdir(new_directory)
    subprocess.run(['Rscript plot.R'], shell=True)
    os.chdir('../../../')
if argv1 == 'step3_LIRVSat6_TE_diff_plot2D_norm':
    R_txt = f'''
          library(ggplot2)
library(dplyr)

# Read data
input_file <- read.table(
  "1_TEdensity_norm",  # Replace with your filename
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  sep = "\t"
)


# Calculate mean and standard deviation for each type
summary_data <- input_file %>%
  group_by(type) %>%
  summarise(
    mean_TEnum = mean(TEnum_norm, na.rm = TRUE),
    mean_TElen = mean(TElen_norm, na.rm = TRUE),
    sd_TEnum = sd(TEnum_norm, na.rm = TRUE),
    sd_TElen = sd(TElen_norm, na.rm = TRUE)
  )
p_summary <- ggplot(summary_data, aes(x = mean_TEnum, y = mean_TElen, color = type)) +
  geom_point(shape = 16, stroke = 0, size = 4) +  # Thicker mean points
  geom_errorbar(
    aes(ymin = mean_TElen - sd_TElen, ymax = mean_TElen + sd_TElen),
    width = 0, alpha = 1
  ) +  # Vertical error bars (y direction)
  geom_errorbarh(
    aes(xmin = mean_TEnum - sd_TEnum, xmax = mean_TEnum + sd_TEnum),
    height = 0, alpha = 1
  ) +  # Horizontal error bars (x direction)
  theme_classic() +
  labs(
    x = "Mean TE Number per Mb",
    y = "Mean TE Length per Mb",
    color = "TE Type"
  ) + coord_equal()+
  theme(
    legend.position = "right",
    text = element_text(size = 12)
  )

# Save as PDF
pdf(file = "1_TEdensity_summary_errorbar_plot2.pdf", width = 20 / 2.54, height = 20 / 2.54)
print(p_summary)
dev.off()



    '''
    with open(f'./chr2LIRinfo/3_LIRVSat6_TE_diff/plot.R', 'w', encoding='utf-8') as f:
        f.write(R_txt)
    new_directory = f'./chr2LIRinfo/3_LIRVSat6_TE_diff/'
    os.chdir(new_directory)
    subprocess.run(['Rscript plot.R'], shell=True)
    os.chdir('../../../')

if argv1 == 'step4_time':
    print('Pairwise comparison')
    if os.path.exists('./chr2LIRinfo/4_compare') == False:
        subprocess.run(["mkdir chr2LIRinfo/4_compare"], shell=True)
    target_list = []
    target_list.append(['V126.hap2', 'Chr2', 14692428, 15030566])
    #target_list.append(['V126.hap2','Chr3',15358845,15777811])   
    #target_list.append(['V126.hap2','Chr4',14056925,14276580])   
    target_list.append(['V126.hap2', 'Chr5', 15178558, 15569106])
    target_list.append(['V126.hap2', 'Chr6', 12110223, 12469057])
    #target_list.append(['V126.hap2','Chr8',7890307,8123339])   
    #target_list.append(['V126.hap2','Chr12',12275803,12799879])   
    #target_list.append(['V126.hap2','Chr14',14994034,15429422])   
    #target_list.append(['V126.hap2','Chr18',16544108,17045390])   
    #target_list.append(['V126.hap2','Chr19',17685851,18107892])   
    
    target_list.append(['PN40024', 'Chr2', 12939361, 13300832])
    #target_list.append(['PN40024_hap1','Chr3',14353397,14553397])
    #target_list.append(['PN40024','Chr4',13291174,13483935])
    target_list.append(['PN40024', 'Chr5', 14008287, 14357682])
    target_list.append(['PN40024', 'Chr6', 10606684, 10966038])
    #target_list.append(['PN40024','Chr8',6679384,7071344])
    #target_list.append(['PN40024','Chr9',15113362,15449399])   
    target_num = len(target_list)
    kk = 0
    jj = 1
    input_list = []
    '''
    while kk<target_num:
        while jj<target_num:
            input_list.append([kk,jj])
            jj+=1
       kk+=1    '''
    # Generate all index pairs for pairwise comparison
    for kk in range(target_num):
        #for jj in range(kk + 1, target_num):  # Avoid duplicate comparisons and self-comparison
        for jj in range(kk, target_num):  # Avoid duplicate comparisons and self-comparison
            input_list.append((kk, jj))
    def run_step1(one):
        kk, jj = one
        if kk == jj:
            return False
        target1 = target_list[kk]
        target2 = target_list[jj]
        name = f"{kk}_{jj}"
        para = f"{target1[0]}:{target1[1]}:{target1[2]}-{target1[3]} {target2[0]}:{target2[1]}:{target2[2]}-{target2[3]} "
        #subprocess.run([f'python /home/lain/aaa_data/run0/new_work_dir/chr2aligntime.py {para} ./chr2LIRinfo//4_compare/{name}  -seq2reverse no'  ], shell=True)       
        
    with Pool(processes=5) as pool:
        # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step1, input_list), start=1):
            # Process results here, e.g., store or print
            progress = (i / len(input_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush()
    with open(f"./chr2LIRinfo/4_compare/all_result", 'w') as f2:
        f2.write(f"name\tkk\tjj\ttarget1\ttarget2\talign_len\tmismatch\tpercent\n")
        for one in input_list:
            kk, jj = one
            target1 = target_list[kk]
            target2 = target_list[jj]
            name = f"{kk}_{jj}"
            if kk == jj:
                align_len = int(target1[3]) - int(target1[2]) + 1
                mismatch = 0
                percent = 1
            else:
                with open(f'./chr2LIRinfo//4_compare/{name}/ALL_SUM', 'r') as f:
                    align_len, mismatch, percent = f.read().strip().split('\t')
            f2.write(f"{name}\t{kk}\t{jj}\t{target1}\t{target2}\t{align_len}\t{mismatch}\t{percent}\n")
        
    R_txt = f'''
            library(ggplot2)
            library(dplyr)
            
            print("")
           
              input_file=read.table('all_result', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')

# Draw percentage stacked bar chart
p <- ggplot()
p <- p + geom_point(data = input_file, aes(x = kk, y = jj,color=percent,size=align_len), shape = 16,stroke = 0)

    p =p+theme_classic() +         
        theme(
          #axis.ticks.y = element_blank(),
          #axis.text.y = element_blank(),
          #legend.position = "none",
          #axis.text.x = element_blank()
        ) 
      
      # Save as PDF
      pdf(file = paste0('all_result', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
      print(p)
      dev.off()


    '''
    with open(f'./chr2LIRinfo/4_compare/plot.R', 'w', encoding='utf-8') as f:
        f.write(R_txt)
    new_directory = f'./chr2LIRinfo/4_compare/'
    os.chdir(new_directory)
    subprocess.run(['Rscript plot.R'], shell=True)
    os.chdir('../../../')

print("\n\n")
print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
time_end = timeit.default_timer()
print('All the running time: %.0f Seconds' % (time_end - time_start))