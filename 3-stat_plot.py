#!/usr/bin/python
# -*- coding: UTF-8 -*-
import sys
argvs=sys.argv
argv_len=len(sys.argv)
if argv_len < 2 or argvs[1] in ("h", "-h", "help", "-help") or not argvs[1].startswith("step") :
    print ("3-stat_plot.py-----help:")
    print ("step0 extracts and organizes chr2blast results for samples_satellite.py. The remaining steps are for additional visualizations.\n\n")
    
    print ("Usage：")
    print ("3-stat_plot.py step0      #Extract and summarize ./new_work_dir/chr2blast/1_blastn/sum and ./samples_satellite/2_good_regions")
    print ("")
    
    print ("3-stat_plot.py step1a/b   Generate some plots")
    print ("3-stat_plot.py step2a     Generate raw data for heatmap")
    
    print ("3-stat_plot.py step4      Determine whether LIR1 and VSat6(cen_383) can be on the same chromosome; LIR1 and VSat6 have high similarity")
    print ("3-stat_plot.py step5      Calculate the proportion of LIRSat2 near VSat6")
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
import re # For regex processing
from multiprocessing import Pool, cpu_count

time_start=timeit.default_timer()
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  

 




if argv1=="stepall" or argv1=="step0" :    
    print("Summarize the 2_good_regions file and try to plot")
    subprocess.run(["mkdir ./stat_plot"], shell=True)
    def chrname_old2new(one_id):
        one_id_lower=one_id.lower()
        if one_id.isdigit():
            one_id="Chr"+one_id
        elif one_id_lower.endswith('chr01') or one_id_lower.endswith('chr1') or one_id_lower=='region_1':   one_id="Chr1";      
        elif one_id_lower.endswith('chr02') or one_id_lower.endswith('chr2') or one_id_lower=='region_2':   one_id="Chr2";      
        elif one_id_lower.endswith('chr03') or one_id_lower.endswith('chr3') or one_id_lower=='region_3':   one_id="Chr3";      
        elif one_id_lower.endswith('chr04') or one_id_lower.endswith('chr4') or one_id_lower=='region_4':   one_id="Chr4";      
        elif one_id_lower.endswith('chr05') or one_id_lower.endswith('chr5') or one_id_lower=='region_5':   one_id="Chr5";      
        elif one_id_lower.endswith('chr06') or one_id_lower.endswith('chr6') or one_id_lower=='region_6':   one_id="Chr6";      
        elif one_id_lower.endswith('chr07') or one_id_lower.endswith('chr7') or one_id_lower=='region_7':   one_id="Chr7";      
        elif one_id_lower.endswith('chr08') or one_id_lower.endswith('chr8') or one_id_lower=='region_8':   one_id="Chr8";      
        elif one_id_lower.endswith('chr09') or one_id_lower.endswith('chr9') or one_id_lower=='region_9':   one_id="Chr9";  
        elif one_id_lower.endswith('chr10') or one_id_lower=='region_10':   one_id="Chr10";      
        elif one_id_lower.endswith('chr11') or one_id_lower=='region_11':   one_id="Chr11";      
        elif one_id_lower.endswith('chr12') or one_id_lower=='region_12':   one_id="Chr12";      
        elif one_id_lower.endswith('chr13') or one_id_lower=='region_13':   one_id="Chr13";      
        elif one_id_lower.endswith('chr14') or one_id_lower=='region_14':   one_id="Chr14";      
        elif one_id_lower.endswith('chr15') or one_id_lower=='region_15':   one_id="Chr15";      
        elif one_id_lower.endswith('chr16') or one_id_lower=='region_16':   one_id="Chr16";      
        elif one_id_lower.endswith('chr17') or one_id_lower=='region_17':   one_id="Chr17";      
        elif one_id_lower.endswith('chr18') or one_id_lower=='region_18':   one_id="Chr18";      
        elif one_id_lower.endswith('chr19') or one_id_lower=='region_19':   one_id="Chr19";      
        elif one_id_lower.endswith('chr20') or one_id_lower=='region_20':   one_id="Chr20";  
        else: one_id='other'
        return one_id
    print('Summarizing')
    print('Count the number and length of blocks without distinguishing chromosomes, by sample')
    dict_sample_type_num1len={}
    
    with open('./stat_plot/0-region2info_raw','w') as f2:
        f2.write(f"sample\tchromosome\tchromosome_new\tcentype\tchr_start\tchr_end\tlength\tstrand\tmatch_percent\n")
        ### BLAST results
        with open('./new_work_dir/chr2blast/1_blastn/sum','r') as f:    
            for line in f:
                eachline_arr=line.strip().split('\t')
                #f2.write("#\tsample\tcenid\tchr\tstrand\tchr_start\tchr_end\tlength\tmatch_len\tmatch_percent\n")
                if len(eachline_arr)!=9:continue
                sample=         eachline_arr[0]
                chromosome=     eachline_arr[2]
                chromosome_new=chrname_old2new(chromosome)
                centype=        eachline_arr[1]
                chr_start=      eachline_arr[4]
                chr_end=        eachline_arr[5]
                length=         eachline_arr[6]
                strand=       eachline_arr[3]
                match_percent=eachline_arr[8]
                f2.write(f"{sample}\t{chromosome}\t{chromosome_new}\t{centype}\t{chr_start}\t{chr_end}\t{length}\t{strand}\t{match_percent}\n")

            
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=8:continue
                if eachline_arr[0]=='sample':continue
                #headline='sample\tregion_name\tregion_pos\tbigblock_chrstart\tbigblock_chrend\tchr_region_length\tstrand\tmatch_percent\n'
                #region_1	V008#1#chr01:1-25278865	18450941	18452876	14611267	14613090	1936	1824	minus	1869	0.965	1514	0.782
                sample=         eachline_arr[0]
                chromosome=     eachline_arr[2].split(':')[0]
                chromosome_new=chrname_old2new(chromosome)
                centype=        'cen_107'
                chr_start=      eachline_arr[3]
                chr_end=        eachline_arr[4]
                length=         eachline_arr[5]
                strand=       eachline_arr[6]
                match_percent=       eachline_arr[7]
                f2.write(f"{sample}\t{chromosome}\t{chromosome_new}\t{centype}\t{chr_start}\t{chr_end}\t{length}\t{strand}\t{match_percent}\n")

        if 1==2:    ###Summary of muscadine grape (Vitis rotundifolia) samples after processing with samples_satellite.py
            with open('./samples_satellite_RO/2_good_regions','r') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if len(eachline_arr)!=8:continue
                    if eachline_arr[0]=='sample':continue
                    #headline='sample\tregion_name\tregion_pos\tbigblock_chrstart\tbigblock_chrend\tchr_region_length\tstrand\tmatch_percent\n'
                    sample=         eachline_arr[0]
                    chromosome=     eachline_arr[2].split(':')[0]
                    chromosome_new=chrname_old2new(chromosome)
                    centype=        'cen_RO21'
                    chr_start=      eachline_arr[3]
                    chr_end=        eachline_arr[4]
                    length=         eachline_arr[5]
                    strand=       eachline_arr[6]
                    match_percent=       eachline_arr[7]
                    f2.write(f"{sample}\t{chromosome}\t{chromosome_new}\t{centype}\t{chr_start}\t{chr_end}\t{length}\t{strand}\t{match_percent}\n")
    subprocess.run([""" 
                        head -n 1 ./stat_plot/0-region2info_raw > ./stat_plot/head
                        tail -n +2 ./stat_plot/0-region2info_raw > ./stat_plot/tail
                        sort  -k1,1 -k3,3V -k5,5n ./stat_plot/tail > ./stat_plot/tail_sort
                        cat ./stat_plot/head ./stat_plot/tail_sort > ./stat_plot/0-region2info
                        rm ./stat_plot/head ./stat_plot/tail ./stat_plot/tail_sort  ./stat_plot/0-region2info_raw
                    """], shell=True)
    
if argv1=="stepall" or argv1=="step1" or argv1=="step1a": 
    subprocess.run(["mkdir ./stat_plot/step1-NoChromosome"], shell=True)
    print('Count the number of each group, without distinguishing chromosomes')
    dict_sample_type_num1len={}
    with open('./stat_plot/0-region2info','r') as f:
        for line in f:
            eachline_arr=line.strip().split('\t')
            if eachline_arr[0]=='sample':continue
            ##{sample}\t{chromosome}\t{chromosome_new}\t{centype}\t{chr_start}\t{chr_end}\t{length}\t{strand}\t{match_percent}
            sample=                 eachline_arr[0]
            chromosome=             eachline_arr[1]
            chromosome_new=         eachline_arr[2]
            centype=                eachline_arr[3]
            chr_start=              eachline_arr[4]
            chr_end=                eachline_arr[5]
            length=                 eachline_arr[6]
            strand=                 eachline_arr[7]
            match_percent=          eachline_arr[8]
            #########
            if sample not in dict_sample_type_num1len:
                dict_sample_type_num1len[sample]={}
            if centype not in dict_sample_type_num1len[sample]:
                dict_sample_type_num1len[sample][centype]={}
                dict_sample_type_num1len[sample]
                dict_sample_type_num1len[sample][centype][0]=0
                dict_sample_type_num1len[sample][centype][1]=0
            dict_sample_type_num1len[sample][centype][0]+=1
            dict_sample_type_num1len[sample][centype][1]+=int(length)                
    with open('./stat_plot/step1-NoChromosome/1-region2centype_stat1','w') as f3:
        f3.write("sample\trepeat_type\tblock_num\tblock_length\n")
        for sample,type_info in dict_sample_type_num1len.items():
            for one_type,one_info in type_info.items():
                f3.write(f"{sample}\t{one_type}\t{str(one_info[0])}\t{str(one_info[1])}\n")
                
if argv1=="stepall" or argv1=="step1" or argv1=="step1b":
    
    R_txt='''library(ggplot2)
library(dplyr)
setwd('./')
input_file1 <- read.table('1-region2centype_stat1', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
input_file1$repeat_type <- factor(input_file1$repeat_type, levels = sort(unique(input_file1$repeat_type)))
p <- ggplot(input_file1, aes(x = repeat_type, y = block_length,fill=repeat_type)) +
    geom_boxplot(outlier.shape = NA) +  
    geom_point(position = position_jitter(width = 0.30, height = 0), alpha = 0.1,size=0.9) +  
    theme_classic()   
p <- p + scale_y_continuous(
  name = "Repeat Length (M)",  # y-axis label
  trans = 'log10',        # Use logarithmic transformation
)  +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Display x-axis labels vertically
  )
  
# Set PDF device size to 10cm x 10cm (convert to inches)
pdf("pic1_grape.pdf", width = 25 / 2.54, height = 15 / 2.54)


# Print the plot
print(p) 
'''
    with open('./stat_plot/step1-NoChromosome/pic_ggplot_step1b.R','w',encoding='utf-8') as f:
        f.write(R_txt)
    new_directory = "./stat_plot/step1-NoChromosome/"
    os.chdir(new_directory)
    subprocess.run(['Rscript pic_ggplot_step1b.R'], shell=True)  
    os.chdir('../../../')
    
if argv1=="stepall" or argv1=="step2" or argv1=="step2a": 
    subprocess.run(["mkdir ./stat_plot/step2-WithChromosome"], shell=True)
    subprocess.run(["mkdir ./stat_plot/step2-WithChromosome/centype"], shell=True)
    print('Count the number of each group, distinguishing chromosomes, each chromosome is counted separately for heatmap plotting')
    dict_type_chr_sample_num1len={}
    with open('./stat_plot/0-region2info','r') as f:
        for line in f:
            eachline_arr=line.strip().split('\t')
            if eachline_arr[0]=='sample':continue
            ##{sample}\t{chromosome}\t{chromosome_new}\t{centype}\t{chr_start}\t{chr_end}\t{length}\t{strand}\t{match_percent}
            sample=                 eachline_arr[0]
            #print(sample)
            if   sample[-5:]=='.hap1' :sample=sample[:-5];  hap='_hap1'
            elif sample[-5:]=='_hap1' :sample=sample[:-5];  hap='_hap1'
            elif sample[-5:]=='.hap2' :sample=sample[:-5];  hap='_hap2'
            elif sample[-5:]=='_hap2' :sample=sample[:-5];  hap='_hap2'
            else:continue
            chromosome=             eachline_arr[1]
            chromosome_new=         eachline_arr[2]
            if chromosome_new=='other':continue
            chromosome_new2=        chromosome_new+hap
            centype=                eachline_arr[3]
            chr_start=              eachline_arr[4]
            chr_end=                eachline_arr[5]
            length=                 eachline_arr[6]
            strand=                 eachline_arr[7]
            match_percent=          eachline_arr[8]
            #########
            if centype not in dict_type_chr_sample_num1len:
                dict_type_chr_sample_num1len[centype]={}
            if chromosome_new2 not in dict_type_chr_sample_num1len[centype]:
                dict_type_chr_sample_num1len[centype][chromosome_new2]={}
            if sample not in dict_type_chr_sample_num1len[centype][chromosome_new2]:
                dict_type_chr_sample_num1len[centype][chromosome_new2][sample]={}
                dict_type_chr_sample_num1len[centype][chromosome_new2][sample][0]=0
                dict_type_chr_sample_num1len[centype][chromosome_new2][sample][1]=0
            dict_type_chr_sample_num1len[centype][chromosome_new2][sample][0]+=1
            dict_type_chr_sample_num1len[centype][chromosome_new2][sample][1]+=int(length)                
    
    #print(dict_type_chr_sample_num1len)
    centype_set=set();    sample_set=set();     chr_set=set()
    for centype,centype_info in dict_type_chr_sample_num1len.items():
        centype_set.add(centype)
        for chromosome,chromosome_info in centype_info.items():
            chr_set.add(chromosome)
            for sample,sample_info in chromosome_info.items():
                sample_set.add(sample)
    centype_list=list(centype_set);     centype_list.sort()
    chr_list=list(chr_set);             chr_list.sort()
    sample_list=list(sample_set);       sample_list.sort()
    print(centype_list)   
    print(sample_list)   
    print(chr_list)  
    #centype_list=['cen_103', 'cen_107', 'cen_145', 'cen_187', 'cen_191', 'cen_355', 'cen_66',  'cen_383', 'cen_RO197', 'cen_RO21', 'cen_RO60', 'cen_RO72', 'LIR_Sat1', 'LIR_Sat2', 'LIR_Sat3']
    centype_list=['LIRSat1', 'LIRSat2', 'LIRSat3']
    #sample_list=['PN40024','Chardonnay','VHP-T2T', 'PinotNoir', 'PinotNoir2',
    #'Baimunage', 'BlackMonukka', 'DavidiiGrape', 'Hongmunage', 'ManicureFinger', 'MuscatHamburg',  'PiasezkiiGrape',  'ShineMuscat', 'ThompsonSeedless', 'ThompsonSeedless2',
     #'V074','V087', 'V088', 'V091', 'V092', 'V093','V005', 'V007', 'V008', 'V012', 'V015', 'V018', 'V019', 'V020', 'V022', 'V023',  'V030', 'V031', 'V032', 'V033', 'V034', 'V036', 'V037', 'V038', 'V040', 'V041', 'V043', 'V048', 'V049', 'V050', 'V051', 'V052', 'V053', 'V055', 'V058', 'V059', 'V060', 'V061', 'V062', 'V063', 'V064', 'V065', 'V066', 'V067', 'V069', 'V070', 'V072', 'V076', 'V077', 'V079', 'V081', 'V096', 'V098', 'V099', 'V100', 'V102', 'V105', 'V106', 'V107', 'V108', 'V112', 'V117', 'V120', 'V123', 'V124', 'V125', 'V126',  'WoollyGrape']

    #Optimized the order
    sample_list=['Baimunage','BlackMonukka','V059','V063','V066','V070','Chardonnay','VHP-T2T','V091','Hongmunage','V074','ManicureFinger','V092','V062','MuscatHamburg','V088','PinotNoir','PinotNoir2','V069','PN40024','V087','V061','V065','V064','ThompsonSeedless' ,'ThompsonSeedless2','V093','V060','V067','V058','V107','V108','V081','V015','V105','V126','V112','DavidiiGrape','V007','V008','V012','V100','V099','V098','V117','V120','PiasezkiiGrape','V123','V124','V125','WoollyGrape','V106','V102','V079','V030','V031','V005','V048','V050','V051','V049','V038','V043','V052','V034','V023','V032','V033','V053','V020','V022','V040','V019','V041','V077','V037','V096','ShineMuscat','V036','V076','V055','V018','V072']    
        
        
    chr_list=['Chr1_hap1', 'Chr1_hap2', 'Chr2_hap1', 'Chr2_hap2', 'Chr3_hap1', 'Chr3_hap2', 'Chr4_hap1', 'Chr4_hap2', 'Chr5_hap1', 'Chr5_hap2', 'Chr6_hap1', 'Chr6_hap2', 'Chr7_hap1', 'Chr7_hap2', 'Chr8_hap1', 'Chr8_hap2', 'Chr9_hap1', 'Chr9_hap2', 'Chr10_hap1', 'Chr10_hap2', 'Chr11_hap1', 'Chr11_hap2', 'Chr12_hap1', 'Chr12_hap2', 'Chr13_hap1', 'Chr13_hap2', 'Chr14_hap1', 'Chr14_hap2', 'Chr15_hap1', 'Chr15_hap2', 'Chr16_hap1', 'Chr16_hap2', 'Chr17_hap1', 'Chr17_hap2', 'Chr18_hap1', 'Chr18_hap2', 'Chr19_hap1', 'Chr19_hap2']#,  'Chr20_hap1','Chr20_hap2'
    
    #sys.exit()   
    #print(dict_type_chr_sample_num1len)
    #sys.exit() 
    for one_centype in centype_list:
        #print(one_centype)
        with open(f'./stat_plot/step2-WithChromosome/centype/{one_centype}','w') as f:
            newline='chr\t'+'\t'.join(sample_list)
            f.write(newline+'\n')
            for one_chr in chr_list:
                #print(one_chr)
                value_list=[]
                for one_sample in sample_list:
                    #print(one_sample)
                    #num= dict_type_chr_sample_num1len[one_centype][sample][chromosome_new2][0] 
                    if one_chr not in dict_type_chr_sample_num1len[one_centype]:
                        length_str= '0'
                    elif one_sample not in dict_type_chr_sample_num1len[one_centype][one_chr]:
                        length_str= '0'
                    else:    
                        length_str= str(dict_type_chr_sample_num1len[one_centype][one_chr][one_sample][1] )
                        #print(dict_type_chr_sample_num1len[one_centype][one_chr][sample])
                    #print(length_str)
                    #sys.exit()    
                    value_list.append(length_str)
                #print(value_list)    
                newline=f"{one_chr}\t"+"\t".join(value_list)
                f.write(newline+'\n')
    
if argv1=="stepall" or argv1=="step3" or argv1=="step3.0":
    subprocess.run(["mkdir ./stat_plot/step3_chr_group"], shell=True)
    dict_samplel_group={}
    with open("./samples_satellite/sample_info",'r') as f:
        for line in f:
            eachline_arr=line.strip().split('\t')
            #Baimunage_hap1	Baimunage	 vinifera	Eurasian	Eurasian_vinifera_Baimunage	Table
            if len(eachline_arr)!=6:continue
            sample=             eachline_arr[0]
            species=            eachline_arr[2]
            species_group=      eachline_arr[3]
            dict_samplel_group[sample]=[species,species_group]
    centype_list= ["cen_66","cen_103","cen_107","cen_191","cen_383","cen_355","LIRSat1","LIRSat2","LIRSat3"]
    list_chr = ["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
    for one_centype in centype_list:
        dict_sample_chromosome_length={}
        
        with open("./stat_plot/0-region2info",'r') as f:
            next(f)
            for line in f:
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                if len(eachline_arr)!=9:continue
                sample=             eachline_arr[0]
                chromosome=         eachline_arr[2]
                centype=            eachline_arr[3]
                length=             eachline_arr[6]
                if centype!=one_centype:continue
                
                if sample not in dict_sample_chromosome_length:
                    dict_sample_chromosome_length[sample]={}
                if chromosome not in dict_sample_chromosome_length[sample]:
                    dict_sample_chromosome_length[sample][chromosome]=0
                dict_sample_chromosome_length[sample][chromosome]+=int(length)                  
        with open(f"./stat_plot/step3_chr_group/{one_centype}",'w') as f2:     
            f2.write("sample\tchromosome\tcentype\tlength\tspecies\tgroup\n")
            for sample,dict_chromosome_length in dict_sample_chromosome_length.items():
                species,group=  dict_samplel_group[sample]
                if group not in ["Eurasian","East_Asia","America"]:continue
                for chromosome in list_chr:
                    if chromosome not in dict_chromosome_length:
                        length=0
                    else:
                        length=dict_chromosome_length[chromosome]
                    newline=f"{sample}\t{chromosome}\t{one_centype}\t{length}\t{species.strip()}\t{group.strip()}\n"
                    f2.write(newline)
                    
if argv1=="stepall" or argv1=="step3" or argv1=="step3.0b":
    print("Calculating mean total length per centype per haplotype...")
    
    # Read sample grouping information
    dict_sample_group = {}
    with open("./samples_satellite/sample_info",'r') as f:
        for line in f:
            eachline_arr = line.strip().split('\t')
            if len(eachline_arr) != 6: continue
            sample = eachline_arr[0]
            species = eachline_arr[2]
            species_group = eachline_arr[3]
            dict_sample_group[sample] = [species, species_group]
    
    centype_list = ["cen_66","cen_103","cen_107","cen_191","cen_383","cen_355","LIRSat1","LIRSat2","LIRSat3"]
    
    # Store results: {centype: [lengths_list]}
    results = {ct: [] for ct in centype_list}
    
    for one_centype in centype_list:
        # Aggregate total length by sample
        sample_total = {}
        
        with open("./stat_plot/0-region2info", 'r') as f:
            next(f)  # skip header
            for line in f:
                eachline_arr = line.strip().split('\t')
                if len(eachline_arr) != 9: continue
                sample = eachline_arr[0]
                centype = eachline_arr[3]
                length = int(eachline_arr[6])
                
                if centype != one_centype: continue
                
                # Only keep required groups
                if sample not in dict_sample_group: continue
                group = dict_sample_group[sample][1]
                #if group not in ["Eurasian","East_Asia","America"]: continue
                if group =="Routundifolia":continue
                sample_total[sample] = sample_total.get(sample, 0) + length
        
        # Collect total lengths for all haplotypes
        if sample_total:
            results[one_centype] = list(sample_total.values())
    
    # Output statistics to file
    with open("./stat_plot/step3_chr_group/centype_mean_length.txt", 'w') as f_out:
        f_out.write("centype\tmean_length_per_hap\tsd\tn_hap\tmin\tmax\n")
        for ct in centype_list:
            lengths = results[ct]
            if lengths:
                mean_len = sum(lengths) / len(lengths)
                sd_len = (sum((x - mean_len)**2 for x in lengths) / len(lengths))**0.5
                min_len = min(lengths)
                max_len = max(lengths)
                f_out.write(f"{ct}\t{mean_len:.2f}\t{sd_len:.2f}\t{len(lengths)}\t{min_len}\t{max_len}\n")
            else:
                f_out.write(f"{ct}\tNA\tNA\t0\tNA\tNA\n")
    
    # Print to console
    print("\n=== Mean total length per centype per haplotype ===")
    for ct in centype_list:
        if results[ct]:
            mean_val = sum(results[ct]) / len(results[ct])
            print(f"{ct}: {mean_val:.2f} bp (n={len(results[ct])} haplotypes)")
        else:
            print(f"{ct}: no data")
    
if argv1=="stepall" or argv1=="step3" or argv1=="step3.0c":
    print("Calculating merged block count (>10kb, gap<10kb merged)...")
    
    # Read sample grouping information
    dict_sample_group = {}
    with open("./samples_satellite/sample_info",'r') as f:
        for line in f:
            eachline_arr = line.strip().split('\t')
            if len(eachline_arr) != 6: continue
            sample = eachline_arr[0]
            species = eachline_arr[2]
            species_group = eachline_arr[3]
            dict_sample_group[sample] = [species, species_group]
    
    centype_list = ["cen_66","cen_103","cen_107","cen_191","cen_383","cen_355","LIRSat1","LIRSat2","LIRSat3"]
    list_chr = ["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
    
    # Store results: {centype: [merged_block_counts_per_hap]}
    results = {ct: [] for ct in centype_list}
    
    for one_centype in centype_list:
        # Store all region coordinates by sample + chromosome
        # Data structure: {sample: {chromosome: [(start, end), ...]}}
        sample_chr_regions = {}
        
        with open("./stat_plot/0-region2info", 'r') as f:
            next(f)  # skip header
            for line in f:
                eachline_arr = line.strip().split('\t')
                if len(eachline_arr) != 9: continue
                sample = eachline_arr[0]
                chromosome = eachline_arr[2]
                centype = eachline_arr[3]
                start = int(eachline_arr[4])
                end = int(eachline_arr[5])
                
                if centype != one_centype: continue
                
                # Only keep required groups
                if sample not in dict_sample_group: continue
                group = dict_sample_group[sample][1]
                if group == "Routundifolia": continue
                if sample=='PN40024':continue
                # Store region coordinates
                if sample not in sample_chr_regions:
                    sample_chr_regions[sample] = {}
                if chromosome not in sample_chr_regions[sample]:
                    sample_chr_regions[sample][chromosome] = []
                sample_chr_regions[sample][chromosome].append((start, end))
        
        # Merge regions for each sample's each chromosome
        for sample, chr_dict in sample_chr_regions.items():
            total_merged_blocks = 0
            
            for chromosome, regions in chr_dict.items():
                if not regions:
                    continue
                
                # Sort by start coordinate
                regions.sort(key=lambda x: x[0])
                
                # Merge adjacent regions (gap < 10000 bp considered same block)
                merged_blocks = []
                current_start, current_end = regions[0]
                
                for start, end in regions[1:]:
                    gap = start - current_end
                    if gap < 20000:  # Gap less than 10kb, merge
                        current_end = max(current_end, end)
                    else:  # Gap >= 10kb, end current block
                        block_len = current_end - current_start + 1
                        if block_len > 3000:  # Only keep blocks > 10kb
                            merged_blocks.append((current_start, current_end, block_len))
                        current_start, current_end = start, end
                
                # Process the last block
                block_len = current_end - current_start + 1
                if block_len > 3000:
                    merged_blocks.append((current_start, current_end, block_len))
                
                # Accumulate block count for this chromosome
                total_merged_blocks += len(merged_blocks)
            
            # Record merged block count for this sample
            results[one_centype].append(total_merged_blocks)
    
    # Output statistics to file
    with open("./stat_plot/step3_chr_group/centype_merged_block_count.txt", 'w') as f_out:
        f_out.write("centype\tmean_merged_block_count_per_hap\tsd\tn_hap\tmin\tmax\tall_block_counts\n")
        for ct in centype_list:
            counts = results[ct]
            if counts:
                mean_cnt = sum(counts) / len(counts)
                sd_cnt = (sum((x - mean_cnt)**2 for x in counts) / len(counts))**0.5
                min_cnt = min(counts)
                max_cnt = max(counts)
                # Join all counts with commas for subsequent analysis
                all_counts_str = ",".join([str(x) for x in counts])
                f_out.write(f"{ct}\t{mean_cnt:.2f}\t{sd_cnt:.2f}\t{len(counts)}\t{min_cnt}\t{max_cnt}\t{all_counts_str}\n")
            else:
                f_out.write(f"{ct}\tNA\tNA\t0\tNA\tNA\tNA\n")
    
    # Print to console
    print("\n=== Merged block count (>10kb, gap<10kb merged) per centype per haplotype ===")
    for ct in centype_list:
        counts = results[ct]
        if counts:
            mean_val = sum(counts) / len(counts)
            print(f"{ct}: {mean_val:.2f} blocks (n={len(counts)} haplotypes)")
        else:
            print(f"{ct}: no data")
    
    print("\nOutput saved to: ./stat_plot/step3_chr_group/centype_merged_block_count.txt")
if argv1 == "stepall" or argv1 == "step3" or argv1 == "step3.0d":
    import sys
    import os
    from multiprocessing import Pool
    import subprocess
    
    print("Calculating GC content for each array per haplotype using coordinates (multiprocessing)...")
    
    # Read sample grouping information
    dict_sample_group = {}
    with open("./samples_satellite/sample_info", 'r') as f:
        for line in f:
            eachline_arr = line.strip().split('\t')
            if len(eachline_arr) != 6:
                continue
            sample = eachline_arr[0]
            species = eachline_arr[2]
            species_group = eachline_arr[3]
            dict_sample_group[sample] = [species, species_group]
    
    # Define satellite types to process
    centype_list = ["cen_66", "cen_103", "cen_107", "cen_191", "cen_383", "cen_355", "LIRSat1", "LIRSat2", "LIRSat3"]
    
    # Base path
    base_dir = "/home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare"
    
    # Read 0-region2info, organize all region information by sample
    # Data structure: [(sample, chromosome, centype, start, end, length), ...]
    all_regions = []
    
    with open("./stat_plot/0-region2info", 'r') as f:
        next(f)  # skip header
        for line in f:
            eachline_arr = line.strip().split('\t')
            if len(eachline_arr) != 9:
                continue
            sample = eachline_arr[0]
            chromosome = eachline_arr[2]  # chromosome_new, corresponding to fasta header
            centype = eachline_arr[3]
            start = int(eachline_arr[4])
            end = int(eachline_arr[5])
            length = int(eachline_arr[6])
            
            # Only process specified centypes
            if centype not in centype_list:
                continue
            
            # Only keep required groups
            if sample not in dict_sample_group:
                continue
            group = dict_sample_group[sample][1]
            if sample=='PN40024':continue
            if group == "Routundifolia":
                continue
            
            all_regions.append((sample, chromosome, centype, start, end, length))
    
    print(f"Loaded {len(all_regions)} regions for {len(set(r[0] for r in all_regions))} samples")
    
    # Group by sample, each sample's all regions as a task unit
    sample_to_regions = {}
    for region in all_regions:
        sample = region[0]
        if sample not in sample_to_regions:
            sample_to_regions[sample] = []
        sample_to_regions[sample].append(region)
    
    # Prepare task list: each element is (sample, regions_list, base_dir)
    task_list = [(sample, regions, base_dir) for sample, regions in sample_to_regions.items()]
    
    print(f"Total tasks: {len(task_list)} samples to process")
    
    def process_sample(args):
        """Process all regions for a single sample, return GC results for this sample (merged by centype then calculated)"""
        sample, regions, base_dir = args
        
        hap_path = os.path.join(base_dir, sample)
        fasta_file = os.path.join(hap_path, f"{sample}.fasta")
        
        if not os.path.exists(fasta_file):
            print(f"Warning: {fasta_file} not found, skipping sample {sample}")
            return []
        
        # Read the fasta file for this hap, store sequences for each chromosome
        chr_sequences = {}
        current_chr = None
        current_seq = []
        
        try:
            with open(fasta_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Save previous sequence
                        if current_chr and current_seq:
                            chr_sequences[current_chr] = ''.join(current_seq)
                        # Parse new chromosome name
                        current_chr = line[1:].split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line)
                # Save the last sequence
                if current_chr and current_seq:
                    chr_sequences[current_chr] = ''.join(current_seq)
        except Exception as e:
            print(f"Error reading {fasta_file}: {e}")
            return []
        
        # Group by centype, merge all sequences of the same type
        centype_sequences = {}  # {centype: [sequence1, sequence2, ...]}
        
        for region in regions:
            sample_name, chromosome, centype, start, end, region_len = region
            
            if chromosome not in chr_sequences:
                # Only print warning for unexpected chromosomes to avoid spam
                if chromosome not in ["other", "Other", "OTHER"]:
                    print(f"  Warning: {chromosome} not found in {sample} fasta")
                continue
            
            seq_full = chr_sequences[chromosome]
            seq_len = len(seq_full)
            
            # Coordinate conversion: fasta sequence is 1-indexed
            if start < 1 or end > seq_len:
                extract_start = max(0, start - 1)
                extract_end = min(seq_len, end)
            else:
                extract_start = start - 1
                extract_end = end
            
            subseq = seq_full[extract_start:extract_end]
            
            # Collect sequences by centype
            if centype not in centype_sequences:
                centype_sequences[centype] = []
            centype_sequences[centype].append(subseq)
        
        # For each centype, merge all sequences and calculate GC content
        results = []
        for centype, seq_list in centype_sequences.items():
            # Merge all sequences
            merged_seq = ''.join(seq_list)
            merged_upper = merged_seq.upper()
            
            gc_count = merged_upper.count('G') + merged_upper.count('C')
            total_bases = len(merged_upper)
            
            if total_bases > 0:
                gc_content = (gc_count / total_bases) * 100
            else:
                gc_content = 0.0
            
            results.append({
                'sample': sample,
                'centype': centype,
                'gc_content': gc_content,
                'total_bases': total_bases,
                'num_arrays': len(seq_list)
            })
        
        return results
    
    # Use multiprocessing
    total_tasks = len(task_list)
    # Store GC values for each centype (one value per sample, not per array)
    results_gc = {ct: [] for ct in centype_list}
    # Store detailed results
    results_detailed = {ct: [] for ct in centype_list}
    
    print(f"Starting multiprocessing with 20 workers...")
    
    with Pool(processes=20) as pool:
        for i, sample_results in enumerate(pool.imap(process_sample, task_list), start=1):
            # Process results returned by each sample
            for result in sample_results:
                centype = result['centype']
                results_gc[centype].append(result['gc_content'])
                results_detailed[centype].append(result)
            
            # Display progress
            progress = (i / total_tasks) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}% ({i}/{total_tasks})")
            sys.stdout.flush()
    
    print("\n\nAll samples processed!")
    
    # Output statistics to file (one value per centype per sample)
    with open("./stat_plot/step3_chr_group/centype_gc_content.txt", 'w') as f_out:
        f_out.write("centype\tmean_gc_percent\tsd\tn_samples\tmin\tmax\tall_gc_values\n")
        for ct in centype_list:
            gc_values = results_gc[ct]
            if gc_values:
                mean_gc = sum(gc_values) / len(gc_values)
                sd_gc = (sum((x - mean_gc)**2 for x in gc_values) / len(gc_values))**0.5
                min_gc = min(gc_values)
                max_gc = max(gc_values)
                all_values_str = ",".join([f"{x:.2f}" for x in gc_values])
                f_out.write(f"{ct}\t{mean_gc:.2f}\t{sd_gc:.2f}\t{len(gc_values)}\t{min_gc:.2f}\t{max_gc:.2f}\t{all_values_str}\n")
            else:
                f_out.write(f"{ct}\tNA\tNA\t0\tNA\tNA\tNA\n")
    
    # Output detailed results (one row per sample per centype)
    with open("./stat_plot/step3_chr_group/centype_gc_content_detailed.txt", 'w') as f_detail:
        f_detail.write("centype\tsample\tgc_percent\ttotal_bases\tnum_arrays\n")
        for ct in centype_list:
            for item in results_detailed[ct]:
                f_detail.write(f"{ct}\t{item['sample']}\t{item['gc_content']:.2f}\t{item['total_bases']}\t{item['num_arrays']}\n")
    
    # Print summary to console
    print("\n=== GC content summary per centype (merged sequences, then calculate) ===")
    for ct in centype_list:
        gc_values = results_gc[ct]
        if gc_values:
            mean_gc = sum(gc_values) / len(gc_values)
            print(f"{ct}: {mean_gc:.2f}% GC (n={len(gc_values)} samples)")
        else:
            print(f"{ct}: no data")
    
    print("\nOutput saved to:")
    print("  - ./stat_plot/step3_chr_group/centype_gc_content.txt (summary)")
    print("  - ./stat_plot/step3_chr_group/centype_gc_content_detailed.txt (per sample)")
if argv1=="stepall" or argv1=="step3" or argv1=="step3.0p":    
    Plot_txt=r"""
library(ggplot2)
library(dplyr)
# Get all command line arguments

### Monomer
print("")
{
  inputfile=read.table('inputfile', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')





  filtered_input <- inputfile %>% 
    filter(group %in% c("Eurasian","East_Asia","America")) 
    
    # Define facet order
list_A <- c("Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap1","V059.hap2","V063.hap1","V063.hap2","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap1","Chardonnay_hap2","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap1","V074.hap2","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap1","V062.hap2","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap1","PinotNoir_hap2","PinotNoir2_hap1","PinotNoir2_hap2","V069.hap1","V069.hap2","PN40024","PN40024_hap1","PN40024_hap2","V087.hap1","V087.hap2","V061.hap1","V061.hap2","V065.hap1","V065.hap2","V064.hap1","V064.hap2","ThompsonSeedless_hap1","ThompsonSeedless_hap2","ThompsonSeedless2_hap1","ThompsonSeedless2_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap1","V067.hap2","V058.hap1","V058.hap2","V107.hap1","V107.hap2","V108.hap1","V108.hap2","V081.hap1","V081.hap2","V015.hap1","V015.hap2","V105.hap1","V105.hap2","V126.hap1","V126.hap2","V112.hap1","V112.hap2","DavidiiGrape_hap1","DavidiiGrape_hap2","V007.hap1","V007.hap2","V008.hap1","V008.hap2","V012.hap1","V012.hap2","V100.hap1","V100.hap2","V099.hap1","V099.hap2","V098.hap1","V098.hap2","V117.hap1","V117.hap2","V120.hap1","V120.hap2","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap1","V124.hap2","V125.hap1","V125.hap2","WoollyGrape_hap1","WoollyGrape_hap2","V106.hap1","V106.hap2","V102.hap1","V102.hap2","V079.hap1","V079.hap2","V030.hap1","V030.hap2","V031.hap1","V031.hap2","V005.hap1","V005.hap2","V048.hap1","V048.hap2","V050.hap1","V050.hap2","V051.hap1","V051.hap2","V049.hap1","V049.hap2","V038.hap1","V038.hap2","V043.hap1","V043.hap2","V052.hap1","V052.hap2","V034.hap1","V034.hap2","V023.hap1","V023.hap2","V032.hap1","V032.hap2","V033.hap1","V033.hap2","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap1","V022.hap2","V040.hap1","V040.hap2","V019.hap1","V019.hap2","V041.hap1","V041.hap2","V077.hap1","V077.hap2","V037.hap1","V037.hap2","V096.hap1","V096.hap2","ShineMuscat_hap1","ShineMuscat_hap2","V036.hap1","V036.hap2","V076.hap1","V076.hap2","V055.hap1","V055.hap2","V018.hap1","V018.hap2","V072.hap1","V072.hap2")

list_chr <- c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19")

 # Convert chromosome column to factor, order according to list_chr
filtered_input$sample <- factor(filtered_input$sample, levels = list_A) 
filtered_input$chromosome <- factor(filtered_input$chromosome, levels = list_chr)
filtered_input$group <- factor(filtered_input$group, levels = c("Eurasian","East_Asia","America"))

  # Define color values
  color_values <- c(
    "Eurasian"="#ff3399",
    "East_Asia"="#0066ff",
    "America"="#33cc33",
    "berlandieri_riparia_labrusca"="#e4d5ab",
    "betulifolia"="#74a0c6",
    "bryoniifolia var. bryoniifolia"="#cc0066",
    "cinerea_riparia"="#148EFF",
    "cinerea_rupestris"="#B2DF8A",
    "cinerea_vinifera"="#008000",
    "davidii"="#FB9A99",
    "heyneana"="#00C08F",
    "labrusca_berlandieri"="#005FAF",
    "labrusca_riparia"="#FFCC00",
    "piasezkii"="#00ccff",
    "pseudoreticulata"="#6A3D9A",
    "retordii"="#FFFF00",
    "riparia"="#A0522D",
    "riparia_rupestris"="#9966cc",
    "rupestris"="#ba03d9",
    "sinocinerea"="#7878ff",
    "vinifera"="#0073CF",
    "wilsoniae"="#FF8C00",
    "wuhanensis"="#ff80df",
    "yunnanensis"="#cc4100"
     #   TRUE ~ "Other"
  )
      
  # Create plot object
  p = ggplot() 
# 1. Add boxplot: x=group, y=length, fill color by group
  p = p+geom_boxplot(
    data = filtered_input,
    aes(x = group, y = length,color=group),
    width = 0.7,linewidth=0.1, fill = 'white',       # Box width
    outlier.shape = NA # Hide outliers in boxplot (avoid duplication with scatter points)
  ) 
  # 2. Add scatter plot: x=group, y=length, fill color by species
  p = p+geom_point(
    data = filtered_input,
    aes(x = group, y = length, color = species),
    position = position_jitter(width = 0.2), # Slight jitter to avoid overlap
    size = 0.5         # Point size
    
  ) 
p <- p + facet_wrap(~ chromosome, ncol = 20)


   p=p+theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      #legend.position = "none",
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      strip.text = element_blank()            # Hide facet titles
    ) +
scale_color_manual(values = color_values, drop = FALSE)
  

  # Save as PDF
  pdf(file = paste0('output', ".pdf"), width = 40 / 2.54, height = 7 / 2.54)
  print(p)
  dev.off()
  
}            
        """
    with open(f'./stat_plot/step3_chr_group/plot.R','w') as f:
        f.write(Plot_txt)
    os.chdir(f'./stat_plot/step3_chr_group/')    
    centype_list= ["cen_107","cen_66","cen_103","cen_191","cen_383","cen_355","LIRSat1","LIRSat2","LIRSat3"]
    for one_centype in centype_list:    
        print(one_centype)
        
        
        subprocess.run([f'cp {one_centype} inputfile'], shell=True)  #
        subprocess.run([f'Rscript plot.R  '], shell=True)  #>null 2>&1 
        #sys.exit()
        subprocess.run([f'mv output.pdf {one_centype}.pdf  '], shell=True) 
        
    os.chdir('../../../')                            
    
        
if argv1=="stepall" or argv1=="step2" or argv1=="step4":    
    subprocess.run(["mkdir ./stat_plot/step4_hap_stat"], shell=True)
    print('Determine whether LIR1 and VSat6(cen_383) can be on the same chromosome; LIR1 and VSat6 have high similarity')
    dict_sample_chromosome_info={}
    centype_list_used=['cen_383','LIRSat1']
    sample_list= [
"Baimunage_hap1",
"Baimunage_hap2",
"BlackMonukka_hap1",
"BlackMonukka_hap2",
"V059.hap1",
"V059.hap2",
"V063.hap1",
"V063.hap2",
"V066.hap1",
"V066.hap2",
"V070.hap1",
"V070.hap2",
"Chardonnay_hap1",
"Chardonnay_hap2",
"VHP-T2T.hap1",
"VHP-T2T.hap2",
"V091.hap1",
"V091.hap2",
"Hongmunage_hap1",
"Hongmunage_hap2",
"V074.hap1",
"V074.hap2",
"ManicureFinger_hap1",
"ManicureFinger_hap2",
"V092.hap1",
"V092.hap2",
"V062.hap1",
"V062.hap2",
"MuscatHamburg_hap1",
"MuscatHamburg_hap2",
"V088.hap1",
"V088.hap2",
"PinotNoir_hap1",
"PinotNoir_hap2",
"PinotNoir2_hap1",
"PinotNoir2_hap2",
"V069.hap1",
"V069.hap2",
"PN40024",
"PN40024_hap1",
"PN40024_hap2",
"V087.hap1",
"V087.hap2",
"V061.hap1",
"V061.hap2",
"V065.hap1",
"V065.hap2",
"V064.hap1",
"V064.hap2",
"ThompsonSeedless_hap1",
"ThompsonSeedless_hap2",
"ThompsonSeedless2_hap1",
"ThompsonSeedless2_hap2",
"V093.hap1",
"V093.hap2",
"V060.hap1",
"V060.hap2",
"V067.hap1",
"V067.hap2",
"V058.hap1",
"V058.hap2",
"V107.hap1",
"V107.hap2",
"V108.hap1",
"V108.hap2",
"V081.hap1",
"V081.hap2",
"V015.hap1",
"V015.hap2",
"V105.hap1",
"V105.hap2",
"V126.hap1",
"V126.hap2",
"V112.hap1",
"V112.hap2",
"DavidiiGrape_hap1",
"DavidiiGrape_hap2",
"V007.hap1",
"V007.hap2",
"V008.hap1",
"V008.hap2",
"V012.hap1",
"V012.hap2",
"V100.hap1",
"V100.hap2",
"V099.hap1",
"V099.hap2",
"V098.hap1",
"V098.hap2",
"V117.hap1",
"V117.hap2",
"V120.hap1",
"V120.hap2",
"PiasezkiiGrape_hap1",
"PiasezkiiGrape_hap2",
"V123.hap1",
"V123.hap2",
"V124.hap1",
"V124.hap2",
"V125.hap1",
"V125.hap2",
"WoollyGrape_hap1",
"WoollyGrape_hap2",
"V106.hap1",
"V106.hap2",
"V102.hap1",
"V102.hap2",
"V079.hap1",
"V079.hap2",
"V030.hap1",
"V030.hap2",
"V031.hap1",
"V031.hap2",
"V005.hap1",
"V005.hap2",
"V048.hap1",
"V048.hap2",
"V050.hap1",
"V050.hap2",
"V051.hap1",
"V051.hap2",
"V049.hap1",
"V049.hap2",
"V038.hap1",
"V038.hap2",
"V043.hap1",
"V043.hap2",
"V052.hap1",
"V052.hap2",
"V034.hap1",
"V034.hap2",
"V023.hap1",
"V023.hap2",
"V032.hap1",
"V032.hap2",
"V033.hap1",
"V033.hap2",
"V053.hap1",
"V053.hap2",
"V020.hap1",
"V020.hap2",
"V022.hap1",
"V022.hap2",
"V040.hap1",
"V040.hap2",
"V019.hap1",
"V019.hap2",
"V041.hap1",
"V041.hap2",
"V077.hap1",
"V077.hap2",
"V037.hap1",
"V037.hap2",
"V096.hap1",
"V096.hap2",
"ShineMuscat_hap1",
"ShineMuscat_hap2",
"V036.hap1",
"V036.hap2",
"V076.hap1",
"V076.hap2",
"V055.hap1",
"V055.hap2",
"V018.hap1",
"V018.hap2",
"V072.hap1",
"V072.hap2"]
    with open('./stat_plot/0-region2info','r') as f: 
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            sample= eachline_arr[0]	
            chromosome_new=eachline_arr[2]	
            centype=eachline_arr[3]	
            if centype not in centype_list_used:continue
            if sample not in sample_list:continue
            if sample not in dict_sample_chromosome_info:dict_sample_chromosome_info[sample]={}
            if chromosome_new not in dict_sample_chromosome_info[sample]:   
                dict_sample_chromosome_info[sample][chromosome_new]={}
                for one_centype in centype_list_used:
                    dict_sample_chromosome_info[sample][chromosome_new][one_centype]=0
            dict_sample_chromosome_info[sample][chromosome_new][centype]+=1
    with open('./stat_plot/step4_hap_stat/0_info','w') as f: 
        head_centype_str='\t'.join(centype_list_used)
        f.write(f'sample\tchromosome\t{head_centype_str}\n')
        for sample,dict_chromosome_info in dict_sample_chromosome_info.items():
            for chromosome,dict_info in dict_chromosome_info.items():
                nums_list=[]
                for one_centype in centype_list_used:
                    one_num=dict_info[one_centype]
                    nums_list.append(f"{one_num}")
                nums_list_str= '\t'.join(nums_list)   
                one_line=f"{sample}\t{chromosome}\t{nums_list_str}\n"
                f.write(f'{one_line}')
    
if argv1=="stepall" or argv1=="step2" or argv1=="step5":
    subprocess.run(["mkdir ./stat_plot/step5_VSat6_LIRSat2"], shell=True)
    print('Calculate the proportion of LIRSat2 near VSat6')    
    merge_distance=10000
    all2_str=''  
    ALL_STR=''
    flank_distance_list=[3000,10000,15000,20000,25000,30000,35000,40000,45000,50000]
    for flank_distance in flank_distance_list:
        #flank_distance=100000
        centype_list_used=['cen_383','LIRSat1','LIRSat2','LIRSat3']
        sample_list= [
    "Baimunage_hap1",
    "Baimunage_hap2",
    "BlackMonukka_hap1",
    "BlackMonukka_hap2",
    "V059.hap1",
    "V059.hap2",
    "V063.hap1",
    "V063.hap2",
    "V066.hap1",
    "V066.hap2",
    "V070.hap1",
    "V070.hap2",
    "Chardonnay_hap1",
    "Chardonnay_hap2",
    "VHP-T2T.hap1",
    "VHP-T2T.hap2",
    "V091.hap1",
    "V091.hap2",
    "Hongmunage_hap1",
    "Hongmunage_hap2",
    "V074.hap1",
    "V074.hap2",
    "ManicureFinger_hap1",
    "ManicureFinger_hap2",
    "V092.hap1",
    "V092.hap2",
    "V062.hap1",
    "V062.hap2",
    "MuscatHamburg_hap1",
    "MuscatHamburg_hap2",
    "V088.hap1",
    "V088.hap2",
    "PinotNoir_hap1",
    "PinotNoir_hap2",
    "PinotNoir2_hap1",
    "PinotNoir2_hap2",
    "V069.hap1",
    "V069.hap2",
    "PN40024",
    "PN40024_hap1",
    "PN40024_hap2",
    "V087.hap1",
    "V087.hap2",
    "V061.hap1",
    "V061.hap2",
    "V065.hap1",
    "V065.hap2",
    "V064.hap1",
    "V064.hap2",
    "ThompsonSeedless_hap1",
    "ThompsonSeedless_hap2",
    "ThompsonSeedless2_hap1",
    "ThompsonSeedless2_hap2",
    "V093.hap1",
    "V093.hap2",
    "V060.hap1",
    "V060.hap2",
    "V067.hap1",
    "V067.hap2",
    "V058.hap1",
    "V058.hap2",
    "V107.hap1",
    "V107.hap2",
    "V108.hap1",
    "V108.hap2",
    "V081.hap1",
    "V081.hap2",
    "V015.hap1",
    "V015.hap2",
    "V105.hap1",
    "V105.hap2",
    "V126.hap1",
    "V126.hap2",
    "V112.hap1",
    "V112.hap2",
    "DavidiiGrape_hap1",
    "DavidiiGrape_hap2",
    "V007.hap1",
    "V007.hap2",
    "V008.hap1",
    "V008.hap2",
    "V012.hap1",
    "V012.hap2",
    "V100.hap1",
    "V100.hap2",
    "V099.hap1",
    "V099.hap2",
    "V098.hap1",
    "V098.hap2",
    "V117.hap1",
    "V117.hap2",
    "V120.hap1",
    "V120.hap2",
    "PiasezkiiGrape_hap1",
    "PiasezkiiGrape_hap2",
    "V123.hap1",
    "V123.hap2",
    "V124.hap1",
    "V124.hap2",
    "V125.hap1",
    "V125.hap2",
    "WoollyGrape_hap1",
    "WoollyGrape_hap2",
    "V106.hap1",
    "V106.hap2",
    "V102.hap1",
    "V102.hap2",
    "V079.hap1",
    "V079.hap2",
    "V030.hap1",
    "V030.hap2",
    "V031.hap1",
    "V031.hap2",
    "V005.hap1",
    "V005.hap2",
    "V048.hap1",
    "V048.hap2",
    "V050.hap1",
    "V050.hap2",
    "V051.hap1",
    "V051.hap2",
    "V049.hap1",
    "V049.hap2",
    "V038.hap1",
    "V038.hap2",
    "V043.hap1",
    "V043.hap2",
    "V052.hap1",
    "V052.hap2",
    "V034.hap1",
    "V034.hap2",
    "V023.hap1",
    "V023.hap2",
    "V032.hap1",
    "V032.hap2",
    "V033.hap1",
    "V033.hap2",
    "V053.hap1",
    "V053.hap2",
    "V020.hap1",
    "V020.hap2",
    "V022.hap1",
    "V022.hap2",
    "V040.hap1",
    "V040.hap2",
    "V019.hap1",
    "V019.hap2",
    "V041.hap1",
    "V041.hap2",
    "V077.hap1",
    "V077.hap2",
    "V037.hap1",
    "V037.hap2",
    "V096.hap1",
    "V096.hap2",
    "ShineMuscat_hap1",
    "ShineMuscat_hap2",
    "V036.hap1",
    "V036.hap2",
    "V076.hap1",
    "V076.hap2",
    "V055.hap1",
    "V055.hap2",
    "V018.hap1",
    "V018.hap2",
    "V072.hap1",
    "V072.hap2"]
        dict_sample_chromosome_info={}
        with open('./stat_plot/0-region2info','r') as f: 
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                sample= eachline_arr[0]	
                chromosome_new=eachline_arr[2]	
                centype=eachline_arr[3]	
                start=int(eachline_arr[4])
                end=int(eachline_arr[5])
                if centype not in centype_list_used:continue
                if sample not in sample_list:continue
                if sample not in dict_sample_chromosome_info:dict_sample_chromosome_info[sample]={}
                if chromosome_new not in dict_sample_chromosome_info[sample]:   
                    dict_sample_chromosome_info[sample][chromosome_new]={}
                    for one_centype in centype_list_used:
                        dict_sample_chromosome_info[sample][chromosome_new][one_centype]=[]
                dict_sample_chromosome_info[sample][chromosome_new][centype].append([start,end])        
        #print(dict_sample_chromosome_info['PN40024']['Chr15']['cen_383'])
        ##To merge intervals in list A_list, with merge condition that the distance between the end of adjacent intervals and the start is less than 10,000, follow the steps below.                
        def merge_intervals(A_list):
            if not A_list:
                return []
            
            # Sort by interval start in ascending order
            A_list.sort(key=lambda x: x[0])
            
            merged = [A_list[0]]
            for current in A_list[1:]:
                last = merged[-1]
                # Check if the current interval's start distance from the previous interval's end is < 10000
                if current[0] - last[1] < merge_distance:
                    # Merge: update the last interval's end to the current interval's end (or other logic)
                    # Here assume merging takes the current interval's end (can be adjusted as needed)
                    merged[-1] = [last[0], max(last[1], current[1])]
                else:
                    merged.append(current)
            return merged                
        
        for sample,dict_chr_info in dict_sample_chromosome_info.items():
            for chromosome,dict_info in dict_chr_info.items():
                for one_centype,A_list in dict_info.items():
                    new_A_list=merge_intervals(A_list)
                    dict_sample_chromosome_info[sample][chromosome][one_centype]=new_A_list
        
        all_num=0 
        yes1_num=0
        yes2_num=0
        yes3_num=0
        min_distance_list=[]
        for sample,dict_chr_info in dict_sample_chromosome_info.items():
            for chromosome,dict_info in dict_chr_info.items():
                cen_383_list=dict_info['cen_383']
                #print(cen_383_list)       
                LIRSat1_list=dict_info['LIRSat1']
                LIRSat2_list=dict_info['LIRSat2']
                LIRSat3_list=dict_info['LIRSat3']
                
                for one_cen_383 in cen_383_list:
                    all_num+=1
                    ##################
                    mark=''
                    for one_LIRSat1 in LIRSat1_list:
                        if not (one_LIRSat1[1] < one_cen_383[0] - flank_distance or one_LIRSat1[0] > one_cen_383[1] + flank_distance): mark='yes';break
                    if mark=='yes':
                        yes1_num+=1
                    ##################
                    mark=''
                    for one_LIRSat2 in LIRSat2_list:
                        if not (one_LIRSat2[1] < one_cen_383[0] - flank_distance or one_LIRSat2[0] > one_cen_383[1] + flank_distance): mark='yes';break
                    if mark=='yes':
                        yes2_num+=1
                        min_distance=min(abs(one_LIRSat2[0]-one_cen_383[0]),abs(one_LIRSat2[0]-one_cen_383[1]),abs(one_LIRSat2[1]-one_cen_383[0]),abs(one_LIRSat2[1]-one_cen_383[1]))
                        min_distance_list.append(min_distance)
                    ##################
                    mark=''
                    for one_LIRSat3 in LIRSat3_list:
                        if not (one_LIRSat3[1] < one_cen_383[0] - flank_distance or one_LIRSat3[0] > one_cen_383[1] + flank_distance): mark='yes';break
                    if mark=='yes':
                        yes3_num+=1                    
        with open('./stat_plot/step5_VSat6_LIRSat2/0_info_mindistance','w')  as f:
            for one in min_distance_list:
                f.write(f'{one}\n')

        with open('./stat_plot/step5_VSat6_LIRSat2/0_info','w')     as f:
            ALL_STR+=(f"flank_distance\t{flank_distance}\n")
            ALL_STR+=(f"Total number of VSat6 arrays (merged if gap <{merge_distance}bp)\t{all_num}\n\n\n")
            ALL_STR+=(f"Number of VSat6 arrays with LIRSat1 within {flank_distance}\t{yes1_num}\n")
            ALL_STR+=(f"Number of VSat6 arrays with LIRSat2 within {flank_distance}\t{yes2_num}\n");            all2_str+=(f"Number of VSat6 arrays with LIRSat2 within {flank_distance}\t{yes2_num}\n")
            ALL_STR+=(f"Number of VSat6 arrays with LIRSat3 within {flank_distance}\t{yes3_num}\n")
            f.write(f"{ALL_STR}\n")
            #print(ALL_STR)
        
        
    print(all2_str)
    

print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))      
time_end=timeit.default_timer()
print('All the running time: %.0f Seconds'%(time_end-time_start))