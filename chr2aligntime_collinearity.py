#!/usr/bin/python
# -*- coding: UTF-8 -*-
import sys
argvs=sys.argv
argv_len=len(sys.argv)
if argv_len < 2 or argvs[1] in ("h", "-h", "help", "-help") :
    print ("chr2aligntime.py-----help:")
    print ("")
    print ("Usage：")
    print ("chr2aligntime_collinearity.py  stepall Manually set, display collinearity of multiple samples")
    print (" ")
    print ("chr2aligntime_collinearity.py  step1 Adjacent pairwise alignment")
    print ("chr2aligntime_collinearity.py  step2 Plotting")
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
import re # Handle regular expressions
from multiprocessing import Pool, cpu_count
import numpy as np
time_start=timeit.default_timer()
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  

if  os.path.exists('./chr2aligntime_collinearity')==False:
    subprocess.run(["mkdir ./chr2aligntime_collinearity"], shell=True)


target_list=[]
'''
target_list.append(["V032.hap2","Chr9",16059125,16258850])
target_list.append(["V036.hap1","Chr9",15852836,16122334])
target_list.append(["V022.hap1","Chr9",12892149,13157886])
target_list.append(["V037.hap1","Chr9",15841557,16135972])
target_list.append(["V051.hap2","Chr9",16109432,16366204])
target_list.append(["V034.hap2","Chr9",16082980,16346901])
target_list.append(["V048.hap1","Chr9",15846112,16075387])

target_list.append(["V061.hap1","Chr9",16735650,17078790])
target_list.append(["V059.hap2","Chr9",16638477,16981625])
target_list.append(["V060.hap2","Chr9",16993605,17360565])

target_list.append(["Hongmunage_hap1","Chr9",17428070,17789105])
target_list.append(["MuscatHamburg_hap1","Chr9",15304624,15653041])
target_list.append(["BlackMonukka_hap2","Chr9",17599949,17960992])
target_list.append(["PinotNoir_hap1","Chr9",17354494,17721462])
target_list.append(["V087.hap1","Chr9",15060064,15408481])

target_list.append(["V108.hap1","Chr9",18211370,18637467])
target_list.append(["V100.hap2","Chr9",17496278,17854348])
target_list.append(["V079.hap2","Chr9",18687495,19097694])
target_list.append(["DavidiiGrape_hap2","Chr9",20415421,20890118])
target_list.append(["WoollyGrape_hap1","Chr9",17693900,17872344])
target_list.append(["V106.hap1","Chr9",19205595,19389522])
target_list.append(["V124.hap2","Chr9",19043760,19198608])
'''
kmer_len=41  ##kmer41 for LIR, kmer200 for VSat1
target_list.append(['PN40024','Chr1',15154005,15435737])
#target_list.append(['PN40024','Chr13',12020240,12180415])
target_list.append(['PN40024','Chr7',13458101,13762988]) 
target_list.append(['PN40024','Chr14',15144034,15480228])
target_list.append(['PN40024','Chr17',15221832,15560758])
target_list.append(['PN40024','Chr2',12939361,13300832])
#target_list.append(['PN40024_hap1','Chr3',14353397,14553397])
target_list.append(['PN40024','Chr4',13291174,13483935])
target_list.append(['PN40024','Chr5',14008287,14357682])
target_list.append(['PN40024','Chr6',10606684,10966038])
target_list.append(['PN40024','Chr8',6679384,7071344])
target_list.append(['PN40024','Chr9',15113362,15449399])
#target_list.append(['PN40024','Chr19',15995936,16195936])
target_list.append(['PN40024','Chr18',16390000,16620000])



##Attached figure, showing homology on both sides of LIR
target_list=[]
kmer_len=201 ##kmer41 for LIR, kmer200 for VSat1
target_list.append(['V077.hap1','Chr13',13155000,13405000])
target_list.append(['V077.hap1','Chr13',13658000,13845000])
##
target_list=[]
kmer_len=201
target_list.append(['V037.hap2','Chr13',13160000,13410000])
target_list.append(['V037.hap2','Chr13',13660000,13835000])
##
target_list=[]
kmer_len=301
target_list.append(['V079.hap1','Chr19',18890000,19135000])
target_list.append(['V079.hap1','Chr19',19420000,19655000])
##
target_list=[]
kmer_len=301
target_list.append(['V105.hap1','Chr1',16560000,17230000])
target_list.append(['V105.hap1','Chr1',17370000,18010000])

##Region near VSat6
target_list=[]
kmer_len=21
target_list.append(['PN40024','Chr10',20040000,20190000])
target_list.append(['PN40024','Chr11',13170000,13320000])
target_list.append(['PN40024','Chr15',8950000,9100000]) 
target_list.append(['PN40024','Chr16',12720000,12870000])      ########These have no results, alignment too short to recognize


target_num=len(target_list)
folder_name=f"{target_list[0][0]}:{target_list[0][1]}.etal___FileNum{target_num}"
if  os.path.exists(f'./chr2aligntime_collinearity/{folder_name}')==False:
    subprocess.run([f"mkdir ./chr2aligntime_collinearity/{folder_name}"], shell=True)
if argv1=='stepall' or argv1=='step1':
    print('Process pairwise comparisons separately, identify highly similar segments')
    if  os.path.exists('./chr2aligntime_collinearity/{folder_name}/0_vsall')==False:
        subprocess.run([f'mkdir ./chr2aligntime_collinearity/{folder_name}/0_vsall'], shell=True) 
    ###

    i=0
    while i< target_num-1:
        one_target1=target_list[i];                 sample1,chr1,start1,end1=one_target1
        one_target2=target_list[i+1];               sample2,chr2,start2,end2=one_target2
        print(f'Progress: {i}/{target_num-1}')
        one_target1_str=f"{sample1}:{chr1}:{start1}-{end1}"
        one_target2_str=f"{sample2}:{chr2}:{start2}-{end2}"
        
        i+=1
        subprocess.run([f'mkdir ./chr2aligntime_collinearity/{folder_name}/0_vsall/{i}'], shell=True) 
        if os.path.exists(f'./chr2aligntime_collinearity/{folder_name}/0_vsall/{i}/kmer_matches.up_lines4')==False:
            cmd=f'python ./chr2aligntime.py {one_target1_str} {one_target2_str} ./chr2aligntime_collinearity/{folder_name}/0_vsall/{i} {kmer_len}'
            print(cmd)
            subprocess.run([cmd], shell=True)  #> ./chr2aligntime_collinearity/0_vsall/{i}/log
    
if argv1=='stepall' or argv1=='step2':
    print('Summarize')
    if  os.path.exists(f'./chr2aligntime_collinearity/{folder_name}/1_sum')==False:
        subprocess.run([f'mkdir ./chr2aligntime_collinearity/{folder_name}/1_sum'], shell=True) 
    i=1;kk=0
    with open(f'./chr2aligntime_collinearity/{folder_name}/1_sum/0_all','w') as f2:
        f2.write(f"index\tindex2\tstart_x\tstart_y\tend_x\tend_y\tpercent\tTime\tplot_y\n")
        while i< target_num: 
            with open(f'./chr2aligntime_collinearity/{folder_name}/0_vsall/{i}/kmer_matches.up_lines4','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    print(eachline_arr)
                    index,start_x,start_y,end_x,end_y,align_length,var_num,percent,Time=eachline_arr
                    kk+=1
                    f2.write(f"{i}-{i+1}\t{kk}\t{start_x}\t{start_y}\t{end_x}\t{end_y}\t{percent}\t{Time}\t{-i}\n")
                    
            i+=1
    #######        
    i=0
    with open(f'./chr2aligntime_collinearity/{folder_name}/1_sum/0_base','w') as f:
        f.write(f"index\tsample\tchr\tstart\tend\n")
        for one in target_list:
            i+=1
            sample1,chr1,start1,end1=one
            f.write(f"{i}\t{sample1}\t{chr1}\t{start1}\t{end1}\n")
            
    ###                
    print('Plotting')
    R_txt=r"""
library(ggplot2)
library(dplyr)
library(tidyr)

input_file=read.table('0_all', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
input_file2=read.table('0_base', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
# Calculate the 4 vertices of the trapezoid
plot_data <- input_file %>%
  rowwise() %>%
  mutate(
    x = list(c(start_x,end_x, end_y, start_y)),
    y = list(c(plot_y + 0.3, plot_y + 0.3, plot_y - 0.3, plot_y - 0.3)),
    group = index2
  ) %>%
  unnest(c(x, y)) 
p=ggplot()
p=p+ geom_rect(data=input_file2, aes(xmin =0 , xmax = end-start, ymin = -index+ 0.45, ymax = -index+ 0.55), fill = 'black')
p=p+ geom_polygon(data=plot_data,aes(x = x,y=y, group = group,fill=percent))

p=p+theme_classic()+
    #scale_color_manual(values = color_values, drop = FALSE)+
            theme(
        #axis.ticks.y = element_blank(),
        #axis.text.y = element_blank(),
        #legend.position = "none",
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank()
        #panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Manually add border
    )

# Save as PDF
pdf(file = paste0('1_all', ".pdf"), width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()        
"""
    with open(f'./chr2aligntime_collinearity/{folder_name}/1_sum/plot.R','w') as f:
        f.write(R_txt)
    os.chdir(f'./chr2aligntime_collinearity/{folder_name}/1_sum/')
    subprocess.run([f'Rscript plot.R  '], shell=True)  #>n
    subprocess.run([f'mv 1_all.pdf {folder_name}.1_all.pdf'], shell=True) 
    os.chdir('../../../')    


    print('Plotting 2')
    R_txt=r"""
library(ggplot2)
library(dplyr)
library(tidyr)

input_file=read.table('0_all', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
input_file2=read.table('0_base', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')

# Calculate the 4 vertices of the trapezoid
plot_data <- input_file %>%
  rowwise() %>%
  mutate(
    x = list(c(start_x,end_x, end_y, start_y)),
    y = list(c(plot_y + 0.3, plot_y + 0.3, plot_y - 0.3, plot_y - 0.3)),
    group = index2
  ) %>%
  mutate(
        category = case_when(
            percent<0.01 ~ "<1",        #2 million years
            percent<=0.02  ~ "1-2",      #4 million years
            percent<=0.03  ~ "2-3",      #6 million years
            percent<=0.04  ~ "3-4",     #8 million years
            percent<=0.08  ~ "4-8",    #16 million years
            percent>=0.08  ~ "8+",   #16 million years+
        ),
      )   %>%
  unnest(c(x, y))
# Define color values
color_values <- c(
        "<1" = "#e03c31",
          "1-2" = "#ff7f41",
          "2-3" = "#f7ea48",
          "3-4" = "#2dc84d",
          "4-8" = "#147bd1",
          "8+" = "#753bbd"
)    
  
p=ggplot()
p=p+ geom_rect(data=input_file2, aes(xmin =0 , xmax = end-start, ymin = -index+ 0.45, ymax = -index+ 0.55), fill = 'black')
p=p+ geom_polygon(data=plot_data,aes(x = x,y=y, group = group,fill=category),alpha=0.5)+scale_fill_manual(values = color_values, drop = FALSE)
  
p=p+theme_classic()+
            theme(
        #axis.ticks.y = element_blank(),
        #axis.text.y = element_blank(),
        #legend.position = "none",
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank()
        #panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Manually add border
    )

# Save as PDF
pdf(file = paste0('2_all', ".pdf"), width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()        
"""
    with open(f'./chr2aligntime_collinearity/{folder_name}/1_sum/plot.R','w') as f:
        f.write(R_txt)
    os.chdir(f'./chr2aligntime_collinearity/{folder_name}/1_sum/')
    subprocess.run([f'Rscript plot.R  '], shell=True)  #>n
    subprocess.run([f'mv 2_all.pdf {folder_name}.2_all.pdf'], shell=True) 
    os.chdir('../../../')   


    print('Plotting 3')
    R_txt=r"""
library(ggplot2)
library(dplyr)
library(tidyr)

input_file=read.table('0_all', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
input_file2=read.table('0_base', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')



# Calculate the 4 vertices of the trapezoid
plot_data <- input_file %>%
  rowwise() %>%
  mutate(
    x = list(c(start_x,end_x, end_y, start_y)),
    y = list(c(plot_y + 0.3, plot_y + 0.3, plot_y - 0.3, plot_y - 0.3)),
    group = index2
  ) %>%
  mutate(
        category = case_when(
            percent<0.01 ~ "<1",        #2 million years
            percent<=0.02  ~ "1-2",      #4 million years
            percent<=0.03  ~ "2-3",      #6 million years
            percent<=0.04  ~ "3-4",     #8 million years
            percent<=0.08  ~ "4-8",    #16 million years
            percent>=0.08  ~ "8+",   #16 million years+
        ),
      )   %>%
  unnest(c(x, y))
# Define color values
color_values <- c(
        "<1" = "#FF0066",
          "1-2" = "#BF008C",
          "2-3" = "#7F00B2",
          "3-4" = "#3F00D8",
          "4-8" = "#0000FF",
          "8+" = "#0000FF"
)    
  
p=ggplot()
p=p+ geom_rect(data=input_file2, aes(xmin =0 , xmax = end-start, ymin = -index+ 0.45, ymax = -index+ 0.55), fill = 'black')
p=p+ geom_polygon(data=plot_data,aes(x = x,y=y, group = group,fill=category),alpha=0.5)+scale_fill_manual(values = color_values, drop = FALSE)
  
p=p+theme_classic()+
            theme(
        #axis.ticks.y = element_blank(),
        #axis.text.y = element_blank(),
        #legend.position = "none",
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank()
        #panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Manually add border
    )

# Save as PDF
pdf(file = paste0('3_all', ".pdf"), width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()        
"""
    with open(f'./chr2aligntime_collinearity/{folder_name}/1_sum/plot.R','w') as f:
        f.write(R_txt)
    os.chdir(f'./chr2aligntime_collinearity/{folder_name}/1_sum/')
    subprocess.run([f'Rscript plot.R  '], shell=True)  #>n
    subprocess.run([f'mv 3_all.pdf {folder_name}.3_all.pdf'], shell=True) 
    os.chdir('../../../')   



























print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))      
time_end=timeit.default_timer()
print('All the running time: %.0f Seconds'%(time_end-time_start))