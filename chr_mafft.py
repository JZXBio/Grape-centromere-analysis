#!/usr/bin/python
# -*- coding: UTF-8 -*-
import sys
argvs=sys.argv
argv_len=len(sys.argv)
if argv_len < 2 or argvs[1] in ("h", "-h", "help", "-help") or not argvs[1].startswith("step") :
    print ("chr_mafft.py-----help:")
    print ("chr_mafft.py is used to identify monomers or subunits. It needs to be placed in the input_folder/sample_name/ directory, with the sample_name.fa file in the same location.")
    print ("\nUsage：")
    print ("chr_mafft.py stepall -i FASTA_file")
    print ("")
    print ("Monomer sequence can be specified using -monomer or modified directly in this script. Default monomer is for grapevine VSat1 (28bp)")
    print ("")
    print ("")
    print ("You can also execute intermediate steps only：")
    print ("chr_mafft.py step0 -i FASTA_file   ##Split sequences by chromosome")
    print ("chr_mafft.py step1/step1r          ##MAFFT alignment of query to frame (1: plus strand, 1r: minus strand)")
    print ("chr_mafft.py step2/step2r/step2s   ##Scoring to form large blocks and identify satellite positions (2: plus strand, 2r: minus strand, 2s: summary)")
    print ("chr_mafft.py step3                 ##SEAT - Place each qualified base into seat positions")
    print ("chr_mafft.py step4                 ##Monomer - Identify circular monomers from seat arrays")
    print ("chr_mafft.py step5                 ##HOR - Primary higher-order repeat identification")
    #################
    print ("")
    
    
    
    print ("-monomer\t\tSpecify monomer sequence, Default monomer is for grapevine VSat1 (28bp)")
    print ("-thread \t\tNumber of threads (default: 50), used for multiprocessing")
    print ("-i      \t\tInput FASTA file (required for step0)")
    print ("-c      \t\tRun from the specified step to the end of the script")   
    print ("-simple \t\tOnly output all monomers(subunits), do not perform statistics")
    print (" ")
    print('Mafft is required')
    print('R, dplyr, and ggplot2 are required for plotting, but not required when using -simple')
    sys.exit()
argv1=argvs[1]

import subprocess
import csv
import os
import math
import time
import timeit 
import copy
import re 
from multiprocessing import Pool, cpu_count


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

time_start=timeit.default_timer()
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  

if "thread" in args_dict:  thread=int(args_dict["thread"])  
else:thread=70

#软件需求
mafft="mafft"   #software #Can change the path



if  os.path.exists('./chr_mafft')==False:
    subprocess.run(["mkdir chr_mafft"], shell=True)


# Generate reverse complementary sequence
def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_complement_seq = ''.join(complement[base] for base in reversed(sequence))
    return reverse_complement_seq
    
monomer="ACTCGCACGGATTCTACCATTTTTCGGT"  #VSat1（28bp）

#VSat2(CEN66):CCTGAGTTGTTCATAGCATCTAGAGCCTTGAGCTCATGACCTAGAGTCGTTCATGTCATCCACAAC
#VSat3(CEN103):TCATGACCGATCGGATCGGGTGCGGTCTATGATGAAAACCAGATAGGACAACGACGTGACCGATCGTATAGGTGTGGTCTACGCCGAAAATGAAACTGAACTG
#VroSat1(Vitis Rotundifolia;21bp)：ACTTTGTTAGAGAGGGTACGA



monomer_r=reverse_complement(monomer)

window_size=math.ceil(len(monomer)*1.1)



if "monomer" in args_dict:  monomer=args_dict["monomer"].upper()   
monomer_len=len(monomer)
frame_repeat_num= math.ceil(1000/monomer_len)+math.floor(500/monomer_len)  

stat_mark=''
if "simple" in args_dict:  stat_mark='yes'

print('Script start')

frame = monomer*frame_repeat_num
frame_r = monomer_r*frame_repeat_num

#step0
if argv1=="stepall" or argv1=="step0":
    def process_fasta(input_file, output_dir):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        info_table = []  #List to store information table
        sequence_data = {}  #Dictionary to store sequence data
        current_id = None
        current_sequence = [] 
        with open(input_file, 'r') as infile:
            for line in infile:
                line = line.strip()
                if line.startswith('>'):
                    #Process new sequence header
                    if current_id and current_sequence:
                        #Add previous sequence data to dictionary and info table
                        cleaned_id = current_id[1:].split(" ")[0]
                        sequence_str = ''.join(current_sequence)  #Remove newlines and concatenate sequence
                        sequence_length = len(sequence_str)
                        #Assume start position is 1, end position is sequence length
                        start_pos = 1
                        end_pos = sequence_length
                        if end_pos>10000000:
                            info_table.append([cleaned_id, start_pos, end_pos])
                            sequence_data[cleaned_id] = sequence_str.upper()
                    #Update current sequence ID and reset sequence list
                    current_id = line
                    current_sequence = []
                else: 
                    #Process sequence line
                    current_sequence.append(line)
        #Process the last sequence if exists
        if current_id and current_sequence:
            cleaned_id = current_id[1:].split(" ")[0]
            sequence_str = ''.join(current_sequence)  #Remove newlines and concatenate sequence
            sequence_str=sequence_str.upper()
            sequence_length = len(sequence_str)
            #Assume start position is 1, end position is sequence length
            start_pos = 1
            end_pos = sequence_length
            if end_pos>10000000:
                info_table.append([cleaned_id, start_pos, end_pos])
                sequence_data[cleaned_id] = sequence_str
        #Write information table to CSV file
        info_table_file = os.path.join(output_dir, 'region_info.all')
        with open(info_table_file, 'w') as outfile:
            outfile.write('chr\tstart\tend\n')
            for row in info_table:
                outfile.write('\t'.join(map(str, row)) + '\n')
        #Write each sequence to a separate file
        for seq_id, seq_str in sequence_data.items():
            output_file = os.path.join(output_dir, f'{seq_id}.fasta')
            with open(output_file, 'w') as outfile:
                outfile.write(f'>{seq_id}\n{seq_str}\n')
    if os.path.exists("./chr_mafft/0_prepare")==False:
        subprocess.run(["mkdir ./chr_mafft/0_prepare"], shell=True)
    print("A coordinate file containing all chromosomes has been generated. Edit this file to narrow down the analysis range")
    if "i"  in args_dict: #print("Missing input fasta file");sys.exit()
        input_fasta = args_dict["i"]
    else:
        file_list = os.listdir('./')
        dir_fasta_file_num=0
        for one in file_list:
            if     one.endswith('.fasta'):input_fasta=one;      dir_fasta_file_num+=1
            elif   one.endswith('.fa'):input_fasta=one;         dir_fasta_file_num+=1
            elif   one.endswith('.fna'):input_fasta=one;        dir_fasta_file_num+=1
        if  dir_fasta_file_num==0:
            print("Missing input fasta file");sys.exit()
        elif  dir_fasta_file_num>1:
            print("Input fasta file not specified");sys.exit()
    process_fasta(input_fasta, "./chr_mafft/0_prepare")
    
    #The following section is for sorting, so that region1 corresponds to chr1
    arr_chr_info=[]
    k=0
    with open("./chr_mafft/0_prepare/region_info.all",'r') as f :
        for line in f:
            k+=1
            if k==1:continue
            eachline_arr=line.strip().split('\t')
            chr_name=eachline_arr[0]
            print(chr_name)
            version_num=''
            for x in chr_name[::-1]:
                #print(x)
                if x.isdigit()==True:
                    version_num=x+version_num
                    #print('version_num:'+version_num)
                else:break   
            if version_num=='':print('Error: Chromosome name does not end with a number.');sys.exit()
            arr_chr_info.append([eachline_arr[0],int(version_num),eachline_arr[1],eachline_arr[2]])
    #print(arr_chr_info)
    arr_chr_info_sort=     sorted(arr_chr_info, key=lambda x: x[1], reverse=False)   
    with open('./chr_mafft/0_prepare/region_info.set', 'w') as f:
        f.write('chr\tstart\tend\n')
        for row in arr_chr_info_sort:
            f.write(f"{row[0]}\t{row[2]}\t{row[3]}\n")    
    subprocess.run(["cp ./chr_mafft/0_prepare/region_info.all ./chr_mafft/0_prepare/region_info.set"], shell=True)

#
start_step=100000
continue_state=''
for one_argv in args_dict:
    if "continue"==one_argv or 'c'==one_argv: continue_state='yes'
if continue_state=='yes':
    print('Run from the specified step to the end of the script.\n')
    if    argv1=='step1':    start_step=1
    elif  argv1=='step1.1':  start_step=2
    elif  argv1=='step1.2':  start_step=3
    elif  argv1=='step1.3':  start_step=4
    elif  argv1=='step1.4':  start_step=5
    elif  argv1=='step1.5':  start_step=6
    elif  argv1=='step1.6':  start_step=7
    elif  argv1=='step1.7':  start_step=8
    elif  argv1=='step1.8':  start_step=9
    elif  argv1=='step1r':    start_step=10
    elif  argv1=='step1.1r':  start_step=11
    elif  argv1=='step1.2r':  start_step=12
    elif  argv1=='step1.3r':  start_step=13
    elif  argv1=='step1.4r':  start_step=14
    elif  argv1=='step1.5r':  start_step=15
    elif  argv1=='step1.6r':  start_step=16
    elif  argv1=='step1.7r':  start_step=17
    elif  argv1=='step1.8r':  start_step=18
    elif  argv1=='step1c':   start_step=19
    elif  argv1=='step2':    start_step=20
    elif  argv1=='step2.0':  start_step=21
    elif  argv1=='step2.1':  start_step=22
    elif  argv1=='step2.2':  start_step=23
    elif  argv1=='step2.3':  start_step=24   
    elif  argv1=='step2r':    start_step=25
    elif  argv1=='step2.0r':  start_step=26
    elif  argv1=='step2.1r':  start_step=27
    elif  argv1=='step2.2r':  start_step=28
    elif  argv1=='step2.3r':  start_step=29  
    elif  argv1=='step2s':    start_step=30
    elif  argv1=='step3':    start_step=31
    elif  argv1=='step3.0':    start_step=32
    elif  argv1=='step3.1':    start_step=33
    elif  argv1=='step3.2':    start_step=34
    elif  argv1=='step3.3':    start_step=35
    elif  argv1=='step3.3s1':    start_step=36    
    elif  argv1=='step3.3s2':    start_step=37 
    elif  argv1=='step3.3p1':    start_step=38 
    elif  argv1=='step3.3p2':    start_step=39 
    elif  argv1=='step3.3p3':    start_step=40 
    elif  argv1=='step3.3p4':    start_step=41 
    elif  argv1=='step3.3p5':    start_step=42  
    elif  argv1=='step4':        start_step=43  
    elif  argv1=='step4.0':      start_step=44 
    elif  argv1=='step4.0_plus': start_step=44 
    elif  argv1=='step4.0p1':    start_step=45  
    elif  argv1=='step4.1':      start_step=46  
    elif  argv1=='step4.1p1':    start_step=47  
    elif  argv1=='step4.2':      start_step=48  
    elif  argv1=='step4.3':      start_step=49 
    elif  argv1=='step5':       start_step=50     
    elif  argv1=='step5.0':    start_step=51     
    elif  argv1=='step5.1':    start_step=52     
    elif  argv1=='step5.1p1':    start_step=53     
   

#step1,mafft，除了mafft预计130s
if argv1=="stepall" or argv1=="step1" or  argv1=="step1.1" or start_step<=1:
    print("Default monomer: "+monomer)
    print("Monomer duplicated to generate alignment template, copy number: "+str(frame_repeat_num))   #Approximately 1500bp
    print('step1.1 ———— Splitting sequences, preparing for MAFFT, then concatenating results')
    if os.path.exists("./chr_mafft/1_mafft")==False:
        subprocess.run(["mkdir ./chr_mafft/1_mafft"], shell=True)
    if not os.path.exists("./chr_mafft/0_prepare/region_info.set"):print("Missing ./chr_mafft/0_prepare/region_info.set, please run step0 first");sys.exit()
    region_list=[];dict_region_info={};i=0
    with open("./chr_mafft/0_prepare/region_info.set",'r') as f:
        for line in f.readlines():
            eachline = line.strip()     
            if len(eachline)==0:continue  
            if eachline[0]=="#":continue
            i+=1
            if i==1:continue
            eachline_arr=eachline.split('\t')
            region_list.append(eachline_arr[0])
            if len(eachline_arr)!=3:continue
            dict_region_info[i-1]=[eachline_arr[0],eachline_arr[1],eachline_arr[2]]
            subprocess.run(["mkdir ./chr_mafft/1_mafft/region_"+str(i-1)], shell=True)  
    region_num=len(region_list)
    print("Number of regions:"+str(region_num))
    ############
    i=0;region_fa_file_list=[];
    with open('./chr_mafft/1_mafft/mafft_result_list','w') as f2:
        f2.write('')
    sum_part_serial=0
    while i<region_num:
        i+=1
        one_region = dict_region_info[i]
        print('Split sequence:',one_region)
        one_region_chr = one_region[0]
        one_region_start = int(one_region[1])
        one_region_end  = int(one_region[2])
    
        one_region_str=one_region_chr+':'+str(one_region_start)+'-'+str(one_region_end)
        one_region_chr_file='./chr_mafft/0_prepare/'+one_region_chr+'.fasta'
        with open(one_region_chr_file,'r') as f:
            for line in f.readlines():
                eachline=line.strip()
                if eachline[0]=='>':continue
                else: one_region_chr_seq=eachline
        ###
        one_start = one_region_start-500
        part_num=0   ; one_mafft_result_list=[]
        end_mark=''
        while one_start < one_region_end-500 or part_num==0:  #Ensure each segment is ~1000bp, last segment 1000-1500bp
            if end_mark=='yes' :break
            one_start+=500;part_num+=1;sum_part_serial+=1
            one_end = one_start+999   
            frame_seq=frame
            if one_region_end-one_end<500: one_end=one_region_end  ;end_mark='yes' ; frame_seq=frame_seq*2;frame_seq=frame_seq[:1900];
            ##
            one_seq = one_region_chr_seq[one_start-1:one_end].upper() 
            one_input='./chr_mafft/1_mafft/region_'+str(i)+'/region_'+str(part_num)
            one_output='./chr_mafft/1_mafft/region_'+str(i)+'/mafft_'+str(part_num)
            region_fa_file_list.append((str(sum_part_serial),one_input,one_output))    
            one_mafft_result_list.append(one_output)
            with open(one_input,'w') as f:
                f.write('>frame\n'+frame_seq+'\n>part\n'+one_seq+'\n')
        one_mafft_result_list_str=','.join(one_mafft_result_list)
        with open('./chr_mafft/1_mafft/mafft_result_list','a') as f2:
            f2.write('region_'+str(i)+'\t'+one_region_str+'\t'+one_mafft_result_list_str+'\n')       
    ############## 
    region_fa_file_list_num=len(region_fa_file_list)
    print('Number of sequences to process with MAFFT:'+str(region_fa_file_list_num))
    print('step1.1 ———— MAFFT task started')
    def run_mafft(region_fa_file_list):
        sum_part_serial, input_file, output_file = region_fa_file_list
        if os.path.exists(output_file)==False: 
            os.system(f"{mafft} --quiet --op 4 --genafpair '{input_file}' > '{output_file}'")
    with Pool(processes=thread) as pool:
        pool.map(run_mafft, region_fa_file_list)     
    ########################Concatenate MAFFT results
if argv1=="stepall" or argv1=="step1" or  argv1=="step1.2" or start_step<=3:
    print('step1.2 ———— Pairing MAFFT results') #Estimated 9s
    with open('./chr_mafft/1_mafft/mafft_result_list2_pair','w') as f2:
        f2.write('')
    kk=1
    with open('./chr_mafft/1_mafft/mafft_result_list','r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=3:continue
            region_name=eachline_arr[0]
            region_pos=eachline_arr[1]
            mafft_result_str=eachline_arr[2]
            region_mafftresult_list=mafft_result_str.split(',')
            region_mafftresult_list_num=len(region_mafftresult_list)
            j=0
            while j< region_mafftresult_list_num-1:
                with open('./chr_mafft/1_mafft/mafft_result_list2_pair','a') as f2:
                    f2.write(str(kk)+'\t'+region_name+'\t'+region_pos+'\t'+region_mafftresult_list[j]+'\t'+region_mafftresult_list[j+1]+'\n')
                j+=1;kk+=1
if argv1=="stepall" or argv1=="step1" or  argv1=="step1.3" or start_step<=4:
    print('step1.3 ———— Calculating good_joint from MAFFT results') #Estimated 42s
    def run_mafft_calc(region_mafftresult_info):
        pair_name,region_name,region_pos,mafft_file1,mafft_file2=region_mafftresult_info
        #Read MAFFT1 alignment results
        frame1_seq='';query1_seq=''
        with open(mafft_file1,'r') as f:
            for line in f.readlines():
                eachline = line.strip()
                if eachline=='>frame':mark=1
                elif eachline=='>part':mark=2
                elif mark==1:frame1_seq+=eachline
                else: query1_seq+=eachline
        if(len(frame1_seq)!=len(query1_seq)): print('error,len||frame1_seq!=query1_seq')
        mafft1_len=len(frame1_seq)
        #Read MAFFT2 alignment results
        frame2_seq='';query2_seq=''
        with open(mafft_file2,'r') as f:
            for line in f.readlines():
                eachline = line.strip()
                if eachline=='>frame':mark=1
                elif eachline=='>part':mark=2
                elif mark==1:frame2_seq+=eachline
                else: query2_seq+=eachline
        if(len(frame2_seq)!=len(query2_seq)): print('error,len||frame2_seq!=query2_seq')
        mafft2_len=len(frame2_seq)

        #Store MAFFT1 results in dictionary
        jj=0;query_base_serial=0;frame_base_serial=0;dict_mafft1={}
        while  jj<mafft1_len:
            current_frame_base=frame1_seq[jj];       
            if current_frame_base!='-':  frame_base_serial+=1
            current_query_base=query1_seq[jj];       
            if current_query_base!='-':  
                query_base_serial+=1
                frame_circle_serial=frame_base_serial%monomer_len
                if current_frame_base!='-':    frame_circle_serial_str=str(frame_circle_serial)
                else:                           frame_circle_serial_str=str(frame_circle_serial)+'-'
                dict_mafft1[str(query_base_serial)]=[jj,frame_circle_serial_str]
            jj+=1
        #Store MAFFT2 results in dictionary
        jj=0;query_base_serial=0;frame_base_serial=0;dict_mafft2={}
        while  jj<mafft2_len:
            current_frame_base=frame2_seq[jj];       
            if current_frame_base!='-':  frame_base_serial+=1
            current_query_base=query2_seq[jj];       
            if current_query_base!='-':  
                query_base_serial+=1
                frame_circle_serial=frame_base_serial%monomer_len
                if current_frame_base!='-':    frame_circle_serial_str=str(frame_circle_serial)
                else:                           frame_circle_serial_str=str(frame_circle_serial)+'-'
                dict_mafft2[str(query_base_serial)]=[jj,frame_circle_serial_str]
            jj+=1            

        ##Identify good_joint positions between adjacent fragments
        i=0;good_joint1='NA';good_joint2='NA';
        while i<245: 
            serial1=str(750-i); serial2=str(250-i)
            if dict_mafft1[serial1][1]==dict_mafft2[serial2][1]: 
                good_joint1=dict_mafft1[serial1][0];     good_joint2=dict_mafft2[serial2][0];        break
            serial1=str(750+i); serial2=str(250+i)
            if dict_mafft1[serial1][1]==dict_mafft2[serial2][1]: 
                good_joint1=dict_mafft1[serial1][0];     good_joint2=dict_mafft2[serial2][0];        break
            i+=1    
    
        ##Save results
        with open('./chr_mafft/1_mafft/mafft_result_list3_pair_tmp','a') as f2:
            f2.write(pair_name+'\t'+region_name+'\t'+region_pos+'\t'+mafft_file1+'\t'+mafft_file2+'\t'+serial1+'\t'+serial2+'\t'+str(good_joint1)+'\t'+str(good_joint2)+'\n')
           

    #Create a process pool
    region_mafftresult_info=[];
    with open('./chr_mafft/1_mafft/mafft_result_list2_pair','r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=5:continue
            pair_name=eachline_arr[0]
            region_name=eachline_arr[1]
            region_pos=eachline_arr[2]
            mafft_file1=eachline_arr[3]
            mafft_file2=eachline_arr[4]
            eachline_tuple=(pair_name,region_name,region_pos,mafft_file1,mafft_file2)
            region_mafftresult_info.append(eachline_tuple)
            
    with open('./chr_mafft/1_mafft/mafft_result_list3_pair_tmp','w') as f2:
        f2.write('')

    with Pool(processes=thread) as pool:
        #Assign tasks to processes in the pool
        pool.map(run_mafft_calc, region_mafftresult_info)     
    subprocess.run(["sort -n -k 1 ./chr_mafft/1_mafft/mafft_result_list3_pair_tmp > ./chr_mafft/1_mafft/mafft_result_list3_pair; rm ./chr_mafft/1_mafft/mafft_result_list3_pair_tmp"], shell=True)  
if argv1=="stepall" or argv1=="step1" or  argv1=="step1.4" or start_step<=5:
    print('step1.4 ———— Annotating MAFFT results from step1.3')#Estimated 40s
    print('Reading')
    dict_mafftfile_joint={};dict_mafftfile_source={}
    mafft_file_set=set()
    with open('./chr_mafft/1_mafft/mafft_result_list4_files_tmp','w') as f2:
        f2.write('')
    with open('./chr_mafft/1_mafft/mafft_result_list3_pair','r') as f:
        mafft_result_list3_pair_line_num=len(f.readlines())
    i=0 ;
    with open('./chr_mafft/1_mafft/mafft_result_list3_pair','r') as f:
        for line in f.readlines():
            i+=1
            eachline_arr = line.strip().split('\t')
            part_serial=eachline_arr[0]
            region_name=eachline_arr[1]
            region_pos=eachline_arr[2]
            mafft_file1=eachline_arr[3]
            mafft_file2=eachline_arr[4]
            mafft_joint1=eachline_arr[7]
            mafft_joint2=eachline_arr[8] 
            dict_mafftfile_source[mafft_file1]=  [region_name,region_pos]         
            dict_mafftfile_source[mafft_file2]=  [region_name,region_pos]   
            #Check and update dict_mafftfile_joint
            if mafft_file1 not in dict_mafftfile_joint:
                dict_mafftfile_joint[mafft_file1] = {}
            if mafft_file2 not in dict_mafftfile_joint:
                dict_mafftfile_joint[mafft_file2] = {}
            dict_mafftfile_joint[mafft_file1]['B']=mafft_joint1
            dict_mafftfile_joint[mafft_file2]['A']=mafft_joint2
            mafft_file_set.update([mafft_file1,mafft_file2])
    print('\nOutput')
    mafft_file_list=list(mafft_file_set);mafft_file_list.sort()
    mafft_file_list_num=len(mafft_file_list)
    i=0
    while i<mafft_file_list_num:                                                    
        i+=1
        mafft_file_current               = mafft_file_list[i-1]   
        region_name                    = dict_mafftfile_source[mafft_file_current][0];     
        region_pos                      = dict_mafftfile_source[mafft_file_current][1];    
        part_serial     =  mafft_file_current.split('/')[-1][6:]
        if 'A' in dict_mafftfile_joint[mafft_file_current]:      mafft_joint_A = dict_mafftfile_joint[mafft_file_current]['A']
        else:                                                    mafft_joint_A = '0'
        if 'B' in dict_mafftfile_joint[mafft_file_current]:      mafft_joint_B = dict_mafftfile_joint[mafft_file_current]['B']
        else:                                                    mafft_joint_B = '-1'         
        #Output
        with open('./chr_mafft/1_mafft/mafft_result_list4_files_tmp','a') as f2:
            f2.write(mafft_file_current+'\t'+region_name+'\t'+region_pos+'\t'+part_serial+'\t'+mafft_joint_A+'\t'+mafft_joint_B+'\n')
        
    print()
    subprocess.run(["sort  -k2,2V -k4,4n ./chr_mafft/1_mafft/mafft_result_list4_files_tmp > ./chr_mafft/1_mafft/mafft_result_list4_files; "], shell=True) 
    subprocess.run(["rm ./chr_mafft/1_mafft/mafft_result_list4_files_tmp "], shell=True)  
if argv1=="stepall" or argv1=="step1" or  argv1=="step1.5" or start_step<=6:
    print('step1.5 ———— Determining output boundaries from MAFFT results')#Estimated 15s
    with open('./chr_mafft/1_mafft/mafft_result_list5_files','w') as f2:
        f2.write('')
    with open('./chr_mafft/1_mafft/mafft_result_list4_files','r') as f:
        mafft_result_list4_files_line_num=len(f.readlines())
    i=0 ;NA_state=''
    with open('./chr_mafft/1_mafft/mafft_result_list4_files','r') as f:
        for line in f.readlines():
            i+=1
            eachline_arr = line.strip().split('\t')
            mafft_file=eachline_arr[0]
            region_name=eachline_arr[1]
            region_pos=eachline_arr[2]
            mafft_joint1=eachline_arr[4]
            mafft_joint2=eachline_arr[5] 
            ##Normal case
            if NA_state==''       and mafft_joint1!='NA' and mafft_joint2!='NA':   new_mafft_joint1=mafft_joint1;    new_mafft_joint2=mafft_joint2
            ##Cannot match the latter, output full segment
            elif NA_state==''     and mafft_joint1!='NA' and mafft_joint2=='NA':   new_mafft_joint1=mafft_joint1;    new_mafft_joint2=''  ;           NA_state='skip_before'           
            ##Do not output this line
            elif NA_state=='skip_before'   :                                          new_mafft_joint1='';                new_mafft_joint2=''  ;           NA_state='skip_after'           
            ##Previous line not output, start from 0 for this line
            elif NA_state=='skip_after'                   and  mafft_joint2!='NA':   new_mafft_joint1='0';              new_mafft_joint2=mafft_joint2 ;    NA_state=''    
            ##Previous line not output, start from 0, and mafft_joint2 remains NA
            elif NA_state=='skip_after'                   and  mafft_joint2=='NA':   new_mafft_joint1='0';              new_mafft_joint2='' ;            NA_state='skip_before'           
            
            else: 
                print('')
                print('NA_state:'+NA_state)
                print('region_name:'+region_name)
                print('region_pos:'+region_pos)
                print('mafft_joint1:'+mafft_joint1)
                print('mafft_joint2:'+mafft_joint2)
                print('')
  
            with open('./chr_mafft/1_mafft/mafft_result_list5_files','a') as f2:
                f2.write(mafft_file+'\t'+region_name+'\t'+region_pos+'\t'+mafft_joint1+'\t'+mafft_joint2+'\t'+new_mafft_joint1+'\t'+new_mafft_joint2+'\t'+'|'+'\n')
if argv1=="stepall" or argv1=="step1" or  argv1=="step1.6" or start_step<=7:
    print('step1.6 ———— Extracting sequences from MAFFT results')#Estimated 5s
    with open('./chr_mafft/1_mafft/mafft_result_list6_files','w') as f2:
        f2.write('')
    i=0 ;
    region_mafftresult_info=[]
    with open('./chr_mafft/1_mafft/mafft_result_list5_files','r') as f:
        for line in f.readlines():
            i+=1
            eachline_arr = line.strip().split('\t')
            if len(eachline_arr)!=8:continue
            mafft_file=eachline_arr[0]
            region_name=eachline_arr[1]
            region_pos=eachline_arr[2]
            mafft_joint1=eachline_arr[5]
            if mafft_joint1=='':continue
            else:                               mafft_joint1=int(mafft_joint1)
            mafft_joint2=eachline_arr[6] ;
            if mafft_joint2!='':                mafft_joint2=int(mafft_joint2)
            eachline_tuple=(mafft_file,region_name,region_pos,mafft_joint1,mafft_joint2)
            region_mafftresult_info.append(eachline_tuple)
    print('Generating')
    def run_mafft_getseq(region_mafftresult_info):
        mafft_file,region_name,region_pos,mafft_joint1,mafft_joint2=region_mafftresult_info
        
        #Read MAFFT alignment results
        frame_seq='';query_seq=''
        with open(mafft_file,'r') as f:
            for line in f.readlines():
                eachline = line.strip()
                if eachline=='>frame':mark=1
                elif eachline=='>part':mark=2
                elif mark==1:frame_seq+=eachline
                else: query_seq+=eachline           
         ##
        if mafft_joint2!='':
            frame_seq_splice = frame_seq[mafft_joint1:mafft_joint2]
            query_seq_splice = query_seq[mafft_joint1:mafft_joint2]
        else:
            frame_seq_splice = frame_seq[mafft_joint1:]
            query_seq_splice = query_seq[mafft_joint1:]
        with open('./chr_mafft/1_mafft/mafft_result_list6_files','a') as f2:
            f2.write(mafft_file+'\t'+region_name+'\t'+region_pos+'\t'+str(mafft_joint1)+'\t'+str(mafft_joint2)+'\t'+frame_seq_splice+'\t'+query_seq_splice+'\n')
    with Pool(processes=thread) as pool:
        #Assign tasks to processes in the pool
        pool.map(run_mafft_getseq, region_mafftresult_info) 
if argv1=="stepall" or argv1=="step1" or  argv1=="step1.7" or start_step<=8:
    print('step1.7 ———— Concatenating MAFFT sequences')#Estimated 7s
    if  os.path.exists('./chr_mafft/1_mafft/sum'):
        subprocess.run(["rm -r ./chr_mafft/1_mafft/sum"], shell=True)
    subprocess.run(["mkdir ./chr_mafft/1_mafft/sum"], shell=True)
    subprocess.run(["mkdir ./chr_mafft/1_mafft/sum/frame"], shell=True)    
    subprocess.run(["mkdir ./chr_mafft/1_mafft/sum/part"], shell=True)    
    dict_mafftfile_seq={}
    with open('./chr_mafft/1_mafft/mafft_result_list6_files','r') as f:
        for line in f.readlines():
            eachline_arr= line.strip().split('\t')
            mafftfile=eachline_arr[0]
            frame_seq=eachline_arr[5]
            part_seq=eachline_arr[6]
            dict_mafftfile_seq[mafftfile]=[frame_seq,part_seq]
    with open('./chr_mafft/1_mafft/mafft_result_list','r') as f:
        for line in f.readlines():
            eachline_arr= line.strip().split('\t')
            region_name=eachline_arr[0]
            region_pos=eachline_arr[1]
            mafft_file_arr=eachline_arr[2].split(',')
            print(region_name)
            region_frame_seq_sum=''
            region_part_seq_sum=''
            k=0
            with open('./chr_mafft/1_mafft/sum/frame/'+region_name,'a') as f2 , open('./chr_mafft/1_mafft/sum/part/'+region_name,'a') as f3:            
                for one_mafft_file in mafft_file_arr:
                    k+=1
                    if one_mafft_file in dict_mafftfile_seq:
                        frame_part_seq=dict_mafftfile_seq[one_mafft_file]
                        frame_seq=frame_part_seq[0]
                        part_seq=frame_part_seq[1]
                        f2.write(frame_seq+'\n')
                        f3.write(part_seq+'\n')
if argv1=="stepall" or argv1=="step1" or  argv1=="step1.8" or start_step<=9:
    print('step1.8 ———— Validating results')#Estimated 9s
    region_name_pos_list=[]
    with open('./chr_mafft/1_mafft/mafft_result_list','r') as f:
        for line in f.readlines():
            eachline_arr= line.strip().split('\t')
            region_name=eachline_arr[0]
            region_pos=eachline_arr[1]
            eachline_tuple=(region_name,region_pos)
            region_name_pos_list.append(eachline_tuple)
    ####

    def run_result_check(region_name_pos_list):
        one_region_name,region_pos=region_name_pos_list
        frame_seq=''
        with open('./chr_mafft/1_mafft/sum/frame/'+one_region_name) as f:
            for line in f.readlines():
                frame_seq+=line.strip()
        query_seq=''
        with open('./chr_mafft/1_mafft/sum/part/'+one_region_name) as f:
            for line in f.readlines():
                query_seq+=line.strip()
        print()
        frame_seq=frame_seq.upper()
        query_seq=query_seq.upper()

        print("one_region_name:"+'\t'+one_region_name+'\t'+str(len(frame_seq))+'\t'+str(len(query_seq))+':\n')
        uu1=0;jj=-1
        for x in frame_seq:
            x=x.strip();
            if x=='':continue
            jj+=1
            if x=='-':continue
            if x!=(monomer[uu1%monomer_len]):
                print('Error in frame sequence');
                print(frame_seq[jj-10:jj]+'_'+frame_seq[jj]+'_'+frame_seq[jj+1:jj+10]);print(x);print(monomer[uu1%monomer_len]);
                print('Error index:'+str(jj)+'\t'+region_name+'\t'+region_pos+'\n');sys.exit();break
            uu1+=1
        uu2=0;part_seq_sum=''
        for x in query_seq:
            if x=='-':continue
            uu2+=1;part_seq_sum+=x
        region_chr=region_pos.split(':')[0]
        region_start_end=region_pos.split(':')[1].split('-')
        region_length=int(region_start_end[1])-int(region_start_end[0])+1
        subprocess.run(["samtools faidx ./chr_mafft/0_prepare/"+region_chr+'.fasta '+region_pos+'>./chr_mafft/1_mafft/'+region_pos], shell=True)
        region_seq=''
        with open('./chr_mafft/1_mafft/'+region_pos,'r') as f:
            for line in f.readlines():
                eachline = line.strip()
                if eachline[0]=='>':continue
                region_seq+=eachline
        subprocess.run(['rm ./chr_mafft/1_mafft/'+region_pos], shell=True)
        if region_seq!=part_seq_sum:
           
            part_seq_sum_len=len(part_seq_sum)
            tt=0;error_mark=''
            while  tt<part_seq_sum_len:
                if region_seq[tt]!=part_seq_sum[tt]:print('Error index:'+str(tt)+'\t'+region_name+'\t'+region_pos);print("Correct sequence:"+region_seq[tt]);print("Wrong sequence:"+part_seq_sum[tt]);error_mark='yes';break
                tt+=1  
            if error_mark=='yes':  print('Error: Sequence mismatch!'+region_name);sys.exit()
            else:  print('Warning: Sequence mismatch, but sequences are consistent. The last part was not included, does not affect analysis.'+region_name)
        else:print('Query sequence: samtools extraction matches concatenated result')
        print(one_region_name,str(len(frame_seq)),str(len(query_seq)))
        
        with open('./chr_mafft/1_mafft/mafft_result_sum','a') as f2:
            f2.write(one_region_name+'\t'+region_pos+'\t'+frame_seq+'\t'+query_seq+'\n')
    with open('./chr_mafft/1_mafft/mafft_result_sum','w') as f2:
        f2.write('')
    with Pool(processes=thread) as pool:
        #Assign tasks to processes in the pool
        pool.map(run_result_check, region_name_pos_list) 
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))) 
#step1r
if argv1=="stepall" or argv1=="step1r" or  argv1=="step1.1r" or start_step<=10:
    print("Default monomer (reverse): "+monomer_r)
    print("Monomer duplicated to generate alignment template, copy number: "+str(frame_repeat_num))   #Approximately 1500bp
    print('step1r ———— Splitting sequences, preparing for MAFFT (reverse strand), then concatenating results')
    if os.path.exists("./chr_mafft/1_mafft_r")==False:
        subprocess.run(["mkdir ./chr_mafft/1_mafft_r"], shell=True)
    if not os.path.exists("./chr_mafft/0_prepare/region_info.set"):print("Missing ./chr_mafft/0_prepare/region_info.set, please run step0 first");sys.exit()
    region_list=[];dict_region_info={};i=0
    with open("./chr_mafft/0_prepare/region_info.set",'r') as f:
        for line in f.readlines():
            eachline = line.strip()     
            if len(eachline)==0:continue  
            if eachline[0]=="#":continue
            i+=1
            if i==1:continue
            eachline_arr=eachline.split('\t')
            region_list.append(eachline_arr[0])
            if len(eachline_arr)!=3:continue
            dict_region_info[i-1]=[eachline_arr[0],eachline_arr[1],eachline_arr[2]]
            subprocess.run(["mkdir ./chr_mafft/1_mafft_r/region_"+str(i-1)], shell=True)  
    region_num=len(region_list)
    print("Number of regions:"+str(region_num))
    ############
    i=0;region_fa_file_list=[];
    with open('./chr_mafft/1_mafft_r/mafft_result_list','w') as f2:
        f2.write('')
    sum_part_serial=0
    while i<region_num:
        i+=1
        one_region = dict_region_info[i]
        print('Split sequence:',one_region)
        one_region_chr = one_region[0]
        one_region_start = int(one_region[1])
        one_region_end  = int(one_region[2])
    
        one_region_str=one_region_chr+':'+str(one_region_start)+'-'+str(one_region_end)
        one_region_chr_file='./chr_mafft/0_prepare/'+one_region_chr+'.fasta'
        with open(one_region_chr_file,'r') as f:
            for line in f.readlines():
                eachline=line.strip()
                if eachline[0]=='>':continue
                else: one_region_chr_seq=eachline
        ###
        one_start = one_region_start-500
        part_num=0   ; one_mafft_result_list=[]
        end_mark=''
        while one_start < one_region_end-500 or part_num==0:  #Ensure each segment is ~1000bp, last segment 1000-1500bp
            if end_mark=='yes' :break
            one_start+=500;part_num+=1;sum_part_serial+=1
            one_end = one_start+999   
            frame_seq=frame_r     
            if one_region_end-one_end<500: one_end=one_region_end  ;end_mark='yes' ; frame_seq=frame_seq*2;frame_seq=frame_seq[:1900];
            ##
            one_seq = one_region_chr_seq[one_start-1:one_end].upper() 
            one_input='./chr_mafft/1_mafft_r/region_'+str(i)+'/region_'+str(part_num)
            one_output='./chr_mafft/1_mafft_r/region_'+str(i)+'/mafft_'+str(part_num)
            region_fa_file_list.append((str(sum_part_serial),one_input,one_output))    
            one_mafft_result_list.append(one_output)
            with open(one_input,'w') as f:
                f.write('>frame\n'+frame_seq+'\n>part\n'+one_seq+'\n')
        one_mafft_result_list_str=','.join(one_mafft_result_list)
        with open('./chr_mafft/1_mafft_r/mafft_result_list','a') as f2:
            f2.write('region_'+str(i)+'\t'+one_region_str+'\t'+one_mafft_result_list_str+'\n')       
    ############## 
    region_fa_file_list_num=len(region_fa_file_list)
    print('Number of sequences to process with MAFFT:'+str(region_fa_file_list_num))
    ########################Multi-threaded MAFFT execution
    #Define the task function to run in parallel
    print('step1r ———— MAFFT task started')
    def run_mafft(region_fa_file_list):
        sum_part_serial, input_file, output_file = region_fa_file_list
        if os.path.exists(output_file)==False: 
            os.system(f"mafft --quiet --op 4 --genafpair '{input_file}' > '{output_file}'")
    #Create a process pool
    with Pool(processes=thread) as pool:
        #Assign tasks to processes in the pool
        pool.map(run_mafft, region_fa_file_list)     
    ########################Concatenate MAFFT results
if argv1=="stepall" or argv1=="step1r" or  argv1=="step1.2r" or start_step<=12:
    print('step1.2r ———— Pairing MAFFT results (reverse strand)')
    with open('./chr_mafft/1_mafft_r/mafft_result_list2_pair','w') as f2:
        f2.write('')
    kk=1
    with open('./chr_mafft/1_mafft_r/mafft_result_list','r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=3:continue
            region_name=eachline_arr[0]
            region_pos=eachline_arr[1]
            mafft_result_str=eachline_arr[2]
            region_mafftresult_list=mafft_result_str.split(',')
            region_mafftresult_list_num=len(region_mafftresult_list)
            j=0
            while j< region_mafftresult_list_num-1:
                with open('./chr_mafft/1_mafft_r/mafft_result_list2_pair','a') as f2:
                    f2.write(str(kk)+'\t'+region_name+'\t'+region_pos+'\t'+region_mafftresult_list[j]+'\t'+region_mafftresult_list[j+1]+'\n')
                j+=1;kk+=1
if argv1=="stepall" or argv1=="step1r" or  argv1=="step1.3r" or start_step<=13:
    print('step1.3r ———— Calculating good_joint from MAFFT results (reverse strand)')
    def run_mafft_calc(region_mafftresult_info):
        pair_name,region_name,region_pos,mafft_file1,mafft_file2=region_mafftresult_info
        #Read MAFFT1 alignment results
        frame1_seq='';query1_seq=''
        with open(mafft_file1,'r') as f:
            for line in f.readlines():
                eachline = line.strip()
                if eachline=='>frame':mark=1
                elif eachline=='>part':mark=2
                elif mark==1:frame1_seq+=eachline
                else: query1_seq+=eachline
        if(len(frame1_seq)!=len(query1_seq)): print('error,len||frame1_seq!=query1_seq')
        mafft1_len=len(frame1_seq)
        #Read MAFFT2 alignment results
        frame2_seq='';query2_seq=''
        with open(mafft_file2,'r') as f:
            for line in f.readlines():
                eachline = line.strip()
                if eachline=='>frame':mark=1
                elif eachline=='>part':mark=2
                elif mark==1:frame2_seq+=eachline
                else: query2_seq+=eachline
        if(len(frame2_seq)!=len(query2_seq)): print('error,len||frame2_seq!=query2_seq')
        mafft2_len=len(frame2_seq)

        #Store MAFFT1 results in dictionary
        jj=0;query_base_serial=0;frame_base_serial=0;dict_mafft1={}
        while  jj<mafft1_len:
            current_frame_base=frame1_seq[jj];       
            if current_frame_base!='-':  frame_base_serial+=1
            current_query_base=query1_seq[jj];       
            if current_query_base!='-':  
                query_base_serial+=1
                frame_circle_serial=frame_base_serial%monomer_len
                if current_frame_base!='-':    frame_circle_serial_str=str(frame_circle_serial)
                else:                           frame_circle_serial_str=str(frame_circle_serial)+'-'
                dict_mafft1[str(query_base_serial)]=[jj,frame_circle_serial_str]
            jj+=1
        #Store MAFFT2 results in dictionary
        jj=0;query_base_serial=0;frame_base_serial=0;dict_mafft2={}
        while  jj<mafft2_len:
            current_frame_base=frame2_seq[jj];       
            if current_frame_base!='-':  frame_base_serial+=1
            current_query_base=query2_seq[jj];       
            if current_query_base!='-':  
                query_base_serial+=1
                frame_circle_serial=frame_base_serial%monomer_len
                if current_frame_base!='-':    frame_circle_serial_str=str(frame_circle_serial)
                else:                           frame_circle_serial_str=str(frame_circle_serial)+'-'
                dict_mafft2[str(query_base_serial)]=[jj,frame_circle_serial_str]
            jj+=1            

        ##Identify good_joint positions between adjacent fragments
        i=0;good_joint1='NA';good_joint2='NA';
        while i<245: 
            serial1=str(750-i); serial2=str(250-i)
            if dict_mafft1[serial1][1]==dict_mafft2[serial2][1]: 
                good_joint1=dict_mafft1[serial1][0];     good_joint2=dict_mafft2[serial2][0];        break
            serial1=str(750+i); serial2=str(250+i)
            if dict_mafft1[serial1][1]==dict_mafft2[serial2][1]: 
                good_joint1=dict_mafft1[serial1][0];     good_joint2=dict_mafft2[serial2][0];        break
            i+=1    
    
        ##Save results
        with open('./chr_mafft/1_mafft_r/mafft_result_list3_pair_tmp','a') as f2:
            f2.write(pair_name+'\t'+region_name+'\t'+region_pos+'\t'+mafft_file1+'\t'+mafft_file2+'\t'+serial1+'\t'+serial2+'\t'+str(good_joint1)+'\t'+str(good_joint2)+'\n')
           

    #Create a process pool
    region_mafftresult_info=[];
    with open('./chr_mafft/1_mafft_r/mafft_result_list2_pair','r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=5:continue
            pair_name=eachline_arr[0]
            region_name=eachline_arr[1]
            region_pos=eachline_arr[2]
            mafft_file1=eachline_arr[3]
            mafft_file2=eachline_arr[4]
            eachline_tuple=(pair_name,region_name,region_pos,mafft_file1,mafft_file2)
            region_mafftresult_info.append(eachline_tuple)
            
    with open('./chr_mafft/1_mafft_r/mafft_result_list3_pair_tmp','w') as f2:
        f2.write('')

    with Pool(processes=thread) as pool:
        #Assign tasks to processes in the pool
        pool.map(run_mafft_calc, region_mafftresult_info)     
    subprocess.run(["sort -n -k 1 ./chr_mafft/1_mafft_r/mafft_result_list3_pair_tmp > ./chr_mafft/1_mafft_r/mafft_result_list3_pair; rm ./chr_mafft/1_mafft_r/mafft_result_list3_pair_tmp"], shell=True)  
if argv1=="stepall" or argv1=="step1r" or  argv1=="step1.4r" or start_step<=14:
    print('step1.4r ———— Annotating MAFFT results from step1.3r (reverse strand)')
    print('Reading')
    dict_mafftfile_joint={};dict_mafftfile_source={}
    mafft_file_set=set()
    with open('./chr_mafft/1_mafft_r/mafft_result_list4_files_tmp','w') as f2:
        f2.write('')
    with open('./chr_mafft/1_mafft_r/mafft_result_list3_pair','r') as f:
        mafft_result_list3_pair_line_num=len(f.readlines())
    i=0 ;
    with open('./chr_mafft/1_mafft_r/mafft_result_list3_pair','r') as f:
        for line in f.readlines():
            i+=1
            eachline_arr = line.strip().split('\t')
            part_serial=eachline_arr[0]
            region_name=eachline_arr[1]
            region_pos=eachline_arr[2]
            mafft_file1=eachline_arr[3]
            mafft_file2=eachline_arr[4]
            mafft_joint1=eachline_arr[7]
            mafft_joint2=eachline_arr[8] 
            dict_mafftfile_source[mafft_file1]=  [region_name,region_pos]         
            dict_mafftfile_source[mafft_file2]=  [region_name,region_pos]   
            #Check and update dict_mafftfile_joint
            if mafft_file1 not in dict_mafftfile_joint:
                dict_mafftfile_joint[mafft_file1] = {}
            if mafft_file2 not in dict_mafftfile_joint:
                dict_mafftfile_joint[mafft_file2] = {}
            dict_mafftfile_joint[mafft_file1]['B']=mafft_joint1
            dict_mafftfile_joint[mafft_file2]['A']=mafft_joint2
            mafft_file_set.update([mafft_file1,mafft_file2])
    print('\nOutput')
    mafft_file_list=list(mafft_file_set);mafft_file_list.sort()
    mafft_file_list_num=len(mafft_file_list)
    i=0
    while i<mafft_file_list_num:                                                    
        i+=1
        mafft_file_current               = mafft_file_list[i-1]   
        region_name                    = dict_mafftfile_source[mafft_file_current][0];     
        region_pos                      = dict_mafftfile_source[mafft_file_current][1];    
        part_serial     =  mafft_file_current.split('/')[-1][6:]
        if 'A' in dict_mafftfile_joint[mafft_file_current]:      mafft_joint_A = dict_mafftfile_joint[mafft_file_current]['A']
        else:                                                    mafft_joint_A = '0'
        if 'B' in dict_mafftfile_joint[mafft_file_current]:      mafft_joint_B = dict_mafftfile_joint[mafft_file_current]['B']
        else:                                                    mafft_joint_B = '-1'         
        #Output
        with open('./chr_mafft/1_mafft_r/mafft_result_list4_files_tmp','a') as f2:
            f2.write(mafft_file_current+'\t'+region_name+'\t'+region_pos+'\t'+part_serial+'\t'+mafft_joint_A+'\t'+mafft_joint_B+'\n')
        
    print()
    subprocess.run(["sort  -k2,2V -k4,4n ./chr_mafft/1_mafft_r/mafft_result_list4_files_tmp > ./chr_mafft/1_mafft_r/mafft_result_list4_files; "], shell=True) 
    subprocess.run(["rm ./chr_mafft/1_mafft_r/mafft_result_list4_files_tmp "], shell=True)  
if argv1=="stepall" or argv1=="step1r" or  argv1=="step1.5r" or start_step<=15:
    print('step1.5r ———— Determining output boundaries from MAFFT results (reverse strand)')
    with open('./chr_mafft/1_mafft_r/mafft_result_list5_files','w') as f2:
        f2.write('')
    with open('./chr_mafft/1_mafft_r/mafft_result_list4_files','r') as f:
        mafft_result_list4_files_line_num=len(f.readlines())
    i=0 ;NA_state=''
    with open('./chr_mafft/1_mafft_r/mafft_result_list4_files','r') as f:
        for line in f.readlines():
            i+=1
            eachline_arr = line.strip().split('\t')
            mafft_file=eachline_arr[0]
            region_name=eachline_arr[1]
            region_pos=eachline_arr[2]
            mafft_joint1=eachline_arr[4]
            mafft_joint2=eachline_arr[5] 
            ##Normal case
            if NA_state==''       and mafft_joint1!='NA' and mafft_joint2!='NA':   new_mafft_joint1=mafft_joint1;    new_mafft_joint2=mafft_joint2
            ##Cannot match the latter, output full segment
            elif NA_state==''     and mafft_joint1!='NA' and mafft_joint2=='NA':   new_mafft_joint1=mafft_joint1;    new_mafft_joint2=''  ;           NA_state='skip_before'           
            ##Do not output this line
            elif NA_state=='skip_before'   :                                          new_mafft_joint1='';                new_mafft_joint2=''  ;           NA_state='skip_after'           
            ##Previous line not output, start from 0 for this line
            elif NA_state=='skip_after'                   and  mafft_joint2!='NA':   new_mafft_joint1='0';              new_mafft_joint2=mafft_joint2 ;    NA_state=''    
            ##Previous line not output, start from 0, and mafft_joint2 remains NA
            elif NA_state=='skip_after'                   and  mafft_joint2=='NA':   new_mafft_joint1='0';              new_mafft_joint2='' ;            NA_state='skip_before'           
            
            else: 
                print('')
                print('NA_state:'+NA_state)
                print('region_name:'+region_name)
                print('region_pos:'+region_pos)
                print('mafft_joint1:'+mafft_joint1)
                print('mafft_joint2:'+mafft_joint2)
                print('')
  
            with open('./chr_mafft/1_mafft_r/mafft_result_list5_files','a') as f2:
                f2.write(mafft_file+'\t'+region_name+'\t'+region_pos+'\t'+mafft_joint1+'\t'+mafft_joint2+'\t'+new_mafft_joint1+'\t'+new_mafft_joint2+'\t'+'|'+'\n')
if argv1=="stepall" or argv1=="step1r" or  argv1=="step1.6r" or start_step<=16:
    print('step1.6r ———— Extracting sequences from MAFFT results (reverse strand)')
    with open('./chr_mafft/1_mafft_r/mafft_result_list6_files','w') as f2:
        f2.write('')
    i=0 ;
    region_mafftresult_info=[]
    with open('./chr_mafft/1_mafft_r/mafft_result_list5_files','r') as f:
        for line in f.readlines():
            i+=1
            eachline_arr = line.strip().split('\t')
            if len(eachline_arr)!=8:continue
            mafft_file=eachline_arr[0]
            region_name=eachline_arr[1]
            region_pos=eachline_arr[2]
            mafft_joint1=eachline_arr[5]
            if mafft_joint1=='':continue
            else:                               mafft_joint1=int(mafft_joint1)
            mafft_joint2=eachline_arr[6] ;
            if mafft_joint2!='':                mafft_joint2=int(mafft_joint2)
            eachline_tuple=(mafft_file,region_name,region_pos,mafft_joint1,mafft_joint2)
            region_mafftresult_info.append(eachline_tuple)
    print('Generating')
    def run_mafft_getseq(region_mafftresult_info):
        mafft_file,region_name,region_pos,mafft_joint1,mafft_joint2=region_mafftresult_info
        
        #Read MAFFT alignment results
        frame_seq='';query_seq=''
        with open(mafft_file,'r') as f:
            for line in f.readlines():
                eachline = line.strip()
                if eachline=='>frame':mark=1
                elif eachline=='>part':mark=2
                elif mark==1:frame_seq+=eachline
                else: query_seq+=eachline           
         ##
        if mafft_joint2!='':
            frame_seq_splice = frame_seq[mafft_joint1:mafft_joint2]
            query_seq_splice = query_seq[mafft_joint1:mafft_joint2]
        else:
            frame_seq_splice = frame_seq[mafft_joint1:]
            query_seq_splice = query_seq[mafft_joint1:]
        with open('./chr_mafft/1_mafft_r/mafft_result_list6_files','a') as f2:
            f2.write(mafft_file+'\t'+region_name+'\t'+region_pos+'\t'+str(mafft_joint1)+'\t'+str(mafft_joint2)+'\t'+frame_seq_splice+'\t'+query_seq_splice+'\n')
    with Pool(processes=thread) as pool:
        #Assign tasks to processes in the pool
        pool.map(run_mafft_getseq, region_mafftresult_info) 
if argv1=="stepall" or argv1=="step1r" or  argv1=="step1.7r" or start_step<=17:
    print('step1.7r ———— Concatenating MAFFT sequences (reverse strand)')
    if  os.path.exists('./chr_mafft/1_mafft_r/sum'):
        subprocess.run(["rm -r ./chr_mafft/1_mafft_r/sum"], shell=True)
    subprocess.run(["mkdir ./chr_mafft/1_mafft_r/sum"], shell=True)
    subprocess.run(["mkdir ./chr_mafft/1_mafft_r/sum/frame"], shell=True)    
    subprocess.run(["mkdir ./chr_mafft/1_mafft_r/sum/part"], shell=True)    
    dict_mafftfile_seq={}
    with open('./chr_mafft/1_mafft_r/mafft_result_list6_files','r') as f:
        for line in f.readlines():
            eachline_arr= line.strip().split('\t')
            mafftfile=eachline_arr[0]
            frame_seq=eachline_arr[5]
            part_seq=eachline_arr[6]
            dict_mafftfile_seq[mafftfile]=[frame_seq,part_seq]

    with open('./chr_mafft/1_mafft_r/mafft_result_list','r') as f:
        for line in f.readlines():
            eachline_arr= line.strip().split('\t')
            region_name=eachline_arr[0]
            region_pos=eachline_arr[1]
            mafft_file_arr=eachline_arr[2].split(',')
            print(region_name)
            region_frame_seq_sum=''
            region_part_seq_sum=''
            k=0
            for one_mafft_file in mafft_file_arr:
                k+=1
                if one_mafft_file in dict_mafftfile_seq:
                    frame_part_seq=dict_mafftfile_seq[one_mafft_file]
                    frame_seq=frame_part_seq[0]
                    part_seq=frame_part_seq[1]
                    with open('./chr_mafft/1_mafft_r/sum/frame/'+region_name,'a') as f2:
                        f2.write(frame_seq+'\n')
                    with open('./chr_mafft/1_mafft_r/sum/part/'+region_name,'a') as f3:
                        f3.write(part_seq+'\n')
if argv1=="stepall" or argv1=="step1r" or  argv1=="step1.8r" or start_step<=18:
    print('step1.8r ———— Validating results (reverse strand)')
    region_name_pos_list=[]
    with open('./chr_mafft/1_mafft_r/mafft_result_list','r') as f:
        for line in f.readlines():
            eachline_arr= line.strip().split('\t')
            region_name=eachline_arr[0]
            region_pos=eachline_arr[1]
            eachline_tuple=(region_name,region_pos)
            region_name_pos_list.append(eachline_tuple)
    ####

    def run_result_check(region_name_pos_list):
        one_region_name,region_pos=region_name_pos_list
        frame_seq=''
        with open('./chr_mafft/1_mafft_r/sum/frame/'+one_region_name) as f:
            for line in f.readlines():
                frame_seq+=line.strip()
        query_seq=''
        with open('./chr_mafft/1_mafft_r/sum/part/'+one_region_name) as f:
            for line in f.readlines():
                query_seq+=line.strip()
        print()
        frame_seq=frame_seq.upper()
        query_seq=query_seq.upper()

        print("one_region_name:"+'\t'+one_region_name+'\t'+str(len(frame_seq))+'\t'+str(len(query_seq))+':\n')
        uu1=0;jj=-1
        for x in frame_seq:
            x=x.strip();
            if x=='':continue
            jj+=1
            if x=='-':continue
            if x!=(monomer_r[uu1%monomer_len]):
                print('Error in frame sequence');
                print(frame_seq[jj-10:jj]+'_'+frame_seq[jj]+'_'+frame_seq[jj+1:jj+10]);print(x);print(monomer_r[uu1%monomer_len]);
                print('Error index:'+str(jj)+'\t'+region_name+'\t'+region_pos+'\n');sys.exit();break
            uu1+=1
        uu2=0;part_seq_sum=''
        for x in query_seq:
            if x=='-':continue
            uu2+=1;part_seq_sum+=x
        region_chr=region_pos.split(':')[0]
        region_start_end=region_pos.split(':')[1].split('-')
        region_length=int(region_start_end[1])-int(region_start_end[0])+1
        subprocess.run(["samtools faidx ./chr_mafft/0_prepare/"+region_chr+'.fasta '+region_pos+'>./chr_mafft/1_mafft_r/'+region_pos], shell=True)
        region_seq=''
        with open('./chr_mafft/1_mafft_r/'+region_pos,'r') as f:
            for line in f.readlines():
                eachline = line.strip()
                if eachline[0]=='>':continue
                region_seq+=eachline
        subprocess.run(['rm ./chr_mafft/1_mafft_r/'+region_pos], shell=True)
        if region_seq!=part_seq_sum:
            part_seq_sum_len=len(part_seq_sum)
            tt=0;error_mark=''
            while  tt<part_seq_sum_len:
                if region_seq[tt]!=part_seq_sum[tt]:print('Error index:'+str(tt)+'\t'+region_name+'\t'+region_pos);print("Correct sequence:"+region_seq[tt]);print("Wrong sequence:"+part_seq_sum[tt]);error_mark='yes';break
                tt+=1  
            if error_mark=='yes':  print('Error: Sequence mismatch!'+region_name);sys.exit()
            else:  print('Warning: Sequence mismatch, but sequences are consistent. The last part was not included, does not affect analysis.'+region_name)
        else:print('Query sequence: samtools extraction matches concatenated result')
        print(one_region_name,str(len(frame_seq)),str(len(query_seq)))
        
        with open('./chr_mafft/1_mafft_r/mafft_result_sum','a') as f2:
            f2.write(one_region_name+'\t'+region_pos+'\t'+frame_seq+'\t'+query_seq+'\n')
    with open('./chr_mafft/1_mafft_r/mafft_result_sum','w') as f2:
        f2.write('')
    with Pool(processes=thread) as pool:
        #Assign tasks to processes in the pool
        pool.map(run_result_check, region_name_pos_list) 
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))) 
#
if argv1=="stepall" or argv1=="step1c" or start_step<=19:
    print('Cleaning up step1 and step1r intermediate files')
    subprocess.run(["rm ./chr_mafft/1_mafft/mafft_result_list* ;rm -r ./chr_mafft/1_mafft/region_*; rm -r ./chr_mafft/1_mafft/sum"], shell=True)
    subprocess.run(["rm ./chr_mafft/1_mafft_r/mafft_result_list* ;rm -r ./chr_mafft/1_mafft_r/region_*; rm -r ./chr_mafft/1_mafft_r/sum"], shell=True)

    
#step2, Scoring and region identification, Estimated 150s*2
if argv1=="stepall" or argv1=="step2" or  argv1=="step2.0" or start_step<=20:
    print('step2.0 ———— Scoring MAFFT alignment results')
    if  os.path.exists('./chr_mafft/2_score')==False:
        subprocess.run(["mkdir ./chr_mafft/2_score"], shell=True)
    if  os.path.exists('./chr_mafft/2_score/2.0'):
        subprocess.run(["rm -r ./chr_mafft/2_score/2.0"], shell=True)
    subprocess.run(["mkdir ./chr_mafft/2_score/2.0"], shell=True)
    print('Storing dict: mafft_serial + query_serial + frame_serial')
    def run_score(run_score_one):
        region_name,region_pos,frame_seq,query_seq=run_score_one
        mafft_result_len=len(frame_seq)
        print('Using list comprehension to create new list')
        pair_list = [x + y for x, y in zip(frame_seq, query_seq)]
        pair_list_str="|".join(pair_list);      
        #Apply substitution rules sequentially
        #print('Applying substitution rules 1')
        list_state_mark= re.sub(r'AA|CC|GG|TT', '1', pair_list_str)
        list_state_mark = re.sub(r'[^-|]-', 'd', list_state_mark)
        list_state_mark = re.sub(r'-[^-|]', 'i', list_state_mark)        
        list_state_mark = re.sub(r'[^-|][^-|]', 'm', list_state_mark)
        list_state_mark = list_state_mark.split('|')
        #print('Applying substitution rules 2')
        list_pair_score = re.sub(r'[^-|]-', '10', pair_list_str)
        list_pair_score = re.sub(r'-[^-|]', '01', list_pair_score)        
        list_pair_score = re.sub(r'[^-|0][^-|0]', '11', list_pair_score)
        #print(list_pair_score);sys.exit()
        list_pair_score = list_pair_score.split('|')
        #Loop
        with open('./chr_mafft/2_score/2.0/' + region_name, 'a') as f_state:
           
            jj=0;match_str=[];frame_serial=0;query_serial=0
            while jj<mafft_result_len:
                state_mark=   list_state_mark[jj];              match_str.append(state_mark)
                pair_score=   list_pair_score[jj]
                frame_serial+=int(pair_score[0])
                query_serial+=int(pair_score[1])
                jj+=1
                #print('Progress: '+str(jj)+'/'+str(mafft_result_len),end='\r')
                if state_mark=='1':
                    f_state.write(f"{jj}\t{query_serial}\t{frame_serial}\n")
        with open('./chr_mafft/2_score/score1', 'a') as f_score:            
            match_str=''.join(match_str)
            f_score.write(f"{region_name}\t{region_pos}\t{frame_seq}\t{query_seq}\t{match_str}\n")
        
    #Create a process pool
    run_score_list=[];
    print('Loading ./chr_mafft/1_mafft/mafft_result_sum')
    with open('./chr_mafft/1_mafft/mafft_result_sum','r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=4:continue
            region_name=eachline_arr[0]
            #if region_name!='region_1':continue
            region_pos=eachline_arr[1]
            frame_seq=eachline_arr[2]
            query_seq=eachline_arr[3]
            eachline_tuple=(region_name,region_pos,frame_seq,query_seq)
            run_score_list.append(eachline_tuple)
            
    with open('./chr_mafft/2_score/score1','w') as f2:
        f2.write('')
    #Assign tasks to processes in the pool
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  
    with Pool(processes=thread) as pool:
        pool.map(run_score, run_score_list)   
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  
if argv1=="stepall" or argv1=="step2" or  argv1=="step2.1" or start_step<=22:##Estimated 73s
    print('step2.1 ———— Cumulative scoring')
    if  os.path.exists('./chr_mafft/2_score/2.1'):
        subprocess.run(["rm -r ./chr_mafft/2_score/2.1"], shell=True)
    subprocess.run(["mkdir ./chr_mafft/2_score/2.1"], shell=True)
    def run_score1(run_score_one):
        region_name,region_pos,frame_seq,query_seq,match_str=run_score_one
        base_all_num=len(match_str)
        print('Scoring')
        #Open file in write mode
        match_str=re.sub(r'[mid]', '0', match_str)
        with open('./chr_mafft/2_score/2.1/'+region_name, "w") as file:
            #Initialize variables
            current_block_char = None  #Current block character ('0' or '1')
            current_block_start = None  #Start index of current block
            #Iterate through the string
            for i, char in enumerate(match_str):
                #print('step2.1 Progress: '+str(i)+'/'+str(base_all_num),end='\r') #
                if current_block_char is None:
                    #Start first block
                    current_block_char = char
                    current_block_start = i
                elif char != current_block_char:
                    #Encounter a new block, write current block info to file
                    if current_block_char=='0':sign='-'
                    else:sign=''
                    file.write(f"{current_block_start+1}\t{sign}{i - current_block_start}\n")
                    #Start new block
                    current_block_char = char
                    current_block_start = i
         
            #Write the last block info to file (won't trigger at end of string)
            if current_block_char is not None:
                if current_block_char=='0':sign='-'
                else:sign=''
                file.write(f"{current_block_start+1}\t{sign}{i - current_block_start}\n")
        
        print('Formatting output')
        #Each line: serial, current match length, subsequent mismatch length
        one_accumulate_left_bad='';
        with open('./chr_mafft/2_score/2.1/'+region_name+'_b','a') as f3:
            try:
                with open('./chr_mafft/2_score/2.1/'+region_name,'r') as f2:
                    for line in f2.readlines():
                        eachline_arr=line.strip().split('\t')
                        serial=eachline_arr[0]
                        accumulate=eachline_arr[1]
                        if one_accumulate_left_bad=='':
                            if accumulate[0]=='-':  one_accumulate_left_bad=accumulate
                        else:
                            if accumulate[0]!='-':  one_accumulate=accumulate;  one_serial=serial
                            else:                   
                                one_accumulate_right_bad=accumulate
                                
                                f3.write(one_serial+'\t'+one_accumulate_left_bad+'\t'+one_accumulate+'\t'+one_accumulate_right_bad+'\n')
                                one_accumulate_left_bad=accumulate
            finally:
                #No need to explicitly close f3, with statement handles it automatically
                pass                    
                
     
    run_score2_list=[]         
    with open('./chr_mafft/2_score/score1','r') as f:
         for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=5:continue
            region_name=eachline_arr[0]
            region_pos=eachline_arr[1]
            frame_seq=eachline_arr[2]
            query_seq=eachline_arr[3]
            match_str=eachline_arr[4]
            eachline_tuple=(region_name,region_pos,frame_seq,query_seq,match_str)
            run_score2_list.append(eachline_tuple)
            
    #Assign tasks to processes in the pool
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  
    with Pool(processes=thread) as pool:
        pool.map(run_score1, run_score2_list)          
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))      
if argv1=="stepall" or argv1=="step2" or  argv1=="step2.2" or start_step<=23:
    print('step2.2 ———— Combining consecutive matches into larger blocks')
    if  os.path.exists('./chr_mafft/2_score/2.2'):
        subprocess.run(["rm -r ./chr_mafft/2_score/2.2"], shell=True)
    subprocess.run(["mkdir ./chr_mafft/2_score/2.2"], shell=True)
    def run_score2(run_score_one):
        region_name,region_pos=run_score_one
        #if region_name!='region_1':return ''
        #
        print('Reading dictionary ./chr_mafft/2_score/2.0/'+region_name)
        dict_mafftserial_queryserial={}
        with open('./chr_mafft/2_score/2.0/'+region_name,'r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                dict_mafftserial_queryserial[eachline_arr[0]]=eachline_arr[1]

        #
        line_serial=0
        ##################
        ###Initialize bases
        #window_size=30              #30bp with 3 mismatches is optimal, 21 matches is secondary
        good_percent=0.7; error_percent=1-good_percent; threshold=window_size*good_percent;  threshold_sub=window_size*error_percent
        #class2_percent=0.7; error2_percent=1-class2_percent; threshold2=window_size*class2_percent;  threshold2_sub=window_size*error2_percent
        print('Reading list ./chr_mafft/2_score/2.1/'+region_name+'_b')
        lineserial_info_list=[]
        with open('./chr_mafft/2_score/2.1/'+region_name+'_b','r') as f:
            for line in f.readlines():
                line_serial+=1
                eachline_arr=line.strip().split('\t')
                lineserial_info_list.append(eachline_arr)
        print("Processing sequentially")
        line_num=len(lineserial_info_list)
        i=0
        block_size=0;error_size=0
        while i<line_num:
            one_info=lineserial_info_list[i]
            one_bad_left_len=       int(one_info[1])
            one_good_len=           int(one_info[2])    
            one_bad_right_len=      int(one_info[3])
            ##Block initialization and update
            if block_size==0:
                if one_good_len+one_bad_right_len<1: i+=1;continue
                else:
                    block_size=one_good_len
                    block_start_index=i
                    block_start_serial=one_info[0]
            else:
                block_size+=(one_good_len-one_bad_left_len)
                error_size-=one_bad_left_len
            ##Block quality control and discard/output
            if block_size>window_size:
                one_percent=error_size/block_size
                if one_percent>error_percent:
                    i=block_start_index+1;block_size=0;error_size=0;continue
                else:
                    with open('./chr_mafft/2_score/2.2/'+region_name,'a') as f2:
                        f2.write(block_start_serial+'\t'+dict_mafftserial_queryserial[block_start_serial]+'\t'+str(block_size)+'\t'+str(round(1-one_percent,3))+'\n')
                    i=block_start_index+1;block_size=0;error_size=0;continue
            else:
                if error_size>threshold_sub: 
                    i=block_start_index+1;block_size=0;error_size=0;continue
                else:i+=1
                #print('Progress: '+str(i)+'/'+str(line_num),end='\r')
        
     
    run_score2_list=[]         
    with open('./chr_mafft/2_score/score1','r') as f:
         for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=5:continue
            region_name=eachline_arr[0]
            region_pos=eachline_arr[1]
            eachline_tuple=(region_name,region_pos)
            run_score2_list.append(eachline_tuple)
        
    #Assign tasks to processes in the pool
    with Pool(processes=thread) as pool:
        pool.map(run_score2, run_score2_list)     
if argv1=="stepall" or argv1=="step2" or  argv1=="step2.3" or start_step<=24:
    print('step2.3 ———— Filtering and summarizing large blocks')
    print('\tResult description:')
    print('\t\t\tregion_length: Length of the region')
    print('\t\t\tmatch_percent: Actual match percentage of the region')
    print('\t\t\tbig_block_size: Large block formed by combining smaller consecutive matches, includes small gaps in between')
    def run_score2(run_score_one):
        region_name,region_pos,frame_seq,query_seq,match_str=run_score_one
        #if region_name!='region_18':return ''
        #
        print('Reading dictionary ./chr_mafft/2_score/2.0/'+region_name)
        dict_mafftserial_queryserial={}
        with open('./chr_mafft/2_score/2.0/'+region_name,'r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                dict_mafftserial_queryserial[eachline_arr[0]]=eachline_arr[1]

        #print("Keeping only consecutive blocks within 3000bp")
        big_block_start1=''
        if  os.path.exists('./chr_mafft/2_score/2.2/'+region_name)==False: return False
        with open('./chr_mafft/2_score/2.2/'+region_name,'r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                one_block_start1=int(eachline_arr[0]);one_block_size=int(eachline_arr[2])
                if big_block_start1=='':
                    big_block_start1    =one_block_start1; 
                    big_block_end1      =one_block_start1+one_block_size-1
                    big_block_size      =one_block_size;
                    last_block_start1   =one_block_start1
                    last_block_size      =big_block_size
                    continue
                ##
                if one_block_start1-last_block_start1>1000:
                    if big_block_size>monomer_len*3:
                        mafft_region_length     =big_block_end1-big_block_start1+1
                        big_block_size_percent  =round(big_block_size/mafft_region_length,3)
                        region_match_num        =match_str[big_block_start1-1:big_block_end1].count('1')
                        region_match_percent    =round(region_match_num/mafft_region_length,3)
                        bigblock_chrstart       =dict_mafftserial_queryserial[str(big_block_start1)]
                        bigblock_chrend         =dict_mafftserial_queryserial[str(big_block_end1)]
                        chr_region_length       =int(bigblock_chrend)-int(bigblock_chrstart)+1
                        with open('./chr_mafft/2_score/score2_tmp','a') as f2:
                            f2.write(region_name+'\t'+region_pos+'\t'+
                                str(big_block_start1)+'\t'+ 
                                str(big_block_end1)+'\t'+   
                                bigblock_chrstart+'\t'+ 
                                bigblock_chrend+'\t'+                                 
                                str(mafft_region_length)+'\t'+str(chr_region_length)+'\t'+'plus'+'\t'+
                                str(big_block_size)+'\t'+str(big_block_size_percent)+'\t'+
                                str(region_match_num)+'\t'+str(region_match_percent)+'\n')                     
                    big_block_start1=one_block_start1
                    big_block_end1=big_block_start1+one_block_size-1
                    big_block_size=one_block_size
                else:
                    #
                    big_block_end1 = one_block_start1+one_block_size-1
                    #
                    last_block_end    = last_block_start1+last_block_size-1
                    last_block_end_shouldmax = one_block_start1-1                    
                    if last_block_end > last_block_end_shouldmax:  ##Most blocks exceed the initial row
                        big_block_size_delta = one_block_size-(last_block_end-last_block_end_shouldmax)
                    else:
                        big_block_size_delta = one_block_size
                    big_block_size+=big_block_size_delta
                    #print('big_block_size_delta:'+str(big_block_size_delta))
                last_block_start1 = one_block_start1
                last_block_size   = one_block_size
     
    run_score2_list=[]         
    with open('./chr_mafft/2_score/score1','r') as f:
         for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=5:continue
            region_name=eachline_arr[0]
            region_pos=eachline_arr[1]
            frame_seq=eachline_arr[2]
            query_seq=eachline_arr[3]
            match_str=eachline_arr[4]
            eachline_tuple=(region_name,region_pos,frame_seq,query_seq,match_str)
            run_score2_list.append(eachline_tuple)
    with open('./chr_mafft/2_score/score2_tmp','w') as f2:
        f2.write('region_name\tregion_pos\tbigblock_mafftstart\tbigblock_mafftend\tbigblock_chrstart\tbigblock_chrend\tmafft_region_length\tchr_region_length\tstrand\tbig_block_size\tbig_block_size_percent\tmatch_num\tmatch_percent\n')            
    #Assign tasks to processes in the pool
    with Pool(processes=thread) as pool:
        pool.map(run_score2, run_score2_list)     
    subprocess.run(["( head -n 1 ./chr_mafft/2_score/score2_tmp && tail -n +2  ./chr_mafft/2_score/score2_tmp | sort -k1,1V -k4,4n ) > ./chr_mafft/2_score/score2;   rm ./chr_mafft/2_score/score2_tmp"], shell=True)
#step2r        
if argv1=="stepall" or argv1=="step2r" or  argv1=="step2.0r" or start_step<=25:
    print('step2.0r ———— Scoring MAFFT alignment results (reverse strand)')
    if  os.path.exists('./chr_mafft/2_score_r')==False:
        subprocess.run(["mkdir ./chr_mafft/2_score_r"], shell=True)
    if  os.path.exists('./chr_mafft/2_score_r/2.0'):
        subprocess.run(["rm -r ./chr_mafft/2_score_r/2.0"], shell=True)
    subprocess.run(["mkdir ./chr_mafft/2_score_r/2.0"], shell=True)
    print('Storing dict: mafft_serial + query_serial + frame_serial')
    def run_score(run_score_one):
        region_name,region_pos,frame_seq,query_seq=run_score_one
        mafft_result_len=len(frame_seq)
        print('Using list comprehension to create new list')
        pair_list = [x + y for x, y in zip(frame_seq, query_seq)]
        pair_list_str="|".join(pair_list);      
        #print('Applying substitution rules 1')
        list_state_mark = re.sub(r'AA|CC|GG|TT', '1', pair_list_str)
        list_state_mark = re.sub(r'[^-|]-', 'd', list_state_mark)
        list_state_mark = re.sub(r'-[^-|]', 'i', list_state_mark)        
        list_state_mark = re.sub(r'[^-|][^-|]', 'm', list_state_mark)
        list_state_mark = list_state_mark.split('|')
        #print('Applying substitution rules 2')
        list_pair_score = re.sub(r'[^-|]-', '10', pair_list_str)
        list_pair_score = re.sub(r'-[^-|]', '01', list_pair_score)        
        list_pair_score = re.sub(r'[^-|0][^-|0]', '11', list_pair_score)
        list_pair_score = list_pair_score.split('|')
        #Loop
        with open('./chr_mafft/2_score_r/2.0/' + region_name, 'a') as f_state:
           
            jj=0;match_str=[];frame_serial=0;query_serial=0
            while jj<mafft_result_len:
                state_mark=   list_state_mark[jj];              match_str.append(state_mark)
                pair_score=   list_pair_score[jj]
                frame_serial+=int(pair_score[0])
                query_serial+=int(pair_score[1])
                jj+=1
                #print('Progress: '+str(jj)+'/'+str(mafft_result_len),end='\r')
                if state_mark=='1':
                    f_state.write(f"{jj}\t{query_serial}\t{frame_serial}\n")
        with open('./chr_mafft/2_score_r/score1', 'a') as f_score:            
            match_str=''.join(match_str)
            f_score.write(f"{region_name}\t{region_pos}\t{frame_seq}\t{query_seq}\t{match_str}\n")
        
    #Create a process pool
    run_score_list=[];
    print('Loading ./chr_mafft/1_mafft_r/mafft_result_sum')
    with open('./chr_mafft/1_mafft_r/mafft_result_sum','r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=4:continue
            region_name=eachline_arr[0]
            region_pos=eachline_arr[1]
            frame_seq=eachline_arr[2]
            query_seq=eachline_arr[3]
            eachline_tuple=(region_name,region_pos,frame_seq,query_seq)
            run_score_list.append(eachline_tuple)
            
    with open('./chr_mafft/2_score_r/score1','w') as f2:
        f2.write('')
    #Assign tasks to processes in the pool
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  
    with Pool(processes=thread) as pool:
        pool.map(run_score, run_score_list)   
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))   
if argv1=="stepall" or argv1=="step2r" or  argv1=="step2.1r" or start_step<=27:##Was too slow (total time 1.5h), need a different approach
    print('step2.1r ———— Cumulative scoring (reverse strand)')
    if  os.path.exists('./chr_mafft/2_score_r/2.1'):
        subprocess.run(["rm -r ./chr_mafft/2_score_r/2.1"], shell=True)
    subprocess.run(["mkdir ./chr_mafft/2_score_r/2.1"], shell=True)
    def run_score1(run_score_one):
        region_name,region_pos,frame_seq,query_seq,match_str=run_score_one
        base_all_num=len(match_str)
        #Open file in write mode
        match_str=re.sub(r'[mid]', '0', match_str)
        with open('./chr_mafft/2_score_r/2.1/'+region_name, "w") as file:
            #Initialize variables
            current_block_char = None  #Current block character ('0' or '1')
            current_block_start = None  #Start index of current block
            #Iterate through the string
            for i, char in enumerate(match_str):
                #print('step2.1 Progress: '+str(i)+'/'+str(base_all_num),end='\r') #
                if current_block_char is None:
                    #Start first block
                    current_block_char = char
                    current_block_start = i
                elif char != current_block_char:
                    #Encounter a new block, write current block info to file
                    if current_block_char=='0':sign='-'
                    else:sign=''
                    file.write(f"{current_block_start+1}\t{sign}{i - current_block_start}\n")
                    #Start new block
                    current_block_char = char
                    current_block_start = i
         
            #Write the last block info to file (won't trigger at end of string)
            if current_block_char is not None:
                if current_block_char=='0':sign='-'
                else:sign=''
                file.write(f"{current_block_start+1}\t{sign}{i - current_block_start}\n")
        
        print('Formatting output')
        #Each line: serial, current match length, subsequent mismatch length
        one_accumulate_left_bad='';
        with open('./chr_mafft/2_score_r/2.1/'+region_name+'_b','a') as f3:
            try:
                with open('./chr_mafft/2_score_r/2.1/'+region_name,'r') as f2:
                    for line in f2.readlines():
                        eachline_arr=line.strip().split('\t')
                        serial=eachline_arr[0]
                        accumulate=eachline_arr[1]
                        if one_accumulate_left_bad=='':
                            if accumulate[0]=='-':  one_accumulate_left_bad=accumulate
                        else:
                            if accumulate[0]!='-':  one_accumulate=accumulate;  one_serial=serial
                            else:                   
                                one_accumulate_right_bad=accumulate
                                
                                f3.write(one_serial+'\t'+one_accumulate_left_bad+'\t'+one_accumulate+'\t'+one_accumulate_right_bad+'\n')
                                one_accumulate_left_bad=accumulate
            finally:
                #No need to explicitly close f3, with statement handles it automatically
                pass                    
                
     
    run_score2_list=[]         
    with open('./chr_mafft/2_score_r/score1','r') as f:
         for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=5:continue
            region_name=eachline_arr[0]
            region_pos=eachline_arr[1]
            frame_seq=eachline_arr[2]
            query_seq=eachline_arr[3]
            match_str=eachline_arr[4]
            eachline_tuple=(region_name,region_pos,frame_seq,query_seq,match_str)
            run_score2_list.append(eachline_tuple)
            
    #Assign tasks to processes in the pool
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  
    with Pool(processes=thread) as pool:
        pool.map(run_score1, run_score2_list)          
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))          
if argv1=="stepall" or argv1=="step2r" or argv1=="step2.2r" or start_step<=28:
    print('step2.2 ———— Combine multiple consecutive matches into large blocks')
    if os.path.exists('./chr_mafft/2_score_r/2.2'):
        subprocess.run(["rm -r ./chr_mafft/2_score_r/2.2"], shell=True)
    subprocess.run(["mkdir ./chr_mafft/2_score_r/2.2"], shell=True)
    def run_score2(run_score_one):
        region_name,region_pos=run_score_one
        #
        print('Reading dict ./chr_mafft/2_score_r/2.0/'+region_name)
        dict_mafftserial_queryserial={}
        with open('./chr_mafft/2_score_r/2.0/'+region_name,'r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                dict_mafftserial_queryserial[eachline_arr[0]]=eachline_arr[1]

        #
        line_serial=0
        ##################
        ### Initialize bases
        #window_size=30              # 30bp with 3 different bases is optimal, 30bp with 21 identical bases is suboptimal
        good_percent=0.7; error_percent=1-good_percent; threshold=window_size*good_percent;  threshold_sub=window_size*error_percent
        #class2_percent=0.7; error2_percent=1-class2_percent; threshold2=window_size*class2_percent;  threshold2_sub=window_size*error2_percent
        print('Reading list ./chr_mafft/2_score_r/2.1/'+region_name+'_b')
        lineserial_info_list=[]
        with open('./chr_mafft/2_score_r/2.1/'+region_name+'_b','r') as f:
            for line in f.readlines():
                line_serial+=1
                eachline_arr=line.strip().split('\t')
                lineserial_info_list.append(eachline_arr)
        print("Start sequential judgment")
        line_num=len(lineserial_info_list)
        i=0
        block_size=0;error_size=0
        while i<line_num:
            one_info=lineserial_info_list[i]
            one_bad_left_len=       int(one_info[1])
            one_good_len=           int(one_info[2])    
            one_bad_right_len=      int(one_info[3])
            ## Block initialization and update
            if block_size==0:
                if one_good_len+one_bad_right_len<1: i+=1;continue
                else:
                    block_size=one_good_len
                    block_start_index=i
                    block_start_serial=one_info[0]
            else:
                block_size+=(one_good_len-one_bad_left_len)
                error_size-=one_bad_left_len
            ## Block quality check and discard/output
            if block_size>window_size:
                one_percent=error_size/block_size
                if one_percent>error_percent:
                    i=block_start_index+1;block_size=0;error_size=0;continue
                else:
                    with open('./chr_mafft/2_score_r/2.2/'+region_name,'a') as f2:
                        f2.write(block_start_serial+'\t'+dict_mafftserial_queryserial[block_start_serial]+'\t'+str(block_size)+'\t'+str(round(1-one_percent,3))+'\n')
                    i=block_start_index+1;block_size=0;error_size=0;continue
            else:
                if error_size>threshold_sub: 
                    i=block_start_index+1;block_size=0;error_size=0;continue
                else:i+=1
                #print('Progress: '+str(i)+'/'+str(line_num),end='\r')
        
     
    run_score2_list=[]         
    with open('./chr_mafft/2_score_r/score1','r') as f:
         for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=5:continue
            region_name=eachline_arr[0]
            region_pos=eachline_arr[1]
            eachline_tuple=(region_name,region_pos)
            run_score2_list.append(eachline_tuple)
        
    # Assign tasks to processes in the process pool
    with Pool(processes=thread) as pool:
        pool.map(run_score2, run_score2_list)     
if argv1=="stepall" or argv1=="step2r" or argv1=="step2.3r" or start_step<=29:
    print('step2.3r ———— Filter and count large blocks')
    print('\tResult description:')
    print('\t\t\tregion_length: length of the region')
    print('\t\t\tmatch_percent: actual match percentage of the region')
    print('\t\t\tbig_block_size: large match block formed by small consecutive matches, including small internal gaps')
    def run_score2(run_score_one):
        region_name,region_pos,frame_seq,query_seq,match_str=run_score_one
        #
        print('Reading dict ./chr_mafft/2_score_r/2.0/'+region_name)
        dict_mafftserial_queryserial={}
        with open('./chr_mafft/2_score_r/2.0/'+region_name,'r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                dict_mafftserial_queryserial[eachline_arr[0]]=eachline_arr[1]

        #print("Keep only consecutive blocks within 3000bp")
        big_block_start1=''
        if os.path.exists('./chr_mafft/2_score_r/2.2/'+region_name)==False: return False
        with open('./chr_mafft/2_score_r/2.2/'+region_name,'r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                one_block_start1=int(eachline_arr[0]);one_block_size=int(eachline_arr[2])
                if big_block_start1=='':
                    big_block_start1    =one_block_start1; 
                    big_block_end1      =one_block_start1+one_block_size-1
                    big_block_size      =one_block_size;
                    last_block_start1   =one_block_start1
                    last_block_size      =big_block_size
                    continue
                ##
                if one_block_start1-last_block_start1>1000:
                    if big_block_size>monomer_len*3:
                        mafft_region_length     =big_block_end1-big_block_start1+1
                        big_block_size_percent  =round(big_block_size/mafft_region_length,3)
                        region_match_num        =match_str[big_block_start1-1:big_block_end1].count('1')
                        region_match_percent    =round(region_match_num/mafft_region_length,3)
                        bigblock_chrstart       =dict_mafftserial_queryserial[str(big_block_start1)]
                        bigblock_chrend         =dict_mafftserial_queryserial[str(big_block_end1)]
                        chr_region_length       =int(bigblock_chrend)-int(bigblock_chrstart)+1
                        with open('./chr_mafft/2_score_r/score2_tmp','a') as f2:
                            f2.write(region_name+'\t'+region_pos+'\t'+
                                str(big_block_start1)+'\t'+ 
                                str(big_block_end1)+'\t'+   
                                bigblock_chrstart+'\t'+ 
                                bigblock_chrend+'\t'+                                 
                                str(mafft_region_length)+'\t'+str(chr_region_length)+'\t'+'minus'+'\t'+
                                str(big_block_size)+'\t'+str(big_block_size_percent)+'\t'+
                                str(region_match_num)+'\t'+str(region_match_percent)+'\n')                     
                    big_block_start1=one_block_start1
                    big_block_end1=big_block_start1+one_block_size-1
                    big_block_size=one_block_size
                else:
                    #
                    big_block_end1 = one_block_start1+one_block_size-1
                    #
                    last_block_end    = last_block_start1+last_block_size-1
                    last_block_end_shouldmax = one_block_start1-1                    
                    if last_block_end > last_block_end_shouldmax:  ## Most blocks exceed the starting line
                        big_block_size_delta = one_block_size-(last_block_end-last_block_end_shouldmax)
                    else:
                        big_block_size_delta = one_block_size
                    big_block_size+=big_block_size_delta
                    #print('big_block_size_delta:'+str(big_block_size_delta))
                last_block_start1 = one_block_start1
                last_block_size   = one_block_size
     
    run_score2_list=[]         
    with open('./chr_mafft/2_score_r/score1','r') as f:
         for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=5:continue
            region_name=eachline_arr[0]
            region_pos=eachline_arr[1]
            frame_seq=eachline_arr[2]
            query_seq=eachline_arr[3]
            match_str=eachline_arr[4]
            eachline_tuple=(region_name,region_pos,frame_seq,query_seq,match_str)
            run_score2_list.append(eachline_tuple)
    with open('./chr_mafft/2_score_r/score2_tmp','w') as f2:
        f2.write('region_name\tregion_pos\tbigblock_mafftstart\tbigblock_mafftend\tbigblock_chrstart\tbigblock_chrend\tmafft_region_length\tchr_region_length\tstrand\tbig_block_size\tbig_block_size_percent\tmatch_num\tmatch_percent\n')            
    # Assign tasks to processes in the process pool
    with Pool(processes=thread) as pool:
        pool.map(run_score2, run_score2_list)     
    subprocess.run(["( head -n 1 ./chr_mafft/2_score_r/score2_tmp && tail -n +2  ./chr_mafft/2_score_r/score2_tmp | sort -k1,1V -k4,4n ) > ./chr_mafft/2_score_r/score2;   rm ./chr_mafft/2_score_r/score2_tmp"], shell=True)

#step2s, summarize positive and negative strands
if argv1=="stepall" or argv1=="step3" or argv1=="step2s" or start_step<=30:
    print('Summarize the results of step2')
    if os.path.exists('./chr_mafft/2_score/score2')==False: print('Need to execute step2'); sys.exit()
    if os.path.exists('./chr_mafft/2_score_r/score2')==False: print('Need to execute step2r'); sys.exit()
    big_block_size_threshold=300
    big_block_size_percent_threshold=0.5        ### Sum the lengths of all blocks greater than 70%
    match_num_threshold=300
    match_percent_threshold=0.5                 ### Regardless of the overall match ratio of the block, both this and block length need to be satisfied simultaneously
    cmd=f'''
# Extract header (assumed to be from file A)
header=$(head -n 1 ./chr_mafft/2_score/score2)
 
# Merge data portion (skip headers of files A and B)
data=$(tail -n +2 ./chr_mafft/2_score/score2 && tail -n +2 ./chr_mafft/2_score_r/score2)
 
# Sort the data portion
sorted_data=$(echo "$data" | sort -k1,1V -k5,5n)

# Filter data (using awk)
filtered_data=$(echo "$sorted_data" | awk '$10 > {big_block_size_threshold} && $11 > {big_block_size_percent_threshold} && $12 > {match_num_threshold} && $13 > {match_percent_threshold}')

# Add header before the sorted data and output to the final result file
echo "$header" > ./chr_mafft/2_good_regions_tmp
echo "$filtered_data" >> ./chr_mafft/2_good_regions_tmp'''
    subprocess.run([cmd], shell=True)     
    with open('./chr_mafft/2_good_regions','w') as f:
        f.write('')
    last_1=''
    with open('./chr_mafft/2_good_regions_tmp','r') as f:
        for line in f.readlines():
            eachline=line.strip()
            eachline_arr=eachline.split('\t')
            if last_1=='':
                with open('./chr_mafft/2_good_regions','a') as f:
                    f.write(eachline+'\n')
            elif eachline_arr[0]!=last_1:
                with open('./chr_mafft/2_good_regions','a') as f:
                    f.write('\n'+eachline+'\n')
            else:
                with open('./chr_mafft/2_good_regions','a') as f:
                    f.write(eachline+'\n') 
            last_1=eachline_arr[0]
    subprocess.run(["rm ./chr_mafft/2_good_regions_tmp"], shell=True)        
    #subprocess.run(["cp ./chr_mafft/2_good_regions ./chr_mafft/2_good_regions.bak"], shell=True)        
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    
#step3, SEAT - Place each base into seat positions, Estimated 120s
if argv1=="stepall" or argv1=="step3"  or  argv1=="step3.0" or start_step<=31:
    print('step3 ———— Analyzing base positions (SEAT)')##3.0 Estimated 10s
    print('step3.0 ———— Extracting operational sequences')
    if  os.path.exists('./chr_mafft/3_seat')==False:
        subprocess.run(["mkdir ./chr_mafft/3_seat"], shell=True)
    if  os.path.exists('./chr_mafft/3_seat/3.0'):
        subprocess.run(["rm -r ./chr_mafft/3_seat/3.0"], shell=True)   
    subprocess.run(["mkdir ./chr_mafft/3_seat/3.0"], shell=True)
    ####
    print('step3.0 Loading dictionaries from ./chr_mafft/2_score (plus and minus strands)')
    dict_plus_region_info={}
    with open('./chr_mafft/2_score/score1','r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=5:continue
            dict_plus_region_info[eachline_arr[0]]=eachline_arr[2:]   #region_pos+'\t'+frame_seq+'\t'+query_seq+'\t'+match_str+
    dict_minus_region_info={}
    with open('./chr_mafft/2_score_r/score1','r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=5:continue
            dict_minus_region_info[eachline_arr[0]]=eachline_arr[2:]             
    ####
    print('step3.0 Extracting step2 results and truncating operational sequences')
    i=0
    with open('./chr_mafft/2_good_regions','r') as f:
        for line in f.readlines():
            i+=1
            if i==1:continue
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=13:
                #print('Row does not have 13 columns: region_name\tregion_pos\tbigblock_mafftstart\tbigblock_mafftend\tbigblock_chrstart\tbigblock_chrend\tmafft_region_length\tchr_region_length\tstrand\tbig_block_size\tbig_block_size_percent\tmatch_num\tmatch_percent\n') #
                #print('Error row number:'+str(i))
                #print(eachline_arr)
                continue
            region_name         =eachline_arr[0]
            bigblock_mafftstart =eachline_arr[2];   one_mafftstart_index=   int(bigblock_mafftstart)-1
            bigblock_mafftend   =eachline_arr[3];   one_mafftend_index=     int(bigblock_mafftend)-1
            bigblock_chrstart   =eachline_arr[4]
            bigblock_chrend     =eachline_arr[5]
            strand              =eachline_arr[8]
            if strand=='plus':  region_info=dict_plus_region_info[region_name];        sign='+'
            if strand=='minus': region_info=dict_minus_region_info[region_name];      sign='-'
            frame_seq=region_info[0];     
            query_seq=region_info[1];
            match_str=region_info[2];    
            one_output_name=region_name+sign+':'+bigblock_chrstart+'-'+bigblock_chrend
            one_frame_seq=frame_seq[one_mafftstart_index:one_mafftend_index+1]
            one_query_seq=query_seq[one_mafftstart_index:one_mafftend_index+1]
            one_match_seq=match_str[one_mafftstart_index:one_mafftend_index+1]
            if one_match_seq[0]!='1':  print('Error 1:'+one_output_name)
            if one_match_seq[-1]!='1': print('Error 2:'+one_output_name)
            with open('./chr_mafft/3_seat/3.0/'+one_output_name,'w') as f2:
                f2.write('>'+one_output_name+'_frame\n')
                f2.write(one_frame_seq+'\n')
                f2.write('>'+one_output_name+'_query\n')
                f2.write(one_query_seq+'\n')
                f2.write('>'+one_output_name+'_match\n')
                f2.write(one_match_seq+'\n')
if argv1=="stepall" or argv1=="step3"  or  argv1=="step3.1" or start_step<=33:
    print('step3.1 ———— Outputting positions')  #Estimated 42s
    if  os.path.exists('./chr_mafft/3_seat/3.1'):
        subprocess.run(["rm -r ./chr_mafft/3_seat/3.1"], shell=True)   
    subprocess.run(["mkdir ./chr_mafft/3_seat/3.1"], shell=True)
    run_seat_list=[]
    i=0
    print('Reading ./chr_mafft/2_good_regions')
    with open('./chr_mafft/2_good_regions','r') as f:
        for line in f.readlines():
            i+=1
            if i==1:continue
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=13:
                #print('Row does not have 13 columns: region_name\tregion_pos\tbigblock_mafftstart\tbigblock_mafftend\tbigblock_chrstart\tbigblock_chrend\tmafft_region_length\tchr_region_length\tstrand\tbig_block_size\tbig_block_size_percent\tmatch_num\tmatch_percent\n') #
                #print('Error row number:'+str(i))
                #print(eachline_arr)
                continue 
            run_seat_list.append(eachline_arr) 
    #print('Main task started\n')        
    #####
    def run_seat(run_seat_list_one):
        #print(run_seat_list_one)
        ##################        
        region_name=        run_seat_list_one[0]
        strand=             run_seat_list_one[8];#print (strand)
        if strand=='plus':  sign='+'
        else:               sign='-'
        #####Load dictionary
        dict_mafftserial_frameserial={}
        if strand=='plus':  dict_file='./chr_mafft/2_score/2.0/'+region_name
        elif strand=='minus':  dict_file='./chr_mafft/2_score_r/2.0/'+region_name
        else: print('error')
        print('3.1 Processing '+dict_file)
        #print('Loading dictionary '+dict_file)
        with open(dict_file,'r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=3:continue
                dict_mafftserial_frameserial[eachline_arr[0]]=eachline_arr[2]
        #######Read sequences   
        bigblock_mafftstart=  run_seat_list_one[2]
        bigblock_mafftend=  run_seat_list_one[3]
        bigblock_chrstart=  run_seat_list_one[4];        bigblock_frame_start=   dict_mafftserial_frameserial[bigblock_mafftstart]
        bigblock_chrend=    run_seat_list_one[5];        bigblock_frame_end=     dict_mafftserial_frameserial[bigblock_mafftend]
        input_name=         region_name+sign+':'+bigblock_chrstart+'-'+bigblock_chrend 
        input_file=         './chr_mafft/3_seat/3.0/'+input_name
        output_file=         './chr_mafft/3_seat/3.1/'+input_name        
        #print('Loading sequences')
        with open (input_file,'r') as f:
            for line in f.readlines():
                eachline=line.strip()
                if eachline=='>'+input_name+'_frame': state='frame';continue
                if eachline=='>'+input_name+'_query': state='query';continue
                if eachline=='>'+input_name+'_match': state='match';continue
                if state=='frame':  frame=eachline
                if state=='query':  query=eachline
                if state=='match':  match=eachline
        ########### 
        #print('Processing sequences')
        mafft_len=len(match)
        bigblock_frame_start_seat=  int(bigblock_frame_start)%monomer_len
        if bigblock_frame_start_seat==0: bigblock_frame_start_seat=monomer_len
        i=0; insert_acc_bases='';deletion_acc_bases=''
        with open (output_file,'a') as f:
            while i<mafft_len:
                one_frame_pos = frame[i]
                one_query_pos = query[i]
                one_match_pos = match[i]
                #mafft_serial=i+1
                ###Determine seat
                if i==0:  
                    seat=bigblock_frame_start_seat
                    chr_serial=int(bigblock_chrstart)
                elif one_frame_pos!='-'  :
                    seat+=1
                    if seat>monomer_len:seat=1
                    if one_query_pos!='-':chr_serial+=1
                ###Output only for real seat positions?
                if      one_match_pos=='i':     insert_acc_bases+=one_query_pos
                elif    one_match_pos=='d':     deletion_acc_bases+=one_frame_pos  ##deletion seems meaningless
                else:
                    f.write(str(chr_serial)+'\t'+one_query_pos+'\t'+one_frame_pos+'\t'+insert_acc_bases+'\t'+str(seat)+'\n')
                    insert_acc_bases=''    
                    deletion_acc_bases=''
                i+=1
        #print('Completed: '+input_name)        
            #print('Progress: '+str(i)+'/'+str(mafft_len),end='\r')
        
        #####
        
                    
    #Assign tasks to processes in the pool
    with Pool(processes=20) as pool:
        pool.map(run_seat, run_seat_list) 
if argv1=="stepall" or argv1=="step3"  or  argv1=="step3.2" or start_step<=34:
    print('step3.2 ———— Format transformation')  ##Estimated 8s
    if  os.path.exists('./chr_mafft/3_seat/3.2'):
        subprocess.run(["rm -r ./chr_mafft/3_seat/3.2"], shell=True)   
    subprocess.run(["mkdir ./chr_mafft/3_seat/3.2"], shell=True)
    run_seat_list=[]
    i=0
    print('step3.2 Reading ./chr_mafft/2_good_regions')
    with open('./chr_mafft/2_good_regions','r') as f:
        for line in f.readlines():
            i+=1
            if i==1:continue
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=13:
                #print('Row does not have 13 columns: region_name\tregion_pos\tbigblock_mafftstart\tbigblock_mafftend\tbigblock_chrstart\tbigblock_chrend\tmafft_region_length\tchr_region_length\tstrand\tbig_block_size\tbig_block_size_percent\tmatch_num\tmatch_percent\n') #
                #print('Error row number:'+str(i))
                #print(eachline_arr)
                continue 
            run_seat_list.append(eachline_arr) 
    #print('Main task started\n')        
    #####
    def run_seat(run_seat_list_one):
        #print(run_seat_list_one)
        ##################        
        region_name=        run_seat_list_one[0]
        
        strand=             run_seat_list_one[8];#print (strand)
        if strand=='plus':  sign='+'
        else:               sign='-'
        bigblock_mafftstart=  run_seat_list_one[2]
        bigblock_mafftend=  run_seat_list_one[3]
        bigblock_chrstart=  run_seat_list_one[4];        
        bigblock_chrend=    run_seat_list_one[5];        
        input_name=         region_name+sign+':'+bigblock_chrstart+'-'+bigblock_chrend 
        #if input_name!='region_4+:12565351-13295334':return ''
        input_file=         './chr_mafft/3_seat/3.1/'+input_name
        output_file=        './chr_mafft/3_seat/3.2/'+input_name     
        ### ##   
        dict_linenum_info={}
        
        i=0
        with open (input_file,'r') as f:
            for line in f.readlines():
                eachline_arr        =line.strip().split('\t');#print(eachline_arr)
                #chr_serial          =eachline_arr[0]
                #one_base            =eachline_arr[1]
                #seat_raw            =eachline_arr[3]   
                #insert_acc_bases    =eachline_arr[4]
                #if len(eachline_arr)!=4 and len(eachline_arr)!=5:continue
                dict_linenum_info[i]=eachline_arr
                i+=1
        line_num=i 
        #print(dict_linenum_info);sys.exit()
        ###################
        with open(output_file,'a') as f:
            i=1
            while i<line_num-1:
                #print('Progress: '+str(i)+'/'+str(line_num),end='\r')
                info1=dict_linenum_info[i-1]
                info2=dict_linenum_info[i]
                info3=dict_linenum_info[i+1]
                chr_serial_1=info1[0]
                chr_serial_2=info2[0]
                chr_serial_3=info3[0]
                base_1=info1[1]
                base_2=info2[1]
                base_3=info3[1]
                mid_1=info2[3]
                mid_2=info3[3]            
                seat_1=info1[4]
                seat_2=info2[4]
                seat_3=info3[4]
                i+=1
                oneline=    chr_serial_1+'//'+chr_serial_2+'//'+chr_serial_3+'\t'+  \
                            base_1+'/'+mid_1+'/'+base_2+'/'+mid_2+'/'+base_3+'\t'+  \
                            seat_1+'//'+seat_2+'//'+seat_3+'\n'
                f.write(oneline)
        print('step3.2 Completed: '+input_name)            
    #Assign tasks to processes in the pool
    with Pool(processes=thread) as pool:
        pool.map(run_seat, run_seat_list) 
if argv1=="stepall" or argv1=="step3"  or  argv1=="step3.3" or start_step<=35:
    print('step3.3 ———— Statistics (combining plus and minus strands)')  #Estimated 20s
    if  os.path.exists('./chr_mafft/3_seat/3.3'):
        subprocess.run(["rm -r ./chr_mafft/3_seat/3.3"], shell=True)   
    subprocess.run(["mkdir ./chr_mafft/3_seat/3.3"], shell=True)    
    subprocess.run(["mkdir ./chr_mafft/3_seat/3.3/stat1"], shell=True)
    subprocess.run(["mkdir ./chr_mafft/3_seat/3.3/stat2"], shell=True)
    i=0
    print('step3.3 Reading ./chr_mafft/2_good_regions and files in ./chr_mafft/3_seat/3.2/')
    
    
    run_seat_list=[]
    with open('./chr_mafft/2_good_regions','r') as f:
        for line in f.readlines():
            i+=1
            if i==1:continue
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=13:
                continue 
            run_seat_list.append(eachline_arr) 
            
    def run_seat(run_seat_list_one):
        #print(run_seat_list_one)
        ##################        
        region_name=        run_seat_list_one[0]
        
        strand=             run_seat_list_one[8];#print (strand)
        if strand=='plus':  sign='+'
        else:               sign='-'
        bigblock_mafftstart=  run_seat_list_one[2]
        bigblock_mafftend=  run_seat_list_one[3]
        bigblock_chrstart=  run_seat_list_one[4];        
        bigblock_chrend=    run_seat_list_one[5];        
        input_name=         region_name+sign+':'+bigblock_chrstart+'-'+bigblock_chrend 
        #if input_name!='region_4+:12565351-13295334':return ''
        input_file=         './chr_mafft/3_seat/3.2/'+input_name
        output_file=        './chr_mafft/3_seat/3.3/'+input_name            
        
        print(input_file)
        dict_seatpath_base={}
        dict_seat3_base={}       
        with open(output_file,'a') as f3:
            with open(input_file,'r') as f2:
                for line in f2.readlines():
                    eachline_arr=line.strip().split('\t')
                    #15059786//15059787//15059788	G/./A/./A	15//16//17
                    #15059787//15059788//15059789	A/./A/./T	16//17//18
                    base_info=eachline_arr[1]
                    seat_info=eachline_arr[2]
                    base_info_arr=base_info.split('/')
                    seat_info_arr=seat_info.split('//')
                    ###
                    #print('eachline_arr:')
                    #print(eachline_arr)
                    if sign=='+':
                        base_info_arr_new=base_info_arr
                        seat_info_arr_new=seat_info_arr
                    else:
                        base_info_arr_r=base_info_arr;base_info_arr_r.reverse()
                        base_info_arr_new = []
                        for one_seq in base_info_arr_r:
                            if one_seq=='.': new_seq='.'
                            else:
                                one_seq_r=one_seq[::-1]
                                replacement_map = {
                                    'A': 'T',
                                    'T': 'A',
                                    'C': 'G',
                                    'G': 'C',
                                    'N':'N'
                                }
                                new_seq = ''.join(replacement_map.get(char, 'N')  for char in one_seq_r)   #replacement_map[char] for char in one_seq_r
                            base_info_arr_new.append(new_seq)
                        
                        ###
                        seat_info_arr_r=seat_info_arr;seat_info_arr_r.reverse()
                        seat_info_arr_new = []
                        
                        for element in seat_info_arr_r:
                            new_element = str(monomer_len-int(element)+1)
                            seat_info_arr_new.append(new_element)                            
                    base_info='/'.join(base_info_arr_new) 
                    seat_info='//'.join(seat_info_arr_new) 
                    ###
                    #First row
                    if dict_seatpath_base=={}:
                        dict_seatpath_base[seat_info_arr_new[0]]={}
                        dict_seatpath_base[seat_info_arr_new[0]][base_info_arr_new[0]]=1
                        dict_seatpath_base[seat_info_arr_new[1]]={}
                        dict_seatpath_base[seat_info_arr_new[1]][base_info_arr_new[2]]=1
                        if base_info_arr_new[1]!='.':
                            seatpath=seat_info_arr_new[0]+'-'+seat_info_arr_new[1]
                            dict_seatpath_base[seatpath]={}
                            dict_seatpath_base[seatpath][base_info_arr_new[1]]=1
                    #All rows
                    if seat_info_arr_new[2] not in dict_seatpath_base:
                        dict_seatpath_base[seat_info_arr_new[2]]={}
                    if base_info_arr_new[4] not in dict_seatpath_base[seat_info_arr_new[2]]:
                        dict_seatpath_base[seat_info_arr_new[2]][base_info_arr_new[4]]=0
                    dict_seatpath_base[seat_info_arr_new[2]][base_info_arr_new[4]]+=1
                    if  base_info_arr_new[3]!='.':
                        seatpath=seat_info_arr_new[1]+'-'+seat_info_arr_new[2]
                        if seatpath not in dict_seatpath_base:
                            dict_seatpath_base[seatpath]={}
                        if base_info_arr_new[3]not in dict_seatpath_base[seatpath]:
                            dict_seatpath_base[seatpath][base_info_arr_new[3]]=0
                        dict_seatpath_base[seatpath][base_info_arr_new[3]]+=1
                    ###
                    
                    f3.write('\t'.join(eachline_arr)+'\t'+sign+"\t"+base_info+'\t'+seat_info+'\n')
                    ###    
                    if seat_info not in dict_seat3_base:
                        dict_seat3_base[seat_info]={}
                    if base_info not in dict_seat3_base[seat_info]:
                        dict_seat3_base[seat_info][base_info]=0
                    dict_seat3_base[seat_info][base_info]+=1
                
        ############################################################################################
        ############################################################################################        
        seat1_list=list(dict_seatpath_base.keys())
        #Define a function to extract the number before the first '-', use sorted() with key parameter for sorting
        def custom_sort(item):
            key_part = item.split('\t')[0]
            if '-' in key_part:
                left, right = key_part.split('-')
                return (int(left), int(right))
            else:
                return (int(key_part),)
         
        #Sort the data
        sorted_seatpath_list = sorted(seat1_list, key=custom_sort)
        with open('./chr_mafft/3_seat/3.3/stat1/'+input_name,'w') as f:
            f.write('')
        for one_seatpath in sorted_seatpath_list:
            with open('./chr_mafft/3_seat/3.3/stat1/'+input_name,'a') as f:
                for base,value in dict_seatpath_base[one_seatpath].items():
                    f.write(one_seatpath+'\t'+base+'\t'+str(value)+'\n')        
        ############################################################################################
        ############################################################################################
        seat3_list = list(dict_seat3_base.keys())
        #Define a function to extract the number before the first '//', use sorted() with key parameter for sorting
        def extract_first_number(s):
            return int(s.split('//')[0])
        sorted_seat_list = sorted(seat3_list, key=extract_first_number)
        with open('./chr_mafft/3_seat/3.3/stat2/'+input_name,'w') as f:
            f.write('')
        for one_seat in sorted_seat_list:
            with open('./chr_mafft/3_seat/3.3/stat2/'+input_name,'a') as f:
                for base,value in dict_seat3_base[one_seat].items():
                    f.write(one_seat+'\t'+base+'\t'+str(value)+'\n')
        print('step3.3 '+input_file+'\tCompleted')
    #Assign tasks to processes in the pool
    with Pool(processes=thread) as pool:
        pool.map(run_seat, run_seat_list) 
if argv1=="stepall" or argv1=="step3"  or  argv1=="step3.3s1" or start_step<=36:
    print('step3.3s1')  #Estimated 20s
    if "dir" in args_dict:      
        dictionary=args_dict["dir"]
    else:                       
        dictionary='./chr_mafft/3_seat/3.3/stat1/';        
    if 'outfile'in args_dict:          
        outfile=args_dict["outfile"]
    else:
        outfile='./chr_mafft/3_seat/3.3_stat1_sum'
    files=os.listdir(dictionary)
    #
    dict_seatpath_base={}
    for onefile in files:
        with open(dictionary+'/'+onefile,'r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                seat=   eachline_arr[0]
                seq=    eachline_arr[1]
                num=    eachline_arr[2]
                if seat not in dict_seatpath_base:
                    dict_seatpath_base[seat]={}
                if seq  not in dict_seatpath_base[seat]:    
                    dict_seatpath_base[seat][seq]=0
                dict_seatpath_base[seat][seq]+=int(num)  
    ##################################################################
    ##################################################################
    seat1_list=list(dict_seatpath_base.keys())
    #Define a function to extract the number before the first '-', use sorted() with key parameter for sorting
    def custom_sort(item):
        key_part = item.split('\t')[0]
        if '-' in key_part:
            left, right = key_part.split('-')
            return (int(left), int(right))
        else:
            return (int(key_part),)
     
    #Sort the data
    sorted_seatpath_list = sorted(seat1_list, key=custom_sort)        
    with open(outfile,'w') as f:
        f.write('')
    for one_seatpath in sorted_seatpath_list:
        with open(outfile,'a') as f:
            for base,value in dict_seatpath_base[one_seatpath].items():
                f.write(one_seatpath+'\t'+base+'\t'+str(value)+'\n')  
    print('Recording the best monomer')      
    last_seat=''
    with open('./chr_mafft/3_seat/3.3_stat1_sum','r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if '-' in eachline_arr[0]:continue
            if last_seat=='': 
                with open('./chr_mafft/3_seat/3.3_stat3_seq','w') as f2:
                    f2.write(eachline_arr[1])
                last_seat= eachline_arr[0];continue   
            if eachline_arr[0]!=last_seat:
                with open('./chr_mafft/3_seat/3.3_stat3_seq','a') as f2:
                    f2.write(eachline_arr[1])
            last_seat= eachline_arr[0];        
if (argv1=="stepall" or argv1=="step3"  or  argv1=="step3.3s2" or start_step<=37) and stat_mark=='':
    print('step3.3s2')
    if "dir" in args_dict:      
        dictionary=args_dict["dir"]
    else:                       
        dictionary='./chr_mafft/3_seat/3.3/stat2/';        
    if 'outfile'in args_dict:          
        outfile=args_dict["outfile"]
    else:
        outfile='./chr_mafft/3_seat/3.3_stat2_sum'
    files=os.listdir(dictionary)
    #
    dict_seatpath_base={}
    for onefile in files:
        with open(dictionary+'/'+onefile,'r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                seat=   eachline_arr[0]
                seq=    eachline_arr[1]
                num=    eachline_arr[2]
                if seat not in dict_seatpath_base:
                    dict_seatpath_base[seat]={}
                if seq  not in dict_seatpath_base[seat]:    
                    dict_seatpath_base[seat][seq]=0
                dict_seatpath_base[seat][seq]+=int(num)  
    ##################################################################
    ##################################################################
    seat3_list=list(dict_seatpath_base.keys())
    #Define a function to extract numbers from '//' separated format, use sorted() with key parameter for sorting
    def custom_sort(item):
        key_part = item.split('\t')[0]
        col1,col2,col3 = key_part.split('//')
        return (int(col1), int(col2), int(col3))

    #Sort the data
    sorted_seatpath_list = sorted(seat3_list, key=custom_sort)        
    with open(outfile,'w') as f:
        f.write('')
    for one_seatpath in sorted_seatpath_list:
        with open(outfile,'a') as f:
            for base,value in dict_seatpath_base[one_seatpath].items():
                f.write(one_seatpath+'\t'+base+'\t'+str(value)+'\n')                       
#step3 Statistics and plotting
#step3 Statistics and plotting
if (argv1=="stepall" or argv1=="step3"  or  argv1=="step3.3p1" or start_step<=38) and stat_mark=='':    
    print('step3.3p1')
    if "ii" in args_dict:      
        input_file=args_dict["ii"]        
    else:
        input_file='./chr_mafft/3_seat/3.3_stat1_sum'
    if "oo" in args_dict:      
        output_file=args_dict["oo"]        
    else:
        output_file='./chr_mafft/3_seat/3.3_plot_1a'
    #output_file1=output_file+''
    #output_file2='./chr_mafft/3_seat/3.3_stat1_sum'+'_tmp2'
    with open(output_file,'w') as f2:
        f2.write('seat\tbase\tseat_x_center\tseat_y_center\tseat_x_centerreal\tseat_y_centerreal\tnum\tx\ty\tradius\n')
    seat_max=0
    num_max=0
    dict_seat_base_num={}
    dict_path_base_num={}
    dict_path_num={}
    with open(input_file,'r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            ##Seat position
            seatpath =eachline_arr[0]; 
            seq      =eachline_arr[1] 
            bad_mark=''
            for one in seq:
                if one not in ['A','C','G','T','-']:bad_mark='yes';break
            if bad_mark=='yes':continue               
            num      =int(eachline_arr[2]) 
            if '-' not in seatpath: 
                seatpath=int(seatpath)
                if num>num_max: num_max=num
                if seatpath>seat_max: seat_max=seatpath
                if seatpath not in dict_seat_base_num:
                    dict_seat_base_num[seatpath]={}
                if seq not in dict_seat_base_num[seatpath]:
                    dict_seat_base_num[seatpath][seq]=num
            else:
                if seatpath not in dict_path_base_num:
                    dict_path_base_num[seatpath]={}
                if  seq not in  dict_path_base_num[seatpath]:
                    dict_path_base_num[seatpath][seq]=num
                ##    
                if seatpath not in dict_path_num:
                    dict_path_num[seatpath]=0
                dict_path_num[seatpath]+=num  
    i=0
    each_angle=360/seat_max
    bigcirc_r_length=1000
    num_max_sqrt=math.sqrt(num_max)
    max_base_width = 2 * math.pi * bigcirc_r_length/seat_max/3
    while i<seat_max:
        seat=(i+1)
        dict_base_num=dict_seat_base_num[seat]
        if 'A' in dict_base_num:    A_num=dict_base_num['A']
        else:                       A_num=0
        if 'C' in dict_base_num:    C_num=dict_base_num['C']
        else:                       C_num=0
        if 'G' in dict_base_num:    G_num=dict_base_num['G']
        else:                       G_num=0
        if 'T' in dict_base_num:    T_num=dict_base_num['T']
        else:                       T_num=0        
        ACGT_num=A_num+C_num+G_num+T_num
        ####
        A_num_width=math.sqrt(A_num)/num_max_sqrt*max_base_width     ##A  C arranged in this order
        C_num_width=math.sqrt(C_num)/num_max_sqrt*max_base_width     ##G  T    
        G_num_width=math.sqrt(G_num)/num_max_sqrt*max_base_width
        T_num_width=math.sqrt(T_num)/num_max_sqrt*max_base_width
        onemax_base_width=max(A_num_width,C_num_width,G_num_width,T_num_width)
        seat_base_pic_width=max(A_num_width+C_num_width,G_num_width+T_num_width)
        seat_base_pic_height=max(A_num_width+G_num_width,C_num_width+T_num_width)
        ###Calculate offset, offset is the distance from new center to original center, i.e., half of (max base box side length minus seat_base_pic_width or seat_base_pic_height)
        offset_x_abs=onemax_base_width-seat_base_pic_width/2
        if      A_num_width==onemax_base_width:    offset_x=offset_x_abs;       offset_y=-offset_x_abs 
        elif    C_num_width==onemax_base_width:    offset_x=-offset_x_abs;      offset_y=-offset_x_abs 
        elif    G_num_width==onemax_base_width:    offset_x=offset_x_abs;       offset_y=offset_x_abs 
        elif    T_num_width==onemax_base_width:    offset_x=-offset_x_abs;      offset_y=offset_x_abs         
        ####Center point coordinates, need to calculate radian first
        angle_radians = math.radians(each_angle*i)    
        cosine_value = math.cos(angle_radians);     one_x = bigcirc_r_length*cosine_value
        sine_value = math.sin(angle_radians);       one_y = bigcirc_r_length*sine_value
        #print(cosine_value)
        ####New center point coordinates
        one_x_real=one_x+offset_x
        one_y_real=one_y+offset_y
        ####Center coordinates for each base
        A_box_x=str(one_x_real-A_num_width/2)
        A_box_y=str(one_y_real+A_num_width/2)
        C_box_x=str(one_x_real+C_num_width/2) 
        C_box_y=str(one_y_real+C_num_width/2)
        G_box_x=str(one_x_real-G_num_width/2) 
        G_box_y=str(one_y_real-G_num_width/2)
        T_box_x=str(one_x_real+T_num_width/2) 
        T_box_y=str(one_y_real-T_num_width/2)
        with open(output_file,'a') as f2:
            f2.write(str(seat)+'\t'+'A'+'\t'+str(one_x)+'\t'+str(one_y)+'\t'+str(one_x_real)+'\t'+str(one_y_real)+'\t'+str(A_num)+'\t'+A_box_x+'\t'+A_box_y+'\t'+str(A_num_width/2)+'\n')
            f2.write(str(seat)+'\t'+'C'+'\t'+str(one_x)+'\t'+str(one_y)+'\t'+str(one_x_real)+'\t'+str(one_y_real)+'\t'+str(C_num)+'\t'+C_box_x+'\t'+C_box_y+'\t'+str(C_num_width/2)+'\n')
            f2.write(str(seat)+'\t'+'G'+'\t'+str(one_x)+'\t'+str(one_y)+'\t'+str(one_x_real)+'\t'+str(one_y_real)+'\t'+str(G_num)+'\t'+G_box_x+'\t'+G_box_y+'\t'+str(G_num_width/2)+'\n')
            f2.write(str(seat)+'\t'+'T'+'\t'+str(one_x)+'\t'+str(one_y)+'\t'+str(one_x_real)+'\t'+str(one_y_real)+'\t'+str(T_num)+'\t'+T_box_x+'\t'+T_box_y+'\t'+str(T_num_width/2)+'\n')
        ####
        i+=1
    ######################################
    #####################################
    path_list=list(dict_path_num.keys())
    #Define a function to extract the number before the first '-', use sorted() with key parameter for sorting
    def custom_sort(item):
        key_part = item.split('\t')[0]
        left, right = key_part.split('-')
        return (int(left), int(right))
    sorted_path_list = sorted(path_list, key=custom_sort)     
    ###
    '''
    with open(output_file2,'w') as f2:
        f2.write('path\t'+'num'+'\n')
    ###
    for one_path in sorted_path_list:
        path_num=dict_path_num[one_path]
        #dict_path_seq=dict_path_base_num[one_path]
        with open(output_file2,'a') as f2:
            f2.write(one_path+'\t'+str(path_num)+'\n')
    '''        
    #print(dict_seat3_base)
    print('Outputting stat results')  

    ############################################################################################
    ############################################################################################
if (argv1=="stepall" or argv1=="step3"  or  argv1=="step3.3p2" or start_step<=39) and stat_mark=='':        
    print('step3.3p2')
    if "ii" in args_dict:      
        input_file=args_dict["ii"]        
    else:
        input_file='./chr_mafft/3_seat/3.3_stat2_sum'
    #
    if "ii2" in args_dict:      
        input_file2=args_dict["ii2"]        
    else:
        input_file2='./chr_mafft/3_seat/3.3_plot_1a'   
    #   
    if "oo" in args_dict:      
        output_file=args_dict["oo"]        
    else:
        output_file='./chr_mafft/3_seat//3.3_plot_1b'
    dict_newpath_num={}
    dict_newpath_seq_num={}
    newpath_seq_num_max=0
    with open(input_file,'r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            path     =eachline_arr[0]
            seq      =eachline_arr[1];  
            bad_mark=''
            for one in seq:
                if one not in ['A','C','G','T','/']:bad_mark='yes';break
            if bad_mark=='yes':continue   
            num      =int(eachline_arr[2]) 
            newpath_seq_num_max+=num
            path_arr=   path.split('//')
            seq_arr=    seq.split('/') 
            new_path=path_arr[1]+seq_arr[2]+'-'+path_arr[2]+seq_arr[4]  #e.g., 1A-3C
            if new_path not  in dict_newpath_num:
                dict_newpath_num[new_path]=0
            dict_newpath_num[new_path]+=num 
            ##
            if new_path not  in dict_newpath_seq_num:
                dict_newpath_seq_num[new_path]={}
            if seq_arr[3]=='':mid_seq='-'
            else: mid_seq=seq_arr[3]
            if mid_seq not in dict_newpath_seq_num[new_path]:
                dict_newpath_seq_num[new_path][mid_seq]=0
            dict_newpath_seq_num[new_path][mid_seq]+=num
      
    ############################################################################################
    ############################################################################################
    newpath_list=list(dict_newpath_num.keys())
    #print(newpath_list)
    #Define a function to extract the number before the first '-', use sorted() with key parameter for sorting
    def custom_sort(item):
        key_part = item.split('\t')[0]
        left, right = key_part.split('-')
        return (int(left[:-1]), int(right[:-1]))
    sorted_newpath_list = sorted(newpath_list, key=custom_sort)     
    ###
    with open(output_file,'w') as f2:
        f2.write('name\tflypast\tpath1\tpath1_x\tpath1_y\tpath2\tpath2_x\tpath2_y\tnum\tline_width\topacity\tinfo\n')
    ###
    dict_reviseseat_pos={}
    seat_max=0
    with open(input_file2,'r') as f3:
        for line in f3.readlines():
            eachline_arr=line.strip().split('\t')
            if eachline_arr[0]=='seat':continue
            if int(eachline_arr[0])>seat_max:seat_max=int(eachline_arr[0])
            reviseseat=eachline_arr[0]+eachline_arr[1]		
            x=eachline_arr[7]
            y=eachline_arr[8]
            dict_reviseseat_pos[reviseseat]=[x,y]
    #1px corresponds to width
    px1_num=newpath_seq_num_max/seat_max/10
    #
    for one_path in sorted_newpath_list:
        one_path_arr=one_path.split('-')
        path1=one_path_arr[0];  path1_x=dict_reviseseat_pos[path1][0];  path1_y=dict_reviseseat_pos[path1][1];
        path2=one_path_arr[1];  path2_x=dict_reviseseat_pos[path2][0];  path2_y=dict_reviseseat_pos[path2][1];
        seat1=int(path1[:-1])
        seat2=int(path2[:-1])
        if seat2>=seat1: flypast=seat2-seat1-1
        else:flypast=seat_max-seat1+seat2-1
        #
        newpath_num=dict_newpath_num[one_path]
        line_width=round(newpath_num/px1_num,1)
        if line_width<1:line_width=1
        #    
        opacity=round(newpath_num/px1_num,3)    
        if opacity>0.8:opacity=0.8
        #    
        dict_path_seq=dict_newpath_seq_num[one_path]
        #print(dict_path_seq)
        sorteddict_path_seq = dict(sorted(dict_path_seq.items(), key=lambda item: item[1], reverse=True)) 
        oneinfo=''
        for key,value in sorteddict_path_seq.items():
            oneinfo+=(key+'|'+str(value)+';')
        with open(output_file,'a') as f2:
            f2.write(one_path+'\t'+str(flypast)+'\t'+path1+'\t'+path1_x+'\t'+path1_y+'\t'+path2+'\t'+path2_x+'\t'+path2_y+'\t'+str(newpath_num)+'\t'+str(line_width)+'\t'+str(opacity)+'\t'+oneinfo+'\n')      
if (argv1=="stepall" or argv1=="step3"  or  argv1=="step3.3p3" or start_step<=40) and stat_mark=='':        
    print('step3.3p3')
    if "ii" in args_dict:      
        input_file=args_dict["ii"]        
    else:
        input_file='./chr_mafft/3_seat/3.3_plot_1a'
    #
    if "ii2" in args_dict:      
        input_file2=args_dict["ii2"]        
    else:
        input_file2='./chr_mafft/3_seat/3.3_plot_1b'   
    #   
    if "oo" in args_dict:      
        output_file=args_dict["oo"]        
    else:
        output_file='./chr_mafft/3_seat/3.3_plot'
    output_file1=output_file+'_2a'      
    output_file2=output_file+'_2b'  
    output_file3=output_file+'_2c'      
    ############################################################################################
    ############################################################################################
    dict_seat_base_num={};seat_max=0;num_max=0
    with open(input_file,'r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if eachline_arr[0]=='seat':continue
            seat=eachline_arr[0]
            base=eachline_arr[1]
            num=int(eachline_arr[6])
            if num>num_max:num_max=num
            if int(seat)>seat_max:seat_max=int(seat)
            if seat not in dict_seat_base_num:
                dict_seat_base_num[seat]={}
            if base not in dict_seat_base_num[seat]:
                dict_seat_base_num[seat][base]=0
            dict_seat_base_num[seat][base]+=num
    with open(output_file1,'w') as f2:
        f2.write('seat\tbase\tseat_x_center\tseat_y_center\tseat_x_centerreal\tseat_y_centerreal\tnum\tx\ty\tradius\n')            
    with open(output_file2,'w') as f3:
        f3.write('seat\tbase\tnum\tx\ty\topacity\n')                
    i=0
    dict_seat_pos={}
    while i<seat_max+1:
        i+=1
        seat=i;seat_show=seat
        if seat>seat_max:seat=1
        dict_base_num=dict_seat_base_num[str(seat)]
        if 'A' in dict_base_num:    A_num=dict_base_num['A']
        else:                       A_num=0
        if 'G' in dict_base_num:    G_num=dict_base_num['G']
        else:                       G_num=0                
        if 'C' in dict_base_num:    C_num=dict_base_num['C']
        else:                       C_num=0        
        if 'T' in dict_base_num:    T_num=dict_base_num['T']
        else:                       T_num=0  
        ACGT_num=A_num+C_num+G_num+T_num
        ####
        max_base_width=40
        num_max_sqrt=math.sqrt(num_max)
        A_num_width=math.sqrt(A_num)/num_max_sqrt*max_base_width     ##A  C arranged in this order
        C_num_width=math.sqrt(C_num)/num_max_sqrt*max_base_width     ##G  T    
        G_num_width=math.sqrt(G_num)/num_max_sqrt*max_base_width
        T_num_width=math.sqrt(T_num)/num_max_sqrt*max_base_width
        onemax_base_width=max(A_num_width,C_num_width,G_num_width,T_num_width)
        seat_base_pic_width=max(A_num_width+C_num_width,G_num_width+T_num_width)
        seat_base_pic_height=max(A_num_width+G_num_width,C_num_width+T_num_width)
        ###Calculate offset, offset is the distance from new center to original center, i.e., half of (max base box side length minus seat_base_pic_width or seat_base_pic_height)
        offset_x_abs=onemax_base_width-seat_base_pic_width/2
        if      A_num_width==onemax_base_width:    offset_x=offset_x_abs;       offset_y=-offset_x_abs 
        elif    C_num_width==onemax_base_width:    offset_x=-offset_x_abs;      offset_y=-offset_x_abs 
        elif    G_num_width==onemax_base_width:    offset_x=offset_x_abs;       offset_y=offset_x_abs 
        elif    T_num_width==onemax_base_width:    offset_x=-offset_x_abs;      offset_y=offset_x_abs 
        ###
        ##Center point
        one_x=70*i
        one_y=180
        ####New center point coordinates
        one_x_real=one_x+offset_x
        one_y_real=one_y+offset_y  
         ####Center coordinates for each base
        A_box_x=str(one_x_real-A_num_width/2)
        A_box_y=str(one_y_real+A_num_width/2)
        C_box_x=str(one_x_real+C_num_width/2) 
        C_box_y=str(one_y_real+C_num_width/2)
        G_box_x=str(one_x_real-G_num_width/2) 
        G_box_y=str(one_y_real-G_num_width/2)
        T_box_x=str(one_x_real+T_num_width/2) 
        T_box_y=str(one_y_real-T_num_width/2)
        with open(output_file1,'a') as f2:
            f2.write(str(seat_show)+'\t'+'A'+'\t'+str(one_x)+'\t'+str(one_y)+'\t'+str(one_x_real)+'\t'+str(one_y_real)+'\t'+str(A_num)+'\t'+A_box_x+'\t'+A_box_y+'\t'+str(A_num_width/2)+'\n')
            f2.write(str(seat_show)+'\t'+'C'+'\t'+str(one_x)+'\t'+str(one_y)+'\t'+str(one_x_real)+'\t'+str(one_y_real)+'\t'+str(C_num)+'\t'+C_box_x+'\t'+C_box_y+'\t'+str(C_num_width/2)+'\n')
            f2.write(str(seat_show)+'\t'+'G'+'\t'+str(one_x)+'\t'+str(one_y)+'\t'+str(one_x_real)+'\t'+str(one_y_real)+'\t'+str(G_num)+'\t'+G_box_x+'\t'+G_box_y+'\t'+str(G_num_width/2)+'\n')
            f2.write(str(seat_show)+'\t'+'T'+'\t'+str(one_x)+'\t'+str(one_y)+'\t'+str(one_x_real)+'\t'+str(one_y_real)+'\t'+str(T_num)+'\t'+T_box_x+'\t'+T_box_y+'\t'+str(T_num_width/2)+'\n')
        ####       
        A_pos=120
        C_pos=90
        G_pos=60
        T_pos=30
        dict_seat_pos[str(seat_show)+'A']=[one_x,A_pos]
        dict_seat_pos[str(seat_show)+'C']=[one_x,C_pos]
        dict_seat_pos[str(seat_show)+'G']=[one_x,G_pos]
        dict_seat_pos[str(seat_show)+'T']=[one_x,T_pos]
        with open(output_file2,'a') as f3:
            f3.write(str(seat_show)+'\t'+'A'+'\t'+str(A_num)+'\t'+str(one_x)+'\t'+str(A_pos)+'\t'+str(A_num/num_max)+'\n')
            f3.write(str(seat_show)+'\t'+'C'+'\t'+str(C_num)+'\t'+str(one_x)+'\t'+str(C_pos)+'\t'+str(C_num/num_max)+'\n')
            f3.write(str(seat_show)+'\t'+'G'+'\t'+str(G_num)+'\t'+str(one_x)+'\t'+str(G_pos)+'\t'+str(G_num/num_max)+'\n')
            f3.write(str(seat_show)+'\t'+'T'+'\t'+str(T_num)+'\t'+str(one_x)+'\t'+str(T_pos)+'\t'+str(T_num/num_max)+'\n')            
    ################################
    ################################
    with open(output_file3,'w') as f3:
        f3.write('name\tflypast\tpath1\tpath1_x\tpath1_y\tpath2\tpath2_x\tpath2_y\tnum\tline_width\topacity\tinfo\n')  
    with open(input_file2,'r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            name	    =eachline_arr[0]
            flypast	    =eachline_arr[1]
            path1	    =eachline_arr[2]
            path2	    =eachline_arr[5]
            num	        =eachline_arr[8]
            line_width  =eachline_arr[9]
            opacity     =eachline_arr[10]
            info        =eachline_arr[11]
            if eachline_arr[1]!='0':continue
            if path2[:-1] =='1': ##Replace the last 1 with 29
                path2=str(seat_max+1)+eachline_arr[5][-1];
                name=path1+'-'+path2
            path1_x=str(dict_seat_pos[path1][0])
            path1_y=str(dict_seat_pos[path1][1])
            path2_x=str(dict_seat_pos[path2][0])
            path2_y=str(dict_seat_pos[path2][1])
            with open(output_file3,'a') as f3:
                f3.write(name+'\t'+flypast+'\t'+path1+'\t'+path1_x+'\t'+path1_y+'\t'+path2+'\t'+path2_x+'\t'+path2_y+'\t'+num+'\t'+line_width+'\t'+opacity+'\t'+info+'\n')
if (argv1=="stepall" or argv1=="step3"  or  argv1=="step3.3p4" or start_step<=41) and stat_mark=='':
    print('step3.3p4')
    R_txt='''library(ggplot2)
library(dplyr)
setwd('./')
#Read data
input_file1 <- read.table('3.3_plot_1a', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
input_file2 <- read.table('3.3_plot_1b', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')

##
input_file1_simple<-unique.data.frame(input_file1[,c(1,3,4)])

sorted_input_file2<-input_file2[order(-input_file2[,9]),]
filtered_input_file2<- sorted_input_file2[1:500, ]

flypast_small<-filtered_input_file2%>% filter(flypast==0)
flypast_big<-filtered_input_file2%>% filter(flypast>0)

##Statistics
sum_num=sum(input_file1$num)
flypast_sum_num=sum(flypast_big$num)
#Initialize counter
count_dash <- 0
count_other <- 0
#Iterate through each element in the info column
for (info in flypast_big$info) {
  elements <- unlist(strsplit(info, ";"))
  for (element in elements) {
    parts <- unlist(strsplit(element, "\\\\|"))
    name <- parts[1]
    quantity <- as.numeric(parts[2])
    if (name == "-") {
      count_dash <- count_dash + quantity
    } else {
      count_other <- count_other + quantity
    }
  }
}
#cat("Sum of Name '-' :", count_dash, "\\n")   ##22173641      
#cat("Sum of Name not '-' :", count_other, "\\n")##2513

flypast_allbig<-input_file2%>% filter(flypast>0)
flypast_allbig_simple<-flypast_allbig[,c(1,2,9)]
flypast_allbig_simple_sorted<-flypast_allbig_simple[order(-flypast_allbig_simple[,3]),]
flypast_allbig_simple_sorted_15<-flypast_allbig_simple_sorted[1:15,]
flypast_allbig_simple_sorted_15a <- flypast_allbig_simple_sorted[16:nrow(flypast_allbig_simple_sorted),]
sum15a=sum(flypast_allbig_simple_sorted_15a$num)
flypast_allbig_simple_sorted_15other <- flypast_allbig_simple_sorted_15 %>%
  add_row(flypast = NA, num = sum15a, name = 'Other')
flypast_allbig_simple_sorted_15other <- flypast_allbig_simple_sorted_15other %>%  
  dplyr::mutate(name_order = row_number())#Add a row number column

#Calculate the four vertex coordinates of the rectangle (only need bottom-left and top-right corners)
input_file1$xmin <- input_file1$x - input_file1$radius
input_file1$ymin <- input_file1$y - input_file1$radius
input_file1$xmax <- input_file1$x + input_file1$radius
input_file1$ymax <- input_file1$y + input_file1$radius

#Define colors
color_A <- '#009e73'
color_C <- '#56b4e9'
color_G <- '#e69f00'
color_T <- '#cc79a7'


#Define an empty plot
p <- ggplot()
#Add line segment layer, but not mapping linewidth directly in aes()
p <- p + geom_segment(data = flypast_big, 
                      aes(x = path1_x, y = path1_y, xend = path2_x, yend = path2_y, alpha = opacity), 
                      linewidth=3,
                      color = "black",  #Set line color to blue
                      arrow = arrow(length = unit(0, "npc"), type = "closed", ends = "last"))
p <- p + geom_segment(data = flypast_small, 
                      aes(x = path1_x, y = path1_y, xend = path2_x, yend = path2_y, alpha = opacity), 
                      linewidth=0.1,
                      color = "black",  #Set line color to blue
                      arrow = arrow(length = unit(0, "npc"), type = "closed", ends = "last"))
#Add rectangle layer
p <- p + geom_rect(data = input_file1, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = base))
#Set fill colors
p <- p + scale_fill_manual(values = c(A = color_A, C = color_C, G = color_G, T = color_T))
#Add text layer
p <- p + geom_text(data = input_file1_simple, aes(x = seat_x_center*1.12, y = seat_y_center*1.12, label = seat),
                   color = "black", size = 3, vjust = 0.5, hjust = 0.5)
#Apply theme
p <- p + theme_classic()
p <- p+coord_equal()
#Add title and axis labels
p <- p + labs(title = "Plot of Rectangles with Segments", x =  paste("Number (All =", sum_num, ")", sep = ""), y = "Y Coordinate")
p <- p + theme(
  #axis.title = element_blank(),  #Hide axis titles
  axis.text = element_blank(),   #Hide axis tick labels
  axis.ticks = element_blank(),  #Hide axis ticks
  axis.line = element_blank()    #Hide axis lines
)


#Set PDF device size to 20cm x 20cm (convert to inches)
pdf("pic_step3.3p4a.pdf", width = 40 / 2.54, height = 40 / 2.54)

#Print graphic to PDF device
print(p)

#Close PDF device
dev.off()



desired_order=flypast_allbig_simple_sorted_15other$name
desired_order=rev(desired_order)
flypast_allbig_simple_sorted_15other$name <- factor(flypast_allbig_simple_sorted_15other$name, levels = desired_order)
#Custom color sequence
custom_colors <- c(
  "1" = "#cc9900",
  "5" = "#4ab1d3",
  "6" = "#ffad99",
  "NA" = "#8c8c8c"
)
flypast_allbig_simple_sorted_15other$flypast <- factor(flypast_allbig_simple_sorted_15other$flypast)
p <- ggplot(data = flypast_allbig_simple_sorted_15other, aes(x = name, y = num,fill=flypast)) +
  geom_col() +  #Use geom_col() to create bar chart
  labs(title = "Bar Chart of num by name", x = "Name", y = paste("Number (All =", flypast_sum_num, ")", sep = "")) +  #Add title and axis labels
  theme_classic() + #Apply theme
  coord_flip()+  #Flip the graph horizontally
  scale_fill_manual(values = custom_colors)#Apply custom color sequence



#Set PDF device size to 10cm x 10cm (convert to inches)
pdf("pic_step3.3p4b.pdf", width = 10 / 2.54, height = 10 / 2.54)

#Print graphic to PDF device
print(p)

#Close PDF device
dev.off()
'''
    with open('./chr_mafft/3_seat/pic_ggplot_step3.3p4.R','w',encoding='utf-8') as f:
        f.write(R_txt)
    new_directory = "./chr_mafft/3_seat/"
    os.chdir(new_directory)
    subprocess.run(['Rscript pic_ggplot_step3.3p4.R'], shell=True)    
    os.chdir('../../')
if (argv1=="stepall" or argv1=="step3"  or  argv1=="step3.3p5" or start_step<=42) and stat_mark=='':
    print('step3.3p5')
    R_txt='''library(ggplot2)
library(dplyr)

setwd('./')
#Read data
input_file1 <- read.table('3.3_plot_2a', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
input_file2 <- read.table('3.3_plot_2b', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
input_file3 <- read.table('3.3_plot_2c', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')

input_file2=input_file2%>%mutate(alpha=0.5*(1-opacity))

sorted_input_file3<-input_file3[order(-input_file3[,9]),]
filtered_input_file3<- sorted_input_file3[1:150, ]
input_file3_big<-filtered_input_file3%>% filter(opacity==0.8)
input_file3_small<-filtered_input_file3%>% filter(opacity<0.8)
#
seat_pos_simple=unique.data.frame(input_file2[,c(1,4)])
#
#For each seat, the base with the highest num and its corresponding x
seat_maxbase <- input_file2[,c(1,2,3,4)] %>%
  group_by(seat) %>%  #Group by seat
  arrange(desc(num)) %>%  #Sort by num descending within each group
  slice(1) %>%  #Select the first row in each group (i.e., row with highest num)
  ungroup()  #Ungroup
#print(seat_maxbase)


#Calculate the four vertex coordinates of the rectangle for input_file1
input_file1$xmin <- input_file1$x - input_file1$radius
input_file1$ymin <- input_file1$y - input_file1$radius
input_file1$xmax <- input_file1$x + input_file1$radius
input_file1$ymax <- input_file1$y + input_file1$radius

#Calculate the four vertex coordinates of the rectangle for input_file2
input_file2$xmin <- input_file2$x - 15
input_file2$ymin <- input_file2$y - 15
input_file2$xmax <- input_file2$x + 15
input_file2$ymax <- input_file2$y + 15

#Define colors
color_A <- '#009e73'
color_C <- '#56b4e9'
color_G <- '#e69f00'
color_T <- '#cc79a7'




#Define an empty plot
p <- ggplot()
#Add line segment layer, but not mapping linewidth directly in aes()

#
p <- p + geom_segment(data = input_file3_small, 
                      aes(x = path1_x, y = path1_y, xend = path2_x, yend = path2_y,alpha=opacity), 
                      linewidth =1,
                      color = "#336699",  #Set line color to blue
)

p <- p + geom_segment(data = input_file3_big, 
                      aes(x = path1_x, y = path1_y, xend = path2_x, yend = path2_y,linewidth=line_width), 
                      alpha = 0.8,
                      color = "#336699",  #Set line color to blue
                      #arrow = arrow(length = unit(0, "npc"), type = "closed", ends = "last")
)
print(p)
#Add rectangle layerp <- p + geom_rect(data = input_file1, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = base))

#p <- p + geom_rect(data = input_file2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white")  #Place an opaque layer at the bottom
p <- p + geom_rect(data = input_file2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = base))
p <- p + geom_rect(data = input_file2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,alpha=alpha), fill = "white") 
##Originally I wanted to add color blocks with different transparencies on a white background, but it affected the line display. So I added white boards with different transparencies on top of the color blocks, which works well. Adjusting the top 150 line segments gives the best display effect

#print(p)
#Set fill colors
p <- p + scale_fill_manual(values = c(A = color_A, C = color_C, G = color_G, T = color_T))
#Add text layer
p <- p + geom_text(data = seat_pos_simple, aes(x = x, y = -50, label = seat),
                   color = "black", size = 3, vjust = 0.5, hjust = 0.5)
p <- p + geom_text(data = seat_maxbase, aes(x = x, y = -10, label = base),
                   color = "black", size = 3, vjust = 0.5, hjust = 0.5)

#Apply theme
p <- p + theme_classic()
p <- p+coord_equal()
#Add title and axis labels
p <- p + labs(title = "Plot of Rectangles with Segments", x = "X Coordinate", y = "Y Coordinate")
p <- p + theme(
  axis.title = element_blank(),  #Hide axis titles
  axis.text = element_blank(),   #Hide axis tick labels
  axis.ticks = element_blank(),  #Hide axis ticks
  axis.line = element_blank()    #Hide axis lines
)

#Set PDF device size to 10cm x 10cm (convert to inches)
pdf("pic_step3.3p5.pdf", width = 25 / 2.54, height = 15 / 2.54)

#Print graphic to PDF device
print(p)

#Close PDF device
dev.off()



'''
    with open('./chr_mafft/3_seat/pic_ggplot_step3.3p5.R','w',encoding='utf-8') as f:
        f.write(R_txt)
    new_directory = "./chr_mafft/3_seat/"
    os.chdir(new_directory)
    subprocess.run(['Rscript pic_ggplot_step3.3p5.R'], shell=True)  
    os.chdir('../../')
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    
#step4 - Assemble SEAT positions into monomers, then UMAP dimensionality reduction and visualization, Estimated 300s
if argv1=="stepall" or argv1=="step4"  or  argv1=="step4.0" or start_step<=43:
    print('step4.0 ———— Analyzing monomers') ###Estimated 94s
    if  os.path.exists('./chr_mafft/4_monomer')==False:
        subprocess.run(["mkdir ./chr_mafft/4_monomer"], shell=True)
    if  os.path.exists('./chr_mafft/4_monomer/4.0'):
        subprocess.run(["rm -r ./chr_mafft/4_monomer/4.0"], shell=True)   
    subprocess.run(["mkdir ./chr_mafft/4_monomer/4.0"], shell=True)   
    #if  os.path.exists('./chr_mafft/4_monomer/monomer_0')==True:
    #    subprocess.run(["rm ./chr_mafft/4_monomer/monomer_0"], shell=True)   
    run_monomer_list=[]
    i=0
    print('Reading ./chr_mafft/2_good_regions')
    with open('./chr_mafft/2_good_regions','r') as f:
        for line in f.readlines():
            i+=1
            if i==1:continue
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=13:continue 
            run_monomer_list.append(eachline_arr) 
    print('Main task started\n')     
    def run_monomer(run_monomer_list_one):
        #print(run_monomer_list_one[1])
        region_name=        run_monomer_list_one[0]
        strand=             run_monomer_list_one[8];#print (strand)
        if strand=='plus':  sign='+'
        else:               sign='-'
        bigblock_mafftstart=  run_monomer_list_one[2]
        bigblock_mafftend=  run_monomer_list_one[3]
        bigblock_chrstart=  run_monomer_list_one[4];        
        bigblock_chrend=    run_monomer_list_one[5];        
        input_name=         region_name+sign+':'+bigblock_chrstart+'-'+bigblock_chrend 
        #if input_name!='region_4+:12565351-13295334':return ''
        input_file=         './chr_mafft/3_seat/3.3/'+input_name
        output_file=        './chr_mafft/4_monomer/4.0/'+input_name
        output_file2=        './chr_mafft/4_monomer/monomer_0'
        #Read file
        with open(input_file,'r') as f:
            lines = f.readlines()    #Read all lines into a list
        if sign=='-': lines.reverse()
        
        #Generate header
        i=0;seat_list=[]
        while i<monomer_len:
            i+=1
            seat_list.append(str(i))
            seat_list.append(str(i+0.5))
        ##Generate dictionary
        #     0        2        4          6    8    10     12     14  16  18 20 
        #14213171//14213172//14213173	A//G//T	26//27//28	-	A//C//T	1//2//3
        #14213172//14213173//14213174	G//T//A	27//28//1	-	T//A//C	28//1//2
        dict_circleserial_seat_seq={}
        dict_onecirc_info={}
        last_seat=''
        circ_serial=0
        chrserial_list=[]
        with open(output_file,'w') as f:
            head='region_name\tcirc_serial\tchrserial_1\tchrserial_2\tstrand\t'+'\t'.join(seat_list)+'\tcirc_seq\n'
            f.write(head)
        for line in lines:
            #print (line)
            #print(dict_onecirc_info)
            if 'N' in line:continue
            eachline_arr=line.strip().split('/')
            chrserial   =eachline_arr[2]
            seat        =eachline_arr[18]
            seat_right  =str(int(seat)+0.5)
            #
            seq         =eachline_arr[14]
            seq_right   =eachline_arr[15]
            #
            if last_seat!='' and int(seat)<=int(last_seat):
                circ_serial+=1
                oneline=region_name+'\t'+'circ_'+str(circ_serial)+'\t'+chrserial_list[0]+'\t'+chrserial_list[-1]+'\t'+sign+'\t'
                one_circ_seq=''
                for one_seat in seat_list:
                    if one_seat in dict_onecirc_info:   oneline+=dict_onecirc_info[one_seat]+'\t';              one_circ_seq+=dict_onecirc_info[one_seat]+'|'
                    else:                               oneline+='\t';                                          one_circ_seq+='|'
                oneline+=one_circ_seq
                with open(output_file,'a') as f:
                    f.write(oneline+'\n')
                #with open(output_file2,'a') as f2:
                #    f2.write(oneline+'\n')                
                dict_onecirc_info={}   
                chrserial_list=[]
            ###########    
            dict_onecirc_info[seat]=seq
            dict_onecirc_info[seat_right]=seq_right
            chrserial_list.append(chrserial)
            #
            last_seat   =seat
        if chrserial_list!=[]:
            circ_serial+=1
            oneline=region_name+'\t'+'circ_'+str(circ_serial)+'\t'+chrserial_list[0]+'\t'+chrserial_list[-1]+'\t'+sign+'\t'
            one_circ_seq=''
            for one_seat in seat_list:
                if one_seat in dict_onecirc_info:   oneline+=dict_onecirc_info[one_seat]+'\t';              one_circ_seq+=dict_onecirc_info[one_seat]+'|'
                else:                               oneline+='\t';                                          one_circ_seq+='|'
            oneline+=one_circ_seq
            with open(output_file,'a') as f:
                f.write(oneline+'\n')            
        
    #Assign tasks to processes in the pool
    with Pool(processes=thread) as pool:
        pool.map(run_monomer, run_monomer_list)     
if argv1=="stepall" or argv1=="step4"  or  argv1=="step4.0_plus" or start_step<=43:
    print('step4.0_plus ———— Preventing generation of overly short monomers') ###Estimated 94s
    if  os.path.exists('./chr_mafft/4_monomer/4.0_plus'):
        subprocess.run(["rm -r ./chr_mafft/4_monomer/4.0_plus"], shell=True)   
    subprocess.run(["mkdir ./chr_mafft/4_monomer/4.0_plus"], shell=True)   
    if  os.path.exists('./chr_mafft/4_monomer/monomer_0')==True:
        subprocess.run(["rm ./chr_mafft/4_monomer/monomer_0"], shell=True)  
    #
    run_monomer_list=[]
    i=0
    print('Reading ./chr_mafft/2_good_regions')
    with open('./chr_mafft/2_good_regions','r') as f:
        for line in f.readlines():
            i+=1
            if i==1:continue
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=13:continue 
            run_monomer_list.append(eachline_arr) 
    print('Main task started\n')     
    def run_monomer(run_monomer_list_one):
        #print(run_monomer_list_one[1])
        region_name=        run_monomer_list_one[0]
        strand=             run_monomer_list_one[8];#print (strand)
        if strand=='plus':  sign='+'
        else:               sign='-'
        bigblock_mafftstart=  run_monomer_list_one[2]
        bigblock_mafftend=  run_monomer_list_one[3]
        bigblock_chrstart=  run_monomer_list_one[4];        
        bigblock_chrend=    run_monomer_list_one[5];        
        input_name=         region_name+sign+':'+bigblock_chrstart+'-'+bigblock_chrend 
        #if input_name!='region_1+:15437180-19063289':return ''
        input_file=         './chr_mafft/4_monomer/4.0/'+input_name
        output_file=        './chr_mafft/4_monomer/4.0_plus/'+input_name
        output_file2=        './chr_mafft/4_monomer/monomer_0'
        i=0
        lastline=''
        circ_serial=0
        with  open (output_file2,'a') as f3:
            with open (output_file,'w') as f2:
                with open (input_file,'r') as f:
                    for line in f.readlines():
                        i+=1
                        #print(i)
                        eachline=line.strip()
                        #if i==40:     sys.exit()
                        if i==1:
                            f2.write(eachline+'\n')
                            continue
                        if lastline=='':lastline=eachline;continue
                        #
                        lastline_arr=lastline.split('\t');      lastline_seq=lastline_arr[-1]
                        eachline_arr=eachline.split('\t');      eachline_seq=eachline_arr[-1]
                        #print(eachline_seq)
                        #
                        part_line_num1=0
                        for k in range(len(lastline_seq) - 1, -1, -1):
                            if lastline_seq[k]=='|':part_line_num1+=1 
                            else:break
                        if part_line_num1%2==1: 
                            circ_serial+=1
                            new_line=lastline_arr[0]+'\tcirc_'+str(circ_serial)+'\t'+'\t'.join(lastline_arr[2:])
                            f2.write(new_line+'\n');f3.write(new_line+'\n');lastline=eachline;continue
                        blank1_num=part_line_num1/2 -1
                        ##
                        part_line_num2=0
                        for x in eachline_seq:
                            if x=='|':part_line_num2+=1
                            else:break
                        if part_line_num2%2==1: 
                            circ_serial+=1
                            new_line=lastline_arr[0]+'\tcirc_'+str(circ_serial)+'\t'+'\t'.join(lastline_arr[2:])
                            f2.write(new_line+'\n');f3.write(new_line+'\n');lastline=eachline;continue                       
                        blank2_num=part_line_num2/2 
                        #print(blank1_num,blank2_num)
                        ##
                        if blank1_num+blank2_num>monomer_len-5:         #Set a merging threshold
                            circ_serial+=1
                            new_info=lastline_arr[0]+'\tcirc_'+str(circ_serial)+'\t'+lastline_arr[2]+'\t'+eachline_arr[3]+'\t'+eachline_arr[4]+'\t'
                            ##||||||||||||||||||||||||||||||||||||||||||||||||||||G||T||
                            ##newseq = left seq + internal extra columns + right seq
                            border_part_line=part_line_num1-1
                            new_seq1=lastline_seq[:-border_part_line] 
                            ###
                            vertilinenum=0;new_seq2='';new_seq3='';change_mark=''
                            for x in eachline_seq[::-1]:
                                if x=="|":vertilinenum+=1
                                if change_mark=='':     new_seq3+=x
                                else:                   new_seq2+=x
                                if vertilinenum==border_part_line: change_mark='yes'
                            new_seq=new_seq1+new_seq2[::-1].replace('|','')+new_seq3[::-1]
                            #print(new_seq3[::-1])
                            ################
                            all_new_line=new_info+'\t'.join(new_seq.split("|"))+new_seq
                            f2.write(all_new_line+'\n');
                            f3.write(all_new_line+'\n');
                            lastline=''
                        else:
                            circ_serial+=1
                            new_line=lastline_arr[0]+'\tcirc_'+str(circ_serial)+'\t'+'\t'.join(lastline_arr[2:])
                            f2.write(new_line+'\n');f3.write(new_line+'\n');lastline=eachline;continue         
                if lastline!='': 
                    lastline_arr=lastline.split('\t')
                    circ_serial+=1
                    new_line=lastline_arr[0]+'\tcirc_'+str(circ_serial)+'\t'+'\t'.join(lastline_arr[2:])
                    f2.write(new_line+'\n');f3.write(new_line+'\n')
            
                    
    #Assign tasks to processes in the pool
    with Pool(processes=thread) as pool:
        pool.map(run_monomer, run_monomer_list)     
        
    print('Statistical analysis of monomers')
    i=0;head_seat_part=''
    while i<monomer_len:
        i+=1
        head_seat_part+=str(i)+'\t'+str(i+0.5)+'\t'
    with open('./chr_mafft/4_monomer/monomer_0_stat','w') as f:
        f.write('monomer_name\t'+head_seat_part+'circ_seq'+'\tin_seat_num\tcirc_len\tnum\tvariant_distance\n')  
    circ_seq_list=[]
    with open('./chr_mafft/4_monomer/monomer_0','r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            circ_seq_list.append('\t'.join(eachline_arr[5:]))
    ###
    print('Calculating frequencies')
    frequency={}
    for item in circ_seq_list:      
        if item not in frequency:
            frequency[item]=1
        else:
            frequency[item]+=1
    #Use sorted() function to sort in descending order by value
    print('Sorting')
    sorted_dict_items = sorted(frequency.items(), key=lambda item: item[1], reverse=True)
    sorted_frequency = dict(sorted_dict_items)#Convert sorted list back to dictionary
    print('Loading genome consensus sequence')
    with open('./chr_mafft/3_seat/3.3_stat3_seq','r') as f:
        for line in f.readlines():
            genome_consensus_monomer=line.strip()
            break

    #Output the sorted dictionary
    print('Calculating variant_distance')
    dict_variantdistance_num={};dict_circlength_num={}
    i=0
    for info,value in sorted_frequency.items():
        i+=1
        info_arr=info.split('\t')
        seat_seq_arr=info_arr[:-1]
        seat_seq_sum_withlines=info_arr[-1]
        kk=0;variant_distance=0;in_seat_num=0
        while kk<monomer_len:
            if seat_seq_arr[2*kk]!=genome_consensus_monomer[kk]:  variant_distance+=1
            if seat_seq_arr[2*kk+1]!='': variant_distance+=len(seat_seq_arr[2*kk+1])
            if seat_seq_arr[2*kk]!='':in_seat_num+=1
            kk+=1
        circ_length=len(seat_seq_sum_withlines.replace('|',''))
        ##
        if str(variant_distance) not in dict_variantdistance_num:dict_variantdistance_num[str(variant_distance)]=0
        #if in_seat_num>monomer_len/2:  dict_variantdistance_num[str(variant_distance)]+=value  ##Some monomers with deletions exceeding half the monomer length are not counted
        dict_variantdistance_num[str(variant_distance)]+=value
        ##
        if str(circ_length) not in dict_circlength_num:dict_circlength_num[str(circ_length)]=0
        #if in_seat_num>monomer_len/2:  dict_circlength_num[str(circ_length)]+=value  ##Some monomers with deletions exceeding half the monomer length are not counted
        dict_circlength_num[str(circ_length)]+=value 
        ##
        with open('./chr_mafft/4_monomer/monomer_0_stat','a') as f:
            f.write('monomer_'+str(i)+'\t'+info+'\t'+str(in_seat_num)+'\t'+str(circ_length)+'\t'+str(value)+'\t'+str(variant_distance)+'\n')
    
    print('Statistical analysis of monomer length')
    with open('./chr_mafft/4_monomer/monomer_0_stat2_Length','w') as f:
        f.write('circ_length\tnum\n')
    sorted_dict_circlength_num = dict(sorted(dict_circlength_num.items(), key=lambda item:int(item[0])))
    for circ_length,num in sorted_dict_circlength_num.items():
        #if int(circ_length)>monomer_len*2:continue
        with open('./chr_mafft/4_monomer/monomer_0_stat2_Length','a') as f:
            f.write(circ_length+'\t'+str(num)+'\n')        
    
    
    print('Statistical analysis of monomer variation')
    with open('./chr_mafft/4_monomer/monomer_0_stat2_Variant','w') as f:
        f.write('variant_distance\tnum\n')
    sorted_dict_variantdistance_num = dict(sorted(dict_variantdistance_num.items(), key=lambda item:int(item[0])))
    for variant_distance,num in sorted_dict_variantdistance_num.items():
        if int(variant_distance)>monomer_len:continue
        with open('./chr_mafft/4_monomer/monomer_0_stat2_Variant','a') as f:
            f.write(variant_distance+'\t'+str(num)+'\n') 
if (argv1=="stepall" or argv1=="step4"  or  argv1=="step4.0p1" or start_step<=45)  and stat_mark=='':
    R_txt='''library(ggplot2)

setwd('./')
#Read data
input_file1 <- read.table('monomer_0_stat2_Length', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
input_file2 <- read.table('monomer_0_stat2_Variant', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')

#Length plot
sum_num=sum(input_file1$num)
p=ggplot(input_file1, aes(x = circ_length, y = num)) +
  geom_bar(stat = "identity")  #+scale_x_continuous(limits = c(1, 28))
p <- p + theme_classic()
p <- p + labs(
  title = "",
  x = paste("Circ length (All =", sum_num, ")", sep = ""),
  y = "Frequency"
)

pdf("pic_step4.0p1a_circlength.pdf", width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()


#Variant distance plot
sum_num=sum(input_file2$num)
p=ggplot(input_file2, aes(x = variant_distance, y = num)) +
  geom_bar(stat = "identity")  #+scale_x_continuous(limits = c(1, 28))
p <- p + theme_classic()
p <- p + labs(
  title = "", 
  x = paste("Variant Distance (All =", sum_num, ")", sep = ""),
  y = "Frequency"
)

pdf("pic_step4.0p1b_variantdistance.pdf", width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()
'''
    with open('./chr_mafft/4_monomer/pic_ggplot_step4.0p1.R','w',encoding='utf-8') as f:
        f.write(R_txt)
    new_directory = "./chr_mafft/4_monomer/"
    os.chdir(new_directory)
    subprocess.run(['Rscript pic_ggplot_step4.0p1.R'], shell=True)    
    os.chdir('../../')
if (argv1=="stepall" or argv1=="step4"  or  argv1=="step4.1" or start_step<=46)   and stat_mark=='':
    print('step4.1 ———— UMAP, reducing 28 bases to 2 dimensions') ###One block of chr1 has about 10,000+ monomers, estimated several thousand after deduplication. Pairwise comparison would be over 100 million pairs - too time-consuming. Distant pairs shouldn't have meaningful relationships anyway.
    print('Loading numpy, umap, and other libraries')
    import numpy as np
    import pandas as pd
    from ast import literal_eval
    import umap
    from umap import UMAP
    
    
    #Manually create a small DataFrame to simulate data
    '''
    #Simulate raw data
    data = {
        'name': ['monomer_1', 'monomer_2', 'monomer_3'],
        'vector_1': ['(1,0,0)', '(0,1,0)', '(0,0,1)'],
        'vector_2': ['(0,1,0)', '(0,0,1)', '(1,0,0)'],
        'vector_3': ['(0,0,1)', '(1,0,0)', '(0,1,0)']
    }
    ''' 
    if  os.path.exists('./chr_mafft/4_monomer/monomer_1_stat_b'):
        subprocess.run(["rm -r ./chr_mafft/4_monomer/monomer_1_stat_b"], shell=True)   
    print('One-Hot Encoding bases to vectors')
    with open('./chr_mafft/4_monomer/monomer_0_stat','r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            newline1=eachline_arr[0]
            if newline1=='monomer_name':continue
            eachline_arr2=eachline_arr[1:monomer_len*2]
            newline2_arr=[]
            for x in eachline_arr2:
                if      x=='A':x="(1,0,0,0,0,0)"
                elif    x=='C':x="(0,1,0,0,0,0)"
                elif    x=='G':x="(0,0,1,0,0,0)"
                elif    x=='T':x="(0,0,0,1,0,0)"
                elif    x=='':x= "(0,0,0,0,1,0)"
                else:         x= "(0,0,0,0,0,1)"   
                newline2_arr.append(x)
            newline2='\t'.join(newline2_arr)    
            with open('./chr_mafft/4_monomer/monomer_1_stat_b','a') as f:
                f.write(newline1+'\t'+newline2+'\n')
         

    print('Reading vector table, converting to numeric vectors')
    df = pd.read_csv('./chr_mafft/4_monomer/monomer_1_stat_b', sep='\t', header=None)

    #Define a function to convert string-represented vectors to numeric vectors
    def parse_vector(vector_str):
        return np.array(literal_eval(vector_str))
     
    #Apply this function to all columns except the first column
    for col in df.columns[1:]:
        df[col] = df[col].apply(parse_vector)

    #Extract feature values, i.e., all columns except the first (name column)
    parsed_df = df.iloc[:, 1:].values
    #print(parsed_df)
    print(parsed_df[0, :]) 
    #print(parsed_df[0, 0]) 
    
    
    print('Data preparation complete')
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))   
    features = np.array([np.concatenate(row) for row in parsed_df])
    #print(features)
     

    #Create UMAP model and adjust n_neighbors parameter
    print('Initializing umap_model')
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))) 
    umap_model = UMAP(n_components=2, n_neighbors=15)#, random_state=42 - setting random_state prevents multi-threading
    #n_neighbors=500, runtime 327s
    #n_neighbors=100, runtime 160s
    #n_neighbors=50, runtime 110s
    #n_neighbors=15, runtime 74s, good
    #n_neighbors=5, runtime 70s
    ###n_neighbors=2, runtime 1030s
    #################### Apply UMAP for dimensionality reduction
    print('Applying UMAP for dimensionality reduction')
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))) 
    try:
        reduced_features = umap_model.fit_transform(features)
        print("Reduced features:")
        print(reduced_features)
        #Merge dimensionality reduction results with name column
        result_df = pd.DataFrame(reduced_features, columns=['x', 'y'])
        result_df['monomer_name'] = df.iloc[:, 0]  #Assuming first column is the name column
        result_df = result_df[['monomer_name', 'x', 'y']]  #Adjust column order
        
        #Write to file
        result_df.to_csv('./chr_mafft/4_monomer/monomer_1_stat_b_umap',sep='\t',  index=False)
        print("Results written to ./chr_mafft/4_monomer/monomer_1_stat_b_umap")        
    except Exception as e:
        print("An error occurred:", e)

    #Read file A and file B
    df_a = pd.read_csv('./chr_mafft/4_monomer/monomer_0_stat',sep='\t',header=0)
    df_b = pd.read_csv('./chr_mafft/4_monomer/monomer_1_stat_b_umap',sep='\t',header=0)
     
    #Assume both file A and file B have a column named 'key' to be used as merge key
    merge_key = 'monomer_name'
     
    #Use merge function to combine DataFrames
    # how='inner' means inner join, only keep rows with matching keys in both DataFrames
    # To keep all rows and fill NaN values, use how='outer'
    merged_df = pd.merge(df_a, df_b, on=merge_key, how='inner')
     
    #Output the merged DataFrame to file c.csv
    merged_df.to_csv('./chr_mafft/4_monomer/monomer_1_stat_b_umap2',sep='\t',   index=False)
     
    #Print the merged DataFrame (optional)
    #Select the first column
    first_column_df = merged_df.iloc[:, [0]]
    #Select the last 7 columns
    last_7_columns_df = merged_df.iloc[:, -7:]
    #Merge the first column and last 7 columns
    result_df = pd.concat([first_column_df, last_7_columns_df], axis=1)
    ##
    result_df.to_csv('./chr_mafft/4_monomer/monomer_1_stat_b_umap2_simple',sep='\t',   index=False)
if (argv1=="stepall" or argv1=="step4"  or argv1=="step4.1"  or  argv1=="step4.1p1" or start_step<=47)  and stat_mark=='':
    print('Generating plots')
    print('1')
    i=0
    dict_circseq_info={};dict_circseq_region_num={}
    with open('./chr_mafft/4_monomer/monomer_1_stat_b_umap2_simple','r') as f2:
        for line in f2.readlines():
            i+=1
            #if i>2001:break
            eachline_arr=line.strip().split('\t')
            circseq=eachline_arr[1]
            if circseq=='circ_seq':continue
            dict_circseq_info[circseq]=line.strip()
            dict_circseq_region_num[circseq]={}
    print('2')
    with open('./chr_mafft/4_monomer/monomer_0','r') as f1:
        for line in f1.readlines():
            eachline_arr=line.strip().split('\t')
            circseq=eachline_arr[-1]
            region_name=eachline_arr[0]
            if circseq not in dict_circseq_region_num: continue
            if region_name not in dict_circseq_region_num[circseq]:
                dict_circseq_region_num[circseq][region_name]=0
            dict_circseq_region_num[circseq][region_name]+=1
    print('3')         
    with open('./chr_mafft/4_monomer/monomer_1_stat_b_umap2_simple_part','w') as f3:
        f3.write('monomer_name\tcirc_seq\tin_seat_num\tcirc_len\tnum\tvariant_distance\tx\ty\tother_region1\tother_region1_num\tother_region2\tother_region2_num\n')
    for circseq,region_num in dict_circseq_region_num.items():
        dict_region_num_sorted = dict(sorted(region_num.items(), key=lambda item: item[1], reverse=True))
        i=0; region_info1='';region_info2=''
        for region,num in dict_region_num_sorted.items():
            i+=1
            if i==1:region_info1=region+'\t'+str(num)
            if i==2:region_info2=region+'\t'+str(num)
            else:break
        if region_info2=='':region_info2='-'+'\t'+'0'
        one_circ_info=dict_circseq_info[circseq]
        with open('./chr_mafft/4_monomer/monomer_1_stat_b_umap2_simple_part','a') as f3:
            f3.write(one_circ_info+'\t'+region_info1+'\t'+region_info2+'\n')
            
            
    print('Generating plots')
    R_txt=f'''library(ggplot2)
library(dplyr)
#install.packages("stringr") 
#library("stringr")

#Set working directory
setwd('./')

#Read data
input_file1 <- read.table('monomer_1_stat_b_umap2_simple_part', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')

input_file2 <- input_file1 %>% filter(circ_len>20 & num >100)
input_file3 <- input_file2 %>%
  arrange(as.numeric(sub("region_", "", other_region1)))
region_order <- unique(input_file3$other_region1)

input_file3 <- input_file3 %>%
  #mutate(other_region1_chr = str_replace(other_region1, "^region_", "Chr"))
  mutate(other_region1_chr = gsub("^region_", "Chr", other_region1))


region_order <- unique(input_file3$other_region1)
chr_order <- unique(input_file3$other_region1_chr)
input_file3$other_region1_chr <- factor(input_file3$other_region1_chr, levels = chr_order)

all_num=sum(input_file3$num)
#Use mutate and case_when to add new column shape
input_file3 <- input_file3 %>%
  mutate(shape = case_when(
    circ_len == 28 ~ 14, # Circle
    circ_len == 23 ~ 15, # Triangle
    TRUE ~ 17          # Square (default for other cases)
  ))

#Create plot
p <- ggplot()
p <- ggplot(data = input_file3, aes(x = x, y = y, size = num,color =other_region1_chr)) +
  geom_point(      aes(shape = factor(shape)),        alpha =0.7, stroke = 0.1)
p <- p + scale_color_discrete(breaks = chr_order)
p <- p + theme_classic()
p <- p + coord_equal() 
#Add title and axis labels
p <- p + labs(title =paste("Umap (All =", all_num, ")", sep = ""), x = "X", y = "Y")
p <- p + theme(
  #axis.title = element_blank(),  #Hide axis titles
  #axis.text = element_blank(),   #Hide axis tick labels
  #axis.ticks = element_blank(),  #Hide axis ticks
  #axis.line = element_blank()    #Hide axis lines
)


pdf("pic_step4.1p1_monomer_umap.pdf", width = 20 / 2.54, height = 20 / 2.54)
print(p)
dev.off()

input_file2 <- input_file1 %>% filter(circ_len>10 & num >10)
input_file3 <- input_file2 %>%
  arrange(as.numeric(sub("region_", "", other_region1)))
region_order <- unique(input_file3$other_region1)

input_file3 <- input_file3 %>%
  #mutate(other_region1_chr = str_replace(other_region1, "^region_", "Chr"))
  mutate(other_region1_chr = gsub("^region_", "Chr", other_region1))
  
region_order <- unique(input_file3$other_region1)
chr_order <- unique(input_file3$other_region1_chr)
input_file3$other_region1_chr <- factor(input_file3$other_region1_chr, levels = chr_order)

all_num=sum(input_file3$num)
#Use mutate and case_when to add new column shape
input_file3 <- input_file3 %>%
  mutate(shape = case_when(
    circ_len == 28 ~ 14, # Circle
    circ_len == 23 ~ 15, # Triangle
    TRUE ~ 17          # Square (default for other cases)
  ))
  
#Create plot
p <- ggplot()
p <- ggplot(data = input_file3, aes(x = x, y = y, size = num,color =other_region1_chr)) +
  geom_point(      aes(shape = factor(shape)),        alpha =0.7, stroke = 0.1)
p <- p + scale_color_discrete(breaks = chr_order)
p <- p + facet_wrap(~ other_region1_chr, ncol = 5)
p <- p + theme_bw()
p <- p + coord_equal() 
#Add title and axis labels
p <- p + labs(title =paste("Umap (All =", all_num, ")", sep = ""), x = "X", y = "Y")
p <- p + theme(
  #axis.title = element_blank(),  #Hide axis titles
  #axis.text = element_blank(),   #Hide axis tick labels
  #axis.ticks = element_blank(),  #Hide axis ticks
  #axis.line = element_blank()    #Hide axis lines
)
pdf("pic_step4.1p1_monomer_umap2.pdf", width = 20 / 2.54, height = 20 / 2.54)
print(p)
dev.off()
'''
    with open('./chr_mafft/4_monomer/pic_ggplot_step4.1.R','w',encoding='utf-8') as f:
        f.write(R_txt)
    new_directory = "./chr_mafft/4_monomer/"
    os.chdir(new_directory)
    subprocess.run(['Rscript pic_ggplot_step4.1.R'], shell=True)    
    os.chdir('../../')   
if (argv1=="stepall" or argv1=="step4"  or  argv1=="step4.2" or start_step<=48)  and stat_mark=='':
    print('step4.2 ———— UMAP, reducing 28 bases to 1 dimension') 
    print('Loading numpy, umap, and other libraries')
    import numpy as np
    import pandas as pd
    from ast import literal_eval
    import umap
    from umap import UMAP



    print('Reading vector table, converting to numeric vectors')
    df = pd.read_csv('./chr_mafft/4_monomer/monomer_1_stat_b', sep='\t', header=None)

    #Define a function to convert string-represented vectors to numeric vectors
    def parse_vector(vector_str):
        return np.array(literal_eval(vector_str))
     
    #Apply this function to all columns except the first column
    for col in df.columns[1:]:
        df[col] = df[col].apply(parse_vector)

    #Extract feature values, i.e., all columns except the first (name column)
    parsed_df = df.iloc[:, 1:].values
    #print(parsed_df)
    print(parsed_df[0, :]) 
    #print(parsed_df[0, 0]) 
    
    
    print('Data preparation complete')
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))   
    features = np.array([np.concatenate(row) for row in parsed_df])
    #print(features)
     

    #Create UMAP model and adjust n_neighbors parameter
    print('Initializing umap_model')
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))) 
    umap_model = UMAP(n_components=1, n_neighbors=15)#, random_state=42 - setting random_state prevents multi-threading

    #################### Apply UMAP for dimensionality reduction
    print('Applying UMAP for dimensionality reduction')
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))) 
    try:
        reduced_features = umap_model.fit_transform(features)
        print("Reduced features:")
        print(reduced_features)
        #Merge dimensionality reduction results with name column
        result_df = pd.DataFrame(reduced_features, columns=['x1'])
        result_df['monomer_name'] = df.iloc[:, 0]  #Assuming first column is the name column
        result_df = result_df[['monomer_name', 'x1']]  #Adjust column order
        
        #Write to file
        result_df.to_csv('./chr_mafft/4_monomer/monomer_2_stat_b_umap',sep='\t',  index=False)
        print("Results written to ./chr_mafft/4_monomer/monomer_2_stat_b_umap")        
    except Exception as e:
        print("An error occurred:", e)

    #Read file A and file B
    df_a = pd.read_csv('./chr_mafft/4_monomer/monomer_0_stat',sep='\t',header=0)
    df_b = pd.read_csv('./chr_mafft/4_monomer/monomer_2_stat_b_umap',sep='\t',header=0)
     
    #Assume both file A and file B have a column named 'key' to be used as merge key
    merge_key = 'monomer_name'
     
    #Use merge function to combine DataFrames
    # how='inner' means inner join, only keep rows with matching keys in both DataFrames
    # To keep all rows and fill NaN values, use how='outer'
    merged_df = pd.merge(df_a, df_b, on=merge_key, how='inner')
     
    #Output the merged DataFrame to file c.csv
    merged_df.to_csv('./chr_mafft/4_monomer/monomer_2_stat_b_umap2',sep='\t',   index=False)
     
    #Print the merged DataFrame (optional)
    #Select the first column
    first_column_df = merged_df.iloc[:, [0]]
    #Select the last 6 columns
    last_6_columns_df = merged_df.iloc[:, -6:]
    #Merge the first column and last 6 columns
    result_df = pd.concat([first_column_df, last_6_columns_df], axis=1)
    ##
    result_df.to_csv('./chr_mafft/4_monomer/monomer_2_stat_b_umap2_simple',sep='\t',   index=False)
if (argv1=="stepall" or argv1=="step4"  or  argv1=="step4.3" or start_step<=49)  and stat_mark=='': 
    print('step4.3 ———— Filling table') 
    if  os.path.exists('./chr_mafft/4_monomer/4.3'):
        subprocess.run(["rm -r ./chr_mafft/4_monomer/4.3"], shell=True)   
    subprocess.run(["mkdir ./chr_mafft/4_monomer/4.3"], shell=True)   

    run_monomer_list=[]
    i=0
    print('Reading ./chr_mafft/2_good_regions')
    with open('./chr_mafft/2_good_regions','r') as f:
        for line in f.readlines():
            i+=1
            if i==1:continue
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=13:continue
            run_monomer_list.append(eachline_arr) 
    print('Main task started\n')     
    def run_monomer(run_monomer_list_one):
        #print(run_monomer_list_one[1])
        region_name=        run_monomer_list_one[0]
        strand=             run_monomer_list_one[8];#print (strand)
        if strand=='plus':  sign='+'
        else:               sign='-'
        bigblock_mafftstart=  run_monomer_list_one[2]
        bigblock_mafftend=  run_monomer_list_one[3]
        bigblock_chrstart=  run_monomer_list_one[4];        
        bigblock_chrend=    run_monomer_list_one[5];        
        input_name=         region_name+sign+':'+bigblock_chrstart+'-'+bigblock_chrend 
        #if input_name!='region_4+:12565351-13295334':return ''
        input_file=         './chr_mafft/4_monomer/4.0_plus/'+input_name
        output_file=        './chr_mafft/4_monomer/4.3/'+input_name
        #print('Loading dictionary 1')
        dict_monomerseq_info={}
        with open('./chr_mafft/4_monomer/monomer_1_stat_b_umap2_simple','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                monomerseq=eachline_arr[1]
                #monomer_56867	A||C||T||C||G||C||A||C||A||G||A||T||T||G||T||A||C||||G||T||T||T||T||C||T||G||G||T||	27	27	1	9	-3.471739	6.086634
                dict_monomerseq_info[monomerseq]='\t'.join(eachline_arr[2:])
        dict_monomerseq_info2={}
        with open('./chr_mafft/4_monomer/monomer_2_stat_b_umap2_simple','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                monomerseq=eachline_arr[1]
                #monomer_1		A||C||T||C||G||C||A||C||C||G||A||T||T||C||T||A||C||C||A||T||T||T||C||C||G||G||G||A||	28	28	18436	8	-18.013376
                dict_monomerseq_info2[monomerseq]=eachline_arr[-1]        
        #print(dict_monomerseq_info2)
        with open(input_file,'r') as f:
            for line in f.readlines():
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                circ_seq=eachline_arr[-1]
                if circ_seq=='circ_seq':
                    addition_cols_str='in_seat_num\tcirc_len\tnum\tvariant_distance\tx\ty\tx1'
                else:
                    addition_cols_str=dict_monomerseq_info[circ_seq]+'\t'+dict_monomerseq_info2[circ_seq]
                #print(eachline)
                #print(addition_cols_str)
                newline=eachline+'\t'+addition_cols_str  
                with open(output_file,'a') as f2:
                    f2.write(newline+'\n')
                
                       
    #Assign tasks to processes in the pool
    with Pool(processes=thread) as pool:
        pool.map(run_monomer, run_monomer_list)
    
#step5 - Preliminary HOR analysis
if (argv1=="stepall" or argv1=="step5"  or  argv1=="step5.0" or start_step<=50) and stat_mark=='':
    print('step5.0 ———— HOR analysis')
    if  os.path.exists('./chr_mafft/5_hor')==False:
        subprocess.run(["mkdir ./chr_mafft/5_hor"], shell=True)
    if  os.path.exists('./chr_mafft/5_hor/5.0'): 
        subprocess.run(["rm -r ./chr_mafft/5_hor/5.0"], shell=True)   
    subprocess.run(["mkdir ./chr_mafft/5_hor/5.0"], shell=True)   
    run_monomer_list=[]
    i=0
    print('Reading ./chr_mafft/2_good_regions')
    with open('./chr_mafft/2_good_regions','r') as f:
        for line in f.readlines():
            i+=1
            if i==1:continue
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=13:continue 
            run_monomer_list.append(eachline_arr) 
    print('Main task started\n')     
    def run_monomer(run_monomer_list_one):
        #print(run_monomer_list_one)
        region_name=        run_monomer_list_one[0]
        strand=             run_monomer_list_one[8];#print (strand)
        if strand=='plus':  sign='+'
        else:               sign='-'
        bigblock_mafftstart=  run_monomer_list_one[2]
        bigblock_mafftend=  run_monomer_list_one[3]
        bigblock_chrstart=  run_monomer_list_one[4];        
        bigblock_chrend=    run_monomer_list_one[5];        
        input_name=         region_name+sign+':'+bigblock_chrstart+'-'+bigblock_chrend 
        #if input_name!='region_4+:12565351-13295334':return ''
        input_file=         './chr_mafft/4_monomer/4.3/'+input_name
        output_file=        './chr_mafft/5_hor/5.0/'+input_name
        #region_name	circ_serial	chrserial_1	chrserial_2	strand	1	1.5	2	2.5	3	3.5	4	4.5	5	5.5	6	6.5	7	7.5	8	8.5	9	9.5	10	10.5	11	11.5	12	12.5	13	13.5	14	14.5	15	15.5	16	16.5	17	17.5	18	18.5	19	19.5	20	20.5	21	21.5	22	22.5	23	23.5	24	24.5	25	25.5	26	26.5	27	27.5	28	28.5	circ_seq	in_seat_num	circ_len	num	variant_distance	x	y	x1
        dict_serial_info={}
        list_x1=[];list_seatlen=[]
        i=-2
        with open(input_file,'r') as f:
            for line in f.readlines():
                i+=1
                if i==-1:headline=line.strip();continue
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                dict_serial_info[i]=eachline
                list_x1.append(eachline_arr[-1])
                list_seatlen.append(eachline_arr[-7])
        with open(output_file,'w') as f2:
            f2.write(headline+'\tHOR_mer_raw\tpenalty_average\n')        
        #
        list_x1_num=len(list_x1)
        shift_max=200
        index_start_1=0  
        while index_start_1<list_x1_num:
            ##########Determine index_start1+shift/2shift based on shift(kk)
            kk=1;good_mark=''
            while kk<shift_max:
                kk+=1
                index_start_2=index_start_1+kk;     
                index_start_3=index_start_2+kk;  
                index_start_4=index_start_3+kk;   #index_start_4 serves as the maximum boundary, HOR requires at least 3 monomers
                if index_start_4>list_x1_num:break
                ##Calculate number of comparisons needed based on shift
                pair_num=2*kk
                penalty_threshold=pair_num#*3    ##A fairly lenient condition, allowing up to 3 difference per comparison in 1D scale
                seatlen_penalty_threshold=pair_num*3
                ############
                jj=0;penalty=0;
                while jj<kk:
                    pos1=float(list_x1[index_start_1+jj])
                    pos2=float(list_x1[index_start_2+jj])
                    pos3=float(list_x1[index_start_3+jj])
                    abs1=abs(pos1-pos2)
                    abs2=abs(pos1-pos3)
                    if abs1>5 or abs2>5:break    #Not too strict for the first pass
                    penalty+=abs1+abs2
                    if penalty>penalty_threshold:break
                    ##
                    seatlen1=int(list_seatlen[index_start_1+jj])
                    seatlen2=int(list_seatlen[index_start_2+jj])
                    seatlen3=int(list_seatlen[index_start_3+jj])
                    seatlen_abs1=abs(seatlen1-seatlen2)
                    seatlen_abs2=abs(seatlen1-seatlen3)
                    if seatlen_abs1>1 or seatlen_abs2>1  :break
                    ##
                    jj+=1
                if jj==kk:  #All position sequences processed
                    good_mark='yes';break
                ###########
            if  good_mark=='yes':
                good_mark=''
                addtion_str=str(kk)+'\t'+str(round(penalty/pair_num,3))
            else:
                addtion_str='.\t.'
            newline=dict_serial_info[index_start_1]+'\t'+addtion_str
            with open(output_file,'a') as f2:
                f2.write(newline+'\n')
            index_start_1+=1 
    #Assign tasks to processes in the pool
    with Pool(processes=thread) as pool:
        pool.map(run_monomer, run_monomer_list)    
if (argv1=="stepall" or argv1=="step5"  or  argv1=="step5.1" or start_step<=52) and stat_mark=='':
    print('step5.1 ———— HOR revision')
    if  os.path.exists('./chr_mafft/5_hor/5.1'): 
        subprocess.run(["rm -r ./chr_mafft/5_hor/5.1"], shell=True)   
    subprocess.run(["mkdir ./chr_mafft/5_hor/5.1"], shell=True)   
    #
    with open ('./chr_mafft/5_hor/hor_1_stat_raw','w') as f:
        f.write('belong_pos\tblock_serial\tHOR_mer_raw\tcirc_start_serial\tcirc_end_serial\tchr_start\tchr_end\tcirc_num\n')
    run_monomer_list=[]
    i=0
    print('Reading ./chr_mafft/2_good_regions')
    with open('./chr_mafft/2_good_regions','r') as f:
        for line in f.readlines():
            i+=1
            if i==1:continue
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=13:continue 
            run_monomer_list.append(eachline_arr) 
    print('Main task started\n')     
    def run_monomer(run_monomer_list_one):
        #print(run_monomer_list_one)
        region_name=        run_monomer_list_one[0]
        strand=             run_monomer_list_one[8];#print (strand)
        if strand=='plus':  sign='+'
        else:               sign='-'
        bigblock_mafftstart=  run_monomer_list_one[2]
        bigblock_mafftend=  run_monomer_list_one[3]
        bigblock_chrstart=  run_monomer_list_one[4];        
        bigblock_chrend=    run_monomer_list_one[5];        
        input_name=         region_name+sign+':'+bigblock_chrstart+'-'+bigblock_chrend 
        #if input_name!='region_1+:15437180-19063289':return ''  #region_4+:12565350-13295346
        input_file=         './chr_mafft/5_hor/5.0/'+input_name
        output_file=        './chr_mafft/5_hor/5.1/'+input_name
        output_file2=        './chr_mafft/5_hor/5.1/'+input_name+'_s'
        #region_name	circ_serial	chrserial_1	chrserial_2	strand	1	1.5	2	2.5	3	3.5	4	4.5	5	5.5	6	6.5	7	7.5	8	8.5	9	9.5	10	10.5	11	11.5	12	12.5	13	13.5	14	14.5	15	15.5	16	16.5	17	17.5	18	18.5	19	19.5	20	20.5	21	21.5	22	22.5	23	23.5	24	24.5	25	25.5	26	26.5	27	27.5	28	28.5	circ_seq	in_seat_num	circ_len	num	variant_distance	x	y	x1    HOR_mer    penalty_average    HOR_mer_revise penalty_average_revise
        dict_serial_info={}
        dict_serial_HORmer={}
        dict_serial_circlen={}
        dict_circname_start={}
        dict_circname_end={}        
        i=-1
        #print('Reading data')
        with open(input_file,'r') as f:
            for line in f.readlines():
                i+=1
                if i==0:headline=line.strip();continue
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                dict_serial_info[i]=eachline
                dict_serial_HORmer[i]=eachline_arr[-2]
                dict_serial_circlen[i]=eachline_arr[-8]
                #
                if eachline_arr[4]!='-':
                    dict_circname_start[eachline_arr[1]]=eachline_arr[2]
                    dict_circname_end[eachline_arr[1]]=eachline_arr[3]
                else:    
                    dict_circname_start[eachline_arr[1]]=eachline_arr[3]
                    dict_circname_end[eachline_arr[1]]=eachline_arr[2]                
        with open(output_file,'w') as f2:
            f2.write(headline+'\tHOR_mer_revise\n')  
        with open(output_file2,'w') as f3:
            f3.write('belong_pos\tblock_serial\tHOR_mer_raw\tcirc_start_serial\tcirc_end_serial\tchr_start\tchr_end\tcirc_num\n')               
        #print('Analyzing data')
        dict_serial_HORmer_num=len(dict_serial_HORmer)
        i=0
        dict_serial_addition={}
        dict_serial_block={}
        dict_serial=0
        current_block_endserial='';current_block_mer='';current_block='';start_bad_mark='';end_result_mark=''
        while i<dict_serial_HORmer_num:
            i+=1 
            one_HOR_mer=dict_serial_HORmer[i]
            if one_HOR_mer=='1':one_HOR_mer='.'  ##Monomers with value 1 are not accurate enough, discard them
            #print(i)
            ##A possible new block start
            if one_HOR_mer!='.' and current_block=='': 
                #print(i)
                
                current_block_mer=int(one_HOR_mer); current_block_startserial=i;  current_block_endserial=i+current_block_mer*3-1;current_block=[current_block_startserial,current_block_endserial]
                #Check if there are smaller values within the range, discard current mer if found
                kk=0;start_bad_mark=''
                while kk<current_block_mer*3-1:
                    kk+=1                        
                    block_inner_serial=i+kk
                    one_HOR_mer_inner=dict_serial_HORmer[block_inner_serial]
                    if one_HOR_mer_inner!='.' and one_HOR_mer_inner!='.' :
                        if int(one_HOR_mer_inner)<current_block_mer : #and int(one_HOR_mer_inner)>1
                            current_block_mer='';current_block_startserial='';current_block_endserial='';current_block='';start_bad_mark='yes';break

                #If current block is valid, check for larger values after it    
                if start_bad_mark!='yes'  :
                    kk=0
                    while kk<current_block_mer*3-1:
                        kk+=1                        
                        block_inner_serial=i+kk
                        one_HOR_mer_inner=dict_serial_HORmer[block_inner_serial]
                        if one_HOR_mer_inner!='.':
                            if int(one_HOR_mer_inner)>current_block_mer:   dict_serial_HORmer[block_inner_serial]='.'

                #If current block is valid, store position
                if start_bad_mark!='yes'  :       
                    dict_serial+=1
                    dict_serial_block[dict_serial]={}
                    dict_serial_block[dict_serial][0]=current_block_mer
                    dict_serial_block[dict_serial][1]=current_block_startserial
                    dict_serial_block[dict_serial][2]=current_block_endserial
            ##Already have a valid block, check if extension is possible        
            elif current_block_endserial!=''  : 
                if current_block_endserial>=i and one_HOR_mer!='.' and current_block!='':
                    #print(1111111111111111)
                    current_block_mer=int(one_HOR_mer);  current_block_endserial=i+current_block_mer*3-1
                    #Check if there are smaller values within the range, discard current mer if found and need to output
                    kk=0;end_result_mark=''
                    while kk<current_block_mer*3-1:
                        kk+=1                        
                        block_inner_serial=i+kk
                        one_HOR_mer_inner=dict_serial_HORmer[block_inner_serial]
                        if one_HOR_mer_inner!='.' and one_HOR_mer_inner!='.' :
                            if int(one_HOR_mer_inner)<current_block_mer and int(one_HOR_mer_inner)>1:
                                end_result_mark='yes';break

                    #If current block is valid, check for larger values after it
                    if end_result_mark!='yes':
                        kk=0
                        while kk<current_block_mer*3-1:
                            kk+=1                            
                            block_inner_serial=i+kk
                            one_HOR_mer_inner=dict_serial_HORmer[block_inner_serial]
                            if one_HOR_mer_inner!='.':
                                if int(one_HOR_mer_inner)>current_block_mer:   dict_serial_HORmer[block_inner_serial]='.'

                    #If current block is valid, update position
                    if end_result_mark!='yes'  :    
                        #print("111111111111111111",dict_serial_block[dict_serial])
                        dict_serial_block[dict_serial][2]=current_block_endserial
            ##Store current mer
            if current_block_mer!='' :          dict_serial_addition[i]=current_block_mer  
            else:                               dict_serial_addition[i]='.'
            
            ##Check if reaching the final boundary    
            if  current_block_endserial!='':
                if current_block_endserial<i+1 :
                    current_block_mer='';current_block_startserial='';current_block_endserial='';current_block='';            
        #print(dict_serial_addition)
        ####Check again, treat very short mers (monomer length less than half) as not monomers
        import copy
        dict_serial_block_copied = copy.deepcopy(dict_serial_block)
        for key,value in dict_serial_block_copied.items():
            mer_num=int(value[0])
            i1=int(value[1])
            i2=int(value[2])
            bad_circ_num=0
            jj=i1
            while jj<=i1+mer_num-1:
                if int(dict_serial_circlen[jj])<monomer_len/2:bad_circ_num+=1;
                jj+=1
            if bad_circ_num>0:
                raw_mer=dict_serial_block[key][0]
                new_mer=raw_mer-bad_circ_num
                dict_serial_block[key][0]=new_mer
                
        ####Output 1 - Detailed
        i=0
        while i<dict_serial_HORmer_num:
            i+=1
            one_old_info=dict_serial_info[i]
            one_new_info=str(dict_serial_addition[i]);       
            with open(output_file,'a') as f2:
                f2.write(one_old_info+'\t'+one_new_info+'\n')
        ####Output 2 - Statistics
        for key,value in dict_serial_block.items():
            circ_start='circ_'+str(value[1]);       chr_start=dict_circname_start[circ_start]
            circ_end='circ_'+str(value[2]);       chr_end=dict_circname_end[circ_end]
            newline=input_name+'\tHOR_'+str(key)+'\t'+str(value[0])+'\t'+circ_start+'\t'+circ_end+'\t'+chr_start+'\t'+chr_end+'\t'+str(value[2]-value[1]+1)
            with open(output_file2,'a') as f3:
                f3.write(newline+'\n')        
            ##Aggregate and output
            with open ('./chr_mafft/5_hor/hor_1_stat_raw','a') as f4:
                f4.write(newline+'\n')
    #Assign tasks to processes in the pool
    with Pool(processes=thread) as pool:
        pool.map(run_monomer, run_monomer_list)    
                
        
    print('Organizing and summarizing all results')
    #Read and output the first line (header), while writing it to output file. Read from the second line, sort by column 1 and column 5, then append results to output file
    subprocess.run(["head -n 1 './chr_mafft/5_hor/hor_1_stat_raw' > './chr_mafft/5_hor/hor_1_stat1'"], shell=True)   
    subprocess.run(["tail -n +2 './chr_mafft/5_hor/hor_1_stat_raw' | sort -k1,1V -k2,2V >> './chr_mafft/5_hor/hor_1_stat1'"], shell=True)   
    subprocess.run(["rm './chr_mafft/5_hor/hor_1_stat_raw'"], shell=True)   
    dict_mer_num={}
    with open('./chr_mafft/5_hor/hor_1_stat1','r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if eachline_arr[0]=='belong_pos':continue
            mer=eachline_arr[2]
            num=int(eachline_arr[7])
            if mer not in dict_mer_num:
                dict_mer_num[mer]=0
            dict_mer_num[mer]+=num    
    #Sort dictionary in ascending order by value    
    dict_mer_num_sorted=dict(sorted(dict_mer_num.items(), key=lambda item:int(item[0])))    
    with open('./chr_mafft/5_hor/hor_1_stat2','w') as f:
        f.write('HOR_mer_revise\tnum\n')
    for key,value in dict_mer_num_sorted.items():
        with open('./chr_mafft/5_hor/hor_1_stat2','a') as f:
            f.write(key+'\t'+str(value)+'\n')        
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))              
if (argv1=="stepall" or argv1=="step5"  or  argv1=="step5.1p1" or start_step<=53) and stat_mark=='':
    R_txt='''library(ggplot2)

setwd('./')
#Read data
input_file1 <- read.table('hor_1_stat2', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')

#HOR
sum_num=sum(input_file1$num)
p=ggplot(input_file1, aes(x = HOR_mer_revise, y = num)) +
  geom_bar(stat = "identity") 
p <- p + theme_classic()
p <- p + labs(
  title = "",
  x = paste("HOR num (All =", sum_num, ")", sep = ""),
  y = "Frequency"
)

pdf("pic_step5.1p1_HOR_mer.pdf", width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()


'''
    with open('./chr_mafft/5_hor/pic_ggplot_step5.1p1.R','w',encoding='utf-8') as f:
        f.write(R_txt)
    new_directory = "./chr_mafft/5_hor/"
    os.chdir(new_directory)
    subprocess.run(['Rscript pic_ggplot_step5.1p1.R'], shell=True)    
    os.chdir('../../')


        
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))      
time_end=timeit.default_timer()
print('All the running time: %.0f Seconds'%(time_end-time_start))






