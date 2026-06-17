#!/usr/bin/python
# -*- coding: UTF-8 -*-
import sys
argvs=sys.argv
argv_len=len(sys.argv)
if argv_len < 2 or argvs[1] in ("h", "-h", "help", "-help") or not argvs[1].startswith("step") :
    print ("samples_satellite.py-----help:")
    print ("samples_satellite.py is used to summarize the chr_mafft.py results from all samples in the input_folder and perform population-level analysis.")
    print ("")
    print ("Usage：")
    print ("step3 -i  <folder>        #Organize and summarize results from chr2mafft script")
    print ("step4 -seed <integer>     #Add seed to stabilize UMAP results, select appropriate distribution for plotting, then proceed")
    
    print ('./samples_satellite/sample_info:')
    print ('\t\tSampleID1	SampleID2	 Species	Class1(Group)	Sum_name	Class2')
    print ('\t\tBaimunage.hap1	Baimunage	 vinifera	Eurasian	Eurasian_vinifera_Baimunage	Table\n\n\n')
    print ("step5                     #UMAP statistics: Grid calculation")
    
    print ("step9                     #Divide monomers using 23bp intervals")
    print ("step10                    #Monomer patterns and statistics")   

    print (" ")
    print ("step11                    #Determine HORs using monomers as the unit")  
    print ("step5s                    #Supplementary display of monomers and HORs in 2D UMAP plots")  
    print (" ")    
    print ("step12                    #Statistics for panel b/c/d of pic3: HOR layers, whether multiple monomer types are included, length distribution of each HOR monomer, total length of each HOR")  
    print ("step12b                   #Monomer statistics, CEN107 and CEN107-like, etc.") 
    print ("step12_var                Based on step11/5sum and step12/1_blocks_stat, calculate the variation of father-class HORs, pushing forward from the first circ until exceeding the HOR range.") 
    
    print (" ")   
    print ("step13                    #HOR and monomer plots for each region_block")  
    print ("")
    print ("step14                    #Chromosome integration, HOR plots, monomer plots. Visualization including genes, TEs, ChIP, etc., Chr18 visualization") 
    print ("step14                    #python 2-samples_satellite9.py step14.22_partall V079.hap1 Chr19 19420000-19655000     Detailed display of a single region segment") 
    
    print ("step15          Perform UMAP with MON135, MON107, MON79")  
    print (" ")  
    print ("step16          MON connection status")      
    print (" ")  
    print ("step17          MON variation status")
    print (" ")  
    print ("-thread \t\tNumber of threads (default: 50), used for multiprocessing in some steps")
    print ("-i      \t\tInput FASTA file (required for step0)")
    print ("-seed <integer>     #(step4/5): Add seed to stabilize UMAP results, select appropriate distribution for plotting, then proceed")
    print ("-simple          #Normal satellite monomer, no subunit. Step9 executes subunit-to-monomer conversion by default. Step10 skip.")
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
import re # Regular expression processing
from multiprocessing import Pool, cpu_count

### Parse/Separate parameters
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

##
time_start=timeit.default_timer()
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  

if "thread" in args_dict:  thread=int(args_dict["thread"])  
else:thread=70


mafft="mafft"   

if  os.path.exists('./samples_satellite')==False:
    subprocess.run(["mkdir samples_satellite"], shell=True)

subunit_mark='yes'  #
if "simple" in args_dict:  subunit_mark='no';




#step3, Collect all samples' monomers and preprocess for UMAP
if argv1=="stepall" or "step3" in argv1:
    if argv1=="stepall" or argv1=="step3"  or argv1=="step3.0":
        print("step3, total 515s, 170 hap genomes")
        print('step3 ———— Collect monomers from all samples and preprocess into UMAP format')
        if  os.path.exists('./samples_satellite/3_samples')==True:
            subprocess.run(["rm -r ./samples_satellite/3_samples"], shell=True)  
        subprocess.run(["mkdir ./samples_satellite/3_samples"], shell=True)    
        subprocess.run(["mkdir ./samples_satellite/3_samples/monomer"], shell=True)    
        if "i"  in args_dict:
            input_dir = args_dict["i"]
            if os.path.exists(input_dir)==False: print('Input folder does not exist');sys.exit()
        else:
            print("Input folder not specified. Format should be: -i folder_name");sys.exit()
        
        def get_first_level_directories(path):
            if not os.path.exists(path):
                raise FileNotFoundError(f"The path '{path}' does not exist.")
            if not os.path.isdir(path):
                raise NotADirectoryError(f"The path '{path}' is not a directory.")
            first_level_dirs = [entry for entry in os.listdir(path) if os.path.isdir(os.path.join(path, entry))]
            return first_level_dirs
     
        directory_path = input_dir
        try:
            dirs = get_first_level_directories(directory_path)
            print("First level directories:")
            print(dirs[:10],'...total '+str(len(dirs))+' samples')
            with open('./samples_satellite/3_samples/3_samples_list','w') as f:
                for x in dirs:
                    f.write(x+'\n')
        except (FileNotFoundError, NotADirectoryError) as e:
            print(e)
        dirs.sort()    
        ###
        print('Loading all monomers')
        dict_oldcircseq_newcircseq={}
        dict_monomer_num={};dict_monomer_info={}
        for x in dirs:
            if os.path.exists(input_dir+'/'+x+'/chr_mafft/4_monomer/monomer_0_stat')==False: print('error');sys.exit()
            k=0
            with open(input_dir+'/'+x+'/chr_mafft/4_monomer/monomer_0_stat') as f2:
                for line in f2.readlines():
                    eachline_arr=line.strip().split('\t')
                    k+=1
                    if k==1:continue
                    circ_seq=           eachline_arr[-5]
                    #####################################################################################################################################
                    # Correct sequence to 23bp positioned at 11-17,,,Only applicable to grapevine, based on subjective judgment. 
                    circ_seq_new=circ_seq.replace('|','')
                    if len(circ_seq_new)==23:
                        circ_seq_new='||'.join(circ_seq_new[:11])+'||||||||||||'+'||'.join(circ_seq_new[11:])+'||'
                        dict_oldcircseq_newcircseq[circ_seq]=circ_seq_new
                        circ_seq=circ_seq_new
                    #####################################################################################################################################
                    in_seat_num=        eachline_arr[-4]
                    circ_len=           eachline_arr[-3]
                    variant_distance=   eachline_arr[-1]
                    if circ_seq not in dict_monomer_info:
                        dict_monomer_info[circ_seq]=[in_seat_num,circ_len,variant_distance]
                    num=int(eachline_arr[-2])
                    if circ_seq not in dict_monomer_num:
                        dict_monomer_num[circ_seq]=0
                    dict_monomer_num[circ_seq]+=num
        sorted_dict_monomer_num = dict(sorted(dict_monomer_num.items(), key=lambda item: item[1], reverse=True))             
                    
                    
        all_monomer_list=list(sorted_dict_monomer_num.keys())
        with open('./samples_satellite/3_samples/3_samples_monomer','w') as f:        
            f.write('circ_seq\tin_seat_num\tcirc_len\tnum\tvariant_distance\n')
            for circ_seq,num in sorted_dict_monomer_num.items():
                infos=dict_monomer_info[circ_seq]
                in_seat_num=        infos[0]
                circ_len=           infos[1]
                variant_distance=   infos[2]
                f.write(circ_seq+'\t'+in_seat_num+'\t'+circ_len+'\t'+str(num)+'\t'+variant_distance+'\n')
                
        print('Generate UMAP format, convert sequence to numerical vectors')           
        with open('./samples_satellite/3_samples/3_samples_monomer_vector','w') as f:  
            for x in all_monomer_list:
                addition_cols_raw_arr=x.split('|')
                circ_seq_len_para_str=str(len(x)/1000)
                addition_cols_arr=[]
                for one in addition_cols_raw_arr:
                    if      one=='A':one=f"(1,0,0,0,0,0,{circ_seq_len_para_str})"
                    elif    one=='C':one=f"(0,1,0,0,0,0,{circ_seq_len_para_str})"
                    elif    one=='G':one=f"(0,0,1,0,0,0,{circ_seq_len_para_str})"
                    elif    one=='T':one=f"(0,0,0,1,0,0,{circ_seq_len_para_str})"
                    elif    one=='':one= f"(0,0,0,0,1,0,{circ_seq_len_para_str})"
                    else:           one= f"(0,0,0,0,0,1,{circ_seq_len_para_str})"   
                    addition_cols_arr.append(one)
                addition_cols='\t'.join(addition_cols_arr) 
                f.write(x+'\t'+addition_cols+'\n')      
                
        print('Loading all regions')     
        sample_region_list=[]
        with open('./samples_satellite/2_good_regions','w') as f:
            headline='sample\tregion_name\tregion_pos\tbigblock_chrstart\tbigblock_chrend\tchr_region_length\tstrand\tmatch_percent\n'
            f.write(headline)
            for one_sample in dirs:
                i=0
                with open(input_dir+'/'+one_sample+'/chr_mafft/2_good_regions','r') as f2:
                    for line in f2.readlines():
                        i+=1
                        if i==1:continue
                        eachline_arr=line.strip().split('\t')
                        if len(eachline_arr)<13:
                            f.write('\n')
                            continue
                        if eachline_arr[8]=='minus':     strand='-'    
                        else:                           strand='+'     
                        one_region_name=eachline_arr[0]+strand+':'+eachline_arr[4]+'-'+eachline_arr[5]
                        #Check that region1 corresponds to chr1
                        region_name2chr=f"Chr{int(eachline_arr[0][7:])}"
                        one_id=eachline_arr[1].split(':')[0]
                        one_id_lower=one_id.lower()
                        if one_id.isdigit():
                            one_id="Chr"+one_id
                        elif one_id_lower.endswith('chr01') or one_id_lower.endswith('chr1'):   one_id="Chr1";      
                        elif one_id_lower.endswith('chr02') or one_id_lower.endswith('chr2'):   one_id="Chr2";      
                        elif one_id_lower.endswith('chr03') or one_id_lower.endswith('chr3'):   one_id="Chr3";      
                        elif one_id_lower.endswith('chr04') or one_id_lower.endswith('chr4'):   one_id="Chr4";      
                        elif one_id_lower.endswith('chr05') or one_id_lower.endswith('chr5'):   one_id="Chr5";      
                        elif one_id_lower.endswith('chr06') or one_id_lower.endswith('chr6'):   one_id="Chr6";      
                        elif one_id_lower.endswith('chr07') or one_id_lower.endswith('chr7'):   one_id="Chr7";      
                        elif one_id_lower.endswith('chr08') or one_id_lower.endswith('chr8'):   one_id="Chr8";      
                        elif one_id_lower.endswith('chr09') or one_id_lower.endswith('chr9'):   one_id="Chr9";      
                        elif one_id_lower.endswith('chr10'):   one_id="Chr10";      
                        elif one_id_lower.endswith('chr11'):   one_id="Chr11";      
                        elif one_id_lower.endswith('chr12'):   one_id="Chr12";      
                        elif one_id_lower.endswith('chr13'):   one_id="Chr13";      
                        elif one_id_lower.endswith('chr14'):   one_id="Chr14";      
                        elif one_id_lower.endswith('chr15'):   one_id="Chr15";      
                        elif one_id_lower.endswith('chr16'):   one_id="Chr16";      
                        elif one_id_lower.endswith('chr17'):   one_id="Chr17";      
                        elif one_id_lower.endswith('chr18'):   one_id="Chr18";      
                        elif one_id_lower.endswith('chr19'):   one_id="Chr19";      
                        elif one_id_lower.endswith('chr20'):   one_id="Chr20"; 
                        else:continue
                        if one_id!=region_name2chr:
                            print(['error:',one_sample,one_id,region_name2chr]);sys.exit()
                        sample_region_list.append([one_sample,one_region_name])
                        newline=one_sample+'\t'+eachline_arr[0]+'\t'+eachline_arr[1]+'\t'+eachline_arr[4]+'\t'+eachline_arr[5]+'\t'+eachline_arr[7]+'\t'+strand+'\t'+eachline_arr[12]+'\n'
                        f.write(newline)
                        
        print('Loading all monomers')
        for one_sample_region in sample_region_list:
            one_sample=     one_sample_region[0]
            one_region_name=     one_sample_region[1]
            new_name=one_sample+':'+one_region_name
            with open('./samples_satellite/3_samples/monomer/'+new_name,'w') as f:
                f.write('region_name\tcirc_serial\tchrserial_1\tchrserial_2\tstrand\tcirc_seq\n')
                i=0
                with open(input_dir+'/'+one_sample+'/chr_mafft/4_monomer/4.0_plus/'+one_region_name,'r') as f2:
                    for line in f2.readlines():
                        i+=1
                        if i==1:continue
                        eachline_arr=line.strip().split('\t')
                        circ_seq=eachline_arr[-1]
                        if circ_seq in dict_oldcircseq_newcircseq:circ_seq=dict_oldcircseq_newcircseq[circ_seq]
                        newline=eachline_arr[0]+'\t'+eachline_arr[1]+'\t'+eachline_arr[2]+'\t'+eachline_arr[3]+'\t'+eachline_arr[4]+'\t'+circ_seq+'\n'
                        f.write(newline)
    if argv1=="stepall" or argv1=="step3"  or argv1=="step3.0_regionshift":                    
        print("2_good_regions has many fragmented regions. For plotting, I want to align the left ends of all significant blocks. Need to analyze a approximate position.")
        extend_small=100000
        extend_big=1000000        
        dict_samplechr_list={}
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if eachline_arr[0]=='':continue
                sample=         eachline_arr[0]
                chromosome=     eachline_arr[1].replace('region_','Chr')
                samplechr=      sample+'|||'+chromosome
                start=          int(eachline_arr[3])
                end=            int(eachline_arr[4])
                length=         int(eachline_arr[5])
                if samplechr not in dict_samplechr_list: dict_samplechr_list[samplechr]=[]
                dict_samplechr_list[samplechr].append([start,end,length])
        with open('./samples_satellite/2_good_regions_main','w') as f2:   
            f2.write(f'sample\tchromosome\tstart\tend\tlength\n')
            for one_samplechr,infos in dict_samplechr_list.items():
                sample=one_samplechr.split('|||')[0]
                chromosome=one_samplechr.split('|||')[1]
                sorted_infos = sorted(infos, key=lambda x: x[2], reverse=True)
                position1= sorted_infos[0]
                kkk_start=position1[0]
                kkk_end=position1[1]
                state_kkk=1
                #Expand the largest region left and right within 100000bp
                while state_kkk>0:
                    state_kkk=0
                    for one in sorted_infos:
                        one_start=one[0]
                        one_end=one[1]
                        if one_start<kkk_start and kkk_start-one_start<extend_small:
                            kkk_start=one_start
                            state_kkk=1
                        if one_end>kkk_end and kkk_end-one_end<extend_small:
                            kkk_end=one_end
                            state_kkk=1                    
                #Find parts within 1M bp that are larger than 10000bp
                for one in sorted_infos:
                    state_kkk=0
                    one_start=one[0]
                    one_end=one[1]
                    if one_start<kkk_start and kkk_start-one_start<extend_big:
                        kkk_start=one_start
                        state_kkk=1
                    if one_end>kkk_end and kkk_end-one_end<extend_big:
                        kkk_end=one_end
                        state_kkk=1     
                    if state_kkk==1:break    
                #After expansion, extend again with small range to get final coordinates
                state_kkk=1
                while state_kkk>0:
                    state_kkk=0
                    for one in sorted_infos:
                        one_start=one[0]
                        one_end=one[1]
                        if one_start<kkk_start and kkk_start-one_start<extend_small:
                            kkk_start=one_start
                            state_kkk=1
                        if one_end>kkk_end and kkk_end-one_end<extend_small:
                            kkk_end=one_end
                            state_kkk=1  
                f2.write(f'{sample}\t{chromosome}\t{kkk_start}\t{kkk_end}\t{kkk_end-kkk_start+1}\n')
                

                    
    if argv1=="stepall" or argv1=="step3"  or argv1=="step3.1" or argv1=="step3.1a":
        print('step3.1a ———— Without multithreading it is very slow (1 hour). Rewritten to take only 2 minutes. To count 18-19-20        C-C-T   11097900')
        region_name_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if eachline_arr[0]=='':continue
                one_region_name=f"{eachline_arr[0]}:{eachline_arr[1]}{eachline_arr[6]}:{eachline_arr[3]}-{eachline_arr[4]}"
                region_name_list.append(one_region_name)
        region_name_list_num=len(region_name_list)        
        print(region_name_list_num)
        with open('./samples_satellite/3_samples/seat_base_stat_tmp','w') as f:
            f.write('')
        def run_convert(one_region_name):
            dict_seat_base={}
            k=0
            seat_base_list=[]
            with open('./samples_satellite/3_samples/monomer/'+one_region_name,'r') as f:
                for line in f:
                    k+=1
                    if k==1:continue
                    circ_seq_seats=line.strip().split('\t')[-1].split('|')
                    j=0.5
                    for one_seat in circ_seq_seats:
                        j+=0.5
                        if one_seat!='':
                            seat_base_list.append([j,one_seat])
            seat_base_list_len=len(seat_base_list)
            s=1
            while s<seat_base_list_len-1:
                seat_base_1=seat_base_list[s-1]
                seat_base_2=seat_base_list[s]
                seat_base_3=seat_base_list[s+1]
                index=f"{seat_base_1[0]}-{seat_base_2[0]}-{seat_base_3[0]}|{seat_base_1[1]}-{seat_base_2[1]}-{seat_base_3[1]}"
                if index not in dict_seat_base:
                    dict_seat_base[index]=0
                dict_seat_base[index]+=1    
                s+=1
            with open('./samples_satellite/3_samples/seat_base_stat_tmp','a') as f:
                for key,num in dict_seat_base.items():
                    keys=key.split('|')
                    f.write(f"{keys[0].replace('.0','')}\t{keys[1]}\t{str(num)}\n")    
        with Pool(processes=thread) as pool:
            for i, result in enumerate(pool.imap(run_convert, region_name_list), start=1):
                progress = (i / len(region_name_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush() 
            
        dict_seat_base={}     
        with open('./samples_satellite/3_samples/seat_base_stat_tmp','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                seat=eachline_arr[0]
                base=eachline_arr[1]
                index=seat+'|'+base
                num=int(eachline_arr[2])
                if index not in dict_seat_base:
                    dict_seat_base[index]=0
                dict_seat_base[index]+=num
        with open('./samples_satellite/3_samples/seat_base_stat','w') as f:
            for key,num in dict_seat_base.items():
                keys=key.split('|')
                f.write(f"{keys[0]}\t{keys[1]}\t{str(num)}\n")
        subprocess.run(["rm ./samples_satellite/3_samples/seat_base_stat_tmp"], shell=True)        
    if argv1=="stepall" or argv1=="step3"  or argv1=="step3.1" or argv1=="step3.1b":    
        print('step3.1b ———— Count seat-base numbers, time 0s')
        dict_seat_base_num={}
        with open('./samples_satellite/3_samples/seat_base_stat','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')       
                middle_seat=eachline_arr[0].split('-')[1]
                middle_base=eachline_arr[1].split('-')[1]
                middle_base_num=int(eachline_arr[2])
                if middle_seat not in dict_seat_base_num:
                    dict_seat_base_num[middle_seat]={}
                if middle_base not in dict_seat_base_num[middle_seat]:
                    dict_seat_base_num[middle_seat][middle_base]=0
                dict_seat_base_num[middle_seat][middle_base]+=middle_base_num   
        print('Count path (seat-insert-seat) numbers')
        dict_path_insert_num={}
        with open('./samples_satellite/3_samples/seat_base_stat','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t') 
                seat_arr=eachline_arr[0].split('-')
                left_seat=seat_arr[0]
                middle_seat=seat_arr[1]
                right_seat=seat_arr[2]
                path_num=int(eachline_arr[2])
                if '.5' not in middle_seat:
                    path=f"{left_seat}-{middle_seat}"
                    insert_base='-'
                else:
                    if '.5' in right_seat:print('error, right_seat contains .5?');sys.exit()
                    path=f"{left_seat}-{right_seat}"
                    insert_base=eachline_arr[1].split('-')[1]              
                if path not in dict_path_insert_num:
                    dict_path_insert_num[path]={}
                if insert_base not in  dict_path_insert_num[path]:
                     dict_path_insert_num[path][insert_base]=0
                dict_path_insert_num[path][insert_base]+=  path_num       
        with open('./samples_satellite/3_samples/seat_base_stat2tmp','w') as f:
            for seat,y in dict_seat_base_num.items():
                for base,num in y.items():
                    f.write(f"{seat}\t{base}\t{str(num)}\n")
            for path,y in dict_path_insert_num.items():
                for insert,num in y.items():
                    f.write(f"{path}\t{insert}\t{str(num)}\n")    
        subprocess.run(["sort -k 1,1n -k 3,3nr ./samples_satellite/3_samples/seat_base_stat2tmp > ./samples_satellite/3_samples/seat_base_stat2"], shell=True)                
        subprocess.run(["rm ./samples_satellite/3_samples/seat_base_stat2tmp "], shell=True)    
        print('Record the best monomer')   
        last_seat=''
        with open('./samples_satellite/3_samples/seat_base_stat2','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if '-' in eachline_arr[0]:continue
                if '.5' in eachline_arr[0]:continue
                if last_seat=='': 
                    with open('./samples_satellite/3_samples/seat_base_stat3_consensusseq','w') as f2:
                        f2.write(eachline_arr[1])
                    last_seat= eachline_arr[0];continue   
                if eachline_arr[0]!=last_seat:
                    with open('./samples_satellite/3_samples/seat_base_stat3_consensusseq','a') as f2:
                        f2.write(eachline_arr[1])
                last_seat= eachline_arr[0];                        
    #step0 statistics and plotting, length and variability           
    if argv1=="stepall" or argv1=="step3"  or argv1=="step3.2" or argv1=="step3.2a":  
        print('step3.2a ———— Continue processing the 3_samples_monomer_vector file, analyze counts in seats, monomer length, and differences from consensus sequence, time 53s')
        with open('./samples_satellite/3_samples/seat_base_stat3_consensusseq') as f:
            consensusseq=f.read().strip()
        print('Multi-sample consensus sequence: '+consensusseq)
        monomer_len=len(consensusseq)
        print('Calculating variant_distance')    
        dict_variantdistance_num={};dict_circlength_num={}
        with open('./samples_satellite/3_samples/3_samples_monomer_plus','w') as f2:
            f2.write('circ_seq\tin_seat_num\tcirc_len\tnum\tvariant_distance2onehap\tvariant_distance2allhaps\n')
            with open('./samples_satellite/3_samples/3_samples_monomer','r') as f:
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    if eachline_arr[0]=='circ_seq':continue
                    seat_seq_sum_withlines=eachline_arr[0]
                    seat_seq_arr=seat_seq_sum_withlines.split('|')
                    value=int(eachline_arr[-2])
                    kk=0;variant_distance=0;in_seat_num=0
                    while kk<monomer_len:
                        if seat_seq_arr[2*kk]!=consensusseq[kk]:  variant_distance+=1
                        if seat_seq_arr[2*kk+1]!='': variant_distance+=len(seat_seq_arr[2*kk+1])
                        if seat_seq_arr[2*kk]!='':in_seat_num+=1
                        kk+=1
                    circ_length=len(seat_seq_sum_withlines.replace('|',''))
                    if str(variant_distance) not in dict_variantdistance_num:dict_variantdistance_num[str(variant_distance)]=0
                    dict_variantdistance_num[str(variant_distance)]+=value
                    if str(circ_length) not in dict_circlength_num:dict_circlength_num[str(circ_length)]=0
                    dict_circlength_num[str(circ_length)]+=value 
                    f2.write(eachline+'\t'+str(variant_distance)+'\n')
            
        print('Monomer length statistics')
        with open('./samples_satellite/3_samples/stat_Length2num','w') as f:
            f.write('circ_length\tnum\n')
            sorted_dict_circlength_num = dict(sorted(dict_circlength_num.items(), key=lambda item:int(item[0])))
            for circ_length,num in sorted_dict_circlength_num.items():
                f.write(circ_length+'\t'+str(num)+'\n')        
            
        print('Monomer variation statistics')
        with open('./samples_satellite/3_samples/stat_Variant2num','w') as f:
            f.write('variant_distance\tnum\n')
            sorted_dict_variantdistance_num = dict(sorted(dict_variantdistance_num.items(), key=lambda item:int(item[0])))
            for variant_distance,num in sorted_dict_variantdistance_num.items():
                if int(variant_distance)>monomer_len:continue
                f.write(variant_distance+'\t'+str(num)+'\n')                
    if argv1=="stepall" or argv1=="step3"  or argv1=="step3.2" or argv1=="step3.2b": 
        print('step3.2b ———— Plotting, time 1s')
        R_txt='''library(ggplot2)
    
    setwd('./')
    input_file1 <- read.table('stat_Length2num', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
    input_file2 <- read.table('stat_Variant2num', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
    
    #Length plot
    sum_num=sum(input_file1$num)
    p=ggplot(input_file1, aes(x = circ_length, y = num)) +
      geom_bar(stat = "identity")  +scale_x_continuous(limits = c(15, 35))
    p <- p + theme_classic()
    p <- p + labs(
      title = "",
      x = paste("Circ length (All =", sum_num, ")", sep = ""),
      y = "Frequency"
    )
    
    pdf("pic1_stat_Length2num.pdf", width = 10 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()
    
    #Length plot percentage
    sum_num=sum(input_file1$num)
    p=ggplot(input_file1, aes(x = circ_length, y = num/sum_num)) +
      geom_bar(stat = "identity",width = 0.7)  +
      scale_x_continuous(limits = c(19, 31))+
      scale_y_continuous(labels = scales::percent)
    p <- p + theme_classic()
    p <- p + labs(
      title = "",
      x = paste("Circ length (All =", sum_num, ")", sep = ""),
      y = "Percentage"
    )
    
    pdf("pic1_stat_Length2num_percent.pdf", width = 4 / 2.54, height = 4 / 2.54)
    print(p)
    dev.off()
    
    #Variation plot
    sum_num=sum(input_file2$num)
    p=ggplot(input_file2, aes(x = variant_distance, y = num)) +
      geom_bar(stat = "identity")  +scale_x_continuous(limits = c(1, 20))
    p <- p + theme_classic()
    p <- p + labs(
      title = "", 
      x = paste("Variant Distance (All =", sum_num, ")", sep = ""),
      y = "Frequency"
    )
    
    pdf("pic2_stat_Variant2num.pdf", width = 10 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()
    
    #Variation plot percentage
    sum_num=sum(input_file2$num)
    p=ggplot(input_file2, aes(x = variant_distance, y = num/sum_num)) +
      geom_bar(stat = "identity",width = 0.7)  +
      scale_x_continuous(limits = c(0, 17))+
      scale_y_continuous(labels = scales::percent)
    p <- p + theme_classic()
    p <- p + labs(
      title = "", 
      x = paste("Variant Distance (All =", sum_num, ")", sep = ""),
      y = "Percentage"
    )
    pdf("pic2_stat_Variant2num_percent.pdf", width = 4 / 2.54, height = 4 / 2.54)
    print(p)
    dev.off()
    '''
        with open('./samples_satellite/3_samples/pic_ggplot_step3.2b.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = "./samples_satellite/3_samples/"
        os.chdir(new_directory)
        subprocess.run(['Rscript pic_ggplot_step3.2b.R'], shell=True)    
        os.chdir('../../')                   
    if argv1=="stepall" or argv1=="step3"  or argv1=="step3.2" or argv1=="step3.2c":          
        print('step3.2c ———— Calculate consensus sequences for 28bp and 23bp separately')
        dict_28seat_base_num={}
        dict_23seat_base_num={}
        with open('./samples_satellite/3_samples/3_samples_monomer_plus','r') as f:
            for line in f:
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                circ_seq=eachline_arr[0].replace('|','')
                if circ_seq=='circ_seq':continue
                circ_len=eachline_arr[2]
                num=int(eachline_arr[3])   
                
                if circ_len=="28":
                    for index,base in enumerate(circ_seq,start=0):
                        if index not in dict_28seat_base_num:dict_28seat_base_num[index]={}
                        if base not in dict_28seat_base_num[index]:dict_28seat_base_num[index][base]=0
                        dict_28seat_base_num[index][base]+=num
                    
                if circ_len=="23":
                    for index,base in enumerate(circ_seq,start=0):
                        if index not in dict_23seat_base_num:dict_23seat_base_num[index]={}
                        if base not in dict_23seat_base_num[index]:dict_23seat_base_num[index][base]=0
                        dict_23seat_base_num[index][base]+=num            
        k=0
        consensusseq28=''
        while k<28:
            dict_base_num=dict_28seat_base_num[k]
            arr_base_num_sorted = list(dict(sorted(dict_base_num.items(), key=lambda item: item[1], reverse=True)))
            consensusseq28+=arr_base_num_sorted[0]
            k+=1
        with open('./samples_satellite/3_samples/seat_base_stat4_consensusseq28','w') as f:
            f.write(consensusseq28)
        print('consensusseq28：'+consensusseq28)
        
        k=0
        consensusseq23=''
        while k<23:
            dict_base_num=dict_23seat_base_num[k]
            arr_base_num_sorted = list(dict(sorted(dict_base_num.items(), key=lambda item: item[1], reverse=True)))
            consensusseq23+=arr_base_num_sorted[0]
            k+=1
        with open('./samples_satellite/3_samples/seat_base_stat4_consensusseq23','w') as f:
            f.write(consensusseq23)        
        print('consensusseq23：'+consensusseq23)                
    if argv1=="stepall" or argv1=="step3"  or argv1=="step3.2" or argv1=="step3.2d":          
        print('step3.2d ———— Calculate consensus sequences for 28bp and 23bp separately, then compute variants')     
    
        with open('./samples_satellite/3_samples/seat_base_stat4_consensusseq28','r') as f:
            consensusseq28=f.read().strip()
        with open('./samples_satellite/3_samples/seat_base_stat4_consensusseq23','r') as f:
            consensusseq23=f.read().strip()        
        with open('./samples_satellite/3_samples/3_samples_monomer28','w') as f2:
            f2.write('circ_seq\tin_seat_num\tcirc_len\tnum\tvariant28_distance\n')
        with open('./samples_satellite/3_samples/3_samples_monomer23','w') as f3:
            f3.write('circ_seq\tin_seat_num\tcirc_len\tnum\tvariant23_distance\n')     
        dict_seq_num_distance={}
        with open('./samples_satellite/3_samples/3_samples_monomer_plus','r') as f:
            for line in f:
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                circ_seq        =eachline_arr[0]
                if circ_seq=='circ_seq':continue
                circ_seq_pure=circ_seq.replace('|','')        
                in_seat_num     =eachline_arr[1]
                circ_len        =eachline_arr[2]
                num             =eachline_arr[3]
                if circ_len=="28":  
                    variant28_distance=0
                    for i,x in enumerate(circ_seq_pure):
                        if x!=consensusseq28[i]:variant28_distance+=1
                    with open('./samples_satellite/3_samples/3_samples_monomer28','a') as f2:
                        f2.write(f'{circ_seq}\t{in_seat_num}\t{circ_len}\t{num}\t{variant28_distance}\n')
                if circ_len=="23":  
                    variant23_distance=5    #Starting from 5
                    for i,x in enumerate(circ_seq_pure):
                        if x!=consensusseq23[i]:variant23_distance+=1
                    circ_seq_pure=circ_seq.replace('|','')
                    if circ_seq_pure not in dict_seq_num_distance:
                        dict_seq_num_distance[circ_seq_pure]={}
                        dict_seq_num_distance[circ_seq_pure][0]=0
                        dict_seq_num_distance[circ_seq_pure][1]=''
                    dict_seq_num_distance[circ_seq_pure][0]+=int(num)
                    dict_seq_num_distance[circ_seq_pure][1]=variant23_distance
                    with open('./samples_satellite/3_samples/3_samples_monomer23','a') as f3:
                        f3.write(f'{circ_seq}\t{in_seat_num}\t{circ_len}\t{num}\t{variant23_distance}\n')       
     
        
        arr_seq_num_distance_sorted = list(dict(sorted(dict_seq_num_distance.items(), key=lambda item: item[1][0], reverse=True)).items())       
        with open('./samples_satellite/3_samples/3_samples_monomer23_revise','w') as f:
            f.write('circ_seq\tin_seat_num\tcirc_len\tnum\tvariant23_distance\n')
            for seq_num1var in  arr_seq_num_distance_sorted:
                seq=                    seq_num1var[0]
                num=                    seq_num1var[1][0]
                variant23_distance=     seq_num1var[1][1]
                f.write(f'{seq}\t{str(num)}\t{variant23_distance}\n')
    if argv1=="stepall" or argv1=="step3"  or argv1=="step3.2" or argv1=="step3.2e":    
        print('step3.2e ———— Merge length-distance statistics for 28bp and 23bp, time 2s')     
        dict28_distance_num={}
        with open('./samples_satellite/3_samples/3_samples_monomer28','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if  eachline_arr[0]=='circ_seq':continue
                distance=   eachline_arr[-1]
                num=        eachline_arr[-2]
                if distance not in dict28_distance_num:
                    dict28_distance_num[distance]=0
                dict28_distance_num[distance]+=int(num)
        dict23_distance_num={}
        with open('./samples_satellite/3_samples/3_samples_monomer23_revise','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if  eachline_arr[0]=='circ_seq':continue
                distance=   eachline_arr[-1]
                num=        eachline_arr[-2]
                if distance not in dict23_distance_num:
                    dict23_distance_num[distance]=0
                dict23_distance_num[distance]+=int(num)
        with open('./samples_satellite/3_samples/stat_Seattype2Variant2num','w') as f:
            f.write('type\tvariant_distance\tnum\n')
            for distance,num in dict28_distance_num.items():
                f.write(f'seat28\t{distance}\t{num}\n')
            for distance,num in dict23_distance_num.items():
                f.write(f'seat23\t{distance}\t{num}\n')         
    if argv1=="stepall" or argv1=="step3"  or argv1=="step3.2" or argv1=="step3.2f":    
        print('step3.2f ———— Plotting, time 1s')
        R_txt='''library(ggplot2)
    
    setwd('./')
    input_file <- read.table('stat_Seattype2Variant2num', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
    input_file2 <- read.table('stat_Variant2num', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
    
    sum_num <- sum(input_file$num)
    sum_num2=sum(input_file2$num)
    
    library(ggplot2)
    
    p <- ggplot() +
        geom_col(data = input_file2, aes(x = variant_distance, y = num), fill = "black")+
        geom_line(data = input_file, aes(x = variant_distance, y = num, color = type, group = type)) +
        geom_point(data = input_file, aes(x = variant_distance, y = num, color = type))
    p <- p + theme_classic()
    p <- p + labs(
      title = "", 
      x = paste("Variant Distance (All =", sum_num2, ")", sep = ""),
      y = "Frequency"
    )
     
    pdf("pic3_stat_Seattype2Variant2num.pdf", width = 10 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()
    
    
    p <- ggplot(input_file, aes(x = variant_distance, y = num, color = type, group = type)) +
      geom_line() +
      geom_point()
     
    p <- p + theme_classic()
    p <- p + labs(
      title = "", 
      x = paste("Variant Distance (All =", sum_num, ")", sep = ""),
      y = "Frequency"
    )
     
    pdf("pic3_stat_Seattype2Variant2num2.pdf", width = 10 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()
    '''
        with open('./samples_satellite/3_samples/pic_ggplot_step3.2f.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = "./samples_satellite/3_samples/"
        os.chdir(new_directory)
        subprocess.run(['Rscript pic_ggplot_step3.2f.R'], shell=True)    
        os.chdir('../../')                   
    #step3 statistics and plotting, boxes, seat connections 
    if argv1=="stepall" or argv1=="step3"  or argv1=="step3.3" or argv1=="step3.3a":
        print('step3.3a')
        input_file='./samples_satellite/3_samples/seat_base_stat2'
        output_file='./samples_satellite/3_samples/seat_base_stat2_plot_1a'    
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
                seatpath =eachline_arr[0]; 
                if '.5' in seatpath:continue
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
            A_num_width=math.sqrt(A_num)/num_max_sqrt*max_base_width
            C_num_width=math.sqrt(C_num)/num_max_sqrt*max_base_width
            G_num_width=math.sqrt(G_num)/num_max_sqrt*max_base_width
            T_num_width=math.sqrt(T_num)/num_max_sqrt*max_base_width
            onemax_base_width=max(A_num_width,C_num_width,G_num_width,T_num_width)
            seat_base_pic_width=max(A_num_width+C_num_width,G_num_width+T_num_width)
            seat_base_pic_height=max(A_num_width+G_num_width,C_num_width+T_num_width)
            offset_x_abs=onemax_base_width-seat_base_pic_width/2
            if      A_num_width==onemax_base_width:    offset_x=offset_x_abs;       offset_y=-offset_x_abs 
            elif    C_num_width==onemax_base_width:    offset_x=-offset_x_abs;      offset_y=-offset_x_abs 
            elif    G_num_width==onemax_base_width:    offset_x=offset_x_abs;       offset_y=offset_x_abs 
            elif    T_num_width==onemax_base_width:    offset_x=-offset_x_abs;      offset_y=offset_x_abs         
            angle_radians = math.radians(each_angle*i)    
            cosine_value = math.cos(angle_radians);     one_x = bigcirc_r_length*cosine_value
            sine_value = math.sin(angle_radians);       one_y = bigcirc_r_length*sine_value
            one_x_real=one_x+offset_x
            one_y_real=one_y+offset_y
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
            i+=1
        path_list=list(dict_path_num.keys())
        def custom_sort(item):
            key_part = item.split('\t')[0]
            left, right = key_part.split('-')
            return (int(left), int(right))
        sorted_path_list = sorted(path_list, key=custom_sort)     
        print('Output stat results')  
    if argv1=="stepall" or argv1=="step3"  or argv1=="step3.3" or argv1=="step3.3b":      
        print('step3.3b')
        input_file='./samples_satellite/3_samples/seat_base_stat'
        input_file2='./samples_satellite/3_samples/seat_base_stat2_plot_1a'   
        output_file='./samples_satellite/3_samples/seat_base_stat2_plot_1b'
        dict_newpath_num={}
        dict_newpath_seq_num={}
        newpath_seq_num_max=0
        with open(input_file,'r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                path     =eachline_arr[0]
                if ".5" in path:continue
                seq      =eachline_arr[1];  
                bad_mark=''
                for one in seq:
                    if one not in ['A','C','G','T','-']:bad_mark='yes';break
                if bad_mark=='yes':continue   
                num      =int(eachline_arr[2]) 
                newpath_seq_num_max+=num
                path_arr=   path.split('-')
                seq_arr=    seq.split('-') 
                new_path=path_arr[1]+seq_arr[1]+'-'+path_arr[2]+seq_arr[2]
                if new_path not  in dict_newpath_num:
                    dict_newpath_num[new_path]=0
                dict_newpath_num[new_path]+=num 
                if new_path not  in dict_newpath_seq_num:
                    dict_newpath_seq_num[new_path]={}
                mid_seq='-'
                if mid_seq not in dict_newpath_seq_num[new_path]:
                    dict_newpath_seq_num[new_path][mid_seq]=0
                dict_newpath_seq_num[new_path][mid_seq]+=num
          
        newpath_list=list(dict_newpath_num.keys())
        def custom_sort(item):
            key_part = item.split('\t')[0]
            left, right = key_part.split('-')
            return (int(left[:-1]), int(right[:-1]))
        sorted_newpath_list = sorted(newpath_list, key=custom_sort)     
        with open(output_file,'w') as f2:
            f2.write('name\tflypast\tpath1\tpath1_x\tpath1_y\tpath2\tpath2_x\tpath2_y\tnum\tline_width\topacity\tinfo\n')
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
        px1_num=newpath_seq_num_max/seat_max/10
        for one_path in sorted_newpath_list:
            one_path_arr=one_path.split('-')
            path1=one_path_arr[0];  path1_x=dict_reviseseat_pos[path1][0];  path1_y=dict_reviseseat_pos[path1][1];
            path2=one_path_arr[1];  path2_x=dict_reviseseat_pos[path2][0];  path2_y=dict_reviseseat_pos[path2][1];
            seat1=int(path1[:-1])
            seat2=int(path2[:-1])
            if seat2>=seat1: flypast=seat2-seat1-1
            else:flypast=seat_max-seat1+seat2-1
            newpath_num=dict_newpath_num[one_path]
            line_width=round(newpath_num/px1_num,1)
            if line_width<1:line_width=1
            opacity=round(newpath_num/px1_num,3)    
            if opacity>0.8:opacity=0.8
            dict_path_seq=dict_newpath_seq_num[one_path]
            sorteddict_path_seq = dict(sorted(dict_path_seq.items(), key=lambda item: item[1], reverse=True)) 
            oneinfo=''
            for key,value in sorteddict_path_seq.items():
                oneinfo+=(key+'|'+str(value)+';')
            with open(output_file,'a') as f2:
                f2.write(one_path+'\t'+str(flypast)+'\t'+path1+'\t'+path1_x+'\t'+path1_y+'\t'+path2+'\t'+path2_x+'\t'+path2_y+'\t'+str(newpath_num)+'\t'+str(line_width)+'\t'+str(opacity)+'\t'+oneinfo+'\n')         
    if argv1=="stepall" or argv1=="step3"  or argv1=="step3.3" or argv1=="step3.3c":        
        print('step3.3c')
        input_file='./samples_satellite/3_samples/seat_base_stat2_plot_1a'
        input_file2='./samples_satellite/3_samples/seat_base_stat2_plot_1b'   
        output_file='./samples_satellite/3_samples/'
        output_file1=output_file+'seat_base_stat2_plot_1c'      
        output_file2=output_file+'seat_base_stat2_plot_1d'  
        output_file3=output_file+'seat_base_stat2_plot_1e'      
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
            max_base_width=40
            num_max_sqrt=math.sqrt(num_max)
            A_num_width=math.sqrt(A_num)/num_max_sqrt*max_base_width
            C_num_width=math.sqrt(C_num)/num_max_sqrt*max_base_width
            G_num_width=math.sqrt(G_num)/num_max_sqrt*max_base_width
            T_num_width=math.sqrt(T_num)/num_max_sqrt*max_base_width
            onemax_base_width=max(A_num_width,C_num_width,G_num_width,T_num_width)
            seat_base_pic_width=max(A_num_width+C_num_width,G_num_width+T_num_width)
            seat_base_pic_height=max(A_num_width+G_num_width,C_num_width+T_num_width)
            offset_x_abs=onemax_base_width-seat_base_pic_width/2
            if      A_num_width==onemax_base_width:    offset_x=offset_x_abs;       offset_y=-offset_x_abs 
            elif    C_num_width==onemax_base_width:    offset_x=-offset_x_abs;      offset_y=-offset_x_abs 
            elif    G_num_width==onemax_base_width:    offset_x=offset_x_abs;       offset_y=offset_x_abs 
            elif    T_num_width==onemax_base_width:    offset_x=-offset_x_abs;      offset_y=offset_x_abs 
            one_x=70*i
            one_y=180
            one_x_real=one_x+offset_x
            one_y_real=one_y+offset_y  
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
                if path2[:-1] =='1':
                    path2=str(seat_max+1)+eachline_arr[5][-1];
                    name=path1+'-'+path2
                path1_x=str(dict_seat_pos[path1][0])
                path1_y=str(dict_seat_pos[path1][1])
                path2_x=str(dict_seat_pos[path2][0])
                path2_y=str(dict_seat_pos[path2][1])
                with open(output_file3,'a') as f3:
                    f3.write(name+'\t'+flypast+'\t'+path1+'\t'+path1_x+'\t'+path1_y+'\t'+path2+'\t'+path2_x+'\t'+path2_y+'\t'+num+'\t'+line_width+'\t'+opacity+'\t'+info+'\n')
    if argv1=="stepall" or argv1=="step3"  or argv1=="step3.3" or argv1=="step3.3d":
        print('step3.3d')
        R_txt='''library(ggplot2)
    library(dplyr)
    setwd('./')
    
    input_file1 <- read.table('seat_base_stat2_plot_1a', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
    input_file2 <- read.table('seat_base_stat2_plot_1b', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
    
    input_file1_simple<-unique.data.frame(input_file1[,c(1,3,4)])
    
    sorted_input_file2<-input_file2[order(-input_file2[,9]),]
    filtered_input_file2<- sorted_input_file2[1:500, ]
    
    flypast_small<-filtered_input_file2%>% filter(flypast==0)
    flypast_big<-filtered_input_file2%>% filter(flypast>0)
    
    sum_num=sum(input_file1$num)
    flypast_sum_num=sum(flypast_big$num)
    count_dash <- 0
    count_other <- 0
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
    
    flypast_allbig<-input_file2%>% filter(flypast>0)
    flypast_allbig_simple<-flypast_allbig[,c(1,2,9)]
    flypast_allbig_simple_sorted<-flypast_allbig_simple[order(-flypast_allbig_simple[,3]),]
    flypast_allbig_simple_sorted_15<-flypast_allbig_simple_sorted[1:15,]
    flypast_allbig_simple_sorted_15a <- flypast_allbig_simple_sorted[16:nrow(flypast_allbig_simple_sorted),]
    sum15a=sum(flypast_allbig_simple_sorted_15a$num)
    flypast_allbig_simple_sorted_15other <- flypast_allbig_simple_sorted_15 %>%
      add_row(flypast = NA, num = sum15a, name = 'Other')
    flypast_allbig_simple_sorted_15other <- flypast_allbig_simple_sorted_15other %>%  
      dplyr::mutate(name_order = row_number())
    
    input_file1$xmin <- input_file1$x - input_file1$radius
    input_file1$ymin <- input_file1$y - input_file1$radius
    input_file1$xmax <- input_file1$x + input_file1$radius
    input_file1$ymax <- input_file1$y + input_file1$radius
    
    color_A <- '#009e73'
    color_C <- '#56b4e9'
    color_G <- '#e69f00'
    color_T <- '#cc79a7'
    
    p <- ggplot()
    p <- p + geom_segment(data = flypast_big, 
                          aes(x = path1_x, y = path1_y, xend = path2_x, yend = path2_y, alpha = opacity), 
                          linewidth=3,
                          color = "black",
                          arrow = arrow(length = unit(0, "npc"), type = "closed", ends = "last"))
    p <- p + geom_segment(data = flypast_small, 
                          aes(x = path1_x, y = path1_y, xend = path2_x, yend = path2_y, alpha = opacity), 
                          linewidth=0.1,
                          color = "black",
                          arrow = arrow(length = unit(0, "npc"), type = "closed", ends = "last"))
    p <- p + geom_rect(data = input_file1, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = base))
    p <- p + scale_fill_manual(values = c(A = color_A, C = color_C, G = color_G, T = color_T))
    p <- p + geom_text(data = input_file1_simple, aes(x = seat_x_center*1.12, y = seat_y_center*1.12, label = seat),
                       color = "black", size = 3, vjust = 0.5, hjust = 0.5)
    p <- p + theme_classic()
    p <- p+coord_equal()
    p <- p + labs(title = "Plot of Rectangles with Segments", x =  paste("Number (All =", sum_num, ")", sep = ""), y = "Y Coordinate")
    p <- p + theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )
    
    pdf("pic4_circle-connect.pdf", width = 10 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()
    
    desired_order=flypast_allbig_simple_sorted_15other$name
    desired_order=rev(desired_order)
    flypast_allbig_simple_sorted_15other$name <- factor(flypast_allbig_simple_sorted_15other$name, levels = desired_order)
    custom_colors <- c(
      "1" = "#cc9900",
      "5" = "#4ab1d3",
      "6" = "#ffad99",
      "NA" = "#8c8c8c"
    )
    flypast_allbig_simple_sorted_15other$flypast <- factor(flypast_allbig_simple_sorted_15other$flypast)
    p <- ggplot(data = flypast_allbig_simple_sorted_15other, aes(x = name, y = num,fill=flypast)) +
      geom_col() +
      labs(title = "Bar Chart of num by name", x = "Name", y = paste("Number (All =", flypast_sum_num, ")", sep = "")) +
      theme_classic() +
      coord_flip()+
      scale_fill_manual(values = custom_colors)
    
    pdf("pic5_circle-connect_flypast.pdf", width = 10 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()
    '''
        with open('./samples_satellite/3_samples/pic_ggplot_step3.3d.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = "./samples_satellite/3_samples/"
        os.chdir(new_directory)
        subprocess.run(['Rscript pic_ggplot_step3.3d.R'], shell=True)    
        os.chdir('../../')
    if argv1=="stepall" or argv1=="step3" or argv1=="step3.3" or argv1=="step3.3e":
        print('step3.3e')
        R_txt='''library(ggplot2)
    library(dplyr)
    
    setwd('./')
    input_file1 <- read.table('seat_base_stat2_plot_1c', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
    input_file2 <- read.table('seat_base_stat2_plot_1d', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
    input_file3 <- read.table('seat_base_stat2_plot_1e', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
    
    input_file2=input_file2%>%mutate(alpha=0.5*(1-opacity))
    
    sorted_input_file3<-input_file3[order(-input_file3[,9]),]
    filtered_input_file3<- sorted_input_file3[1:500, ]
    input_file3_big<-filtered_input_file3%>% filter(opacity==0.8)
    input_file3_small<-filtered_input_file3%>% filter(opacity<0.8)
    
    seat_pos_simple=unique.data.frame(input_file2[,c(1,4)])
    
    seat_maxbase <- input_file2[,c(1,2,3,4)] %>%
      group_by(seat) %>%
      arrange(desc(num)) %>%
      slice(1) %>%
      ungroup()
    
    input_file1$xmin <- input_file1$x - input_file1$radius
    input_file1$ymin <- input_file1$y - input_file1$radius
    input_file1$xmax <- input_file1$x + input_file1$radius
    input_file1$ymax <- input_file1$y + input_file1$radius
    
    input_file2$xmin <- input_file2$x - 15
    input_file2$ymin <- input_file2$y - 15
    input_file2$xmax <- input_file2$x + 15
    input_file2$ymax <- input_file2$y + 15
    
    color_A <- '#009e73'
    color_C <- '#56b4e9'
    color_G <- '#e69f00'
    color_T <- '#cc79a7'
    
    p <- ggplot()
    p <- p + geom_segment(data = input_file3_small, 
                          aes(x = path1_x, y = path1_y, xend = path2_x, yend = path2_y,alpha=opacity), 
                          linewidth =1,
                          color = "#336699",
    )
    
    p <- p + geom_segment(data = input_file3_big, 
                          aes(x = path1_x, y = path1_y, xend = path2_x, yend = path2_y,linewidth=line_width), 
                          alpha = 0.8,
                          color = "#336699",
    )
    
    p <- p + geom_rect(data = input_file1, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = base))
    p <- p + geom_rect(data = input_file2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = base))
    p <- p + geom_rect(data = input_file2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,alpha=alpha), fill = "white") 
    p <- p + scale_fill_manual(values = c(A = color_A, C = color_C, G = color_G, T = color_T))
    p <- p + geom_text(data = seat_pos_simple, aes(x = x, y = -50, label = seat),
                       color = "black", size = 3, vjust = 0.5, hjust = 0.5)
    p <- p + geom_text(data = seat_maxbase, aes(x = x, y = -10, label = base),
                       color = "black", size = 3, vjust = 0.5, hjust = 0.5)
    p <- p + theme_classic()
    p <- p+coord_equal()
    p <- p + labs(title = "Plot of Rectangles with Segments", x = "X Coordinate", y = "Y Coordinate")
    p <- p + theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )
    
    pdf("pic6_base-matrixs-connections.pdf", width = 60 / 2.54, height = 15 / 2.54)
    print(p)
    dev.off()
    '''
        with open('./samples_satellite/3_samples/pic_ggplot_step3.3e.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = "./samples_satellite/3_samples/"
        os.chdir(new_directory)
        subprocess.run(['Rscript pic_ggplot_step3.3e.R'], shell=True)  
        os.chdir('../../')
        print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))

"""
2_good_regions
       ↓
step3.0_plot_centromere_pos → sum_for_plot.pdf
       ↓
step3.0_regionshift2 → 2_good_regions_main_66_103
       ↓
step3.0_block_type → 2_good_regions_interarray
       ↓
step3.0_block_type_LIR → 2_good_regions_interarray_LIR 
       ↓
step3.0_block_type_LIR_strand → 2_good_regions_interarray_LIR_strand
"""
if argv1=="step3.0_plot_centromere_pos" :
    print('Plot centromere positions')
    if  os.path.exists('./samples_satellite/plot_centromere_pos')==False:
        subprocess.run(["mkdir ./samples_satellite/plot_centromere_pos"], shell=True)  
    with open ("./samples_satellite/plot_centromere_pos/input",'w')     as f2:
        f2.write("sample\tchromosome\tstart\tend\tlength\tchr_length\n")
        with open ("./samples_satellite/2_good_regions_main",'r')     as f:
            next(f)
            for line in f:
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                sample=eachline_arr[0]
                chromosome=eachline_arr[1]
                length2=0
                with open(f'/home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample}/{sample}.fasta.fai','r') as ff:
                    for line2 in ff:
                        eachline_arr2=line2.split('\t')
                        chromosome2,length2=eachline_arr2[:2]
                        if chromosome2==chromosome:
                            length=int(length2);break
                if length2==0:print('error')       ;sys.exit() 
                f2.write(f"{eachline}\t{length2}\n")
                
    with open ("./samples_satellite/plot_centromere_pos/input2",'w')     as f2:
        f2.write("sample\tchromosome\tstart\tend\tlength\tchr_length\n")
        with open ("./samples_satellite/2_good_regions",'r')     as f:
            next(f)
            for line in f:
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                if len(eachline_arr)!=8:continue
                sample=eachline_arr[0]
                chromosome=eachline_arr[1].replace('region_','Chr')
                start=eachline_arr[3]
                end=eachline_arr[4]
                length=eachline_arr[5]
                length2=0
                with open(f'/home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample}/{sample}.fasta.fai','r') as ff:
                    for line2 in ff:
                        eachline_arr2=line2.split('\t')
                        chromosome2,length2=eachline_arr2[:2]
                        if chromosome2==chromosome:
                            length=int(length2);break
                if length2==0:print('error')       ;sys.exit() 
                f2.write(f"{sample}\t{chromosome}\t{start}\t{end}\t{length}\t{length2}\n")    
    Plot_txt=r"""
    library(ggplot2)
    library(dplyr)
    #Get all command line arguments
    
    ###Monomer
    print("")
    {
      input_file=read.table('input2', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')

      input_file$chromosome <- factor(
          input_file$chromosome,
          levels = c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10','Chr11','Chr12','Chr13','Chr14','Chr15','Chr16','Chr17','Chr18','Chr19')
        )  
        selected_cols <- input_file[, c(1, 2, 6)]
        unique_data <- unique(selected_cols)
        
      color_values=c(

         )
          
      p = ggplot() 
    
    p =p+  geom_rect(data = unique_data,aes(xmin = 0,xmax = chr_length/1000000, ymin = 1,ymax =  2),fill = 'grey')   
    p =p+  geom_rect(data = input_file,aes(xmin = start/1000000,xmax = end/1000000, ymin = 1,ymax =  2),fill = 'blue')   
    

      ##
      ###################################
    p=p+facet_grid(sample ~ chromosome)
                       
        
       p=p+theme_classic() +         
        theme(
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          #legend.position = "none",
          #axis.text.x = element_blank()
        )# +
       # scale_fill_manual(values = color_values, drop = FALSE)
      
      pdf(file = paste0('sum_for_plot', ".pdf"), width = 240 / 2.54, height = 80 / 2.54)
      print(p)
      dev.off()
      
    }            
    """
    with open(f'./samples_satellite/plot_centromere_pos/plot.R','w') as f:
        f.write(Plot_txt)
    os.chdir(f'./samples_satellite/plot_centromere_pos/')
    subprocess.run([f'Rscript plot.R  '], shell=True)
    os.chdir('../') 
if argv1=="step3.0_regionshift2":                    
    print("2_good_regions has many fragmented regions. For plotting, I want to align the left ends of all significant blocks. Need to analyze a approximate position.")
    block_size=10000
    extend_small=100000
    extend_big=1000000        
    dict_samplechr_list={}
    with open('./stat_plot/0-region2info','r') as f:
        #sample	chromosome	chromosome_new	centype	chr_start	chr_end	length	strand	match_percent
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            if eachline_arr[0]=='':continue
            sample=         eachline_arr[0]
            chromosome=     eachline_arr[2]
            centype=        eachline_arr[3]
            if centype not in ['cen_107','cen_66','cen_103']: continue
            samplechr=      sample+'|||'+chromosome
            start=          int(eachline_arr[4])
            end=            int(eachline_arr[5])
            length=         int(eachline_arr[6])
            if samplechr not in dict_samplechr_list: dict_samplechr_list[samplechr]=[]
            dict_samplechr_list[samplechr].append([start,end,length])
    with open('./samples_satellite/2_good_regions_main_66_103','w') as f2:   
        f2.write(f'sample\tchromosome\tstart\tend\tlength\n')
        for one_samplechr,infos in dict_samplechr_list.items():
            sample=one_samplechr.split('|||')[0]
            chromosome=one_samplechr.split('|||')[1]
            sorted_infos = sorted(infos, key=lambda x: x[2], reverse=True)
            position1= sorted_infos[0]
            kkk_start=position1[0]
            kkk_end=position1[1]
            state_kkk=1
            #Expand the largest region left and right within 100000bp
            while state_kkk>0:
                state_kkk=0
                for one in sorted_infos:
                    one_start=one[0]
                    one_end=one[1]
                    if one_start<kkk_start and kkk_start-one_start<extend_small:
                        kkk_start=one_start
                        state_kkk=1
                    if one_end>kkk_end and kkk_end-one_end<extend_small:
                        kkk_end=one_end
                        state_kkk=1  
            n=10
            while n>0:
                n-=1
                #Find parts within 1M bp that are larger than 10000bp
                for one in sorted_infos:
                    state_kkk=0
                    one_start=one[0]
                    one_end=one[1]
                    if one_start<kkk_start and kkk_start-one_start<extend_big:
                        kkk_start=one_start
                        state_kkk=1
                    if one_end>kkk_end and kkk_end-one_end<extend_big:
                        kkk_end=one_end
                        state_kkk=1     
                    if state_kkk==1:break    
                #After expansion, extend again with small range to get final coordinates
                state_kkk=1
                while state_kkk>0:
                    state_kkk=0
                    for one in sorted_infos:
                        one_start=one[0]
                        one_end=one[1]
                        if one_start<kkk_start and kkk_start-one_start<extend_small and abs(one_end-one_start)>block_size:
                            kkk_start=one_start
                            state_kkk=1
                        if one_end>kkk_end and kkk_end-one_end<extend_small and abs(one_end-one_start)>block_size:
                            kkk_end=one_end
                            state_kkk=1  
                     
            f2.write(f'{sample}\t{chromosome}\t{kkk_start}\t{kkk_end}\t{kkk_end-kkk_start+1}\n')                
if argv1=="step3.0_block_type":                    
    print("VSat1 and array in 2_good_regions. Need to analyze the start points of VSat1 on a chromosome, internal points, and inter-array positions.")
    dict_samplechr_list={}
    with open('./samples_satellite/2_good_regions','r') as f:
        for line in f:
            eachline_arr=line.strip().split('\t')
            if eachline_arr[0]=='sample':continue
            if eachline_arr[0]=='':continue
            sample=         eachline_arr[0]
            chromosome=     eachline_arr[1].replace('region_','Chr')
            samplechr=      sample+'|||'+chromosome
            start=          int(eachline_arr[3])
            end=            int(eachline_arr[4])
            length=         int(eachline_arr[5])
            if samplechr not in dict_samplechr_list: dict_samplechr_list[samplechr]=[]
            dict_samplechr_list[samplechr].append([start,end,length])  
    print('Processing data')        
    with open('./samples_satellite/2_good_regions_interarray','w') as ff:
        ff.write(f"sample\tchromosome\tbigblock_start\tbigblock_end\tall_VSat1_length\tother_VSat1_len\tcore_VSat1_length\tinter_array_num\tinter_array_length\tinter_arrays\tlargest_part_start\tlargest_part_end\tbigblock_part_str\n")
        for  one_samplechr,infos in dict_samplechr_list.items():
            sample,chromosome=one_samplechr.split('|||')
            print(sample,chromosome)
            
            sorted_infos_pos=     sorted(infos, key=lambda x: x[0], reverse=False)  
            
            sorted_infos_length = sorted(infos, key=lambda x: x[2], reverse=True)
            max_array_start=sorted_infos_length[0][0]
            max_array_end=sorted_infos_length[0][1]
            VSat1_num=len(sorted_infos_pos)
            iii=0;bad_iii_set=set()
            good_iii=''
            while iii<VSat1_num:
                start0,end0,len0=sorted_infos_pos[iii]
                if start0-max_array_end>0:  direction='VSat1 is to the right of the largest block'
                elif max_array_start-end0>0:  direction='VSat1 is to the left of the largest block'
                else: direction='Current is the largest block';good_iii=iii
                bad_mark=''
                if direction=='VSat1 is to the right of the largest block':
                    for one_info in sorted_infos_pos:
                        tmp_start,tmp_end,tmp_length=one_info
                        if tmp_start>max_array_end  and tmp_end< start0 and start0-tmp_end>1000000:bad_mark='yes';break
                if direction=='VSat1 is to the left of the largest block':
                    for one_info in sorted_infos_pos:
                        tmp_start,tmp_end,tmp_length=one_info
                        if tmp_start<max_array_start  and tmp_start> end0 and tmp_start-end0>1000000:bad_mark='yes';break
                if bad_mark=='yes':bad_iii_set.add(iii)
                iii+=1

            print('Initial screening completed')    
            iii=0
            bigblock_start=max_array_start
            bigblock_end=max_array_end
            other_VSat1_len=0
            all_VSat1_length=0
            while iii<VSat1_num:
                start0,end0,len0=sorted_infos_length[iii]
                all_VSat1_length+=len0
                if iii in bad_iii_set:
                    iii+=1;
                    other_VSat1_len+=one_length;
                    continue;
                iii+=1
                #Based on experience, discard those more than 5M away from the largest array
                if start0-max_array_end>5000000:other_VSat1_len+=one_length;continue
                if max_array_start-end0>5000000:other_VSat1_len+=one_length;continue
                if len0>10000:
                    if start0<bigblock_start:           bigblock_start=start0
                    if end0>bigblock_end:               bigblock_end=end0
            break_mark='no'
            while break_mark=='no':
                break_mark=''
                for one_info in sorted_infos_pos:
                    one_start,one_end,one_length=one_info
                    if bigblock_start>one_end and bigblock_start-one_end<50000: bigblock_start=one_start; break_mark='no'
                    elif bigblock_end<one_start and one_start-bigblock_end<50000: bigblock_end=one_end; break_mark='no'
                
            core_VSat1_length=0
            for one_info in sorted_infos_pos:
                one_start,one_end,one_length=one_info
                if one_end<bigblock_start:      continue
                if one_start>bigblock_end:      continue
                core_VSat1_length+=one_length
            
            inter_array_num=0;inter_array_length=0;inter_array_pos=[];inter_array_str='NA'
            if VSat1_num>=2:
                iii=0
                while iii<VSat1_num-1:
                    start0,end0,len0=sorted_infos_pos[iii]
                    start1,end1,len1=sorted_infos_pos[iii+1]
                    if start0>=bigblock_start and start0<=bigblock_end and start1>=bigblock_start and start1<=bigblock_end:
                        if start1-end0>50000:inter_array_num+=1;inter_array_pos.append(f"{end0}-{start1}");inter_array_length+=(start1-end0+1)
                    iii+=1
            if len(inter_array_pos)>=1:
                inter_array_str='|'.join(inter_array_pos)
            largest_part_start=bigblock_start;largest_part_end=bigblock_end;bigblock_part_str=f'{largest_part_start}-{largest_part_end}'
            if inter_array_str!='NA':
                inter_arrays=inter_array_pos
                inter_arrays_len=len(inter_array_pos)
                bigblock_part_list=[];bigblock_part_str_list=[]
                kk=0
                while kk<inter_arrays_len:
                    
                    if kk==0:                           
                        A=bigblock_start
                        B=int(inter_arrays[kk].split('-')[0])
                        bigblock_part_list.append([A,B,B-A+1])
                        bigblock_part_str_list.append(f"{A}-{B}")
                    if kk+1<inter_arrays_len:           
                        A=int(inter_arrays[kk].split('-')[1])
                        B=int(inter_arrays[kk+1].split('-')[0])
                        bigblock_part_list.append([A,B,B-A+1])
                        bigblock_part_str_list.append(f"{A}-{B}")
                    if kk==inter_arrays_len-1:            
                        A=int(inter_arrays[kk].split('-')[1])
                        B=bigblock_end
                        bigblock_part_list.append([A,B,B-A+1]) 
                        bigblock_part_str_list.append(f"{A}-{B}")
                    kk+=1
                sorted_bigblock_part_list = sorted(bigblock_part_list, key=lambda x: x[2], reverse=True)    
                largest_part_start=sorted_bigblock_part_list[0][0]
                largest_part_end=sorted_bigblock_part_list[0][1]    
                bigblock_part_str='|'.join(bigblock_part_str_list)
                
            ff.write(f"{sample}\t{chromosome}\t{bigblock_start}\t{bigblock_end}\t{all_VSat1_length}\t{other_VSat1_len}\t{core_VSat1_length}\t{inter_array_num}\t{inter_array_length}\t{inter_array_str}\t{largest_part_start}\t{largest_part_end}\t{bigblock_part_str}\n")
if argv1=="step3.0_block_type_LIR":                  
    print('Based on LIR positions from moddotplot, filter LIR-related interarrays')
    dict_samplechr_list={}
    with open('./new_work_dir/chr2moddotplot3/4_VSat1out_down_sum/3_lines_filter','r') as f:
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            sample=         eachline_arr[0]
            chromosome=     eachline_arr[1]	
            start_x=        int(eachline_arr[3])
            start_y=        int(eachline_arr[4])
            sample_chr=sample+'___'+chromosome
            if sample_chr not in dict_samplechr_list:dict_samplechr_list[sample_chr]=[]
            dict_samplechr_list[sample_chr].append([start_x,start_y])
    with open('./samples_satellite/2_good_regions_interarray_LIR','w') as ff:
        ff.write(f"sample\tchromosome\tbigblock_start\tbigblock_end\tall_VSat1_length\tother_VSat1_len\tcore_VSat1_length\tinter_array_num\tinter_array_length\tinter_arrays\tlargest_part_start\tlargest_part_end\tbigblock_part_str\tLIRs_interarray\n")    
        with open('./samples_satellite/2_good_regions_interarray','r') as f:
            next(f)
            for line in f:
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                sample=         eachline_arr[0]
                chromosome=     eachline_arr[1]	
                sample_chr=sample+'___'+chromosome
                inter_array=eachline_arr[9]
                if inter_array=='NA':LIRs_interarray='NA'
                elif sample_chr not in dict_samplechr_list: LIRs_interarray='NA'
                else:
                    LIR_list=dict_samplechr_list[sample_chr]
                    inter_arrays=inter_array.split('|')
                    LIR_interarray_list=[]
                    for one_interarray in inter_arrays:
                        start,end=one_interarray.split('-')
                        start=int(start)
                        end=int(end)
                        mark=''
                        for one_LIR in LIR_list:
                            one_LIR_start,one_LIR_end=one_LIR
                            if start>one_LIR_start and end<one_LIR_end: mark='yes';break
                            if start<one_LIR_start and end>one_LIR_start: mark='yes';break                        
                            if start<one_LIR_end   and end>one_LIR_end: mark='yes';break    
                        if mark=="yes":
                            LIR_interarray_list.append(f"{start}-{end}")
                    if len(LIR_interarray_list)>0:
                        LIRs_interarray='|'.join(LIR_interarray_list)
                    else:        
                        LIRs_interarray='NA' 
                    ff.write(f"{eachline}\t{LIRs_interarray}\n")       
if argv1=="step3.0_block_type_LIR_strand": 
    dict_samplechr_pos_strand={}
    with open('./samples_satellite/2_good_regions','r') as f:
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=8:continue
            sample=eachline_arr[0]
            chromosome=eachline_arr[1].replace('region_','Chr')
            start=eachline_arr[3]
            end=eachline_arr[4]
            strand=eachline_arr[6]
            samplechr=sample+'___'+chromosome
            if samplechr not in dict_samplechr_pos_strand:      dict_samplechr_pos_strand[samplechr]=[]
            dict_samplechr_pos_strand[samplechr].append([start,end,strand])                
    with open('./samples_satellite/2_good_regions_interarray_LIR_strand','w') as f2:        
        f2.write(f"sample\tchromosome\tone_interarray\tleft_array_start\tleft_array_end\tleft_array_strand\tright_array_start\tright_array_end\tright_array_strand\n")
        with open('./samples_satellite/2_good_regions_interarray_LIR','r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                sample,chromosome=eachline_arr[:2]
                samplechr=sample+'___'+chromosome
                inter_arrays=eachline_arr[-1].split('|')
                maybe_array= dict_samplechr_pos_strand[samplechr]
                for one_interarray in inter_arrays:
                    if one_interarray=="NA":continue
                    start,end=one_interarray.split('-')
                    left_info=''
                    right_info=''
                    for one_maybe in maybe_array:
                        one_maybe_start,one_maybe_end,strand = one_maybe
                        if one_maybe_end==start: left_info=[one_maybe_start,one_maybe_end,strand]
                        if one_maybe_start==end:right_info=[one_maybe_start,one_maybe_end,strand]
                    if left_info=='' or right_info=='': print('error');continue
                    f2.write(f"{sample}\t{chromosome}\t{one_interarray}\t{left_info[0]}\t{left_info[1]}\t{left_info[2]}\t{right_info[0]}\t{right_info[1]}\t{right_info[2]}\n")

#step4，umap
if argv1=="stepall" or "step4" in argv1:
    if 'seed'=='seed':
        import random
        random_seed=7309 
        if "seed"  in args_dict: #print("Missing input fasta file");sys.exit()
            random_seed = int(args_dict["seed"])
        elif 'random_seed' in locals():
            print(f'umap manually set seed to {random_seed}')  
        else:
            print("Random seed")
            random_seed = random.randint(1, 10000)  # You can choose any integer as the random seed
            #random_seed=355
        seed_str=f'_seed{random_seed}'
    if argv1=="stepall" or argv1=="step4" or argv1=="step4.0":   
        if  os.path.exists(f'./samples_satellite/4_umap{seed_str}')==False:
            subprocess.run([f"mkdir ./samples_satellite/4_umap{seed_str}"], shell=True) 
        
         
        print('step4.0 ———— UMAP dimensionality reduction to 1D and 2D')
        print(f'umap seed is {seed_str}')
        import numpy as np
        import pandas as pd
        from ast import literal_eval
        import umap
        from umap import UMAP          
        umap_list=  [
                        [1,15,'./samples_satellite/3_samples/3_samples_monomer_vector',f'./samples_satellite/4_umap{seed_str}/monomer_vector_umap1v'],
                        [2,15,'./samples_satellite/3_samples/3_samples_monomer_vector',f'./samples_satellite/4_umap{seed_str}/monomer_vector_umap2v']
                    ]
        def run_umap(umap_list_one):
            print(umap_list_one)
            umap_component_num=     umap_list_one[0]
            n_neighbors=            umap_list_one[1]
            input_file=             umap_list_one[2]
            output_file=            umap_list_one[3]
            if os.path.exists(output_file)==True:
                print(f"UMAP dimensionality reduction completed: {output_file}") 
                return False
            ##
            df = pd.read_csv(input_file, sep='\t', header=None)
            def parse_vector(vector_str):
                return np.array(literal_eval(vector_str))
            for col in df.columns[1:]:
                df[col] = df[col].apply(parse_vector)
            parsed_df = df.iloc[:, 1:].values
            #print(parsed_df[0, :]) 
            features = np.array([np.concatenate(row) for row in parsed_df])
            
            # Set random seed
            #                                                                                  ###Set to 15 in umap_list       ###42,
            umap_model = UMAP(n_components=umap_component_num, n_neighbors=n_neighbors, random_state=random_seed)
    
            try:
                reduced_features = umap_model.fit_transform(features)
                print("Reduced features:")
                print(reduced_features)
                if umap_component_num==1:
                    result_df = pd.DataFrame(reduced_features, columns=['umap1v_x'])
                    result_df['id'] = df.iloc[:, 0]  # Assuming the first column is the name column
                    result_df = result_df[['id', 'umap1v_x']]  # Adjust column order
                else:
                    result_df = pd.DataFrame(reduced_features, columns=['umap2v_x', 'umap2v_y'])
                    result_df['id'] = df.iloc[:, 0]  # Assuming the first column is the name column
                    result_df = result_df[['id', 'umap2v_x', 'umap2v_y']]  # Adjust column order                
                result_df.to_csv(output_file ,sep='\t',  index=False)
                print(f"UMAP dimensionality reduction completed: {output_file}")        
            except Exception as e:
                print("An error occurred:", e)
          
    
        with Pool(processes=thread) as pool:
            pool.map(run_umap, umap_list)  
    if argv1=="step4.01":     
        print('BSCAN clustering, total 2.57 million points, each cluster requires more than 10,000')
        #import pandas as pd
        #import hdbscan
        #import matplotlib.pyplot as plt
        
        import pandas as pd
        from sklearn.cluster import DBSCAN
         
        # 1. Read tab-separated data file
        file_path = f'./samples_satellite/4_umap{seed_str}/monomer_vector_umap2v'  # Replace with your file path
        data = pd.read_csv(file_path, sep='\t')  # Specify separator as tab
         
        # 2. Extract UMAP coordinates (assuming column names are 'UMAP1' and 'UMAP2')
        umap_coords = data[['umap2v_x', 'umap2v_y']]
        print(umap_coords.head()) # Check data
        # 3. Use DBSCAN for clustering
        dbscan = DBSCAN(eps=0.5, min_samples=10000)  # Adjust parameters based on data distribution
        clusters = dbscan.fit_predict(umap_coords)
         
        # 4. Add clustering results to original data
        #data['Cluster'] = clusters
         
        # 5. Visualize clustering results (commented out)
        #plt.figure(figsize=(10, 8))
        #plt.scatter(umap_coords['UMAP1'], umap_coords['UMAP2'], c=clusters, cmap='Spectral', s=50)
        #plt.colorbar(label='Cluster')
        #plt.title('DBSCAN Clustering on UMAP')
        #plt.xlabel('UMAP1')
        #plt.ylabel('UMAP2')
        #plt.show()
         
        # 6. Save results to new file (optional)
        output_file = f'./samples_satellite/4_umap{seed_str}/monomer_vector_umap2v_dbscan'
        data.to_csv(output_file, sep='\t', index=False)  # Save as tab-separated file
        print(f"Clustering results saved to {output_file}")
    if argv1=="stepall" or argv1=="step4" or argv1=="step4.1":   
        print('step4.1 ———— Combine 1D and 2D UMAP dimensionality reduction results, 11s')     
        dict_circseq_umap1vx={} 
        with open(f'./samples_satellite/4_umap{seed_str}/monomer_vector_umap1v','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                dict_circseq_umap1vx[eachline_arr[0]]=eachline_arr[1]
        dict_circseq_umap2vxy={} 
        with open(f'./samples_satellite/4_umap{seed_str}/monomer_vector_umap2v','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                dict_circseq_umap2vxy[eachline_arr[0]]=[eachline_arr[1] ,eachline_arr[2] ]
        print('Converting')               
        with open(f'samples_satellite/4_umap{seed_str}/monomer_umap','w') as f:
            with open('samples_satellite/3_samples/3_samples_monomer','r') as f2:
                for line in f2.readlines():
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    circ_seq=eachline_arr[0]
                    if circ_seq=='circ_seq':
                        f.write('circ_seq\tin_seat_num\tcirc_len\tnum\tvariant_distance\tumap1v_x\tumap2v_x\tumap2v_y\n')
                    else:    
                        umap1vx=dict_circseq_umap1vx[circ_seq]      
                        umap2vxy=dict_circseq_umap2vxy[circ_seq]          
                        newline=eachline+'\t'+umap1vx+'\t'+umap2vxy[0]+'\t'+umap2vxy[1]+'\n'
                        f.write(newline)
    if argv1=="stepall" or argv1=="step4"  or argv1=="step4.2": 
        print('step4.2 ———— Load monomer length information etc. and UMAP data into each sample, time 324s')    
        
        if  os.path.exists(f'./samples_satellite/4_umap{seed_str}/monomer')==True:
            subprocess.run([f"rm -r ./samples_satellite/4_umap{seed_str}/monomer"], shell=True)  
        subprocess.run([f"mkdir ./samples_satellite/4_umap{seed_str}/monomer"], shell=True) 
        
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)
        print('Loading dictionary')
        #circ_seq        in_seat_num     circ_len        num     variant_distance2onehap variant_distance2allhaps
        dict_circseq_infos={}
        with open ("./samples_satellite/3_samples/3_samples_monomer_plus",'r') as f :
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                dict_circseq_infos[eachline_arr[0]]=[eachline_arr[1] ,eachline_arr[2],eachline_arr[3] ,eachline_arr[4],eachline_arr[5]]
        
        print('Loading dictionary')        
        dict_circseq_umap1vx={} 
        with open(f'./samples_satellite/4_umap{seed_str}/monomer_vector_umap1v','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                dict_circseq_umap1vx[eachline_arr[0]]=eachline_arr[1]
        dict_circseq_umap2vxy={} 
        with open(f'./samples_satellite/4_umap{seed_str}/monomer_vector_umap2v','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                dict_circseq_umap2vxy[eachline_arr[0]]=[eachline_arr[1] ,eachline_arr[2] ]
        print('Converting')    
        region_list.sort()
        region_list_len=len(region_list)
        print("region_list_len:"+str(region_list_len))
        i=0
        for one in region_list:
            i+=1
            print('Progress: '+str(i)+'/'+str(region_list_len),end='\r')
            input_file=     './samples_satellite/3_samples/monomer/'+one
            output_file=    f'./samples_satellite/4_umap{seed_str}/monomer/'+one
            with open(output_file,'w') as f:
                f.write('region_name\tcirc_serial\tchrserial_1\tchrserial_2\tstrand\tcirc_seq\tumap1v_x\tumap2v_x\tumap2v_y\tin_seat_num\tcirc_len\tnum\tvariant_distance2onehap\tvariant_distance2allhaps\n')
                k=0
                with open(input_file,'r') as f2:
                    for line in f2.readlines():
                        k+=1
                        if k==1:continue
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        circ_seq=   eachline_arr[5]
                        umap1v_x=   dict_circseq_umap1vx[circ_seq]
                        umap2v_xy=  dict_circseq_umap2vxy[circ_seq]
                        umap2v_x=   umap2v_xy[0]
                        umap2v_y=   umap2v_xy[1]
                        infos   =   dict_circseq_infos[circ_seq]
                        in_seat_num=                infos[0]
                        circ_len=                   infos[1]        
                        num=                        infos[2] 
                        variant_distance2onehap=    infos[3] 
                        variant_distance2allhaps=   infos[4] 
                        f.write(eachline+'\t'+umap1v_x+'\t'+umap2v_x+'\t'+umap2v_y+'\t'+in_seat_num+'\t'+circ_len+'\t'+num+'\t'+variant_distance2onehap+'\t'+variant_distance2allhaps+'\n')
        #####
    if argv1=="stepall" or argv1=="step4"  or argv1=="step4.3":
        print('Point plotting, 216s')
        #Plot: chromosome ———— point ———— count
        if  os.path.exists(f'./samples_satellite/4_umap{seed_str}/stat')==True:
            subprocess.run([f"rm -r ./samples_satellite/4_umap{seed_str}/stat"], shell=True)  
        subprocess.run([f"mkdir ./samples_satellite/4_umap{seed_str}/stat"], shell=True)     
        #Load all regions
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)    
        region_list.sort()
        region_list_len=len(region_list)
        print("region_list_len:"+str(region_list_len))
        #
    
        #Retrieve all regions
        print('Loading monomer data for all regions')
        dict_chr1circseq_num={}
        dict_circseq_infos={}
        i=0
        for one in region_list:
            i+=1
            #if i>10:break
            print('Progress: '+str(i)+'/'+str(region_list_len),end='\r')
            input_file=f'./samples_satellite/4_umap{seed_str}/monomer/'+one
            with open(input_file,'r') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    #('region_name\tcirc_serial\tchrserial_1\tchrserial_2\tstrand\tcirc_seq\tumap1v_x\tumap2v_x\tumap2v_y\tin_seat_num\tcirc_len\tnum\tvariant_distance2onehap\tvariant_distance2allhaps\n')
                    region_name=eachline_arr[0];    
                    if region_name=='region_name':continue
                    circ_seq=eachline_arr[5]
                    umap2v_x=eachline_arr[7]
                    umap2v_y=eachline_arr[8]
                    circ_len=eachline_arr[10]
                    variant_distance2allhaps=eachline_arr[13]
                    dict_circseq_infos[circ_seq]=[umap2v_x,umap2v_y,circ_len,variant_distance2allhaps]
                    chr1circseq=region_name+'_||_'+circ_seq
                    if chr1circseq not in dict_chr1circseq_num:
                        dict_chr1circseq_num[chr1circseq]=0
                    dict_chr1circseq_num[chr1circseq]+=1 
        #print(dict_chr1circseq_num)
        # Sort dictionary by value in descending order
        dict_chr1circseq_num_sorted = dict(sorted(dict_chr1circseq_num.items(), key=lambda item: item[1], reverse=True))
        with open (f'./samples_satellite/4_umap{seed_str}/stat/4.3-chr_circseq_num_info','w') as f:
            f.write(f"chr\tcircseq\tnum\tumap2v_x\tumap2v_y\tcirc_len\tvariant_distance2allhaps\n")
            for chr1circseq,num in dict_chr1circseq_num_sorted.items():
                chr1circseq_arr=chr1circseq.split('_||_')
                chromosome_name     =chr1circseq_arr[0]
                circseq             =chr1circseq_arr[1]
                infos               =dict_circseq_infos[circseq]
                f.write(f"{chromosome_name}\t{circseq}\t{num}\t{infos[0]}\t{infos[1]}\t{infos[2]}\t{infos[3]}\n")
    if argv1=="stepall" or argv1=="step4"  or argv1=="step4.3" or argv1=="step4.3p":
        print('Point plotting, 55s')            
        R_txt=f'''library(ggplot2)
    library(dplyr)
    #install.packages("stringr") 
    library("stringr")
    
    # Set working directory
    setwd('./')
    
    # Read data
    input_file1 <- read.table('4.3-chr_circseq_num_info', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
    
    input_file2 <- input_file1 %>% filter(circ_len>20 & num >10)
    input_file3 <- input_file2 %>%
      arrange(as.numeric(sub("region_", "", chr)))
    region_order <- unique(input_file3$chr)
    
    input_file3 <- input_file3 %>%
      mutate(chr_standard = str_replace(chr, "^region_", "Chr"))
      
    input_file3 <- input_file3 %>%
          mutate(
            category = case_when(
              #percent < 0.5 ~ "Mixed",
              chr == "Chr1" ~ "Chr1",
              chr == "Chr2" ~ "Chr2",
              chr == "Chr3" ~ "Chr3",
              chr == "Chr4" ~ "Chr4",
              chr == "Chr5" ~ "Chr5",
              chr == "Chr6" ~ "Chr6",
              chr == "Chr7" ~ "Chr7",
              chr == "Chr8" ~ "Chr8",
              chr == "Chr9" ~ "Chr9",
              chr == "Chr10" ~ "Chr10",
              chr == "Chr11" ~ "Chr11",
              chr == "Chr12" ~ "Chr12",
              chr == "Chr13" ~ "Chr13",
              chr == "Chr14" ~ "Chr14",
              chr == "Chr15" ~ "Chr15",
              chr == "Chr16" ~ "Chr16",
              chr == "Chr17" ~ "Chr17",
              chr == "Chr18" ~ "Chr18",
              chr == "Chr19" ~ "Chr19",
              TRUE ~ "Other"
            ),
          )

    region_order <- unique(input_file3$chr)
    chr_order <- unique(input_file3$chr_standard)
    input_file3$chr_standard <- factor(input_file3$chr_standard, levels = chr_order)
    
    all_num=sum(input_file3$num)
    # Use mutate and case_when to add new column shape
    input_file3 <- input_file3 %>%
      mutate(shape = case_when(
        circ_len == 28 ~ 14, # Circle
        circ_len == 23 ~ 15, # Triangle
        TRUE ~ 17          # Square (default for other cases)
      ))
    
    # Create plot
    p <- ggplot()
    p <- ggplot(data = input_file3, aes(x = umap2v_x, y = umap2v_y, size = num,color =chr_standard)) +
      geom_point(      aes(shape = factor(shape)),        alpha =0.3, stroke = 0.1)
    p <- p + scale_color_discrete(breaks = chr_order)
    p <- p + theme_classic()
    p <- p + coord_equal() 
    # Add title and axis labels
    p <- p + labs(title =paste("Umap (All =", all_num, ")", sep = ""), x = "X", y = "Y")
    p <- p + theme(
      #axis.title = element_blank(),  # Hide axis titles
      #axis.text = element_blank(),   # Hide axis tick labels
      #axis.ticks = element_blank(),  # Hide axis ticks
      #axis.line = element_blank()    # Hide axis lines
    )+scale_color_manual(name = "Variant Mode", 
                           values = c("Mixed" = "grey80", 
                                      "Chr1" = "red",  
                                      "Chr2" = "#8CC63F",  
                                      "Chr3" = "#CC0066",  
                                      "Chr4" = "#148EFF",  
                                      "Chr5" = "#B2DF8A",   
                                      "Chr6" = "#008000",   
                                      "Chr7" = "#FB9A99",   
                                      "Chr8" = "#00C08F",   
                                      "Chr9" = "#005FAF",  
                                      "Chr10" = "#FFCC00",  
                                      "Chr11" = "#00ccff",  
                                      "Chr12" = "#6A3D9A",  
                                      "Chr13" = "#FFFF00",  
                                      "Chr14" = "#A0522D",  
                                      "Chr15" = "#9966CC",  
                                      "Chr16" = "#ff4dff",  
                                      "Chr17" = "#ff6699", 
                                      "Chr18" = "#0073CF",  
                                      "Chr19" = "#FF8C00",  
                                      "Other" = "#CCEBC5")) 
    
    
    pdf("4.3-pic_p1_monomer_umap.pdf", width = 60 / 2.54, height = 60 / 2.54)
    print(p)
    dev.off()
    
    input_file2 <- input_file1 %>% filter(circ_len>10 & num >10)
    input_file3 <- input_file2 %>%
      arrange(as.numeric(sub("region_", "", chr)))
    region_order <- unique(input_file3$chr)
    
    input_file3 <- input_file3 %>%
      mutate(chr_standard = str_replace(chr, "^region_", "Chr"))
    
    region_order <- unique(input_file3$chr)
    chr_order <- unique(input_file3$chr_standard)
    input_file3$chr_standard <- factor(input_file3$chr_standard, levels = chr_order)
    
    all_num=sum(input_file3$num)
    # Use mutate and case_when to add new column shape
    input_file3 <- input_file3 %>%
      mutate(shape = case_when(
        circ_len == 28 ~ 14, # Circle
        circ_len == 23 ~ 15, # Triangle
        TRUE ~ 17          # Square (default for other cases)
      ))
    
    # Create plot
    p <- ggplot()
    p <- ggplot(data = input_file3, aes(x = umap2v_x, y = umap2v_y, size = num,color =chr_standard)) +
      geom_point(      aes(shape = factor(shape)),        alpha =0.3, stroke = 0.1)
    p <- p + scale_color_discrete(breaks = chr_order)
    p <- p + facet_wrap(~ chr_standard, ncol = 4)
    p <- p + theme_bw()
    p <- p + coord_equal() 
    # Add title and axis labels
    p <- p + labs(title =paste("Umap (All =", all_num, ")", sep = ""), x = "X", y = "Y")
    p <- p + theme(
      #axis.title = element_blank(),  # Hide axis titles
      #axis.text = element_blank(),   # Hide axis tick labels
      #axis.ticks = element_blank(),  # Hide axis ticks
      #axis.line = element_blank()    # Hide axis lines
    )+scale_color_manual(name = "Variant Mode", 
                           values = c("Mixed" = "grey80", 
                                      "Chr1" = "red",  
                                      "Chr2" = "#8CC63F",  
                                      "Chr3" = "#CC0066",  
                                      "Chr4" = "#148EFF",  
                                      "Chr5" = "#B2DF8A",   
                                      "Chr6" = "#008000",   
                                      "Chr7" = "#FB9A99",   
                                      "Chr8" = "#00C08F",   
                                      "Chr9" = "#005FAF",  
                                      "Chr10" = "#FFCC00",  
                                      "Chr11" = "#00ccff",  
                                      "Chr12" = "#6A3D9A",  
                                      "Chr13" = "#FFFF00",  
                                      "Chr14" = "#A0522D",  
                                      "Chr15" = "#9966CC",  
                                      "Chr16" = "#ff4dff",  
                                      "Chr17" = "#ff6699", 
                                      "Chr18" = "#0073CF",  
                                      "Chr19" = "#FF8C00",  
                                      "Other" = "#CCEBC5")) 
    pdf("4.3-pic_p1_monomer_umap2.pdf", width = 100 / 2.54, height = 100 / 2.54)
    print(p)
    dev.off()
    
    '''
        with open(f'./samples_satellite/4_umap{seed_str}/stat/4.3-pic_ggplot.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = f"./samples_satellite/4_umap{seed_str}/stat/"
        os.chdir(new_directory)
        subprocess.run(['Rscript 4.3-pic_ggplot.R'], shell=True)    
        os.chdir('../../../')    
    
    ##Change according to the actual situation.
    if argv1=="step4.6_riparia":
        print('Based on the step4.3 plot, extract specific parts for statistics')
        umap2v_x_region=[4.7,6.2]
        umap2v_y_region=[6.7,8]
        square_name=f"{umap2v_x_region[0]}-{umap2v_x_region[1]}_{umap2v_y_region[0]}-{umap2v_y_region[1]}"
        output_dir=f'./samples_satellite/4_umap{seed_str}/stat/{square_name}'
    
        if os.path.exists(output_dir)==True:subprocess.run([f"rm -r {output_dir}"], shell=True)      
        subprocess.run([f"mkdir {output_dir}"], shell=True)      
        subprocess.run([f"mkdir {output_dir}/raw"], shell=True)      
        ###########
        info_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                info_list.append([new_name,square_name])    
        ##########
        def run_step4_6(one_info):
            one_region,square_name=one_info
            square_name_arr=square_name.split('_')
            square_x_arr=square_name_arr[0].split('-')
            square_y_arr=square_name_arr[1].split('-')
            square_x1=float(square_x_arr[0])
            square_x2=float(square_x_arr[1])
            square_y1=float(square_y_arr[0])
            square_y2=float(square_y_arr[1])       
            input_file=f'./samples_satellite/4_umap{seed_str}/monomer/{one_region}'
            output_dir=f'./samples_satellite/4_umap{seed_str}/stat/{square_name}/'
            output_stat=f'./samples_satellite/4_umap{seed_str}/stat/{square_name}/4.6-stat'
            
            one_region_arr=one_region.split(':')
            one_region_sample=one_region_arr[0]
            one_region_chr=one_region_arr[1][:-1]  # Remove +/- sign
            one_region_pos=one_region_arr[2]
            target_monomer_num=0
            with open(input_file,'r') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if eachline_arr[7]=='umap2v_x':continue
                    umap2v_x=float(eachline_arr[7])
                    umap2v_y=float(eachline_arr[8])
                    if umap2v_x>square_x1 and umap2v_x<square_x2 and umap2v_y>square_y1 and umap2v_y<square_y2:
                        target_monomer_num+=1
                        #######   
                        output_sampleanno=f'./samples_satellite/4_umap{seed_str}/stat/{square_name}/4.6-sampleanno'
                        with open(output_sampleanno,'a') as f:
                            if one_region_sample in ['V050.hap1','V050.hap2','V048.hap1','V048.hap2','V051.hap2','V051.hap1']:
                                f.write(f"{line.strip()}\trip\n") 
                            elif one_region_sample in ['V038.hap1','V038.hap2','V049.hap1','V049.hap2']:
                                f.write(f"{line.strip()}\trup\n") 
                            elif one_region_sample in ['V036.hap2','V032.hap2','V053.hap1','V052.hap2','V043.hap1','V052.hap1','V053.hap2','V023.hap2','V033.hap1','V033.hap2','V023.hap1','V034.hap2','V034.hap1','V032.hap1','V043.hap2','V019.hap1','V020.hap2','V019.hap2','V040.hap2','V041.hap1','V040.hap1','V022.hap2','V041.hap2','V022.hap1','V020.hap1','V037.hap1','V037.hap2','V018.hap1','V018.hap2','V055.hap1','V055.hap2']:
                                f.write(f"{line.strip()}\trip_rup\n") 
                            else:
                                f.write(f"{line.strip()}\tother\n") 
                        ######
            if  target_monomer_num>0:     
                #print(target_monomer_num)
                with open(output_stat,'a') as f:
                    f.write(f"{one_region_sample}\t{one_region_chr}\t{one_region_pos}\t{target_monomer_num}\n") 
                subprocess.run([f"cp {input_file} {output_dir}/raw/{one_region}"], shell=True)          
                
                    
                    
        # Assign tasks to processes in the process pool    
        with Pool(processes=30) as pool:
             # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step4_6, info_list), start=1):
                # Process results here, e.g., store or print
                progress = (i / len(info_list)) * 100
                #print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))      
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()  
        
        dict_sample_num={}
        dict_chr_num={}
        with open(f'{output_dir}/4.6-stat','r') as f:        
            for line in f:
                eachline_arr=line.strip().split('\t')
                one_sample=eachline_arr[0]
                one_chr=eachline_arr[1]
                one_num=int(eachline_arr[3])
                if one_sample not in dict_sample_num:
                    dict_sample_num[one_sample]=0
                dict_sample_num[one_sample]+=one_num
                if one_chr not in dict_chr_num:
                    dict_chr_num[one_chr]=0                
                dict_chr_num[one_chr]+=one_num
        print('\nOrganizing results')        
        dict_sample_num_sorted = dict(sorted(dict_sample_num.items(), key=lambda item: item[1], reverse=True)) 
        dict_chr_num_sorted = dict(sorted(dict_chr_num.items(), key=lambda item: item[1], reverse=True))
        with open(f'{output_dir}/4.6-stat2_sample_num','w') as f:
            for one_sample,num in dict_sample_num_sorted.items():
                f.write(f'{one_sample}\t{num}\n')       
        with open(f'{output_dir}/4.6-stat2_chr_num','w') as f:     
            for onechr,num in dict_chr_num_sorted.items():
                f.write(f'{onechr}\t{num}\n')
        subprocess.run([f"sort -n -k 4 -r {output_dir}/4.6-stat >{output_dir}/4.6-stat_sorted "], shell=True)      
    
        dict_sampletpye_circseq_num={}
        dict_circseq_infos={}
        with open(f'./samples_satellite/4_umap{seed_str}/stat/{square_name}/4.6-sampleanno','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                circseq=eachline_arr[5]
                umap2v_x=eachline_arr[7]
                umap2v_y=eachline_arr[8]
                sample_type=eachline_arr[-1]
                dict_circseq_infos[circseq]=[umap2v_x,umap2v_y]
                if sample_type not in dict_sampletpye_circseq_num:
                    dict_sampletpye_circseq_num[sample_type]={}
                if circseq not in dict_sampletpye_circseq_num[sample_type]:
                    dict_sampletpye_circseq_num[sample_type][circseq]=0
                dict_sampletpye_circseq_num[sample_type][circseq]+=1
        with open(f'./samples_satellite/4_umap{seed_str}/stat/{square_name}/4.6-sampleanno_simple','w') as f:
            f.write(f"sampletype\tcircseq\tumap2v_x\tumap2v_y\tnum\n")
            for sampletype,circseqnum in dict_sampletpye_circseq_num.items():
                for circseq,num in circseqnum.items():
                    circseq_infos=dict_circseq_infos[circseq]
                    f.write(f"{sampletype}\t{circseq}\t{circseq_infos[0]}\t{circseq_infos[1]}\t{num}\n")
    if argv1=="step4.6p_riparia":        
        print('Point plotting')            
        R_txt=f'''library(ggplot2)
    library(dplyr)
    library("stringr")
    
    # Set working directory
    setwd('./')
    
    # Read data
    input_file1 <- read.table('4.6-sampleanno_simple', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    input_file2 <- input_file1 %>% filter(num >=1)
    
    all_num=sum(input_file2$num)
    
    # Create plot
    p <- ggplot()
    p <- ggplot(data = input_file2, aes(x = umap2v_x, y = umap2v_y, size = num,color =sampletype)) +
      geom_point(      aes(),        alpha =0.5, stroke = 0.1)+
      
      scale_color_manual(name = "", 
                           values = c(
                                      "rip" = "#26547c",  
                                      "rup" = "#48e5c2",  
                                      "rip_rup" = "#ef476f",  
                                      "other" = "grey")) 
    p <- p + theme_classic()
    #p <- p + coord_equal() 
    p <- p + theme(
      #axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      #axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    ) 
    # Add title and axis labels
    p <- p + labs(title =paste("Umap (All =", all_num, ")", sep = ""), x = "X", y = "Y")
    
    print(p)
    
    
    pdf("4.6-sampleanno_simple.pdf", width = 20 / 2.54, height = 20 / 2.54)
    print(p)
    dev.off()
    
    '''
        umap2v_x_region=[4.7,6.2]
        umap2v_y_region=[6.7,8]
        square_name=f"{umap2v_x_region[0]}-{umap2v_x_region[1]}_{umap2v_y_region[0]}-{umap2v_y_region[1]}"
        new_directory =f'./samples_satellite/4_umap{seed_str}/stat/{square_name}'
        os.chdir(new_directory)
        with open('./4.6-pic_ggplot.R','w',encoding='utf-8') as f:       f.write(R_txt)
        subprocess.run(['Rscript 4.6-pic_ggplot.R'], shell=True)    
        os.chdir('../../../../')         

    if argv1=="step4.7":
        with open('./samples_satellite/3_samples/seat_base_stat3_consensusseq','r') as f:
            seat_base_stat3_consensusseq=f.read().strip()
        with open('./samples_satellite/3_samples/seat_base_stat4_consensusseq23','r') as f:
            seat_base_stat4_consensusseq23=f.read().strip()        
        with open('./samples_satellite/3_samples/seat_base_stat4_consensusseq28','r') as f:
            seat_base_stat4_consensusseq28=f.read().strip()        
        seat_base_stat3_consensusseq_new='||'.join(seat_base_stat3_consensusseq)+'||'
        seat_base_stat4_consensusseq23_new='||'.join(seat_base_stat4_consensusseq23)+'||';  
        seat_base_stat4_consensusseq23_new=seat_base_stat4_consensusseq23_new[:33]+'|||||||||'+seat_base_stat4_consensusseq23_new[-37:]
        seat_base_stat4_consensusseq28_new='||'.join(seat_base_stat4_consensusseq28)+'||'
        ##
        dict_circseq_infos={}
        with open(f'./samples_satellite/4_umap{seed_str}/monomer_vector_umap2v','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                dict_circseq_infos[eachline_arr[0]]=[eachline_arr[1],eachline_arr[2]]
        seat_base_stat3_consensusseq_new_info       =dict_circseq_infos[seat_base_stat3_consensusseq_new]
        seat_base_stat4_consensusseq23_new_info     =dict_circseq_infos[seat_base_stat4_consensusseq23_new]
        seat_base_stat4_consensusseq28_new_info     =dict_circseq_infos[seat_base_stat4_consensusseq28_new]
        with open(f'./samples_satellite/4_umap{seed_str}/stat/4.7-core_site','w') as f:
            f.write(f"seat_base_stat3_consensusseq\t{seat_base_stat3_consensusseq_new}\t{seat_base_stat3_consensusseq_new_info[0]}\t{seat_base_stat3_consensusseq_new_info[1]}\n")
            f.write(f"seat_base_stat4_consensusseq23\t{seat_base_stat4_consensusseq23_new}\t{seat_base_stat4_consensusseq23_new_info[0]}\t{seat_base_stat4_consensusseq23_new_info[1]}\n")
            f.write(f"seat_base_stat4_consensusseq28\t{seat_base_stat4_consensusseq28_new}\t{seat_base_stat4_consensusseq28_new_info[0]}\t{seat_base_stat4_consensusseq28_new_info[1]}\n")
    
    if argv1=="step4_getpart_seq":
        if os.path.exists(f'./samples_satellite/4_umap{seed_str}/getpart')==False:
            subprocess.run([f"mkdir ./samples_satellite/4_umap{seed_str}/getpart"], shell=True) 
        def nearest_point_on_line(m, k, X, Y):
            """
            Calculate the foot of the perpendicular from point (X, Y) to line y = m*x + k.
            
            Parameters:
                m (float): Slope of the line.
                k (float): Intercept of the line.
                X (float): x-coordinate of the point.
                Y (float): y-coordinate of the point.
            
            Returns:
                tuple: Foot coordinates (x_foot, y_foot).
            """
            if m == 0:
                # Special case: horizontal line, foot is (X, k)
                return (X, k)
            else:
                # General case: non-zero slope
                x_foot = (X + m * (Y - k)) / (m**2 + 1)
                y_foot = m * x_foot + k
                return (x_foot, y_foot)
        
        ## Manually modify the current run parameters        
        x_min=4.7;x_max=6.2;y_min=6.7;y_max=8    ;      m=-0.8789; k=12.056    
        x_min=6.5;x_max=8.8;y_min=1;y_max=2.9;            m=0.625;k=-2.875
      
        print('Remove points within range') 
        result_list=[]
        with open(f'./samples_satellite/4_umap{seed_str}/monomer_vector_umap2v','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                seq=eachline_arr[0]
                if seq=='id':continue
                v2x=float(eachline_arr[1])
                v2y=float(eachline_arr[2])
                if v2x>x_min and v2x<x_max and v2y>y_min and v2y<y_max:
                    result_list.append([seq,v2x,v2y])
        with open(f'./samples_satellite/4_umap{seed_str}/getpart/{x_min}_{x_max}_{y_min}_{y_max}','w') as f:
            for one_info in result_list:
                f.write(f"{one_info[0]}\t{one_info[1]}\t{one_info[2]}\n")
                
        print('Remove irrelevant points not on the trend line')   
        with open(f'./samples_satellite/4_umap{seed_str}/getpart/{x_min}_{x_max}_{y_min}_{y_max}_filter','w') as f:
            for one_info in result_list:
                if abs(m*one_info[1] +k -one_info[2])<0.1:  # Calculate trend line using Excel
                    # Calculate foot of perpendicular
                    X=one_info[1]
                    Y=one_info[2]
                    x_foot, y_foot=nearest_point_on_line(m, k,X , Y)
                    #     
                    seq=one_info[0]
                    seq_arr=seq.split('|')
                    seq_arr_len=len(seq_arr)
                    ii=0;new_seq_arr=[]
                    while ii< seq_arr_len:
                        one_base=seq_arr[ii]
                        if ii%2==0:
                            new_seq_arr.append(one_base)
                        ii+=1    
                    new_seq='\t'.join(new_seq_arr)        
                    f.write(f"{seq}\t{X}\t{Y}\t{x_foot}\t{y_foot}\t{new_seq}\n")       
        print('Plotting')
        R_txt=r'''
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    
    # Get command line arguments
    args <- commandArgs(TRUE)
     
    # Check number of arguments
    if (length(args) < 1) {
      stop("Usage: Rscript script.R <input_file>")
    }
     
    # Input file path
    inputfile <- args[1]
     
    # Output file path (append '.pdf' to input file path)
    outputfile1 <- paste0(inputfile, "_real.pdf")
    outputfile2 <- paste0(inputfile, "_simple.pdf")
    
    # Set working directory
    setwd("./")
    # Read data
    df <- read.table(inputfile, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    
    
    # Assume df is your data frame, and column 4 is named "y_pos"
    # Extract column 4 as y_pos, and convert the letter columns from wide to long format
    df_long <- df %>%
      select(y_pos = V5, everything()) %>%  # Extract column 4 and name it y_pos, keep other columns
      # Add a sorting column, sort by y_pos in descending order
      arrange((y_pos)) %>%  # Sort by y_pos descending #desc
      mutate(y_order = row_number()) %>%  # Assign a sequence number to each y_pos group
      ##
      pivot_longer(cols = 6:ncol(df), names_to = "Position", values_to = "Base") %>%
      filter(!is.na(Base) & Base != "") %>%  # Remove rows where Base is NA or empty string
      mutate(Position = as.numeric(factor(Position, levels = unique(Position)))) %>%  # Convert Position to integer starting from 1
      #
      select(y_pos,y_order, Position,Base)  # Keep only y_pos, Position, Base columns
    
    
    
    # Define color mapping
    # Define colors
    color_A <- '#009e73'
    color_C <- '#56b4e9'
    color_G <- '#e69f00'
    color_T <- '#cc79a7'
    #color_map <- c("A" = "green", "C" = "blue", "G" = "orange", "T" = "#ff00ff")
    color_map <- c("A" = color_A, "C" = color_C, "G" = color_G, "T" = color_T)
    
    # Create plot (horizontal segments)
    p=ggplot(df_long, aes(x = Position, xend = Position + 1, y = y_pos, yend = y_pos)) +
      geom_segment(aes(color = Base), linewidth = 0.1) +  # Adjust linewidth to control segment thickness
      scale_color_manual(values = color_map) +
      labs(x = "Position", y = "Y Position (Column 4)", title = "Base Distribution (Horizontal Segments)") +
      theme_classic() 
      
    pdf(outputfile1, width = 20 / 2.54, height = 20 / 2.54)
    print(p)
    dev.off()
    
    # Create plot (horizontal segments)
    p=ggplot(df_long, aes(x = Position, xend = Position + 1, y = y_order, yend = y_order)) +
      geom_segment(aes(color = Base), linewidth = 0.1) +  # Adjust linewidth to control segment thickness
      scale_color_manual(values = color_map) +
      labs(x = "Position", y = "Y Position (Column 4)", title = "Base Distribution (Horizontal Segments)") +
      theme_classic() 
    pdf(outputfile2, width = 20 / 2.54, height = 20 / 2.54)
    print(p)
    dev.off()
    '''
        new_directory = f"./samples_satellite/4_umap{seed_str}/getpart/"
        os.chdir(new_directory)
        with open('base_plot.R','w',encoding='utf-8') as f:     f.write(R_txt)    
        subprocess.run([f'Rscript base_plot.R {x_min}_{x_max}_{y_min}_{y_max}_filter'], shell=True)    
        os.chdir('../../../../')   
    

            
###UMAP statistics: Grid calculation
if argv1=="stepall" or "step5" in argv1:
    if 'seed'=='seed':
        import random
        random_seed=7309 
        if "seed"  in args_dict: #print("Missing input fasta file");sys.exit()
            random_seed = int(args_dict["seed"])
        elif 'random_seed' in locals():
            print(f'umap manually set seed to {random_seed}')  
        else:
            print("Random seed")
            random_seed = random.randint(1, 10000)  # You can choose any integer as the random seed
            #random_seed=355
        seed_str=f'_seed{random_seed}'
    if argv1=="stepall" or argv1=="step5" or argv1=="step5.0":
        print('Calculate and count all units for each sample. Based on UMAP results, divide into 50 grids horizontally and vertically, each grid size 0.1, step size 25')
        if os.path.exists(f'./samples_satellite/5_umap_grid')==False:
            subprocess.run([f"mkdir ./samples_satellite/5_umap_grid"], shell=True)     
        
        
        import math 
        def umap2center_0_5(number):
            #return 0.2 * math.floor(number / 0.2)+0.1   # Round down to specified decimal places, e.g., 3.54 becomes square 3.4-3.6, center point 3.5
            return 0.1 * math.floor(number / 0.1)+0.05
        print('Round down to specified decimal places, 0.1 per grid cell, e.g., 3.49 becomes square 3.4-3.5, center point 3.45')
        #circ_seq        in_seat_num     circ_len        num     variant_distance        umap1v_x        umap2v_x        umap2v_y
        #A||C||T||C||G||C||A||C||G||G||A||T||T||C||T||A||C||C||T||T||T||T||C||C||C||G||G||T||    28      28      1562284 9       24.036808       14.522831       15.895154
    
    
        with open(f'./samples_satellite/5_umap_grid/monomer_umap_grid','w') as f2:
            f2.write('circ_seq\tin_seat_num\tcirc_len\tnum\tvariant_distance\tumap1v_x\tumap2v_x\tumap2v_y\tumap2v_x_grid\tumap2v_y_grid\n')
            i=0
            with open(f'./samples_satellite/4_umap{seed_str}/monomer_umap','r') as f:
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    i+=1
                    if i==1: continue
                    umap2v_x_grid=umap2center_0_5(float(eachline_arr[-2]))
                    umap2v_y_grid=umap2center_0_5(float(eachline_arr[-1]))
                    f2.write(f"{eachline}\t{umap2v_x_grid:.2f}\t{umap2v_y_grid:.2f}\n")     
                    #:.1f will be formatted without decimal places; if float, keep one decimal place.
    if argv1=="stepall" or argv1=="step5" or argv1=="step5.1":    
        if os.path.exists(f'./samples_satellite/5_umap_grid/monomer')==False:
            subprocess.run([f"mkdir ./samples_satellite/5_umap_grid/monomer"], shell=True)         
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)
        print('Loading dictionary')        
        dict_circseq_gridinfos={} 
        with open(f'./samples_satellite/5_umap_grid/monomer_umap_grid','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                dict_circseq_gridinfos[eachline_arr[0]]=[eachline_arr[-2],eachline_arr[-1]]
        
        print('Converting')    
        region_list.sort()
        region_list_len=len(region_list)
        print("region_list_len:"+str(region_list_len))
        i=0
        for one in region_list:
            i+=1
            print('Progress: '+str(i)+'/'+str(region_list_len),end='\r')
            input_file=     f'./samples_satellite/4_umap{seed_str}/monomer/'+one
            output_file=    './samples_satellite/5_umap_grid/monomer/'+one
            with open(output_file,'w') as f:
                f.write('region_name\tcirc_serial\tchrserial_1\tchrserial_2\tstrand\tcirc_seq\tumap1v_x\tumap2v_x\tumap2v_y\tin_seat_num\tcirc_len\tnum\tvariant_distance2onehap\tvariant_distance2allhaps\tumap2v_x_grid\tumap2v_y_grid\n')
                k=0
                with open(input_file,'r') as f2:
                    for line in f2.readlines():
                        k+=1
                        if k==1:continue
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        circ_seq=   eachline_arr[5]
                        gridinfos=dict_circseq_gridinfos[circ_seq]
                        umap2v_x_grid=gridinfos[0]
                        umap2v_y_grid=gridinfos[1]
                        f.write(eachline+'\t'+umap2v_x_grid+'\t'+umap2v_y_grid+'\n')
    if argv1=="stepall" or argv1=="step5" or argv1=="step5.2":       
        print('Statistical analysis of different intervals')
        if os.path.exists(f'./samples_satellite/5_umap_grid/monomer_stat')==False:
            subprocess.run([f"mkdir ./samples_satellite/5_umap_grid/monomer_stat"], shell=True) 
        subprocess.run([f"mkdir ./samples_satellite/5_umap_grid/monomer_stat/1-dot_density"], shell=True) 
        subprocess.run([f"mkdir ./samples_satellite/5_umap_grid/monomer_stat/2-circlen_density"], shell=True)  
        subprocess.run([f"mkdir ./samples_satellite/5_umap_grid/monomer_stat/3-variant_density"], shell=True)     
        subprocess.run([f"mkdir ./samples_satellite/5_umap_grid/monomer_stat/4-blocksize_density"], shell=True)
        subprocess.run([f"mkdir ./samples_satellite/5_umap_grid/monomer_stat/5-direction_density"], shell=True) 
        subprocess.run([f"mkdir ./samples_satellite/5_umap_grid/monomer_stat/6-chromosome_density"], shell=True) 
        subprocess.run([f"mkdir ./samples_satellite/5_umap_grid/monomer_stat/7-sample_density"], shell=True)
        subprocess.run([f"mkdir ./samples_satellite/5_umap_grid/monomer_stat/8-order_density"], shell=True)
        
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)
        def run_step5_2(one_region):
            input_file=f"./samples_satellite/5_umap_grid/monomer/{one_region}"
            output_file1=f"./samples_satellite/5_umap_grid/monomer_stat/1-dot_density/{one_region}"
            output_file2=f"./samples_satellite/5_umap_grid/monomer_stat/2-circlen_density/{one_region}"
            output_file3=f"./samples_satellite/5_umap_grid/monomer_stat/3-variant_density/{one_region}"
            output_file4=f"./samples_satellite/5_umap_grid/monomer_stat/4-blocksize_density/{one_region}"
            output_file5=f"./samples_satellite/5_umap_grid/monomer_stat/5-direction_density/{one_region}"
            output_file6=f"./samples_satellite/5_umap_grid/monomer_stat/6-chromosome_density/{one_region}"
            output_file7=f"./samples_satellite/5_umap_grid/monomer_stat/7-sample_density/{one_region}"
            output_file8=f"./samples_satellite/5_umap_grid/monomer_stat/8-order_density/{one_region}"
            one_region_arr=one_region.split(':')
            sample=     one_region_arr[0]
            chrname=    one_region_arr[1].replace('region_','Chr')[:-1]
            block_pos_arr= one_region_arr[2].split('-')
            block_size=int(block_pos_arr[1])-int(block_pos_arr[0])
            i=0;order=''
            dict_pos_infos={}
            with open(input_file,'r') as f:
                for line in f:
                    i+=1
                    if i==1:continue
                    eachline_arr=line.strip().split('\t')
                    circ_len=                   eachline_arr[-6]
                    variant_distance2allhaps=   eachline_arr[-3]
                    umap2v_x_grid=              eachline_arr[-2]
                    umap2v_y_grid=              eachline_arr[-1]
                    pos=umap2v_x_grid+'_'+umap2v_y_grid
                    if i==2:    
                        last_pos=pos
                        continue
                    if pos not in dict_pos_infos:
                        dict_pos_infos[pos]={}
                        dict_pos_infos[pos]['dot']=0
                        dict_pos_infos[pos]['circlen']=[]
                        dict_pos_infos[pos]['variant']=[]
                        dict_pos_infos[pos]['blocksize']=[]
                        dict_pos_infos[pos]['direction']=[]
                        dict_pos_infos[pos]['chromosome']=[]
                        dict_pos_infos[pos]['sample']=[]
                        dict_pos_infos[pos]['order']=[]
                    dict_pos_infos[pos]['dot']+=1
                    dict_pos_infos[pos]['circlen'].append(circ_len)
                    dict_pos_infos[pos]['variant'].append(variant_distance2allhaps)
                    dict_pos_infos[pos]['blocksize'].append(block_size)
                    dict_pos_infos[pos]['chromosome'].append(chrname)
                    dict_pos_infos[pos]['sample'].append(sample)
                    if i>3:
                        dict_pos_infos[last_pos]['direction'].append(pos)
                    last_pos=pos
                    #####
                    if order!='':
                        order+=1
                        dict_pos_infos[pos]['order'].append(order)      
                    if circ_len=='23':order=1
                    #####                    
            with open(output_file1,'w')  as f1:
                with open(output_file2,'w')  as f2:
                    with open(output_file3,'w')  as f3:
                        with open(output_file4,'w')  as f4:
                            with open(output_file5,'w')  as f5:
                                with open(output_file6,'w')  as f6:
                                    with open(output_file7,'w')  as f7:
                                        with open(output_file8,'w')  as f8:
                                            for pos,infos in dict_pos_infos.items():
                                                one_dot_num=            infos['dot']
                                                one_circlen_list=        infos['circlen']
                                                one_variant_list=        infos['variant']
                                                one_blocksize_list=      infos['blocksize']
                                                one_direction_list=      infos['direction']
                                                one_chromosome_list=     infos['chromosome']
                                                one_sample_list=         infos['sample']
                                                one_order_list=          infos['order']
                                                f1.write(f"{pos}\t{one_dot_num}\n")
                                                #
                                                for one_circlen in one_circlen_list:
                                                    f2.write(f"{pos}\t{one_circlen}\n")
                                                #
                                                for one_variant in one_variant_list:
                                                    f3.write(f"{pos}\t{one_variant}\n")
                                                #
                                                for one_blocksize in one_blocksize_list:
                                                    f4.write(f"{pos}\t{one_blocksize}\n")
                                                #
                                                for one_direction in one_direction_list:
                                                    f5.write(f"{pos}\t{one_direction}\n")  
                                                #
                                                for one_chromosome in one_chromosome_list:
                                                    f6.write(f"{pos}\t{one_chromosome}\n")    
                                                #
                                                for one_sample in one_sample_list:
                                                    f7.write(f"{pos}\t{one_sample}\n")   
                                                #
                                                for one_order in one_order_list:
                                                    f8.write(f"{pos}\t{one_order}\n")                                                   
    
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step5_2, region_list), start=1):
                # Process results here, e.g., store or print
                progress = (i / len(region_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()    
    if argv1=="stepall" or argv1=="step5" or argv1=="step5.2s":    
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)
        ### Find mode and percentage of the list
        from collections import Counter
        def calculate_mode_and_percentage(num_list):
            counts = Counter(num_list)
            max_count = max(counts.values())
            mode = next(num for num, count in counts.items() if count == max_count)
            num_list_num=len(num_list)
            mode_percentage = round((max_count / num_list_num),3)
            return mode, num_list_num, mode_percentage
        ##############################################    
        print('\nStatistics 1-dot_density')
        print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
        if os.path.exists(f'./samples_satellite/5_umap_grid/monomer_stat/1-stat_dot_density')==False:
            dict_pos_dotnum={}
            for one_region in region_list:
                input_file=f"./samples_satellite/5_umap_grid/monomer_stat/1-dot_density/{one_region}"
                with open(input_file,'r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        one_pos=    eachline_arr[0]
                        num=        int(eachline_arr[1])
                        if one_pos not in dict_pos_dotnum:   dict_pos_dotnum[one_pos]=0
                        dict_pos_dotnum[one_pos]+=num
            with open('./samples_satellite/5_umap_grid/monomer_stat/1-stat_dot_density','w') as f :
                f.write(f"pos\tpos_x\tpos_y\tdotnum\n")
                for one_pos,num in dict_pos_dotnum.items():
                    arr=one_pos.split('_')
                    f.write(f"{one_pos}\t{arr[0]}\t{arr[1]}\t{num}\n")
            
        print('\nStatistics 2-circlen_density')
        print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
        if os.path.exists(f'./samples_satellite/5_umap_grid/monomer_stat/2-stat_circlen_density')==False:
            dict_pos_circlenlist={}
            for one_region in region_list:
                input_file=f"./samples_satellite/5_umap_grid/monomer_stat/2-circlen_density/{one_region}"
                with open(input_file,'r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        one_pos=    eachline_arr[0]
                        circlen=    int(eachline_arr[1])
                        if one_pos not in dict_pos_circlenlist:   dict_pos_circlenlist[one_pos]=[]
                        dict_pos_circlenlist[one_pos].append(circlen)
            with open('./samples_satellite/5_umap_grid/monomer_stat/2-stat_circlen_density','w') as f :
                f.write(f"pos\tpos_x\tpos_y\tcirclen_mode\tsum_num\tpercent\n")    
                for one_pos,circlenlist in dict_pos_circlenlist.items():
                    arr=one_pos.split('_')
                    mode, sum_num, percent=calculate_mode_and_percentage(circlenlist)
                    f.write(f"{one_pos}\t{arr[0]}\t{arr[1]}\t{mode}\t{sum_num}\t{percent}\n")    
    
        print('\nStatistics 3-variant_density')
        print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
        if os.path.exists(f'./samples_satellite/5_umap_grid/monomer_stat/3-stat_variant_density')==False:
            dict_pos_variantlist={}
            for one_region in region_list:
                input_file=f"./samples_satellite/5_umap_grid/monomer_stat/3-variant_density/{one_region}"
                with open(input_file,'r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        one_pos=    eachline_arr[0]
                        variant=    int(eachline_arr[1])
                        if one_pos not in dict_pos_variantlist:   dict_pos_variantlist[one_pos]=[]
                        dict_pos_variantlist[one_pos].append(variant)
            with open('./samples_satellite/5_umap_grid/monomer_stat/3-stat_variant_density','w') as f :
                f.write(f"pos\tpos_x\tpos_y\tvariant_mode\tsum_num\tpercent\n")    
                for one_pos,variantlist in dict_pos_variantlist.items():
                    arr=one_pos.split('_')
                    mode, sum_num, percent=calculate_mode_and_percentage(variantlist)
                    f.write(f"{one_pos}\t{arr[0]}\t{arr[1]}\t{mode}\t{sum_num}\t{percent}\n")      
        
        print('\nStatistics 4-blocksize_density')
        print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
        if os.path.exists(f'./samples_satellite/5_umap_grid/monomer_stat/4-stat_blocksize_density')==False:
            dict_pos_blocksizelist={}
            for one_region in region_list:
                input_file=f"./samples_satellite/5_umap_grid/monomer_stat/4-blocksize_density/{one_region}"
                with open(input_file,'r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        one_pos=    eachline_arr[0]
                        blocksize=    int(eachline_arr[1])
                        blocksize_standard= math.log10(blocksize)
                        if blocksize_standard>6:blocksize_standard2='1M+'
                        elif blocksize_standard>5:blocksize_standard2='0.1M-1M'
                        elif  blocksize_standard>4:blocksize_standard2='10k-0.1M'
                        else:blocksize_standard2='<10k'
                        if one_pos not in dict_pos_blocksizelist:   dict_pos_blocksizelist[one_pos]=[]
                        dict_pos_blocksizelist[one_pos].append(blocksize_standard2)
            with open('./samples_satellite/5_umap_grid/monomer_stat/4-stat_blocksize_density','w') as f :
                f.write(f"pos\tpos_x\tpos_y\tblocksize_mode\tsum_num\tpercent\n")    
                for one_pos,blocksizelist in dict_pos_blocksizelist.items():
                    arr=one_pos.split('_')
                    mode, sum_num, percent=calculate_mode_and_percentage(blocksizelist)
                    f.write(f"{one_pos}\t{arr[0]}\t{arr[1]}\t{mode}\t{sum_num}\t{percent}\n")         
        
        print('\nStatistics 5-direction_density')
        print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
        if os.path.exists(f'./samples_satellite/5_umap_grid/monomer_stat/5-stat_direction_density')==False:
            dict_pos_directionlist={}
            for one_region in region_list:
                input_file=f"./samples_satellite/5_umap_grid/monomer_stat/5-direction_density/{one_region}"
                with open(input_file,'r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        one_pos=    eachline_arr[0]
                        direction=    eachline_arr[1]
                        if one_pos not in dict_pos_directionlist:   dict_pos_directionlist[one_pos]=[]
                        dict_pos_directionlist[one_pos].append(direction)
            with open('./samples_satellite/5_umap_grid/monomer_stat/5-stat_direction_density','w') as f :
                f.write(f"pos\tpos_x\tpos_y\tdirection_mode\tpos2_x\tpos2_y\tsum_num\tpercent\n")    
                for one_pos,directionlist in dict_pos_directionlist.items():
                    arr1=one_pos.split('_')
                    mode, sum_num, percent=calculate_mode_and_percentage(directionlist)
                    arr2=mode.split('_')
                    f.write(f"{one_pos}\t{arr1[0]}\t{arr1[1]}\t{mode}\t{arr2[0]}\t{arr2[1]}\t{sum_num}\t{percent}\n")        
                    
        print('\nStatistics 6-chromosome_density')
        print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
        if os.path.exists(f'./samples_satellite/5_umap_grid/monomer_stat/6-stat_chromosome_density')==False:
            dict_pos_chromosomelist={}
            for one_region in region_list:
                input_file=f"./samples_satellite/5_umap_grid/monomer_stat/6-chromosome_density/{one_region}"
                with open(input_file,'r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        one_pos=    eachline_arr[0]
                        chromosome=   eachline_arr[1]
                        if one_pos not in dict_pos_chromosomelist:   dict_pos_chromosomelist[one_pos]=[]
                        dict_pos_chromosomelist[one_pos].append(chromosome)
            with open('./samples_satellite/5_umap_grid/monomer_stat/6-stat_chromosome_density','w') as f :
                f.write(f"pos\tpos_x\tpos_y\tchromosome_mode\tsum_num\tpercent\n")    
                for one_pos,chromosomelist in dict_pos_chromosomelist.items():
                    arr=one_pos.split('_')
                    mode, sum_num, percent=calculate_mode_and_percentage(chromosomelist)
                    f.write(f"{one_pos}\t{arr[0]}\t{arr[1]}\t{mode}\t{sum_num}\t{percent}\n") 
                    
                    
        print('\nStatistics 7-sample_density')
        print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
        if os.path.exists(f'./samples_satellite/5_umap_grid/monomer_stat/7-stat_sample_density')==False:
            ##
            print('Loading sample information table ./samples_satellite/sample_info')
            dict_sample_type={}
            with open('./samples_satellite/sample_info') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if len(eachline_arr)!=6:continue
                    sample_name=eachline_arr[0]
                    sample_type=eachline_arr[3]
                    dict_sample_type[sample_name]=sample_type
            ##  
            print('Continuing')
            dict_pos_samplelist={}
            for one_region in region_list:
                input_file=f"./samples_satellite/5_umap_grid/monomer_stat/7-sample_density/{one_region}"
                with open(input_file,'r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        one_pos=    eachline_arr[0]
                        sample_name=   eachline_arr[1]
                        if sample_name not in dict_sample_type:continue
                        sample_type=   dict_sample_type[sample_name]
                        if one_pos not in dict_pos_samplelist:   dict_pos_samplelist[one_pos]=[]
                        dict_pos_samplelist[one_pos].append(sample_type)
            with open('./samples_satellite/5_umap_grid/monomer_stat/7-stat_sample_density','w') as f :
                f.write(f"pos\tpos_x\tpos_y\tsample_mode\tsum_num\tpercent\n")    
                for one_pos,samplelist in dict_pos_samplelist.items():
                    arr=one_pos.split('_')
                    mode, sum_num, percent=calculate_mode_and_percentage(samplelist)
                    f.write(f"{one_pos}\t{arr[0]}\t{arr[1]}\t{mode}\t{sum_num}\t{percent}\n")        
    
        print('\nStatistics 8-order_density')
        print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
        if os.path.exists(f'./samples_satellite/5_umap_grid/monomer_stat/8-stat_order_density')==False:
            dict_pos_orderlist={}
            for one_region in region_list:
                input_file=f"./samples_satellite/5_umap_grid/monomer_stat/8-order_density/{one_region}"
                with open(input_file,'r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        one_pos=    eachline_arr[0]
                        order=   eachline_arr[1]
                        if one_pos not in dict_pos_orderlist:   dict_pos_orderlist[one_pos]=[]
                        dict_pos_orderlist[one_pos].append(order)
            with open('./samples_satellite/5_umap_grid/monomer_stat/8-stat_order_density','w') as f :
                f.write(f"pos\tpos_x\tpos_y\torder_mode\tsum_num\tpercent\n")    
                for one_pos,orderlist in dict_pos_orderlist.items():
                    arr=one_pos.split('_')
                    mode, sum_num, percent=calculate_mode_and_percentage(orderlist)
                    f.write(f"{one_pos}\t{arr[0]}\t{arr[1]}\t{mode}\t{sum_num}\t{percent}\n")                 
    if argv1=="stepall" or argv1=="step5" or argv1=="step5.2p":             
        R_txt=r'''
    library(ggplot2)
    library(dplyr)
    
    setwd("./")
    
    print("1-stat_dot_density")
    {
      # Read data
      input_file1 <- read.table('1-stat_dot_density', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
      filtered_input_file1 <- input_file1 %>%  filter(dotnum > 100)
      
      filtered_input_file1 <- filtered_input_file1 %>%
        mutate(dotnum_category = cut(dotnum, 
                                     breaks = c(0,100, 1000, 5000, 10000,50000, 100000, 2000000),   ## 6 boundaries correspond to 5 intervals
                                     labels = c("0","1","2", "3", "4", "5", "6"), 
                                     right = FALSE))    #right = FALSE means interval [0, 10): includes 0 but excludes 10.
        
      # Use ggplot2 to draw scatter plot, dotnum for color mapping
      p=ggplot(filtered_input_file1, aes(x = pos_x, y = pos_y, color = dotnum_category)) +
        geom_point(size = 0.1,shape = 15) +
        scale_color_manual(values = c("0" = "#b3b3cc",
                                      "1" = "#666699",
                                      "2" = "#0099cc", 
                                      "3" = "#009999", 
                                      "4" = "#99cc00", 
                                      "5" = "#ffcc00", 
                                      "6" = "#ff5050"),name = "") +
        theme_bw()+ 
        coord_equal()+
        theme(
          panel.grid.major = element_blank(),  # Hide major grid lines
          panel.grid.minor = element_blank(),  # Hide minor grid lines
          #legend.position = "none",  # Hide legend
          legend.title = element_text(size = 14),  # Legend title font size
          legend.text = element_text(size = 12),   # Legend label font size
        )+ guides(color = guide_legend(override.aes = list(size = 5)))  # Ensure legend points are large enough to display as blocks
      
      pdf("pic_1-stat_dot_density.pdf", width = 14 / 2.54, height = 10 / 2.54)
      print(p)
      dev.off()
      print(p)
    }
    
    ###
    print("2-stat_circlen_density")
    {
      input_file2 <- read.table('2-stat_circlen_density', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
      filtered_input_file2 <- input_file2 %>%  filter(sum_num > 100)
      # Convert circlen_mode to factor and ensure factor level order is correct (optional)
      filtered_input_file2 <- filtered_input_file2 %>%
        mutate(circlen_mode_factor = factor(
          circlen_mode,
          levels = c(30, 29, 28, 27, 26, 25, 24, 23, 22, 21)  # Specify order as needed
        ))
      # Use ggplot2 to draw scatter plot, dotnum for color mapping
      # Use ggplot2 to draw scatter plot, define colors via scale_color_manual
      p <- ggplot(filtered_input_file2, aes(
        x = pos_x, 
        y = pos_y, 
        color = circlen_mode_factor  # Use factor variable
      )) +
        geom_point(size = 0.1, shape = 15) +
        scale_color_manual(
          values = c(
            "30" = "#0099FF",
            "29" = "#0099FF",
            "28" = "#003399",
            "27" = "#0099FF",
            "26" = "#0099FF",
            "25" = "#A5A559",
            "24" = "#A5A559",
            "23" = "#FF9900",
            "22" = "#cccc00",
            "21" = "#cccc00"
          ),
          name = "circlen_mode"  # Legend title
        ) +
        theme_bw() +
        coord_equal() +
        theme(
          panel.grid.major = element_blank(),  # Hide major grid lines
          panel.grid.minor = element_blank(),  # Hide minor grid lines
          #legend.position = "none",  # Hide legend
          legend.title = element_text(size = 10),  # Legend title font size
          legend.text = element_text(size = 10)   # Legend label font size
        ) +
        guides(color = guide_legend(override.aes = list(size = 3)))  # Ensure legend points are large enough to display as blocks
      
      pdf("pic_2-stat_circlen_density.pdf", width = 14 / 2.54, height = 10 / 2.54)
      print(p)
      dev.off()
      print(p)
    }
    
    ###3 Variation plot, binned legend
    print("3-stat_variant_density")
    {
      input_file3 <- read.table('3-stat_variant_density', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
      filtered_input_file3 <- input_file3 %>%  filter(sum_num > 100 & variant_mode<25)
      # Filter data
      filtered_input_file3 <- input_file3 %>% 
        filter(sum_num > 100 & variant_mode < 25) %>%
        mutate(variant_mode_bin = cut(variant_mode, 
                                      breaks = c(-Inf, 5, 10, 15, Inf), 
                                      labels = c("0-5", "5-10", "10-15","15+"), 
                                      right = FALSE))
      
      # Use ggplot2 to draw scatter plot
      p <- ggplot(filtered_input_file3, aes(x = pos_x, y = pos_y, color = variant_mode_bin)) +
        geom_point(size = 0.1, shape = 15) +
        scale_color_manual(values = c("0-5" = "#00cc99", "5-10" = "#006666", "10-15" = "#999966", "15+" = "#993333"), 
                           name = "") +
        theme_bw() + 
        coord_equal() +
        theme(
          panel.grid.major = element_blank(),  # Hide major grid lines
          panel.grid.minor = element_blank(),  # Hide minor grid lines
          #legend.position = "none",  # Hide legend
          legend.title = element_text(size = 10),  # Legend title font size
          legend.text = element_text(size = 10)   # Legend label font size
        )+guides(color = guide_legend(override.aes = list(size = 3)))  # Ensure legend points are large enough to display as blocks
      
      # Save plot to PDF
      pdf("pic_3-stat_variant_density.pdf", width = 14 / 2.54, height = 10 / 2.54)
      print(p)
      dev.off()
      print(p)
    }
    
    
    
    ### Binned legend blocksize
    print("4-stat_blocksize_density")
    {
      # Read data
      input_file4 <- read.table('4-stat_blocksize_density', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
      filtered_input_file4 <- input_file4 %>% filter(sum_num > 100)
      
      # Create categorical variable and color mapping
      filtered_input_file4 <- filtered_input_file4 %>%
        mutate(
          # Create categorical column
          category = case_when(
            percent < 0.5 ~ "Low Percent",
            blocksize_mode == "<10k" ~ "<10k",
            blocksize_mode == "10k-0.1M" ~ "10k-0.1M",
            blocksize_mode == "0.1M-1M" ~ "0.1M-1M",
            blocksize_mode == "1M+" ~ "1M+",
            TRUE ~ "Other"
          ),
        )
      
      # Use ggplot2 to draw scatter plot
      p <- ggplot(filtered_input_file4, aes(x = pos_x, y = pos_y, color = category)) +  # Use category for aesthetic mapping
        geom_point(size = 0.1, shape = 15) +
        theme_bw() +
        coord_equal() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #legend.position = "none",  # Hide legend
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)
        ) +
        scale_color_manual(name = "Variant Mode", 
                           values = c("Low Percent" = "grey80", 
                                      "<10k" = "black", 
                                      "10k-0.1M" = "blue", 
                                      "0.1M-1M" = "#ffb74d", 
                                      "1M+" = "red", 
                                      "Other" = "grey")) +
        guides(color = guide_legend(override.aes = list(size = 3), title.position = "top", title.hjust = 0.5))
      
      # Save as PDF
      pdf("pic_4-stat_blocksize_density2.pdf", width = 14 / 2.54, height = 10 / 2.54)
      print(p)
      dev.off()
      
      # Print plot (optional)
      print(p)
      
    }

    
    
    ####
    print("5-stat_direction_density")
    {
      library(grid)  # Provides arrow() function
      
      # Read data
      input_file5 <- read.table('5-stat_direction_density', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
      
      # Filter data
      filtered_input_file5 <- input_file5 %>%
        filter(sum_num > 10000 & percent > 0.2)
      
      # Add a new column to indicate the direction color of the line segment
      filtered_input_file5 <- filtered_input_file5 %>%
        mutate(
          direction_color = case_when(
            pos2_x > pos_x & pos2_y > pos_y ~ "red",    # From bottom-left to top-right
            pos2_x > pos_x & pos2_y < pos_y ~ "#000066",   # From top-left to bottom-right
            pos2_x < pos_x & pos2_y < pos_y ~ "#339966",  # From top-right to bottom-left (assuming y decreases when x decreases)
            pos2_x < pos_x & pos2_y > pos_y ~ "orange"  # From bottom-right to top-left
          )
        )
      
      # Draw line segment plot with arrows, set color based on direction
      p <- ggplot(filtered_input_file5) +
        geom_segment(aes(x = pos_x, y = pos_y, xend = pos2_x, yend = pos2_y, color = direction_color), 
                     size = 0.1, 
                     arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +  # Add arrow
        scale_color_identity() +  # Use identity scale to directly use our defined colors
        theme_bw() + 
        coord_equal() +
        theme(
          panel.grid.major = element_blank(),  # Hide major grid lines
          panel.grid.minor = element_blank()#  # Hide minor grid lines
          #legend.position = "none",  # Hide legend
        )
      
      # Save to PDF file
      pdf("pic_5-stat_direction_density.pdf", width = 14 / 2.54, height = 10 / 2.54)
      print(p)
      dev.off()
      
      # Print plot (optional)
      print(p)
    }
    
    
    
    
    ### Binned legend chromosome
    print("6-stat_chromosome_density")
    {
        input_file6 <- read.table('6-stat_chromosome_density', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
        filtered_input_file6 <- input_file6 %>%  filter(sum_num > 100)
        # Create categorical column
        filtered_input_file6 <- filtered_input_file6 %>%
          mutate(
            category = case_when(
              percent < 0.5 ~ "Mixed",
              chromosome_mode == "Chr1" ~ "Chr1",
              chromosome_mode == "Chr2" ~ "Chr2",
              chromosome_mode == "Chr3" ~ "Chr3",
              chromosome_mode == "Chr4" ~ "Chr4",
              chromosome_mode == "Chr5" ~ "Chr5",
              chromosome_mode == "Chr6" ~ "Chr6",
              chromosome_mode == "Chr7" ~ "Chr7",
              chromosome_mode == "Chr8" ~ "Chr8",
              chromosome_mode == "Chr9" ~ "Chr9",
              chromosome_mode == "Chr10" ~ "Chr10",
              chromosome_mode == "Chr11" ~ "Chr11",
              chromosome_mode == "Chr12" ~ "Chr12",
              chromosome_mode == "Chr13" ~ "Chr13",
              chromosome_mode == "Chr14" ~ "Chr14",
              chromosome_mode == "Chr15" ~ "Chr15",
              chromosome_mode == "Chr16" ~ "Chr16",
              chromosome_mode == "Chr17" ~ "Chr17",
              chromosome_mode == "Chr18" ~ "Chr18",
              chromosome_mode == "Chr19" ~ "Chr19",
              TRUE ~ "Other"
            ),
          )
        
        # Use ggplot2 to draw scatter plot
        p <- ggplot(filtered_input_file6, aes(x = pos_x, y = pos_y, color = category)) +  # Use category for aesthetic mapping
          geom_point(size = 0.1, shape = 15) +
          theme_bw() +
          coord_equal() +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            #legend.position = "none",  # Hide legend
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12)
          ) +
          scale_color_manual(name = "Variant Mode", 
                           values = c("Mixed" = "grey80", 
                                      "Chr1" = "red",  
                                      "Chr2" = "#8CC63F",  
                                      "Chr3" = "#CC0066",  
                                      "Chr4" = "#148EFF",  
                                      "Chr5" = "#B2DF8A",   
                                      "Chr6" = "#008000",   
                                      "Chr7" = "#FB9A99",   
                                      "Chr8" = "#00C08F",   
                                      "Chr9" = "#005FAF",  
                                      "Chr10" = "#FFCC00",  
                                      "Chr11" = "#00ccff",  
                                      "Chr12" = "#6A3D9A",  
                                      "Chr13" = "#FFFF00",  
                                      "Chr14" = "#A0522D",  
                                      "Chr15" = "#9966CC",  
                                      "Chr16" = "#ff4dff",  
                                      "Chr17" = "#ff6699", 
                                      "Chr18" = "#0073CF",  
                                      "Chr19" = "#FF8C00",  
                                      "Other" = "#CCEBC5")) +
          guides(color = guide_legend(override.aes = list(size = 1), title.position = "top", title.hjust = 0.5,keywidth  = unit(0.1, "cm"),keyheight = unit(0.1, "cm")))
     
        # Save as PDF
        pdf("pic_6-stat_chromosome_density.pdf", width = 14 / 2.54, height = 10 / 2.54)
        print(p)
        dev.off()
        print(p)
    }
    
    ##
    print('7-stat_sample_density')
    {
      input_file7 <- read.table('7-stat_sample_density', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
      filtered_input_file7 <- input_file7 %>%  filter(sum_num > 100)
      filtered_input_file7 <- filtered_input_file7 %>%
        mutate(
          category = case_when(
            percent < 0.5 ~ "Mixed",
            sample_mode == "Eurasian" ~ "Eurasian",
            sample_mode == "East_Asia" ~ "East_Asia",
            sample_mode == "America" ~ "America",
            sample_mode == "Routundifolia" ~ "Routundifolia",
            TRUE ~ "Other"
          ),
        )
      # Use ggplot2 to draw scatter plot, dotnum for color mapping
      # Use ggplot2 to draw scatter plot, define colors via scale_color_manual
      p <- ggplot(filtered_input_file7, aes(x = pos_x, y = pos_y, color = category)) +
        geom_point(size = 0.1, shape = 15) +
        scale_color_manual(
          values = c(
            'Mixed'="grey80",
            "Eurasian" = "#ff3399",
            "East_Asia" = "#0066ff",
            "America" = "#33cc33",
            "Routundifolia" = "black"
          ),
          name = "" 
        ) +
        theme_bw() +
        coord_equal() +
        theme(
          panel.grid.major = element_blank(),  # Hide major grid lines
          panel.grid.minor = element_blank(),  # Hide minor grid lines
          #legend.position = "none",  # Hide legend
          legend.title = element_text(size = 10),  # Legend title font size
          legend.text = element_text(size = 10)   # Legend label font size
        ) +
        guides(color = guide_legend(override.aes = list(size = 3)))  # Ensure legend points are large enough to display as blocks
      
      pdf("pic_7-stat_sample_density.pdf", width = 14 / 2.54, height = 10 / 2.54)
      print(p)
      dev.off()
      print(p)
    }
    
    
    ### Binned legend order
    print("8-stat_order_density")
    {
        input_file8 <- read.table('8-stat_order_density', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
        filtered_input_file8 <- input_file8 %>%  filter(sum_num > 100)
        # Create categorical column
        filtered_input_file8 <- filtered_input_file8 %>%
          mutate(
            category = case_when(
              percent < 0.5 ~ "Mixed",
              order_mode == "1" ~ "1",
              order_mode == "2" ~ "2",
              order_mode == "3" ~ "3",
              order_mode == "4" ~ "4",
              order_mode == "5" ~ "5",
              order_mode == "6" ~ "6",
              order_mode == "7" ~ "7",
              order_mode == "8" ~ "8",
              order_mode == "9" ~ "9",
              order_mode == "10" ~ "10",
              order_mode == "11" ~ "11",
              order_mode == "12" ~ "12",
              TRUE ~ "Other"
            ),
          )
        
        # Use ggplot2 to draw scatter plot
        p <- ggplot(filtered_input_file8, aes(x = pos_x, y = pos_y, color = category)) +  # Use category for aesthetic mapping
          geom_point(size = 0.1, shape = 15) +
          theme_bw() +
          coord_equal() +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            #legend.position = "none",  # Hide legend
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12)
          ) +
          scale_color_manual(name = "Variant Mode", 
                           values = c("Mixed" = "grey80", 
                                      "1" = "red",  
                                      "2" = "#8CC63F",  
                                      "3" = "#CC0066",  
                                      "4" = "#148EFF",  
                                      "5" = "#B2DF8A",   
                                      "6" = "#008000",   
                                      "7" = "#FB9A99",   
                                      "8" = "#00C08F",   
                                      "9" = "#005FAF",  
                                      "10" = "#FFCC00",  
                                      "11" = "#00ccff",  
                                      "12" = "#6A3D9A",  
                                      "Other" = "CCEBC5")) +
          guides(color = guide_legend(override.aes = list(size = 1), title.position = "top", title.hjust = 0.5,keywidth  = unit(0.1, "cm"),keyheight = unit(0.1, "cm")))
     
        # Save as PDF
        pdf("pic_8-stat_order_density.pdf", width = 14 / 2.54, height = 10 / 2.54)
        print(p)
        dev.off()
        print(p)
    }
    '''
        new_directory = "./samples_satellite/5_umap_grid/monomer_stat/"
        os.chdir(new_directory)
        with open('./pic_ggplot_step5_2.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        subprocess.run(['Rscript pic_ggplot_step5_2.R'], shell=True)    
        os.chdir('../../../../')

#subunit to mononmer
if (argv1=="stepall" or "step9" in argv1) and subunit_mark=='yes':
    if 'seed'=='seed':
        import random
        random_seed=7309 
        if "seed"  in args_dict: #print("Missing input fasta file");sys.exit()
            random_seed = int(args_dict["seed"])
        elif 'random_seed' in locals():
            print(f'umap manually set seed to {random_seed}')  
        else:
            print("Random seed")
            random_seed = random.randint(1, 10000)  # You can choose any integer as the random seed
            #random_seed=355
        seed_str=f'_seed{random_seed}'
    ##9.0a uses 23 as the marker, 23------23
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.0a" or  argv1=="step9.0a_s0":
        print('step9.0 ———— Determine preliminary HOR based on the appearance of 21-25bp monomers')
        if  os.path.exists('./samples_satellite/9_circlenarr')==False:
            subprocess.run(["mkdir ./samples_satellite/9_circlenarr"], shell=True)
        if  os.path.exists('./samples_satellite/9_circlenarr/a_monomer'): 
            subprocess.run(["rm -r ./samples_satellite/9_circlenarr/a_monomer"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/9_circlenarr/a_monomer"], shell=True)  
    
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)   
        print('Main task started\n')     
        def run_step9_1(one_region):    
            input_file=f"./samples_satellite/4_umap{seed_str}/monomer/{one_region}"
            output=f"./samples_satellite/9_circlenarr/a_monomer/{one_region}"
            i=0
            start_mark=''
            with open(input_file,'r') as f:
                with open(output,'w') as f2:
                    for line in f:
                        i+=1
                        if i==1:
                            f2.write(f"region_name\tcirclen_arr\tstart\tend\tstrand\tseq_arr\tumap2v_arr\tdistance_arr\tarr_serial\tlength\n")
                            continue
                        eachline_arr=line.strip().split('\t')
                        circ_len=int(eachline_arr[10])
                        circ_seq=eachline_arr[5]
                        umap2v=f"{eachline_arr[7]}_{eachline_arr[8]}"
                        distance=eachline_arr[13]
                        end=eachline_arr[3]
                        if circ_len>=21 and circ_len<=25: 
                            if start_mark=='yes' and arr_serial>0:
                                start_mark=''
                                arr_serial+=1
                                circlen_arr.append(str(circ_len))
                                seq_arr.append(circ_seq)
                                umap2v_arr.append(umap2v)
                                distance_arr.append(distance)
                                circlen_arr_str='-'.join(circlen_arr)
                                seq_arr_str='_'.join(seq_arr)
                                umap2v_arr_str='_'.join(umap2v_arr)
                                distance_arr_str='_'.join(distance_arr)
                                region_name=eachline_arr[0]
                                strand=eachline_arr[4]
                                length=abs(int(end)-int(start))+1
                                f2.write(f"{region_name}\t{circlen_arr_str}\t{start}\t{end}\t{strand}\t{seq_arr_str}\t{umap2v_arr_str}\t{distance_arr_str}\t{arr_serial}\t{length}\n")
                            start_mark='yes'
                            arr_serial=0
                            circlen_arr=[]
                            seq_arr=[]
                            umap2v_arr=[]
                            distance_arr=[]
                            continue
                        elif start_mark=='yes':
                            arr_serial+=1
                            if arr_serial==1:start=eachline_arr[2]
                            circlen_arr.append(str(circ_len))
                            seq_arr.append(circ_seq)
                            umap2v_arr.append(umap2v)
                            distance_arr.append(distance)
        
        with Pool(processes=thread) as pool:
            for i, result in enumerate(pool.imap(run_step9_1, region_list), start=1):
                progress = (i / len(region_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush() 
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.0a" or  argv1=="step9.0a_s1": 
        print('step9.0s1 ———— Calculate the length of all regions')
        if  os.path.exists('./samples_satellite/9_circlenarr/a_stat')==False:
            subprocess.run(["mkdir ./samples_satellite/9_circlenarr/a_stat"], shell=True)    
        print('Loading region name list')
        region_all_len=0
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_all_len+= abs(int(eachline_arr[3])-int(eachline_arr[4]))+1
        with open('./samples_satellite/9_circlenarr/a_stat/all_region_length','w') as f:        
            f.write(str(region_all_len))
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.0a" or  argv1=="step9.0a_s2":        
        print('step9.0 ———— Determine preliminary HOR based on 21-25bp monomers, statistics')
        print('Loading all_region_length')
        with open('./samples_satellite/9_circlenarr/a_stat/all_region_length','r') as f:        
            all_region_length=int(f.read().strip())
        print('Loading region name list')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)  
        region_list_len=len(region_list)        
        print('Loading each region')
        dict_circlenarr_num={}
        i=0
        for one_region in region_list:
            i+=1
            print(f'Progress: {i}/{region_list_len}',end='\r')
            input_file=f"./samples_satellite/9_circlenarr/a_monomer/{one_region}"
            k=0
            with open(input_file,'r') as f:
                for line in f:
                    k+=1
                    if k==1:continue
                    eachline_arr=line.strip().split('\t')
                    circlen_arr_str=    eachline_arr[1]
                    if circlen_arr_str not in dict_circlenarr_num :dict_circlenarr_num[circlen_arr_str]=0
                    dict_circlenarr_num[circlen_arr_str]+=1
        sorted_dict_circlenarr_num = dict(sorted(dict_circlenarr_num.items(), key=lambda item:item[1],reverse=True))
    
        with open('./samples_satellite/9_circlenarr/a_stat/0-circlenarr_num','w') as f :
            f.write(f"circlenarr\tone_parts\tone_length\tnum\tall_length\tpercent\n")    
            for circlenarr,num in sorted_dict_circlenarr_num.items():
                one_length= sum(int(num) for num in circlenarr.split('-'))
                all_length=one_length*num
                percent=round(all_length/all_region_length,10)
                one_parts=len(circlenarr.split('-'))
                f.write(f"{circlenarr}\t{one_parts}\t{one_length}\t{num}\t{all_length}\t{percent}\n")          
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.0a" or  argv1=="step9.0a_s3":
        circlenarr_list=[]
        circlenarr_list.append([28,23])
        circlenarr_list.append([28,28,23])
        circlenarr_list.append([28,28,28,23])
        circlenarr_list.append([28,28,28,28,23])
        circlenarr_list.append([28,28,28,28,28,23])
        circlenarr_list.append([28,28,28,28,28,28,23])
        circlenarr_list.append([28,28,28,28,28,28,28,23])
        circlenarr_list.append([28,28,28,28,28,28,28,28,23])
        circlenarr_list.append([28,28,28,28,28,28,28,28,28,23])
        circlenarr_list.append([28,28,28,28,28,28,28,28,28,28,23])
        circlenarr_list_num=len(circlenarr_list)
        dict_k_info={}
        i=0
        with open('./samples_satellite/9_circlenarr/a_stat/0-circlenarr_num2other','w') as f2:    
            with open('./samples_satellite/9_circlenarr/a_stat/0-circlenarr_num','r') as f:
                for line in f:
                    i+=1
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    if len(eachline_arr)!=6:continue
                    circlenarr=         eachline_arr[0]
                    if circlenarr=='circlenarr':continue
                    circlenarr_arr=         circlenarr.split('-')
                    circlenarr_arr_len=     len(circlenarr_arr)
                    k=0;good_state=''
                    while k<circlenarr_list_num:
                        one_circlenarr=circlenarr_list[k]
                        if len(one_circlenarr)==circlenarr_arr_len:
                            kk=0;good_mark=0
                            while kk<circlenarr_arr_len:
                                delta=abs(one_circlenarr[kk]-int(circlenarr_arr[kk]))
                                kk+=1
                                if delta<=2:good_mark+=1
                            if good_mark==circlenarr_arr_len:
                                if k not in dict_k_info:dict_k_info[k]=[]
                                dict_k_info[k].append(eachline)
                                good_state='yes'
                        if good_state=='yes':break
                        k+=1    
                    if good_state=='':
                        f2.write(f"{eachline}\n")
            f2.write(f"\n")            
        k=0
        with open('./samples_satellite/9_circlenarr/a_stat/0-circlenarr_num2','w') as f:    
            while k<circlenarr_list_num:
                headstr="-".join([str(num) for num in circlenarr_list[k]])
                f.write(f"\n##################{headstr}\n")
                for oneinfo in dict_k_info[k]:
                    f.write(f"{oneinfo}\n")
                k+=1
            f.write(f"\n")      
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.0a" or  argv1=="step9.0a_s4":
        start_mark=''
        with open('./samples_satellite/9_circlenarr/a_stat/0-circlenarr_num2s','w') as f2: 
            with open('./samples_satellite/9_circlenarr/a_stat/0-circlenarr_num2','r') as f: 
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    if len(eachline)==0 and start_mark=='yes':
                        f2.write(f'#\t#\t#\t{num_acc}\t{length_acc}\t{percent_acc}\n')
                        num_acc=0;length_acc=0;percent_acc=0;start_mark=''            
                    if len(eachline_arr)==1: 
                        f2.write(eachline+'\n')  
                        num_acc=0;length_acc=0;percent_acc=0;start_mark='yes'
                    elif start_mark=='yes': 
                        num_acc+=int(eachline_arr[3])
                        length_acc+=int(eachline_arr[4])
                        percent_acc+=float(eachline_arr[5])
    ##9.0b uses 28 as the marker, 28------28            
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.0b" or  argv1=="step9.0b_s0":    
        print('step9.0 ———— Determine preliminary HOR based on the appearance of 21-25bp monomers')
        if  os.path.exists('./samples_satellite/9_circlenarr')==False:
            subprocess.run(["mkdir ./samples_satellite/9_circlenarr"], shell=True)
        if  os.path.exists('./samples_satellite/9_circlenarr/b_monomer'): 
            subprocess.run(["rm -r ./samples_satellite/9_circlenarr/b_monomer"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/9_circlenarr/b_monomer"], shell=True)   

        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)   
        print('Main task started\n')     
        def run_step9_1(one_region):    
            input_file=f"./samples_satellite/4_umap{seed_str}/monomer/{one_region}"
            output=f"./samples_satellite/9_circlenarr/b_monomer/{one_region}"
            i=0
            start_mark=''
            with open(input_file,'r') as f:
                with open(output,'w') as f2:
                    for line in f:
                        i+=1
                        if i==1:
                            f2.write(f"region_name\tcirclen_arr\tstart\tend\tstrand\tseq_arr\tumap2v_arr\tdistance_arr\tarr_serial\tlength\n")
                            continue
                        eachline_arr=line.strip().split('\t')
                        circ_len=int(eachline_arr[10])
                        circ_seq=eachline_arr[5]
                        umap2v=f"{eachline_arr[7]}_{eachline_arr[8]}"
                        distance=eachline_arr[13]
                        end=eachline_arr[3]
                        if circ_len>=26 and circ_len<=30: 
                            if start_mark=='yes' and arr_serial>0:
                                start_mark=''
                                arr_serial+=1
                                circlen_arr.append(str(circ_len))
                                seq_arr.append(circ_seq)
                                umap2v_arr.append(umap2v)
                                distance_arr.append(distance)
                                circlen_arr_str='-'.join(circlen_arr)
                                seq_arr_str='_'.join(seq_arr)
                                umap2v_arr_str='_'.join(umap2v_arr)
                                distance_arr_str='_'.join(distance_arr)
                                region_name=eachline_arr[0]
                                strand=eachline_arr[4]
                                length=abs(int(end)-int(start))+1
                                f2.write(f"{region_name}\t{circlen_arr_str}\t{start}\t{end}\t{strand}\t{seq_arr_str}\t{umap2v_arr_str}\t{distance_arr_str}\t{arr_serial}\t{length}\n")
                            start_mark='yes'
                            arr_serial=0
                            circlen_arr=[]
                            seq_arr=[]
                            umap2v_arr=[]
                            distance_arr=[]
                            continue
                        elif start_mark=='yes':
                            arr_serial+=1
                            if arr_serial==1:start=eachline_arr[2]
                            circlen_arr.append(str(circ_len))
                            seq_arr.append(circ_seq)
                            umap2v_arr.append(umap2v)
                            distance_arr.append(distance)
        
        with Pool(processes=thread) as pool:
            for i, result in enumerate(pool.imap(run_step9_1, region_list), start=1):
                progress = (i / len(region_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush() 
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.0b" or  argv1=="step9.0b_s1": 
        print('step9.0s1 ———— Calculate the length of all regions')
        if  os.path.exists('./samples_satellite/9_circlenarr/b_stat')==False:
            subprocess.run(["mkdir ./samples_satellite/9_circlenarr/b_stat"], shell=True)    
        print('Loading region name list')
        region_all_len=0
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_all_len+= abs(int(eachline_arr[3])-int(eachline_arr[4]))+1
        with open('./samples_satellite/9_circlenarr/b_stat/all_region_length','w') as f:        
            f.write(str(region_all_len))
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.0b" or  argv1=="step9.0b_s2":        
        print('step9.0 ———— Determine preliminary HOR based on 21-25bp monomers, statistics')
        print('Loading all_region_length')
        with open('./samples_satellite/9_circlenarr/b_stat/all_region_length','r') as f:        
            all_region_length=int(f.read().strip())
        print('Loading region name list')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)  
        region_list_len=len(region_list)        
        print('Loading each region')
        dict_circlenarr_num={}
        i=0
        for one_region in region_list:
            i+=1
            print(f'Progress: {i}/{region_list_len}',end='\r')
            input_file=f"./samples_satellite/9_circlenarr/b_monomer/{one_region}"
            k=0
            with open(input_file,'r') as f:
                for line in f:
                    k+=1
                    if k==1:continue
                    eachline_arr=line.strip().split('\t')
                    circlen_arr_str=    eachline_arr[1]
                    if circlen_arr_str not in dict_circlenarr_num :dict_circlenarr_num[circlen_arr_str]=0
                    dict_circlenarr_num[circlen_arr_str]+=1
        sorted_dict_circlenarr_num = dict(sorted(dict_circlenarr_num.items(), key=lambda item:item[1],reverse=True))
    
        with open('./samples_satellite/9_circlenarr/b_stat/0-circlenarr_num','w') as f :
            f.write(f"circlenarr\tone_parts\tone_length\tnum\tall_length\tpercent\n")    
            for circlenarr,num in sorted_dict_circlenarr_num.items():
                one_length= sum(int(num) for num in circlenarr.split('-'))
                all_length=one_length*num
                percent=round(all_length/all_region_length,10)
                one_parts=len(circlenarr.split('-'))
                f.write(f"{circlenarr}\t{one_parts}\t{one_length}\t{num}\t{all_length}\t{percent}\n")          
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.0b" or  argv1=="step9.0b_s3":
        circlenarr_list=[]
        circlenarr_list.append([23,28])
        circlenarr_list.append([23,23,28])
        circlenarr_list.append([23,23,23,28])
        circlenarr_list.append([23,23,23,23,28])
        circlenarr_list.append([23,23,23,23,23,28])
        circlenarr_list.append([23,23,23,23,23,23,28])
        circlenarr_list.append([23,23,23,23,23,23,23,28])
        circlenarr_list.append([23,23,23,23,23,23,23,23,28])
        circlenarr_list.append([23,23,23,23,23,23,23,23,23,28])
        circlenarr_list.append([23,23,23,23,23,23,23,23,23,23,28])
        circlenarr_list_num=len(circlenarr_list)
        dict_k_info={}
        i=0
        with open('./samples_satellite/9_circlenarr/b_stat/0-circlenarr_num2other','w') as f2:    
            with open('./samples_satellite/9_circlenarr/b_stat/0-circlenarr_num','r') as f:
                for line in f:
                    i+=1
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    if len(eachline_arr)!=6:continue
                    circlenarr=         eachline_arr[0]
                    if circlenarr=='circlenarr':continue
                    circlenarr_arr=         circlenarr.split('-')
                    circlenarr_arr_len=     len(circlenarr_arr)
                    k=0;good_state=''
                    while k<circlenarr_list_num:
                        one_circlenarr=circlenarr_list[k]
                        if len(one_circlenarr)==circlenarr_arr_len:
                            kk=0;good_mark=0
                            while kk<circlenarr_arr_len:
                                delta=abs(one_circlenarr[kk]-int(circlenarr_arr[kk]))
                                kk+=1
                                if delta<=2:good_mark+=1
                            if good_mark==circlenarr_arr_len:
                                if k not in dict_k_info:dict_k_info[k]=[]
                                dict_k_info[k].append(eachline)
                                good_state='yes'
                        if good_state=='yes':break
                        k+=1    
                    if good_state=='':
                        f2.write(f"{eachline}\n")
            f2.write(f"\n")             
        k=0
        with open('./samples_satellite/9_circlenarr/b_stat/0-circlenarr_num2','w') as f:    
            while k<circlenarr_list_num:
                headstr="-".join([str(num) for num in circlenarr_list[k]])
                f.write(f"\n##################{headstr}\n")
                for oneinfo in dict_k_info[k]:
                    f.write(f"{oneinfo}\n")
                k+=1
            f.write(f"\n")     
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.0b" or  argv1=="step9.0b_s4":
        start_mark=''
        with open('./samples_satellite/9_circlenarr/b_stat/0-circlenarr_num2s','w') as f2: 
            with open('./samples_satellite/9_circlenarr/b_stat/0-circlenarr_num2','r') as f: 
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    if len(eachline)==0 and start_mark=='yes':
                        f2.write(f'#\t#\t#\t{num_acc}\t{length_acc}\t{percent_acc}\n')
                        num_acc=0;length_acc=0;percent_acc=0;start_mark=''            
                    if len(eachline_arr)==1: 
                        f2.write(eachline+'\n')  
                        num_acc=0;length_acc=0;percent_acc=0;start_mark='yes'
                    elif start_mark=='yes': 
                        num_acc+=int(eachline_arr[3])
                        length_acc+=int(eachline_arr[4])
                        percent_acc+=float(eachline_arr[5])    
    ##step9 is Pic2 monomers composed of subunits     
    if argv1=="stepall" or argv1=="step9_readme":
        print('Print instructions')
        with open('./samples_satellite/9_circlenarr/readme','w') as f:
            txt=r'''
            a_monomer and b_monomer divide monomers using 23 and 28 as boundaries respectively
            sum_monomer1 merges the above two, adding four columns A1A2 and B1B2
            sum_monomer2 converts the four columns into two additional columns, divided into types such as 28-28...-23 and 23-28... etc., totaling about 20 groups plus "other", effectively classifying monomer types
            sum_monomer3 is a statistic; each file contains counts for all types
            sum_monomer6 converts subunit patterns to monomer patterns
            '''
            f.write(txt)
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.1":
        print('Merge the results of 23----23 and 28----28 patterns. After UMAP analysis, add four columns A1A2 and B1B2, indicating the pattern type and current ordinal of each subunit in both patterns')
        if  os.path.exists('./samples_satellite/9_circlenarr/sum_monomer1'): 
            subprocess.run(["rm -r ./samples_satellite/9_circlenarr/sum_monomer1"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/9_circlenarr/sum_monomer1"], shell=True)      
    
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)       
    
        print('Main task started\n')     
        def run_step9_1(one_region):    
            input_file=f"./samples_satellite/4_umap{seed_str}/monomer/{one_region}"
            input_file1=f"./samples_satellite/9_circlenarr/a_monomer/{one_region}"
            input_file2=f"./samples_satellite/9_circlenarr/b_monomer/{one_region}"
            output=f"./samples_satellite/9_circlenarr/sum_monomer1/{one_region}"
            arr1=[]
            with open(input_file1) as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if eachline_arr[0]=='region_name':continue
                    start=          eachline_arr[2]
                    end=            eachline_arr[3]
                    circlenarr=     eachline_arr[1]
                    arr1.append([start,end,circlenarr])
            arr1_len=len(arr1)                
            arr2=[]
            with open(input_file2) as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if eachline_arr[0]=='region_name':continue
                    start=          eachline_arr[2]
                    end=            eachline_arr[3]
                    circlenarr=     eachline_arr[1]
                    arr2.append([start,end,circlenarr])       
            arr2_len=len(arr2)                
            k1=0
            k2=0
            i=0
            serial_1=0;serial_2=0
            A1='';A2=''
            B1='';B2=''
            with open (output,'w') as f2:
                f2.write("region_name\tcirc_serial\tchrserial_1\tchrserial_2\tstrand\tcirc_seq\tumap1v_x\tumap2v_x\tumap2v_y\tin_seat_num\tcirc_len\tnum\tvariant_distance2onehap\tvariant_distance2allhaps\tsubunit1\tserial1\tsubunit2\tserial2\n")
                with open (input_file,'r') as f:
                    for line in f:
                        i+=1
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        if i==1:continue
                        if i==2:
                            if arr1_len>0:
                                arr1_one=arr1[k1]
                                arr1_one_start=     arr1_one[0]
                                arr1_one_end=       arr1_one[1]
                                arr1_one_name=      arr1_one[2]
                            else:
                                arr1_one_start=''
                                arr1_one_end=''
                                arr1_one_name=''
                            if arr2_len>0:    
                                arr2_one=arr2[k2]
                                arr2_one_start=     arr2_one[0]
                                arr2_one_end=       arr2_one[1]
                                arr2_one_name=      arr2_one[2]    
                            else:
                                arr2_one_start=''
                                arr2_one_end=''
                                arr2_one_name=''                            
                        start=      eachline_arr[2]
                        end=        eachline_arr[3]   
                        if start==arr1_one_start:
                            serial_1=1
                            A1=arr1_one_name
                        elif serial_1!=0:
                            serial_1+=1
                        if serial_1>0:
                            A2=serial_1
                        if start==arr2_one_start:
                            serial_2=1
                            B1=arr2_one_name
                        elif serial_2!=0:
                            serial_2+=1
                        if serial_2>0:
                            B2=serial_2
                        if A1=='':A1='.'
                        if A2=='':A2='.'
                        if B1=='':B1='.'
                        if B2=='':B2='.'
                        newline=f'{eachline}\t{A1}\t{A2}\t{B1}\t{B2}\n'
                        f2.write(newline)
                        if end==arr1_one_end:  
                            serial_1=0
                            A1='';A2=''
                            k1+=1
                            if k1>=arr1_len:continue
                            arr1_one=arr1[k1]
                            arr1_one_start=     arr1_one[0]
                            arr1_one_end=       arr1_one[1]
                            arr1_one_name=      arr1_one[2]                        
                        if end==arr2_one_end: 
                            serial_2=0                    
                            B1='';B2=''
                            k2+=1
                            if k2>=arr2_len:continue
                            arr2_one=arr2[k2]
                            arr2_one_start=     arr2_one[0]
                            arr2_one_end=       arr2_one[1]
                            arr2_one_name=      arr2_one[2]  
        
        with Pool(processes=thread) as pool:
            for i, result in enumerate(pool.imap(run_step9_1, region_list), start=1):
                progress = (i / len(region_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush() 
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.2":
        print('Add two more columns. Use the four added columns to determine a standard name. I only count 28*10 patterns')
        if os.path.exists('./samples_satellite/9_circlenarr/sum_monomer2'): 
            subprocess.run(["rm -r ./samples_satellite/9_circlenarr/sum_monomer2"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/9_circlenarr/sum_monomer2"], shell=True)      
    
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)       
    
        print('Main task started\n')     
        def run_step9_2(one_region):    
            input_file=f"./samples_satellite/9_circlenarr/sum_monomer1/{one_region}"
            output_file=f"./samples_satellite/9_circlenarr/sum_monomer2/{one_region}"
            i=0
            with open(output_file,'w') as f2:
                f2.write("region_name\tcirc_serial\tchrserial_1\tchrserial_2\tstrand\tcirc_seq\tumap1v_x\tumap2v_x\tumap2v_y\tin_seat_num\tcirc_len\tnum\tvariant_distance2onehap\tvariant_distance2allhaps\tsubunit1\tserial1\tsubunit2\tserial2\tsubunit_revised_type\tsubunit_revised_serial\n")
                with open(input_file,'r') as f:
                    for line in f:
                        i+=1
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        if i==1:continue
                        subunit1=       eachline_arr[-4];subunit1_arr=subunit1.split('-')	;subunit1_len=len(subunit1_arr)
                        serial1=        eachline_arr[-3]	
                        subunit2=       eachline_arr[-2];subunit2_arr=subunit2.split('-')	;subunit2_len=len(subunit2_arr)	
                        serial2=        eachline_arr[-1]
                        new_type='.';new_serial='.'
                        if subunit1!='.':
                            if subunit1_len>0  :
                                circlenarr_list=[]
                                circlenarr_list.append([28,23])
                                circlenarr_list.append([28,28,23])
                                circlenarr_list.append([28,28,28,23])
                                circlenarr_list.append([28,28,28,28,23])
                                circlenarr_list.append([28,28,28,28,28,23])
                                circlenarr_list.append([28,28,28,28,28,28,23])
                                circlenarr_list.append([28,28,28,28,28,28,28,23])
                                circlenarr_list.append([28,28,28,28,28,28,28,28,23])
                                circlenarr_list.append([28,28,28,28,28,28,28,28,28,23])
                                circlenarr_list.append([28,28,28,28,28,28,28,28,28,28,23])
                                circlenarr_list_num=len(circlenarr_list)                              
                                k=0
                                while k<circlenarr_list_num:
                                    circlenarr_one=circlenarr_list[k]
                                    k+=1
                                    circlenarr_one_len=len(circlenarr_one)
                                    if circlenarr_one_len==subunit1_len:
                                        kk=0;good_num=0
                                        while kk<circlenarr_one_len:
                                            if abs(circlenarr_one[kk]-int(subunit1_arr[kk]))<=2:good_num+=1
                                            kk+=1
                                        if good_num==circlenarr_one_len:        new_type= '-'.join(list(map(str, circlenarr_one)))  ;new_serial=serial1
                                        else:                                   new_type='other1'
                                    if new_type!='other1' and new_type!='.':break
                        else:            
                            if  subunit2_len>=3 and serial2!='.' and int(serial2)>=2:   
                                circlenarr_list=[]
                                circlenarr_list.append([23,23,28])
                                circlenarr_list.append([23,23,23,28])
                                circlenarr_list.append([23,23,23,23,28])
                                circlenarr_list.append([23,23,23,23,23,28])
                                circlenarr_list.append([23,23,23,23,23,23,28])
                                circlenarr_list.append([23,23,23,23,23,23,23,28])
                                circlenarr_list.append([23,23,23,23,23,23,23,23,28])
                                circlenarr_list.append([23,23,23,23,23,23,23,23,23,28])
                                circlenarr_list.append([23,23,23,23,23,23,23,23,23,23,28])
                                circlenarr_list_num=len(circlenarr_list)                    
                                k=0
                                while k<circlenarr_list_num:
                                    circlenarr_one=circlenarr_list[k]
                                    k+=1
                                    circlenarr_one_len=len(circlenarr_one)
                                    if circlenarr_one_len==subunit2_len:
                                        kk=0;good_num=0
                                        while kk<circlenarr_one_len:
                                            if abs(circlenarr_one[kk]-int(subunit2_arr[kk]))<=2:good_num+=1
                                            kk+=1
                                        if good_num==circlenarr_one_len:        
                                            new_type= '-'+'-'.join(list(map(str, circlenarr_one[1:-1])))+'-'  ;
                                            new_serial=int(serial2)-1
                                        else:                                   new_type='other2'
                                    if new_type!='other2' and new_type!='.':break                            
                        newline=f"{eachline}\t{new_type}\t{new_serial}"
                        f2.write(newline+'\n')         
        with Pool(processes=thread) as pool:
            for i, result in enumerate(pool.imap(run_step9_2, region_list), start=1):
                progress = (i / len(region_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()     
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.3":
        print('Add two more columns. Use the four added columns to determine a standard name. I only count 28*10 patterns... statistics')
        if os.path.exists('./samples_satellite/9_circlenarr/sum_monomer3'): 
            subprocess.run(["rm -r ./samples_satellite/9_circlenarr/sum_monomer3"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/9_circlenarr/sum_monomer3"], shell=True)      
    
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)   
        def  run_step9_3(one_region):
            input_file=f"./samples_satellite/9_circlenarr/sum_monomer2/{one_region}"
            output_file=f"./samples_satellite/9_circlenarr/sum_monomer3/{one_region}"
            arr=one_region.split(':')[-1].split('-')
            region_length= int(arr[1])-int(arr[0])+1
            i=0
            dict_subunit_num={}
            with open(input_file) as f:
                for line in f:
                    i+=1
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    if i==1:continue
                    subunit_type=   eachline_arr[-2]
                    subunit_num=    eachline_arr[-1]
                    circlen=        int(eachline_arr[10])
                    if subunit_type not in dict_subunit_num:
                        dict_subunit_num[subunit_type]={}
                        dict_subunit_num[subunit_type]['num']=0
                        dict_subunit_num[subunit_type]['length']=0
                    dict_subunit_num[subunit_type]['num'] +=1
                    dict_subunit_num[subunit_type]['length']+=circlen
            with open(output_file,'w') as f:
                for subunit,dict_num1len in dict_subunit_num.items():
                    num=        dict_num1len['num']
                    length=        dict_num1len['length']
                    percent=length/region_length
                    f.write(f"{one_region}\t{subunit}\t{num}\t{length}\t{percent}\n")
        with Pool(processes=thread) as pool:
            for i, result in enumerate(pool.imap(run_step9_3, region_list), start=1):
                progress = (i / len(region_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()              
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.4":
        print('Add two more columns. Use the four added columns to determine a standard name. I only count 28*10 patterns... statistics')
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)
        subunit_set=set()            
        dict_block_subunit_length={}
        for one_region in region_list:
            input_file=f"./samples_satellite/9_circlenarr/sum_monomer3/{one_region}"
            with open(input_file,'r') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    block           =eachline_arr[0]
                    subunit         =eachline_arr[1]
                    num             =eachline_arr[2]
                    length          =eachline_arr[3]
                    percent         =eachline_arr[4]
                    if block not in  dict_block_subunit_length:dict_block_subunit_length[block]={}
                    dict_block_subunit_length[block][subunit]=length
                    subunit_set.add(subunit)
        subunit_list=list(subunit_set)        
        subunit_list.sort()  
        print('Manually output order')
        subunit_list=[  '28-23',
                        '28-28-23',
                        '28-28-28-23',
                        '28-28-28-28-23',
                        '28-28-28-28-28-23',
                        '28-28-28-28-28-28-23',
                        '28-28-28-28-28-28-28-23',
                        '28-28-28-28-28-28-28-28-23',
                        '28-28-28-28-28-28-28-28-28-23',
                        '28-28-28-28-28-28-28-28-28-28-23',
                        '-23-',
                        '-23-23-',
                        '-23-23-23-',
                        '-23-23-23-23-',
                        '-23-23-23-23-23-',
                        '-23-23-23-23-23-23-',
                        '-23-23-23-23-23-23-23-',
                        '-23-23-23-23-23-23-23-23-',
                        '-23-23-23-23-23-23-23-23-23-',
                        'other1',
                        'other2',
                        '.']
        with open('./samples_satellite/9_circlenarr/sum_monomer4_circlenarr_list','w') as f:
            for one in subunit_list:
                f.write(one+'\n')
        print(subunit_list)
        
        with open("./samples_satellite/9_circlenarr/sum_monomer4_allinfo",'w') as f2:
            i=0
            with open('./samples_satellite/2_good_regions','r') as f:
                for line in f.readlines():
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    i+=1
                    if i==1:
                        newline=eachline+'\t'+'\t'.join(subunit_list)+'\n'
                    elif len(eachline_arr)!=8:
                        newline='\n'
                    else:        
                        block=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                        length_list=[]
                        for one in subunit_list:
                            if block not in dict_block_subunit_length:
                                length='0'
                            elif one not in  dict_block_subunit_length[block]:
                                length='0'
                            else:    
                                length= dict_block_subunit_length[block][one]
                            length_list.append(length)
                        newline=eachline+'\t'+'\t'.join(length_list)+'\n'
                    f2.write(newline)
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.5s1":               
        print('Sample statistics — length percentage')
        dict_hap_info={}
        dict_sample_type={}
        with open('./samples_satellite/sample_info','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=6:continue
                sample_name_hap=    eachline_arr[0]
                sample_name=        eachline_arr[1]
                sample_type=        eachline_arr[3]
                dict_hap_info[sample_name_hap]=[sample_name,sample_type]
                dict_sample_type[sample_name]=sample_type
        
        subunitgroup_list=[]
        i=0
        dict_sample_subunit_length={}
        dict_sample_sumlength={}
        with open('./samples_satellite/9_circlenarr/sum_monomer4_allinfo','r') as f:
            for line in f:
                i+=1
                eachline_arr=line.strip().split('\t')
                eachline_arr_len=len(eachline_arr)
                if eachline_arr_len<9:continue
                if i==1:
                    k=8
                    while k<eachline_arr_len :
                        subunitgroup=eachline_arr[k]
                        subunitgroup_list.append(subunitgroup)
                        k+=1
                else:
                    sample_hap=eachline_arr[0]
                    if sample_hap=='PN40024':continue
                    sample_name=dict_hap_info[sample_hap][0]
                    sample_type=dict_hap_info[sample_hap][1]
                    if sample_type=="Routundifolia":continue
                    if sample_name not in dict_sample_sumlength:
                        dict_sample_sumlength[sample_name]=0
                    dict_sample_sumlength[sample_name]+=int(eachline_arr[5])    
                    k=8
                    while k<eachline_arr_len :
                        subunitgroup=subunitgroup_list[k-8]
                        if sample_name not in dict_sample_subunit_length:
                            dict_sample_subunit_length[sample_name]={}
                        if subunitgroup not in dict_sample_subunit_length[sample_name]:
                            dict_sample_subunit_length[sample_name][subunitgroup]=0
                        dict_sample_subunit_length[sample_name][subunitgroup]+=int(eachline_arr[k])
                        k+=1
        with open('./samples_satellite/9_circlenarr/sum_monomer4_allinfo_stat1','w') as f:   
            f.write("sample\tsample_type\tsubunitgroup\tlength\tpercent\n")
            for sample,subunitgroup_length in dict_sample_subunit_length.items():
                sumlength=dict_sample_sumlength[sample]
                sample_type=dict_sample_type[sample]
                if sample_type=="Routundifolia":continue
                if sample_type=="hybrid":continue
                for subunitgroup,length in subunitgroup_length.items():
                    f.write(f"{sample}\t{sample_type}\t{subunitgroup}\t{length}\t{length/sumlength}\n")
                    
        print('Simplify the plot')
        dict_raw_new={
            '28-23':'(23)-28-23', 
            '28-28-23':'(23)-28x2-23', 
            '28-28-28-23':'(23)-28x3-23', 
            '28-28-28-28-23':'(23)-28x4-23', 
            '28-28-28-28-28-23':'(23)-28x5-23', 
            '28-28-28-28-28-28-23':'(23)-28x6-23', 
            '28-28-28-28-28-28-28-23':'(23)-28x7-23', 
            '28-28-28-28-28-28-28-28-23':'(23)-28x8-23', 
            '28-28-28-28-28-28-28-28-28-23':'(23)-28x9-23', 
            '28-28-28-28-28-28-28-28-28-28-23':'(23)-28x10-23', 
            '-23-':'(28-23)-23xn-(28)', 
            '-23-23-':'(28-23)-23xn-(28)', 
            '-23-23-23-':'(28-23)-23xn-(28)', 
            '-23-23-23-23-':'(28-23)-23xn-(28)', 
            '-23-23-23-23-23-':'(28-23)-23xn-(28)', 
            '-23-23-23-23-23-23-':'(28-23)-23xn-(28)', 
            '-23-23-23-23-23-23-23-':'(28-23)-23xn-(28)', 
            '-23-23-23-23-23-23-23-23-':'(28-23)-23xn-(28)', 
            '-23-23-23-23-23-23-23-23-23-':'(28-23)-23xn-(28)', 
            'other1':'other', 
            'other2':'other', 
            '.':'other'
        }
        print(dict_raw_new)
        dict_sample_subunit_length={}
        with open('./samples_satellite/9_circlenarr/sum_monomer4_allinfo_stat1','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)<5:continue
                sample	        =eachline_arr[0]
                sample_type	    =eachline_arr[1]
                subunitgroup	=eachline_arr[2]  
                subunitgroup_new =dict_raw_new[subunitgroup]
                length	        =eachline_arr[3]
                percent         =eachline_arr[4]
                if sample not in dict_sample_subunit_length:dict_sample_subunit_length[sample]={}
                if subunitgroup_new not in dict_sample_subunit_length[sample]:
                    dict_sample_subunit_length[sample][subunitgroup_new]={}
                    dict_sample_subunit_length[sample][subunitgroup_new]=0
                dict_sample_subunit_length[sample][subunitgroup_new]+=int(length)
    
        with open('./samples_satellite/9_circlenarr/sum_monomer4_allinfo_stat2','w') as f:   
            f.write("sample\tsample_type\tsubunitgroup\tlength\tpercent\n")
            for sample,subunitgroup_length in dict_sample_subunit_length.items():
                sumlength=dict_sample_sumlength[sample]
                sample_type=dict_sample_type[sample]
                for subunitgroup,length in subunitgroup_length.items():
                    f.write(f"{sample}\t{sample_type}\t{subunitgroup}\t{length}\t{length/sumlength}\n")
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.5p1":  
        print('step9.5p — plotting')
        R_txt=r'''library(ggplot2)
    library(dplyr)
    
    setwd('./')
    input_file1 <- read.table('sum_monomer4_allinfo_stat2', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    
    subunitgroups <- c("(23)-28-23", "(23)-28x2-23", "(23)-28x3-23", "(23)-28x4-23", 
                       "(23)-28x5-23", "(23)-28x6-23", "(23)-28x7-23", "(23)-28x8-23", 
                       "(23)-28x9-23", "(23)-28x10-23", "(28-23)-23xn-(28)", "other")
    reversed_subunitgroups <- rev(subunitgroups)
    input_file1$subunitgroup <- factor(input_file1$subunitgroup, 
          levels = reversed_subunitgroups)
    
    p=ggplot(input_file1, aes(x = subunitgroup, y = percent)) +
      geom_boxplot(
        alpha = 0.5,
        outlier.shape = NA,
        width = 0.4
      ) +
      geom_jitter(
        aes(color = sample_type),
        width = 0.1,
        height = 0,
        alpha = 0.7,
        size = 0.5
      ) +
      coord_flip() +
      ylim(0.7, 0.9) + 
      labs(
        title = "Distribution of percent by Subunit Group",
        x = "Percent Value",
        y = "Subunit Group"
      ) +
      theme_classic() 
    
    pdf("sum_monomer4_allinfo_stat2.pdf", width = 50 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()
    '''
        new_directory = "./samples_satellite/9_circlenarr/"
        os.chdir(new_directory)
        with open('pic_ggplot_step9.5p1.R','w',encoding='utf-8') as f:
            f.write(R_txt)    
        subprocess.run(['Rscript pic_ggplot_step9.5p1.R'], shell=True)    
        os.chdir('../../')                                           
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.5s2":
        print('Chromosome monomer statistics')
        dict_hap_info={}
        dict_sample_type={}
        with open('./samples_satellite/sample_info','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=6:continue
                sample_name_hap=    eachline_arr[0]
                sample_name=        eachline_arr[1]
                sample_type=        eachline_arr[3]
                dict_hap_info[sample_name_hap]=[sample_name,sample_type]
                dict_sample_type[sample_name]=sample_type
        print('Simplify the plot')
        dict_raw_new={
            '28-23':'(23)-28-23', 
            '28-28-23':'(23)-28x2-23', 
            '28-28-28-23':'(23)-28x3-23', 
            '28-28-28-28-23':'(23)-28x4-23', 
            '28-28-28-28-28-23':'(23)-28x5-23', 
            '28-28-28-28-28-28-23':'(23)-28x6-23', 
            '28-28-28-28-28-28-28-23':'(23)-28x7-23', 
            '28-28-28-28-28-28-28-28-23':'(23)-28x8-23', 
            '28-28-28-28-28-28-28-28-28-23':'(23)-28x9-23', 
            '28-28-28-28-28-28-28-28-28-28-23':'(23)-28x10-23', 
            '-23-':'(28-23)-23xn-(28)', 
            '-23-23-':'(28-23)-23xn-(28)', 
            '-23-23-23-':'(28-23)-23xn-(28)', 
            '-23-23-23-23-':'(28-23)-23xn-(28)', 
            '-23-23-23-23-23-':'(28-23)-23xn-(28)', 
            '-23-23-23-23-23-23-':'(28-23)-23xn-(28)', 
            '-23-23-23-23-23-23-23-':'(28-23)-23xn-(28)', 
            '-23-23-23-23-23-23-23-23-':'(28-23)-23xn-(28)', 
            '-23-23-23-23-23-23-23-23-23-':'(28-23)-23xn-(28)', 
            'other1':'other', 
            'other2':'other', 
            '.':'other'
        }
        subunitgroup_list=[]
        i=0
        chr_list=[]
        dict_subunit_chr_length={}
        dict_subunit_sumlength={}
        subunitgroup_new_list=[]
        dict_newsubunit_chr_percent={}
        with open('./samples_satellite/9_circlenarr/sum_monomer4_allinfo','r') as f:
            for line in f:
                i+=1
                eachline_arr=line.strip().split('\t')
                eachline_arr_len=len(eachline_arr)
                if eachline_arr_len<9:continue
                chromosome=eachline_arr[1].replace('region_','Chr')
                if i==1:
                    k=8
                    while k<eachline_arr_len :
                        subunitgroup=eachline_arr[k]
                        subunitgroup_list.append(subunitgroup)
                        k+=1
                else:
                    if chromosome not in chr_list:chr_list.append(chromosome)
                    sample_hap=eachline_arr[0]
                    if sample_hap=='PN40024':continue
                    sample_name=dict_hap_info[sample_hap][0]
                    sample_type=dict_hap_info[sample_hap][1]
                    if sample_type=="Routundifolia":continue
                    k=8
                    while k<eachline_arr_len :
                        subunitgroup=subunitgroup_list[k-8]
                        subunitgroup_new=dict_raw_new[subunitgroup]
                        if subunitgroup_new not in subunitgroup_new_list:subunitgroup_new_list.append(subunitgroup_new)
                        if subunitgroup_new not in dict_subunit_sumlength:
                            dict_subunit_sumlength[subunitgroup_new]=0
                        dict_subunit_sumlength[subunitgroup_new]+=int(eachline_arr[k])    
                        if subunitgroup_new not in dict_subunit_chr_length:
                            dict_subunit_chr_length[subunitgroup_new]={}
                        if chromosome not in dict_subunit_chr_length[subunitgroup_new]:
                            dict_subunit_chr_length[subunitgroup_new][chromosome]=0
                        dict_subunit_chr_length[subunitgroup_new][chromosome]+=int(eachline_arr[k])
                        k+=1  
        with open('./samples_satellite/9_circlenarr/sum_monomer4_allinfo_stat3','w') as f:   
            f.write("subunitgroup\tchromosome\tlength\tpercent\n")
            for subunitgroup_new,chromosome_length in dict_subunit_chr_length.items():
                sumlength=dict_subunit_sumlength[subunitgroup_new]
                for chromosome,length in chromosome_length.items():
                    if subunitgroup_new not in dict_newsubunit_chr_percent:dict_newsubunit_chr_percent[subunitgroup_new]={}
                    dict_newsubunit_chr_percent[subunitgroup_new][chromosome]=length/sumlength
                    f.write(f"{subunitgroup_new}\t{chromosome}\t{length}\t{length/sumlength}\n")            
       
        print('Convert to heatmap data format')                
        with open('./samples_satellite/9_circlenarr/sum_monomer4_allinfo_stat4','w') as f:   
            newline='subunit'
            for one_chr in chr_list:
                if one_chr=="Chr20":continue
                newline+='\t'+one_chr
            f.write(newline+'\n')
            for subunitgroup_new_one in subunitgroup_new_list:
                newline=subunitgroup_new_one
                for one_chr in chr_list:
                    if one_chr=="Chr20":continue
                    newline+="\t"+str(dict_newsubunit_chr_percent[subunitgroup_new_one][one_chr])
                f.write(newline+'\n')    
    if argv1=="stepall" or argv1=="step9"  or  argv1=="step9.5p2":  
        print('step9.5p — plotting')
        R_txt=r'''library(ggplot2)
    library(pheatmap)
    
    setwd('./')
    input_file1 <- read.table('sum_monomer4_allinfo_stat4', skip = 0, header = TRUE, stringsAsFactors = FALSE,row.names=1, check.names = FALSE, sep = '\t')
    
    pdf("sum_monomer4_allinfo_stat4.pdf", width = 10 / 2.54, height = 10 / 2.54)
    
    pheatmap(input_file1, 
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             scale = "none",
             display_numbers = FALSE,
             color = colorRampPalette(c("white", "red"))(50)
    )   
    dev.off()
    
    input_file2 <- read.table('sum_monomer4_allinfo_stat3', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    
    subunitgroups <- c("(23)-28-23", "(23)-28x2-23", "(23)-28x3-23", "(23)-28x4-23", 
                       "(23)-28x5-23", "(23)-28x6-23", "(23)-28x7-23", "(23)-28x8-23", 
                       "(23)-28x9-23", "(23)-28x10-23", "(28-23)-23xn-(28)", "other")
    reversed_subunitgroups <- rev(subunitgroups)
    input_file2$subunitgroup <- factor(input_file2$subunitgroup, 
                                       levels = reversed_subunitgroups)
    
    p=ggplot(input_file2, aes(x =subunitgroup, y = percent , fill = percent)) +
      geom_bar(stat = "identity", position = "stack",color='white') +
      theme_classic()+
      coord_flip() 
    print(p)
    
    pdf("sum_monomer4_allinfo_stat3.pdf", width = 10 / 2.54, height = 10 / 2.54)
    print(p)
    dev.off()  
    '''
        new_directory = "./samples_satellite/9_circlenarr/"
        os.chdir(new_directory)
        with open('pic_ggplot_step9.5p2.R','w',encoding='utf-8') as f:
            f.write(R_txt)    
        subprocess.run(['Rscript pic_ggplot_step9.5p2.R'], shell=True)    
        os.chdir('../../')                                                           
if (argv1=="stepall" or "step9" in argv1) and subunit_mark=='no':     
    if 'seed'=='seed':
        import random
        random_seed=7309 
        if "seed"  in args_dict: #print("Missing input fasta file");sys.exit()
            random_seed = int(args_dict["seed"])
        elif 'random_seed' in locals():
            print(f'umap manually set seed to {random_seed}')  
        else:
            print("Random seed")
            random_seed = random.randint(1, 10000)  # You can choose any integer as the random seed
            #random_seed=355
        seed_str=f'_seed{random_seed}'    
    print("Normal satellite monomer, no subunit structure. Generate simplified sum_monomer2 for step11 compatibility")
    
    # Create necessary directories
    if os.path.exists('./samples_satellite/9_circlenarr')==False:
        subprocess.run(["mkdir ./samples_satellite/9_circlenarr"], shell=True)
    
    # Create sum_monomer2 directory (required by step11)
    if os.path.exists('./samples_satellite/9_circlenarr/sum_monomer2'): 
        subprocess.run(["rm -r ./samples_satellite/9_circlenarr/sum_monomer2"], shell=True)   
    subprocess.run(["mkdir ./samples_satellite/9_circlenarr/sum_monomer2"], shell=True)
    
    # Also create sum_monomer1 and sum_monomer3 for consistency
    if os.path.exists('./samples_satellite/9_circlenarr/sum_monomer1'): 
        subprocess.run(["rm -r ./samples_satellite/9_circlenarr/sum_monomer1"], shell=True)   
    subprocess.run(["mkdir ./samples_satellite/9_circlenarr/sum_monomer1"], shell=True)
    
    if os.path.exists('./samples_satellite/9_circlenarr/sum_monomer3'): 
        subprocess.run(["rm -r ./samples_satellite/9_circlenarr/sum_monomer3"], shell=True)   
    subprocess.run(["mkdir ./samples_satellite/9_circlenarr/sum_monomer3"], shell=True)
    
    print('Loading regions')
    region_list=[]
    with open('./samples_satellite/2_good_regions','r') as f:
        for line in f.readlines():
            eachline_arr=line.strip().split('\t')
            if eachline_arr[0]=='sample':continue
            if len(eachline_arr)!=8:continue
            new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
            region_list.append(new_name)
    
    print('Main task started\n')
    print('Generating simplified sum_monomer2 files (each monomer treated as a single subunit)')
    
    def run_step9_simple(one_region):
        """Generate sum_monomer2 format for normal satellites where each monomer is a single subunit"""
        input_file = f"./samples_satellite/4_umap{seed_str}/monomer/{one_region}"
        output_file = f"./samples_satellite/9_circlenarr/sum_monomer2/{one_region}"
        
        with open(output_file, 'w') as f2:
            # Write header (matching original sum_monomer2 format)
            f2.write("region_name\tcirc_serial\tchrserial_1\tchrserial_2\tstrand\tcirc_seq\tumap1v_x\tumap2v_x\tumap2v_y\tin_seat_num\tcirc_len\tnum\tvariant_distance2onehap\tvariant_distance2allhaps\tsubunit1\tserial1\tsubunit2\tserial2\tsubunit_revised_type\tsubunit_revised_serial\n")
            
            k=0
            with open(input_file, 'r') as f:
                for line in f:
                    k+=1
                    if k==1: continue  # Skip header
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    
                    # Extract basic information
                    circ_seq = eachline_arr[5]
                    circ_len = eachline_arr[10]
                    
                    # For normal satellites, each monomer is itself a subunit
                    # subunit_revised_type = monomer length (e.g., "171" for human alpha satellite)
                    # subunit_revised_serial = always "1" (no internal subunits)
                    newline = f"{eachline}\t{circ_len}\t1\t.\t.\t{circ_len}\t1\n"
                    f2.write(newline)
    
    with Pool(processes=thread) as pool:
        for i, result in enumerate(pool.imap(run_step9_simple, region_list), start=1):
            progress = (i / len(region_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush()
    
    print('\nGenerating simplified sum_monomer1 files')
    
    def run_step9_simple_sum1(one_region):
        """Generate sum_monomer1 format (simplified, with placeholder columns)"""
        input_file = f"./samples_satellite/4_umap{seed_str}/monomer/{one_region}"
        output_file = f"./samples_satellite/9_circlenarr/sum_monomer1/{one_region}"
        
        with open(output_file, 'w') as f2:
            f2.write("region_name\tcirc_serial\tchrserial_1\tchrserial_2\tstrand\tcirc_seq\tumap1v_x\tumap2v_x\tumap2v_y\tin_seat_num\tcirc_len\tnum\tvariant_distance2onehap\tvariant_distance2allhaps\tsubunit1\tserial1\tsubunit2\tserial2\n")
            
            k=0
            with open(input_file, 'r') as f:
                for line in f:
                    k+=1
                    if k==1: continue
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    
                    circ_len = eachline_arr[10]
                    # Each monomer is a single subunit, no A/B pattern distinction
                    newline = f"{eachline}\t{circ_len}\t1\t.\t.\n"
                    f2.write(newline)
    
    with Pool(processes=thread) as pool:
        for i, result in enumerate(pool.imap(run_step9_simple_sum1, region_list), start=1):
            progress = (i / len(region_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush()
    
    print('\nGenerating simplified sum_monomer3 statistics')
    
    def run_step9_simple_sum3(one_region):
        """Generate simplified statistics"""
        input_file = f"./samples_satellite/9_circlenarr/sum_monomer2/{one_region}"
        output_file = f"./samples_satellite/9_circlenarr/sum_monomer3/{one_region}"
        
        arr=one_region.split(':')[-1].split('-')
        region_length= int(arr[1])-int(arr[0])+1
        
        dict_subunit_num={}
        with open(input_file, 'r') as f:
            next(f)  # Skip header
            for line in f:
                eachline_arr=line.strip().split('\t')
                subunit_type = eachline_arr[-2]
                circlen = int(eachline_arr[10])
                if subunit_type not in dict_subunit_num:
                    dict_subunit_num[subunit_type] = {'num': 0, 'length': 0}
                dict_subunit_num[subunit_type]['num'] += 1
                dict_subunit_num[subunit_type]['length'] += circlen
        
        with open(output_file, 'w') as f:
            for subunit, dict_num1len in dict_subunit_num.items():
                num = dict_num1len['num']
                length = dict_num1len['length']
                percent = length/region_length
                f.write(f"{one_region}\t{subunit}\t{num}\t{length}\t{percent}\n")
    
    with Pool(processes=thread) as pool:
        for i, result in enumerate(pool.imap(run_step9_simple_sum3, region_list), start=1):
            progress = (i / len(region_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush()
    
    print('\nstep9 (simple mode) completed. Files are ready for step11.')
    print('Note: This is the simplified mode for normal satellites. No subunit pattern analysis was performed.')
        

# Analyze the evolutionary history of monomer subunits
if (argv1=="stepall" or "step10" in argv1)  and subunit_mark=='yes':
    if 'seed'=='seed':
        import random
        random_seed=7309 
        if "seed"  in args_dict: #print("Missing input fasta file");sys.exit()
            random_seed = int(args_dict["seed"])
        elif 'random_seed' in locals():
            print(f'umap manually set seed to {random_seed}')  
        else:
            print("Random seed")
            random_seed = random.randint(1, 10000)  # You can choose any integer as the random seed
            #random_seed=355
        seed_str=f'_seed{random_seed}'
    if argv1=="stepall" or argv1=="step10_readme":
        print('Print instructions')
        with open('./samples_satellite/10_monomer_model/readme','w') as f:
            txt=r'''
            10.1 ———— Preset 10 monomer types, extract corresponding monomers from the population (0_raw folder)
            Deprecated — 10.2_raw1, calculate the number of base variations between each pair of the four subunits of 28-28-28-28-23
            Deprecated — 10.2_raw2, seems useless, such results are too few
            Deprecated — 10.2 ————
            Deprecated — 10.2s ————
            10.3
            Deprecated — 10.3s_a1/10.3s_a2 ———— Calculate averages, the average results are meaningless
            10.3s_b1/10.3s_b2 ———— Good, find the mode, then classify by simple colors for plotting, simple and clear
            '''
            f.write(txt)
    if argv1=="stepall" or   argv1=="step10"  or   argv1=="step10.1":     
        print('Analyze the evolutionary history of monomer subunits. Note that unlike step9, each subunit allows 2bp length variation. Step10 analysis does not allow length variation')
        target_list=[]
        target_list.append("28-23")
        target_list.append("28-28-23")
        target_list.append("28-28-28-23")
        target_list.append("28-28-28-28-23")
        target_list.append("28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-28-28-23")
        
        
        if os.path.exists('./samples_satellite/10_monomer_model')==False: 
            subprocess.run(["mkdir ./samples_satellite/10_monomer_model"], shell=True)   
        
        for target in target_list:
        
            print(f'\nAnalyzing {target}')    
            subprocess.run([f"mkdir ./samples_satellite/10_monomer_model/{target}"], shell=True)          
            subprocess.run([f"mkdir ./samples_satellite/10_monomer_model/{target}/monomer"], shell=True) 
            
            
            print('Loading regions')
            info_list=[]
            with open('./samples_satellite/2_good_regions','r') as f:
                for line in f.readlines():
                    eachline_arr=line.strip().split('\t')
                    if eachline_arr[0]=='sample':continue
                    if len(eachline_arr)!=8:continue
                    new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                    info_list.append([new_name,target])    
            
            def run_step10_1(one_info):
                one_region,target=one_info
                input_file=f"./samples_satellite/9_circlenarr/a_monomer/{one_region}"
                output_file=f"./samples_satellite/10_monomer_model/{target}/monomer/{one_region}"
                
                with open(output_file,'w') as f2:
                    with open(input_file) as f:
                        for line in f:
                            eachline=line.strip()
                            eachline_arr=eachline.split('\t')
                            if eachline_arr[1]==target:
                                f2.write(one_region+'\t'+eachline+'\n')
        
            with Pool(processes=thread) as pool:
                for i, result in enumerate(pool.imap(run_step10_1, info_list), start=1):
                    progress = (i / len(info_list)) * 100
                    sys.stdout.write(f"\rProgress: {progress:.2f}%")
                    sys.stdout.flush()                     
            subprocess.run([f"find ./samples_satellite/10_monomer_model/{target}/monomer/ -type f -print0 | xargs -0 cat > ./samples_satellite/10_monomer_model/{target}/0_raw"], shell=True)                     
    if  argv1=="step10.2_raw1":                      
        target="28-28-28-28-23"
        ### Calculate the difference between two equal-length sequences
        def count_differences(str1, str2):
            difference_count = 0
            for char1, char2 in zip(str1, str2):
                if char1 != char2:
                    difference_count += 1
            return difference_count    
        i=0
        with open(f"./samples_satellite/10_monomer_model/{target}/1_vs",'w') as f2:
            with open(f"./samples_satellite/10_monomer_model/{target}/0_raw",'r') as f:
                for line in f:
                    i+=1
                    print(f"{i}",end='\r')
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    seq_arr=eachline_arr[6].split('_')
                    seq1=seq_arr[0]
                    seq2=seq_arr[1]
                    seq3=seq_arr[2]
                    seq4=seq_arr[3]
                    seq1_2=count_differences(seq1,seq2)
                    seq1_3=count_differences(seq1,seq3)
                    seq1_4=count_differences(seq1,seq4)
                    seq2_3=count_differences(seq2,seq3)
                    seq2_4=count_differences(seq2,seq4)
                    seq3_4=count_differences(seq3,seq4)
                    newline=f"{eachline}\t{seq1_2}\t{seq1_3}\t{seq1_4}\t{seq2_3}\t{seq2_4}\t{seq3_4}\n"
                    f2.write(newline)
    if  argv1=="step10.2_raw2":     
        print('Pairwise comparison of subunits is not very effective, with few results, because only CEN135 has two consecutive subunits co-migrating and replicating')
        target_list=[]
        #target_list.append("28-28-23")
        #target_list.append("28-28-28-23")
        target_list.append("28-28-28-28-23")
        #target_list.append("28-28-28-28-28-23")
        #target_list.append("28-28-28-28-28-28-23")
        #target_list.append("28-28-28-28-28-28-28-23")
        #target_list.append("28-28-28-28-28-28-28-28-23")
        #target_list.append("28-28-28-28-28-28-28-28-28-23")
        #target_list.append("28-28-28-28-28-28-28-28-28-28-23")
        ####
        def find_most_similar_pair(strings):
            n = len(strings)
            min_diff = float('inf')
            most_similar_pair = None
            for i in range(n):
                for j in range(i + 1, n):
                    diff = count_differences(strings[i], strings[j])
                    if diff < min_diff:
                        min_diff = diff
                        most_similar_pair = (i, j)
            if most_similar_pair is not None:
                i, j = most_similar_pair
                all_pairs_satisfy = True
                for k in range(n):
                    for l in range(k + 1, n):
                        if (k == i and l == j) or (k == j and l == i):
                            continue
                        current_diff = count_differences(strings[k], strings[l])
                        if current_diff <= 6:
                            all_pairs_satisfy = False
                            break
                    if not all_pairs_satisfy:
                        break
         
                if all_pairs_satisfy and min_diff < 3:
                    return f"{i + 1}-{j + 1}"  # Return 1-based index
                else:
                    return 'unknown'
            else:
                return 'unknown'
        def count_differences(str1, str2):
            if len(str1) != len(str2):
                raise ValueError("The strings must be of the same length.")
            difference_count = 0
            for char1, char2 in zip(str1, str2):
                if char1 != char2:
                    difference_count += 1
            return difference_count
        
        ####
        for target in target_list:
            print(f'\nAnalyzing {target}')    
    
            input_file=f"./samples_satellite/10_monomer_model/{target}/0_raw"
            output_file=f"./samples_satellite/10_monomer_model/{target}/1_vs"
            
            with open(output_file,'w') as f2:
                with open(input_file) as f:
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        seq_arr=eachline_arr[6].split('_')[:-1]                        
                        result=find_most_similar_pair(seq_arr)
                        newline=eachline+'\t'+result+'\n'
                        f2.write(newline)
    if argv1=="stepall" or   argv1=="step10"  or   argv1=="step10.2":                    
        print('Analyze all pairwise monomers, add headers, calculate the number of differing bases between each monomer')
        target_list=[]
        target_list.append("28-28-23")
        target_list.append("28-28-28-23")
        target_list.append("28-28-28-28-23")
        target_list.append("28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-28-28-23")    
        ### Calculate the difference between two equal-length sequences
        def count_differences(str1, str2):
            difference_count = 0
            for char1, char2 in zip(str1, str2):
                if char1 != char2:
                    difference_count += 1
            return difference_count  
            
        def generate_pairs_manual(num):   #num=4#[(1, 2), (2, 3), (3, 4), (1, 3), (2, 4), (1, 4)]
            pairs = []
            for i in range(1, num + 1):
                for j in range(i + 1, num + 1):
                    pairs.append((i, j))
            pairs.sort(key=lambda x: x[1] - x[0])
            return pairs
        ####
        for target in target_list:
            print(f'\nAnalyzing {target}')    
            pairs_list=generate_pairs_manual(len(target.split('-'))-1)
            input_file=f"./samples_satellite/10_monomer_model/{target}/0_raw"
            output_file=f"./samples_satellite/10_monomer_model/{target}/1_vs"
            
            with open(output_file,'w') as f2:
                head_addition=''
                for one_pair in pairs_list:
                    head_addition+=f"\t#{one_pair[0]}-{one_pair[1]}"
                f2.write(f"block\tchr\tcirclen_arr\tstart\tend\tstrand\tseq_arr\tumap2v_arr\tdistance_arr\tarr_serial\tlength{head_addition}\n")
                with open(input_file,'r') as f:
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        seq_arr=eachline_arr[6].split('_')[:-1]                        
                        addition_str=''
                        for one_pair in pairs_list:
                            index1=one_pair[0]-1
                            index2=one_pair[1]-1
                            difference_count=count_differences(seq_arr[index1],seq_arr[index2])
                            addition_str+=f'\t{difference_count}'
                        newline=eachline+addition_str+'\n'
                        f2.write(newline)                        
    if argv1=="stepall" or   argv1=="step10"  or   argv1=="step10.2s":
        print('Statistics')    
        target_list=[]
        #target_list.append("28-23")
        target_list.append("28-28-23")
        target_list.append("28-28-28-23")
        target_list.append("28-28-28-28-23")
        target_list.append("28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-28-28-23") 
        ## Remove the top and bottom 10% and then calculate the average
        def calculate_trimmed_average(lst, trim_percentage=0.1):
            if not lst:  # Check if the list is empty
                return 0  # Or raise an exception, depending on requirements
            
            n = len(lst)
            trim_count = int(n * trim_percentage)  # Calculate the number of elements to remove
            
            # Remove the first and last trim_count elements
            trimmed_list = lst[trim_count:n - trim_count]
            
            # Calculate the average of the remaining elements
            return sum(trimmed_list) / len(trimmed_list)
        
        with open(f'./samples_satellite/10_monomer_model/monomer_vs','w') as f2:
            for target in target_list:  
                print(target)
                f2.write(target+'\n')
                #################
                dict_twopos_list={}
                i=0        
                with open(f'./samples_satellite/10_monomer_model/{target}/1_vs','r') as f:
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=line.split('\t')
                        i+=1
                        if i==1:
                            eachline_arr_num=len(eachline_arr)
                            twopos_list=eachline_arr[11:]
                            current_twopos_list=[] 
                            for one in twopos_list:
                                if one not in dict_twopos_list:
                                    dict_twopos_list[one]=[]
                                current_twopos_list.append(one)
                        else:        
                            k=11
                            while k< eachline_arr_num:
                                dict_twopos_list[current_twopos_list[k-11]].append(int(eachline_arr[k]))
                                k+=1
                ####
                str1=''
                str2=''
                for twopos,var_list in dict_twopos_list.items():
                    Trimmed_Mean=calculate_trimmed_average(var_list)        
                    str1+=f'\t{twopos}'
                    str2+=f'\t{Trimmed_Mean}'
                f2.write(f"{str1.strip()}\n{str2.strip()}\n") 
                f2.write(f"\n")     
    if argv1=="stepall" or   argv1=="step10"  or   argv1=="step10.3":   
        print('Calculating pairwise base distances is too complicated; using umap1v to assign colors may be more intuitive')
        target_list=[]
        target_list.append("28-23")
        target_list.append("28-28-23")
        target_list.append("28-28-28-23")
        target_list.append("28-28-28-28-23")
        target_list.append("28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-28-28-23")
        
        print('Loading monomer_vector_umap1v')
        dict_circseq_umap1vx={}
        with open(f"./samples_satellite/4_umap{seed_str}/monomer_vector_umap1v",'r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                dict_circseq_umap1vx[eachline_arr[0]]=eachline_arr[1]
        for target in target_list:
            print(target)
            with open(f'./samples_satellite/10_monomer_model/{target}/1_umap1v','w') as f2:
                f2.write(f"block\tchr\tcirclen_arr\tstart\tend\tstrand\tseq_arr\tumap2v_arr\tdistance_arr\tarr_serial\tlength\tumap1vx_list\n")
                with open(f'./samples_satellite/10_monomer_model/{target}/0_raw','r') as f:
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=line.split('\t')
                        circseq_arr=    eachline_arr[6].split('_')
                        umap1vx_list=[]
                        for one in circseq_arr:
                            umap1vx_list.append(dict_circseq_umap1vx[one])
                        umap1vx_list_str='|'.join(umap1vx_list)    
                        newline=eachline+'\t'+umap1vx_list_str+'\n'
                        f2.write(newline)
    if argv1=="stepall" or   argv1=="step10"  or   argv1=="step10.3s_a1":   
        print('Statistics, #### This can be optimized further later. Taking trimmed averages is not very accurate. After averaging, it may not represent any real sequence...')
        target_list=[]
        target_list.append("28-23")
        target_list.append("28-28-23")
        target_list.append("28-28-28-23")
        target_list.append("28-28-28-28-23")
        target_list.append("28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-28-28-23")    
        ## Remove the top and bottom 10% and then calculate the average
        def calculate_trimmed_average(lst, trim_percentage=0.1):
            if not lst:  # Check if the list is empty
                return 0  # Or raise an exception, depending on requirements
            
            n = len(lst)
            trim_count = int(n * trim_percentage)  # Calculate the number of elements to remove
            
            # Remove the first and last trim_count elements
            trimmed_list = lst[trim_count:n - trim_count]
            
            # Calculate the average of the remaining elements
            return sum(trimmed_list) / len(trimmed_list)
        ##
        with open(f'./samples_satellite/10_monomer_model/monomer_umap1v','w') as f2:
            for target in target_list:
                print(target)                        
                f2.write(f"target\tumap1vx_list\n")
                target_len=len(target.split('-'))
                k=0;dict_index_list={}
                while k<target_len:
                    dict_index_list[k]=[]
                    k+=1
                
                
                with open(f'./samples_satellite/10_monomer_model/{target}/1_umap1v','r') as f:
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=line.split('\t')
                        if eachline_arr[0]=='block':continue
                        circseq_arr=    eachline_arr[-1].split('|')
                        k=0
                        while k< target_len:
                            dict_index_list[k].append(float(circseq_arr[k]))
                            k+=1
                addition_list=[]
                for index,umap1vx_list in dict_index_list.items():
                    Trimmed_Mean=calculate_trimmed_average(umap1vx_list)
                    addition_list.append(str(Trimmed_Mean))    
                addition_str='|'.join(addition_list)        
                f2.write(f"{target}\t{addition_str}\n")        
    if argv1=="stepall" or   argv1=="step10"  or   argv1=="step10.3s_a2": 
        print('Statistics, the proportion of all monomers matching the average')
        with open (f'./samples_satellite/10_monomer_model/monomer_umap1v_stat','w') as f3:
            with open (f'./samples_satellite/10_monomer_model/monomer_umap1v','r') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if eachline_arr[0]=='target':continue
                    target          =eachline_arr[0]
                    umap1vx_list    =eachline_arr[1].split('|')
                    #####
                    all_num=0;good_num=0
                    with open(f'./samples_satellite/10_monomer_model/{target}/1_umap1v','r') as f2:
                        for line in f2:
                            eachline_arr=line.strip().split('\t')
                            if eachline_arr[0]=='block':continue
                            all_num+=1
                            print(f"{target}\t{all_num}",end='\r')
                            one_umap1vx_list=eachline_arr[-1].split('|')
                            mark=''
                            for one1,one2 in zip(one_umap1vx_list, umap1vx_list): 
                                if abs(float(one1)-float(one2))>5:
                                    mark='bad'
                            if mark=='':good_num+=1
                    f3.write(f"{target}\t{good_num/all_num}\n"        )
    if argv1=="stepall" or   argv1=="step10"  or   argv1=="step10.3s_b1":
        print('10.3 Calculating the average may not be representative. Directly calculate the mode with divisions like 0-3, 1-4, 2-5, etc. Make three classifications, count numbers, 123s')
        def subunit_classA(umap1vx):
            umap1vx=float(umap1vx)
            if umap1vx<3:       A=1
            elif umap1vx<6:     A=2
            elif umap1vx<9:     A=3
            elif umap1vx<12:     A=4
            elif umap1vx<15:     A=5
            elif umap1vx<18:     A=6
            elif umap1vx<21:     A=7
            elif umap1vx<24:     A=8
            else:                A=9
            return A
        def subunit_classB(umap1vx):
            umap1vx=float(umap1vx)
            if umap1vx<4:       B=1
            elif umap1vx<7:     B=2
            elif umap1vx<10:     B=3
            elif umap1vx<13:     B=4
            elif umap1vx<16:     B=5
            elif umap1vx<19:     B=6
            elif umap1vx<22:     B=7
            elif umap1vx<25:     B=8
            else:                B=9
            return B
        def subunit_classC(umap1vx):
            umap1vx=float(umap1vx)
            if umap1vx<5:       C=1
            elif umap1vx<8:     C=2
            elif umap1vx<11:     C=3
            elif umap1vx<14:     C=4
            elif umap1vx<17:     C=5
            elif umap1vx<20:     C=6
            elif umap1vx<23:     C=7
            elif umap1vx<26:     C=8
            else:                C=9   
            return C
        def subunit_listclassA(umap1vx_list):
            A_str=''
            for one in umap1vx_list:
                A_str+=str(subunit_classA(one))
            return f"A{A_str}"    
        def subunit_listclassB(umap1vx_list):
            B_str=''
            for one in umap1vx_list:
                B_str+=str(subunit_classB(one))
            return f"B{B_str}" 
        def subunit_listclassC(umap1vx_list):
            C_str=''
            for one in umap1vx_list:
                C_str+=str(subunit_classC(one))
            return f"C{C_str}"     
        ##    
        target_list=[]
        target_list.append("28-23")
        target_list.append("28-28-23")
        target_list.append("28-28-28-23")
        target_list.append("28-28-28-28-23")
        target_list.append("28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-28-28-23")      
        def run_step(target):
            with open(f'./samples_satellite/10_monomer_model/{target}/1_umap1v_classABC','w') as f2:
                f2.write("block\tchr\tcirclen_arr\tstart\tend\tstrand\tseq_arr\tumap2v_arr\tdistance_arr\tarr_serial\tlength\tumap1vx_list\tclassA\tclassB\tclassC\n")
                with open(f'./samples_satellite/10_monomer_model/{target}/1_umap1v','r') as f:
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=line.split('\t')
                        if eachline_arr[0]=='block':continue
                        umap1vx_list=    eachline_arr[-1].split('|')[:-1]
                        classA  =subunit_listclassA(umap1vx_list)
                        classB  =subunit_listclassB(umap1vx_list)
                        classC  =subunit_listclassC(umap1vx_list)
                        newline=f"{eachline}\t{classA}\t{classB}\t{classC}\n"
                        f2.write(newline)
        with Pool(processes=thread) as pool:
            for i, result in enumerate(pool.imap(run_step, target_list), start=1):
                progress = (i / len(target_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush() 
    if argv1=="stepall" or   argv1=="step10"  or   argv1=="step10.3s_b2":
        print('10.3s_b2, 67s')
        target_list=[]
        target_list.append("28-23")
        target_list.append("28-28-23")
        target_list.append("28-28-28-23")
        target_list.append("28-28-28-28-23")
        target_list.append("28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-28-23")
        target_list.append("28-28-28-28-28-28-28-28-28-28-23")
        print()
        # For a list of strings, find the mode and its proportion
        from collections import Counter
        def find_mode_and_proportion(string_list):
            counter = Counter(string_list)
            mode = counter.most_common(1)[0][0]  # Mode string
            mode_count = counter.most_common(1)[0][1]  # Mode occurrence count
            total_count = len(string_list)
            proportion = mode_count / total_count
            return mode, proportion
            
        # For a list of strings, find the top five modes and their proportions
        def find_top_five_modes_and_proportions(string_list):
            counter = Counter(string_list)
            top_five = counter.most_common(5)  # Get the top five modes and their occurrence counts
            total_count = len(string_list)
            # Calculate the proportion for each mode, store as a list of (mode, proportion)
            result = [(mode, count / total_count) for mode, count in top_five]
            return result   
            
        with open(f'./samples_satellite/10_monomer_model/monomer_umap1vmode','w') as f2:    
            for target in target_list:
                print(target)
                listA=[]
                listB=[]
                listC=[]
                with open(f'./samples_satellite/10_monomer_model/{target}/1_umap1v_classABC','r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        listA.append(eachline_arr[-3])
                        listB.append(eachline_arr[-2])
                        listC.append(eachline_arr[-1])
                A_mode, A_proportion=find_mode_and_proportion(listA)
                B_mode, B_proportion=find_mode_and_proportion(listB)
                C_mode, C_proportion=find_mode_and_proportion(listC)
                f2.write(f"{target}\t{A_mode}\t{A_proportion}\t{B_mode}\t{B_proportion}\t{C_mode}\t{C_proportion}\n")    
        with open(f'./samples_satellite/10_monomer_model/monomer_umap1vmode_top5','w') as f2:    
            for target in target_list:
                print(target)
                listA=[]
                listB=[]
                listC=[]
                with open(f'./samples_satellite/10_monomer_model/{target}/1_umap1v_classABC','r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        listA.append(eachline_arr[-3])
                        listB.append(eachline_arr[-2])
                        listC.append(eachline_arr[-1])
                A_result_arr=find_top_five_modes_and_proportions(listA)
                B_result_arr=find_top_five_modes_and_proportions(listB)
                C_result_arr=find_top_five_modes_and_proportions(listC)
                min_length = min(len(A_result_arr), len(B_result_arr), len(C_result_arr))
                
                for i in range(min_length):
                    A_mode, A_proportion = A_result_arr[i]
                    B_mode, B_proportion = B_result_arr[i]
                    C_mode, C_proportion = C_result_arr[i]
                
                    f2.write(f"{target}\t{A_mode}\t{A_proportion}\t{B_mode}\t{B_proportion}\t{C_mode}\t{C_proportion}\n")    
                
    if argv1=="stepall" or   argv1=="step10"  or   argv1=="step10.5":
        print('Find typical examples of several typical monomer types, i.e., cases that are exactly duplicated, 5s')
        target_list=[]
        #target_list.append("28-23")
        #target_list.append("28-28-23")
        #target_list.append("28-28-28-23")
        target_list.append(["28-28-28-28-23","23","C7442","AB2CX"])
        target_list.append(["28-28-28-28-23","12","C8852","A2BCX"])
        
        target_list.append(["28-28-28-28-28-23","24|35","C85252","A(BC)2X"])
        target_list.append(["28-28-28-28-28-23","13|24","C86852","(AB)2CX"])
        target_list.append(["28-28-28-28-28-23","13|24","C85853","(AB)2CX"])
        
        target_list.append(["28-28-28-28-28-28-23","14|25|36","C853853","(ABC)2X"])
        target_list.append(["28-28-28-28-28-28-23","14|25|36","C852852","(ABC)2X"])
        target_list.append(["28-28-28-28-28-28-23","14|25|36","C863863","(ABC)2X"])
        target_list.append(["28-28-28-28-28-28-23","14|25|36","C853753","(ABC)2X"])
        
        target_list.append(["28-28-28-28-28-28-28-23","14|25|36","C8148142","(ACB)2CX",])
        target_list.append(["28-28-28-28-28-28-28-23","24|35|26|37","C8525252","A(BC)3X",])
        target_list.append(["28-28-28-28-28-28-28-23","24|35|26|37","C8535353","A(BC)3X",])
        target_list.append(["28-28-28-28-28-28-28-23","15|26|37|34","C8533853","ABC2ABCX",])
        
        target_list.append(["28-28-28-28-28-28-28-28-23","15|26|37|48|23","C74427442","(AB2C)2X"])
        target_list.append(["28-28-28-28-28-28-28-28-23","15|26|37|48|23","C74425442","(AB2C)2X"])
        target_list.append(["28-28-28-28-28-28-28-28-23","15|26|37|48|23","C71427142","(ACBC)2X"])
        target_list.append(["28-28-28-28-28-28-28-28-23","15|26|37|48|23","C75427542","(ABBC)2X"])
        
        target_list.append(["28-28-28-28-28-28-28-28-28-23","14|25|36|17|28|39","C853853853","(ABC)3X"])
        target_list.append(["28-28-28-28-28-28-28-28-28-23","14|25|36|17|28|39","C852852852","(ABC)3X"])
        target_list.append(["28-28-28-28-28-28-28-28-28-23","14|25|36|17|28|39","C953953953","(ABC)3X"])
        
        target_list.append(["28-28-28-28-28-28-28-28-28-28-23","14|25|36|17|28|39","C7147147142","(ACB)3CX"])  
        for one_target in target_list:
            one_target, pairs, C3, name= one_target
            pair_arr=pairs.split("|")
            input_file=f"./samples_satellite/10_monomer_model/{one_target}/1_umap1v_classABC"
            output_file=f"./samples_satellite/10_monomer_model/{one_target}/1_umap1v_classABC_example_{name}_{C3}"
            with open (output_file,'w') as f2:
                f2.write(f"block\tchr\tcirclen_arr\tstart\tend\tstrand\tseq_arr\tumap2v_arr\tdistance_arr\tarr_serial\tlength\tumap1vx_list\tclassA\tclassB\tclassC\n")
                with open (input_file,'r') as f:
                    next(f)
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        umap1vx_list=eachline_arr[11].split('|')
                        bad_mark=''
                        if eachline_arr[-1]!=C3:continue
                        for one_pair in pair_arr:
                            if umap1vx_list[int(one_pair[0])-1]!=umap1vx_list[int(one_pair[1])-1]:
                                bad_mark='yes'
                                break
                        if bad_mark=='':
                            f2.write(eachline+'\n')
            

if argv1=="stepall" or "step11" in argv1:
    if argv1=="stepall" or argv1=="step11_readme":
        print('Print instructions')
        with open('./samples_satellite/11_hor/readme','w') as f:
            txt=r'''
            11.0 ———— Keep all 23-separated ones, treat other2 from 28-separated as '.' and keep subunit format, add two columns: mark and mark_serial
            11.1 ———— Analyze HOR. shift=300
            11.2 ———— Perform a complex operation: transpose all results, sort by penalty in ascending order. If the lower one conflicts with the upper one (i.e., exceeds a certain range of the upper one), discard it
            11.3 ———— Determine HOR
            11.4 ———— Determine HOR within HOR,,, loop until no internal HOR remains
            11.5 ———— Summarize
            #11.6 ———— Remove parts with repeats less than 3 to get sum_repeat3........ Feels troublesome, not counting. Repeat count of 2 is actually correct, no need to require repeat count > 3
            '''
            f.write(txt)    
    if argv1=="stepall" or   argv1=="step11" or  argv1=="step11.0":
        print('step11.0 ———— Re-analyze HOR')
        if os.path.exists('./samples_satellite/11_hor')==False:
            subprocess.run(["mkdir ./samples_satellite/11_hor"], shell=True)
        if os.path.exists('./samples_satellite/11_hor/0_markserial_a'): 
            subprocess.run(["rm -r ./samples_satellite/11_hor/0_markserial_a"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/11_hor/0_markserial_a"], shell=True) 
        if os.path.exists('./samples_satellite/11_hor/0_markserial_b'): 
            subprocess.run(["rm -r ./samples_satellite/11_hor/0_markserial_b"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/11_hor/0_markserial_b"], shell=True)         
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)
               
        def run_step(one_region):
            input_file=f"./samples_satellite/9_circlenarr/sum_monomer2/{one_region}"    
            output_file=f"./samples_satellite/11_hor/0_markserial_a/{one_region}"   
            dict_markserial_circarr={}
            with open(output_file,'w') as f2:
                k=0
                serial=0
                with open(input_file) as f:
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        k+=1
                        if k==1:
                            newline=eachline+'\tmark\tmark_serilal'
                        else:
                            circ_serial=                eachline_arr[1]
                            subunit_revised_type=	    eachline_arr[-2]
                            subunit_revised_serial=     eachline_arr[-1]
                            subunit1=                   eachline_arr[-6]
                            serial1=                    eachline_arr[-5]
                            subunit2=                   eachline_arr[-4]
                            serial2=                    eachline_arr[-3]
                            ## Basically follow 20 groups. For other1 from the 23-interval classification method, take monomer; for other2 from the 28-segmentation method, ignore and treat as subunit format (non-monomer)
                            if subunit_revised_type=='.':
                                mark='.'
                                serial+=1
                            elif subunit_revised_type=='other2':             
                                mark='.'  
                                serial+=1
                            elif subunit_revised_type=='other1':
                                if serial1=='1':
                                    mark=subunit1
                                    serial+=1
                                else:    
                                    mark='delete'
                            elif subunit_revised_serial=='1':
                                mark=subunit_revised_type
                                serial+=1
                            else:
                                mark='delete'
                            ##
                            if serial not in dict_markserial_circarr:
                                dict_markserial_circarr[serial]=[]
                            dict_markserial_circarr[serial].append(circ_serial)    
    
                            newline=    eachline+'\t'+mark+'\t'+str(serial)
                        f2.write(newline+'\n')
            output_file=f"./samples_satellite/11_hor/0_markserial_b/{one_region}"               
            with open (output_file,'w') as f:
                f.write('markserial\tcirc_start\tcirc_end\tindex_start\tindex_end\n')
                for markserial,circarr in dict_markserial_circarr.items():
                    circ_start=circarr[0]
                    circ_end=circarr[-1]
                    index_start=int(circ_start[5:])-1
                    index_end=int(circ_end[5:])-1
                    f.write(f"{markserial}\t{circ_start}\t{circ_end}\t{index_start}\t{index_end}\n")
                            
        with Pool(processes=thread) as pool:
            for i, result in enumerate(pool.imap(run_step, region_list), start=1):
                progress = (i / len(region_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()     
    if argv1=="stepall" or   argv1=="step11" or  argv1=="step11.1":
        print('step11.1 ———— Analyze HOR. shift=300')
        if os.path.exists('./samples_satellite/11_hor/1_shiftscore')==False: 
            subprocess.run(["mkdir ./samples_satellite/11_hor/1_shiftscore"], shell=True)  
        
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)
        
        print('Main task started\n')     
        def run_step(one_region):
            ####
            input_file=f"./samples_satellite/11_hor/0_markserial_a/{one_region}"    
            list_x2=[];list_y2=[];list_seatlen=[]
            dict_serial_info={}
            with open(input_file,'r') as f:
                next(f)
                for line in f.readlines():
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    list_x2.append(eachline_arr[7])
                    list_y2.append(eachline_arr[8])
                    list_seatlen.append(eachline_arr[10])
            ##
            input_file=f"./samples_satellite/11_hor/0_markserial_b/{one_region}"    
            dict_markserial_info={}
            markserial_num=0
            with open(input_file,'r') as f:
                next(f)
                for line in f.readlines():
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    dict_markserial_info[int(eachline_arr[0])]=[int(eachline_arr[-2]),int(eachline_arr[-1])]       
                    markserial_num+=1
            #
            shift_max=300
            hor_start=1 
            dict_markserial_result={}
            while hor_start<markserial_num:
                ## Determine index_start1+shift/2shift based on shift(kk)
                kk=0;good_mark='';addtion_str_list=[]         ####20250606 I originally started with kk=1, meaning monomers of length 1 are not considered by default, but actually should be considered
                while kk<shift_max:
                    kk+=1   ##kk is the length of a monomer
                    hor_start_2=hor_start+kk;     
                    if hor_start_2+kk-1>markserial_num:break
                    ## Actual circ start and end
                    hor1_startserial_infos=     dict_markserial_info[hor_start]
                    hor1_endserial_infos=       dict_markserial_info[hor_start+kk-1]
                    hor2_startserial_infos=     dict_markserial_info[hor_start_2]
                    hor2_endserial_infos=       dict_markserial_info[hor_start_2+kk-1]
                    hor1_index_start=   hor1_startserial_infos[0]
                    hor1_index_end=     hor1_endserial_infos[1]
                    hor2_index_start=   hor2_startserial_infos[0]
                    hor2_index_end=     hor2_endserial_infos[1]       
                    hor1_circ_num=  hor1_index_end-hor1_index_start+1
                    hor2_circ_num=  hor2_index_end-hor2_index_start+1
                    if hor1_circ_num!=hor2_circ_num:continue
                    hor_circ_num=hor1_circ_num
                    ###
                    jj=0;penalty1=0;penalty2=0   ##
                    while jj<hor_circ_num:
                        pos1=float(list_x2[hor1_index_start+jj])
                        pos2=float(list_x2[hor2_index_start+jj])
                        abs1=abs(pos1-pos2)
                        if abs1>5 :break
                        penalty1+=abs1
                        ##
                        pos1=float(list_y2[hor1_index_start+jj])
                        pos2=float(list_y2[hor2_index_start+jj])
                        abs1=abs(pos1-pos2)
                        if abs1>5:break
                        penalty2+=abs1
                        ###
                        sum_penalty=penalty1+penalty2
                        if sum_penalty>hor_circ_num*2:break     #Allow a difference of 1 per circ
                        ##
                        seatlen1=int(list_seatlen[hor1_index_start+jj])
                        seatlen2=int(list_seatlen[hor2_index_start+jj])
                        seatlen_abs1=abs(seatlen1-seatlen2)
                        if seatlen_abs1>min(2,kk/5)   :break
                        ##
                        jj+=1
                    if jj==hor_circ_num:  #All position sequences processed
                        good_mark='yes';addtion_str_list.append(f"{kk}:{round((penalty1+penalty2)/(hor_circ_num),3)}")
                    ###########
                if  good_mark=='yes':
                    good_mark=''
                    addtion_str='|'.join(addtion_str_list)
                else:
                    addtion_str='.'
                dict_markserial_result[hor_start]=addtion_str
                hor_start+=1 
            input_file=f"./samples_satellite/11_hor/0_markserial_a/{one_region}"
            output_file=f"./samples_satellite/11_hor/1_shiftscore/{one_region}"           
            with open(output_file,'w') as f2:
                f2.write('region_name\tcirc_serial\tchrserial_1\tchrserial_2\tstrand\tcirc_seq\tumap1v_x\tumap2v_x\tumap2v_y\tin_seat_num\tcirc_len\tnum\tvariant_distance2onehap\tvariant_distance2allhaps\tsubunit1\tserial1\tsubunit2\tserial2\tsubunit_revised_type\tsubunit_revised_serial\tmark\tmark_serilal\tHOR_mer_raw|penalty_average\n')        
                with open (input_file,'r')  as f:
                    next(f)
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        mark_serilal=eachline_arr[-1]
                        if int(mark_serilal) not in dict_markserial_result:
                            result='.'
                        else:    
                            result=dict_markserial_result[int(mark_serilal)]  
                        newline=f"{eachline}\t{result}\n"
                        f2.write(newline)
        with Pool(processes=thread) as pool:
            for i, result in enumerate(pool.imap(run_step, region_list), start=1):
                progress = (i / len(region_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()                 
    if argv1=="stepall" or   argv1=="step11" or   argv1=="step11.2":    
        print('step11.2 ———— Perform a complex operation: transpose all results, sort by penalty in ascending order. If the lower one conflicts with the upper one (i.e., exceeds a certain range of the upper one), discard it, 286s')
        if os.path.exists('./samples_satellite/11_hor/2_cascade')==True:     
            subprocess.run(["rm -r ./samples_satellite/11_hor/2_cascade"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/11_hor/2_cascade"], shell=True)     
    
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)    
                
        def run_step(one_region):
            input_file=f"./samples_satellite/11_hor/1_shiftscore/{one_region}"  
            output_file=f"./samples_satellite/11_hor/2_cascade/{one_region}"  
            ###
            info_arr=[]
            with open(input_file,'r') as f:
                next(f)
                for line in f: 
                    eachline=line.strip() 
                    eachline_arr=eachline.split('\t')
                    mark=               eachline_arr[-3]
                    if mark=='delete':continue
                    mark_serilal=       int(eachline_arr[-2])
                    mer_penalty=        eachline_arr[-1]
                    mer_penalty_arr=    mer_penalty.split('|')                        
                    if mer_penalty!='.':
                        for one in mer_penalty_arr:
                            one_arr=one.split(':')
                            one_mer=    int(one_arr[0])
                            one_penalty=float(one_arr[1])
                            end_serial=    mark_serilal+one_mer*2 -1  #20250607 Originally it was 3
                            info_arr.append([one_mer,one_penalty,mark_serilal,end_serial])
            # Sort by the second element of each sublist in ascending order
            sorted_info_arr = sorted(info_arr, key=lambda x: x[1])
            ###
            with open(output_file,'w') as f2:
                f2.write('mer\tpenalty\tstart_serial\tend_serial\n')
                for one in sorted_info_arr:
                    f2.write(f'{one[0]}\t{one[1]:.3f}\t{one[2]}\t{one[3]}\n')
            
        with Pool(processes=thread) as pool:
            for i, result in enumerate(pool.imap(run_step, region_list), start=1):
                progress = (i / len(region_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()                
    if argv1=="stepall" or   argv1=="step11" or   argv1=="step11.3":
        print('step11.3 ———— Determine HOR')
        if os.path.exists('./samples_satellite/11_hor/3_cover')==True:     
            subprocess.run(["rm -r ./samples_satellite/11_hor/3_cover"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/11_hor/3_cover"], shell=True)     
    
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)    
                
        def run_step(one_region):
            input_file=f"./samples_satellite/11_hor/2_cascade/{one_region}"  
            ###
            dict_markserial_mer={}      ### Good results
            dict_markserial_penalty={}
            i=0
            with open(input_file,'r') as f:
                next(f)
                for line in f: 
                    i+=1
                    eachline=line.strip() 
                    eachline_arr=eachline.split('\t')
                    mer=	        int(eachline_arr[0])
                    penalty=	    eachline_arr[1]
                    start_serial=	int(eachline_arr[2])
                    end_serial=	    int(eachline_arr[3])
                    if i==1:
                        k=start_serial
                        while k<=end_serial:
                            if k not in dict_markserial_mer:
                                dict_markserial_mer[k]=mer
                                dict_markserial_penalty[k]=penalty
                            k+=1
                    if    end_serial not  in dict_markserial_mer and start_serial not in dict_markserial_mer:     
                        k=start_serial;bad_mark=''
                        while k<=end_serial:
                            if k in dict_markserial_mer:
                                bad_mark='yes';break
                            k+=1 
                        if bad_mark=='yes':continue                        
                        k=start_serial
                        while k<=end_serial:
                            if k not in dict_markserial_mer:
                                dict_markserial_mer[k]=mer
                                dict_markserial_penalty[k]=penalty
                            k+=1
                    if  (start_serial in dict_markserial_mer and end_serial not in dict_markserial_mer) or (end_serial in dict_markserial_mer and start_serial not in dict_markserial_mer):
                        k=start_serial;bad_mark='';k_arr=[]
                        while k<=end_serial:
                            if k in dict_markserial_mer:
                                if dict_markserial_mer[k]!=mer: bad_mark='yes';break
                            else:
                                k_arr.append(k)
                            k+=1 
                        if bad_mark=='yes':continue
                        for k in k_arr:
                            dict_markserial_mer[k]=mer
                            dict_markserial_penalty[k]=penalty
                        
            sorted_keys = sorted(dict_markserial_mer.keys())
            output_file=f"./samples_satellite/11_hor/3_cover/{one_region}"  
            
            mer_last='';markserial_last=''
            with open(output_file,'w') as f2:  
                f2.write('Hor_name\tmark_serilal\tmer\tpenalty\n')
                for markserial in sorted_keys:
                    mer=      dict_markserial_mer[markserial]
                    penalty=    dict_markserial_penalty[markserial]
                    if mer_last=='':
                        serial=1;
                    elif  mer_last!=mer or markserial!=markserial_last+1:
                        serial+=1
                    mer_last=mer;markserial_last=markserial
                    f2.write(f"Hor-{serial}\t{markserial}\t{mer}\t{penalty}\n")    
 
        with Pool(processes=thread) as pool:
            for i, result in enumerate(pool.imap(run_step, region_list), start=1):
                progress = (i / len(region_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()   
    if argv1=="stepall" or   argv1=="step11" or   argv1=="step11.4":
        print('step11.4 ———— Determine HOR within HOR, 200s')
        last_dir,current_dir=['3_cover','4_inner']
        run_list=['3_cover','4_inner1','4_inner2','4_inner3','4_inner4','4_inner5','4_inner6','4_inner7','4_inner8','4_inner9','4_inner10','4_inner11','4_inner12','4_inner13','4_inner14','4_inner15']
        run_num=len(run_list)
        iii=0
        with open('./samples_satellite/11_hor/4_inner.log','w') as fiii:
            while iii<run_num-1:
                last_dir,current_dir=run_list[iii],run_list[iii+1]
                iii+=1
                if os.path.exists(f'./samples_satellite/11_hor/{current_dir}')==True:     
                    subprocess.run([f"rm -r ./samples_satellite/11_hor/{current_dir}"], shell=True)   
                subprocess.run([f"mkdir ./samples_satellite/11_hor/{current_dir}"], shell=True)     
            
                print('Loading regions')
                region_list=[]
                with open('./samples_satellite/2_good_regions','r') as f:
                    for line in f.readlines():
                        eachline_arr=line.strip().split('\t')
                        if eachline_arr[0]=='sample':continue
                        if len(eachline_arr)!=8:continue
                        new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                        region_list.append(new_name)                            
                    
                def run_step(one_region):
                    input_file=f"./samples_satellite/11_hor/{last_dir}/{one_region}"  
                    dict_markserial_horname_used={}
                    dict_markserial_mer_used={}
                    with open(input_file,'r') as f:
                        next(f)
                        for line in f: 
                            eachline=line.strip() 
                            eachline_arr=eachline.split('\t')
                            hor_name=       eachline_arr[0]
                            mark_serial=    int(eachline_arr[1])
                            HOR_mer=        int(eachline_arr[2])
                            dict_markserial_horname_used[mark_serial]=hor_name
                            dict_markserial_mer_used[mark_serial]=HOR_mer
                            
        
                    input_file=f"./samples_satellite/11_hor/2_cascade/{one_region}"  
                    info_arr=[]
                    dict_markserial_mer={}
                    dict_markserial_penalty={}       
                    i=0
                    with open(input_file,'r') as f:
                        next(f)
                        for line in f: 
                            
                            eachline=line.strip() 
                            eachline_arr=eachline.split('\t')
                            mer=	        int(eachline_arr[0])
                            penalty=	    eachline_arr[1]
                            start_serial=	int(eachline_arr[2])
                            end_serial=	    int(eachline_arr[3])              
                            ## Start and end are within the same HOR
                            if start_serial in dict_markserial_horname_used  and end_serial in dict_markserial_horname_used  \
                                and  dict_markserial_horname_used [start_serial]==dict_markserial_horname_used [end_serial]:
                                ### The child HOR's mer is less than the parent HOR's mer
                                if mer<=dict_markserial_mer_used[start_serial]/2:  ######20250607, originally it was 3
                                    i+=1
                                    if i==1:
                                        k=start_serial
                                        while k<=end_serial:
                                            if k not in dict_markserial_mer:
                                                dict_markserial_mer[k]=mer
                                                dict_markserial_penalty[k]=penalty
                                            k+=1
                                    if    end_serial not  in dict_markserial_mer and start_serial not in dict_markserial_mer:     
                                        k=start_serial;bad_mark=''
                                        while k<=end_serial:
                                            if k in dict_markserial_mer:
                                                bad_mark='yes';break
                                            k+=1 
                                        if bad_mark=='yes':continue                        
                                        k=start_serial
                                        while k<=end_serial:
                                            if k not in dict_markserial_mer:
                                                dict_markserial_mer[k]=mer
                                                dict_markserial_penalty[k]=penalty
                                            k+=1
                                    if  (start_serial in dict_markserial_mer and end_serial not in dict_markserial_mer) or (end_serial in dict_markserial_mer and start_serial not in dict_markserial_mer):
                                        k=start_serial;bad_mark='';k_arr=[]
                                        while k<=end_serial:
                                            if k in dict_markserial_mer:
                                                if dict_markserial_mer[k]!=mer: bad_mark='yes';break
                                            else:
                                                k_arr.append(k)
                                            k+=1 
                                        if bad_mark=='yes':continue
                                        for k in k_arr:
                                            dict_markserial_mer[k]=mer
                                            dict_markserial_penalty[k]=penalty
                    sorted_keys = sorted(dict_markserial_mer.keys())
                    output_file=f"./samples_satellite/11_hor/{current_dir}/{one_region}"  
                    
                    mer_last='';markserial_last=''
                    with open(output_file,'w') as f2:  
                        f2.write('Hor_name\tmark_serilal\tmer\tpenalty\n')
                        for markserial in sorted_keys:
                            mer=      dict_markserial_mer[markserial]
                            penalty=    dict_markserial_penalty[markserial]
                            if mer_last=='':
                                serial=1;
                            elif  mer_last!=mer or markserial!=markserial_last+1:
                                serial+=1
                            mer_last=mer;markserial_last=markserial
                            f2.write(f"Hor{'Hor'*iii}-{serial}\t{markserial}\t{mer}\t{penalty}\n")  
                        
                with Pool(processes=thread) as pool:
                    for i, result in enumerate(pool.imap(run_step, region_list), start=1):
                        progress = (i / len(region_list)) * 100
                        sys.stdout.write(f"\rProgress: {progress:.2f}%")
                        sys.stdout.flush()     
                ## Statistics
                if iii==1:
                    all_num=0
                    for one_region in region_list:
                        with open(f"./samples_satellite/11_hor/{last_dir}/{one_region}",'r') as f:
                            next(f)  
                            for line in f:
                                eachline=line.strip()
                                if len(eachline)>1:
                                    all_num+=1    
                    print('Total results count: '+str(all_num))
                    fiii.write(f'{last_dir}\t{all_num}\n')
                ###    
                all_num=0
                for one_region in region_list:
                    with open(f"./samples_satellite/11_hor/{current_dir}/{one_region}",'r') as f:
                        next(f)  
                        for line in f:
                            eachline=line.strip()
                            if len(eachline)>1:
                                all_num+=1    
                print('Total results count: '+str(all_num))
                fiii.write(f'{current_dir}\t{all_num}\n')
                if all_num==0:
                    break
    if argv1=="stepall" or   argv1=="step11" or   argv1=="step11.5":    
        print('step11.5 ———— Summarize, 33s')
        if os.path.exists('./samples_satellite/11_hor/5_sum')==True:     
            subprocess.run(["rm -r ./samples_satellite/11_hor/5_sum"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/11_hor/5_sum"], shell=True)     
    
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)                            
            
        def run_step(one_region):
            input_file=f"./samples_satellite/11_hor/3_cover/{one_region}"  
            dict_layer_markserial_info={}
            dict_layer_markserial_info[1]={}
            with open(input_file,'r') as f:
                next(f)
                for line in f: 
                    eachline=line.strip() 
                    eachline_arr=eachline.split('\t')
                    hor_name=	        eachline_arr[0]  
                    mark_serilal=       eachline_arr[1]
                    mer	=               eachline_arr[2]
                    penalty	=           eachline_arr[3]
                    dict_layer_markserial_info[1][mark_serilal]=[hor_name,mer,penalty]
            ####    
            inner_result_dir=[]
            with open(f"./samples_satellite/11_hor/4_inner.log",'r') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if eachline_arr[0]=='3_cover':continue
                    if int(eachline_arr[1])>0:
                        inner_result_dir.append(eachline_arr[0])
            
            layer=1
            for one_innerresult_dir in inner_result_dir:
                layer+=1
                dict_layer_markserial_info[layer]={}
                input_file=f"./samples_satellite/11_hor/{one_innerresult_dir}/{one_region}"  
                
                with open(input_file,'r') as f:
                    next(f)
                    for line in f: 
                        eachline=line.strip() 
                        eachline_arr=eachline.split('\t')
                        hor_name=	        eachline_arr[0]  
                        mark_serilal=       eachline_arr[1]
                        mer	=               eachline_arr[2]
                        penalty	=           eachline_arr[3]
                        dict_layer_markserial_info[layer][mark_serilal]=[hor_name,mer,penalty]
            intput_file=f"./samples_satellite/11_hor/1_shiftscore/{one_region}"
            output_file=f"./samples_satellite/11_hor/5_sum/{one_region}" 
            with open(output_file,'w') as f2:
                head_addition=''
                iii=0
                while iii<layer:
                    iii+=1
                    head_addition+=f"\thor_{iii}\tmer_{iii}\tpenalty_{iii}"
                f2.write(f"region_name\tcirc_serial\tchrserial_1\tchrserial_2\tstrand\tcirc_seq\tumap1v_x\tumap2v_x\tumap2v_y\tin_seat_num\tcirc_len\tnum\tvariant_distance2onehap\tvariant_distance2allhaps\tsubunit1\tserial1\tsubunit2\tserial2\tsubunit_revised_type\tsubunit_revised_serial\tmark\tmark_serilal\tHOR_mer_raw|penalty_average{head_addition}\n")
                with open(intput_file,'r') as f:
                    next(f)
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        mark_serilal=eachline_arr[-2]
                        iii=0
                        addtion_str=''
                        while iii<layer:
                            iii+=1  
                            dict_markserial_info=dict_layer_markserial_info[iii]
                            if mark_serilal in dict_markserial_info:
                                hor_name,mer,penalty=dict_markserial_info[mark_serilal]
                            else:
                                hor_name,mer,penalty=['.','.','.']
                            addtion_str+=f"\t{hor_name}\t{mer}\t{penalty}"
                        newline=eachline+addtion_str
                        f2.write(f'{newline}\n')

                    
        with Pool(processes=thread) as pool:
            for i, result in enumerate(pool.imap(run_step, region_list), start=1):
                progress = (i / len(region_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()
       

## Grid statistics for monomers and HORs, continue from step5
if argv1=="stepall" or "step5s" in argv1:    
    if argv1=="stepall" or argv1=="step5s"  or argv1=="step5s_step5.2":       
        print('Statistical analysis of different intervals')
        if  os.path.exists(f'./samples_satellite/5_umap_grid/monomer_stat')==False:
            subprocess.run([f"mkdir ./samples_satellite/5_umap_grid/monomer_stat"], shell=True) 
        subprocess.run([f"mkdir ./samples_satellite/5_umap_grid/monomer_stat/9-monomer_density"], shell=True)
        subprocess.run([f"mkdir ./samples_satellite/5_umap_grid/monomer_stat/10-HORhor_density"], shell=True)
        
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)
        def run_step5_2(one_region):
            input_file=f"./samples_satellite/5_umap_grid/monomer/{one_region}"
            output_file9=f"./samples_satellite/5_umap_grid/monomer_stat/9-monomer_density/{one_region}"
            output_file10=f"./samples_satellite/5_umap_grid/monomer_stat/10-HORhor_density/{one_region}"
            one_region_arr=one_region.split(':')
            sample=     one_region_arr[0]
            chrname=    one_region_arr[1].replace('region_','Chr')[:-1]
            block_pos_arr= one_region_arr[2].split('-')
            block_size=int(block_pos_arr[1])-int(block_pos_arr[0])
            i=0;order=''
            dict_circserial_pos={}
            with open(input_file,'r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    circserial=                 eachline_arr[1]
                    umap2v_x_grid=	            eachline_arr[-2]
                    umap2v_y_grid=              eachline_arr[-1]
                    pos=umap2v_x_grid+'_'+umap2v_y_grid
                    dict_circserial_pos[circserial]=pos
            dict_pos_infos={}
            input_file9=f"./samples_satellite/11_hor/5_sum/{one_region}"
            with open(input_file9,'r') as f:
                next(f)
                last=''
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    circserial=                 eachline_arr[1]   
                    mark=                       eachline_arr[20]   
                    ##
                    if last=='':
                        last=mark
                    elif  mark=='delete':   mark=last
                    else: last=mark
                    ###
                    pos=dict_circserial_pos[circserial]
                    monomer=mark
                    if pos not in dict_pos_infos:
                        dict_pos_infos[pos]={}
                    if 'monomer' not in dict_pos_infos[pos]:
                        dict_pos_infos[pos]['monomer']=[]
                    dict_pos_infos[pos]['monomer'].append(monomer)
            #
            input_file10=f"./samples_satellite/11_hor/5_sum/{one_region}"
            with open(input_file10,'r') as f:
                next(f)
                last=''
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    circserial=                 eachline_arr[1]   
                    ##
                    mer1=                 eachline_arr[24]   
                    mer2=                 eachline_arr[27]   
                    mer3=                 eachline_arr[30]   
                    mer4=                 eachline_arr[33]   
                    mer5=                 eachline_arr[36]   
                    mer6=                 eachline_arr[39]   
                    HORhor_list=[]
                    if mer1!='.' and mer1!='1':     HORhor_list.append(mer1)
                    if mer2!='.' and mer2!='1':     HORhor_list.append(mer2) 
                    if mer3!='.' and mer3!='1':     HORhor_list.append(mer3) 
                    if mer4!='.' and mer4!='1':     HORhor_list.append(mer4) 
                    if mer5!='.' and mer5!='1':     HORhor_list.append(mer5) 
                    if mer6!='.' and mer6!='1':     HORhor_list.append(mer6) 
                
                    if len(HORhor_list)==0: HORhor_list=['noHor']
                    ###
                    pos=dict_circserial_pos[circserial]
                    if pos not in dict_pos_infos:
                        dict_pos_infos[pos]={}
                    if 'HORhor' not in  dict_pos_infos[pos]:   
                        dict_pos_infos[pos]['HORhor']=[]
                    for one in HORhor_list:
                        dict_pos_infos[pos]['HORhor'].append(one)         
            ####            
            with open(output_file9,'w')  as f9:
                with open(output_file10,'w')  as f10:
                    for pos,infos in dict_pos_infos.items():
                        one_monomer_list=         infos['monomer']
                        one_HORhor_list=          infos['HORhor']
                        #
                        for one_monomer in one_monomer_list:
                            f9.write(f"{pos}\t{one_monomer}\n")
                        #
                        for one_HORhor in one_HORhor_list:
                            f10.write(f"{pos}\t{one_HORhor}\n")
        with Pool(processes=thread) as pool:
            for i, result in enumerate(pool.imap(run_step5_2, region_list), start=1):
                progress = (i / len(region_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()      
    
    if argv1=="stepall" or argv1=="step5s"  or argv1=="step5s_step5.2s":    
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)
        ### Find mode and percentage of the list
        from collections import Counter
        def calculate_mode_and_percentage(num_list):
            counts = Counter(num_list)
            max_count = max(counts.values())
            mode = next(num for num, count in counts.items() if count == max_count)
            num_list_num=len(num_list)
            mode_percentage = round((max_count / num_list_num),3)
            return mode, num_list_num, mode_percentage
        ##############################################    
        print('\nStatistics 9-monomer_density')
        print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
        if  os.path.exists(f'./samples_satellite/5_umap_grid/monomer_stat/9-stat_monomer_density')==False:
            dict_pos_monomerlist={}
            for one_region in region_list:
                input_file=f"./samples_satellite/5_umap_grid/monomer_stat/9-monomer_density/{one_region}"
                with open(input_file,'r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        one_pos=    eachline_arr[0]
                        monomer=   eachline_arr[1]
                        if monomer=="28-28-28-23":monomer='CEN107'
                        elif monomer=="28-28-28-28-23":monomer='CEN135'
                        elif monomer=="28-23":monomer='other28n-23'
                        elif monomer=="28-28-23":monomer='other28n-23'
                        elif monomer=="28-28-28-28-28-23":monomer='other28n-23'
                        elif monomer=="28-28-28-28-28-28-23":monomer='other28n-23'
                        elif monomer=="28-28-28-28-28-28-28-23":monomer='other28n-23'
                        elif monomer=="28-28-28-28-28-28-28-28-23":monomer='other28n-23'
                        elif monomer=="28-28-28-28-28-28-28-28-28-23":monomer='other28n-23'
                        elif monomer=="28-28-28-28-28-28-28-28-28-28-23":monomer='other28n-23'
                        elif monomer=="28-28-28-28-28-28-28-28-28-28-28-23":monomer='other28n-23'
                        else: monomer='Other'
                        
                        if one_pos not in dict_pos_monomerlist:   dict_pos_monomerlist[one_pos]=[]
                        dict_pos_monomerlist[one_pos].append(monomer)
            with open('./samples_satellite/5_umap_grid/monomer_stat/9-stat_monomer_density','w') as f :
                f.write(f"pos\tpos_x\tpos_y\tmonomer_mode\tsum_num\tpercent\n")    
                for one_pos,monomerlist in dict_pos_monomerlist.items():
                    arr=one_pos.split('_')
                    mode, sum_num, percent=calculate_mode_and_percentage(monomerlist)
                    f.write(f"{one_pos}\t{arr[0]}\t{arr[1]}\t{mode}\t{sum_num}\t{percent}\n") 
        ##############################################    
        print('\nStatistics 10-HORhor_density')
        print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
        if  os.path.exists(f'./samples_satellite/5_umap_grid/monomer_stat/10-stat_HORhor_density')==False:
            dict_pos_HORhorlist={}
            for one_region in region_list:
                input_file=f"./samples_satellite/5_umap_grid/monomer_stat/10-HORhor_density/{one_region}"
                with open(input_file,'r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        one_pos=    eachline_arr[0]
                        HORhor=   eachline_arr[1]
                        if HORhor=='noHor':HORhor='noHOR'
                        elif int(HORhor)>=2 and int(HORhor)<=4 :HORhor='mer2-4'
                        elif int(HORhor)>=5 and int(HORhor)>=10:HORhor='mer5-10'
                        else:HORhor='mer10plus'
                        if one_pos not in dict_pos_HORhorlist:   dict_pos_HORhorlist[one_pos]=[]
                        dict_pos_HORhorlist[one_pos].append(HORhor)
            with open('./samples_satellite/5_umap_grid/monomer_stat/10-stat_HORhor_density','w') as f :
                f.write(f"pos\tpos_x\tpos_y\tHORhor_mode\tsum_num\tpercent\n")    
                for one_pos,HORhorlist in dict_pos_HORhorlist.items():
                    arr=one_pos.split('_')
                    mode, sum_num, percent=calculate_mode_and_percentage(HORhorlist)
                    f.write(f"{one_pos}\t{arr[0]}\t{arr[1]}\t{mode}\t{sum_num}\t{percent}\n")                     
                 
    if argv1=="stepall" or argv1=="step5s"  or argv1=="step5s_step5.2p":             
        R_txt=r'''
library(ggplot2)
library(dplyr)
### Monomer
print("9-stat_monomer_density")
{
  input_file9 <- read.table('9-stat_monomer_density', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  filtered_input_file9 <- input_file9 %>%  filter(sum_num > 100)
  # Create categorical column
  filtered_input_file9 <- filtered_input_file9 %>%
    mutate(
      category = case_when(
        percent < 0.5 ~ "Mixed",
        monomer_mode == "CEN107" ~ "CEN107",
        monomer_mode == "CEN135" ~ "CEN135",
        monomer_mode == "other28n-23" ~ "28n-23",
        monomer_mode == "Other" ~ "Other",
        TRUE ~ "Other"
      ),
    )
  
  # Use ggplot2 to draw scatter plot
  p <- ggplot(filtered_input_file9, aes(x = pos_x, y = pos_y, color = category)) +
    geom_point(size = 0.1, shape = 15) +
    theme_bw() +
    coord_equal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    ) +
    scale_color_manual(name = "monomer Mode", 
                       values = c("Mixed" = "grey80", 
                                  "CEN107" = "#336699",  
                                  "CEN135" = "#8CC63F",  
                                  "28n-23" = "#CC0066",  
                                  "Other" = "#7f7f7f"))  +
    guides(color = guide_legend(override.aes = list(size = 4), title.position = "top", title.hjust = 0.5,keywidth  = unit(0.1, "cm"),keyheight = unit(0.1, "cm")))
  
  
  # Save as PDF
  pdf("pic_9-stat_monomer_density2.pdf", width = 14 / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
  print(p)
}

### HORhor
print("10-stat_HORhor_density")
{
  input_file10 <- read.table('10-stat_HORhor_density', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  filtered_input_file10 <- input_file10 %>%  filter(sum_num > 100)
  # Create categorical column
  filtered_input_file10 <- filtered_input_file10 %>%
    mutate(
      category = case_when(
        percent < 0.5 ~ "Mixed",
        HORhor_mode == "noHOR" ~ "noHOR",
        HORhor_mode == "mer2-4" ~ "mer2-4",
        HORhor_mode == "mer5-10" ~ "mer5-10",
        HORhor_mode =="mer10plus" ~"mer10plus",
        TRUE ~ "Other"
      ),
    )
  
  # Use ggplot2 to draw scatter plot
  p <- ggplot(filtered_input_file10, aes(x = pos_x, y = pos_y, color = category)) +
    geom_point(size = 0.1, shape = 15) +
    theme_bw() +
    coord_equal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    ) +
    scale_color_manual(name = "HORhor Mode", 
                       values = c("Mixed" = "grey80", 
                                  "noHOR" = "black",  
                                  "mer2-4" = "#3366cc",  
                                  "mer5-10" = "#00b386",  
                                  "mer10plus" = "#ff5050",
                                  "Other" = "#ff5050"))  +
    guides(color = guide_legend(override.aes = list(size = 1), title.position = "top", title.hjust = 0.5,keywidth  = unit(0.1, "cm"),keyheight = unit(0.1, "cm")))
  
  
  # Save as PDF
  pdf("pic_10-stat_HORhor_density.pdf", width = 14 / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
  print(p)
}
    '''
        new_directory = "./samples_satellite/5_umap_grid/monomer_stat/"
        os.chdir(new_directory)
        with open('./pic_ggplot_step5_2.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        subprocess.run(['Rscript pic_ggplot_step5_2.R'], shell=True)    
        os.chdir('../../../../')

if argv1=="stepall" or "step12" in argv1:    
    print ("step12          Statistics for panel b/c of pic3, regarding HOR layers and whether multiple monomer types are included, filtering out HORs with repeat count less than 3")  
    if  os.path.exists(f'./samples_satellite/12_Hor_stat')==False:
        subprocess.run([f"mkdir ./samples_satellite/12_Hor_stat"], shell=True) 
    if argv1=="stepall" or argv1=="step12_readme":
        print('Print instructions')
        with open('./samples_satellite/12_Hor_stat/readme','w') as f:
            txt=r'''
            12.1 ———— 1_blocks_stat folder, analyze each block separately, includes MT/ST column, layer / father_layer, HOR repeat count
            12.2 ———— Length of different mers, Figure 3a
            12.3 ———— Figure 3d
            12.4 ———— Table, mer lengths of different segments
            12.5 ———— Figure 3c, lengths of different layers
            12.6 ———— Find the longest HOR, see how many mers it can extend
            '''
            f.write(txt) 
    if argv1=="stepall" or argv1=="step12"  or argv1=="step12.1":
        print('Statistical analysis of different mers')
        ##28-28-28-23|28-28-23|28-28-28-28-23 HOR monomer arrangement, but it is equivalent to 28-28-28-28-23|28-28-28-23|28-28-23, i.e., rotation, translocation. How to obtain a unique index
        def normalize_hor_pattern_rotation_only(sequence, separator='|'):
            """
            Normalize HOR monomer sequence, only handle rotation, not considering reversal
            
            For example: "28-28-28-23|28-28-23|28-28-28-28-23"
               -> normalized to unique string
            
            Strategy: Find the smallest lexicographical order among all rotations
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
            
            # 3. Return the smallest lexicographical string
            return min(candidates) 
            
        if os.path.exists('./samples_satellite/12_Hor_stat/1_blocks_stat')==True:     
            subprocess.run(["rm -r ./samples_satellite/12_Hor_stat/1_blocks_stat"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/12_Hor_stat/1_blocks_stat"], shell=True)    
        
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)
                
        def run_step(one_region):
            input_file=f"./samples_satellite/11_hor/5_sum/{one_region}" 
            output_file=f"./samples_satellite/12_Hor_stat/1_blocks_stat/{one_region}" 
            dict_Hor_info={}
            with open(input_file,'r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t')    
                        hor1=               eachline_arr[23]     
                        mer1=               eachline_arr[24]   
                        penalty1=	        eachline_arr[25]  
                        hor2=               eachline_arr[26]     
                        mer2=               eachline_arr[27]   
                        penalty2=	        eachline_arr[28]
                        hor3=               eachline_arr[29]     
                        mer3=               eachline_arr[30]   
                        penalty3=	        eachline_arr[31]
                        hor4=               eachline_arr[32]     
                        mer4=               eachline_arr[33]   
                        penalty4=	        eachline_arr[34]
                        hor5=               eachline_arr[35]     
                        mer5=               eachline_arr[36]   
                        penalty5=	        eachline_arr[37]
                        hor6=               eachline_arr[38]     
                        mer6=               eachline_arr[39]   
                        penalty6=	        eachline_arr[40]   
                        ###
                        circ_serial	=       eachline_arr[1] 
                        chrserial_1=	    eachline_arr[2]    
                        chrserial_2=        eachline_arr[3] 
                        circ_len=           eachline_arr[10]
                        
                        mark=               eachline_arr[20]    
                        markserial=         eachline_arr[21]    
                        if mer1 != '.' and mer1!='1':  
                            if hor1 not in dict_Hor_info:       
                                dict_Hor_info[hor1]={}
                                dict_Hor_info[hor1]['mer']=mer1
                                dict_Hor_info[hor1]['segment-end']=int(markserial)+int(mer1)-1 
                                dict_Hor_info[hor1]['segment']=mark
                                dict_Hor_info[hor1]['segment_length']=int(circ_len)
                                dict_Hor_info[hor1]['circ-start']=circ_serial
                                dict_Hor_info[hor1]['circ-end']=circ_serial
                                dict_Hor_info[hor1]['position-start']=chrserial_1
                                dict_Hor_info[hor1]['position-end']=chrserial_2
                                dict_Hor_info[hor1]['markserial-start']=markserial
                                dict_Hor_info[hor1]['markserial-end']=markserial                                
                                dict_Hor_info[hor1]['father_layer']=1
                                dict_Hor_info[hor1]['layer']=1
                            else: 
                                if mark!='delete' and dict_Hor_info[hor1]['segment-end']>=int(markserial):
                                    dict_Hor_info[hor1]['segment']+=f"|{mark}"
                                if dict_Hor_info[hor1]['segment-end']>=int(markserial):   
                                    dict_Hor_info[hor1]['segment_length']+=int(circ_len)
                            dict_Hor_info[hor1]['circ-end']=circ_serial
                            dict_Hor_info[hor1]['position-end']=chrserial_2
                            dict_Hor_info[hor1]['markserial-end']=markserial 
                            if    dict_Hor_info[hor1]['father_layer']<6  and   mer6 != '.' and mer6!='1':dict_Hor_info[hor1]['father_layer']=6
                            elif  dict_Hor_info[hor1]['father_layer']<5  and   mer5 != '.' and mer5!='1':dict_Hor_info[hor1]['father_layer']=5
                            elif  dict_Hor_info[hor1]['father_layer']<4  and   mer4 != '.' and mer4!='1':dict_Hor_info[hor1]['father_layer']=4
                            elif  dict_Hor_info[hor1]['father_layer']<3  and   mer3 != '.' and mer3!='1':dict_Hor_info[hor1]['father_layer']=3
                            elif  dict_Hor_info[hor1]['father_layer']<2  and   mer2 != '.' and mer2!='1':dict_Hor_info[hor1]['father_layer']=2
                        if mer2 != '.' and mer2!='1':  
                            if hor2 not in dict_Hor_info:       
                                dict_Hor_info[hor2]={}
                                dict_Hor_info[hor2]['mer']=mer2
                                dict_Hor_info[hor2]['segment-end']=int(markserial)+int(mer2)-1 
                                dict_Hor_info[hor2]['segment']=mark
                                dict_Hor_info[hor2]['segment_length']=int(circ_len)                                
                                dict_Hor_info[hor2]['circ-start']=circ_serial
                                dict_Hor_info[hor2]['circ-end']=circ_serial
                                dict_Hor_info[hor2]['position-start']=chrserial_1
                                dict_Hor_info[hor2]['position-end']=chrserial_2
                                dict_Hor_info[hor2]['markserial-start']=markserial
                                dict_Hor_info[hor2]['markserial-end']=markserial                                    
                                dict_Hor_info[hor2]['father_layer']=1
                                dict_Hor_info[hor2]['layer']=2
                            else: 
                                if mark!='delete' and dict_Hor_info[hor2]['segment-end']>=int(markserial):
                                    dict_Hor_info[hor2]['segment']+=f"|{mark}"
                                if dict_Hor_info[hor2]['segment-end']>=int(markserial):   
                                    dict_Hor_info[hor2]['segment_length']+=int(circ_len)                                    
                            dict_Hor_info[hor2]['circ-end']=circ_serial
                            dict_Hor_info[hor2]['position-end']=chrserial_2
                            dict_Hor_info[hor2]['markserial-end']=markserial 
                            if    dict_Hor_info[hor2]['father_layer']<5  and   mer6 != '.' and mer6!='1':dict_Hor_info[hor2]['father_layer']=5
                            elif  dict_Hor_info[hor2]['father_layer']<4  and   mer5 != '.' and mer5!='1':dict_Hor_info[hor2]['father_layer']=4
                            elif  dict_Hor_info[hor2]['father_layer']<3  and   mer4 != '.' and mer4!='1':dict_Hor_info[hor2]['father_layer']=3
                            elif  dict_Hor_info[hor2]['father_layer']<2  and   mer3 != '.' and mer3!='1':dict_Hor_info[hor2]['father_layer']=2
                        if mer3 != '.' and mer3!='1':  
                            if hor3 not in dict_Hor_info:       
                                dict_Hor_info[hor3]={}
                                dict_Hor_info[hor3]['mer']=mer3
                                dict_Hor_info[hor3]['segment-end']=int(markserial)+int(mer3)-1 
                                dict_Hor_info[hor3]['segment']=mark
                                dict_Hor_info[hor3]['segment_length']=int(circ_len)
                                dict_Hor_info[hor3]['circ-start']=circ_serial
                                dict_Hor_info[hor3]['circ-end']=circ_serial
                                dict_Hor_info[hor3]['position-start']=chrserial_1
                                dict_Hor_info[hor3]['position-end']=chrserial_2
                                dict_Hor_info[hor3]['markserial-start']=markserial
                                dict_Hor_info[hor3]['markserial-end']=markserial                                    
                                dict_Hor_info[hor3]['father_layer']=1
                                dict_Hor_info[hor3]['layer']=3
                            else: 
                                if mark!='delete' and dict_Hor_info[hor3]['segment-end']>=int(markserial):
                                    dict_Hor_info[hor3]['segment']+=f"|{mark}"
                                if dict_Hor_info[hor3]['segment-end']>=int(markserial):   
                                    dict_Hor_info[hor3]['segment_length']+=int(circ_len)                                    
                            dict_Hor_info[hor3]['circ-end']=circ_serial
                            dict_Hor_info[hor3]['position-end']=chrserial_2
                            dict_Hor_info[hor3]['markserial-end']=markserial 
                            if    dict_Hor_info[hor3]['father_layer']<4  and   mer6 != '.' and mer6!='1':dict_Hor_info[hor3]['father_layer']=4
                            elif  dict_Hor_info[hor3]['father_layer']<3  and   mer5 != '.' and mer5!='1':dict_Hor_info[hor3]['father_layer']=3
                            elif  dict_Hor_info[hor3]['father_layer']<2  and   mer4 != '.' and mer4!='1':dict_Hor_info[hor3]['father_layer']=2
                        if mer4 != '.' and mer4!='1':  
                            if hor4 not in dict_Hor_info:       
                                dict_Hor_info[hor4]={}
                                dict_Hor_info[hor4]['mer']=mer4
                                dict_Hor_info[hor4]['segment-end']=int(markserial)+int(mer4)-1 
                                dict_Hor_info[hor4]['segment']=mark
                                dict_Hor_info[hor4]['segment_length']=int(circ_len)
                                dict_Hor_info[hor4]['circ-start']=circ_serial
                                dict_Hor_info[hor4]['circ-end']=circ_serial
                                dict_Hor_info[hor4]['position-start']=chrserial_1
                                dict_Hor_info[hor4]['position-end']=chrserial_2
                                dict_Hor_info[hor4]['markserial-start']=markserial
                                dict_Hor_info[hor4]['markserial-end']=markserial                                    
                                dict_Hor_info[hor4]['father_layer']=1
                                dict_Hor_info[hor4]['layer']=4
                            else: 
                                if mark!='delete' and dict_Hor_info[hor4]['segment-end']>=int(markserial):
                                    dict_Hor_info[hor4]['segment']+=f"|{mark}"
                                if dict_Hor_info[hor4]['segment-end']>=int(markserial):   
                                    dict_Hor_info[hor4]['segment_length']+=int(circ_len)                                    
                            dict_Hor_info[hor4]['circ-end']=circ_serial
                            dict_Hor_info[hor4]['position-end']=chrserial_2
                            dict_Hor_info[hor4]['markserial-end']=markserial 
                            if    dict_Hor_info[hor4]['father_layer']<3  and   mer6 != '.' and mer6!='1':dict_Hor_info[hor4]['father_layer']=3
                            elif  dict_Hor_info[hor4]['father_layer']<2  and   mer5 != '.' and mer5!='1':dict_Hor_info[hor4]['father_layer']=2
                        if mer5 != '.' and mer5!='1':  
                            if hor5 not in dict_Hor_info:       
                                dict_Hor_info[hor5]={}
                                dict_Hor_info[hor5]['mer']=mer5
                                dict_Hor_info[hor5]['segment-end']=int(markserial)+int(mer5)-1 
                                dict_Hor_info[hor5]['segment']=mark
                                dict_Hor_info[hor5]['segment_length']=int(circ_len)
                                dict_Hor_info[hor5]['circ-start']=circ_serial
                                dict_Hor_info[hor5]['circ-end']=circ_serial
                                dict_Hor_info[hor5]['position-start']=chrserial_1
                                dict_Hor_info[hor5]['position-end']=chrserial_2
                                dict_Hor_info[hor5]['markserial-start']=markserial
                                dict_Hor_info[hor5]['markserial-end']=markserial                                    
                                dict_Hor_info[hor5]['father_layer']=1
                                dict_Hor_info[hor5]['layer']=5
                            else: 
                                if mark!='delete' and dict_Hor_info[hor5]['segment-end']>=int(markserial):
                                    dict_Hor_info[hor5]['segment']+=f"|{mark}"
                                if dict_Hor_info[hor5]['segment-end']>=int(markserial):   
                                    dict_Hor_info[hor5]['segment_length']+=int(circ_len)                                    
                            dict_Hor_info[hor5]['circ-end']=circ_serial
                            dict_Hor_info[hor5]['position-end']=chrserial_2
                            dict_Hor_info[hor5]['markserial-end']=markserial 
                            if    dict_Hor_info[hor5]['father_layer']<2  and   mer6 != '.' and mer6!='1':dict_Hor_info[hor5]['father_layer']=2
                        if mer6 != '.' and mer6!='1':  
                            if hor6 not in dict_Hor_info:       
                                dict_Hor_info[hor6]={}
                                dict_Hor_info[hor6]['mer']=mer6
                                dict_Hor_info[hor6]['segment-end']=int(markserial)+int(mer6)-1 
                                dict_Hor_info[hor6]['segment']=mark
                                dict_Hor_info[hor6]['segment_length']=int(circ_len)
                                dict_Hor_info[hor6]['circ-start']=circ_serial
                                dict_Hor_info[hor6]['circ-end']=circ_serial
                                dict_Hor_info[hor6]['position-start']=chrserial_1
                                dict_Hor_info[hor6]['position-end']=chrserial_2
                                dict_Hor_info[hor6]['markserial-start']=markserial
                                dict_Hor_info[hor6]['markserial-end']=markserial                                    
                                dict_Hor_info[hor6]['father_layer']=1
                                dict_Hor_info[hor6]['layer']=6
                            else: 
                                if mark!='delete' and dict_Hor_info[hor6]['segment-end']>=int(markserial):
                                    dict_Hor_info[hor6]['segment']+=f"|{mark}"
                                if dict_Hor_info[hor6]['segment-end']>=int(markserial):   
                                    dict_Hor_info[hor6]['segment_length']+=int(circ_len)                                    
                            dict_Hor_info[hor6]['circ-end']=circ_serial
                            dict_Hor_info[hor6]['position-end']=chrserial_2
                            dict_Hor_info[hor6]['markserial-end']=markserial 
            with open(output_file,'w') as f2:
                f2.write(f"mer\tsegment\tsegment_type\tsegment_type_num\tsegment_length\tHOR_type\tcirc_start\tcirc_end\tpos_start\tpos_end\tmarkserial_start\tmarkserial_end\tlayer\tfather_layer\tpos_length\tmarkserial_num\tHOR_repeat_num\n")
                for Hor,info in dict_Hor_info.items():
                    mer=                    int(info['mer'])
                    segment=                info['segment']
                    segment_length=         info['segment_length']
                    circ_start=             info['circ-start']
                    circ_end=               info['circ-end']
                    pos_start=              int(info['position-start'])
                    pos_end=                int(info['position-end'])
                    markserial_start=       int(info['markserial-start'])
                    markserial_end=         int(info['markserial-end'])    

                    layer=                  info['layer']
                    father_layer=           info['father_layer']
                    
                    pos_length=     abs(pos_end-pos_start)+1
                    markserial_num= markserial_end-markserial_start+1
                    HOR_repeat_num= markserial_num/mer
                    
                    segment_arr=    segment.split('|')
                    segment_set=set()
                    for one in segment_arr:
                        segment_set.add(one)
                    segment_list=list(segment_set)
                    segment_list.sort()
                    segment_type=    ';'.join(segment_list)
                    segment_type_num= len(segment_list)
                    if segment_type_num>1:HOR_type='MT-HOR'
                    else:HOR_type='ST-HOR'
                    segment=normalize_hor_pattern_rotation_only(segment)
                    f2.write(f"{mer}\t{segment}\t{segment_type}\t{segment_type_num}\t{segment_length}\t{HOR_type}\t{circ_start}\t{circ_end}\t{pos_start}\t{pos_end}\t{markserial_start}\t{markserial_end}\t{layer}\t{father_layer}\t{pos_length}\t{markserial_num}\t{HOR_repeat_num}\n")          
                        
                
        with Pool(processes=thread) as pool:
            for i, result in enumerate(pool.imap(run_step, region_list), start=1):
                progress = (i / len(region_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush() 
    if argv1=="stepall" or argv1=="step12"  or argv1=="step12.2":    
            
        if os.path.exists('./samples_satellite/12_Hor_stat/2_plot_mer_type')==True:     
            subprocess.run(["rm -r ./samples_satellite/12_Hor_stat/2_plot_mer_type"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/12_Hor_stat/2_plot_mer_type"], shell=True)   
        
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)  
        dict_mer_type_length={}      
        ALL_length=0;ALL_length_pure=0
        for one_region in region_list:
            input_file=f"./samples_satellite/12_Hor_stat/1_blocks_stat/{one_region}" 
            
            with open(input_file,'r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t') 
                        
                        mer=                    int(eachline_arr[0])
                        HOR_type=               eachline_arr[5]
                        pos_length=             int(eachline_arr[-3])
                        layer=                  int(eachline_arr[12])
                        if mer not in dict_mer_type_length:         dict_mer_type_length[mer]={}
                        if HOR_type not in dict_mer_type_length[mer]:           dict_mer_type_length[mer][HOR_type]=0
                        dict_mer_type_length[mer][HOR_type]+=pos_length
                        ALL_length+=pos_length
                        if layer==1:ALL_length_pure+=pos_length

        print(f'Total length {ALL_length}\t\t\t\tNon-overlapping/pure HOR length {ALL_length_pure}')                
        with open(f"./samples_satellite/12_Hor_stat/2_plot_mer_type/0_HOR_ALL_length" ,'w')  as f:           
            f.write(f"{ALL_length}\t##This total length includes all layers, which is longer than pure layer1\n")
            f.write(f"{ALL_length_pure}\t##This total length includes layer=1 HORs, which is longer than pure HOR length\n")
                    
        with open(f"./samples_satellite/12_Hor_stat/2_plot_mer_type/1_mer_type" ,'w')  as f1:
            f1.write(f"mer\tHOR_type\tlength\tpercent\n")
            with open(f"./samples_satellite/12_Hor_stat/2_plot_mer_type/1_mer_MT2ST" ,'w')  as f2:
                f2.write(f"mer\tST-length\tMT-length\tsum_length\tSTMT_percent\tMT_vs_ST\n")
                for mer,dict_type_len in dict_mer_type_length.items():
                    if 'ST-HOR' in    dict_type_len:
                        ST_length=  dict_type_len['ST-HOR'];        ST_percent=ST_length/ALL_length_pure
                    else:
                        ST_length=0;ST_percent=0
                    if 'MT-HOR' in    dict_type_len:    
                        MT_length=  dict_type_len['MT-HOR'];        MT_percent=MT_length/ALL_length_pure
                    else:
                        MT_length=0;MT_percent=0    
                    STMT_length= ST_length+MT_length;        STMT_percent=STMT_length/ALL_length_pure
                    if ST_length!=0:
                        MT_vs_ST= MT_length/ST_length
                    else:
                        MT_vs_ST='.'
                    f1.write(f"{mer}\tST-HOR\t{ST_length}\t{ST_percent}\n")
                    f1.write(f"{mer}\tMT-HOR\t{MT_length}\t{MT_percent}\n")
                    f2.write(f"{mer}\t{ST_length}\t{MT_length}\t{STMT_length}\t{STMT_percent}\t{MT_vs_ST}\n")
    if argv1=="stepall" or argv1=="step12"  or argv1=="step12.2p":               
        R_txt=r'''
library(ggplot2)
library(dplyr)

###HORhor
print("step12.2p")
{
  input_file <- read.table('1_mer_type', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file$HOR_type <- factor(input_file$HOR_type, levels = c("ST-HOR", "MT-HOR"))
  # Use ggplot2 to draw scatter plot
  p <- ggplot(input_file, aes(x = mer, y = percent, fill = HOR_type)) +
    geom_bar(stat = "identity", position = "stack",width=0.65) +
    labs(title = "Basic Stacked Bar Chart", x = "Category", y = "Value", fill = "Subcategory")+
    theme_classic() +scale_x_continuous(name = "MER", limits = c(1, 43))+   scale_y_continuous(name = "MER")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Hide legend
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  
  # Save as PDF
  pdf("1_mer_type.pdf", width = 6 / 2.54, height = 5 / 2.54)
  print(p)
  dev.off()


}
    
        '''
        new_directory = "./samples_satellite/12_Hor_stat/2_plot_mer_type"
        os.chdir(new_directory)
        with open('./pic_ggplot_step12_2.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        subprocess.run(['Rscript pic_ggplot_step12_2.R'], shell=True)    
        os.chdir('../../../')                    
    if argv1=="stepall" or argv1=="step12"  or argv1=="step12.3":    
        print('step12.3_Statistics of mer-length-num')
        if os.path.exists('./samples_satellite/12_Hor_stat/3_plot_mer_length')==True:     
            subprocess.run(["rm -r ./samples_satellite/12_Hor_stat/3_plot_mer_length"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/12_Hor_stat/3_plot_mer_length"], shell=True)   
        
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)  
        dict_monomer_mer_length_info={}   
        HOR_len_all_repeatnum3=0
        for one_region in region_list:
            input_file=f"./samples_satellite/12_Hor_stat/1_blocks_stat/{one_region}" 
            
            with open(input_file,'r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t') 
                        mer=                    int(eachline_arr[0])
                        segment_type=           eachline_arr[2]
                        segment_length=         int(eachline_arr[4])
                        pos_length=             int(eachline_arr[-3])
                        HOR_repeat_num=         float(eachline_arr[-1])
                        if HOR_repeat_num<3:continue
                        HOR_len_all_repeatnum3+=pos_length
                        if segment_type=='28-28-28-23':monomer='107n'
                        elif segment_type=='28-28-28-28-23':monomer='135n'
                        else: monomer='other'
                        if monomer not in dict_monomer_mer_length_info:     dict_monomer_mer_length_info[monomer]={}
                        if mer not in dict_monomer_mer_length_info[monomer]:         dict_monomer_mer_length_info[monomer][mer]={}
                        if segment_length not in dict_monomer_mer_length_info[monomer][mer]:           
                            dict_monomer_mer_length_info[monomer][mer][segment_length]={}
                            dict_monomer_mer_length_info[monomer][mer][segment_length]['num']=0
                            dict_monomer_mer_length_info[monomer][mer][segment_length]['sum_length']=0
                        dict_monomer_mer_length_info[monomer][mer][segment_length]['num']+=1
                        dict_monomer_mer_length_info[monomer][mer][segment_length]['sum_length']+=pos_length
        with open ("./samples_satellite/12_Hor_stat/3_plot_mer_length/0_HOR_ALL_length",'w') as f:
            f.write(f"{HOR_len_all_repeatnum3}\t#This count excludes HORs with repeat count less than 3\n")
        print(f"HOR_len_all_repeatnum3\t\t{HOR_len_all_repeatnum3}\t\t\n")
        with open(f"./samples_satellite/12_Hor_stat/3_plot_mer_length/1_mer_length" ,'w')  as f:
            f.write(f"monomer\tmer\tunit_length\tnum\tsum_length\tsum_length_percent\tsum_length_lg10\n")
            for monomer,dict_mer_length_info in dict_monomer_mer_length_info.items():
                for mer,dict_length_info in dict_mer_length_info.items():
                    for length,info in dict_length_info.items():
                        num=            info['num']
                        sum_length=     int(info['sum_length'])
                        sum_length_percent=sum_length/HOR_len_all_repeatnum3
                        sum_length_lg10=    math.log10(sum_length)
                        f.write(f"{monomer}\t{mer}\t{length}\t{num}\t{sum_length}\t{sum_length_percent}\t{sum_length_lg10}\n")
    if argv1=="stepall" or argv1=="step12"  or argv1=="step12.3" or argv1=="step12.3p":               
        R_txt=r'''
library(ggplot2)
library(dplyr)

###HORhor
print("step12.3p")
{
  input_file <- read.table('1_mer_length', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    filtered_input_file<- input_file %>%  filter(num>100 & sum_length > 100)
  # Use ggplot2 to draw scatter plot
p=ggplot(filtered_input_file) +
  geom_point(aes(x = mer, y = unit_length, size = sum_length,alpha=sum_length_lg10,color=monomer)) +
  labs(x = "mer", y = "unit_length", color = "sum_length") +
  theme_classic()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #legend.position = "none",  # Hide legend
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +guides(color = guide_legend(override.aes = list(size = 3), title.position = "top", title.hjust = 0.5,keywidth  = unit(0.3, "cm"),keyheight = unit(0.3, "cm")))
  
  # Save as PDF
  pdf("1_mer_length_dot.pdf", width = 15 / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()


}
{
  # Use ggplot2 to draw bar chart
p=ggplot(filtered_input_file) +
  geom_col(aes(
    x = unit_length,      # x-axis variable
    y = sum_length_percent,       # Bar height
    fill = monomer        # Fill color by monomer
  ),
  width = 20,          # Widen bars (default 0.9, range 0~1)
  ) +
  labs(
    x = "Unit Length",    # x-axis label
    y = "Sum Length percent",     # y-axis label
    fill = "Monomer"      # Legend title
  )  +
   scale_fill_manual(
    values = c("#006699", "#33cc33", "#757575"),  # Manually specify color order
    name = "Monomer"      # Legend title (optional)
  ) +
  theme_classic()+   coord_cartesian(ylim = c(0, 0.6))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #legend.position = "none",  # Hide legend
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +guides(color = guide_legend(override.aes = list(size = 3), title.position = "top", title.hjust = 0.5,keywidth  = unit(0.3, "cm"),keyheight = unit(0.3, "cm")))
  
  # Save as PDF
  pdf("1_mer_length_col.pdf", width = 15 / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
}

{
  # Use ggplot2 to draw bar chart
p=ggplot(filtered_input_file) +
  geom_col(aes(
    x = unit_length,      # x-axis variable
    y = sum_length_percent,       # Bar height
    fill = monomer        # Fill color by monomer
  ),
  width = 1,          # Widen bars (default 0.9, range 0~1)
  ) +
  labs(
    x = "Unit Length",    # x-axis label
    y = "Sum Length percent",     # y-axis label
    fill = "Monomer"      # Legend title
  )  +
   scale_fill_manual(
    values = c("#006699", "#33cc33", "#757575"),  # Manually specify color order
    name = "Monomer"      # Legend title (optional)
  ) +
  theme_classic()+
  scale_x_continuous(name = "Unit Length", limits = c(100, 501))+   coord_cartesian(ylim = c(0, 0.013)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #legend.position = "none",  # Hide legend
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +guides(color = guide_legend(override.aes = list(size = 3), title.position = "top", title.hjust = 0.5,keywidth  = unit(0.3, "cm"),keyheight = unit(0.3, "cm")))
  
  # Save as PDF
  pdf("1_mer_length_col_part.pdf", width = 15 / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
}

        '''
        new_directory = "./samples_satellite/12_Hor_stat/3_plot_mer_length"
        os.chdir(new_directory)
        with open('./pic_ggplot_step12_3.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        subprocess.run(['Rscript pic_ggplot_step12_3.R'], shell=True)    
        os.chdir('../../../')                    
    if argv1=="stepall" or argv1=="step12"  or argv1=="step12.4":    
        print('step12.4_Statistics of mer-length-num')
        if os.path.exists('./samples_satellite/12_Hor_stat/4_segmenttype_mer_length_stat')==True:     
            subprocess.run(["rm -r ./samples_satellite/12_Hor_stat/4_segmenttype_mer_length_stat"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/12_Hor_stat/4_segmenttype_mer_length_stat"], shell=True)   
        
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)  
        dict_segmenttype_mer_length_info={}      
        for one_region in region_list:
            input_file=f"./samples_satellite/12_Hor_stat/1_blocks_stat/{one_region}" 
            
            with open(input_file,'r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t') 
                        mer=                    int(eachline_arr[0])
                        segment_type=           eachline_arr[2]
                        segment_length=         int(eachline_arr[4])
                        pos_length=             int(eachline_arr[-3])
                        HOR_repeat_num=         float(eachline_arr[-1])
                        if HOR_repeat_num<3:continue
                        if segment_type not in dict_segmenttype_mer_length_info:     dict_segmenttype_mer_length_info[segment_type]={}
                        if mer not in dict_segmenttype_mer_length_info[segment_type]:         dict_segmenttype_mer_length_info[segment_type][mer]={}
                        if segment_length not in dict_segmenttype_mer_length_info[segment_type][mer]:           
                            dict_segmenttype_mer_length_info[segment_type][mer][segment_length]={}
                            dict_segmenttype_mer_length_info[segment_type][mer][segment_length]['num']=0
                            dict_segmenttype_mer_length_info[segment_type][mer][segment_length]['sum_length']=0
                        dict_segmenttype_mer_length_info[segment_type][mer][segment_length]['num']+=1
                        dict_segmenttype_mer_length_info[segment_type][mer][segment_length]['sum_length']+=pos_length
                    
        with open(f"./samples_satellite/12_Hor_stat/4_segmenttype_mer_length_stat/1_segmenttype_mer_length_stat" ,'w')  as f:
            f.write(f"segment_type\tmer\tunit_length\tnum\tsum_length\tsum_length_lg10\n")
            for segment_type,dict_mer_length_info in dict_segmenttype_mer_length_info.items():
                for mer,dict_length_info in dict_mer_length_info.items():
                    for length,info in dict_length_info.items():
                        num=            info['num']
                        sum_length=     int(info['sum_length'])
                        sum_length_lg10=    math.log10(sum_length)
                        f.write(f"{segment_type}\t{mer}\t{length}\t{num}\t{sum_length}\t{sum_length_lg10}\n")            
    if argv1=="stepall" or argv1=="step12"  or argv1=="step12.5":    
        print('step12.5_Statistics of mer-length-num')
        if os.path.exists('./samples_satellite/12_Hor_stat/5_layer')==True:     
            subprocess.run(["rm -r ./samples_satellite/12_Hor_stat/5_layer"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/12_Hor_stat/5_layer"], shell=True)   
        
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)  
        dict_layer_info={}     
        dict_fatherlayer_info={} 
        HOR_len_all_repeatnum3=0
        for one_region in region_list:
            input_file=f"./samples_satellite/12_Hor_stat/1_blocks_stat/{one_region}" 
            
            with open(input_file,'r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t') 
                        mer=                    int(eachline_arr[0])
                        layer=                  int(eachline_arr[12])
                        father_layer=           int(eachline_arr[13])
                        pos_length=             int(eachline_arr[-3])
                        HOR_repeat_num=         float(eachline_arr[-1])
                        if HOR_repeat_num<3:continue
                        HOR_len_all_repeatnum3+=pos_length
                        ###
                        if layer not in dict_layer_info:     
                            dict_layer_info[layer]={}
                            dict_layer_info[layer]['num']=0
                            dict_layer_info[layer]['len']=0
                        dict_layer_info[layer]['num']+=1    
                        dict_layer_info[layer]['len']+=pos_length
                        ###
                        if father_layer not in dict_fatherlayer_info:     
                            dict_fatherlayer_info[father_layer]={}
                            dict_fatherlayer_info[father_layer]['num']=0
                            dict_fatherlayer_info[father_layer]['len']=0
                        dict_fatherlayer_info[father_layer]['num']+=1    
                        dict_fatherlayer_info[father_layer]['len']+=pos_length                        
        with open(f"./samples_satellite/12_Hor_stat/5_layer/1_layer_stat" ,'w')  as f:
            f.write(f"layer\tnum\tlength\tpercentage\n")
            for layer,info in dict_layer_info.items():
                num=info['num']
                length=info['len']
                f.write(f"L{layer}\t{num}\t{length}\t{length/HOR_len_all_repeatnum3}\n")                       
            for layer,info in dict_fatherlayer_info.items():
                num=info['num']
                length=info['len']
                f.write(f"B{layer}\t{num}\t{length}\t{length/HOR_len_all_repeatnum3}\n")    
    if argv1=="stepall" or argv1=="step12"  or argv1=="step12.6":    
        print('step12.6_Statistics of mer-length-num')
        if os.path.exists('./samples_satellite/12_Hor_stat/6_HOR_size')==True:     
            subprocess.run(["rm -r ./samples_satellite/12_Hor_stat/6_HOR_size"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/12_Hor_stat/6_HOR_size"], shell=True)   
        
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)  
        output_file=f"./samples_satellite/12_Hor_stat/6_HOR_size/0_sum"         
        with open (output_file,'w') as f2:
            for one_region in region_list:
                input_file=f"./samples_satellite/12_Hor_stat/1_blocks_stat/{one_region}"
                with open(input_file,'r') as f:
                    next(f)
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t') 
                        HOR_repeat_num=         float(eachline_arr[-1])
                        if HOR_repeat_num<3:continue
                        sample=one_region.split(':')[0]
                        chromosome=one_region.split(':')[1][:-1].replace('region_','Chr')
                        f2.write(f"{sample}\t{chromosome}\t{one_region}\t{eachline}\n")
        subprocess.run(["sort -k 18,18nr ./samples_satellite/12_Hor_stat/6_HOR_size/0_sum > ./samples_satellite/12_Hor_stat/6_HOR_size/0_sum2"], shell=True) 
        ###
        with open ("./samples_satellite/12_Hor_stat/6_HOR_size/0_sum",'w') as f2:
            f2.write("sample\tchromosome\tregion\tmer\tsegment\tsegment_type\tsegment_type_num\tsegment_length\tHOR_type\tcirc_start\tcirc_end\tpos_start\tpos_end\tmarkserial_start\tmarkserial_end\tlayer\tfather_layer\tpos_length\tmarkserial_num\tHOR_repeat_num\n")
            with open ("./samples_satellite/12_Hor_stat/6_HOR_size/0_sum2",'r') as f:
                for line in f:
                    eachline=line.strip()
                    f2.write(eachline+'\n')
        subprocess.run(["rm ./samples_satellite/12_Hor_stat/6_HOR_size/0_sum2"], shell=True)         
    if argv1=="stepall" or argv1=="step12"  or argv1=="step12.6p1":                  
        R_txt=r'''
library(ggplot2)
library(dplyr)

print("")
{
  input_file=read.table('0_sum', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
# Filter data
filtered_input <- input_file %>% 
  #filter(HOR_repeat_num >= 3 ) %>%
  mutate(
    #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
    category = case_when(
      mer == 2 ~ "2",
      mer == 3 ~ "3",
      mer == 4 ~ "4",
      mer == 5 ~ "5",
      mer == 6 ~ "6",
      mer == 7 ~ "7",
      mer == 8 ~ "8",
      mer == 9 ~ "9",
      mer == 10 ~ "10",
      mer == 11 ~ "11",
      mer == 12 ~ "12",
      mer == 13 ~ "13",
      mer == 14 ~ "14",
      mer == 15 ~ "15",
      mer == 16 ~ "16",
      mer == 17 ~ "17",
      mer == 18 ~ "18",
      mer == 19 ~ "19",
      mer == 20 ~ "20",
      mer>20 & mer<=30 ~ "21-30",
      mer>30 ~ "30+",
      TRUE ~ "Other"
    )
  )%>%
  mutate(
    dot_size = pos_length/30000*3

  )
# Define color values
color_values <- c(
  "2" = "#b5e48c",
  "3" = "#ff8566",
  "4" = "#ECD678",
  "5" = "#ade8f4", 
  "6" = "#cc99ff",
  "7" = "#7f9fcc",  
  "8" = "#ccccff"  ,
  "9" = "#ffccff",
  "10" = "#BAD694",   
  "11" = "#bcdcd6",  
  "12" = "#f49cbb",   
  "13" = "#DDDD7A",   
  "14" = "#bbdefb", 
  "15" = "#339933",  
  "16" = "#ffb4a2",   
  "17" = "#3385ff",  
  "18" = "#00b359",  
  "19" = "#45ABC3", 
  "20" = "#e5e600",
  "21-30"= "#ff80df",
  "30+"= "#cc3300",
  "Other"="black"
)

 
# Create plot object
p <- ggplot() +
    geom_point(data = filtered_input, aes(x = HOR_repeat_num, y = pos_length, color = category,size=dot_size)) +
theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_color_manual(values = color_values, drop = FALSE)+
     scale_size_identity()  ## Key: disable automatic scaling of size, directly use original values of dot_size
 
  
  # Save as PDF
  pdf(file = paste0('0_sum_mer_unit_length', ".pdf"), width = 24/ 2.54, height = 20 / 2.54)
  print(p)
  dev.off()
  
 
  # Create plot object
  filtered_input2<- filtered_input %>% filter(sample=="PN40024")
p <- ggplot() +
    geom_point(data = filtered_input2, aes(x = HOR_repeat_num, y = pos_length, color = category,size=dot_size)) +
theme_classic() +         
    scale_color_manual(values = color_values, drop = FALSE)+
     scale_size_identity()  ## Key: disable automatic scaling of size, directly use original values of dot_size
  pdf(file = paste0('0_sum_mer_unit_length_PN40024', ".pdf"), width = 24/ 2.54, height = 20 / 2.54)
  print(p)
  dev.off()
}
    '''
        with open('./samples_satellite/12_Hor_stat/6_HOR_size/0_sum_mer_unit_length.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './samples_satellite/12_Hor_stat/6_HOR_size/'
        os.chdir(new_directory)
        subprocess.run(['Rscript 0_sum_mer_unit_length.R'], shell=True)    
        os.chdir('../../../')                                                                    
    if argv1=="stepall" or argv1=="step12"  or argv1=="step12.6p2":                  
        R_txt=r'''
library(ggplot2)
library(dplyr)

print("")
{
  input_file=read.table('0_sum', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
# Filter data
    filtered_input<- input_file %>%
      mutate(
        category = case_when(
          chromosome == "Chr1" ~ "Chr1",
          chromosome == "Chr2" ~ "Chr2",
          chromosome == "Chr3" ~ "Chr3",
          chromosome == "Chr4" ~ "Chr4",
          chromosome == "Chr5" ~ "Chr5",
          chromosome == "Chr6" ~ "Chr6",
          chromosome == "Chr7" ~ "Chr7",
          chromosome == "Chr8" ~ "Chr8",
          chromosome == "Chr9" ~ "Chr9",
          chromosome == "Chr10" ~ "Chr10",
          chromosome == "Chr11" ~ "Chr11",
          chromosome == "Chr12" ~ "Chr12",
          chromosome == "Chr13" ~ "Chr13",
          chromosome == "Chr14" ~ "Chr14",
          chromosome == "Chr15" ~ "Chr15",
          chromosome == "Chr16" ~ "Chr16",
          chromosome == "Chr17" ~ "Chr17",
          chromosome == "Chr18" ~ "Chr18",
          chromosome == "Chr19" ~ "Chr19",
          TRUE ~ "Other"
        ),
      )%>%
  mutate(
    dot_size = pos_length/30000*3
  )
# Define color values
color_values <- c(
                                      "Chr1" = "red",  
                                      "Chr2" = "#8CC63F",  
                                      "Chr3" = "#CC0066",  
                                      "Chr4" = "#148EFF",  
                                      "Chr5" = "#B2DF8A",   
                                      "Chr6" = "#008000",   
                                      "Chr7" = "#FB9A99",   
                                      "Chr8" = "#00C08F",   
                                      "Chr9" = "#005FAF",  
                                      "Chr10" = "#FFCC00",  
                                      "Chr11" = "#00ccff",  
                                      "Chr12" = "#6A3D9A",  
                                      "Chr13" = "#FFFF00",  
                                      "Chr14" = "#A0522D",  
                                      "Chr15" = "#9966CC",  
                                      "Chr16" = "#ff4dff",  
                                      "Chr17" = "#ff6699", 
                                      "Chr18" = "#0073CF",  
                                      "Chr19" = "#FF8C00",  
                                      "Other" = "#CCEBC5"
)

 
# Create plot object
p <- ggplot() +
    geom_point(data = filtered_input, aes(x = HOR_repeat_num, y = pos_length, color = category,size=dot_size)) +
theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_color_manual(values = color_values, drop = FALSE)+
     scale_size_identity()  ## Key: disable automatic scaling of size, directly use original values of dot_size
 
  
  # Save as PDF
  pdf(file = paste0('0_sum_chr_unit_length', ".pdf"), width = 24/ 2.54, height = 20 / 2.54)
  print(p)
  dev.off()
  
  
  
 # Create plot object 
 filtered_input2<- filtered_input %>% filter(sample=="PN40024")
p <- ggplot() +
    geom_point(data = filtered_input2, aes(x = HOR_repeat_num, y = pos_length, color = category,size=dot_size)) +
theme_classic() +         
    scale_color_manual(values = color_values, drop = FALSE)+
     scale_size_identity()  ## Key: disable automatic scaling of size, directly use original values of dot_size
  pdf(file = paste0('0_sum_chr_unit_length_PN40024', ".pdf"), width = 24/ 2.54, height = 20 / 2.54)
  print(p)
  dev.off()
  
}
    '''
        with open('./samples_satellite/12_Hor_stat/6_HOR_size/0_sum_chr_unit_length.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './samples_satellite/12_Hor_stat/6_HOR_size/'
        os.chdir(new_directory)
        subprocess.run(['Rscript 0_sum_chr_unit_length.R'], shell=True)    
        os.chdir('../../../../')

## Statistics of different monomer types
if argv1=="stepall" or "step12b" in argv1: 
    print ("step12b          Statistics of monomers, CEN107 and CEN107-like, etc.") 
    if argv1=="stepall" or argv1=="step12b"  or argv1=="step12b.1":
        print('Statistics of different mers')
        if  os.path.exists(f'./samples_satellite/12_monomer_stat')==False:
            subprocess.run([f"mkdir ./samples_satellite/12_monomer_stat"], shell=True) 
            
        if os.path.exists('./samples_satellite/12_monomer_stat/1_blocks_stat')==True:     
            subprocess.run(["rm -r ./samples_satellite/12_monomer_stat/1_blocks_stat"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/12_monomer_stat/1_blocks_stat"], shell=True)    
        
        print('Loading regions')
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)
                
        def run_step(one_region):
            input_file=f"./samples_satellite/11_hor/5_sum/{one_region}" 
            output_file=f"./samples_satellite/12_monomer_stat/1_blocks_stat/{one_region}" 
            dict_serial_info={};current_monomer_type=''
            with open(input_file,'r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t')    
                        ###
                        circ_serial	=       eachline_arr[1] 
                        chrserial_1=	    eachline_arr[2]    
                        chrserial_2=        eachline_arr[3] 
                        circ_len=           eachline_arr[10]
                        subunit1=                   eachline_arr[14]
                        subunit_revised_type=       eachline_arr[18]
                        if subunit_revised_type=='28-28-28-23':
                            if subunit_revised_type==subunit1:   monomer_type='CEN107'
                            else:                               monomer_type='CEN107-like'
                        elif subunit_revised_type=='28-28-28-28-23':
                            if subunit_revised_type==subunit1:   monomer_type='CEN135'
                            else:                               monomer_type='CEN135-like'  
                        elif subunit_revised_type=='28-28-23':
                            if subunit_revised_type==subunit1:   monomer_type='CEN79'
                            else:                               monomer_type='CEN79-like'   
                        elif subunit_revised_type=='28-23':
                            if subunit_revised_type==subunit1:   monomer_type='CEN51'
                            else:                               monomer_type='CEN51-like'                             
                        else:  monomer_type='other'
                        #mark=               eachline_arr[20]    
                        markserial=         eachline_arr[21]    
                        if current_monomer_type=='':    
                            serial=1
                            current_monomer_type=     monomer_type
                            dict_serial_info[serial]={}
                            current_monomer_start=      chrserial_1
                            current_monomer_end=        chrserial_2
                            dict_serial_info[serial]['circ_start']=circ_serial
                            dict_serial_info[serial]['circ_end']=circ_serial
                            dict_serial_info[serial]['position_start']=chrserial_1
                            dict_serial_info[serial]['position_end']=chrserial_2
                            dict_serial_info[serial]['monomer_type']=monomer_type
                        elif current_monomer_type!=monomer_type:  
                            serial+=1
                            dict_serial_info[serial]={}
                            current_monomer_type=monomer_type
                            current_monomer_start=      chrserial_1
                            current_monomer_end=        chrserial_2
                            dict_serial_info[serial]['circ_start']=circ_serial
                            dict_serial_info[serial]['circ_end']=circ_serial                            
                            dict_serial_info[serial]['position_start']=chrserial_1
                            dict_serial_info[serial]['position_end']=chrserial_2
                            dict_serial_info[serial]['monomer_type']=monomer_type
                        else:       
                            dict_serial_info[serial]['circ_end']=circ_serial
                            dict_serial_info[serial]['position_end']=chrserial_2
            with open(output_file,'w') as f2:
                f2.write(f"serial\tcirc_start\tcirc_end\tpos_start\tpos_end\tlength\tmonomer_type\n")
                for serial,info in dict_serial_info.items():
                    monomer_type=         info['monomer_type']
                    circ_start=             info['circ_start']
                    circ_end=               info['circ_end']
                    position_start=             int(info['position_start'])
                    position_end=               int(info['position_end'] )                   
                    length=abs(position_end-position_start)+1
                    f2.write(f"{serial}\t{circ_start}\t{circ_end}\t{position_start}\t{position_end}\t{length}\t{monomer_type}\n")          
                        
                
        with Pool(processes=thread) as pool:
            for i, result in enumerate(pool.imap(run_step, region_list), start=1):
                progress = (i / len(region_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()

if argv1=="stepall" or "step12_var" in argv1:
    subprocess.run(["mkdir -p ./samples_satellite/12_Hor_var"], shell=True)
    subprocess.run([f"mkdir -p ./samples_satellite/12_Hor_var/1_blocks_stat/"], shell=True)
    
    if argv1=="stepall" or argv1=="step12_var.1":    
        all_items = os.listdir('./samples_satellite/12_Hor_stat/1_blocks_stat')
        def run_step (one):
            with open(f"./samples_satellite/12_Hor_var/1_blocks_stat/{one}",'w') as f2:
                f2.write(f"mer\tsegment\tsegment_type\tsegment_type_num\tsegment_length\tHOR_type\tcirc_start\tcirc_end\tpos_start\tpos_end\tmarkserial_start\tmarkserial_end\tlayer\tfather_layer\tpos_length\tmarkserial_num\tHOR_repeat_num\tvar_percent\n")
            #print(one)
            dict_circ_info={}
            dict_serial_circ={}
            with open(f"./samples_satellite/11_hor/5_sum/{one}",'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=line.split('\t')        
                    circ_serial=eachline_arr[1]
                    circ_serial_num=int(circ_serial[5:])
                    circ_seq=eachline_arr[5]
                    mark=eachline_arr[20]
                    mark_serilal=int(eachline_arr[21])
                    dict_circ_info[circ_serial_num]=[circ_seq,mark_serilal]
                    if mark!='delete':dict_serial_circ[mark_serilal]=circ_serial_num
            if dict_serial_circ=={} or dict_circ_info=={}:return False
            #print(dict_serial_circ)
            
            with open(f"./samples_satellite/12_Hor_stat/1_blocks_stat/{one}",'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=line.split('\t')
                    mer,segment,segment_type,segment_type_num,segment_length,HOR_type,circ_start,circ_end,pos_start,pos_end,markserial_start,markserial_end,layer,father_layer,pos_length,markserial_num,HOR_repeat_num=eachline_arr
                    if '.' in segment_type:continue
                    circ_start_num=int(circ_start[5:])
                    circ_end_num=int(circ_end[5:])
                    if int(father_layer)>1:
                        #print(circ_start)
                        serial_start=int(markserial_start)
                        serial_start_delta=serial_start+int(mer)
                        circ_start_delta=dict_serial_circ[serial_start_delta]-circ_start_num
                        seqlen_analysis=int(pos_length)*2-int(segment_length)   #####Sliding window, the final number of bases compared is twice the interval size - unit size
                        #if one=="V112.hap1:region_19-:17001902-17217884" and circ_start=="circ_734": print(f"\n\n{circ_start_delta}\n\n")
                        kk=circ_start_num
                        mismatch=0
                        while kk+circ_start_delta<=circ_end_num:
                            seq1=dict_circ_info[kk][0]
                            seq2=dict_circ_info[kk+circ_start_delta][0]
                            seq1_arr=seq1.split('|')
                            seq2_arr=seq2.split('|')
                            for base1,base2 in zip(seq1_arr,seq2_arr):
                                if base1!=base2:mismatch+=1
                            kk+=1
                        var_percent=mismatch/seqlen_analysis*100
                        with open(f"./samples_satellite/12_Hor_var/1_blocks_stat/{one}",'a') as f2:
                            f2.write(f"{eachline}\t{var_percent:.2f}\n")
                        
                        
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, all_items), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(all_items)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()    
    if argv1=="stepall" or argv1=="step12_var.1s":     
        all_items = os.listdir('./samples_satellite/12_Hor_stat/1_blocks_stat')
        with open(f"./samples_satellite/12_Hor_var/1_blocks_stat_sum",'w') as f3:
            f3.write(f"file\tsample\tchromosme\tsamplechr\tmer\tsegment\tsegment_type\tsegment_type_num\tsegment_length\tHOR_type\tcirc_start\tcirc_end\tpos_start\tpos_end\tmarkserial_start\tmarkserial_end\tlayer\tfather_layer\tpos_length\tmarkserial_num\tHOR_repeat_num\tvar_percent\n")
            for one in all_items:
                sample=one.split(':')[0]
                chromosome=one.split(':')[1][:-1].replace('region_','Chr')
                samplechr=sample+'|'+chromosome
                with open(f"./samples_satellite/12_Hor_var/1_blocks_stat/{one}",'r') as f:
                    next(f)
                    for line in f:
                        eachline=line.strip()
                        f3.write(f"{one}\t{sample}\t{chromosome}\t{samplechr}\t{eachline}\n")
    if argv1=="stepall" or argv1=="step12_var.2":       
        subprocess.run([f"mkdir -p ./samples_satellite/12_Hor_var/2_plot/"], shell=True)
        ##
        good_samples=["Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap1","V059.hap2","V063.hap1","V063.hap2","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap1","Chardonnay_hap2","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap1","V074.hap2","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap1","V062.hap2","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap1","PinotNoir_hap2","PinotNoir2_hap1","PinotNoir2_hap2","V069.hap1","V069.hap2","PN40024","PN40024_hap1","PN40024_hap2","V087.hap1","V087.hap2","V061.hap1","V061.hap2","V065.hap1","V065.hap2","V064.hap1","V064.hap2","ThompsonSeedless_hap1","ThompsonSeedless_hap2","ThompsonSeedless2_hap1","ThompsonSeedless2_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap1","V067.hap2","V058.hap1","V058.hap2","V107.hap1","V107.hap2","V108.hap1","V108.hap2","V081.hap1","V081.hap2","V015.hap1","V015.hap2","V105.hap1","V105.hap2","V126.hap1","V126.hap2","V112.hap1","V112.hap2","DavidiiGrape_hap1","DavidiiGrape_hap2","V007.hap1","V007.hap2","V008.hap1","V008.hap2","V012.hap1","V012.hap2","V100.hap1","V100.hap2","V099.hap1","V099.hap2","V098.hap1","V098.hap2","V117.hap1","V117.hap2","V120.hap1","V120.hap2","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap1","V124.hap2","V125.hap1","V125.hap2","WoollyGrape_hap1","WoollyGrape_hap2","V106.hap1","V106.hap2","V102.hap1","V102.hap2","V079.hap1","V079.hap2","V030.hap1","V030.hap2","V031.hap1","V031.hap2","V005.hap1","V005.hap2","V048.hap1","V048.hap2","V050.hap1","V050.hap2","V051.hap1","V051.hap2","V049.hap1","V049.hap2","V038.hap1","V038.hap2","V043.hap1","V043.hap2","V052.hap1","V052.hap2","V034.hap1","V034.hap2","V023.hap1","V023.hap2","V032.hap1","V032.hap2","V033.hap1","V033.hap2","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap1","V022.hap2","V040.hap1","V040.hap2","V019.hap1","V019.hap2","V041.hap1","V041.hap2","V077.hap1","V077.hap2","V037.hap1","V037.hap2","V096.hap1","V096.hap2","ShineMuscat_hap1","ShineMuscat_hap2","V036.hap1","V036.hap2","V076.hap1","V076.hap2","V055.hap1","V055.hap2","V018.hap1","V018.hap2","V072.hap1","V072.hap2"]
        
        ##
        dict_samplechr={}
        with open("./samples_satellite/2_good_regions_main",'r') as f:
            for line in f:
                eachline=line.strip()
                eachline_arr=eachline.split('\t')                
                sample,chromosome,start,end,length=eachline_arr
                if sample not in good_samples:continue
                samplechr=sample+'|'+chromosome
                dict_samplechr[samplechr]=[int(start),int(end)]
        
        with open('./samples_satellite/12_Hor_var/2_plot/0_Vsat1_pos','w') as f2:
            f2.write("sample\tchromosome\trevise_start\trevise_end\n")
            with open('./samples_satellite/2_good_regions') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    if eachline=='':continue
                    eachline_arr=eachline.split('\t')   
                    if len(eachline_arr)!=8:continue
                    sample,region_name,region_pos,bigblock_chrstart,bigblock_chrend,chr_region_length,strand,match_percent=eachline_arr
                    if sample not in good_samples:continue
                    chromosome=region_name.replace('region_','Chr')
                    samplechr=sample+'|'+chromosome
                    border1,border2=dict_samplechr[samplechr]
                    #print(chromosome)
                    if int(bigblock_chrstart)<border1 or int(bigblock_chrstart)>border2 :continue
                    if int(bigblock_chrend)<border1 or int(bigblock_chrend)>border2 :continue
                    #print(border1,border2)
                    revise_start=int(bigblock_chrstart)-border1+1
                    revise_end=int(bigblock_chrend)-border1+1
                    f2.write(f"{sample}\t{chromosome}\t{revise_start}\t{revise_end}\n")
                
        ##
        with open(f"./samples_satellite/12_Hor_var/2_plot/0_horvar",'w') as f2:
            f2.write("sample\tchromosome\tpos_start\tpos_end\tpos_start_revise\tpos_end_revise\tvar_percent\n")
            with open(f"./samples_satellite/12_Hor_var/1_blocks_stat_sum",'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    sample=eachline_arr[1]
                    chromosome=eachline_arr[2]
                    if sample not in good_samples:continue
                    samplechr=eachline_arr[3]
                    chrVSat1_start=dict_samplechr[samplechr][0]

                    pos_start	=int(eachline_arr[12])
                    pos_end=int(eachline_arr[13])
                    #print(pos_start,pos_end,chrVSat1_start)
                    pos_start_revise= pos_start-chrVSat1_start+1
                    pos_end_revise= pos_end-chrVSat1_start+1
                    var_percent=eachline_arr[-1]
                    f2.write(f"{sample}\t{chromosome}\t{pos_start}\t{pos_end}\t{pos_start_revise}\t{pos_end_revise}\t{var_percent}\n")
                
    if argv1=="stepall" or argv1=="step12_var.2p":                  
        Plot_txt=r"""
        library(ggplot2)
        library(dplyr)
        library(scales) 
        # Get all command line arguments
        
        ### Monomer
        print("")
        {
          input_file0=read.table('0_Vsat1_pos', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
          input_file0$chromosome <- factor(
              input_file0$chromosome,
              levels = c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10','Chr11','Chr12','Chr13','Chr14','Chr15','Chr16','Chr17','Chr18','Chr19')
            )  
        input_file0$sample <- factor(input_file0$sample, levels = c(
        "Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap1","V059.hap2","V063.hap1","V063.hap2","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap1","Chardonnay_hap2","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap1","V074.hap2","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap1","V062.hap2","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap1","PinotNoir_hap2","PinotNoir2_hap1","PinotNoir2_hap2","V069.hap1","V069.hap2","PN40024","PN40024_hap1","PN40024_hap2","V087.hap1","V087.hap2","V061.hap1","V061.hap2","V065.hap1","V065.hap2","V064.hap1","V064.hap2","ThompsonSeedless_hap1","ThompsonSeedless_hap2","ThompsonSeedless2_hap1","ThompsonSeedless2_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap1","V067.hap2","V058.hap1","V058.hap2","V107.hap1","V107.hap2","V108.hap1","V108.hap2","V081.hap1","V081.hap2","V015.hap1","V015.hap2","V105.hap1","V105.hap2","V126.hap1","V126.hap2","V112.hap1","V112.hap2","DavidiiGrape_hap1","DavidiiGrape_hap2","V007.hap1","V007.hap2","V008.hap1","V008.hap2","V012.hap1","V012.hap2","V100.hap1","V100.hap2","V099.hap1","V099.hap2","V098.hap1","V098.hap2","V117.hap1","V117.hap2","V120.hap1","V120.hap2","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap1","V124.hap2","V125.hap1","V125.hap2","WoollyGrape_hap1","WoollyGrape_hap2","V106.hap1","V106.hap2","V102.hap1","V102.hap2","V079.hap1","V079.hap2","V030.hap1","V030.hap2","V031.hap1","V031.hap2","V005.hap1","V005.hap2","V048.hap1","V048.hap2","V050.hap1","V050.hap2","V051.hap1","V051.hap2","V049.hap1","V049.hap2","V038.hap1","V038.hap2","V043.hap1","V043.hap2","V052.hap1","V052.hap2","V034.hap1","V034.hap2","V023.hap1","V023.hap2","V032.hap1","V032.hap2","V033.hap1","V033.hap2","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap1","V022.hap2","V040.hap1","V040.hap2","V019.hap1","V019.hap2","V041.hap1","V041.hap2","V077.hap1","V077.hap2","V037.hap1","V037.hap2","V096.hap1","V096.hap2","ShineMuscat_hap1","ShineMuscat_hap2","V036.hap1","V036.hap2","V076.hap1","V076.hap2","V055.hap1","V055.hap2","V018.hap1","V018.hap2","V072.hap1","V072.hap2"))
         # Check data and handle missing values
        input_file0 <- input_file0 %>%
          group_by(sample, chromosome) %>%
          mutate(
            ymin = as.numeric(factor(sample)) - 0.4,  # Assign y-axis position for each sample
            ymax = as.numeric(factor(sample)) + 0.4
          )          
          
        input_file=read.table('0_horvar', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')    
          input_file$chromosome <- factor(
              input_file$chromosome,
              levels = c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10','Chr11','Chr12','Chr13','Chr14','Chr15','Chr16','Chr17','Chr18','Chr19')
            )  
        input_file$sample <- factor(input_file$sample, levels = c(
        "Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap1","V059.hap2","V063.hap1","V063.hap2","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap1","Chardonnay_hap2","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap1","V074.hap2","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap1","V062.hap2","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap1","PinotNoir_hap2","PinotNoir2_hap1","PinotNoir2_hap2","V069.hap1","V069.hap2","PN40024","PN40024_hap1","PN40024_hap2","V087.hap1","V087.hap2","V061.hap1","V061.hap2","V065.hap1","V065.hap2","V064.hap1","V064.hap2","ThompsonSeedless_hap1","ThompsonSeedless_hap2","ThompsonSeedless2_hap1","ThompsonSeedless2_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap1","V067.hap2","V058.hap1","V058.hap2","V107.hap1","V107.hap2","V108.hap1","V108.hap2","V081.hap1","V081.hap2","V015.hap1","V015.hap2","V105.hap1","V105.hap2","V126.hap1","V126.hap2","V112.hap1","V112.hap2","DavidiiGrape_hap1","DavidiiGrape_hap2","V007.hap1","V007.hap2","V008.hap1","V008.hap2","V012.hap1","V012.hap2","V100.hap1","V100.hap2","V099.hap1","V099.hap2","V098.hap1","V098.hap2","V117.hap1","V117.hap2","V120.hap1","V120.hap2","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap1","V124.hap2","V125.hap1","V125.hap2","WoollyGrape_hap1","WoollyGrape_hap2","V106.hap1","V106.hap2","V102.hap1","V102.hap2","V079.hap1","V079.hap2","V030.hap1","V030.hap2","V031.hap1","V031.hap2","V005.hap1","V005.hap2","V048.hap1","V048.hap2","V050.hap1","V050.hap2","V051.hap1","V051.hap2","V049.hap1","V049.hap2","V038.hap1","V038.hap2","V043.hap1","V043.hap2","V052.hap1","V052.hap2","V034.hap1","V034.hap2","V023.hap1","V023.hap2","V032.hap1","V032.hap2","V033.hap1","V033.hap2","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap1","V022.hap2","V040.hap1","V040.hap2","V019.hap1","V019.hap2","V041.hap1","V041.hap2","V077.hap1","V077.hap2","V037.hap1","V037.hap2","V096.hap1","V096.hap2","ShineMuscat_hap1","ShineMuscat_hap2","V036.hap1","V036.hap2","V076.hap1","V076.hap2","V055.hap1","V055.hap2","V018.hap1","V018.hap2","V072.hap1","V072.hap2"))
         # Check data and handle missing values
        input_file <- input_file %>%
          #filter(!is.na(start), !is.na(end), !is.na(var_percent)) %>%  # Filter invalid rows
          group_by(sample, chromosome) %>%
          mutate(
            ymin = as.numeric(factor(sample)) - 0.4,  # Assign y-axis position for each sample
            ymax = as.numeric(factor(sample)) + 0.4
          )
        
        # Create plot object
        p <- ggplot() +
          geom_rect(
            data = input_file0,	
            aes(
              xmin = revise_start / 1000000,  # Convert to Mb
              xmax = revise_end / 1000000,
              ymin = ymin,
              ymax = ymax,
              
            ),fill = 'black'  # Fill with black
          )+
          geom_rect(
            data = input_file,
            aes(
              xmin = pos_start_revise / 1000000,  # Convert to Mb
              xmax = pos_end_revise / 1000000,
              ymin = ymin,
              ymax = ymax,
              fill = var_percent  # Fill with var_percent
            )
          ) +
          scale_fill_gradientn(
            colors = c("blue", "white", "red"),  # Color gradient
            values = scales::rescale(c(0, 3, 3.0001)),  # Key point: 3 is the cutoff, values >3 will be mapped to red
            limits = c(0, 3),  # Force display range to 0-3
            oob = squish,  # Values outside the range (>3) are compressed to the maximum value (3)
            name = "Variation (%)"
          ) +
          facet_grid(sample ~ chromosome,  drop = FALSE) +
          theme_classic() +
          theme(
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            strip.text.y = element_text(angle = 0)  # Display sample names horizontally
          )
        
        # Save as PDF
        pdf(file = "sum_for_plot.pdf", width = 400 / 2.54, height = 150 / 2.54)
        print(p)
        dev.off()
        }            
        """
        with open(f'./samples_satellite/12_Hor_var/2_plot/plot.R','w') as f:
            f.write(Plot_txt)
        os.chdir(f'./samples_satellite/12_Hor_var/2_plot/')
        subprocess.run([f'Rscript plot.R  '], shell=True)  #>null 2>&1 
        os.chdir('../../../')
        
## Visualize HORs in different blocks
if argv1=="stepall" or "step13" in argv1:
    if argv1=="stepall" or argv1=="step13"  or argv1=="step13.1":       
        print('Visualize HORs, 45min')
        if  os.path.exists(f'./samples_satellite/13_image/')==False:
            subprocess.run([f"mkdir ./samples_satellite/13_image/"], shell=True) 
            
        if  os.path.exists('./samples_satellite/13_image/1_Hor'): 
            subprocess.run(["rm -r ./samples_satellite/13_image/1_Hor"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/13_image/1_Hor"], shell=True)  
        subprocess.run(["mkdir ./samples_satellite/13_image/1_Hor/blocks"], shell=True)
        
        print('Loading regions')
        info_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                info_list.append([new_name,eachline_arr[3],eachline_arr[4]])

        R_txt=r'''
library(ggplot2)
library(dplyr)
# Get all command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Print arguments (for debugging)
print(args)
input_file_path=args[1]
chrserial_start=as.numeric(args[2])
chrserial_end=as.numeric(args[3])
print(input_file_path)
print(chrserial_start)
print(chrserial_end)
### Monomer
print("")
{input_file=read.table(input_file_path, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')

# Filter data
filtered_input <- input_file %>% 
  filter(pos_length > 300 ) %>%
  mutate(
    #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
    category = case_when(
      mer == 2 ~ "2",
      mer == 3 ~ "3",
      mer == 4 ~ "4",
      mer == 5 ~ "5",
      mer == 6 ~ "6",
      mer == 7 ~ "7",
      mer == 8 ~ "8",
      mer == 9 ~ "9",
      mer == 10 ~ "10",
      mer == 11 ~ "11",
      mer == 12 ~ "12",
      mer == 13 ~ "13",
      mer == 14 ~ "14",
      mer == 15 ~ "15",
      mer == 16 ~ "16",
      mer == 17 ~ "17",
      mer == 18 ~ "18",
      mer == 19 ~ "19",
      mer == 20 ~ "20",
      mer>20 & mer<=30 ~ "21-30",
      mer>30 ~ "30+",
      TRUE ~ "Other"
    )
  )


filtered_input_father1 <- filtered_input %>% filter(father_layer == 1)
filtered_input_father2 <- filtered_input %>% filter(father_layer == 2)
filtered_input_father3 <- filtered_input %>% filter(father_layer == 3)


filtered_input_father1_hor3 <- filtered_input_father1 %>% filter(HOR_repeat_num >= 3) 
filtered_input_father2_hor3 <- filtered_input_father2 %>% filter(HOR_repeat_num >= 3) 
filtered_input_father3_hor3 <- filtered_input_father3 %>% filter(HOR_repeat_num >= 3) 

filtered_input_layer1 <- filtered_input %>% filter(layer == 1)
filtered_input_layer2 <- filtered_input %>% filter(layer == 2)
filtered_input_layer3 <- filtered_input %>% filter(layer == 3)


filtered_input_layer1_hor3 <- filtered_input_layer1 %>% filter(HOR_repeat_num >= 3) 
filtered_input_layer2_hor3 <- filtered_input_layer2 %>% filter(HOR_repeat_num >= 3) 
filtered_input_layer3_hor3 <- filtered_input_layer3 %>% filter(HOR_repeat_num >= 3) 
# Define color values
color_values <- c(
  "2" = "#b5e48c",
  "3" = "#ff8566",
  "4" = "#ECD678",
  "5" = "#ade8f4", 
  "6" = "#cc99ff",
  "7" = "#7f9fcc",  
  "8" = "#ccccff"  ,
  "9" = "#ffccff",
  "10" = "#BAD694",   
  "11" = "#bcdcd6",  
  "12" = "#f49cbb",   
  "13" = "#DDDD7A",   
  "14" = "#bbdefb", 
  "15" = "#339933",  
  "16" = "#ffb4a2",   
  "17" = "#3385ff",  
  "18" = "#00b359",  
  "19" = "#45ABC3", 
  "20" = "#e5e600",
  "21-30"= "#ff80df",
  "30+"= "#cc3300",
  "Other"="black"
)

# Create plot object
p = ggplot() 
p =p+  geom_rect(aes(xmin = chrserial_start,xmax = chrserial_end,ymin = 95,ymax = 100),fill = '#6A356C',) 
p =p+  geom_rect(data = filtered_input_father3,
                 aes(xmin = pos_start,xmax = pos_end, ymin = 95,ymax = 100,fill = category)) 
##
p =p+  geom_rect(aes(xmin = chrserial_start,xmax = chrserial_end,ymin = 85,ymax = 90),fill = '#6A356C',) 
p =p+  geom_rect(data = filtered_input_father2,
                 aes(xmin = pos_start,xmax = pos_end, ymin = 85,ymax = 90,fill = category)) 
##
p =p+  geom_rect(aes(xmin = chrserial_start,xmax = chrserial_end,ymin = 75,ymax = 80),fill = '#6A356C',) 
p =p+  geom_rect(data = filtered_input_father1,
                 aes(xmin = pos_start,xmax = pos_end, ymin = 75,ymax = 80,fill = category)) 
p =p+  annotate("text", x = (chrserial_start+chrserial_end)/2, y = 65, label = "B1--B2--B3 ",size=4)
##
p =p+  geom_rect(aes(xmin = chrserial_start,xmax = chrserial_end,ymin = 45,ymax = 50),fill = '#6A356C',) 
p =p+  geom_rect(data = filtered_input_father3_hor3,
                 aes(xmin = pos_start,xmax = pos_end, ymin = 45,ymax = 50,fill = category)) 
##
p =p+  geom_rect(aes(xmin = chrserial_start,xmax = chrserial_end,ymin = 35,ymax = 40),fill = '#6A356C',) 
p =p+  geom_rect(data = filtered_input_father2_hor3,
                 aes(xmin = pos_start,xmax = pos_end, ymin = 35,ymax = 40,fill = category)) 
##
p =p+  geom_rect(aes(xmin = chrserial_start,xmax = chrserial_end,ymin = 25,ymax = 30),fill = '#6A356C',) 
p =p+  geom_rect(data = filtered_input_father1_hor3,
                 aes(xmin = pos_start,xmax = pos_end, ymin = 25,ymax = 30,fill = category)) 
p =p+  annotate("text", x = (chrserial_start+chrserial_end)/2, y = 15, label = "B1--B2--B3(repeat>=3) ",size=4)
###################################
p =p+  geom_rect(aes(xmin = chrserial_start,xmax = chrserial_end,ymin = 195,ymax = 200),fill = '#6A356C',) 
p =p+  geom_rect(data = filtered_input_layer3,
                 aes(xmin = pos_start,xmax = pos_end, ymin = 195,ymax = 200,fill = category)) 
##
p =p+  geom_rect(aes(xmin = chrserial_start,xmax = chrserial_end,ymin = 185,ymax = 190),fill = '#6A356C',) 
p =p+  geom_rect(data = filtered_input_layer2,
                 aes(xmin = pos_start,xmax = pos_end, ymin = 185,ymax = 190,fill = category)) 
##
p =p+  geom_rect(aes(xmin = chrserial_start,xmax = chrserial_end,ymin = 175,ymax = 180),fill = '#6A356C',) 
p =p+  geom_rect(data = filtered_input_layer1,
                 aes(xmin = pos_start,xmax = pos_end, ymin = 175,ymax = 180,fill = category)) 
p =p+  annotate("text", x = (chrserial_start+chrserial_end)/2, y = 165, label = "L1--L2--L3 ",size=4)
##
p =p+  geom_rect(aes(xmin = chrserial_start,xmax = chrserial_end,ymin = 145,ymax = 150),fill = '#6A356C',) 
p =p+  geom_rect(data = filtered_input_layer3_hor3,
                 aes(xmin = pos_start,xmax = pos_end, ymin = 145,ymax = 150,fill = category)) 
##
p =p+  geom_rect(aes(xmin = chrserial_start,xmax = chrserial_end,ymin = 135,ymax = 140),fill = '#6A356C',) 
p =p+  geom_rect(data = filtered_input_layer2_hor3,
                 aes(xmin = pos_start,xmax = pos_end, ymin = 135,ymax = 140,fill = category)) 
##
p =p+  geom_rect(aes(xmin = chrserial_start,xmax = chrserial_end,ymin = 125,ymax = 130),fill = '#6A356C',) 
p =p+  geom_rect(data = filtered_input_layer1_hor3,
                 aes(xmin = pos_start,xmax = pos_end, ymin = 125,ymax = 130,fill = category)) 
p =p+  annotate("text", x = (chrserial_start+chrserial_end)/2, y = 115, label = "L1--L2--L3(repeat>=3) ",size=4)
###################################

  #geom_vline(
  #  xintercept = c(chrserial_start , chrserial_end ),
  #  linetype = "dashed",
  #  color = "red",
  #  linewidth = 0.5
  #) 

p=p+  xlim(chrserial_start, chrserial_end) + 
  ylim(0, 200) + 
  labs(
    x = paste0(chrserial_start, '-', chrserial_end),
    y = "",
    fill = ""
  ) +
  theme_classic() +         
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_blank()
  ) +
  scale_fill_manual(values = color_values, drop = FALSE)
  
  
  # Save as PDF
  pdf(file = paste0(input_file_path, ".pdf"), width = 10 / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()

}
    '''
        with open('./samples_satellite/13_image/1_Hor/pic_ggplot_13.1.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        '''Which Science color scheme is this?
        "2" = "#CD6953", 
            "3" = "#8D4323",  
            "4" = "#ECD678",  
            "5" = "#7BB9AD",  
            "6" = "#A83E6E",  
            "7" = "#7EAE8E",   
            "8" = "#823D82",   
            "9" = "#6061A0",   
            "10" = "#96C086",   
            "11" = "#C679B0",  
            "12" = "#3A8FBA",  
            "13" = "#827DBA",  
            "14" = "#87B7CF",  
            "15" = "#066656",  
            "16" = "#45ABC3",  
            "17" = "#AAA3D0",  
            "18" = "#9F6FAD",  
            "19" = "#3C6298", 
            "20" = "#BAD694",  
            "21" = "#A85A81",  
            "22" = "#DDDD7A",  
            "23" = "#339933",  
            "24" = "#E4935C"
        '''
        R_legend_txt=r'''
        # Load required packages
library(ggplot2)
 
# Define color values
values <- c( 
"2" = "#b5e48c",
  "3" = "#ff8566",
  "4" = "#ECD678",
  "5" = "#ade8f4", 
  "6" = "#cc99ff",
  "7" = "#7f9fcc",  
  "8" = "#ccccff"  ,
  "9" = "#ffccff",
  "10" = "#BAD694",   
  "11" = "#bcdcd6",  
  "12" = "#f49cbb",   
  "13" = "#DDDD7A",   
  "14" = "#bbdefb", 
  "15" = "#339933",  
  "16" = "#ffb4a2",   
  "17" = "#3385ff",  
  "18" = "#00b359",  
  "19" = "#45ABC3", 
  "20" = "#e5e600",
  "21-30"= "#ff80df",
  "30+"= "#cc3300",
  "Other"="black"
)
 
# Create a dummy data frame for plotting
df <- data.frame(
  category = factor(names(values), levels = names(values)),  # Ensure categories are factors and in order
  dummy = 1  # Create a dummy column for plotting
)
 
# Create a plot object, only for displaying the legend
p <- ggplot(df, aes(x = dummy, y = dummy, fill = category)) +
  geom_point(shape = 22, size = 3, color = "black") +  # Use square points to display the legend
  scale_fill_manual(values = values, name = "Categories") +  # Set colors and legend title
  scale_shape_manual(values = 22) +  # Ensure all points are squares
  theme_void() +  # Hide all plot elements
  theme(legend.position = "bottom", legend.key.size = unit(1, "cm"))  # Set legend position and size
 
 
# Print the legend
  # Save as PDF
  pdf("./samples_satellite/13_image/1_Hor/legend.pdf", width = 20 / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
        '''
        with open('./samples_satellite/13_image/1_Hor/pic_ggplot_13.1_lengend.R','w',encoding='utf-8') as f:
            f.write(R_legend_txt)
        subprocess.run([f'Rscript ./samples_satellite/13_image/1_Hor/pic_ggplot_13.1_lengend.R '], shell=True)  
        
        
        import hashlib
        
        #sys.exit()
        def get_file_md5(file_path):
            """Calculate the MD5 hash of a file"""
            hash_md5 = hashlib.md5()
            with open(file_path, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    hash_md5.update(chunk)
            return hash_md5.hexdigest()
        
        def run_step(one_info):
            one_region,chrserial_start,chrserial_end=one_info
            length=abs(int(chrserial_end)-int(chrserial_start))+1
            input_file=f"./samples_satellite/12_Hor_stat/1_blocks_stat/{one_region}"
            num=0
            with open (input_file,'r') as f:
                next(f)
                for line in f:
                    num+=1
            if num<2:return False        
            input_file_hash=get_file_md5(input_file)
            subprocess.run([f'cp {input_file} ./samples_satellite/13_image/1_Hor/{input_file_hash}' ], shell=True)    
            os.chdir('./samples_satellite/13_image/1_Hor/')
            subprocess.run([f'Rscript pic_ggplot_13.1.R  {input_file_hash}  {chrserial_start} {chrserial_end} >null 2>&1 '], shell=True)  
            subprocess.run([f'rm {input_file_hash} ' ], shell=True)    
            subprocess.run([f'mv  {input_file_hash}.pdf  ./blocks/{one_region}____{length}.pdf'], shell=True)  
            os.chdir('../../../')                   

            #subprocess.run([f'Rscript ./samples_satellite/13_image/1_Hor/pic_ggplot_13.1.R {input_file}  {output_file} {chrserial_start} {chrserial_end}'], shell=True)    
            
            
        with Pool(processes=10) as pool:
            for i, result in enumerate(pool.imap(run_step, info_list), start=1):
                progress = (i / len(info_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()      
    if argv1=="stepall" or argv1=="step13"  or argv1=="step13.2":       
        print('Visualize monomers')
        if  os.path.exists(f'./samples_satellite/13_image/')==False:
            subprocess.run([f"mkdir ./samples_satellite/13_image/"], shell=True) 
            
        if  os.path.exists('./samples_satellite/13_image/2_monomer'): 
            subprocess.run(["rm -r ./samples_satellite/13_image/2_monomer"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/13_image/2_monomer"], shell=True)  
        subprocess.run(["mkdir ./samples_satellite/13_image/2_monomer/blocks"], shell=True)
        
        print('Loading regions')
        info_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                info_list.append([new_name,eachline_arr[3],eachline_arr[4]])

        R_txt=r'''
library(ggplot2)
library(dplyr)
# Get all command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Print arguments (for debugging)
print(args)
input_file_path=args[1]
chrserial_start=as.numeric(args[2])
chrserial_end=as.numeric(args[3])
print(input_file_path)
print(chrserial_start)
print(chrserial_end)
### Monomer
print("")
{input_file=read.table(input_file_path, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')

# Filter data
filtered_input <- input_file %>% 
  #filter(pos_length > 300 ) %>%
  mutate(
    #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
    category = case_when(
      monomer_type == 'CEN107' ~ "CEN107",
      monomer_type == 'CEN107-like' ~ "CEN107-like",
      monomer_type == 'CEN135' ~ "CEN135",
      monomer_type == 'CEN135-like' ~ "CEN135-like",
      monomer_type == 'CEN79' ~ "CEN79",
      monomer_type == 'CEN79-like' ~ "CEN79-like",
      monomer_type == 'CEN51' ~ "CEN51",
      monomer_type == 'CEN51-like' ~ "CEN51-like",      
      monomer_type == 'other' ~ "other_monomer",
      TRUE ~ "other_monomer"
    )
  )
  
# Define color values
color_values <- c(
  "CEN107" = "#066292",
  "CEN107-like" = "#85d2fa",
  "CEN135" = "#51B33F",
  "CEN135-like" = "#bbe4b4", 
  "CEN79" = "#ff9900",
  "CEN79-like" = "#ffd699", 
  "CEN51" = "#ff0000",
  "CEN51-like" = "#ffb3b3",   
  "other_monomer" = "black" 
)

# Create plot object
p = ggplot() 
p =p+  geom_rect(aes(xmin = chrserial_start,xmax = chrserial_end,ymin = 5,ymax = 10),fill = '#6A356C',) 
p =p+  geom_rect(data = filtered_input,
                 aes(xmin = pos_start,xmax = pos_end, ymin = 5,ymax = 10,fill = category)) 

###################################

  #geom_vline(
  #  xintercept = c(chrserial_start , chrserial_end ),
  #  linetype = "dashed",
  #  color = "red",
  #  linewidth = 0.5
  #) 

p=p+  xlim(chrserial_start, chrserial_end) + 
  ylim(0, 200) + 
  labs(
    x = paste0(chrserial_start, '-', chrserial_end),
    y = "",
    fill = ""
  ) +
  theme_classic() +         
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_blank()
  ) +
  scale_fill_manual(values = color_values, drop = FALSE)
  
  
  # Save as PDF
  pdf(file = paste0(input_file_path, ".pdf"), width = 10 / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()

}
    '''
        with open('./samples_satellite/13_image/2_monomer/pic_ggplot_13.1.R','w',encoding='utf-8') as f:
            f.write(R_txt)

        R_legend_txt=r'''
        # Load required packages
library(ggplot2)
 
# Define color values
values <- c( 
  "CEN107" = "#066292",
  "CEN107-like" = "#85d2fa",
  "CEN135" = "#51B33F",
  "CEN135-like" = "#bbe4b4", 
  "CEN79" = "#ff9900",
  "CEN79-like" = "#ffd699",  
  "CEN51" = "#ff0000",
  "CEN51-like" = "#ffb3b3",   
  "other_monomer" = "black" 
)
 
# Create a dummy data frame for plotting
df <- data.frame(
  category = factor(names(values), levels = names(values)),  # Ensure categories are factors and in order
  dummy = 1  # Create a dummy column for plotting
)
 
# Create a plot object, only for displaying the legend
p <- ggplot(df, aes(x = dummy, y = dummy, fill = category)) +
  geom_point(shape = 22, size = 3, color = "black") +  # Use square points to display the legend
  scale_fill_manual(values = values, name = "Categories") +  # Set colors and legend title
  scale_shape_manual(values = 22) +  # Ensure all points are squares
  theme_void() +  # Hide all plot elements
  theme(legend.position = "bottom", legend.key.size = unit(1, "cm"))  # Set legend position and size
 
 
# Print the legend
  # Save as PDF
  pdf("./samples_satellite/13_image/2_monomer/legend.pdf", width = 20 / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
        '''
        with open('./samples_satellite/13_image/1_Hor/pic_ggplot_13.1_lengend.R','w',encoding='utf-8') as f:
            f.write(R_legend_txt)
        subprocess.run([f'Rscript ./samples_satellite/13_image/1_Hor/pic_ggplot_13.1_lengend.R '], shell=True)  
        
        
        import hashlib
        
        #sys.exit()
        def get_file_md5(file_path):
            """Calculate the MD5 hash of a file"""
            hash_md5 = hashlib.md5()
            with open(file_path, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    hash_md5.update(chunk)
            return hash_md5.hexdigest()
        
        def run_step(one_info):
            one_region,chrserial_start,chrserial_end=one_info
            length=abs(int(chrserial_end)-int(chrserial_start))+1
            input_file=f"./samples_satellite/12_monomer_stat/1_blocks_stat/{one_region}"
            num=0
            with open (input_file,'r') as f:
                next(f)
                for line in f:
                    num+=1
            if num<2:return False        
            input_file_hash=get_file_md5(input_file)
            subprocess.run([f'cp {input_file} ./samples_satellite/13_image/2_monomer/{input_file_hash}' ], shell=True)    
            os.chdir('./samples_satellite/13_image/2_monomer/')
            subprocess.run([f'Rscript pic_ggplot_13.1.R  {input_file_hash}  {chrserial_start} {chrserial_end} >null 2>&1 '], shell=True)  #>null 2>&1
            subprocess.run([f'rm {input_file_hash} ' ], shell=True)    
            subprocess.run([f'mv  {input_file_hash}.pdf  ./blocks/{one_region}____{length}.pdf'], shell=True)  
            os.chdir('../../../')                   

            #subprocess.run([f'Rscript ./samples_satellite/13_image/1_Hor/pic_ggplot_13.1.R {input_file}  {output_file} {chrserial_start} {chrserial_end}'], shell=True)    
            
            
        with Pool(processes=10) as pool:
            for i, result in enumerate(pool.imap(run_step, info_list), start=1):
                progress = (i / len(info_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()

###Pictures
if argv1=="stepall" or "step14" in argv1:
    subprocess.run(["mkdir -p ./samples_satellite/14_chr_image"], shell=True)     
    if argv1=="stepall" or argv1=="step14_readme":
        print('Print instructions')
        with open('./samples_satellite/14_chr_image/readme','w') as f:
            txt=r'''
            14.0————1_hor_chr_stat, 2_monomer_chr_stat, 3_block_chr_stat
            14.02————4_otherblock_chr_stat
            
            14.9————Visualize HORs and monomers!!!
            14.9a_hor2————All, including repeats with count == 2
            14.9b_hor3————Filtered for HOR-repeat < 3, commonly used!!!!
            14.9c_part23————Specify partial regions (step14.9c sample chromosome position simple/normal)
            
            14.10————Visualize genes and TEs in PN40024
            
            step14.11_PN40024————ChIP-seq visualization
            
            step14.12————Starting from CENH3 peak regions on each chromosome, calculate the nearest Vsat1-6 on all other chromosomes
            
            14.13_Chr18————Visualize Chr18 (pic1f, without HOR) ———— p1p2 respectively show all and selected parts of Chr18
            
            14.14————Visualize chromosome range, using PN40024 as an example, including ChIP and gene density
            
            14.15_ChrXX————Whole chromosome with HOR and monomer information. Although vector graphics, overly detailed features may still be inexplicably compressed, less clear than 14.9b_hor3
            step14.15p_ChrXX_unit2, plotting, unit>=2
            step14.15p_ChrXX_unit3, plotting, unit>=3
            
            14.16————Visualize chromosome range, using PN40024 as an example, without ChIP, only gene density
            
            14.17————Based on 14.15, perform chromosome-level statistics by HOR type
            
            14.18————step14.18_HOR_type_ChrXX, plot results from step14.17 on all chromosomes, i.e., from step14.15p_ChrXX_unit3 files, plus classification from 14.17
            step14.18_part_HORclass————Partial display from 14.8
            
            14.19————Statistic all HORclass, based on 14.15
            
            step14.22_partall————Display all parts at precise positions: python 2-samples_satellite9.py step14.22_partall V079.hap1 Chr19 19420000-19655000
            '''
            f.write(txt) 
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.0":       
        print('Summarize individual chromosomes')
        if  os.path.exists(f'./samples_satellite/14_image/')==False:
            subprocess.run([f"mkdir ./samples_satellite/14_chr_image/"], shell=True) 
            
        if os.path.exists('./samples_satellite/14_chr_image/1_hor_chr_stat')==True:     
            subprocess.run(["rm -r ./samples_satellite/14_chr_image/1_hor_chr_stat"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/14_chr_image/1_hor_chr_stat"], shell=True) 
        
        if os.path.exists('./samples_satellite/14_chr_image/2_monomer_chr_stat')==True:     
            subprocess.run(["rm -r ./samples_satellite/14_chr_image/2_monomer_chr_stat"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/14_chr_image/2_monomer_chr_stat"], shell=True)      
        
        if os.path.exists('./samples_satellite/14_chr_image/3_block_chr_stat')==True:     
            subprocess.run(["rm -r ./samples_satellite/14_chr_image/3_block_chr_stat"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/14_chr_image/3_block_chr_stat"], shell=True) 
        
    
        print('Loading regions')
        dict_hapchr_regionlist={}
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                hapchr=eachline_arr[0]+'___'+eachline_arr[1].replace('region_','Chr')
                if hapchr not in dict_hapchr_regionlist:
                    dict_hapchr_regionlist[hapchr]=[]
                dict_hapchr_regionlist[hapchr].append(new_name)    
        hapchr_list=list(dict_hapchr_regionlist.keys())
        hapchr_list.sort()
        with open('./samples_satellite/14_chr_image/hapchr','w') as f:
            for one in hapchr_list:
                f.write(f"{one}\n")

        
        print('Merge HOR files')
        for hapchr,regionlist in dict_hapchr_regionlist.items():
            output_file=f"./samples_satellite/14_chr_image/1_hor_chr_stat/{hapchr}"
            with open (output_file,'w') as f2:
                f2.write("mer\tsegment\tsegment_type\tsegment_type_num\tsegment_length\tHOR_type\tcirc_start\tcirc_end\tpos_start\tpos_end\tmarkserial_start\tmarkserial_end\tlayer\tfather_layer\tpos_length\tmarkserial_num\tHOR_repeat_num\tblock_name\tblock_start\tblock_end\n")
                for one_region in regionlist:
                    input_file=f"./samples_satellite/12_Hor_stat/1_blocks_stat/{one_region}"
                    block_start=one_region.split(":")[-1].split('-')[0]
                    block_end=one_region.split(":")[-1].split('-')[1]
                    with open (input_file,'r') as f:
                        next(f)
                        for line in f:
                            f2.write(f"{line.strip()}\t{one_region}\t{block_start}\t{block_end}\n")
         
        print('Merge monomer files')
        for hapchr,regionlist in dict_hapchr_regionlist.items():
            output_file=f"./samples_satellite/14_chr_image/2_monomer_chr_stat/{hapchr}"
            with open (output_file,'w') as f2:
                f2.write("serial\tcirc_start\tcirc_end\tpos_start\tpos_end\tlength\tmonomer_type\tblock_name\tblock_start\tblock_end\n")
                for one_region in regionlist:
                    input_file=f"./samples_satellite/12_monomer_stat/1_blocks_stat/{one_region}"
                    block_start=one_region.split(":")[-1].split('-')[0]
                    block_end=one_region.split(":")[-1].split('-')[1]                    
                    with open (input_file,'r') as f:
                        next(f)
                        for line in f:
                            f2.write(f"{line.strip()}\t{one_region}\t{block_start}\t{block_end}\n")

        print('Generate block position files for each hap___chr')
        for hapchr,regionlist in dict_hapchr_regionlist.items():
            output_file=f"./samples_satellite/14_chr_image/3_block_chr_stat/{hapchr}"
            with open (output_file,'w') as f2:
                f2.write("hapchr\tblock_name\tblock_start\tblock_end\n")
                for one_region in regionlist:
                    block_start=one_region.split(":")[-1].split('-')[0]
                    block_end=one_region.split(":")[-1].split('-')[1]                    
                    f2.write(f"{hapchr}\t{one_region}\t{block_start}\t{block_end}\n")
                    
        print('Generate block position files for other satellites')
        for hapchr,regionlist in dict_hapchr_regionlist.items():
            output_file=f"./samples_satellite/14_chr_image/3_block_chr_stat/{hapchr}"
            with open (output_file,'w') as f2:
                f2.write("hapchr\tblock_name\tblock_start\tblock_end\tstrand\n")
                for one_region in regionlist:
                    block_start=one_region.split(":")[-1].split('-')[0]
                    block_end=one_region.split(":")[-1].split('-')[1]         
                    strand=one_region.split(":")[-2][-1]
                    f2.write(f"{hapchr}\t{one_region}\t{block_start}\t{block_end}\t{strand}\n")   
                    
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.02":         
        print('Read other satellite sequences (non-Vsat1) from ./new_work_dir/chr2blast/1_blastn/sum')                     
        if os.path.exists('./samples_satellite/14_chr_image/4_otherblock_chr_stat')==True:     
            subprocess.run(["rm -r ./samples_satellite/14_chr_image/4_otherblock_chr_stat"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/14_chr_image/4_otherblock_chr_stat"], shell=True)     
        with open('./stat_plot/0-region2info','r') as f:
            dict_hapchr_regionlist={}
            #sample	chromosome	chromosome_new	centype	chr_start	chr_end	length	strand	match_percent
            for line in f:
                eachline_arr=line.strip().split('\t')
                hapchr= eachline_arr[0]+'___'+eachline_arr[2]
                if hapchr not in dict_hapchr_regionlist: dict_hapchr_regionlist[hapchr]=[]
                dict_hapchr_regionlist[hapchr].append(eachline_arr)
        
        hapchr_used=[]        
        with open('./samples_satellite/14_chr_image/hapchr','r') as f:
            for line in f:
                eachline=line.strip()
                if len(eachline)>0:hapchr_used.append(eachline)
                
                
        for hapchr in hapchr_used:
            if hapchr not in dict_hapchr_regionlist: print(f'error\t{hapchr}');sys.exit()
            regionlist = dict_hapchr_regionlist[hapchr]
            output_file=f"./samples_satellite/14_chr_image/4_otherblock_chr_stat/{hapchr}"
            with open (output_file,'w') as f2:
                f2.write("hapchr\tcentype\tchr_start\tchr_end\tlength\tstrand\tmatch_percent\n")
                for one_region_infos in regionlist:
                    sample,chromosome,chromosome_new,centype,chr_start,chr_end,length,strand,match_percent=one_region_infos
                    f2.write(f"{hapchr}\t{centype}\t{chr_start}\t{chr_end}\t{length}\t{strand}\t{match_percent}\n")
    
    
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.9" or argv1=="step14.9a_hor2":       
        print('Visualizing SatHOR/monomer/strand')
        if  os.path.exists('./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand'): 
            subprocess.run(["rm -r ./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand"], shell=True)  
        subprocess.run(["mkdir ./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand/blocks"], shell=True)

        
        print('Loading regions')
        dict_hapchr_regionlist={}
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                #if new_name!='PN40024:region_3-:13509958-14079218':continue  #PN40024:region_1-:15125963-15154005
                hapchr=eachline_arr[0]+'___'+eachline_arr[1].replace('region_','Chr')
                #print(hapchr)
                #if hapchr!='PN40024___Chr3':continue
                if hapchr not in dict_hapchr_regionlist:
                    dict_hapchr_regionlist[hapchr]=[]
                dict_hapchr_regionlist[hapchr].append(new_name)  
        hapchr_list=list(dict_hapchr_regionlist.keys() )
        hapchr_list.sort()
        print("hapchr_list:",len(hapchr_list))
        
        
        #############
        R_legend_txt=r'''
        # Load required packages
library(ggplot2)
 
# Define color values
values <- c( 
    "2" = "#b5e48c",
    "3" = "#ff8566",
    "4" = "#ECD678",
    "5" = "#ade8f4",
    "6" = "#cc99ff",
    "7" = "#7f9fcc",  
    "8" = "#ccccff",
    "9" = "#ffccff",
    "10" = "#BAD694",   
    "11" = "#bcdcd6",  
    "12" = "#f49cbb",   
    "13" = "#DDDD7A",   
    "14" = "#bbdefb", 
    "15" = "#339933",  
    "16" = "#ffb4a2",   
    "17" = "#3385ff",  
    "18" = "#00b359",  
    "19" = "#45ABC3", 
    "20" = "#e5e600",
    "21-30"= "#ff80df",
    "30+"= "#cc3300",
    "Other"="white",
    "plus"="red",
    "minus"="blue",
    "cen_66"="#ff6666",
    "cen_103"="#6666ff",
    "cen_107"="#006666",
    "cen_191"="black",
    "other_satellite"="black",
        "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699",  
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3",     
    "other_monomer" = "black"
)
 
# Create a dummy data frame for plotting
df <- data.frame(
  category = factor(names(values), levels = names(values)),  # Ensure categories are factors and ordered
  dummy = 1  # Create a dummy column for plotting
)
 
# Create a plot object, only for displaying the legend
p <- ggplot(df, aes(x = dummy, y = dummy, fill = category)) +
  geom_point(shape = 22, size = 3, color = "black") +  # Use square points to display the legend
  scale_fill_manual(values = values, name = "Categories") +  # Set colors and legend title
  scale_shape_manual(values = 22) +  # Ensure all points are square
  theme_void() +  # Hide all plot elements
  theme(legend.position = "bottom", legend.key.size = unit(1, "cm"))  # Set legend position and size
 
 
# Print the legend
  # Save as PDF
  pdf("./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand/legend.pdf", width = 20 / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
        '''
        with open('./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand/pic_ggplot_14.9_lengend.R','w',encoding='utf-8') as f:
            f.write(R_legend_txt)
        subprocess.run([f'Rscript ./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand/pic_ggplot_14.9_lengend.R '], shell=True)  
        ###########
        R_txt=r'''
library(ggplot2)
library(dplyr)
# Get all command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Print arguments (for debugging)
print(args)
input_file_name=args[1]
input_file1=paste0(input_file_name,'.hor')
input_file2=paste0(input_file_name,'.monomer')
input_file3=paste0(input_file_name,'.block')
input_file4=paste0(input_file_name,'.otherblock')
print(input_file1)
print(input_file2)
print(input_file3)
print(input_file4)


### Monomer
print("")
{
  input_file1=read.table(input_file1, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file2=read.table(input_file2, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file3=read.table(input_file3, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file4=read.table(input_file4, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  
  min3 <- min(min(input_file3$block_start, na.rm = TRUE),min(input_file3$block_end, na.rm = TRUE))
  max3 <- max(max(input_file3$block_start, na.rm = TRUE),max(input_file3$block_end, na.rm = TRUE))
  min4 <- min(min(input_file4$chr_start, na.rm = TRUE),min(input_file4$chr_end, na.rm = TRUE))
  max4 <- max(max(input_file4$chr_start, na.rm = TRUE),max(input_file4$chr_end, na.rm = TRUE))
  min_all=min(min3,min4)
  max_all=max(max3,max4)
  print(min_all)
  print(max_all)
  
  # Filter data
  filtered_input1 <- input_file1 %>% 
    filter(pos_length > 300 ) %>%
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        mer == 2 ~ "2",
        mer == 3 ~ "3",
        mer == 4 ~ "4",
        mer == 5 ~ "5",
        mer == 6 ~ "6",
        mer == 7 ~ "7",
        mer == 8 ~ "8",
        mer == 9 ~ "9",
        mer == 10 ~ "10",
        mer == 11 ~ "11",
        mer == 12 ~ "12",
        mer == 13 ~ "13",
        mer == 14 ~ "14",
        mer == 15 ~ "15",
        mer == 16 ~ "16",
        mer == 17 ~ "17",
        mer == 18 ~ "18",
        mer == 19 ~ "19",
        mer == 20 ~ "20",
        mer>20 & mer<=30 ~ "21-30",
        mer>30 ~ "30+",
        TRUE ~ "Other"
      )
    )
  
  
  B1 <- filtered_input1 %>% filter(father_layer == 1)
  B2 <- filtered_input1 %>% filter(father_layer == 2)
  B3 <- filtered_input1 %>% filter(father_layer == 3)
  
  
  
  # Define color values
  color_values <- c(
    "2" = "#b5e48c",
    "3" = "#ff8566",
    "4" = "#ECD678",
    "5" = "#ade8f4",
    "6" = "#cc99ff",
    "7" = "#7f9fcc",  
    "8" = "#ccccff",
    "9" = "#ffccff",
    "10" = "#BAD694",   
    "11" = "#bcdcd6",  
    "12" = "#f49cbb",   
    "13" = "#DDDD7A",   
    "14" = "#bbdefb", 
    "15" = "#339933",  
    "16" = "#ffb4a2",   
    "17" = "#3385ff",  
    "18" = "#00b359",  
    "19" = "#45ABC3", 
    "20" = "#e5e600",
    "21-30"= "#ff80df",
    "30+"= "#cc3300",
    "Other"="white",
    "plus"="red",
    "minus"="blue",
    "cen_66"="#ff6666",
    "cen_103"="#6666ff",
    "cen_107"="#006666",
    "cen_191"="black",
    "other_satellite"="black",
    
    "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699",  
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3",       
    "other_monomer" = "black"
  )
    filtered_input2 <- input_file2 %>% 
    mutate(
      category = case_when(
        monomer_type == 'CEN107' ~ "CEN107",
        monomer_type == 'CEN107-like' ~ "CEN107-like",
        monomer_type == 'CEN135' ~ "CEN135",
        monomer_type == 'CEN135-like' ~ "CEN135-like",
        monomer_type == 'CEN79' ~ "CEN79",
        monomer_type == 'CEN79-like' ~ "CEN79-like",
        monomer_type == 'CEN51' ~ "CEN51",
        monomer_type == 'CEN51-like' ~ "CEN51-like",        
        monomer_type == 'other' ~ "other_monomer",
        TRUE ~ "other_monomer"
      )
    )
    
  filtered_input3 <- input_file3 %>% 
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        strand == '+' ~ "plus",
        strand == '-' ~ "minus",
        TRUE ~ "Other"
      )
    )

  filtered_input4 <- input_file4 %>% 
    #filter(match_percent > 0.5 ) %>%
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        centype == 'cen_66' ~ "cen_66",
        centype == 'cen_103' ~ "cen_103",
        centype == 'cen_107' ~ "cen_107",
        centype == 'cen_191' ~ "cen_191",
        TRUE ~ "other_satellite"
      )
    )  
  filtered_input4 <- filtered_input4 %>% 
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category2 = case_when(
        strand == '+' ~ "plus",
        strand == '-' ~ "minus",
        TRUE ~ "Other"
      )
    )  
  # Create plot object
  p = ggplot() 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = 95,ymax = 100),fill = '#6A356C',) 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = 85,ymax = 90),fill = '#6A356C',) 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = 75,ymax = 80),fill = '#6A356C',) 
  p =p+  geom_rect(data = B3,
                   aes(xmin = pos_start,xmax = pos_end, ymin = 95,ymax = 100,fill = category)) 
  p =p+  geom_rect(data = B2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = 85,ymax = 90,fill = category)) 
  p =p+  geom_rect(data = B1,
                   aes(xmin = pos_start,xmax = pos_end, ymin = 75,ymax = 80,fill = category)) 

 ##monomer
  p =p+  geom_rect(data = filtered_input2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = 65,ymax = 70,fill = category))   
  
  ##OTHER satellite
  p =p+  geom_rect(data = filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end, ymin = 55,ymax = 60,fill = category)) 
  
  ##strand B1
  p =p+  geom_rect(data = filtered_input3,
                   aes(xmin = block_start,xmax = block_end, ymin = 45,ymax = 50,fill = category)) 
  ##strand other
  p =p+  geom_rect(data = filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end, ymin = 45,ymax = 50,fill = category2)) 
  ##name
  p =p+  annotate("text", x = (min_all+max_all)/2, y = 35, label = "strand---satellite----B1--B2--B3 ",size=4)

  ##
  ###################################
  ##Strandness
  
  
  p=p+  xlim(min_all, max_all) + 
    ylim(0, 200) + 
    labs(
      x = paste0(min_all, '-', max_all),
      y = "",
      fill = ""
    ) +
    theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  
  width_revise=max(10,(max_all-min_all)/1000000*10)
  # Save as PDF
  pdf(file = paste0(input_file_name, ".pdf"), width = width_revise / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
  
}
    '''
        with open('./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand/pic_ggplot_14.9.R','w',encoding='utf-8') as f:
            f.write(R_txt)

        

        def run_step(one_hapchr):
            input_file1=f"./samples_satellite/14_chr_image/1_hor_chr_stat/{one_hapchr}"
            input_file2=f"./samples_satellite/14_chr_image/2_monomer_chr_stat/{one_hapchr}"
            input_file3=f"./samples_satellite/14_chr_image/3_block_chr_stat/{one_hapchr}"
            input_file4=f"./samples_satellite/14_chr_image/4_otherblock_chr_stat/{one_hapchr}"
            output_file1=f'./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand/{one_hapchr}.hor'
            output_file2=f'./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand/{one_hapchr}.monomer'
            output_file3=f'./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand/{one_hapchr}.block'
            output_file4=f'./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand/{one_hapchr}.otherblock'
            subprocess.run([f'cp {input_file1} {output_file1}' ], shell=True)   
            subprocess.run([f'cp {input_file2} {output_file2}' ], shell=True) 
            subprocess.run([f'cp {input_file3} {output_file3}' ], shell=True)   
            subprocess.run([f'cp {input_file4} {output_file4}' ], shell=True)             
            os.chdir('./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand/')
            subprocess.run([f'Rscript pic_ggplot_14.9.R  {one_hapchr}  >null 2>&1 '], shell=True)  #>null 2>&1 
            subprocess.run([f'rm {one_hapchr}.hor {one_hapchr}.monomer {one_hapchr}.block {one_hapchr}.otherblock' ], shell=True)    
            subprocess.run([f'mv  {one_hapchr}.pdf  ./blocks/{one_hapchr}.pdf'], shell=True)  
            os.chdir('../../../')                   

            #subprocess.run([f'Rscript ./samples_satellite/13_image/1_Hor/pic_ggplot_13.1.R {input_file}  {output_file} {chrserial_start} {chrserial_end}'], shell=True)    
            
            
        with Pool(processes=10) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, hapchr_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(hapchr_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()              
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.9" or argv1=="step14.9b_hor3":       
        print('Visualizing SatHOR/monomer/strand')
        if  os.path.exists('./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand_hor3repeat'): 
            subprocess.run(["rm -r ./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand_hor3repeat"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand_hor3repeat"], shell=True)  
        subprocess.run(["mkdir ./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand_hor3repeat/blocks"], shell=True)

        
        print('Loading regions')
        dict_hapchr_regionlist={}
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                #if new_name!='PN40024:region_3-:13509958-14079218':continue  #PN40024:region_1-:15125963-15154005
                hapchr=eachline_arr[0]+'___'+eachline_arr[1].replace('region_','Chr')
                #print(hapchr)
                #if hapchr!='PN40024___Chr3':continue
                #if 'PN40024' !=eachline_arr[0]:continue
                if hapchr not in dict_hapchr_regionlist:
                    dict_hapchr_regionlist[hapchr]=[]
                dict_hapchr_regionlist[hapchr].append(new_name)  
        hapchr_list=list(dict_hapchr_regionlist.keys() )
        hapchr_list.sort()
        print("hapchr_list:",len(hapchr_list))
        
        
        #############
        R_legend_txt=r'''
        # Load required packages
library(ggplot2)
 
# Define color values
values <- c( 
    "2" = "#b5e48c",
    "3" = "#ff8566",
    "4" = "#ECD678",
    "5" = "#ade8f4",
    "6" = "#cc99ff",
    "7" = "#7f9fcc",  
    "8" = "#ccccff",
    "9" = "#ffccff",
    "10" = "#BAD694",   
    "11" = "#bcdcd6",  
    "12" = "#f49cbb",   
    "13" = "#DDDD7A",   
    "14" = "#bbdefb", 
    "15" = "#339933",  
    "16" = "#ffb4a2",   
    "17" = "#3385ff",  
    "18" = "#00b359",  
    "19" = "#45ABC3", 
    "20" = "#e5e600",
    "21-30"= "#ff80df",
    "30+"= "#cc3300",
    "Other"="white",
    "plus"="red",
    "minus"="blue",
    "cen_66"="#00ccff",
    "cen_103"="#ff579c",
    "cen_107"="#199200",
    "cen_191"="#9933ff",
    "cen_355"="#995c17",
    "cen_383"="#ffe600",
    #"other_satellite"="black",
    
    "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699", 
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3",    
    "other_monomer" = "black"
)
 
# Create a dummy data frame for plotting
df <- data.frame(
  category = factor(names(values), levels = names(values)),  # Ensure categories are factors and ordered
  dummy = 1  # Create a dummy column for plotting
)
 
# Create a plot object, only for displaying the legend
p <- ggplot(df, aes(x = dummy, y = dummy, fill = category)) +
  geom_point(shape = 22, size = 3, color = "black") +  # Use square points to display the legend
  scale_fill_manual(values = values, name = "Categories") +  # Set colors and legend title
  scale_shape_manual(values = 22) +  # Ensure all points are square
  theme_void() +  # Hide all plot elements
  theme(legend.position = "bottom", legend.key.size = unit(1, "cm"))  # Set legend position and size
 
 
# Print the legend
  # Save as PDF
  pdf("./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand_hor3repeat/legend.pdf", width = 20 / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
        '''
        with open('./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand_hor3repeat/pic_ggplot_14.9_lengend.R','w',encoding='utf-8') as f:
            f.write(R_legend_txt)
        subprocess.run([f'Rscript ./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand_hor3repeat/pic_ggplot_14.9_lengend.R '], shell=True)  
        ###########
        R_txt=r'''
library(ggplot2)
library(dplyr)
# Get all command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Print arguments (for debugging)
print(args)
input_file_name=args[1]
input_file1=paste0(input_file_name,'.hor')
input_file2=paste0(input_file_name,'.monomer')
input_file3=paste0(input_file_name,'.block')
input_file4=paste0(input_file_name,'.otherblock')
print(input_file1)
print(input_file2)
print(input_file3)
print(input_file4)


### Monomer
print("")
{
  input_file1=read.table(input_file1, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file2=read.table(input_file2, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file3=read.table(input_file3, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file4=read.table(input_file4, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  
  min3 <- min(min(input_file3$block_start, na.rm = TRUE),min(input_file3$block_end, na.rm = TRUE))
  max3 <- max(max(input_file3$block_start, na.rm = TRUE),max(input_file3$block_end, na.rm = TRUE))
  min4 <- min(min(input_file4$chr_start, na.rm = TRUE),min(input_file4$chr_end, na.rm = TRUE))
  max4 <- max(max(input_file4$chr_start, na.rm = TRUE),max(input_file4$chr_end, na.rm = TRUE))
  min_all=min(min3,min4)
  max_all=max(max3,max4)
  print(min_all)
  print(max_all)
  
  # Filter data
  filtered_input1 <- input_file1 %>% 
    filter(pos_length > 300 & HOR_repeat_num>=3) %>%     ############################################################## The difference from step9 is adding this filter condition in HOR identification
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        mer == 2 ~ "2",
        mer == 3 ~ "3",
        mer == 4 ~ "4",
        mer == 5 ~ "5",
        mer == 6 ~ "6",
        mer == 7 ~ "7",
        mer == 8 ~ "8",
        mer == 9 ~ "9",
        mer == 10 ~ "10",
        mer == 11 ~ "11",
        mer == 12 ~ "12",
        mer == 13 ~ "13",
        mer == 14 ~ "14",
        mer == 15 ~ "15",
        mer == 16 ~ "16",
        mer == 17 ~ "17",
        mer == 18 ~ "18",
        mer == 19 ~ "19",
        mer == 20 ~ "20",
        mer>20 & mer<=30 ~ "21-30",
        mer>30 ~ "30+",
        TRUE ~ "Other"
      )
    )
  
  
  B1 <- filtered_input1 %>% filter(father_layer == 1)
  B2 <- filtered_input1 %>% filter(father_layer == 2)
  B3 <- filtered_input1 %>% filter(father_layer == 3)
  
  
  
  # Define color values
  color_values <- c(
    "2" = "#b5e48c",
    "3" = "#ff8566",
    "4" = "#ECD678",
    "5" = "#ade8f4",
    "6" = "#cc99ff",
    "7" = "#7f9fcc",  
    "8" = "#ccccff",
    "9" = "#ffccff",
    "10" = "#BAD694",   
    "11" = "#bcdcd6",  
    "12" = "#f49cbb",   
    "13" = "#DDDD7A",   
    "14" = "#bbdefb", 
    "15" = "#339933",  
    "16" = "#ffb4a2",   
    "17" = "#3385ff",  
    "18" = "#00b359",  
    "19" = "#45ABC3", 
    "20" = "#e5e600",
    "21-30"= "#ff80df",
    "30+"= "#cc3300",
    "Other"="white",
    "plus"="red",
    "minus"="blue",
    "cen_66"="#00ccff",
    "cen_103"="#ff579c",
    "cen_107"="#199200",
    "cen_191"="#9933ff",
    "cen_355"="#995c17",
    "cen_383"="#ffe600",
    #"cen_677"="black",
    #"other_satellite"="black",
    
    "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699",  
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3",    
    "other_monomer" = "black"
  )
filtered_input2 <- input_file2 %>% 
    mutate(
      category = case_when(
        monomer_type == 'CEN107' ~ "CEN107",
        monomer_type == 'CEN107-like' ~ "CEN107-like",
        monomer_type == 'CEN135' ~ "CEN135",
        monomer_type == 'CEN135-like' ~ "CEN135-like",
        monomer_type == 'CEN79' ~ "CEN79",
        monomer_type == 'CEN79-like' ~ "CEN79-like",
        monomer_type == 'CEN51' ~ "CEN51",
        monomer_type == 'CEN51-like' ~ "CEN51-like",        
        monomer_type == 'other' ~ "other_monomer",
        TRUE ~ "other_monomer"
      )
    )
    
  filtered_input3 <- input_file3 %>% 
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        strand == '+' ~ "plus",
        strand == '-' ~ "minus",
        TRUE ~ "Other"
      )
    )

  filtered_input4 <- input_file4 %>% 
    filter(centype %in% c("cen_66","cen_103","cen_107","cen_191","cen_355","cen_383")) %>%
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        centype == 'cen_66' ~ "cen_66",
        centype == 'cen_103' ~ "cen_103",
        centype == 'cen_107' ~ "cen_107",
        centype == 'cen_191' ~ "cen_191",
        centype == 'cen_355' ~ "cen_355",
        centype == 'cen_383' ~ "cen_383",
        #TRUE ~ "other_satellite"
      )
    )  
  filtered_input4 <- filtered_input4 %>% 
    filter( centype %in% c("cen_66","cen_103","cen_107","cen_191","cen_355","cen_383") ) %>%
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category2 = case_when(
        strand == '+' ~ "plus",
        strand == '-' ~ "minus",
        TRUE ~ "Other"
      )
    )  
  # Create plot object
  p = ggplot() 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = 95,ymax = 100),fill = '#6A356C',) 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = 85,ymax = 90),fill = '#6A356C',) 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = 75,ymax = 80),fill = '#6A356C',) 
  p =p+  geom_rect(data = B3,
                   aes(xmin = pos_start,xmax = pos_end, ymin = 95,ymax = 100,fill = category)) 
  p =p+  geom_rect(data = B2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = 85,ymax = 90,fill = category)) 
  p =p+  geom_rect(data = B1,
                   aes(xmin = pos_start,xmax = pos_end, ymin = 75,ymax = 80,fill = category)) 

##monomer
  p =p+  geom_rect(data = filtered_input2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = 65,ymax = 70,fill = category))   
  
  ##OTHER satellite
  p =p+  geom_rect(data = filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end, ymin = 55,ymax = 60,fill = category)) 
  
  ##strand B1
  p =p+  geom_rect(data = filtered_input3,
                   aes(xmin = block_start,xmax = block_end, ymin = 45,ymax = 50,fill = category)) 
  ##strand other
  p =p+  geom_rect(data = filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end, ymin = 45,ymax = 50,fill = category2)) 
  ##name
  p =p+  annotate("text", x = (min_all+max_all)/2, y = 35, label = "strand---satellite----B1--B2--B3 ",size=4)

  ##
  ###################################
  ##Strandness
  
  
  p=p+  xlim(min_all, max_all) + 
    ylim(0, 200) + 
    labs(
      x = paste0(min_all, '-', max_all),
      y = "",
      fill = ""
    ) +
    theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  
  width_revise=max(10,(max_all-min_all)/1000000*10)
  # Save as PDF
  pdf(file = paste0(input_file_name, ".pdf"), width = width_revise / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
  
}
    '''
        with open('./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand_hor3repeat/pic_ggplot_14.9.R','w',encoding='utf-8') as f:
            f.write(R_txt)

        

        def run_step(one_hapchr):
            input_file1=f"./samples_satellite/14_chr_image/1_hor_chr_stat/{one_hapchr}"
            input_file2=f"./samples_satellite/14_chr_image/2_monomer_chr_stat/{one_hapchr}"
            input_file3=f"./samples_satellite/14_chr_image/3_block_chr_stat/{one_hapchr}"
            input_file4=f"./samples_satellite/14_chr_image/4_otherblock_chr_stat/{one_hapchr}"
            output_file1=f'./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand_hor3repeat/{one_hapchr}.hor'
            output_file2=f'./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand_hor3repeat/{one_hapchr}.monomer'
            output_file3=f'./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand_hor3repeat/{one_hapchr}.block'
            output_file4=f'./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand_hor3repeat/{one_hapchr}.otherblock'
            subprocess.run([f'cp {input_file1} {output_file1}' ], shell=True)   
            subprocess.run([f'cp {input_file2} {output_file2}' ], shell=True) 
            subprocess.run([f'cp {input_file3} {output_file3}' ], shell=True)   
            subprocess.run([f'cp {input_file4} {output_file4}' ], shell=True)             
            os.chdir('./samples_satellite/14_chr_image/9_sat_Hor_monomer_strand_hor3repeat/')
            subprocess.run([f'Rscript pic_ggplot_14.9.R  {one_hapchr} >null 2>&1  '], shell=True)  #>null 2>&1 
            subprocess.run([f'rm {one_hapchr}.hor {one_hapchr}.monomer {one_hapchr}.block {one_hapchr}.otherblock' ], shell=True)    
            subprocess.run([f'mv  {one_hapchr}.pdf  ./blocks/{one_hapchr}.pdf'], shell=True)  
            os.chdir('../../../')                   

            #subprocess.run([f'Rscript ./samples_satellite/13_image/1_Hor/pic_ggplot_13.1.R {input_file}  {output_file} {chrserial_start} {chrserial_end}'], shell=True)    
            
            
        with Pool(processes=10) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, hapchr_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(hapchr_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()                              
    if   argv1=="step14.9c_part23":       
        print('step14.9c sample chromosome position mode             (missing mode parameter, default mode=simple, repeat>=3)')
        print('step14.9c sample chromosome position normal(repeat>=2)')
        print('step14.9c sample chromosome position simple(repeat>=3)')
        print('Visualizing SatHOR/monomer/strand')    
        if  os.path.exists('./samples_satellite/14_chr_image/9_mono_hor_part')==False: 
            subprocess.run(["mkdir ./samples_satellite/14_chr_image/9_mono_hor_part"], shell=True)  
        sample=sys.argv[2]
        chromosome= sys.argv[3]
        pos= sys.argv[4]
        pos1=pos.split('-')[0]
        pos2=pos.split('-')[1]
        if len(sys.argv)==6:
            mode=sys.argv[5]
        else:mode='simple'    
        if mode not in ['simple','normal']:print("step14.9c sample chromosome position normal(repeat>=2)/simple(repeat>=3)");sys.exit()
        if mode=='normal': mode_num=2
        elif mode=='simple': mode_num=3
        print(f"{sample}\t{chromosome}\t{pos}\t{mode}\t")

        one_hapchr=f"{sample}___{chromosome}"
        R_txt=r'''
library(ggplot2)
library(dplyr)
# Get all command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Print arguments (for debugging)
print(args)
input_file_name=args[1]
input_file1=paste0(input_file_name,'.hor')
input_file2=paste0(input_file_name,'.monomer')
input_file3=paste0(input_file_name,'.block')
input_file4=paste0(input_file_name,'.otherblock')
print(input_file1)
print(input_file2)
print(input_file3)
print(input_file4)
pos1=as.numeric(args[2])
pos2=as.numeric(args[3])
mode_num=as.numeric(args[4])
print(pos1)
print(pos2)
### Monomer
print("")
{
  input_file1=read.table(input_file1, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file2=read.table(input_file2, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file3=read.table(input_file3, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file4=read.table(input_file4, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  
  min3 <- min(min(input_file3$block_start, na.rm = TRUE),min(input_file3$block_end, na.rm = TRUE))
  max3 <- max(max(input_file3$block_start, na.rm = TRUE),max(input_file3$block_end, na.rm = TRUE))
  min4 <- min(min(input_file4$chr_start, na.rm = TRUE),min(input_file4$chr_end, na.rm = TRUE))
  max4 <- max(max(input_file4$chr_start, na.rm = TRUE),max(input_file4$chr_end, na.rm = TRUE))
  min_all=min(min3,min4)
  max_all=max(max3,max4)
  #print(min_all)
  #print(max_all)
  
  # Filter data
  filtered_input1 <- input_file1 %>% 
  filter(pos_length > 300) %>%
    filter(HOR_repeat_num>=mode_num) %>%     ###########
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        mer == 2 ~ "2",
        mer == 3 ~ "3",
        mer == 4 ~ "4",
        mer == 5 ~ "5",
        mer == 6 ~ "6",
        mer == 7 ~ "7",
        mer == 8 ~ "8",
        mer == 9 ~ "9",
        mer == 10 ~ "10",
        mer == 11 ~ "11",
        mer == 12 ~ "12",
        mer == 13 ~ "13",
        mer == 14 ~ "14",
        mer == 15 ~ "15",
        mer == 16 ~ "16",
        mer == 17 ~ "17",
        mer == 18 ~ "18",
        mer == 19 ~ "19",
        mer == 20 ~ "20",
        mer>20 & mer<=30 ~ "21-30",
        mer>30 ~ "30+",
        TRUE ~ "Other"
      )
    )
  
  
  B1 <- filtered_input1 %>% filter(father_layer == 1)
  B2 <- filtered_input1 %>% filter(father_layer == 2)
  B3 <- filtered_input1 %>% filter(father_layer == 3)
  
  
  
  # Define color values
  color_values <- c(
    "2" = "#b5e48c",
    "3" = "#ff8566",
    "4" = "#ECD678",
    "5" = "#ade8f4",
    "6" = "#cc99ff",
    "7" = "#7f9fcc",  
    "8" = "#ccccff",
    "9" = "#ffccff",
    "10" = "#BAD694",   
    "11" = "#bcdcd6",  
    "12" = "#f49cbb",   
    "13" = "#DDDD7A",   
    "14" = "#bbdefb", 
    "15" = "#339933",  
    "16" = "#ffb4a2",   
    "17" = "#3385ff",  
    "18" = "#00b359",  
    "19" = "#45ABC3", 
    "20" = "#e5e600",
    "21-30"= "#ff80df",
    "30+"= "#cc3300",
    "Other"="white",
    "plus"="red",
    "minus"="blue",
    "cen_66"="#00ccff",
    "cen_103"="#ff579c",
    "cen_107"="#199200",
    "cen_191"="#9933ff",
    "cen_355"="#995c17",
    "cen_383"="#ffe600",
    #"cen_677"="black",
    #"other_satellite"="black",
    
    "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699",  
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3",        
    "other_monomer" = "black"
  )
filtered_input2 <- input_file2 %>% 
    mutate(
      category = case_when(
        monomer_type == 'CEN107' ~ "CEN107",
        monomer_type == 'CEN107-like' ~ "CEN107-like",
        monomer_type == 'CEN135' ~ "CEN135",
        monomer_type == 'CEN135-like' ~ "CEN135-like",
        monomer_type == 'CEN79' ~ "CEN79",
        monomer_type == 'CEN79-like' ~ "CEN79-like",
        monomer_type == 'CEN51' ~ "CEN51",
        monomer_type == 'CEN51-like' ~ "CEN51-like",        
        monomer_type == 'other' ~ "other_monomer",
        TRUE ~ "other_monomer"
      )
    )
    
  filtered_input3 <- input_file3 %>% 
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        strand == '+' ~ "plus",
        strand == '-' ~ "minus",
        TRUE ~ "Other"
      )
    )

  filtered_input4 <- input_file4 %>% 
    filter(centype %in% c("cen_66","cen_103","cen_107","cen_191","cen_355","cen_383")) %>%
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        centype == 'cen_66' ~ "cen_66",
        centype == 'cen_103' ~ "cen_103",
        centype == 'cen_107' ~ "cen_107",
        centype == 'cen_191' ~ "cen_191",
        centype == 'cen_355' ~ "cen_355",
        centype == 'cen_383' ~ "cen_383",
        #TRUE ~ "other_satellite"
      )
    )  
  filtered_input4 <- filtered_input4 %>% 
    filter( centype %in% c("cen_66","cen_103","cen_107","cen_191","cen_355","cen_383") ) %>%
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category2 = case_when(
        strand == '+' ~ "plus",
        strand == '-' ~ "minus",
        TRUE ~ "Other"
      )
    )  
  # Create plot object
  p = ggplot() 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = 95,ymax = 100),fill = '#6A356C',) 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = 85,ymax = 90),fill = '#6A356C',) 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = 75,ymax = 80),fill = '#6A356C',) 
  p =p+  geom_rect(data = B3,
                   aes(xmin = pos_start,xmax = pos_end, ymin = 95,ymax = 100,fill = category)) 
  p =p+  geom_rect(data = B2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = 85,ymax = 90,fill = category)) 
  p =p+  geom_rect(data = B1,
                   aes(xmin = pos_start,xmax = pos_end, ymin = 75,ymax = 80,fill = category)) 

##monomer
  p =p+  geom_rect(data = filtered_input2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = 65,ymax = 70,fill = category))   
  
  ##OTHER satellite
  p =p+  geom_rect(data = filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end, ymin = 55,ymax = 60,fill = category)) 
  
  ##strand B1
  p =p+  geom_rect(data = filtered_input3,
                   aes(xmin = block_start,xmax = block_end, ymin = 45,ymax = 50,fill = category)) 
  ##strand other
  p =p+  geom_rect(data = filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end, ymin = 45,ymax = 50,fill = category2)) 
  ##name
  #p =p+  annotate("text", x = (min_all+max_all)/2, y = 35, label = "strand---satellite----B1--B2--B3 ",size=4)

p =p+ geom_vline(
    xintercept = c(pos1, pos2),
    linetype = "dashed",  # Dashed line type
    color = "red",        # Line color
    size = 1              # Line thickness
  ) 
  ##
  ###################################
  ##Strandness
  
  
  p=p+  #xlim(pos1, pos2) + 
  coord_cartesian(xlim = c(pos1, pos2), expand = FALSE)+           ##xlim hides elements not included but spanning across; coord_cartesian displays everything
    ylim(30, 110) + 
    labs(
      x = paste0(pos1, '-', pos2),
      y = "",
      fill = ""
    ) +
    theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  
  width_revise=max(10,(pos2-pos1)/1000000*10)
  # Save as PDF
  pdf(file = paste0(input_file_name, ".pdf"), width = width_revise / 2.54+3, height = 10 / 2.54)
  print(p)
  dev.off()
  
}
    '''
        with open('./samples_satellite/14_chr_image/9_mono_hor_part/pic_ggplot_14.9c.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        input_file1=f"./samples_satellite/14_chr_image/1_hor_chr_stat/{one_hapchr}"  # {one_hapchr}
        input_file2=f"./samples_satellite/14_chr_image/2_monomer_chr_stat/{one_hapchr}"
        input_file3=f"./samples_satellite/14_chr_image/3_block_chr_stat/{one_hapchr}"
        input_file4=f"./samples_satellite/14_chr_image/4_otherblock_chr_stat/{one_hapchr}"
        output_file1=f'./samples_satellite/14_chr_image/9_mono_hor_part/{sample}.hor'
        output_file2=f'./samples_satellite/14_chr_image/9_mono_hor_part/{sample}.monomer'
        output_file3=f'./samples_satellite/14_chr_image/9_mono_hor_part/{sample}.block'
        output_file4=f'./samples_satellite/14_chr_image/9_mono_hor_part/{sample}.otherblock'
        subprocess.run([f'cp {input_file1} {output_file1}' ], shell=True)   
        subprocess.run([f'cp {input_file2} {output_file2}' ], shell=True) 
        subprocess.run([f'cp {input_file3} {output_file3}' ], shell=True)   
        subprocess.run([f'cp {input_file4} {output_file4}' ], shell=True)             
        os.chdir('./samples_satellite/14_chr_image/9_mono_hor_part/')
        subprocess.run([f'Rscript pic_ggplot_14.9c.R  {sample} {pos1} {pos2} {mode_num}'], shell=True)        
        os.chdir('../../../')
        subprocess.run([f'rm ./samples_satellite/14_chr_image/9_mono_hor_part/{sample}.hor ./samples_satellite/14_chr_image/9_mono_hor_part/{sample}.monomer ./samples_satellite/14_chr_image/9_mono_hor_part/{sample}.block ./samples_satellite/14_chr_image/9_mono_hor_part/{sample}.otherblock'], shell=True)        
        
    #The code for VHP is similar.           
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.10" or argv1=="step14.10_PN40024":
        print('Visualizing genes and TEs of PN40024')
        if  os.path.exists('./samples_satellite/14_chr_image/10_geneTE_PN40024'): 
            subprocess.run(["rm -r ./samples_satellite/14_chr_image/10_geneTE_PN40024"], shell=True)   
        subprocess.run(["mkdir ./samples_satellite/14_chr_image/10_geneTE_PN40024"], shell=True)          
        sample='PN40024'
        dict_sample_cen_pos={}
        dict_sample_cen_pos['PN40024']={}
        dict_sample_cen_pos['PN40024']['Chr1']=[13593437,18778696]
        dict_sample_cen_pos['PN40024']['Chr2']=[12314344,14598337]
        dict_sample_cen_pos['PN40024']['Chr3']=[13359195,14694588]
        dict_sample_cen_pos['PN40024']['Chr4']=[12497123,13875971]
        dict_sample_cen_pos['PN40024']['Chr5']=[12982258,15234075]
        dict_sample_cen_pos['PN40024']['Chr6']=[9480469,13216847]
        dict_sample_cen_pos['PN40024']['Chr7']=[12546375,14159463]
        dict_sample_cen_pos['PN40024']['Chr8']=[6040738,8153004]
        dict_sample_cen_pos['PN40024']['Chr9']=[14971644,16213361]
        dict_sample_cen_pos['PN40024']['Chr10']=[20055902,22365281]
        dict_sample_cen_pos['PN40024']['Chr11']=[13226414,14996086]
        dict_sample_cen_pos['PN40024']['Chr12']=[10170306,11694530]
        dict_sample_cen_pos['PN40024']['Chr13']=[11384695,12926588]
        dict_sample_cen_pos['PN40024']['Chr14']=[14913798,15899857]
        dict_sample_cen_pos['PN40024']['Chr15']=[7025767,9336095]
        dict_sample_cen_pos['PN40024']['Chr16']=[7953854,14385543]
        dict_sample_cen_pos['PN40024']['Chr17']=[14676650,15762428]
        dict_sample_cen_pos['PN40024']['Chr18']=[14676674,18327642]
        dict_sample_cen_pos['PN40024']['Chr19']=[15528645,19536633]
        #EDTA_result_GFF=f'./new_work_dir/chr2EDTA/1_EDTA/{sample}/{sample}.fasta.mod.EDTA.TEanno.gff3'
        EDTA_result_GFF=f'./new_work_dir/chr2EDTA/1_EDTA/{sample}/{sample}.fasta.mod.EDTA.intact.gff3'
        Gene_GFF=f'./orthofinder/input/{sample}_T2T2_agat.gff'
        with open('./samples_satellite/14_chr_image/10_geneTE_PN40024/0_TEgene','w') as f2:
            f2.write(f"chromosome\ttype\tstart\tend\n")
            ##
            for chromosome,borders in dict_sample_cen_pos['PN40024'].items():
                border1,border2=borders
                f2.write(f"{chromosome}\t'border'\t{border1-1000000}\t{border1}\n")
                f2.write(f"{chromosome}\t'border'\t{border2}\t{border2+1000000}\n")
            ##
            dict_type_num={}
            with open(EDTA_result_GFF,'r') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if len(eachline_arr)!=9:continue
                    chromosome=     eachline_arr[0]
                    centromere_border=dict_sample_cen_pos[sample][chromosome]
                    border1,border2=centromere_border                
                    part_type=        eachline_arr[2]
                    part_start=       int(eachline_arr[3])
                    part_end=         int(eachline_arr[4])
                    if part_end-part_start>30000:continue
                    if part_type=="LTR_retrotransposon":                      part_type="RNA-LTR-normal"
                    elif part_type=="Gypsy_LTR_retrotransposon":              part_type="RNA-LTR-Gypsy"
                    elif part_type=="Copia_LTR_retrotransposon":              part_type="RNA-LTR-Copia"
                    #   
                    elif part_type=="Rpart_LINE_retrotransposon":             part_type="RNA-LINE-RTE"
                    elif part_type=="L1_LINE_retrotransposon":                part_type="RNA-LINE-L1"
                    #
                    elif part_type=="Mutator_TIR_transposon":                 part_type="DNA-TIR-Mutator"
                    elif part_type=="PIF_Harbinger_TIR_transposon":           part_type="DNA-TIR-PIFHarbinger"
                    elif part_type=="CACTA_TIR_transposon":                   part_type="DNA-TIR-CACTA"
                    elif part_type=="Tc1_Mariner_TIR_transposon":             part_type="DNA-TIR-Tc1Mariner"
                    elif part_type=="hAT_TIR_transposon":                     part_type="DNA-TIR-hAT"
                    #
                    elif part_type=="helitron":                               part_type="helitron"
                    #
                    elif part_type=="pararetrovirus":                         part_type="pararetrovirus"
                    #
                    elif part_type=="low_complexity":                         part_type="Simple_repeat"    
                    else:continue
                    #
                    #if border1<=part_start and part_end <= border2:
                    if part_type not in dict_type_num:
                        dict_type_num[part_type]={}
                        dict_type_num[part_type]['length']=0
                        dict_type_num[part_type]['num']=0
                    dict_type_num[part_type]['num']+=1
                    dict_type_num[part_type]['length']+=(part_end-part_start+1)
                    f2.write(f"{chromosome}\t{part_type}\t{part_start}\t{part_end}\n")

            with open(Gene_GFF,'r') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if len(eachline_arr)!=9:continue
                    chromosome=     eachline_arr[0]
                    centromere_border=dict_sample_cen_pos[sample][chromosome]
                    border1,border2=centromere_border
                    part_type=        eachline_arr[2]
                    part_start=       int(eachline_arr[3])
                    part_end=         int(eachline_arr[4])    
                    if part_type!='gene':continue
                    if part_end-part_start>30000:continue
                    #if border1<=part_start and part_end <= border2:
                    f2.write(f"{chromosome}\t{part_type}\t{part_start}\t{part_end}\n")
            
            with open("./samples_satellite/14_chr_image/10_geneTE_PN40024/0_TEgene_stat",'w')        as f:
                for parttype,dict_info in dict_type_num.items():
                    f.write(f"{parttype}\t{dict_info['num']}\t{dict_info['length']}\n")
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.10" or argv1=="step14.10p_PN40024":            
        R_txt=r'''
library(ggplot2)
library(dplyr)

print("")
{
  input_file5=read.table('0_TEgene', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  


  # Create plot object
  p = ggplot() 
  p =p+  geom_rect(data=input_file5,
                   aes(xmin = start,xmax = end,ymin = 95,ymax = 100,fill=type)) 
  # Add faceting
  p = p + facet_grid(chromosome ~ ., scales = "free", space = "free")
  
  # Save as PDF
  pdf(file = paste0('0_TEgene', ".pdf"), width = 20 / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
  
  ###########################################################################################################################################################################################################################################################################################################################################################################
    filtered_input5 <- input_file5 %>% 
    mutate(
      category = case_when(
        type == 'RNA-LTR-normal' ~ "RNA",
        type == 'RNA-LTR-Gypsy' ~ "RNA",
        type == 'RNA-LTR-Copia' ~ "RNA",
        type == 'RNA-LINE-RTE' ~ "RNA",
        type == 'RNA-LINE-L1' ~ "RNA",
        type == 'DNA-TIR-Mutator' ~ "DNA",
        type == 'DNA-TIR-PIFHarbinger' ~ "DNA",
        type == 'DNA-TIR-CACTA' ~ "DNA",
        type == 'DNA-TIR-Tc1Mariner' ~ "DNA",
        type == 'DNA-TIR-hAT' ~ "DNA",
        type == 'helitron' ~ "other",
        type == 'pararetrovirus' ~ "other",
        type == 'Simple_repeat' ~ "other",
        type == 'border' ~ "border",
        type == 'gene' ~ "gene",
        TRUE ~ "other"
      )
    ) %>%filter(category!="other")
   # Define color values
  color_values <- c(
    "RNA" = "#b5e48c",
    "DNA" = "#ff8566",
    "border" = "#ECD678",
    "gene" = "blue",
    "other" = "#ade8f4"
  )
    # Create plot object
  p = ggplot() 
  p =p+  geom_rect(data=filtered_input5,
                   aes(xmin = start,xmax = end,ymin = 95,ymax = 100,fill=category)) 
  # Add faceting
  p = p + facet_grid(chromosome ~ ., scales = "free", space = "free") +
    scale_fill_manual(values = color_values, drop = FALSE)
  
  # Save as PDF
  pdf(file = paste0('0_TEgene_simple', ".pdf"), width = 20 / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
  
  
  
  ###########################################################################################################################################################################################################################################################################################################################################################################
  filtered_input5 <- input_file5 %>% 
    mutate(
      category = case_when(
        type == 'RNA-LTR-Gypsy' ~ "Gypsy",
        type == 'RNA-LTR-Copia' ~ "Copia",
        type == 'border' ~ "border",
        TRUE ~ "other"
      )
    ) %>%filter(category!="other")
   # Define color values
  color_values <- c(
    "Gypsy" = "#47476b",
    "Copia" = "#b30059",
    "border" = "#ECD678",
    "gene" = "#669999",
    "other" = "#ade8f4"
  )
    filtered_input5_2 <- input_file5 %>% 
    mutate(
      category = case_when(
        type == 'gene' ~ "gene",
        TRUE ~ "other"
      )
    ) %>%filter(category!="other")
  # Create plot object
  p = ggplot() 
  p =p+  geom_rect(data=filtered_input5,
                   aes(xmin = start,xmax = end,ymin = 95,ymax = 100,fill=category)) 
  p =p+  geom_rect(data=filtered_input5_2,
                   aes(xmin = start,xmax = end,ymin = 97,ymax = 98,fill=category))                    
  # Add faceting
  p = p + facet_grid(chromosome ~ ., scales = "free", space = "free") +theme(strip.text.y = element_text(angle = 90, size = 4))  +
    scale_fill_manual(values = color_values, drop = FALSE)    +theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    ) + 
          
  
  # Save as PDF
  pdf(file = paste0('0_TEgene_simple_only2LTR', ".pdf"), width = 20 / 2.54, height = 20 / 2.54)
  print(p)
  dev.off()
}
    '''
        with open('./samples_satellite/14_chr_image/10_geneTE_PN40024/pic_ggplot_14.10.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = "./samples_satellite/14_chr_image/10_geneTE_PN40024"
        os.chdir(new_directory)
        subprocess.run(['Rscript pic_ggplot_14.10.R'], shell=True)    
        os.chdir('../../../')
        
    #The code for VHP is similar.    
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.11" or argv1=="step14.11_PN40024":    
        if  os.path.exists('./samples_satellite/14_chr_image/11_Chip')==False: 
            subprocess.run(["mkdir ./samples_satellite/14_chr_image/11_Chip"], shell=True)      
        print('Calculating distances and sizes of the 3 closest Vsat1-6 to ChIP-seq peaks')
        print('macs3 callpeak -t {chipseq_ip_bam} -c {chipseq_in_bam} -f BAM -g 500000000  -n name -B  --broad   --outdir {step_dir} --nomodel')
        
        print('Not marking asterisk single points in the 3D plot, consider line charts etc.')

        ##Save a boundary        
        sample='PN40024'
        dict_sample_cen_pos={}
        dict_sample_cen_pos['PN40024']={}
        dict_sample_cen_pos['PN40024']['Chr1']=[13593437,18778696]
        dict_sample_cen_pos['PN40024']['Chr2']=[12314344,14598337]
        dict_sample_cen_pos['PN40024']['Chr3']=[13359195,14694588]
        dict_sample_cen_pos['PN40024']['Chr4']=[12497123,13875971]
        dict_sample_cen_pos['PN40024']['Chr5']=[12982258,15234075]
        dict_sample_cen_pos['PN40024']['Chr6']=[9480469,13216847]
        dict_sample_cen_pos['PN40024']['Chr7']=[12546375,14159463]
        dict_sample_cen_pos['PN40024']['Chr8']=[6040738,8153004]
        dict_sample_cen_pos['PN40024']['Chr9']=[14971644,16213361]
        dict_sample_cen_pos['PN40024']['Chr10']=[20055902,22365281]
        dict_sample_cen_pos['PN40024']['Chr11']=[13226414,14996086]
        dict_sample_cen_pos['PN40024']['Chr12']=[10170306,11694530]
        dict_sample_cen_pos['PN40024']['Chr13']=[11384695,12926588]
        dict_sample_cen_pos['PN40024']['Chr14']=[14913798,15899857]
        dict_sample_cen_pos['PN40024']['Chr15']=[7025767,9336095]
        dict_sample_cen_pos['PN40024']['Chr16']=[7953854,14385543]
        dict_sample_cen_pos['PN40024']['Chr17']=[14676650,15762428]
        dict_sample_cen_pos['PN40024']['Chr18']=[14676674,18327642]
        dict_sample_cen_pos['PN40024']['Chr19']=[15528645,19536633]  
        ##
        with open('./samples_satellite/14_chr_image/11_Chip/0_PN40024','w') as f2:
            f2.write(f"chromosome\tstart\tend\tfoldchange\tqvalue\n")
            ##
            for chromosome,borders in dict_sample_cen_pos['PN40024'].items():
                border1,border2=borders
                f2.write(f"{chromosome}\t{border1-1000000}\t{border1}\t1\t0\n")
                f2.write(f"{chromosome}\t{border2}\t{border2+1000000}\t1\t0\n")       
            
            with open('./samples_satellite/14_chr_image/11_Chip/PN40024.name_peaks','r') as f:
                next(f)
                for line in  f:
                    eachline_arr=line.strip().split('\t')
                    if len(eachline_arr)!=9:continue
                    chromosome=         eachline_arr[0]
                    start=	            int(eachline_arr[1])
                    end=                int(eachline_arr[2])
                    length=             int(eachline_arr[3])
                    foldchange=         float(eachline_arr[6])
                    qvalue=             float(eachline_arr[7])
                    #if qvalue<10:continue
                    #if length<100:continue
                
                    centromere_border=dict_sample_cen_pos[sample][chromosome]
                    border1,border2=centromere_border                
                    #if border1<=start and end <= border2:
                    f2.write(f"{chromosome}\t{start}\t{end}\t{foldchange}\t{qvalue}\n")
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.11" or argv1=="step14.11p_PN40024":            
        R_txt=r'''
library(ggplot2)
library(dplyr)

print("")
{
  input_file5=read.table('0_PN40024', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  filtered_input5 <- input_file5 %>% 
    mutate(
        category = case_when(
          qvalue < 40 ~ "value1",
          qvalue >= 40 & qvalue < 100 ~ "value10",
          qvalue >= 100 ~ "value100",
          TRUE ~ "other"
      )
    )
color_values <- c(
    "value1" = "#666699",
    "value10" = "#0066cc",
    "value100" = "#cc6699"
  )

  # Create plot object
  p = ggplot() 
  p =p+ geom_rect(data = filtered_input5,
            aes(xmin = start, xmax = end, ymin = 95, ymax = 95 + foldchange, fill = category)) 
  # Add faceting
  p = p + facet_grid(chromosome ~ ., scales = "free", space = "free")+
    theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    )+
    scale_fill_manual(values = color_values, drop = FALSE)
  
  # Save as PDF
  pdf(file = paste0('0_PN40024', ".pdf"), width = 20 / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
  
 
}
    '''
        with open('./samples_satellite/14_chr_image/11_Chip/0_PN40024.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = "./samples_satellite/14_chr_image/11_Chip"
        os.chdir(new_directory)
        subprocess.run(['Rscript 0_PN40024.R'], shell=True)    
        os.chdir('../../../')

    ######################    
    ######################    
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.12" or argv1=="step14.12":    
        if  os.path.exists('./samples_satellite/14_chr_image/12_Chip_distance')==False: 
            subprocess.run(["mkdir ./samples_satellite/14_chr_image/12_Chip_distance"], shell=True)  
        print('Starting from the CENH3 peak region on each chromosome, calculate the nearest Vsat1-6 on all other chromosomes')
        with open('./samples_satellite/14_chr_image/12_Chip_distance/1_PN40024_VHPhap1','w') as f2:
            with open('./samples_satellite/14_chr_image/12_Chip_distance/1_PN40024_VHPhap1_format','w') as f5:
                f2.write(f"sample\tchromosome\tstart\tend\tfoldchange\tqvalue\tcen66\tcen103\tcen107\tcen191\tcen355\tcen383\n")
                f5.write(f"sample\tchromosome\ttype\tdistance\n")
                with open('./samples_satellite/14_chr_image/11_Chip/PN40024.hap1.name_peaks','r') as f:
                    next(f)
                    for line in  f:
                        eachline_arr=line.strip().split('\t')
                        if len(eachline_arr)!=9:continue
                        chromosome=         eachline_arr[0].replace('chr','Chr')
                        start=	            int(eachline_arr[1])
                        end=                int(eachline_arr[2])
                        length=             int(eachline_arr[3])
                        foldchange=         float(eachline_arr[6])
                        qvalue=             float(eachline_arr[7])
                        if qvalue<100:continue
                        if length<10000:continue
                        dict_type_distance={}
                        dict_type_distance['cen_66']=10000000
                        dict_type_distance['cen_103']=10000000
                        dict_type_distance['cen_107']=10000000
                        dict_type_distance['cen_191']=10000000
                        dict_type_distance['cen_355']=10000000
                        dict_type_distance['cen_383']=10000000
                        with open(f"./samples_satellite/14_chr_image/4_otherblock_chr_stat/PN40024___{chromosome}",'r') as f3:
                            next(f3)
                            for line in f3:
                                eachline_arr=line.strip().split('\t')
                                centype	=   eachline_arr[1]
                                chr_start	=   int(eachline_arr[2])
                                chr_end=   int(eachline_arr[3])
                                one_distance=min(abs(chr_start-start),abs(chr_start-end),abs(chr_end-start),abs(chr_end-end))
                                if centype not in dict_type_distance:continue
                                #if dict_type_distance[centype]=='':dict_type_distance[centype]=one_distance
                                if one_distance<dict_type_distance[centype]:dict_type_distance[centype]=one_distance
                        distance_cen_66=    dict_type_distance['cen_66']
                        distance_cen_103=    dict_type_distance['cen_103']
                        distance_cen_107=    dict_type_distance['cen_107']
                        distance_cen_191=    dict_type_distance['cen_191']
                        distance_cen_355=    dict_type_distance['cen_355']
                        distance_cen_383=    dict_type_distance['cen_383']
                        f2.write(f"PN40024\t{chromosome}\t{start}\t{end}\t{foldchange}\t{qvalue}\t{distance_cen_66}\t{distance_cen_103}\t{distance_cen_107}\t{distance_cen_191}\t{distance_cen_355}\t{distance_cen_383}\n")
                        f5.write(f"PN40024\t{chromosome}\tcen66\t{distance_cen_66}\n")
                        f5.write(f"PN40024\t{chromosome}\tcen103\t{distance_cen_103}\n")
                        f5.write(f"PN40024\t{chromosome}\tcen107\t{distance_cen_107}\n")
                        f5.write(f"PN40024\t{chromosome}\tcen191\t{distance_cen_191}\n")
                        f5.write(f"PN40024\t{chromosome}\tcen355\t{distance_cen_355}\n")
                        f5.write(f"PN40024\t{chromosome}\tcen383\t{distance_cen_383}\n")
                with open('./samples_satellite/14_chr_image/11_Chip/VHP-hap1.name_peaks','r') as f:
                    next(f)
                    for line in  f:
                        eachline_arr=line.strip().split('\t')
                        if len(eachline_arr)!=9:continue
                        chromosome=         eachline_arr[0].replace('chr','Chr')
                        start=	            int(eachline_arr[1])
                        end=                int(eachline_arr[2])
                        length=             int(eachline_arr[3])
                        foldchange=         float(eachline_arr[6])
                        qvalue=             float(eachline_arr[7])
                        dict_type_distance={}
                        dict_type_distance['cen_66']=10000000
                        dict_type_distance['cen_103']=10000000
                        dict_type_distance['cen_107']=10000000
                        dict_type_distance['cen_191']=10000000
                        dict_type_distance['cen_355']=10000000
                        dict_type_distance['cen_383']=10000000
                        with open(f"./samples_satellite/14_chr_image/4_otherblock_chr_stat/VHP-T2T.hap1___{chromosome}",'r') as f3:
                            next(f3)
                            for line in f3:
                                eachline_arr=line.strip().split('\t')
                                centype	=   eachline_arr[1]
                                chr_start	=   int(eachline_arr[2])
                                chr_end=   int(eachline_arr[3])
                                one_distance=min(abs(chr_start-start),abs(chr_start-end),abs(chr_end-start),abs(chr_end-end))
                                if centype not in dict_type_distance:continue
                                #if dict_type_distance[centype]=='':dict_type_distance[centype]=one_distance
                                if one_distance<dict_type_distance[centype]:dict_type_distance[centype]=one_distance
                        distance_cen_66=    dict_type_distance['cen_66']
                        distance_cen_103=    dict_type_distance['cen_103']
                        distance_cen_107=    dict_type_distance['cen_107']
                        distance_cen_191=    dict_type_distance['cen_191']
                        distance_cen_355=    dict_type_distance['cen_355']
                        distance_cen_383=    dict_type_distance['cen_383']                    
                        if qvalue<10:continue
                        if length<10000:continue
                        f2.write(f"VHP-T2T.hap1\t{chromosome}\t{start}\t{end}\t{foldchange}\t{qvalue}\t{distance_cen_66}\t{distance_cen_103}\t{distance_cen_107}\t{distance_cen_191}\t{distance_cen_355}\t{distance_cen_383}\n")
                        f5.write(f"VHP-T2T.hap1\t{chromosome}\tcen66\t{distance_cen_66}\n")
                        f5.write(f"VHP-T2T.hap1\t{chromosome}\tcen103\t{distance_cen_103}\n")
                        f5.write(f"VHP-T2T.hap1\t{chromosome}\tcen107\t{distance_cen_107}\n")
                        f5.write(f"VHP-T2T.hap1\t{chromosome}\tcen191\t{distance_cen_191}\n")
                        f5.write(f"VHP-T2T.hap1\t{chromosome}\tcen355\t{distance_cen_355}\n")
                        f5.write(f"VHP-T2T.hap1\t{chromosome}\tcen383\t{distance_cen_383}\n")
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.12" or argv1=="step14.12p":            
        R_txt=r'''
library(ggplot2)
library(dplyr)

print("")
{
  input_file=read.table('1_PN40024_VHPhap1_format', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  # Filter data
filtered_input_normal <- input_file %>% 
  mutate(
    category = case_when(
      distance == 10000000 ~ "NA",
      TRUE ~ "normal"
    )
  ) %>% filter(category != "NA")
 
filtered_input_NA <- input_file %>% 
  mutate(
    category = case_when(
      distance == 10000000 ~ "NA",
      TRUE ~ "normal"
    )
  ) %>% filter(category == "NA")
 
# Define color values
color_values <- c(
  "PN40024" = "#cc3399",
  "VHP-T2T.hap1" = "#009999"
)
 
# Define x-axis order
x_order <- factor(c("cen383", "cen355", "cen191", "cen103", "cen66", "cen107"), levels = c("cen383", "cen355", "cen191", "cen103", "cen66", "cen107"))
 
# Create plot object
p <- ggplot() +
  geom_boxplot(data = filtered_input_normal, aes(x = type, y = distance), outlier.shape = NA) +  # Draw boxplot, hide outliers
  geom_jitter(data = filtered_input_normal %>% filter(sample == "PN40024"), 
              aes(x = as.numeric(factor(type, levels = levels(x_order))) + 0.15, y = distance, color = sample), 
              width = 0.1, height = 0.2, size = 1, shape = 19) +  # Add jitter points for PN40024 (left side)
  geom_jitter(data = filtered_input_normal %>% filter(sample == "VHP-T2T.hap1"), 
              aes(x = as.numeric(factor(type, levels = levels(x_order))) - 0.15, y = distance, color = sample), 
              width = 0.1, height = 0.2, size = 1, shape = 19) +  # Add jitter points for VHP-T2T.hap1 (right side)
  geom_jitter(data = filtered_input_NA %>% filter(sample == "PN40024"), 
              aes(x = as.numeric(factor(type, levels = levels(x_order))) + 0.15, y = 6000000, color = sample), 
              width = 0.1, height = 0.2, size = 1, shape = 19) +  # Add jitter points for PN40024 (left side)
  geom_jitter(data = filtered_input_NA %>% filter(sample == "VHP-T2T.hap1"), 
              aes(x = as.numeric(factor(type, levels = levels(x_order))) - 0.15, y = 6000000, color = sample), 
              width = 0.1, height = 0.2, size = 1, shape = 19) +  # Add jitter points for VHP-T2T.hap1 (right side)
 
  scale_x_discrete(limits = levels(x_order)) +  # Fix x-axis order
  scale_fill_manual(values = color_values) +  # Set fill colors
  scale_color_manual(values = color_values) +  # Set border colors
  coord_flip() +  # Flip the plot horizontally
  theme_classic() +
  labs(title = "Boxplot with Jittered Points", x = "Type", y = "Distance", color = "Sample")
 
  
  # Save as PDF
  pdf(file = paste0('1_PN40024_VHPhap1_format', ".pdf"), width = 17.5 / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
  
 
}
    '''
        with open('./samples_satellite/14_chr_image/12_Chip_distance/1_PN40024_VHPhap1_format.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = "./samples_satellite/14_chr_image/12_Chip_distance"
        os.chdir(new_directory)
        subprocess.run(['Rscript 1_PN40024_VHPhap1_format.R'], shell=True)    
        os.chdir('../../../../')                                                            
    if argv1=="stepall" or argv1=="step14b":
        if os.path.exists('./samples_satellite/14_chr_image/12b_Chip_distance_multiPeaks') == False:
            subprocess.run(["mkdir ./samples_satellite/14_chr_image/12b_Chip_distance_multiPeaks"], shell=True)
        
        print('Starting from the CENH3 peak region on each chromosome (allowing multiple peaks), calculate the shortest distance to each VSat1-6, satellite size must be >10000bp')
        
        with open('./samples_satellite/14_chr_image/12b_Chip_distance_multiPeaks/1_PN40024_VHPhap1', 'w') as f2:
            with open('./samples_satellite/14_chr_image/12b_Chip_distance_multiPeaks/1_PN40024_VHPhap1_format', 'w') as f5:
                f2.write(f"sample\tchromosome\tpeak_start\tpeak_end\tfoldchange\tqvalue\tminDist_cen66\tminDist_cen103\tminDist_cen107\tminDist_cen191\tminDist_cen355\tminDist_cen383\n")
                f5.write(f"sample\tchromosome\ttype\tdistance\n")
                
                # Process PN40024
                with open('./samples_satellite/14_chr_image/11_Chip/PN40024.hap1.name_peaks', 'r') as f:
                    next(f)
                    # Collect all peaks by chromosome
                    peaks_by_chr = {}
                    for line in f:
                        eachline_arr = line.strip().split('\t')
                        if len(eachline_arr) != 9:
                            continue
                        chromosome = eachline_arr[0].replace('chr', 'Chr')
                        start = int(eachline_arr[1])
                        end = int(eachline_arr[2])
                        length = int(eachline_arr[3])
                        foldchange = float(eachline_arr[6])
                        qvalue = float(eachline_arr[7])
                        
                        if qvalue < 50:
                            continue
                        if length < 5000:
                            continue
                        
                        if chromosome not in peaks_by_chr:
                            peaks_by_chr[chromosome] = []
                        peaks_by_chr[chromosome].append({
                            'start': start,
                            'end': end,
                            'foldchange': foldchange,
                            'qvalue': qvalue
                        })
                    
                    # For each chromosome, read satellite data and calculate minimum distance
                    for chromosome, peaks in peaks_by_chr.items():
                        # Read satellite data for this chromosome
                        satellite_file = f"./samples_satellite/14_chr_image/4_otherblock_chr_stat/PN40024___{chromosome}"
                        if not os.path.exists(satellite_file):
                            print(f"Warning: {satellite_file} not found, skipping {chromosome}")
                            continue
                        
                        # Collect all satellites (size > 10000bp)
                        satellites = []
                        with open(satellite_file, 'r') as f3:
                            next(f3)
                            for line in f3:
                                eachline_arr = line.strip().split('\t')
                                if len(eachline_arr) < 4:
                                    continue
                                centype = eachline_arr[1]
                                chr_start = int(eachline_arr[2])
                                chr_end = int(eachline_arr[3])
                                sat_length = chr_end - chr_start
                                if sat_length < 10000:  # Satellite size must be >10000bp
                                    continue
                                satellites.append({
                                    'type': centype,
                                    'start': chr_start,
                                    'end': chr_end
                                })
                        
                        # For each satellite type, calculate the minimum distance from all peaks to this type of satellite
                        sat_types = ['cen_66', 'cen_103', 'cen_107', 'cen_191', 'cen_355', 'cen_383']
                        min_dist_by_type = {sat_type: 10000000 for sat_type in sat_types}
                        
                        for sat in satellites:
                            sat_type = sat['type']
                            if sat_type not in min_dist_by_type:
                                continue
                            sat_start = sat['start']
                            sat_end = sat['end']
                            
                            # Calculate the minimum distance from this satellite to all peaks
                            for peak in peaks:
                                peak_start = peak['start']
                                peak_end = peak['end']
                                one_distance = min(
                                    abs(sat_start - peak_start),
                                    abs(sat_start - peak_end),
                                    abs(sat_end - peak_start),
                                    abs(sat_end - peak_end)
                                )
                                if one_distance < min_dist_by_type[sat_type]:
                                    min_dist_by_type[sat_type] = one_distance
                        
                        # Output the results for this chromosome (one summary row)
                        # Take the maximum foldchange among all peaks as representative (or take average, here take max)
                        max_foldchange = max([p['foldchange'] for p in peaks]) if peaks else 0
                        max_qvalue = max([p['qvalue'] for p in peaks]) if peaks else 0
                        # Record the peak interval range
                        all_starts = [p['start'] for p in peaks]
                        all_ends = [p['end'] for p in peaks]
                        
                        f2.write(f"PN40024\t{chromosome}\t{min(all_starts)}\t{max(all_ends)}\t{max_foldchange}\t{max_qvalue}\t"
                                f"{min_dist_by_type['cen_66']}\t{min_dist_by_type['cen_103']}\t{min_dist_by_type['cen_107']}\t"
                                f"{min_dist_by_type['cen_191']}\t{min_dist_by_type['cen_355']}\t{min_dist_by_type['cen_383']}\n")
                        
                        # Output format file (one row per satellite type)
                        for sat_type in sat_types:
                            dist = min_dist_by_type[sat_type]
                            type_short = sat_type.replace('cen_', 'cen')
                            f5.write(f"PN40024\t{chromosome}\t{type_short}\t{dist}\n")
                
                # Process VHP-T2T.hap1
                with open('./samples_satellite/14_chr_image/11_Chip/VHP-hap1.name_peaks', 'r') as f:
                    next(f)
                    # Collect all peaks by chromosome
                    peaks_by_chr = {}
                    for line in f:
                        eachline_arr = line.strip().split('\t')
                        if len(eachline_arr) != 9:
                            continue
                        chromosome = eachline_arr[0].replace('chr', 'Chr')
                        start = int(eachline_arr[1])
                        end = int(eachline_arr[2])
                        length = int(eachline_arr[3])
                        foldchange = float(eachline_arr[6])
                        qvalue = float(eachline_arr[7])
                        
                        if qvalue < 10:  # Lower threshold for VHP
                            continue
                        if length < 10000:
                            continue
                        
                        if chromosome not in peaks_by_chr:
                            peaks_by_chr[chromosome] = []
                        peaks_by_chr[chromosome].append({
                            'start': start,
                            'end': end,
                            'foldchange': foldchange,
                            'qvalue': qvalue
                        })
                    
                    # For each chromosome, read satellite data and calculate minimum distance
                    for chromosome, peaks in peaks_by_chr.items():
                        satellite_file = f"./samples_satellite/14_chr_image/4_otherblock_chr_stat/VHP-T2T.hap1___{chromosome}"
                        if not os.path.exists(satellite_file):
                            print(f"Warning: {satellite_file} not found, skipping {chromosome}")
                            continue
                        
                        # Collect all satellites (size > 10000bp)
                        satellites = []
                        with open(satellite_file, 'r') as f3:
                            next(f3)
                            for line in f3:
                                eachline_arr = line.strip().split('\t')
                                if len(eachline_arr) < 4:
                                    continue
                                centype = eachline_arr[1]
                                chr_start = int(eachline_arr[2])
                                chr_end = int(eachline_arr[3])
                                sat_length = chr_end - chr_start
                                if sat_length < 10000:
                                    continue
                                satellites.append({
                                    'type': centype,
                                    'start': chr_start,
                                    'end': chr_end
                                })
                        
                        # For each satellite type, calculate the minimum distance from all peaks to this type of satellite
                        sat_types = ['cen_66', 'cen_103', 'cen_107', 'cen_191', 'cen_355', 'cen_383']
                        min_dist_by_type = {sat_type: 10000000 for sat_type in sat_types}
                        
                        for sat in satellites:
                            sat_type = sat['type']
                            if sat_type not in min_dist_by_type:
                                continue
                            sat_start = sat['start']
                            sat_end = sat['end']
                            
                            for peak in peaks:
                                peak_start = peak['start']
                                peak_end = peak['end']
                                one_distance = min(
                                    abs(sat_start - peak_start),
                                    abs(sat_start - peak_end),
                                    abs(sat_end - peak_start),
                                    abs(sat_end - peak_end)
                                )
                                if one_distance < min_dist_by_type[sat_type]:
                                    min_dist_by_type[sat_type] = one_distance
                        
                        # Output the results for this chromosome
                        max_foldchange = max([p['foldchange'] for p in peaks]) if peaks else 0
                        max_qvalue = max([p['qvalue'] for p in peaks]) if peaks else 0
                        all_starts = [p['start'] for p in peaks]
                        all_ends = [p['end'] for p in peaks]
                        
                        f2.write(f"VHP-T2T.hap1\t{chromosome}\t{min(all_starts)}\t{max(all_ends)}\t{max_foldchange}\t{max_qvalue}\t"
                                f"{min_dist_by_type['cen_66']}\t{min_dist_by_type['cen_103']}\t{min_dist_by_type['cen_107']}\t"
                                f"{min_dist_by_type['cen_191']}\t{min_dist_by_type['cen_355']}\t{min_dist_by_type['cen_383']}\n")
                        
                        for sat_type in sat_types:
                            dist = min_dist_by_type[sat_type]
                            type_short = sat_type.replace('cen_', 'cen')
                            f5.write(f"VHP-T2T.hap1\t{chromosome}\t{type_short}\t{dist}\n")
    if argv1=="stepall" or argv1=="step14b" or argv1 == "step14b_plot":
        R_txt = r'''
    library(ggplot2)
    library(dplyr)
    
    print("Step12b: Plotting multi-peak CENH3 to satellite distances")
    
    {
      input_file = read.table('1_PN40024_VHPhap1_format', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
      
      # Filter data
      filtered_input_normal <- input_file %>% 
        mutate(
          category = case_when(
            distance == 10000000 ~ "NA",
            TRUE ~ "normal"
          )
        ) %>% filter(category != "NA")
      
      filtered_input_NA <- input_file %>% 
        mutate(
          category = case_when(
            distance == 10000000 ~ "NA",
            TRUE ~ "normal"
          )
        ) %>% filter(category == "NA")
      
      # Define color values
      color_values <- c(
        "PN40024" = "#cc3399",
        "VHP-T2T.hap1" = "#009999"
      )
      
      # Define x-axis order
      x_order <- factor(c("cen383", "cen355", "cen191", "cen103", "cen66", "cen107"), 
                        levels = c("cen383", "cen355", "cen191", "cen103", "cen66", "cen107"))
      
      # Create plot object
      p <- ggplot() +
        geom_boxplot(data = filtered_input_normal, aes(x = type, y = distance), outlier.shape = NA) +
        geom_jitter(data = filtered_input_normal %>% filter(sample == "PN40024"), 
                    aes(x = as.numeric(factor(type, levels = levels(x_order))) + 0.15, y = distance, color = sample), 
                    width = 0.2, height = 0, size = 2, shape = 16,stroke=0) +
        geom_jitter(data = filtered_input_normal %>% filter(sample == "VHP-T2T.hap1"), 
                    aes(x = as.numeric(factor(type, levels = levels(x_order))) - 0.15, y = distance, color = sample), 
                    width = 0.2, height = 0, size = 2, shape = 16,stroke=0) +
        geom_jitter(data = filtered_input_NA %>% filter(sample == "PN40024"), 
                    aes(x = as.numeric(factor(type, levels = levels(x_order))) + 0.15, y = 6000000, color = sample), 
                    width = 0.2, height = 0, size = 2, shape = 16,stroke=0) +
        geom_jitter(data = filtered_input_NA %>% filter(sample == "VHP-T2T.hap1"), 
                    aes(x = as.numeric(factor(type, levels = levels(x_order))) - 0.15, y = 6000000, color = sample), 
                    width = 0.2, height = 0, size = 2, shape = 16,stroke=0) +
        scale_x_discrete(limits = levels(x_order)) +
        scale_fill_manual(values = color_values) +
        scale_color_manual(values = color_values) +
        coord_flip() +
        theme_classic() +
        labs(title = "Multi-peak CENH3 to Satellite Distance (min distance across all peaks per chromosome)", 
             x = "Satellite Type", y = "Distance (bp)", color = "Sample")
      
      # Save as PDF
      pdf(file = paste0('1_PN40024_VHPhap1_format_multiPeak', ".pdf"), width = 17.5 / 2.54, height = 10 / 2.54)
      print(p)
      dev.off()
    }
        '''
        
        with open('./samples_satellite/14_chr_image/12b_Chip_distance_multiPeaks/1_PN40024_VHPhap1_format.R', 'w', encoding='utf-8') as f:
            f.write(R_txt)
        
        new_directory = "./samples_satellite/14_chr_image/12b_Chip_distance_multiPeaks"
        os.chdir(new_directory)
        subprocess.run(['Rscript 1_PN40024_VHPhap1_format.R'], shell=True)
        os.chdir('../../../../')
        
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.13" or argv1=="step14.13_Chr18": 
        ##
        used_dict=    {
        'V107.hap1':100,
        'V107.hap2':99,
        'V108.hap1':98,
        'V108.hap2':97,
        'V081.hap1':96,
        'V081.hap2':95,
        'V015.hap1':94,
        'V015.hap2':93,
        'V112.hap1':90,
        'V112.hap2':89,
        'DavidiiGrape_hap1':88,
        'DavidiiGrape_hap2':87,  
        'V117.hap1':84,
        'V117.hap2':83,  
        'V120.hap1':82,
        'V120.hap2':81,   
        'PiasezkiiGrape_hap1':80,
        'PiasezkiiGrape_hap2':79,
        'V079.hap1':76,
        'V079.hap2':75,

        'V105.hap1':50,
        'V105.hap2':49,
        'V126.hap1':48,
        'V126.hap2':47,
        'WoollyGrape_hap1':44,
        'WoollyGrape_hap2':43,
        'V106.hap1':40,
        'V106.hap2':39,
        'V102.hap1':38,
        'V102.hap2':37,
        'V030.hap1':34,
        'V030.hap2':33,
        'V031.hap1':32,
        'V031.hap2':31,
        'V005.hap1':28,
        'V005.hap2':27}

        if  os.path.exists('./samples_satellite/14_chr_image/13_Chr18')==False: 
            subprocess.run(["mkdir ./samples_satellite/14_chr_image/13_Chr18"], shell=True)          
        print('Plotting Chr18 for all samples')
        print('Loading regions')
        region_list=[]
        sample_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                sample=eachline_arr[0]
                if len(eachline_arr)!=8:continue
                #new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                #region_list.append(new_name)
                if sample not in sample_list:
                    sample_list.append(sample)  
        dict_sample_info={}            
        with open('./samples_satellite/sample_info','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                #if len(eachline_arr)!=5:continue
                #Baimunage_hap1	Baimunage	 vinifera	Eurasian	Eurasian_vinifera_Baimunage	Table
                sample=     eachline_arr[0]
                #species=    eachline_arr[2]
                #group=      eachline_arr[3]
                dict_sample_info[sample]=eachline_arr
        #print(dict_sample_info)
        dict_sample_info_2={}
        with open('./samples_satellite/14_chr_image/13_Chr18/0_sum','w') as  f2:
            f2.write("serial\tsample\tgroup\tspecies\tcultivar\thapchr\tcentype\tchr_start\tchr_end\tlength\tstrand\tmatch_percent\n")
            serial=0
            for one_sample in sample_list:
                
                if one_sample not in dict_sample_info:print(f'skip {one_sample}');continue
                if os.path.exists(f'./samples_satellite/14_chr_image/4_otherblock_chr_stat/{one_sample}___Chr18')==False: print(f'skip {one_sample}');continue
                serial+=1
                if one_sample not in dict_sample_info_2: dict_sample_info_2[one_sample]=[]
                with open(f'./samples_satellite/14_chr_image/4_otherblock_chr_stat/{one_sample}___Chr18','r') as f:
                    next(f)
                    
                    for line in f:
                        eachline=line.strip()
                        
                        species=        dict_sample_info[one_sample][2]
                        group=          dict_sample_info[one_sample][3]
                        cultivar=       dict_sample_info[one_sample][4]
                        newline=f"{serial}\t{one_sample}\t{group}\t{species}\t{cultivar}\t{eachline}\n"
                        dict_sample_info_2[one_sample].append(newline)
                        f2.write(newline)
        with open('./samples_satellite/14_chr_image/13_Chr18/0_sum_used','w') as  f2:
            f2.write("pos\tserial\tsample\tgroup\tspecies\tcultivar\thapchr\tcentype\tchr_start\tchr_end\tlength\tstrand\tmatch_percent\n")
            for one_sample,pos in used_dict.items(): 
                lines=dict_sample_info_2[one_sample]
                for oneline in lines:
                    f2.write(f"{pos}\t{oneline}")
                    
        with open(f'./samples_satellite/14_chr_image/13_Chr18/0_chromosome','w') as f2:
            f2.write(f"one_sample\tpos\tstart\tend\n")
            for one_sample,pos in used_dict.items(): 
                with open(f'./new_work_dir/chr2EDTA/0_prepare/{one_sample}/{one_sample}.fasta.fai','r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        if eachline_arr[0]!='Chr18':continue
                        chr18_length=int(eachline_arr[1])
                        f2.write(f"{one_sample}\t{pos}\t1\t{chr18_length}\n")                
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.13" or argv1=="step14.13p1_Chr18":            
        R_txt=r'''
library(ggplot2)
library(dplyr)

print("")
{
  input_file4=read.table('0_sum', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  filtered_input4 <- input_file4 %>% 
    filter(centype %in% c("cen_66","cen_103","cen_107","cen_191","cen_355","cen_383","cen_6","cen_68")) %>%
    mutate(
      category = case_when(
        centype == 'cen_66' ~ "cen_66",
        centype == 'cen_103' ~ "cen_103",
        centype == 'cen_107' ~ "cen_107",
        centype == 'cen_191' ~ "cen_191",
        centype == 'cen_355' ~ "cen_355",
        centype == 'cen_383' ~ "cen_383",
        #TRUE ~ "other_satellite"
      )
    ) 
  
  
    # Extract columns 1 and 3
    name_file <- filtered_input4[, c(1, 2)]
  
  # Remove duplicate rows
  name_file <- unique(name_file)
  
  
  # Define color values
  color_values <- c(
    "cen_66"="#00ccff",
    "cen_103"="#ff579c",
    "cen_107"="#199200",
    "cen_191"="#9933ff",
    "cen_355"="#995c17",
    "cen_383"="#ffe600",
    "cen_6"="#372c97",
    "cen_68"="#eab215",    
    "other" = "black"
  )
  
  # Create plot object
  p = ggplot() 
  p =p+  geom_rect(data=filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end,ymin = serial,ymax = serial+0.5,fill=centype)) +
    
    theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  ## Labels
  p =p+   geom_text(data = name_file, aes(x = 0, y = serial, label = sample), 
                    size = 4, vjust = -0.5,hjust = 0)
  
  # Save as PDF
  pdf(file = paste0('0_sum', ".pdf"), width = 20 / 2.54, height = 100 / 2.54)
  print(p)
  dev.off()
  
}
    '''
        with open('./samples_satellite/14_chr_image/13_Chr18/pic_ggplot_14.13.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = "./samples_satellite/14_chr_image/13_Chr18"
        os.chdir(new_directory)
        subprocess.run(['Rscript pic_ggplot_14.13.R'], shell=True)    
        os.chdir('../../../')                                        
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.13" or argv1=="step14.13p2_Chr18":            
        R_txt=r'''
library(ggplot2)
library(dplyr)

input_file0=read.table('0_chromosome', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
print("")
{
  input_file4=read.table('0_sum_used', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  filtered_input4 <- input_file4 %>% 
    filter(centype %in% c("cen_66","cen_103","cen_107","cen_191","cen_355","cen_383","cen_6","cen_68")) %>%
    mutate(
      category = case_when(
        centype == 'cen_66' ~ 0.5,
        centype == 'cen_103' ~  0.5,
        centype == 'cen_107' ~  0.5,
        centype == 'cen_191' ~  0.5,
        centype == 'cen_355' ~  0.5,
        centype == 'cen_383' ~  0.5,
        centype == 'cen_6' ~  0.7,
        centype == 'cen_68' ~  0.7,
        #TRUE ~ "other_satellite"
      )
    ) 
  
  
    # Extract columns 1, 3, and 5
    name_file <- filtered_input4[, c(1,3, 5)]
  
  # Remove duplicate rows
  name_file <- unique(name_file)
  
  
  # Define color values
  color_values <- c(
    "cen_66"="#00ccff",
    "cen_103"="#ff579c",
    "cen_107"="#199200",
    "cen_191"="#9933ff",
    "cen_355"="#995c17",
    "cen_383"="#ffe600",
    "cen_6"="#372c97",
    "cen_68"="#eab215",   
    "other" = "black"
  )
  
  # Create plot object
  p = ggplot() 
  p =p+  geom_rect(data=filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end,ymin = pos,ymax = pos+category,fill=centype)) +
    
    theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  ## Labels
  #p =p+   geom_text(data = name_file, aes(x = 0, y = pos, label = species), 
  #                  size = 4, vjust = -0.5,hjust = 0)
                    
  p <- p + geom_rect(data = input_file0,
                   aes(xmin = start, xmax = end, ymin =pos, ymax = pos),
                   color = 'black',linewidth=0.1,fill='#efefef')
  # Save as PDF
  pdf(file = paste0('0_sum_used', ".pdf"), width = 10 / 2.54, height = 20 / 2.54)
  print(p)
  dev.off()
  
}
    '''
        with open('./samples_satellite/14_chr_image/13_Chr18/pic_ggplot_14.13b.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = "./samples_satellite/14_chr_image/13_Chr18"
        os.chdir(new_directory)
        subprocess.run(['Rscript pic_ggplot_14.13b.R'], shell=True)    
        os.chdir('../../../')
                    
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.14" or argv1=="step14.14_sample":       
        print('Visualizing SatHOR/monomer/strand')
        if  os.path.exists('./samples_satellite/14_chr_image/14_chr_complete')==False: 
            subprocess.run(["mkdir ./samples_satellite/14_chr_image/14_chr_complete"], shell=True)  
        
        sample_list=['PN40024']
        
        def run_step(one_sample):
            subprocess.run([f"mkdir ./samples_satellite/14_chr_image/14_chr_complete/{one_sample}"], shell=True) 
            chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
            dict_sheet_info={}
            dict_sheet_info['1_hor_chr_stat']=[]
            dict_sheet_info['2_monomer_chr_stat']=[]
            dict_sheet_info['3_block_chr_stat']=[]
            dict_sheet_info['4_otherblock_chr_stat']=[]
            for one_chr in chr_list:
                one_hapchr=one_sample+'___'+one_chr
                input_file1=f"./samples_satellite/14_chr_image/1_hor_chr_stat/{one_hapchr}"
                input_file2=f"./samples_satellite/14_chr_image/2_monomer_chr_stat/{one_hapchr}"
                input_file3=f"./samples_satellite/14_chr_image/3_block_chr_stat/{one_hapchr}"
                input_file4=f"./samples_satellite/14_chr_image/4_otherblock_chr_stat/{one_hapchr}"
                k=0
                with open(input_file1,'r') as f:
                    for line in f:
                        eachline=line.strip()
                        k+=1
                        if k==1:hor_chr_stat_head=eachline
                        else:
                            dict_sheet_info['1_hor_chr_stat'].append(f"{eachline}\t{one_chr}")
                k=0
                with open(input_file2,'r') as f:
                    for line in f:
                        eachline=line.strip()
                        k+=1
                        if k==1:monomer_chr_stat_head=eachline
                        else:
                            dict_sheet_info['2_monomer_chr_stat'].append(f"{eachline}\t{one_chr}")
                k=0
                with open(input_file3,'r') as f:
                    for line in f:
                        eachline=line.strip()
                        k+=1
                        if k==1:block_chr_stat_head=eachline
                        else:
                            dict_sheet_info['3_block_chr_stat'].append(f"{eachline}\t{one_chr}")
                k=0
                with open(input_file4,'r') as f:
                    for line in f:
                        eachline=line.strip()
                        k+=1
                        if k==1:otherblock_chr_head=eachline
                        else:
                            dict_sheet_info['4_otherblock_chr_stat'].append(f"{eachline}\t{one_chr}")
            #################################################################################                
            with open(f'./samples_satellite/14_chr_image/14_chr_complete/{one_sample}/0_hor','w') as f:
                f.write(f"{hor_chr_stat_head}\tchr_note\n")
                for one in dict_sheet_info['1_hor_chr_stat']:
                    f.write(f"{one}\n")
            with open(f'./samples_satellite/14_chr_image/14_chr_complete/{one_sample}/0_monomer','w') as f:
                f.write(f"{monomer_chr_stat_head}\tchr_note\n")
                for one in dict_sheet_info['2_monomer_chr_stat']:
                    f.write(f"{one}\n")
            with open(f'./samples_satellite/14_chr_image/14_chr_complete/{one_sample}/0_block','w') as f:
                f.write(f"{block_chr_stat_head}\tchr_note\n")
                for one in dict_sheet_info['3_block_chr_stat']:
                    f.write(f"{one}\n")
            with open(f'./samples_satellite/14_chr_image/14_chr_complete/{one_sample}/0_otherblock','w') as f:
                f.write(f"{otherblock_chr_head}\tchr_note\n")
                for one in dict_sheet_info['4_otherblock_chr_stat']:
                    f.write(f"{one}\n")
            
            info_list=[]
            with open(f'./samples_satellite/14_chr_image/14_chr_complete/{one_sample}/0_TEgene','w') as f2:
                f2.write(f"chr_note\ttype\tstart\tend\n")
                with open('./samples_satellite/14_chr_image/10_geneTE_PN40024/0_TEgene','r') as f:
                    next(f)
                    for line in f:
                        f2.write(line)
                        #######
                        eachline_arr=line.strip().split('\t')
                        if 'border' in eachline_arr[1]:continue
                        chromosome=eachline_arr[0]
                        one_type=eachline_arr[1]
                        start=int(eachline_arr[2])
                        end=int(eachline_arr[3])
                        if abs(start-end)>100000:continue
                        info_list.append([chromosome,one_type,start,end])
            
            with open(f'./samples_satellite/14_chr_image/14_chr_complete/{one_sample}/0_Chip','w') as f2:
                f2.write(f"chr_note\tstart\tend\tfoldchange\tqvalue\n")
                with open('./samples_satellite/14_chr_image/11_Chip/0_PN40024','r') as f:
                    next(f)
                    for line in f:
                        f2.write(line) 
            dict_chr_length={}            
            with open(f'./samples_satellite/14_chr_image/14_chr_complete/{one_sample}/0_chromosome','w') as f2:
                f2.write(f"chr_note\tstart\tend\n")
                with open(f'./new_work_dir/chr2EDTA/0_prepare/{one_sample}/{one_sample}.fasta.fai','r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        dict_chr_length[eachline_arr[0]]=int(eachline_arr[1])
                        f2.write(f"{eachline_arr[0]}\t1\t{eachline_arr[1]}\n")    
            #print(dict_chr_length);sys.exit()
            dict_region_info={}
            for chromosome,length in dict_chr_length.items():
                kk=1
                print(chromosome)
                while kk+499999<length:
                    kk_end=kk+499999
                    region=f"{chromosome}-{kk}-{kk_end}"
                    dict_region_info[region]={}
                    dict_region_info[region]['gene']={}
                    dict_region_info[region]['gypsy']={}
                    dict_region_info[region]['copia']={}
                    #
                    dict_region_info[region]['gene']['length']=0
                    dict_region_info[region]['gene']['num']=0
                    #
                    dict_region_info[region]['gypsy']['length']=0
                    dict_region_info[region]['gypsy']['num']=0
                    #
                    dict_region_info[region]['copia']['length']=0
                    dict_region_info[region]['copia']['num']=0                    
                    #
                    #iii=0
                    for one_info in info_list:
                        
                        #iii+=1
                        #if iii>4344:sys.exit()
                        one_chr,one_type,one_start,one_end=one_info
                        #print(one_info)
                        if one_chr!=chromosome:continue
                        delta_num=1
                        if one_type=='gene':
                            if one_start>=kk   and one_end<=kk_end:
                                delta_length=abs(one_end-one_start)+1
                            elif one_start<kk and one_end<=kk_end and one_end>=kk:                         
                                delta_length=abs(one_end-kk)+1
                                #print(delta_length)
                            elif one_start>=kk and one_end>kk_end and one_start<=kk_end:                         
                                delta_length=abs(kk_end-one_start)+1
                                
                            elif one_start<kk and one_end>kk_end:   
                                delta_length=500000
                            else:
                                delta_length=0
                                delta_num=0
                            dict_region_info[region]['gene']['length']+=delta_length    
                            dict_region_info[region]['gene']['num']+=delta_num   
                        elif one_type=='RNA-LTR-Gypsy':
                            if one_start>=kk   and one_end<=kk_end:
                                delta_length=abs(one_end-one_start)+1
                            elif one_start<kk and one_end<=kk_end and one_end>=kk:            
                                delta_length=abs(one_end-kk)+1
                            elif one_start>=kk and one_end>kk_end and one_start<=kk_end:                  
                                delta_length=abs(kk_end-one_start)+1
                            elif one_start<kk and one_end>kk_end:   
                                delta_length=500000
                            else:
                                delta_length=0
                                delta_num=0                                
                            #print(delta_length)
                            dict_region_info[region]['gypsy']['length']+=delta_length    
                            dict_region_info[region]['gypsy']['num']+=delta_num
                        elif one_type=='RNA-LTR-Copia':  
                            if one_start>=kk   and one_end<=kk_end:
                                delta_length=abs(one_end-one_start)+1
                            elif one_start<kk and one_end<=kk_end and one_end>=kk:                     
                                delta_length=abs(one_end-kk)+1
                            elif one_start>=kk and one_end>kk_end and one_start<=kk_end:                            
                                delta_length=abs(kk_end-one_start)+1
                            elif one_start<kk and one_end>kk_end:   
                                delta_length=500000
                            else:
                                delta_length=0
                                delta_num=0                                
                            dict_region_info[region]['copia']['length']+=delta_length    
                            dict_region_info[region]['copia']['num']+=delta_num 
                        
                    kk+=500000
            #print(dict_region_info)  
            #sys.exit()
            with open(f'./samples_satellite/14_chr_image/14_chr_complete/{one_sample}/0_gene_stat','w') as f2:
                f2.write('chr_note\tstart\tend\tlength\tnum\n')
                for region,dict_infos in dict_region_info.items():
                    chromosome=region.split('-')[0]
                    start=region.split('-')[1]
                    end=region.split('-')[2]
                    length=dict_infos['gene']['length']
                    number=dict_infos['gene']['num']
                    f2.write(f"{chromosome}\t{start}\t{end}\t{length}\t{number}\n")    
            with open(f'./samples_satellite/14_chr_image/14_chr_complete/{one_sample}/0_TE_stat','w') as f2:
                f2.write('chr_note\tstart\tend\tlength\tnum\n')
                for region,dict_infos in dict_region_info.items():
                    chromosome=region.split('-')[0]
                    start=region.split('-')[1]
                    end=region.split('-')[2]
                    length=dict_infos['gypsy']['length']+dict_infos['copia']['length']
                    number=dict_infos['gypsy']['num']+dict_infos['copia']['num']
                    f2.write(f"{chromosome}\t{start}\t{end}\t{length}\t{number}\n")    
                    
                                
            Plot_txt=r"""
library(ggplot2)
library(dplyr)
# Get all command line arguments

### Monomer
print("")
{
  input_file0=read.table('0_chromosome', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file1=read.table('0_hor', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file2=read.table('0_monomer', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file3=read.table('0_block', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file4=read.table('0_otherblock', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  
  input_file5=read.table('0_Chip', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file6=read.table('0_TEgene', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  
  input_file7=read.table('0_TE_stat', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file8=read.table('0_gene_stat', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  
  
  foldchange_max=max(input_file5$foldchange, na.rm = TRUE)
  input_file7_maxlen=max(input_file7$length, na.rm = TRUE)
  input_file8_maxlen=max(input_file8$length, na.rm = TRUE)
  input_file7_maxnum=max(input_file7$num, na.rm = TRUE)
  input_file8_maxnum=max(input_file8$num, na.rm = TRUE)
  
  
  min3 <- min(min(input_file3$block_start, na.rm = TRUE),min(input_file3$block_end, na.rm = TRUE))
  max3 <- max(max(input_file3$block_start, na.rm = TRUE),max(input_file3$block_end, na.rm = TRUE))
  min4 <- min(min(input_file4$chr_start, na.rm = TRUE),min(input_file4$chr_end, na.rm = TRUE))
  max4 <- max(max(input_file4$chr_start, na.rm = TRUE),max(input_file4$chr_end, na.rm = TRUE))
  min_all=min(min3,min4)
  max_all=max(max3,max4)
  print(min_all)
  print(max_all)
  
  # Filter data
  filtered_input1 <- input_file1 %>% 
    filter(pos_length > 300 & HOR_repeat_num>=3) %>%     ############################################################## The difference from step9 is adding this filter condition in HOR identification
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        mer == 2 ~ "2",
        mer == 3 ~ "3",
        mer == 4 ~ "4",
        mer == 5 ~ "5",
        mer == 6 ~ "6",
        mer == 7 ~ "7",
        mer == 8 ~ "8",
        mer == 9 ~ "9",
        mer == 10 ~ "10",
        mer == 11 ~ "11",
        mer == 12 ~ "12",
        mer == 13 ~ "13",
        mer == 14 ~ "14",
        mer == 15 ~ "15",
        mer == 16 ~ "16",
        mer == 17 ~ "17",
        mer == 18 ~ "18",
        mer == 19 ~ "19",
        mer == 20 ~ "20",
        mer>20 & mer<=30 ~ "21-30",
        mer>30 ~ "30+",
        TRUE ~ "Other"
      )
    )
  
  
  B1 <- filtered_input1 %>% filter(father_layer == 1)
  B2 <- filtered_input1 %>% filter(father_layer == 2)
  B3 <- filtered_input1 %>% filter(father_layer == 3)
  
  
  
  # Define color values
  color_values <- c(
    "2" = "#b5e48c",
    "3" = "#ff8566",
    "4" = "#ECD678",
    "5" = "#ade8f4",
    "6" = "#cc99ff",
    "7" = "#7f9fcc",  
    "8" = "#ccccff",
    "9" = "#ffccff",
    "10" = "#BAD694",   
    "11" = "#bcdcd6",  
    "12" = "#f49cbb",   
    "13" = "#DDDD7A",   
    "14" = "#bbdefb", 
    "15" = "#339933",  
    "16" = "#ffb4a2",   
    "17" = "#3385ff",  
    "18" = "#00b359",  
    "19" = "#45ABC3", 
    "20" = "#e5e600",
    "21-30"= "#ff80df",
    "30+"= "#cc3300",
    "Other"="white",
    "plus"="red",
    "minus"="blue",
    "cen_66"="#00ccff",
    "cen_103"="#ff579c",
    "cen_107"="#199200",
    "cen_191"="#9933ff",
    "cen_355"="#995c17",
    "cen_383"="#ffe600",
    #"cen_677"="black",
    #"other_satellite"="black",
    
    "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699",  
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3",  
    "other_monomer" = "black",
    "Gypsy" = "#47476b",
    "Copia" = "#b30059",
    "border" = "#ECD678",
    "gene" = "#669999",
    "other" = "#ade8f4",
        "value1" = "#666699",
    "value10" = "#0066cc",
    "value100" = "#cc6699"
  )
  


 
filtered_input2 <- input_file2 %>% 
    mutate(
      category = case_when(
        monomer_type == 'CEN107' ~ "CEN107",
        monomer_type == 'CEN107-like' ~ "CEN107-like",
        monomer_type == 'CEN135' ~ "CEN135",
        monomer_type == 'CEN135-like' ~ "CEN135-like",
        monomer_type == 'CEN79' ~ "CEN79",
        monomer_type == 'CEN79-like' ~ "CEN79-like",
        monomer_type == 'CEN51' ~ "CEN51",
        monomer_type == 'CEN51-like' ~ "CEN51-like",        
        monomer_type == 'other' ~ "other_monomer",
        TRUE ~ "other_monomer"
      )
    )
    
  filtered_input3 <- input_file3 %>% 
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        strand == '+' ~ "plus",
        strand == '-' ~ "minus",
        TRUE ~ "Other"
      )
    )

  filtered_input4 <- input_file4 %>% 
    filter(centype %in% c("cen_66","cen_103","cen_107","cen_191","cen_355","cen_383")) %>%
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        centype == 'cen_66' ~ "cen_66",
        centype == 'cen_103' ~ "cen_103",
        centype == 'cen_107' ~ "cen_107",
        centype == 'cen_191' ~ "cen_191",
        centype == 'cen_355' ~ "cen_355",
        centype == 'cen_383' ~ "cen_383",
        #TRUE ~ "other_satellite"
      )
    )  
  filtered_input4 <- filtered_input4 %>% 
    filter( centype %in% c("cen_66","cen_103","cen_107","cen_191","cen_355","cen_383") ) %>%
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category2 = case_when(
        strand == '+' ~ "plus",
        strand == '-' ~ "minus",
        TRUE ~ "Other"
      )
    )  
    
    

  filtered_input5 <- input_file5 %>% 
    mutate(
        category = case_when(
          qvalue < 40 ~ "value1",
          qvalue >= 40 & qvalue < 100 ~ "value10",
          qvalue >= 100 ~ "value100",
          TRUE ~ "other"
      )
    )

    filtered_input6 <- input_file6 %>% 
    mutate(
      category = case_when(
        type == 'RNA-LTR-normal' ~ "RNA",
        type == 'RNA-LTR-Gypsy' ~ "RNA",
        type == 'RNA-LTR-Copia' ~ "RNA",
        type == 'RNA-LINE-RTE' ~ "RNA",
        type == 'RNA-LINE-L1' ~ "RNA",
        type == 'DNA-TIR-Mutator' ~ "DNA",
        type == 'DNA-TIR-PIFHarbinger' ~ "DNA",
        type == 'DNA-TIR-CACTA' ~ "DNA",
        type == 'DNA-TIR-Tc1Mariner' ~ "DNA",
        type == 'DNA-TIR-hAT' ~ "DNA",
        type == 'helitron' ~ "other",
        type == 'pararetrovirus' ~ "other",
        type == 'Simple_repeat' ~ "other",
        type == 'border' ~ "border",
        type == 'gene' ~ "gene",
        TRUE ~ "other"
      )
    ) %>%filter(category!="other")
    filtered_input6_2 <- input_file6 %>% 
    mutate(
      category = case_when(
        type == 'gene' ~ "gene",
        TRUE ~ "other"
      )
    ) %>%filter(category!="other")
  # Create plot object
      
  # Create plot object
  p = ggplot() 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = 95,ymax = 100),fill = '#6A356C',) 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = 85,ymax = 90),fill = '#6A356C',) 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = 75,ymax = 80),fill = '#6A356C',) 
  p =p+  geom_rect(data = B3,
                   aes(xmin = pos_start,xmax = pos_end, ymin = 95,ymax = 100,fill = category)) 
  p =p+  geom_rect(data = B2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = 85,ymax = 90,fill = category)) 
  p =p+  geom_rect(data = B1,
                   aes(xmin = pos_start,xmax = pos_end, ymin = 75,ymax = 80,fill = category)) 

##monomer
  p =p+  geom_rect(data = filtered_input2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = 65,ymax = 70,fill = category))   
  
  ##OTHER satellite
  p =p+  geom_rect(data = filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end, ymin = 55,ymax = 60,fill = category)) 
  
  ##strand B1
  p =p+  geom_rect(data = filtered_input3,
                   aes(xmin = block_start,xmax = block_end, ymin = 45,ymax = 50,fill = category)) 
  ##strand other
  p =p+  geom_rect(data = filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end, ymin = 45,ymax = 50,fill = category2)) 
  ##name
  #p =p+  annotate("text", x = (min_all+max_all)/2, y = 35, label = "strand---satellite----B1--B2--B3 ",size=4)

  ##
  ###################################
  ##Strandness

  p =p+ geom_rect(data = filtered_input5,
            aes(xmin = start, xmax = end, ymin = 105, ymax = 105 + foldchange/foldchange_max*20, fill = category)) 
            
  p =p+  geom_rect(data=filtered_input6,
                   aes(xmin = start,xmax = end,ymin = 25,ymax = 30,fill=category)) 
  p =p+  geom_rect(data=filtered_input6_2,
                   aes(xmin = start,xmax = end,ymin = 27,ymax = 28,fill=category))      
  p <- p + geom_rect(data = input_file0,
                   aes(xmin = start, xmax = end, ymin = 15, ymax = 20),
                   fill = 'black', color = 'black')
                   
  p <- p + geom_rect(data = input_file7,
                   aes(xmin = start, xmax = end, ymin = -50, ymax = -50+num/input_file7_maxnum*50),
                   fill = 'grey', color = 'black')
                   
  p <- p + geom_rect(data = input_file8,
                   aes(xmin = start, xmax = end, ymin = -100, ymax = -100+length/input_file8_maxlen*50),
                   fill = '#cceeff', color = 'black')                   
  #p=p+  xlim(min_all, max_all) + 
   p=p+   ylim(-150, 150) + 
    labs(
      x = paste0(min_all, '-', max_all),
      y = "",
      fill = ""
    ) +
    facet_grid(chr_note ~ ., scales = "free", space = "free")+
    theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  
  #width_revise=max(10,(max_all-min_all)/1000000*10)
  # Save as PDF
  pdf(file = paste0('Plot1', ".pdf"), width = 200 / 2.54, height = 200 / 2.54)
  print(p)
  dev.off()
  
}            
            """
            with open(f'./samples_satellite/14_chr_image/14_chr_complete/{one_sample}/Plot1.R','w') as f:
                f.write(Plot_txt)
            os.chdir(f'./samples_satellite/14_chr_image/14_chr_complete/{one_sample}')
            subprocess.run([f'Rscript Plot1.R  '], shell=True)  #>null 2>&1 
            
            
            Plot_txt=r"""
library(ggplot2)
library(dplyr)
# Get all command line arguments

### Monomer
print("")
{
  input_file0=read.table('0_chromosome', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file1=read.table('0_hor', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file2=read.table('0_monomer', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file3=read.table('0_block', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file4=read.table('0_otherblock', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  
  input_file5=read.table('0_Chip', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file6=read.table('0_TEgene', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  
  input_file7=read.table('0_TE_stat', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file8=read.table('0_gene_stat', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  
  
  foldchange_max=max(input_file5$foldchange, na.rm = TRUE)
  input_file7_maxlen=max(input_file7$length, na.rm = TRUE)
  input_file8_maxlen=max(input_file8$length, na.rm = TRUE)
  input_file7_maxnum=max(input_file7$num, na.rm = TRUE)
  input_file8_maxnum=max(input_file8$num, na.rm = TRUE)
  
  
  min3 <- min(min(input_file3$block_start, na.rm = TRUE),min(input_file3$block_end, na.rm = TRUE))
  max3 <- max(max(input_file3$block_start, na.rm = TRUE),max(input_file3$block_end, na.rm = TRUE))
  min4 <- min(min(input_file4$chr_start, na.rm = TRUE),min(input_file4$chr_end, na.rm = TRUE))
  max4 <- max(max(input_file4$chr_start, na.rm = TRUE),max(input_file4$chr_end, na.rm = TRUE))
  min_all=min(min3,min4)
  max_all=max(max3,max4)
  print(min_all)
  print(max_all)
  
  # Filter data
  filtered_input1 <- input_file1 %>% 
    filter(pos_length > 300 & HOR_repeat_num>=3) %>%     ############################################################## The difference from step9 is adding this filter condition in HOR identification
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        mer == 2 ~ "2",
        mer == 3 ~ "3",
        mer == 4 ~ "4",
        mer == 5 ~ "5",
        mer == 6 ~ "6",
        mer == 7 ~ "7",
        mer == 8 ~ "8",
        mer == 9 ~ "9",
        mer == 10 ~ "10",
        mer == 11 ~ "11",
        mer == 12 ~ "12",
        mer == 13 ~ "13",
        mer == 14 ~ "14",
        mer == 15 ~ "15",
        mer == 16 ~ "16",
        mer == 17 ~ "17",
        mer == 18 ~ "18",
        mer == 19 ~ "19",
        mer == 20 ~ "20",
        mer>20 & mer<=30 ~ "21-30",
        mer>30 ~ "30+",
        TRUE ~ "Other"
      )
    )
  
  
  B1 <- filtered_input1 %>% filter(father_layer == 1)
  B2 <- filtered_input1 %>% filter(father_layer == 2)
  B3 <- filtered_input1 %>% filter(father_layer == 3)
  
  
  
  # Define color values
  color_values <- c(
    "2" = "#b5e48c",
    "3" = "#ff8566",
    "4" = "#ECD678",
    "5" = "#ade8f4",
    "6" = "#cc99ff",
    "7" = "#7f9fcc",  
    "8" = "#ccccff",
    "9" = "#ffccff",
    "10" = "#BAD694",   
    "11" = "#bcdcd6",  
    "12" = "#f49cbb",   
    "13" = "#DDDD7A",   
    "14" = "#bbdefb", 
    "15" = "#339933",  
    "16" = "#ffb4a2",   
    "17" = "#3385ff",  
    "18" = "#00b359",  
    "19" = "#45ABC3", 
    "20" = "#e5e600",
    "21-30"= "#ff80df",
    "30+"= "#cc3300",
    "Other"="white",
    "plus"="red",
    "minus"="blue",
    "cen_66"="#00ccff",
    "cen_103"="#ff579c",
    "cen_107"="#199200",
    "cen_191"="#9933ff",
    "cen_355"="#995c17",
    "cen_383"="#ffe600",
    #"cen_677"="black",
    #"other_satellite"="black",
    
    "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699",  
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3",      
    "other_monomer" = "black",
    "Gypsy" = "#47476b",
    "Copia" = "#b30059",
    "border" = "#ECD678",
    "gene" = "#669999",
    "other" = "#ade8f4",
        "value1" = "#666699",
    "value10" = "#0066cc",
    "value100" = "#cc6699"
  )
  


 
filtered_input2 <- input_file2 %>% 
    mutate(
      category = case_when(
        monomer_type == 'CEN107' ~ "CEN107",
        monomer_type == 'CEN107-like' ~ "CEN107-like",
        monomer_type == 'CEN135' ~ "CEN135",
        monomer_type == 'CEN135-like' ~ "CEN135-like",
        monomer_type == 'CEN79' ~ "CEN79",
        monomer_type == 'CEN79-like' ~ "CEN79-like",
        monomer_type == 'CEN51' ~ "CEN51",
        monomer_type == 'CEN51-like' ~ "CEN51-like",        
        monomer_type == 'other' ~ "other_monomer",
        TRUE ~ "other_monomer"
      )
    )
    
  filtered_input3 <- input_file3 %>% 
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        strand == '+' ~ "plus",
        strand == '-' ~ "minus",
        TRUE ~ "Other"
      )
    )

  filtered_input4 <- input_file4 %>% 
    filter(centype %in% c("cen_66","cen_103","cen_107","cen_191","cen_355","cen_383")) %>%
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        centype == 'cen_66' ~ "cen_66",
        centype == 'cen_103' ~ "cen_103",
        centype == 'cen_107' ~ "cen_107",
        centype == 'cen_191' ~ "cen_191",
        centype == 'cen_355' ~ "cen_355",
        centype == 'cen_383' ~ "cen_383",
        #TRUE ~ "other_satellite"
      )
    )  
  filtered_input4 <- filtered_input4 %>% 
    filter( centype %in% c("cen_66","cen_103","cen_107","cen_191","cen_355","cen_383") ) %>%
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category2 = case_when(
        strand == '+' ~ "plus",
        strand == '-' ~ "minus",
        TRUE ~ "Other"
      )
    )  
    
    

  filtered_input5 <- input_file5 %>% 
    mutate(
        category = case_when(
          qvalue < 40 ~ "value1",
          qvalue >= 40 & qvalue < 100 ~ "value10",
          qvalue >= 100 ~ "value100",
          TRUE ~ "other"
      )
    )

    filtered_input6 <- input_file6 %>% 
    mutate(
      category = case_when(
        type == 'RNA-LTR-normal' ~ "RNA",
        type == 'RNA-LTR-Gypsy' ~ "RNA",
        type == 'RNA-LTR-Copia' ~ "RNA",
        type == 'RNA-LINE-RTE' ~ "RNA",
        type == 'RNA-LINE-L1' ~ "RNA",
        type == 'DNA-TIR-Mutator' ~ "DNA",
        type == 'DNA-TIR-PIFHarbinger' ~ "DNA",
        type == 'DNA-TIR-CACTA' ~ "DNA",
        type == 'DNA-TIR-Tc1Mariner' ~ "DNA",
        type == 'DNA-TIR-hAT' ~ "DNA",
        type == 'helitron' ~ "other",
        type == 'pararetrovirus' ~ "other",
        type == 'Simple_repeat' ~ "other",
        type == 'border' ~ "border",
        type == 'gene' ~ "gene",
        TRUE ~ "other"
      )
    ) %>%filter(category!="other")
    filtered_input6_2 <- input_file6 %>% 
    mutate(
      category = case_when(
        type == 'gene' ~ "gene",
        TRUE ~ "other"
      )
    ) %>%filter(category!="other")
  # Create plot object
      
  # Create plot object
  p = ggplot() 
  
##All satellite OTHER satellite
  p =p+  geom_rect(data = filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end, ymin = -15,ymax =0 ,fill = category)) 
  ##strand B1
  p =p+  geom_rect(data = filtered_input3,
                   aes(xmin = block_start,xmax = block_end, ymin = -20,ymax = -18,fill = category)) 
  ##strand other
  p =p+  geom_rect(data = filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end, ymin = -20,ymax = -18,fill = category2)) 
                   
  ##monomer
  p =p+  geom_rect(data = filtered_input2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = -33,ymax = -25,fill = category))                      
  ##hor                 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = -63,ymax = -58),fill = '#6A356C',) 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = -55,ymax = -50),fill = '#6A356C',) 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = -47,ymax = -42),fill = '#6A356C',) 
  p =p+  geom_rect(data = B3,
                   aes(xmin = pos_start,xmax = pos_end, ymin = -63,ymax = -58,fill = category)) 
  p =p+  geom_rect(data = B2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = -55,ymax = -50,fill = category)) 
  p =p+  geom_rect(data = B1,
                   aes(xmin = pos_start,xmax = pos_end, ymin = -47,ymax = -42,fill = category)) 



#chip
  p =p+ geom_rect(data = filtered_input5,
            aes(xmin = start, xmax = end, ymin = 90, ymax = 90 + foldchange/foldchange_max*10, fill = category)) 
            
  p =p+  geom_rect(data=filtered_input6,
                   aes(xmin = start,xmax = end,ymin = 5,ymax = 10,fill=category)) 
  p =p+  geom_rect(data=filtered_input6_2,
                   aes(xmin = start,xmax = end,ymin = 7,ymax = 8,fill=category))     
                   
  p <- p + geom_rect(data = input_file0,
                   aes(xmin = start, xmax = end, ymin = 2, ymax = 3),
                   fill = 'black', color = 'black')
                   
  p <- p + geom_rect(data = input_file7,
                   aes(xmin = start, xmax = end, ymin = 70, ymax = 70+num/input_file7_maxnum*10),
                   fill = 'grey', color = 'black')
                   
  p <- p + geom_rect(data = input_file8,
                   aes(xmin = start, xmax = end, ymin = 50, ymax = 50+length/input_file8_maxlen*10),
                   fill = '#cceeff', color = 'black')                   
  #p=p+  xlim(min_all, max_all) + 
   p=p+   ylim(-150, 500) + 
    labs(
      x = paste0(min_all, '-', max_all),
      y = "",
      fill = ""
    ) +
    facet_grid(chr_note ~ ., scales = "free", space = "free")+
    theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  
  #width_revise=max(10,(max_all-min_all)/1000000*10)
  # Save as PDF
  pdf(file = paste0('Plot2', ".pdf"), width = 200 / 2.54, height = 500 / 2.54)
  print(p)
  dev.off()
  
}            
            """
            with open(f'Plot2.R','w') as f:
                f.write(Plot_txt)
            subprocess.run([f'Rscript Plot2.R  '], shell=True)
            
            os.chdir('../../../../')                   

        
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, sample_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(sample_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()  
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.14" or argv1=="step14.14_sample_simple":       
        print('Visualizing Vsat1')
        if  os.path.exists('./samples_satellite/14_chr_image/14_chr_complete')==False: 
            subprocess.run(["mkdir ./samples_satellite/14_chr_image/14_chr_complete"], shell=True)  
        
        sample_list=['PN40024']
        
        def run_step(one_sample):
                                
            Plot_txt=r"""
library(ggplot2)
library(dplyr)
# Get all command line arguments

### Monomer
print("")
{
  input_file0=read.table('0_chromosome', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')

  input_file3=read.table('0_block', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file4=read.table('0_otherblock', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  
  
  min3 <- min(min(input_file3$block_start, na.rm = TRUE),min(input_file3$block_end, na.rm = TRUE))
  max3 <- max(max(input_file3$block_start, na.rm = TRUE),max(input_file3$block_end, na.rm = TRUE))
  min4 <- min(min(input_file4$chr_start, na.rm = TRUE),min(input_file4$chr_end, na.rm = TRUE))
  max4 <- max(max(input_file4$chr_start, na.rm = TRUE),max(input_file4$chr_end, na.rm = TRUE))
  min_all=min(min3,min4)
  max_all=max(max3,max4)
  print(min_all)
  print(max_all)
  
  # Filter data

  
  # Define color values
  color_values <- c(
    "2" = "#b5e48c",
    "3" = "#ff8566",
    "4" = "#ECD678",
    "5" = "#ade8f4",
    "6" = "#cc99ff",
    "7" = "#7f9fcc",  
    "8" = "#ccccff",
    "9" = "#ffccff",
    "10" = "#BAD694",   
    "11" = "#bcdcd6",  
    "12" = "#f49cbb",   
    "13" = "#DDDD7A",   
    "14" = "#bbdefb", 
    "15" = "#339933",  
    "16" = "#ffb4a2",   
    "17" = "#3385ff",  
    "18" = "#00b359",  
    "19" = "#45ABC3", 
    "20" = "#e5e600",
    "21-30"= "#ff80df",
    "30+"= "#cc3300",
    "Other"="white",
    "plus"="red",
    "minus"="blue",
    "cen_66"="#00ccff",
    "cen_103"="#ff579c",
    "cen_107"="#199200",
    "cen_191"="#9933ff",
    "cen_355"="#995c17",
    "cen_383"="#ffe600",
    #"cen_677"="black",
    #"other_satellite"="black",
    
    "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699",  
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3",   
    "other_monomer" = "black",
    "Gypsy" = "#47476b",
    "Copia" = "#b30059",
    "border" = "#ECD678",
    "gene" = "#669999",
    "other" = "#ade8f4",
        "value1" = "#666699",
    "value10" = "#0066cc",
    "value100" = "#cc6699",
            "cen_RO21"="#8585ad",    #Purple-gray
    "cen_RO60"="#2784ff",    #Aqua
    "cen_RO72"="#ff6236"    #Orange-yellow
  )
  

  filtered_input3 <- input_file3 %>% 
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        strand == '+' ~ "plus",
        strand == '-' ~ "minus",
        TRUE ~ "Other"
      )
    )

  filtered_input4 <- input_file4 %>% 
    filter(centype %in% c("cen_107","cen_66","cen_103","cen_RO21","cen_RO60","cen_RO72")) %>%
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        centype == 'cen_66' ~ "cen_66",
        centype == 'cen_103' ~ "cen_103",
        centype == 'cen_107' ~ "cen_107",
        centype == 'cen_191' ~ "cen_191",
        centype == 'cen_355' ~ "cen_355",
        centype == 'cen_383' ~ "cen_383",
        centype == 'cen_RO21' ~ "cen_RO21",
        centype == 'cen_RO60' ~ "cen_RO60",
        centype == 'cen_RO72' ~ "cen_RO72",
        #TRUE ~ "other_satellite"
      )
    )  
  filtered_input4 <- filtered_input4 %>% 
    filter( centype %in% c("cen_107","cen_66","cen_103","cen_RO21","cen_RO60","cen_RO72") ) %>%
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category2 = case_when(
        strand == '+' ~ "plus",
        strand == '-' ~ "minus",
        TRUE ~ "Other"
      )
    )  
    
    
      
  # Create plot object
  p = ggplot() 

    p <- p + geom_rect(data = input_file0,
                   aes(xmin = start, xmax = end, ymin =0, ymax = 10),
                   color = 'black',linewidth=0.1,fill='#efefef')
  ##OTHER satellite
  p =p+  geom_rect(data = filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end, ymin = 0,ymax = 10,fill = category)) 
  
  ##strand B1
  #p =p+  geom_rect(data = filtered_input3,
  #                 aes(xmin = block_start,xmax = block_end, ymin = 45,ymax = 50,fill = category)) 
  ##strand other
  #p =p+  geom_rect(data = filtered_input4,
  #                 aes(xmin = chr_start,xmax = chr_end, ymin = 45,ymax = 50,fill = category2)) 
  ##name
  #p =p+  annotate("text", x = (min_all+max_all)/2, y = 35, label = "strand---satellite----B1--B2--B3 ",size=4)

  ##
  ###################################

                   
         
  #p=p+  xlim(min_all, max_all) + 
   p=p+   ylim(-1, 20) + 
    labs(
      x = paste0(min_all, '-', max_all),
      y = "",
      fill = ""
    ) +
    facet_grid(chr_note ~ ., scales = "free", space = "free")+
    theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  
  #width_revise=max(10,(max_all-min_all)/1000000*10)
  # Save as PDF
  pdf(file = paste0('Plot3_simple', ".pdf"), width = 10 / 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
  
}            
            """
            with open(f'./samples_satellite/14_chr_image/14_chr_complete/{one_sample}/Plot3.R','w') as f:
                f.write(Plot_txt)
            os.chdir(f'./samples_satellite/14_chr_image/14_chr_complete/{one_sample}')
            subprocess.run([f'Rscript Plot3.R  '], shell=True)  #>null 2>&1 
        
            
            os.chdir('../../../../')                   

        
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, sample_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(sample_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()
                
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.15" or argv1=="step14.15_ChrXX":
        ##
        if  os.path.exists('./samples_satellite/14_chr_image/15_ChrXX')==False: 
            subprocess.run(["mkdir ./samples_satellite/14_chr_image/15_ChrXX"], shell=True)     
        used_samples= [
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
        chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
        
        print('Plot ChrXX for all samples')
        for one_chr in chr_list:
            iii=0
            subprocess.run([f"mkdir ./samples_satellite/14_chr_image/15_ChrXX/{one_chr}"], shell=True)     
            with open(f"./samples_satellite/14_chr_image/15_ChrXX/{one_chr}/0_hor",'w') as f2:
                f2.write("mer\tsegment\tsegment_type\tsegment_type_num\tsegment_length\tHOR_type\tcirc_start\tcirc_end\tpos_start\tpos_end\tmarkserial_start\tmarkserial_end\tlayer\tfather_layer\tpos_length\tmarkserial_num\tHOR_repeat_num\tblock_name\tblock_start\tblock_end\tsample\tplot_y\n")
                for one_sample in used_samples:
                    iii-=1
                    if os.path.exists(f'./samples_satellite/14_chr_image/1_hor_chr_stat/{one_sample}___{one_chr}')==False:continue
                    with open(f"./samples_satellite/14_chr_image/1_hor_chr_stat/{one_sample}___{one_chr}",'r') as f:
                        next(f)
                        for line in f:
                            eachline=line.strip()
                            f2.write(f"{eachline}\t{one_sample}\t{iii}\n")
            iii=0                
            with open(f"./samples_satellite/14_chr_image/15_ChrXX/{one_chr}/0_monomer",'w') as f2:
                f2.write("serial\tcirc_start\tcirc_end\tpos_start\tpos_end\tlength\tmonomer_type\tblock_name\tblock_start\tblock_end\tsample\tplot_y\n")
                for one_sample in used_samples:
                    iii-=1
                    if os.path.exists(f'./samples_satellite/14_chr_image/2_monomer_chr_stat/{one_sample}___{one_chr}')==False:continue
                    with open(f"./samples_satellite/14_chr_image/2_monomer_chr_stat/{one_sample}___{one_chr}",'r') as f:
                        next(f)
                        for line in f:
                            eachline=line.strip()
                            f2.write(f"{eachline}\t{one_sample}\t{iii}\n")                            
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.15" or argv1=="step14.15p_ChrXX_unit2":        
        chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
        for one_chr in chr_list:
            Plot_txt=r"""
library(ggplot2)
library(dplyr)
# Get all command line parameters

### Monomer
print("")
{
  input_file3=read.table('0_hor', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    min3 <- min(min(input_file3$block_start, na.rm = TRUE),min(input_file3$block_end, na.rm = TRUE))  
    max3 <- max(max(input_file3$block_start, na.rm = TRUE),max(input_file3$block_end, na.rm = TRUE))
  sample_num=-min(input_file3$plot_y, na.rm = TRUE)
  #print(sample_num)
  

  input_file2=read.table('0_monomer', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
filtered_input2 <- input_file2 %>% 
    mutate(
      category = case_when(
        monomer_type == 'CEN107' ~ "CEN107",
        monomer_type == 'CEN107-like' ~ "CEN107-like",
        monomer_type == 'CEN135' ~ "CEN135",
        monomer_type == 'CEN135-like' ~ "CEN135-like",
        monomer_type == 'CEN79' ~ "CEN79",
        monomer_type == 'CEN79-like' ~ "CEN79-like",
        monomer_type == 'CEN51' ~ "CEN51",
        monomer_type == 'CEN51-like' ~ "CEN51-like",        
        monomer_type == 'other' ~ "other_monomer",
        TRUE ~ "other_monomer"
      )
    )
  
  
  
  filtered_input_file3 <- input_file3 %>% 
    filter(pos_length > 300 & HOR_repeat_num>=2) %>%     
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        mer == 2 ~ "2",
        mer == 3 ~ "3",
        mer == 4 ~ "4",
        mer == 5 ~ "5",
        mer == 6 ~ "6",
        mer == 7 ~ "7",
        mer == 8 ~ "8",
        mer == 9 ~ "9",
        mer == 10 ~ "10",
        mer == 11 ~ "11",
        mer == 12 ~ "12",
        mer == 13 ~ "13",
        mer == 14 ~ "14",
        mer == 15 ~ "15",
        mer == 16 ~ "16",
        mer == 17 ~ "17",
        mer == 18 ~ "18",
        mer == 19 ~ "19",
        mer == 20 ~ "20",
        mer>20 & mer<=30 ~ "21-30",
        mer>30 ~ "30+",
        TRUE ~ "Other"
      )
    )
    
    
  B1 <- filtered_input_file3 %>% filter(father_layer == 1)
  B2 <- filtered_input_file3 %>% filter(father_layer == 2)
  B3 <- filtered_input_file3 %>% filter(father_layer == 3)
 
 # Extract sub-table with deduplication
  block_start_block_end_plot_y <- filtered_input_file3 %>%
  select(block_start,block_end,plot_y) %>%       
  distinct()  
  
  sample_plot_y <- filtered_input_file3 %>%
  select(sample,plot_y) %>%      
  distinct()  
  
  
  
  # Define color values
  color_values <- c(
    "2" = "#b5e48c",
    "3" = "#ff8566",
    "4" = "#ECD678",
    "5" = "#ade8f4",
    "6" = "#cc99ff",
    "7" = "#7f9fcc",  
    "8" = "#ccccff",
    "9" = "#ffccff",
    "10" = "#BAD694",   
    "11" = "#bcdcd6",  
    "12" = "#f49cbb",   
    "13" = "#DDDD7A",   
    "14" = "#bbdefb", 
    "15" = "#339933",  
    "16" = "#ffb4a2",   
    "17" = "#3385ff",  
    "18" = "#00b359",  
    "19" = "#45ABC3", 
    "20" = "#e5e600",
    "21-30"= "#ff80df",
    "30+"= "#cc3300",
    "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699",  
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3",      
    "other_monomer" = "black"
     #   TRUE ~ "Other"
  )

    
      
  # Create plot object
  p = ggplot() 


p =p+  geom_rect(data = filtered_input2,
                   aes(xmin = pos_start,xmax = pos_end, ymin =  plot_y+0.1,ymax =  plot_y+0.2,fill = category))   

  p =p+  geom_rect(data=block_start_block_end_plot_y,
                   aes(xmin = block_start,xmax = block_end,ymin = plot_y-0.2,ymax = plot_y),fill = '#6A356C',) 
  p =p+  geom_rect(data=block_start_block_end_plot_y,
                   aes(xmin = block_start,xmax = block_end,ymin = plot_y-0.45,ymax = plot_y-0.25),fill = '#6A356C',) 
  p =p+  geom_rect(data=block_start_block_end_plot_y,
                   aes(xmin = block_start,xmax = block_end,ymin = plot_y-0.7,ymax = plot_y-0.5),fill = '#6A356C',) 
  p =p+  geom_rect(data = B1,
                   aes(xmin = pos_start,xmax = pos_end, ymin = plot_y-0.2,ymax = plot_y,fill = category)) 
  p =p+  geom_rect(data = B2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = plot_y-0.45,ymax = plot_y-0.25,fill = category)) 
  p =p+  geom_rect(data = B3,
                   aes(xmin = pos_start,xmax = pos_end, ymin = plot_y-0.7,ymax = plot_y-0.5,fill = category)) 
  
  ## Names
p <- p +  geom_text(data = sample_plot_y, aes(x = min3, y = plot_y-0.25, label = sample), size =10)


  ##
  ###################################

                   
    
   p=p+theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  
 height_revise=max(10,sample_num+1)
  # Save as PDF
  pdf(file = paste0('hor_unit2', ".pdf"), width = 100 / 2.54, height = height_revise / 2.54)
  print(p)
  dev.off()
  
}            
        """
            with open(f'./samples_satellite//14_chr_image/15_ChrXX/{one_chr}/hor.R','w') as f:
                f.write(Plot_txt)
            os.chdir(f'./samples_satellite//14_chr_image/15_ChrXX/{one_chr}')
            subprocess.run([f'Rscript hor.R  '], shell=True)  #>null 2>&1 
            subprocess.run([f'mv hor_unit2.pdf {one_chr}_unit2.pdf  '], shell=True) 
            
            os.chdir('../../../../')                   
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.15" or argv1=="step14.15p_ChrXX_unit3":        
        chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
        for one_chr in chr_list:
            Plot_txt=r"""
library(ggplot2)
library(dplyr)
# Get all command line parameters

### Monomer
print("")
{
  input_file3=read.table('0_hor', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    min3 <- min(min(input_file3$block_start, na.rm = TRUE),min(input_file3$block_end, na.rm = TRUE))  
    max3 <- max(max(input_file3$block_start, na.rm = TRUE),max(input_file3$block_end, na.rm = TRUE))
  sample_num=-min(input_file3$plot_y, na.rm = TRUE)
  #print(sample_num)
  

  input_file2=read.table('0_monomer', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
filtered_input2 <- input_file2 %>% 
    mutate(
      category = case_when(
        monomer_type == 'CEN107' ~ "CEN107",
        monomer_type == 'CEN107-like' ~ "CEN107-like",
        monomer_type == 'CEN135' ~ "CEN135",
        monomer_type == 'CEN135-like' ~ "CEN135-like",
        monomer_type == 'CEN79' ~ "CEN79",
        monomer_type == 'CEN79-like' ~ "CEN79-like",
        monomer_type == 'CEN51' ~ "CEN51",
        monomer_type == 'CEN51-like' ~ "CEN51-like",        
        monomer_type == 'other' ~ "other_monomer",
        TRUE ~ "other_monomer"
      )
    )
  
  
  
  filtered_input_file3 <- input_file3 %>% 
    filter(pos_length > 300 & HOR_repeat_num>=3) %>%     
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        mer == 2 ~ "2",
        mer == 3 ~ "3",
        mer == 4 ~ "4",
        mer == 5 ~ "5",
        mer == 6 ~ "6",
        mer == 7 ~ "7",
        mer == 8 ~ "8",
        mer == 9 ~ "9",
        mer == 10 ~ "10",
        mer == 11 ~ "11",
        mer == 12 ~ "12",
        mer == 13 ~ "13",
        mer == 14 ~ "14",
        mer == 15 ~ "15",
        mer == 16 ~ "16",
        mer == 17 ~ "17",
        mer == 18 ~ "18",
        mer == 19 ~ "19",
        mer == 20 ~ "20",
        mer>20 & mer<=30 ~ "21-30",
        mer>30 ~ "30+",
    
      )
    )
    
    
  B1 <- filtered_input_file3 %>% filter(father_layer == 1)
  B2 <- filtered_input_file3 %>% filter(father_layer == 2)
  B3 <- filtered_input_file3 %>% filter(father_layer == 3)
 
 # Extract sub-table with deduplication
  block_start_block_end_plot_y <- input_file2 %>%
  select(block_start,block_end,plot_y) %>%       
  distinct()  
  
  sample_plot_y <- input_file2 %>%
  select(sample,plot_y) %>%      
  distinct()  
  
  
  
  # Define color values
  color_values <- c(
    "2" = "#b5e48c",
    "3" = "#ff8566",
    "4" = "#ECD678",
    "5" = "#ade8f4",
    "6" = "#cc99ff",
    "7" = "#7f9fcc",  
    "8" = "#ccccff",
    "9" = "#ffccff",
    "10" = "#BAD694",   
    "11" = "#bcdcd6",  
    "12" = "#f49cbb",   
    "13" = "#DDDD7A",   
    "14" = "#bbdefb", 
    "15" = "#339933",  
    "16" = "#ffb4a2",   
    "17" = "#3385ff",  
    "18" = "#00b359",  
    "19" = "#45ABC3", 
    "20" = "#e5e600",
    "21-30"= "#ff80df",
    "30+"= "#cc3300",
    "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699", 
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3",      
    "other_monomer" = "black"
     #   TRUE ~ "Other"
  )

    
      
  # Create plot object
  p = ggplot() 


p =p+  geom_rect(data = filtered_input2,
                   aes(xmin = pos_start,xmax = pos_end, ymin =  plot_y+0.1,ymax =  plot_y+0.2,fill = category))   


  p =p+  geom_rect(data=block_start_block_end_plot_y,
                   aes(xmin = block_start,xmax = block_end,ymin = plot_y-0.2,ymax = plot_y),fill = '#6A356C',) 
  p =p+  geom_rect(data=block_start_block_end_plot_y,
                   aes(xmin = block_start,xmax = block_end,ymin = plot_y-0.45,ymax = plot_y-0.25),fill = '#6A356C',) 
  p =p+  geom_rect(data=block_start_block_end_plot_y,
                   aes(xmin = block_start,xmax = block_end,ymin = plot_y-0.7,ymax = plot_y-0.5),fill = '#6A356C',) 
  p =p+  geom_rect(data = B1,
                   aes(xmin = pos_start,xmax = pos_end, ymin = plot_y-0.2,ymax = plot_y,fill = category)) 
  p =p+  geom_rect(data = B2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = plot_y-0.45,ymax = plot_y-0.25,fill = category)) 
  p =p+  geom_rect(data = B3,
                   aes(xmin = pos_start,xmax = pos_end, ymin = plot_y-0.7,ymax = plot_y-0.5,fill = category)) 
  
  ## Names
p <- p +  geom_text(data = sample_plot_y, aes(x = min3, y = plot_y-0.25, label = sample), size =10)


  ##
  ###################################

                   
    
   p=p+theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  
 height_revise=max(10,sample_num)
  # Save as PDF
  pdf(file = paste0('hor_unit3', ".pdf"), width = 100 / 2.54, height = height_revise / 2.54)
  print(p)
  dev.off()
  
}            
        """
            with open(f'./samples_satellite//14_chr_image/15_ChrXX/{one_chr}/hor.R','w') as f:
                f.write(Plot_txt)
            os.chdir(f'./samples_satellite//14_chr_image/15_ChrXX/{one_chr}')
            subprocess.run([f'Rscript hor.R  '], shell=True)  #>null 2>&1 
            subprocess.run([f'mv hor_unit3.pdf {one_chr}_unit3.pdf  '], shell=True) 
            
            os.chdir('../../../../')                                   
            
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.16" or argv1=="step14.16_sample":       
        print('Visualize SatHOR/monomer/strand')
        if  os.path.exists('./samples_satellite/14_chr_image/16_chr_complete')==False: 
            subprocess.run(["mkdir ./samples_satellite/14_chr_image/16_chr_complete"], shell=True)  
        
        sample_list=['PN40024']
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
        def run_step(one_sample):
            subprocess.run([f"mkdir ./samples_satellite/14_chr_image/16_chr_complete/{one_sample}"], shell=True) 
            chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
            dict_sheet_info={}
            dict_sheet_info['1_hor_chr_stat']=[]
            dict_sheet_info['2_monomer_chr_stat']=[]
            dict_sheet_info['3_block_chr_stat']=[]
            dict_sheet_info['4_otherblock_chr_stat']=[]
            for one_chr in chr_list:
                one_hapchr=one_sample+'___'+one_chr
                input_file1=f"./samples_satellite/14_chr_image/1_hor_chr_stat/{one_hapchr}"
                input_file2=f"./samples_satellite/14_chr_image/2_monomer_chr_stat/{one_hapchr}"
                input_file3=f"./samples_satellite/14_chr_image/3_block_chr_stat/{one_hapchr}"
                input_file4=f"./samples_satellite/14_chr_image/4_otherblock_chr_stat/{one_hapchr}"
                if  os.path.exists(input_file1)==True: 
                    k=0
                    with open(input_file1,'r') as f:
                        for line in f:
                            eachline=line.strip()
                            k+=1
                            if k==1:hor_chr_stat_head=eachline
                            else:
                                dict_sheet_info['1_hor_chr_stat'].append(f"{eachline}\t{one_chr}")
                if  os.path.exists(input_file2)==True:                 
                    k=0
                    with open(input_file2,'r') as f:
                        for line in f:
                            eachline=line.strip()
                            k+=1
                            if k==1:monomer_chr_stat_head=eachline
                            else:
                                dict_sheet_info['2_monomer_chr_stat'].append(f"{eachline}\t{one_chr}")
                if  os.path.exists(input_file3)==True:                 
                    k=0
                    with open(input_file3,'r') as f:
                        for line in f:
                            eachline=line.strip()
                            k+=1
                            if k==1:block_chr_stat_head=eachline
                            else:
                                dict_sheet_info['3_block_chr_stat'].append(f"{eachline}\t{one_chr}")
                if  os.path.exists(input_file4)==True:                 
                    k=0
                    with open(input_file4,'r') as f:
                        for line in f:
                            eachline=line.strip()
                            k+=1
                            if k==1:otherblock_chr_head=eachline
                            else:
                                dict_sheet_info['4_otherblock_chr_stat'].append(f"{eachline}\t{one_chr}")
            #########################################################################                
            with open(f'./samples_satellite/14_chr_image/16_chr_complete/{one_sample}/0_hor','w') as f:
                f.write(f"{hor_chr_stat_head}\tchr_note\n")
                for one in dict_sheet_info['1_hor_chr_stat']:
                    f.write(f"{one}\n")
            with open(f'./samples_satellite/14_chr_image/16_chr_complete/{one_sample}/0_monomer','w') as f:
                f.write(f"{monomer_chr_stat_head}\tchr_note\n")
                for one in dict_sheet_info['2_monomer_chr_stat']:
                    f.write(f"{one}\n")
            with open(f'./samples_satellite/14_chr_image/16_chr_complete/{one_sample}/0_block','w') as f:
                f.write(f"{block_chr_stat_head}\tchr_note\n")
                for one in dict_sheet_info['3_block_chr_stat']:
                    f.write(f"{one}\n")
            with open(f'./samples_satellite/14_chr_image/16_chr_complete/{one_sample}/0_otherblock','w') as f:
                f.write(f"{otherblock_chr_head}\tchr_note\n")
                for one in dict_sheet_info['4_otherblock_chr_stat']:
                    f.write(f"{one}\n")
                    
            with open(f'./samples_satellite/14_chr_image/16_chr_complete/{one_sample}/0_chromosome','w') as f2:
                f2.write(f"chr_note\tstart\tend\n")
                with open(f'./new_work_dir/chr2EDTA/0_prepare/{one_sample}/{one_sample}.fasta.fai','r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        f2.write(f"{eachline_arr[0]}\t1\t{eachline_arr[1]}\n")  

            os.chdir(f'./samples_satellite/14_chr_image/16_chr_complete/{one_sample}')
            #subprocess.run([f'RScript Plot1.R  '], shell=True)  #>null 2>&1 
            
            
            Plot_txt=r"""
library(ggplot2)
library(dplyr)
# Get all command line parameters

### Monomer
print("")
{
  input_file0=read.table('0_chromosome', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file1=read.table('0_hor', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file2=read.table('0_monomer', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file3=read.table('0_block', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file4=read.table('0_otherblock', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  

  
  min3 <- min(min(input_file3$block_start, na.rm = TRUE),min(input_file3$block_end, na.rm = TRUE))
  max3 <- max(max(input_file3$block_start, na.rm = TRUE),max(input_file3$block_end, na.rm = TRUE))
  min4 <- min(min(input_file4$chr_start, na.rm = TRUE),min(input_file4$chr_end, na.rm = TRUE))
  max4 <- max(max(input_file4$chr_start, na.rm = TRUE),max(input_file4$chr_end, na.rm = TRUE))
  min_all=min(min3,min4)
  max_all=max(max3,max4)
  print(min_all)
  print(max_all)
  
  # Filter data
  filtered_input1 <- input_file1 %>% 
    filter(pos_length > 300 & HOR_repeat_num>=3) %>%     ############################################## Difference from step9 is the filter in HOR identification
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        mer == 2 ~ "2",
        mer == 3 ~ "3",
        mer == 4 ~ "4",
        mer == 5 ~ "5",
        mer == 6 ~ "6",
        mer == 7 ~ "7",
        mer == 8 ~ "8",
        mer == 9 ~ "9",
        mer == 10 ~ "10",
        mer == 11 ~ "11",
        mer == 12 ~ "12",
        mer == 13 ~ "13",
        mer == 14 ~ "14",
        mer == 15 ~ "15",
        mer == 16 ~ "16",
        mer == 17 ~ "17",
        mer == 18 ~ "18",
        mer == 19 ~ "19",
        mer == 20 ~ "20",
        mer>20 & mer<=30 ~ "21-30",
        mer>30 ~ "30+",
        TRUE ~ "Other"
      )
    )
  
  
  B1 <- filtered_input1 %>% filter(father_layer == 1)
  B2 <- filtered_input1 %>% filter(father_layer == 2)
  B3 <- filtered_input1 %>% filter(father_layer == 3)
  
  
  
  # Define color values
  color_values <- c(
    "2" = "#b5e48c",
    "3" = "#ff8566",
    "4" = "#ECD678",
    "5" = "#ade8f4",
    "6" = "#cc99ff",
    "7" = "#7f9fcc",  
    "8" = "#ccccff",
    "9" = "#ffccff",
    "10" = "#BAD694",   
    "11" = "#bcdcd6",  
    "12" = "#f49cbb",   
    "13" = "#DDDD7A",   
    "14" = "#bbdefb", 
    "15" = "#339933",  
    "16" = "#ffb4a2",   
    "17" = "#3385ff",  
    "18" = "#00b359",  
    "19" = "#45ABC3", 
    "20" = "#e5e600",
    "21-30"= "#ff80df",
    "30+"= "#cc3300",
    "Other"="white",
    "plus"="red",
    "minus"="blue",
    "cen_66"="#00ccff",
    "cen_103"="#ff579c",
    "cen_107"="#199200",
    "cen_191"="#9933ff",
    "cen_355"="#995c17",
    "cen_383"="#ffe600",
    #"cen_677"="black",
    #"other_satellite"="black",
    
    "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699",  
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3",     
    "other_monomer" = "black",
    "Gypsy" = "#47476b",
    "Copia" = "#b30059",
    "border" = "#ECD678",
    "gene" = "#669999",
    "other" = "#ade8f4",
        "value1" = "#666699",
    "value10" = "#0066cc",
    "value100" = "#cc6699"
  )
  


 
filtered_input2 <- input_file2 %>% 
    mutate(
      category = case_when(
        monomer_type == 'CEN107' ~ "CEN107",
        monomer_type == 'CEN107-like' ~ "CEN107-like",
        monomer_type == 'CEN135' ~ "CEN135",
        monomer_type == 'CEN135-like' ~ "CEN135-like",
        monomer_type == 'CEN79' ~ "CEN79",
        monomer_type == 'CEN79-like' ~ "CEN79-like",
        monomer_type == 'CEN51' ~ "CEN51",
        monomer_type == 'CEN51-like' ~ "CEN51-like",        
        monomer_type == 'other' ~ "other_monomer",
        TRUE ~ "other_monomer"
      )
    )
    
  filtered_input3 <- input_file3 %>% 
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        strand == '+' ~ "plus",
        strand == '-' ~ "minus",
        TRUE ~ "Other"
      )
    )

  filtered_input4 <- input_file4 %>% 
    filter(centype %in% c("cen_66","cen_103","cen_107","cen_191","cen_355","cen_383")) %>%
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        centype == 'cen_66' ~ "cen_66",
        centype == 'cen_103' ~ "cen_103",
        centype == 'cen_107' ~ "cen_107",
        centype == 'cen_191' ~ "cen_191",
        centype == 'cen_355' ~ "cen_355",
        centype == 'cen_383' ~ "cen_383",
        #TRUE ~ "other_satellite"
      )
    )  
  filtered_input4 <- filtered_input4 %>% 
    filter( centype %in% c("cen_66","cen_103","cen_107","cen_191","cen_355","cen_383") ) %>%
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category2 = case_when(
        strand == '+' ~ "plus",
        strand == '-' ~ "minus",
        TRUE ~ "Other"
      )
    )  
    
    

      
  # Create plot object
  p = ggplot() 
  
## All satellite OTHER satellite
  p =p+  geom_rect(data = filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end, ymin = -15,ymax =0 ,fill = category)) 
  ##strand B1
  p =p+  geom_rect(data = filtered_input3,
                   aes(xmin = block_start,xmax = block_end, ymin = -20,ymax = -18,fill = category)) 
  ##strand other
  p =p+  geom_rect(data = filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end, ymin = -20,ymax = -18,fill = category2)) 
                   
  ##monomer
  p =p+  geom_rect(data = filtered_input2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = -33,ymax = -25,fill = category))                      
  ##hor                 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = -63,ymax = -58),fill = '#6A356C',) 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = -55,ymax = -50),fill = '#6A356C',) 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = -47,ymax = -42),fill = '#6A356C',) 
  p =p+  geom_rect(data = B3,
                   aes(xmin = pos_start,xmax = pos_end, ymin = -63,ymax = -58,fill = category)) 
  p =p+  geom_rect(data = B2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = -55,ymax = -50,fill = category)) 
  p =p+  geom_rect(data = B1,
                   aes(xmin = pos_start,xmax = pos_end, ymin = -47,ymax = -42,fill = category)) 

  input_file0=read.table('0_chromosome', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  p <- p + geom_rect(data = input_file0,
                   aes(xmin = start, xmax = end, ymin = 2, ymax = 3),
                   fill = 'black', color = 'black')
   p=p+   ylim(-65, 5) + 
    labs(
      x = paste0(min_all, '-', max_all),
      y = "",
      fill = ""
    ) +
    facet_grid(chr_note ~ ., scales = "free", space = "free")+
    theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  
  # Save as PDF
  pdf(file = paste0('Plot2', ".pdf"), width = 50 / 2.54, height = 50 / 2.54)
  print(p)
  dev.off()
  
}            
        """
            with open(f'Plot2.R','w') as f:
                f.write(Plot_txt)
            subprocess.run([f'RScript Plot2.R  '], shell=True)
            
            os.chdir('../../../../')                   

        
        with Pool(processes=10) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, sample_list), start=1):
                # Here you can process the results, such as store or print
                progress = (i / len(sample_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()       
     
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.17" or argv1=="step14.17_HOR_type":
        print("step14.17_HOR_type,17s")
        if  os.path.exists('./samples_satellite/14_chr_image/17_HOR_type')==False: 
            subprocess.run(["mkdir ./samples_satellite/14_chr_image/17_HOR_type"], shell=True)    
        chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
        
        ##28-28-28-23|28-28-23|28-28-28-28-23 HOR monomer arrangement, but it is equivalent to 28-28-28-28-23|28-28-28-23|28-28-23, i.e., rotation, translocation. How to obtain a unique index
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
        dict_sample_chr_type_length={}    
        for one_chr in chr_list:
            input_file=f"./samples_satellite/14_chr_image/15_ChrXX/{one_chr}/0_hor"
            with open(input_file,'r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split()
                    mer=	                eachline_arr[0]
                    segment=	            eachline_arr[1]
                    segment=normalize_hor_pattern_rotation_only(segment)
                    segment_type=           eachline_arr[2]
                    segment_rough_length=   0
                    if "." in segment:continue
                    for one_subunit_length in segment.replace('-','|').replace('||','|').split('|'):
                        if len(one_subunit_length)==0:continue
                        segment_rough_length+=int(one_subunit_length)
                    HOR_class=f"{mer}-mer|{segment_rough_length}|{segment_type}|||{segment}"
                    length=         int(eachline_arr[14])
                    chromosome=     eachline_arr[17].split(':')[1][:-1].replace('region_','Chr')
                    sample=         eachline_arr[20]
                    #
                    if sample not in dict_sample_chr_type_length:
                        dict_sample_chr_type_length[sample]={}
                    if chromosome not in dict_sample_chr_type_length[sample]:
                        dict_sample_chr_type_length[sample][chromosome]={}
                    if HOR_class not in dict_sample_chr_type_length[sample][chromosome]:
                        dict_sample_chr_type_length[sample][chromosome][HOR_class]=0
                    dict_sample_chr_type_length[sample][chromosome][HOR_class]+=length
                    
        dict_sample_type_length={} 
        dict_type_info={}
        with open('./samples_satellite/14_chr_image/17_HOR_type/0_sample_chr_type_length','w') as f:
            f.write(f"sample\tchromosome\tHOR_class\tlength\n")
            for one_sample,dict_chr_type_length in dict_sample_chr_type_length.items():
                for one_chr, dict_type_length in dict_chr_type_length.items():
                    for one_type, length in dict_type_length.items():
                        f.write(f"{one_sample}\t{one_chr}\t{one_type}\t{length}\n")
                        ###
                        if sample not in dict_sample_type_length:
                            dict_sample_type_length[sample]={}
                        if one_type not in dict_sample_type_length[sample]:
                            dict_sample_type_length[sample][one_type]=0
                        dict_sample_type_length[sample][one_type]+=   int(length) 
                        ##
                        if one_type not in dict_type_info:
                            dict_type_info[one_type]={}
                            dict_type_info[one_type]['length']=0
                            dict_type_info[one_type]['list']=[]
                        dict_type_info[one_type]['length']+=   int(length)    
                        dict_type_info[one_type]['list'].append(int(length))
        ## Without distinguishing chromosomes
        with open('./samples_satellite/14_chr_image/17_HOR_type/0_sample_type_length','w') as f:
            f.write(f"sample\tHOR_class\tlength\n")
            for one_sample,dict_type_length in dict_sample_type_length.items():
                for one_type, length in dict_type_length.items():
                    f.write(f"{one_sample}\t{one_type}\t{length}\n")
        ## Without distinguishing chromosomes and samples
        with open('./samples_satellite/14_chr_image/17_HOR_type/0_type_length','w') as f:
            f.write(f"HOR_class\tlength\tmax_one_chr\n")
            for  one_type, info in dict_type_info.items():
                length=info['length']
                length_list=info['list']
                max_one_chr=max(length_list)
                f.write(f"{one_type}\t{length}\t{max_one_chr}\n")
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.17" or argv1=="step14.17s":   
        print('Ranking, e.g., find the types ranked 1, 2, 3, 4, 5 on a certain chromosome, length must be greater than 200k otherwise ignore')
        info_list=[]
        with open('./samples_satellite/14_chr_image/17_HOR_type/0_sample_chr_type_length','r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                sample	        =eachline_arr[0]
                chromosome	    =eachline_arr[1]
                HOR_class	    =eachline_arr[2]
                length          =int(eachline_arr[3])
                if length<200000: continue
                info_list.append([sample,chromosome,HOR_class,length])
        sorted_list = sorted(info_list, key=lambda x: x[3], reverse=True)
        # print(sorted_list)
        #
        with open('./samples_satellite/14_chr_image/17_HOR_type/0_sample_chr_type_length_order','w') as f:
            f.write(f"sample\tchromosome\tHOR_class\tlength\torder\n")
            dict_samplechr_order={}
            for info in sorted_list:
                sample,chromosome,HOR_class,length=info
                
                samplechr=f"{sample}|{chromosome}"
                if samplechr not in dict_samplechr_order:
                    dict_samplechr_order[samplechr]=0
                dict_samplechr_order[samplechr]+=1
                order=dict_samplechr_order[samplechr]            
                f.write(f"{sample}\t{chromosome}\t{HOR_class}\t{length}\t{order}\n")
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.17" or argv1=="step14.17s2":    
        print('Importance ranking')
        def format_mon_list(mon_list):
            """
            Convert MON_LIST to a compact string
            
            Example:
                ["MON79", "MON79", "MON107"]  -> "(MON79)2-MON107"
                ["MON107", "MON107"]          -> "(MON107)2"
                ["MON51", "MON107"]           -> "MON51-MON107"
                ["MON79", "MON79", "MON79"]   -> "(MON79)3"
                ["MON79", "MON79", "MON107", "MON107", "MON107"] -> "(MON79)2-(MON107)3"
            """
            if not mon_list:
                return ""
            
            result_parts = []
            current = mon_list[0]
            count = 1
            
            for mon in mon_list[1:]:
                if mon == current:
                    count += 1
                else:
                    # Output current unit
                    if count == 1:
                        result_parts.append(current)
                    else:
                        result_parts.append(f"({current}){count}")
                    current = mon
                    count = 1
            
            # Process the last unit
            if count == 1:
                result_parts.append(current)
            else:
                result_parts.append(f"({current}){count}")
            
            return "-".join(result_parts)
        dict_horclass_len={}
        top1_set=set()
        top2_set=set()
        top3_set=set()        
        dict_one_bestorder={}
        dict_HORname_desc={}
        with  open('./samples_satellite/14_chr_image/17_HOR_type/0_sample_chr_type_length_order','r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                sample,chromosome,HOR_class,length,order=eachline_arr
                if '-23-' in HOR_class:continue
                HOR_ARR=HOR_class.split('|||')[-1].strip('-').split('|')
                MON_LIST=[]
                for one in HOR_ARR:
                    one_arr=one.split('-')
                    one_len=0
                    #print(HOR_class)
                    for one_part in one_arr:
                        one_len+=int(one_part)
                    MON_LIST.append(f"MON{one_len}")
                MON_name=format_mon_list(MON_LIST) 
                dict_HORname_desc[MON_name]=HOR_class
                if MON_name not in dict_horclass_len:dict_horclass_len[MON_name]=0
                dict_horclass_len[MON_name]+=int(length)
                ###
                order=int(order)
                if MON_name not in dict_one_bestorder:dict_one_bestorder[MON_name]=order
                elif order< dict_one_bestorder[MON_name]: dict_one_bestorder[MON_name]=order
                if int(length)<200000:continue
                if order<=1: top1_set.add(MON_name)
                if order<=2: top2_set.add(MON_name)
                if order<=3: top3_set.add(MON_name)
        top1_list=list(top1_set)       ;top1_list.sort()
        top2_list=list(top2_set)       ;top2_list.sort()
        top3_list=list(top3_set)       ;top3_list.sort()

        with open('./samples_satellite/14_chr_image/17_HOR_type/0_sample_chr_type_length_order_stat_top3','w')        as f:    
            f.write(f"\n\n#####################\n")    
            for one in top1_list:
                f.write(f"Top1\t{one}\t{dict_horclass_len[one]}\n")
                
            f.write(f"\n\n#####################\n")    
            for one in top2_list:
                f.write(f"Top2\t{one}\t{dict_horclass_len[one]}\n")
     
            f.write(f"\n\n#####################\n")    
            for one in top3_list:
                f.write(f"Top3\t{one}\t{dict_horclass_len[one]}\n")    
                
        with open('./samples_satellite/14_chr_image/17_HOR_type/0_sample_chr_type_length_order_stat_alllength','w')        as f:
            f.write(f"MON_name\tHOR_desc\tlength\tbestorder\n")
            # Sort by length in descending order
            for MON_name, length in sorted(dict_horclass_len.items(), key=lambda x: x[1], reverse=True):
                if MON_name not in dict_one_bestorder: bestorder='NA'
                else:bestorder=dict_one_bestorder[MON_name]
                HOR_class=dict_HORname_desc[MON_name]
                f.write(f"{MON_name}\t{HOR_class}\t{length}\t{bestorder}\n")
        
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.17" or argv1=="step14.17p":
        if  os.path.exists('./samples_satellite/14_chr_image/17_HOR_type/1_chr_plot')==False: 
            subprocess.run(["mkdir ./samples_satellite/14_chr_image/17_HOR_type/1_chr_plot/"], shell=True)        
        print('Filter, keep top 3 per chromosome, length > 200k, manually convert names here')
        '''dict_change={
            "2-mer|158|28-23;28-28-28-23|||28-23|28-28-28-23":"MON51-MON107",
            "2-mer|158|28-23;28-28-28-23|||28-28-28-23|28-23":"MON51-MON107",
            "2-mer|186|28-28-23;28-28-28-23|||28-28-23|28-28-28-23":"MON79-MON107",
            "2-mer|186|28-28-23;28-28-28-23|||28-28-28-23|28-28-23":"MON79-MON107",
            "2-mer|214|28-28-28-23|||28-28-28-23|28-28-28-23":"(MON107)2",
            "2-mer|242|28-28-28-23;28-28-28-28-23|||28-28-28-23|28-28-28-28-23":"MON107-MON135",
            "2-mer|242|28-28-28-23;28-28-28-28-23|||28-28-28-28-23|28-28-28-23":"MON107-MON135",
            "2-mer|270|28-28-28-23;28-28-28-28-28-23|||28-28-28-23|28-28-28-28-28-23":"MON107-MON163",
            "2-mer|270|28-28-28-28-23|||28-28-28-28-23|28-28-28-28-23":"(MON135)2",
            "3-mer|265|28-23;28-28-28-23|||28-28-28-23|28-28-28-23|28-23":"MON51-(MON107)2",
            "3-mer|265|28-28-23;28-28-28-23|||28-28-23|28-28-23|28-28-28-23":"(MON79)2-MON107",
            "3-mer|265|28-28-23;28-28-28-23|||28-28-23|28-28-28-23|28-28-23":"(MON79)2-MON107",
            "3-mer|265|28-28-23;28-28-28-23|||28-28-28-23|28-28-23|28-28-23":"(MON79)2-MON107",
            "3-mer|293|28-28-23;28-28-28-23|||28-28-23|28-28-28-23|28-28-28-23":"MON79-(MON107)2",
            "3-mer|293|28-28-23;28-28-28-23|||28-28-28-23|28-28-28-23|28-28-23":"MON79-(MON107)2",
            "3-mer|293|28-28-23;28-28-28-23|||28-28-28-23|28-28-23|28-28-28-23":"MON79-(MON107)2",
            "3-mer|321|28-28-23;28-28-28-23;28-28-28-28-23|||28-28-23|28-28-28-23|28-28-28-28-23":"MON79-MON107-MON135",
            "3-mer|321|28-28-28-23|||28-28-28-23|28-28-28-23|28-28-28-23":"(MON107)3",
            "3-mer|405|28-28-28-28-23|||28-28-28-28-23|28-28-28-28-23|28-28-28-28-23":"(MON135)3",
            "4-mer|316|28-23;28-28-23;28-28-28-23|||28-28-23|28-28-28-23|28-28-23|28-23":"MON51-(MON79)2-MON107",
            "4-mer|428|28-28-28-23|||28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23":"(MON107)4",
            "5-mer|535|28-28-28-23|||28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23":"(MON107)5",
            "6-mer|642|28-28-28-23|||28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23":"(MON107)6",
            "10-mer|1070|28-28-28-23|||28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23":"(MON107)10"
        }'''
        dict_change = {
            "2-mer|158|28-23;28-28-28-23|||28-23|28-28-28-23": "MON51-MON107",
            "2-mer|186|28-28-23;28-28-28-23|||28-28-23|28-28-28-23": "MON79-MON107",
            "2-mer|214|28-28-28-23|||28-28-28-23|28-28-28-23": "(MON107)2",
            "2-mer|242|28-28-28-23;28-28-28-28-23|||28-28-28-23|28-28-28-28-23": "MON107-MON135",
            "2-mer|270|28-28-28-23;28-28-28-28-28-23|||28-28-28-23|28-28-28-28-28-23": "MON107-MON163",
            "2-mer|270|28-28-28-28-23|||28-28-28-28-23|28-28-28-28-23": "(MON135)2",
            "3-mer|265|28-23;28-28-28-23|||28-23|28-28-28-23|28-28-28-23": "MON51-(MON107)2",           # Normalized
            "3-mer|265|28-28-23;28-28-28-23|||28-28-23|28-28-23|28-28-28-23": "(MON79)2-MON107",      # Normalized
            "3-mer|293|28-28-23;28-28-28-23|||28-28-23|28-28-28-23|28-28-28-23": "MON79-(MON107)2",      # Normalized (note different from above)
            "3-mer|321|28-28-23;28-28-28-23;28-28-28-28-23|||28-28-23|28-28-28-23|28-28-28-28-23": "MON79-MON107-MON135",
            "3-mer|321|28-28-28-23|||28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)3",
            "3-mer|405|28-28-28-28-23|||28-28-28-28-23|28-28-28-28-23|28-28-28-28-23": "(MON135)3",
            "4-mer|316|28-23;28-28-23;28-28-28-23|||28-23|28-28-23|28-28-28-23|28-28-23": "MON51-(MON79)2-MON107",  # Normalized
            "4-mer|428|28-28-28-23|||28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)4",
            "5-mer|535|28-28-28-23|||28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)5",
            "6-mer|642|28-28-28-23|||28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)6",
            "10-mer|1070|28-28-28-23|||28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)10"
        }
        dict_sample_type={}
        with open('./samples_satellite/sample_info') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=6:continue
                sample_name=eachline_arr[0]
                sample_type=eachline_arr[3]
                dict_sample_type[sample_name]=sample_type
        dict_type_group_num={}        
        dict_type_num={} 
        with open ('./samples_satellite/14_chr_image/17_HOR_type/1_chr_plot/0_sample_chr_type_length_order','w') as f2:   
            f2.write(f"sample\tchromosome\tHOR_class\tHOR_class_new\tlength\torder\tsample_type\n")
            with open ('./samples_satellite/14_chr_image/17_HOR_type/0_sample_chr_type_length_order','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    sample,chromosome,HOR_class,length,order=eachline_arr
                    if HOR_class not in dict_change:
                        HOR_class_new='other'
                    else:    
                        HOR_class_new=dict_change[HOR_class]
                    if int(order)>3: continue 
                    if int(length)<200000:continue
                    sample_type=dict_sample_type[sample]
                    if sample_type not in ['East_Asia','Eurasian','America']:continue
                    f2.write(f"{sample}\t{chromosome}\t{HOR_class}\t{HOR_class_new}\t{length}\t{order}\t{sample_type}\n")
                    #
                    if HOR_class_new not in dict_type_group_num:
                        dict_type_group_num[HOR_class_new]={}
                    if sample_type not in dict_type_group_num[HOR_class_new]:    
                        dict_type_group_num[HOR_class_new][sample_type]=0
                    dict_type_group_num[HOR_class_new][sample_type]+=1
                    #
                    if HOR_class_new not in dict_type_num:    
                        dict_type_num[HOR_class_new]=0
                    dict_type_num[HOR_class_new]+=1
        with open ('./samples_satellite/14_chr_image/17_HOR_type/1_chr_plot/1_HORclass_groups_num','w') as f2:  
            f2.write(f"HOR_class_new\tgroup\tnum\n")
            for HOR_class_new,dict_group_num in dict_type_group_num.items():
                for group,num in dict_group_num.items():
                    f2.write(f"{HOR_class_new}\t{group}\t{num}\n")
        with open ('./samples_satellite/14_chr_image/17_HOR_type/1_chr_plot/1_HORclass_num','w') as f2:  
            f2.write(f"HOR_class_new\tnum\n")
            for HOR_class_new,num in dict_type_num.items():
                f2.write(f"{HOR_class_new}\t{num}\n")                    
                                    
        Plot_txt=r"""
    library(ggplot2)
    library(dplyr)
    # Get all command line parameters
    
    ### Monomer
    print("")
    {
      input_file=read.table('0_sample_chr_type_length_order', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
      
    
      # Define color values
    color_values <- c(
        "(MON107)2" = "#c2c2d6",    #1050
        "(MON107)3" = "#3333ff",    #272
        "(MON107)4"= "#00ccff",     #177
        "(MON107)5"= "#0099cc",     #38
        "(MON107)6"= "#9999ff",
        "(MON107)10"='#000066',
        
        "(MON135)2" = "#00cc99",    
        "(MON135)3" = "#99ff33",    
        
        "MON51-MON107" = "#ffff66",
        "MON79-MON107" = "#ffcc00",
        "MON107-MON135" = "#ff9900",
        "MON107-MON163" = "#ff9999",
        "MON51-(MON107)2" = "#ff99ff",
        "(MON79)2-MON107" = "#ff66cc",
        "MON79-(MON107)2" = "#ff3399",
        "MON79-MON107-MON135" = "#cc00cc",
        "MON51-(MON79)2-MON107" = "red",
        "other"="#cc6600",
        "America"="#33CC33",
        "Eurasian" = "#ff3399",
        "East_Asia" = "#0066ff"
      )
    
    
    
    # Fix chromosome order
    input_file$chromosome <- factor(input_file$chromosome, levels = c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"))
    # Fix sample order    
    input_file$sample <- factor(input_file$sample, levels = c("Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap1","V059.hap2","V063.hap1","V063.hap2","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap1","Chardonnay_hap2","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap1","V074.hap2","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap1","V062.hap2","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap1","PinotNoir_hap2","PinotNoir2_hap1","PinotNoir2_hap2","V069.hap1","V069.hap2","PN40024","PN40024_hap1","PN40024_hap2","V087.hap1","V087.hap2","V061.hap1","V061.hap2","V065.hap1","V065.hap2","V064.hap1","V064.hap2","ThompsonSeedless_hap1","ThompsonSeedless_hap2","ThompsonSeedless2_hap1","ThompsonSeedless2_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap1","V067.hap2","V058.hap1","V058.hap2","V107.hap1","V107.hap2","V108.hap1","V108.hap2","V081.hap1","V081.hap2","V015.hap1","V015.hap2","V105.hap1","V105.hap2","V126.hap1","V126.hap2","V112.hap1","V112.hap2","DavidiiGrape_hap1","DavidiiGrape_hap2","V007.hap1","V007.hap2","V008.hap1","V008.hap2","V012.hap1","V012.hap2","V100.hap1","V100.hap2","V099.hap1","V099.hap2","V098.hap1","V098.hap2","V117.hap1","V117.hap2","V120.hap1","V120.hap2","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap1","V124.hap2","V125.hap1","V125.hap2","WoollyGrape_hap1","WoollyGrape_hap2","V106.hap1","V106.hap2","V102.hap1","V102.hap2","V079.hap1","V079.hap2","V030.hap1","V030.hap2","V031.hap1","V031.hap2","V005.hap1","V005.hap2","V048.hap1","V048.hap2","V050.hap1","V050.hap2","V051.hap1","V051.hap2","V049.hap1","V049.hap2","V038.hap1","V038.hap2","V043.hap1","V043.hap2","V052.hap1","V052.hap2","V034.hap1","V034.hap2","V023.hap1","V023.hap2","V032.hap1","V032.hap2","V033.hap1","V033.hap2","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap1","V022.hap2","V040.hap1","V040.hap2","V019.hap1","V019.hap2","V041.hap1","V041.hap2","V077.hap1","V077.hap2","V037.hap1","V037.hap2","V096.hap1","V096.hap2","ShineMuscat_hap1","ShineMuscat_hap2","V036.hap1","V036.hap2","V076.hap1","V076.hap2","V055.hap1","V055.hap2","V018.hap1","V018.hap2","V072.hap1","V072.hap2"))
    
      # Create plot object
      p = ggplot() 
        # 1. Draw background (background rectangles)
    p=p+  geom_point(data=input_file,
            aes(x=sample,y=log(length,10),color = HOR_class_new),
            size =2.5,  # Point size
            shape=15,stroke = 0
          ) 
    p=p+  geom_point(data=input_file,
            aes(x=sample,y=8,color = sample_type),
            size =1.5,  # Point size
            shape=15,stroke = 0
          ) 
      
                       
    p=p+facet_grid(. ~ chromosome, scales = "free_x", space = "free_x")     
       p=p+theme_classic() +         
        theme(
        panel.spacing.x = unit(0.5, "cm"), 
          #axis.ticks.y = element_blank(),
          #axis.text.y = element_blank(),
          legend.position = "none",
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()
        ) + scale_y_continuous(
        breaks = log10(c(200000, 300000, 500000,  1000000, 2000000, 4000000, 6000000, 8000000, 10000000)),  # Specify tick positions (after log10)
        labels = c("200k", "300k", "500k",  "1M", "2M", "4M", "6M", "8M", "10M"),  # Custom labels
        name = "Value (log10 scale)"  # Y-axis title
      ) +
        scale_color_manual(values = color_values, drop = FALSE)
        #scale_fill_manual(values = color_values, drop = FALSE)
      
      pdf(file = paste0('0_sample_chr_type_length_order', ".pdf"), width = 50 / 2.54, height = 30 / 2.54)
      print(p)
      dev.off()
      }
      
    # Create plot object
    input_file=read.table('1_HORclass_groups', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
      # Draw faceted stacked bar chart (instead of pie chart facets)
    p <- ggplot(data = input_file, aes(x = HOR_class_new, y = num, fill = group)) +
      geom_bar(stat = "identity", position = "stack", width = 0.8) +  # Stacked bar chart
      labs(title = "Faceted stacked bar chart") +
      theme_classic() +
      theme(
        #axis.ticks.y = element_blank(),
        #axis.text.y = element_blank(),
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
      ) +
      facet_grid(. ~ HOR_class_new, scales = "free_x", space = "free_x")  # Facets
     
    # Save as PDF
    pdf(file = "1_HORclass_groups_bar.pdf", width = 25/2.54, height = 6/2.54)
    print(p)
    dev.off()
    
    
    
    values= color_values
     #print(values)
    # Create a dummy data frame for plotting
    df <- data.frame(
      category = factor(names(values), levels = names(values)),  # Ensure category is factor and in order
      dummy = 1  # Create a dummy column for plotting
    )
    
    # Create a plot object, only for displaying legend
    p <- ggplot(df, aes(x = dummy, y = dummy, fill = category)) +
      geom_point(shape = 22, size = 3, color = "black") +  # Use square points to display legend
      scale_fill_manual(values = values, name = "Categories") +  # Set colors and legend title
      scale_shape_manual(values = 22) +  # Ensure all points are square
      theme_void() +  # Hide all plot elements
      theme(legend.position = "bottom", legend.key.size = unit(1, "cm"))  # Set legend position and size
     
    # Print legend
      # Save as PDF
      pdf("legend.pdf", width = 20 / 2.54, height = 10 / 2.54)
      print(p)
      dev.off()
           
            """
        with open(f'./samples_satellite//14_chr_image/17_HOR_type/1_chr_plot/0_sample_chr_type_length_order.R','w') as f:
            f.write(Plot_txt)
        os.chdir(f'./samples_satellite//14_chr_image/17_HOR_type/1_chr_plot/')
        subprocess.run([f'Rscript 0_sample_chr_type_length_order.R  '], shell=True)  #>null 2>&1 
        
        
        os.chdir('../../../../')

        
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.18" or argv1=="step14.18_ChrXX_HORclass":        
        print('Plotting the results from step14.17 on all chromosomes, i.e., the file from step14.15p_ChrXX_unit3, plus the step14.17 classification')
        if  os.path.exists('./samples_satellite/14_chr_image/18_ChrXX')==False: 
            subprocess.run(["mkdir ./samples_satellite/14_chr_image/18_ChrXX"], shell=True)             
        chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
        
        print('Plotting ChrXX for all samples')
        for one_chr in chr_list:
            subprocess.run([f"mkdir ./samples_satellite/14_chr_image/18_ChrXX/{one_chr}"], shell=True)    
            
            subprocess.run([f"cp ./samples_satellite/14_chr_image/15_ChrXX/{one_chr}/0_monomer ./samples_satellite/14_chr_image/18_ChrXX/{one_chr}/0_monomer "], shell=True)    
            with open(f"./samples_satellite/14_chr_image/18_ChrXX/{one_chr}/0_hor",'w') as f2:
                f2.write("mer\tsegment\tsegment_type\tsegment_type_num\tsegment_length\tHOR_type\tcirc_start\tcirc_end\tpos_start\tpos_end\tmarkserial_start\tmarkserial_end\tlayer\tfather_layer\tpos_length\tmarkserial_num\tHOR_repeat_num\tblock_name\tblock_start\tblock_end\tsample\tplot_y\tHOR_class_new\n")
                with open(f"./samples_satellite/14_chr_image/15_ChrXX/{one_chr}/0_hor",'r') as f:
                    next(f)
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        ##
                        mer=	                eachline_arr[0]
                        segment=	            eachline_arr[1]
                        segment_type=           eachline_arr[2]
                        if "." not in segment:
                            segment_rough_length=0
                            for one_subunit_length in segment.replace('-','|').replace('||','|').split('|'):
                                if len(one_subunit_length)==0:continue
                                segment_rough_length+=int(one_subunit_length)
                            HOR_class=f"{mer}-mer|{segment_rough_length}|{segment_type}|||{segment}"
                            length=         int(eachline_arr[14])
                        ##
                        dict_change = {
                            "2-mer|158|28-23;28-28-28-23|||28-23|28-28-28-23": "MON51-MON107",
                            "2-mer|186|28-28-23;28-28-28-23|||28-28-23|28-28-28-23": "MON79-MON107",
                            "2-mer|214|28-28-28-23|||28-28-28-23|28-28-28-23": "(MON107)2",
                            "2-mer|242|28-28-28-23;28-28-28-28-23|||28-28-28-23|28-28-28-28-23": "MON107-MON135",
                            "2-mer|270|28-28-28-23;28-28-28-28-28-23|||28-28-28-23|28-28-28-28-28-23": "MON107-MON163",
                            "2-mer|270|28-28-28-28-23|||28-28-28-28-23|28-28-28-28-23": "(MON135)2",
                            "3-mer|265|28-23;28-28-28-23|||28-23|28-28-28-23|28-28-28-23": "MON51-(MON107)2",           # After normalization
                            "3-mer|265|28-28-23;28-28-28-23|||28-28-23|28-28-23|28-28-28-23": "(MON79)2-MON107",      # After normalization
                            "3-mer|293|28-28-23;28-28-28-23|||28-28-23|28-28-28-23|28-28-28-23": "MON79-(MON107)2",      # After normalization (note difference from above)
                            "3-mer|321|28-28-23;28-28-28-23;28-28-28-28-23|||28-28-23|28-28-28-23|28-28-28-28-23": "MON79-MON107-MON135",
                            "3-mer|321|28-28-28-23|||28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)3",
                            "3-mer|405|28-28-28-28-23|||28-28-28-28-23|28-28-28-28-23|28-28-28-28-23": "(MON135)3",
                            "4-mer|316|28-23;28-28-23;28-28-28-23|||28-23|28-28-23|28-28-28-23|28-28-23": "MON51-(MON79)2-MON107",  # After normalization
                            "4-mer|428|28-28-28-23|||28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)4",
                            "5-mer|535|28-28-28-23|||28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)5",
                            "6-mer|642|28-28-28-23|||28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)6",
                            "10-mer|1070|28-28-28-23|||28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23|28-28-28-23": "(MON107)10"
                        }
                        if HOR_class in dict_change:
                            HOR_class_new=dict_change[HOR_class]
                        else:
                            continue
                            #HOR_class_new='.'
                        f2.write(f"{eachline}\t{HOR_class_new}\n")
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.18" or argv1=="step14.18p_ChrXX_HORclass":        
        chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
        for one_chr in chr_list:
            Plot_txt=r"""
library(ggplot2)
library(dplyr)
# Get all command line arguments

### Monomer
print("")
{
  input_file3=read.table('0_hor', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    min3 <- min(min(input_file3$block_start, na.rm = TRUE),min(input_file3$block_end, na.rm = TRUE))  
    max3 <- max(max(input_file3$block_start, na.rm = TRUE),max(input_file3$block_end, na.rm = TRUE))
  sample_num=-min(input_file3$plot_y, na.rm = TRUE)
  #print(sample_num)
  

  input_file2=read.table('0_monomer', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
filtered_input2 <- input_file2 %>% 
    mutate(
      category = case_when(
        monomer_type == 'CEN107' ~ "CEN107",
        monomer_type == 'CEN107-like' ~ "CEN107-like",
        monomer_type == 'CEN135' ~ "CEN135",
        monomer_type == 'CEN135-like' ~ "CEN135-like",
        monomer_type == 'CEN79' ~ "CEN79",
        monomer_type == 'CEN79-like' ~ "CEN79-like",
        monomer_type == 'CEN51' ~ "CEN51",
        monomer_type == 'CEN51-like' ~ "CEN51-like",        
        monomer_type == 'other' ~ "other_monomer",
        TRUE ~ "other_monomer"
      )
    )
  
  
  
  filtered_input_file3 <- input_file3 %>% 
    filter(pos_length > 300 & HOR_repeat_num>=3)
    
    
  B1 <- filtered_input_file3 %>% filter(father_layer == 1)
  B2 <- filtered_input_file3 %>% filter(father_layer == 2)
  B3 <- filtered_input_file3 %>% filter(father_layer == 3)
 
 # Extract sub-table to remove duplicates
  block_start_block_end_plot_y <- input_file2 %>%
  select(block_start,block_end,plot_y) %>%       
  distinct()  
  
  sample_plot_y <- input_file2 %>%
  select(sample,plot_y) %>%      
  distinct()  
  
  
  
  # Define color values
color_values <- c(
    "(MON107)2" = "#c2c2d6",    #1050
    "(MON107)3" = "#3333ff",    #272
    "(MON107)4"= "#00ccff",     #177
    "(MON107)5"= "#0099cc",     #38
    "(MON107)6"= "#9999ff",
    "(MON107)10"='#000066',
    
    "(MON135)2" = "#00cc99",    
    "(MON135)3" = "#99ff33",    
    
    "MON51-MON107" = "#ffff66",
    "MON79-MON107" = "#ffcc00",
    "MON107-MON135" = "#ff9900",
    "MON107-MON163" = "#ff9999",
    "MON51-(MON107)2" = "#ff99ff",
    "(MON79)2-MON107" = "#ff66cc",
    "MON79-(MON107)2" = "#ff3399",
    "MON79-MON107-MON135" = "#cc00cc",
    "MON51-(MON79)2-MON107" = "red",
   
    "America"="#33CC33",
    "Eurasian" = "#ff3399",
    "East_Asia" = "#0066ff",

    "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699",  
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3", 
    "other_monomer" = "black"
     #   TRUE ~ "Other"
  )

    
      
  # Create plot object
  p = ggplot() 


p =p+  geom_rect(data = filtered_input2,
                   aes(xmin = pos_start,xmax = pos_end, ymin =  plot_y+0.1,ymax =  plot_y+0.2,fill = category))       ######## Monomer


  p =p+  geom_rect(data=block_start_block_end_plot_y,
                   aes(xmin = block_start,xmax = block_end,ymin = plot_y-0.2,ymax = plot_y),fill = '#f2f2f2',) 
  p =p+  geom_rect(data=block_start_block_end_plot_y,
                   aes(xmin = block_start,xmax = block_end,ymin = plot_y-0.45,ymax = plot_y-0.25),fill = '#f2f2f2',) 
  p =p+  geom_rect(data=block_start_block_end_plot_y,
                   aes(xmin = block_start,xmax = block_end,ymin = plot_y-0.7,ymax = plot_y-0.5),fill = '#f2f2f2',) #6A356C
  p =p+  geom_rect(data = B1,
                   aes(xmin = pos_start,xmax = pos_end, ymin = plot_y-0.2,ymax = plot_y,fill = HOR_class_new)) 
  p =p+  geom_rect(data = B2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = plot_y-0.45,ymax = plot_y-0.25,fill = HOR_class_new)) 
  p =p+  geom_rect(data = B3,
                   aes(xmin = pos_start,xmax = pos_end, ymin = plot_y-0.7,ymax = plot_y-0.5,fill = HOR_class_new)) 
  
  ## Labels
p <- p +  geom_text(data = sample_plot_y, aes(x = min3, y = plot_y-0.25, label = sample), size =10)


  ##
  ###################################

                   
    
   p=p+theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  
 height_revise=max(10,sample_num)
  # Save as PDF
  pdf(file = paste0('hor_unit3', ".pdf"), width = 100 / 2.54, height = height_revise / 2.54)
  print(p)
  dev.off()
  
}            
        """
            with open(f'./samples_satellite//14_chr_image/18_ChrXX/{one_chr}/hor.R','w') as f:
                f.write(Plot_txt)
            os.chdir(f'./samples_satellite//14_chr_image/18_ChrXX/{one_chr}')
            subprocess.run([f'Rscript hor.R  '], shell=True)  #>null 2>&1 
            subprocess.run([f'mv hor_unit3.pdf {one_chr}_unit3_HORclass.pdf  '], shell=True) 
            
            os.chdir('../../../../')                                   
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.18" or argv1=="step14.18p_ChrXX_HORclass_compress":
        print('Compress three layers into one layer for display')
        chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
        for one_chr in chr_list:
            Plot_txt=r"""
library(ggplot2)
library(dplyr)
# Get all command line arguments

### Monomer
print("")
{
  input_file3=read.table('0_hor', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    min3 <- min(min(input_file3$block_start, na.rm = TRUE),min(input_file3$block_end, na.rm = TRUE))  
    max3 <- max(max(input_file3$block_start, na.rm = TRUE),max(input_file3$block_end, na.rm = TRUE))
  sample_num=-min(input_file3$plot_y, na.rm = TRUE)
  #print(sample_num)
  

  input_file2=read.table('0_monomer', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
filtered_input2 <- input_file2 %>% 
    mutate(
      category = case_when(
        monomer_type == 'CEN107' ~ "CEN107",
        monomer_type == 'CEN107-like' ~ "CEN107-like",
        monomer_type == 'CEN135' ~ "CEN135",
        monomer_type == 'CEN135-like' ~ "CEN135-like",
        monomer_type == 'CEN79' ~ "CEN79",
        monomer_type == 'CEN79-like' ~ "CEN79-like",
        monomer_type == 'CEN51' ~ "CEN51",
        monomer_type == 'CEN51-like' ~ "CEN51-like",        
        monomer_type == 'other' ~ "other_monomer",
        TRUE ~ "other_monomer"
      )
    )
  
  
  
  filtered_input_file3 <- input_file3 %>% 
    filter(pos_length > 300 & HOR_repeat_num>=3)
    
    
  B1 <- filtered_input_file3 %>% filter(father_layer == 1)
  B2 <- filtered_input_file3 %>% filter(father_layer == 2)
  B3 <- filtered_input_file3 %>% filter(father_layer == 3)
 
 # Extract sub-table to remove duplicates
  block_start_block_end_plot_y <- input_file2 %>%
  select(block_start,block_end,plot_y) %>%       
  distinct()  
  
  sample_plot_y <- input_file2 %>%
  select(sample,plot_y) %>%      
  distinct()  
  
  
  
  # Define color values
color_values <- c(
    "(MON107)2" = "#c2c2d6",    #1050
    "(MON107)3" = "#3333ff",    #272
    "(MON107)4"= "#00ccff",     #177
    "(MON107)5"= "#0099cc",     #38
    "(MON107)6"= "#9999ff",
    "(MON107)10"='#000066',
    
    "(MON135)2" = "#00cc99",    
    "(MON135)3" = "#99ff33",    
    
    "MON51-MON107" = "#ffff66",
    "MON79-MON107" = "#ffcc00",
    "MON107-MON135" = "#ff9900",
    "MON107-MON163" = "#ff9999",
    "MON51-(MON107)2" = "#ff99ff",
    "(MON79)2-MON107" = "#ff66cc",
    "MON79-(MON107)2" = "#ff3399",
    "MON79-MON107-MON135" = "#cc00cc",
    "MON51-(MON79)2-MON107" = "red",
   
    "America"="#33CC33",
    "Eurasian" = "#ff3399",
    "East_Asia" = "#0066ff",

    "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699",  
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3", 
    "other_monomer" = "black"
     #   TRUE ~ "Other"
  )

    
      
  # Create plot object
  p = ggplot() 


#p =p+  geom_rect(data = filtered_input2,
#                   aes(xmin = pos_start,xmax = pos_end, ymin =  plot_y+0.1,ymax =  plot_y+0.2,fill = category))       ######## Monomer


  p =p+  geom_rect(data=block_start_block_end_plot_y,
                   aes(xmin = block_start,xmax = block_end,ymin = plot_y-0.8,ymax = plot_y),fill = '#f2f2f2',) 
  
  p =p+  geom_rect(data = B1,
                   aes(xmin = pos_start,xmax = pos_end, ymin = plot_y-0.8,ymax = plot_y,fill = HOR_class_new)) 
  p =p+  geom_rect(data = B2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = plot_y-0.8,ymax = plot_y,fill = HOR_class_new)) 
  p =p+  geom_rect(data = B3,
                   aes(xmin = pos_start,xmax = pos_end, ymin = plot_y-0.8,ymax = plot_y,fill = HOR_class_new)) 
  
  ## Labels
p <- p +  geom_text(data = sample_plot_y, aes(x = min3, y = plot_y-0.25, label = sample), size =10)


  ##
  ###################################

                   
    
   p=p+theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  
 height_revise=max(10,sample_num)
  # Save as PDF
  pdf(file = paste0('hor_unit3', ".pdf"), width = 100 / 2.54, height = height_revise / 2.54)
  print(p)
  dev.off()
  
}            
        """
            with open(f'./samples_satellite//14_chr_image/18_ChrXX/{one_chr}/hor.R','w') as f:
                f.write(Plot_txt)
            os.chdir(f'./samples_satellite//14_chr_image/18_ChrXX/{one_chr}')
            subprocess.run([f'Rscript hor.R  '], shell=True)  #>null 2>&1 
            subprocess.run([f'mv hor_unit3.pdf {one_chr}_unit3_HORclass_compress.pdf  '], shell=True) 
            
            os.chdir('../../../../')                                           
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.18" or argv1=="step14.18_ChrXX_HORclass_compress_shift":
        print('Using the revised starting point as the plot starting point')
        dict_samplechr_shift={}
        with open('./samples_satellite/2_good_regions_main','r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split("\t")
                dict_samplechr_shift[f"{eachline_arr[0]}|||{eachline_arr[1]}"]=eachline_arr[2]
        chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
        for one_chr in chr_list:
            with open(f'./samples_satellite//14_chr_image/18_ChrXX/{one_chr}/0_hor_shift','w') as f2:
                f2.write("mer\tsegment\tsegment_type\tsegment_type_num\tsegment_length\tHOR_type\tcirc_start\tcirc_end\tpos_start\tpos_end\tmarkserial_start\tmarkserial_end\tlayer\tfather_layer\tpos_length\tmarkserial_num\tHOR_repeat_num\tblock_name\tblock_start\tblock_end\tsample\tplot_y\tHOR_class_new\tshift_delta\n")
                with open(f'./samples_satellite//14_chr_image/18_ChrXX/{one_chr}/0_hor','r') as f:
                    next(f)
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=line.split("\t")
                        sample=     eachline_arr[20]
              
                        chromosome= eachline_arr[17].split(':')[1][:-1].replace('region_','Chr')
                        shift=dict_samplechr_shift[f"{sample}|||{chromosome}"]
                        f2.write(f'{eachline}\t{shift}\n')
            with open(f'./samples_satellite//14_chr_image/18_ChrXX/{one_chr}/0_monomer_shift','w') as f2:
                f2.write("serial\tcirc_start\tcirc_end\tpos_start\tpos_end\tlength\tmonomer_type\tblock_name\tblock_start\tblock_end\tsample\tplot_y\tshift_delta\n")
                with open(f'./samples_satellite//14_chr_image/18_ChrXX/{one_chr}/0_monomer','r') as f:
                    next(f)
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=line.split("\t")
                        sample=     eachline_arr[10]
              
                        chromosome= eachline_arr[7].split(':')[1][:-1].replace('region_','Chr')
                        shift=dict_samplechr_shift[f"{sample}|||{chromosome}"]
                        f2.write(f'{eachline}\t{shift}\n')            
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.18" or argv1=="step14.18p_ChrXX_HORclass_compress_shift":   
        print('Compress three layers into one layer for display')
        chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
        for one_chr in chr_list:
            Plot_txt=r"""
library(ggplot2)
library(dplyr)
# Get all command line arguments

### Monomer
print("")
{
  input_file3=read.table('0_hor_shift', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    min3 <- min(min(input_file3$block_start, na.rm = TRUE),min(input_file3$block_end, na.rm = TRUE))  
    max3 <- max(max(input_file3$block_start, na.rm = TRUE),max(input_file3$block_end, na.rm = TRUE))
  sample_num=-min(input_file3$plot_y, na.rm = TRUE)
  #print(sample_num)
  

  input_file2=read.table('0_monomer_shift', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
filtered_input2 <- input_file2 %>% 
    mutate(
      category = case_when(
        monomer_type == 'CEN107' ~ "CEN107",
        monomer_type == 'CEN107-like' ~ "CEN107-like",
        monomer_type == 'CEN135' ~ "CEN135",
        monomer_type == 'CEN135-like' ~ "CEN135-like",
        monomer_type == 'CEN79' ~ "CEN79",
        monomer_type == 'CEN79-like' ~ "CEN79-like",
        monomer_type == 'CEN51' ~ "CEN51",
        monomer_type == 'CEN51-like' ~ "CEN51-like",        
        monomer_type == 'other' ~ "other_monomer",
        TRUE ~ "other_monomer"
      )
    )
  
  
  
  filtered_input_file3 <- input_file3 %>% 
    filter(pos_length > 300 & HOR_repeat_num>=3)
    
    
  B1 <- filtered_input_file3 %>% filter(father_layer == 1)
  B2 <- filtered_input_file3 %>% filter(father_layer == 2)
  B3 <- filtered_input_file3 %>% filter(father_layer == 3)
 
 # Extract sub-table to remove duplicates
  block_start_block_end_plot_y <- input_file2 %>%
  select(block_start,block_end,plot_y,shift_delta) %>%       
  distinct()  
  
  sample_plot_y <- input_file2 %>%
  select(sample,plot_y,shift_delta) %>%      
  distinct()  
  
  
  
  # Define color values
color_values <- c(
    "(MON107)2" = "#c2c2d6",    #1050
    "(MON107)3" = "#3333ff",    #272
    "(MON107)4"= "#00ccff",     #177
    "(MON107)5"= "#0099cc",     #38
    "(MON107)6"= "#9999ff",
    "(MON107)10"='#000066',
    
    "(MON135)2" = "#00cc99",    
    "(MON135)3" = "#99ff33",    
    
    "MON51-MON107" = "#ffff66",
    "MON79-MON107" = "#ffcc00",
    "MON107-MON135" = "#ff9900",
    "MON107-MON163" = "#ff9999",
    "MON51-(MON107)2" = "#ff99ff",
    "(MON79)2-MON107" = "#ff66cc",
    "MON79-(MON107)2" = "#ff3399",
    "MON79-MON107-MON135" = "#cc00cc",
    "MON51-(MON79)2-MON107" = "red",
   
    "America"="#33CC33",
    "Eurasian" = "#ff3399",
    "East_Asia" = "#0066ff",

    "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699",  
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3", 
    "other_monomer" = "black"
     #   TRUE ~ "Other"
  )

    
      
  # Create plot object
  p = ggplot() 


#p =p+  geom_rect(data = filtered_input2,
#                   aes(xmin = pos_start-shift_delta,xmax = pos_end-shift_delta, ymin =  plot_y+0.1,ymax =  plot_y+0.2,fill = category))       ######## Monomer


  p =p+  geom_rect(data=block_start_block_end_plot_y,
                   aes(xmin = block_start-shift_delta,xmax = block_end-shift_delta,ymin = plot_y-0.8,ymax = plot_y),fill = '#f2f2f2',) 
  
  p =p+  geom_rect(data = B1,
                   aes(xmin = pos_start-shift_delta,xmax = pos_end-shift_delta, ymin = plot_y-0.8,ymax = plot_y,fill = HOR_class_new)) 
  p =p+  geom_rect(data = B2,
                   aes(xmin = pos_start-shift_delta,xmax = pos_end-shift_delta, ymin = plot_y-0.8,ymax = plot_y,fill = HOR_class_new)) 
  p =p+  geom_rect(data = B3,
                   aes(xmin = pos_start-shift_delta,xmax = pos_end-shift_delta, ymin = plot_y-0.8,ymax = plot_y,fill = HOR_class_new)) 
  
  ## Labels
p <- p +  geom_text(data = sample_plot_y, aes(x = min3-shift_delta, y = plot_y-0.25, label = sample), size =10)


  ##
  ###################################

                   
    
   p=p+theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  
 height_revise=max(10,sample_num)
  # Save as PDF
  pdf(file = paste0('hor_unit3', ".pdf"), width = 100 / 2.54, height = 50 / 2.54)
  print(p)
  dev.off()
  
}            
        """
            with open(f'./samples_satellite//14_chr_image/18_ChrXX/{one_chr}/hor.R','w') as f:
                f.write(Plot_txt)
            os.chdir(f'./samples_satellite//14_chr_image/18_ChrXX/{one_chr}')
            subprocess.run([f'Rscript hor.R  '], shell=True)  #>null 2>&1 
            subprocess.run([f'mv hor_unit3.pdf {one_chr}_unit3_HORclass_compress_shift.pdf  '], shell=True) 
            
            os.chdir('../../../../')
            
    #The code for other Chr is similar.        
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.18" or argv1=="step14.18p_ChrXX_HORclass_compress_shift_DIY_plus_Chr3":   
        print('Compress three layers into one layer for display, customize output samples and chromosomes')
        chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
        chr_list=["Chr3"]
        for one_chr in chr_list:
            SAMPLE_used=["V077.hap2","V037.hap1","V036.hap2","V043.hap1","V041.hap1","V040.hap2","V022.hap2","V020.hap1","V023.hap1","V032.hap2","V019.hap1","V053.hap2","V052.hap2","V038.hap1","V049.hap1","V020.hap2","V040.hap1","V038.hap2","V022.hap1","V049.hap2","V041.hap2","V055.hap2","V051.hap1","V050.hap2","V051.hap2","V037.hap2","V050.hap1","V033.hap2","V043.hap2","V019.hap2","V055.hap1","V033.hap1","V048.hap2","V032.hap1","V053.hap1","V052.hap1","V034.hap2","V048.hap1","V023.hap2","V007.hap2","V007.hap1","V008.hap2","V008.hap1","V012.hap2","V012.hap1","V108.hap2","V108.hap1","V120.hap1","V081.hap2","V117.hap2","V120.hap2","V117.hap1","V081.hap1","V031.hap2","V031.hap1","WoollyGrape_hap2","WoollyGrape_hap1","V100.hap1","V098.hap2","V100.hap2","V098.hap1","V079.hap2","V079.hap1","DavidiiGrape_hap2","DavidiiGrape_hap1","V126.hap2","V106.hap2","V106.hap1","V102.hap2","V102.hap1","V126.hap1","V105.hap2","V105.hap1","V030.hap2","V030.hap1","V005.hap2","V005.hap1","V124.hap2","V124.hap1","V123.hap2","V123.hap1","V062.hap2","V062.hap1","V065.hap2","V065.hap1","V087.hap1","V092.hap2","V087.hap2","V092.hap1","V064.hap2","V059.hap2","V063.hap1","V061.hap1","V059.hap1","V061.hap2","V063.hap2","V067.hap1","V064.hap1","Chardonnay_hap1","PN40024_hap2","PN40024_hap1","PinotNoir_hap1","V067.hap2","V060.hap2","PinotNoir_hap2","V060.hap1","Chardonnay_hap2","MuscatHamburg_hap1","V091.hap1","BlackMonukka_hap2","ManicureFinger_hap2","BlackMonukka_hap1","Hongmunage_hap2","Baimunage_hap1","V091.hap2","MuscatHamburg_hap2","ManicureFinger_hap1","Hongmunage_hap1","Baimunage_hap2"]
            SAMPLE_used_str=', '.join(f'"{item}"' for item in SAMPLE_used) 
            with open(f'./samples_satellite//14_chr_image/18_ChrXX/{one_chr}/0_hor_shift_DIY','w') as f2:
                with open(f'./samples_satellite//14_chr_image/18_ChrXX/{one_chr}/0_hor_shift','r') as f:
                    kk=0
                    for line in f:
                        kk+=1
                        eachline=line.strip()
                        if kk==1:f2.write(f"{eachline}\n")
                        eachline_arr=eachline.split('\t')
                        sample=eachline_arr[20]
                        if sample not in SAMPLE_used:continue
                        plot_y=-SAMPLE_used.index(sample)-1
                        newline='\t'.join(eachline_arr[:21])+f'\t{plot_y}\t'+  '\t'.join(eachline_arr[22:])  +'\n'                    
                        f2.write(f"{newline}")
            with open(f'./samples_satellite//14_chr_image/18_ChrXX/{one_chr}/0_monomer_shift_DIY','w') as f2:
                with open(f'./samples_satellite//14_chr_image/18_ChrXX/{one_chr}/0_monomer_shift','r') as f:
                    kk=0
                    for line in f:
                        kk+=1
                        eachline=line.strip()
                        if kk==1:f2.write(f"{eachline}\n")
                        eachline_arr=eachline.split('\t')
                        sample=eachline_arr[10]
                        if sample not in SAMPLE_used:continue
                        plot_y=-SAMPLE_used.index(sample)-1
                        newline='\t'.join(eachline_arr[:11])+f'\t{plot_y}\t'+  eachline_arr[12]  +'\n'                    
                        f2.write(f"{newline}")                        
            Plot_txt=f"""
library(ggplot2)
library(dplyr)
# Get all command line arguments

### Monomer
print("")

  input_file3=read.table('0_hor_shift_DIY', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
    ###
  list_A <- c({SAMPLE_used_str})

# Convert sample to factor and specify levels
input_file3 <- input_file3 %>%
  mutate(sample = factor(sample, levels = list_A))
  
  
    min3 <- min(min(input_file3$block_start, na.rm = TRUE),min(input_file3$block_end, na.rm = TRUE))  
    max3 <- max(max(input_file3$block_start, na.rm = TRUE),max(input_file3$block_end, na.rm = TRUE))
  sample_num=-min(input_file3$plot_y, na.rm = TRUE)
  #print(sample_num)

  input_file2=read.table('0_monomer_shift_DIY', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
filtered_input2 <- input_file2 %>% 
    mutate(
      category = case_when(
        monomer_type == 'CEN107' ~ "CEN107",
        monomer_type == 'CEN107-like' ~ "CEN107-like",
        monomer_type == 'CEN135' ~ "CEN135",
        monomer_type == 'CEN135-like' ~ "CEN135-like",
        monomer_type == 'CEN79' ~ "CEN79",
        monomer_type == 'CEN79-like' ~ "CEN79-like",
        monomer_type == 'CEN51' ~ "CEN51",
        monomer_type == 'CEN51-like' ~ "CEN51-like",        
        monomer_type == 'other' ~ "other_monomer",
        TRUE ~ "other_monomer"
      )
    )

  
  filtered_input_file3 <- input_file3 %>% 
    filter(pos_length > 300 & HOR_repeat_num>=3)
    
    #print(input_file2)
    
  B1 <- filtered_input_file3 %>% filter(father_layer == 1)
  B2 <- filtered_input_file3 %>% filter(father_layer == 2)
  B3 <- filtered_input_file3 %>% filter(father_layer == 3)
 
 # Extract sub-table to remove duplicates
  block_start_block_end_plot_y <- input_file2 %>%
  select(block_start,block_end,plot_y,shift_delta,sample) %>%       
  distinct()  %>% 
   mutate(sample = factor(sample, levels = list_A))   # Add this line
  
  sample_plot_y <- input_file2 %>%
  select(sample,plot_y,shift_delta) %>%      
  distinct()   %>%
   mutate(sample = factor(sample, levels = list_A))   # Add this line
  

  
  
  # Define color values
color_values <- c(
    "(MON107)2" = "#c2c2d6",    #1050
    "(MON107)3" = "#3333ff",    #272
    "(MON107)4"= "#00ccff",     #177
    "(MON107)5"= "#0099cc",     #38
    "(MON107)6"= "#9999ff",
    "(MON107)10"='#000066',
    
    "(MON135)2" = "#00cc99",    
    "(MON135)3" = "#99ff33",    
    
    "MON51-MON107" = "#ffff66",
    "MON79-MON107" = "#ffcc00",
    "MON107-MON135" = "#ff9900",
    "MON107-MON163" = "#ff9999",
    "MON51-(MON107)2" = "#ff99ff",
    "(MON79)2-MON107" = "#ff66cc",
    "MON79-(MON107)2" = "#ff3399",
    "MON79-MON107-MON135" = "#cc00cc",
    "MON51-(MON79)2-MON107" = "red",
   
    "America"="#33CC33",
    "Eurasian" = "#ff3399",
    "East_Asia" = "#0066ff",

    "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699",  
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3", 
    "other_monomer" = "black"
     #   TRUE ~ "Other"
  )

    
     
  # Create plot object
  p = ggplot() 


#p =p+  geom_rect(data = filtered_input2,
#                   aes(xmin = pos_start-shift_delta,xmax = pos_end-shift_delta, ymin =  plot_y+0.1,ymax =  plot_y+0.2,fill = category))       ######## Monomer


 p =p+  geom_rect(data=block_start_block_end_plot_y,
                   aes(xmin = block_start-shift_delta,xmax = block_end-shift_delta),fill = '#f2f2f2',ymin = 0,ymax = 1) 
  
  p =p+  geom_rect(data = B1,
                   aes(xmin = pos_start-shift_delta,xmax = pos_end-shift_delta, fill = HOR_class_new),ymin = 0,ymax = 1) 
  p =p+  geom_rect(data = B2,
                   aes(xmin = pos_start-shift_delta,xmax = pos_end-shift_delta, fill = HOR_class_new),ymin = 0,ymax = 1) 
  p =p+  geom_rect(data = B3,
                   aes(xmin = pos_start-shift_delta,xmax = pos_end-shift_delta, fill = HOR_class_new),ymin = 0,ymax = 1) 
  
  ## Labels
p <- p +  geom_text(data = sample_plot_y, aes(x = min3-shift_delta,  label = sample), y=0.4,size =3)
p <- p + facet_wrap(~ sample, ncol = 1, strip.position = "left")

  ##
  ###################################

                   
    
   p=p+theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  
 height_revise=max(10,sample_num)
  # Save as PDF
  pdf(file = paste0('hor_unit3', ".pdf"), width = 100 / 2.54, height = 100 / 2.54)
  print(p)
  dev.off()
  
         
        """
            with open(f'./samples_satellite//14_chr_image/18_ChrXX/{one_chr}/hor.R','w') as f:
                f.write(Plot_txt)
            os.chdir(f'./samples_satellite//14_chr_image/18_ChrXX/{one_chr}')
            subprocess.run([f'Rscript hor.R  '], shell=True)  #>null 2>&1 
            subprocess.run([f'mv hor_unit3.pdf {one_chr}_a_unit3_HORclass_compress_shift_DIY.pdf  '], shell=True) 
            
            os.chdir('../../../../')                                      
    if argv1=="step14.18_part_HORclass":
        print("14.18 image is too large, output a portion")
        sample=sys.argv[2]
        chromosome=sys.argv[3]
        if  os.path.exists('./samples_satellite/14_chr_image/18_ChrXX_part')==False: 
            subprocess.run(["mkdir ./samples_satellite/14_chr_image/18_ChrXX_part"], shell=True)         
        with open(f"./samples_satellite/14_chr_image/18_ChrXX_part/0_hor",'w') as f2:
            f2.write("mer\tsegment\tsegment_type\tsegment_type_num\tsegment_length\tHOR_type\tcirc_start\tcirc_end\tpos_start\tpos_end\tmarkserial_start\tmarkserial_end\tlayer\tfather_layer\tpos_length\tmarkserial_num\tHOR_repeat_num\tblock_name\tblock_start\tblock_end\tsample\tplot_y\tHOR_class_new\n")
            with open(f"./samples_satellite/14_chr_image/18_ChrXX/{chromosome}/0_hor",'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    #mer	segment	segment_type	segment_type_num	segment_length	HOR_type	circ_start	circ_end	pos_start	pos_end	markserial_start	markserial_end	layer	father_layer	pos_length	markserial_num	HOR_repeat_num	block_name	block_start	block_end	sample	plot_y	HOR_class_new
                    one_sample=eachline_arr[-3]
                    if sample!=one_sample:continue
                    f2.write(eachline+'\n')
        with open(f"./samples_satellite/14_chr_image/18_ChrXX_part/0_monomer",'w') as f2:
            f2.write("serial\tcirc_start\tcirc_end\tpos_start\tpos_end\tlength\tmonomer_type\tblock_name\tblock_start\tblock_end\tsample\tplot_y\n")
            with open(f"./samples_satellite/14_chr_image/18_ChrXX/{chromosome}/0_monomer",'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    #serial	circ_start	circ_end	pos_start	pos_end	length	monomer_type	block_name	block_start	block_end	sample	plot_y
                    one_sample=eachline_arr[-2]
                    if sample!=one_sample:continue
                    f2.write(eachline+'\n')                    
        ###
        Plot_txt=r"""
library(ggplot2)
library(dplyr)
# Get all command line arguments

### Monomer
print("")
{
  input_file3=read.table('0_hor', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    min3 <- min(min(input_file3$block_start, na.rm = TRUE),min(input_file3$block_end, na.rm = TRUE))  
    max3 <- max(max(input_file3$block_start, na.rm = TRUE),max(input_file3$block_end, na.rm = TRUE))
  sample_num=-min(input_file3$plot_y, na.rm = TRUE)
  #print(sample_num)
  

  input_file2=read.table('0_monomer', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
filtered_input2 <- input_file2 %>% 
    mutate(
      category = case_when(
        monomer_type == 'CEN107' ~ "CEN107",
        monomer_type == 'CEN107-like' ~ "CEN107-like",
        monomer_type == 'CEN135' ~ "CEN135",
        monomer_type == 'CEN135-like' ~ "CEN135-like",
        monomer_type == 'CEN79' ~ "CEN79",
        monomer_type == 'CEN79-like' ~ "CEN79-like",
        monomer_type == 'CEN51' ~ "CEN51",
        monomer_type == 'CEN51-like' ~ "CEN51-like",        
        monomer_type == 'other' ~ "other_monomer",
        TRUE ~ "other_monomer"
      )
    )
  
  
  
  filtered_input_file3 <- input_file3 %>% 
    filter(pos_length > 300 & HOR_repeat_num>=3)
    
    
  B1 <- filtered_input_file3 %>% filter(father_layer == 1)
  B2 <- filtered_input_file3 %>% filter(father_layer == 2)
  B3 <- filtered_input_file3 %>% filter(father_layer == 3)
 
 # Extract sub-table to remove duplicates
  block_start_block_end_plot_y <- input_file2 %>%
  select(block_start,block_end,plot_y) %>%       
  distinct()  
  
  sample_plot_y <- input_file2 %>%
  select(sample,plot_y) %>%      
  distinct()  
  
  
  
  # Define color values
color_values <- c(
    "(MON107)2" = "#c2c2d6",    #1050
    "(MON107)3" = "#3333ff",    #272
    "(MON107)4"= "#00ccff",     #177
    "(MON107)5"= "#0099cc",     #38
    "(MON107)6"= "#9999ff",
    "(MON107)10"='#000066',
    
    "(MON135)2" = "#00cc99",    
    "(MON135)3" = "#99ff33",    
    
    "MON51-MON107" = "#ffff66",
    "MON79-MON107" = "#ffcc00",
    "MON107-MON135" = "#ff9900",
    "MON107-MON163" = "#ff9999",
    "MON51-(MON107)2" = "#ff99ff",
    "(MON79)2-MON107" = "#ff66cc",
    "MON79-(MON107)2" = "#ff3399",
    "MON79-MON107-MON135" = "#cc00cc",
    "MON51-(MON79)2-MON107" = "red",

    "America"="#33CC33",
    "Eurasian" = "#ff3399",
    "East_Asia" = "#0066ff",

    "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699",
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3",     
    "other_monomer" = "black"
     #   TRUE ~ "Other"
  )

    
      
  # Create plot object
  p = ggplot() 


p =p+  geom_rect(data = filtered_input2,
                   aes(xmin = pos_start,xmax = pos_end, ymin =  plot_y+0.1,ymax =  plot_y+0.2,fill = category))   


  p =p+  geom_rect(data=block_start_block_end_plot_y,
                   aes(xmin = block_start,xmax = block_end,ymin = plot_y-0.2,ymax = plot_y),fill = '#f2f2f2',) 
  p =p+  geom_rect(data=block_start_block_end_plot_y,
                   aes(xmin = block_start,xmax = block_end,ymin = plot_y-0.45,ymax = plot_y-0.25),fill = '#f2f2f2',) 
  p =p+  geom_rect(data=block_start_block_end_plot_y,
                   aes(xmin = block_start,xmax = block_end,ymin = plot_y-0.7,ymax = plot_y-0.5),fill = '#f2f2f2',) #6A356C
  p =p+  geom_rect(data = B1,
                   aes(xmin = pos_start,xmax = pos_end, ymin = plot_y-0.2,ymax = plot_y,fill = HOR_class_new)) 
  p =p+  geom_rect(data = B2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = plot_y-0.45,ymax = plot_y-0.25,fill = HOR_class_new)) 
  p =p+  geom_rect(data = B3,
                   aes(xmin = pos_start,xmax = pos_end, ymin = plot_y-0.7,ymax = plot_y-0.5,fill = HOR_class_new)) 
  
  ## Labels
p <- p +  geom_text(data = sample_plot_y, aes(x = min3, y = plot_y-0.25, label = sample), size =10)


  ##
  ###################################

                   
    
   p=p+theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  
 height_revise=max(10,sample_num)
  # Save as PDF
  pdf(file = paste0('hor_unit3', ".pdf"), width = 100 / 2.54, height = height_revise / 2.54)
  print(p)
  dev.off()
  
}            
        """
        with open(f'./samples_satellite//14_chr_image/18_ChrXX_part/hor.R','w') as f:
            f.write(Plot_txt)
        os.chdir(f'./samples_satellite//14_chr_image/18_ChrXX_part/')
        subprocess.run([f'Rscript hor.R  '], shell=True)  #>null 2>&1 
        subprocess.run([f'mv hor_unit3.pdf {sample}:{chromosome}_unit3_HORclass.pdf  '], shell=True) 
        
        os.chdir('../../../../')
     
    ## 
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.19" or argv1=="step14.19_ChrXX_strand":        
        print('Plot step14.17 results on all chromosomes, i.e., files from step14.15p_ChrXX_unit3, with a 14.17 classification')
        if  os.path.exists('./samples_satellite/14_chr_image/19_ChrXX_strand')==False: 
            subprocess.run(["mkdir ./samples_satellite/14_chr_image/19_ChrXX_strand"], shell=True)             
        chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
        used_samples= ["Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap1","V059.hap2","V063.hap1","V063.hap2","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap1","Chardonnay_hap2","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap1","V074.hap2","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap1","V062.hap2","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap1","PinotNoir_hap2","PinotNoir2_hap1","PinotNoir2_hap2","V069.hap1","V069.hap2","PN40024","PN40024_hap1","PN40024_hap2","V087.hap1","V087.hap2","V061.hap1","V061.hap2","V065.hap1","V065.hap2","V064.hap1","V064.hap2","ThompsonSeedless_hap1","ThompsonSeedless_hap2","ThompsonSeedless2_hap1","ThompsonSeedless2_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap1","V067.hap2","V058.hap1","V058.hap2","V107.hap1","V107.hap2","V108.hap1","V108.hap2","V081.hap1","V081.hap2","V015.hap1","V015.hap2","V105.hap1","V105.hap2","V126.hap1","V126.hap2","V112.hap1","V112.hap2","DavidiiGrape_hap1","DavidiiGrape_hap2","V007.hap1","V007.hap2","V008.hap1","V008.hap2","V012.hap1","V012.hap2","V100.hap1","V100.hap2","V099.hap1","V099.hap2","V098.hap1","V098.hap2","V117.hap1","V117.hap2","V120.hap1","V120.hap2","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap1","V124.hap2","V125.hap1","V125.hap2","WoollyGrape_hap1","WoollyGrape_hap2","V106.hap1","V106.hap2","V102.hap1","V102.hap2","V079.hap1","V079.hap2","V030.hap1","V030.hap2","V031.hap1","V031.hap2","V005.hap1","V005.hap2","V048.hap1","V048.hap2","V050.hap1","V050.hap2","V051.hap1","V051.hap2","V049.hap1","V049.hap2","V038.hap1","V038.hap2","V043.hap1","V043.hap2","V052.hap1","V052.hap2","V034.hap1","V034.hap2","V023.hap1","V023.hap2","V032.hap1","V032.hap2","V033.hap1","V033.hap2","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap1","V022.hap2","V040.hap1","V040.hap2","V019.hap1","V019.hap2","V041.hap1","V041.hap2","V077.hap1","V077.hap2","V037.hap1","V037.hap2","V096.hap1","V096.hap2","ShineMuscat_hap1","ShineMuscat_hap2","V036.hap1","V036.hap2","V076.hap1","V076.hap2","V055.hap1","V055.hap2","V018.hap1","V018.hap2","V072.hap1","V072.hap2"]
        print('Plot ChrXX for all samples')
        for one_chr in chr_list:
            subprocess.run([f"mkdir ./samples_satellite/14_chr_image/19_ChrXX_strand/{one_chr}"], shell=True)    
            with open(f"./samples_satellite/14_chr_image/19_ChrXX_strand/{one_chr}/0_blockstrand",'w') as f2:
                f2.write("hapchr\tcentype\tchr_start\tchr_end\tlength\tstrand\tmatch_percent\tsample\tplot_y\n")
                iii=0
                for one_sample in used_samples:
                    iii-=1
                    if os.path.exists(f'./samples_satellite/14_chr_image/4_otherblock_chr_stat/{one_sample}___{one_chr}')==False:continue
                    with open(f"./samples_satellite/14_chr_image/4_otherblock_chr_stat/{one_sample}___{one_chr}",'r') as f:
                        next(f)
                        for line in f:
                            eachline=line.strip()
                            f2.write(f"{eachline}\t{one_sample}\t{iii}\n")
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.19" or argv1=="step14.19_ChrXX_strand_shift":
        print('Use revised start position as plotting start position')
        dict_samplechr_shift={}
        with open('./samples_satellite/2_good_regions_main_66_103','r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split("\t")
                dict_samplechr_shift[f"{eachline_arr[0]}|||{eachline_arr[1]}"]=[eachline_arr[2],eachline_arr[3]]
                
        chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
        for one_chr in chr_list:
            with open(f'./samples_satellite//14_chr_image/19_ChrXX_strand/{one_chr}/0_blockstrand_shift','w') as f2:
                f2.write("hapchr\tcentype\tchr_start\tchr_end\tlength\tstrand\tmatch_percent\tsample\tplot_y\tcore_start\tcore_end\tmark\n")
                with open(f'./samples_satellite//14_chr_image/19_ChrXX_strand/{one_chr}/0_blockstrand','r') as f:
                    next(f)
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=line.split("\t")
                        chr_start=      int(eachline_arr[2])
                        chr_end=        int(eachline_arr[3])
                        sample=         eachline_arr[7]
                        chromosome=     eachline_arr[0].split('___')[1]
                        info=           dict_samplechr_shift[f"{sample}|||{chromosome}"]
                        core_start=int(info[0])
                        core_end=int(info[1])
                        if core_start<=chr_start and chr_end<=core_end:mark='core'
                        else: 
                            mark='no_core'
                            
                        f2.write(f'{eachline}\t{core_start}\t{core_end}\t{mark}\n')                    
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.19" or argv1=="step14.19p_ChrXX_strand_shift":     
        print('Facet')
        chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
        #chr_list=["Chr18"]
        for one_chr in chr_list:
            Plot_txt=r"""
library(ggplot2)
library(dplyr)
# Get all command line parameters

### Monomer
print("")
{
  input_file4=read.table('0_blockstrand_shift', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
sample_plot_y <- input_file4 %>%
  select(sample,plot_y) %>%      
  distinct()  




  filtered_input4 <- input_file4 %>% 
    filter(centype %in% c("cen_66","cen_103","cen_107","cen_191","cen_355","cen_383","cen_6","cen_68")) %>%    #"cen_6","cen_68",
    mutate( 
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        centype == 'cen_6' & mark=='core'~ "cen_6",
        centype == 'cen_68' & mark=='core'~ "cen_68",
        centype == 'cen_66' & mark=='core'~ "cen_66",
        centype == 'cen_103'& mark=='core' ~ "cen_103",
        centype == 'cen_107'& mark=='core' ~ "cen_107",
        centype == 'cen_191'& mark=='core' ~ "cen_191",
        centype == 'cen_355'& mark=='core' ~ "cen_355",
        centype == 'cen_383'& mark=='core' ~ "cen_383",
        centype == 'cen_6' & mark!='core'~ "cen_6_nocore",
        centype == 'cen_68' & mark!='core'~ "cen_68_nocore",        
        centype == 'cen_66' & mark!='core'~ "cen_66_nocore",
        centype == 'cen_103'& mark!='core' ~ "cen_103_nocore",
        centype == 'cen_107'& mark!='core' ~ "cen_107_nocore",
        centype == 'cen_191'& mark!='core' ~ "cen_191_nocore",
        centype == 'cen_355'& mark!='core' ~ "cen_355_nocore",
        centype == 'cen_383'& mark!='core' ~ "cen_383_nocore"        
        #TRUE ~ "other_satellite"
      )
    )  
    
    filtered_input4_high<- filtered_input4 %>% filter(centype %in% c("cen_6","cen_68"))
    
    filtered_input4 <- filtered_input4 %>% filter(centype %in% c("cen_66","cen_103","cen_107","cen_191","cen_355","cen_383")) 
    
    
    # Define facet order
list_A <- c("Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap1","V059.hap2","V063.hap1","V063.hap2","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap1","Chardonnay_hap2","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap1","V074.hap2","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap1","V062.hap2","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap1","PinotNoir_hap2","PinotNoir2_hap1","PinotNoir2_hap2","V069.hap1","V069.hap2","PN40024","PN40024_hap1","PN40024_hap2","V087.hap1","V087.hap2","V061.hap1","V061.hap2","V065.hap1","V065.hap2","V064.hap1","V064.hap2","ThompsonSeedless_hap1","ThompsonSeedless_hap2","ThompsonSeedless2_hap1","ThompsonSeedless2_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap1","V067.hap2","V058.hap1","V058.hap2","V107.hap1","V107.hap2","V108.hap1","V108.hap2","V081.hap1","V081.hap2","V015.hap1","V015.hap2","V105.hap1","V105.hap2","V126.hap1","V126.hap2","V112.hap1","V112.hap2","DavidiiGrape_hap1","DavidiiGrape_hap2","V007.hap1","V007.hap2","V008.hap1","V008.hap2","V012.hap1","V012.hap2","V100.hap1","V100.hap2","V099.hap1","V099.hap2","V098.hap1","V098.hap2","V117.hap1","V117.hap2","V120.hap1","V120.hap2","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap1","V124.hap2","V125.hap1","V125.hap2","WoollyGrape_hap1","WoollyGrape_hap2","V106.hap1","V106.hap2","V102.hap1","V102.hap2","V079.hap1","V079.hap2","V030.hap1","V030.hap2","V031.hap1","V031.hap2","V005.hap1","V005.hap2","V048.hap1","V048.hap2","V050.hap1","V050.hap2","V051.hap1","V051.hap2","V049.hap1","V049.hap2","V038.hap1","V038.hap2","V043.hap1","V043.hap2","V052.hap1","V052.hap2","V034.hap1","V034.hap2","V023.hap1","V023.hap2","V032.hap1","V032.hap2","V033.hap1","V033.hap2","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap1","V022.hap2","V040.hap1","V040.hap2","V019.hap1","V019.hap2","V041.hap1","V041.hap2","V077.hap1","V077.hap2","V037.hap1","V037.hap2","V096.hap1","V096.hap2","ShineMuscat_hap1","ShineMuscat_hap2","V036.hap1","V036.hap2","V076.hap1","V076.hap2","V055.hap1","V055.hap2","V018.hap1","V018.hap2","V072.hap1","V072.hap2")
 
# Convert sample to factor, and specify levels
filtered_input4 <- filtered_input4 %>%
  mutate(sample = factor(sample, levels = list_A))

filtered_input4_high <- filtered_input4_high %>%
  mutate(sample = factor(sample, levels = list_A))
  
df_plus <- filtered_input4 %>% filter(strand == '+')
df_minus <- filtered_input4 %>% filter(strand == '-')

df_plus_high <- filtered_input4_high %>% filter(strand == '+')
df_minus_high <- filtered_input4_high %>% filter(strand == '-')


  # Define color values
  color_values <- c(
    #"plus"="red",
    #"minus"="blue",
    "cen_6"="#372c97",
    "cen_68"="#eab215",
    "cen_66"="#00ccff",
    "cen_103"="#ff579c",
    "cen_107"="#199200",
    "cen_191"="#9933ff",
    "cen_355"="#995c17",
    "cen_383"="#ffe600",
    "cen_6_nocore"="#372c97",
    "cen_68_nocore"="#eab215",    
    "cen_66_nocore"="#00ccff",
    "cen_103_nocore"="#ff579c",
    "cen_107_nocore"="#199200",
    "cen_191_nocore"="#9933ff",
    "cen_355_nocore"="#995c17",
    "cen_383_nocore"="#ffe600"
     #   TRUE ~ "Other"
  )

        min3 <- min(min(filtered_input4$chr_start, na.rm = TRUE),min(filtered_input4$chr_end, na.rm = TRUE)) 
      
  # Create plot object
  p = ggplot() 

  p =p+  geom_rect(data = df_plus, aes(xmin = chr_start-core_start,xmax = chr_end-core_start, ymin = -0.5,ymax =-0.1 ,fill = category)) 
  p =p+  geom_rect(data = df_minus, aes(xmin = chr_start-core_start,xmax = chr_end-core_start, ymin = -0.9,ymax =-0.5 ,fill = category)) 
  p =p+  geom_rect(data = df_plus_high, aes(xmin = chr_start-core_start,xmax = chr_end-core_start, ymin = -0.5,ymax =-0.05 ,fill = category)) 
  p =p+  geom_rect(data = df_minus_high, aes(xmin = chr_start-core,xmax = chr_end-core_start, ymin = -0.95,ymax =-0.5 ,fill = category)) 
  
  p <- p + annotate("segment", x = 0, xend = 3e7, y = -0.5, yend = -0.5, color = "black", linewidth = 0.001 )
  ##strand other
  #p =p+  geom_rect(data = filtered_input4,
                   #aes(xmin = chr_start,xmax = chr_end, ymin = plot_y-0.9,ymax =plot_y-0.7,fill = category2)) 
  
  ##names
#p <- p +  geom_text(data = sample_plot_y, aes(x = min3, y = plot_y-0.25, label = sample), size =0.3)

p <- p + facet_wrap(~ sample, ncol = 1, strip.position = "left")


   p=p+theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank(),
      #strip.text = element_blank()            #Hide facet titles
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  

  # Save as PDF
  pdf(file = paste0('blockstrand', ".pdf"), width = 50 / 2.54, height = 100 / 2.54)
  print(p)
  dev.off()
  
}            
        """
            with open(f'./samples_satellite/14_chr_image/19_ChrXX_strand/{one_chr}/blockstrand.R','w') as f:
                f.write(Plot_txt)
            print(one_chr)
            os.chdir(f'./samples_satellite//14_chr_image/19_ChrXX_strand/{one_chr}')
            subprocess.run([f'Rscript blockstrand.R  '], shell=True)  #>null 2>&1 
            subprocess.run([f'mv blockstrand.pdf {one_chr}_blockstrand_facet.pdf  '], shell=True) 
            
            os.chdir('../../../../')                                                                               

    
    #The code for other Chr is similar.
    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.19" or argv1=="step14.19p_ChrXX_strand_shift_DIY_Chr3":     
        print('Faceting')
        chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
        chr_list=["Chr3"]
        #chr_list=["Chr18"]
        for one_chr in chr_list:
            Plot_txt=r"""
library(ggplot2)
library(dplyr)
# Get all command line arguments

### Monomer
print("")

  input_file4=read.table('0_blockstrand_shift', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
sample_plot_y <- input_file4 %>%
  select(sample,plot_y) %>%      
  distinct()  




  filtered_input4 <- input_file4 %>% 
    filter(centype %in% c("cen_66","cen_103","cen_107","cen_191","cen_355","cen_383","cen_6","cen_68")) %>%    #"cen_6","cen_68",
    mutate( 
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        centype == 'cen_6' & mark=='core'~ "cen_6",
        centype == 'cen_68' & mark=='core'~ "cen_68",
        centype == 'cen_66' & mark=='core'~ "cen_66",
        centype == 'cen_103'& mark=='core' ~ "cen_103",
        centype == 'cen_107'& mark=='core' ~ "cen_107",
        centype == 'cen_191'& mark=='core' ~ "cen_191",
        centype == 'cen_355'& mark=='core' ~ "cen_355",
        centype == 'cen_383'& mark=='core' ~ "cen_383",
        centype == 'cen_6' & mark!='core'~ "cen_6_nocore",
        centype == 'cen_68' & mark!='core'~ "cen_68_nocore",        
        centype == 'cen_66' & mark!='core'~ "cen_66_nocore",
        centype == 'cen_103'& mark!='core' ~ "cen_103_nocore",
        centype == 'cen_107'& mark!='core' ~ "cen_107_nocore",
        centype == 'cen_191'& mark!='core' ~ "cen_191_nocore",
        centype == 'cen_355'& mark!='core' ~ "cen_355_nocore",
        centype == 'cen_383'& mark!='core' ~ "cen_383_nocore"        
        #TRUE ~ "other_satellite"
      )
    )  
    
    filtered_input4_high<- filtered_input4 %>% filter(centype %in% c("cen_6","cen_68"))
    
    filtered_input4 <- filtered_input4 %>% filter(centype %in% c("cen_66","cen_103","cen_107","cen_191","cen_355","cen_383")) 
    
    
    # Define faceting order
list_A <- c("V077.hap2","V037.hap1","V036.hap2","V043.hap1","V041.hap1","V040.hap2","V022.hap2","V020.hap1","V023.hap1","V032.hap2","V019.hap1","V053.hap2","V052.hap2","V038.hap1","V049.hap1","V020.hap2","V040.hap1","V038.hap2","V022.hap1","V049.hap2","V041.hap2","V055.hap2","V051.hap1","V050.hap2","V051.hap2","V037.hap2","V050.hap1","V033.hap2","V043.hap2","V019.hap2","V055.hap1","V033.hap1","V048.hap2","V032.hap1","V053.hap1","V052.hap1","V034.hap2","V048.hap1","V023.hap2","V007.hap2","V007.hap1","V008.hap2","V008.hap1","V012.hap2","V012.hap1","V108.hap2","V108.hap1","V120.hap1","V081.hap2","V117.hap2","V120.hap2","V117.hap1","V081.hap1","V031.hap2","V031.hap1","WoollyGrape_hap2","WoollyGrape_hap1","V100.hap1","V098.hap2","V100.hap2","V098.hap1","V079.hap2","V079.hap1","DavidiiGrape_hap2","DavidiiGrape_hap1","V126.hap2","V106.hap2","V106.hap1","V102.hap2","V102.hap1","V126.hap1","V105.hap2","V105.hap1","V030.hap2","V030.hap1","V005.hap2","V005.hap1","V124.hap2","V124.hap1","V123.hap2","V123.hap1","V062.hap2","V062.hap1","V065.hap2","V065.hap1","V087.hap1","V092.hap2","V087.hap2","V092.hap1","V064.hap2","V059.hap2","V063.hap1","V061.hap1","V059.hap1","V061.hap2","V063.hap2","V067.hap1","V064.hap1","Chardonnay_hap1","PN40024_hap2","PN40024_hap1","PinotNoir_hap1","V067.hap2","V060.hap2","PinotNoir_hap2","V060.hap1","Chardonnay_hap2","MuscatHamburg_hap1","V091.hap1","BlackMonukka_hap2","ManicureFinger_hap2","BlackMonukka_hap1","Hongmunage_hap2","Baimunage_hap1","V091.hap2","MuscatHamburg_hap2","ManicureFinger_hap1","Hongmunage_hap1","Baimunage_hap2")
filtered_input4 <- filtered_input4[filtered_input4$sample %in% list_A, ]

filtered_input4_high <- filtered_input4_high[filtered_input4_high$sample %in% list_A, ]
# Convert sample to factor and specify levels
filtered_input4 <- filtered_input4 %>%
  mutate(sample = factor(sample, levels = list_A))

filtered_input4_high <- filtered_input4_high %>%
  mutate(sample = factor(sample, levels = list_A))
  
df_plus <- filtered_input4 %>% filter(strand == '+')
df_minus <- filtered_input4 %>% filter(strand == '-')

df_plus_high <- filtered_input4_high %>% filter(strand == '+')
df_minus_high <- filtered_input4_high %>% filter(strand == '-')


  # Define color values
  color_values <- c(
    #"plus"="red",
    #"minus"="blue",
    "cen_6"="#372c97",
    "cen_68"="#eab215",
    "cen_66"="#00ccff",
    "cen_103"="#ff579c",
    "cen_107"="#199200",
    "cen_191"="#9933ff",
    "cen_355"="#995c17",
    "cen_383"="#ffe600",
    "cen_6_nocore"="#372c97",
    "cen_68_nocore"="#eab215",    
    "cen_66_nocore"="#00ccff",
    "cen_103_nocore"="#ff579c",
    "cen_107_nocore"="#199200",
    "cen_191_nocore"="#9933ff",
    "cen_355_nocore"="#995c17",
    "cen_383_nocore"="#ffe600",
    "cen_650"='#ff99ff'
     #   TRUE ~ "Other"
  )

        min3 <- min(min(filtered_input4$chr_start, na.rm = TRUE),min(filtered_input4$chr_end, na.rm = TRUE)) 
      
  # Create plot object
  p = ggplot() 

  p =p+  geom_rect(data = df_plus, aes(xmin = chr_start-core_start,xmax = chr_end-core_start, ymin = -0.5,ymax =-0.1 ,fill = category)) 
  p =p+  geom_rect(data = df_minus, aes(xmin = chr_start-core_start,xmax = chr_end-core_start, ymin = -0.9,ymax =-0.5 ,fill = category)) 
  #p =p+  geom_rect(data = df_plus_high, aes(xmin = chr_start-core_start,xmax = chr_end-core_start, ymin = -0.5,ymax =-0.05 ,fill = category)) 
  #p =p+  geom_rect(data = df_minus_high, aes(xmin = chr_start-core_start,xmax = chr_end-core_start, ymin = -0.95,ymax =-0.5 ,fill = category)) 
  
  p <- p + annotate("segment", x = 0, xend = 3e7, y = -0.5, yend = -0.5, color = "black", linewidth = 0.001 )
  ##strand other
  #p =p+  geom_rect(data = filtered_input4,
                   #aes(xmin = chr_start,xmax = chr_end, ymin = plot_y-0.9,ymax =plot_y-0.7,fill = category2)) 
  
  ## Labels
#p <- p +  geom_text(data = sample_plot_y, aes(x = min3, y = plot_y-0.25, label = sample), size =0.3)

p <- p + facet_wrap(~ sample, ncol = 1, strip.position = "left")


   p=p+theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank(),
      #strip.text = element_blank()            # Hide facet titles
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  

  # Save as PDF
  pdf(file = paste0('blockstrand', ".pdf"), width = 50 / 2.54, height = 100 / 2.54)
  print(p)
  dev.off()
          
        """
            with open(f'./samples_satellite/14_chr_image/19_ChrXX_strand/{one_chr}/blockstrand.R','w') as f:
                f.write(Plot_txt)
            print(one_chr)
            os.chdir(f'./samples_satellite//14_chr_image/19_ChrXX_strand/{one_chr}')
            subprocess.run([f'Rscript blockstrand.R  '], shell=True)  #>null 2>&1 
            subprocess.run([f'mv blockstrand.pdf {one_chr}_a_blockstrand_facet_DIY.pdf  '], shell=True) 
            
            os.chdir('../../../../')                                                                                           

    if argv1=="stepall" or argv1=="step14"  or argv1=="step14.20" or argv1=="step14.20_HORclass_stat":  
        print("step14.20")
        if  os.path.exists('./samples_satellite/14_chr_image/20_HORclass_stat')==False: 
            subprocess.run(["mkdir ./samples_satellite/14_chr_image/20_HORclass_stat"], shell=True)    
            

        for one_chr in chr_list:
            with open(f"./samples_satellite/14_chr_image/19_ChrXX_strand/{one_chr}/0_blockstrand",'r') as f2:
                f2.write("hapchr\tcentype\tchr_start\tchr_end\tlength\tstrand\tmatch_percent\tsample\tplot_y\n")
                iii=0
                for one_sample in used_samples:
                    iii-=1
                    if os.path.exists(f'./samples_satellite/14_chr_image/4_otherblock_chr_stat/{one_sample}___{one_chr}')==False:continue
                    with open(f"./samples_satellite/14_chr_image/4_otherblock_chr_stat/{one_sample}___{one_chr}",'r') as f:
                        next(f)
                        for line in f:
                            eachline=line.strip()
                            f2.write(f"{eachline}\t{one_sample}\t{iii}\n")        
        
    if argv1=="step14.22_partall":        
        sample=sys.argv[2]
        chromosome=sys.argv[3]
        pos=sys.argv[4]
        start=int(pos.split('-')[0])
        end=int(pos.split('-')[1])
        one_hapchr = f"{sample}___{chromosome}"
        print("step14.22_partall")
        if  os.path.exists('./samples_satellite/14_chr_image/22_partall')==False: 
            subprocess.run(["mkdir ./samples_satellite/14_chr_image/22_partall"], shell=True)
        file_name=f"{sample}:{chromosome}:{pos}"
        subprocess.run([f"mkdir ./samples_satellite/14_chr_image/22_partall/{file_name}"], shell=True)      
        ######   
        input_file1=f"./samples_satellite/14_chr_image/1_hor_chr_stat/{one_hapchr}"  # {one_hapchr}
        input_file2=f"./samples_satellite/14_chr_image/2_monomer_chr_stat/{one_hapchr}"
        input_file3=f"./samples_satellite/14_chr_image/3_block_chr_stat/{one_hapchr}"
        input_file4=f"./samples_satellite/14_chr_image/4_otherblock_chr_stat/{one_hapchr}"    
        input_file5=f"./samples_satellite/14_chr_image/18_ChrXX/{chromosome}/0_hor"
        ###
        with open(f"./samples_satellite/14_chr_image/22_partall/{file_name}/1_hor",'w') as f2:
            f2.write("mer\tsegment\tsegment_type\tsegment_type_num\tsegment_length\tHOR_type\tcirc_start\tcirc_end\tpos_start\tpos_end\tmarkserial_start\tmarkserial_end\tlayer\tfather_layer\tpos_length\tmarkserial_num\tHOR_repeat_num\tblock_name\tblock_start\tblock_end\n")
            with open(input_file1,'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    pos_start=int(eachline_arr[8])
                    pos_end=int(eachline_arr[9])
                    if (pos_start>=start and pos_end<=end) or (pos_start<=start and pos_end>=start) or (pos_start<=end and pos_end>=end):
                        f2.write(eachline+'\n')  
        ###
        with open(f"./samples_satellite/14_chr_image/22_partall/{file_name}/2_monomer",'w') as f2:
            f2.write("serial\tcirc_start\tcirc_end\tpos_start\tpos_end\tlength\tmonomer_type\tblock_name\tblock_start\tblock_end\n")
            with open(input_file2,'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    pos_start=int(eachline_arr[3])
                    pos_end=int(eachline_arr[4])
                    if (pos_start>=start and pos_end<=end) or (pos_start<=start and pos_end>=start) or (pos_start<=end and pos_end>=end):
                        f2.write(eachline+'\n')          
        ###
        with open(f"./samples_satellite/14_chr_image/22_partall/{file_name}/3_block",'w') as f2:
            f2.write("hapchr\tblock_name\tblock_start\tblock_end\tstrand\n")
            with open(input_file3,'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    pos_start=int(eachline_arr[2])
                    pos_end=int(eachline_arr[3])
                    if (pos_start>=start and pos_end<=end) or (pos_start<=start and pos_end>=start) or (pos_start<=end and pos_end>=end):
                        f2.write(eachline+'\n')          
        ###
        with open(f"./samples_satellite/14_chr_image/22_partall/{file_name}/4_otherblock",'w') as f2:
            f2.write("hapchr\tcentype\tchr_start\tchr_end\tlength\tstrand\tmatch_percent\n")
            with open(input_file4,'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    pos_start=int(eachline_arr[2])
                    pos_end=int(eachline_arr[3])
                    if (pos_start>=start and pos_end<=end) or (pos_start<=start and pos_end>=start) or (pos_start<=end and pos_end>=end):
                        f2.write(eachline+'\n') 
        #   
        with open(f"./samples_satellite/14_chr_image/22_partall/{file_name}/5_hor",'w') as f2:
            f2.write("mer\tsegment\tsegment_type\tsegment_type_num\tsegment_length\tHOR_type\tcirc_start\tcirc_end\tpos_start\tpos_end\tmarkserial_start\tmarkserial_end\tlayer\tfather_layer\tpos_length\tmarkserial_num\tHOR_repeat_num\tblock_name\tblock_start\tblock_end\tsample\tplot_y\tHOR_class_new\n")
            with open(input_file5,'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    one_sample=eachline_arr[-3]
                    if sample!=one_sample:continue
                    pos_start=int(eachline_arr[8])
                    pos_end=int(eachline_arr[9])
                    if (pos_start>=start and pos_end<=end) or (pos_start<=start and pos_end>=start) or (pos_start<=end and pos_end>=end):
                        f2.write(eachline+'\n')  
        R_txt=r"""
        library(ggplot2)
library(dplyr)
# Get all command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Print arguments (for debugging)
print(args)
input_file_name=args[1]
input_file1=paste0('1_hor')
input_file2=paste0('2_monomer')
input_file3=paste0('3_block')
input_file4=paste0('4_otherblock')
input_file5=paste0('5_hor')
print(input_file1)
print(input_file2)
print(input_file3)
print(input_file4)
print(input_file5)
pos1=as.numeric(args[2])
pos2=as.numeric(args[3])
mode_num=as.numeric(args[4])
print(pos1)
print(pos2)
### Monomer
print("")
{
  input_file1=read.table(input_file1, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file2=read.table(input_file2, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file3=read.table(input_file3, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file4=read.table(input_file4, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  input_file5=read.table(input_file5, skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
  
  min3 <- min(min(input_file3$block_start, na.rm = TRUE),min(input_file3$block_end, na.rm = TRUE))
  max3 <- max(max(input_file3$block_start, na.rm = TRUE),max(input_file3$block_end, na.rm = TRUE))
  min4 <- min(min(input_file4$chr_start, na.rm = TRUE),min(input_file4$chr_end, na.rm = TRUE))
  max4 <- max(max(input_file4$chr_start, na.rm = TRUE),max(input_file4$chr_end, na.rm = TRUE))
  min_all=min(min3,min4)
  max_all=max(max3,max4)
  #print(min_all)
  #print(max_all)
  
  filtered_input_file5 <- input_file5 %>% 
    filter(pos_length > 300 & HOR_repeat_num>=3)

  BB1 <- filtered_input_file5 %>% filter(father_layer == 1)
  BB2 <- filtered_input_file5 %>% filter(father_layer == 2)
  BB3 <- filtered_input_file5 %>% filter(father_layer == 3)
  
  # Filter data
  filtered_input1 <- input_file1 %>% 
  filter(pos_length > 300) %>%
    filter(HOR_repeat_num>=mode_num) %>%     ###########
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        mer == 2 ~ "2",
        mer == 3 ~ "3",
        mer == 4 ~ "4",
        mer == 5 ~ "5",
        mer == 6 ~ "6",
        mer == 7 ~ "7",
        mer == 8 ~ "8",
        mer == 9 ~ "9",
        mer == 10 ~ "10",
        mer == 11 ~ "11",
        mer == 12 ~ "12",
        mer == 13 ~ "13",
        mer == 14 ~ "14",
        mer == 15 ~ "15",
        mer == 16 ~ "16",
        mer == 17 ~ "17",
        mer == 18 ~ "18",
        mer == 19 ~ "19",
        mer == 20 ~ "20",
        mer>20 & mer<=30 ~ "21-30",
        mer>30 ~ "30+",
        TRUE ~ "Other"
      )
    )
  
  
  B1 <- filtered_input1 %>% filter(father_layer == 1)
  B2 <- filtered_input1 %>% filter(father_layer == 2)
  B3 <- filtered_input1 %>% filter(father_layer == 3)
  
  
  
  # Define color values
  color_values <- c(
    "2" = "#b5e48c",
    "3" = "#ff8566",
    "4" = "#ECD678",
    "5" = "#ade8f4",
    "6" = "#cc99ff",
    "7" = "#7f9fcc",  
    "8" = "#ccccff",
    "9" = "#ffccff",
    "10" = "#BAD694",   
    "11" = "#bcdcd6",  
    "12" = "#f49cbb",   
    "13" = "#DDDD7A",   
    "14" = "#bbdefb", 
    "15" = "#339933",  
    "16" = "#ffb4a2",   
    "17" = "#3385ff",  
    "18" = "#00b359",  
    "19" = "#45ABC3", 
    "20" = "#e5e600",
    "21-30"= "#ff80df",
    "30+"= "#cc3300",
    "Other"="white",
    "plus"="red",
    "minus"="blue",
    "cen_6"="#372c97",
    "cen_68"="#eab215",
    "cen_66"="#00ccff",
    "cen_103"="#ff579c",
    "cen_107"="#199200",
    "cen_191"="#9933ff",
    "cen_355"="#995c17",
    "cen_383"="#ffe600",
    #"cen_677"="black",
    #"other_satellite"="black",
    
    "CEN107" = "#066292",
    "CEN107-like" = "#85d2fa",
    "CEN135" = "#51B33F",
    "CEN135-like" = "#bbe4b4", 
    "CEN79" = "#ff9900",
    "CEN79-like" = "#ffd699",  
      "CEN51" = "#ff0000",
      "CEN51-like" = "#ffb3b3",        
    "other_monomer" = "black",
     "(MON107)2" = "#c2c2d6",    #1050
    "(MON107)3" = "#3333ff",    #272
    "(MON107)4"= "#00ccff",     #177
    "(MON107)5"= "#0099cc",     #38
    "(MON107)6"= "#9999ff",
    "(MON107)10"='#000066',
    
    "(MON135)2" = "#00cc99",    
    "(MON135)3" = "#99ff33",    
    
    "MON51-MON107" = "#ffff66",
    "MON79-MON107" = "#ffcc00",
    "MON107-MON135" = "#ff9900",
    "MON107-MON163" = "#ff9999",
    "MON51-(MON107)2" = "#ff99ff",
    "(MON79)2-MON107" = "#ff66cc",
    "MON79-(MON107)2" = "#ff3399",
    "MON79-MON107-MON135" = "#cc00cc",
    "MON51-(MON79)2-MON107" = "red"
  )

filtered_input2 <- input_file2 %>% 
    mutate(
      category = case_when(
        monomer_type == 'CEN107' ~ "CEN107",
        monomer_type == 'CEN107-like' ~ "CEN107-like",
        monomer_type == 'CEN135' ~ "CEN135",
        monomer_type == 'CEN135-like' ~ "CEN135-like",
        monomer_type == 'CEN79' ~ "CEN79",
        monomer_type == 'CEN79-like' ~ "CEN79-like",
        monomer_type == 'CEN51' ~ "CEN51",
        monomer_type == 'CEN51-like' ~ "CEN51-like",        
        monomer_type == 'other' ~ "other_monomer",
        TRUE ~ "other_monomer"
      )
    )
    
  filtered_input3 <- input_file3 %>% 
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        strand == '+' ~ "plus",
        strand == '-' ~ "minus",
        TRUE ~ "Other"
      )
    )

  filtered_input4 <- input_file4 %>% 
    filter(centype %in% c("cen_6","cen_68","cen_66","cen_103","cen_107","cen_191","cen_355","cen_383")) %>%
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category = case_when(
        centype == 'cen_6' ~ "cen_6",
        centype == 'cen_68' ~ "cen_68",      
        centype == 'cen_66' ~ "cen_66",
        centype == 'cen_103' ~ "cen_103",
        centype == 'cen_107' ~ "cen_107",
        centype == 'cen_191' ~ "cen_191",
        centype == 'cen_355' ~ "cen_355",
        centype == 'cen_383' ~ "cen_383",
        #TRUE ~ "other_satellite"
      )
    )  
  filtered_input4 <- filtered_input4 %>% 
    filter( centype %in% c("cen_6","cen_68","cen_66","cen_103","cen_107","cen_191","cen_355","cen_383") ) %>%
    mutate(
      #mer = as.numeric(mer)#factor(mer),  # Ensure mer is a factor
      category2 = case_when(
        strand == '+' ~ "plus",
        strand == '-' ~ "minus",
        TRUE ~ "Other"
      )
    )  
  
  df_plus <- filtered_input4 %>% filter(strand == '+')
  df_minus <- filtered_input4 %>% filter(strand == '-')
  
  # Create plot object
  p = ggplot()
  p =p+  geom_rect(data = df_plus, aes(xmin = chr_start,xmax = chr_end, ymin = 1,ymax =1.5 ,fill = category)) 
  p =p+  geom_rect(data = df_minus, aes(xmin = chr_start,xmax = chr_end, ymin = 0.5,ymax =1 ,fill = category)) 
  
  
  ##OTHER satellite
  p =p+  geom_rect(data = filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end, ymin = -0.7,ymax = 0,fill = category)) 
                   
  ##BLOCK strand B1
  p =p+  geom_rect(data = filtered_input3,
                   aes(xmin = block_start,xmax = block_end, ymin = -1.7,ymax = -1,fill = category)) 
  ##BLOCK strand other
  p =p+  geom_rect(data = filtered_input4,
                   aes(xmin = chr_start,xmax = chr_end, ymin = -1.7,ymax = -1,fill = category2)) 

                   
##monomer
  p =p+  geom_rect(data = filtered_input2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = -2.7,ymax = -2,fill = category)) 
                   
###General HOR
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = -3.7,ymax = -3),fill = '#6A356C',) 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = -4.7,ymax = -4),fill = '#6A356C',) 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = -5.7,ymax = -5),fill = '#6A356C',) 
  p =p+  geom_rect(data = B1,
                   aes(xmin = pos_start,xmax = pos_end, ymin = -3.7,ymax = -3,fill = category))   
  p =p+  geom_rect(data = B2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = -4.7,ymax = -4,fill = category))                    
  p =p+  geom_rect(data = B3,
                   aes(xmin = pos_start,xmax = pos_end, ymin = -5.7,ymax = -5,fill = category)) 


  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = -6.7,ymax = -6),fill = '#f2f2f2',) 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = -7.7,ymax = -7),fill = '#f2f2f2',) 
  p =p+  geom_rect(data=input_file3,
                   aes(xmin = block_start,xmax = block_end,ymin = -8.7,ymax = -8),fill = '#f2f2f2',) #6A356C
  p =p+  geom_rect(data = BB1,
                   aes(xmin = pos_start,xmax = pos_end, ymin = -6.7,ymax = -6,fill = HOR_class_new)) 
  p =p+  geom_rect(data = BB2,
                   aes(xmin = pos_start,xmax = pos_end, ymin = -7.7,ymax = -7,fill = HOR_class_new)) 
  p =p+  geom_rect(data = BB3,
                   aes(xmin = pos_start,xmax = pos_end, ymin = -8.7,ymax = -8,fill = HOR_class_new)) 
 
  

  



  
  
p =p+ geom_vline(
    xintercept = c(pos1, pos2),
    linetype = "dashed",  # Dashed line type
    color = "red",        # Line color
    size = 1              # Line thickness
  ) 
  ##
  ###################################
  ##Strandness
  
  
  p=p+  #xlim(pos1, pos2) + 
  coord_cartesian(xlim = c(pos1, pos2), expand = FALSE)+           ##xlim hides elements not included but spanning across; coord_cartesian displays everything
    ylim(-10, 2) + 
    labs(
      x = paste0(pos1, '-', pos2),
      y = "",
      fill = ""
    ) +
    theme_classic() +         
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_fill_manual(values = color_values, drop = FALSE)
  
  width_revise=max(10,(pos2-pos1)/1000000*10)
  # Save as PDF
  pdf(file = paste0(input_file_name, ".pdf"), width = width_revise / 2.54+3, height = 10 / 2.54)
  print(p)
  dev.off()
  
}
        """
        with open(f"./samples_satellite/14_chr_image/22_partall/{file_name}/plot.R",'w') as f:
             f.write(R_txt)
        os.chdir(f'./samples_satellite/14_chr_image/22_partall/{file_name}')
        mode_num=3 # Filter by HOR_repeat>=3
        subprocess.run([f'Rscript plot.R  {file_name} {start} {end} {mode_num}'], shell=True)        
        os.chdir('../../../../')

### Using monomers as units, such as 79, 107, 135
if argv1=="stepall" or "step15" in argv1:
    if  os.path.exists('./samples_satellite/15_MON_umap')==False:
        subprocess.run(["mkdir ./samples_satellite/15_MON_umap"], shell=True)
    if argv1=="stepall" or argv1=="step15_readme":
        print('Printing instructions')
        with open('./samples_satellite/15_MON_umap/readme','w') as f:
            txt=r'''
            15.0 — Generate format required for UMAP, convert sequences to numerical vectors
            15.1 — Slow, MON107 has many sequences, may take 1 day; other MON79 and MON135 may take 3 hours
            15.2 — Backfill to the original file
            15.3 — Plot UMAP 2D scatter plot
            '''
            f.write(txt)     
    if argv1=="stepall" or argv1=="step15.0":    
        print('Summarizing MON79, MON107, MON135')
        region_name_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if eachline_arr[0]=='':continue
                one_region_name=f"{eachline_arr[0]}:{eachline_arr[1]}{eachline_arr[6]}:{eachline_arr[3]}-{eachline_arr[4]}"
                region_name_list.append(one_region_name)
        print(len(region_name_list)) 
        
        MON_list=['28-28-23','28-28-28-28-23','28-28-28-23']
        for one_MON in MON_list:
            if  os.path.exists(f'./samples_satellite/15_MON_umap/{one_MON}')==False:
                subprocess.run([f"mkdir ./samples_satellite/15_MON_umap/{one_MON}"], shell=True)
            with open(f"./samples_satellite/15_MON_umap/{one_MON}/sum_MON","w") as f2:
                for one_region in region_name_list:
                    with open(f"./samples_satellite/10_monomer_model/{one_MON}/monomer/{one_region}",'r') as f:
                        for line in f:
                            eachline_arr=line.strip().split('\t')
                            seq=eachline_arr[6]
                            f2.write(f"{seq}\n")
            subprocess.run([f"sort ./samples_satellite/15_MON_umap/{one_MON}/sum_MON | uniq > ./samples_satellite/15_MON_umap/{one_MON}/sum_MON.uniq"], shell=True)    
            subprocess.run([f"rm ./samples_satellite/15_MON_umap/{one_MON}/sum_MON "], shell=True)   
            
            
            print('Generate format required for UMAP, convert sequences to numerical vectors')           
            with open(f'./samples_satellite/15_MON_umap/{one_MON}/sum_MON.vector','w') as f2:  
                with open(f"./samples_satellite/15_MON_umap/{one_MON}/sum_MON.uniq",'r') as f:
                    for line in f:
                        seq=line.strip()
                        seq_arr=seq.replace('_',"|").split('|')
                        circ_seq_len_para_str=str(len(seq)/1000)
                        addition_cols_arr=[]
                        for one in seq_arr:
                            if      one=='A':one=f"(1,0,0,0,0,0,{circ_seq_len_para_str})"
                            elif    one=='C':one=f"(0,1,0,0,0,0,{circ_seq_len_para_str})"
                            elif    one=='G':one=f"(0,0,1,0,0,0,{circ_seq_len_para_str})"
                            elif    one=='T':one=f"(0,0,0,1,0,0,{circ_seq_len_para_str})"
                            elif    one=='':one= f"(0,0,0,0,1,0,{circ_seq_len_para_str})"
                            else:           one= f"(0,0,0,0,0,1,{circ_seq_len_para_str})"   
                            addition_cols_arr.append(one)
                        addition_cols='\t'.join(addition_cols_arr) 
                        f2.write(seq+'\t'+addition_cols+'\n')          
    if argv1=="stepall" or argv1=="step15.1":
        print('step4.0 — UMAP dimensionality reduction to 1D and 2D')
        import numpy as np
        import pandas as pd
        from ast import literal_eval
        import umap
        from umap import UMAP    
        ####
        import numba
        numba.set_num_threads(70)   # UMAP uses Numba for acceleration under the hood; enable multi-threading via environment variables or Numba configuration
        #
        umap_list=  [
                        [1,15,f'./samples_satellite/15_MON_umap/28-28-23/sum_MON.vector',f'./samples_satellite/15_MON_umap/28-28-23/sum_MON.umap1v'],
                        [2,15,f'./samples_satellite/15_MON_umap/28-28-23/sum_MON.vector',f'./samples_satellite/15_MON_umap/28-28-23/sum_MON.umap2v'],
                        [1,15,f'./samples_satellite/15_MON_umap/28-28-28-28-23/sum_MON.vector',f'./samples_satellite/15_MON_umap/28-28-28-28-23/sum_MON.umap1v'],
                        [2,15,f'./samples_satellite/15_MON_umap/28-28-28-28-23/sum_MON.vector',f'./samples_satellite/15_MON_umap/28-28-28-28-23/sum_MON.umap2v'],
                        [1,15,f'./samples_satellite/15_MON_umap/28-28-28-23/sum_MON.vector',f'./samples_satellite/15_MON_umap/28-28-28-23/sum_MON.umap1v'],
                        [2,15,f'./samples_satellite/15_MON_umap/28-28-28-23/sum_MON.vector',f'./samples_satellite/15_MON_umap/28-28-28-23/sum_MON.umap2v']
                    ]
        #def run_umap(umap_list_one):
        for umap_list_one in umap_list:
            print(umap_list_one)
            umap_component_num=     umap_list_one[0]
            n_neighbors=            umap_list_one[1]
            input_file=             umap_list_one[2]
            output_file=            umap_list_one[3]
            if os.path.exists(output_file)==True:
                print(f"UMAP dimensionality reduction completed: {output_file}") 
                #return False
            ##
            print('1 — Read input file')
            print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  
            df = pd.read_csv(input_file, sep='\t', header=None)
            
            print('2 — Parse vector strings')
            print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  
            def parse_vector(vector_str):
                return np.array(literal_eval(vector_str))
            for col in df.columns[1:]:
                df[col] = df[col].apply(parse_vector)
                
            print("3 — Extract feature matrix")   
            print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  
            parsed_df = df.iloc[:, 1:].values
            features = np.array([np.concatenate(row) for row in parsed_df])
            
            # Set random seed
            print("4 — Set UMAP model") 
            print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  ### Set to 15 in umap_list       ###42,
            umap_model = UMAP(n_components=umap_component_num, n_neighbors=n_neighbors, n_jobs=-1) #n_jobs=-1 uses all cores; removed random_state=random_seed to allow parallelism
    
            try:
                print("5 — Execute UMAP dimensionality reduction")  
                print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  
                reduced_features = umap_model.fit_transform(features)
                print("Reduced features:")
                print(reduced_features)
                print('6 — Save results:')
                print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  
                if umap_component_num==1:
                    result_df = pd.DataFrame(reduced_features, columns=['umap1v_x'])
                    result_df['id'] = df.iloc[:, 0]  # Assuming the first column is the name column
                    result_df = result_df[['id', 'umap1v_x']]  # Adjust column order
                else:
                    result_df = pd.DataFrame(reduced_features, columns=['umap2v_x', 'umap2v_y'])
                    result_df['id'] = df.iloc[:, 0]  # Assuming the first column is the name column
                    result_df = result_df[['id', 'umap2v_x', 'umap2v_y']]  # Adjust column order                
                result_df.to_csv(output_file ,sep='\t',  index=False)
                print(f"UMAP dimensionality reduction completed: {output_file}")        
            except Exception as e:
                print("An error occurred:", e)
          
    
        #with Pool(processes=10) as pool:
        #    pool.map(run_umap, umap_list)          
    if argv1=="stepall" or argv1=="step15.1" or argv1=="step15.1a":
        print('Merge and generate index')
        MON_list=['28-28-23','28-28-28-28-23','28-28-28-23']
        #MON_list=['28-28-28-23']
        for one_MON in MON_list:
            dict_seq_v1={}
            with open(f'./samples_satellite/15_MON_umap/{one_MON}/sum_MON.umap1v','r')  as f:
                next(f)
                for line in f:
                    
                    eachline_arr=line.strip().split('\t')
                    dict_seq_v1[eachline_arr[0]]=eachline_arr[1]
            dict_seq_v2={}
            with open(f'./samples_satellite/15_MON_umap/{one_MON}/sum_MON.umap2v','r')  as f:
                next(f)
                for line in f:
                    
                    eachline_arr=line.strip().split('\t')
                    dict_seq_v2[eachline_arr[0]]=[eachline_arr[1],eachline_arr[2]]                
                
            with open(f"./samples_satellite/15_MON_umap/{one_MON}/sum_MON.umap1v2v","w") as f2:    
                f2.write("seq\tumap1v_x\tumap2v_x\tumap2v_y\n")
                with open(f"./samples_satellite/15_MON_umap/{one_MON}/sum_MON.uniq","r") as f:
                    for line in f:
                        seq=line.strip()
                        umap1v_x=dict_seq_v1[seq]
                        umap2v_x,umap2v_y=dict_seq_v2[seq]
                        f2.write(f"{seq}\t{umap1v_x}\t{umap2v_x}\t{umap2v_y}\n")
    if argv1=="stepall" or argv1=="step15.2":     
        print('Backfill to the original file')
        region_name_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if eachline_arr[0]=='':continue
                one_region_name=f"{eachline_arr[0]}:{eachline_arr[1]}{eachline_arr[6]}:{eachline_arr[3]}-{eachline_arr[4]}"
                region_name_list.append(one_region_name)  
      
        MON_list=['28-28-23','28-28-28-28-23','28-28-28-23']
        #MON_list=['28-28-28-23']
        for one_MON in MON_list:
            if  os.path.exists(f'./samples_satellite/15_MON_umap/{one_MON}/all')==False:
                subprocess.run([f"mkdir ./samples_satellite/15_MON_umap/{one_MON}/all"], shell=True)
            dict_seq_info={}        
            with open(f"./samples_satellite/15_MON_umap/{one_MON}/sum_MON.umap1v2v","r") as f:     
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    dict_seq_info[eachline_arr[0]]=eachline_arr[1:]
            for one_region in region_name_list:    
                with open(f"./samples_satellite/15_MON_umap/{one_MON}/all/{one_region}","w") as f2:
                    with open(f"./samples_satellite/10_monomer_model/{one_MON}/monomer/{one_region}",'r') as f:
                        for line in f:
                            eachline=line.strip()
                            eachline_arr=eachline.split('\t')
                            seq=eachline_arr[6]
                            info =dict_seq_info[seq]
                            f2.write(f"{eachline}\t{info[0]}\t{info[1]}\t{info[2]}\n")
    if argv1=="stepall" or argv1=="step15.3": 
        print('Scatter plot')
        # Load all regions
        region_list=[]
        with open('./samples_satellite/2_good_regions','r') as f:
            for line in f.readlines():
                eachline_arr=line.strip().split('\t')
                if eachline_arr[0]=='sample':continue
                if len(eachline_arr)!=8:continue
                new_name=eachline_arr[0]+':'+eachline_arr[1]+eachline_arr[6]+':'+eachline_arr[3]+'-'+eachline_arr[4]
                region_list.append(new_name)    
        region_list.sort()
        region_list_len=len(region_list)
        print("region_list_len:"+str(region_list_len))        
        ## Plot: chromosome — point — count
        MON_list=['28-28-23','28-28-28-28-23','28-28-28-23']
  
        for one_MON in MON_list:
            if  os.path.exists(f'./samples_satellite/15_MON_umap/{one_MON}/stat')==True:
                subprocess.run([f"rm -r ./samples_satellite/15_MON_umap/{one_MON}/stat"], shell=True)  
            subprocess.run([f"mkdir ./samples_satellite/15_MON_umap/{one_MON}/stat"], shell=True)     
        
        #
    
            # Search all regions
            print('Loading monomer data for all regions')
            dict_chr1MONseq_num={}
            dict_MONseq_infos={}
            i=0
            for one in region_list:
                i+=1
                #if i>10:break
                print('Progress: '+str(i)+'/'+str(region_list_len),end='\r')
                input_file=f'./samples_satellite/15_MON_umap/{one_MON}/all/'+one
                with open(input_file,'r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        region_name=eachline_arr[1];    
                        MONseq=eachline_arr[6]
                        umap2v_x=eachline_arr[12]
                        umap2v_y=eachline_arr[13]
                        dict_MONseq_infos[MONseq]=[umap2v_x,umap2v_y]
                        chr1MONseq=region_name+'_||_'+MONseq
                        if chr1MONseq not in dict_chr1MONseq_num:
                            dict_chr1MONseq_num[chr1MONseq]=0
                        dict_chr1MONseq_num[chr1MONseq]+=1 
            #print(dict_chr1MONseq_num)
            # Sort dictionary by value in descending order
            dict_chr1MONseq_num_sorted = dict(sorted(dict_chr1MONseq_num.items(), key=lambda item: item[1], reverse=True))
            with open (f'./samples_satellite/15_MON_umap/{one_MON}/stat/chr_MONseq_num_info','w') as f:
                f.write(f"chr\tMONseq\tnum\tumap2v_x\tumap2v_y\n")
                for chr1MONseq,num in dict_chr1MONseq_num_sorted.items():
                    chr1MONseq_arr=chr1MONseq.split('_||_')
                    chromosome_name     =chr1MONseq_arr[0].replace('region_','Chr')
                    MONseq             =chr1MONseq_arr[1]
                    infos               =dict_MONseq_infos[MONseq]
                    f.write(f"{chromosome_name}\t{MONseq}\t{num}\t{infos[0]}\t{infos[1]}\n")
    if argv1=="stepall" or argv1=="step15.3" or argv1=="step15.3p":
        MON_list=['28-28-23','28-28-28-28-23','28-28-28-23']

        for one_MON in MON_list:        
            print('Scatter plot, 55s')            
            R_txt=f'''library(ggplot2)
    library(dplyr)
    #install.packages("stringr") 
    library("stringr")
    
    # Set working directory
    setwd('./')
    
    # Read data
    input_file1 <- read.table('chr_MONseq_num_info', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
    input_file1_filter <- input_file1 %>% filter( num >0)  ## Or not filter
    # Sort
    input_file1_filter <- input_file1_filter %>%arrange(as.numeric(sub("Chr", "", chr)))
    # Get order and assign to factor
    chr_order <- unique(input_file1_filter$chr)
    input_file1_filter$chr <- factor(input_file1_filter$chr, levels = chr_order)

    # Create plot
    p <- ggplot()
    p <- p+geom_point(data = input_file1_filter, aes(x = umap2v_x, y = umap2v_y, size = num,color =chr)) 
    p <- p + theme_classic()
    p <- p + coord_equal() 
    # Add title and axis labels
    p <- p + theme(
      #axis.title = element_blank(),  # Hide axis titles
      #axis.text = element_blank(),   # Hide axis tick text
      #axis.ticks = element_blank(),  # Hide axis ticks
      #axis.line = element_blank()    # Hide axis lines
    )+scale_color_manual(name = "Variant Mode", 
                           values = c("Mixed" = "grey80", 
                                      "Chr1" = "red",  
                                      "Chr2" = "#8CC63F",  
                                      "Chr3" = "#CC0066",  
                                      "Chr4" = "#148EFF",  
                                      "Chr5" = "#B2DF8A",   
                                      "Chr6" = "#008000",   
                                      "Chr7" = "#FB9A99",   
                                      "Chr8" = "#00C08F",   
                                      "Chr9" = "#005FAF",  
                                      "Chr10" = "#FFCC00",  
                                      "Chr11" = "#00ccff",  
                                      "Chr12" = "#6A3D9A",  
                                      "Chr13" = "#FFFF00",  
                                      "Chr14" = "#A0522D",  
                                      "Chr15" = "#9966CC",  
                                      "Chr16" = "#ff4dff",  
                                      "Chr17" = "#ff6699", 
                                      "Chr18" = "#0073CF",  
                                      "Chr19" = "#FF8C00",  
                                      "Other" = "#CCEBC5")) 
    
    
    pdf("chr_MONseq_num_info.pdf", width = 60 / 2.54, height = 60 / 2.54)
    print(p)
    dev.off()
    
    # Create plot
    p <- ggplot()
    p <- p + geom_point(data = input_file1_filter, aes(x = umap2v_x, y = umap2v_y, size = num,color =chr))
    p <- p + facet_wrap(~ chr, ncol = 4)
    p <- p + theme_classic()
    p <- p + coord_equal() 
    # Add title and axis labels
    p <- p + theme(
      #axis.title = element_blank(),  # Hide axis titles
      #axis.text = element_blank(),   # Hide axis tick text
      #axis.ticks = element_blank(),  # Hide axis ticks
      #axis.line = element_blank()    # Hide axis lines
    )+scale_color_manual(name = "Variant Mode", 
                           values = c("Mixed" = "grey80", 
                                      "Chr1" = "red",  
                                      "Chr2" = "#8CC63F",  
                                      "Chr3" = "#CC0066",  
                                      "Chr4" = "#148EFF",  
                                      "Chr5" = "#B2DF8A",   
                                      "Chr6" = "#008000",   
                                      "Chr7" = "#FB9A99",   
                                      "Chr8" = "#00C08F",   
                                      "Chr9" = "#005FAF",  
                                      "Chr10" = "#FFCC00",  
                                      "Chr11" = "#00ccff",  
                                      "Chr12" = "#6A3D9A",                          
                                      "Chr13" = "#FFFF00",  
                                      "Chr14" = "#A0522D",  
                                      "Chr15" = "#9966CC",  
                                      "Chr16" = "#ff4dff",  
                                      "Chr17" = "#ff6699", 
                                      "Chr18" = "#0073CF",  
                                      "Chr19" = "#FF8C00",  
                                      "Other" = "#CCEBC5")) 
    pdf("chr_MONseq_num_info2.pdf", width = 100 / 2.54, height = 100 / 2.54)
    print(p)
    dev.off()
    
    '''
            with open(f'./samples_satellite/15_MON_umap/{one_MON}/stat/plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)
            new_directory = f"./samples_satellite/15_MON_umap/{one_MON}/stat/"
            os.chdir(new_directory)
            subprocess.run(['Rscript plot.R'], shell=True)    
            os.chdir('../../../../')            

### Explore MON connections, only take blocks exceeding 10000bp 
if argv1=="stepall" or "step16" in argv1: 
    if  os.path.exists('./samples_satellite/16_MONlink')==False:
        subprocess.run(["mkdir ./samples_satellite/16_MONlink"], shell=True)
    if argv1=="stepall" or argv1=="step16_readme":
        print('Printing instructions')
        with open('./samples_satellite/16_MONlink/readme','w') as f:
            txt=r'''
            16.0 —
            '''
            f.write(txt)     
    if argv1=="stepall" or argv1=="step16" or argv1=="step16.0":  
        print('step16.0_Analyze MON links, 50s')    
        if  os.path.exists('./samples_satellite/16_MONlink/0_allstat')==False:
            subprocess.run(["mkdir ./samples_satellite/16_MONlink/0_allstat"], shell=True)        
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
        ##
        
        def run_step(one_region):
            dict_chr_MON_link={}
            with open(f"./samples_satellite/11_hor/5_sum/{one_region}",'r') as f:
                A='';B=''
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    chromosome=eachline_arr[0].replace('region_','Chr')
                    MON_type=eachline_arr[20]
                    if MON_type=='delete':continue
                    if MON_type=='.': A='';B='';
                    if A=='':A=MON_type;continue
                    B=MON_type
                    if chromosome not in dict_chr_MON_link:  dict_chr_MON_link[chromosome]={}
                    if A not in dict_chr_MON_link[chromosome]:dict_chr_MON_link[chromosome][A]={}
                    if B not in dict_chr_MON_link[chromosome][A]:dict_chr_MON_link[chromosome][A][B]=0
                    dict_chr_MON_link[chromosome][A][B]+=1
                    A=B
            with open(f'./samples_satellite/16_MONlink/0_allstat/{one_region}','w') as f:
                for chromosome,dict_MON_link in dict_chr_MON_link.items():
                    for MON,dict_link_num in dict_MON_link.items():
                        for link,num in dict_link_num.items(): 
                            f.write(f"{chromosome}\t{MON}\t{link}\t{num}\n")
                        
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, region_name_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(region_name_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()                      
                
        print('Summarizing')
        dict_chr_MONlink={}
        dict_chr_sum={}
        for one_region in region_name_list: 
            with open(f'./samples_satellite/16_MONlink/0_allstat/{one_region}','r') as f:
                pos_arr=one_region.split(':')[-1].split('-')
                region_length=int(pos_arr[1])-int(pos_arr[0])
                if region_length<10000:continue
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    chromosome,MON,link,num=eachline_arr
                    if chromosome not in dict_chr_MONlink:dict_chr_MONlink[chromosome]={}
                    if chromosome not in dict_chr_sum:dict_chr_sum[chromosome]=0
                    dict_chr_sum[chromosome]+=int(num)
                    MON_link=f"{MON}___{link}"
                    if MON_link not in dict_chr_MONlink[chromosome]:dict_chr_MONlink[chromosome][MON_link]=0
                    dict_chr_MONlink[chromosome][MON_link]+=int(num)
                    
        with open(f'./samples_satellite/16_MONlink/0_sum_region10000bp','w') as f:
            for chromosome,dict_MONlink in dict_chr_MONlink.items():  
                chr_sum=dict_chr_sum[chromosome]
                for MON_link,num in dict_MONlink.items(): 
                    MON,link=MON_link.split('___')
                    f.write(f"{chromosome}\t{MON}\t{link}\t{num}\t{chr_sum}\t{round(int(num)/int(chr_sum)*100,5)}\n")
                    
                          
        print('Loading a sample information table ./samples_satellite/sample_info')
        dict_sample_type={}
        with open('./samples_satellite/sample_info') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=6:continue
                sample_name=eachline_arr[0]
                sample_type=eachline_arr[3]
                dict_sample_type[sample_name]=sample_type                         
        print('Summarizing Eurasian species')
        dict_chr_MONlink={}
        dict_chr_sum={}
        for one_region in region_name_list: 
            with open(f'./samples_satellite/16_MONlink/0_allstat/{one_region}','r') as f:
                pos_arr=one_region.split(':')[-1].split('-')
                sample=one_region.split(':')[0]
                sample_type=dict_sample_type[sample]
                if sample_type !="Eurasian":continue
                region_length=int(pos_arr[1])-int(pos_arr[0])
                if region_length<10000:continue
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    chromosome,MON,link,num=eachline_arr
                    if chromosome not in dict_chr_MONlink:dict_chr_MONlink[chromosome]={}
                    if chromosome not in dict_chr_sum:dict_chr_sum[chromosome]=0
                    dict_chr_sum[chromosome]+=int(num)
                    MON_link=f"{MON}___{link}"
                    if MON_link not in dict_chr_MONlink[chromosome]:dict_chr_MONlink[chromosome][MON_link]=0
                    dict_chr_MONlink[chromosome][MON_link]+=int(num)
        print('Summarizing East Asian species')
        dict_chr_MONlink={}
        dict_chr_sum={}
        for one_region in region_name_list: 
            with open(f'./samples_satellite/16_MONlink/0_allstat/{one_region}','r') as f:
                pos_arr=one_region.split(':')[-1].split('-')
                sample=one_region.split(':')[0]
                sample_type=dict_sample_type[sample]
                if sample_type !="East_Asia":continue
                region_length=int(pos_arr[1])-int(pos_arr[0])
                if region_length<10000:continue
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    chromosome,MON,link,num=eachline_arr
                    if chromosome not in dict_chr_MONlink:dict_chr_MONlink[chromosome]={}
                    if chromosome not in dict_chr_sum:dict_chr_sum[chromosome]=0
                    dict_chr_sum[chromosome]+=int(num)
                    MON_link=f"{MON}___{link}"
                    if MON_link not in dict_chr_MONlink[chromosome]:dict_chr_MONlink[chromosome][MON_link]=0
                    dict_chr_MONlink[chromosome][MON_link]+=int(num)
                                        
        with open(f'./samples_satellite/16_MONlink/0_sum_region10000bp_Asian','w') as f:
            for chromosome,dict_MONlink in dict_chr_MONlink.items():  
                chr_sum=dict_chr_sum[chromosome]
                for MON_link,num in dict_MONlink.items(): 
                    MON,link=MON_link.split('___')
                    f.write(f"{chromosome}\t{MON}\t{link}\t{num}\t{chr_sum}\t{round(int(num)/int(chr_sum)*100,5)}\n")            
        print('Summarizing American species')
        dict_chr_MONlink={}
        dict_chr_sum={}
        for one_region in region_name_list: 
            with open(f'./samples_satellite/16_MONlink/0_allstat/{one_region}','r') as f:
                pos_arr=one_region.split(':')[-1].split('-')
                sample=one_region.split(':')[0]
                sample_type=dict_sample_type[sample]
                if sample_type !="America":continue
                region_length=int(pos_arr[1])-int(pos_arr[0])
                if region_length<10000:continue
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    chromosome,MON,link,num=eachline_arr
                    if chromosome not in dict_chr_MONlink:dict_chr_MONlink[chromosome]={}
                    if chromosome not in dict_chr_sum:dict_chr_sum[chromosome]=0
                    dict_chr_sum[chromosome]+=int(num)
                    MON_link=f"{MON}___{link}"
                    if MON_link not in dict_chr_MONlink[chromosome]:dict_chr_MONlink[chromosome][MON_link]=0
                    dict_chr_MONlink[chromosome][MON_link]+=int(num)
                                        
        with open(f'./samples_satellite/16_MONlink/0_sum_region10000bp_American','w') as f:
            for chromosome,dict_MONlink in dict_chr_MONlink.items():  
                chr_sum=dict_chr_sum[chromosome]
                for MON_link,num in dict_MONlink.items(): 
                    MON,link=MON_link.split('___')
                    f.write(f"{chromosome}\t{MON}\t{link}\t{num}\t{chr_sum}\t{round(int(num)/int(chr_sum)*100,5)}\n")                      
    if argv1=="stepall" or argv1=="step16" or argv1=="step16.1": 
        print('Loading file')
        dict_chr_result_list={} 
        with open (f'./samples_satellite/16_MONlink/0_sum_region10000bp','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                chromosome,MON,link,num,chrall,percent=eachline_arr
                if chromosome not in dict_chr_result_list: dict_chr_result_list[chromosome]=[]
                dict_chr_result_list[chromosome].append([MON,link,int(num),chrall,percent])
        
        print('Output') 
        dict_chr_pair_num={}
        #
        with open (f'./samples_satellite/16_MONlink/1_sum_filter','w') as f2:
            f2.write("Chromosome\tFrom\tTo\tNum\tChrAll\tPercentage\n")
            for chromosome,result_list in dict_chr_result_list.items():
                k=0
                sorted_result_list = sorted(result_list, key=lambda x: x[2], reverse=True)
                #print(sorted_result_list)
                for one_result in sorted_result_list:
                    MON,link,num,chrall,percent=one_result
                    k+=1
                    #if k>10:break
                    LIST=['28-23','28-28-23','28-28-28-23','28-28-28-28-23','28-28-28-28-28-23','28-28-28-28-28-28-23']
                    if MON not in LIST or link not in LIST:continue
                    f2.write(f"{chromosome}\t{MON}\t{link}\t{num}\t{chrall}\t{percent}\n")
                    ########
                    if len(link)>len(MON):  pair=f"{MON}___{link}"
                    else:                   pair=f"{link}___{MON}"
                    if chromosome not in dict_chr_pair_num:dict_chr_pair_num[chromosome]={}
                    if pair not in dict_chr_pair_num[chromosome]:
                        dict_chr_pair_num[chromosome][pair]={}
                        dict_chr_pair_num[chromosome][pair]['num']=0
                        dict_chr_pair_num[chromosome][pair]['percent']=0
                        dict_chr_pair_num[chromosome][pair]['chrall']=0
                    dict_chr_pair_num[chromosome][pair]['num']+=int(num)
                    dict_chr_pair_num[chromosome][pair]['percent']+=float(percent)
                    dict_chr_pair_num[chromosome][pair]['chrall']=chrall
        with open (f'./samples_satellite/16_MONlink/1_sum_filter_pair','w') as f2:
            f2.write("Chromosome\tpair\tNum\tChrAll\tPercentage\n")
            for chromosome,dict_pair_info in dict_chr_pair_num.items():
                for pair,infos in dict_pair_info.items():
                    num=infos['num']
                    percent=infos['percent']
                    chrall=infos['chrall']
                    f2.write(f"{chromosome}\t{pair}\t{num}\t{chrall}\t{percent}\n")
    if argv1=="stepall" or argv1=="step16" or argv1=="step16.1a": 
        print('Loading file')
        dict_chr_result_list={} 
        with open (f'./samples_satellite/16_MONlink/0_sum_region10000bp_Eurasian','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                chromosome,MON,link,num,chrall,percent=eachline_arr
                if chromosome not in dict_chr_result_list: dict_chr_result_list[chromosome]=[]
                dict_chr_result_list[chromosome].append([MON,link,int(num),chrall,percent])

        print('Output') 
        dict_chr_pair_num={}
        #
        with open (f'./samples_satellite/16_MONlink/1_sum_filter_Eurasian','w') as f2:
            f2.write("Chromosome\tFrom\tTo\tNum\tChrAll\tPercentage\n")
            for chromosome,result_list in dict_chr_result_list.items():
                k=0
                sorted_result_list = sorted(result_list, key=lambda x: x[2], reverse=True)
                #print(sorted_result_list)
                for one_result in sorted_result_list:
                    MON,link,num,chrall,percent=one_result
                    k+=1
                    #if k>10:break
                    LIST=['28-23','28-28-23','28-28-28-23','28-28-28-28-23','28-28-28-28-28-23','28-28-28-28-28-28-23']
                    if MON not in LIST or link not in LIST:continue
                    f2.write(f"{chromosome}\t{MON}\t{link}\t{num}\t{chrall}\t{percent}\n")
                    ########
                    if len(link)>len(MON):  pair=f"{MON}___{link}"
                    else:                   pair=f"{link}___{MON}"
                    if chromosome not in dict_chr_pair_num:dict_chr_pair_num[chromosome]={}
                    if pair not in dict_chr_pair_num[chromosome]:
                        dict_chr_pair_num[chromosome][pair]={}
                        dict_chr_pair_num[chromosome][pair]['num']=0
                        dict_chr_pair_num[chromosome][pair]['percent']=0
                        dict_chr_pair_num[chromosome][pair]['chrall']=0
                    dict_chr_pair_num[chromosome][pair]['num']+=int(num)
                    dict_chr_pair_num[chromosome][pair]['percent']+=float(percent)
                    dict_chr_pair_num[chromosome][pair]['chrall']=chrall
        with open (f'./samples_satellite/16_MONlink/1_sum_filter_pair_Eurasian','w') as f2:
            f2.write("Chromosome\tpair\tNum\tChrAll\tPercentage\n")
            for chromosome,dict_pair_info in dict_chr_pair_num.items():
                for pair,infos in dict_pair_info.items():
                    num=infos['num']
                    percent=infos['percent']
                    chrall=infos['chrall']
                    f2.write(f"{chromosome}\t{pair}\t{num}\t{chrall}\t{percent}\n")                    
    if argv1=="stepall" or argv1=="step16" or argv1=="step16.1b": 
        print('Loading file')
        dict_chr_result_list={} 
        with open (f'./samples_satellite/16_MONlink/0_sum_region10000bp_Asian','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                chromosome,MON,link,num,chrall,percent=eachline_arr
                if chromosome not in dict_chr_result_list: dict_chr_result_list[chromosome]=[]
                dict_chr_result_list[chromosome].append([MON,link,int(num),chrall,percent])

        print('Output') 
        dict_chr_pair_num={}
        #
        with open (f'./samples_satellite/16_MONlink/1_sum_filter_Asian','w') as f2:
            f2.write("Chromosome\tFrom\tTo\tNum\tChrAll\tPercentage\n")
            for chromosome,result_list in dict_chr_result_list.items():
                k=0
                sorted_result_list = sorted(result_list, key=lambda x: x[2], reverse=True)
                #print(sorted_result_list)
                for one_result in sorted_result_list:
                    MON,link,num,chrall,percent=one_result
                    k+=1
                    #if k>10:break
                    LIST=['28-23','28-28-23','28-28-28-23','28-28-28-28-23','28-28-28-28-28-23','28-28-28-28-28-28-23']
                    if MON not in LIST or link not in LIST:continue
                    f2.write(f"{chromosome}\t{MON}\t{link}\t{num}\t{chrall}\t{percent}\n")
                    ########
                    if len(link)>len(MON):  pair=f"{MON}___{link}"
                    else:                   pair=f"{link}___{MON}"
                    if chromosome not in dict_chr_pair_num:dict_chr_pair_num[chromosome]={}
                    if pair not in dict_chr_pair_num[chromosome]:
                        dict_chr_pair_num[chromosome][pair]={}
                        dict_chr_pair_num[chromosome][pair]['num']=0
                        dict_chr_pair_num[chromosome][pair]['percent']=0
                        dict_chr_pair_num[chromosome][pair]['chrall']=0
                    dict_chr_pair_num[chromosome][pair]['num']+=int(num)
                    dict_chr_pair_num[chromosome][pair]['percent']+=float(percent)
                    dict_chr_pair_num[chromosome][pair]['chrall']=chrall
        with open (f'./samples_satellite/16_MONlink/1_sum_filter_pair_Asian','w') as f2:
            f2.write("Chromosome\tpair\tNum\tChrAll\tPercentage\n")
            for chromosome,dict_pair_info in dict_chr_pair_num.items():
                for pair,infos in dict_pair_info.items():
                    num=infos['num']
                    percent=infos['percent']
                    chrall=infos['chrall']
                    f2.write(f"{chromosome}\t{pair}\t{num}\t{chrall}\t{percent}\n")          
    if argv1=="stepall" or argv1=="step16" or argv1=="step16.1c": 
        print('Loading file')
        dict_chr_result_list={} 
        with open (f'./samples_satellite/16_MONlink/0_sum_region10000bp_American','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                chromosome,MON,link,num,chrall,percent=eachline_arr
                if chromosome not in dict_chr_result_list: dict_chr_result_list[chromosome]=[]
                dict_chr_result_list[chromosome].append([MON,link,int(num),chrall,percent])

        print('Output') 
        dict_chr_pair_num={}
        #
        with open (f'./samples_satellite/16_MONlink/1_sum_filter_American','w') as f2:
            f2.write("Chromosome\tFrom\tTo\tNum\tChrAll\tPercentage\n")
            for chromosome,result_list in dict_chr_result_list.items():
                k=0
                sorted_result_list = sorted(result_list, key=lambda x: x[2], reverse=True)
                #print(sorted_result_list)
                for one_result in sorted_result_list:
                    MON,link,num,chrall,percent=one_result
                    k+=1
                    #if k>10:break
                    LIST=['28-23','28-28-23','28-28-28-23','28-28-28-28-23','28-28-28-28-28-23','28-28-28-28-28-28-23']
                    if MON not in LIST or link not in LIST:continue
                    f2.write(f"{chromosome}\t{MON}\t{link}\t{num}\t{chrall}\t{percent}\n")
                    ########
                    if len(link)>len(MON):  pair=f"{MON}___{link}"
                    else:                   pair=f"{link}___{MON}"
                    if chromosome not in dict_chr_pair_num:dict_chr_pair_num[chromosome]={}
                    if pair not in dict_chr_pair_num[chromosome]:
                        dict_chr_pair_num[chromosome][pair]={}
                        dict_chr_pair_num[chromosome][pair]['num']=0
                        dict_chr_pair_num[chromosome][pair]['percent']=0
                        dict_chr_pair_num[chromosome][pair]['chrall']=0
                    dict_chr_pair_num[chromosome][pair]['num']+=int(num)
                    dict_chr_pair_num[chromosome][pair]['percent']+=float(percent)
                    dict_chr_pair_num[chromosome][pair]['chrall']=chrall
        with open (f'./samples_satellite/16_MONlink/1_sum_filter_pair_American','w') as f2:
            f2.write("Chromosome\tpair\tNum\tChrAll\tPercentage\n")
            for chromosome,dict_pair_info in dict_chr_pair_num.items():
                for pair,infos in dict_pair_info.items():
                    num=infos['num']
                    percent=infos['percent']
                    chrall=infos['chrall']
                    f2.write(f"{chromosome}\t{pair}\t{num}\t{chrall}\t{percent}\n")          
                   
### Explore MON variation, only take blocks exceeding 10000bp
if argv1=="stepall" or "step17" in argv1:
    if  os.path.exists('./samples_satellite/17_MONvar')==False:
        subprocess.run(["mkdir ./samples_satellite/17_MONvar"], shell=True)
    if argv1=="stepall" or argv1=="step17_readme":
        print('Printing instructions')
        with open('./samples_satellite/17_MONvar/readme','w') as f:
            txt=r'''
            17.0 — Summarize all 28-28-28-23, requiring block > 10000bp, precise position is 107bp
            17.1 — Summarize all 28-28-28-28-23, requiring block > 10000bp, precise position is 135bp
            17.2 — Summarize all 28-28-23, requiring block > 10000bp, precise position is 79bp
            
            b — Calculate frequencies
            c — Output images
            '''
            f.write(txt)         
    if argv1=="stepall" or argv1=="step17" or argv1=="step17.0":  
        print('step17.0_Statistics for CEN107')        
        if  os.path.exists('./samples_satellite/17_MONvar/0_allstat')==False:
            subprocess.run(["mkdir ./samples_satellite/17_MONvar/0_allstat"], shell=True)     
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
        ##
        
        def run_step(one_region):
            pos_arr=one_region.split(':')[-1].split('-')
            region_length=int(pos_arr[1])-int(pos_arr[0])
            if region_length<10000:
                with open(f'./samples_satellite/17_MONvar/0_allstat/{one_region}','w') as f:
                    f.write('')
                return False
            dict_chr_MON_link={}
            with open(f"./samples_satellite/11_hor/5_sum/{one_region}",'r') as f:
                next(f)
                MON_tmp_list=[]
                MON_seq_list=[]
                mark_serilal_tmp=''
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    MON_type=eachline_arr[20]
                    circ_seq=eachline_arr[5]
                    if MON_tmp_list!=[] and eachline_arr[21]!=mark_serilal_tmp and mark_serilal_tmp!='':
                        MON_seq='_'.join(MON_tmp_list)
                        MON_seq_list.append(MON_seq)
                    if MON_type=='28-28-28-23':
                        mark_serilal_tmp=eachline_arr[21]
                        MON_tmp_list=[circ_seq]
                        continue
                    if eachline_arr[21]==mark_serilal_tmp:
                        MON_tmp_list.append(circ_seq)
                        
            with open(f'./samples_satellite/17_MONvar/0_allstat/{one_region}','w') as f:
                for one in MON_seq_list:
                    f.write(f"{one}\n")
                        
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, region_name_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(region_name_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()                              
        i=0
        with open(f'./samples_satellite/17_MONvar/0_all_seq_MON107','w') as f2:
            for one_region in region_name_list: 
                i+=1
                print(i,end='\r')
                with open(f'./samples_satellite/17_MONvar/0_allstat/{one_region}','r') as f:
                    for line in f:
                        seq=line.strip()
                        seq_pure=seq.replace("||||||||||||",'').replace('||','').replace('_','')
                        if len(seq_pure)!=107:continue   #'28-28-28-23'
                        if len(seq_pure.replace('|',''))!=107:continue
                        f2.write(seq_pure+'\n')
    if argv1=="stepall" or argv1=="step17" or argv1=="step17.1":  
        print('step17.1_Statistics for CEN135')        
        if  os.path.exists('./samples_satellite/17_MONvar/1_allstat')==False:
            subprocess.run(["mkdir ./samples_satellite/17_MONvar/1_allstat"], shell=True)     
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
        ##
        
        def run_step(one_region):
            pos_arr=one_region.split(':')[-1].split('-')
            region_length=int(pos_arr[1])-int(pos_arr[0])
            if region_length<10000:
                with open(f'./samples_satellite/17_MONvar/1_allstat/{one_region}','w') as f:
                    f.write('')
                return False
            dict_chr_MON_link={}
            with open(f"./samples_satellite/11_hor/5_sum/{one_region}",'r') as f:
                next(f)
                MON_tmp_list=[]
                MON_seq_list=[]
                mark_serilal_tmp=''
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    MON_type=eachline_arr[20]
                    circ_seq=eachline_arr[5]
                    if MON_tmp_list!=[] and eachline_arr[21]!=mark_serilal_tmp and mark_serilal_tmp!='':
                        MON_seq='_'.join(MON_tmp_list)
                        MON_seq_list.append(MON_seq)
                    if MON_type=='28-28-28-28-23':
                        mark_serilal_tmp=eachline_arr[21]
                        MON_tmp_list=[circ_seq]
                        continue
                    if eachline_arr[21]==mark_serilal_tmp:
                        MON_tmp_list.append(circ_seq)
                        
            with open(f'./samples_satellite/17_MONvar/1_allstat/{one_region}','w') as f:
                for one in MON_seq_list:
                    f.write(f"{one}\n")
                        
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, region_name_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(region_name_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()                              
        i=0
        with open(f'./samples_satellite/17_MONvar/1_all_seq_MON135','w') as f2:
            for one_region in region_name_list: 
                i+=1
                print(i,end='\r')
                with open(f'./samples_satellite/17_MONvar/1_allstat/{one_region}','r') as f:
                    for line in f:
                        seq=line.strip()
                        seq_pure=seq.replace("||||||||||||",'').replace('||','').replace('_','')
                        seq_pure_len=len(seq_pure)
                        #print(seq_pure_len)
                        if seq_pure_len!=135:continue   #'28-28-28-23'
                        if len(seq_pure.replace('|',''))!=135:continue
                        f2.write(seq_pure+'\n')                    
    if argv1=="stepall" or argv1=="step17" or argv1=="step17.2":  
        print('step17.2_Statistics for CEN79')        
        if  os.path.exists('./samples_satellite/17_MONvar/2_allstat')==False:
            subprocess.run(["mkdir ./samples_satellite/17_MONvar/2_allstat"], shell=True)     
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
        ##
        
        def run_step(one_region):
            pos_arr=one_region.split(':')[-1].split('-')
            region_length=int(pos_arr[1])-int(pos_arr[0])
            if region_length<10000:
                with open(f'./samples_satellite/17_MONvar/2_allstat/{one_region}','w') as f:
                    f.write('')
                return False
            dict_chr_MON_link={}
            with open(f"./samples_satellite/11_hor/5_sum/{one_region}",'r') as f:
                next(f)
                MON_tmp_list=[]
                MON_seq_list=[]
                mark_serilal_tmp=''
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    MON_type=eachline_arr[20]
                    circ_seq=eachline_arr[5]
                    if MON_tmp_list!=[] and eachline_arr[21]!=mark_serilal_tmp and mark_serilal_tmp!='':
                        MON_seq='_'.join(MON_tmp_list)
                        MON_seq_list.append(MON_seq)
                    if MON_type=='28-28-23':
                        mark_serilal_tmp=eachline_arr[21]
                        MON_tmp_list=[circ_seq]
                        continue
                    if eachline_arr[21]==mark_serilal_tmp:
                        MON_tmp_list.append(circ_seq)
                        
            with open(f'./samples_satellite/17_MONvar/2_allstat/{one_region}','w') as f:
                for one in MON_seq_list:
                    f.write(f"{one}\n")
                        
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, region_name_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(region_name_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()                              
        i=0
        with open(f'./samples_satellite/17_MONvar/2_all_seq_MON79','w') as f2:
            for one_region in region_name_list: 
                i+=1
                print(i,end='\r')
                with open(f'./samples_satellite/17_MONvar/2_allstat/{one_region}','r') as f:
                    for line in f:
                        seq=line.strip()
                        seq_pure=seq.replace("||||||||||||",'').replace('||','').replace('_','')
                        seq_pure_len=len(seq_pure)
                        #print(seq_pure_len)
                        if seq_pure_len!=79:continue   #'28-28-28-23'
                        if len(seq_pure.replace('|',''))!=79:continue
                        f2.write(seq_pure+'\n')                                        
    #
    if argv1=="stepall" or argv1=="step17" or argv1=="step17.0b":  
        print('Extracting data')              
        # Initialize statistical array: 107 positions, 4 bases each (A/C/G/T)
        counts = [[0] * 4 for _ in range(107)]  # 107x4 nested list
        base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}  # Base to index mapping
         
        with open('./samples_satellite/17_MONvar/0_all_seq_MON107', 'r') as f:
            # Read all data and remove newline characters
            data = f.read().replace('\n', '')  # Concatenate directly into a single line string
            total_length = len(data)
            if total_length % 107 != 0:
                print(f"Warning: Total length {total_length} is not a multiple of 107, data may be incomplete")
         
            # Count by column: the base sequence at each position is data[pos::107]
            for pos in range(107):
                print(f"Processing position {pos + 1}/107 ...", end='\r')  # \r moves cursor back to the beginning of the line  
                column = data[pos::107]  # Extract the same position from all rows (step size 107)
                counts[pos][0] = column.count('A')  # Count A
                counts[pos][1] = column.count('C')  # Count C
                counts[pos][2] = column.count('G')  # Count G
                counts[pos][3] = column.count('T')  # Count T
         
        # Output results (same as before)
        with open('./samples_satellite/17_MONvar/0_all_seq_MON107_stat', 'w') as f:
            f.write('POS\tA\tC\tG\tT\tsum\tA_percent\tC_percent\tG_percent\tT_percent\tMain\tMain_num\tMain_percent\n')
            for pos in range(107):
                a, c, g, t = counts[pos]
                total = a + c + g + t
                a_p = round(a / total, 3) if total > 0 else 0
                c_p = round(c / total, 3) if total > 0 else 0
                g_p = round(g / total, 3) if total > 0 else 0
                t_p = round(t / total, 3) if total > 0 else 0
                main_idx = max(range(4), key=lambda i: counts[pos][i])
                main_base = ['A', 'C', 'G', 'T'][main_idx]
                main_num = counts[pos][main_idx]
                main_percent = round(main_num / total, 3) if total > 0 else 0
                f.write(f"{pos+1}\t{a}\t{c}\t{g}\t{t}\t{total}\t{a_p}\t{c_p}\t{g_p}\t{t_p}\t{main_base}\t{main_num}\t{main_percent}\n")    
    if argv1=="stepall" or argv1=="step17" or argv1=="step17.1b":  
        print('Extracting data')              
        # Initialize statistical array: 135 positions, 4 bases each (A/C/G/T)
        counts = [[0] * 4 for _ in range(135)]  # 135x4 nested list
        base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}  # Base to index mapping
         
        with open('./samples_satellite/17_MONvar/1_all_seq_MON135', 'r') as f:
            # Read all data and remove newline characters
            data = f.read().replace('\n', '')  # Concatenate directly into a single line string
            total_length = len(data)
            if total_length % 135 != 0:
                print(f"Warning: Total length {total_length} is not a multiple of 135, data may be incomplete")
         
            # Count by column: the base sequence at each position is data[pos::135]
            for pos in range(135):
                print(f"Processing position {pos + 1}/135 ...", end='\r')  # \r moves cursor back to the beginning of the line  
                column = data[pos::135]  # Extract the same position from all rows (step size 135)
                counts[pos][0] = column.count('A')  # Count A
                counts[pos][1] = column.count('C')  # Count C
                counts[pos][2] = column.count('G')  # Count G
                counts[pos][3] = column.count('T')  # Count T
         
        # Output results (same as before)
        with open('./samples_satellite/17_MONvar/1_all_seq_MON135_stat', 'w') as f:
            f.write('POS\tA\tC\tG\tT\tsum\tA_percent\tC_percent\tG_percent\tT_percent\tMain\tMain_num\tMain_percent\n')
            for pos in range(135):
                a, c, g, t = counts[pos]
                total = a + c + g + t
                a_p = round(a / total, 3) if total > 0 else 0
                c_p = round(c / total, 3) if total > 0 else 0
                g_p = round(g / total, 3) if total > 0 else 0
                t_p = round(t / total, 3) if total > 0 else 0
                main_idx = max(range(4), key=lambda i: counts[pos][i])
                main_base = ['A', 'C', 'G', 'T'][main_idx]
                main_num = counts[pos][main_idx]
                main_percent = round(main_num / total, 3) if total > 0 else 0
                f.write(f"{pos+1}\t{a}\t{c}\t{g}\t{t}\t{total}\t{a_p}\t{c_p}\t{g_p}\t{t_p}\t{main_base}\t{main_num}\t{main_percent}\n")
    if argv1=="stepall" or argv1=="step17" or argv1=="step17.2b":  
        print('Extracting data')              
        # Initialize statistical array: 79 positions, 4 bases each (A/C/G/T)
        counts = [[0] * 4 for _ in range(79)]  # 79x4 nested list
        base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}  # Base to index mapping
         
        with open('./samples_satellite/17_MONvar/2_all_seq_MON79', 'r') as f:
            # Read all data and remove newline characters
            data = f.read().replace('\n', '')  # Concatenate directly into a single line string
            total_length = len(data)
            if total_length % 79 != 0:
                print(f"Warning: Total length {total_length} is not a multiple of 79, data may be incomplete")
         
            # Count by column: the base sequence at each position is data[pos::79]
            for pos in range(79):
                print(f"Processing position {pos + 1}/79 ...", end='\r')  # \r moves cursor back to the beginning of the line  
                column = data[pos::79]  # Extract the same position from all rows (step size 79)
                counts[pos][0] = column.count('A')  # Count A
                counts[pos][1] = column.count('C')  # Count C
                counts[pos][2] = column.count('G')  # Count G
                counts[pos][3] = column.count('T')  # Count T
         
        # Output results (same as before)
        with open('./samples_satellite/17_MONvar/2_all_seq_MON79_stat', 'w') as f:
            f.write('POS\tA\tC\tG\tT\tsum\tA_percent\tC_percent\tG_percent\tT_percent\tMain\tMain_num\tMain_percent\n')
            for pos in range(79):
                a, c, g, t = counts[pos]
                total = a + c + g + t
                a_p = round(a / total, 3) if total > 0 else 0
                c_p = round(c / total, 3) if total > 0 else 0
                g_p = round(g / total, 3) if total > 0 else 0
                t_p = round(t / total, 3) if total > 0 else 0
                main_idx = max(range(4), key=lambda i: counts[pos][i])
                main_base = ['A', 'C', 'G', 'T'][main_idx]
                main_num = counts[pos][main_idx]
                main_percent = round(main_num / total, 3) if total > 0 else 0
                f.write(f"{pos+1}\t{a}\t{c}\t{g}\t{t}\t{total}\t{a_p}\t{c_p}\t{g_p}\t{t_p}\t{main_base}\t{main_num}\t{main_percent}\n")
                
    if argv1=="stepall" or argv1=="step17" or argv1=="step17.0c":  
        R_txt=f'''library(ggplot2)
# Set working directory
setwd('./')

# Read data
input_file1 <- read.table('0_all_seq_MON107_stat', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')

# Create plot
p <- ggplot()
p <- p+geom_col(data = input_file1, aes(x = POS, y = 1-Main_percent,fill=Main),color='black',size=0.1) 
p <- p + theme_classic()
# Add title and axis labels
p <- p + theme(
  axis.title = element_blank(),  # Hide axis titles
  axis.text = element_blank(),   # Hide axis tick text
  legend.position = "none",
  axis.ticks = element_blank(),  # Hide axis ticks
  axis.line = element_blank()    # Hide axis lines
)+scale_fill_manual(values = c("A" = "#009e73", 
                                      "C" = "#56b4e9",  
                                      "G" = "#e69f00",  
                                      "T" = "#cc79a7")) 


pdf("0_all_seq_MON107_stat.col.pdf", width = 10.7 / 2.54, height = 6 / 2.54)
print(p)
dev.off()
'''
        with open(f'./samples_satellite/17_MONvar/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = f"./samples_satellite/17_MONvar/"
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')
    if argv1=="stepall" or argv1=="step17" or argv1=="step17.1c":  
        R_txt=f'''library(ggplot2)
# Set working directory
setwd('./')

# Read data
input_file1 <- read.table('1_all_seq_MON135_stat', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')

# Create plot
p <- ggplot()
p <- p+geom_col(data = input_file1, aes(x = POS, y = 1-Main_percent,fill=Main),color='black',size=0.1) 
p <- p + theme_classic()
# Add title and axis labels
p <- p + theme(
  axis.title = element_blank(),  # Hide axis titles
  axis.text = element_blank(),   # Hide axis tick text
  legend.position = "none",
  axis.ticks = element_blank(),  # Hide axis ticks
  axis.line = element_blank()    # Hide axis lines
)+scale_fill_manual(values = c("A" = "#009e73", 
                                      "C" = "#56b4e9",  
                                      "G" = "#e69f00",  
                                      "T" = "#cc79a7")) 


pdf("1_all_seq_MON135_stat.col.pdf", width = 13.5 / 2.54, height = 6 / 2.54)
print(p)
dev.off()
'''
        with open(f'./samples_satellite/17_MONvar/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = f"./samples_satellite/17_MONvar/"
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')            
    if argv1=="stepall" or argv1=="step17" or argv1=="step17.2c":  
        R_txt=f'''library(ggplot2)
# Set working directory
setwd('./')

# Read data
input_file1 <- read.table('2_all_seq_MON79_stat', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')

# Create plot
p <- ggplot()
p <- p+geom_col(data = input_file1, aes(x = POS, y = 1-Main_percent,fill=Main),color='black',size=0.1) 
p <- p + theme_classic()
# Add title and axis labels
p <- p + theme(
  axis.title = element_blank(),  # Hide axis titles
  axis.text = element_blank(),   # Hide axis tick text
  legend.position = "none",
  axis.ticks = element_blank(),  # Hide axis ticks
  axis.line = element_blank()    # Hide axis lines
)+scale_fill_manual(values = c("A" = "#009e73", 
                                      "C" = "#56b4e9",  
                                      "G" = "#e69f00",  
                                      "T" = "#cc79a7")) 


pdf("2_all_seq_MON79_stat.col.pdf", width = 7.9 / 2.54, height = 6 / 2.54)
print(p)
dev.off()
'''
        with open(f'./samples_satellite/17_MONvar/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = f"./samples_satellite/17_MONvar/"
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')
        
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))      
time_end=timeit.default_timer()
print('All the running time: %.0f Seconds'%(time_end-time_start))






















