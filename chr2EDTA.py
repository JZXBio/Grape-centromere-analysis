#!/usr/bin/python
# -*- coding: UTF-8 -*-
import sys
argvs=sys.argv
argv_len=len(sys.argv)
if not argvs[1].startswith("part"):
    if argv_len < 2 or argvs[1] in ("h", "-h", "help", "-help") or not argvs[1].startswith("step"):
        print ("chr2EDTA.py-----help:")
        print ("")
        print ("Usage：")
        print ("chr2EDTA.py step0 -i  Specify a genome folder    ## Very fast, redo takes 1min")
        print ("chr2EDTA.py step1   Run EDTA")
        print ("chr2EDTA.py part2   Run TEsorter")
        print ("chr2EDTA.py step2_TEsorter_plant   Identify additional sequences")
        
        print ("chr2EDTA.py step3   EDTA,,, statistics")
        print ("chr2EDTA.py part4   TEsorter,,, statistics")     
        
        print ("chr2EDTA.py part5   Build phylogenetic tree")   
        txt=r'''
            A total of over 60,000 Athila sequences were detected,
            99.9% are GAG|Athila PROT|Athila RT|Athila RH|Athila INT|Athila
            5.0 - Get info
            5.0s - Add additional columns: whether in VSat1, peri, and interarray
            ##
            step5.0s2a — Use EDTA's judgment to determine LTR sequences, calculate similarity and infer insertion time
            step5.0s3 — JZX manual analysis of LTR regions
            
            5.1_cdhit deduplication
            5.1s_ Process to get cdhit.index
            5.2-mafft167s
            5.3-trimAl
            5.3s_ Format fasta
            5.4-iqtree
            5.5 - Generate an empty file, output the iTOL tree in Newick format, copy into this file
            5.6 - Generate serial file
            5.7 - Following step5.0s info, add additional columns and sort according to phylogenetic tree
            5.8 - Plot
            
            step5.10 — Select a portion to redo phylogenetic tree ./chr2EDTA/5_tree/{TE_Class}/10_ClassII/
            '''
        print(txt)
        
        print ("chr2EDTA.py part6   Analyze whether a certain position is associated with a trait")    
        print("s")
        
        print ("-thread \t\tNumber of threads (default 10), some steps use multiprocessing")
        print ("-i      \t\tRequired for step0, input fasta file")
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


EDTA="EDTA"
Interval='S'*100
terminal='X' 

    
    
if  os.path.exists('./chr2EDTA')==False:
    subprocess.run(["mkdir chr2EDTA"], shell=True)

#step0
if argv1=="stepall" or argv1=="step0":
    if  os.path.exists('./chr2EDTA/0_prepare')==True:
        subprocess.run(["rm -r ./chr2EDTA/0_prepare"], shell=True)  
    subprocess.run(["mkdir ./chr2EDTA/0_prepare"], shell=True)  
    
    if  os.path.exists('./chr2EDTA/error')==True:subprocess.run(["rm ./chr2EDTA/error"], shell=True)  
    
    
    
    if "i"  in args_dict: #print("缺少输入的fasta文件");sys.exit()
        input_dictionary = args_dict["i"]
        with open('./chr2EDTA/0_prepare/sample_source','w') as f:
            f.write(input_dictionary)
    else:
        print("缺少输入的fasta文件");sys.exit()
    
    dir_fasta_file_num=0
    dir_file_name_list=[]
    files=os.listdir(input_dictionary)
    with open('./chr2EDTA/0_prepare/sample_list','w') as f:
        for one in files:
            if     one.endswith('.fasta'):one_name=one[:-6]
            elif   one.endswith('.fa'):one_name=one[:-3]
            elif   one.endswith('.fna'):one_name=one[:-4] 
            else :continue
            dir_fasta_file_num+=1
            dir_file_name_list.append([input_dictionary,one,one_name])
            f.write(one_name+'\n')

    print('基因组文件数量：'+str(dir_fasta_file_num))

    def run_step0(dir_file):
        input_dir=dir_file[0]
        file_name=dir_file[1]
        sample_name=dir_file[2]
        input_file=input_dir+'/'+file_name
        output_dir="./chr2EDTA/0_prepare/"+sample_name+"/"
        #if  os.path.exists(output_dir)==True:return False
        output_file=output_dir+sample_name
        subprocess.run(["mkdir "+output_dir], shell=True)  
        ###取长度1000000限制
        chr_num=0
        with open(f'{output_dir}/{sample_name}.chrname_old2new','w') as f3:
            with open(f'{output_dir}/{sample_name}.fasta','w') as f2:
                one_seq=''
                with open(input_file,'r') as f:
                    for line in f.readlines():
                        eachline_arr=line.strip().split(' ')
                        eachline=eachline_arr[0]
                        if eachline[0]=='>':
                            ##
                            if len(one_seq)>1000000:
                                one_id_old=one_id
                                one_id_lower=one_id.lower()
                                good_mark=''
                                if one_id.isdigit():
                                    one_id="Chr"+one_id
                                elif one_id_lower.endswith('chr01') or one_id_lower.endswith('chr1'):   one_id="Chr1";      good_mark='yes'
                                elif one_id_lower.endswith('chr02') or one_id_lower.endswith('chr2'):   one_id="Chr2";      good_mark='yes'
                                elif one_id_lower.endswith('chr03') or one_id_lower.endswith('chr3'):   one_id="Chr3";      good_mark='yes'
                                elif one_id_lower.endswith('chr04') or one_id_lower.endswith('chr4'):   one_id="Chr4";      good_mark='yes'
                                elif one_id_lower.endswith('chr05') or one_id_lower.endswith('chr5'):   one_id="Chr5";      good_mark='yes'
                                elif one_id_lower.endswith('chr06') or one_id_lower.endswith('chr6'):   one_id="Chr6";      good_mark='yes'
                                elif one_id_lower.endswith('chr07') or one_id_lower.endswith('chr7'):   one_id="Chr7";      good_mark='yes'
                                elif one_id_lower.endswith('chr08') or one_id_lower.endswith('chr8'):   one_id="Chr8";      good_mark='yes'
                                elif one_id_lower.endswith('chr09') or one_id_lower.endswith('chr9'):   one_id="Chr9";      good_mark='yes'
                                elif one_id_lower.endswith('chr10'):   one_id="Chr10";      good_mark='yes'
                                elif one_id_lower.endswith('chr11'):   one_id="Chr11";      good_mark='yes'
                                elif one_id_lower.endswith('chr12'):   one_id="Chr12";      good_mark='yes'
                                elif one_id_lower.endswith('chr13'):   one_id="Chr13";      good_mark='yes'
                                elif one_id_lower.endswith('chr14'):   one_id="Chr14";      good_mark='yes'
                                elif one_id_lower.endswith('chr15'):   one_id="Chr15";      good_mark='yes'
                                elif one_id_lower.endswith('chr16'):   one_id="Chr16";      good_mark='yes'
                                elif one_id_lower.endswith('chr17'):   one_id="Chr17";      good_mark='yes'
                                elif one_id_lower.endswith('chr18'):   one_id="Chr18";      good_mark='yes'
                                elif one_id_lower.endswith('chr19'):   one_id="Chr19";      good_mark='yes'
                                elif one_id_lower.endswith('chr20'):   one_id="Chr20";      good_mark='yes'
                                if good_mark=='yes':
                                    chr_num+=1
                                    f2.write(f'>{one_id}\n{one_seq}\n')
                                    f3.write(f"{sample_name}\t{one_id_old}\t{one_id}\n")
                                    good_mark=''    
                            ##
                            one_id=eachline[1:]
                            one_seq=''
                        else:
                            one_seq+=eachline
                if len(one_seq)>1000000:
                    one_id_old=one_id
                    one_id_lower=one_id.lower()
                    if one_id.isdigit():
                        one_id="Chr"+one_id
                    elif one_id_lower.endswith('chr01') or one_id_lower.endswith('chr1'):   one_id="Chr1";      good_mark='yes'
                    elif one_id_lower.endswith('chr02') or one_id_lower.endswith('chr2'):   one_id="Chr2";      good_mark='yes'
                    elif one_id_lower.endswith('chr03') or one_id_lower.endswith('chr3'):   one_id="Chr3";      good_mark='yes'
                    elif one_id_lower.endswith('chr04') or one_id_lower.endswith('chr4'):   one_id="Chr4";      good_mark='yes'
                    elif one_id_lower.endswith('chr05') or one_id_lower.endswith('chr5'):   one_id="Chr5";      good_mark='yes'
                    elif one_id_lower.endswith('chr06') or one_id_lower.endswith('chr6'):   one_id="Chr6";      good_mark='yes'
                    elif one_id_lower.endswith('chr07') or one_id_lower.endswith('chr7'):   one_id="Chr7";      good_mark='yes'
                    elif one_id_lower.endswith('chr08') or one_id_lower.endswith('chr8'):   one_id="Chr8";      good_mark='yes'
                    elif one_id_lower.endswith('chr09') or one_id_lower.endswith('chr9'):   one_id="Chr9";      good_mark='yes'
                    elif one_id_lower.endswith('chr10'):   one_id="Chr10";      good_mark='yes'
                    elif one_id_lower.endswith('chr11'):   one_id="Chr11";      good_mark='yes'
                    elif one_id_lower.endswith('chr12'):   one_id="Chr12";      good_mark='yes'
                    elif one_id_lower.endswith('chr13'):   one_id="Chr13";      good_mark='yes'
                    elif one_id_lower.endswith('chr14'):   one_id="Chr14";      good_mark='yes'
                    elif one_id_lower.endswith('chr15'):   one_id="Chr15";      good_mark='yes'
                    elif one_id_lower.endswith('chr16'):   one_id="Chr16";      good_mark='yes'
                    elif one_id_lower.endswith('chr17'):   one_id="Chr17";      good_mark='yes'
                    elif one_id_lower.endswith('chr18'):   one_id="Chr18";      good_mark='yes'
                    elif one_id_lower.endswith('chr19'):   one_id="Chr19";      good_mark='yes'
                    elif one_id_lower.endswith('chr20'):   one_id="Chr20";      good_mark='yes'
                    if good_mark=='yes':
                        chr_num+=1
                        f2.write(f'>{one_id}\n{one_seq}\n')
                        f3.write(f"{sample_name}\t{one_id_old}\t{one_id}\n")
                        good_mark=''    
        if chr_num!=20 and chr_num!=19:
            with open("./chr2EDTA/error",'a') as f:
                f.write(f"{sample_name}的染色体数量不为19或20, 为{str(chr_num)}\n")
                ###
        subprocess.run([f"samtools faidx {output_dir}/{sample_name}.fasta "], shell=True)     
    
    # 将任务分配给进程池中的进程    
    with Pool(processes=thread) as pool:
        # 使用 imap 逐个获取结果
        for i, result in enumerate(pool.imap(run_step0, dir_file_name_list), start=1):
            # 这里可以处理结果，例如存储或打印
            progress = (i / len(dir_file_name_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush() 
    dir_file_name_list.sort()
    with open('./chr2EDTA/chrname_old2new','w') as f3:
        f3.write(f"sample_name\tone_id_old\tone_id\n")    
    for infos in dir_file_name_list:
        input_dir,file_name,sample_name=infos
        
        subprocess.run([f"cat ./chr2EDTA/0_prepare/{sample_name}/{sample_name}.chrname_old2new >> ./chr2EDTA/chrname_old2new"], shell=True) 
        

    
    
    
#step0
if argv1=="stepall" or argv1=="step0":
    if  os.path.exists('./chr2EDTA/0_prepare')==True:
        subprocess.run(["rm -r ./chr2EDTA/0_prepare"], shell=True)  
    subprocess.run(["mkdir ./chr2EDTA/0_prepare"], shell=True)  
    
    if  os.path.exists('./chr2EDTA/error')==True:subprocess.run(["rm ./chr2EDTA/error"], shell=True)  
    
    
    
    if "i"  in args_dict: #print("Missing input fasta file");sys.exit()
        input_dictionary = args_dict["i"]
        with open('./chr2EDTA/0_prepare/sample_source','w') as f:
            f.write(input_dictionary)
    else:
        print("Missing input fasta file");sys.exit()
    
    dir_fasta_file_num=0
    dir_file_name_list=[]
    files=os.listdir(input_dictionary)
    with open('./chr2EDTA/0_prepare/sample_list','w') as f:
        for one in files:
            if     one.endswith('.fasta'):one_name=one[:-6]
            elif   one.endswith('.fa'):one_name=one[:-3]
            elif   one.endswith('.fna'):one_name=one[:-4] 
            else :continue
            dir_fasta_file_num+=1
            dir_file_name_list.append([input_dictionary,one,one_name])
            f.write(one_name+'\n')

    print('Number of genome files: '+str(dir_fasta_file_num))

    def run_step0(dir_file):
        input_dir=dir_file[0]
        file_name=dir_file[1]
        sample_name=dir_file[2]
        input_file=input_dir+'/'+file_name
        output_dir="./chr2EDTA/0_prepare/"+sample_name+"/"
        #if  os.path.exists(output_dir)==True:return False
        output_file=output_dir+sample_name
        subprocess.run(["mkdir "+output_dir], shell=True)  
        ### Limit to length 1000000
        chr_num=0
        with open(f'{output_dir}/{sample_name}.chrname_old2new','w') as f3:
            with open(f'{output_dir}/{sample_name}.fasta','w') as f2:
                one_seq=''
                with open(input_file,'r') as f:
                    for line in f.readlines():
                        eachline_arr=line.strip().split(' ')
                        eachline=eachline_arr[0]
                        if eachline[0]=='>':
                            ##
                            if len(one_seq)>1000000:
                                one_id_old=one_id
                                one_id_lower=one_id.lower()
                                good_mark=''
                                if one_id.isdigit():
                                    one_id="Chr"+one_id
                                elif one_id_lower.endswith('chr01') or one_id_lower.endswith('chr1'):   one_id="Chr1";      good_mark='yes'
                                elif one_id_lower.endswith('chr02') or one_id_lower.endswith('chr2'):   one_id="Chr2";      good_mark='yes'
                                elif one_id_lower.endswith('chr03') or one_id_lower.endswith('chr3'):   one_id="Chr3";      good_mark='yes'
                                elif one_id_lower.endswith('chr04') or one_id_lower.endswith('chr4'):   one_id="Chr4";      good_mark='yes'
                                elif one_id_lower.endswith('chr05') or one_id_lower.endswith('chr5'):   one_id="Chr5";      good_mark='yes'
                                elif one_id_lower.endswith('chr06') or one_id_lower.endswith('chr6'):   one_id="Chr6";      good_mark='yes'
                                elif one_id_lower.endswith('chr07') or one_id_lower.endswith('chr7'):   one_id="Chr7";      good_mark='yes'
                                elif one_id_lower.endswith('chr08') or one_id_lower.endswith('chr8'):   one_id="Chr8";      good_mark='yes'
                                elif one_id_lower.endswith('chr09') or one_id_lower.endswith('chr9'):   one_id="Chr9";      good_mark='yes'
                                elif one_id_lower.endswith('chr10'):   one_id="Chr10";      good_mark='yes'
                                elif one_id_lower.endswith('chr11'):   one_id="Chr11";      good_mark='yes'
                                elif one_id_lower.endswith('chr12'):   one_id="Chr12";      good_mark='yes'
                                elif one_id_lower.endswith('chr13'):   one_id="Chr13";      good_mark='yes'
                                elif one_id_lower.endswith('chr14'):   one_id="Chr14";      good_mark='yes'
                                elif one_id_lower.endswith('chr15'):   one_id="Chr15";      good_mark='yes'
                                elif one_id_lower.endswith('chr16'):   one_id="Chr16";      good_mark='yes'
                                elif one_id_lower.endswith('chr17'):   one_id="Chr17";      good_mark='yes'
                                elif one_id_lower.endswith('chr18'):   one_id="Chr18";      good_mark='yes'
                                elif one_id_lower.endswith('chr19'):   one_id="Chr19";      good_mark='yes'
                                elif one_id_lower.endswith('chr20'):   one_id="Chr20";      good_mark='yes'
                                if good_mark=='yes':
                                    chr_num+=1
                                    f2.write(f'>{one_id}\n{one_seq}\n')
                                    f3.write(f"{sample_name}\t{one_id_old}\t{one_id}\n")
                                    good_mark=''    
                            ##
                            one_id=eachline[1:]
                            one_seq=''
                        else:
                            one_seq+=eachline
                if len(one_seq)>1000000:
                    one_id_old=one_id
                    one_id_lower=one_id.lower()
                    if one_id.isdigit():
                        one_id="Chr"+one_id
                    elif one_id_lower.endswith('chr01') or one_id_lower.endswith('chr1'):   one_id="Chr1";      good_mark='yes'
                    elif one_id_lower.endswith('chr02') or one_id_lower.endswith('chr2'):   one_id="Chr2";      good_mark='yes'
                    elif one_id_lower.endswith('chr03') or one_id_lower.endswith('chr3'):   one_id="Chr3";      good_mark='yes'
                    elif one_id_lower.endswith('chr04') or one_id_lower.endswith('chr4'):   one_id="Chr4";      good_mark='yes'
                    elif one_id_lower.endswith('chr05') or one_id_lower.endswith('chr5'):   one_id="Chr5";      good_mark='yes'
                    elif one_id_lower.endswith('chr06') or one_id_lower.endswith('chr6'):   one_id="Chr6";      good_mark='yes'
                    elif one_id_lower.endswith('chr07') or one_id_lower.endswith('chr7'):   one_id="Chr7";      good_mark='yes'
                    elif one_id_lower.endswith('chr08') or one_id_lower.endswith('chr8'):   one_id="Chr8";      good_mark='yes'
                    elif one_id_lower.endswith('chr09') or one_id_lower.endswith('chr9'):   one_id="Chr9";      good_mark='yes'
                    elif one_id_lower.endswith('chr10'):   one_id="Chr10";      good_mark='yes'
                    elif one_id_lower.endswith('chr11'):   one_id="Chr11";      good_mark='yes'
                    elif one_id_lower.endswith('chr12'):   one_id="Chr12";      good_mark='yes'
                    elif one_id_lower.endswith('chr13'):   one_id="Chr13";      good_mark='yes'
                    elif one_id_lower.endswith('chr14'):   one_id="Chr14";      good_mark='yes'
                    elif one_id_lower.endswith('chr15'):   one_id="Chr15";      good_mark='yes'
                    elif one_id_lower.endswith('chr16'):   one_id="Chr16";      good_mark='yes'
                    elif one_id_lower.endswith('chr17'):   one_id="Chr17";      good_mark='yes'
                    elif one_id_lower.endswith('chr18'):   one_id="Chr18";      good_mark='yes'
                    elif one_id_lower.endswith('chr19'):   one_id="Chr19";      good_mark='yes'
                    elif one_id_lower.endswith('chr20'):   one_id="Chr20";      good_mark='yes'
                    if good_mark=='yes':
                        chr_num+=1
                        f2.write(f'>{one_id}\n{one_seq}\n')
                        f3.write(f"{sample_name}\t{one_id_old}\t{one_id}\n")
                        good_mark=''    
        if chr_num!=20 and chr_num!=19:
            with open("./chr2EDTA/error",'a') as f:
                f.write(f"Chromosome count for {sample_name} is not 19 or 20, it is {str(chr_num)}\n")
                ###
        subprocess.run([f"samtools faidx {output_dir}/{sample_name}.fasta "], shell=True)     
    
    # Assign tasks to processes in the process pool    
    with Pool(processes=thread) as pool:
        # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step0, dir_file_name_list), start=1):
            # Results can be processed here, e.g., stored or printed
            progress = (i / len(dir_file_name_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush() 
    dir_file_name_list.sort()
    with open('./chr2EDTA/chrname_old2new','w') as f3:
        f3.write(f"sample_name\tone_id_old\tone_id\n")    
    for infos in dir_file_name_list:
        input_dir,file_name,sample_name=infos
        
        subprocess.run([f"cat ./chr2EDTA/0_prepare/{sample_name}/{sample_name}.chrname_old2new >> ./chr2EDTA/chrname_old2new"], shell=True) 
        

    
    
    
if argv1=="stepall" or argv1=="step1":    
    if  os.path.exists('./chr2EDTA/1_EDTA')!=True:
        subprocess.run(["mkdir ./chr2EDTA/1_EDTA"], shell=True)  

    
    with open('./chr2EDTA/0_prepare/sample_list','r') as f:
        sample_list=f.read().split('\n')
    sample_list.sort()
    
    
    #sample_list=["VHP-T2T.hap1","VHP-T2T.hap2"]+sample_list        # Run some first
    
    
    print(sample_list)
    def run_step1(sample_name):
        #if sample_name!='PN40024_hap1': return False
        if len(sample_name)==0:return False
        
        #input_file=f"./chr2EDTA/0_prepare/{sample_name}/{sample_name}.fa"
        output_dir="./chr2EDTA/1_EDTA/"+sample_name+"/"
        #if  os.path.exists(f"{output_dir}/{sample_name}.fa.pass.list")==True:return False   
        if  os.path.exists(f"{output_dir}/{sample_name}.fasta.mod.EDTA.TEanno.gff3")==True:return False   
        print(sample_name)
        subprocess.run(["mkdir "+output_dir], shell=True) 
        os.chdir(output_dir)
        subprocess.run([f"source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate EDTA && EDTA.pl --genome ../../../chr2EDTA/0_prepare/{sample_name}/{sample_name}.fasta -species others --step all --threads 20 --u 2.5e-9 --sensitive 1 --anno 1"], shell=True) 
        ####2025.3.19, added --anno 1, seems to get all incomplete TEs
        os.chdir('../../../')
        
    # Assign tasks to processes in the process pool    
    with Pool(processes=3) as pool:
        # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step1, sample_list), start=1):
            # Results can be processed here, e.g., stored or printed
            progress = (i / len(sample_list)) * 100
            sys.stdout.write(f"\nProgress: {progress:.2f}% {str(i)} / {str(len(sample_list))}#####################################################################\n")
            sys.stdout.flush()    
if argv1=="stepall" or argv1=="step_delete":   
    with open('./chr2EDTA/0_prepare/sample_list','r') as f:
        sample_list=f.read().split('\n')
    sample_list.sort()    
    print(sample_list)
    for one_sample in sample_list:
        subprocess.run([f"rm  ./chr2EDTA/1_EDTA/{one_sample}/{one_sample}.fasta.mod.EDTA.raw"], shell=True) 
   
if argv1=="stepall" or argv1=="step2":    
    if  os.path.exists('./chr2EDTA/2_TEsorter ')!=True:
        subprocess.run(["mkdir ./chr2EDTA/2_TEsorter"], shell=True)  

    with open('./chr2EDTA/0_prepare/sample_list','r') as f:
        sample_list=f.read().split('\n')
    sample_list.sort()   
   
    print(sample_list)
    def run_step1(sample_name):
        #if sample_name!='PN40024_hap1': return False
        if len(sample_name)==0:return False
        
        #input_file=f"./chr2EDTA/0_prepare/{sample_name}/{sample_name}.fa"
        output_dir="./chr2EDTA/2_TEsorter/"+sample_name+"/"
        #if  os.path.exists(f"{output_dir}/{sample_name}.fa.pass.list")==True:return False   
        #if  os.path.exists(f"{output_dir}/{sample_name}.fasta.mod.EDTA.TEanno.gff3")==True:return False   
        print(sample_name)
        subprocess.run(["mkdir "+output_dir], shell=True) 
        subprocess.run([f"cp ./chr2EDTA/1_EDTA/{sample_name}/{sample_name}.fasta.mod.EDTA.intact.fa {output_dir}/EDTA.fa"], shell=True) 
        
        os.chdir(output_dir)
        
        subprocess.run([f"source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate EDTA && TEsorter -p 20 -db rexdb-plant EDTA.fa"], shell=True) 
        ####2025.3.19, added --anno 1, seems to get all incomplete TEs
        os.chdir('../../../')
        
    # Assign tasks to processes in the process pool    
    with Pool(processes=3) as pool:
        # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step1, sample_list), start=1):
            # Results can be processed here, e.g., stored or printed
            progress = (i / len(sample_list)) * 100
            sys.stdout.write(f"\nProgress: {progress:.2f}% {str(i)} / {str(len(sample_list))}#####################################################################\n")
            sys.stdout.flush()      
            
if argv1=="stepall" or argv1=="step2_TEsorter_plant":
    if  os.path.exists('./chr2EDTA/2_TEsorter_protinput/ ')!=True:
        subprocess.run(["mkdir ./chr2EDTA/2_TEsorter_protinput"], shell=True)      
    input_file=sys.argv[2]
    os.chdir("./chr2EDTA/2_TEsorter_protinput")
    subprocess.run([f"source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate EDTA && TEsorter -p 20 -db rexdb-plant ../../{input_file}"], shell=True) 
if argv1=="stepall" or argv1=="step2_TEsorter_common":
    if  os.path.exists('./chr2EDTA/2_TEsorter_protinput/ ')==True:
        subprocess.run(["rm -r ./chr2EDTA/2_TEsorter_protinput"], shell=True)    
    subprocess.run(["mkdir ./chr2EDTA/2_TEsorter_protinput"], shell=True)      
    input_file=sys.argv[2]
    os.chdir("./chr2EDTA/2_TEsorter_protinput")
    subprocess.run([f"source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate EDTA && TEsorter -p 20 -db rexdb ../../{input_file}"], shell=True)     
        
    
if argv1=="stepall" or "step3" in argv1:  
    print('stat')
    if  os.path.exists('./chr2EDTA/3_EDTA_stat')!=True:
        subprocess.run(["mkdir ./chr2EDTA/3_EDTA_stat"], shell=True)    
        
    if argv1=="stepall" or argv1=="step3_readme":
        print('Printing instructions')
        with open('./chr2EDTA/3_EDTA_stat/readme','w') as f:
            txt=r'''
            3.0————
            0_VSat1         Inside/outside VSat1                         0_VSat1_innerTE records specific positions
            0_interarray    Inside/outside interarray
            0_coreVSat1     Inside/outside coreVSat1, core extended by 1M
            
            '''
            f.write(txt)     
    if argv1=="stepall"  or argv1=="step3" or argv1=="step3.0":    
        print('')  
                          
        print('Loading a sample information table ../samples_satellite/sample_info')
        dict_sample_type={}
        with open('../samples_satellite/sample_info') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=6:continue
                sample_name=eachline_arr[0]
                sample_type=eachline_arr[3]
                if sample_type not in ["Eurasian","East_Asia","America"]:continue
                dict_sample_type[sample_name]=sample_type      
        sample_list=list(dict_sample_type.keys())
        sample_list.sort()
        #sample_list=['PN40024']
        print('Loading a sample information table ../samples_satellite/2_good_regions_interarray')
        dict_samplechr_interarrays={}
        dict_samplechr_interarray={}
        dict_samplechr_bigblock_expand={}
        with open('../samples_satellite/2_good_regions_interarray') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                sample_name=eachline_arr[0]
                chromosome=eachline_arr[1]
                samplechr=sample_name+'___'+chromosome
                inter_array_length= int(eachline_arr[8])
                ########
                dict_samplechr_interarray[samplechr]=inter_array_length
                ##########
                if eachline_arr[9]!="NA"    :
                    if samplechr not in dict_samplechr_interarrays:
                        dict_samplechr_interarrays[samplechr]=[]                    
                    interarrays=eachline_arr[9].split('|')
                    for one_inter_array in interarrays:
                        pos_arr=one_inter_array.split('-')
                        dict_samplechr_interarrays[samplechr].append([         int(pos_arr[0])       ,     int(pos_arr[1])   ])
                ####
                bigblock_start=     int(eachline_arr[2])
                bigblock_end=       int(eachline_arr[3])
                bigblock_len=       int(eachline_arr[4])
                dict_samplechr_bigblock_expand[samplechr]=[bigblock_start-1000000,bigblock_end+1000000,bigblock_len]

        print('Loading a sample information table ../samples_satellite/2_good_regions')
        dict_samplechr_allblocks={}
        with open('../samples_satellite/2_good_regions') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=8:continue
                sample_name=eachline_arr[0]
                chromosome=eachline_arr[1].replace('region_','Chr')
                samplechr=sample_name+'___'+chromosome
                start=          int(eachline_arr[3])
                end=            int(eachline_arr[4])
                block_length= int(eachline_arr[5])
                ########
                if block_length>10000:
                    if samplechr not in dict_samplechr_allblocks:dict_samplechr_allblocks[samplechr]=[]
                    dict_samplechr_allblocks[samplechr].append([         start      ,     end   ])                    
                
        print('Counting complete TEs')      
        with open(f'./chr2EDTA/3_EDTA_stat/0_VSat1','w') as f2:         ## Use VSat1 >10000bp from 2_good_regions
            f2.write(f"sample\tsample_type\tClass\tin_num\tout_num\n")
        with open(f'./chr2EDTA/3_EDTA_stat/0_coreVSat1','w') as f2:     
            f2.write(f"sample\tsample_type\tClass\tin_num\tout_num\tin_region\tout_region\tin_density\tout_density\tdelta_density\n")            
        with open(f'./chr2EDTA/3_EDTA_stat/0_interarray','w') as f2:     
            f2.write(f"sample\tsample_type\tClass\tin_num\tout_num\tin_region\tout_region\tin_density\tout_density\tdelta_density\n")            
        with open(f'./chr2EDTA/3_EDTA_stat/0_VSat1_innerTE','w') as ff:   ### These are specific positions
            ff.write(f"\n")            
        def run_step(one_sample):
            #
            dict_Sat1_inall={}
            dict_CoreVSat1_inall={}
            dict_LIR_inall={}
            ##
            ## Load the total length of a chromosome
            dict_samplechr_length={}
            with open(f"./chr2EDTA/0_prepare/{one_sample}/{one_sample}.fasta.fai",'r') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    samplechr=one_sample+'___'+eachline_arr[0]
                    dict_samplechr_length[samplechr]=int(eachline_arr[1])
            ##
            dict_sample_CoreVSat1_inout={}
            for samplechr,bigblock_expand_info in dict_samplechr_bigblock_expand.items():
                sample,chromosome=samplechr.split('___')
                bigblock_expand_start,bigblock_expand_end,bigblock_len=bigblock_expand_info
                if sample!=one_sample:continue
                if sample not in dict_sample_CoreVSat1_inout:
                    dict_sample_CoreVSat1_inout[sample]={}
                    dict_sample_CoreVSat1_inout[sample]['in']=0
                    dict_sample_CoreVSat1_inout[sample]['all']=0
                dict_sample_CoreVSat1_inout[sample]['all']+=dict_samplechr_length[samplechr]
                dict_sample_CoreVSat1_inout[sample]['in']+=bigblock_len
            ##
            dict_sample_LIR_inout={}
            for samplechr,interarrays in dict_samplechr_interarrays.items():
                sample,chromosome=samplechr.split('___')
                if sample!=one_sample:continue
                if sample not in dict_sample_LIR_inout:
                    dict_sample_LIR_inout[sample]={}
                    dict_sample_LIR_inout[sample]['in']=0
                    dict_sample_LIR_inout[sample]['all']=0
                dict_sample_LIR_inout[sample]['all']+=dict_samplechr_length[samplechr]
                for one_inter_array in interarrays:
                    one_inter_array_start,one_inter_array_end=one_inter_array
                    dict_sample_LIR_inout[sample]['in']+=one_inter_array_end-one_inter_array_start
            ##
            dict_Sat1_type_num={}
            dict_CoreVSat1_type_num={}
            dict_LIR_type_num={}
            with open(f"./chr2EDTA/1_EDTA/{one_sample}/{one_sample}.fasta.mod.EDTA.TEanno.gff3",'r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if len(eachline_arr)!=9:continue
                    chromosome= eachline_arr[0]
                    EDTA_class=       eachline_arr[2] 
                    start=      int(eachline_arr[3])
                    end=        int(eachline_arr[4])
                    #if one_type not in ["Gypsy_LTR_retrotransposon","LTR_retrotransposon","Copia_LTR_retrotransposon"]:continue
                    if EDTA_class not in ["Gypsy_LTR_retrotransposon","LTR_retrotransposon","PIF_Harbinger_TIR_transposon","Copia_LTR_retrotransposon","CACTA_TIR_transposon","hAT_TIR_transposon","helitron","Mutator_TIR_transposon","Tc1_Mariner_TIR_transposon"]:continue
                    if EDTA_class=="Gypsy_LTR_retrotransposon":         EDTA_class='Gypsy'
                    elif EDTA_class=="LTR_retrotransposon":         EDTA_class='OtherLTR'
                    elif EDTA_class=="PIF_Harbinger_TIR_transposon":         EDTA_class='PIF'
                    elif EDTA_class=="Copia_LTR_retrotransposon":         EDTA_class='Copia'
                    elif EDTA_class=="CACTA_TIR_transposon":         EDTA_class='CACTA'
                    elif EDTA_class=="hAT_TIR_transposon":         EDTA_class='hTA'
                    elif EDTA_class=="helitron":         EDTA_class='helitron'
                    elif EDTA_class=="Mutator_TIR_transposon":         EDTA_class='Mutator'
                    elif EDTA_class=="Tc1_Mariner_TIR_transposon":         EDTA_class='Tc1'
                    else: continue
                    ######
                    samplechr=one_sample+'___'+chromosome
                    ## Count inside and outside VSat1 >10000bp
                    if samplechr in dict_samplechr_allblocks:
                        VSat1_list=dict_samplechr_allblocks[samplechr]
                        VSat1_mark=''
                        for one_VSat1 in VSat1_list:
                            VSat1_start,VSat1_end=one_VSat1
                            if start>VSat1_start and end<VSSat1_end: VSat1_mark='yes';break
                        ###
                        ###############
                        if EDTA_class not in dict_Sat1_type_num:
                            dict_Sat1_type_num[EDTA_class]={}
                            dict_Sat1_type_num[EDTA_class]['in']=0
                            dict_Sat1_type_num[EDTA_class]['out']=0
                        if VSat1_mark=='yes':      
                            dict_Sat1_type_num[EDTA_class]['in']+=1    
                            with open(f'./chr2EDTA/3_EDTA_stat/0_VSat1_innerTE','a') as ff:
                                ff.write(f"{samplechr}\t{line.strip()}\n")
                        else:       dict_Sat1_type_num[EDTA_class]['out']+=1        
                            
                    ## Count near CoreVsat1
                    if samplechr in dict_samplechr_bigblock_expand:
                        core_start,core_end,bigblock_len=dict_samplechr_bigblock_expand[samplechr]
                        if bigblock_len>100000:
                            coreVSat1_mark=''
                            if start>core_start and end<core_end:           coreVSat1_mark='yes'
                            ###
                            ###############
                            if EDTA_class not in dict_CoreVSat1_type_num:
                                dict_CoreVSat1_type_num[EDTA_class]={}
                                dict_CoreVSat1_type_num[EDTA_class]['in']=0
                                dict_CoreVSat1_type_num[EDTA_class]['out']=0
                            if coreVSat1_mark=='yes':      dict_CoreVSat1_type_num[EDTA_class]['in']+=1    
                            else:                           dict_CoreVSat1_type_num[EDTA_class]['out']+=1                           
                    ## Only keep those with interarray
                    if samplechr in dict_samplechr_interarray:
                        if dict_samplechr_interarray[samplechr]>150000 and dict_samplechr_interarray[samplechr]<350000:
                            ###
                            interarrays=    dict_samplechr_interarrays[samplechr]
                            interarray_mark=''
                            for one_inter_array in interarrays:
                                inter_array_start,inter_array_end=one_inter_array
                                if start>inter_array_start and end<inter_array_end: interarray_mark='yes';break
                            ###
                            ###############
                            if EDTA_class not in dict_LIR_type_num:
                                dict_LIR_type_num[EDTA_class]={}      
                                dict_LIR_type_num[EDTA_class]['in']=0
                                dict_LIR_type_num[EDTA_class]['out']=0
                            if interarray_mark=='yes':      dict_LIR_type_num[EDTA_class]['in']+=1    
                            else:                           dict_LIR_type_num[EDTA_class]['out']+=1        
            #print('Output')
            with open(f'./chr2EDTA/3_EDTA_stat/0_VSat1','a') as f2:    
                #with open(f'./chr2EDTA/4_TEsorter_stat/{one_sample}.coreVSat1_stat','w') as f:
                    #f.write(f"sample\tsample_type\tClass\tClass1\tClass2\tClass3\tinter_array_num\tout_array_num\n")
                for EDTA_class,dict_info in dict_Sat1_type_num.items():
                    in_num=     dict_info['in']
                    out_num=    dict_info['out']
                    sample_type=dict_sample_type[one_sample]
                    f2.write(f"{one_sample}\t{sample_type}\t{EDTA_class}\t{in_num}\t{out_num}\n")
                    #f.write(f"{one_sample}\t{sample_type}\t{TEsorter_class1}|{TEsorter_class2}|{TEsorter_class3}\t{TEsorter_class1}\t{TEsorter_class2}\t{TEsorter_class3}\t{in_num}\t{out_num}\n")
                           
            with open(f'./chr2EDTA/3_EDTA_stat/0_coreVSat1','a') as f2:    
                #with open(f'./chr2EDTA/4_TEsorter_stat/{one_sample}.coreVSat1_stat','w') as f:
                    #f.write(f"sample\tsample_type\tClass\tClass1\tClass2\tClass3\tinter_array_num\tout_array_num\n")
                for EDTA_class,dict_info in dict_CoreVSat1_type_num.items():

                    in_num=     dict_info['in']
                    out_num=    dict_info['out']
                    sample_type=dict_sample_type[one_sample]
                    #
                    allin=dict_sample_CoreVSat1_inout[one_sample]['in']/1000000
                    alllength=dict_sample_CoreVSat1_inout[one_sample]['all']/1000000
                    allout=alllength-allin
                    in_density=round(in_num/allin,4)
                    out_density=round(out_num/allout,4)
                    delta_density=in_density-out_density
                    #
                    f2.write(f"{one_sample}\t{sample_type}\t{EDTA_class}\t{in_num}\t{out_num}\t{allin:.2f}\t{allout:.2f}\t{in_density}\t{out_density}\t{delta_density}\n")
                    #f.write(f"{one_sample}\t{sample_type}\t{TEsorter_class1}|{TEsorter_class2}|{TEsorter_class3}\t{TEsorter_class1}\t{TEsorter_class2}\t{TEsorter_class3}\t{in_num}\t{out_num}\n")
                            
            with open(f'./chr2EDTA/3_EDTA_stat/0_interarray','a') as f2:   
                #with open(f'./chr2EDTA/4_TEsorter_stat/{one_sample}.interarray_stat','w') as f:
                    #f.write(f"sample\tsample_type\tClass\tClass1\tClass2\tClass3\tinter_array_num\tout_array_num\n")
                for EDTA_class,dict_info in dict_LIR_type_num.items():
                    in_num=     dict_info['in']
                    out_num=    dict_info['out']
                    sample_type=dict_sample_type[one_sample]
                    #
                    allin=dict_sample_LIR_inout[one_sample]['in']/1000000
                    alllength=dict_sample_LIR_inout[one_sample]['all']/1000000
                    allout=alllength-allin
                    in_density=round(in_num/allin,4)
                    out_density=round(out_num/allout,4)
                    delta_density=in_density-out_density
                    #
                    f2.write(f"{one_sample}\t{sample_type}\t{EDTA_class}\t{in_num}\t{out_num}\t{allin:.2f}\t{allout:.2f}\t{in_density}\t{out_density}\t{delta_density}\n")
                        #f.write(f"{one_sample}\t{sample_type}\t{TEsorter_class1}|{TEsorter_class2}|{TEsorter_class3}\t{TEsorter_class1}\t{TEsorter_class2}\t{TEsorter_class3}\t{in_num}\t{out_num}\n")
                                
            
                    
                    
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, sample_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(sample_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()          
    if argv1=="stepall"  or argv1=="step3" or argv1=="step3.0" or argv1=="step3.0p":                                
        print('Plotting')            
        R_txt=f'''library(ggplot2)
library(dplyr)
# Set working directory
setwd('./')

# Read data
input_file1 <- read.table('0_VSat1', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class, y =in_num/(in_num+out_num), color = sample_type),position = position_jitter(width = 0.3)) 
p <- p + theme_bw()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_VSat1.pdf", width = 20 / 2.54, height = 50 / 2.54)
print(p)
dev.off()


# Read data
input_file1 <- read.table('0_coreVSat1', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class, y =in_num/(in_num+out_num), color = sample_type),position = position_jitter(width = 0.3)) 
p <- p + theme_bw()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_corVSat1.pdf", width = 20 / 2.54, height = 50 / 2.54)
print(p)
dev.off()

# Read data
input_file1 <- read.table('0_coreVSat1', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')


# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class, y =delta_density, color = delta_density),position = position_jitter(width = 0.3)) 
p <- p + theme_bw()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_corVSat1_density.pdf", width = 20 / 2.54, height = 50 / 2.54)
print(p)
dev.off()

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class, y =delta_density, color = sample_type),position = position_jitter(width = 0.3)) 
p <- p + theme_bw()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_corVSat1_density_sample_type.pdf", width = 20 / 2.54, height = 50 / 2.54)
print(p)
dev.off()

# Read data
input_file1 <- read.table('0_interarray', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class, y =in_num/(in_num+out_num), color = sample_type),position = position_jitter(width = 0.3)) 
p <- p + theme_bw()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_interarray.pdf", width = 20 / 2.54, height = 50 / 2.54)
print(p)
dev.off()

# Read data
input_file1 <- read.table('0_interarray', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class, y =delta_density, color = delta_density),position = position_jitter(width = 0.3)) 
p <- p + theme_bw()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_interarray_density.pdf", width = 20 / 2.54, height = 50 / 2.54)
print(p)
dev.off()

# Read data
input_file1 <- read.table('0_interarray', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class, y =delta_density, color = sample_type),position = position_jitter(width = 0.3)) 
p <- p + theme_bw()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_interarray_density_sample_type.pdf", width = 20 / 2.54, height = 50 / 2.54)
print(p)
dev.off()
'''
        with open(f'./chr2EDTA/3_EDTA_stat//plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = f"./chr2EDTA/3_EDTA_stat/"
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')                                       
   
if argv1=="stepall" or "step4" in argv1:  
    print('stat, complete TEs')
    if  os.path.exists('./chr2EDTA/4_TEsorter_stat')!=True:
        subprocess.run(["mkdir ./chr2EDTA/4_TEsorter_stat"], shell=True)    
        
    if argv1=="stepall" or argv1=="step4_readme":
        print('Printing instructions')
        with open('./chr2EDTA/4_TEsorter_stat/readme','w') as f:
            txt=r'''
            4.0————
            0_VSat1         Inside/outside VSat1                         0_VSat1_innerTE records specific positions
            0_interarray    Inside/outside interarray (all interarrays in step4 are LIRs)
            0_coreVSat1     Inside/outside coreVSat1, core extended by 1M
            4.1———— Plotting peri-VSat1 and LIR together in one figure doesn't look good
            4.2———— Plot peri-VSat1 and LIR in two figures, using facet makes it barely acceptable
            4.3———— Merge identical samples, e.g., hap1, hap2, etc.
            '''
            f.write(txt)     
    if argv1=="stepall"  or argv1=="step4" or argv1=="step4.0":    
        print('')  
                          
        print('Loading a sample information table ../samples_satellite/sample_info')
        dict_sample_type={}
        with open('../samples_satellite/sample_info') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=6:continue
                sample_name=eachline_arr[0]
                sample_type=eachline_arr[3]
                #if sample_type not in ["Eurasian","East_Asia","America"]:continue
                dict_sample_type[sample_name]=sample_type      
        sample_list=list(dict_sample_type.keys())
        sample_list.sort()
        #sample_list=['PN40024']
        print('Loading a sample information table ../samples_satellite/2_good_regions_interarray_LIR')
        dict_samplechr_interarrays={}
        dict_samplechr_interarray={}
        dict_samplechr_bigblock_expand={}
        with open('../samples_satellite/2_good_regions_interarray_LIR') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                sample_name=eachline_arr[0]
                chromosome=eachline_arr[1]
                ###if chromosome not in ['Chr7','Chr14','Chr17']:continue
                samplechr=sample_name+'___'+chromosome
                inter_array_length= int(eachline_arr[8])
                ########
                dict_samplechr_interarray[samplechr]=inter_array_length
                ##########
                if eachline_arr[-1]!="NA"    :     ####### Originally [9], changed to -1, i.e., LIR determined by moddotplot
                    if samplechr not in dict_samplechr_interarrays:
                        dict_samplechr_interarrays[samplechr]=[]                    
                    interarrays=eachline_arr[-1].split('|')
                    for one_inter_array in interarrays:
                        pos_arr=one_inter_array.split('-')
                        dict_samplechr_interarrays[samplechr].append([         int(pos_arr[0])       ,     int(pos_arr[1])   ])
                ####
                bigblock_start=     int(eachline_arr[2])
                bigblock_end=       int(eachline_arr[3])
                bigblock_len=       int(eachline_arr[4])
                dict_samplechr_bigblock_expand[samplechr]=[bigblock_start-1000000,bigblock_end+1000000,bigblock_len]
        
        print('Loading a sample information table ../samples_satellite/2_good_regions')
        dict_samplechr_allblocks={}
        with open('../samples_satellite/2_good_regions','r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=8:continue
                sample_name=eachline_arr[0]
                chromosome=eachline_arr[1].replace('region_','Chr')
                samplechr=sample_name+'___'+chromosome
                start=          int(eachline_arr[3])
                end=            int(eachline_arr[4])
                block_length= int(eachline_arr[5])
                ########
                if block_length>10000:
                    if samplechr not in dict_samplechr_allblocks:dict_samplechr_allblocks[samplechr]=[]
                    dict_samplechr_allblocks[samplechr].append([         start      ,     end   ])          

        if 1==2:
            print('In this case, the VSat1 output below is actually the region between 66 and 103, but the result is not better; maybe the 2_good_regions_main_66_103 file is a bit rough')
            print('Loading a sample information table ../samples_satellite/2_good_regions_main_66_103')
            dict_samplechr_allblocks={}
            with open('../samples_satellite/2_good_regions_main_66_103','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if len(eachline_arr)!=5:continue
                    sample_name=eachline_arr[0]
                    chromosome=eachline_arr[1]
                    samplechr=sample_name+'___'+chromosome
                    start=          int(eachline_arr[2])
                    end=            int(eachline_arr[3])
                    block_length= int(eachline_arr[4])
                    ########
                    if block_length>10000:
                        if samplechr not in dict_samplechr_allblocks:dict_samplechr_allblocks[samplechr]=[]
                        dict_samplechr_allblocks[samplechr].append([         start      ,     end   ])  
                    
        print('Loading a sample information table /home/lain/aaa_data/run0/stat_plot/0-region2info') 
        dict_samplechr_type_allblocks={}
        dict_blocktype_num={}
        with open('/home/lain/aaa_data/run0/stat_plot/0-region2info','r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=9:continue
                sample_name=eachline_arr[0]
                chromosome=eachline_arr[2]
                centype=eachline_arr[3]
                strand=eachline_arr[7]
                samplechr=sample_name+'___'+chromosome
                start=          int(eachline_arr[4])
                end=            int(eachline_arr[5])
                block_length= int(eachline_arr[6])
                ########
                if block_length>10000:
                    
                    if samplechr not in dict_samplechr_type_allblocks:dict_samplechr_type_allblocks[samplechr]={}
                    if centype not in dict_samplechr_type_allblocks[samplechr]:dict_samplechr_type_allblocks[samplechr][centype]=[]
                    dict_samplechr_type_allblocks[samplechr][centype].append([         start      ,     end  ,strand ])
                    ## Store a searched total count of CENTYPE
                    if centype not in dict_blocktype_num:dict_blocktype_num[centype]=0
                    dict_blocktype_num[centype]+=1
        #print(dict_samplechr_type_allblocks)
        #sys.exit()
        print('Counting complete TEs')      
        with open(f'./chr2EDTA/4_TEsorter_stat/0_VSat1','w') as f2:       ## Use VSat1 >10000bp from 2_good_regions
            f2.write(f"sample\tsample_type\tClass\tClass1\tClass2\tClass3\tin_num\tout_num\n")
        with open(f'./chr2EDTA/4_TEsorter_stat/0_cen_stat','w') as f2:       ## Use VSat1 >10000bp from 2_good_regions
            f2.write(f"centype\tsample\tsample_type\tClass\tClass1\tClass2\tClass3\tpart1_num\tpart2_num\tpart3_num\tpart4_num\tpart5_num\tpart_out_num\tpart_all_num\tpart1_percent\tpart2_percent\tpart3_percent\tpart4_percent\tpart5_percent\tpart_out_percent\n")

        ############
        with open(f'./chr2EDTA/4_TEsorter_stat/0_coreVSat1','w') as f2:     
            f2.write(f"sample\tsample_type\tClass\tClass1\tClass2\tClass3\tin_num\tout_num\tin_region\tout_region\tin_density\tout_density\tdelta_density\n")            
        with open(f'./chr2EDTA/4_TEsorter_stat/0_interarray','w') as f2:     
            f2.write(f"sample\tsample_type\tClass\tClass1\tClass2\tClass3\tin_num\tout_num\tin_region\tout_region\tin_density\tout_density\tdelta_density\n")            
        with open(f'./chr2EDTA/4_TEsorter_stat/0_VSat1_innerTE','w') as ff:      ## Specific positions
            ff.write(f"#TE\tOrder\tSuperfamily\tClade\tComplete\tStrand\tDomains\n")          
            
            
        def run_step(one_sample):
            #
            dict_Sat1_inall={}
            dict_CoreVSat1_inall={}
            dict_LIR_inall={}
            ##
            ## Load the total length of a chromosome
            dict_samplechr_length={}
            with open(f"./chr2EDTA/0_prepare/{one_sample}/{one_sample}.fasta.fai",'r') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    samplechr=one_sample+'___'+eachline_arr[0]
                    dict_samplechr_length[samplechr]=int(eachline_arr[1])
            ##
            dict_sample_CoreVSat1_inout={}
            for samplechr,bigblock_expand_info in dict_samplechr_bigblock_expand.items():
                sample,chromosome=samplechr.split('___')
                bigblock_expand_start,bigblock_expand_end,bigblock_len=bigblock_expand_info
                if sample!=one_sample:continue
                if sample not in dict_sample_CoreVSat1_inout:
                    dict_sample_CoreVSat1_inout[sample]={}
                    dict_sample_CoreVSat1_inout[sample]['in']=0
                    dict_sample_CoreVSat1_inout[sample]['all']=0
                dict_sample_CoreVSat1_inout[sample]['all']+=dict_samplechr_length[samplechr]
                dict_sample_CoreVSat1_inout[sample]['in']+=bigblock_expand_end-bigblock_expand_start
            ##
            dict_sample_LIR_inout={}
            for samplechr,interarrays in dict_samplechr_interarrays.items():
                sample,chromosome=samplechr.split('___')
                if sample!=one_sample:continue
                if sample not in dict_sample_LIR_inout:
                    dict_sample_LIR_inout[sample]={}
                    dict_sample_LIR_inout[sample]['in']=0
                    dict_sample_LIR_inout[sample]['all']=0
                dict_sample_LIR_inout[sample]['all']+=dict_samplechr_length[samplechr]
                for one_inter_array in interarrays:
                    one_inter_array_start,one_inter_array_end=one_inter_array
                    dict_sample_LIR_inout[sample]['in']+=one_inter_array_end-one_inter_array_start
            ##
            dict_Sat1_type_num={}
            dict_CoreVSat1_type_num={}
            dict_LIR_type_num={}
            
            dict_centype_type_num={}

            
            
            with open (f'./chr2EDTA/2_TEsorter/{one_sample}/EDTA.fa.rexdb-plant.cls.tsv','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    TEsorter_class1=eachline_arr[1]
                    TEsorter_class2=eachline_arr[2]
                    TEsorter_class3=eachline_arr[3]
                    Complete=       eachline_arr[4]
                    Strand=         eachline_arr[5]	
                    Domains=        eachline_arr[6]
                    if Complete!='yes':continue
                    ######
                    #print(eachline_arr[0])
                    pos_type_arr=eachline_arr[0].split("_Chr")[1].split('#')
                    pos_arr=pos_type_arr[0].split(":")
                    chromosome='Chr'+pos_arr[0]
                    pos_arr2=pos_arr[1].split('..')
                    start=int(pos_arr2[0])
                    end=int(pos_arr2[1])
                    type_arr=pos_type_arr[1].split('/')
                    EDTA_class1=type_arr[0]
                    EDTA_class2=type_arr[1]   
                    if TEsorter_class2!=EDTA_class2:continue
                    samplechr=one_sample+'___'+chromosome
                    ## Count inside and outside VSat1 >10000bp
                    if samplechr in dict_samplechr_allblocks:
                        VSat1_list=dict_samplechr_allblocks[samplechr]
                        VSat1_mark=''
                        for one_VSat1 in VSat1_list:
                            VSat1_start,VSat1_end=one_VSat1
                            if start>VSat1_start and end<VSat1_end: VSat1_mark='yes';break
                        ###
                        ###############
                        
                        if TEsorter_class1 not in dict_Sat1_type_num:
                            dict_Sat1_type_num[TEsorter_class1]={}
                        if TEsorter_class2 not in dict_Sat1_type_num[TEsorter_class1]:
                            dict_Sat1_type_num[TEsorter_class1][TEsorter_class2]={}
                        if TEsorter_class3 not in dict_Sat1_type_num[TEsorter_class1][TEsorter_class2]:
                            dict_Sat1_type_num[TEsorter_class1][TEsorter_class2][TEsorter_class3]={}      
                            dict_Sat1_type_num[TEsorter_class1][TEsorter_class2][TEsorter_class3]['in']=0
                            dict_Sat1_type_num[TEsorter_class1][TEsorter_class2][TEsorter_class3]['out']=0
                        if VSat1_mark=='yes':      
                            dict_Sat1_type_num[TEsorter_class1][TEsorter_class2][TEsorter_class3]['in']+=1    
                            with open(f'./chr2EDTA/4_TEsorter_stat/0_VSat1_innerTE','a') as ff:
                                ff.write(f"{samplechr}\t{line.strip()}\n")
                        else:                           dict_Sat1_type_num[TEsorter_class1][TEsorter_class2][TEsorter_class3]['out']+=1        
                            
                    ## Count near CoreVsat1
                    if samplechr in dict_samplechr_bigblock_expand:
                        core_start,core_end,bigblock_len=dict_samplechr_bigblock_expand[samplechr]
                        if bigblock_len>100000:
                            coreVSat1_mark=''
                            if start>core_start and end<core_end:           coreVSat1_mark='yes'
                            ###
                            ###############
                            if TEsorter_class1 not in dict_CoreVSat1_type_num:
                                dict_CoreVSat1_type_num[TEsorter_class1]={}
                            if TEsorter_class2 not in dict_CoreVSat1_type_num[TEsorter_class1]:
                                dict_CoreVSat1_type_num[TEsorter_class1][TEsorter_class2]={}
                            if TEsorter_class3 not in dict_CoreVSat1_type_num[TEsorter_class1][TEsorter_class2]:
                                dict_CoreVSat1_type_num[TEsorter_class1][TEsorter_class2][TEsorter_class3]={}      
                                dict_CoreVSat1_type_num[TEsorter_class1][TEsorter_class2][TEsorter_class3]['in']=0
                                dict_CoreVSat1_type_num[TEsorter_class1][TEsorter_class2][TEsorter_class3]['out']=0
                            if coreVSat1_mark=='yes':      dict_CoreVSat1_type_num[TEsorter_class1][TEsorter_class2][TEsorter_class3]['in']+=1    
                            else:                           dict_CoreVSat1_type_num[TEsorter_class1][TEsorter_class2][TEsorter_class3]['out']+=1                           
                    ## Only keep those with interarray
                    #if samplechr in dict_samplechr_interarray:
                        #if dict_samplechr_interarray[samplechr]>150000 and dict_samplechr_interarray[samplechr]<350000:
                            ###
                    if samplechr in    dict_samplechr_interarrays:     
                        interarrays=    dict_samplechr_interarrays[samplechr]
                        interarray_mark=''
                        for one_inter_array in interarrays:
                            inter_array_start,inter_array_end=one_inter_array
                            if start>inter_array_start and end<inter_array_end: interarray_mark='yes';break
                        ###
                        ###############
                        if TEsorter_class1 not in dict_LIR_type_num:
                            dict_LIR_type_num[TEsorter_class1]={}
                        if TEsorter_class2 not in dict_LIR_type_num[TEsorter_class1]:
                            dict_LIR_type_num[TEsorter_class1][TEsorter_class2]={}
                        if TEsorter_class3 not in dict_LIR_type_num[TEsorter_class1][TEsorter_class2]:
                            dict_LIR_type_num[TEsorter_class1][TEsorter_class2][TEsorter_class3]={}      
                            dict_LIR_type_num[TEsorter_class1][TEsorter_class2][TEsorter_class3]['in']=0
                            dict_LIR_type_num[TEsorter_class1][TEsorter_class2][TEsorter_class3]['out']=0
                        if interarray_mark=='yes':      dict_LIR_type_num[TEsorter_class1][TEsorter_class2][TEsorter_class3]['in']+=1    
                        else:                           dict_LIR_type_num[TEsorter_class1][TEsorter_class2][TEsorter_class3]['out']+=1    
                        
                    ### Count blocks    

                    if samplechr in    dict_samplechr_type_allblocks: 

                        dict_type_list=dict_samplechr_type_allblocks[samplechr]
                        centype_list=['cen_107','cen_66','cen_103','cen_191','cen_355','cen_383']
                        for one_centype in centype_list:
                            if one_centype in dict_type_list:         ###### This causes if a certain sample's chr doesn't have cen_383 detected, the following statistics won't run. Taking Athila as an example, the results in Figure 5a all consider chromosomes that have that type
                                array=dict_type_list[one_centype]  
                                array_mark='';mark_type=''
                                for one in array:
                                    one_start,one_end,one_strand=one
                                    #region1= one_start-40000;      region2=one_start-20000;        region3=one_end+20000;        region4=one_end+40000;
                                    region1= one_start-60000;      region2=one_start-30000;        region3=one_end+30000;        region4=one_end+60000;
                                    #region1= one_start-80000;      region2=one_start-40000;        region3=one_end+40000;        region4=one_end+80000;
                                    if one_strand=='+':
                                        if start<region1 and end>region1: array_mark='yes';mark_type='part1';break
                                        elif start>region1 and end<region2: array_mark='yes';mark_type='part1';break
                                        ##
                                        elif start<region2 and end>region2: array_mark='yes';mark_type='part2';break
                                        elif start>region2 and end<one_start: array_mark='yes';mark_type='part2';break
                                        elif start<one_start and end>one_start: array_mark='yes';mark_type='part2';break
                                        ##
                                        elif start>one_start and end<one_end: array_mark='yes';mark_type='part3';break
                                        ##
                                        elif start<one_end and end>one_end: array_mark='yes';mark_type='part4';break
                                        elif start>one_end and end<region3: array_mark='yes';mark_type='part4';break
                                        elif start<region3 and end>region3: array_mark='yes';mark_type='part4';break  
                                        ##
                                        elif start>region3 and end<region4: array_mark='yes';mark_type='part5';break
                                        elif start<region4 and end>region4: array_mark='yes';mark_type='part5';break                                    
                                    else:
                                        if start<region1 and end>region1: array_mark='yes';mark_type='part5';break
                                        elif start>region1 and end<region2: array_mark='yes';mark_type='part5';break
                                        ##
                                        elif start<region2 and end>region2: array_mark='yes';mark_type='part4';break
                                        elif start>region2 and end<one_start: array_mark='yes';mark_type='part4';break
                                        elif start<one_start and end>one_start: array_mark='yes';mark_type='part4';break
                                        ##
                                        elif start>one_start and end<one_end: array_mark='yes';mark_type='part3';break
                                        ##
                                        elif start<one_end and end>one_end: array_mark='yes';mark_type='part2';break
                                        elif start>one_end and end<region3: array_mark='yes';mark_type='part2';break
                                        elif start<region3 and end>region3: array_mark='yes';mark_type='part2';break  
                                        ##
                                        elif start>region3 and end<region4: array_mark='yes';mark_type='part1';break
                                        elif start<region4 and end>region4: array_mark='yes';mark_type='part1';break 
                                        ##
                            
                                #########
                                if    one_centype not in  dict_centype_type_num:
                                    dict_centype_type_num[one_centype]={}
                                if TEsorter_class3 not in dict_centype_type_num[one_centype]:
                                    dict_centype_type_num[one_centype][TEsorter_class3]={} 
                                
                                    dict_centype_type_num[one_centype][TEsorter_class3]['part1']=0
                                    dict_centype_type_num[one_centype][TEsorter_class3]['part2']=0
                                    dict_centype_type_num[one_centype][TEsorter_class3]['part3']=0
                                    dict_centype_type_num[one_centype][TEsorter_class3]['part4']=0
                                    dict_centype_type_num[one_centype][TEsorter_class3]['part5']=0
                                    dict_centype_type_num[one_centype][TEsorter_class3]['part_out']=0
                                if array_mark=='yes':      
                                    if mark_type=='part1':          dict_centype_type_num[one_centype][TEsorter_class3]['part1']+=1   
                                    if mark_type=='part2':          dict_centype_type_num[one_centype][TEsorter_class3]['part2']+=1   
                                    if mark_type=='part3':          dict_centype_type_num[one_centype][TEsorter_class3]['part3']+=1   
                                    if mark_type=='part4':          dict_centype_type_num[one_centype][TEsorter_class3]['part4']+=1   
                                    if mark_type=='part5':          dict_centype_type_num[one_centype][TEsorter_class3]['part5']+=1                                   
                                else:                               dict_centype_type_num[one_centype][TEsorter_class3]['part_out']+=1  
                            
                        
                        
            print('Output')
            with open(f'./chr2EDTA/4_TEsorter_stat/0_VSat1','a') as f2:    
                #with open(f'./chr2EDTA/4_TEsorter_stat/{one_sample}.coreVSat1_stat','w') as f:
                    #f.write(f"sample\tsample_type\tClass\tClass1\tClass2\tClass3\tinter_array_num\tout_array_num\n")
                for TEsorter_class1,dict_Class2_info in dict_Sat1_type_num.items():
                    for TEsorter_class2,dict_Class3_info in dict_Class2_info.items():
                        for TEsorter_class3,dict_info   in dict_Class3_info.items():
                            in_num=     dict_info['in']
                            out_num=    dict_info['out']
                            sample_type=dict_sample_type[one_sample]
                            f2.write(f"{one_sample}\t{sample_type}\t{TEsorter_class1}|{TEsorter_class2}|{TEsorter_class3}\t{TEsorter_class1}\t{TEsorter_class2}\t{TEsorter_class3}\t{in_num}\t{out_num}\n")
                            #f.write(f"{one_sample}\t{sample_type}\t{TEsorter_class1}|{TEsorter_class2}|{TEsorter_class3}\t{TEsorter_class1}\t{TEsorter_class2}\t{TEsorter_class3}\t{in_num}\t{out_num}\n")
                           
            with open(f'./chr2EDTA/4_TEsorter_stat/0_coreVSat1','a') as f2:    
                #with open(f'./chr2EDTA/4_TEsorter_stat/{one_sample}.coreVSat1_stat','w') as f:
                    #f.write(f"sample\tsample_type\tClass\tClass1\tClass2\tClass3\tinter_array_num\tout_array_num\n")
                for TEsorter_class1,dict_Class2_info in dict_CoreVSat1_type_num.items():
                    for TEsorter_class2,dict_Class3_info in dict_Class2_info.items():
                        for TEsorter_class3,dict_info   in dict_Class3_info.items():
                            in_num=     dict_info['in']
                            out_num=    dict_info['out']
                            sample_type=dict_sample_type[one_sample]
                            #
                            allin=dict_sample_CoreVSat1_inout[one_sample]['in']/1000000
                            alllength=dict_sample_CoreVSat1_inout[one_sample]['all']/1000000
                            allout=alllength-allin
                            in_density=round(in_num/allin,4)
                            out_density=round(out_num/allout,4)
                            delta_density=in_density-out_density
                            #
                            f2.write(f"{one_sample}\t{sample_type}\t{TEsorter_class1}|{TEsorter_class2}|{TEsorter_class3}\t{TEsorter_class1}\t{TEsorter_class2}\t{TEsorter_class3}\t{in_num}\t{out_num}\t{allin:.2f}\t{allout:.2f}\t{in_density}\t{out_density}\t{delta_density}\n")
                            #f.write(f"{one_sample}\t{sample_type}\t{TEsorter_class1}|{TEsorter_class2}|{TEsorter_class3}\t{TEsorter_class1}\t{TEsorter_class2}\t{TEsorter_class3}\t{in_num}\t{out_num}\n")
                            
            with open(f'./chr2EDTA/4_TEsorter_stat/0_interarray','a') as f2:   
                #with open(f'./chr2EDTA/4_TEsorter_stat/{one_sample}.interarray_stat','w') as f:
                    #f.write(f"sample\tsample_type\tClass\tClass1\tClass2\tClass3\tinter_array_num\tout_array_num\n")
                for TEsorter_class1,dict_Class2_info in dict_LIR_type_num.items():
                    for TEsorter_class2,dict_Class3_info in dict_Class2_info.items():
                        for TEsorter_class3,dict_info   in dict_Class3_info.items():
                            in_num=     dict_info['in']
                            out_num=    dict_info['out']
                            sample_type=dict_sample_type[one_sample]
                            #
                            allin=dict_sample_LIR_inout[one_sample]['in']/1000000
                            alllength=dict_sample_LIR_inout[one_sample]['all']/1000000
                            allout=alllength-allin
                            in_density=round(in_num/allin,4)
                            out_density=round(out_num/allout,4)
                            delta_density=in_density-out_density
                            #
                            f2.write(f"{one_sample}\t{sample_type}\t{TEsorter_class1}|{TEsorter_class2}|{TEsorter_class3}\t{TEsorter_class1}\t{TEsorter_class2}\t{TEsorter_class3}\t{in_num}\t{out_num}\t{allin:.2f}\t{allout:.2f}\t{in_density}\t{out_density}\t{delta_density}\n")
                                #f.write(f"{one_sample}\t{sample_type}\t{TEsorter_class1}|{TEsorter_class2}|{TEsorter_class3}\t{TEsorter_class1}\t{TEsorter_class2}\t{TEsorter_class3}\t{in_num}\t{out_num}\n")
            with open(f'./chr2EDTA/4_TEsorter_stat/0_cen_stat','a') as f2:   
                #with open(f'./chr2EDTA/4_TEsorter_stat/{one_sample}.interarray_stat','w') as f:
                    #f.write(f"sample\tsample_type\tClass\tClass1\tClass2\tClass3\tinter_array_num\tout_array_num\n")
                for one_centype,dict_type_num in dict_centype_type_num.items():
                    for TEsorter_class3,dict_info   in dict_type_num.items():                                
                        part1_num=     int(dict_info['part1'])
                        part2_num=     int(dict_info['part2'])
                        part3_num=     int(dict_info['part3'])
                        part4_num=     int(dict_info['part4'])
                        part5_num=     int(dict_info['part5'])
                        part_out_num=    int(dict_info['part_out'] )
                        part_all_num=part1_num+part2_num+part3_num+part4_num+part5_num+part_out_num
                        ###
                        part1_percent= round(part1_num/part_all_num*100,3)
                        part2_percent= round(part2_num/part_all_num*100,3)
                        part3_percent= round(part3_num/part_all_num*100,3)
                        part4_percent= round(part4_num/part_all_num*100,3)
                        part5_percent= round(part5_num/part_all_num*100,3)
                        part_out_percent= round(part_out_num/part_all_num*100,3)
                        f2.write(f"{one_centype}\t{one_sample}\t{sample_type}\t{TEsorter_class1}|{TEsorter_class2}|{TEsorter_class3}\t{TEsorter_class1}\t{TEsorter_class2}\t{TEsorter_class3}\t{part1_num}\t{part2_num}\t{part3_num}\t{part4_num}\t{part5_num}\t{part_out_num}\t{part_all_num}\t{part1_percent}\t{part2_percent}\t{part3_percent}\t{part4_percent}\t{part5_percent}\t{part_out_percent}\n")
                    
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, sample_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(sample_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()          
        dict_centype_TEtype_info={}        
        with open(f'./chr2EDTA/4_TEsorter_stat/0_cen_stat','r') as f: 
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                one_centype=eachline_arr[0]
                TEsorter_class3=eachline_arr[6]
                if one_centype not in dict_centype_TEtype_info:
                    dict_centype_TEtype_info[one_centype]={}
                if TEsorter_class3 not in    dict_centype_TEtype_info[one_centype]:
                    dict_centype_TEtype_info[one_centype][TEsorter_class3]={}
                    dict_centype_TEtype_info[one_centype][TEsorter_class3]['part1'] =0
                    dict_centype_TEtype_info[one_centype][TEsorter_class3]['part2'] =0
                    dict_centype_TEtype_info[one_centype][TEsorter_class3]['part3'] =0
                    dict_centype_TEtype_info[one_centype][TEsorter_class3]['part4'] =0
                    dict_centype_TEtype_info[one_centype][TEsorter_class3]['part5'] =0
                    dict_centype_TEtype_info[one_centype][TEsorter_class3]['part_all'] =0
                dict_centype_TEtype_info[one_centype][TEsorter_class3]['part1']+=int(eachline_arr[7])
                dict_centype_TEtype_info[one_centype][TEsorter_class3]['part2']+=int(eachline_arr[8])
                dict_centype_TEtype_info[one_centype][TEsorter_class3]['part3']+=int(eachline_arr[9])
                dict_centype_TEtype_info[one_centype][TEsorter_class3]['part4']+=int(eachline_arr[10])
                dict_centype_TEtype_info[one_centype][TEsorter_class3]['part5']+=int(eachline_arr[11])
                dict_centype_TEtype_info[one_centype][TEsorter_class3]['part_all']+=int(eachline_arr[13])
        with open(f'./chr2EDTA/4_TEsorter_stat/0_cen_stat_sum','w') as f:   
            f.write(f"centype\tTEtype\tpart\tnum\tpercent_by_TE\tpercent_by_statelliteblocknum\n")
            for one_centype,dict_TEtype_info in dict_centype_TEtype_info.items():
                part_all_sum=0
                for TEtype, dict_info in dict_TEtype_info.items():
                    part_all=int(dict_info['part_all'])
                    part_all_sum+=part_all
                for TEtype, dict_info in dict_TEtype_info.items():
                    part1=int(dict_info['part1'])
                    part2=int(dict_info['part2'])
                    part3=int(dict_info['part3'])
                    part4=int(dict_info['part4'])
                    part5=int(dict_info['part5'])
                    
                    part1_percent=round(part1/part_all_sum*100,5)
                    part2_percent=round(part2/part_all_sum*100,5)
                    part3_percent=round(part3/part_all_sum*100,5)
                    part4_percent=round(part4/part_all_sum*100,5)
                    part5_percent=round(part5/part_all_sum*100,5)
                    ### part1_percent is divided by the total number of TEs; I might need to divide by the total number of block events
                    blocks_checked_num=dict_blocktype_num[one_centype]
                    part1_percent_byblock=round(part1/blocks_checked_num*100,5)
                    part2_percent_byblock=round(part2/blocks_checked_num*100,5)
                    part3_percent_byblock=round(part3/blocks_checked_num*100,5)
                    part4_percent_byblock=round(part4/blocks_checked_num*100,5)
                    part5_percent_byblock=round(part5/blocks_checked_num*100,5)
                    ##########
                    f.write(f"{one_centype}\t{TEtype}\tP1\t{part1}\t{part1_percent}\t{part1_percent_byblock}\n")
                    f.write(f"{one_centype}\t{TEtype}\tP2\t{part2}\t{part2_percent}\t{part2_percent_byblock}\n")
                    #f.write(f"{one_centype}\t{TEtype}\t3\t{part3}\n")
                    f.write(f"{one_centype}\t{TEtype}\tP4\t{part4}\t{part4_percent}\t{part4_percent_byblock}\n")
                    f.write(f"{one_centype}\t{TEtype}\tP5\t{part5}\t{part5_percent}\t{part5_percent_byblock}\n")
                    
                    
                    
    if argv1=="stepall"  or argv1=="step4" or argv1=="step4.0" or argv1=="step4.0p":                                
        print('Plotting')            
        R_txt=f'''library(ggplot2)
library(dplyr)
# Set working directory
setwd('./')

# Read data
input_file1 <- read.table('0_VSat1', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class3, y =in_num/(in_num+out_num), color = sample_type),position = position_jitter(width = 0.3)) 
p <- p + theme_bw()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_VSat1.pdf", width = 20 / 2.54, height = 50 / 2.54)
print(p)
dev.off()


# Read data
input_file1 <- read.table('0_coreVSat1', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class3, y =in_num/(in_num+out_num), color = sample_type),position = position_jitter(width = 0.3)) 
p <- p + theme_bw()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_corVSat1.pdf", width = 20 / 2.54, height = 50 / 2.54)
print(p)
dev.off()

# Read data
input_file1 <- read.table('0_coreVSat1', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class3, y =delta_density, color = delta_density),position = position_jitter(width = 0.3)) 
p <- p + theme_bw()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_corVSat1_density.pdf", width = 20 / 2.54, height = 50 / 2.54)
print(p)
dev.off()

# Read data
input_file1 <- read.table('0_coreVSat1', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class3, y =delta_density), color = 'darkblue',size=0.1,position = position_jitter(width = 0.3)) 
p <- p + theme_classic()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_corVSat1_density_simple.pdf", width = 20 / 2.54, height = 8 / 2.54)
print(p)
dev.off()


# Read data
input_file1 <- read.table('0_coreVSat1', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class3, y =delta_density, color = sample_type),position = position_jitter(width = 0.3)) 
p <- p + theme_bw()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_corVSat1_density_sample_type.pdf", width = 20 / 2.54, height = 50 / 2.54)
print(p)
dev.off()

# Read data
input_file1 <- read.table('0_interarray', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class3, y =in_num/(in_num+out_num), color = sample_type),position = position_jitter(width = 0.3)) 
p <- p + theme_bw()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_interarray.pdf", width = 20 / 2.54, height = 50 / 2.54)
print(p)
dev.off()

# Read data
input_file1 <- read.table('0_interarray', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class3, y =delta_density, color = delta_density),position = position_jitter(width = 0.3)) 
p <- p + theme_bw()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_interarray_density.pdf", width = 20 / 2.54, height = 50 / 2.54)
print(p)
dev.off()

# Read data
input_file1 <- read.table('0_interarray', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class3, y =delta_density), color = 'darkblue',size=0.1,position = position_jitter(width = 0.3)) 
p <- p + theme_classic()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_interarray_density_simple.pdf", width = 20 / 2.54, height = 8 / 2.54)
print(p)
dev.off()

# Read data
input_file1 <- read.table('0_interarray', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class3, y =delta_density, color = sample_type),position = position_jitter(width = 0.3)) 
p <- p + theme_bw()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_interarray_density_sample_type.pdf", width = 20 / 2.54, height = 50 / 2.54)
print(p)
dev.off()
'''
        with open(f'./chr2EDTA/4_TEsorter_stat//plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = f"./chr2EDTA/4_TEsorter_stat/"
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')                                       
    if argv1=="stepall"  or argv1=="step4" or argv1=="step4.0" or argv1=="step4.0p2":                                
        print('Plotting')            
        R_txt=f'''library(ggplot2)
library(dplyr)
# Set working directory
setwd('./')

# Read data
input_file1 <- read.table('0_cen_stat_sum', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')
input_file1$centype <- factor(input_file1$centype, levels = c("cen_107", "cen_66", "cen_103", "cen_191", "cen_355", "cen_383"))
input_file1$TEtype <- factor(input_file1$TEtype, levels = c("Ale","Alesia","Angela","Bianca","Ikeros","Ivana","SIRE","TAR","Tork","Athila","CRM","Galadriel","Ogre","Reina","Retand","Tekay"))

color_values=c(
    "Ale"="#ffca33",
    "Alesia"="#a57a00",
    "Angela"="#a57a00",
    "Bianca"="#a57a00",
    "Ikeros"="blue",
    "Ivana"="#a57a00",
    "SIRE"="#a57a00",
    "TAR"="#a57a00",
    "Athila"="#00aac1",
    "CRM"="#891555",
    "Galadriel"="#ff3300",
    "Ogre"="#666699",
    "Reina"="#ff9999",
    "Retand"="#00cc66",
    "Tekay"="#c900c4"
    )
# Create plot
p <- ggplot()
p <-p+ geom_bar(data=input_file1,aes(x=part,y=percent_by_TE,fill=TEtype),stat = "identity")
p <- p + theme_classic()+scale_fill_manual(values = color_values, drop = FALSE)
p <- p + facet_wrap(~ centype, ncol = 10, scales = "free_y")
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_VSats_sum_by_TE.pdf", width = 20 / 2.54, height = 10 / 2.54)
print(p)
dev.off()
# Create plot
p <- ggplot()
p <-p+ geom_bar(data=input_file1,aes(x=part,y=percent_by_statelliteblocknum,fill=TEtype),stat = "identity")
p <- p + theme_classic()+scale_fill_manual(values = color_values, drop = FALSE)
p <- p + facet_wrap(~ centype, ncol = 10, scales = "free_y")
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_VSats_sum_by_statelliteblocknum.pdf", width = 20 / 2.54, height = 10 / 2.54)
print(p)
dev.off()

'''
        with open(f'./chr2EDTA/4_TEsorter_stat//plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = f"./chr2EDTA/4_TEsorter_stat/"
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')                                       
    
    if argv1=="stepall"  or argv1=="step4" or argv1=="step4.0_backgroud_stat":
        print('Statistics on insertion probability in regions outside VSat1-6. Not very keen to analyze; feels like the first row of Figure 5 is a bit crowded')
        
        centype_list=['cen_107','cen_66','cen_103','cen_191','cen_355','cen_383']            
        
        print('Loading a sample information table ../samples_satellite/sample_info')
        dict_sample_type={}
        with open('../samples_satellite/sample_info') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=6:continue
                sample_name=eachline_arr[0]
                sample_type=eachline_arr[3]
                if sample_type not in ["Eurasian","East_Asia","America"]:continue
                dict_sample_type[sample_name]=sample_type      
        sample_list=list(dict_sample_type.keys())
        sample_list.sort()
        
        
        print('Loading a sample information table /home/lain/aaa_data/run0/stat_plot/0-region2info') 
        dict_samplechr_allblocks={}
        with open('/home/lain/aaa_data/run0/stat_plot/0-region2info','r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=9:continue
                sample_name=eachline_arr[0]
                chromosome=eachline_arr[2]
                centype=eachline_arr[3]
                strand=eachline_arr[7]
                samplechr=sample_name+'___'+chromosome
                start=          int(eachline_arr[4])
                end=            int(eachline_arr[5])
                block_length= int(eachline_arr[6])
                ########
                if block_length>1000:
                    if samplechr not in dict_samplechr_allblocks:dict_samplechr_allblocks[samplechr]=[]
                    if centype not in centype_list:continue
                    dict_samplechr_allblocks[samplechr].append((start,end))   
        
            
        def merge_intervals(intervals):
            if not intervals:
                return []
            sorted_intervals = sorted(intervals, key=lambda x: x[0])
            merged = [sorted_intervals[0]]
            for current in sorted_intervals[1:]:
                last = merged[-1]
                if current[0] <= last[1]:  # Overlapping or adjacent
                    new_start = last[0]
                    new_end = max(last[1], current[1])
                    merged[-1] = (new_start, new_end)
                else:
                    merged.append(current)
            return merged            
        
        dict_samplechr_mergedparts={}  
        all_parts_length=0
        for samplechr,allblocks in dict_samplechr_allblocks.items():
            sample=samplechr.split('___')[0]
            if sample not in sample_list:continue
            merged_parts = merge_intervals(allblocks)
            dict_samplechr_mergedparts[samplechr]=merged_parts
            for one_part in merged_parts:
                one_part_start,one_part_end=one_part
                all_parts_length+=one_part_end-one_part_start+1
 
        

        
        

        
        def run_step(one_sample):
            dict_type_num={}
            TE_type_list=["Ale","Alesia","Angela","Bianca","Ikeros","Ivana","SIRE","TAR","Tork","Athila","CRM","Galadriel","Ogre","Reina","Retand","Tekay"]
            for one_TE_type in TE_type_list:
                dict_type_num[one_TE_type]=0            
            with open (f'./chr2EDTA/2_TEsorter/{one_sample}/EDTA.fa.rexdb-plant.cls.tsv','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    TEsorter_class1=eachline_arr[1]
                    TEsorter_class2=eachline_arr[2]
                    TEsorter_class3=eachline_arr[3]
                    if TEsorter_class3 not in TE_type_list:continue
                    Complete=       eachline_arr[4]
                    Strand=         eachline_arr[5]	
                    Domains=        eachline_arr[6]
                    if Complete!='yes':continue
                    ######
                    #print(eachline_arr[0])
                    pos_type_arr=eachline_arr[0].split("_Chr")[1].split('#')
                    pos_arr=pos_type_arr[0].split(":")
                    chromosome='Chr'+pos_arr[0]
                    pos_arr2=pos_arr[1].split('..')
                    start=int(pos_arr2[0])
                    end=int(pos_arr2[1])
                    type_arr=pos_type_arr[1].split('/')
                    EDTA_class1=type_arr[0]
                    EDTA_class2=type_arr[1]   
                    samplechr=one_sample+'___'+chromosome
                    satellite_out_mark='yes'
                    mergedparts=dict_samplechr_mergedparts[samplechr]
                    #print(mergedparts)
                    for one_part in mergedparts:
                        one_part_start      ,     one_part_end   =one_part
                        if start>one_part_end or end<one_part_start:continue
                        else: satellite_out_mark='';break
                    if satellite_out_mark=='yes':
                        dict_type_num[TEsorter_class3]+=1
            return(dict_type_num)
        result_list=[]            
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, sample_list), start=1):
                result_list.append(result)
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(sample_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()  
        #print(result_list)   
        dict_type_num={}
        TE_type_list=["Ale","Alesia","Angela","Bianca","Ikeros","Ivana","SIRE","TAR","Tork","Athila","CRM","Galadriel","Ogre","Reina","Retand","Tekay"]
        for one_TE_type in TE_type_list:
            dict_type_num[one_TE_type]=0   
             
        for one_result in result_list:
            for one_type,num in one_result.items():
                dict_type_num[one_type]+=num
            
            
        all_length=0    
        for one_sample in sample_list:
            with open(f'./chr2EDTA/0_prepare/{one_sample}/{one_sample}.fasta.fai','r') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    chromosome=eachline_arr[0]
                    length=int(eachline_arr[1])
                    all_length+=length
        
        satellite_out_length=all_length-all_parts_length
        with open(f'./chr2EDTA/4_TEsorter_stat/0_backgroud_stat','w') as f:
            f.write('The analysis feels off. My Figure 5a is based on the probability of a TE landing in a given region vs outside, divided by the total number of TEs')
            f.write('step4.0_backgroud_stat is based on chromosome-level analysis, i.e., the probability of having this TE in any random 30kbp region')
            f.write(f"satellite_out_length\t{satellite_out_length}\n")
            target_unit_num=round(satellite_out_length/30000,0)
            for one_type,num in dict_type_num.items():
                f.write(f"{one_type}\t{num}\t{num/target_unit_num:.5f}\n")
                
        
        
        
        
    if argv1=="stepall"  or argv1=="step4" or argv1=="step4.1" :
        print('Only count chromosomes containing Interarray, then calculate into the core region.')
        if  os.path.exists('./chr2EDTA/4_TEsorter_stat/1_combine/')!=True:
            subprocess.run(["mkdir ./chr2EDTA/4_TEsorter_stat/1_combine/"], shell=True)    
                   
        dict_class_interarraynum={}            
        with open(f'./chr2EDTA/4_TEsorter_stat/0_interarray',"r") as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                sample_class=eachline_arr[0]+'___'+eachline_arr[2]
                dict_class_interarraynum[sample_class]=eachline_arr[6]
        with open(f'./chr2EDTA/4_TEsorter_stat/1_combine/0_coreVSat1_interarray',"w") as f2:     
            f2.write(f"sample\tsample_type\tClass\tClass1\tClass2\tClass3\tin_num\tout_num\tin_region\tout_region\tin_density\tout_density\tdelta_density\tinterarray_num\tinterarry_percent\n")            
            with open(f'./chr2EDTA/4_TEsorter_stat/0_coreVSat1',"r") as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    sample_class=eachline_arr[0]+'___'+eachline_arr[2]
                    if sample_class not in dict_class_interarraynum:continue
                    interarray_num=int(dict_class_interarraynum[sample_class])
                    in_num=int(eachline_arr[6])
                    if in_num==0:
                        interarry_percent=0
                    else:    
                        interarry_percent=  round(interarray_num/ in_num,4)
                    f2.write(f"{eachline}\t{interarray_num}\t{interarry_percent}\n")    
                    
    if argv1=="stepall"  or argv1=="step4" or argv1=="step4.1" or argv1=="step4.1p":                                
        print('Plotting')            
        R_txt=f'''library(ggplot2)
library(dplyr)
# Set working directory
setwd('./')

# Read data
input_file1 <- read.table('0_coreVSat1_interarray', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')

# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class3, y =delta_density, color = interarry_percent),position = position_jitter(width = 0.3)) 
p <- p + theme_bw()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_interarray_density.pdf", width = 20 / 2.54, height = 50 / 2.54)
print(p)
dev.off()
'''
        with open(f'./chr2EDTA/4_TEsorter_stat/1_combine/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = f"./chr2EDTA/4_TEsorter_stat/1_combine"
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../../')                                                  
            

    if argv1=="stepall"  or argv1=="step4" or argv1=="step4.2":
        print('Plot interarray and VSat1Core together')
        if  os.path.exists('./chr2EDTA/4_TEsorter_stat/2_combine/')!=True:
            subprocess.run(["mkdir ./chr2EDTA/4_TEsorter_stat/2_combine/"], shell=True)    
                   
        with open(f'./chr2EDTA/4_TEsorter_stat/2_combine/0_coreVSat1_interarray',"w") as f2:     
            f2.write(f"sample\tsample_type\tClass\tClass1\tClass2\tClass3\tin_num\tout_num\tin_region\tout_region\tin_density\tout_density\tdelta_density\tgroup\n")                 
            with open(f'./chr2EDTA/4_TEsorter_stat/0_interarray',"r") as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    f2.write(f"{eachline}\tinterarray\n")
            with open(f'./chr2EDTA/4_TEsorter_stat/0_coreVSat1',"r") as f:
                next(f)
                for line in f:
                    eachline=line.strip()   
                    f2.write(f"{eachline}\tcoreVSat1\n")
                    
    if argv1=="stepall"  or argv1=="step4" or argv1=="step4.2" or argv1=="step4.2p":                                
        print('Plotting')            
        R_txt=f'''library(ggplot2)
library(dplyr)
# Set working directory
setwd('./')

# Read data
input_file1 <- read.table('0_coreVSat1_interarray', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')
input_file1$Class3 <- factor(input_file1$Class3, levels = c("Ale","Alesia","Angela","Bianca","Ikeros","Ivana","SIRE","TAR","Tork","Athila","CRM","Galadriel","Ogre","Reina","Retand","Tekay"))


# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class3, y =delta_density, color = group,group=group),
    position = position_jitter(width = 0.3),
    #position = position_jitterdodge(
    #  dodge.width = 0.5,  # dodge width (same as position_dodge)
    #  jitter.width = 0.2, # horizontal jitter range (recommended to be smaller than dodge.width)
    #  jitter.height = 0   # vertical jitter range (0 means no jitter)
    #),
    size=0.1) 
p <- p + facet_wrap(~ group, ncol = 1,scales = "free_y")
p <- p + theme_classic()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_interarray_density.pdf", width = 15 / 2.54, height = 10 / 2.54)
print(p)
dev.off()
'''
        with open(f'./chr2EDTA/4_TEsorter_stat/2_combine/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = f"./chr2EDTA/4_TEsorter_stat/2_combine"
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../../')                                                  
  
    if argv1=="stepall"  or argv1=="step4" or argv1=="step4.3":
        print('Plot interarray and VSat1Core together')
        if  os.path.exists('./chr2EDTA/4_TEsorter_stat/3_combine/')!=True:
            subprocess.run(["mkdir ./chr2EDTA/4_TEsorter_stat/3_combine/"], shell=True)    
        dict_group_sample_type_info={}     
        dict_sample_sampletype={}
        with open(f'./chr2EDTA/4_TEsorter_stat/2_combine/0_coreVSat1_interarray',"r") as f:     
            #f2.write(f"sample\tsample_type\tClass\tClass1\tClass2\tClass3\tin_num\tout_num\tin_region\tout_region\tin_density\tout_density\tdelta_density\tgroup\n")                 
            next(f)
            for line in f:
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                sample,sample_type,Class,Class1,Class2,Class3,in_num,out_num,in_region,out_region,in_density,out_density,delta_density,group=eachline_arr
                
                sample_real=sample.replace('_hap1','').replace('.hap1','').replace('_hap2','').replace('.hap2','')
                dict_sample_sampletype[sample_real]=sample_type
                if group not in dict_group_sample_type_info:
                    dict_group_sample_type_info[group]={}
                if sample_real not in dict_group_sample_type_info[group]:
                    dict_group_sample_type_info[group][sample_real]={}
                if Class not in  dict_group_sample_type_info[group][sample_real]:
                    dict_group_sample_type_info[group][sample_real][Class]={}
                    dict_group_sample_type_info[group][sample_real][Class]['in_num']=0
                    dict_group_sample_type_info[group][sample_real][Class]['out_num']=0
                    dict_group_sample_type_info[group][sample_real][Class]['in_region']=0
                    dict_group_sample_type_info[group][sample_real][Class]['out_region']=0
                dict_group_sample_type_info[group][sample_real][Class]['in_num']+=int(in_num)
                dict_group_sample_type_info[group][sample_real][Class]['out_num']+=int(out_num)
                dict_group_sample_type_info[group][sample_real][Class]['in_region']+=float(in_region)
                dict_group_sample_type_info[group][sample_real][Class]['out_region']+=float(out_region)
        
        with open(f'./chr2EDTA/4_TEsorter_stat/3_combine/0_coreVSat1_interarray',"w") as f:     
            f.write(f"sample\tsample_type\tClass\tClass1\tClass2\tClass3\tin_num\tout_num\tin_region\tout_region\tin_density\tout_density\tdelta_density\tgroup\n")                 
            for group,dict_sample_type_info in dict_group_sample_type_info.items():
                for sample,dict_type_info in dict_sample_type_info.items():
                    sample_type=dict_sample_sampletype[sample]
                    for one_type, dict_info in dict_type_info.items():
                        Class1,Class2,Class3=one_type.split('|')
                        in_num=dict_info['in_num']
                        out_num=dict_info['out_num']
                        in_region=dict_info['in_region']
                        out_region=dict_info['out_region']
                        in_density=round(in_num/in_region,4)
                        out_density=round(out_num/out_region,4)
                        delta_density=in_density-out_density
                        f.write(f"{sample}\t{sample_type}\t{Class}\t{Class1}\t{Class2}\t{Class3}\t{in_num}\t{out_num}\t{in_region}\t{out_region}\t{in_density}\t{out_density}\t{delta_density}\t{group}\n")                 
    if argv1=="stepall"  or argv1=="step4" or argv1=="step4.3" or argv1=="step4.3p":                                
        print('Plotting')            
        R_txt=f'''library(ggplot2)
library(dplyr)
# Set working directory
setwd('./')

# Read data
input_file1 <- read.table('0_coreVSat1_interarray', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,row.names = NULL  , sep = '\\t')
input_file1$Class3 <- factor(input_file1$Class3, levels = c("Ale","Alesia","Angela","Bianca","Ikeros","Ivana","SIRE","TAR","Tork","Athila","CRM","Galadriel","Ogre","Reina","Retand","Tekay"))


# Create plot
p <- ggplot()
p <- p + geom_point(data = input_file1,aes(x = Class3, y =delta_density, color = group,group=group),
    position = position_jitter(width = 0.3),
    #position = position_jitterdodge(
    #  dodge.width = 0.5,  # dodge width (same as position_dodge)
    #  jitter.width = 0.2, # horizontal jitter range (recommended to be smaller than dodge.width)
    #  jitter.height = 0   # vertical jitter range (0 means no jitter)
    #),
    size=0.1) 
p <- p + facet_wrap(~ group, ncol = 1,scales = "free_y")
p <- p + theme_classic()
#p <- p + facet_wrap(~ Class, ncol = 1)
# Add title and axis labels
p <- p + theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  #axis.title = element_blank(),  # Hide axis titles
  #axis.text = element_blank(),   # Hide axis tick text
  #legend.position = "none",
  #axis.ticks = element_blank(),  # Hide axis ticks
  #axis.line = element_blank()    # Hide axis lines
)+
pdf("1_interarray_density.pdf", width = 15 / 2.54, height = 10 / 2.54)
print(p)
dev.off()
'''
        with open(f'./chr2EDTA/4_TEsorter_stat/3_combine/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)
        new_directory = f"./chr2EDTA/4_TEsorter_stat/3_combine"
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../../')
        
if argv1=="stepall" or "step4c" in argv1:
    print('stat, including incomplete ones')
    if  os.path.exists('./chr2EDTA/4_TEsorter_stat3_Chr')!=True:
        subprocess.run(["mkdir ./chr2EDTA/4_TEsorter_stat3_Chr"], shell=True)        
    
     
    TE_type_list=["Ale","Alesia","Angela","Bianca","Ikeros","Ivana","SIRE","TAR","Tork","Athila","CRM","Galadriel","Ogre","Reina","Retand","Tekay"]
    TE_used_list=["Athila","CRM","Galadriel","Ogre","Reina","Retand","Tekay"]
    print('Loading a sample information table ../samples_satellite/2_good_regions_interarray_LIR')
    dict_samplechr_interarrays={}
    dict_chr_interarray_num={}
    with open('../samples_satellite/2_good_regions_interarray_LIR') as f:
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            sample_name=eachline_arr[0]
            chromosome=eachline_arr[1]
            samplechr=sample_name+'___'+chromosome
            ########
            ##########
            if chromosome not in dict_chr_interarray_num:dict_chr_interarray_num[chromosome]=0
            dict_chr_interarray_num[chromosome]+=1
            if eachline_arr[-1]!="NA"    :     ####### Originally [9], changed to -1, i.e., LIR determined by moddotplot
                if samplechr not in dict_samplechr_interarrays:
                    dict_samplechr_interarrays[samplechr]=[]                    
                interarrays=eachline_arr[-1].split('|')
                for one_inter_array in interarrays:
                    pos_arr=one_inter_array.split('-')
                    dict_samplechr_interarrays[samplechr].append([         int(pos_arr[0])       ,     int(pos_arr[1])   ])
                    
    print('Loading a sample information table ../samples_satellite/sample_info')
    dict_sample_type={}
    with open('../samples_satellite/sample_info') as f:
        for line in f:
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=6:continue
            sample_name=eachline_arr[0]
            sample_type=eachline_arr[3]
            #if sample_type not in ["Eurasian","East_Asia","America"]:continue
            dict_sample_type[sample_name]=sample_type      
    sample_list=list(dict_sample_type.keys())
    sample_list.sort()   
    dict_chr_type_sample_num={}

    for one_sample in sample_list:
        with open (f'./chr2EDTA/2_TEsorter/{one_sample}/EDTA.fa.rexdb-plant.cls.tsv','r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                TEsorter_class1=eachline_arr[1]
                TEsorter_class2=eachline_arr[2]
                TEsorter_class3=eachline_arr[3]
                if TEsorter_class3 not in TE_used_list:continue
                #print(TEsorter_class3)
                #if TEsorter_class3 !='Athila' and TEsorter_class3 !='Tekay':continue
                Complete=       eachline_arr[4]
                Strand=         eachline_arr[5]	
                Domains=        eachline_arr[6]
                if Complete!='yes':continue
                ######
                #print(eachline_arr[0])
                pos_type_arr=eachline_arr[0].split("_Chr")[1].split('#')
                pos_arr=pos_type_arr[0].split(":")
                chromosome='Chr'+pos_arr[0]
                pos_arr2=pos_arr[1].split('..')
                start=int(pos_arr2[0])
                end=int(pos_arr2[1])
                type_arr=pos_type_arr[1].split('/')
                EDTA_class1=type_arr[0]
                EDTA_class2=type_arr[1]   
                samplechr=one_sample+'___'+chromosome
                if samplechr in    dict_samplechr_interarrays:     
                    interarrays=    dict_samplechr_interarrays[samplechr]
                    interarray_mark=''
                    for one_inter_array in interarrays:
                        inter_array_start,inter_array_end=one_inter_array
                        if start>inter_array_start and end<inter_array_end: interarray_mark='yes';break
                    ###
                    if interarray_mark=='yes':      
                        if chromosome not in dict_chr_type_sample_num:dict_chr_type_sample_num[chromosome]={}
                        if TEsorter_class3 not in dict_chr_type_sample_num[chromosome]:  dict_chr_type_sample_num[chromosome][TEsorter_class3]={}
                        if one_sample not in dict_chr_type_sample_num[chromosome][TEsorter_class3]:dict_chr_type_sample_num[chromosome][TEsorter_class3][one_sample]=0
                        dict_chr_type_sample_num[chromosome][TEsorter_class3][one_sample]+=1  
    with open('./chr2EDTA/4_TEsorter_stat3_Chr/0_allchr','w') as f: 
        f.write(f"chromosome\tchromosome_num\ttype\tall_num\tnum_onesamplenum1\tLIRnum\taverage_num\tinsert_probability\n")
        for chromosome,dict_type_sample_num in dict_chr_type_sample_num.items():
            chromosome_num=chromosome.replace('Chr','')
            chr_interarray_num=dict_chr_interarray_num[chromosome]
            
            for one_type,dict_sample_num in dict_type_sample_num.items():
                num_over1=0;all_num=0
                for sample,num in dict_sample_num.items():
                    num_over1+=1
                    all_num+=num
                percent_by_LTRnum=round(num_over1/chr_interarray_num,5)
                average_num=round(all_num/chr_interarray_num,5)
                f.write(f"{chromosome}\t{chromosome_num}\t{one_type}\t{all_num}\t{num_over1}\t{chr_interarray_num}\t{average_num}\t{percent_by_LTRnum}\n")
    if 1==1:
        dict_sample_cen6_LIR={}
        dict_sample_cen107_LIR={}
        cen6_LTR_NUM=0
        cen107_LTR_NUM=0
        chr18_LIR='/home/lain/aaa_data/run0/samples_satellite/chr18_LIR'
        with open(chr18_LIR,'r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                sample,start,end,mark=eachline_arr
                if mark=='.':continue
                start=  10000000+float(start)*1000000
                end=    10000000+float(end)*1000000
                if sample not in dict_sample_cen6_LIR:dict_sample_cen6_LIR[sample]=[]
                if sample not in dict_sample_cen107_LIR:dict_sample_cen107_LIR[sample]=[]
                if mark=='cen6_LIR':        dict_sample_cen6_LIR[sample].append([start,end]);cen6_LTR_NUM+=1
                if mark=='cen107_LIR':      dict_sample_cen107_LIR[sample].append([start,end]);cen107_LTR_NUM+=1
        dict_LTRtype_LTRNUM={}
        dict_LTRtype_LTRNUM['cen6_LIR']=cen6_LTR_NUM
        dict_LTRtype_LTRNUM['cen107_LIR']=cen107_LTR_NUM
        ##############
        dict_LTRtype_type_sample_num={}
        dict_LTRtype_type_sample_num['cen6_LIR']={}
        dict_LTRtype_type_sample_num['cen107_LIR']={}
        for one_sample in sample_list:
            with open (f'./chr2EDTA/2_TEsorter/{one_sample}/EDTA.fa.rexdb-plant.cls.tsv','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    TEsorter_class1=eachline_arr[1]
                    TEsorter_class2=eachline_arr[2]
                    TEsorter_class3=eachline_arr[3]
                    #print(TEsorter_class3)
                    #if TEsorter_class3 !='Athila' and TEsorter_class3 !='Tekay':continue
                    Complete=       eachline_arr[4]
                    Strand=         eachline_arr[5]	
                    Domains=        eachline_arr[6]
                    if Complete!='yes':continue
                    ######
                    #print(eachline_arr[0])
                    pos_type_arr=eachline_arr[0].split("_Chr")[1].split('#')
                    pos_arr=pos_type_arr[0].split(":")
                    chromosome='Chr'+pos_arr[0]
                    if chromosome !='Chr18':continue
                    pos_arr2=pos_arr[1].split('..')
                    start=int(pos_arr2[0])
                    end=int(pos_arr2[1])
                    type_arr=pos_type_arr[1].split('/')
                    EDTA_class1=type_arr[0]
                    EDTA_class2=type_arr[1]   
                    samplechr=one_sample+'___'+chromosome
                    #####
                    cen6_mark=''
                    if one_sample in dict_sample_cen6_LIR:
                        cen6_LIR_list=dict_sample_cen6_LIR[one_sample]
                        
                        for one in cen6_LIR_list:
                            one_start,one_end=one
                            if start>=one_start and end<one_end:cen6_mark='yes';break
                    cen107_mark=''    
                    if one_sample in dict_sample_cen6_LIR:    
                        cen107_LIR_list=dict_sample_cen107_LIR[one_sample]    
                        
                        for one in cen107_LIR_list:
                            one_start,one_end=one
                            if start>=one_start and end<one_end:cen107_mark='yes';break
                    if  cen6_mark==''  and    cen107_mark=='':continue
                    if TEsorter_class3 not in dict_LTRtype_type_sample_num['cen6_LIR']:dict_LTRtype_type_sample_num['cen6_LIR'][TEsorter_class3]={}
                    if TEsorter_class3 not in dict_LTRtype_type_sample_num['cen107_LIR']:dict_LTRtype_type_sample_num['cen107_LIR'][TEsorter_class3]={}
                    if one_sample not in dict_LTRtype_type_sample_num['cen6_LIR'][TEsorter_class3]:dict_LTRtype_type_sample_num['cen6_LIR'][TEsorter_class3][one_sample]=0  
                    if one_sample not in dict_LTRtype_type_sample_num['cen107_LIR'][TEsorter_class3]:dict_LTRtype_type_sample_num['cen107_LIR'][TEsorter_class3][one_sample]=0   
                    if cen6_mark=='yes'    :
                        dict_LTRtype_type_sample_num['cen6_LIR'][TEsorter_class3][one_sample]+=1
                    if cen107_mark=='yes'    :
                        dict_LTRtype_type_sample_num['cen107_LIR'][TEsorter_class3][one_sample]+=1
        #print(dict_LTRtype_type_sample_num)
             
        with open('./chr2EDTA/4_TEsorter_stat3_Chr/1_chr18_only','w') as f2: 
            f2.write(f"LTRtype\tsample\tchromosome_num\ttype\tnum\n")
            with open('./chr2EDTA/4_TEsorter_stat3_Chr/1_chr18_only_sum','w') as f: 
                f.write(f"LTRtype\tchromosome_num\ttype\tall_num\tnum_onesamplenum1\tLIRnum\taverage_num\tinsert_probability\n")
                for one_LTRtype,dict_type_sample_num in dict_LTRtype_type_sample_num.items():
                    LTRNUM=dict_LTRtype_LTRNUM[one_LTRtype]
                    if one_LTRtype=='cen6_LIR':chromosome_num=22
                    if one_LTRtype=='cen107_LIR':chromosome_num=21                
                    for one_type,dict_sample_num in dict_type_sample_num.items():
                        if one_type not in TE_used_list:continue
                        num_over1=0;all_num=0
                        for sample,num in dict_sample_num.items():
                            if num>0:
                                num_over1+=1
                                #print(num)
                                all_num+=num
                            #######
                            f2.write(f"{one_LTRtype}\t{sample}\t{chromosome_num}\t{one_type}\t{num}\n")
                        percent_by_LTRnum=round(num_over1/LTRNUM,5)
                        average_num=round(all_num/LTRNUM,5)
                        
                        f.write(f"{one_LTRtype}\t{chromosome_num}\t{one_type}\t{all_num}\t{num_over1}\t{LTRNUM}\t{average_num}\t{percent_by_LTRnum}\n")
                        
    R_txt=f"""
        library(ggplot2)
        library(dplyr)
          input_file=read.table('0_allchr', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
            input_file2=read.table('1_chr18_only_sum', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
            combined_table <- bind_rows(input_file, input_file2)
            
            combined_table$type <- factor(combined_table$type, levels = c("Ale","Alesia","Angela","Bianca","Ikeros","Ivana","SIRE","TAR","Tork","Athila","CRM","Galadriel","Ogre","Reina","Retand","Tekay"))
            
            color_values=c(
                "Ale"="#ffca33",
                "Alesia"="#a57a00",
                "Angela"="#a57a00",
                "Bianca"="#a57a00",
                "Ikeros"="blue",
                "Ivana"="#a57a00",
                "SIRE"="#a57a00",
                "TAR"="green",
                "Athila"="#00aac1",
                "CRM"="#891555",
                "Galadriel"="#ff3300",
                "Ogre"="#666699",
                "Reina"="#ff9999",
                "Retand"="#00cc66",
                "Tekay"="#c900c4"
                )
            
        # Create plot object
        p <- ggplot() +
            geom_bar(data=combined_table,aes(x = chromosome_num, y = insert_probability, fill = type),stat = "identity", position = "stack") +
            theme_classic()+scale_fill_manual(values = color_values, drop = FALSE) #+  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position = "none",axis.text.x = element_blank()) 
         
          
          # Save as PDF
          pdf(file = paste0('1_allchr_onechrfor1', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
          print(p)
          dev.off()
          
        # Create plot object
        p <- ggplot() +
            geom_bar(data=combined_table,aes(x = chromosome_num, y = average_num, fill = type),stat = "identity", position = "stack") +
            theme_classic()+scale_fill_manual(values = color_values, drop = FALSE) #+  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position = "none",axis.text.x = element_blank()) 
         
          
          # Save as PDF
          pdf(file = paste0('1_allchr_average_insertnum', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
          print(p)
          dev.off()           
        """
    with open(f'./chr2EDTA/4_TEsorter_stat3_Chr//plot.R','w',encoding='utf-8') as f:
        f.write(R_txt)  
    new_directory = f'./chr2EDTA/4_TEsorter_stat3_Chr/'
    os.chdir(new_directory)
    subprocess.run(['Rscript plot.R'], shell=True)    
    os.chdir('../../')
                        
if argv1=="stepall" or "step5" in argv1:            
    print('step5')
    if  os.path.exists('./chr2EDTA/5_tree')!=True:
        subprocess.run(["mkdir ./chr2EDTA/5_tree"], shell=True)    

    if argv1=="stepall" or argv1=="step5_readme":
        print('Printing instructions')
        with open('./chr2EDTA/5_tree/readme','w') as f:
            txt=r'''
            A total of over 60,000 athila sequences were detected,
            99.9% are GAG|Athila PROT|Athila RT|Athila RH|Athila INT|Athila
            5.0 - Get info
            5.0s - Add additional columns: whether in VSat1, peri, and interarray
            ##
            step5.0s2a — Use EDTA's judgment to determine LTR sequences, calculate similarity and estimated insertion time
            step5.0s3 — JZX manual analysis of LTR regions
            
            5.1_cdhit deduplication
            5.1s_Process to get cdhit.index
            5.2-mafft167s
            5.3-triaml
            5.3s_Format fasta
            5.4-iqtree
            5.5 - Generate an empty file, output the iTOL tree in Newick format, copy into this file
            5.6 - Generate serial file
            5.7 - Following step5.0s info, add additional columns and sort according to phylogenetic tree
            5.8 - Plot
            
            step5.10 — Select a portion to redo phylogenetic tree ./chr2EDTA/5_tree/{TE_Class}/10_ClassII/
            '''
            f.write(txt) 
            
    task_list=[]      #TE_Class      ###domin                                               #dom sequence name                               #Name used for phylogenetic tree, trimal parameters  -gt0.9 -cons0.6 -resoverlap0.4 
    
    #task_list.append(["Ale","GAG|Ale PROT|Ale INT|Ale RT|Ale RH|Ale",                       "Ale-GAG|Ale-PROT|Ale-RT|Ale-RH|Ale-INT",                       "Ale-RT|Ale-INT",                                                      "-gt 0.6 -cons 0.5 -resoverlap 0.4 -seqoverlap 0.5"])
    #task_list.append(["Ale","GAG|Ale PROT|Ale INT|Ale RT|Ale RH|Ale",                       "Ale-GAG|Ale-PROT|Ale-RT|Ale-RH|Ale-INT",                       "Ale-RT|Ale-INT",             "-automated1"]) 
    task_list.append(["Athila","GAG|Athila PROT|Athila RT|Athila RH|Athila INT|Athila",     "Athila-GAG|Athila-PROT|Athila-RT|Athila-RH|Athila-INT",        "Athila-RT|Athila-INT",       "-automated1"])  
    #task_list.append(["CRM","GAG|CRM PROT|CRM RT|CRM RH|CRM INT|CRM CHD|CRM",               "CRM-GAG|CRM-PROT|CRM-RT|CRM-RH|CRM-INT",        "CRM-RT|CRM-INT|CRM-CHD",       "-automated1"])  
    #task_list.append(["Galadriel","GAG|Galadriel PROT|Galadriel RT|Galadriel RH|Galadriel INT|Galadriel CHD|Galadriel",   "Galadriel-GAG|Galadriel-PROT|Galadriel-RT|Galadriel-RH|Galadriel-INT|Galadriel-CHD", "Galadriel-RT|Galadriel-INT|Galadriel-CHD","-automated1"])  
    #task_list.append(["Ogre","GAG|Ogre PROT|Ogre RT|Ogre RH|Ogre INT|Ogre",                       "Ogre-GAG|Ogre-PROT|Ogre-RT|Ogre-RH|Ogre-INT",        "Ogre-RT|Ogre-INT",       "-automated1"])  
    #task_list.append(["Reina","GAG|Reina PROT|Reina RT|Reina RH|Reina INT|Reina CHD|Reina",       "Reina-GAG|Reina-PROT|Reina-RT|Reina-RH|Reina-INT|Reina-CHD",     "Reina-RT|Reina-INT|Reina-CHD",       "-automated1"])  
    #task_list.append(["Retand","GAG|Retand PROT|Retand RT|Retand RH|Retand INT|Retand",           "Retand-GAG|Retand-PROT|Retand-RT|Retand-RH|Retand-INT", "Retand-RT|Retand-INT",       "-automated1"])  
    task_list.append(["Tekay","GAG|Tekay PROT|Tekay RT|Tekay RH|Tekay INT|Tekay CHD|Tekay", "Tekay-GAG|Tekay-PROT|Tekay-RT|Tekay-RH|Tekay-INT|Tekay-CHD",   "Tekay-RT|Tekay-INT|Tekay-CHD","-automated1"])


    cdhid_c=0.95  ## Default 0.95, cluster at this similarity level

    print(f"Number of TE types: {len(task_list)}")
    print('Loading a sample information table ../samples_satellite/sample_info')
    dict_sample_type={}
    with open('../samples_satellite/sample_info') as f:
        for line in f:
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=6:continue
            sample_name=eachline_arr[0]
            sample_type=eachline_arr[3]
            #if sample_type not in ["Eurasian","East_Asia","America"]:continue
            dict_sample_type[sample_name]=sample_type      
    sample_list=list(dict_sample_type.keys())
    sample_list.sort()  
    sample_list_len=len(sample_list)    
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.0":
        
        i=0
        for one_task in task_list:
            TE_Class,domain,domain_seq_name,tree_used_domain,  trimal_para=one_task
            print("\n\n")
            i+=1
            print(one_task)
            if  os.path.exists(f'./chr2EDTA/5_tree/{TE_Class}')!=True:
                subprocess.run([f"mkdir ./chr2EDTA/5_tree/{TE_Class}"], shell=True)    
            with open(f'./chr2EDTA/5_tree/{TE_Class}/0_info_tmp','w') as f2:     
                f2.write(f"sample\tsample_type\tchromosome\tstart\tend\tstrand\tTE_name\tClass\tClass1\tClass2\tClass3\tComplete\tDomains\tdom_names\tseqs\tseqs_used\n")   
            with open(f'./chr2EDTA/5_tree/{TE_Class}/0_index.fa','w') as f:     
                f.write('')
            ii=0  
            kkk=0
            seq_set=set()
            for one_sample in  sample_list:
                ii+=1
                print(f"{i}/{len(task_list)}-----{ii}/{sample_list_len}",end='\r')
                sample_type=dict_sample_type[one_sample]
                dict_TE_domain_seq={}
                with open (f'./chr2EDTA/2_TEsorter/{one_sample}/EDTA.fa.rexdb-plant.dom.faa','r') as f:
                    for line in f:
                        eachline=line.strip()
                        
                        if eachline[0]=='>':
                            TEname=eachline.split('|')[0][1:]
                            domain_name=eachline.split('|')[-1].split(';')[1][5:]
                            if TEname not in dict_TE_domain_seq:dict_TE_domain_seq[TEname]={}
                            if domain_name not in dict_TE_domain_seq[TEname]:dict_TE_domain_seq[TEname][domain_name]={}
                        else:
                            dict_TE_domain_seq[TEname][domain_name]=eachline
                            TEname='';domain_name=''

                ###############                
                
                with open (f'./chr2EDTA/2_TEsorter/{one_sample}/EDTA.fa.rexdb-plant.cls.tsv','r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        TE_name=        eachline_arr[0]
                        TEsorter_class1=eachline_arr[1]
                        TEsorter_class2=eachline_arr[2]
                        TEsorter_class3=eachline_arr[3]
                        if TEsorter_class3!=TE_Class:continue
                        Complete=       eachline_arr[4]
                        ##Strand=         eachline_arr[5]	this strand is incorrect because everything is from the positive strand of the library
                        Domains=        eachline_arr[6]
                        if Complete!='yes':continue
                        if Domains!=domain:continue
                        ######
                        #print(eachline_arr[0])
                        pos_type_arr=eachline_arr[0].split("_Chr")[1].split('#')
                        pos_arr=pos_type_arr[0].split(":")
                        chromosome='Chr'+pos_arr[0]
                        pos_arr2=pos_arr[1].split('..')
                        start_tmp=int(pos_arr2[0])
                        end_tmp=int(pos_arr2[1])
                        if start_tmp>end_tmp:   start=end_tmp;      end=start_tmp;      Strand='-'
                        else:                   start=start_tmp;    end=end_tmp;        Strand='+'
                        type_arr=pos_type_arr[1].split('/')
                        EDTA_class1=type_arr[0]
                        EDTA_class2=type_arr[1]   
                        samplechr=one_sample+'___'+chromosome
                        ##all_seqs
                        seqs=dict_TE_domain_seq[TE_name]
                        seqs_all=[]
                        for one in domain_seq_name.split("|"):
                            seqs_all.append(seqs[one])
                        seqs_all_str='|'.join(seqs_all)
                        ##used seq
                        seqs_all=[]
                        for one in tree_used_domain.split("|"):
                            seqs_all.append(seqs[one])
                        seqs_used_str='|'.join(seqs_all)    
                        ##
                        if TEsorter_class1!=EDTA_class1:continue
                        
                        with open(f'./chr2EDTA/5_tree/{TE_Class}/0_info_tmp','a') as f2:     
                            f2.write(f"{one_sample}\t{sample_type}\t{chromosome}\t{start}\t{end}\t{Strand}\t{TE_name}\t{TEsorter_class1}|{TEsorter_class2}|{TEsorter_class3}\t{TEsorter_class1}\t{TEsorter_class2}\t{TEsorter_class3}\t{Complete}\t{Domains}\t{domain_seq_name}\t{seqs_all_str}\t{seqs_used_str}\n")
                        ##
                        seq_set.add(seqs_used_str)
            seq_set_list=list(seq_set)
            seq_set_list.sort()
            kk=0
            dict_one_index={}
            with open(f'./chr2EDTA/5_tree/{TE_Class}/0_index.fa','a') as f2:     
                #with open(f'./chr2EDTA/5_tree/{TE_Class}/0_index','a') as f:  
                for one in seq_set_list:
                    kk+=1
                    one_mafftinput=one.replace('|',Interval).replace("*",terminal)
                    #f.write(f"{kk}\t{one}\t{one_mafftinput}\n")
                    f2.write(f">{kk}\n{one_mafftinput}\n")
                    dict_one_index[one]=[kk,one_mafftinput]
            with  open(f'./chr2EDTA/5_tree/{TE_Class}/0_info','w') as f2:     
                f2.write(f"sample\tsample_type\tchromosome\tstart\tend\tstrand\tTE_name\tClass\tClass1\tClass2\tClass3\tComplete\tDomains\tdom_names\tseqs\tseqs_used\tmafftinput\tseq_index\n")   
                with open(f'./chr2EDTA/5_tree/{TE_Class}/0_info_tmp','r') as f:    
                    next(f)
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        seq_used=eachline_arr[15]
                        index,one_mafftinput=dict_one_index[seq_used]
                        f2.write(f"{eachline}\t{one_mafftinput}\t{index}\n")
            subprocess.run([f"rm ./chr2EDTA/5_tree/{TE_Class}/0_info_tmp"], shell=True)                
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.0" or argv1=="step5.0s":  
        print('Determine whether in centromere region, determine whether')
        print('Loading a sample information table ../samples_satellite/2_good_regions_interarray_LIR')
        dict_samplechr_interarrays={}
        dict_samplechr_interarray={}
        dict_samplechr_bigblock_expand={}
        with open('../samples_satellite/2_good_regions_interarray_LIR') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                sample_name=eachline_arr[0]
                chromosome=eachline_arr[1]
                samplechr=sample_name+'___'+chromosome
                inter_array_length= int(eachline_arr[8])
                ########
                dict_samplechr_interarray[samplechr]=inter_array_length
                ##########
                if eachline_arr[-1]!="NA"    :     ######## Originally [9], changed to -1, i.e., LIR determined by moddotplot
                    if samplechr not in dict_samplechr_interarrays:
                        dict_samplechr_interarrays[samplechr]=[]                    
                    interarrays=eachline_arr[-1].split('|')
                    for one_inter_array in interarrays:
                        pos_arr=one_inter_array.split('-')
                        dict_samplechr_interarrays[samplechr].append([         int(pos_arr[0])       ,     int(pos_arr[1])   ])
                ####
                bigblock_start=     int(eachline_arr[2])
                bigblock_end=       int(eachline_arr[3])
                bigblock_len=       int(eachline_arr[4])
                dict_samplechr_bigblock_expand[samplechr]=[bigblock_start-1000000,bigblock_end+1000000,bigblock_len]
       
        if 1==2:   
            print('Using step5.0s——————step5.7——————step5.8, the effect seems not very good, there might be some issues')
            dict_samplechr_bigblock_expand={}
            with open('../samples_satellite/2_good_regions_main','r') as f: #2_good_regions_main _66_103
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    sample,chromosome,start,end,length=eachline_arr
                    samplechr=sample_name+'___'+chromosome
                    dict_samplechr_bigblock_expand[samplechr]=[int(start),int(end),int(length)]
                    
        #print(dict_samplechr_bigblock_expand['ManicureFinger_hap2___Chr11'])   
        print('Loading a sample information table ../samples_satellite/2_good_regions')
        dict_samplechr_allblocks={}
        with open('../samples_satellite/2_good_regions','r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=8:continue
                sample_name=eachline_arr[0]
                chromosome=eachline_arr[1].replace('region_','Chr')
                samplechr=sample_name+'___'+chromosome
                start=          int(eachline_arr[3])
                end=            int(eachline_arr[4])
                block_length= int(eachline_arr[5])
                ########
                if block_length>10000:
                    if samplechr not in dict_samplechr_allblocks:dict_samplechr_allblocks[samplechr]=[]
                    dict_samplechr_allblocks[samplechr].append([         start      ,     end   ])           
        #print(dict_samplechr_allblocks)            
                 
        for one_task in task_list:  
            TE_Class=one_task[0]
            with open(f'./chr2EDTA/5_tree/{TE_Class}/0_info_more','w') as f2:   
                f2.write(f"sample\tsample_type\tchromosome\tstart\tend\tstrand\tTE_name\tClass\tClass1\tClass2\tClass3\tComplete\tDomains\tdom_names\tseqs\tseqs_used\tmafftinput\tseq_index\tin_VSat1\tin_Peri\tin_interarray\n")   
                with open(f'./chr2EDTA/5_tree/{TE_Class}/0_info','r') as f: 
                    next(f)
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        sample,chromosome,chromosome,start,end=eachline_arr[:5]
                        start=int(start)
                        end=int(end)
                        samplechr=sample+'___'+chromosome
                        ### Determine if inside VSat1
                        if samplechr not in dict_samplechr_allblocks:
                            in_VSat1='.'
                        else:
                            mark=''
                            for one in dict_samplechr_allblocks[samplechr]:
                                #print(one[0],start,end,one[1])
                                if one[0]<start and end<one[1]:     mark='yes';break
                                elif start<one[0] and end>one[0]:     mark='yes';break
                                elif start<one[1] and end>one[1]:     mark='yes';break
                            if mark=='yes':   in_VSat1='1';#print(11111111)
                            else: in_VSat1='0'
                        ### Determine if near CoreVSat1
                        #print(dict_samplechr_bigblock_expand[samplechr])
                        if samplechr not in dict_samplechr_bigblock_expand:
                            in_Peri='.'
                        else:
                            core_regions= dict_samplechr_bigblock_expand[samplechr]
                            
                            if core_regions[0]<start and end<core_regions[1]:       in_Peri='1'
                            else:                                                   in_Peri='0'
                        ## Determine if in LIR
                        if samplechr not in dict_samplechr_interarrays or chromosome in ['Chr10','Chr11','Chr15','Chr16'] or (chromosome=="Chr18" and sample_type !="East_Asia"):
                            in_interarray='.'
                        else:
                            mark=''
                            for one in dict_samplechr_interarrays[samplechr]:
                                if one[0]<start and end<one[1]:     mark='yes';break
                            if mark=='yes':   in_interarray='1'
                            else: in_interarray='0'
                        f2.write(f"{eachline}\t{in_VSat1}\t{in_Peri}\t{in_interarray}\n")    
                        #sys.exit()
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.0" or argv1=="step5.0s2a": 
        ####
        print('Relying on EDTA\'s judgment to determine LTR sequences')
        # Generate reverse complement sequence
        def reverse_complement(sequence):
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','N':'N'}
            reverse_complement_seq = ''.join(complement[base] for base in reversed(sequence))
            return reverse_complement_seq  
        from Bio import pairwise2
        from Bio.Seq import Seq
        from Bio import Align


        print('Loading LTR sequences, etc.')
        for one_task in task_list:
            TE_Class,domain,domain_seq_name,tree_used_domain,  trimal_para=one_task
            with open(f'./chr2EDTA/5_tree/{TE_Class}/0_info_LTR','w') as f3:
                f3.write('sample\tsample_type\tchromosome\tstart\tend\tstrand\tTE_name\tClass\tClass1\tClass2\tClass3\tComplete\tDomains\tdom_names\tseqs\tseqs_used\tmafftinput\tseq_index\tLTR1_start\tLTR1_end\tLTR1_len\tseq1\tLTR1_TSD\tLTR2_start\tLTR2_end\tLTR2_len\tseq2\tLTR2_TSD\tsimilarity\ttime_years_Ma\n')
                
            result_list=[]
            with open(f'./chr2EDTA/5_tree/{TE_Class}/0_info','r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=line.split('\t')
                    #sample,sample_type,chromosome,start,end,strand,TE_name,Class,Class1,Class2,Class3,Complete,Domains,dom_names,seqs,seqs_used,mafftinput,seq_index=eachline_arr
                    result_list.append(eachline_arr)
            #result_list=result_list[0]
            def run_step(one_result):
                sample,sample_type,chromosome,start,end,strand,TE_name,Class,Class1,Class2,Class3,Complete,Domains,dom_names,seqs,seqs_used,mafftinput,seq_index=one_result
                #if [sample,sample_type,chromosome,start,end,strand,TE_name]!=['Baimunage_hap1', 'Eurasian', 'Chr4', '21578485', '21583211', '-', 'TE_00001339_Chr4:21583211..21578485#LTR/Copia']:return False
                eachline='\t'.join(one_result).strip()
                with open(f'./chr2EDTA/5_tree/{TE_Class}/0_info_LTR','a') as f3:
                    TE_name_simple=TE_name.split("_Chr")[0]
                    good_info_list=[]
                    record_mark=''
                    left_LTR=''
                    right_LTR=''
                    with open(f"./chr2EDTA/1_EDTA/{sample}/{sample}.fasta.mod.EDTA.intact.gff3",'r') as f2:
                        for line2 in f2:
                            eachline_arr2=line2.strip().split('\t')
                            if len(eachline_arr2)!=9:continue
                            chromosome2=eachline_arr2[0]
                            if chromosome2 != chromosome:continue
                            type2=eachline_arr2[2]
                            
                            if type2=='long_terminal_repeat' or type2=='repeat_region':
                                #
                                infos=eachline_arr2[8].split(';')
                                if  type2=='repeat_region':
                                    TE_name2=infos[1][5:]   
                                    if  TE_name2==TE_name_simple:
                                        start2=int(eachline_arr2[3])
                                        end2=int(eachline_arr2[4])
                                        if start2<=int(start) and end2>=int(end):
                                            record_mark='yes'
                                            continue
                                #print(eachline_arr2)
                                if record_mark=='yes':
                                    #print(eachline_arr2)
                                    #print(111111111)
                                    #print(eachline_arr2)
                                    TE_name2=infos[2][5:]   
                                    if TE_name2==TE_name_simple and  type2=='long_terminal_repeat':
                                        start2=int(eachline_arr2[3])
                                        end2=int(eachline_arr2[4])
                                        length2=end2-start2+1
                                        TSD=infos[-1][4:]  #TSD
                                        if left_LTR=='' and start2==int(start) and end2<int(end):      ###### It seems to be exactly stuck at the boundary
                                            good_info_list.append([start2,end2,length2,TSD]);left_LTR='yes'
                                        elif left_LTR=='yes'  and end2==int(end):      ###### It seems to be exactly stuck at the boundary
                                            good_info_list.append([start2,end2,length2,TSD]);right_LTR='yes'                                           
                                        if len(good_info_list)==2:record_mark='';left_LTR='';right_LTR='';break
                    #print(eachline_arr)                
                    #print(good_info_list)
                    LTR1_start=         good_info_list[0][0]   
                    LTR1_end=           good_info_list[0][1]
                    LTR1_length=        good_info_list[0][2]
                    LTR1_TSD=           good_info_list[0][3]
                    LTR2_start=         good_info_list[1][0]   
                    LTR2_end=           good_info_list[1][1]
                    LTR2_length=        good_info_list[1][2]
                    LTR2_TSD=           good_info_list[1][3]   
                    pos1=f"{LTR1_start}-{LTR1_end}"
                    pos2=f"{LTR2_start}-{LTR2_end}"
                    if LTR1_start>LTR2_start:
                        print([sample,sample_type,chromosome,start,end,strand,TE_name]);
                        print(start,end)
                        print(f"{LTR1_start}\t{LTR2_start}");
                        print(good_info_list)
                    #sys.exit()
                    CMD = ["samtools","faidx",f"/home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample}/{sample}.fasta",f"{chromosome}:{pos1}"  ]
                    #print(CMD)
                    seq1=subprocess.run(CMD, capture_output=True, text=True).stdout.strip()
                    #print(seq1)
                    seq1=''.join(seq1.split('\n')[1:]).strip().upper()
                    CMD = ["samtools","faidx",f"/home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample}/{sample}.fasta",f"{chromosome}:{pos2}"  ]
                    seq2=subprocess.run(CMD, capture_output=True, text=True).stdout.strip()
                    seq2=''.join(seq2.split('\n')[1:])
                    #seq2=reverse_complement(seq2).strip().upper()
                  
                    ##

                    # Global alignment (gap_open=-2, gap_extend=-1)
                    aligner = Align.PairwiseAligner()
                    aligner.mode = "global"
                    aligner.match_score = 1
                    aligner.mismatch_score = -0.5
                    aligner.open_gap_score = -5
                    aligner.extend_gap_score = -2
                    
                    # Perform alignment and get the best alignment result
                    alignments = aligner.align(seq1.upper(), seq2.upper())
                    best_alignment = alignments[0]  # Take the first (best) alignment result
                    
                    # Extract aligned sequences (including gaps)
                    aligned_seq1, aligned_seq2 = str(best_alignment.sequences[0]), str(best_alignment.sequences[1])
                    '''
                    # Calculate number of matches (only count positions without gaps)
                    matches = 0
                    total_positions = 0  # Total length of non-gap alignment
                    
                    for a, b in zip(aligned_seq1, aligned_seq2):
                        if a != '-' and b != '-':  # Only count non-gap positions
                            total_positions += 1
                            if a == b:
                                matches += 1
                    # Avoid division by zero
                    similarity = (matches / total_positions) * 100 if total_positions > 0 else 0.0                                
                    '''
                    #### Only compare positions without indels
                    def calculate_snp_similarity(aligned_seq1, aligned_seq2):
                        matches = 0
                        total_non_gap_positions = 0
                        
                        for a, b in zip(aligned_seq1, aligned_seq2):
                            if a != '-' and b != '-':  # Only count columns where both have no gaps
                                total_non_gap_positions += 1
                                if a == b:
                                    matches += 1
                        similarity = (matches / total_non_gap_positions) * 100 if total_non_gap_positions > 0 else 0.0
                        return similarity            
                    similarity = calculate_snp_similarity(aligned_seq1, aligned_seq2)
                    #similarity=1
                    #LTR1 and LTR2 sequence divergence (number of substitutions per base, i.e., Kimura 2-parameter distance)
                    def calculate_insertion_time(similarity_percent, ltr_length, mutation_rate):
                        p = 1 - similarity_percent / 100.0
                        if p >= 1:
                            return float('inf')  # Completely different, cannot calculate
                        K = -0.5 * math.log(1 - (2 * p) / (1 + p))
                        time_years = K / (2 * mutation_rate)
                        return time_years / 1e6  # Convert to millions of years
                    time_years=  calculate_insertion_time(similarity,(LTR1_length+LTR2_length)/2,2.5e-9)  
                    #print(seq1)
                    newline=f"{eachline}"
                    newline+=f"\t{LTR1_start}\t{LTR1_end}\t{LTR1_length}\t{seq1}\t{LTR1_TSD}"
                    newline+=f"\t{LTR2_start}\t{LTR2_end}\t{LTR2_length}\t{seq2}\t{LTR2_TSD}\t{similarity:.2f}\t{time_years}\n"
                    f3.write(newline)
     
            with Pool(processes=thread) as pool:
                # Use imap to get results one by one
                for i, result in enumerate(pool.imap(run_step, result_list), start=1):
                    # Results can be processed here, e.g., stored or printed
                    progress = (i / len(result_list)) * 100
                    sys.stdout.write(f"\rProgress: {progress:.2f}%")
                    sys.stdout.flush()           
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.0" or argv1=="step5.0s2b":
        TE_Class='Tekay'
        if  os.path.exists(f'./chr2EDTA/5_tree/{TE_Class}/LTRtime_plot')!=True:
            subprocess.run([f"mkdir ./chr2EDTA/5_tree/{TE_Class}/LTRtime_plot"], shell=True)     
        
        with open(f"./chr2EDTA/5_tree/{TE_Class}/LTRtime_plot/0_info_LTR",'w') as f2:
            f2.write(f"sample\tsample_type\tchromosome\tsimilarity\ttime_years_Ma\n")
            with open(f"./chr2EDTA/5_tree/{TE_Class}/0_info_LTR",'r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    sample,sample_type,chromosome=eachline_arr[:3]
                    similarity,time_years_Ma=eachline_arr[-2:]
                    f2.write(f"{sample}\t{sample_type}\t{chromosome}\t{similarity}\t{time_years_Ma}\n")
                
        #subprocess.run([f"cp ./chr2EDTA/5_tree/{TE_Class}/0_info_LTR ./chr2EDTA/5_tree/{TE_Class}/LTRtime_plot/0_info_LTR "], shell=True)         
        R_txt=f"""
        library(ggplot2)
        library(dplyr)
          input_file=read.table('0_info_LTR', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
          
          list_chr <- c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19")
            input_file$chromosome <- factor(input_file$chromosome, levels = list_chr)
            ##time_years_Ma
        # Create plot object
        p <- ggplot() +
            geom_point(data=input_file,aes(x = chromosome, y = similarity	, color = sample_type)) 
         
          
          # Save as PDF
          pdf(file = paste0('plot', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
          print(p)
          dev.off()
         
        """
        with open(f'./chr2EDTA/5_tree/{TE_Class}/LTRtime_plot/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = f'./chr2EDTA/5_tree/{TE_Class}/LTRtime_plot'
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../../../')  
            
            
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.0" or argv1=="step5.0s3":
        # Generate reverse complement sequence
        def reverse_complement(sequence):
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','N':'N'}
            reverse_complement_seq = ''.join(complement[base] for base in reversed(sequence))
            return reverse_complement_seq  
        from Bio import pairwise2
        from Bio.Seq import Seq
        from Bio import Align
        print('JZX manual analysis of LTR regions')
        for one_task in task_list:
            TE_Class,domain,domain_seq_name,tree_used_domain,  trimal_para=one_task
            with open(f'./chr2EDTA/5_tree/{TE_Class}/0_info_LTR_JZX','w') as f3:
                f3.write('sample\tsample_type\tchromosome\tstart\tend\tstrand\tTE_name\tClass\tClass1\tClass2\tClass3\tComplete\tDomains\tdom_names\tseqs\tseqs_used\tmafftinput\tseq_index\tLTR1\tLTR2\tsimilarity\tmatches\ttotal_aligned\tdirection\n')
                
            result_list=[]
            with open(f'./chr2EDTA/5_tree/{TE_Class}/0_info','r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=line.split('\t')
                    #sample,sample_type,chromosome,start,end,strand,TE_name,Class,Class1,Class2,Class3,Complete,Domains,dom_names,seqs,seqs_used,mafftinput,seq_index=eachline_arr
                    result_list.append(eachline_arr)
            #result_list=result_list[0]
            def run_step(one_result):
                sample,sample_type,chromosome,start,end,strand,TE_name,Class,Class1,Class2,Class3,Complete,Domains,dom_names,seqs,seqs_used,mafftinput,seq_index=one_result
                #if [sample,sample_type,chromosome,start,end,strand,TE_name]!=['Baimunage_hap1', 'Eurasian', 'Chr4', '21578485', '21583211', '-', 'TE_00001339_Chr4:21583211..21578485#LTR/Copia']:return False
                eachline='\t'.join(one_result).strip()  
                start, end = int(start), int(end)

                # Define analysis windows at both ends of the LTR (e.g., TE boundary ±100bp, adjustable as needed)
                ltr_window = 5000  
                pos1 = f"{start - ltr_window}-{start}"  # Region near the left end of LTR
                pos2 = f"{end}-{end + ltr_window}"      # Region near the right end of LTR

                CMD = ["samtools","faidx",f"/home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample}/{sample}.fasta",f"{chromosome}:{pos1}"  ]
                seq1=subprocess.run(CMD, capture_output=True, text=True).stdout.strip()
                seq1=''.join(seq1.split('\n')[1:]).strip().upper()
                seq1_rc=reverse_complement(seq1)
                CMD = ["samtools","faidx",f"/home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample}/{sample}.fasta",f"{chromosome}:{pos2}"  ]
                seq2=subprocess.run(CMD, capture_output=True, text=True).stdout.strip()
                seq2=''.join(seq2.split('\n')[1:]).strip().upper()
                seq2_rc =reverse_complement(seq2)

                # Local alignment (gap_open=-2, gap_extend=-1)
                aligner = Align.PairwiseAligner()
                aligner.mode = "local"
                aligner.match_score = 1
                aligner.mismatch_score = -0.5
                aligner.open_gap_score = -5
                aligner.extend_gap_score = -2
                
                # Perform alignment
                alignments_forward = aligner.align(seq1, seq2)  # Forward strand alignment
                alignments_reverse = aligner.align(seq1, seq2_rc)  # Reverse strand alignment                
                
                # Determine the best alignment direction (choose the one with higher score)
                if alignments_reverse.score > alignments_forward.score:
                    best_alignment = alignments_reverse[0]  # Reverse alignment is better
                    alignment_type = 'reverse'
                else:
                    best_alignment = alignments_forward[0]  # Forward alignment is better (or equal)
                    alignment_type = 'forward'
                
                # Parse alignment results
                aligned_seq1 = str(best_alignment.sequences[0])  # First aligned sequence (seq1)
                aligned_seq2 = str(best_alignment.sequences[1])  # Second aligned sequence (seq2 or seq2_rc)
                
                # Get alignment coordinates (shape (2, n), n is the length after alignment)
                #coords = best_alignment.coordinates  # Correct: call on a single PairwiseAlignment object
                #print("Alignment Coordinates:\n", coords)
                
                # Optional: print aligned sequences
                #print("Aligned seq1:", aligned_seq1)
                #print("Aligned seq2:", aligned_seq2)
                
                # Calculate alignment length and percentage
                matches = 0
                total_aligned_positions = 0
                
                for a, b in zip(aligned_seq1, aligned_seq2):
                    if a != "-" and b != "-":  # Exclude gaps
                        total_aligned_positions += 1
                        if a == b:
                            matches += 1
                
                identity_percentage = (matches / total_aligned_positions) * 100 if total_aligned_positions > 0 else 0
                
                # Get alignment regions (length without gaps)
                aligned_segments = best_alignment.aligned
                seq1_start, seq1_end = aligned_segments[0][0]
                aligned_length = seq1_end - seq1_start  # Alignment length without gaps
                with open(f'./chr2EDTA/5_tree/{TE_Class}/0_info_LTR_JZX','a') as f3:
                    f3.write(f"{eachline}\t{alignment_type}\t{matches}\t{total_aligned_positions}\t{identity_percentage}\n")
                    
                    
                    
                    
                    
                    
                    
                    
            with Pool(processes=thread) as pool:
                # Use imap to get results one by one
                for i, result in enumerate(pool.imap(run_step, result_list), start=1):
                    # Results can be processed here, e.g., stored or printed
                    progress = (i / len(result_list)) * 100
                    sys.stdout.write(f"\rProgress: {progress:.2f}%")
                    sys.stdout.flush()           
        
        
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.1":    
        print('Estimated 10s')
        def run_step(one_task):
            TE_Class,domain,domain_seq_name,tree_used_domain,  trimal_para=one_task
            env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
            subprocess.run([f"{env_source_str} cd-hit -M 30000 -T 20 -d 100 -c {cdhid_c} -n 5 -i ./chr2EDTA/5_tree/{TE_Class}/0_index.fa -o ./chr2EDTA/5_tree/{TE_Class}/1_cdhit "], shell=True)
            
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, task_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(task_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()    
                
        
        ##############        
        
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.1" or argv1=="step5.1s": 
        def run_step(one_task):
            TE_Class,domain,domain_seq_name,tree_used_domain,  trimal_para=one_task
       
            def parse_cdhit_clstr(clstr_file):
                dict_core_others = {}
                current_cluster = None
                core_seq = None
                other_seqs = []
            
                with open(clstr_file, "r") as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith(">Cluster"):
                            # When encountering a new cluster, save the previous cluster's information (if any)
                            if current_cluster is not None:
                                dict_core_others[core_seq] = other_seqs.copy()
                                other_seqs.clear()
                            current_cluster = line
                        else:
                            # Parse sequence line
                            parts = line.split(">")
                            seq_name = parts[1].split("...")[0]  # Extract sequence name (remove trailing "...")
                            if line.endswith("*"):
                                core_seq = seq_name  # Representative sequence
                            else:
                                other_seqs.append(seq_name)  # Other sequences
                    
                    # Process the last cluster
                    if current_cluster is not None:
                        dict_core_others[core_seq] = other_seqs
            
                return dict_core_others
            
            # Example call
            clstr_file = f"./chr2EDTA/5_tree/{TE_Class}/1_cdhit.clstr"
            dict_core_others = parse_cdhit_clstr(clstr_file)
            ######
            ii=0
            dict_seq_serial={}
            with open(f"./chr2EDTA/5_tree/{TE_Class}/1_cdhit.index",'w') as f:
                f.write(f"group1\tgroup2\tcdhit_index\tseq_index\n")
                for core,others in dict_core_others.items():
                    ii+=1
                    f.write(f'{ii}\t0\t{core}\t{core}\n')
                    kk=0
                    for one in others:
                        kk+=1
                        f.write(f"{ii}\t{kk}\t{core}\t{one}\n")
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, task_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(task_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()                      
        
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.2":  
        print('There was an issue with this step previously. Should decompose the fasta by delimiter, perform mafft on each part separately, then merge the mafft results and run trimal')
                
        
        #print('Estimated 167s')
        def run_step(one_task):
            TE_Class,domain,domain_seq_name,tree_used_domain,  trimal_para=one_task
            tree_used_domain_num=len(tree_used_domain.split('|'))
            
            index_list = list(range(tree_used_domain_num))
            kk=0
            while kk<tree_used_domain_num:
                with open(f"./chr2EDTA/5_tree/{TE_Class}/1_cdhit.{kk}",'w') as f:
                    f.write('')
                kk+=1    
            
            seq_name_seq_list=[]
            with open(f"./chr2EDTA/5_tree/{TE_Class}/1_cdhit",'r') as f:
                for line in f:
                    eachline=line.strip()
                    if eachline[0]=='>': name=eachline[1:]
                    else: seq=eachline;     seq_name_seq_list.append([name,seq])
            name_arr=[]        
            for one in seq_name_seq_list:
                name,seq=one
                name_arr.append(name)
                seq_arr=seq.split(Interval)
                if len(seq_arr)!=tree_used_domain_num:print('error'); sys.exit()
                kk=0
                while kk<tree_used_domain_num:
                    with open(f"./chr2EDTA/5_tree/{TE_Class}/1_cdhit.{kk}",'a') as f:
                        f.write(f">{name}\n{seq_arr[kk]}\n")
                    kk+=1 
                    
            ###########
            
            ##########
            kk=0
            while kk<tree_used_domain_num:
                env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
                subprocess.run([f"{env_source_str} mafft --auto ./chr2EDTA/5_tree/{TE_Class}/1_cdhit.{kk} > ./chr2EDTA/5_tree/{TE_Class}/2_mafft.{kk}.fasta"], shell=True)
                kk+=1

            
            dict_kk_name_mafftseq={}
            dict_kk_len={}  
            kk=0
            while kk<tree_used_domain_num:
                with open(f"./chr2EDTA/5_tree/{TE_Class}/2_mafft.{kk}.fasta",'r') as f:
                    for line in f:
                        eachline=line.strip()
                        if eachline[0]=='>': 
                            name=eachline[1:]
                            if name not in dict_kk_name_mafftseq:
                                dict_kk_name_mafftseq[name]={}
                            dict_kk_name_mafftseq[name][kk]=''
                            dict_kk_len[kk]=0
                        else: 
                            seq=eachline
                            dict_kk_name_mafftseq[name][kk]+=seq
                            dict_kk_len[kk]+=len(seq)
                kk+=1
            
            kk=0  
            with open(f"./chr2EDTA/5_tree/{TE_Class}/2_mafft.part_len",'w') as f:    
                f.write(f"part_index\tlength\n")
                while kk<tree_used_domain_num:
                    f.write(f"{kk}\t{dict_kk_len[kk]}\n")
                    kk+=1
            with open(f"./chr2EDTA/5_tree/{TE_Class}/2_mafft.fasta",'w') as f:
                for one_name in name_arr:
                    dict_kk_seq=dict_kk_name_mafftseq[one_name]
                    #print(dict_kk_seq)
                    kk=0
                    mafft_all=[]
                    while kk<tree_used_domain_num:
                        mafft_part=dict_kk_seq[kk]
                        mafft_all.append(mafft_part)
                        kk+=1
                    mafft_all_str=''.join(mafft_all)    
                    f.write(f">{one_name}\n{mafft_all_str}\n")
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, task_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(task_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush() 
                
                
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.3": 
        print('trimal')
        def run_step(one_task):
            TE_Class,domain,domain_seq_name,tree_used_domain,  trimal_para=one_task
            env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate orthofinder && "
            cmd=f"{env_source_str} trimal {trimal_para} -in  ./chr2EDTA/5_tree/{TE_Class}/2_mafft.fasta -out  ./chr2EDTA/5_tree/{TE_Class}/3_trimAl.fasta "
            print(cmd)
            subprocess.run([cmd], shell=True)
            print(f'trimal completed: {TE_Class}')
            
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, task_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(task_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()     
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.3s":
        print('Formatting')
        def run_step(one_task):
            TE_Class,domain,domain_seq_name,tree_used_domain,  trimal_para=one_task
            env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "

            subprocess.run([f'{env_source_str} fasta_formatter -i ./chr2EDTA/5_tree/{TE_Class}/3_trimAl.fasta -o ./chr2EDTA/5_tree/{TE_Class}/3_trimAl.fa'], shell=True)
            subprocess.run([f'mv ./chr2EDTA/5_tree/{TE_Class}/3_trimAl.fa ./chr2EDTA/5_tree/{TE_Class}/3_trimAl.fasta'], shell=True)
            print(f'trimal completed: {TE_Class}')
            
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, task_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(task_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()          
        
        
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.4":  
        def run_step(one_task):
            TE_Class,domain,domain_seq_name,tree_used_domain,  trimal_para=one_task
            env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
            subprocess.run([f"{env_source_str} iqtree2 -s ./chr2EDTA/5_tree/{TE_Class}/3_trimAl.fasta -m TEST -bb 1000 -nt AUTO -mem 0.8 -pre ./chr2EDTA/5_tree/{TE_Class}/4_iqtree.fasta "], shell=True)
            if  os.path.exists(f'./chr2EDTA/5_tree/{TE_Class}/5_iTOL.Newick')!=True:
                subprocess.run([f"touch ./chr2EDTA/5_tree/{TE_Class}/5_iTOL.Newick"], shell=True)
            ### -m TESTNEW tests more models, -m TEST, but I'll stick with the default for now
            ## -wsl is a parameter for outputting site log-likelihoods, which writes each site's contribution to the current tree to a file (.sitelh) in a specific format for deeper tree tests (e.g., topology comparison) by downstream tools like CONSEL. Let's remove this advanced parameter since it's not usable.
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, task_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(task_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush() 
        # -o specifies outgroup
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.4" or argv1=="step5.5":   
        for one_task in task_list:
            TE_Class,domain,domain_seq_name,tree_used_domain,  trimal_para=one_task
            if  os.path.exists(f'./chr2EDTA/5_tree/{TE_Class}/5_iTOL.Newick')!=True:
                subprocess.run([f"touch ./chr2EDTA/5_tree/{TE_Class}/5_iTOL.Newick"], shell=True)
            print('After visualizing the tree obtained from iqtree in iTOL, to get the sample order, click export to Newick format, copy into 5_iTOL.Newick')
            sys.exit()
        
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.6":          
        print('After visualizing the tree obtained from iqtree in iTOL, to get the sample order, click export to Newick format, copy into 5_iTOL.Newick')
        ##task_list=[['Athila']]
        #task_list=[['Tekay']]
        for one_task in task_list:
            TE_Class=one_task[0]
            import re
            with open(f"./chr2EDTA/5_tree/{TE_Class}/5_iTOL.Newick2",'w') as f2:
                with open(f"./chr2EDTA/5_tree/{TE_Class}/5_iTOL.Newick",'r') as f:
                    text=f.read().strip()
                cleaned_text = re.sub(r'\[.*?\]', '', text)
                f2.write(cleaned_text)
            from ete3 import Tree
            # Read Newick file
            tree = Tree(f"./chr2EDTA/5_tree/{TE_Class}/5_iTOL.Newick2")
            # Get sample order (tip labels)
            tip_order = [leaf.name for leaf in tree.get_leaves()]
            
            
            with open(f"./chr2EDTA/5_tree/{TE_Class}/6_serial",'w') as f:
                f.write('seq_index\tserial\n')
                ii=0
                for one in tip_order:
                    ii+=1
                    f.write(f"{one}\t{ii}\n")
            
            subprocess.run([f"rm ./chr2EDTA/5_tree/{TE_Class}/5_iTOL.Newick2"], shell=True)
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.7":                        
        for one_task in task_list:
            TE_Class=one_task[0]
            dict_seqindex_serial={} 
            with open(f"./chr2EDTA/5_tree/{TE_Class}/6_serial",'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    seqindex,serial=eachline_arr
                    dict_seqindex_serial[seqindex]=serial
            
            dict_seqindex_cdhitindex={}
            with open(f"./chr2EDTA/5_tree/{TE_Class}/1_cdhit.index",'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    _,_,cdhit_index,seq_index=eachline_arr
                    dict_seqindex_cdhitindex[seq_index]=cdhit_index
            result_list=[]        
            with open(f"./chr2EDTA/5_tree/{TE_Class}/0_info_more",'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    seq_index=eachline_arr[17]
                    cdhit_index=dict_seqindex_cdhitindex[seq_index]
                    serial=dict_seqindex_serial[cdhit_index]
                    result_list.append([eachline,seq_index,cdhit_index,serial])
                    
                    
            sorted_result=sorted(result_list, key=lambda x:  int(x[3]))        
            with open(f"./chr2EDTA/5_tree/{TE_Class}/7_info_onlycdhitresult",'w') as f2:
                f2.write(f"sample\tsample_type\tchromosome\tstart\tend\tstrand\tTE_name\tClass\tClass1\tClass2\tClass3\tComplete\tDomains\tdom_names\tseqs\tseqs_used\tmafftinput\tseq_index\tin_VSat1\tin_Peri\tin_interarray\tcdhit_index\tserial\n")
                with open(f"./chr2EDTA/5_tree/{TE_Class}/7_info",'w') as f:
                    f.write(f"sample\tsample_type\tchromosome\tstart\tend\tstrand\tTE_name\tClass\tClass1\tClass2\tClass3\tComplete\tDomains\tdom_names\tseqs\tseqs_used\tmafftinput\tseq_index\tin_VSat1\tin_Peri\tin_interarray\tcdhit_index\tserial\n")
                    for one in sorted_result:
                        eachline,seq_index,cdhit_index,serial=one
                        f.write(f"{eachline}\t{cdhit_index}\t{serial}\n")
                        
                        if seq_index==cdhit_index:
                            f2.write(f"{eachline}\t{cdhit_index}\t{serial}\n")
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.7_addition":
        print('Go to the EDTA structure to find the LTR similarity column, append it')
        sample_set=set()
        for one_task in task_list:
            TE_Class=one_task[0]        
            ###
            with open(f"./chr2EDTA/5_tree/{TE_Class}/7_info",'r') as f:
                i=0
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    i+=1
                    if i==1:head=line ;continue
                    sample=eachline_arr[0]
                    chromosome=eachline_arr[2]
                    TE_name=eachline_arr[6]
                    sample_set.add(sample)
            sample_list=list(sample_set)    
        sample_list_len=len(sample_list)    
        print('Loading EDTA results. Store as dictionary')    
        dict_sample_TE_LTRidentity={}
        kk=0
        for one_sample in sample_list:
            #if one_sample!='V005.hap1':continue
            kk+=1
            print(f'Progress: {kk}/{sample_list_len}',end='\r')
            with open(f"./chr2EDTA/1_EDTA/{one_sample}/{one_sample}.fasta.mod.EDTA.intact.gff3",'r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t') 
                    if len(eachline_arr)!=9:continue
                    posmin= int(eachline_arr[3])
                    posmax= int(eachline_arr[4])
                    chromosome=eachline_arr[0]
                    strand=eachline_arr[6]
                    if strand=='+':pos1=posmin;pos2=posmax
                    else:pos1=posmax;pos2=posmin
                    
                    one_info= eachline_arr[-1]
                    #if "TE_0000065" in line:print(eachline_arr)        
                    #if "ltr_identity" not in one_info:continue
                    one_info_arr=one_info.split(';')
                    TE_name='';ltr_identity=''
                        
                    for one in one_info_arr:
                        if one[:5]=='Name=':TE_name=f"{one[5:]}_{chromosome}:{pos1}..{pos2}"
                        if one[:13]=='ltr_identity=':ltr_identity=one[13:];break
                    #if  TE_name=='' or ltr_identity==''  :print('error') 
                    if one_sample not in dict_sample_TE_LTRidentity:
                        dict_sample_TE_LTRidentity[one_sample]={}
                    dict_sample_TE_LTRidentity[one_sample][TE_name]    =ltr_identity
                    #if "TE_00000658" in TE_name:print(TE_name,one_sample,ltr_identity)            
        print('\n\nOutput results')                
        for one_task in task_list:
            TE_Class=one_task[0]        
            ###
            with open(f"./chr2EDTA/5_tree/{TE_Class}/7_info_EDTA_LTRidentity",'w') as f2:
                with open(f"./chr2EDTA/5_tree/{TE_Class}/7_info",'r') as f:
                    i=0
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        i+=1
                        if i==1:
                            f2.write(f"{eachline}\tEDTA_ltr_identity\n");continue
                        one_sample=eachline_arr[0]
                        #if one_sample!='V005.hap1':continue
                        #chromosome=eachline_arr[2]
                        TE_name=eachline_arr[6].split('#')[0]
                        #if TE_name[:3]=='Chr':      TE_name_revise=TE_name.split('_')[0]
                        #else:                       TE_name_revise='_'.join(TE_name.split('_')[:2]) 
                        if TE_name not in dict_sample_TE_LTRidentity[one_sample]:
                            ltr_identity='.';
                            #if "TE_00000658" in TE_name:print('2:',one_sample,TE_name)
                        else:ltr_identity=dict_sample_TE_LTRidentity[one_sample][TE_name] 
                        f2.write(f"{eachline}\t{ltr_identity}\n")
                     

    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.8":        
        for one_task in task_list:
            TE_Class=one_task[0]
            
            dict_cdhitindex_info={}                   
            with open(f"./chr2EDTA/5_tree/{TE_Class}/7_info",'r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    in_VSat1,in_Peri,in_interarray,cdhit_index,serial=eachline_arr[-5:]
                    if cdhit_index not in dict_cdhitindex_info: 
                        dict_cdhitindex_info[cdhit_index]={}
                        dict_cdhitindex_info[cdhit_index]['in_VSat1']={}        ;dict_cdhitindex_info[cdhit_index]['in_VSat1']['yes']=0;dict_cdhitindex_info[cdhit_index]['in_VSat1']['no']=0
                        dict_cdhitindex_info[cdhit_index]['in_Peri']={}         ;dict_cdhitindex_info[cdhit_index]['in_Peri']['yes']=0;dict_cdhitindex_info[cdhit_index]['in_Peri']['no']=0
                        dict_cdhitindex_info[cdhit_index]['in_interarray']={}   ;dict_cdhitindex_info[cdhit_index]['in_interarray']['yes']=0;dict_cdhitindex_info[cdhit_index]['in_interarray']['no']=0
                    if in_VSat1=='1':   dict_cdhitindex_info[cdhit_index]['in_VSat1']['yes']+=1
                    elif in_VSat1=='0':   dict_cdhitindex_info[cdhit_index]['in_VSat1']['no']+=1
                    ##
                    if in_Peri=='1':   dict_cdhitindex_info[cdhit_index]['in_Peri']['yes']+=1
                    elif in_Peri=='0':   dict_cdhitindex_info[cdhit_index]['in_Peri']['no']+=1
                    ##
                    if in_interarray=='1':   dict_cdhitindex_info[cdhit_index]['in_interarray']['yes']+=1
                    elif in_interarray=='0':   dict_cdhitindex_info[cdhit_index]['in_interarray']['no']+=1     
                    
            cdhit_index_num=len(list(dict_cdhitindex_info.keys()))     
            valid_size=90   ##350 degrees
            addition_blank=int(round(cdhit_index_num/valid_size *(360-valid_size),0))
            
            with open(f"./chr2EDTA/5_tree/{TE_Class}/8_info_bar_percent_for_iTOL",'w') as f2:
                f2.write("serial\tcdhit_index\tVSat1_yes\tVSat1_no\tVSat1_yes_percent\tVSat1_no_percent\tPeri_yes\tPeri_no\tPeri_yes_percent\tPeri_no_percent\tinterarray_yes\tinterarray_no\tinterarray_yes_percent\tinterarray_no_percent\n")
                with open(f"./chr2EDTA/5_tree/{TE_Class}/8_info_bar",'w') as f:
                    f.write('serial\tcdhit_index\tVSat1_type\tVSat1_num\tPeri_type\tPeri_num\tinterarray_type\tinterarray_num\n')
                    ii=0
                    for cdhit_index,dict_infos in dict_cdhitindex_info.items():
                        ii+=1
                        dict_VSat1=dict_infos['in_VSat1'] 
                        VSat1_yes=int(dict_VSat1['yes'])
                        VSat1_no=int(dict_VSat1['no'])
                        VSat1_yesno=VSat1_yes+VSat1_no
                        VSat1_yes_percent=round(VSat1_yes/VSat1_yesno,3)if VSat1_yesno != 0 else 0
                        VSat1_no_percent=1-VSat1_yes_percent
                        #
                        dict_Peri=dict_infos['in_Peri'] 
                        Peri_yes=int(dict_Peri['yes'])
                        Peri_no=int(dict_Peri['no'])
                        Peri_yesno=Peri_yes+Peri_no
                        Peri_yes_percent=round(Peri_yes/Peri_yesno,3)if Peri_yesno != 0 else 0
                        Peri_no_percent=1-Peri_yes_percent                        
                        ##
                        dict_interarray=dict_infos['in_interarray'] 
                        interarray_yes=int(dict_interarray['yes'])
                        interarray_no=int(dict_interarray['no'])
                        interarray_yesno=interarray_yes+interarray_no
                        interarray_yes_percent=round(interarray_yes/interarray_yesno,3)if interarray_yesno != 0 else 0
                        interarray_no_percent=1-interarray_yes_percent                           
                        ##
                        f.write(f"{ii}\t{cdhit_index}\tyes\t{VSat1_yes}\tyes\t{Peri_yes}\tyes\t{interarray_yes}\n")
                        f.write(f"{ii}\t{cdhit_index}\tno\t{VSat1_no}\tno\t{Peri_no}\tno\t{interarray_no}\n")
            ################
                        txt=  f"{VSat1_yes}\t{VSat1_no}\t{VSat1_yes_percent}\t{VSat1_no_percent}\t"
                        txt+= f"{Peri_yes}\t{Peri_no}\t{Peri_yes_percent}\t{Peri_no_percent}\t"
                        txt+= f"{interarray_yes}\t{interarray_no}\t{interarray_yes_percent}\t{interarray_no_percent}"
                        f2.write(f"{ii}\t{cdhit_index}\t{txt}\n")
                    
                    kk=0
                    while kk<addition_blank-1:
                        kk+=1
                        ii+=1 
                        f.write(f"{ii}\t.\tyes\t100\tyes\t100\tyes\t100\n")
                        f.write(f"{ii}\t.\tno\t100\tno\t100\tno\t100\n")
            
    if argv1=="stepall"  or argv1=="step5"or argv1=="step5.8" or argv1=="step5.8p":        
        for one_task in task_list:
            TE_Class=one_task[0]
            
            R_txt=f"""
            library(ggplot2)
            library(dplyr)
              input_file=read.table('8_info_bar', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
              # Convert cdhit to a factor and sort by integer from large to small
            input_file <- input_file %>%
              mutate(serial2 = factor(serial, levels = sort(unique(serial))))
            # Create plot object
            p <- ggplot() +
                geom_bar(data=input_file,aes(x = serial2, y = VSat1_num, fill = VSat1_type),stat = "identity", position = "stack") +
                theme_classic() +  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position = "none",axis.text.x = element_blank()) 
             
              
              # Save as PDF
              pdf(file = paste0('8_plot_VSat1', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
              print(p)
              dev.off()
            
            # Create plot object
            p <- ggplot() +
                geom_bar(data=input_file,aes(x = serial2, y = Peri_num, fill = Peri_type),stat = "identity", position = "stack") +
                theme_classic() +  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position = "none",axis.text.x = element_blank()) 
             
              
              # Save as PDF
              pdf(file = paste0('8_plot_Peri', ".pdf"), width = 200/ 2.54, height = 50 / 2.54)
              print(p)
              dev.off()

            
            # Create plot object
            p <- ggplot() +
                geom_bar(data=input_file,aes(x = serial2, y = Peri_num, fill = Peri_type),stat = "identity", position = "stack") +
                theme_classic() +  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position = "none",axis.text.x = element_blank()) 
            p <- p+ coord_polar()+ylim(-2000,NA)   ##+ scale_y_continuous(-500,100) 

              # Save as PDF
              pdf(file = paste0('8_plot_Peri_circle', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
              print(p)
              dev.off()                
            
            # Create plot object
            p <- ggplot() +
                geom_bar(data=input_file,aes(x = serial2, y = interarray_num, fill = interarray_type),stat = "identity", position = "stack") +
                theme_classic() + # theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position = "none",axis.text.x = element_blank()) 
              
              # Save as PDF
              pdf(file = paste0('8_plot_interarray', ".pdf"), width = 200/ 2.54, height = 10 / 2.54)
              print(p)
              dev.off()   
              
            
            # Create plot object
            p <- ggplot() +
                geom_bar(data=input_file,aes(x = serial2, y = interarray_num, fill = interarray_type),stat = "identity", position = "stack") +
                theme_classic() +  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position = "none",axis.text.x = element_blank()) 
            p <- p+ coord_polar()+ylim(-300,NA)   ##+ scale_y_continuous(-500,100) 

              # Save as PDF
              pdf(file = paste0('8_plot_interarray_circle', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
              print(p)
              dev.off()                
                """
            with open(f'./chr2EDTA/5_tree/{TE_Class}/8_plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)  
            new_directory = f'./chr2EDTA/5_tree/{TE_Class}'
            os.chdir(new_directory)
            subprocess.run(['Rscript 8_plot.R'], shell=True)    
            os.chdir('../../../')                             
           
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.10": 
        print('Select a portion to redo phylogenetic tree')
        index_list=[741,765,444,848,778,390,478,477,799,785,479,479,472,480,899,459,800,800,947,160,359,161,907,1379,520,520,858,939,545,3165,664,839,452,452,488,942,463,330,3134,3134,3134,329,527,609,270,949,700,746,13,940,156,918,918,661,920,918,920,920,918,919,919,919,919,918,918,918,1392,900,919,662,484,260,788,486,934,361,482,921,788,484,260,908,507,805,775,495,759,451,387,534,531,607,837,608,604,637,534,634,486,260,484,934,788,361,788,457,808,486,934,788,457,484,921,289,487,495,487,495,495,1393,170,385,829,516,613,3096,814,668,827,487,805,775,495,908,507,487,908,507,805,775,495,439,788,457,934,788,457,482,260,788,460,486,934,486,934,260,460,788,560,535,491,157,523,552,441,610,753,395,494,722,611,766,652,1409,690,740,708,706,709,579,365,579,928,553,591,677,726,738,713,691,579,758,731,639,665,718,714,933,728,696,623,714,844,3144,676,65,576,716,748,705,691,691,579,697,653,553,707,652,652,719,536,726,616,714,652,648,594,682,652,539,652,579,471,652,787,733,644,574,550,718,592,629,716,738,652,548,652,727,651,511,271,768,503,692,562,715,712,654,620,579,652,725,652,721,614,471,579,471,652,760,655,655,724,593,714,631,582,548,649,697,653,758,577,744,718,620,658,657,650,515,721,575,501,58,678,719,720,458,522,677,537,716,647,581,3066,581,652,652,652,745,565,769,768,502,652,652,620,656,652,766,575,501,58,561,631,582,738,713,691,579,578,698,577,719,726,566,652,620,543,830,827,602,935,827,567,169,590,944,446,542,542,761,564,772,938,567,3108,945,445,845,737,618,925,641,865,877,897,869,877,849,856,782,851,902,877,461,896,877,877,877,51,496,894,866,861,51,870,877,877,864,888,877,877,875,877,791,877,496,455,777,877,877,861,3189,877,790,877,499,870,877,455,877,860,888,877,877,892,3189,922,877,875,777,877,312,872,51,893,861,866,533,652,684,158,916,369,688,685,526,695,3096,1382,693,688,685,946,597,588,466,589,685,764,636,551,3070,666,685,683,687,685,485,540,601,532,3100,549,301,813,680,757,153,685,917,492,711,3071,679,685,583,563,598,525,570,917,736,1376,596,154,1378,752,456,872,890,930,861,872,877,781,929,932,932,872,877,893,930,877,877,500,877,877,881,877,887,871,487,876,805,775,889,877,934,843,805,775,891,877,930,877,870,877,862,877,877,312,877,888,877,877,875,877,854,455,777,877,864,892,3189,791,877,870,877,777,877,877,312,860,892,3189,877,922,875,877,499,455,888,877,877,854,499,931,870,877,862,868,877,877,312,860,864,892,3189,877,922,875,857,885,877,882,855,877,850,300,879,853,303,884,874,932,872,877,893,930,854,791,877,878,862,868,877,931,877,312,877,864,888,877,877,3189,875,877,870,877,455,777,867,877,877,783,932,859,932,930,872,51,877,893,861,51,877,893,861,872,880,852,863,877,877,886,877,883,877,3081,873,877,498,497,870,508,895,862,864,331,788,457,762,826,921,547,934,908,507,805,775,495,487,484,807,776,773,438,926,361,486,774,846,846,484,803,484,934,805,775,901,484,333,843,912,906,910,795,934,361,773,934,923,805,775,901,908,507,926,775,934,803,487,923,805,775,495,934,843,798,805,775,934,908,507,361,484,805,775,843,803,482,486,368,367,826,767,877,843,506,3068,816,482,805,775,934,460,366,568,930,675,569,826,877,843,660,518,646,474,779,571,448,909,786,521,620,509,3067,554,530,620,640,606,3072,467,750,674,802,671,673,827,334,642,528,635,672,827,827,543,827,817,702,166,618,3168,835,796,832,735,556,486,484,921,788,457,289,486,788,361,486,788,457,332,332,809,362,486,934,260,460,482,788,439,934,788,457,788,457,294,734,612,295,453,780,3098,723,903,847,789,3099,669,632,810,1410,1394,630,659,579,905,505,840,3069,464,730,580,294,734,295,453,612,780,754,755,595,685,686,685,489,517,685,685,364,470,3077,751,749,584,490,469,363,363,792,941,8,743,620,624,927,386,948,667,620,710,546,572,620,572,159,620,450,605,605,620,620,948,619,1422,573,572,538,820,937,913,389,842,827,833,914,541,558,904,838,620,811,827,559,842,770,822,821,633,702,475,599,819,600,468,771,834,793,642,542,704,815,671,823,898,645,476,557,465,642,299,544,812,449,473,824,827,529,818,3097,739,821,828,831,670,557,168,671,671,825,911,483,943,836,841,454,703,756,504,756,524,586,11,620,587,717,462,442,481,514,763,934,628,806,804,519,484,512,699,493,732,613,685,302,625,681,694,167,827,1411,805,775,908,507,495,742,701,620,620,617,729,621,620,797,620,603,622,620,626,620,924,620,627,615,620,620,617,620,915,3175,620,797,620,638,620,924,620,555,620,585,1413,1413,1413,663,801,1413,784,447,784,784,447,784,784,447,784,784,447,784,784,447,643,513,689,3083,3082,747,510,936,794,155,1096]
        print('Outgroup index is 1096')
        outgroup_list=[1096]


        TE_Class='Tekay'
        tree_used_domain_num=3
        if  os.path.exists(f'./chr2EDTA/5_tree/{TE_Class}/10_ClassII')!=True:
            subprocess.run([f"mkdir ./chr2EDTA/5_tree/{TE_Class}/10_ClassII"], shell=True)  
        
        ii=0
        dict_newindex_oldindex={}
        new_index_list=[]
        outgroup_newindex_list=[]
        with open(f"./chr2EDTA/5_tree/{TE_Class}/10_ClassII/0_info",'w') as f2:
            f2.write("sample\tsample_type\tchromosome\tstart\tend\tstrand\tTE_name\tClass\tClass1\tClass2\tClass3\tComplete\tDomains\tdom_names\tseqs\tseqs_used\tmafftinput\tseq_index\tin_VSat1\tin_Peri\tin_interarray\tcdhit_index\tserial\tnew_index\n")
            with open(f"./chr2EDTA/5_tree/{TE_Class}/7_info",'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    seq_index=int(eachline_arr[17])
                    if seq_index in index_list:
                        ii+=1
                        if seq_index in outgroup_list:
                            newindex=f"outgroup_{ii}"
                            outgroup_newindex_list.append(newindex)
                        else:
                            newindex=f"newindex_{ii}"
                        dict_newindex_oldindex[newindex]=seq_index
                        f2.write(f"{eachline}\t{newindex}\n")
                        new_index_list.append(f"{newindex}")
        outgroup_newindex_list_str=','.join(outgroup_newindex_list)                
            
        dict_oldindex_seq={}
        with open(f"./chr2EDTA/5_tree/{TE_Class}/0_index.fa",'r') as f:
            for line in f:
                eachline=line.strip()
                if eachline[0]=='>':name=int(eachline[1:])
                else:seq=eachline;dict_oldindex_seq[name]=seq
                
        #print(dict_oldindex_seq)        
        ############
        kk=0
        while kk<tree_used_domain_num:
            with open(f"./chr2EDTA/5_tree/{TE_Class}/10_ClassII/0_index.{kk}",'w') as f:
                f.write('')
            kk+=1   
        ##############
        with open(f"./chr2EDTA/5_tree/{TE_Class}/10_ClassII/0_index.fa",'w') as f:
            for one_index in new_index_list:
                old_index=dict_newindex_oldindex[one_index]
                seq=dict_oldindex_seq[old_index]
                seq_arr=seq.split(Interval)
                if len(seq_arr)!=tree_used_domain_num: print('error');sys.exit()
                f.write(f">{one_index}\n{seq}\n")
                ##########
                kk=0
                while kk<tree_used_domain_num:
                    with open(f"./chr2EDTA/5_tree/{TE_Class}/10_ClassII/0_index.{kk}",'a') as f2:
                        f2.write(f">{one_index}\n{seq_arr[kk]}\n")
                    kk+=1     
                ###########

        kk=0
        while kk<tree_used_domain_num:
            env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
            subprocess.run([f"{env_source_str} mafft --auto ./chr2EDTA/5_tree/{TE_Class}/10_ClassII/0_index.{kk} > ./chr2EDTA/5_tree/{TE_Class}/10_ClassII/1_mafft.{kk}.fasta"], shell=True)
            kk+=1

    
        dict_kk_name_mafftseq={}
        dict_kk_len={}  
        kk=0
        while kk<tree_used_domain_num:
            with open(f"./chr2EDTA/5_tree/{TE_Class}/10_ClassII/1_mafft.{kk}.fasta",'r') as f:
                for line in f:
                    eachline=line.strip()
                    if eachline[0]=='>': 
                        name=eachline[1:]
                        if name not in dict_kk_name_mafftseq:
                            dict_kk_name_mafftseq[name]={}
                        dict_kk_name_mafftseq[name][kk]=''
                        dict_kk_len[kk]=0
                    else: 
                        seq=eachline
                        dict_kk_name_mafftseq[name][kk]+=seq
                        dict_kk_len[kk]+=len(seq)
            kk+=1
        
        kk=0  
        with open(f"./chr2EDTA/5_tree/{TE_Class}/10_ClassII/1_mafft.part_len",'w') as f:    
            f.write(f"part_index\tlength\n")
            while kk<tree_used_domain_num:
                f.write(f"{kk}\t{dict_kk_len[kk]}\n")
                kk+=1
        with open(f"./chr2EDTA/5_tree/{TE_Class}/10_ClassII/1_mafft.fasta",'w') as f:
            for one_name in new_index_list:
                dict_kk_seq=dict_kk_name_mafftseq[one_name]
                #print(dict_kk_seq)
                kk=0
                mafft_all=[]
                while kk<tree_used_domain_num:
                    mafft_part=dict_kk_seq[kk]
                    mafft_all.append(mafft_part)
                    kk+=1
                mafft_all_str=''.join(mafft_all)    
                f.write(f">{one_name}\n{mafft_all_str}\n")            
        
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate orthofinder && "
        cmd=f"{env_source_str} trimal -automated1 -in  ./chr2EDTA/5_tree/{TE_Class}/10_ClassII/1_mafft.fasta -out  ./chr2EDTA/5_tree/{TE_Class}/10_ClassII/2_trimAl.fasta "
        subprocess.run([cmd], shell=True)
        
        #####
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
        subprocess.run([f'{env_source_str} fasta_formatter -i ./chr2EDTA/5_tree/{TE_Class}//10_ClassII/2_trimAl.fasta -o ./chr2EDTA/5_tree/{TE_Class}/10_ClassII/2_trimAl.fa'], shell=True)
        subprocess.run([f'mv ./chr2EDTA/5_tree/{TE_Class}/10_ClassII/2_trimAl.fa ./chr2EDTA/5_tree/{TE_Class}//10_ClassII/2_trimAl.fasta'], shell=True)
            
        
        
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
        subprocess.run([f"{env_source_str} iqtree2 -s ./chr2EDTA/5_tree/{TE_Class}/10_ClassII/2_trimAl.fasta -m TEST -bb 1000 -nt AUTO -mem 0.8 -pre ./chr2EDTA/5_tree/{TE_Class}/10_ClassII/3_iqtree.fasta  -o {outgroup_newindex_list_str}"], shell=True)
        if  os.path.exists(f'./chr2EDTA/5_tree/{TE_Class}/10_ClassII/4_iTOL.Newick')!=True:
            subprocess.run([f"touch ./chr2EDTA/5_tree/{TE_Class}/10_ClassII/4_iTOL.Newick"], shell=True)        

        
    if argv1=="stepall"  or argv1=="step5" or argv1=="step5.11": 
        print('Select a portion to redo phylogenetic tree, but use an unrooted tree instead; the branch lengths of rooted trees are too extreme')
        index_list=[741,765,444,848,778,390,478,477,799,785,479,479,472,480,899,459,800,800,947,160,359,161,907,1379,520,520,858,939,545,3165,664,839,452,452,488,942,463,330,3134,3134,3134,329,527,609,270,949,700,746,13,940,156,918,918,661,920,918,920,920,918,919,919,919,919,918,918,918,1392,900,919,662,484,260,788,486,934,361,482,921,788,484,260,908,507,805,775,495,759,451,387,534,531,607,837,608,604,637,534,634,486,260,484,934,788,361,788,457,808,486,934,788,457,484,921,289,487,495,487,495,495,1393,170,385,829,516,613,3096,814,668,827,487,805,775,495,908,507,487,908,507,805,775,495,439,788,457,934,788,457,482,260,788,460,486,934,486,934,260,460,788,560,535,491,157,523,552,441,610,753,395,494,722,611,766,652,1409,690,740,708,706,709,579,365,579,928,553,591,677,726,738,713,691,579,758,731,639,665,718,714,933,728,696,623,714,844,3144,676,65,576,716,748,705,691,691,579,697,653,553,707,652,652,719,536,726,616,714,652,648,594,682,652,539,652,579,471,652,787,733,644,574,550,718,592,629,716,738,652,548,652,727,651,511,271,768,503,692,562,715,712,654,620,579,652,725,652,721,614,471,579,471,652,760,655,655,724,593,714,631,582,548,649,697,653,758,577,744,718,620,658,657,650,515,721,575,501,58,678,719,720,458,522,677,537,716,647,581,3066,581,652,652,652,745,565,769,768,502,652,652,620,656,652,766,575,501,58,561,631,582,738,713,691,579,578,698,577,719,726,566,652,620,543,830,827,602,935,827,567,169,590,944,446,542,542,761,564,772,938,567,3108,945,445,845,737,618,925,641,865,877,897,869,877,849,856,782,851,902,877,461,896,877,877,877,51,496,894,866,861,51,870,877,877,864,888,877,877,875,877,791,877,496,455,777,877,877,861,3189,877,790,877,499,870,877,455,877,860,888,877,877,892,3189,922,877,875,777,877,312,872,51,893,861,866,533,652,684,158,916,369,688,685,526,695,3096,1382,693,688,685,946,597,588,466,589,685,764,636,551,3070,666,685,683,687,685,485,540,601,532,3100,549,301,813,680,757,153,685,917,492,711,3071,679,685,583,563,598,525,570,917,736,1376,596,154,1378,752,456,872,890,930,861,872,877,781,929,932,932,872,877,893,930,877,877,500,877,877,881,877,887,871,487,876,805,775,889,877,934,843,805,775,891,877,930,877,870,877,862,877,877,312,877,888,877,877,875,877,854,455,777,877,864,892,3189,791,877,870,877,777,877,877,312,860,892,3189,877,922,875,877,499,455,888,877,877,854,499,931,870,877,862,868,877,877,312,860,864,892,3189,877,922,875,857,885,877,882,855,877,850,300,879,853,303,884,874,932,872,877,893,930,854,791,877,878,862,868,877,931,877,312,877,864,888,877,877,3189,875,877,870,877,455,777,867,877,877,783,932,859,932,930,872,51,877,893,861,51,877,893,861,872,880,852,863,877,877,886,877,883,877,3081,873,877,498,497,870,508,895,862,864,331,788,457,762,826,921,547,934,908,507,805,775,495,487,484,807,776,773,438,926,361,486,774,846,846,484,803,484,934,805,775,901,484,333,843,912,906,910,795,934,361,773,934,923,805,775,901,908,507,926,775,934,803,487,923,805,775,495,934,843,798,805,775,934,908,507,361,484,805,775,843,803,482,486,368,367,826,767,877,843,506,3068,816,482,805,775,934,460,366,568,930,675,569,826,877,843,660,518,646,474,779,571,448,909,786,521,620,509,3067,554,530,620,640,606,3072,467,750,674,802,671,673,827,334,642,528,635,672,827,827,543,827,817,702,166,618,3168,835,796,832,735,556,486,484,921,788,457,289,486,788,361,486,788,457,332,332,809,362,486,934,260,460,482,788,439,934,788,457,788,457,294,734,612,295,453,780,3098,723,903,847,789,3099,669,632,810,1410,1394,630,659,579,905,505,840,3069,464,730,580,294,734,295,453,612,780,754,755,595,685,686,685,489,517,685,685,364,470,3077,751,749,584,490,469,363,363,792,941,8,743,620,624,927,386,948,667,620,710,546,572,620,572,159,620,450,605,605,620,620,948,619,1422,573,572,538,820,937,913,389,842,827,833,914,541,558,904,838,620,811,827,559,842,770,822,821,633,702,475,599,819,600,468,771,834,793,642,542,704,815,671,823,898,645,476,557,465,642,299,544,812,449,473,824,827,529,818,3097,739,821,828,831,670,557,168,671,671,825,911,483,943,836,841,454,703,756,504,756,524,586,11,620,587,717,462,442,481,514,763,934,628,806,804,519,484,512,699,493,732,613,685,302,625,681,694,167,827,1411,805,775,908,507,495,742,701,620,620,617,729,621,620,797,620,603,622,620,626,620,924,620,627,615,620,620,617,620,915,3175,620,797,620,638,620,924,620,555,620,585,1413,1413,1413,663,801,1413,784,447,784,784,447,784,784,447,784,784,447,784,784,447,643,513,689,3083,3082,747,510,936,794,155]
        print('Outgroup index 1096 is disabled; no outgroup set')
        outgroup_list=[]


        TE_Class='Tekay'
        tree_used_domain_num=3
        if  os.path.exists(f'./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot')!=True:
            subprocess.run([f"mkdir ./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot"], shell=True)  
        
        ii=0
        dict_newindex_oldindex={}
        new_index_list=[]
        outgroup_newindex_list=[]
        with open(f"./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/0_info",'w') as f2:
            f2.write("sample\tsample_type\tchromosome\tstart\tend\tstrand\tTE_name\tClass\tClass1\tClass2\tClass3\tComplete\tDomains\tdom_names\tseqs\tseqs_used\tmafftinput\tseq_index\tin_VSat1\tin_Peri\tin_interarray\tcdhit_index\tserial\tnew_index\n")
            with open(f"./chr2EDTA/5_tree/{TE_Class}/7_info",'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    seq_index=int(eachline_arr[17])
                    if seq_index in index_list:
                        ii+=1
                        if seq_index in outgroup_list:
                            newindex=f"outgroup_{ii}"
                            outgroup_newindex_list.append(newindex)
                        else:
                            newindex=f"newindex_{ii}"
                        dict_newindex_oldindex[newindex]=seq_index
                        f2.write(f"{eachline}\t{newindex}\n")
                        new_index_list.append(f"{newindex}")
        outgroup_newindex_list_str=','.join(outgroup_newindex_list)                
            
        dict_oldindex_seq={}
        with open(f"./chr2EDTA/5_tree/{TE_Class}/0_index.fa",'r') as f:
            for line in f:
                eachline=line.strip()
                if eachline[0]=='>':name=int(eachline[1:])
                else:seq=eachline;dict_oldindex_seq[name]=seq
                
        #print(dict_oldindex_seq)        
        ############
        kk=0
        while kk<tree_used_domain_num:
            with open(f"./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/0_index.{kk}",'w') as f:
                f.write('')
            kk+=1   
        ##############
        with open(f"./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/0_index.fa",'w') as f:
            for one_index in new_index_list:
                old_index=dict_newindex_oldindex[one_index]
                seq=dict_oldindex_seq[old_index]
                seq_arr=seq.split(Interval)
                if len(seq_arr)!=tree_used_domain_num: print('error');sys.exit()
                f.write(f">{one_index}\n{seq}\n")
                ##########
                kk=0
                while kk<tree_used_domain_num:
                    with open(f"./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/0_index.{kk}",'a') as f2:
                        f2.write(f">{one_index}\n{seq_arr[kk]}\n")
                    kk+=1     
                ###########

        kk=0
        while kk<tree_used_domain_num:
            env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
            subprocess.run([f"{env_source_str} mafft --auto ./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/0_index.{kk} > ./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/1_mafft.{kk}.fasta"], shell=True)
            kk+=1

    
        dict_kk_name_mafftseq={}
        dict_kk_len={}  
        kk=0
        while kk<tree_used_domain_num:
            with open(f"./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/1_mafft.{kk}.fasta",'r') as f:
                for line in f:
                    eachline=line.strip()
                    if eachline[0]=='>': 
                        name=eachline[1:]
                        if name not in dict_kk_name_mafftseq:
                            dict_kk_name_mafftseq[name]={}
                        dict_kk_name_mafftseq[name][kk]=''
                        dict_kk_len[kk]=0
                    else: 
                        seq=eachline
                        dict_kk_name_mafftseq[name][kk]+=seq
                        dict_kk_len[kk]+=len(seq)
            kk+=1
        
        kk=0  
        with open(f"./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/1_mafft.part_len",'w') as f:    
            f.write(f"part_index\tlength\n")
            while kk<tree_used_domain_num:
                f.write(f"{kk}\t{dict_kk_len[kk]}\n")
                kk+=1
        with open(f"./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/1_mafft.fasta",'w') as f:
            for one_name in new_index_list:
                dict_kk_seq=dict_kk_name_mafftseq[one_name]
                #print(dict_kk_seq)
                kk=0
                mafft_all=[]
                while kk<tree_used_domain_num:
                    mafft_part=dict_kk_seq[kk]
                    mafft_all.append(mafft_part)
                    kk+=1
                mafft_all_str=''.join(mafft_all)    
                f.write(f">{one_name}\n{mafft_all_str}\n")            
        
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate orthofinder && "
        cmd=f"{env_source_str} trimal -automated1 -in  ./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/1_mafft.fasta -out  ./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/2_trimAl.fasta "
        subprocess.run([cmd], shell=True)
        
        #####
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
        subprocess.run([f'{env_source_str} fasta_formatter -i ./chr2EDTA/5_tree/{TE_Class}//11_ClassII_noroot/2_trimAl.fasta -o ./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/2_trimAl.fa'], shell=True)
        subprocess.run([f'mv ./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/2_trimAl.fa ./chr2EDTA/5_tree/{TE_Class}//11_ClassII_noroot/2_trimAl.fasta'], shell=True)
            
        
        
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
        subprocess.run([f"{env_source_str} iqtree2 -s ./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/2_trimAl.fasta -m MFP+MERGE  --bnni -bb 1000 -nt AUTO -mem 0.8 -pre ./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/3_iqtree.fasta  "], shell=True)
        if  os.path.exists(f'./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/4_iTOL.Newick')!=True:
            subprocess.run([f"touch ./chr2EDTA/5_tree/{TE_Class}/11_ClassII_noroot/4_iTOL.Newick"], shell=True)                    
            
                
if argv1=="stepall" or "step6" in argv1:  
    print('step6, Calculate correlation')
    if  os.path.exists('./chr2EDTA/6_POS2AA')!=True:
        subprocess.run(["mkdir ./chr2EDTA/6_POS2AA"], shell=True)    
    
    task_list=[]
    #task_list.append(["Ale",2])
    task_list.append(["Athila",2])  ### Name and number of sequence parts, mainly differences in RT, INT, and CHD, primarily whether CHD is present or not
    task_list.append(["CRM",3])
    task_list.append(["Galadriel",3])
    task_list.append(["Ogre",2])
    task_list.append(["Reina",3])
    task_list.append(["Retand",2])
    task_list.append(["Tekay",3])
    if argv1=="stepall" or argv1=="step6_readme":
        print('Printing instructions')
        with open('./chr2EDTA/6_POS2AA/readme','w') as f:
            txt=r'''
            step6.1 — Separate sequences, perform mafft on each, then combine into a large mafft alignment
            step6.2 — Merge info table
            step6.4 — Calculate correlation
            
            The approach in 6.4 is too rough because the entire Class may have huge differences
            
            6.9 — Only analyze Class2 of Tekay, with serial159.
            '''
            f.write(txt)     
    if argv1=="stepall"  or argv1=="step6" or argv1=="step6.1":  
        print('There was an issue with this step before. The fasta should be split by delimiter, each part aligned with mafft separately, then the mafft results merged and processed with trimal')
        #print('Estimated 167s')
        def run_step(one_task):
            TE_Class,tree_used_domain_num=one_task
            
            if  os.path.exists(f'./chr2EDTA/6_POS2AA/{TE_Class}')!=True:
                subprocess.run([f"mkdir ./chr2EDTA/6_POS2AA/{TE_Class}"], shell=True)    
            
            index_list = list(range(tree_used_domain_num))
            kk=0
            while kk<tree_used_domain_num:
                with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/0_index.{kk}",'w') as f:
                    f.write('')
                kk+=1    
            
            seq_name_seq_list=[]
            with open(f"./chr2EDTA/5_tree/{TE_Class}/0_index.fa",'r') as f:
                for line in f:
                    eachline=line.strip()
                    if eachline[0]=='>': name=eachline[1:]
                    else: seq=eachline;     seq_name_seq_list.append([name,seq])
            name_arr=[]        
            for one in seq_name_seq_list:
                name,seq=one
                name_arr.append(name)
                seq_arr=seq.split(Interval)
                if len(seq_arr)!=tree_used_domain_num:print('error'); print(one_task);sys.exit()
                kk=0
                while kk<tree_used_domain_num:
                    with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/0_index.{kk}",'a') as f:
                        f.write(f">{name}\n{seq_arr[kk]}\n")
                    kk+=1 
                    
            ###########
            
            ##########
            kk=0
            while kk<tree_used_domain_num:
                env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
                subprocess.run([f"{env_source_str} mafft --auto ./chr2EDTA/6_POS2AA/{TE_Class}/0_index.{kk} > ./chr2EDTA/6_POS2AA/{TE_Class}/1_mafft.{kk}.fasta"], shell=True)
                kk+=1

            
            dict_kk_name_mafftseq={}
            dict_kk_len={}  
            kk=0
            while kk<tree_used_domain_num:
                with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/1_mafft.{kk}.fasta",'r') as f:
                    for line in f:
                        eachline=line.strip()
                        if eachline[0]=='>': 
                            name=eachline[1:]
                            if name not in dict_kk_name_mafftseq:
                                dict_kk_name_mafftseq[name]={}
                            dict_kk_name_mafftseq[name][kk]=''
                            dict_kk_len[kk]=0
                        else: 
                            seq=eachline
                            dict_kk_name_mafftseq[name][kk]+=seq
                            dict_kk_len[kk]+=len(seq)
                kk+=1
            
            kk=0  
            with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/1_mafft.part_len",'w') as f:    
                f.write(f"part_index\tlength\n")
                while kk<tree_used_domain_num:
                    f.write(f"{kk}\t{dict_kk_len[kk]}\n")
                    kk+=1
            with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/1_mafft.fasta",'w') as f:
                for one_name in name_arr:
                    dict_kk_seq=dict_kk_name_mafftseq[one_name]
                    #print(dict_kk_seq)
                    kk=0
                    mafft_all=[]
                    while kk<tree_used_domain_num:
                        mafft_part=dict_kk_seq[kk]
                        mafft_all.append(mafft_part)
                        kk+=1
                    mafft_all_str=''.join(mafft_all)    
                    f.write(f">{one_name}\n{mafft_all_str}\n")        

        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, task_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(task_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()  
                
                
    if argv1=="stepall"  or argv1=="step6" or argv1=="step6.2":          
        for one_task in task_list:
            TE_Class=one_task[0]
            print(TE_Class)
            ##########
            dict_index_mafftseq={}
            with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/1_mafft.fasta",'r') as f:
                for line in f:
                    eachline=line.strip()
                    if eachline[0]==">":index=eachline[1:]; dict_index_mafftseq[index]=''
                    else: dict_index_mafftseq[index]+= eachline

            with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/2_sum",'w') as f2:
                f2.write("sample\tsample_type\tchromosome\tstart\tend\tstrand\tTE_name\tClass\tClass1\tClass2\tClass3\tComplete\tDomains\tdom_names\tseqs\tseqs_used\tmafftinput\tseq_index\tin_VSat1\tin_Peri\tin_interarray\tmafftresult\n")
                with open(f"./chr2EDTA/5_tree/{TE_Class}/0_info_more",'r') as f:
                    next(f)
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        seq_index=eachline_arr[17]
                        mafftresult=dict_index_mafftseq[seq_index]
                        f2.write(f"{eachline}\t{mafftresult}\n")

    if argv1=="stepall"  or argv1=="step6" or argv1=="step6.3" and 1==2:     
        import pandas as pd
        from scipy.stats import chi2_contingency
        for one_task in task_list:
            TE_Class=one_task[0]
            print(TE_Class)
            in_VSat1_yes_list=[];           in_VSat1_no_list=[]
            in_Peri_yes_list=[];            in_Peri_no_list=[]
            in_interarray_yes_list=[];      in_interarray_no_list=[]
            with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/2_sum",'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    in_VSat1=           eachline_arr[16]
                    in_Peri=            eachline_arr[17]
                    in_interarray=      eachline_arr[18]
                    mafftseq=eachline_arr[19]
                    #print(mafftseq);sys.exit()    
                    if in_VSat1!='1':       in_VSat1_yes_list.append(mafftseq)
                    else:                   in_VSat1_no_list.append(mafftseq)
                    if in_Peri!='1':        in_Peri_yes_list.append(mafftseq)
                    else:                   in_Peri_no_list.append(mafftseq)                        
                    if in_interarray!='1':  in_interarray_yes_list.append(mafftseq)
                    else:                   in_interarray_no_list.append(mafftseq)                    
            #print(in_VSat1_yes_list)    ;sys.exit()    
            mafftseq_len=len(in_VSat1_yes_list[0])        
            with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/3_AA2other",'w') as f2:
                f2.write(f"pos\tanimo\tall\tA_yes\tA_no\tB_yes\tB_no\tchi2\tPvalue\n")
                kk=0
                while kk<mafftseq_len:
                    yes_list=[]
                    for one_mafft in in_Peri_yes_list:
                        base=one_mafft[kk]
                        if base=="-":continue
                        yes_list.append(base)
                    no_list=[]
                    for one_mafft in in_Peri_no_list:
                        base=one_mafft[kk]
                        if base=="-":continue
                        no_list.append(base)                    
                    kk+=1
                    #print(len(yes_list),len(no_list))
                    # Count frequency of each amino acid in yes and no
                    # All possible amino acids (converted to list)
                    all_amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
                    # Initialize frequency count dictionary
                    yes_counts = {aa: 0 for aa in all_amino_acids}
                    no_counts = {aa: 0 for aa in all_amino_acids}
                    
                    # Actually count amino acids in yes_list and no_list
                    for base in yes_list:
                        if base in yes_counts:
                            yes_counts[base] += 1
                    
                    for base in no_list:
                        if base in no_counts:
                            no_counts[base] += 1
                    
                    # Fill in statistical results
                    for aa in all_amino_acids:
                        A_yes = yes_counts[aa]  # Count of current amino acid in yes list
                        A_no = no_counts[aa]   # Count of current amino acid in no list
                        B_yes = len(yes_list) - A_yes  # Total count of other amino acids in yes list
                        B_no = len(no_list) - A_no   # Total count of other amino acids in no list
                        total_observations = A_yes + A_no + B_yes + B_no

                        # Chi-square test (use underscores to ignore the two return values if not needed)
                        if total_observations<10:continue
                        if A_yes+A_no==0:continue
                        if A_yes+B_yes==0:continue
                        if B_yes+B_no==0:continue
                        if A_no+B_no==0:continue
                        chi2_statistic, p_value, dof, expected = chi2_contingency([[A_yes, A_no], [B_yes, B_no]])

         
                        f2.write(f"{kk}\t{aa}\t{total_observations}\t{A_yes}\t{A_no}\t{B_yes}\t{B_no}\t{round(chi2_statistic, 4)}\t{p_value}\n")
                    
    if argv1=="stepall"  or argv1=="step6" or argv1=="step6.4":     
        import pandas as pd
        from scipy.stats import chi2_contingency
        import itertools
        ######
        all_amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        amino_acid_pairs = list(itertools.combinations(all_amino_acids, 2))
        #print(amino_acid_pairs)
        #####
        for one_task in task_list:
            TE_Class=one_task[0]
            print()
            print()
            print(TE_Class+'#######################')
            in_VSat1_yes_list=[];           in_VSat1_no_list=[]
            in_Peri_yes_list=[];            in_Peri_no_list=[]
            in_interarray_yes_list=[];      in_interarray_no_list=[]
            with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/2_sum",'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    in_VSat1=           eachline_arr[18]
                    in_Peri=            eachline_arr[19]
                    in_interarray=      eachline_arr[20]
                    mafftresult=eachline_arr[21]
                    #print(mafftresult);sys.exit()    
                    if in_VSat1=='1':               in_VSat1_yes_list.append(mafftresult)
                    elif in_VSat1=='0':             in_VSat1_no_list.append(mafftresult)
                    if in_Peri=='1':                in_Peri_yes_list.append(mafftresult)
                    elif in_Peri=='0':              in_Peri_no_list.append(mafftresult)                        
                    if in_interarray=='1':          in_interarray_yes_list.append(mafftresult)
                    elif in_interarray=='0':        in_interarray_no_list.append(mafftresult)                    
            #print(in_VSat1_yes_list)    ;sys.exit()    
            type_list=[]
            type_list.append(['VSat1',in_VSat1_yes_list,in_VSat1_no_list])
            type_list.append(['Peri',in_Peri_yes_list,in_Peri_no_list])
            type_list.append(['Interarray',in_interarray_yes_list,in_interarray_no_list])
            
            
            for one_type in type_list:
                print(one_type[0])
                type_name,type_yes_list,type_no_list=one_type
                type_yes_list_len=len(type_yes_list)
                type_no_list_len=len(type_no_list)
                if len(type_yes_list)==0:print('No TEs located in the region');continue
                #print(f"{len(type_yes_list)}")
                mafftresult_len=len(type_no_list[0])   
                index_list = list(range(mafftresult_len))
                with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/4_AA2BB.{TE_Class}.{type_name}",'w') as f2:
                    f2.write(f"pos/{mafftresult_len}\tyes_num\tno_num\tyes_amino\tno_amino\tused\taminoIN/aminoOUT\tin_aa\tout_aa\tin_yes\tin_no\tout_yes\tout_no\tchi2\tPvalue\n")
                #kk=0
                #while kk<mafftresult_len:
                def run_step(kk):
                    yes_list=[]
                    for one_mafft in type_yes_list:
                        base=one_mafft[kk]
                        if base=="-":continue
                        yes_list.append(base)
                    no_list=[]
                    for one_mafft in type_no_list:
                        base=one_mafft[kk]
                        if base=="-":continue
                        no_list.append(base)                    
                    #kk+=1
                    yes_list_len=len(yes_list)
                    no_list_len=len(no_list)

                    #print(len(yes_list),len(no_list))
                    # Count frequency of each amino acid in yes and no
                    # All possible amino acids (converted to list)
                    
                    #def run_step(pair):
                    for pair in amino_acid_pairs:
                        A,B=pair
                        A_yes = yes_list.count(A)
                        A_no =  no_list.count(A)
                        B_yes = yes_list.count(B)
                        B_no =  no_list.count(B)          
                    
                        total_observations = A_yes + A_no + B_yes + B_no
                        
                        # Chi-square test (use underscores to ignore the two return values if not needed)

                        if not(total_observations<10 or A_yes+A_no==0 or A_yes+B_yes==0 or B_yes+B_no==0 or A_no+B_no==0):
                            chi2_statistic, p_value, dof, expected = chi2_contingency([[A_yes, A_no], [B_yes, B_no]])
                            if p_value<0.05:
                                ##
                                if      A_no==0:            in_aa=A;out_aa=B
                                elif    B_no==0:            in_aa=B;out_aa=A
                                elif A_yes/A_no>B_yes/B_no: in_aa=A;out_aa=B
                                else:                       in_aa=B;out_aa=A
                                if in_aa==A:        in_yes,in_no,out_yes,out_no=A_yes,A_no,B_yes,B_no
                                else:               in_yes,in_no,out_yes,out_no=B_yes,B_no,A_yes,A_no
                                ##
                                with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/4_AA2BB.{TE_Class}.{type_name}",'a') as f2:
                                    f2.write(f"{kk+1}\t{type_yes_list_len}\t{type_no_list_len}\t{yes_list_len}\t{no_list_len}\t{total_observations}\t{in_aa}/{out_aa}\t{in_aa}\t{out_aa}\t{in_yes}\t{in_no}\t{out_yes}\t{out_no}\t{round(chi2_statistic, 4)}\t{p_value}\n")
                                
            
                with Pool(processes=thread) as pool:
                    # Use imap to get results one by one
                    for i, result in enumerate(pool.imap(run_step, index_list), start=1):
                        # Results can be processed here, e.g., stored or printed
                        progress = (i / len(index_list)) * 100
                        sys.stdout.write(f"\rProgress: {progress:.2f}%")
                        sys.stdout.flush() 
                print()
                print()
        #if argv1=="stepall"  or argv1=="step6" or argv1=="step6.4s":
        for one_task in task_list:
            TE_Class=one_task[0]    
            type_list=['VSat1','Peri','Interarray']
            for one_type in type_list:
                print(one_type)
                if os.path.exists(f"./chr2EDTA/6_POS2AA/{TE_Class}/4_AA2BB.{TE_Class}.{one_type}")==False:continue
                result_list=[]
                i=0
                with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/4_AA2BB.{TE_Class}.{one_type}",'r') as f:
                    for line in f:
                        i+=1
                        eachline=line.strip()
                        if i==1:head=eachline;continue
                        eachline_arr=eachline.split('\t')
                        chi2=float(eachline_arr[-2])
                        result_list.append([eachline,chi2])
                sorted_list = sorted(result_list, key=lambda x: x[1],reverse=True)
                        
            
                with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/4_AA2BB.{TE_Class}.{one_type}.sort",'w') as f2:
                    f2.write(f"{head}\n")   
                    for one in sorted_list:
                        f2.write(one[0]+'\n')
            
                subprocess.run([f"mv ./chr2EDTA/6_POS2AA/{TE_Class}/4_AA2BB.{TE_Class}.{one_type}.sort ./chr2EDTA/6_POS2AA/{TE_Class}/4_AA2BB.{TE_Class}.{one_type}"], shell=True)
            
    if argv1=="stepall"  or argv1=="step6" or argv1=="step6.9":     
        current_TE='Tekay'   ## Only analyze this one
        current_seqindex=[741,765,444,848,778,390,478,477,799,785,479,479,472,480,899,459,800,800,947,160,359,161,907,1379,520,520,858,939,545,3165,664,839,452,452,488,942,463,330,3134,3134,3134,329,527,609,270,949,700,746,13,940,156,918,918,661,920,918,920,920,918,919,919,919,919,918,918,918,1392,900,919,662,484,260,788,486,934,361,482,921,788,484,260,908,507,805,775,495,759,451,387,534,531,607,837,608,604,637,534,634,486,260,484,934,788,361,788,457,808,486,934,788,457,484,921,289,487,495,487,495,495,1393,170,385,829,516,613,3096,814,668,827,487,805,775,495,908,507,487,908,507,805,775,495,439,788,457,934,788,457,482,260,788,460,486,934,486,934,260,460,788,560,535,491,157,523,552,441,610,753,395,494,722,611,766,652,1409,690,740,708,706,709,579,365,579,928,553,591,677,726,738,713,691,579,758,731,639,665,718,714,933,728,696,623,714,844,3144,676,65,576,716,748,705,691,691,579,697,653,553,707,652,652,719,536,726,616,714,652,648,594,682,652,539,652,579,471,652,787,733,644,574,550,718,592,629,716,738,652,548,652,727,651,511,271,768,503,692,562,715,712,654,620,579,652,725,652,721,614,471,579,471,652,760,655,655,724,593,714,631,582,548,649,697,653,758,577,744,718,620,658,657,650,515,721,575,501,58,678,719,720,458,522,677,537,716,647,581,3066,581,652,652,652,745,565,769,768,502,652,652,620,656,652,766,575,501,58,561,631,582,738,713,691,579,578,698,577,719,726,566,652,620,543,830,827,602,935,827,567,169,590,944,446,542,542,761,564,772,938,567,3108,945,445,845,737,618,925,641,865,877,897,869,877,849,856,782,851,902,877,461,896,877,877,877,51,496,894,866,861,51,870,877,877,864,888,877,877,875,877,791,877,496,455,777,877,877,861,3189,877,790,877,499,870,877,455,877,860,888,877,877,892,3189,922,877,875,777,877,312,872,51,893,861,866,533,652,684,158,916,369,688,685,526,695,3096,1382,693,688,685,946,597,588,466,589,685,764,636,551,3070,666,685,683,687,685,485,540,601,532,3100,549,301,813,680,757,153,685,917,492,711,3071,679,685,583,563,598,525,570,917,736,1376,596,154,1378,752,456,872,890,930,861,872,877,781,929,932,932,872,877,893,930,877,877,500,877,877,881,877,887,871,487,876,805,775,889,877,934,843,805,775,891,877,930,877,870,877,862,877,877,312,877,888,877,877,875,877,854,455,777,877,864,892,3189,791,877,870,877,777,877,877,312,860,892,3189,877,922,875,877,499,455,888,877,877,854,499,931,870,877,862,868,877,877,312,860,864,892,3189,877,922,875,857,885,877,882,855,877,850,300,879,853,303,884,874,932,872,877,893,930,854,791,877,878,862,868,877,931,877,312,877,864,888,877,877,3189,875,877,870,877,455,777,867,877,877,783,932,859,932,930,872,51,877,893,861,51,877,893,861,872,880,852,863,877,877,886,877,883,877,3081,873,877,498,497,870,508,895,862,864,331,788,457,762,826,921,547,934,908,507,805,775,495,487,484,807,776,773,438,926,361,486,774,846,846,484,803,484,934,805,775,901,484,333,843,912,906,910,795,934,361,773,934,923,805,775,901,908,507,926,775,934,803,487,923,805,775,495,934,843,798,805,775,934,908,507,361,484,805,775,843,803,482,486,368,367,826,767,877,843,506,3068,816,482,805,775,934,460,366,568,930,675,569,826,877,843,660,518,646,474,779,571,448,909,786,521,620,509,3067,554,530,620,640,606,3072,467,750,674,802,671,673,827,334,642,528,635,672,827,827,543,827,817,702,166,618,3168,835,796,832,735,556,486,484,921,788,457,289,486,788,361,486,788,457,332,332,809,362,486,934,260,460,482,788,439,934,788,457,788,457,294,734,612,295,453,780,3098,723,903,847,789,3099,669,632,810,1410,1394,630,659,579,905,505,840,3069,464,730,580,294,734,295,453,612,780,754,755,595,685,686,685,489,517,685,685,364,470,3077,751,749,584,490,469,363,363,792,941,8,743,620,624,927,386,948,667,620,710,546,572,620,572,159,620,450,605,605,620,620,948,619,1422,573,572,538,820,937,913,389,842,827,833,914,541,558,904,838,620,811,827,559,842,770,822,821,633,702,475,599,819,600,468,771,834,793,642,542,704,815,671,823,898,645,476,557,465,642,299,544,812,449,473,824,827,529,818,3097,739,821,828,831,670,557,168,671,671,825,911,483,943,836,841,454,703,756,504,756,524,586,11,620,587,717,462,442,481,514,763,934,628,806,804,519,484,512,699,493,732,613,685,302,625,681,694,167,827,1411,805,775,908,507,495,742,701,620,620,617,729,621,620,797,620,603,622,620,626,620,924,620,627,615,620,620,617,620,915,3175,620,797,620,638,620,924,620,555,620,585,1413,1413,1413,663,801,1413,784,447,784,784,447,784,784,447,784,784,447,784,784,447,643,513,689]
        import pandas as pd
        from scipy.stats import chi2_contingency
        import itertools
        ######
        all_amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V','-']
        amino_acid_pairs = list(itertools.combinations(all_amino_acids, 2))
        #print(amino_acid_pairs)
        #####
        for one_task in task_list:
            TE_Class=one_task[0]
            if TE_Class!=current_TE:continue
            subprocess.run([f"mkdir ./chr2EDTA/6_POS2AA/{TE_Class}/9_AA"], shell=True)            
        
            print()
            print()
            print(TE_Class+'#######################')
            in_VSat1_yes_list=[];           in_VSat1_no_list=[]
            in_Peri_yes_list=[];            in_Peri_no_list=[]
            in_interarray_yes_list=[];      in_interarray_no_list=[]
            with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/2_sum",'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    seq_index=  int(eachline_arr[17])
                    if seq_index not in current_seqindex:continue
                    in_VSat1=           eachline_arr[18]
                    in_Peri=            eachline_arr[19]
                    in_interarray=      eachline_arr[20]
                    mafftresult=    eachline_arr[21]
                    
                    #print(mafftresult);sys.exit()    
                    if in_VSat1=='1':               in_VSat1_yes_list.append(mafftresult)
                    elif in_VSat1=='0':             in_VSat1_no_list.append(mafftresult)
                    if in_Peri=='1':                in_Peri_yes_list.append(mafftresult)
                    elif in_Peri=='0':              in_Peri_no_list.append(mafftresult)                        
                    if in_interarray=='1':          in_interarray_yes_list.append(mafftresult)
                    elif in_interarray=='0':        in_interarray_no_list.append(mafftresult)                    
            #print(in_VSat1_yes_list)    ;sys.exit()    
            type_list=[]
            type_list.append(['VSat1',in_VSat1_yes_list,in_VSat1_no_list])
            type_list.append(['Peri',in_Peri_yes_list,in_Peri_no_list])
            type_list.append(['Interarray',in_interarray_yes_list,in_interarray_no_list])
            
            
            for one_type in type_list:
                print(one_type[0])
                type_name,type_yes_list,type_no_list=one_type
                type_yes_list_len=len(type_yes_list)
                type_no_list_len=len(type_no_list)
                if len(type_yes_list)==0:print('No TEs located in the region');continue
                #print(f"{len(type_yes_list)}")
                mafftresult_len=len(type_no_list[0])   
                index_list = list(range(mafftresult_len))
                with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/9_AA/0_AA2BB.{TE_Class}.{type_name}",'w') as f2:
                    f2.write(f"pos/{mafftresult_len}\tyes_num\tno_num\tyes_amino\tno_amino\tused\taminoIN/aminoOUT\tin_aa\tout_aa\tin_yes\tin_no\tout_yes\tout_no\tchi2\tPvalue\n")
                #kk=0
                #while kk<mafftresult_len:
                def run_step(kk):
                    yes_list=[]
                    for one_mafft in type_yes_list:
                        base=one_mafft[kk]
                        #if base=="-":continue
                        yes_list.append(base)
                    no_list=[]
                    for one_mafft in type_no_list:
                        base=one_mafft[kk]
                        #if base=="-":continue
                        no_list.append(base)                    
                    #kk+=1
                    yes_list_len=len(yes_list)
                    no_list_len=len(no_list)

                    #print(len(yes_list),len(no_list))
                    # Count frequency of each amino acid in yes and no
                    # All possible amino acids (converted to list)
                    
                    #def run_step(pair):
                    for pair in amino_acid_pairs:
                        A,B=pair
                        A_yes = yes_list.count(A)
                        A_no =  no_list.count(A)
                        B_yes = yes_list.count(B)
                        B_no =  no_list.count(B)          
                    
                        total_observations = A_yes + A_no + B_yes + B_no
                        
                        # Chi-square test (use underscores to ignore the two return values if not needed)

                        if not(total_observations<10 or A_yes+A_no==0 or A_yes+B_yes==0 or B_yes+B_no==0 or A_no+B_no==0):
                            chi2_statistic, p_value, dof, expected = chi2_contingency([[A_yes, A_no], [B_yes, B_no]])
                            if p_value<0.05:
                                ##
                                if      A_no==0:            in_aa=A;out_aa=B
                                elif    B_no==0:            in_aa=B;out_aa=A
                                elif A_yes/A_no>B_yes/B_no: in_aa=A;out_aa=B
                                else:                       in_aa=B;out_aa=A
                                if in_aa==A:        in_yes,in_no,out_yes,out_no=A_yes,A_no,B_yes,B_no
                                else:               in_yes,in_no,out_yes,out_no=B_yes,B_no,A_yes,A_no
                                ##
                                with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/9_AA/0_AA2BB.{TE_Class}.{type_name}",'a') as f2:
                                    f2.write(f"{kk+1}\t{type_yes_list_len}\t{type_no_list_len}\t{yes_list_len}\t{no_list_len}\t{total_observations}\t{in_aa}/{out_aa}\t{in_aa}\t{out_aa}\t{in_yes}\t{in_no}\t{out_yes}\t{out_no}\t{round(chi2_statistic, 4)}\t{p_value}\n")
                                
            
                with Pool(processes=thread) as pool:
                    # Use imap to get results one by one
                    for i, result in enumerate(pool.imap(run_step, index_list), start=1):
                        # Results can be processed here, e.g., stored or printed
                        progress = (i / len(index_list)) * 100
                        sys.stdout.write(f"\rProgress: {progress:.2f}%")
                        sys.stdout.flush() 
                print()
                print()
        #if argv1=="stepall"  or argv1=="step6" or argv1=="step6.4s":
        for one_task in task_list:
            TE_Class=one_task[0]    
            if TE_Class!=current_TE:continue
            
            type_list=['VSat1','Peri','Interarray']
            for one_type in type_list:
                
                print(one_type)
                if os.path.exists(f"./chr2EDTA/6_POS2AA/{TE_Class}/9_AA/0_AA2BB.{TE_Class}.{one_type}")==False:continue
                result_list=[]
                i=0
                with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/9_AA/0_AA2BB.{TE_Class}.{one_type}",'r') as f:
                    for line in f:
                        i+=1
                        eachline=line.strip()
                        if i==1:head=eachline;continue
                        eachline_arr=eachline.split('\t')
                        chi2=float(eachline_arr[-2])
                        result_list.append([eachline,chi2])
                sorted_list = sorted(result_list, key=lambda x: x[1],reverse=True)
     
            
                with open(f"./chr2EDTA/6_POS2AA/{TE_Class}/9_AA/0_AA2BB.{TE_Class}.{one_type}.sort",'w') as f2:
                    f2.write(f"{head}\n")   
                    for one in sorted_list:
                        f2.write(one[0]+'\n')
            
                subprocess.run([f"mv ./chr2EDTA/6_POS2AA/{TE_Class}/9_AA/0_AA2BB.{TE_Class}.{one_type}.sort ./chr2EDTA/6_POS2AA/{TE_Class}/9_AA/0_AA2BB.{TE_Class}.{one_type}"], shell=True)
    
        
if argv1=="stepall" or "step7" in argv1:  
    print('step6, Calculate correlation')
    if  os.path.exists('./chr2EDTA/7_stat')!=True:
        subprocess.run(["mkdir ./chr2EDTA/7_stat"], shell=True)    
    
    task_list=[]
    #task_list.append(["Ale",2])
    task_list.append(["Athila"])  ### Names and the number of sequence parts, mainly differences in RT, INT, and CHD, primarily whether CHD is present or not
    task_list.append(["CRM"])
    task_list.append(["Galadriel"])
    task_list.append(["Ogre"])
    task_list.append(["Reina"])
    task_list.append(["Retand"])
    task_list.append(["Tekay"])
    if argv1=="stepall" or argv1=="step7_readme":
        print('Printing instructions')
        with open('./chr2EDTA/7_stat/readme','w') as f:
            txt=r'''

            '''
            f.write(txt)     
    if argv1=="stepall"  or argv1=="step7" or argv1=="step7.1":  
        print('Count the number of different chromosomes')
        #print('Estimated 167s')
        def run_step(one_task):
            TE_Class=one_task[0]
            
            if  os.path.exists(f'./chr2EDTA/7_stat/{TE_Class}')!=True:
                subprocess.run([f"mkdir ./chr2EDTA/7_stat/{TE_Class}"], shell=True)    
            
            dict_sample_type={}
            dict_chr_sample_num ={}
            with open(f'./chr2EDTA/5_tree/{TE_Class}/0_info_more','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    sample,sample_type,chromosome=eachline_arr[:3]
                    if sample_type=='Routundifolia':continue
                    dict_sample_type[sample]=sample_type
                    samplechr=sample+'___'+chromosome
                    Peri_num=eachline_arr[19]
                    #if Peri_num!='1':continue
                    if chromosome not in dict_chr_sample_num:
                        dict_chr_sample_num[chromosome]={}
                    if sample not in    dict_chr_sample_num[chromosome]:
                        dict_chr_sample_num[chromosome][sample]=0
                    dict_chr_sample_num[chromosome][sample]+=  1
            with open(f'./chr2EDTA/7_stat/{TE_Class}/0_info','w') as f:   
                f.write(f"chromosome\tchromosome_num\tsample\tsample_type\tnum\n")
                for   chromosome,dict_sample_num in dict_chr_sample_num.items():
                    chromosome_num=chromosome.replace('Chr','')
                    for sample , num in dict_sample_num.items():
                        sample_type=dict_sample_type[sample]
                        f.write(f"{chromosome}\t{chromosome_num}\t{sample}\t{sample_type}\t{num}\n")
            R_txt=r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file=read.table('0_info', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    
        # Create plot object
        p <- ggplot() +
            geom_boxplot(
                data = input_file,
                aes(x = as.factor(chromosome_num), y = num),
                width = 0.5,
                alpha = 0.5,outlier.shape = NA
              ) +
              geom_jitter(
                data = input_file,
                aes(x = as.factor(chromosome_num), y = num,color=sample_type),
                size=0.1,
                width = 0.2,  # Jitter width
                alpha = 0.7
              ) 

            p <- p +scale_color_manual(name = "", 
                           values = c("Eurasian" = "#ff3399",   
                                        "Table" = "#E5E54C",   
                                      "Wine" = "#7e318e",  
                                      "East_Asia" = "#0066ff",  
                                      "America" = "#29a329",
                                      "hybrid" = "#d9d9d9",
                                      "Routundifolia"='black'
                                      ))  +
        theme_classic() +         
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
         
          
          # Save as PDF
          pdf(file = paste0('chr', ".pdf"), width = 20/ 2.54, height = 15 / 2.54)
          print(p)
          dev.off()
        
        }
            '''
            with open(f'./chr2EDTA/7_stat/{TE_Class}/plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)  
            new_directory = f'./chr2EDTA/7_stat/{TE_Class}/'
            os.chdir(new_directory)
            subprocess.run(['Rscript plot.R'], shell=True)    
            os.chdir('../../../')                      
        
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, task_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(task_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()
        
        

def merge_intervals(intervals):
    """Merge overlapping intervals"""
    if not intervals:
        return []
    # Sort by start position
    intervals.sort(key=lambda x: x[0])
    merged = [list(intervals[0])]
    for current_start, current_end in intervals[1:]:
        last_start, last_end = merged[-1]
        if current_start <= last_end:  # Overlap exists
            merged[-1][1] = max(last_end, current_end)
        else:
            merged.append([current_start, current_end])
    return merged
    
#PN40024，LIR
if "part1" in argv1:
    if argv1 == "part1":     
        if os.path.exists('./chr2EDTA/part1') == False:
            subprocess.run(["mkdir ./chr2EDTA/part1"], shell=True)
    
        target_list = []
        
        target_list.append(['PN40024','Chr1',15154005,15435737])
        target_list.append(['PN40024','Chr13',12020240,12180415])
        target_list.append(['PN40024','Chr7',13458101,13762988]) 
        target_list.append(['PN40024','Chr14',15144034,15480228])
        target_list.append(['PN40024','Chr17',15221832,15560758])
        target_list.append(['PN40024','Chr2',12939361,13300832])
        target_list.append(['PN40024_hap1','Chr3',14353397,14553397])
        target_list.append(['PN40024','Chr4',13291174,13483935])
        target_list.append(['PN40024','Chr5',14008287,14357682])
        target_list.append(['PN40024','Chr6',10606684,10966038])
        target_list.append(['PN40024','Chr8',6679384,7071344])
        target_list.append(['PN40024','Chr9',15113362,15449399])
        target_list.append(['PN40024','Chr19',15995936,16195936])
        target_list.append(['PN40024','Chr18',16390000,16620000])    

        i = 0
        with open('./chr2EDTA/part1/background_length', 'w') as f:
            f.write(f"i\tlength\tplot_y\n")     
            for one in target_list:
                length = one[3] - one[2] + 1
                f.write(f"{i}\t{length}\t{-i}\n")        
                i += 1
        
        dict_i_type_len = {}
        dict_i_pos_info = {}
        dict_i_deltastart = {}
        
        with open("./chr2EDTA/part1/sum", 'w') as f2:
            f2.write(f"serial\tsample\tchromosome\ttype\tstart\tend\tstrand\tdelta_start\tdelta_end\tplot_y\tintact_mark\n")
            i = 0
            for one in target_list:
                print(one)
                sample, chromosome, start, end = one
                dict_i_deltastart[i] = start
                
                length = end - start
                unit = 40000
                
                dict_i_pos_info[i] = {}
                current_start = 0
                while current_start < length:
                    region_start = current_start
                    region_end = min(current_start + unit, length)
                    dict_i_pos_info[i][f"{region_start}-{region_end}"] = {
                        'te_intervals': []  # 存储窗口内的TE区间用于合并
                    }
                    current_start += unit
                
                # 读取intact标记
                dict_intact_pos = {}
                try:
                    with open(f"./chr2EDTA/1_EDTA/{sample}/{sample}.fasta.mod.EDTA.intact.gff3", 'r') as f:
                        next(f)
                        for line in f:
                            eachline_arr = line.strip().split('\t') 
                            if len(eachline_arr) != 9: continue
                            one_chromosome = eachline_arr[0]
                            one_type = eachline_arr[2] 
                            one_start = int(eachline_arr[3])
                            one_end = int(eachline_arr[4])
                            if one_chromosome not in dict_intact_pos:
                                dict_intact_pos[one_chromosome] = set()
                            dict_intact_pos[one_chromosome].add(f"{one_type}_{one_start}_{one_end}")
                except FileNotFoundError:
                    print(f"Warning: intact file not found for {sample}")
                
                # 收集当前区域内所有TE的相对坐标区间
                all_te_intervals = []  # 每个元素: [rel_start, rel_end, te_key]
                
                try:
                    with open(f"./chr2EDTA/1_EDTA/{sample}/{sample}.fasta.mod.EDTA.TEanno.gff3", 'r') as f:
                        next(f)
                        for line in f:
                            eachline_arr = line.strip().split('\t')
                            if len(eachline_arr) != 9: continue
                            one_chromosome = eachline_arr[0]
                            one_type = eachline_arr[2] 
                            one_start = int(eachline_arr[3])
                            one_end = int(eachline_arr[4])
                            if one_end-one_start>30000:continue  
                            one_strand = eachline_arr[5]
                            
                            if one_type not in ["Gypsy_LTR_retrotransposon", "LTR_retrotransposon", 
                                                "PIF_Harbinger_TIR_transposon", "Copia_LTR_retrotransposon",
                                                "CACTA_TIR_transposon", "hAT_TIR_transposon", "helitron",
                                                "Mutator_TIR_transposon", "Tc1_Mariner_TIR_transposon"]:
                                continue
                            
                            if one_chromosome == chromosome:
                                # 计算与目标区域的重叠
                                overlap_start = max(one_start, start)
                                overlap_end = min(one_end, end)
                                
                                if overlap_start < overlap_end:
                                    overlap_length = overlap_end - overlap_start + 1
                                    total_te_length = one_end - one_start + 1
                                    overlap_ratio = overlap_length / total_te_length
                                    
                                    te_key = f"{one_type}_{one_start}_{one_end}"
                                    
                                    # 判断intact状态
                                    intact_mark = 'intact' if te_key in dict_intact_pos.get(one_chromosome, set()) else 'broken'
                                    
                                    # 写入sum文件（保持原有格式）
                                    delta_start = one_start - start
                                    delta_end = one_end - start
                                    f2.write(f"{i}\t{sample}\t{chromosome}\t{one_type}\t{one_start}\t{one_end}\t{one_strand}\t{delta_start}\t{delta_end}\t{-i}\t{intact_mark}\n")
                                    
                                    # 统计每种类型的总长度（按重叠长度计算）
                                    if i not in dict_i_type_len:
                                        dict_i_type_len[i] = {}
                                    if one_type not in dict_i_type_len[i]:
                                        dict_i_type_len[i][one_type] = 0
                                    dict_i_type_len[i][one_type] += overlap_length
                                    
                                    # 存储相对坐标区间（用于窗口合并计算）
                                    rel_start = overlap_start - start
                                    rel_end = overlap_end - start
                                    all_te_intervals.append([rel_start, rel_end])
                except FileNotFoundError:
                    print(f"Warning: TEanno file not found for {sample}")
                
                # 按窗口处理，合并重叠区间后计算每个窗口的覆盖长度
                for pos, dict_info in dict_i_pos_info[i].items():
                    region_range = pos.split('-')
                    region_start = int(region_range[0])
                    region_end = int(region_range[1])
                    
                    # 找出当前窗口内的所有TE区间
                    window_intervals = []
                    for te_rel_start, te_rel_end in all_te_intervals:
                        overlap_start = max(te_rel_start, region_start)
                        overlap_end = min(te_rel_end, region_end)
                        if overlap_start < overlap_end:
                            window_intervals.append([overlap_start, overlap_end])
                    
                    # 合并重叠区间
                    merged = merge_intervals(window_intervals)
                    
                    # 计算合并后的总长度
                    total_len = sum(m_end - m_start for m_start, m_end in merged)
                    
                    # 存储结果（保持原有格式的键名）
                    dict_info['num'] = len(window_intervals)  # TE片段数量
                    dict_info['len'] = total_len  # 合并后的覆盖长度（不会超过窗口大小）
                
                i += 1
        
        print('输出')                        
        with open("./chr2EDTA/part1/sum_stat", 'w') as f:
            f.write(f"i\ttype\tlength\n")
            for i, dict_type_len in dict_i_type_len.items():
                for one_type, length in dict_type_len.items():
                    f.write(f"{i}\t{one_type}\t{length}\n")
        
        print(dict_i_deltastart)
        
        max_num=0
        max_len=0
        # 保持原有输出格式
        with open("./chr2EDTA/part1/sum_density", 'w') as f:  
            f.write(f"i\tregion_start\tregion_end\tnum\tlength\n")
            for i, dict_pos_info in dict_i_pos_info.items():
                for pos, info in dict_pos_info.items():
                    region = pos.split('-')
                    region_start = int(region[0])
                    region_end = int(region[1])
                    num = info['num']
                    length = info['len']
                    if length == 0:
                        continue
                    f.write(f"{i}\t{region_start}\t{region_end}\t{num}\t{length}\n")
                    if num>max_num:max_num=num
                    if length>max_len:max_len=length
        with open("./chr2EDTA/part1/sum_density_info",'w') as f:
            f.write(f'window:\t{unit}\n')
            f.write(f'max_num:\t{max_num}\n')
            f.write(f'max_length:\t{max_len}\n')
    if argv1=="part1" or argv1=="part1p1":
        R_txt=r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file=read.table('sum', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
          input_file2=read.table('background_length', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
        # 过滤数据
        intact_part=input_file%>%filter(intact_mark=="intact")
    
        
        # 定义颜色值
        color_values <- c(
            "Gypsy_LTR_retrotransposon"="#f7776c",
            "Copia_LTR_retrotransposon"="#cc9702",
            "LTR_retrotransposon"="#7eb008",
    
            "Tc1_Mariner_TIR_transposon"="#00bf67",
            "Mutator_TIR_transposon"="#01c0c3",
            "hAT_TIR_transposon"="#6c89c0",
            "helitron"="#02aafd",
            "CACTA_TIR_transposon"="#c87eff",
            
            "PIF_Harbinger_TIR_transposon"="#fe63cb"
        )
         
        # 创建绘图对象
        p <- ggplot() +
            geom_rect(data = input_file2, aes(xmin = 0, xmax = length, ymin = plot_y-0.1, ymax = plot_y+0.1), fill = 'black')+
            geom_rect(data = input_file, aes(xmin = delta_start, xmax = delta_end, ymin = plot_y-0.1, ymax = plot_y+0.1, fill = type))+
            
            geom_rect(data = intact_part, aes(xmin = delta_start, xmax = delta_end, ymin = plot_y-0.3, ymax = plot_y-0.2), fill = 'red')+
            
            scale_fill_manual(values = color_values, drop = FALSE)+
        theme_classic() +         
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
         
          
          # 保存为 PDF
          pdf(file = paste0('plot', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
          print(p)
          dev.off()
        
        }
            '''
        with open('./chr2EDTA/part1/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2EDTA/part1//'
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')                             
    if argv1=="part1" or argv1=="part1p2":
        R_txt=r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file=read.table('sum_stat', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
        
        # 过滤数据
    
    
        
        # 定义颜色值
        color_values <- c(
            "Gypsy_LTR_retrotransposon"="#f7776c",
            "Copia_LTR_retrotransposon"="#cc9702",
            "LTR_retrotransposon"="#7eb008",
    
            "Tc1_Mariner_TIR_transposon"="#00bf67",
            "Mutator_TIR_transposon"="#01c0c3",
            "hAT_TIR_transposon"="#6c89c0",
            "helitron"="#02aafd",
            "CACTA_TIR_transposon"="#c87eff",
            
            "PIF_Harbinger_TIR_transposon"="#fe63cb"
        )
         
        # 创建绘图对象
        p <- ggplot() +
            geom_col(data = input_file, aes(x=i, y = length, fill = type),position = "stack")+
            scale_fill_manual(values = color_values, drop = FALSE)+
        theme_classic() +         
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
         
          
          # 保存为 PDF
          pdf(file = paste0('plot_stat', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
          print(p)
          dev.off()
        
        }
            '''
        with open('./chr2EDTA/part1/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2EDTA/part1//'
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')                                                     
    if argv1=="part1" or argv1=="part1p3":
        R_txt=r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file=read.table('sum_density', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    
        # 创建绘图对象
        p <- ggplot() 
            p <- p + geom_rect(data = input_file, aes(xmin = region_start, xmax = region_end, ymin = 0, ymax = num, fill = length))
            p <- p + facet_wrap(~ i, ncol = 100)+
           # scale_fill_manual(values = color_values, drop = FALSE)+
scale_fill_gradientn(
  colors = c("#0D4E5C", "#2E7A6B", "#D98A3C", "#D98A3C"),
  values = c(0, 0.4, 0.8, 1)
)+
        theme_classic() +         
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
         
          
          # 保存为 PDF
          pdf(file = paste0('plot_density', ".pdf"), width = 20/ 2.54, height = 5 / 2.54)
          print(p)
          dev.off()
        
        }
            '''
        with open('./chr2EDTA/part1/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2EDTA/part1//'
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')  
    if argv1=="part1" or argv1=="part1p4" and 1==2:
        print('试图融合两个指标，效果一般，算了')
        import pandas as pd
        # 读取你的数据
        df = pd.read_csv("./chr2EDTA/part1/sum_density", sep='\t')            
        df['num_norm'] = (df['num'] - df['num'].min()) / (df['num'].max() - df['num'].min())
        df['length_norm'] = (df['length'] - df['length'].min()) / (df['length'].max() - df['length'].min())
        df['combined_norm'] = df['num_norm'] * df['length_norm']
        # 保存包含新指标的表
        df.to_csv("./chr2EDTA/part1/sum_density_combined", sep='\t', index=False)           
        R_txt=r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file=read.table('sum_density_combined', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    	
        # 创建绘图对象
        p <- ggplot() 
            p <- p + geom_rect(data = input_file, aes(xmin = region_start, xmax = region_end, ymin = 0, ymax = combined_norm, fill = combined_norm))
            p <- p + facet_wrap(~ i, ncol = 100)+
           # scale_fill_manual(values = color_values, drop = FALSE)+
scale_fill_gradientn(
  colors = c("#1A237E", "#4A148C", "#AD1457", "#AD1457"),
  values = c(0, 0.2, 0.5, 1)
)+
        theme_classic() +         
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
         
          
          # 保存为 PDF
          pdf(file = paste0('plot_density_norm', ".pdf"), width = 20/ 2.54, height = 5 / 2.54)
          print(p)
          dev.off()
        
        }
            '''
        with open('./chr2EDTA/part1/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2EDTA/part1//'
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')  

#PN40024,Chr10,11,15,16        
if "part2" in argv1:   
    if argv1=="part2":     
        if  os.path.exists('./chr2EDTA/part2')==False:
            subprocess.run(["mkdir ./chr2EDTA/part2"], shell=True)
    
    
        target_list=[]
        target_list.append(['PN40024','Chr10',19282777,22365281])        
        target_list.append(['PN40024','Chr11',12682161,14996086])   
        target_list.append(['PN40024','Chr15',7025767,9336095])       ########默认宽度太窄了
        #target_list.append(['PN40024','Chr15',6500000,10000000]) 
        target_list.append(['PN40024','Chr16',7953854,14385543])  
        
        '''target_list.append(['PN40024','Chr10',19282777,20196420])        
        target_list.append(['PN40024','Chr11',12682161,13446851])   
        target_list.append(['PN40024','Chr15',8965378,9336095])       
        target_list.append(['PN40024','Chr16',12780139,14385543])  '''
     
        i = 0
        with open('./chr2EDTA/part2/background_length', 'w') as f:
            f.write(f"i\tlength\tplot_y\n")     
            for one in target_list:
                length = one[3] - one[2] + 1
                f.write(f"{i}\t{length}\t{-i}\n")        
                i += 1
        
        dict_i_type_len = {}
        dict_i_pos_info = {}
        dict_i_deltastart = {}
        
        with open("./chr2EDTA/part2/sum", 'w') as f2:
            f2.write(f"serial\tsample\tchromosome\ttype\tstart\tend\tstrand\tdelta_start\tdelta_end\tplot_y\tintact_mark\n")
            i = 0
            for one in target_list:
                print(one)
                sample, chromosome, start, end = one
                dict_i_deltastart[i] = start
                
                length = end - start
                unit = 40000
                
                dict_i_pos_info[i] = {}
                current_start = 0
                while current_start < length:
                    region_start = current_start
                    region_end = min(current_start + unit, length)
                    dict_i_pos_info[i][f"{region_start}-{region_end}"] = {
                        'te_intervals': []  # 存储窗口内的TE区间用于合并
                    }
                    current_start += unit
                
                # 读取intact标记
                dict_intact_pos = {}
                try:
                    with open(f"./chr2EDTA/1_EDTA/{sample}/{sample}.fasta.mod.EDTA.intact.gff3", 'r') as f:
                        next(f)
                        for line in f:
                            eachline_arr = line.strip().split('\t') 
                            if len(eachline_arr) != 9: continue
                            one_chromosome = eachline_arr[0]
                            one_type = eachline_arr[2] 
                            one_start = int(eachline_arr[3])
                            one_end = int(eachline_arr[4])
                            if one_chromosome not in dict_intact_pos:
                                dict_intact_pos[one_chromosome] = set()
                            dict_intact_pos[one_chromosome].add(f"{one_type}_{one_start}_{one_end}")
                except FileNotFoundError:
                    print(f"Warning: intact file not found for {sample}")
                
                # 收集当前区域内所有TE的相对坐标区间
                all_te_intervals = []  # 每个元素: [rel_start, rel_end, te_key]
                
                try:
                    with open(f"./chr2EDTA/1_EDTA/{sample}/{sample}.fasta.mod.EDTA.TEanno.gff3", 'r') as f:
                        next(f)
                        for line in f:
                            eachline_arr = line.strip().split('\t')
                            if len(eachline_arr) != 9: continue
                            one_chromosome = eachline_arr[0]
                            one_type = eachline_arr[2] 
                            one_start = int(eachline_arr[3])
                            one_end = int(eachline_arr[4])
                            one_strand = eachline_arr[5]
                            if one_end-one_start>30000:continue  
                            if one_type not in ["Gypsy_LTR_retrotransposon", "LTR_retrotransposon", 
                                                "PIF_Harbinger_TIR_transposon", "Copia_LTR_retrotransposon",
                                                "CACTA_TIR_transposon", "hAT_TIR_transposon", "helitron",
                                                "Mutator_TIR_transposon", "Tc1_Mariner_TIR_transposon"]:
                                continue
                            
                            if one_chromosome == chromosome:
                                # 计算与目标区域的重叠
                                overlap_start = max(one_start, start)
                                overlap_end = min(one_end, end)
                                
                                if overlap_start < overlap_end:
                                    overlap_length = overlap_end - overlap_start + 1
                                    total_te_length = one_end - one_start + 1
                                    overlap_ratio = overlap_length / total_te_length
                                    
                                    te_key = f"{one_type}_{one_start}_{one_end}"
                                    
                                    # 判断intact状态
                                    intact_mark = 'intact' if te_key in dict_intact_pos.get(one_chromosome, set()) else 'broken'
                                    
                                    # 写入sum文件（保持原有格式）
                                    delta_start = one_start - start
                                    delta_end = one_end - start
                                    f2.write(f"{i}\t{sample}\t{chromosome}\t{one_type}\t{one_start}\t{one_end}\t{one_strand}\t{delta_start}\t{delta_end}\t{-i}\t{intact_mark}\n")
                                    
                                    # 统计每种类型的总长度（按重叠长度计算）
                                    if i not in dict_i_type_len:
                                        dict_i_type_len[i] = {}
                                    if one_type not in dict_i_type_len[i]:
                                        dict_i_type_len[i][one_type] = 0
                                    dict_i_type_len[i][one_type] += overlap_length
                                    
                                    # 存储相对坐标区间（用于窗口合并计算）
                                    rel_start = overlap_start - start
                                    rel_end = overlap_end - start
                                    all_te_intervals.append([rel_start, rel_end])
                except FileNotFoundError:
                    print(f"Warning: TEanno file not found for {sample}")
                
                # 按窗口处理，合并重叠区间后计算每个窗口的覆盖长度
                for pos, dict_info in dict_i_pos_info[i].items():
                    region_range = pos.split('-')
                    region_start = int(region_range[0])
                    region_end = int(region_range[1])
                    
                    # 找出当前窗口内的所有TE区间
                    window_intervals = []
                    for te_rel_start, te_rel_end in all_te_intervals:
                        overlap_start = max(te_rel_start, region_start)
                        overlap_end = min(te_rel_end, region_end)
                        if overlap_start < overlap_end:
                            window_intervals.append([overlap_start, overlap_end])
                    
                    # 合并重叠区间
                    merged = merge_intervals(window_intervals)
                    
                    # 计算合并后的总长度
                    total_len = sum(m_end - m_start for m_start, m_end in merged)
                    
                    # 存储结果（保持原有格式的键名）
                    dict_info['num'] = len(window_intervals)  # TE片段数量
                    dict_info['len'] = total_len  # 合并后的覆盖长度（不会超过窗口大小）
                
                i += 1
        
        print('输出')                        
        with open("./chr2EDTA/part2/sum_stat", 'w') as f:
            f.write(f"i\ttype\tlength\n")
            for i, dict_type_len in dict_i_type_len.items():
                for one_type, length in dict_type_len.items():
                    f.write(f"{i}\t{one_type}\t{length}\n")
        
        print(dict_i_deltastart)
        
        max_num=0
        max_len=0
        # 保持原有输出格式
        with open("./chr2EDTA/part2/sum_density", 'w') as f:  
            f.write(f"i\tregion_start\tregion_end\tnum\tlength\n")
            for i, dict_pos_info in dict_i_pos_info.items():
                for pos, info in dict_pos_info.items():
                    region = pos.split('-')
                    region_start = int(region[0])
                    region_end = int(region[1])
                    num = info['num']
                    length = info['len']
                    if length == 0:
                        continue
                    f.write(f"{i}\t{region_start}\t{region_end}\t{num}\t{length}\n")
                    if num>max_num:max_num=num
                    if length>max_len:max_len=length
        with open("./chr2EDTA/part2/sum_density_info",'w') as f:
            f.write(f'window:\t{unit}\n')
            f.write(f'max_num:\t{max_num}\n')
            f.write(f'max_length:\t{max_len}\n')
    if argv1=="part2" or argv1=="part2p1":
        R_txt=r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file=read.table('sum', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
          input_file2=read.table('background_length', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
        # 过滤数据
        intact_part=input_file%>%filter(intact_mark=="intact")
    
        
        # 定义颜色值
        color_values <- c(
            "Gypsy_LTR_retrotransposon"="#f7776c",
            "Copia_LTR_retrotransposon"="#cc9702",
            "LTR_retrotransposon"="#7eb008",
    
            "Tc1_Mariner_TIR_transposon"="#00bf67",
            "Mutator_TIR_transposon"="#01c0c3",
            "hAT_TIR_transposon"="#6c89c0",
            "helitron"="#02aafd",
            "CACTA_TIR_transposon"="#c87eff",
            
            "PIF_Harbinger_TIR_transposon"="#fe63cb"
        )
         
        # 创建绘图对象
        p <- ggplot() +
            geom_rect(data = input_file2, aes(xmin = 0, xmax = length, ymin = plot_y-0.1, ymax = plot_y+0.1), fill = 'black')+
            geom_rect(data = input_file, aes(xmin = delta_start, xmax = delta_end, ymin = plot_y-0.1, ymax = plot_y+0.1, fill = type))+
            
            geom_rect(data = intact_part, aes(xmin = delta_start, xmax = delta_end, ymin = plot_y-0.3, ymax = plot_y-0.2), fill = 'red')+
            
            scale_fill_manual(values = color_values, drop = FALSE)+
        theme_classic() +         
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
         
          
          # 保存为 PDF
          pdf(file = paste0('plot', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
          print(p)
          dev.off()
        
        }
            '''
        with open('./chr2EDTA/part2/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2EDTA/part2//'
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')                             
    if argv1=="part2" or argv1=="part2p2":
        R_txt=r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file=read.table('sum_stat', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
        
        # 过滤数据
    
    
        
        # 定义颜色值
        color_values <- c(
            "Gypsy_LTR_retrotransposon"="#f7776c",
            "Copia_LTR_retrotransposon"="#cc9702",
            "LTR_retrotransposon"="#7eb008",
    
            "Tc1_Mariner_TIR_transposon"="#00bf67",
            "Mutator_TIR_transposon"="#01c0c3",
            "hAT_TIR_transposon"="#6c89c0",
            "helitron"="#02aafd",
            "CACTA_TIR_transposon"="#c87eff",
            
            "PIF_Harbinger_TIR_transposon"="#fe63cb"
        )
         
        # 创建绘图对象
        p <- ggplot() +
            geom_col(data = input_file, aes(x=i, y = length, fill = type),position = "stack")+
            scale_fill_manual(values = color_values, drop = FALSE)+
        theme_classic() +         
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
         
          
          # 保存为 PDF
          pdf(file = paste0('plot_stat', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
          print(p)
          dev.off()
        
        }
            '''
        with open('./chr2EDTA/part2/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2EDTA/part2//'
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')                                                     
    if argv1=="part2" or argv1=="part2p3":
        R_txt=r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file=read.table('sum_density', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    
        # 创建绘图对象
        p <- ggplot() 
            p <- p + geom_rect(data = input_file, aes(xmin = region_start, xmax = region_end, ymin = 0, ymax = num, fill = length))
            p <- p + facet_wrap(~ i, ncol = 100)+
            #scale_fill_manual(values = color_values, drop = FALSE)+
        theme_classic() +         
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
         
          
          # 保存为 PDF
          pdf(file = paste0('plot_density', ".pdf"), width = 20/ 2.54, height = 5 / 2.54)
          print(p)
          dev.off()
        
        # 创建绘图对象
        p <- ggplot() 
        p <- p + geom_rect(data = input_file, aes(xmin = region_start, xmax = region_end, ymin = 0, ymax = num, fill = length))
            p <- p + facet_wrap(~ i, ncol = 100)+

scale_fill_gradientn(
  colors = c("#0D4E5C", "#2E7A6B", "#D98A3C", "#D98A3C"),
  values = c(0, 0.4, 0.8, 1)
)+
theme_classic() +         
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
         
          
          # 保存为 PDF
          pdf(file = paste0('plot_density_2', ".pdf"), width = 20/ 2.54, height = 5 / 2.54)
          print(p)
          dev.off()
        }
            '''
        with open('./chr2EDTA/part2/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2EDTA/part2//'
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')                                                     
    if argv1=="part2"  or argv1=="part2.2":
        dict_i_length={}
        with open ('./chr2EDTA/part2/sum_stat','r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                i,one_type,length=eachline_arr
                if i not in dict_i_length:dict_i_length[i]=0
                dict_i_length[i]+=int(length)
        with open ('./chr2EDTA/part2/sum_stat100','w') as f2:  
            f2.write(f"i\ttype\tlength\tpercent\n")
            with open ('./chr2EDTA/part2/sum_stat','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    i,one_type,length=eachline_arr
                    percent=round(int(length)/dict_i_length[i],3)
                    f2.write(f"{i}\t{one_type}\t{length}\t{percent}\n")
        R_txt=r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file=read.table('sum_stat100', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
        
        # 过滤数据
    
    
        
        # 定义颜色值
        color_values <- c(
            "Gypsy_LTR_retrotransposon"="#f7776c",
            "Copia_LTR_retrotransposon"="#cc9702",
            "LTR_retrotransposon"="#7eb008",
    
            "Tc1_Mariner_TIR_transposon"="#00bf67",
            "Mutator_TIR_transposon"="#01c0c3",
            "hAT_TIR_transposon"="#6c89c0",
            "helitron"="#02aafd",
            "CACTA_TIR_transposon"="#c87eff",
            
            "PIF_Harbinger_TIR_transposon"="#fe63cb"
        )
        input_file$type <- factor(
          input_file$type,
          levels = c("Gypsy_LTR_retrotransposon", "Copia_LTR_retrotransposon", "LTR_retrotransposon","Tc1_Mariner_TIR_transposon","Mutator_TIR_transposon","hAT_TIR_transposon","helitron","CACTA_TIR_transposon", "PIF_Harbinger_TIR_transposon")
        ) 
        # 创建绘图对象
        p <- ggplot() +
            geom_col(data = input_file, aes(x=i, y = percent, fill = type),width=0.5,position = "stack")+
            scale_fill_manual(values = color_values, drop = FALSE)+
        theme_classic() +         
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
         
          
          # 保存为 PDF
          pdf(file = paste0('plot_stat100', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
          print(p)
          dev.off()
        
        }
            '''
        with open('./chr2EDTA/part2/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2EDTA/part2//'
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')                               
   
        
##chr9  
if "part4" in argv1:
    if argv1 == "part4":     
        if os.path.exists('./chr2EDTA/part4') == False:
            subprocess.run(["mkdir ./chr2EDTA/part4"], shell=True)
    
        target_list = []
        
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
        
        i = 0
        with open('./chr2EDTA/part4/background_length', 'w') as f:
            f.write(f"i\tlength\tplot_y\n")     
            for one in target_list:
                length = one[3] - one[2] + 1
                f.write(f"{i}\t{length}\t{-i}\n")        
                i += 1
        
        dict_i_type_len = {}
        dict_i_pos_info = {}
        dict_i_deltastart = {}
        
        with open("./chr2EDTA/part4/sum", 'w') as f2:
            f2.write(f"serial\tsample\tchromosome\ttype\tstart\tend\tstrand\tdelta_start\tdelta_end\tplot_y\tintact_mark\n")
            i = 0
            for one in target_list:
                print(one)
                sample, chromosome, start, end = one
                dict_i_deltastart[i] = start
                
                length = end - start
                unit = 40000
                
                dict_i_pos_info[i] = {}
                current_start = 0
                while current_start < length:
                    region_start = current_start
                    region_end = min(current_start + unit, length)
                    dict_i_pos_info[i][f"{region_start}-{region_end}"] = {
                        'te_intervals': []  
                    }
                    current_start += unit
                
                dict_intact_pos = {}
                try:
                    with open(f"./chr2EDTA/1_EDTA/{sample}/{sample}.fasta.mod.EDTA.intact.gff3", 'r') as f:
                        next(f)
                        for line in f:
                            eachline_arr = line.strip().split('\t') 
                            if len(eachline_arr) != 9: continue
                            one_chromosome = eachline_arr[0]
                            one_type = eachline_arr[2] 
                            one_start = int(eachline_arr[3])
                            one_end = int(eachline_arr[4])
                            
                            if one_chromosome not in dict_intact_pos:
                                dict_intact_pos[one_chromosome] = set()
                            dict_intact_pos[one_chromosome].add(f"{one_type}_{one_start}_{one_end}")
                except FileNotFoundError:
                    print(f"Warning: intact file not found for {sample}")
                

                all_te_intervals = []  
                
                try:
                    with open(f"./chr2EDTA/1_EDTA/{sample}/{sample}.fasta.mod.EDTA.TEanno.gff3", 'r') as f:
                        next(f)
                        for line in f:
                            eachline_arr = line.strip().split('\t')
                            if len(eachline_arr) != 9: continue
                            one_chromosome = eachline_arr[0]
                            one_type = eachline_arr[2] 
                            one_start = int(eachline_arr[3])
                            one_end = int(eachline_arr[4])
                            if one_end-one_start>30000:continue  
                            one_strand = eachline_arr[5]
                            
                            if one_type not in ["Gypsy_LTR_retrotransposon", "LTR_retrotransposon", 
                                                "PIF_Harbinger_TIR_transposon", "Copia_LTR_retrotransposon",
                                                "CACTA_TIR_transposon", "hAT_TIR_transposon", "helitron",
                                                "Mutator_TIR_transposon", "Tc1_Mariner_TIR_transposon"]:
                                continue
                            
                            if one_chromosome == chromosome:
                                # 计算与目标区域的重叠
                                overlap_start = max(one_start, start)
                                overlap_end = min(one_end, end)
                                
                                if overlap_start < overlap_end:
                                    overlap_length = overlap_end - overlap_start + 1
                                    total_te_length = one_end - one_start + 1
                                    overlap_ratio = overlap_length / total_te_length
                                    
                                    te_key = f"{one_type}_{one_start}_{one_end}"
                                    
                                    intact_mark = 'intact' if te_key in dict_intact_pos.get(one_chromosome, set()) else 'broken'
                                    
                                    delta_start = one_start - start
                                    delta_end = one_end - start
                                    f2.write(f"{i}\t{sample}\t{chromosome}\t{one_type}\t{one_start}\t{one_end}\t{one_strand}\t{delta_start}\t{delta_end}\t{-i}\t{intact_mark}\n")
                                    
                                    if i not in dict_i_type_len:
                                        dict_i_type_len[i] = {}
                                    if one_type not in dict_i_type_len[i]:
                                        dict_i_type_len[i][one_type] = 0
                                    dict_i_type_len[i][one_type] += overlap_length

                                    rel_start = overlap_start - start
                                    rel_end = overlap_end - start
                                    all_te_intervals.append([rel_start, rel_end])
                except FileNotFoundError:
                    print(f"Warning: TEanno file not found for {sample}")

                for pos, dict_info in dict_i_pos_info[i].items():
                    region_range = pos.split('-')
                    region_start = int(region_range[0])
                    region_end = int(region_range[1])
                    
                    window_intervals = []
                    for te_rel_start, te_rel_end in all_te_intervals:
                        overlap_start = max(te_rel_start, region_start)
                        overlap_end = min(te_rel_end, region_end)
                        if overlap_start < overlap_end:
                            window_intervals.append([overlap_start, overlap_end])
                    
                    merged = merge_intervals(window_intervals)
                    
                    total_len = sum(m_end - m_start for m_start, m_end in merged)
                    
                    dict_info['num'] = len(window_intervals)  
                    dict_info['len'] = total_len  
                
                i += 1
        
                     
        with open("./chr2EDTA/part4/sum_stat", 'w') as f:
            f.write(f"i\ttype\tlength\n")
            for i, dict_type_len in dict_i_type_len.items():
                for one_type, length in dict_type_len.items():
                    f.write(f"{i}\t{one_type}\t{length}\n")
        
        print(dict_i_deltastart)
        
        max_num=0
        max_len=0
        # 保持原有输出格式
        with open("./chr2EDTA/part4/sum_density", 'w') as f:  
            f.write(f"i\tregion_start\tregion_end\tnum\tlength\n")
            for i, dict_pos_info in dict_i_pos_info.items():
                for pos, info in dict_pos_info.items():
                    region = pos.split('-')
                    region_start = int(region[0])
                    region_end = int(region[1])
                    num = info['num']
                    length = info['len']
                    if length == 0:
                        continue
                    f.write(f"{i}\t{region_start}\t{region_end}\t{num}\t{length}\n")
                    if num>max_num:max_num=num
                    if length>max_len:max_len=length
        with open("./chr2EDTA/part4/sum_density_info",'w') as f:
            f.write(f'window:\t{unit}\n')
            f.write(f'max_num:\t{max_num}\n')
            f.write(f'max_length:\t{max_len}\n')
    if argv1=="part4" or argv1=="part4p1":
        R_txt=r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file=read.table('sum', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
          input_file2=read.table('background_length', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
        # 过滤数据
        intact_part=input_file%>%filter(intact_mark=="intact")
    
        
        # 定义颜色值
        color_values <- c(
            "Gypsy_LTR_retrotransposon"="#f7776c",
            "Copia_LTR_retrotransposon"="#cc9702",
            "LTR_retrotransposon"="#7eb008",
    
            "Tc1_Mariner_TIR_transposon"="#00bf67",
            "Mutator_TIR_transposon"="#01c0c3",
            "hAT_TIR_transposon"="#6c89c0",
            "helitron"="#02aafd",
            "CACTA_TIR_transposon"="#c87eff",
            
            "PIF_Harbinger_TIR_transposon"="#fe63cb"
        )
         
        # 创建绘图对象
        p <- ggplot() +
            geom_rect(data = input_file2, aes(xmin = 0, xmax = length, ymin = plot_y-0.1, ymax = plot_y+0.1), fill = 'black')+
            geom_rect(data = input_file, aes(xmin = delta_start, xmax = delta_end, ymin = plot_y-0.1, ymax = plot_y+0.1, fill = type))+
            
            geom_rect(data = intact_part, aes(xmin = delta_start, xmax = delta_end, ymin = plot_y-0.3, ymax = plot_y-0.2), fill = 'red')+
            
            scale_fill_manual(values = color_values, drop = FALSE)+
        theme_classic() +         
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
         
          
          # 保存为 PDF
          pdf(file = paste0('plot', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
          print(p)
          dev.off()
        
        }
            '''
        with open('./chr2EDTA/part4/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2EDTA/part4//'
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')                             
    if argv1=="part4" or argv1=="part4p2":
        R_txt=r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file=read.table('sum_stat', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
        
        # 过滤数据
    # 按照color_values的顺序设置type为因子
type_order <- c(
    "Gypsy_LTR_retrotransposon",
    "Copia_LTR_retrotransposon",
    "LTR_retrotransposon",
    "Tc1_Mariner_TIR_transposon",
    "Mutator_TIR_transposon",
    "hAT_TIR_transposon",
    "helitron",
    "CACTA_TIR_transposon",
    "PIF_Harbinger_TIR_transposon"
)

input_file$type <- factor(input_file$type, levels = type_order)
    
        
        # 定义颜色值
        color_values <- c(
            "Gypsy_LTR_retrotransposon"="#f7776c",
            "Copia_LTR_retrotransposon"="#cc9702",
            "LTR_retrotransposon"="#7eb008",
    
            "Tc1_Mariner_TIR_transposon"="#00bf67",
            "Mutator_TIR_transposon"="#01c0c3",
            "hAT_TIR_transposon"="#6c89c0",
            "helitron"="#02aafd",
            "CACTA_TIR_transposon"="#c87eff",
            
            "PIF_Harbinger_TIR_transposon"="#fe63cb"
        )
         
        # 创建绘图对象
        p <- ggplot() +
            geom_col(data = input_file, aes(x=i, y = length, fill = type),position = "stack")+
            scale_fill_manual(values = color_values, drop = FALSE)+
        theme_classic() +         
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
         
          
          # 保存为 PDF
          pdf(file = paste0('plot_stat', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
          print(p)
          dev.off()
        
        
# 计算每个i的总长度
library(dplyr)
input_file <- input_file %>%
  group_by(i) %>%
  mutate(total_length = sum(length)) %>%
  ungroup() %>%
  mutate(percentage = length / total_length * 100)

# 定义颜色值
color_values <- c(
    "Gypsy_LTR_retrotransposon"="#f7776c",
    "Copia_LTR_retrotransposon"="#cc9702",
    "LTR_retrotransposon"="#7eb008",

    "Tc1_Mariner_TIR_transposon"="#00bf67",
    "Mutator_TIR_transposon"="#01c0c3",
    "hAT_TIR_transposon"="#6c89c0",
    "helitron"="#02aafd",
    "CACTA_TIR_transposon"="#c87eff",
    
    "PIF_Harbinger_TIR_transposon"="#fe63cb"
)

# 创建百分比堆积图
p <- ggplot() +
    geom_col(data = input_file, aes(x = factor(i), y = percentage, fill = type), position = "stack",width=0.5) +
    scale_fill_manual(values = color_values, drop = FALSE) +
    labs(x = "i", y = "Percentage (%)", title = "TE类型百分比堆积图") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
    )+coord_flip()   # 翻转坐标轴，使图形竖过来

# 保存为 PDF
pdf(file = paste0('plot_stat_percentage', ".pdf"), width = 20 / 2.54, height = 20 / 2.54)
print(p)
dev.off()
        }
            '''
        with open('./chr2EDTA/part4/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2EDTA/part4//'
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')                                                     
    if argv1=="part4" or argv1=="part4p3":
        R_txt=r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file=read.table('sum_density', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    
        # 创建绘图对象
        p <- ggplot() 
            p <- p + geom_rect(data = input_file, aes(xmin = region_start, xmax = region_end, ymin = 0, ymax = num, fill = length))
            p <- p + facet_wrap(~ i, ncol = 100)+
           # scale_fill_manual(values = color_values, drop = FALSE)+
scale_fill_gradientn(
  colors = c("#0D4E5C", "#2E7A6B", "#D98A3C", "#D98A3C"),
  values = c(0, 0.4, 0.8, 1)
)+
        theme_classic() +         
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
         
          
          # 保存为 PDF
          pdf(file = paste0('plot_density', ".pdf"), width = 20/ 2.54, height = 5 / 2.54)
          print(p)
          dev.off()
        
        }
            '''
        with open('./chr2EDTA/part4/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2EDTA/part4//'
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')  

if argv1=="part_allaverage":
    sample_list=[]
    with open('/home/lain/aaa_data/run0/samples_satellite/sample_info','r') as f:
        for line in f:
            eachline_arr=line.strip().split()
            if len(eachline_arr)<3:continue
            sample=eachline_arr[0]
            group=eachline_arr[3]
            if group=='Routundifolia':continue            
            sample_list.append(sample)
    dict_samplechr={}       
    
    
    ##获得染色体臂的位置
    with open('/home/lain/aaa_data/run0/samples_satellite/2_good_regions_main_66_103','r') as f:
        #sample	chromosome	start	end	length
        next(f)
        for line in f:
            eachline_arr=line.strip().split()
            if len(eachline_arr)<5:continue
            sample,chromosome,start,end,length=eachline_arr
            ##着丝粒区域为start-end
            
            
 
    
    def run_step1(one_sample):
        # 存储染色体总长的字典
        chr_length_dict = {}
        with open(f'./chr2EDTA/0_prepare/{one_sample}/{one_sample}.fasta.fai','r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t')
                chr_name = eachline_arr[0]
                chr_len = int(eachline_arr[1])
                chr_length_dict[chr_name] = chr_len
        
        # 存储每个染色体的 TE 区间（过滤后）
        from collections import defaultdict
        chr_intervals = defaultdict(list)        
        
        with open(f"./chr2EDTA/1_EDTA/{one_sample}/{one_sample}.fasta.mod.EDTA.TEanno.gff3",'r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=9:continue
                one_chromosome= eachline_arr[0]
                one_type=       eachline_arr[2] 
                one_start=      int(eachline_arr[3])
                one_end=        int(eachline_arr[4])
                if one_end-one_start>30000:continue
                if one_type not in ["Gypsy_LTR_retrotransposon","LTR_retrotransposon","PIF_Harbinger_TIR_transposon","Copia_LTR_retrotransposon","CACTA_TIR_transposon","hAT_TIR_transposon","helitron","Mutator_TIR_transposon","Tc1_Mariner_TIR_transposon"]:continue            
        
                # 只记录需要统计的 TE
                chr_intervals[one_chromosome].append((one_start, one_end))
        
        # 分染色体统计结果列表
        result_list = []
        
        for chrom, intervals in chr_intervals.items():
            if not intervals:
                # 没有TE的染色体也要记录
                chr_total_len = chr_length_dict.get(chrom, 0)
                result_list.append([one_sample, chrom, chr_total_len, 0, 0])
                continue
        
            # 按起始位置排序
            intervals.sort()
            merged = []
            cur_start, cur_end = intervals[0]
        
            for s, e in intervals[1:]:
                if s <= cur_end:
                    # 重叠或相邻（若相邻不算重叠则改为 s <= cur_end + 1）
                    cur_end = max(cur_end, e)
                else:
                    merged.append((cur_start, cur_end))
                    cur_start, cur_end = s, e
            merged.append((cur_start, cur_end))
        
            # 统计该染色体的TE长度和数量
            te_total_len = 0
            te_total_count = 0
            for s, e in merged:
                te_total_len += (e - s + 1)   # 按区间长度（包含两端）
                te_total_count += 1
            
            # 获取染色体总长
            chr_total_len = chr_length_dict.get(chrom, 0)
            
            # 添加到结果列表
            result_list.append([one_sample, chrom, chr_total_len, te_total_count, te_total_len])
        
        return result_list
        
    with open('./chr2EDTA/allaverage','w')as f:
        f.write(f"one_sample\tchrom\tchr_total_len\tte_total_count\tte_total_len\n")
        with Pool(processes=50) as pool:

            for i, result in enumerate(pool.imap(run_step1, sample_list), start=1):
                for one in result:
                    
                    one_sample, chrom, chr_total_len, te_total_count, te_total_len=one
                    f.write(f"{one_sample}\t{chrom}\t{chr_total_len}\t{te_total_count}\t{te_total_len}\n")
                progress = (i / len(sample_list)) * 100
                sys.stdout.write(f"\nProgress: {progress:.2f}% {str(i)} / {str(len(sample_list))}#####################################################################\n")
                sys.stdout.flush()              
        
if argv1 == "part_allaverage_nocen":
    sample_list = []
    with open('/home/lain/aaa_data/run0/samples_satellite/sample_info', 'r') as f:
        for line in f:
            eachline_arr = line.strip().split()
            if len(eachline_arr) < 3:
                continue
            sample = eachline_arr[0]
            group = eachline_arr[3]
            if group == 'Routundifolia':
                continue
            sample_list.append(sample)

    ## Obtain centromere region positions (centromere interval for each chromosome of each sample)
    dict_cen_region = {}  # key: f"{sample}|||{chromosome}", value: [start, end]
    with open('/home/lain/aaa_data/run0/samples_satellite/2_good_regions_main_66_103', 'r') as f:
        # sample    chromosome    start    end    length
        next(f)
        for line in f:
            eachline_arr = line.strip().split()
            if len(eachline_arr) < 5:
                continue
            sample, chromosome, start, end, length = eachline_arr
            sample_chr = f"{sample}|||{chromosome}"
            dict_cen_region[sample_chr] = [int(start), int(end)]

    def run_step1(one_sample):
        # Dictionary to store chromosome total lengths
        chr_length_dict = {}
        with open(f'./chr2EDTA/0_prepare/{one_sample}/{one_sample}.fasta.fai', 'r') as f:
            for line in f:
                eachline_arr = line.strip().split('\t')
                chr_name = eachline_arr[0]
                chr_len = int(eachline_arr[1])
                chr_length_dict[chr_name] = chr_len

        # Store TE intervals for each chromosome (after filtering)
        from collections import defaultdict
        chr_intervals = defaultdict(list)

        with open(f"./chr2EDTA/1_EDTA/{one_sample}/{one_sample}.fasta.mod.EDTA.TEanno.gff3", 'r') as f:
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

                # Only record TEs that need to be counted
                chr_intervals[one_chromosome].append((one_start, one_end))

        # Result lists for different region types
        result_list = []  # Whole chromosome
        result_list_arm = []  # Chromosome arms (after removing centromere)
        result_list_cen = []  # Centromere region

        for chrom, intervals in chr_intervals.items():
            # Get chromosome total length
            chr_total_len = chr_length_dict.get(chrom, 0)

            # ========== 1. Count whole chromosome ==========
            if not intervals:
                result_list.append([one_sample, chrom, chr_total_len, 0, 0])
            else:
                # Sort by start position
                intervals.sort()
                merged = []
                cur_start, cur_end = intervals[0]
                for s, e in intervals[1:]:
                    if s <= cur_end:
                        cur_end = max(cur_end, e)
                    else:
                        merged.append((cur_start, cur_end))
                        cur_start, cur_end = s, e
                merged.append((cur_start, cur_end))

                te_total_len = 0
                te_total_count = 0
                for s, e in merged:
                    te_total_len += (e - s + 1)
                    te_total_count += 1

                result_list.append([one_sample, chrom, chr_total_len, te_total_count, te_total_len])

            # ========== 2. Count chromosome arms (after removing centromere) ==========
            sample_chr_key = f"{one_sample}|||{chrom}"
            cen_region = dict_cen_region.get(sample_chr_key, None)

            if cen_region is None:
                # If no centromere information, skip or record as empty
                result_list_arm.append([one_sample, chrom, chr_total_len, 0, 0, "no_cen_info"])
                result_list_cen.append([one_sample, chrom, 0, 0, 0, "no_cen_info"])  # Record centromere region as empty as well
                continue

            cen_start, cen_end = cen_region

            # Chromosome arm regions = whole chromosome - centromere region
            arm_regions = []
            if cen_start > 0:
                arm_regions.append((0, cen_start))  # Left arm (or short arm)
            if cen_end < chr_total_len:
                arm_regions.append((cen_end, chr_total_len))  # Right arm (or long arm)

            # Calculate total length of chromosome arms (excluding centromere)
            arm_total_len = sum(end - start for start, end in arm_regions)

            # Find TE intervals that fall within the chromosome arm regions
            if not intervals:
                result_list_arm.append([one_sample, chrom, arm_total_len, 0, 0])
            else:
                arm_te_intervals = []
                for te_start, te_end in intervals:
                    # Check if this TE overlaps with any arm region
                    for arm_start, arm_end in arm_regions:
                        overlap_start = max(te_start, arm_start)
                        overlap_end = min(te_end, arm_end)
                        if overlap_start < overlap_end:
                            # Only record the overlapping portion (if TE spans the centromere, only count the part on the arm)
                            arm_te_intervals.append((overlap_start, overlap_end))

                # Merge TE intervals on chromosome arms
                if arm_te_intervals:
                    arm_te_intervals.sort()
                    arm_merged = []
                    cur_s, cur_e = arm_te_intervals[0]
                    for s, e in arm_te_intervals[1:]:
                        if s <= cur_e:
                            cur_e = max(cur_e, e)
                        else:
                            arm_merged.append((cur_s, cur_e))
                            cur_s, cur_e = s, e
                    arm_merged.append((cur_s, cur_e))

                    arm_te_total_len = 0
                    arm_te_total_count = 0
                    for s, e in arm_merged:
                        arm_te_total_len += (e - s + 1)
                        arm_te_total_count += 1
                else:
                    arm_te_total_len = 0
                    arm_te_total_count = 0

                result_list_arm.append([one_sample, chrom, arm_total_len, arm_te_total_count, arm_te_total_len])

            # ========== 3. Count TEs in centromere region ==========
            if not intervals:
                result_list_cen.append([one_sample, chrom, (cen_end - cen_start), 0, 0])
            else:
                cen_te_intervals = []
                for te_start, te_end in intervals:
                    # Check if this TE overlaps with the centromere region
                    overlap_start = max(te_start, cen_start)
                    overlap_end = min(te_end, cen_end)
                    if overlap_start < overlap_end:
                        # Only record the overlapping portion
                        cen_te_intervals.append((overlap_start, overlap_end))

                # Merge TE intervals on the centromere region
                if cen_te_intervals:
                    cen_te_intervals.sort()
                    cen_merged = []
                    cur_s, cur_e = cen_te_intervals[0]
                    for s, e in cen_te_intervals[1:]:
                        if s <= cur_e:
                            cur_e = max(cur_e, e)
                        else:
                            cen_merged.append((cur_s, cur_e))
                            cur_s, cur_e = s, e
                    cen_merged.append((cur_s, cur_e))

                    cen_te_total_len = 0
                    cen_te_total_count = 0
                    for s, e in cen_merged:
                        cen_te_total_len += (e - s + 1)
                        cen_te_total_count += 1
                else:
                    cen_te_total_len = 0
                    cen_te_total_count = 0

                result_list_cen.append([one_sample, chrom, (cen_end - cen_start), cen_te_total_count, cen_te_total_len])

        return result_list, result_list_arm, result_list_cen

    # Assign tasks to processes in the process pool
    with open('./chr2EDTA/allaverage', 'w') as f:
        f.write(f"one_sample\tchrom\tchr_total_len\tte_total_count\tte_total_len\tregion_type\n")

    with open('./chr2EDTA/allaverage_arm', 'w') as f_arm:
        f_arm.write(f"one_sample\tchrom\tarm_total_len\tte_total_count\tte_total_len\tregion_type\n")

    with open('./chr2EDTA/allaverage_cen', 'w') as f_cen:
        f_cen.write(f"one_sample\tchrom\tcen_total_len\tte_total_count\tte_total_len\tregion_type\n")

    with Pool(processes=50) as pool:
        # Use imap to get results one by one
        for i, (result, result_arm, result_cen) in enumerate(pool.imap(run_step1, sample_list), start=1):
            # Write whole chromosome results
            for one in result:
                one_sample, chrom, chr_total_len, te_total_count, te_total_len = one
                with open('./chr2EDTA/allaverage', 'a') as f:
                    f.write(f"{one_sample}\t{chrom}\t{chr_total_len}\t{te_total_count}\t{te_total_len}\twhole_chromosome\n")

            # Write results after removing centromere
            for one in result_arm:
                if len(one) == 6 and one[5] == "no_cen_info":
                    # Skip or specially mark those without centromere information
                    continue
                one_sample, chrom, arm_total_len, te_total_count, te_total_len = one[:5]
                with open('./chr2EDTA/allaverage_arm', 'a') as f_arm:
                    f_arm.write(f"{one_sample}\t{chrom}\t{arm_total_len}\t{te_total_count}\t{te_total_len}\tchromosome_arm\n")

            # Write centromere region results
            for one in result_cen:
                if len(one) == 6 and one[5] == "no_cen_info":
                    # Skip or specially mark those without centromere information
                    continue
                one_sample, chrom, cen_total_len, te_total_count, te_total_len = one[:5]
                with open('./chr2EDTA/allaverage_cen', 'a') as f_cen:
                    f_cen.write(f"{one_sample}\t{chrom}\t{cen_total_len}\t{te_total_count}\t{te_total_len}\tcentromere\n")

            progress = (i / len(sample_list)) * 100
            sys.stdout.write(f"\nProgress: {progress:.2f}% {str(i)} / {str(len(sample_list))}#####################################################################\n")
            sys.stdout.flush()

        
        
        
        
        

print("\n\n")        
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))      
time_end=timeit.default_timer()
print('All the running time: %.0f Seconds'%(time_end-time_start))
    
    
