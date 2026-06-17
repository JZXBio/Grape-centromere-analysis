#!/usr/bin/python
# -*- coding: UTF-8 -*-
import sys
argvs=sys.argv
argv_len=len(sys.argv)
if argv_len < 2 or argvs[1] in ("h", "-h", "help", "-help") :
    print ("chr2moddotplot3.py-----help:")
    print ("")
    print ("Usage：")
    print ("chr2moddotplot3.py step0 -i  Specify a genome folder")
    print ("chr2moddotplot3.py step1      Run the main program")
    print ("chr2moddotplot3.py step1.5    Classify into in and out based on relative position to VSat1")    
    print ("chr2moddotplot3.py step2      2a/2b, Perform line identification analysis for in and out separately")
    print ("chr2moddotplot3.py step2_combine   step2a(b)_combine, Merge in and out separately into a single file")
    print ("chr2moddotplot3.py step3a    Find plotting methods for up lines within VSat1")
    print ("chr2moddotplot3.py step3b    For down lines within VSat1, very few, not very useful")
    
    
    print ( "python chr2moddotplot3.py part1/part2 PN40024 Chr1 12000000-20000000      #####part1 static/part2 interactive visualization")
    print ( "python chr2moddotplot3.py vs1/vs2 PN40024 Chr1 12000000-20000000 V024.hap1 Chr1 12000000-20000000  Static/interactive visualization of fixed positions on two chromosomes respectively, requires 7 parameters, pos can be 'all'")
    
    print ( "python chr2moddotplot3.py vs_all V024.hap1 PN40024  Static visualization of all pairwise chromosome comparisons, generating multiple subplots statically")
    
    print ( "python chr2moddotplot3.py select  Manually enter a position, then combine into one sequence, compare with itself")
    
    print ( "python chr2moddotplot3.py fa2fa  sequence file 1 sequence file 2   ### Input two fasta sequences, dynamic alignment")
    
    print ( "python chr2moddotplot3.py Chrall   Chr??       # Output all Chr? core regions")
    
    print ( "python chr2moddotplot3.py Chrall2one   Chr2 PN40024:Chr2:12000000-15000000      ")
    print ( "python chr2moddotplot3.py Chrall2one   Chr2 PN40024:Chr17:14000000-16000000     ")
    print ( "python chr2moddotplot3.py Chrall2one   Chr2 V033.hap2:Chr2:11500000 13500000     ")
    print ( "python chr2moddotplot3.py Chrall2one   Chr3  ShineMuscat_hap1:Chr3:13000000-15000000       # Output comparison chart of all Chr? core regions and specified regions  Chr3(or all) ")
    print ( "python chr2moddotplot3.py chrall2one   Chr4 PN40024:Chr4:12000000-14000000      ")
    print ( "python chr2moddotplot3.py Chrall2one   Chr7 PN40024:Chr7:12000000-15000000      ")
    print ( "python chr2moddotplot3.py Chrall2one   Chr10  V043.hap1:Chr10:19800000-20200000    ")
    print ( "python chr2moddotplot3.py Chrall2one   Chr13  PN40024:Chr13:11800000-12200000 ")
    print ( "python chr2moddotplot3.py Chrall2one   Chr13  V043.hap2:Chr13:13000000-14500000 ") 
    print ( "python chr2moddotplot3.py Chrall2one   Chr13  PN40024:Chr9:15000000-15500000 ") 
    print ( "python chr2moddotplot3.py Chrall2one   Chr17  PN40024:Chr17:14000000-16000000    ")
    
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
import re # Handle regular expressions
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
#print(args_dict)
##
time_start=timeit.default_timer()
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  

if "thread" in args_dict:  thread=int(args_dict["thread"])  
else:thread=50   


moddotplot="moddotplot"


    
    
if  os.path.exists('./chr2moddotplot3')==False:
    subprocess.run(["mkdir chr2moddotplot3"], shell=True)

#step0
#step0
if argv1=="stepall" or argv1=="step0":
    if  os.path.exists('./chr2moddotplot3/0_prepare')==True:
        subprocess.run(["rm -r ./chr2moddotplot3/0_prepare"], shell=True)  
    subprocess.run(["mkdir ./chr2moddotplot3/0_prepare"], shell=True)  
    
    if  os.path.exists('./chr2moddotplot3/error')==True:subprocess.run(["rm ./chr2moddotplot3/error"], shell=True)  

    if "i"  in args_dict: #print("Missing input fasta file");sys.exit()
        input_dictionary = args_dict["i"]
        with open('./chr2moddotplot3/0_prepare/sample_source','w') as f:
            f.write(input_dictionary)
    else:
        print("Missing input fasta file");sys.exit()
    
    dir_fasta_file_num=0
    dir_file_name_list=[]
    files=os.listdir(input_dictionary)
    with open('./chr2moddotplot3/0_prepare/sample_list','w') as f:
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
        output_dir="./chr2moddotplot3/0_prepare/"+sample_name+"/"
        if  os.path.exists(output_dir)==True:return False
        output_file=output_dir+sample_name
        subprocess.run(["mkdir "+output_dir], shell=True)  
        ### Limit length to 1000000
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
            with open("./chr2moddotplot3/error",'a') as f:
                f.write(f"Chromosome count for {sample_name} is not 19 or 20, it is {str(chr_num)}\n")
                ###
            
    
    # Assign tasks to processes in the process pool    
    with Pool(processes=thread) as pool:
        pool.map(run_step0, dir_file_name_list)      
    
    dir_file_name_list.sort()
    with open('./chr2moddotplot3/chrname_old2new','w') as f3:
        f3.write(f"sample_name\tone_id_old\tone_id\n")    
    for infos in dir_file_name_list:
        input_dir,file_name,sample_name=infos
        
        subprocess.run([f"cat ./chr2moddotplot3/0_prepare/{sample_name}/{sample_name}.chrname_old2new >> ./chr2moddotplot3/chrname_old2new"], shell=True) 
             

    
     
if argv1=="stepall" or argv1=="step1":    
    print("Run main program moddotplot")
    if  os.path.exists('./chr2moddotplot3/1_moddotplot')==False:
        subprocess.run(["mkdir ./chr2moddotplot3/1_moddotplot"], shell=True)  
    
    with open('./chr2moddotplot3/0_prepare/sample_list','r') as f:
        sample_list=f.read().split('\n')

    def run_step0(sample_name):
        #if sample_name!='V072.hap1': return False   #V108.hap1
        if len(sample_name)==0: return False
        input_file=f"./chr2moddotplot3/0_prepare/{sample_name}/{sample_name}.fasta"
        output_dir="./chr2moddotplot3/1_moddotplot/"+sample_name+"/"
        if  os.path.exists(output_dir)==True:return False
        print(sample_name)
        subprocess.run(["mkdir "+output_dir], shell=True) 
        ##ulimit -u  8192
        moddotplot_software='ulimit -u 100000 && source /home/lain/software/venv/ModDotPlot/venv/bin/activate &&  /home/lain/software/venv/ModDotPlot/venv/bin/moddotplot '
        subprocess.run([f"{moddotplot_software} static  --fasta {input_file} --kmer 21 --output-dir {output_dir} --identity 80 --dpi 3000 --resolution 5000 "], shell=True)  
    
    # Assign tasks to processes in the process pool    
    with Pool(processes=10) as pool:
         # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step0, sample_list), start=1):
            # Results can be processed here, e.g., stored or printed
            progress = (i / len(sample_list)) * 100
            print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))      
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush() 
if argv1=="stepall" or argv1=="step1.5":
    print('Extract points inside and outside VSat1')
    subprocess.run(["mkdir ./chr2moddotplot3/1_VSat12split"], shell=True)  
    
    print('Loading input file list')
    with open('./chr2moddotplot3/0_prepare/sample_list','r') as f:
        sample_list=f.read().split('\n')
    info_list=[]
    with open('./chr2moddotplot3/0_prepare/sample_list','r') as f:sample_list=f.read().split('\n')
    for sample_name in sample_list:
        input_dir=f"./chr2moddotplot3/1_moddotplot/{sample_name}"
        output_dir=f"./chr2moddotplot3/1_VSat12split/{sample_name}" 
        info_list.append([sample_name,input_dir,output_dir])   
        
    def run_step2(one_info):
        sample_name,input_dir,output_dir=one_info
        subprocess.run([f"mkdir {output_dir}"], shell=True)  
        print(sample_name)
        if len(sample_name)==0: return False    
        #### Remove points related to VSat1
        files=os.listdir(input_dir)
        chr_name_list=[]
        for x in files:
            if x[-4:]=='.bed':
                chr_name=x[:-4]
                chr_name_list.append(chr_name)
        chr_name_list.sort()  
        ####
        
        ####
        print('Processing bed')
        for one_chr in chr_name_list:   
            print('Generating dictionary')
            dict_VSat1={}
            with open ("/home/lain/aaa_data/run0/samples_satellite/2_good_regions",'r') as f:
                next(f)
                for line in f:
                    eachline_arr =line.strip().split('\t')
                    if len(eachline_arr)!=8:continue
                    #sample	region_name	region_pos	bigblock_chrstart	bigblock_chrend	chr_region_length	strand	match_percent
                    sample= eachline_arr[0]
                    chromosome=eachline_arr[1].replace("region_",'Chr') 
                    if sample!=sample_name: continue
                    if chromosome!=one_chr: continue
                    
                    chrstart=       int(eachline_arr[3])
                    chrend=         int(eachline_arr[4])
                    kkk=chrstart
                    while kkk<=chrend:
                        dict_VSat1[kkk]=1
                        kkk+=1
            with open(f"{output_dir}/{one_chr}.bed.VSat1_in",'w') as f2:
                f2.write("#query_name\tquery_start\tquery_end\treference_name\treference_start\treference_end\tperID_by_events\n")
                with open(f"{output_dir}/{one_chr}.bed.VSat1_out",'w') as f3:
                    f3.write("#query_name\tquery_start\tquery_end\treference_name\treference_start\treference_end\tperID_by_events\n")
                    with open(f"{output_dir}/{one_chr}.bed.VSat1_cross",'w') as f4:
                        f4.write("#query_name\tquery_start\tquery_end\treference_name\treference_start\treference_end\tperID_by_events\n")                    
                        with open(f"{input_dir}/{one_chr}.bed",'r') as f:
                            next(f)
                            for line in f:
                                eachline=line.strip()
                                eachline_arr=eachline.split('\t')
                                query_name,query_start,query_end,reference_name,reference_start,reference_end,perID_by_events=eachline_arr
                                if int(query_start) not in dict_VSat1 and \
                                     int(query_end) not in dict_VSat1 and \
                                     int(reference_start) not in dict_VSat1 and \
                                     int(reference_end) not in dict_VSat1 :
                                    f3.write(eachline+'\n')
                                elif  int(query_start) in dict_VSat1 and \
                                     int(query_end) in dict_VSat1 and \
                                     int(reference_start)  in dict_VSat1 and \
                                     int(reference_end)  in dict_VSat1 :
                                    f2.write(eachline+'\n')
                                else:
                                    f4.write(eachline+'\n')

    with Pool(processes=thread) as pool:
         # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step2, info_list), start=1):
            # Results can be processed here, e.g., stored or printed
            progress = (i / len(info_list)) * 100
            sys.stdout.write(f"\nProgress: {progress:.2f}% {str(i)} / {str(len(sample_list))}#####################################################################\n")
            sys.stdout.flush()   
                
if argv1=="stepall" or argv1=="step2" or argv1=="step2a": 
    print('70s')
    import numpy as np
    #from scipy import stats
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    #from scipy.signal import find_peaks
    #from scipy.interpolate import interp1d
    #import matplotlib.pyplot as plt  
    
    print("Analyzing moddotplot bed files, extracting reverse sequences")
    if  os.path.exists('./chr2moddotplot3/2_slashin')==False:
        subprocess.run(["mkdir ./chr2moddotplot3/2_slashin"], shell=True)     

    ### Set parameters
    parafile='./chr2moddotplot3/2_slashin/parameter'
    match_identity=80   # Preliminary filtering of some points seems unnecessary since later we take the largest 10 anyway
    x_ydot_max_num=10        # Filter the largest 10 values for each x position
    x_ydot_max_add_bynum=0.01    # If >100, filter the largest 11 values, etc.
    x_ydot_max_add_byvalue=0.3   # identity, to prevent repetitive sequences
    x_ydot_region_max=5000 # Prevent overly distant matches, max is 5000
    #
    x_ydot_neighbor=5    # At an x position, allowed gap tolerance is 5, meaning the minimum distance between two y points on the same x will not be less than 5
    x_gap_max=5         # x_ydot_neighbor seems related; after x moves right by +1, y can only have one point within ±5,      # Horizontal search
    y_gap_max=15                                                                                    # Vertical search
    extend_tolerance=30
    #
    with open('./chr2moddotplot3/2_slashin/parameter','w') as f:
        f.write(f"match_identity:{match_identity}\n")   
        f.write(f"x_ydot_max_num:{x_ydot_max_num}\n")
        f.write(f"x_ydot_max_add_bynum:{x_ydot_max_add_bynum}\n")
        f.write(f"x_ydot_max_add_byvalue:{x_ydot_max_add_byvalue}\n")
        f.write(f"x_ydot_region_max:{x_ydot_region_max}\n")  
        f.write(f"x_ydot_neighbor:{x_ydot_neighbor}\n")  
        f.write(f"x_gap_max:{x_gap_max}\n")     
        f.write(f"y_gap_max:{y_gap_max}\n")    
        f.write(f"extend_tolerance:{extend_tolerance}\n")    
    
    info_list=[]
    with open('./chr2moddotplot3/0_prepare/sample_list','r') as f:sample_list=f.read().split('\n')
    for sample_name in sample_list:
        input_dir=f"./chr2moddotplot3/1_VSat12split/{sample_name}"
        output_dir=f"./chr2moddotplot3/2_slashin/{sample_name}"        
        info_list.append([sample_name,input_dir,output_dir,parafile])        
    

        
    # Assign tasks to processes in the process pool    
    def run_step2(one_info):
        sample_name,input_dir,output_dir,parafile=one_info
        if len(sample_name)==0: return False
        #if sample_name!='PN40024_hap1': return False
        print(sample_name)
        
        subprocess.run([f"mkdir {output_dir}"], shell=True)     
        files=os.listdir(input_dir)
        
        ## Load parameters
        with open(parafile,'r') as f:
            parameters=f.read().split('\n')
            for one_parameter in parameters:
                if "match_identity:" in one_parameter: match_identity=float(one_parameter.split('match_identity:')[1])
                if "x_ydot_max_num:" in one_parameter: x_ydot_max_num=int(one_parameter.split('x_ydot_max_num:')[1])
                if "x_ydot_max_add_bynum:" in one_parameter: x_ydot_max_add_bynum=float(one_parameter.split('x_ydot_max_add_bynum:')[1])
                if "x_ydot_max_add_byvalue:" in one_parameter: x_ydot_max_add_byvalue=float(one_parameter.split('x_ydot_max_add_byvalue:')[1])
                if "x_ydot_neighbor:" in one_parameter: x_ydot_neighbor=int(one_parameter.split('x_ydot_neighbor:')[1])
                if "x_ydot_region_max:" in one_parameter: x_ydot_region_max=int(one_parameter.split('x_ydot_region_max:')[1])
                if "x_gap_max:" in one_parameter: x_gap_max=int(one_parameter.split('x_gap_max:')[1])
                if "y_gap_max:" in one_parameter: y_gap_max=int(one_parameter.split('y_gap_max:')[1])
                if "extend_tolerance:" in one_parameter: extend_tolerance=int(one_parameter.split('extend_tolerance:')[1])
                
        #print(f"match_identity:{str(match_identity)}")        
        #print(f"x_ydot_max_num:{str(x_ydot_max_num)}")
        #print(f"x_ydot_max_add_bynum:{str(x_ydot_max_add_bynum)}")
        #print(f"x_ydot_max_add_byvalue:{str(x_ydot_max_add_byvalue)}")
        #print(f"x_ydot_neighbor:{str(x_ydot_neighbor)}")
        #print(f"x_ydot_region_max:{str(x_ydot_region_max)}")
        #print(f"x_gap_max:{str(x_gap_max)}")
        #print(f"y_gap_max:{str(y_gap_max)}")
        #print(f"extend_tolerance:{str(extend_tolerance)}")
        
        
        #########################
        chr_name_list=[]
        for x in files:
            if x[-13:]=='.bed.VSat1_in':
                chr_name=x[:-13]
                chr_name_list.append(chr_name)
        chr_name_list.sort()        
        #print(chr_name_list)
        ########################
        
        for one_chr in chr_name_list:   
            #print(one_chr)
            ### Load bed file into dictionary dict_x_y_z
            dict_x_y_z={}
            x_list=[]
            with open(f'{output_dir}/{one_chr}.1_newbed','w') as ff:
                #ff.write(f"pos1\tpos2\tperID_by_events\n")   
                unit_len=''
                with open(f"{input_dir}/{one_chr}.bed.VSat1_in") as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        if eachline_arr[0]=='#query_name':continue
                        if unit_len=='':  
                            unit_len=int(eachline_arr[2])-int(eachline_arr[1])+1
                            ff.write(f"#unit_len={unit_len}\n")
  
                        ##query_name	query_start	query_end	reference_name	reference_start	reference_end	perID_by_events
                        ##BMNGT2T_hap1_chr02	1	21454	BMNGT2T_hap1_chr02	1	21454	100
                        new_pos1 = int((int(eachline_arr[1])-1)/unit_len);          new_pos1_str=str(new_pos1);              x_list.append(new_pos1)
                        new_pos2 = int((int(eachline_arr[4])-1)/unit_len);          new_pos2_str=str(new_pos2);
                        percent_str  = eachline_arr[6];                             percent_float = float(percent_str)
                        #if new_pos2-new_pos1>500:continue    #20000*200 roughly 4M
                        ff.write(f"{new_pos1_str}\t{new_pos2_str}\t{percent_str}\n")
                        
                        if new_pos2-new_pos1>x_ydot_region_max:continue
                        if percent_float<match_identity:continue 
                        if new_pos1_str not in dict_x_y_z: dict_x_y_z[new_pos1_str]={}
                        dict_x_y_z[new_pos1_str][new_pos2_str]=percent_float
            
            ### Filter out some points, keep the top 10 (variable) points for each x position
            dict2_x_y_z={}
            with open(f'{output_dir}/{one_chr}.2_simplebed','w') as fff: 
                fff.write(f"#unit_len={unit_len}\n")
                for new_pos1_str,new_pos2_str_dict in dict_x_y_z.items():
                    # Sort dictionary values in descending order
                    x2y_valuelist = sorted(new_pos2_str_dict.values(), key=lambda item: float(item), reverse=True)
                    
                    # Fixed number + (total count / 100)
                    keep_num=x_ydot_max_num+round(len(x2y_valuelist) * float(x_ydot_max_add_bynum))
                    # Find the minimum value, subtract x_ydot_max_add_byvalue
                    x2y_valuelist_min=min(x2y_valuelist[:keep_num-1])-x_ydot_max_add_byvalue

                    # Create a new dictionary with values greater than x2y_valuelist_min
                    filtered_dict = {key: value for key, value in new_pos2_str_dict.items() if float(value) > x2y_valuelist_min}
                    
                    for new_pos2_str,z in filtered_dict.items():
                        fff.write(f"{new_pos1_str}\t{new_pos2_str}\t{z}\n")
                    #
                    dict2_x_y_z[new_pos1_str]=filtered_dict
                    
            ### Filter out some points more. For the remaining points, merge by position. Points with neighbor distance <5 (variable) are considered the same peak, keep the max value
            with open(f'{output_dir}/{one_chr}.3_simplebed','w') as ffff: 
                ffff.write(f"#unit_len={unit_len}\n")
                for new_pos1_str,new_pos2_str_dict in dict2_x_y_z.items():
                    # Get keys as list sorted by value in descending order
                    sorted_keys = sorted(new_pos2_str_dict.keys(),key=lambda key: int(key), reverse=True  )
                    sorted_keys_list = list(sorted_keys)
                    #print(sorted_keys_list);
                    goodkey='';key_edge='';goodvalue='';
                    for new_pos2_str in sorted_keys_list:
                        value=float(new_pos2_str_dict[new_pos2_str])
                        #print(value)
                        if goodkey=='': goodkey=new_pos2_str;     key_edge=new_pos2_str;      goodvalue=value
                        #print(new_pos2_str,int(key_edge)+x_ydot_neighbor)
                        if int(new_pos2_str)>=int(key_edge)-x_ydot_neighbor:   
                            key_edge=new_pos2_str
                            if goodvalue<=value:     goodkey=new_pos2_str;    goodvalue=value;           ## Why have =? Mainly because near the diagonal there may be two 100s, ensuring 100 is the lowest one, i.e., the diagonal
                        else:    
                            if new_pos1_str!=goodkey:
                                ffff.write(f"{new_pos1_str}\t{goodkey}\t{str(goodvalue)}\n")
                            #print(value)
                            goodkey=new_pos2_str;key_edge=new_pos2_str;goodvalue=value;
                    if goodkey!='':   
                        if new_pos1_str!=goodkey:
                            ffff.write(f"{new_pos1_str}\t{goodkey}\t{str(goodvalue)}\n")
                    goodkey='';key_edge='';goodvalue='';
                    
            subprocess.run([f"rm {output_dir}/{one_chr}.1_newbed"], shell=True)              
            subprocess.run([f"rm {output_dir}/{one_chr}.2_simplebed"], shell=True) 
            
            
            # Point pool 
            dict3_x_y_z={}  ## This is a numeric dictionary xy
            with open(f'{output_dir}/{one_chr}.3_simplebed','r') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if len(eachline_arr)<3:continue
                    pos1,pos2,value=eachline_arr
                    pos1=int(pos1)
                    pos2=int(pos2)
                    if pos1 not in dict3_x_y_z: dict3_x_y_z[pos1]={}
                    dict3_x_y_z[pos1][pos2]=value
            x_list=sorted(dict3_x_y_z.keys())       

            ### Calculate diagonally upward
            with open(f'{output_dir}/{one_chr}.4_line2up','w') as fffff:
                fffff.write(f"#unit_len={unit_len}\n")
                fffff.write(f"#start_x\tstart_y\tend_x\tend_y\tseq1_len\tseq2_len\tdot_num\tdotskip_num\tdot_used\tline_slope\tline_R2\tback_accuracy1\tback_accuracy2\tback_distance1\tback_distance2\n")  #\tdots_string
                tmp_dot_list=[]
                
                for x in x_list:
                    y_list=sorted(dict3_x_y_z[x].keys(), reverse=True)     ####y_list=list(dict3_x_y_z[x].keys);  y_list.sort()   
                    #print(y_list)
                    for y in y_list:
                        # Current point is (x,y)
                        x_start=x
                        y_start=y
                        end_mark=''
                        tolerance=0
                        #print(x,y,type(dict3_x_y_z))
                        serial=0;
                        while end_mark=='' and tolerance<extend_tolerance:
                            ## Check if the current start point has been used in previous loops
                            if x_start not in dict3_x_y_z:              end_mark='yes';break
                            if y_start not in dict3_x_y_z[x_start]:     end_mark='yes';break
                            xy_distance=y_start-x_start
                            ##
                            if tmp_dot_list==[]:tmp_dot_list.append([x,y]);serial+=1
                            # Collect the next point
                            i1=0;good_mark=''
                            while i1<x_gap_max:
                                i1+=1
                                
                                if good_mark=='yes':break
                                if serial<3:i2=1   # Slope of the first three points cannot be horizontal
                                elif xy_distance<20:    i2=1        # The smallest duplication might be 100k, i.e., 20 points
                                elif xy_distance<40:    i2=0
                                elif xy_distance<90:    i2=-1
                                else :                  i2=-2
                                while i2<y_gap_max:
                                    x_tmp=x_start+i1
                                    y_tmp=y_start+i2  ### Move x one step right, y first stays (+0) then moves upward (++1)
                                    i2+=1
                                    if x_tmp in dict3_x_y_z:
                                        if y_tmp in dict3_x_y_z[x_tmp]:
                                            # New coordinates are (x_tmp,y_tmp)
                                            good_mark='yes';break
                            ## If the next point is collected, continue with it as the starting point            
                            if good_mark=='yes': 
                                tmp_dot_list.append([x_tmp,y_tmp])
                                x_start=x_tmp
                                y_start=y_tmp
                                serial+=1
                            ## If no next point is collected, end the current line segment    
                            else:
                                tolerance+=1
                                x_start+=1
                                y_start+=1  ##
                                
                      
                        ### Determine the number of points in the currently collected point sequence, slope, and fit. [[x1,y1],[x2,y2],...], calculate the slope and R² of the fitted line for all points
                        tmp_dot_list_len=len(tmp_dot_list)
                        if tmp_dot_list_len>5:
                            points = np.array(tmp_dot_list)
                            np_x = points[:, 0]    # Extract x and y values
                            np_y = points[:, 1]
                            slope, intercept = np.polyfit(np_x, np_y, 1)  # Use numpy's polyfit function for linear regression, 1 indicates first-order linear fit
                            #print(f"Fitted line slope: {slope}")
                            #print(f"Fitted line intercept: {intercept}")
                            y_pred = slope * np_x + intercept  # Calculate R² goodness of fit
                            ss_total = np.sum((np_y - np.mean(np_y))**2)
                            ss_res = np.sum((np_y - y_pred)**2)
                            if ss_total == 0:       r_squared = 0  # or other default value
                            else:                   r_squared = 1 - (ss_res / ss_total)
                            #print(f"R² goodness of fit: {r_squared}")
                            if r_squared>0.5 and abs(slope)>0.7 and abs(slope)<1.42 and tmp_dot_list_len>5:
                                all_start_x=tmp_dot_list[0][0]*unit_len+1   ;               seq1_start=all_start_x     #seq1_start
                                all_start_y=tmp_dot_list[0][1]*unit_len+1+(unit_len-1)   ;  seq2_start=all_start_y     #seq2_start     #+(unit_len-1) because the simplified start position was used
                                all_end_x=tmp_dot_list[-1][0]*unit_len+1    ;               seq1_end=all_end_x      #seq1_end
                                all_end_y=tmp_dot_list[-1][1]*unit_len+1+(unit_len-1)    ;  seq2_end=all_end_y      #seq2_end
                                seq1_len=abs(seq1_end-seq1_start)+1
                                seq2_len=abs(seq2_end-seq2_start)+1
                                dots_string = '|'.join(f'{x},{y}' for x, y in tmp_dot_list)
                                #if dots_string=="2777,2815|2781,2817|2785,2817|2786,2817|2787,2817|2788,2817|2791,2820|2795,2827|2800,2831|2801,2839|2802,2839|2803,2839|2804,2848|2805,2849|2806,2849|2810,2849|2811,2849|2815,2849|2816,2849|2818,2849|2820,2849|2821,2849|2822,2849":print(points)
                                
                                all_x_num=tmp_dot_list[-1][0]-tmp_dot_list[0][0]+1
                                good_ratio=tmp_dot_list_len/all_x_num
                                
                                delete_dot_mark='yes'
                                if tmp_dot_list_len<3000 and good_ratio>0.7:
                                    #### Calculate the proportion of the line and the surrounding circle, -1~+1 width 3, -3~-5 and +3~+5 width 3,
                                    dict_tmp_dot_list = {int(item[0]): int(item[1]) for item in tmp_dot_list}
                                    line_start_serial=tmp_dot_list[0][0]
                                    line_end_serial=tmp_dot_list[-1][0]
                                    j=line_start_serial
                                    middle_list=[];surround_up_list=[];surround_down_list=[]
                                    middle_valuelist=[];surround_up_valuelist=[];surround_down_valuelist=[]
                                    #print(line_start_serial,line_end_serial)
                                    while j<line_end_serial:
                                        if str(j) in dict_x_y_z:        dict_yz=dict_x_y_z[str(j)];##print(dict_yz)
                                        
                                        if j in dict_tmp_dot_list:
                                            j_y=dict_tmp_dot_list[j]
                                        else:
                                            j_y+=1
                                        ###    
                                        if str(j_y) in dict_yz:     value=float(dict_yz[str(j_y)]);       middle_list.append([j,value]);         middle_valuelist.append(value)
                                        if str(j_y-1) in dict_yz:   value=float(dict_yz[str(j_y-1)]);     middle_list.append([j,value]);         middle_valuelist.append(value)
                                        if str(j_y+1) in dict_yz:   value=float(dict_yz[str(j_y+1)]);     middle_list.append([j,value]);         middle_valuelist.append(value)
                                        if str(j_y+4) in dict_yz:   value=float(dict_yz[str(j_y+4)]);     surround_up_list.append([j,value]);         surround_up_valuelist.append(value)
                                        if str(j_y+5) in dict_yz:   value=float(dict_yz[str(j_y+5)]);     surround_up_list.append([j,value]);         surround_up_valuelist.append(value)
                                        if str(j_y+6) in dict_yz:   value=float(dict_yz[str(j_y+6)]);     surround_up_list.append([j,value]);         surround_up_valuelist.append(value)
                                        if str(j_y-4) in dict_yz:   value=float(dict_yz[str(j_y-4)]);     surround_down_list.append([j,value]);         surround_down_valuelist.append(value)
                                        if str(j_y-5) in dict_yz:   value=float(dict_yz[str(j_y-5)]);     surround_down_list.append([j,value]);         surround_down_valuelist.append(value)
                                        if str(j_y-6) in dict_yz:   value=float(dict_yz[str(j_y-6)]);     surround_down_list.append([j,value]);         surround_down_valuelist.append(value)  
                                        ###
                                        j+=1
                                    # UP    
                                    if len(middle_list)<5 or len(surround_up_list)<5:
                                        accuracy1="NA"
                                    else:    
                                        line_A = np.array(middle_list)
                                        line_B = np.array(surround_up_list)
                                        X = np.vstack((line_A, line_B))      # Merge data and create labels
                                        y = np.hstack((np.zeros(len(line_A)), np.ones(len(line_B))))
                                        lda = LinearDiscriminantAnalysis()# Train LDA model
                                        lda.fit(X, y)
                                        predictions = lda.predict(X) # Check model accuracy
                                        accuracy1 = np.mean(predictions == y)
                                        #print("Coefficients:", lda.coef_)  # Coefficients of the linear boundary
                                        #print("Intercept:", lda.intercept_)  # Intercept of the linear boundary
                                        #print("Accuracy:", accuracy)  # Classification accuracy
                                    # DOWN
                                    if len(middle_list)<5 or len(surround_down_list)<5:
                                        accuracy2="NA"
                                    else:                                            
                                        line_A = np.array(middle_list)
                                        line_B = np.array(surround_down_list)
                                        X = np.vstack((line_A, line_B))      # Merge data and create labels
                                        y = np.hstack((np.zeros(len(line_A)), np.ones(len(line_B))))
                                        lda = LinearDiscriminantAnalysis()# Train LDA model
                                        lda.fit(X, y)
                                        predictions = lda.predict(X) # Check model accuracy
                                        accuracy2 = np.mean(predictions == y)
                                    ####
                                    middle_mean = sum(middle_valuelist) / len(middle_valuelist)
                                    if len(surround_up_valuelist)==0:        up_distance_mean="NA"
                                    else:                               up_mean=sum(surround_up_valuelist) / len(surround_up_valuelist);          up_distance_mean= round(middle_mean-up_mean,2)
                                    if len(surround_down_valuelist)==0:      down_distance_mean="NA"
                                    else:                               down_mean=sum(surround_down_valuelist) / len(surround_down_valuelist);          down_distance_mean= round(middle_mean-down_mean,2)
                                    ####
                                    if accuracy1!='NA' and accuracy1<0.6:                          delete_dot_mark='no'
                                    elif accuracy2!='NA' and accuracy2<0.6:                          delete_dot_mark='no'
                                    elif up_distance_mean!='NA' and up_distance_mean<0.2:       delete_dot_mark='no'
                                    elif down_distance_mean!='NA' and down_distance_mean<0.2:   delete_dot_mark='no' 
                                    else:    
                                        fffff.write(f"{str(all_start_x)}\t{str(all_start_y)}\t{str(all_end_x)}\t{str(all_end_y)}\t{str(seq1_len)}\t{str(seq2_len)}\t{str(tmp_dot_list_len)}\t{str(all_x_num-tmp_dot_list_len)}\t{str(round(good_ratio,3))}\t{str(round(slope,2))}\t{str(round(r_squared,2))}\t{accuracy1}\t{accuracy2}\t{str(up_distance_mean)}\t{str(down_distance_mean)}\n")   #\t{dots_string}
           
                                        
                                ## Remove these points from the dictionary
                                if delete_dot_mark=='yes':
                                    for one_result in tmp_dot_list:
                                        one_x,one_y=one_result
                                        del dict3_x_y_z[one_x][one_y]
                        ## Initialize
                        tmp_dot_list=[]
                        #return False                
                        
                ### Calculate diagonally downward
            
            ### Calculate diagonally downward
            with open(f'{output_dir}/{one_chr}.4_line2down','w') as fffff:
                fffff.write(f"#unit_len={unit_len}\n")
                fffff.write(f"#start_x\tstart_y\tend_x\tend_y\tseq1_len\tseq2_len\tdot_num\tdotskip_num\tdot_used\tline_slope\tline_R2\tback_accuracy1\tback_accuracy2\tback_distance1\tback_distance2\n") #\tdots_string
                tmp_dot_list=[]
                
                for x in x_list:
                    y_list=sorted(dict3_x_y_z[x].keys(), reverse=True)     ####y_list=list(dict3_x_y_z[x].keys);  y_list.sort()   
                    #print(y_list)
                    for y in y_list:
                        # Current point is (x,y)
                        x_start=x
                        y_start=y
                        end_mark=''
                        tolerance=0
                        #print(x,y,type(dict3_x_y_z))
                        serial=0;
                        while end_mark=='' and tolerance<extend_tolerance:
                            ## Check if the current start point has been used in previous loops
                            if x_start not in dict3_x_y_z:              end_mark='yes';break
                            if y_start not in dict3_x_y_z[x_start]:     end_mark='yes';break
                            xy_distance=y_start-x_start
                            ##
                            if tmp_dot_list==[]:tmp_dot_list.append([x,y]);serial+=1
                            # Collect the next point
                            i1=0;good_mark=''
                            while i1<x_gap_max:
                                i1+=1
                                
                                if good_mark=='yes':break
                                if serial<3:i2=1   # Slope of the first three points cannot be horizontal
                                elif xy_distance<20:    i2=1        # The smallest duplication might be 100k, i.e., 20 points
                                elif xy_distance<40:    i2=0
                                elif xy_distance<90:    i2=-1
                                else :                  i2=-2
                                while i2<y_gap_max:
                                    x_tmp=x_start+i1
                                    y_tmp=y_start-i2  ### Move x one step right, y first stays (+0) then moves upward (++1)             !!!!!!!!!!!!!!!!!!!!!!!!! Change direction here
                                    i2+=1
                                    if x_tmp in dict3_x_y_z:
                                        if y_tmp in dict3_x_y_z[x_tmp]:
                                            # New coordinates are (x_tmp,y_tmp)
                                            good_mark='yes';break
                            ## If the next point is collected, continue with it as the starting point            
                            if good_mark=='yes': 
                                tmp_dot_list.append([x_tmp,y_tmp])
                                x_start=x_tmp
                                y_start=y_tmp
                                serial+=1
                            ## If no next point is collected, end the current line segment    
                            else:
                                tolerance+=1
                                x_start+=1
                                y_start-=1  ##                                                           !!!!!!!!!!!!!!!!!!!!!!!!! Change direction here
                                
                      
                        ### Determine the number of points in the currently collected point sequence, slope, and fit. [[x1,y1],[x2,y2],...], calculate the slope and R² of the fitted line for all points
                        tmp_dot_list_len=len(tmp_dot_list)
                        if tmp_dot_list_len>5:
                            points = np.array(tmp_dot_list)
                            np_x = points[:, 0]    # Extract x and y values
                            np_y = points[:, 1]
                            slope, intercept = np.polyfit(np_x, np_y, 1)  # Use numpy's polyfit function for linear regression, 1 indicates first-order linear fit
                            #print(f"Fitted line slope: {slope}")
                            #print(f"Fitted line intercept: {intercept}")
                            y_pred = slope * np_x + intercept  # Calculate R² goodness of fit
                            ss_total = np.sum((np_y - np.mean(np_y))**2)
                            ss_res = np.sum((np_y - y_pred)**2)
                            if ss_total == 0:       r_squared = 0  # or other default value
                            else:                   r_squared = 1 - (ss_res / ss_total)
                            #print(f"R² goodness of fit: {r_squared}")
                            if r_squared>0.5 and abs(slope)>0.7 and abs(slope)<1.42 and tmp_dot_list_len>5:
                                all_start_x=tmp_dot_list[0][0]*unit_len+1   ;               seq1_start=all_start_x     #seq1_start
                                all_start_y=tmp_dot_list[0][1]*unit_len+1+(unit_len-1)   ;  seq2_start=all_start_y     #seq2_start     #+(unit_len-1) because the simplified start position was used
                                all_end_x=tmp_dot_list[-1][0]*unit_len+1    ;               seq1_end=all_end_x      #seq1_end
                                all_end_y=tmp_dot_list[-1][1]*unit_len+1+(unit_len-1)    ;  seq2_end=all_end_y      #seq2_end
                                seq1_len=abs(seq1_end-seq1_start)+1
                                seq2_len=abs(seq2_end-seq2_start)+1
                                dots_string = '|'.join(f'{x},{y}' for x, y in tmp_dot_list)
                                #if dots_string=="2777,2815|2781,2817|2785,2817|2786,2817|2787,2817|2788,2817|2791,2820|2795,2827|2800,2831|2801,2839|2802,2839|2803,2839|2804,2848|2805,2849|2806,2849|2810,2849|2811,2849|2815,2849|2816,2849|2818,2849|2820,2849|2821,2849|2822,2849":print(points)
                                
                                all_x_num=tmp_dot_list[-1][0]-tmp_dot_list[0][0]+1
                                good_ratio=tmp_dot_list_len/all_x_num
                                
                                delete_dot_mark='yes'
                                if tmp_dot_list_len<3000 and good_ratio>0.7:
                                    #### Calculate the proportion of the line and the surrounding circle, -1~+1 width 3, -3~-5 and +3~+5 width 3,
                                    dict_tmp_dot_list = {int(item[0]): int(item[1]) for item in tmp_dot_list}
                                    line_start_serial=tmp_dot_list[0][0]
                                    line_end_serial=tmp_dot_list[-1][0]
                                    j=line_start_serial
                                    middle_list=[];surround_up_list=[];surround_down_list=[]
                                    middle_valuelist=[];surround_up_valuelist=[];surround_down_valuelist=[]
                                    #print(line_start_serial,line_end_serial)
                                    while j<line_end_serial:
                                        if str(j) in dict_x_y_z:        dict_yz=dict_x_y_z[str(j)];##print(dict_yz)
                                        
                                        if j in dict_tmp_dot_list:
                                            j_y=dict_tmp_dot_list[j]
                                        else:
                                            j_y+=1
                                        ###    
                                        if str(j_y) in dict_yz:     value=float(dict_yz[str(j_y)]);       middle_list.append([j,value]);         middle_valuelist.append(value)
                                        if str(j_y-1) in dict_yz:   value=float(dict_yz[str(j_y-1)]);     middle_list.append([j,value]);         middle_valuelist.append(value)
                                        if str(j_y+1) in dict_yz:   value=float(dict_yz[str(j_y+1)]);     middle_list.append([j,value]);         middle_valuelist.append(value)
                                        if str(j_y+4) in dict_yz:   value=float(dict_yz[str(j_y+4)]);     surround_up_list.append([j,value]);         surround_up_valuelist.append(value)
                                        if str(j_y+5) in dict_yz:   value=float(dict_yz[str(j_y+5)]);     surround_up_list.append([j,value]);         surround_up_valuelist.append(value)
                                        if str(j_y+6) in dict_yz:   value=float(dict_yz[str(j_y+6)]);     surround_up_list.append([j,value]);         surround_up_valuelist.append(value)
                                        if str(j_y-4) in dict_yz:   value=float(dict_yz[str(j_y-4)]);     surround_down_list.append([j,value]);         surround_down_valuelist.append(value)
                                        if str(j_y-5) in dict_yz:   value=float(dict_yz[str(j_y-5)]);     surround_down_list.append([j,value]);         surround_down_valuelist.append(value)
                                        if str(j_y-6) in dict_yz:   value=float(dict_yz[str(j_y-6)]);     surround_down_list.append([j,value]);         surround_down_valuelist.append(value)  
                                        ###
                                        j+=1
                                    # UP    
                                    if len(middle_list)<5 or len(surround_up_list)<5:
                                        accuracy1="NA"
                                    else:    
                                        line_A = np.array(middle_list)
                                        line_B = np.array(surround_up_list)
                                        X = np.vstack((line_A, line_B))      # Merge data and create labels
                                        y = np.hstack((np.zeros(len(line_A)), np.ones(len(line_B))))
                                        lda = LinearDiscriminantAnalysis()# Train LDA model
                                        lda.fit(X, y)
                                        predictions = lda.predict(X) # Check model accuracy
                                        accuracy1 = np.mean(predictions == y)
                                        #print("Coefficients:", lda.coef_)  # Coefficients of the linear boundary
                                        #print("Intercept:", lda.intercept_)  # Intercept of the linear boundary
                                        #print("Accuracy:", accuracy)  # Classification accuracy
                                    # DOWN
                                    if len(middle_list)<5 or len(surround_down_list)<5:
                                        accuracy2="NA"
                                    else:                                            
                                        line_A = np.array(middle_list)
                                        line_B = np.array(surround_down_list)
                                        X = np.vstack((line_A, line_B))      # Merge data and create labels
                                        y = np.hstack((np.zeros(len(line_A)), np.ones(len(line_B))))
                                        lda = LinearDiscriminantAnalysis()# Train LDA model
                                        lda.fit(X, y)
                                        predictions = lda.predict(X) # Check model accuracy
                                        accuracy2 = np.mean(predictions == y)
                                    ####
                                    middle_mean = sum(middle_valuelist) / len(middle_valuelist)
                                    if len(surround_up_valuelist)==0:        up_distance_mean="NA"
                                    else:                               up_mean=sum(surround_up_valuelist) / len(surround_up_valuelist);          up_distance_mean= round(middle_mean-up_mean,2)
                                    if len(surround_down_valuelist)==0:      down_distance_mean="NA"
                                    else:                               down_mean=sum(surround_down_valuelist) / len(surround_down_valuelist);          down_distance_mean= round(middle_mean-down_mean,2)
                                    ####
                                    if accuracy1!='NA' and accuracy1<0.6:                          delete_dot_mark='no'
                                    elif accuracy2!='NA' and accuracy2<0.6:                          delete_dot_mark='no'
                                    elif up_distance_mean!='NA' and up_distance_mean<0.2:       delete_dot_mark='no'
                                    elif down_distance_mean!='NA' and down_distance_mean<0.2:   delete_dot_mark='no' 
                                    else:    
                                        fffff.write(f"{str(all_start_x)}\t{str(all_start_y)}\t{str(all_end_x)}\t{str(all_end_y)}\t{str(seq1_len)}\t{str(seq2_len)}\t{str(tmp_dot_list_len)}\t{str(all_x_num-tmp_dot_list_len)}\t{str(round(good_ratio,3))}\t{str(round(slope,2))}\t{str(round(r_squared,2))}\t{accuracy1}\t{accuracy2}\t{str(up_distance_mean)}\t{str(down_distance_mean)}\n")   #\t{dots_string}
           
                                        
                                ## Remove these points from the dictionary
                                if delete_dot_mark=='yes':
                                    for one_result in tmp_dot_list:
                                        one_x,one_y=one_result
                                        del dict3_x_y_z[one_x][one_y]
                        ## Initialize
                        tmp_dot_list=[]
                        #return False                
                        
                ### Calculate diagonally downward
            
            
    
    with Pool(processes=thread) as pool:
         # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step2, info_list), start=1):
            # Results can be processed here, e.g., stored or printed
            progress = (i / len(info_list)) * 100
            sys.stdout.write(f"\nProgress: {progress:.2f}% {str(i)} / {str(len(sample_list))}#####################################################################\n")
            sys.stdout.flush()    
             
    print("Summarizing step2 results")
    with open('./chr2moddotplot3/0_prepare/sample_list','r') as f:
        sample_list=f.read().split('\n')
        def chrname_old2new(one_id):
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
            else: one_id='other'
            return one_id
    with open("./chr2moddotplot3/2_slashin_updown",'w') as fout:
        fout.write(f"#sample\tchr\ttype\tstart_x\tstart_y\tend_x\tend_y\tseq1_len\tseq2_len\tdot_num\tdotskip_num\tdot_used\tline_slope\tline_R2\tback_accuracy1\tback_accuracy2\tback_distance1\tback_distance2\n")
        for one_sample in sample_list:
            input_dir=f"./chr2moddotplot3/2_slashin/{one_sample}/"
            files=os.listdir(input_dir)
            chr_list=[]
            for onefile in files: 
                if onefile[-12:]==".3_simplebed":
                    chr_name=onefile[:-12]
                    chr_list.append(chr_name)    
            for one_chr in chr_list:
                with open(f"{input_dir}/{one_chr}.4_line2up") as f2:  
                    one_chr_new=chrname_old2new(one_chr)
                    if one_chr_new=='other': continue
                    for line in f2:
                        if line[0]=='#':continue
                        eachline=line.strip()
                        fout.write(f"{one_sample}\t{one_chr_new}\tup\t{eachline}\n")
                        
                with open(f"{input_dir}/{one_chr}.4_line2down") as f3:
                    one_chr_new=chrname_old2new(one_chr)
                    if one_chr_new=='other': continue
                    for line in f3:
                        if line[0]=='#':continue
                        eachline=line.strip()
                        fout.write(f"{one_sample}\t{one_chr_new}\tdown\t{eachline}\n")
if argv1=="stepall" or argv1=="step2" or argv1=="step2b": 
    print('70s')
    import numpy as np
    #from scipy import stats
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    #from scipy.signal import find_peaks
    #from scipy.interpolate import interp1d
    #import matplotlib.pyplot as plt  
    
    print("Analyzing moddotplot bed files, extracting reverse sequences")
    if  os.path.exists('./chr2moddotplot3/2_slashout')==False:
        subprocess.run(["mkdir ./chr2moddotplot3/2_slashout"], shell=True)     

    ### Set parameters
    parafile='./chr2moddotplot3/2_slashout/parameter'
    match_identity=80   # Preliminary filtering of some points seems unnecessary since later we take the largest 10 anyway
    x_ydot_max_num=10        # Filter the largest 10 values for each x position
    x_ydot_max_add_bynum=0.01    # If >100, filter the largest 11 values, etc.
    x_ydot_max_add_byvalue=0.3   # identity, to prevent repetitive sequences
    x_ydot_region_max=5000 # Prevent overly distant matches, max is 5000
    #
    x_ydot_neighbor=5    # At an x position, allowed gap tolerance is 5, meaning the minimum distance between two y points on the same x will not be less than 5
    x_gap_max=5         # x_ydot_neighbor seems related; after x moves right by +1, y can only have one point within ±5,      # Horizontal search
    y_gap_max=15                                                                                    # Vertical search
    extend_tolerance=30
    #
    with open('./chr2moddotplot3/2_slashout/parameter','w') as f:
        f.write(f"match_identity:{match_identity}\n")   
        f.write(f"x_ydot_max_num:{x_ydot_max_num}\n")
        f.write(f"x_ydot_max_add_bynum:{x_ydot_max_add_bynum}\n")
        f.write(f"x_ydot_max_add_byvalue:{x_ydot_max_add_byvalue}\n")
        f.write(f"x_ydot_region_max:{x_ydot_region_max}\n")  
        f.write(f"x_ydot_neighbor:{x_ydot_neighbor}\n")  
        f.write(f"x_gap_max:{x_gap_max}\n")     
        f.write(f"y_gap_max:{y_gap_max}\n")    
        f.write(f"extend_tolerance:{extend_tolerance}\n")    
    
    info_list=[]
    with open('./chr2moddotplot3/0_prepare/sample_list','r') as f:sample_list=f.read().split('\n')
    for sample_name in sample_list:
        input_dir=f"./chr2moddotplot3/1_VSat12split/{sample_name}"
        output_dir=f"./chr2moddotplot3/2_slashout/{sample_name}"        
        info_list.append([sample_name,input_dir,output_dir,parafile])        
    

        
    # Assign tasks to processes in the process pool    
    def run_step2(one_info):
        sample_name,input_dir,output_dir,parafile=one_info
        if len(sample_name)==0: return False
        #if sample_name!='PN40024_hap1': return False
        print(sample_name)
        
        subprocess.run([f"mkdir {output_dir}"], shell=True)     
        files=os.listdir(input_dir)
        
        ## Load parameters
        with open(parafile,'r') as f:
            parameters=f.read().split('\n')
            for one_parameter in parameters:
                if "match_identity:" in one_parameter: match_identity=float(one_parameter.split('match_identity:')[1])
                if "x_ydot_max_num:" in one_parameter: x_ydot_max_num=int(one_parameter.split('x_ydot_max_num:')[1])
                if "x_ydot_max_add_bynum:" in one_parameter: x_ydot_max_add_bynum=float(one_parameter.split('x_ydot_max_add_bynum:')[1])
                if "x_ydot_max_add_byvalue:" in one_parameter: x_ydot_max_add_byvalue=float(one_parameter.split('x_ydot_max_add_byvalue:')[1])
                if "x_ydot_neighbor:" in one_parameter: x_ydot_neighbor=int(one_parameter.split('x_ydot_neighbor:')[1])
                if "x_ydot_region_max:" in one_parameter: x_ydot_region_max=int(one_parameter.split('x_ydot_region_max:')[1])
                if "x_gap_max:" in one_parameter: x_gap_max=int(one_parameter.split('x_gap_max:')[1])
                if "y_gap_max:" in one_parameter: y_gap_max=int(one_parameter.split('y_gap_max:')[1])
                if "extend_tolerance:" in one_parameter: extend_tolerance=int(one_parameter.split('extend_tolerance:')[1])
                
        #print(f"match_identity:{str(match_identity)}")        
        #print(f"x_ydot_max_num:{str(x_ydot_max_num)}")
        #print(f"x_ydot_max_add_bynum:{str(x_ydot_max_add_bynum)}")
        #print(f"x_ydot_max_add_byvalue:{str(x_ydot_max_add_byvalue)}")
        #print(f"x_ydot_neighbor:{str(x_ydot_neighbor)}")
        #print(f"x_ydot_region_max:{str(x_ydot_region_max)}")
        #print(f"x_gap_max:{str(x_gap_max)}")
        #print(f"y_gap_max:{str(y_gap_max)}")
        #print(f"extend_tolerance:{str(extend_tolerance)}")
        
        
        #########################
        chr_name_list=[]
        for x in files:
            if x[-14:]=='.bed.VSat1_out':
                chr_name=x[:-14]
                chr_name_list.append(chr_name)
        chr_name_list.sort()        
        #print(chr_name_list)
        ########################
        
        for one_chr in chr_name_list:   
            #print(one_chr)
            ### Load bed file into dictionary dict_x_y_z
            dict_x_y_z={}
            x_list=[]
            with open(f'{output_dir}/{one_chr}.1_newbed','w') as ff:
                #ff.write(f"pos1\tpos2\tperID_by_events\n")   
                unit_len=''
                with open(f"{input_dir}/{one_chr}.bed.VSat1_out") as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        if eachline_arr[0]=='#query_name':continue
                        if unit_len=='':  
                            unit_len=int(eachline_arr[2])-int(eachline_arr[1])+1
                            ff.write(f"#unit_len={unit_len}\n")
  
                        ##query_name	query_start	query_end	reference_name	reference_start	reference_end	perID_by_events
                        ##BMNGT2T_hap1_chr02	1	21454	BMNGT2T_hap1_chr02	1	21454	100
                        new_pos1 = int((int(eachline_arr[1])-1)/unit_len);          new_pos1_str=str(new_pos1);              x_list.append(new_pos1)
                        new_pos2 = int((int(eachline_arr[4])-1)/unit_len);          new_pos2_str=str(new_pos2);
                        percent_str  = eachline_arr[6];                             percent_float = float(percent_str)
                        #if new_pos2-new_pos1>500:continue    #20000*200 roughly 4M
                        ff.write(f"{new_pos1_str}\t{new_pos2_str}\t{percent_str}\n")
                        
                        if new_pos2-new_pos1>x_ydot_region_max:continue
                        if percent_float<match_identity:continue 
                        if new_pos1_str not in dict_x_y_z: dict_x_y_z[new_pos1_str]={}
                        dict_x_y_z[new_pos1_str][new_pos2_str]=percent_float
            
            ### Filter out some points, keep the top 10 (variable) points for each x position
            dict2_x_y_z={}
            with open(f'{output_dir}/{one_chr}.2_simplebed','w') as fff: 
                fff.write(f"#unit_len={unit_len}\n")
                for new_pos1_str,new_pos2_str_dict in dict_x_y_z.items():
                    # Sort dictionary values in descending order
                    x2y_valuelist = sorted(new_pos2_str_dict.values(), key=lambda item: float(item), reverse=True)
                    
                    # Fixed number + (total count / 100)
                    keep_num=x_ydot_max_num+round(len(x2y_valuelist) * float(x_ydot_max_add_bynum))
                    # Find the minimum value, subtract x_ydot_max_add_byvalue
                    x2y_valuelist_min=min(x2y_valuelist[:keep_num-1])-x_ydot_max_add_byvalue

                    # Create a new dictionary with values greater than x2y_valuelist_min
                    filtered_dict = {key: value for key, value in new_pos2_str_dict.items() if float(value) > x2y_valuelist_min}
                    
                    for new_pos2_str,z in filtered_dict.items():
                        fff.write(f"{new_pos1_str}\t{new_pos2_str}\t{z}\n")
                    #
                    dict2_x_y_z[new_pos1_str]=filtered_dict
                    
            ### Filter out some points more. For the remaining points, merge by position. Points with neighbor distance <5 (variable) are considered the same peak, keep the max value
            with open(f'{output_dir}/{one_chr}.3_simplebed','w') as ffff: 
                ffff.write(f"#unit_len={unit_len}\n")
                for new_pos1_str,new_pos2_str_dict in dict2_x_y_z.items():
                    # Get keys as list sorted by value in descending order
                    sorted_keys = sorted(new_pos2_str_dict.keys(),key=lambda key: int(key), reverse=True  )
                    sorted_keys_list = list(sorted_keys)
                    #print(sorted_keys_list);
                    goodkey='';key_edge='';goodvalue='';
                    for new_pos2_str in sorted_keys_list:
                        value=float(new_pos2_str_dict[new_pos2_str])
                        #print(value)
                        if goodkey=='': goodkey=new_pos2_str;     key_edge=new_pos2_str;      goodvalue=value
                        #print(new_pos2_str,int(key_edge)+x_ydot_neighbor)
                        if int(new_pos2_str)>=int(key_edge)-x_ydot_neighbor:   
                            key_edge=new_pos2_str
                            if goodvalue<=value:     goodkey=new_pos2_str;    goodvalue=value;           ## Why have =? Mainly because near the diagonal there may be two 100s, ensuring 100 is the lowest one, i.e., the diagonal
                        else:    
                            if new_pos1_str!=goodkey:
                                ffff.write(f"{new_pos1_str}\t{goodkey}\t{str(goodvalue)}\n")
                            #print(value)
                            goodkey=new_pos2_str;key_edge=new_pos2_str;goodvalue=value;
                    if goodkey!='':   
                        if new_pos1_str!=goodkey:
                            ffff.write(f"{new_pos1_str}\t{goodkey}\t{str(goodvalue)}\n")
                    goodkey='';key_edge='';goodvalue='';
                    
            subprocess.run([f"rm {output_dir}/{one_chr}.1_newbed"], shell=True)              
            subprocess.run([f"rm {output_dir}/{one_chr}.2_simplebed"], shell=True) 
            
            
            # Point pool 
            dict3_x_y_z={}  ## This is a numeric dictionary xy
            with open(f'{output_dir}/{one_chr}.3_simplebed','r') as f:
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if len(eachline_arr)<3:continue
                    pos1,pos2,value=eachline_arr
                    pos1=int(pos1)
                    pos2=int(pos2)
                    if pos1 not in dict3_x_y_z: dict3_x_y_z[pos1]={}
                    dict3_x_y_z[pos1][pos2]=value
            x_list=sorted(dict3_x_y_z.keys())       

            ### Calculate diagonally upward
            with open(f'{output_dir}/{one_chr}.4_line2up','w') as fffff:
                fffff.write(f"#unit_len={unit_len}\n")
                fffff.write(f"#start_x\tstart_y\tend_x\tend_y\tseq1_len\tseq2_len\tdot_num\tdotskip_num\tdot_used\tline_slope\tline_R2\tback_accuracy1\tback_accuracy2\tback_distance1\tback_distance2\n")  #\tdots_string
                tmp_dot_list=[]
                
                for x in x_list:
                    y_list=sorted(dict3_x_y_z[x].keys(), reverse=True)     ####y_list=list(dict3_x_y_z[x].keys);  y_list.sort()   
                    #print(y_list)
                    for y in y_list:
                        # Current point is (x,y)
                        x_start=x
                        y_start=y
                        end_mark=''
                        tolerance=0
                        #print(x,y,type(dict3_x_y_z))
                        serial=0;
                        while end_mark=='' and tolerance<extend_tolerance:
                            ## Check if the current start point has been used in previous loops
                            if x_start not in dict3_x_y_z:              end_mark='yes';break
                            if y_start not in dict3_x_y_z[x_start]:     end_mark='yes';break
                            xy_distance=y_start-x_start
                            ##
                            if tmp_dot_list==[]:tmp_dot_list.append([x,y]);serial+=1
                            # Collect the next point
                            i1=0;good_mark=''
                            while i1<x_gap_max:
                                i1+=1
                                
                                if good_mark=='yes':break
                                if serial<3:i2=1   # Slope of the first three points cannot be horizontal
                                elif xy_distance<20:    i2=1        # The smallest duplication might be 100k, i.e., 20 points
                                elif xy_distance<40:    i2=0
                                elif xy_distance<90:    i2=-1
                                else :                  i2=-2
                                while i2<y_gap_max:
                                    x_tmp=x_start+i1
                                    y_tmp=y_start+i2  ### Move x one step right, y first stays (+0) then moves upward (++1)
                                    i2+=1
                                    if x_tmp in dict3_x_y_z:
                                        if y_tmp in dict3_x_y_z[x_tmp]:
                                            # New coordinates are (x_tmp,y_tmp)
                                            good_mark='yes';break
                            ## If the next point is collected, continue with it as the starting point            
                            if good_mark=='yes': 
                                tmp_dot_list.append([x_tmp,y_tmp])
                                x_start=x_tmp
                                y_start=y_tmp
                                serial+=1
                            ## If no next point is collected, end the current line segment    
                            else:
                                tolerance+=1
                                x_start+=1
                                y_start+=1  ##
                                
                      
                        ### Determine the number of points in the currently collected point sequence, slope, and fit. [[x1,y1],[x2,y2],...], calculate the slope and R² of the fitted line for all points
                        tmp_dot_list_len=len(tmp_dot_list)
                        if tmp_dot_list_len>5:
                            points = np.array(tmp_dot_list)
                            np_x = points[:, 0]    # Extract x and y values
                            np_y = points[:, 1]
                            slope, intercept = np.polyfit(np_x, np_y, 1)  # Use numpy's polyfit function for linear regression, 1 indicates first-order linear fit
                            #print(f"Fitted line slope: {slope}")
                            #print(f"Fitted line intercept: {intercept}")
                            y_pred = slope * np_x + intercept  # Calculate R² goodness of fit
                            ss_total = np.sum((np_y - np.mean(np_y))**2)
                            ss_res = np.sum((np_y - y_pred)**2)
                            if ss_total == 0:       r_squared = 0  # or other default value
                            else:                   r_squared = 1 - (ss_res / ss_total)
                            #print(f"R² goodness of fit: {r_squared}")
                            if r_squared>0.5 and abs(slope)>0.7 and abs(slope)<1.42 and tmp_dot_list_len>5:
                                all_start_x=tmp_dot_list[0][0]*unit_len+1   ;               seq1_start=all_start_x     #seq1_start
                                all_start_y=tmp_dot_list[0][1]*unit_len+1+(unit_len-1)   ;  seq2_start=all_start_y     #seq2_start     #+(unit_len-1) because the simplified start position was used
                                all_end_x=tmp_dot_list[-1][0]*unit_len+1    ;               seq1_end=all_end_x      #seq1_end
                                all_end_y=tmp_dot_list[-1][1]*unit_len+1+(unit_len-1)    ;  seq2_end=all_end_y      #seq2_end
                                seq1_len=abs(seq1_end-seq1_start)+1
                                seq2_len=abs(seq2_end-seq2_start)+1
                                dots_string = '|'.join(f'{x},{y}' for x, y in tmp_dot_list)
                                #if dots_string=="2777,2815|2781,2817|2785,2817|2786,2817|2787,2817|2788,2817|2791,2820|2795,2827|2800,2831|2801,2839|2802,2839|2803,2839|2804,2848|2805,2849|2806,2849|2810,2849|2811,2849|2815,2849|2816,2849|2818,2849|2820,2849|2821,2849|2822,2849":print(points)
                                
                                all_x_num=tmp_dot_list[-1][0]-tmp_dot_list[0][0]+1
                                good_ratio=tmp_dot_list_len/all_x_num
                                
                                delete_dot_mark='yes'
                                if tmp_dot_list_len<3000 and good_ratio>0.7:
                                    #### Calculate the proportion of the line and the surrounding circle, -1~+1 width 3, -3~-5 and +3~+5 width 3,
                                    dict_tmp_dot_list = {int(item[0]): int(item[1]) for item in tmp_dot_list}
                                    line_start_serial=tmp_dot_list[0][0]
                                    line_end_serial=tmp_dot_list[-1][0]
                                    j=line_start_serial
                                    middle_list=[];surround_up_list=[];surround_down_list=[]
                                    middle_valuelist=[];surround_up_valuelist=[];surround_down_valuelist=[]
                                    #print(line_start_serial,line_end_serial)
                                    while j<line_end_serial:
                                        if str(j) in dict_x_y_z:        dict_yz=dict_x_y_z[str(j)];##print(dict_yz)
                                        
                                        if j in dict_tmp_dot_list:
                                            j_y=dict_tmp_dot_list[j]
                                        else:
                                            j_y+=1
                                        ###    
                                        if str(j_y) in dict_yz:     value=float(dict_yz[str(j_y)]);       middle_list.append([j,value]);         middle_valuelist.append(value)
                                        if str(j_y-1) in dict_yz:   value=float(dict_yz[str(j_y-1)]);     middle_list.append([j,value]);         middle_valuelist.append(value)
                                        if str(j_y+1) in dict_yz:   value=float(dict_yz[str(j_y+1)]);     middle_list.append([j,value]);         middle_valuelist.append(value)
                                        if str(j_y+4) in dict_yz:   value=float(dict_yz[str(j_y+4)]);     surround_up_list.append([j,value]);         surround_up_valuelist.append(value)
                                        if str(j_y+5) in dict_yz:   value=float(dict_yz[str(j_y+5)]);     surround_up_list.append([j,value]);         surround_up_valuelist.append(value)
                                        if str(j_y+6) in dict_yz:   value=float(dict_yz[str(j_y+6)]);     surround_up_list.append([j,value]);         surround_up_valuelist.append(value)
                                        if str(j_y-4) in dict_yz:   value=float(dict_yz[str(j_y-4)]);     surround_down_list.append([j,value]);         surround_down_valuelist.append(value)
                                        if str(j_y-5) in dict_yz:   value=float(dict_yz[str(j_y-5)]);     surround_down_list.append([j,value]);         surround_down_valuelist.append(value)
                                        if str(j_y-6) in dict_yz:   value=float(dict_yz[str(j_y-6)]);     surround_down_list.append([j,value]);         surround_down_valuelist.append(value)  
                                        ###
                                        j+=1
                                    # UP    
                                    if len(middle_list)<5 or len(surround_up_list)<5:
                                        accuracy1="NA"
                                    else:    
                                        line_A = np.array(middle_list)
                                        line_B = np.array(surround_up_list)
                                        X = np.vstack((line_A, line_B))      # Merge data and create labels
                                        y = np.hstack((np.zeros(len(line_A)), np.ones(len(line_B))))
                                        lda = LinearDiscriminantAnalysis()# Train LDA model
                                        lda.fit(X, y)
                                        predictions = lda.predict(X) # Check model accuracy
                                        accuracy1 = np.mean(predictions == y)
                                        #print("Coefficients:", lda.coef_)  # Coefficients of the linear boundary
                                        #print("Intercept:", lda.intercept_)  # Intercept of the linear boundary
                                        #print("Accuracy:", accuracy)  # Classification accuracy
                                    # DOWN
                                    if len(middle_list)<5 or len(surround_down_list)<5:
                                        accuracy2="NA"
                                    else:                                            
                                        line_A = np.array(middle_list)
                                        line_B = np.array(surround_down_list)
                                        X = np.vstack((line_A, line_B))      # Merge data and create labels
                                        y = np.hstack((np.zeros(len(line_A)), np.ones(len(line_B))))
                                        lda = LinearDiscriminantAnalysis()# Train LDA model
                                        lda.fit(X, y)
                                        predictions = lda.predict(X) # Check model accuracy
                                        accuracy2 = np.mean(predictions == y)
                                    ####
                                    middle_mean = sum(middle_valuelist) / len(middle_valuelist)
                                    if len(surround_up_valuelist)==0:        up_distance_mean="NA"
                                    else:                               up_mean=sum(surround_up_valuelist) / len(surround_up_valuelist);          up_distance_mean= round(middle_mean-up_mean,2)
                                    if len(surround_down_valuelist)==0:      down_distance_mean="NA"
                                    else:                               down_mean=sum(surround_down_valuelist) / len(surround_down_valuelist);          down_distance_mean= round(middle_mean-down_mean,2)
                                    ####
                                    if accuracy1!='NA' and accuracy1<0.6:                          delete_dot_mark='no'
                                    elif accuracy2!='NA' and accuracy2<0.6:                          delete_dot_mark='no'
                                    elif up_distance_mean!='NA' and up_distance_mean<0.2:       delete_dot_mark='no'
                                    elif down_distance_mean!='NA' and down_distance_mean<0.2:   delete_dot_mark='no' 
                                    else:    
                                        fffff.write(f"{str(all_start_x)}\t{str(all_start_y)}\t{str(all_end_x)}\t{str(all_end_y)}\t{str(seq1_len)}\t{str(seq2_len)}\t{str(tmp_dot_list_len)}\t{str(all_x_num-tmp_dot_list_len)}\t{str(round(good_ratio,3))}\t{str(round(slope,2))}\t{str(round(r_squared,2))}\t{accuracy1}\t{accuracy2}\t{str(up_distance_mean)}\t{str(down_distance_mean)}\n")   #\t{dots_string}
           
                                        
                                ## Remove these points from the dictionary
                                if delete_dot_mark=='yes':
                                    for one_result in tmp_dot_list:
                                        one_x,one_y=one_result
                                        del dict3_x_y_z[one_x][one_y]
                        ## Initialize
                        tmp_dot_list=[]
                        #return False                
                        
                ### Calculate diagonally downward
            
            ### Calculate diagonally downward
            with open(f'{output_dir}/{one_chr}.4_line2down','w') as fffff:
                fffff.write(f"#unit_len={unit_len}\n")
                fffff.write(f"#start_x\tstart_y\tend_x\tend_y\tseq1_len\tseq2_len\tdot_num\tdotskip_num\tdot_used\tline_slope\tline_R2\tback_accuracy1\tback_accuracy2\tback_distance1\tback_distance2\n") #\tdots_string
                tmp_dot_list=[]
                
                for x in x_list:
                    y_list=sorted(dict3_x_y_z[x].keys(), reverse=True)     ####y_list=list(dict3_x_y_z[x].keys);  y_list.sort()   
                    #print(y_list)
                    for y in y_list:
                        # Current point is (x,y)
                        x_start=x
                        y_start=y
                        end_mark=''
                        tolerance=0
                        #print(x,y,type(dict3_x_y_z))
                        serial=0;
                        while end_mark=='' and tolerance<extend_tolerance:
                            ## Check if the current start point has been used in previous loops
                            if x_start not in dict3_x_y_z:              end_mark='yes';break
                            if y_start not in dict3_x_y_z[x_start]:     end_mark='yes';break
                            xy_distance=y_start-x_start
                            ##
                            if tmp_dot_list==[]:tmp_dot_list.append([x,y]);serial+=1
                            # Collect the next point
                            i1=0;good_mark=''
                            while i1<x_gap_max:
                                i1+=1
                                
                                if good_mark=='yes':break
                                if serial<3:i2=1   # Slope of the first three points cannot be horizontal
                                elif xy_distance<20:    i2=1        # The smallest duplication might be 100k, i.e., 20 points
                                elif xy_distance<40:    i2=0
                                elif xy_distance<90:    i2=-1
                                else :                  i2=-2
                                #i2=1
                                while i2<y_gap_max:
                                    x_tmp=x_start+i1
                                    y_tmp=y_start-i2  ### Move x one step right, y first stays (+0) then moves upward (++1)             !!!!!!!!!!!!!!!!!!!!!!!!! Change direction here
                                    i2+=1
                                    if x_tmp in dict3_x_y_z:
                                        if y_tmp in dict3_x_y_z[x_tmp]:
                                            # New coordinates are (x_tmp,y_tmp)
                                            good_mark='yes';break
                            ## If the next point is collected, continue with it as the starting point            
                            if good_mark=='yes': 
                                tmp_dot_list.append([x_tmp,y_tmp])
                                x_start=x_tmp
                                y_start=y_tmp
                                serial+=1
                            ## If no next point is collected, end the current line segment    
                            else:
                                tolerance+=1
                                x_start+=1
                                y_start-=1  ##                                                           !!!!!!!!!!!!!!!!!!!!!!!!! Change direction here
                                
                      
                        ### Determine the number of points in the currently collected point sequence, slope, and fit. [[x1,y1],[x2,y2],...], calculate the slope and R² of the fitted line for all points
                        tmp_dot_list_len=len(tmp_dot_list)
                        if tmp_dot_list_len>5:
                            points = np.array(tmp_dot_list)
                            np_x = points[:, 0]    # Extract x and y values
                            np_y = points[:, 1]
                            slope, intercept = np.polyfit(np_x, np_y, 1)  # Use numpy's polyfit function for linear regression, 1 indicates first-order linear fit
                            #print(f"Fitted line slope: {slope}")
                            #print(f"Fitted line intercept: {intercept}")
                            y_pred = slope * np_x + intercept  # Calculate R² goodness of fit
                            ss_total = np.sum((np_y - np.mean(np_y))**2)
                            ss_res = np.sum((np_y - y_pred)**2)
                            if ss_total == 0:       r_squared = 0  # or other default value
                            else:                   r_squared = 1 - (ss_res / ss_total)
                            #print(f"R² goodness of fit: {r_squared}")
                            if r_squared>0.5 and abs(slope)>0.7 and abs(slope)<1.42 and tmp_dot_list_len>5:
                                all_start_x=tmp_dot_list[0][0]*unit_len+1   ;               seq1_start=all_start_x     #seq1_start
                                all_start_y=tmp_dot_list[0][1]*unit_len+1+(unit_len-1)   ;  seq2_start=all_start_y     #seq2_start     #+(unit_len-1) because the simplified start position was used
                                all_end_x=tmp_dot_list[-1][0]*unit_len+1    ;               seq1_end=all_end_x      #seq1_end
                                all_end_y=tmp_dot_list[-1][1]*unit_len+1+(unit_len-1)    ;  seq2_end=all_end_y      #seq2_end
                                seq1_len=abs(seq1_end-seq1_start)+1
                                seq2_len=abs(seq2_end-seq2_start)+1
                                dots_string = '|'.join(f'{x},{y}' for x, y in tmp_dot_list)
                                #if dots_string=="2777,2815|2781,2817|2785,2817|2786,2817|2787,2817|2788,2817|2791,2820|2795,2827|2800,2831|2801,2839|2802,2839|2803,2839|2804,2848|2805,2849|2806,2849|2810,2849|2811,2849|2815,2849|2816,2849|2818,2849|2820,2849|2821,2849|2822,2849":print(points)
                                
                                all_x_num=tmp_dot_list[-1][0]-tmp_dot_list[0][0]+1
                                good_ratio=tmp_dot_list_len/all_x_num
                                
                                delete_dot_mark='yes'
                                if tmp_dot_list_len<3000 and good_ratio>0.7:
                                    #### Calculate the proportion of the line and the surrounding circle, -1~+1 width 3, -3~-5 and +3~+5 width 3,
                                    dict_tmp_dot_list = {int(item[0]): int(item[1]) for item in tmp_dot_list}
                                    line_start_serial=tmp_dot_list[0][0]
                                    line_end_serial=tmp_dot_list[-1][0]
                                    j=line_start_serial
                                    middle_list=[];surround_up_list=[];surround_down_list=[]
                                    middle_valuelist=[];surround_up_valuelist=[];surround_down_valuelist=[]
                                    #print(line_start_serial,line_end_serial)
                                    while j<line_end_serial:
                                        if str(j) in dict_x_y_z:        dict_yz=dict_x_y_z[str(j)];##print(dict_yz)
                                        
                                        if j in dict_tmp_dot_list:
                                            j_y=dict_tmp_dot_list[j]
                                        else:
                                            j_y+=1
                                        ###    
                                        if str(j_y) in dict_yz:     value=float(dict_yz[str(j_y)]);       middle_list.append([j,value]);         middle_valuelist.append(value)
                                        if str(j_y-1) in dict_yz:   value=float(dict_yz[str(j_y-1)]);     middle_list.append([j,value]);         middle_valuelist.append(value)
                                        if str(j_y+1) in dict_yz:   value=float(dict_yz[str(j_y+1)]);     middle_list.append([j,value]);         middle_valuelist.append(value)
                                        if str(j_y+4) in dict_yz:   value=float(dict_yz[str(j_y+4)]);     surround_up_list.append([j,value]);         surround_up_valuelist.append(value)
                                        if str(j_y+5) in dict_yz:   value=float(dict_yz[str(j_y+5)]);     surround_up_list.append([j,value]);         surround_up_valuelist.append(value)
                                        if str(j_y+6) in dict_yz:   value=float(dict_yz[str(j_y+6)]);     surround_up_list.append([j,value]);         surround_up_valuelist.append(value)
                                        if str(j_y-4) in dict_yz:   value=float(dict_yz[str(j_y-4)]);     surround_down_list.append([j,value]);         surround_down_valuelist.append(value)
                                        if str(j_y-5) in dict_yz:   value=float(dict_yz[str(j_y-5)]);     surround_down_list.append([j,value]);         surround_down_valuelist.append(value)
                                        if str(j_y-6) in dict_yz:   value=float(dict_yz[str(j_y-6)]);     surround_down_list.append([j,value]);         surround_down_valuelist.append(value)  
                                        ###
                                        j+=1
                                    # UP    
                                    if len(middle_list)<5 or len(surround_up_list)<5:
                                        accuracy1="NA"
                                    else:    
                                        line_A = np.array(middle_list)
                                        line_B = np.array(surround_up_list)
                                        X = np.vstack((line_A, line_B))      # Merge data and create labels
                                        y = np.hstack((np.zeros(len(line_A)), np.ones(len(line_B))))
                                        lda = LinearDiscriminantAnalysis()# Train LDA model
                                        lda.fit(X, y)
                                        predictions = lda.predict(X) # Check model accuracy
                                        accuracy1 = np.mean(predictions == y)
                                        #print("Coefficients:", lda.coef_)  # Coefficients of the linear boundary
                                        #print("Intercept:", lda.intercept_)  # Intercept of the linear boundary
                                        #print("Accuracy:", accuracy)  # Classification accuracy
                                    # DOWN
                                    if len(middle_list)<5 or len(surround_down_list)<5:
                                        accuracy2="NA"
                                    else:                                            
                                        line_A = np.array(middle_list)
                                        line_B = np.array(surround_down_list)
                                        X = np.vstack((line_A, line_B))      # Merge data and create labels
                                        y = np.hstack((np.zeros(len(line_A)), np.ones(len(line_B))))
                                        lda = LinearDiscriminantAnalysis()# Train LDA model
                                        lda.fit(X, y)
                                        predictions = lda.predict(X) # Check model accuracy
                                        accuracy2 = np.mean(predictions == y)
                                    ####
                                    middle_mean = sum(middle_valuelist) / len(middle_valuelist)
                                    if len(surround_up_valuelist)==0:        up_distance_mean="NA"
                                    else:                               up_mean=sum(surround_up_valuelist) / len(surround_up_valuelist);          up_distance_mean= round(middle_mean-up_mean,2)
                                    if len(surround_down_valuelist)==0:      down_distance_mean="NA"
                                    else:                               down_mean=sum(surround_down_valuelist) / len(surround_down_valuelist);          down_distance_mean= round(middle_mean-down_mean,2)
                                    ####
                                    if accuracy1!='NA' and accuracy1<0.6:                          delete_dot_mark='no'
                                    elif accuracy2!='NA' and accuracy2<0.6:                          delete_dot_mark='no'
                                    elif up_distance_mean!='NA' and up_distance_mean<0.2:       delete_dot_mark='no'
                                    elif down_distance_mean!='NA' and down_distance_mean<0.2:   delete_dot_mark='no' 
                                    else:    
                                        fffff.write(f"{str(all_start_x)}\t{str(all_start_y)}\t{str(all_end_x)}\t{str(all_end_y)}\t{str(seq1_len)}\t{str(seq2_len)}\t{str(tmp_dot_list_len)}\t{str(all_x_num-tmp_dot_list_len)}\t{str(round(good_ratio,3))}\t{str(round(slope,2))}\t{str(round(r_squared,2))}\t{accuracy1}\t{accuracy2}\t{str(up_distance_mean)}\t{str(down_distance_mean)}\n")   #\t{dots_string}
           
                                        
                                ## Remove these points from the dictionary
                                if delete_dot_mark=='yes':
                                    for one_result in tmp_dot_list:
                                        one_x,one_y=one_result
                                        del dict3_x_y_z[one_x][one_y]
                        ## Initialize
                        tmp_dot_list=[]
                        #return False                
                        
                ### Calculate diagonally downward
            
            
    
    with Pool(processes=thread) as pool:
         # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step2, info_list), start=1):
            # Results can be processed here, e.g., stored or printed
            progress = (i / len(info_list)) * 100
            sys.stdout.write(f"\nProgress: {progress:.2f}% {str(i)} / {str(len(sample_list))}#####################################################################\n")
            sys.stdout.flush()    
               
    print("Summarizing step2 results")
    with open('./chr2moddotplot3/0_prepare/sample_list','r') as f:
        sample_list=f.read().split('\n')
        def chrname_old2new(one_id):
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
            else: one_id='other'
            return one_id
    with open("./chr2moddotplot3/2_slashout_updown",'w') as fout:
        fout.write(f"#sample\tchr\ttype\tstart_x\tstart_y\tend_x\tend_y\tseq1_len\tseq2_len\tdot_num\tdotskip_num\tdot_used\tline_slope\tline_R2\tback_accuracy1\tback_accuracy2\tback_distance1\tback_distance2\n")
        for one_sample in sample_list:
            input_dir=f"./chr2moddotplot3/2_slashout/{one_sample}/"
            files=os.listdir(input_dir)
            chr_list=[]
            for onefile in files: 
                if onefile[-12:]==".3_simplebed":
                    chr_name=onefile[:-12]
                    chr_list.append(chr_name)    
            for one_chr in chr_list:
                with open(f"{input_dir}/{one_chr}.4_line2up") as f2:  
                    one_chr_new=chrname_old2new(one_chr)
                    if one_chr_new=='other': continue
                    for line in f2:
                        if line[0]=='#':continue
                        eachline=line.strip()
                        fout.write(f"{one_sample}\t{one_chr_new}\tup\t{eachline}\n")
                        
                with open(f"{input_dir}/{one_chr}.4_line2down") as f3:
                    one_chr_new=chrname_old2new(one_chr)
                    if one_chr_new=='other': continue
                    for line in f3:
                        if line[0]=='#':continue
                        eachline=line.strip()
                        fout.write(f"{one_sample}\t{one_chr_new}\tdown\t{eachline}\n")                        
if argv1=="stepall" or argv1=="step2_combine" or argv1=="step2a_combine":  
    print('Merge and extend sequences')
    print("Calculate the internal length of a single chromosome, i.e., the length of all alignment lines, only count sequences inside VSat1)")
    ## Simply load the VSat1 length for each chromosome
    dict_samplechr_VSat1len={}
    with open("/home/lain/aaa_data/run0/samples_satellite/2_good_regions_interarray",'r') as f:
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            samplechr=eachline_arr[0]+"_"+eachline_arr[1]
            dict_samplechr_VSat1len[samplechr]=int(eachline_arr[4])
    ##
    # Load sample information
    dict_sample_speciesgroup={}
    with open("/home/lain/aaa_data/run0/samples_satellite/sample_info",'r') as f:
        #Baimunage_hap2	Baimunage	 vinifera	Eurasian	Eurasian_vinifera_Baimunage	Table
        for line in f:
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=6:continue
            sample=             eachline_arr[0]
            species_group=      eachline_arr[3]
            dict_sample_speciesgroup[sample]=species_group
    ##  
    dict_sample_chr_type_info={}
    with open("./chr2moddotplot3/2_slashin_updown",'r') as f:   
        #("sample\tchromosome\ttype\tstart_1(x1)\tstart_2(y1)\tend_1(x2)\tend_2(y2)\tseq1_len\tseq2_len\tdot_num\tdotskip_num\tdot_used\tline_slope\tline_R2\tback_accuracy1\tback_accuracy2\tback_distance1\tback_distance2\tdistance\tinclude_length\tclass2VSat1\tunit_length\tunit_number\tclass2VSat1_revised\tmark\n")
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            sample=             eachline_arr[0]
            chromosome=         eachline_arr[1]
            one_type=         eachline_arr[2]
            alignmentline_length=         int(eachline_arr[7])+int(eachline_arr[8])  # Sum of the mapped genomic lengths of the alignment line, i.e., sum of seq1 and seq2
            if sample not in dict_sample_chr_type_info:                  dict_sample_chr_type_info[sample]={}
            if chromosome not in dict_sample_chr_type_info[sample]:      dict_sample_chr_type_info[sample][chromosome]={}
            if one_type not in dict_sample_chr_type_info[sample][chromosome]:
                dict_sample_chr_type_info[sample][chromosome][one_type]={}
                dict_sample_chr_type_info[sample][chromosome][one_type]['length']=0
                dict_sample_chr_type_info[sample][chromosome][one_type]['num']=0 
                dict_sample_chr_type_info[sample][chromosome][one_type]['pos']=[]
                dict_sample_chr_type_info[sample][chromosome][one_type]['pos_str']=[]
            dict_sample_chr_type_info[sample][chromosome][one_type]['length']+=alignmentline_length
            dict_sample_chr_type_info[sample][chromosome][one_type]['num']+=1  
            dict_sample_chr_type_info[sample][chromosome][one_type]['pos'].append([int(eachline_arr[3]),int(eachline_arr[4]),int(eachline_arr[5]),int(eachline_arr[6])])
            dict_sample_chr_type_info[sample][chromosome][one_type]['pos_str'].append(f"{eachline_arr[3]}-{eachline_arr[4]}-{eachline_arr[5]}-{eachline_arr[6]}")
    used_samples= ["Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap2","V059.hap1","V063.hap2","V063.hap1","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap2","Chardonnay_hap1","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap2","V074.hap1","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap2","V062.hap1","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap2","PinotNoir_hap1","PinotNoir2_hap2","PinotNoir2_hap1","V069.hap2","V069.hap1","PN40024","PN40024_hap1","PN40024_hap2","V087.hap2","V087.hap1","V061.hap2","V061.hap1","V065.hap2","V065.hap1","V064.hap2","V064.hap1","ThompsonSeedless_hap1","ThompsonSeedless_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap2","V067.hap1","V058.hap1","V058.hap2","V107.hap2","V107.hap1","V108.hap2","V108.hap1","V081.hap2","V081.hap1","V015.hap2","V015.hap1","V105.hap2","V105.hap1","V126.hap2","V126.hap1","V112.hap2","V112.hap1","DavidiiGrape_hap2","DavidiiGrape_hap1","V007.hap1","V007.hap2","V008.hap2","V008.hap1","V012.hap1","V012.hap2","V100.hap2","V100.hap1","V099.hap1","V099.hap2","V098.hap2","V098.hap1","V117.hap2","V117.hap1","V120.hap2","V120.hap1","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap2","V124.hap1","V125.hap2","V125.hap1","WoollyGrape_hap2","WoollyGrape_hap1","V106.hap1","V106.hap2","V102.hap2","V102.hap1","V079.hap1","V079.hap2","V030.hap2","V030.hap1","V031.hap1","V031.hap2","V005.hap2","V005.hap1","V048.hap1","V048.hap2","V050.hap2","V050.hap1","V051.hap1","V051.hap2","V049.hap2","V049.hap1","V038.hap2","V038.hap1","V043.hap1","V043.hap2","V052.hap2","V052.hap1","V034.hap2","V034.hap1","V023.hap2","V023.hap1","V032.hap2","V032.hap1","V033.hap2","V033.hap1","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap2","V022.hap1","V040.hap1","V040.hap2","V019.hap2","V019.hap1","V041.hap2","V041.hap1","V077.hap1","V077.hap2","V037.hap2","V037.hap1","V096.hap1","V096.hap2","ShineMuscat_hap2","ShineMuscat_hap1","V036.hap2","V036.hap1","V076.hap1","V076.hap2","V055.hap2","V055.hap1","V018.hap1","V018.hap2","V072.hap1","V072.hap2"]
    chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]            
    
    result_list=[]
    for one_sample in used_samples:
        for one_chr in chr_list:
            if one_sample not in dict_sample_chr_type_info: num=0;length=0;pos_str='.'
            else:
                if one_chr not in dict_sample_chr_type_info[one_sample]: num=0;length=0;pos_str='.'
                else: 
                    for one_type,dict_info in dict_sample_chr_type_info[one_sample][one_chr].items():
                        result_list.append([one_sample,one_chr,one_type,dict_info])
                
    with open("./chr2moddotplot3/2_slashin_updown_sum",'w') as f:
        f.write(f"sample\tchromosome\ttype\tnum\tlength\tpure_length\tpos_str_raw\tpos_str_new\tchr_VSat1len\tspeciesgroup\n")                        
    def run_step(one_info):
        one_sample,one_chr,one_type,dict_info=one_info
        num=dict_info['num']
        length=dict_info['length']
        pos_str='|'.join(dict_info['pos_str'])
        pos_set=set()
        for one_pos in dict_info['pos']:
            pos1,pos2,pos3,pos4=one_pos
            kk=pos1
            while kk<=pos3:
                pos_set.add(kk)
                kk+=1
            ##
            if pos2<pos4:kk=pos2;kk_end=pos4
            else: kk=pos4;kk_end=pos2
            while kk<=kk_end:
                pos_set.add(kk)
                kk+=1
        pos_list=list(pos_set)        
        pos_list.sort()
        pure_length=len(pos_list)
        tmp_pos='';pos_result_list=[];pos_result_list_str=''
        for one_pos in pos_list:
            if tmp_pos=='': tmp_pos=one_pos;one_pos_start=tmp_pos
            else: 
                if one_pos-tmp_pos==1:
                    tmp_pos=one_pos
                else:
                    pos_result_list.append(f"{one_pos_start}-{tmp_pos}")
                    tmp_pos=one_pos;one_pos_start=tmp_pos
        pos_result_list.append(f"{one_pos_start}-{tmp_pos}")
        pos_result_list_str="|".join(pos_result_list)
        ####
        sample_chr=f"{one_sample}_{one_chr}"        
        if sample_chr not in dict_samplechr_VSat1len: samplechr_VSat1len=0
        else: samplechr_VSat1len=dict_samplechr_VSat1len[sample_chr]
        speciesgroup=dict_sample_speciesgroup[one_sample] 
        with open("./chr2moddotplot3/2_slashin_updown_sum",'a') as f:
            f.write(f"{one_sample}\t{one_chr}\t{one_type}\t{num}\t{length}\t{pure_length}\t{pos_str}\t{pos_result_list_str}\t{samplechr_VSat1len}\t{speciesgroup}\n")      

        
    with Pool(processes=70) as pool:
         # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step, result_list), start=1):
            # Results can be processed here, e.g., stored or printed
            progress = (i / len(result_list)) * 100
            sys.stdout.write(f"\nProgress: {progress:.2f}% {str(i)} / {str(len(result_list))}#####################################################################\n")
            sys.stdout.flush()                               
if argv1=="stepall" or argv1=="step2_combine" or argv1=="step2b_combine":
    print('Merge and extend sequences')
    print("Calculate the internal length of a single chromosome, i.e., the length of all alignment lines, only count sequences outside VSat1)")
    ## Simply load the VSat1 length for each chromosome
    dict_samplechr_VSat1len={}
    with open("/home/lain/aaa_data/run0/samples_satellite/2_good_regions_interarray",'r') as f:
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            samplechr=eachline_arr[0]+"_"+eachline_arr[1]
            dict_samplechr_VSat1len[samplechr]=int(eachline_arr[4])
    ##
    # Load sample information
    dict_sample_speciesgroup={}
    with open("/home/lain/aaa_data/run0/samples_satellite/sample_info",'r') as f:
        #Baimunage_hap2	Baimunage	 vinifera	Eurasian	Eurasian_vinifera_Baimunage	Table
        for line in f:
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=6:continue
            sample=             eachline_arr[0]
            species_group=      eachline_arr[3]
            dict_sample_speciesgroup[sample]=species_group
    ##  
    dict_sample_chr_type_info={}
    with open("./chr2moddotplot3/2_slashout_updown",'r') as f:   
        #("sample\tchromosome\ttype\tstart_1(x1)\tstart_2(y1)\tend_1(x2)\tend_2(y2)\tseq1_len\tseq2_len\tdot_num\tdotskip_num\tdot_used\tline_slope\tline_R2\tback_accuracy1\tback_accuracy2\tback_distance1\tback_distance2\tdistance\tinclude_length\tclass2VSat1\tunit_length\tunit_number\tclass2VSat1_revised\tmark\n")
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            sample=             eachline_arr[0]
            chromosome=         eachline_arr[1]
            one_type=         eachline_arr[2]
            alignmentline_length=         int(eachline_arr[7])+int(eachline_arr[8])  # Sum of the mapped genomic lengths of the alignment line, i.e., sum of seq1 and seq2
            if sample not in dict_sample_chr_type_info:                  dict_sample_chr_type_info[sample]={}
            if chromosome not in dict_sample_chr_type_info[sample]:      dict_sample_chr_type_info[sample][chromosome]={}
            if one_type not in dict_sample_chr_type_info[sample][chromosome]:
                dict_sample_chr_type_info[sample][chromosome][one_type]={}
                dict_sample_chr_type_info[sample][chromosome][one_type]['length']=0
                dict_sample_chr_type_info[sample][chromosome][one_type]['num']=0 
                dict_sample_chr_type_info[sample][chromosome][one_type]['pos']=[]
                dict_sample_chr_type_info[sample][chromosome][one_type]['pos_str']=[]
            dict_sample_chr_type_info[sample][chromosome][one_type]['length']+=alignmentline_length
            dict_sample_chr_type_info[sample][chromosome][one_type]['num']+=1  
            dict_sample_chr_type_info[sample][chromosome][one_type]['pos'].append([int(eachline_arr[3]),int(eachline_arr[4]),int(eachline_arr[5]),int(eachline_arr[6])])
            dict_sample_chr_type_info[sample][chromosome][one_type]['pos_str'].append(f"{eachline_arr[3]}-{eachline_arr[4]}-{eachline_arr[5]}-{eachline_arr[6]}")
    used_samples= ["Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap2","V059.hap1","V063.hap2","V063.hap1","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap2","Chardonnay_hap1","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap2","V074.hap1","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap2","V062.hap1","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap2","PinotNoir_hap1","PinotNoir2_hap2","PinotNoir2_hap1","V069.hap2","V069.hap1","PN40024","PN40024_hap1","PN40024_hap2","V087.hap2","V087.hap1","V061.hap2","V061.hap1","V065.hap2","V065.hap1","V064.hap2","V064.hap1","ThompsonSeedless_hap1","ThompsonSeedless_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap2","V067.hap1","V058.hap1","V058.hap2","V107.hap2","V107.hap1","V108.hap2","V108.hap1","V081.hap2","V081.hap1","V015.hap2","V015.hap1","V105.hap2","V105.hap1","V126.hap2","V126.hap1","V112.hap2","V112.hap1","DavidiiGrape_hap2","DavidiiGrape_hap1","V007.hap1","V007.hap2","V008.hap2","V008.hap1","V012.hap1","V012.hap2","V100.hap2","V100.hap1","V099.hap1","V099.hap2","V098.hap2","V098.hap1","V117.hap2","V117.hap1","V120.hap2","V120.hap1","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap2","V124.hap1","V125.hap2","V125.hap1","WoollyGrape_hap2","WoollyGrape_hap1","V106.hap1","V106.hap2","V102.hap2","V102.hap1","V079.hap1","V079.hap2","V030.hap2","V030.hap1","V031.hap1","V031.hap2","V005.hap2","V005.hap1","V048.hap1","V048.hap2","V050.hap2","V050.hap1","V051.hap1","V051.hap2","V049.hap2","V049.hap1","V038.hap2","V038.hap1","V043.hap1","V043.hap2","V052.hap2","V052.hap1","V034.hap2","V034.hap1","V023.hap2","V023.hap1","V032.hap2","V032.hap1","V033.hap2","V033.hap1","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap2","V022.hap1","V040.hap1","V040.hap2","V019.hap2","V019.hap1","V041.hap2","V041.hap1","V077.hap1","V077.hap2","V037.hap2","V037.hap1","V096.hap1","V096.hap2","ShineMuscat_hap2","ShineMuscat_hap1","V036.hap2","V036.hap1","V076.hap1","V076.hap2","V055.hap2","V055.hap1","V018.hap1","V018.hap2","V072.hap1","V072.hap2"]
    chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]            
    
    result_list=[]
    for one_sample in used_samples:
        for one_chr in chr_list:
            if one_sample not in dict_sample_chr_type_info: num=0;length=0;pos_str='.'
            else:
                if one_chr not in dict_sample_chr_type_info[one_sample]: num=0;length=0;pos_str='.'
                else: 
                    for one_type,dict_info in dict_sample_chr_type_info[one_sample][one_chr].items():
                        result_list.append([one_sample,one_chr,one_type,dict_info])
                
    with open("./chr2moddotplot3/2_slashout_updown_sum",'w') as f:
        f.write(f"sample\tchromosome\ttype\tnum\tlength\tpure_length\tpos_str_raw\tpos_str_new\tchr_VSat1len\tspeciesgroup\n")                        
    def run_step(one_info):
        one_sample,one_chr,one_type,dict_info=one_info
        num=dict_info['num']
        length=dict_info['length']
        pos_str='|'.join(dict_info['pos_str'])
        pos_set=set()
        for one_pos in dict_info['pos']:
            pos1,pos2,pos3,pos4=one_pos
            kk=pos1
            while kk<=pos3:
                pos_set.add(kk)
                kk+=1
            ##
            if pos2<pos4:kk=pos2;kk_end=pos4
            else: kk=pos4;kk_end=pos2
            while kk<=kk_end:
                pos_set.add(kk)
                kk+=1
        pos_list=list(pos_set)        
        pos_list.sort()
        pure_length=len(pos_list)
        tmp_pos='';pos_result_list=[];pos_result_list_str=''
        for one_pos in pos_list:
            if tmp_pos=='': tmp_pos=one_pos;one_pos_start=tmp_pos
            else: 
                if one_pos-tmp_pos==1:
                    tmp_pos=one_pos
                else:
                    pos_result_list.append(f"{one_pos_start}-{tmp_pos}")
                    tmp_pos=one_pos;one_pos_start=tmp_pos
        pos_result_list.append(f"{one_pos_start}-{tmp_pos}")
        pos_result_list_str="|".join(pos_result_list)
        ####
        sample_chr=f"{one_sample}_{one_chr}"        
        if sample_chr not in dict_samplechr_VSat1len: samplechr_VSat1len=0
        else: samplechr_VSat1len=dict_samplechr_VSat1len[sample_chr]
        speciesgroup=dict_sample_speciesgroup[one_sample] 
        with open("./chr2moddotplot3/2_slashout_updown_sum",'a') as f:
            f.write(f"{one_sample}\t{one_chr}\t{one_type}\t{num}\t{length}\t{pure_length}\t{pos_str}\t{pos_result_list_str}\t{samplechr_VSat1len}\t{speciesgroup}\n")      

        
    with Pool(processes=70) as pool:
         # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step, result_list), start=1):
            # Results can be processed here, e.g., stored or printed
            progress = (i / len(result_list)) * 100
            sys.stdout.write(f"\nProgress: {progress:.2f}% {str(i)} / {str(len(result_list))}#####################################################################\n")
            sys.stdout.flush()

#step3a Mainly analyze large fragment duplications within VSat1 (forward duplications)
if argv1=="stepall" or argv1=="step3" or "step3a" in argv1:
    if argv1=="stepall" or argv1=="step3" or argv1=="step3a" or argv1=="step3a.0":  
        if  os.path.exists('./chr2moddotplot3/3_VSat1in_up_sum')==False:
            subprocess.run(["mkdir ./chr2moddotplot3/3_VSat1in_up_sum"], shell=True) 
        print('Process up class, merging method: if Apos and Bpos intersect, merge AB')    
        ##Store in dictionary
        dict_samplechr_lines={}
        with open ('./chr2moddotplot3/2_slashin_updown','r') as f:
            next(f)
            for line in f:
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                ##
                one_sample=	        eachline_arr[0]
                #if one_sample!='PN40024':continue
                one_chr	=	        eachline_arr[1]
                one_type=	        eachline_arr[2]	
                start_x	=	        int(eachline_arr[3])
                start_y	=	        int(eachline_arr[4])
                end_x	=	        int(eachline_arr[5])
                end_y=	            int(eachline_arr[6])
                if one_type=='down':continue
                samplechr=one_sample+'|||'+one_chr
                if samplechr not in dict_samplechr_lines:
                    dict_samplechr_lines[samplechr]=[]
                dict_samplechr_lines[samplechr].append([start_x,start_y,end_x,end_y]) 
        result_list= []
        for samplechr,lines in dict_samplechr_lines.items():
            result_list.append([samplechr,lines])
        ##
        def run_step(one_result):
            samplechr,lines=one_result
            used_kkk=set()
            lines_num=len(lines)
            one_sample,one_chr=samplechr.split("|||")
            ###Re-plan the idea: put ball 1 into frames 2, 3, 4 sequentially; if satisfied, expand frames 2/3/4
            ###First set all positions as frame sets
            i=0
            dict_index_set={}
            while i<lines_num:
                pos1,pos2,pos3,pos4=lines[i]
                dict_index_set[i]={}
                dict_index_set[i]['pos']=set()
                dict_index_set[i]['index']=set()
                dict_index_set[i]['index'].add(i)
                jj=pos1
                while jj<=pos3:
                    dict_index_set[i]['pos'].add(jj)
                    jj+=1
                jj=pos2
                while jj<=pos4:
                    dict_index_set[i]['pos'].add(jj)
                    jj+=1
                i+=1
            ##Start throwing balls, if thrown then record in bad_list
            bad_list=[]
            i=0 
            while i<lines_num:
                current_pos_set=dict_index_set[i]['pos']
                current_index_set=dict_index_set[i]['index']
                j=i+1
                while j<lines_num:
                    if j in bad_list:j+=1;continue
                    next_set=dict_index_set[j]['pos']
                    if not current_pos_set.isdisjoint(next_set):  # If there is intersection
                        dict_index_set[j]['pos'].update(current_pos_set)
                        bad_list.append(i)
                        ##
                        dict_index_set[j]['index'].update(current_index_set)
                        #dict_index_set[j]['index'].remove(i)
                        break
                    j+=1
                i+=1
            ##Final inventory
            with open ('./chr2moddotplot3/3_VSat1in_up_sum/0_class','a') as f:
                i=0 
                while i<lines_num:
                    if i in bad_list:i+=1;continue
                    cluster_set=dict_index_set[i]['index']
                    for one_index in cluster_set:
                        start_x,start_y,end_x,end_y=lines[one_index]
                        f.write(f"{one_sample}\t{one_chr}\t{i}\t{start_x}\t{start_y}\t{end_x}\t{end_y}\n")
                    i+=1
                
     
                        
        with open ('./chr2moddotplot3/3_VSat1in_up_sum/0_class','w') as f:
                f.write("sample\tchromosome\tclass\tstart_x\tstart_y\tend_x\tend_y\n")        
            
        with Pool(processes=70) as pool:
             # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, result_list), start=1):
                # Process results here, e.g., store or print
                progress = (i / len(result_list)) * 100
                sys.stdout.write(f"\nProgress: {progress:.2f}% {str(i)} / {str(len(result_list))}#####################################################################\n")
                sys.stdout.flush()            
    if argv1=="stepall" or argv1=="step3" or argv1=="step3a" or argv1=="step3a.1":                
        print('Complete information')
        result_list=[]
        with open ('./chr2moddotplot3/3_VSat1in_up_sum/0_class','r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                #one_sample,one_chr,start_x,start_y,end_x,end_y=eachline_arr
                result_list.append(eachline_arr)
        sort_result_list = sorted(result_list, key=lambda x: x[4], reverse=False)
        sort_result_list = sorted(result_list, key=lambda x: (x[0],x[1], x[2], x[3], x[4]), reverse=False)
        print('Load original table')
        dict_index_info={}
        with open ('./chr2moddotplot3/2_slashin_updown','r') as f:
            next(f)
            for line in f:
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                ##
                one_sample=	        eachline_arr[0]
                one_chr	=	        eachline_arr[1]
                one_type=	        eachline_arr[2]
                start_x	=	        int(eachline_arr[3])
                start_y	=	        int(eachline_arr[4])
                end_x	=	        int(eachline_arr[5])
                end_y=	            int(eachline_arr[6])
                dict_index_info[f"{one_sample}|{one_chr}|{start_x}|{start_y}|{end_x}|{end_y}"]=eachline_arr[7:]
        with open ('./chr2moddotplot3/3_VSat1in_up_sum/1_class','w') as f:      
            f.write("sample\tchromosome\ttype\tclass\tstart_x\tstart_y\tend_x\tend_y\tseq1_len\tseq2_len\tdot_num\tdotskip_num\tdot_used\tline_slope\tline_R2\tback_accuracy1\tback_accuracy2\tback_distance1\tback_distance2\n")
            for one_result in sort_result_list:
                index=f"{one_result[0]}|{one_result[1]}|{one_result[3]}|{one_result[4]}|{one_result[5]}|{one_result[6]}"
                newline=f"{one_result[0]}\t{one_result[1]}\tup\t{one_result[2]}\t{one_result[3]}\t{one_result[4]}\t{one_result[5]}\t{one_result[6]}\t"+'\t'.join(dict_index_info[index])+'\n'
                f.write(newline)
    if argv1=="stepall" or argv1=="step3" or argv1=="step3a" or argv1=="step3a.2":
        print('Determine the unit length and extension length from the longest alignment line in each class, calculate the length of a single line and the number of extensions, then estimate the length and number of extensions using n(N-1) method')
        dict_samplechr_VSat1len={}
        with open("/home/lain/aaa_data/run0/samples_satellite/2_good_regions_interarray",'r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                samplechr=eachline_arr[0]+"_"+eachline_arr[1]
                dict_samplechr_VSat1len[samplechr]=int(eachline_arr[4])
        ##
        #Load sample information
        dict_sample_speciesgroup={}
        with open("/home/lain/aaa_data/run0/samples_satellite/sample_info",'r') as f:
            #Baimunage_hap2	Baimunage	 vinifera	Eurasian	Eurasian_vinifera_Baimunage	Table
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=6:continue
                sample=             eachline_arr[0]
                species_group=      eachline_arr[3]
                dict_sample_speciesgroup[sample]=species_group     
        ##########        
        dict_sample_chr_type_info={}
        with open("./chr2moddotplot3/3_VSat1in_up_sum/1_class",'r') as f:   
            #("sample\tchromosome\ttype\tstart_1(x1)\tstart_2(y1)\tend_1(x2)\tend_2(y2)\tseq1_len\tseq2_len\tdot_num\tdotskip_num\tdot_used\tline_slope\tline_R2\tback_accuracy1\tback_accuracy2\tback_distance1\tback_distance2\tdistance\tinclude_length\tclass2VSat1\tunit_length\tunit_number\tclass2VSat1_revised\tmark\n")
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                sample=             eachline_arr[0]
                chromosome=         eachline_arr[1]
                one_class=         eachline_arr[3]
                alignmentline_length=         int(eachline_arr[8])+int(eachline_arr[9])  #Sum of genomic mapping lengths corresponding to alignment lines, i.e., sum of seq1 and seq2
                if sample not in dict_sample_chr_type_info:                  dict_sample_chr_type_info[sample]={}
                if chromosome not in dict_sample_chr_type_info[sample]:      dict_sample_chr_type_info[sample][chromosome]={}
                if one_class not in dict_sample_chr_type_info[sample][chromosome]:
                    dict_sample_chr_type_info[sample][chromosome][one_class]={}
                    dict_sample_chr_type_info[sample][chromosome][one_class]['length']=0
                    dict_sample_chr_type_info[sample][chromosome][one_class]['num']=0 
                    dict_sample_chr_type_info[sample][chromosome][one_class]['pos']=[]
                    dict_sample_chr_type_info[sample][chromosome][one_class]['pos_str']=[]
                dict_sample_chr_type_info[sample][chromosome][one_class]['length']+=alignmentline_length
                dict_sample_chr_type_info[sample][chromosome][one_class]['num']+=1  
                dict_sample_chr_type_info[sample][chromosome][one_class]['pos'].append([int(eachline_arr[4]),int(eachline_arr[5]),int(eachline_arr[6]),int(eachline_arr[7])])
                dict_sample_chr_type_info[sample][chromosome][one_class]['pos_str'].append(f"{eachline_arr[4]}-{eachline_arr[5]}-{eachline_arr[6]}-{eachline_arr[7]}")
        used_samples= ["Baimunage_hap1","Baimunage_hap2","BlackMonukka_hap1","BlackMonukka_hap2","V059.hap2","V059.hap1","V063.hap2","V063.hap1","V066.hap1","V066.hap2","V070.hap1","V070.hap2","Chardonnay_hap2","Chardonnay_hap1","VHP-T2T.hap1","VHP-T2T.hap2","V091.hap1","V091.hap2","Hongmunage_hap1","Hongmunage_hap2","V074.hap2","V074.hap1","ManicureFinger_hap1","ManicureFinger_hap2","V092.hap1","V092.hap2","V062.hap2","V062.hap1","MuscatHamburg_hap1","MuscatHamburg_hap2","V088.hap1","V088.hap2","PinotNoir_hap2","PinotNoir_hap1","PinotNoir2_hap2","PinotNoir2_hap1","V069.hap2","V069.hap1","PN40024","PN40024_hap1","PN40024_hap2","V087.hap2","V087.hap1","V061.hap2","V061.hap1","V065.hap2","V065.hap1","V064.hap2","V064.hap1","ThompsonSeedless_hap1","ThompsonSeedless_hap2","V093.hap1","V093.hap2","V060.hap1","V060.hap2","V067.hap2","V067.hap1","V058.hap1","V058.hap2","V107.hap2","V107.hap1","V108.hap2","V108.hap1","V081.hap2","V081.hap1","V015.hap2","V015.hap1","V105.hap2","V105.hap1","V126.hap2","V126.hap1","V112.hap2","V112.hap1","DavidiiGrape_hap2","DavidiiGrape_hap1","V007.hap1","V007.hap2","V008.hap2","V008.hap1","V012.hap1","V012.hap2","V100.hap2","V100.hap1","V099.hap1","V099.hap2","V098.hap2","V098.hap1","V117.hap2","V117.hap1","V120.hap2","V120.hap1","PiasezkiiGrape_hap1","PiasezkiiGrape_hap2","V123.hap1","V123.hap2","V124.hap2","V124.hap1","V125.hap2","V125.hap1","WoollyGrape_hap2","WoollyGrape_hap1","V106.hap1","V106.hap2","V102.hap2","V102.hap1","V079.hap1","V079.hap2","V030.hap2","V030.hap1","V031.hap1","V031.hap2","V005.hap2","V005.hap1","V048.hap1","V048.hap2","V050.hap2","V050.hap1","V051.hap1","V051.hap2","V049.hap2","V049.hap1","V038.hap2","V038.hap1","V043.hap1","V043.hap2","V052.hap2","V052.hap1","V034.hap2","V034.hap1","V023.hap2","V023.hap1","V032.hap2","V032.hap1","V033.hap2","V033.hap1","V053.hap1","V053.hap2","V020.hap1","V020.hap2","V022.hap2","V022.hap1","V040.hap1","V040.hap2","V019.hap2","V019.hap1","V041.hap2","V041.hap1","V077.hap1","V077.hap2","V037.hap2","V037.hap1","V096.hap1","V096.hap2","ShineMuscat_hap2","ShineMuscat_hap1","V036.hap2","V036.hap1","V076.hap1","V076.hap2","V055.hap2","V055.hap1","V018.hap1","V018.hap2","V072.hap1","V072.hap2"]
        chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]            
        
        result_list=[]
        for one_sample in used_samples:
            for one_chr in chr_list:
                if one_sample not in dict_sample_chr_type_info: num=0;length=0;pos_str='.'
                else:
                    if one_chr not in dict_sample_chr_type_info[one_sample]: num=0;length=0;pos_str='.'
                    else: 
                        for one_class,dict_info in dict_sample_chr_type_info[one_sample][one_chr].items():
                            result_list.append([one_sample,one_chr,one_class,dict_info])
                    
        with open("./chr2moddotplot3/3_VSat1in_up_sum/2_estimate",'w') as f:
            f.write(f"sample\tchromosome\tclass\tnum\tlength\tpure_length\tpos_str_new\tchr_VSat1len\tspeciesgroup\tunit_number\tunit_length\tmaxline_pos\tmax_linelen\tmaxline_unit_length\tmaxline_unit_number\tCascade_model_state\tunit_num_type\n")    #pos_str_raw\t                     
        def run_step(one_info):
            one_sample,one_chr,one_class,dict_info=one_info
            num=dict_info['num']
            length=dict_info['length']
            pos_str='|'.join(dict_info['pos_str'])
            #####Find the position of the longest alignment line
            max_pos1,max_pos2,max_pos3,max_pos4=dict_info['pos'][0]
            max_linelen=abs(max_pos3-max_pos1)+abs(max_pos4-max_pos2)
            for one_pos in dict_info['pos']:
                tmp_pos1,tmp_pos2,tmp_pos3,tmp_pos4=one_pos
                tmp_delta=abs(tmp_pos3-tmp_pos1)+abs(tmp_pos4-tmp_pos2)
                if tmp_delta>max_linelen:
                    max_pos1,max_pos2,max_pos3,max_pos4=one_pos
                    max_linelen=tmp_delta
            maxline_pos=f"{max_pos1}-{max_pos2}-{max_pos3}-{max_pos4}"        
            maxline_unit_length=   round(((max_pos2-max_pos1)+(max_pos4-max_pos3))/2,0)
            maxline_unit_number=round(max_linelen/maxline_unit_length,3)
            #####
            pos_set=set()
            for one_pos in dict_info['pos']:
                pos1,pos2,pos3,pos4=one_pos
                kk=pos1
                while kk<=pos3:
                    pos_set.add(kk)
                    kk+=1
                ##
                if pos2<pos4:kk=pos2;kk_end=pos4
                else: kk=pos4;kk_end=pos2
                while kk<=kk_end:
                    pos_set.add(kk)
                    kk+=1
            pos_list=list(pos_set)        
            pos_list.sort()
            pure_length=len(pos_list)
            tmp_pos='';pos_result_list=[];pos_result_list_str=''
            for one_pos in pos_list:
                if tmp_pos=='': tmp_pos=one_pos;one_pos_start=tmp_pos
                else: 
                    if one_pos-tmp_pos==1:
                        tmp_pos=one_pos
                    else:
                        pos_result_list.append(f"{one_pos_start}-{tmp_pos}")
                        tmp_pos=one_pos;one_pos_start=tmp_pos
            pos_result_list.append(f"{one_pos_start}-{tmp_pos}")
            pos_result_list_str="|".join(pos_result_list)
            ####
            sample_chr=f"{one_sample}_{one_chr}"        
            if sample_chr not in dict_samplechr_VSat1len: samplechr_VSat1len=0
            else: samplechr_VSat1len=dict_samplechr_VSat1len[sample_chr]
            speciesgroup=dict_sample_speciesgroup[one_sample] 
            #
            unit_num=           length/pure_length+1
            unit_length=        pure_length/unit_num
            #
            if abs((unit_length-maxline_unit_length)/maxline_unit_number)<0.3 or abs((unit_length-maxline_unit_length)/unit_length)<0.3 :normal_or_complex='normal'
            else :
                normal_or_complex='complex'
            ##
            if maxline_unit_number >1.8 and maxline_unit_number <2.2 and unit_num>1.8 and unit_num <2.2:
                unit_num_type='copy1'
            else:
                unit_num_type='other'
            with open("./chr2moddotplot3/3_VSat1in_up_sum/2_estimate",'a') as f:
                f.write(f"{one_sample}\t{one_chr}\t{one_class}\t{num}\t{length}\t{pure_length}\t{pos_result_list_str}\t{samplechr_VSat1len}\t{speciesgroup}\t{unit_num}\t{unit_length}\t{maxline_pos}\t{max_linelen}\t{maxline_unit_length}\t{maxline_unit_number}\t{normal_or_complex}\t{unit_num_type}\n")      ##pos_str
        



            
        with Pool(processes=70) as pool:
             # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, result_list), start=1):
                # Process results here, e.g., store or print
                progress = (i / len(result_list)) * 100
                sys.stdout.write(f"\nProgress: {progress:.2f}% {str(i)} / {str(len(result_list))}#####################################################################\n")
                sys.stdout.flush()                   
    if argv1=="stepall" or argv1=="step3" or argv1=="step3a" or argv1=="step3a.3p":
        if  os.path.exists('./chr2moddotplot3/3_VSat1in_up_sum/3_plot')==False:
                subprocess.run(["mkdir ./chr2moddotplot3/3_VSat1in_up_sum/3_plot"], shell=True) 
        subprocess.run(["cp ./chr2moddotplot3/3_VSat1in_up_sum/2_estimate ./chr2moddotplot3/3_VSat1in_up_sum/3_plot/inputfile"], shell=True)     
        R_txt=r'''
library(ggplot2)
library(dplyr)

print("")
{
  input_file=read.table('inputfile', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
print(colnames(input_file))
 input_file$unit_number <- as.numeric(input_file$unit_number)
# Filter data
filtered_input <- input_file %>% 
  filter(maxline_unit_number >1.8 ) %>%
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
          chromosome == "Chr19" ~ "Chr19"
          #TRUE ~ "Other"
        ),
      )      %>%
        mutate(
          category2 = case_when(
            speciesgroup == "Eurasian" ~ "Eurasian",
            speciesgroup == "East_Asia" ~ "East_Asia",
            speciesgroup == "America" ~ "America"#,
            #speciesgroup == "Routundifolia" ~ "Routundifolia"#,
            #TRUE ~ "Other"
          ),
        ) %>%
        mutate(
          category3 = case_when(
            num == "1" ~ "1",
            TRUE ~ "16"
          ),
        )%>%
        mutate(
          category4 = case_when(
            unit_number <2.1 ~ 1,
            TRUE ~ 16
          ),
        )%>%
        mutate(
          category5 = case_when(
            Cascade_model_state =='normal' ~ 16,
            TRUE ~ 1
          ),
        )%>%
        mutate(
          category6 = case_when(
            unit_num_type=='copy1'~ 1,
            TRUE ~ 16
          ),
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
                                      "Other" = "#CCEBC5",
              "Eurasian" = "#ff3399",
                "East_Asia" = "#0066ff",
                "America" = "#33cc33"#,
                #"Routundifolia" = "black"
)


# Create plot object
p <- ggplot() +
    geom_point(
        data = filtered_input,
        aes(x = chr_VSat1len	, y = length,color=category2)
    ) +coord_equal() +

theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank()
    ) #+
    #scale_color_manual(values = color_values, drop = FALSE)
    
 filtered_input$category3 <- as.numeric(filtered_input$category3)
  
  # Save as PDF
  pdf(file = paste0('up_chr_length', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
  print(p)
  dev.off()


  
  # Create plot object
p <- ggplot() +
    geom_boxplot(
        data = filtered_input,
        aes(x = category2	, y = length,color=category2),
        width = 0.8,          # Reduce boxplot width
        outlier.shape = NA   # Hide boxplot outliers (since scatter plot already shows them)
    ) +
    geom_jitter(
        data = filtered_input,
        aes(x = category2	, y = length,color=category2,shape=category6,size=max_linelen/500000),
        width = 0.2,         # Jitter amplitude in x-axis direction
        #size = 1.5,          # Point size
        alpha = 0.6          # Transparency
        #shape = 16,       # Solid circle (no border)
        #stroke = 0
    ) +scale_shape_identity()+  scale_size_identity()+


theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_x_discrete(limits = c("Eurasian", "East_Asia", "America")) +
    scale_color_manual(values = color_values, drop = FALSE)
 
 
  
  # Save as PDF
  pdf(file = paste0('up_chr_length3', ".pdf"), width = 6/ 2.54, height = 10 / 2.54)
  print(p)
  dev.off() 
  
 
}
    '''
        with open('./chr2moddotplot3/3_VSat1in_up_sum/3_plot/up_chr_length.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2moddotplot3//3_VSat1in_up_sum/3_plot/'
        os.chdir(new_directory)
        subprocess.run(['Rscript up_chr_length.R'], shell=True)    
        os.chdir('../../../')       
    if argv1=="stepall" or argv1=="step3" or argv1=="step3a" or argv1=="step3a.5":
        print('For each chromosome, identify only the longest alignment line, then count the lengths of its related alignment lines. If there are none, it is single; if there are, all are counted as point sizes')
        dict_index_info={}
        i=0
        with open ('./chr2moddotplot3/3_VSat1in_up_sum/2_estimate','r') as f:
            for line in f:
                eachline=line.strip()
                i+=1
                if i==1: head=eachline;continue
                eachline_arr=eachline.split('\t')
                maxline_unit_length=int(eachline_arr[12])
                sample_chr=eachline_arr[0]+"___"+eachline_arr[1]
                if sample_chr not in dict_index_info:
                    dict_index_info[sample_chr]={}
                    dict_index_info[sample_chr]['maxline_len']=0
                    dict_index_info[sample_chr]['class_info']=''
                if maxline_unit_length>dict_index_info[sample_chr]['maxline_len']:
                    dict_index_info[sample_chr]['maxline_len']=maxline_unit_length
                    dict_index_info[sample_chr]['class_info']=eachline
        with open ('./chr2moddotplot3/3_VSat1in_up_sum/5_estimate_chrmaxonly','w') as f:
            f.write(head+'\n')
            for sample_chr,info in dict_index_info.items():
                line=info['class_info']
                f.write(line+'\n')
    if argv1=="stepall" or argv1=="step3" or argv1=="step3a" or argv1=="step3a.6p1":
        if  os.path.exists('./chr2moddotplot3/3_VSat1in_up_sum/6_plot')==False:
                subprocess.run(["mkdir ./chr2moddotplot3/3_VSat1in_up_sum/6_plot"], shell=True) 
        subprocess.run(["cp ./chr2moddotplot3/3_VSat1in_up_sum/5_estimate_chrmaxonly ./chr2moddotplot3/3_VSat1in_up_sum/6_plot/inputfile"], shell=True)     
        R_txt=r'''
library(ggplot2)
library(dplyr)

print("")
{
  input_file=read.table('inputfile', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
print(colnames(input_file))
 input_file$unit_number <- as.numeric(input_file$unit_number)
# Filter data
filtered_input <- input_file %>% 
  filter(maxline_unit_number >1.8 ) %>%
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
          chromosome == "Chr19" ~ "Chr19"
          #TRUE ~ "Other"
        ),
      )      %>%
        mutate(
          category2 = case_when(
            speciesgroup == "Eurasian" ~ "Eurasian",
            speciesgroup == "East_Asia" ~ "East_Asia",
            speciesgroup == "America" ~ "America"#,
            #speciesgroup == "Routundifolia" ~ "Routundifolia"#,
            #TRUE ~ "Other"
          ),
        ) %>%
        mutate(
          category3 = case_when(
            num == "1" ~ "1",
            TRUE ~ "16"
          ),
        )%>%
        mutate(
          category4 = case_when(
            unit_number <2.1 ~ 1,
            TRUE ~ 16
          ),
        )%>%
        mutate(
          category5 = case_when(
            Cascade_model_state =='normal' ~ 16,
            TRUE ~ 1
          ),
        )%>%
        mutate(
          category6 = case_when(
            unit_num_type=='copy1'~ 1,
            TRUE ~ 16
          ),
        ) %>%
        mutate(
          category7 = case_when(
            maxline_unit_number <2.5 ~ 1,
            TRUE ~ 16
          ),
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
                                      "Other" = "#CCEBC5",
              "Eurasian" = "#ff3399",
                "East_Asia" = "#0066ff",
                "America" = "#33cc33"#,
                #"Routundifolia" = "black"
)


# Create plot object
p <- ggplot() +
    geom_point(
        data = filtered_input,
        aes(x = chr_VSat1len	, y = length,color=category2)
    ) +coord_equal() +

theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank()
    ) #+
    #scale_color_manual(values = color_values, drop = FALSE)
    
 filtered_input$category3 <- as.numeric(filtered_input$category3)
  
  # Save as PDF
  pdf(file = paste0('up_chr_length', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
  print(p)
  dev.off()


  
  # Create plot object
p <- ggplot() +
    geom_boxplot(
        data = filtered_input,
        aes(x = category2	, y = max_linelen,color=category2),
        width = 0.8,          # Reduce boxplot width
        outlier.shape = NA   # Hide boxplot outliers (since scatter plot already shows them)
    ) +
    geom_jitter(
        data = filtered_input,
        aes(x = category2	, y = max_linelen,color=category2,shape=category7,size=length/1000000),
        width = 0.2,         # Jitter amplitude in x-axis direction
        #size = 1.5,          # Point size
        alpha = 0.6          # Transparency
        #shape = 16,       # Solid circle (no border)
        #stroke = 0
    ) +scale_shape_identity()+  scale_size_identity()+


theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_x_discrete(limits = c("Eurasian", "East_Asia", "America")) +
    scale_color_manual(values = color_values, drop = FALSE)
 
 
  
  # Save as PDF
  pdf(file = paste0('up_chr_length3', ".pdf"), width = 6/ 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
}
    '''
        with open('./chr2moddotplot3/3_VSat1in_up_sum/6_plot/up_chr_length.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2moddotplot3//3_VSat1in_up_sum/6_plot/'
        os.chdir(new_directory)
        subprocess.run(['Rscript up_chr_length.R'], shell=True)    
        os.chdir('../../../')                   
    if argv1=="stepall" or argv1=="step3" or argv1=="step3a" or argv1=="step3a.6p2":
        if  os.path.exists('./chr2moddotplot3/3_VSat1in_up_sum/6_plot')==False:
                subprocess.run(["mkdir ./chr2moddotplot3/3_VSat1in_up_sum/6_plot"], shell=True) 
        subprocess.run(["cp ./chr2moddotplot3/3_VSat1in_up_sum/5_estimate_chrmaxonly ./chr2moddotplot3/3_VSat1in_up_sum/6_plot/inputfile"], shell=True)     
        R_txt=r'''
library(ggplot2)
library(dplyr)

print("")
{
  input_file=read.table('inputfile', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
print(colnames(input_file))
 input_file$unit_number <- as.numeric(input_file$unit_number)
# Filter data
filtered_input <- input_file %>% 
  filter(maxline_unit_number >1.8 ) %>%
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
          chromosome == "Chr19" ~ "Chr19"
          #TRUE ~ "Other"
        ),
      )      %>%
    mutate(
          category2 = case_when(
            speciesgroup == "Eurasian" ~ "Eurasian",
            speciesgroup == "East_Asia" ~ "East_Asia",
            speciesgroup == "America" ~ "America"#,
            #speciesgroup == "Routundifolia" ~ "Routundifolia"#,
            #TRUE ~ "Other"
          ),
        ) %>%
    mutate(        
        category3 = case_when(
          Cascade_model_state=="complex" ~ "complex",
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
          chromosome == "Chr19" ~ "Chr19"
          #TRUE ~ "Other"
        ),
      )   
      
# Define color values
color_values <- c(
                                      "complex"="grey",
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
                                      "Other" = "#CCEBC5",
              "Eurasian" = "#ff3399",
                "East_Asia" = "#0066ff",
                "America" = "#33cc33"#,
                #"Routundifolia" = "black"
)
  
 # Create plot object
    B1 <- filtered_input %>% filter(maxline_unit_number < 2.5)
    B2 <- filtered_input %>% filter(maxline_unit_number >=2.5)   

    p <- ggplot() +
        # Boxplot (in Chr1-Chr19 order)
        geom_boxplot(
            data = filtered_input,
            aes(x =category, y = max_linelen),
            outlier.shape = NA,  # Hide boxplot outliers (since scatter plot will show them)
            width = 0.5,         # Adjust boxplot width
            alpha = 0.5          # Semi-transparent effect
        ) +
        # Scatter plot (jitter to avoid overlap)
        geom_jitter(
            data = B2,
            aes(x = category, y = max_linelen, color = category3,size=length/1000000),
            width = 0.2,         # Adjust scatter horizontal jitter range
            shape=16,
            #size = 2,            # Adjust scatter point size
            alpha = 0.9          # Semi-transparent effect
        ) +   
        geom_jitter(
            data = B1,
            aes(x = category, y = max_linelen, color = category3,size=length/1000000),
            width = 0.2,         # Adjust scatter horizontal jitter range
            shape=1,
            #size = 2,            # Adjust scatter point size
            alpha = 0.9          # Semi-transparent effect
        ) +        
        # Set x-axis label order
        scale_x_discrete(limits = paste0("Chr", 1:19)) +
    theme_classic() +         
        theme(
          #axis.ticks.y = element_blank(),
          #axis.text.y = element_blank(),
          #legend.position = "none",
          #axis.text.x = element_blank()
        ) +
        scale_color_manual(values = color_values, drop = FALSE) #+  scale_size_identity()

      
      # Save as PDF
      pdf(file = paste0('chr_bar', ".pdf"), width = 24/ 2.54, height = 13 / 2.54)
      print(p)
      dev.off()
}
    '''
        with open('./chr2moddotplot3/3_VSat1in_up_sum/6_plot/up_chr_length.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2moddotplot3//3_VSat1in_up_sum/6_plot/'
        os.chdir(new_directory)
        subprocess.run(['Rscript up_chr_length.R'], shell=True)    
        os.chdir('../../../')                  

if argv1=="stepall" or "step3_rough" in argv1:
    if  os.path.exists('./chr2moddotplot3/3_VSat1in_up_rough')==False:
        subprocess.run(["mkdir ./chr2moddotplot3/3_VSat1in_up_rough"], shell=True) 
    #########Analyze the longest alignment line
    dict_samplechr_len={}
    with open('./chr2moddotplot3/2_slashin_updown','r') as f:
        next(f)
        for line in f:
            eachline=line.strip()
            eachline_arr=eachline.split('\t')   
            if eachline_arr[2]!='up':continue
            seq_len=int(eachline_arr[7])+int(eachline_arr[8])
            samplechr=eachline_arr[0]+'___'+eachline_arr[1]
            if samplechr not in dict_samplechr_len:
                dict_samplechr_len[samplechr]={}
                dict_samplechr_len[samplechr]['maxline_unit_length']=0
                dict_samplechr_len[samplechr]['max_linelen']=0
            if seq_len>dict_samplechr_len[samplechr]['maxline_unit_length']:
                dict_samplechr_len[samplechr]['maxline_unit_length']=seq_len
                dict_samplechr_len[samplechr]['max_linelen']=int(eachline_arr[4])-int(eachline_arr[3])
    #########
    if argv1=="stepall" or argv1=="step3_rough" or argv1=="step3.0_rough":
        with open ('./chr2moddotplot3/3_VSat1in_up_rough/slashin_updown_sum','w') as f2:
            f2.write("sample\tchromosome\ttype\tnum\tlength\tpure_length\tpos_str_raw\tpos_str_new\tchr_VSat1len\tspeciesgroup\tunit_number\tunit_length\tmaxline_len\tmaxline_unitlen\tmaxline_unitnum\n")
            with open ('./chr2moddotplot3/2_slashin_updown_sum','r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    onetype=            eachline_arr[2]
                    if onetype=='down':continue
                    samplechr=eachline_arr[0]+'___'+eachline_arr[1]
                    length=             int(eachline_arr[4])
                    pure_length=         int(eachline_arr[5])
                    unit_num=           length/pure_length+1
                    unit_length=        pure_length/unit_num
                    #
                    maxline_unit_length=    dict_samplechr_len[samplechr]['maxline_unit_length']
                    max_linelen=            dict_samplechr_len[samplechr]['max_linelen']
                    maxline_unitnum=        round(maxline_unit_length/max_linelen,3)
                    
                    f2.write(f"{eachline}\t{unit_num}\t{unit_length}\t{max_linelen}\t{maxline_unit_length}\t{maxline_unitnum}\n")
    if argv1=="stepall" or argv1=="step3_rough" or argv1=="step3.1_rough":
        print('Statistics of alignment line lengths per chromosome')
        if  os.path.exists('./chr2moddotplot3/3_VSat1in_up_rough/1_Chr_len')==False:
            subprocess.run(["mkdir ./chr2moddotplot3/3_VSat1in_up_rough/1_Chr_len"], shell=True) 
        subprocess.run(["cp ./chr2moddotplot3/3_VSat1in_up_rough/slashin_updown_sum ./chr2moddotplot3/3_VSat1in_up_rough/1_Chr_len/slashin_updown_sum"], shell=True)     
        R_txt=r'''
    library(ggplot2)
    library(dplyr)
    
    print("")
    {
      input_file=read.table('slashin_updown_sum', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
    print(colnames(input_file))
    up_ratio <- mean(input_file$unit_number>= 2.5, na.rm = TRUE)  # Calculate proportion
    up_percentage <- round(up_ratio * 100, 2.5)
    # Filter data
    filtered_input <- input_file %>% 
      filter(type == "up" ) %>%
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
              chromosome == "Chr19" ~ "Chr19"
              #TRUE ~ "Other"
            ),
          )
        B1 <- filtered_input %>% filter(unit_number < 2.5)
        B2 <- filtered_input %>% filter(unit_number >=2.5)      
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
                                          #'up'='#006699',
                                          #'down'='red'
    )
    
    #filtered_input$category <- factor(filtered_input$category, levels = paste0("Chr", 1:19))
    
     
    # Create plot object
    p <- ggplot() +
        # Boxplot (in Chr1-Chr19 order)
        geom_boxplot(
            data = filtered_input,
            aes(x =category, y = unit_length),
            outlier.shape = NA,  # Hide boxplot outliers (since scatter plot will show them)
            width = 0.5,         # Adjust boxplot width
            alpha = 0.5          # Semi-transparent effect
        ) +
        # Scatter plot (jitter to avoid overlap)
        geom_jitter(
            data = B1,
            aes(x = category, y = unit_length, color = category,size=pure_length),
            width = 0.2,         # Adjust scatter horizontal jitter range
            shape=1,
            #size = 2,            # Adjust scatter point size
            alpha = 0.9          # Semi-transparent effect
        ) +
        geom_jitter(
            data = B2,
            aes(x = category, y = unit_length, color = category,size=pure_length),
            width = 0.2,         # Adjust scatter horizontal jitter range
            shape=16,
            #size = 2,            # Adjust scatter point size
            alpha = 0.9          # Semi-transparent effect
        ) +        
        # Set x-axis label order
        scale_x_discrete(limits = paste0("Chr", 1:19)) +
    theme_classic() +         
        theme(
          #axis.ticks.y = element_blank(),
          #axis.text.y = element_blank(),
          #legend.position = "none",
          #axis.text.x = element_blank()
        ) +
        scale_color_manual(values = color_values, drop = FALSE)+
        labs(
            title = paste("Proportion of 'up':", up_percentage, "%"),
            #x = "Type",
            #y = "Count"
          )
     
      
      # Save as PDF
      pdf(file = paste0('slashin_updown_sum', ".pdf"), width = 24/ 2.54, height = 13 / 2.54)
      print(p)
      dev.off()
    
    }
        '''
        with open('./chr2moddotplot3/3_VSat1in_up_rough/1_Chr_len/slashin_updown.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2moddotplot3/3_VSat1in_up_rough/1_Chr_len/'
        os.chdir(new_directory)
        subprocess.run(['Rscript slashin_updown.R'], shell=True)    
        os.chdir('../../../')     
    if argv1=="stepall" or argv1=="step3_rough" or argv1=="step3.2_rough":
        if  os.path.exists('./chr2moddotplot3/3_VSat1in_up_rough/2_species_len')==False:
            subprocess.run(["mkdir ./chr2moddotplot3/3_VSat1in_up_rough/2_species_len"], shell=True) 
        subprocess.run(["cp ./chr2moddotplot3/3_VSat1in_up_rough/slashin_updown_sum ./chr2moddotplot3/3_VSat1in_up_rough/2_species_len/slashin_updown_sum"], shell=True)     
        R_txt=r'''
library(ggplot2)
library(dplyr)

print("")
{
  input_file=read.table('slashin_updown_sum', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
print(colnames(input_file))
 input_file$unit_number <- as.numeric(input_file$unit_number)
# Filter data
filtered_input <- input_file %>% 
  filter(type == "up" ) %>%
  filter(maxline_unitnum>1.8 & unit_number >1.8 ) %>%
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
          chromosome == "Chr19" ~ "Chr19"
          #TRUE ~ "Other"
        ),
      )      %>%
        mutate(
          category2 = case_when(
            speciesgroup == "Eurasian" ~ "Eurasian",
            speciesgroup == "East_Asia" ~ "East_Asia",
            speciesgroup == "America" ~ "America"#,
            #speciesgroup == "Routundifolia" ~ "Routundifolia"#,
            #TRUE ~ "Other"
          ),
        ) %>%
        mutate(
          category3 = case_when(
            num == "1" ~ "1",
            TRUE ~ "16"
          ),
        )%>%
        mutate(
          category4 = case_when(
            unit_number <2.1 ~ 1,
            TRUE ~ 16
          ),
        )%>%
        mutate(
          category5 = case_when(
            maxline_unitnum <2.2 & unit_number<2.2 ~ 1,
            TRUE ~ 16
          ),
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
                                      "Other" = "#CCEBC5",
              "Eurasian" = "#ff3399",
                "East_Asia" = "#0066ff",
                "America" = "#33cc33"#,
                #"Routundifolia" = "black"
)


# Create plot object
p <- ggplot() +
    geom_point(
        data = filtered_input,
        aes(x = chr_VSat1len, y = length, color = category2)
    ) + coord_equal() +

theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank()
    ) #+
    #scale_color_manual(values = color_values, drop = FALSE)
    
 filtered_input$category3 <- as.numeric(filtered_input$category3)
  
  # Save as PDF
  pdf(file = paste0('up_chr_length', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
  print(p)
  dev.off()

# Create plot object
p <- ggplot() +
    geom_boxplot(
        data = filtered_input,
        aes(x = category2, y = pure_length, color = category2),
        width = 0.8,          # Reduce boxplot width
        outlier.shape = NA,   # Hide boxplot outliers (since scatter plot shows them)
    ) +
    geom_jitter(
        data = filtered_input,
        aes(x = category2, y = pure_length, color = category2, shape = category3),
        width = 0.2,         # Jitter range on x-axis
        size = 1.5,          # Point size
        alpha = 0.6,         # Transparency
        #shape = 16,       # Solid circle (no border)
        #stroke = 0
    ) + scale_shape_identity() +


theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_x_discrete(limits = c("Eurasian", "East_Asia", "America")) +
    scale_color_manual(values = color_values, drop = FALSE)
 
 
  
  # Save as PDF
  pdf(file = paste0('up_chr_length2', ".pdf"), width = 6/ 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
  
  
  # Create plot object
p <- ggplot() +
    geom_boxplot(
        data = filtered_input,
        aes(x = category2, y = unit_length, color = category2),
        width = 0.8,          # Reduce boxplot width
        outlier.shape = NA   # Hide boxplot outliers (since scatter plot shows them)
    ) +
    geom_jitter(
        data = filtered_input,
        aes(x = category2, y = unit_length, color = category2, shape = category3, size = length/1000000),
        width = 0.2,         # Jitter range on x-axis
        #size = 1.5,          # Point size
        alpha = 0.6          # Transparency
        #shape = 16,       # Solid circle (no border)
        #stroke = 0
    ) + scale_shape_identity() + scale_size_identity() +


theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_x_discrete(limits = c("Eurasian", "East_Asia", "America")) +
    scale_color_manual(values = color_values, drop = FALSE)
 
 
  
  # Save as PDF
  pdf(file = paste0('up_chr_length3', ".pdf"), width = 6/ 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
  
  
    # Create plot object
p <- ggplot() +
    geom_boxplot(
        data = filtered_input,
        aes(x = category2, y = length/pure_length, color = category2),
        width = 0.8,          # Reduce boxplot width
        outlier.shape = NA   # Hide boxplot outliers (since scatter plot shows them)
    ) +
    geom_jitter(
        data = filtered_input,
        aes(x = category2, y = length/pure_length, color = category2, shape = category4),  #, size = length/1000000
        width = 0.2,         # Jitter range on x-axis
        #size = 1.5,          # Point size
        alpha = 0.6          # Transparency
        #shape = 16,       # Solid circle (no border)
        #stroke = 0
    ) + scale_shape_identity() +  #scale_size_identity() +


theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_x_discrete(limits = c("Eurasian", "East_Asia", "America")) +
    scale_color_manual(values = color_values, drop = FALSE)
 
 
  
  # Save as PDF
  pdf(file = paste0('up_chr_length4', ".pdf"), width = 6/ 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
  
      # Create plot object
p <- ggplot() +
    geom_boxplot(
        data = filtered_input,
        aes(x = category2, y = length, color = category2),
        width = 0.8,          # Reduce boxplot width
        outlier.shape = NA   # Hide boxplot outliers (since scatter plot shows them)
    ) +
    geom_jitter(
        data = filtered_input,
        aes(x = category2, y = length, color = category2, shape = category4),  #, size = length/1000000
        width = 0.2,         # Jitter range on x-axis
        #size = 1.5,          # Point size
        alpha = 0.6          # Transparency
        #shape = 16,       # Solid circle (no border)
        #stroke = 0
    ) + scale_shape_identity() +  #scale_size_identity() +


theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_x_discrete(limits = c("Eurasian", "East_Asia", "America")) +
    scale_color_manual(values = color_values, drop = FALSE)
 
# Create plot object
p <- ggplot() +
    geom_boxplot(
        data = filtered_input,
        aes(x = category2, y = pure_length, color = category2),
        width = 0.8,          # Reduce boxplot width
        outlier.shape = NA   # Hide boxplot outliers (since scatter plot shows them)
    ) +
    geom_jitter(
        data = filtered_input,
        aes(x = category2, y = pure_length, color = category2, shape = category5, size = length/1000000),  #, size = length/1000000
        width = 0.2,         # Jitter range on x-axis
        #size = 1.5,          # Point size
        alpha = 0.6          # Transparency
        #shape = 16,       # Solid circle (no border)
        #stroke = 0
    ) + scale_shape_identity() + scale_size_identity() +


theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank()
    ) +
    scale_x_discrete(limits = c("Eurasian", "East_Asia", "America")) +
    scale_color_manual(values = color_values, drop = FALSE) 
  
  # Save as PDF
  pdf(file = paste0('up_chr_length6', ".pdf"), width = 6/ 2.54, height = 10 / 2.54)
  print(p)
  dev.off()
}
    '''
        with open('./chr2moddotplot3//3_VSat1in_up_rough/2_species_len/up_chr_length.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2moddotplot3//3_VSat1in_up_rough/2_species_len/'
        os.chdir(new_directory)
        subprocess.run(['Rscript up_chr_length.R'], shell=True)    
        os.chdir('../../../')
    
if argv1=="stepall" or argv1=="step3" or "step3b" in argv1:
    if argv1=="stepall" or argv1=="step3" or argv1=="step3b" or argv1=="step3b.0":  
        if  os.path.exists('./chr2moddotplot3/3_VSat1in_down_sum')==False:
            subprocess.run(["mkdir ./chr2moddotplot3/3_VSat1in_down_sum"], shell=True) 
        print('Process up class, merge method, for all positions, if Apos and Bpos intersect, merge AB')    
        ##Store in dictionary
        dict_samplechr_lines={}
        with open ('./chr2moddotplot3/2_slashin_updown','r') as f:
            next(f)
            for line in f:
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                ##
                one_sample=	        eachline_arr[0]
                #if one_sample!='PN40024':continue
                one_chr	=	        eachline_arr[1]
                one_type=	        eachline_arr[2]	
                start_x	=	        int(eachline_arr[3])
                start_y	=	        int(eachline_arr[4])
                end_x	=	        int(eachline_arr[5])
                end_y=	            int(eachline_arr[6])
                if one_type!='down':continue
                samplechr=one_sample+'|||'+one_chr
                if samplechr not in dict_samplechr_lines:
                    dict_samplechr_lines[samplechr]=[]
                dict_samplechr_lines[samplechr].append([start_x,start_y,end_x,end_y]) 
        result_list= []
        for samplechr,lines in dict_samplechr_lines.items():
            result_list.append([samplechr,lines])
        ##
        def run_step(one_result):
            samplechr,lines=one_result
            used_kkk=set()
            lines_num=len(lines)
            one_sample,one_chr=samplechr.split("|||")
            ###Re-plan the approach: throw ball 1 into frames 2, 3, 4 sequentially, if condition is met, expand frames 2/3/4
            ###First set all positions as frame sets
            i=0
            dict_index_set={}
            while i<lines_num:
                pos1,pos2,pos3,pos4=lines[i]
                dict_index_set[i]={}
                dict_index_set[i]['pos']=set()
                dict_index_set[i]['index']=set()
                dict_index_set[i]['index'].add(i)
                jj=pos1
                while jj<=pos3:
                    dict_index_set[i]['pos'].add(jj)
                    jj+=1
                jj=pos4                 ##############Difference between up and down is here
                while jj<=pos2:
                    dict_index_set[i]['pos'].add(jj)
                    jj+=1
                i+=1
            ##Start throwing balls, if thrown, record in bad_list
            bad_list=[]
            i=0 
            while i<lines_num:
                current_pos_set=dict_index_set[i]['pos']
                current_index_set=dict_index_set[i]['index']
                j=i+1
                while j<lines_num:
                    if j in bad_list:j+=1;continue
                    next_set=dict_index_set[j]['pos']
                    if not current_pos_set.isdisjoint(next_set):  # If there is an intersection
                        dict_index_set[j]['pos'].update(current_pos_set)
                        bad_list.append(i)
                        ##
                        dict_index_set[j]['index'].update(current_index_set)
                        #dict_index_set[j]['index'].remove(i)
                        break
                    j+=1
                i+=1
            ##Final count
            with open ('./chr2moddotplot3/3_VSat1in_down_sum/0_class','a') as f:
                i=0 
                while i<lines_num:
                    if i in bad_list:i+=1;continue
                    cluster_set=dict_index_set[i]['index']
                    for one_index in cluster_set:
                        start_x,start_y,end_x,end_y=lines[one_index]
                        f.write(f"{one_sample}\t{one_chr}\t{i}\t{start_x}\t{start_y}\t{end_x}\t{end_y}\n")
                    i+=1
                
     
                        
        with open ('./chr2moddotplot3/3_VSat1in_down_sum/0_class','w') as f:
                f.write("sample\tchromosome\tclass\tstart_x\tstart_y\tend_x\tend_y\n")        
            
        with Pool(processes=70) as pool:
             # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, result_list), start=1):
                # Process results here, e.g., store or print
                progress = (i / len(result_list)) * 100
                sys.stdout.write(f"\nProgress: {progress:.2f}% {str(i)} / {str(len(result_list))}#####################################################################\n")
                sys.stdout.flush()            
    if argv1=="stepall" or argv1=="step3" or argv1=="step3b" or argv1=="step3b.1":                
        print('Complete information')
        result_list=[]
        with open ('./chr2moddotplot3/3_VSat1in_down_sum/0_class','r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                #one_sample,one_chr,start_x,start_y,end_x,end_y=eachline_arr
                result_list.append(eachline_arr)
        sort_result_list = sorted(result_list, key=lambda x: x[4], reverse=False)
        sort_result_list = sorted(result_list, key=lambda x: (x[0],x[1], x[2], x[3], x[4]), reverse=False)
        print('Load original table')
        dict_index_info={}
        with open ('./chr2moddotplot3/2_slashin_updown','r') as f:
            next(f)
            for line in f:
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                ##
                one_sample=	        eachline_arr[0]
                one_chr	=	        eachline_arr[1]
                one_type=	        eachline_arr[2]
                start_x	=	        int(eachline_arr[3])
                start_y	=	        int(eachline_arr[4])
                end_x	=	        int(eachline_arr[5])
                end_y=	            int(eachline_arr[6])
                dict_index_info[f"{one_sample}|{one_chr}|{start_x}|{start_y}|{end_x}|{end_y}"]=eachline_arr[7:]
        with open ('./chr2moddotplot3/3_VSat1in_down_sum/1_class','w') as f:      
            f.write("sample\tchromosome\ttype\tclass\tstart_x\tstart_y\tend_x\tend_y\tseq1_len\tseq2_len\tdot_num\tdotskip_num\tdot_used\tline_slope\tline_R2\tback_accuracy1\tback_accuracy2\tback_distance1\tback_distance2\n")
            for one_result in sort_result_list:
                index=f"{one_result[0]}|{one_result[1]}|{one_result[3]}|{one_result[4]}|{one_result[5]}|{one_result[6]}"
                newline=f"{one_result[0]}\t{one_result[1]}\tdown\t{one_result[2]}\t{one_result[3]}\t{one_result[4]}\t{one_result[5]}\t{one_result[6]}\t"+'\t'.join(dict_index_info[index])+'\n'
                f.write(newline)                
    if argv1=="stepall" or argv1=="step3" or argv1=="step3b" or argv1=="step3b.2":
        print('20s')
        with open ('./chr2moddotplot3/3_VSat1in_down_sum/2_class','w') as f2:
            i=0
            with open ('./chr2moddotplot3/3_VSat1in_down_sum/1_class','r') as f:
                for line in f:
                    eachline=line.strip()
                    i+=1
                    if i==1: 
                        head=f"{eachline}\tregion_length\tregion_valid_length\tblank_length\tvalid_percent\n"
                        f2.write(head)
                        continue
                    eachline_arr=eachline.split('\t')
                    sample,chromosome,onetype,oneclass,start_x,start_y,end_x,end_y,seq1_len,seq2_len,dot_num,dotskip_num,dot_used,line_slope,line_R2,back_accuracy1,back_accuracy2,back_distance1,back_distance2=eachline_arr
                    start_x=int(start_x)   
                    start_y=int(start_y)   
                    end_x=int(end_x)   
                    end_y=int(end_y)   
                    region_length=start_y-start_x+1
                    region_valid_length=end_x-start_x+1+start_y-end_y+1
                    blank_length=end_y-end_x+1
                    valid_percent=round(1-blank_length/region_length,3)   
                    f2.write(f"{eachline}\t{region_length}\t{region_valid_length}\t{blank_length}\t{valid_percent}\n")
    if argv1=="stepall" or argv1=="step3" or argv1=="step3b" or argv1=="step3b.3":
        print('20s')
        dict_samplechr_info={}
        i=0
        with open ('./chr2moddotplot3/3_VSat1in_down_sum/2_class','r') as f:
            for line in f:
                eachline=line.strip()
                i+=1
                if i==1: head=eachline;continue
                eachline_arr=eachline.split('\t')
                sample,chromosome,onetype,oneclass,start_x,start_y,end_x,end_y,seq1_len,seq2_len,dot_num,dotskip_num,dot_used,line_slope,line_R2,back_accuracy1,back_accuracy2,back_distance1,back_distance2,region_length,region_valid_length,blank_length,valid_percent=eachline_arr
                samplechr=f"{sample}___{chromosome}"
                if samplechr not in dict_samplechr_info:
                    dict_samplechr_info[samplechr]=[]
                dict_samplechr_info[samplechr].append([oneclass,int(start_x),int(start_y),int(end_x),int(end_y),int(dot_num),dot_used,region_length,int(region_valid_length),blank_length,float(valid_percent)])
        print("#Load sample information")
        dict_sample_speciesgroup={}
        with open("/home/lain/aaa_data/run0/samples_satellite/sample_info",'r') as f:
            #Baimunage_hap2	Baimunage	 vinifera	Eurasian	Eurasian_vinifera_Baimunage	Table
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=6:continue
                sample=             eachline_arr[0]
                species_group=      eachline_arr[3]
                dict_sample_speciesgroup[sample]=species_group     
        print('Output')
        with open ('./chr2moddotplot3/3_VSat1in_down_sum/3_chrsum','w') as f:
            f.write(f"sample\tspecies_group\tchromosome\tclass\tstart_x\tstart_y\tend_x\tend_y\tdot_num\tdot_used\tregion_length\tregion_valid_length\tblank_length\tvalid_percent\tchr_IR_num\tchr_allIR_region_length\n")
            for samplechr,info_list in dict_samplechr_info.items():
                #if samplechr!="WoollyGrape_hap2___Chr7":continue
                sample,chromosome=samplechr.split('___')
                #sorted_items = sorted(info_list.items(), key=lambda x: x[1][4], reverse=True)        
                sorted_list = sorted(info_list, key=lambda x: x[8], reverse=True)
                oneclass,start_x,start_y,end_x,end_y,dot_num,dot_used,region_length,region_valid_length,blank_length,valid_percent=sorted_list[0]
                ###Calculate total region
                pos_set=set()
                chr_IR_num=0
                for one_info in info_list:
                    if valid_percent<0.2:continue
                    if region_valid_length<100000:continue
                    if dot_num<10:continue
                    pos_set|=set(range(one_info[1], one_info[3] + 1))
                    pos_set|=set(range(one_info[4], one_info[2] + 1))
                    chr_IR_num+=1
                chr_allIR_region_length=len(pos_set)   
                
                if chr_IR_num==0:continue
                chr_IR_num=f"({chr_IR_num}/{len(info_list)})"
                #print(sorted_list)
                #print(chr_allIR_region_length)
                ###Calculate ratio
                ###
                species_group=dict_sample_speciesgroup[sample]
                f.write(f"{sample}\t{species_group}\t{chromosome}\t{oneclass}\t{start_x}\t{start_y}\t{end_x}\t{end_y}\t{dot_num}\t{dot_used}\t{region_length}\t{region_valid_length}\t{blank_length}\t{valid_percent}\t{chr_IR_num}\t{chr_allIR_region_length}\n")
            
        print('')
    if argv1=="stepall" or argv1=="step3" or argv1=="step3b" or argv1=="step3b.4":        
        print('Plotting IR within VSat1')
        if  os.path.exists('./chr2moddotplot3/3_VSat1in_down_sum/4_plot')==False:
            subprocess.run(["mkdir ./chr2moddotplot3/3_VSat1in_down_sum/4_plot"], shell=True) 
        subprocess.run(["cp ./chr2moddotplot3/3_VSat1in_down_sum/3_chrsum ./chr2moddotplot3/3_VSat1in_down_sum/4_plot/inputfile"], shell=True)     
        R_txt=r'''
library(ggplot2)
library(dplyr)

print("")
{
  input_file=read.table('inputfile', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')




filtered_input <- input_file %>% 
  mutate(
        category = case_when(
          species_group == "Eurasian" ~ "Eurasian",
          species_group == "East_Asia" ~ "East_Asia",
          species_group == "America" ~ "America",
          #TRUE ~ "Other"
        ),
      )    
# Define color values
color_values <- c(
    "Eurasian" = "#ff3399",
    "East_Asia" = "#0066ff",
    "America" = "#33cc33"
)

filtered_input$category <- factor(filtered_input$category, levels = c("Eurasian", "East_Asia", "America"))

# Create plot object
p <- ggplot() +
    geom_boxplot(
        data = filtered_input,
        aes(x = category, y = region_valid_length, color = category),  # Optional: fill by group
        width = 0.5,               # Boxplot width
        outlier.shape = NA        # Hide boxplot outliers (since scatter plot shows them)
        
    ) +

    geom_jitter(
        data = filtered_input,
        aes(x = category, y = region_valid_length, color = category, size = chr_allIR_region_length / 150000, alpha = valid_percent),
        width = 0.2,  # Control horizontal jitter range
        shape=16, stroke = 0
        #alpha = 0.6  # Adjust transparency to avoid overlap
    )


p <- p+theme_classic() + scale_size_identity() +     
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    ) + scale_y_continuous(limits = c(0, NA), breaks = seq(0, 10000000, by = 250000)) + 
    scale_color_manual(values = color_values, drop = FALSE)
  
  # Save as PDF
  pdf(file = paste0('down_chr_length1', ".pdf"), width = 20/ 2.54, height = 17 / 2.54)
  print(p)
  dev.off()


p <- ggplot() +
    geom_point(
        data = input_file,
        aes(x = blank_length, y = region_valid_length, color = chr_allIR_region_length)
    ) + scale_size_identity()+

theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
    ) #+
    #scale_color_manual(values = color_values, drop = FALSE)
  
  # Save as PDF
  pdf(file = paste0('down_chr_length2', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
  print(p)
  dev.off()
}
    '''
        with open('./chr2moddotplot3/3_VSat1in_down_sum/4_plot/down_chr_length.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2moddotplot3/3_VSat1in_down_sum/4_plot/'
        os.chdir(new_directory)
        subprocess.run(['Rscript down_chr_length.R'], shell=True)    
        os.chdir('../../../')                          
    if argv1=="stepall" or argv1=="step3" or argv1=="step3b" or argv1=="step3b.5":        
        print('VSat1 internal, plot by position')
        if  os.path.exists('./chr2moddotplot3/3_VSat1in_down_sum/5_plot')==False:
            subprocess.run(["mkdir ./chr2moddotplot3/3_VSat1in_down_sum/5_plot"], shell=True)                 
        print('Load sample information table ../samples_satellite/2_good_regions')
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
                block_length=   int(eachline_arr[5])
                ########
                if samplechr not in dict_samplechr_allblocks:dict_samplechr_allblocks[samplechr]=[]
                dict_samplechr_allblocks[samplechr].append([         start      ,     end   ])          
        with open("./chr2moddotplot3/3_VSat1in_down_sum/5_plot/input_file",'w') as f2:    
            f2.write(f"sample\tchromosome\tcenter_revise\tIR_length\tIR_length_revise\n")
            with open("./chr2moddotplot3/3_VSat1in_down_sum/2_class",'r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    samplechr=eachline_arr[0]+'___'+eachline_arr[1]
                    center=int(eachline_arr[5])+int(round((int(eachline_arr[5])-int(eachline_arr[4]))/2,0))
                    MARK=''
                    for one_blocks in dict_samplechr_allblocks[samplechr]:
                        block_start=one_blocks[0]
                        block_end=one_blocks[1]
                        if center>block_start and center<block_end:
                            MARK='YES';break
                    if MARK=='':
                        #print('ERROR');
                        continue
                    center_revise=(center-block_start)/(block_end-block_start)    
                    seq1_len=int(eachline_arr[8])
                    seq2_len=int(eachline_arr[9])
                    IR_length=seq1_len+seq2_len
                    if IR_length>600000:    IR_length_revise=600000
                    else:                   IR_length_revise=IR_length
                    f2.write(f"{eachline_arr[0]}\t{eachline_arr[1]}\t{center_revise}\t{IR_length}\t{IR_length_revise}\n")
        R_txt=r'''
library(ggplot2)
library(dplyr)

print("")
{
input_file=read.table('input_file', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
input_file2 <- input_file %>% arrange(as.numeric(sub("Chr", "", chromosome)))
chr_order <- unique(input_file2$chromosome)
input_file2$chromosome <- factor(input_file2$chromosome, levels = chr_order)
        
# Create plot object
p <- ggplot() +
    geom_point(data = input_file2, aes(x = center_revise, y = IR_length_revise), size=2, shape=16, stroke = 0, alpha=0.7, color="#527a7a")     
p = p + ylim(0, 600000)    
p <- p + theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
      #panel.border = element_rect(color = "black", fill = NA, size = 0.4)
    )
p <- p + 
  facet_grid( chromosome ~ ., 
             scales = "free_x",  # If x-axis needs independent scaling across facets
             space = "free_x",   # If x-axis space needs to adapt
             switch = "x")   
    

  
  # Save as PDF
  pdf(file = paste0('IR', ".pdf"), width = 7.5/ 2.54, height = 20 / 2.54)
  print(p)
  dev.off()
}
    '''
        with open('./chr2moddotplot3/3_VSat1in_down_sum/5_plot/IR.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2moddotplot3//3_VSat1in_down_sum/5_plot/'
        os.chdir(new_directory)
        subprocess.run(['Rscript IR.R'], shell=True)    
        os.chdir('../../../')           


if argv1=="stepall" or argv1=="step4" or "step4a" in argv1:                
    if argv1=="stepall" or argv1=="step4" or argv1=="step4a" or argv1=="step4a.2":    
        print('4a.1 results are too complex, need simplification')
        print('Identify the most important position of VSat1, ultimately normalize to zero')
        print('Identify IR positions, no overlap possible, otherwise delete, left or right')
        print('Positions of VSat1 at the left and right edges of IR')
        if  os.path.exists('./chr2moddotplot3/4_VSat1out_down_sum')==False:
            subprocess.run(["mkdir ./chr2moddotplot3/4_VSat1out_down_sum"], shell=True) 
        print('Loading inter-VSat1')    
        dict_samplechr_VSat1_interarray={}
        dict_samplechr_VSat1_coreborder={}
        dict_samplechr_VSat1_largest_border={}
        dict_samplechr_VSat1_parts={}
        with open("/home/lain/aaa_data/run0/samples_satellite/2_good_regions_interarray",'r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                samplechr=eachline_arr[0]+"___"+eachline_arr[1]
                all_VSat1_length=int(eachline_arr[4])
                if all_VSat1_length<100000:continue
                inter_arrays_str=eachline_arr[9]
                bigblock_start=     int(eachline_arr[2])
                bigblock_end=       int(eachline_arr[3])
                dict_samplechr_VSat1_interarray[samplechr]=inter_arrays_str
                inter_arrays=inter_arrays_str.split('|')
                largest_part_pos1=int(eachline_arr[10])
                largest_part_pos2=int(eachline_arr[11])
                dict_samplechr_VSat1_largest_border[samplechr]=[largest_part_pos1,largest_part_pos2]
                dict_samplechr_VSat1_parts[samplechr]=eachline_arr[12]
                ##
                dict_samplechr_VSat1_coreborder[samplechr]=[bigblock_start,bigblock_end]
        print('Loading All-VSat1, some bigblocks may break the IR in the inter-array, load all VSat1 directly to determine left and right distances')    
        dict_samplechr_VSat1_array={}
        with open("/home/lain/aaa_data/run0/samples_satellite/2_good_regions",'r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=8:continue
                sample= eachline_arr[0]
                chromosome=eachline_arr[1].replace("region_",'Chr') 
                
                samplechr=sample+"___"+chromosome
                if samplechr not in dict_samplechr_VSat1_array:dict_samplechr_VSat1_array[samplechr]=[]
                dict_samplechr_VSat1_array[samplechr].append([int(eachline_arr[3]),int(eachline_arr[4])])
        
        #Given boundaries A (POS1, POS2) and B (POS3, POS4), find the length of overlap between A and B
        def calculate_overlap(POS1, POS2, POS3, POS4):
            # Calculate the start and end of the overlapping interval
            overlap_start = max(POS1, POS3)
            overlap_end = min(POS2, POS4)
            
            # Determine if there is overlap
            if overlap_start > overlap_end:
                return 0  # No overlap
            else:
                return overlap_end - overlap_start + 1  # Closed interval length (+1)
                
                
        print('Processing down class')  
        headline=f"sample\tchromosome\ttype\tstart_x\tstart_y\tend_x\tend_y\tseq1_len\tseq2_len\tdot_num\tdotskip_num\tdot_used\tline_slope\tline_R2\tback_accuracy1\tback_accuracy2\tback_distance1\tback_distance2\tcenter_pos\tcenter_blank\tIR_type0\tblank_percent\tvalid_type\tlargest_part_pos1\tlargest_part_pos2\tpos_type\tcore1\tcore2\tinterarray_percent\tIR_type\tdistance2part1\tVSat1core_left_neibour\tVSat1core_right_neibour\n"
        with open("./chr2moddotplot3/4_VSat1out_down_sum/2_lines_filter",'w') as f3:
            f3.write(headline)
            with open("./chr2moddotplot3/4_VSat1out_down_sum/2_lines",'w') as f2:
                f2.write(headline)
                bad_num=0;other_num=0;i=0
                with open("./chr2moddotplot3/2_slashout_updown",'r') as f:
                    next(f)
                    for line in f:
                        i+=1
                        print(f'Progress:{i}',end='\r')
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        samplechr=eachline_arr[0]+"___"+eachline_arr[1]
                        one_type=eachline_arr[2]
                        if one_type=='up':continue
                        if samplechr not in dict_samplechr_VSat1_interarray:bad_num+=1;continue   ##Effectively removes those with all_VSat1_length<20000, requiring presence of VSat1
                        
                        #if samplechr!="V022.hap1___Chr15":continue
                        
                        start_x=        int(eachline_arr[3])
                        start_y=        int(eachline_arr[4])
                        end_x=          int(eachline_arr[5])
                        end_y=          int(eachline_arr[6])
                        region_length=  start_y-start_x+1
                        #Calculate IR center point, blank ratio
                        center_pos=end_x+(end_y-end_x+1)/2
                        center_blank=end_y-end_x+1
                        blank_percent=center_blank/(start_y-start_x+1)
                        if center_blank>500000:IR_type0='distant'
                        else:IR_type0='continue'
                        if blank_percent<0.2: valid_type='0.8+'
                        elif  blank_percent<0.5: valid_type='0.5-0.8'
                        elif  blank_percent<0.8: valid_type='0.2-0.5'
                        else: valid_type='0.2-'
                        ##Determine if it is inter-array. Subtract boundaries directly, if the length within the boundary exceeds 50%, it is inter-array
                        core1,core2=dict_samplechr_VSat1_coreborder[samplechr]
                        overlap_length = calculate_overlap(start_x, start_y, core1, core2)
                        interarray_percent=     round(overlap_length/region_length,3)
                        if interarray_percent==1:  
                            IR_type='inter-array'
                        else:
                            IR_type='out-array'                    
                        ##Determine distance to the largest VSat1 (part, not original VSat1 block, directly connected within 50000bp)
                        largest_part_pos1,largest_part_pos2=dict_samplechr_VSat1_largest_border[samplechr]
                        largest_part_center=largest_part_pos1+(largest_part_pos2-largest_part_pos1)/2
                        distance2part1='.'
                        if  center_pos < largest_part_center:    
                            pos_type='left'
                            ##If any edge of IR is greater than largest_part_pos1, mark for deletion
                            if start_y>largest_part_pos1: IR_type='delete'
                            distance2part1=largest_part_pos1-start_y
                        elif center_pos > largest_part_center:       
                            pos_type='right'
                            ##If any edge of IR is less than largest_part_pos2, mark for deletion
                            if start_x<largest_part_pos2: IR_type='delete'
                            distance2part1=start_x-largest_part_pos2
                            #print(final_left,final_right)
                        
                        ##Find possible flanking VSat1
                        VSat1_parts=dict_samplechr_VSat1_parts[samplechr].split('|')
                        ###
                        Vsat1_blocks=dict_samplechr_VSat1_array[samplechr]
                        tmp_left=0;tmp_right=99999999
                        for one_block in Vsat1_blocks:
                            one1,one2=one_block
                            if one2<start_x and one2>tmp_left:tmp_left=one2
                            if one1>start_y and one1<tmp_right:tmp_right=one1
                        if tmp_left==0 or  start_x-tmp_left>1000000 or tmp_left<core1 or tmp_left>core2:tmp_left='.'
                        if tmp_right==99999999 or  tmp_right-start_y>1000000   or tmp_right<core1 or tmp_right>core2: tmp_right='.'
                        
                        if IR_type=='inter-array':
                            if tmp_left=='.' or tmp_right=='.':IR_type=='out-array'

                        newline=f"{eachline}\t{center_pos}\t{IR_type0}\t{center_blank}\t{blank_percent}\t{valid_type}\t{largest_part_pos1}\t{largest_part_pos2}\t{pos_type}\t{core1}\t{core2}\t{interarray_percent}\t{IR_type}\t{distance2part1}\t{tmp_left}\t{tmp_right}\n"
                        
                        f2.write(newline)
                        if IR_type!='delete' and blank_percent<0.5 and IR_type0=='continue':f3.write(newline)
                            
            print(f"Skipped {bad_num}, possibly no VSat1")   
            print(f"other_num={other_num}, may have crossed VSat1")     
    if argv1=="stepall" or argv1=="step4" or argv1=="step4a" or argv1=="step4a.3":        
        print('Loading VSat5')
        dict_sample_chr_poslist={}
        with open("/home/lain/aaa_data/run0/stat_plot/0-region2info",'r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')                
                sample,chromosome,chromosome_new,centype,chr_start,chr_end,length,strand,match_percent=eachline_arr
                if centype!='cen_355':continue
                sample_chr=f'{sample}___{chromosome_new}'
                if sample_chr not in dict_sample_chr_poslist:dict_sample_chr_poslist[sample_chr]=[]
                dict_sample_chr_poslist[sample_chr].append([int(chr_start),int(chr_end)])
                
        print('#Determine if VSat5 is contained')
        with open("./chr2moddotplot3/4_VSat1out_down_sum/3_lines_filter",'w') as f2:
            with open("./chr2moddotplot3/4_VSat1out_down_sum/2_lines_filter",'r') as f:
                i=0
                for line in f:
                    i+=1
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    if i==1:
                        head=eachline
                        f2.write(f"{eachline}\tVSat5_state\n")
                    else:
                        sample_chr=f"{eachline_arr[0]}___{eachline_arr[1]}"
                        mark='no-VSat5'
                        if sample_chr not in dict_sample_chr_poslist:
                            mark='no-VSat5'
                        else:    
                            pos_list=dict_sample_chr_poslist[sample_chr]
                            start=int(eachline_arr[3])
                            end=int(eachline_arr[4])
                            for one_pos in pos_list:
                                tmp_pos1,tmp_pos2=one_pos
                                if start<tmp_pos1 and end>tmp_pos1:   mark='VSat5_related';break
                                elif start<tmp_pos2  and end>tmp_pos2:    mark='VSat5_related';break
                                elif start<tmp_pos1 and end>tmp_pos2:   mark='VSat5_related';break
                                elif start>tmp_pos1 and end<tmp_pos2:   mark='VSat5_related';break
                        f2.write(f"{eachline}\t{mark}\n")
        #############################################################      20260427,,, Supplement a judgment to determine if LIRSat is contained          
    if argv1=="stepall" or argv1=="step4" or argv1=="step4a" or argv1=="step4a.4":
        print('step4a.3')
        with open ("./chr2moddotplot3/4_VSat1out_down_sum/4_plot_input_neibour",'w') as f3:
            f3.write(f"sample\tchromosome\tIR_center_percent\tIR_y\tpos_type\tIR_left_percent\tIR_right_percent\tVSat1core_neibour_space\tVSat5_state\n")
            with open ("./chr2moddotplot3/4_VSat1out_down_sum/4_plot_input",'w') as f2:
                f2.write(f"sample\tchromosome\tIR_center\tIR_y\tIR_left_border\tIR_right_border\tIR_type\tIR_type2\tVSat5_state\n")            
                with open ("./chr2moddotplot3/4_VSat1out_down_sum/3_lines_filter",'r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        sample,chromosome,onetype,start_x,start_y,end_x,end_y,seq1_len,seq2_len,dot_num,dotskip_num,dot_used,line_slope,line_R2,back_accuracy1,back_accuracy2,back_distance1,back_distance2,center_pos,center_blank,IR_type0,blank_percent,valid_type,largest_part_pos1,largest_part_pos2,pos_type,core1,core2,interarray_percent,IR_type,distance2part1,VSat1core_left_neibour,VSat1core_right_neibour,VSat5_state=eachline_arr
                        if pos_type=='left':        IR_center=-int(distance2part1)
                        else:                       IR_center=int(distance2part1)
                        IR_y=(int(seq1_len)+int(seq2_len))
                        IR_left_border=IR_center-(float(center_pos)-int(start_x))
                        IR_right_border=IR_center+(int(start_y)-float(center_pos))
                        if IR_type=='inter-array':IR_type2='inter-array'
                        else:IR_type2=IR_type+'_'+eachline_arr[25]
                        newline=f"{sample}\t{chromosome}\t{IR_center}\t{IR_y}\t{IR_left_border}\t{IR_right_border}\t{IR_type}\t{IR_type2}\t{VSat5_state}\n"
                        f2.write(newline)
                        ###########
                        if IR_type=='inter-array':
                            if VSat1core_left_neibour=='.' or VSat1core_right_neibour=='.':continue
                            VSat1core_neibour_space=int(VSat1core_right_neibour)-int(VSat1core_left_neibour)+1
                            IR_center_percent=round((float(center_pos)-int(VSat1core_left_neibour))/VSat1core_neibour_space,3)
                            IR_left_percent=round((float(center_pos)-int(IR_left_border))/VSat1core_neibour_space,3)
                            IR_right_percent=round((float(center_pos)-int(IR_right_border))/VSat1core_neibour_space,3)
                            newline=f"{sample}\t{chromosome}\t{IR_center_percent}\t{IR_y}\t{pos_type}\t{IR_left_percent}\t{IR_right_percent}\t{VSat1core_neibour_space}\t{VSat5_state}\n"
                            f3.write(newline)    
                      
    if argv1=="stepall" or argv1=="step4" or argv1=="step4a" or argv1=="step4a.5" and 1==2:      
        if  os.path.exists('./chr2moddotplot3/4_VSat1out_down_sum/5_plot')==False:
            subprocess.run(["mkdir ./chr2moddotplot3/4_VSat1out_down_sum/5_plot"], shell=True) 
        subprocess.run(["cp ./chr2moddotplot3/4_VSat1out_down_sum/4_plot_input ./chr2moddotplot3/4_VSat1out_down_sum/5_plot/input_file"], shell=True)     
        R_txt=r'''
library(ggplot2)
library(dplyr)

print("")
{
input_file=read.table('input_file', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
input_file$chromosome <- factor(input_file$chromosome, levels = c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"))

input_file$IR_type2 <- factor(input_file$IR_type2, levels = c("out-array_left","inter-array","out-array_right"))           
input_file2 <- input_file %>% filter(IR_type2 != "inter-array" ) 
             
# Define color values
color_values <- c(
                "out-array"="grey",
                "inter-array"="blue"
)

# Create plot object
p <- ggplot() +
    geom_point(data = input_file2,aes(x = IR_center, y = IR_y, color = VSat5_state), size=0.5) #, size=IR_y/200000 + scale_size_identity()#+  color=IR_type,
    #geom_segment(data = input_file2, aes(x = IR_center + IR_left, xend = IR_center + IR_right, y = IR_y, yend = IR_y, color = IR_type), size = 0.5 )
    #coord_equal() +
p = p + ylim(0, 600000)     
p <- p + 
  facet_grid( chromosome ~ IR_type2, 
             scales = "free_x",  # If x-axis needs independent scaling across facets
             space = "free_x",   # If x-axis space needs to adapt
             switch = "x")   
p <- p + theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      legend.position = "none",
      #axis.text.x = element_blank(),
      #panel.border = element_rect(color = "black", fill = NA, size = 0.4)
    ) #+
    #scale_color_manual(values = color_values, drop = FALSE)
    

  
  # Save as PDF
  pdf(file = paste0('IR1', ".pdf"), width = 10/ 2.54, height = 20 / 2.54)
  print(p)
  dev.off()
}
    '''
        with open('./chr2moddotplot3/4_VSat1out_down_sum/5_plot/IR.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2moddotplot3//4_VSat1out_down_sum/5_plot/'
        os.chdir(new_directory)
        subprocess.run(['Rscript IR.R'], shell=True)    
        os.chdir('../../../')       
    if argv1=="stepall" or argv1=="step4" or argv1=="step4a" or argv1=="step4a.5" and 1==2:      
        subprocess.run(["cp ./chr2moddotplot3/4_VSat1out_down_sum/4_plot_input_neibour ./chr2moddotplot3/4_VSat1out_down_sum/5_plot/input_file2"], shell=True)  
        R_txt=r'''
library(ggplot2)
library(dplyr)

print("")
{
input_file=read.table('input_file2', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
input_file$chromosome <- factor(input_file$chromosome, levels = c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"))
        
        
# Define color values
color_values <- c(
                #"out-array"="grey",
                #"inter-array"="blue"
                "left"='#ff9900',
                "right"='#666699'
)

# Create plot object
p <- ggplot() +
    geom_point(data = input_file, aes(x = IR_center_percent, y = IR_y, color = pos_type), size=0.5)# , size=IR_y/200000+ scale_size_identity()#+
    #geom_segment(data = input_file, aes(x = IR_center_percent + IR_left_percent, xend = IR_center_percent + IR_right_percent, y = IR_y, yend = IR_y), size = 0.5, alpha=0.5 )
    #coord_equal() +
#p = p + ylim(0, 600000)    
p = p + xlim(0, 1) 
p <- p + facet_wrap(~ chromosome, ncol = 1, strip.position = "left") 
p <- p + theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
      #panel.border = element_rect(color = "black", fill = NA, size = 0.4)
    ) +
    scale_color_manual(values = color_values, drop = FALSE)
    

  
  # Save as PDF
  pdf(file = paste0('IR2', ".pdf"), width = 10/ 2.54, height = 20 / 2.54)
  print(p)
  dev.off()
}
    '''
        with open('./chr2moddotplot3/4_VSat1out_down_sum/5_plot/IR.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2moddotplot3//4_VSat1out_down_sum/5_plot/'
        os.chdir(new_directory)
        subprocess.run(['Rscript IR.R'], shell=True)    
        os.chdir('../../../')           
    
    if argv1=="stepall" or argv1=="step4" or argv1=="step4a" or argv1=="step4a.6":
        with open ('./chr2moddotplot3/4_VSat1out_down_sum/6_plot_input','w') as f2:
            f2.write(f"sample\tchromosome\tcenter_revise\tIR_y\tIR_type2\tVSat5_state\n")
            with open ('./chr2moddotplot3/4_VSat1out_down_sum/4_plot_input','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split("\t")
                    sample,chromosome,IR_center,IR_y,IR_left_border,IR_right_border,IR_type,IR_type2,VSat5_state=eachline_arr
                    if IR_type!='out-array':continue
                    IR_center=float(IR_center)
                    if IR_type2=='out-array_left':
                        if IR_center<-20000000:continue
                        center_revise=1-(-IR_center)/30000000+0     ##This 3000000 adjusts the proportion of the three parts of the image
                    if IR_type2=='out-array_right':
                        if IR_center>20000000:continue   
                        center_revise=(IR_center)/30000000+4
                    f2.write(f"{sample}\t{chromosome}\t{center_revise}\t{IR_y}\t{IR_type2}\t{VSat5_state}\n")
            with open ('./chr2moddotplot3/4_VSat1out_down_sum/4_plot_input_neibour','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split("\t")
                    sample,chromosome,IR_center_percent,IR_y,pos_type,IR_left_percent,IR_right_percent,VSat1core_neibour_space,VSat5_state=eachline_arr
                    IR_center_percent=float(IR_center_percent)
                    center_revise=IR_center_percent+2
                    f2.write(f"{sample}\t{chromosome}\t{center_revise}\t{IR_y}\tinter-array\t{VSat5_state}\n")    
    if argv1=="stepall" or argv1=="step4" or argv1=="step4a" or argv1=="step4a.7":   
        if  os.path.exists('./chr2moddotplot3/4_VSat1out_down_sum/7_plot')==False:
            subprocess.run(["mkdir ./chr2moddotplot3/4_VSat1out_down_sum/7_plot"], shell=True) 
        subprocess.run(["cp ./chr2moddotplot3/4_VSat1out_down_sum/6_plot_input ./chr2moddotplot3/4_VSat1out_down_sum/7_plot/input_file"], shell=True)  
        R_txt=r'''
library(ggplot2)
library(dplyr)

print("")
{
input_file=read.table('input_file', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
input_file$chromosome <- factor(input_file$chromosome, levels = c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"))
input_file$IR_type2 <- factor(input_file$IR_type2, levels = c("out-array_left","inter-array","out-array_right"))      
       
filtered_input <- input_file %>% 
  #filter(maxline_unit_number >1.8 ) %>%
  mutate(
        category = case_when(
            VSat5_state == "VSat5_related" ~ "VSat5_related",
          IR_type2 == "inter-array" ~ "inter",
          
          IR_type2 == "out-array_left" ~ "left",
          
          IR_type2 == "out-array_right" ~ "right"
          #TRUE ~ "Other"
        ),
      )         
        
# Define color values
color_values <- c(
                "left"='#527a7a',
                "right"='#527a7a',
                "inter"='#00b300',
                "VSat5_related"='#cc6600'
)

# Create plot object
p <- ggplot() +
    geom_point(data = filtered_input, aes(x = center_revise, y = IR_y, color = category), size=2, shape=16, stroke=0, alpha=0.7)            # , size=IR_y/200000+ scale_size_identity()#+
    #geom_segment(data = filtered_input, aes(x = IR_center_percent + IR_left_percent, xend = IR_center_percent + IR_right_percent, y = IR_y, yend = IR_y), size = 0.5, alpha=0.5 )
    #coord_equal() +
#p = p + ylim(0, 600000)    
#p = p + xlim(0, 1) 
p <- p + 
  facet_grid( chromosome ~ IR_type2, 
             scales = "free_x",  # If x-axis needs independent scaling across facets
             space = "free_x",   # If x-axis space needs to adapt
             switch = "x")  
p <- p + theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
      #panel.border = element_rect(color = "black", fill = NA, size = 0.4)
    ) +
    scale_color_manual(values = color_values, drop = FALSE)
    

  
  # Save as PDF
  pdf(file = paste0('IR', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
  print(p)
  dev.off()
}
    '''
        with open('./chr2moddotplot3/4_VSat1out_down_sum/7_plot/IR.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2moddotplot3//4_VSat1out_down_sum/7_plot/'
        os.chdir(new_directory)
        subprocess.run(['Rscript IR.R'], shell=True)    
        os.chdir('../../../')           
        
    #Determine if LIRSat123 is contained   
    if argv1=="stepall" or argv1=="step4" or argv1=="step4a" or argv1=="step4a.3_b":  
        print('#Determine if LIRSat123 is contained')
        dict_samplechr_pos={}
        with open("/home/lain/aaa_data/run0/stat_plot/0-region2info",'r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                sample,chromosome,chromosome_new,centype,chr_start,chr_end,length,strand,match_percent=eachline_arr
                samplechr=f"{sample}___{chromosome_new}"
                if samplechr not in dict_samplechr_pos:dict_samplechr_pos[samplechr]=[]
                if (abs(int(chr_start)-int(chr_end))<5000):continue
                dict_samplechr_pos[samplechr].append([int(chr_start),int(chr_end)])
        with open("./chr2moddotplot3/4_VSat1out_down_sum/3_lines_filter_addition",'w') as f2:
            with open("./chr2moddotplot3/4_VSat1out_down_sum/3_lines_filter",'r') as f:
                i=0
                for line in f:
                    i+=1
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    if i==1:
                        head=eachline
                        f2.write(f"{eachline}\tLIRSat_state\n")
                    else:
                        sample_chr=f"{eachline_arr[0]}___{eachline_arr[1]}"
                        mark='no-LIRSat'
                        if sample_chr not in dict_samplechr_pos:
                            mark='no-LIRSat'
                        else:    
                            pos_list=dict_samplechr_pos[sample_chr]
                            start=int(eachline_arr[3])
                            end=int(eachline_arr[4])
                            for one_pos in pos_list:
                                tmp_pos1,tmp_pos2=one_pos
                                if start<tmp_pos1 and end>tmp_pos1:   mark='LIRSat_have';break
                                elif start<tmp_pos2  and end>tmp_pos2:    mark='LIRSat_have';break
                                elif start<tmp_pos1 and end>tmp_pos2:   mark='LIRSat_have';break
                                elif start>tmp_pos1 and end<tmp_pos2:   mark='LIRSat_have';break
                        f2.write(f"{eachline}\t{mark}\n")                            
    if argv1=="stepall" or argv1=="step4" or argv1=="step4a" or argv1=="step4a.4_b":
        print('step4a.3')
        with open ("./chr2moddotplot3/4_VSat1out_down_sum/4_plot_input_neibour_addition",'w') as f3:
            f3.write(f"sample\tchromosome\tIR_center_percent\tIR_y\tpos_type\tIR_left_percent\tIR_right_percent\tVSat1core_neibour_space\tVSat5_state\tLIRSat_state\n")
            with open ("./chr2moddotplot3/4_VSat1out_down_sum/4_plot_input_addition",'w') as f2:
                f2.write(f"sample\tchromosome\tIR_center\tIR_y\tIR_left_border\tIR_right_border\tIR_type\tIR_type2\tVSat5_state\tLIRSat_state\n")            
                with open ("./chr2moddotplot3/4_VSat1out_down_sum/3_lines_filter_addition",'r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        sample,chromosome,onetype,start_x,start_y,end_x,end_y,seq1_len,seq2_len,dot_num,dotskip_num,dot_used,line_slope,line_R2,back_accuracy1,back_accuracy2,back_distance1,back_distance2,center_pos,center_blank,IR_type0,blank_percent,valid_type,largest_part_pos1,largest_part_pos2,pos_type,core1,core2,interarray_percent,IR_type,distance2part1,VSat1core_left_neibour,VSat1core_right_neibour,VSat5_state,LIRSat_state=eachline_arr
                        if pos_type=='left':        IR_center=-int(distance2part1)
                        else:                       IR_center=int(distance2part1)
                        IR_y=(int(seq1_len)+int(seq2_len))
                        IR_left_border=IR_center-(float(center_pos)-int(start_x))
                        IR_right_border=IR_center+(int(start_y)-float(center_pos))
                        if IR_type=='inter-array':IR_type2='inter-array'
                        else:IR_type2=IR_type+'_'+eachline_arr[25]
                        newline=f"{sample}\t{chromosome}\t{IR_center}\t{IR_y}\t{IR_left_border}\t{IR_right_border}\t{IR_type}\t{IR_type2}\t{VSat5_state}\t{LIRSat_state}\n"
                        f2.write(newline)
                        ###########
                        if IR_type=='inter-array':
                            if VSat1core_left_neibour=='.' or VSat1core_right_neibour=='.':continue
                            VSat1core_neibour_space=int(VSat1core_right_neibour)-int(VSat1core_left_neibour)+1
                            IR_center_percent=round((float(center_pos)-int(VSat1core_left_neibour))/VSat1core_neibour_space,3)
                            IR_left_percent=round((float(center_pos)-int(IR_left_border))/VSat1core_neibour_space,3)
                            IR_right_percent=round((float(center_pos)-int(IR_right_border))/VSat1core_neibour_space,3)
                            newline=f"{sample}\t{chromosome}\t{IR_center_percent}\t{IR_y}\t{pos_type}\t{IR_left_percent}\t{IR_right_percent}\t{VSat1core_neibour_space}\t{VSat5_state}\t{LIRSat_state}\n"
                            f3.write(newline)                    
    if argv1=="stepall" or argv1=="step4" or argv1=="step4a" or argv1=="step4a.6_b":
        with open ('./chr2moddotplot3/4_VSat1out_down_sum/6_plot_input_addition','w') as f2:
            f2.write(f"sample\tchromosome\tcenter_revise\tIR_y\tIR_type2\tVSat5_state\tLIRSat_state\n")
            with open ('./chr2moddotplot3/4_VSat1out_down_sum/4_plot_input_addition','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split("\t")
                    sample,chromosome,IR_center,IR_y,IR_left_border,IR_right_border,IR_type,IR_type2,VSat5_state,LIRSat_state=eachline_arr
                    if IR_type!='out-array':continue
                    IR_center=float(IR_center)
                    if IR_type2=='out-array_left':
                        if IR_center<-20000000:continue
                        center_revise=1-(-IR_center)/30000000+0     ##This 3000000 adjusts the proportion of the three parts of the image
                    if IR_type2=='out-array_right':
                        if IR_center>20000000:continue   
                        center_revise=(IR_center)/30000000+4
                    f2.write(f"{sample}\t{chromosome}\t{center_revise}\t{IR_y}\t{IR_type2}\t{VSat5_state}\t{LIRSat_state}\n")
            with open ('./chr2moddotplot3/4_VSat1out_down_sum/4_plot_input_neibour_addition','r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split("\t")
                    sample,chromosome,IR_center_percent,IR_y,pos_type,IR_left_percent,IR_right_percent,VSat1core_neibour_space,VSat5_state,LIRSat_state=eachline_arr
                    IR_center_percent=float(IR_center_percent)
                    center_revise=IR_center_percent+2
                    f2.write(f"{sample}\t{chromosome}\t{center_revise}\t{IR_y}\tinter-array\t{VSat5_state}\t{LIRSat_state}\n")                            
    if argv1=="stepall" or argv1=="step4" or argv1=="step4a" or argv1=="step4a.7_b":   
        if  os.path.exists('./chr2moddotplot3/4_VSat1out_down_sum/7_plot_b')==False:
            subprocess.run(["mkdir ./chr2moddotplot3/4_VSat1out_down_sum/7_plot_b"], shell=True) 
        subprocess.run(["cp ./chr2moddotplot3/4_VSat1out_down_sum/6_plot_input_addition ./chr2moddotplot3/4_VSat1out_down_sum/7_plot_b/input_file"], shell=True)  
        R_txt=r'''
library(ggplot2)
library(dplyr)

print("")
{
input_file=read.table('input_file', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
input_file$chromosome <- factor(input_file$chromosome, levels = c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"))
input_file$IR_type2 <- factor(input_file$IR_type2, levels = c("out-array_left","inter-array","out-array_right"))      
       
filtered_input <- input_file %>% 
  #filter(maxline_unit_number >1.8 ) %>%
  mutate(
        category = case_when(
            
            LIRSat_state == "no-LIRSat" ~ "no-LIRSat",
            LIRSat_state == "LIRSat_have" ~ "LIRSat_have",
            VSat5_state == "VSat5_related" ~ "VSat5_related"
          #IR_type2 == "inter-array" ~ "inter",
          
          #IR_type2 == "out-array_left" ~ "left",
          
          #IR_type2 == "out-array_right" ~ "right"
          #TRUE ~ "Other"
        ),
      )         
        
# Define color values
color_values <- c(
                "left"='#527a7a',
                "right"='#527a7a',
                "inter"='#00b300',
                "VSat5_related"='#cc6600',
                "no-LIRSat"='grey',
                "LIRSat_have"='#00b300'
)

# Create plot object
p <- ggplot() +
    geom_point(data = filtered_input, aes(x = center_revise, y = IR_y, color = category), size=2, shape=16, stroke=0, alpha=0.7)            # , size=IR_y/200000+ scale_size_identity()#+
    #geom_segment(data = filtered_input, aes(x = IR_center_percent + IR_left_percent, xend = IR_center_percent + IR_right_percent, y = IR_y, yend = IR_y), size = 0.5, alpha=0.5 )
    #coord_equal() +
#p = p + ylim(0, 600000)    
#p = p + xlim(0, 1) 
p <- p + 
  facet_grid( chromosome ~ IR_type2, 
             scales = "free_x",  # If x-axis needs independent scaling across facets
             space = "free_x",   # If x-axis space needs to adapt
             switch = "x")  
p <- p + theme_classic() +         
    theme(
      #axis.ticks.y = element_blank(),
      #axis.text.y = element_blank(),
      #legend.position = "none",
      #axis.text.x = element_blank()
      #panel.border = element_rect(color = "black", fill = NA, size = 0.4)
    ) +
    scale_color_manual(values = color_values, drop = FALSE)
    

  
  # Save as PDF
  pdf(file = paste0('IR', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
  print(p)
  dev.off()
}
    '''
        with open('./chr2moddotplot3/4_VSat1out_down_sum/7_plot_b/IR.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = './chr2moddotplot3//4_VSat1out_down_sum/7_plot_b/'
        os.chdir(new_directory)
        subprocess.run(['Rscript IR.R'], shell=True)    
        os.chdir('../../../')      
        
        
if argv1=="part1" or  argv1=="part2" :
    if  argv1=="part1": mode='static'
    elif argv1=="part2":mode='interactive'
    sample=sys.argv[2]
    chromosome=sys.argv[3]
    pos=sys.argv[4]
    print(f'Input parameters: {sample}\t\t{chromosome}\t\t{pos}')
    if  os.path.exists('./chr2moddotplot3/part')==False: 
        subprocess.run(["mkdir ./chr2moddotplot3/part"], shell=True)   

    
    input_name=f"{sample}:{chromosome}:{pos}"
    input_file=f"./chr2moddotplot3/part/{input_name}/{input_name}.fa"
    subprocess.run([f"mkdir ./chr2moddotplot3/part/{input_name}"], shell=True) 
    subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample}/{sample}.fasta {chromosome}:{pos} > {input_file}  "], shell=True)  

    output_dir=f"./chr2moddotplot3/part/{input_name}"

    moddotplot_software='ulimit -u 100000 && source /home/lain/software/venv/ModDotPlot/venv/bin/activate && /home/lain/software/venv/ModDotPlot/venv/bin/moddotplot '
    
    if mode=='static':
        subprocess.run([f"{moddotplot_software} static  --fasta {input_file} --kmer 21 --output-dir {output_dir} --identity 80  --no-bed --no-hist"], shell=True)     #--dpi 3000 --resolution 5000 
    elif mode=='interactive':
        import random
        port=8050+random.randint(0,300)
        print(f"port:{port}")
        # 1. Find processes occupying the port   # 2. Parse output to get PID (assuming output format matches expectations) # 3. Terminate all related processes
        result = subprocess.run(["lsof", "-i", f":{port}"], capture_output=True, text=True)
        output = result.stdout
        pids = []
        for line in output.split('\n'):
            if 'PID' in line:  # Skip header
                continue
            parts = line.split()
            if len(parts) > 1:
                pids.append(parts[1])  # Assuming PID is the second column
        for pid in pids:
            subprocess.run(["kill", "-9", pid])

        subprocess.run([f"{moddotplot_software} interactive --fasta {input_file} --identity 80 --port {port}"], shell=True)   

if argv1=="part1" or  argv1=="part3" :
    mode='interactive'
    sample=sys.argv[2]
    chromosome=sys.argv[3]
    pos=sys.argv[4]
    print(f'Input parameters: {sample}\t\t{chromosome}\t\t{pos}')
    if  os.path.exists('./chr2moddotplot3/part')==False: 
        subprocess.run(["mkdir ./chr2moddotplot3/part"], shell=True)   

    
    input_name=f"{sample}:{chromosome}:{pos}"
    input_file=f"./chr2moddotplot3/part/{input_name}/{input_name}.fa"
    subprocess.run([f"mkdir ./chr2moddotplot3/part/{input_name}"], shell=True) 
    subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample}/{sample}.fasta {chromosome}:{pos} > {input_file}  "], shell=True)  

    output_dir=f"./chr2moddotplot3/part/{input_name}"

    moddotplot_software='ulimit -u 100000 && source /home/lain/software/venv/ModDotPlot/venv/bin/activate && /home/lain/software/venv/ModDotPlot/venv/bin/moddotplot '
    
    if mode=='static':
        subprocess.run([f"{moddotplot_software} static  --fasta {input_file} --kmer 21 --output-dir {output_dir} --identity 80  --resolution 500  --no-bed --no-hist"], shell=True)     #--dpi 3000 --resolution 5000 
    elif mode=='interactive':
        import random
        port=8050+random.randint(0,300)
        print(f"port:{port}")
        # 1. Find processes occupying the port   # 2. Parse output to get PID (assuming output format matches expectations) # 3. Terminate all related processes
        result = subprocess.run(["lsof", "-i", f":{port}"], capture_output=True, text=True)
        output = result.stdout
        pids = []
        for line in output.split('\n'):
            if 'PID' in line:  # Skip header
                continue
            parts = line.split()
            if len(parts) > 1:
                pids.append(parts[1])  # Assuming PID is the second column
        for pid in pids:
            subprocess.run(["kill", "-9", pid])

        subprocess.run([f"{moddotplot_software} interactive --fasta {input_file} --identity 80 --resolution 300  --port {port}"], shell=True)   
        
if argv1=="vs1" or  argv1=="vs2" :
    if  argv1=="vs1": mode='static'
    elif argv1=="vs2":mode='interactive'
    sample1=sys.argv[2]
    chromosome1=sys.argv[3]
    pos1=sys.argv[4]
    sample2=sys.argv[5]
    chromosome2=sys.argv[6]
    pos2=sys.argv[7]
    print(f'Input parameters: {sample1}\t\t{chromosome1}\t\t{pos1}')
    print(f'Input parameters: {sample2}\t\t{chromosome2}\t\t{pos2}')
    if  os.path.exists('./chr2moddotplot3/vs')==False: 
        subprocess.run(["mkdir ./chr2moddotplot3/vs"], shell=True)   

    
    input_name1=f"{sample1}:{chromosome1}:{pos1}"
    input_file1=f"./chr2moddotplot3/vs/{input_name1}/{input_name1}.fa"
    subprocess.run([f"mkdir ./chr2moddotplot3/vs/{input_name1}"], shell=True) 
    if pos1=='all':
        subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample1}/{sample1}.fasta {chromosome1} > {input_file1}  "], shell=True)  
    else:    
        subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample1}/{sample1}.fasta {chromosome1}:{pos1} > {input_file1}  "], shell=True)  
    #
    input_name2=f"{sample2}:{chromosome2}:{pos2}"
    input_file2=f"./chr2moddotplot3/vs/{input_name2}/{input_name2}.fa"
    subprocess.run([f"mkdir ./chr2moddotplot3/vs/{input_name2}"], shell=True) 
    if pos2=='all':
        subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample2}/{sample2}.fasta {chromosome2} > {input_file2}  "], shell=True)  
    else:    
        subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample2}/{sample2}.fasta {chromosome2}:{pos2} > {input_file2}  "], shell=True)  
    
    
    output_name=f"{input_name1}_{input_name2}"
    
    output_dir=f"./chr2moddotplot3/vs/{output_name}"
    output_file=f"{output_dir}/{output_name}.fa"
    subprocess.run([f"mkdir {output_dir}"], shell=True) 
    seq1='';seq2=''
    with open(input_file1,'r') as f:
        next(f)
        for line in f:
            seq1+=line.strip()
    with open(input_file2,'r') as f:
        next(f)
        for line in f:
            seq2+=line.strip()            
    with open(f"{output_file}",'w') as f:
        f.write(f">{input_name1}\n{seq1}\n>{input_name2}\n{seq2}\n")
            
    moddotplot_software='ulimit -u 100000 && source /home/lain/software/venv/ModDotPlot/venv/bin/activate &&  /home/lain/software/venv/ModDotPlot/venv/bin/moddotplot '

    if mode=='static':
        subprocess.run([f"{moddotplot_software} static  --fasta {output_file} --kmer 21 -w 5000 --output-dir {output_dir} --identity 80  --no-bed --no-hist --compare-only"], shell=True)     #--dpi 3000 --resolution 5000 
    elif mode=='interactive':
        # 1. Find processes occupying port 8050   # 2. Parse output to get PID (assuming output format matches expectations) # 3. Terminate all related processes
        result = subprocess.run(["lsof", "-i", ":8050"], capture_output=True, text=True)
        output = result.stdout
        pids = []
        for line in output.split('\n'):
            if 'PID' in line:  # Skip header
                continue
            parts = line.split()
            if len(parts) > 1:
                pids.append(parts[1])  # Assuming PID is the second column
        for pid in pids: 
            subprocess.run(["kill", "-9", pid])
        
        
        subprocess.run([f"{moddotplot_software} interactive --fasta {output_file} --identity 80 --compare-only --port 8050"], shell=True)   
     
if argv1=="vs_all"  :
    argv1=="vs1"
    mode='static'
    chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19","Chr20"]
    for one_chr in chr_list:
        sample1=sys.argv[2]
        sample2=sys.argv[3]
        chromosome1=one_chr
        chromosome2=one_chr
        #if chromosome2!='Chr20':continue
        #if chromosome2=='Chr20':chromosome2="Chr7"
        pos1='all'
        pos2='all'
        pos1='10000000-20000000'
        pos2='10000000-20000000'        
        print(f'Input parameters: {sample1}\t\t{chromosome1}\t\t{pos1}')
        print(f'Input parameters: {sample2}\t\t{chromosome2}\t\t{pos2}')
        if  os.path.exists('./chr2moddotplot3/vs_all')==False: 
            subprocess.run(["mkdir ./chr2moddotplot3/vs_all"], shell=True)   
    
        
        input_name1=f"{sample1}:{chromosome1}:{pos1}"
        input_file1=f"./chr2moddotplot3/vs_all/{input_name1}/{input_name1}.fa"
        subprocess.run([f"mkdir ./chr2moddotplot3/vs_all/{input_name1}"], shell=True) 
        if pos1=='all':
            subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample1}/{sample1}.fasta {chromosome1} > {input_file1}  "], shell=True)  
        else:    
            subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample1}/{sample1}.fasta {chromosome1}:{pos1} > {input_file1}  "], shell=True)  
        #
        input_name2=f"{sample2}:{chromosome2}:{pos2}"
        input_file2=f"./chr2moddotplot3/vs_all/{input_name2}/{input_name2}.fa"
        subprocess.run([f"mkdir ./chr2moddotplot3/vs_all/{input_name2}"], shell=True) 
        if pos2=='all':
            subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample2}/{sample2}.fasta {chromosome2} > {input_file2}  "], shell=True)  
        else:    
            subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample2}/{sample2}.fasta {chromosome2}:{pos2} > {input_file2}  "], shell=True)  
        
        
        output_name=f"{input_name1}_{input_name2}"
        
        output_dir=f"./chr2moddotplot3/vs_all/{output_name}"
        output_file=f"{output_dir}/{output_name}.fa"
        subprocess.run([f"mkdir {output_dir}"], shell=True) 
        seq1='';seq2=''
        with open(input_file1,'r') as f:
            next(f)
            for line in f:
                seq1+=line.strip()
        with open(input_file2,'r') as f:
            next(f)
            for line in f:
                seq2+=line.strip()            
        with open(f"{output_file}",'w') as f:
            f.write(f">{input_name1}\n{seq1}\n>{input_name2}\n{seq2}\n")
                
        moddotplot_software='ulimit -u 100000 && source /home/lain/software/venv/ModDotPlot/venv/bin/activate &&  /home/lain/software/venv/ModDotPlot/venv/bin/moddotplot '
    
        if mode=='static':
            subprocess.run([f"{moddotplot_software} static  --fasta {output_file} --kmer 21 --output-dir {output_dir} --identity 80  --no-bed --no-hist --compare-only"], shell=True)     #--dpi 3000 --resolution 5000 
        elif mode=='interactive':
            # 1. Find processes occupying port 8050   # 2. Parse output to get PID (assuming output format matches expectations) # 3. Terminate all related processes
            result = subprocess.run(["lsof", "-i", ":8050"], capture_output=True, text=True)
            output = result.stdout
            pids = []
            for line in output.split('\n'):
                if 'PID' in line:  # Skip header
                    continue
                parts = line.split()
                if len(parts) > 1:
                    pids.append(parts[1])  # Assuming PID is the second column
            for pid in pids: 
                subprocess.run(["kill", "-9", pid])
            
            
            subprocess.run([f"{moddotplot_software} interactive --fasta {output_file} --identity 80 --compare-only --port 8050"], shell=True)   
        
                     
            
            
# Display all regions within the maximum centromere boundary, including VSat1-3            
if argv1=="select1":
    print('Display all regions within the maximum centromere boundary, including VSat1-3')
    if  os.path.exists('./chr2moddotplot3/select1')==False:
        subprocess.run(["mkdir ./chr2moddotplot3/select1"], shell=True)
    if  os.path.exists('./chr2moddotplot3/select1/tmp.fa')==True:    
        subprocess.run([f"rm  ./chr2moddotplot3/select1/tmp.fa ./chr2moddotplot3/select1/seq_all.fa ./chr2moddotplot3/select1/seq_all.fa.fai "], shell=True)    
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
 
 
    interval='S'*1000000
    chr_list=["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19"]
    #chr_list=["Chr1","Chr2"]
    seq_all_list=[]
    for one_chr in chr_list:
        pos_info=dict_sample_cen_pos['PN40024'][one_chr]
        if  one_chr not in ['Chr7','Chr12','Chr14','Chr17']:continue
        pos=f"{pos_info[0]}-{pos_info[1]}"
        subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample}/{sample}.fasta {one_chr}:{pos} > ./chr2moddotplot3/select1/tmp.fa  "], shell=True)
        seq_tmp=''
        with open('./chr2moddotplot3/select1/tmp.fa','r') as f:
            next(f)
            for line in f:
                seq_tmp+=line.strip()
        seq_all_list.append(seq_tmp)
    seq_all=interval.join(seq_all_list)
    with open('./chr2moddotplot3/select1/seq_all.fa','w') as f:
        f.write(f">seq_all_1\n{seq_all}\n>seq_all_2\n{seq_all}")
 
 
    # 1. Find processes occupying port 8050   # 2. Parse output to get PID (assuming output format matches expectations) # 3. Terminate all related processes
    result = subprocess.run(["lsof", "-i", ":8050"], capture_output=True, text=True)
    output = result.stdout
    pids = []
    for line in output.split('\n'):
        if 'PID' in line:  # Skip header
            continue
        parts = line.split()
        if len(parts) > 1:
            pids.append(parts[1])  # Assuming PID is the second column
    for pid in pids:
        subprocess.run(["kill", "-9", pid])
        
    moddotplot_software='ulimit -u 100000 && source /home/lain/software/venv/ModDotPlot/venv/bin/activate &&  /home/lain/software/venv/ModDotPlot/venv/bin/moddotplot '
    input_file='./chr2moddotplot3/select1/seq_all.fa'
    CMD=f"{moddotplot_software} interactive --fasta {input_file} --identity 80 --port 8050"
    print(CMD)
    subprocess.run([CMD], shell=True)    
    
if argv1=="select17":
    print('Display inverted repeat regions in the Chr18 centromere region')
    if  os.path.exists('./chr2moddotplot3/select17')==False:
        subprocess.run(["mkdir ./chr2moddotplot3/select17"], shell=True)
    if  os.path.exists('./chr2moddotplot3/select17/tmp.fa')==True:    
        subprocess.run([f"rm  ./chr2moddotplot3/select17/tmp.fa ./chr2moddotplot3/select17/seq_all.fa ./chr2moddotplot3/select17/seq_all.fa.fai "], shell=True)    

    target_list=[]
    
    target_list.append(['V041.hap2','Chr2',13300000,13850000])  ##No VSat1, only LIR
    target_list.append(['V049.hap2','Chr2',12570000,12810000])  ##
    target_list.append(['V108.hap2','Chr2',13380000,13640000])
    #target_list.append(['V020.hap1','Chr2',12180000,12240000])  ##Too small, and appears to be a VSat2-related inversion
    ###
    #target_list.append(['V053.hap1','Chr18',15800000,16190000])   ####C2
    #target_list.append(['V100.hap1','Chr18',15080000,15330000])   ##C1
    #target_list.append(['V106.hap2','Chr18',15580000,16110000])##C1#
    #target_list.append(['PN40024_hap1','Chr18',16320000,16680000])   ##C3
    interval='S'*100000
    
    #chr_list=["Chr1","Chr2"]
    seq_all_list=[]
    i=0
    for one_target in target_list:
        i+=1
        sample,one_chr,start,end=one_target
        pos=f"{start}-{end}"
        
        subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample}/{sample}.fasta {one_chr}:{pos} > ./chr2moddotplot3/select17/tmp.fa"], shell=True)
        seq_tmp=''
        with open('./chr2moddotplot3/select17/tmp.fa','r') as f:
            next(f)
            for line in f:
                seq_tmp+=line.strip()
        seq_all_list.append(seq_tmp)
        seq_all_tmp=interval.join(seq_all_list)
        print(f'{i}\t{one_chr}:{start}-{end}\tlength:{int(end)-int(start)}\t\t\t{len(seq_all_tmp)}')
    seq_all=interval.join(seq_all_list)
    with open('./chr2moddotplot3/select17/seq_all.fa','w') as f:
        f.write(f">seq_all_1\n{seq_all}\n>seq_all_2\n{seq_all}")
    print(f'Synthesized sequence length (including spacer sequences): {len(seq_all)}')
 
    
    moddotplot_software='ulimit -u 100000 && source /home/lain/software/venv/ModDotPlot/venv/bin/activate &&  /home/lain/software/venv/ModDotPlot/venv/bin/moddotplot '
    input_file='./chr2moddotplot3/select17/seq_all.fa'
    if 1==1:
        # 1. Find processes occupying port 8050   # 2. Parse output to get PID (assuming output format matches expectations) # 3. Terminate all related processes
        result = subprocess.run(["lsof", "-i", ":8050"], capture_output=True, text=True)
        output = result.stdout
        pids = []
        for line in output.split('\n'):
            if 'PID' in line:  # Skip header
                continue
            parts = line.split()
            if len(parts) > 1:
                pids.append(parts[1])  # Assuming PID is the second column
        for pid in pids:
            subprocess.run(["kill", "-9", pid])
            
        
        input_file='./chr2moddotplot3/select17/seq_all.fa'
        CMD=f"{moddotplot_software} interactive --fasta {input_file} --identity 80 --port 8050"
        print(CMD)
        subprocess.run([CMD], shell=True)     
    else:      
        os.chdir(f'./chr2moddotplot3/select17/')
        subprocess.run([f"{moddotplot_software} static  --fasta seq_all.fa --kmer 21 --output-dir output_dir --identity 80  --dpi 3000 --resolution 5000"], shell=True)  
        
if argv1=="Chrall":   
    print ( "Output all core regions ")
    Chromosome_list=['Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10','Chr11','Chr12','Chr13','Chr14','Chr15','Chr16','Chr17','Chr18','Chr19']
    
    
    input_list=[]
    with open("/home/lain/aaa_data/run0/samples_satellite/2_good_regions_main",'r') as f:
        next(f)
        for line in f:
            eachline_arr=line.strip().split()
            sample,chromosome,start,end,length=eachline_arr

            #if chromosome!=Chromosome:continue
            input_list.append(eachline_arr)
            
            
            
    
    def run_step(one_input):
        sample,chromosome,start,end,length=one_input
        output_dir=f"./chr2moddotplot3/Chrall/{chromosome}/{sample}"
        subprocess.run([f'mkdir -p {output_dir}'], shell=True)  
        start_1M=int(start)-1000000
        end_1M=int(end)+1000000
        pos=f"{start_1M}-{end_1M}"
        input_file=f"{output_dir}/input.fa"
        subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample}/{sample}.fasta {chromosome}:{pos} > {input_file}  "], shell=True) 
        if  os.path.exists(f'{output_dir}/{chromosome}:{start_1M}-{end_1M}_FULL.pdf')==True: return False  
        moddotplot_software='ulimit -u 100000 && source /home/lain/software/venv/ModDotPlot/venv/bin/activate &&  /home/lain/software/venv/ModDotPlot/venv/bin/moddotplot '
        subprocess.run([f"{moddotplot_software} static  --fasta {input_file} --kmer 21 --palette Blues_9 --output-dir {output_dir} --identity 80 --dpi 300 --resolution 1000  "], shell=True)  
        
                    ###View the colors used
                   ##from palettable.colorbrewer import sequential, diverging, qualitative
                    # Print all sequential palettes
                    ##print("Sequential Palettes:")
                    ##print([p for p in dir(sequential) if not p.startswith("_")])
                    
                    ###
    with Pool(processes=4) as pool:
        # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step, input_list), start=1):
            # Process results here, e.g., store or print
            progress = (i / len(input_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush() 
        
if argv1=="Chrall2one":   
    print ( "Output all core regions ")
    Chromosome_list=['Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10','Chr11','Chr12','Chr13','Chr14','Chr15','Chr16','Chr17','Chr18','Chr19']
    one=sys.argv[2]
    one_region=sys.argv[3]
    #ShineMuscat_hap1:Chr3:13908633-14381843
    input_list=[]
    with open("/home/lain/aaa_data/run0/samples_satellite/2_good_regions_main",'r') as f:
        next(f)
        for line in f:
            eachline_arr=line.strip().split()
            sample,chromosome,start,end,length=eachline_arr
            if one!='all':
                if chromosome!=one:continue
            #if chromosome!=Chromosome:continue
            input_list.append(eachline_arr)
    
    TARGET_sample,TARGET_chr,TARGET_pos=one_region.split(':')        
    subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{TARGET_sample}/{TARGET_sample}.fasta {TARGET_chr}:{TARGET_pos} > ./chr2moddotplot3/Chrall2one/tmp.fa  "], shell=True)         
   
    
    def run_step(one_input):
        sample,chromosome,start,end,length=one_input
        output_dir=f"./chr2moddotplot3/Chrall2one/{chromosome}/{sample}"
        subprocess.run([f'mkdir -p {output_dir}'], shell=True)  
        start_1M=int(start)-1000000
        end_1M=int(end)+1000000
        pos=f"{start_1M}-{end_1M}"
        input_file=f"{output_dir}/input.fa"
        subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample}/{sample}.fasta {chromosome}:{pos} > {input_file}  "], shell=True) 
        subprocess.run([f"cat ./chr2moddotplot3/Chrall2one/tmp.fa >> {input_file}  "], shell=True) 
        if  os.path.exists(f'{output_dir}/{chromosome}:{start_1M}-{end_1M}_FULL.pdf')==True: return False  
        moddotplot_software='ulimit -u 100000 && source /home/lain/software/venv/ModDotPlot/venv/bin/activate &&  /home/lain/software/venv/ModDotPlot/venv/bin/moddotplot '
        subprocess.run([f"{moddotplot_software} static  --fasta {input_file} -w 5000 --compare-only  --kmer 21 --palette Blues_9 --output-dir {output_dir} --identity 80 --dpi 300 "], shell=True)  
        
                    ###View the colors used
                   ##from palettable.colorbrewer import sequential, diverging, qualitative
                    # Print all sequential palettes
                    ##print("Sequential Palettes:")
                    ##print([p for p in dir(sequential) if not p.startswith("_")])
                    
                    ###
    with Pool(processes=4) as pool:
        # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step, input_list), start=1):
            # Process results here, e.g., store or print
            progress = (i / len(input_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush()         
        
        
        
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))      
time_end=timeit.default_timer()
print('All the running time: %.0f Seconds'%(time_end-time_start))






            
            
            
            