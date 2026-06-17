#!/usr/bin/python
# -*- coding: UTF-8 -*-
import sys
argvs=sys.argv
argv_len=len(sys.argv)
if argv_len < 2 or argvs[1] in ("h", "-h", "help", "-help") or not argvs[1].startswith("step") :
    print ("chr2blast.py-----help:")
    print ("Roughly estimate the satellite positions and ranges. For precise identification when the monomer length is <50 bp, chr_mafft.py should be used.")
    print ("")
    print ("Usage:")
    print ("chr2blast.py step0 -i  <genome_folder>")
    print ("chr2blast.py step1  all/cen_66...        ##")
    print ("chr2blast.py step1s   Automatically summarize step1 results")
    
    print ("chr2blast.py step_generate_oneseq   Extract block regions to form a large sequence for NTRprism monomer analysis")
    
   
    print ("-thread \t\tNumber of threads (default: 70), used for multiprocessing in some steps")
    print ("-i      \t\tInput FASTA file (required for step0)")
    print("\nManually modify the script, or specify a name and monomer sequence:")
    print ("-cenid  name   -mononmer unit_sequence")
    print("unit_sequence supports multiple reference sequences, using || as the separator.")

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

### Parse parameters
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
else:thread=30

# Software requirements
blast="blast"


def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_complement_seq = ''.join(complement[base] for base in reversed(sequence))
    return reverse_complement_seq
    
    
    
if  os.path.exists('./chr2blast')==False:
    subprocess.run(["mkdir chr2blast"], shell=True)

#step0
if argv1=="step0":
    print('Building BLAST index, this step is fast even if completely redone')
    if  os.path.exists('./chr2blast/0_samples')==True:
        subprocess.run(["rm -r ./chr2blast/0_index"], shell=True)  
    subprocess.run(["mkdir ./chr2blast/0_index"], shell=True)  
    
    if "i"  in args_dict:
        input_dictionary = args_dict["i"]
        with open('./chr2blast/0_index/sample_source','w') as f:
            f.write(input_dictionary)
    else:
        print("Missing input fasta file");sys.exit()
    
    dir_fasta_file_num=0
    dir_file_name_list=[]
    files=os.listdir(input_dictionary)
    with open('./chr2blast/0_index/sample_list','w') as f:
        for one in files:
            if     one.endswith('.fasta'):one_name=one[:-6]
            elif   one.endswith('.fa'):one_name=one[:-3]
            elif   one.endswith('.fna'):one_name=one[:-4] 
            else :continue
            dir_file_name_list.append([input_dictionary,one,one_name]);      dir_fasta_file_num+=1
            f.write(one_name+'\n')
    
    print('Number of genome files: '+str(dir_fasta_file_num))

    def run_step0(dir_file):
        input_dir=dir_file[0]
        file_name=dir_file[1]
        sample_name=dir_file[2]
        input_file=input_dir+'/'+file_name
        output_dir="./chr2blast/0_index/"+sample_name+"/"
        if  os.path.exists(output_dir)==True:return False
        output_file=output_dir+sample_name
        subprocess.run(["mkdir "+output_dir], shell=True)  
        ### Limit by length 1,000,000
        with open(f'{output_dir}/{sample_name}.fa','w') as f2:
            one_seq=''
            with open(input_file,'r') as f:
                for line in f.readlines():
                    eachline_arr=line.strip().split(' ')
                    eachline=eachline_arr[0]
                    if eachline[0]=='>':
                        ##
                        if len(one_seq)>1000000:
                            f2.write(f'>{one_id}\n{one_seq}\n')
                        ##
                        one_id=eachline[1:]
                        one_seq=''
                    else:
                        one_seq+=eachline
            if len(one_seq)>1000000:
                f2.write(f'>{one_id}\n{one_seq}\n')                
        ###        
        subprocess.run([f"makeblastdb -dbtype nucl -in {output_dir}/{sample_name}.fa -out {output_file}"], shell=True)  
    
    # Assign tasks to processes in the process pool    
    with Pool(processes=thread) as pool:
        pool.map(run_step0, dir_file_name_list)      
    

if argv1=="stepall" or argv1=="step1":
    if  os.path.exists('./chr2blast/1_blastn')==False:
        subprocess.run(["mkdir ./chr2blast/1_blastn"], shell=True)  
    argv2=sys.argv[2]
    info_list=[]
    
    
    if argv2=="cen_6" :
        cenid="cen_6"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer='GGTGACAGTAACGGTTATGGTAAC'   #GGTAAC
        info_list.append([cenid,output_dir,monomer])
        
    if argv2=="cen_25" :
        cenid="cen_25"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer='GCGTTCAGAGCCACCATTTTTATAG'
        info_list.append([cenid,output_dir,monomer])         

        
    if argv2=="cen_26" :
        cenid="cen_26"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer='AAGTGGCTATAAACGCCTAACTCAGA'
        info_list.append([cenid,output_dir,monomer])        
        
        #28BP#ACTCGCACGGATACCACCATTTTTCGGT
        
    if argv2=="cen_31" :
        cenid="cen_31"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer='AGCCATCATCATTCTCAGTTTTGACGTTCAG'
        info_list.append([cenid,output_dir,monomer])            
        
    if argv2=="cen_43" :
        cenid="cen_43"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer='CAAGCAACGCCCCGTGGCACCTTGGCTTGGCTCGCGGGTGGGC'
        info_list.append([cenid,output_dir,monomer])
        
    '''cen53 belongs to cen103    
    if argv2=="cen_53" or argv2=="all":
        cenid="cen_53"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer='TTCGAGTCAGGACCGCACCTATATCGATCGGTCATCATCCTATCAGTTTGAGT'
        info_list.append([cenid,output_dir,monomer])        
    '''    

    if argv2=="cen_60" :
        cenid="cen_60"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer='AAGTTTCTTCTTTGTTTCTTTTGTTTCTCATCATTGGGGCACCTTCCTAAGTTTTATTTT'
        info_list.append([cenid,output_dir,monomer])        
        
        
    if argv2=="cen_66" or argv2=="all":    #VSat2
        cenid="cen_66"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer='CCTGAGTTGTTCATAGCATCTAGAGCCTTGAGCTCATGACCTAGAGTCGTTCATGTCATCCACAAC'
        info_list.append([cenid,output_dir,monomer])

    if argv2=="cen_68":                   #VSat8
        cenid="cen_68"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer='GAGTAGTTTTTCTCATTTATGACCTAGAGTCATTTATTTTCCATTTAGAGCCCTAAGCTCACGACCCA'
        info_list.append([cenid,output_dir,monomer])        


    if argv2=="cen_101" :
        cenid="cen_101"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer='AATCACATCGTTCTCCTGTCTGCTATTCATGATACATAGTATCCTGTTCGATCAGTCATGATAGTTCAGTTTCGTTTTTGGCATAGACCGCACTTGTATAA'
        info_list.append([cenid,output_dir,monomer])        

    if argv2=="cen_103" or argv2=="all":   #VSat3
        cenid="cen_103"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="TCATGACCGATCGGATCGGGTGCGGTCTATGATGAAAACCAGATAGGACAACGACGTGACCGATCGTATAGGTGTGGTCTACGCCGAAAATGAAACTGAACTG"
        info_list.append([cenid,output_dir,monomer])
        
    #The localization of VSat1 does not use this script; instead, chr_mafft.py is used for greater precision.    
    if argv2=="VSat1" :           ##VSat1. unit_sequence supports multiple reference sequences, using || as the separator.
        cenid="VSat1"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="ACTCGCACGGATTCTACCTTTTCCCGGTACTCGCCCGGATTCTACCATTTTTCGGTACTCACACTGATTCTACCCTTTTCCGGTACTCGCCCGGACAGTTTTTCGGT||ACTCGCACGGATTCTACCTTTTCCCGGTACTCGCCCGGATTCTACCTTTTCCCGGTACTCGCCCGGATTCTACCATTTTTCGGTACTCACACTGATTCTACCCTTTTCCGGTACTCGTCCGGACAGTTTTTCGAT||ACTCGCACGGATTCTACCATTTTCCGGTACTCACACTGATTCTACCATTTTTCGGTACTCGCCCGGACAGTTTTTCGAT"
        info_list.append([cenid,output_dir,monomer])        
         
    if argv2=="cen_145" :
        cenid="cen_145"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="TGGGTTTAGCGTTTATAGAGTTCATGGTCTTTATGGTTTAGGGTTAGGTTTTATGGTTTAGGGTTTAAGGTGTTAGGGTTGATGGTATCGGGTTTCAGGGTATAGGGCTTCCAGATTTGGGTTTAAGGTTTTGGGTTTTAAGGTT"   
        info_list.append([cenid,output_dir,monomer])
    if argv2=="cen_168" :
        cenid="cen_168"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="AATTTTGTGAAGGATGTACCAAGGGTTCCAAAGTGCCATCTGAAGTGGGGATCGACCCTAAATCACTTTGGTATGGGTTTCTAAATGAAATTATGACCCATATATGCTAATATGGGGGTCCCGAAAGCACTAGGCTCGGATGTACCAAATTTGAAGAAAGATGTACCA"   
        info_list.append([cenid,output_dir,monomer])
        

    if argv2=="cen_187":
        cenid="cen_187"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="AAGTGGAGAATTGGGAGGGTACGAAGAATGTGAGGATTCCTCACATTCTTCGTACCCTTGTGAATTTTTTAGTATTTTTGGATTTCAATTTTTCCGGATGTCATTTGACCACCTTCAAAGTGTAACGAACCAACTTAAACAACATCAAAGTCATATGATGCCCCAAAAATTTTTAAAAAGTCCAAAG"   
        info_list.append([cenid,output_dir,monomer])
    if argv2=="cen_191" or argv2=="all":       #VSat4
        cenid="cen_191"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="CCGGCTGAACCGGTTTGGAACCGGTTCAACCGGTTCAGCACCATAGCTGGCTGCCTCTGTCAGCCATGAGGAACACCTGCCTTCCTTTGTCCGGTTCTCACTCATCCGAGCTCGGAATTGAGAACCGCTTCTTTTTTTGTCTTCCTTAATCACTTAGGAACTTTCTGACCAAATTAGGGATTTTTTACCAA"
        info_list.append([cenid,output_dir,monomer])
    if argv2=="cen_355" or argv2=="all":      #VSat5
        cenid="cen_355"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="AAGGGTGGAATTTTTCAGGTTTGCCGGGAACAAGAGAAGTTTATTTAGAAATAGTTGTTGGGTAAGTTTCTTCCGTGTTCCCTATTCTCTCAAATTTCTCATTGTTCTTCATTTTTCTCTTGTATGTTTGCGTGTTATGTGGTAAATGTATGCGTGTACTAGCAAAATTATTTTTCATTATGTTGATTTCTTTGTTTTCTTGTGATCTGGTTGTTTCGTTGTGTATTACTCATTTATTTATTTTTGTTTGAAAATACTTTGATGTCCATTCAACAACATGGTGCCTTAATAGGTAATTTTCTTATTCTCGTTTTTACAATTGCTCTCTACTTGTTTTCCTATCTTCTATAAGCTT"                
        info_list.append([cenid,output_dir,monomer])
    '''383 same as 370
    if argv2=="cen_370" or argv2=="all":
        cenid="cen_370"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="CATTGGGACATGATCCTTTGGTTTCGCATCGTTTGGACACCATCGAGGTTTGATGGATGTTTATTCTGTCAAAATTTTGGACGACATTGTAGGTTCTTTCTAAGGTCATGCATTTTTTACTTTGGTTTGTCATCATGGGGACGCCCCCTAGGCGCTATCCAAGCTTAATTGTTTGTTCCTTTAATTTCTCACCGTTTAGACACCATCCTAGGTTTTGCGTATGTTTATACTTCATCGATTTTGTTTCCTATCATTGGAACACCATCTAGGTTTATTCTGAGTTCATTTGTCTATTCCAGTCGTTTTTTATCATTGGGACAGCTTTCTAGGTTCTTTCTACGCTCATTCCAGTTTTCCATTTGTAACTAAT"
        info_list.append([cenid,output_dir,monomer])   '''   
    
    if argv2=="cen_383" or argv2=="all":     #VSat6
        cenid="cen_383"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="CATTGGGACATGATCCTTTGGTTTCTCATCGTTTGGACACCATCTAGGTTTTATGGTAGTTTCTTCTTTCTCCATTCGGTTTCTCATTTTGGGAGGACATTGCAGGTTTTTTCTAAGGTCTTGCATTCTTTACTTTGGTTTGTCATCATTGGGACTCCCCCTAGGCCCTATCTAAGCTTAACAGTTTGTTCCTTTAATTTCTCATCGTTTAGGCACCATCCTAGGTTCTGTCTATGTTCATACTTCATCCATTTTGTTTCCTATGATTGGAACACCATCTAGGTTTAAACTATGTTCATTTGTCTTTTCCACTCGGTTTTCATCATTGGGACAGCTTTCTAGGTTCTATCTATGCTGATTCCTTCTTTCCATTTGTAACTAAT"
        info_list.append([cenid,output_dir,monomer]) 
    
    if argv2=="cen_650":  #VSat5b
        cenid="cen_650"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="TCGTTCGTTATGTCTGTTATGTTTGGGTATTATTTTTTGACTTTTTGTAAGCGTTATAGAATTTATTTAGAATATTTTTATTTCTGTGTTTGGTTGTGGGTTTTATGTGTGACTTTTTATTTGTCATTTATTTGATTGTTTGTGACAAGTATTTGATTTCCTTTGAATTGTGTTTTTCCATGGGATGTTATTTTGATTTTCTTGTCTTTTCAACTGATCTCGGTGTCTTGTCCCACCATTAATAAGCATATAAGTGGATTTTTCAGGATTCTCGGGAGCAGGAGAAGTGTGTGAGGAAAGTCGTCTTAGATAAGTTTTCTTCGTATTCCCTGTTCTTTCATATACCGTAGTGTTCAATATATCTTTTATGTTTGGGTATGATTTGTTGACTTGATCTATACATTATTGAATTCATTTCGAATATGTTCATTATGTTTTGGTTTGTGCGCTTTACGTTTCATTTTTTATTTCTCATTTGTTTCTTTGTGGGTTAGAATAGTGTAATTTCCATTTCTGGTGTTTTCATCTGATCTCGGCTTCTTTGGCAAGCATTTCTAAGCATATCACTGGATTTTTCAGGTTTCTCGGAAGCAAGAGAAGTTTAAGAGGAAGTACTTCTTATCTAAGTTTACTTCGTATTCCCTTTTGTC"
        info_list.append([cenid,output_dir,monomer])           
        
    
    if argv2=="cen_677":  ##same as 354
        cenid="cen_677"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="AAGGATGGAATTTTTCAGGTTTCTCAGGGACAAGAGAAGATTCGTCGGAGATACCTCTGTAGTATGTTTATTTCGTATTCCCTATTATCCTCCATTTCGCACTGTTCTTCATTCTTCTGTTTTATCTTTGGGTGTTGTTTGTTGACTGCAAATATGTATTAGAGAATTAATTTTGAATCCTTTGATTTATTTGTTGGTTTCTGAGCTTCATGCTGCGTTTTCTAATACTCGTTTCATTCGTTTCGTATCAGAATAATATTACTTGCTTTTAATGAGATTGTGCCATAAGTTGTAATTTTCTGATTCTCGTTCTTACATTTGATCTCTGTTCTTTTCTATACCTTCTATAGGCGTAAGAGTGGAATTTCAAGGTTTATCGGGGACAAGAGAAGTTTGATTGGAAATAGGTCTTAGGTAAGATTTCTCCGTATTCCCTATTCTTTAGAATTTCTCGTTGTACTTTAATTTTCTCTTGTATGCTTGGGTGTTATATGTCGAATGCATGTGTGTAGTAGGAAAATTAATTTTCAATATATTGATTTCTTTGTTTTCTTCTAATGTTGATGTTCTGTTGTCTATTTTAATTCCATTGAACGACATTGTTCATTAAGAGGTAATTTTCTTGTTCTTGTTTTTTACGTTTGATCTCTACTCCCTTCACTACCATCTTTAAATGT"
        info_list.append([cenid,output_dir,monomer])         
        
        

    if argv2=="cen_783" :
        cenid="cen_783"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="ACTGCGGAGTTCTGATGGGATCCGGTGCATTAGTGCTGGTATGATCGCACCCGTTGTGTTATTGCCAACGAGTCCTCTTATGTGGTATGATTGAGGTGCAACTTCTAGTGTGTTGGATAGAAGGAAATAGGAGGAAAAGCTTGTTTTGACAGATGTATTGTTATGCGAGTATCTGGGGATGTTCTTCCAGATCGGACGGTGTGCCCCTCTTAAAAATTTGCACATTCTTTCGATCCGTACGGTTGGATTCGCGTTTCGGTGCATGGATTGGGGTTTCATCGTCACTTTAAAGGGTGGATGGGGGAAATGTGGAAACGAGAGTTTACCAGGTGGTTTTTGGGCCCTTCGGGTAGAGTTTGGGGACCCTTTTTTTTTTCAAACCATGAGTCTCGACGCGAGGATCGCGATGGTGGGCTCGGATCGTCCTTTCGACGGACGGATGGGGAGACTAGGGAGGAAAGGCGTTTTGGGCCATAACTTGGGCCTCGGAGCTCGAAAAAGGTTAGCGTCTGCGCCGATCGTCCTAAAATGGGGGAAACTTTCACGTGTGAGAAAACGACAAAAATCGGTGACGATCGATTTCGACGGCTGATAGCTCACGATCGGGTGGTCCGATCGACGATCGGCCGGTCCGGATGGAAATGGGAAGGGCGGAACCTATTGGCACGACAAGTGCGGACATAAAATCGTAGAAATTCCGTCGTAAAAGGGAAAAAGGGGTGCAACACGAGGACTTCCCAGGAGGTCACCCATCCTAGTACTACTCTCGCCCAAGCACGCTTA"
        info_list.append([cenid,output_dir,monomer])         

    if argv2=="LIRSat1" or argv2=="all" :######cen_1305    
        cenid="LIRSat1"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="TTGGACCATAATTCTTGGTTCTATCTTTGTTTATTCTTAGTTCCTTTAGTTACCCATCATTGGGTTGCTATCCTAGTTTTAATATGAGTTTATTCGTTTGTTCCCTCTTATTTCTCATCACTAGGATGTCATCCTAGGTTCTATATCCAATTTGAGTTAATTTTGTGTTTGATCGTCTACATACCATCCTAGGTTCTACCTAAGATTAGTTGATTGTTTGTTATGTTACTCGTCATTGTGATATCATCCTAATTTCTATCTAGGTTTGTTCTTAGTGCCTTTGGGTTCAAACATTAGGATCTAATCCTAGGTTTGGCCTATGTTCATTGTTTGTTCCTTCGATTTCCCATCATTGTGAACTATCTAAGCTTATTGATATGTTCCTTTGGTTTCTGATAATTAGGACGCCATCCTATGAGTTTTTCTTAAGCGTATGCACTATTCCTTGGGTACGTTCAGGACGTCCTCCTTCCTTGAGTTGTTTTTTGGCATTAGTTTGTTTTGGGCATTCCCCTTCATTAGGTTTATTTAGGGATCCTCTTTCCTTGAGTGTTATAAGGGTATTCCCATTCCACTACTTTTTTAGGGCATTCCTGTTCATTGGGCTATCAAACTACCTTGATTTACTTTAGTCCATCATGCTTCCTTCAGTTTCTTAGGTCATCCATCTTCATGGAGATTTTTAATTCATTCCCTTTCCTAGAGTTTGTTTACACCAACCAACTTCCTTCAATTCTTTTAGGGCATACCCTTTCCTAGGGATTACTTAGCAGACTCCCACTTCCTTTAGTGTAGTAGGACATCCCCATCTAGTAGTTTGTTTAGCGCAACTCCATTTCTTGAGGTTTTTTTAGCTCAACCGTTTTTGTTGAGATTGGGCATCCCCCTTGGTAGAGTTTGTTTGGGGCATGCCCTTTCTATGAGTCTTATTAGTGAATCCCCCTTCTTCATGTTCATTTACAGCGTTCTACTTCCATTGGCTTTTTAAGGGCATCCCACTTACTCGAATTTGCTTAGAGCATCCCCCTTCCTCGACTTTGTTTGGAGAATTCCCTTCATTAAGTTTATTAAGGGCATCCCACTTCCTTAAGTTGGTTTAGGGTATCGTTCGTCCTTCACCCCATAATCAGGATGTCATCCTATGTTCTATCTAAACATAATCATGTGTTCCCGAGATTCCTCATCATTGGAACACCATTCATAGTTCTATTCATTTTCATTCTTAGTTTATTTAGTTTCACGTTGTTGAGATGCCATCCTAGGTTCTATCTAAGCTAATTCATTTGTTCTTTTGATTGCTCATCA"
        info_list.append([cenid,output_dir,monomer]) 
    
    if argv2=="LIRSat2" or argv2=="all" :##########cen_185 
        cenid="LIRSat2"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="CTTAGTTTAGTTTTTAAGCATTATTAGTTTGCTTACTTTCATGATATTATACTAAGTTTAGTTTCTTTCAAAAAGTCAGTTAATCTCATGTTTAGTTTAGTTTTGAAAAATTATTATTTTGGTTTCTTTGTTGATATTATGCTTACTTTAGTTTTTGTTTTTTGTTTTTTTTAAAAAAATCATTAGTTTGGTCATTCCTTTCACCTCATGTTTAGTTTAGTTTATAAAAATTATTAGTTTGGTTGCTTTGTTGATATTATACTTAGCTTAGTTTCTTTTAAAAATCATTAGTTATGTGGTTATTCCATTAATCTCATGTTTAGTTTAGTTTATTTTATAAAAAATTATTAATTTGGGTGCTTTGTTGATATTGTGCTTAGTTTAGT"
        info_list.append([cenid,output_dir,monomer])   
    
    if argv2=="LIRSat3" or argv2=="all" :########cen_650    
        cenid="LIRSat3"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="TCCAAGTCAATTTGTTAATAAATAAATCATTTTTCAAAATAAAAAGCTAGTTTGCCATCGAGTTGAAGCGCATGAGCCGAAATGCGGGGTCCACAAGTATTGTAATTGATATTTAAGCTAAAAATTTAATCTCCATGTACCCAAAGATTGGGTTGTGCCAAATTTATTTTTAACATGAATCCGGCCCATTGGCTTGCCCAACTTAAATTTTGGGCCTCACCTCAAATAAAAGCGAAGCTCGAACCTTAACCACTTAAAAATATTTGGGCTACGTTGGGGTAAGAAATAAATTCCTCATAATAGCCATCCAACCAAGGGAATGAAAGTACCCAATGGGCTTTCATTTTAGTTTGGCTTTTTCGGTCTTTTGTTTGTTGGAATGGAATTACCCAATGGGCTTGTATTTTTAGTTGGACCCTCCAATTATTTTTTTAACCTCATCCTAAACACCCCAAGGACCCAAAAGAAAAAAAAAATTGTTTTCTTTGAGAAGTATATCCCAACGCCTCCTAAATAAATTGTTAATATTTTACTATAAACTTTTATAAATAGAAATTTCAATATTATAATAAATAAATAAGGAGTCACACTTAGGGTTTTCTTTCTTATTTTGTTTACCCTTTAAAAATAAAACAAAAATAAGTGACGAC"
        info_list.append([cenid,output_dir,monomer])        
      

    
    if argv2=='cen_240' or argv2=="all":   ###The second most abundant satellite in the LIR of Chr17
        cenid="cen_240"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="ATGAGAAACCAAAGGAACAAACAAATAAGGTTATATAGAACATAAGATGACGTCCTAATGATGAGATGAGAAATAAAAGGACCAAAGAGTAAACTTAGATAGAACCTATGATGTTGTCCCAATGAAAAGAATCCAAAGGACCAAATGATTAAGCTTAGATAGAATCGAGGATGGTGTCCCCATGATGAGAAACTAAATGAACAAAGAAGATACTTAGAACTTAGGATGATGTCCCATTGATGAGAAACTAAATGAACAAAGAAGATACTTAGAACTTAGGATGATGTCCCATTGATGAGAAACCAAAGGAACAAATGATTAAGCTTAGAGAGAACCTAGAATGGTATCCCATTGATGAGAAATGAAATGAAAAAGAATGAAGTAGGAGAAAACCTAGGGTGGTGTCCTAATGATGAGAAACTAAATGAAGAAAGAAGAAACTTAGATAGAGCCTAGGATGGTGTCCTAAAGATGAGAAACCAAAGGAACAAACAATTAAACTTAGATAGAACCTAGGACGGCATTTCAATGATAAGAAACCAAGGGAAAAAAAACGAATATTCTTAGATAGTACCTAAGAGTGTGTCCCAATGATGATAAATGAATGGAACTAAGAATAAACTTATATAGAAGTTATGATGGAGTTCCAATGATGAGAAAACAAAGGAACAAAGAATAAACTTAAATAGAACCTAGAATTGTAAGCCAATG"
        info_list.append([cenid,output_dir,monomer])   
        
    if argv2=='cen_TekayLTR' or argv2=="all":   ##
        cenid="cen_TekayLTR"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="TGTAACACCTAGATTATTTTAGTACTTAACTTAAGCTTGCTTAATTAGATATTAAAACTAGTTTAGTGCTAATCATGTTTAATTATGAATTAAACTTTGATTAAGGTTAAGTGGGATTAGTTAATGATCTAAGCATGCTTATTAGGTTCTTAGGGACTAATTAAGGTTATGAAAGCCTAAAGGACTAATGGGCAATTATGGGTTTCATTTTTGATAACTTGAAGGACTAAAGTGCAAATTTAGAAACTTGAGTTAGTGGAATGCCACCTCCTACTTGGTGGATTGGTGGCAGCCATGTGGGCTGCCACCTCATGCTTGGGTGCTTGGAGACCCCTTTATAAGGGGGGCTGCAGCTCGTTTTGCACCATTTCACCTGCAGCAACTTGGGTTAGAGAGAGAGACTAGAGAGAGAAAGTGTAGATTTGAGGTAAGCAAGTGTTTAATTTTGTAATTTTCGTTTATTTTGGTTCCTAATGGTATGTTGAGAAGGATTAATGGTAAAAATTGTGTAAAAGTTGAGTATAATTAGCTATAAATCTGTTTTAGTTAAAAATCACAATTTATTAAAGATTAAAGTATTTTAAGTTAAGTATTTTATTAAGGGTAAGATCAGCATAAATCTTTGTTTAATTTGGTTTGGTTGGAAAATAGATTAGTAAACATGTGATGAATTGGGTTATTTGGACACTTAGGATTTGTTCATAGGTTAGGCTTTAAGAAAATCAAAAGAATTAGTGTTTGGTAATTAATTAGGGAAATTAATATGCATTGTAAATATTGTTTAATTCAATTCATATTTGATTTGAAAACATGTAGGATATAGATTAGAGAATGTGAGAATTAAAAATTTGCAACAATAAACTAAAACCTAAATTATGTTGTTTCATTTCAGTTCCTCGTAATAAGCAAAGATCGCCTACTTAGAAGAAACACAAAAAGGAAGTCGGTGTAAGGCAGGGAATTTTATGCTATTCTATGGTGTATCTGATTTCTTTCCTTGAATTATTTGTTGGAATTTTATGTTGTTTCTAGGATTTCATAAAATCTTTTAAAGTGTCATTGCATCACGGTTTATTGCATTGAAATGTTACCATGATGTTGGAACGTATATTGGTATGAAGATTCGTGGTTTTGTGTTGCCGAAATGGTTTATGGTGTGATTGCATTGCAAGAAGGATATGGAATTATCACAGATCATATTCAAAGAATTTATTGGTGGCGAAGTGCTTATGTGATGCTACAGGCCTTATCATCTGATATATTGTTTATGTGGGACCGTGGGTCCTTGGGTGAAAGTCCCTAAAGCCCTCGAGGAAGACACTCCGATGTGGGTGTATGGGGTTGGTCCTGCCCCTGGGTTTTAGTCCCTAAAGACACGGGGTTGGTCCTGCCCTTGGGTATGAGTCCCAAAAGACACGGGGTTGGTCCTGCCCTTGGGTATGAGTCCCTAAAGACACGGGGTATGACTTGCCCTTGGGTGATGTCCCAGAATAGTCATTATTATATTGATCAAGTATCTTGTGGAGCATTGATATCTGCTGTGTACATTGTCGGGAGTCAGCAATTTATCAGTTGTGAAAGGATGGATGGGGAAAAGCAATGATGAAACCAACAAACACATGCATACATGTTTGACATGATTACACATTTGTTAATGATTATGAAATATTTATCATGCATGTTATTAAACGTTTGATTCAAGGTTTTTAAGGGTATACTTAGGATGGTTATAAACTTTCCTACTGAGTTGTGAACTCACCCTATCCCTTCCACCTCTAGATGCAGGTCAGAAGTCCTATGCAGGAAAGAATGCTTGAGCAATGCTTTACTGTGATTGGATGCTGTTTGCTCATATTGAAAGTCTTTTGAGGCCATGCCTTTTGTATAAAGACTAAACCTTTTGTTGGAATGATGTAATGTACATATTTGTCAGATACACTTGTAGTATGGGCTACCCGTAGTAAACAATGGGGTTTGGTATATGTAAATTAATGTTGTACAAGTTTCTTTTGGGAAACAAAACATATTTTGTTCACTTCATAGTCACAATGTTTAATACAGGATGTACTCTTTATAACAGGTTTTCCAAATCCTCTTTTTGAACAAATTCAAGTAGGAACTCAACATTTAACTTTTGGAAAGCACTTATACTGTATTGAGTATCAATTGCTTAGTTGAAGAATAAATAATAATAATAATAAGAAAAAAAAAATAATGGGGTGTGACA"
        info_list.append([cenid,output_dir,monomer])          
        
    if argv2=='cen_CRMLTR' or argv2=="all":  
        cenid="cen_CRMLTR"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'        
        monomer="TGATGGGGGAAAGCCTTTTAAAGTTTATTCAAGGAAGAAAGGGAAGACCTAAGTAAACTTTAAGTAGGATTAACATTAATTAGAATTGAATAGGGATTGGTATTTGTTGGAGTCTAATTAGGATTAGCTAGTTGAGTTATAATATAATTAGAATTTGAGTCATGATAGGTTATTAGAGTCCTAGTAGACTTTGGATTTCTTAGAGAAGCCTATAAATAGGCTAATCAATGTAAATCAAAGGAAGGAATTTTGATGAATAATATTATACTTTCTTTCATTGCAAGGTTGCAACCCTCAATGGTGAGACTCCATTGATTTTCTTCTAGGTGAGACTCCTAGAAGGCCTTAGTGAGACTCTAAGGTTTTCCATCTTTTCTTCATTGTTTCTTCTTTCTCTCTATTTCTTATCTTATAATTTTCTCTACCATAAAGTTTATTCCTTGTTCCTCTCCTTACACCCTAAAAATAAAACCCTAGCCCACCTACCCTTGAGTGTAGGCAACCATCCTAGGGTTGTCCTACATCA"
        info_list.append([cenid,output_dir,monomer])          

    if argv2=="cen_RO60" or argv2=="all":      #VroSat2
        cenid="cen_RO60"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="AATTCGTTTCTCATCTTTGGGATAACTTCCTAAGTTTTAGTTTAGCTTAATCAATTAATC"
        info_list.append([cenid,output_dir,monomer])

                 
    if argv2=="cen_RO72" or argv2=="all":      #VroSat3
        cenid="cen_RO72"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="ACGGGGAATCCAATTCCGCTGAGTACCGAAGAAGAGATCGTGTCCGGTTGAGCACCAGAGAATACGAGACAG"
        info_list.append([cenid,output_dir,monomer])
    if argv2=="cen_RO124" or argv2=="all":
        cenid="cen_RO124"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="CATGAATTGGTTTCTCATCTTAAGGACAACTTCCTAAGTTTTAGTTTAGCTTAATCAATTAGTCAATTGGTTTCTTACCATTGGGACACCGTCCTAACCTTTAGTTTTACTTAATCAATTAATC"
        info_list.append([cenid,output_dir,monomer])        
    if argv2=="cen_RO197" or argv2=="all":
        cenid="cen_RO197"
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer="TCTTTTGGTCTTGGAATCGCCGGTACTCGGACGAATTCGATTCCCCGTCTCTACCGTACTCTTCGGTACTCGCCTGGATACCACCCTTTCTGTCTTGGATTCTCCGATACTCGCACGAATTTCATTCCCCGTTTGTCTTGTATTCTCCATTGCTCGACCGGATACGATATCTTCTCCGGTACTCGCGAGGATGCCAT"        
        info_list.append([cenid,output_dir,monomer])
    if "mononmer"  in args_dict and "cenid"  in args_dict: 
        cenid=args_dict["cenid"]
        output_dir="./chr2blast/1_blastn/"+cenid+'/'
        monomer=args_dict["mononmer"] 
        info_list.append([cenid,output_dir,monomer])

    if info_list==[]:
        print("Option 1: Missing parameter argv2, can be 'all' or 'cen_66', etc.");
        print("Option 2: Missing parameters -mononmer or -cenid (custom name like cen_66)");sys.exit() 
    info_list_num=len(info_list)
    process_serial=0
    for one in info_list:
        process_serial+=1
        print("Progress:"+str(process_serial)+'/'+str(info_list_num))
        cenid,output_dir,monomer=one
        ##########
        monomer_len=len(monomer) 

        #print('\t\tmonomer:'+monomer)    
        print('\t\tcenid:'+cenid) 
        print('\t\tMonomer length:'+str(monomer_len)) 
        print('\t\tOutput folder:'+output_dir)
        subprocess.run(["mkdir "+output_dir], shell=True)  
        if "||" not in monomer:
            subprocess.run([f"echo '>{cenid}\n{monomer}\n' > {output_dir}{cenid}"], shell=True)
            multi_mark=''
        else:
            multi_mark='yes'
            SEQs=monomer.split("||")
            kk=0;sum_str=''
            while kk< len(SEQs):
                SEQs_one=SEQs[kk]
                kk+=1
                sum_str+=f">{cenid}.....{kk}\n{SEQs_one}\n"
            with open(f"{output_dir}{cenid}",'w') as f:
                f.write(sum_str)
                        
        def run_step1(outdir_name):
            #print(outdir_name)
            output_dir=outdir_name[0]
            one_name=outdir_name[1]
            #if one_name!="PN40024_hap1":return False
            cenid=outdir_name[2]
            monomer_file=output_dir+cenid
            input_file1=f"./chr2blast/0_index/{one_name}/{one_name}"
            output_file1=output_dir+one_name+'.outfmt6'
            #print(output_file1);sys.exit()
            if  os.path.exists(output_file1+'.stat')==True:
                #print('skip'+output_file1);
                return False
            #-word_size 4/5 test    
            subprocess.run([f"blastn -outfmt 6 -evalue 1 -task blastn-short  -query {monomer_file} -db {input_file1} -out {output_file1}"], shell=True) 
            
            ### Sort
            input_file2=output_file1
            output_file2=output_file1+'.sort'
            subprocess.run([f"sort  -k2,2V -k9,9n {input_file2} > {output_file2}"], shell=True)  
            subprocess.run([f"rm {input_file2}"], shell=True)  
            ### Determine intervals
            input_file3=output_file2
            output_file3=output_file1+'.stat'
            current_monomer='';#current_chr='';current_match=0;current_chr_start='';current_chr_end='';current_strand=''
            current_match=0
            with open(output_file3,'w') as f2:
                f2.write('chr\tstrand\tchr_start\tchr_end\tlength\tmatch_len\tmatch_percent\n')
                with open(input_file3,'r') as f:
                    for line in f.readlines():
                        eachline_arr=line.strip().split('\t')
                        if multi_mark=='yes':
                            one_monomer=eachline_arr[0].split('.....')[0]
                        else:    
                            one_monomer=eachline_arr[0]
                        one_chr=eachline_arr[1]
                        one_match=int(eachline_arr[3])
                        one_start=int(eachline_arr[8])
                        one_end=int(eachline_arr[9])
                        #print(one_match,one_start,one_end)
                        if one_start<one_end:
                            one_strand='+'
                            one_min=one_start
                            one_max=one_end
                        else:    
                            one_strand='-'
                            one_min=one_end
                            one_max=one_start      
                        ######################    
                        if current_monomer=='':
                            current_monomer=one_monomer
                            current_monomer_num=1
                            current_chr=one_chr 
                            current_match_set=set()
                            kkk=one_min
                            while kkk<=one_max:
                                current_match_set.add(kkk)
                                kkk+=1
                            current_match=one_match                                  
                            current_strand=one_strand
                            current_min=one_min
                            current_max=one_max
                        elif one_monomer==current_monomer and one_chr==current_chr and one_min-current_max<1000:
                            #current_match+=one_match
                            kkk=one_min
                            while kkk<=one_max:
                                current_match_set.add(kkk)
                                kkk+=1
                            current_match=len(current_match_set)    
                            current_monomer_num+=1
                            if one_max>current_max:
                                current_max=one_max
                            #if current_chr=='Chr20':print(one_min)
                            #print(current_match)
                        else:
                            if current_match>3*monomer_len:#current_monomer_num>=3 and 
                                #print(current_monomer_num)
                                current_length=current_max-current_min+1
                                current_match_percent=round(current_match/current_length,3)
                                f2.write(current_chr+'\t'+current_strand+'\t'+str(current_min)+'\t'+str(current_max)+'\t'+str(current_length)+'\t'+str(current_match)+'\t'+str(current_match_percent)+'\n')
                            current_monomer=one_monomer
                            current_monomer_num=1
                            current_chr=one_chr 
                            current_match_set=set()
                            kkk=one_min
                            while kkk<=one_max:
                                current_match_set.add(kkk)
                                kkk+=1
                            current_match=one_match    
                            current_strand=one_strand
                            current_min=one_min
                            current_max=one_max        
                    if current_match>3*monomer_len:#current_monomer_num>=3 and 
                        #print(current_monomer_num)
                        current_length=current_max-current_min+1
                        current_match_percent=round(current_match/current_length,3)
                        f2.write(current_chr+'\t'+current_strand+'\t'+str(current_min)+'\t'+str(current_max)+'\t'+str(current_length)+'\t'+str(current_match)+'\t'+str(current_match_percent)+'\n')                    
                    
                    
        # Assign tasks to processes in the process pool    
        outdir_name_list=[]
        with open('./chr2blast/0_index/sample_list','r') as f:
            for one in f.readlines():
                one_name=one.strip()
                #if one_name!="V028.hap1":continue
                outdir_name_list.append([output_dir,one_name,cenid])
        #with Pool(processes=thread) as pool:
        #    pool.map(run_step1, outdir_name_list) 
            
        with Pool(processes=thread) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step1, outdir_name_list), start=1):
                progress = (i / len(outdir_name_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()           
        
        
if argv1=="stepall" or argv1=="step1" or argv1=="step1s":    


    cenid_list = []
    for file in os.listdir('./chr2blast/1_blastn/'):
        if os.path.isdir(os.path.join('./chr2blast/1_blastn/', file)):
            cenid_list.append(file)
            
    sample_list=[]
    #with open('./chr2blast/0_index/sample_list','r') as f:  sample_list = f.read().strip().split('\n')
    sample_list=[]
    with open('../samples_satellite/sample_info','r') as f:  
        for line in f:
            eachline_arr=line.strip().split('\t')
            sample=eachline_arr[0]
            if len(sample)==0:continue
            sample_list.append(sample) 
            
    print('Number of samples:'+str(len(sample_list)))
    print(cenid_list)
    ###
    i=0
    with open('./chr2blast/1_blastn/sum_tmp','w') as f2:
        f2.write("#\tsample\tcenid\tchr\tstrand\tchr_start\tchr_end\tlength\tmatch_len\tmatch_percent\n")
        for one_cenid in cenid_list:
            for one_sample in sample_list:
                input_file = f'./chr2blast/1_blastn/{one_cenid}/{one_sample}.outfmt6.stat'
                k=0
                with open(input_file,'r') as f:
                    for line in f.readlines():
                        eachline=line.strip()
                        k+=1
                        if k==1:continue
                        i+=1
                        f2.write(f"\t{one_sample}\t{one_cenid}\t{eachline}\n")
    subprocess.run(["sort -k1,1  -k3,3V  -k5,5V -k12,12n ./chr2blast/1_blastn/sum_tmp > ./chr2blast/1_blastn/sum; rm ./chr2blast/1_blastn/sum_tmp"], shell=True)               
    print('Number of result blocks:'+str(i))
    
if argv1=="stepall" or argv1=="step1s2":
    print('Organize the results')
    dict_sample_type_list={}
    with open("./chr2blast/1_blastn/sum",'r') as f:
        for line in f:
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=9:continue
            sample,cenid,chromosome,strand,chr_start,chr_end,length,match_len,match_percent=eachline_arr
            if sample not in dict_sample_type_list:dict_sample_type_list[sample]={}
            if cenid not in dict_sample_type_list[sample]: dict_sample_type_list[sample][cenid]=[]
            dict_sample_type_list[sample][cenid].append([int(chr_start),int(chr_end),strand])
    with open("./chr2blast/1_blastn/sum2",'w') as f:
        f.write(f"sample\tcenid\tnumber\tlength\tpositions\n")
        for sample,dict_cenid_list in dict_sample_type_list.items():
            for cenid ,one_list in dict_cenid_list.items():
                number=len(one_list)
                length=0
                str_list=[]
                for one in one_list:
                    one_start,one_end,strand=one
                    length+=abs(one_end-one_start)+1
                    str_list.append(f"[{one_start}-{one_end}]({strand})")
                str_all=','.join(str_list)    
                f.write(f"{sample}\t{cenid}\t{number}\t{length}\t{str_all}\n")    

if argv1=="stepall" or argv1=="step1s3":
    print('Organize results into wide format, average per haplotype, statistics by sample')    
    dict_sample_type_info={}
    with open("./chr2blast/1_blastn/sum2",'r') as f:    
        next(f)        
        for line in f:
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=5:continue
            samplehap,centype,num,length,info=eachline_arr
            if samplehap[-4:] not in ['hap1','hap2','hap3','hap4']:print('error, non-hap1234 exists');#sys.exit
            sample=samplehap[:-5]
            #print(sample)
            if sample not in dict_sample_type_info:
                dict_sample_type_info[sample]={}
            if centype not in dict_sample_type_info[sample]:  
                dict_sample_type_info[sample][centype]={}
                dict_sample_type_info[sample][centype]['length']=[]
                dict_sample_type_info[sample][centype]['num']=[]
            dict_sample_type_info[sample][centype]['length'].append(int(length))    
            dict_sample_type_info[sample][centype]['num'].append(int(num))   
    with open("./chr2blast/1_blastn/sum3",'w') as f:        
        f.write(f"sample\tcentype\tnum_str\tlength_str\tlength_average\n")       
        for sample,dict_centype_info in  dict_sample_type_info.items():
            for centype,dict_info in dict_centype_info.items():
                num_list= dict_info['num']
                length_list= dict_info['length']
                num_str='|'.join(str(num) for num in num_list)
                length_str='|'.join(str(num) for num in length_list)
                length_average = sum(length_list) / len(length_list)    
                f.write(f"{sample}\t{centype}\t{num_str}\t{length_str}\t{length_average}\n")
if argv1=="stepall" or argv1=="step1s4":        
    dict_sample_type_len={}
    with open("./chr2blast/1_blastn/sum3",'r') as f:    
        next(f)        
        for line in f:
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=5:continue   
            sample,centype,num_str,length_str,length_average=eachline_arr
            if sample not in dict_sample_type_len:dict_sample_type_len[sample]={}
            dict_sample_type_len[sample][centype]=int(round(float(length_average),0))
            
    with open("./chr2blast/1_blastn/sum4",'w') as f:        
        f.write(f"sample\tVSat1\tVSat2\tVSat3\tVSat4\tVSat5\tVSat6\n")           
        for sample,dict_centype_len in dict_sample_type_len.items():
            VSat1 = dict_centype_len.get('cen_107', 0)
            VSat2 = dict_centype_len.get('cen_66', 0)
            VSat3 = dict_centype_len.get('cen_103', 0)
            VSat4 = dict_centype_len.get('cen_191', 0)
            VSat5 = dict_centype_len.get('cen_355', 0)
            VSat6 = dict_centype_len.get('cen_383', 0)
            f.write(f"{sample}\t{VSat1}\t{VSat2}\t{VSat3}\t{VSat4}\t{VSat5}\t{VSat6}\n")           
            
if argv1=="step_generate_oneseq":
    if  os.path.exists('./chr2blast/ZZZ_____generate_oneseq')==False:
        subprocess.run(["mkdir ./chr2blast/ZZZ_____generate_oneseq"], shell=True)  
    print('Concatenate all scattered sequences into a continuous sequence')
    sample_list=[]
    with open('/home/lain/aaa_data/run0/samples_satellite/sample_info','r') as f:  
        for line in f:
            eachline_arr=line.strip().split('\t')
            sample=eachline_arr[0]
            if len(sample)==0:continue
            sample_list.append(sample) 
            
            
    dict_sample_inputlist={}
    with open('./chr2blast/1_blastn/sum','r') as f:
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            sample = eachline_arr[0]
            if sample not in sample_list:continue
            cenid=  eachline_arr[1]
        
            chromosome= eachline_arr[2]
            strand=eachline_arr[3]
            start=eachline_arr[4]
            end=eachline_arr[5]
            if sample not in dict_sample_inputlist: dict_sample_inputlist[sample]=[]
            dict_sample_inputlist[sample].append([sample,cenid,chromosome,strand,start,end])
    sample_used_list=(dict_sample_inputlist.keys())   
    #print('sample_used_list:')
    #print(sample_used_list)
    # Generate reverse complementary sequence
    def reverse_complement(sequence):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','N':'N'}
        reverse_complement_seq = ''.join(complement[base] for base in reversed(sequence))
        return reverse_complement_seq             
    
    #input_list_Len=len(input_list)
    #i=0
    #dict_centype_seqlist={}
    #for one_input in  input_list:
    def run_step(one_sample):
        input_list=dict_sample_inputlist[one_sample]
        print(one_sample+"\t\tNum:"+str(len(input_list)))
        dict_centype_seqlist={}
        for one_input in  input_list:
            #i+=1
            #print(f"{i}/{input_list_Len}",end='\r')
            sample,cenid,chromosome,strand,start,end=one_input
            pos=f"{chromosome}:{start}-{end}"
            #if chromosome !='BMNGT2T_hap1_chr01':return False
            #print(pos)
            CMD = ["samtools","faidx",f"./chr2blast/0_index/{sample}/{sample}.fa",f"{pos}"  ]
            # Execute command and check errors
            result = subprocess.run(CMD, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"Error: {sample} samtools failed - {result.stderr.strip()}")
                #continue#
                continue
            seq=result.stdout.strip()    
            #print(seq)
            seq=''.join(seq.split('\n')[1:]).strip().upper()
            if strand=='-':seq=reverse_complement(seq)
            #print([cenid,seq])
            
            if cenid not in dict_centype_seqlist: dict_centype_seqlist[cenid]=[]
            dict_centype_seqlist[cenid].append(seq)
        return dict_centype_seqlist
        
    dict_result_list={}
    with Pool(processes=thread) as pool:
        # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step, sample_used_list), start=1):
            #print(result)
            for cenid,seq_list in result.items():
                if cenid not in dict_result_list: dict_result_list[cenid]=[]
                for one_seq in seq_list:
                    dict_result_list[cenid].append(one_seq)
            progress = (i / len(sample_used_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush()  
            
    interval='N'*50000
    for cenid,seq_list in dict_result_list.items():
        seq_str=interval.join(seq_list)
        seq_num=len(seq_list)
        with open(f'./chr2blast/ZZZ_____generate_oneseq/{cenid}__num{seq_num}','w') as f:
            f.write(f'>{cenid}\n{seq_str}')
if argv1=="step_generate_oneseq_ultra_cen383":
    print('Extract given cenid and specific chromosomes, step_generate_oneseq_ultra cen_383 Chr15,Chr16,Chr10,Chr11')
    if  os.path.exists('./chr2blast/ZZZ_____generate_oneseq')==False:
        subprocess.run(["mkdir ./chr2blast/ZZZ_____generate_oneseq"], shell=True)  
    print('Concatenate all scattered sequences into a continuous sequence')
    sample_list=[]
    with open('/home/lain/aaa_data/run0/samples_satellite/sample_info','r') as f:  
        for line in f:
            eachline_arr=line.strip().split('\t')
            sample=eachline_arr[0]
            if len(sample)==0:continue
            sample_list.append(sample) 
            
            
    dict_sample_inputlist={}
    with open('./chr2blast/1_blastn/sum','r') as f:
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            sample = eachline_arr[0]
            if sample not in sample_list:continue
            cenid=  eachline_arr[1]
            if cenid!='cen_383':continue
            chromosome= eachline_arr[2]
            #print(chromosome)
            chr_used_list=["Chr15","Chr16","Chr10","Chr11"]
            bad_mark='yes'
            for onechr in chr_used_list:
                if onechr.lower() in chromosome.lower():bad_mark='';break
            if bad_mark=='yes':continue
            strand=eachline_arr[3]
            start=eachline_arr[4]
            end=eachline_arr[5]
            if sample not in dict_sample_inputlist: dict_sample_inputlist[sample]=[]
            dict_sample_inputlist[sample].append([sample,cenid,chromosome,strand,start,end])
    sample_used_list=(dict_sample_inputlist.keys())   
    #print(sample_used_list)
    #print('sample_used_list:')
    #print(sample_used_list)
    # Generate reverse complementary sequence
    def reverse_complement(sequence):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','N':'N'}
        reverse_complement_seq = ''.join(complement[base] for base in reversed(sequence))
        return reverse_complement_seq             
    
    #input_list_Len=len(input_list)
    #i=0
    #dict_centype_seqlist={}
    #for one_input in  input_list:
    with open("./chr2blast/ZZZ_____generate_oneseq/VSat6_seqs.fa",'w') as ff:
        ff.write("")    
    def run_step(one_sample):
        input_list=dict_sample_inputlist[one_sample]
        print(one_sample+"\t\tNum:"+str(len(input_list)))
        dict_centype_seqlist={}
        for one_input in  input_list:
            #i+=1
            #print(f"{i}/{input_list_Len}",end='\r')
            sample,cenid,chromosome,strand,start,end=one_input
            pos=f"{chromosome}:{start}-{end}"
            #if chromosome !='BMNGT2T_hap1_chr01':return False
            #print(pos)
            
            CMD = ["samtools","faidx",f"./chr2blast/0_index/{sample}/{sample}.fa",f"{pos}"  ]
            # Execute command and check errors
            result = subprocess.run(CMD, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"Error: {sample} samtools failed - {result.stderr.strip()}")
                #continue#
                continue
            seq=result.stdout.strip()    
            #print(seq)
            seq=''.join(seq.split('\n')[1:]).strip().upper()
            if strand=='-':seq=reverse_complement(seq)
            #print([cenid,seq])
            with open("./chr2blast/ZZZ_____generate_oneseq/VSat6_seqs.fa",'a') as ff:
                ff.write(f">{sample}_{chromosome}:{pos}\n{seq}\n")
            if cenid not in dict_centype_seqlist: dict_centype_seqlist[cenid]=[]
            dict_centype_seqlist[cenid].append(seq)
        return dict_centype_seqlist
        
    dict_result_list={}
    with Pool(processes=thread) as pool:
        # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step, sample_used_list), start=1):
            #print(result)
            for cenid,seq_list in result.items():
                if cenid not in dict_result_list: dict_result_list[cenid]=[]
                for one_seq in seq_list:
                    dict_result_list[cenid].append(one_seq)
            progress = (i / len(sample_used_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush()  
            
    interval='N'*50000
    for cenid,seq_list in dict_result_list.items():
        seq_str=interval.join(seq_list)
        seq_num=len(seq_list)
        with open(f'./chr2blast/ZZZ_____generate_oneseq/{cenid}__num{seq_num}_Chr10,11,15,16','w') as f:
            f.write(f'>{cenid}\n{seq_str}')
if argv1=="step_generate_oneseq_ultra_LIRSat1":
    print('Extract given cenid and specific chromosomes, step_generate_oneseq_ultra LIRSat1, only analyze chromosomes without VSat6')
    if  os.path.exists('./chr2blast/ZZZ_____generate_oneseq')==False:
        subprocess.run(["mkdir ./chr2blast/ZZZ_____generate_oneseq"], shell=True)  
    print('Concatenate all scattered sequences into a continuous sequence')
    sample_list=[]
    with open('/home/lain/aaa_data/run0/samples_satellite/sample_info','r') as f:  
        for line in f:
            eachline_arr=line.strip().split('\t')
            sample=eachline_arr[0]
            if len(sample)==0:continue
            sample_list.append(sample) 
    
    #Find chromosomes without VSat6
    dict_samplechr_len={}
    with open('./chr2blast/1_blastn/sum','r') as f:
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            sample,cenid,chromsome,strand,chr_start,chr_end,length,match_len,match_percent=eachline_arr
            if cenid!='cen_383':continue
            samplechr=sample+"|"+chromsome
            if samplechr not in dict_samplechr_len:dict_samplechr_len[samplechr]=0
            dict_samplechr_len[samplechr]+=int(length)
            
            
                        
    dict_sample_inputlist={}
    with open('./chr2blast/1_blastn/sum','r') as f:
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            sample = eachline_arr[0]
            if sample not in sample_list:continue
            cenid=  eachline_arr[1]
            if cenid!='LIRSat1':continue
            chromosome= eachline_arr[2]
            samplechr=sample+"|"+chromsome
            if samplechr in dict_samplechr_len:continue
            #print(chromosome)
            chr_used_list=["Chr15","Chr16","Chr10","Chr11"]
            bad_mark=''
            for onechr in chr_used_list:
                if onechr.lower() in chromosome.lower():bad_mark='yes';break
            if bad_mark=='yes':continue
            strand=eachline_arr[3]
            start=eachline_arr[4]
            end=eachline_arr[5]
            if sample not in dict_sample_inputlist: dict_sample_inputlist[sample]=[]
            dict_sample_inputlist[sample].append([sample,cenid,chromosome,strand,start,end])
    sample_used_list=(dict_sample_inputlist.keys())   
    #print(sample_used_list)
    #print('sample_used_list:')
    #print(sample_used_list)
    # Generate reverse complementary sequence
    def reverse_complement(sequence):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','N':'N'}
        reverse_complement_seq = ''.join(complement[base] for base in reversed(sequence))
        return reverse_complement_seq             
    
    #input_list_Len=len(input_list)
    #i=0
    #dict_centype_seqlist={}
    #for one_input in  input_list:
    with open("./chr2blast/ZZZ_____generate_oneseq/LIRSat1_seqs.fa",'w') as ff:
        ff.write("")
    def run_step(one_sample):
        input_list=dict_sample_inputlist[one_sample]
        print(one_sample+"\t\tNum:"+str(len(input_list)))
        dict_centype_seqlist={}
        for one_input in  input_list:
            #i+=1
            #print(f"{i}/{input_list_Len}",end='\r')
            sample,cenid,chromosome,strand,start,end=one_input
            pos=f"{chromosome}:{start}-{end}"
            #if chromosome !='BMNGT2T_hap1_chr01':return False
            #print(pos)
            
            CMD = ["samtools","faidx",f"./chr2blast/0_index/{sample}/{sample}.fa",f"{pos}"  ]
            # Execute command and check errors
            result = subprocess.run(CMD, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"Error: {sample} samtools failed - {result.stderr.strip()}")
                #continue#
                continue
            seq=result.stdout.strip()    

            #print(seq)
            seq=''.join(seq.split('\n')[1:]).strip().upper()
            if strand=='-':seq=reverse_complement(seq)
            #print([cenid,seq])
            with open("./chr2blast/ZZZ_____generate_oneseq/LIRSat1_seqs.fa",'a') as ff:
                ff.write(f">{sample}_{chromosome}:{pos}\n{seq}\n")
                
            if cenid not in dict_centype_seqlist: dict_centype_seqlist[cenid]=[]
            dict_centype_seqlist[cenid].append(seq)
        return dict_centype_seqlist
        
    dict_result_list={}
    with Pool(processes=thread) as pool:
        # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step, sample_used_list), start=1):
            #print(result)
            for cenid,seq_list in result.items():
                if cenid not in dict_result_list: dict_result_list[cenid]=[]
                for one_seq in seq_list:
                    dict_result_list[cenid].append(one_seq)
            progress = (i / len(sample_used_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush()  
            
    interval='N'*50000
    for cenid,seq_list in dict_result_list.items():
        seq_str=interval.join(seq_list)
        seq_num=len(seq_list)
        with open(f'./chr2blast/ZZZ_____generate_oneseq/{cenid}__num{seq_num}_noVSat6','w') as f:
            f.write(f'>{cenid}\n{seq_str}')            

    
    
    
    
    
    
    
    
    
    
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))      
time_end=timeit.default_timer()
print('All the running time: %.0f Seconds'%(time_end-time_start))