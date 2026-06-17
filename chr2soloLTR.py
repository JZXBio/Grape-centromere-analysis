#!/usr/bin/python
# -*- coding: UTF-8 -*-
import sys
argvs=sys.argv
argv_len=len(sys.argv)
if not argvs[1].startswith("part"):
    if argv_len < 2 or argvs[1] in ("h", "-h", "help", "-help") or not argvs[1].startswith("step"):
        print ("chr2soloLTRs.py-----help:")
        print ("")
        print ("Usage：")
        print ("chr2soloLTRs.py   ## Directly use the genome organized by chr2EDTA step0")
        print ("chr2soloLTRs.py   step1        #makeblastdb build index for each sample")
        print ("chr2soloLTRs.py   step2        #Take about 10 examples each of typical ClassI, II, II (based on amino acid tree), BLAST to find all (solo)-LTRs (300k)")
        print ("chr2soloLTRs.py   step3        #All cdhit0.8 build tree to find Tekay ClassII, then extract only ClassII————————————————————————step results 3_filter.all (100k) extract V126, V008 for separate tree building")
        print ("chr2soloLTRs.py   step4        #Build tree only for ClassII")        
        print (" ")
        print ("chr2soloLTRs.py   step14        #Re-analyze, based on step3, now I just want to distinguish a few ancestral sequences of ClassII")        
        print ("chr2soloLTRs.py   part        #Build separate tree for specific samples using cdhit")        
        print ("chr2soloLTRs.py   part_nocdhit        ##Build separate tree for specific samples without cdhit")   
        print ("chr2soloLTRs.py   part_pos2cluster              #For specific samples without building tree, classify into clusters based on step14")
        print ("chr2soloLTRs.py   part_pos2cluster_plus         #For specific samples (specific positions) without building tree, classify into clusters based on step14")
        
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
else:thread=70    #### Too many threads may cause crashes

if  os.path.exists('./chr2soloLTRs')!=True:
        subprocess.run(["mkdir ./chr2soloLTRs"], shell=True)  
if argv1=="stepall" or argv1=="step0":    
    subprocess.run(["cp -r ./chr2EDTA/0_prepare ./chr2soloLTRs/0_prepare"], shell=True)  


with open('./chr2EDTA/0_prepare/sample_list','r') as f:
    sample_list=f.read().split('\n')
    sample_list = [x for x in sample_list if x or x == 0 or x is False]
    sample_list.sort()
print('sample_list:')
print(len(sample_list))

if argv1=="stepall" or argv1=="step1"  or argv1=="step1.0": 
    print('makeblastdb build index')
    if  os.path.exists('./chr2soloLTRs/1_index')!=True:
        subprocess.run(["mkdir ./chr2soloLTRs/1_index"], shell=True)      
    def run_step(one):
        subprocess.run([f"makeblastdb -in ./chr2soloLTRs/0_prepare/{one}/{one}.fasta -dbtype nucl -out ./chr2soloLTRs/1_index/{one}/{one}"], shell=True)  
                    
    with Pool(processes=thread) as pool:
        # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step, sample_list), start=1):
            # Results can be processed here, e.g., stored or printed
            progress = (i / len(sample_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush()    


if argv1=="stepall" or argv1=="step1"  or argv1=="step1s":  
    print('Use the copied information table to extract sequences')
    '''
    soloLTR="TGTAACACCTAGATTATTTTAGTACTTAACTTAAGCTTGCTTAATTAGATATTAAAACTAGTTTAGTGCTAATCATGTTTAATTATGAATTAAACTTTGATTAAGGTTAAGTGGGATTAGTTAATGATCTAAGCATGCTTATTAGGTTCTTAGGGACTAATTAAGGTTATGAAAGCCTAAAGGACTAATGGGCAATTATGGGTTTCATTTTTGATAACTTGAAGGACTAAAGTGCAAATTTAGAAACTTGAGTTAGTGGAATGCCACCTCCTACTTGGTGGATTGGTGGCAGCCATGTGGGCTGCCACCTCATGCTTGGGTGCTTGGAGACCCCTTTATAAGGGGGGCTGCAGCTCGTTTTGCACCATTTCACCTGCAGCAACTTGGGTTAGAGAGAGAGACTAGAGAGAGAAAGTGTAGATTTGAGGTAAGCAAGTGTTTAATTTTGTAATTTTCGTTTATTTTGGTTCCTAATGGTATGTTGAGAAGGATTAATGGTAAAAATTGTGTAAAAGTTGAGTATAATTAGCTATAAATCTGTTTTAGTTAAAAATCACAATTTATTAAAGATTAAAGTATTTTAAGTTAAGTATTTTATTAAGGGTAAGATCAGCATAAATCTTTGTTTAATTTGGTTTGGTTGGAAAATAGATTAGTAAACATGTGATGAATTGGGTTATTTGGACACTTAGGATTTGTTCATAGGTTAGGCTTTAAGAAAATCAAAAGAATTAGTGTTTGGTAATTAATTAGGGAAATTAATATGCATTGTAAATATTGTTTAATTCAATTCATATTTGATTTGAAAACATGTAGGATATAGATTAGAGAATGTGAGAATTAAAAATTTGCAACAATAAACTAAAACCTAAATTATGTTGTTTCATTTCAGTTCCTCGTAATAAGCAAAGATCGCCTACTTAGAAGAAACACAAAAAGGAAGTCGGTGTAAGGCAGGGAATTTTATGCTATTCTATGGTGTATCTGATTTCTTTCCTTGAATTATTTGTTGGAATTTTATGTTGTTTCTAGGATTTCATAAAATCTTTTAAAGTGTCATTGCATCACGGTTTATTGCATTGAAATGTTACCATGATGTTGGAACGTATATTGGTATGAAGATTCGTGGTTTTGTGTTGCCGAAATGGTTTATGGTGTGATTGCATTGCAAGAAGGATATGGAATTATCACAGATCATATTCAAAGAATTTATTGGTGGCGAAGTGCTTATGTGATGCTACAGGCCTTATCATCTGATATATTGTTTATGTGGGACCGTGGGTCCTTGGGTGAAAGTCCCTAAAGCCCTCGAGGAAGACACTCCGATGTGGGTGTATGGGGTTGGTCCTGCCCCTGGGTTTTAGTCCCTAAAGACACGGGGTTGGTCCTGCCCTTGGGTATGAGTCCCAAAAGACACGGGGTTGGTCCTGCCCTTGGGTATGAGTCCCTAAAGACACGGGGTATGACTTGCCCTTGGGTGATGTCCCAGAATAGTCATTATTATATTGATCAAGTATCTTGTGGAGCATTGATATCTGCTGTGTACATTGTCGGGAGTCAGCAATTTATCAGTTGTGAAAGGATGGATGGGGAAAAGCAATGATGAAACCAACAAACACATGCATACATGTTTGACATGATTACACATTTGTTAATGATTATGAAATATTTATCATGCATGTTATTAAACGTTTGATTCAAGGTTTTTAAGGGTATACTTAGGATGGTTATAAACTTTCCTACTGAGTTGTGAACTCACCCTATCCCTTCCACCTCTAGATGCAGGTCAGAAGTCCTATGCAGGAAAGAATGCTTGAGCAATGCTTTACTGTGATTGGATGCTGTTTGCTCATATTGAAAGTCTTTTGAGGCCATGCCTTTTGTATAAAGACTAAACCTTTTGTTGGAATGATGTAATGTACATATTTGTCAGATACACTTGTAGTATGGGCTACCCGTAGTAAACAATGGGGTTTGGTATATGTAAATTAATGTTGTACAAGTTTCTTTTGGGAAACAAAACATATTTTGTTCACTTCATAGTCACAATGTTTAATACAGGATGTACTCTTTATAACAGGTTTTCCAAATCCTCTTTTTGAACAAATTCAAGTAGGAACTCAACATTTAACTTTTGGAAAGCACTTATACTGTATTGAGTATCAATTGCTTAGTTGAAGAATAAATAATAATAATAATAAGAAAAAAAAAATAATGGGGTGTGACA"
    with open('./chr2soloLTRs/soloLTR_example.fa','w') as f:
        f.write(f">soloLTR_example\n{soloLTR}\n")'''
        
    # I want to take a few examples from different Tekay types
    bed='''
    V030.hap2	Chr7	-	21336456	21338597	Class1
    MuscatHamburg_hap2	Chr11	-	7841109	7843211	Class1
    V032.hap2	Chr9	-	13932161	13934272	Class1
    PinotNoir2_hap2	Chr10	-	2597383	2599541	Class1
    V019.hap1	Chr6	+	14616707	14618973	Class1
    V040.hap1	Chr9	+	4289189	4291438	Class1
    V043.hap2	Chr8	+	17885281	17887542	Class1
    V067.hap1	Chr3	+	19689785	19692033	Class1
    V076.hap2	Chr8	-	1956741	1958763	Class1
    V031.hap2	Chr5	-	4807246	4809231	Class1
    V065.hap2	Chr3	-	14655393	14657567	Class2
    Hongmunage_hap2	Chr7	+	9909843	9912057	Class2
    V012.hap2	Chr5	+	14217242	14219490	Class2
    V019.hap1	Chr13	+	21779876	21782089	Class2
    V037.hap2	Chr16	+	10325347	10327550	Class2
    V041.hap2	Chr18	-	15960934	15963172	Class2
    V050.hap1	Chr16	+	9623531	9625732	Class2
    V061.hap2	Chr7	+	14572511	14574768	Class2
    V081.hap1	Chr12	+	13063078	13065325	Class2
    V092.hap1	Chr7	+	9768987	9771201	Class2
    V106.hap2	Chr16	+	6625426	6627602	Class2
    Baimunage_hap1	Chr5	+	16633058	16635314	Class3
    PN40024_hap2	Chr13	+	29083661	29085880	Class3
    ShineMuscat_hap2	Chr10	-	1947884	1950101	Class3
    V019.hap1	Chr16	+	12462595	12464823	Class3
    V050.hap2	Chr12	+	4864418	4866628	Class3
    V062.hap1	Chr11	+	17791755	17793983	Class3
    V066.hap1	Chr18	+	1539522	1541783	Class3
    V087.hap2	Chr5	+	16523232	16525488	Class3
    V108.hap1	Chr2	+	22572598	22574824	Class3
    VHP-T2T.hap1	Chr13	+	28755649	28757868	Class3
    '''
    bed_lines=bed.strip().split('\n')
    iii=0
    with open('./chr2soloLTRs/soloLTR_examples.fa','w') as f:
        for oneline in bed_lines:
            iii+=1
            oneline=oneline.strip().split('\t')
            sample,chromosome,strand,start,end,mark=oneline
            addition_para=''
            if strand=='-':
                CMD = [
                    "samtools", "faidx",
                    f"/home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample}/{sample}.fasta",
                    "--reverse-complement",  # Correct: separate option
                    "--mark-strand", "sign",  # Correct: specify type (e.g., "sign" for (+) / (-))
                    f"{chromosome}:{start}-{end}"
                ]
            else:
                CMD = ["samtools","faidx",f"/home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample}/{sample}.fasta",f"{chromosome}:{start}-{end}"  ]
            result = subprocess.run(CMD, capture_output=True, text=True)

            if result.returncode != 0:
                print(f"Error: {result.stderr}")
            info_seq = result.stdout.strip()
            seq="".join(info_seq.split('\n')[1:])
            f.write(f">{mark}_{iii}\n{seq}\n")

if argv1=="stepall" or argv1=="step2"  or argv1=="step2.0":  
    print('31 example sequences, BLAST to all genomes, 22s')
    if  os.path.exists('./chr2soloLTRs/2_blastn')!=True:
        subprocess.run(["mkdir ./chr2soloLTRs/2_blastn"], shell=True)      
    def run_step(one):
        #subprocess.run([f"blastn -query ./chr2soloLTRs/soloLTR_example.fa -db  ./chr2soloLTRs/1_index/{one}/{one}  -out ./chr2soloLTRs/2_blastn/{one}.blastn -outfmt 6"], shell=True)  
        subprocess.run([f"blastn -query ./chr2soloLTRs/soloLTR_examples.fa -db  ./chr2soloLTRs/1_index/{one}/{one}  -out ./chr2soloLTRs/2_blastn/{one}.blastn -outfmt 6"], shell=True)  
                    
    with Pool(processes=thread) as pool:
        # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step, sample_list), start=1):
            # Results can be processed here, e.g., stored or printed
            progress = (i / len(sample_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush()   


if argv1=="stepall" or argv1=="step2"  or argv1=="step2.1":  
    print('Remove overlaps. Multiple sequences BLASTed together will definitely have duplicates. Sort by score from high to low, create an empty dictionary, then check each one for 80% similarity with recorded entries. If not, store it. Finally, remove those shorter than 500bp, 23s')
    if  os.path.exists('./chr2soloLTRs/2_blastn')!=True:
        subprocess.run(["mkdir ./chr2soloLTRs/2_blastn"], shell=True)      
    
    MIN_length=500   # soloLTR needs length > 500bp    
        
    def run_step(one):
        result_arr=[]
        with open(f"./chr2soloLTRs/2_blastn/{one}.blastn",'r') as f:
            for line in f:
                eachline_arr=line.strip().split('\t') 
                query_name=eachline_arr[0]
                chromosome=eachline_arr[1]
                start=int(eachline_arr[8])
                end=int(eachline_arr[9])
                percent=eachline_arr[2]
                pvalue=eachline_arr[-2]
                strand = '-' if end - start < 0 else '+'  # Determine strand direction
                bitscore=float(eachline_arr[-1])
                # Ensure interval is in [min, max] format
                interval_start, interval_end = min(start, end), max(start, end)                
                result_arr.append([query_name, chromosome,strand, interval_start, interval_end, percent,pvalue, bitscore])
        #print(result_arr[0])        
        # Sort descending
        result_arr_sorted = sorted(result_arr, key=lambda x: x[-1], reverse=True)        
        ###
        def is_80_percent_overlap(A, B, C, D):
            len_AB = B - A
            len_CD = D - C
            target_overlap = 0.8 * len_CD  # 80% of CD length
        
            # Calculate actual overlap
            overlap_start = max(A, C)
            overlap_end = min(B, D)
            overlap = overlap_end - overlap_start
        
            # Determine if condition is met
            return overlap >= target_overlap
        
        # Example test
        #A, B = 14101367, 14103569  # Original interval [10, 20]
        #C, D = 14101367, 14102824  # New interval [12, 18]
        #print(is_80_percent_overlap(A, B, C, D))  # Output True or False
        ###
        dict_chromosome_pos={}
        with open(f"./chr2soloLTRs/2_blastn/{one}.filter",'w') as f:
            for one_result in result_arr_sorted:
                query_name, chromosome,strand, start, end, percent,pvalue, bitscore=one_result
                if chromosome not in dict_chromosome_pos:
                    dict_chromosome_pos[chromosome]=[]
                bad_mark = False
                
                for x in dict_chromosome_pos[chromosome]:
                    x_query_name, x_chromosome, x_strand, x_start, x_end, x_percent, x_pvalue, x_bitscore = x
                    # No need to consider strand direction
                    if is_80_percent_overlap(x_start, x_end, start, end):
                        bad_mark = True
                        break
                if not bad_mark:
                    dict_chromosome_pos[chromosome].append([ query_name, chromosome,strand, start, end, percent,pvalue, bitscore])
                    if end-start>MIN_length:
                        ####
                        if strand=='-':
                            CMD = [
                                "samtools", "faidx",
                                f"/home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{one}/{one}.fasta",
                                "--reverse-complement",  # Correct: separate option
                                "--mark-strand", "sign",  # Correct: specify type (e.g., "sign" for (+) / (-))
                                f"{chromosome}:{start}-{end}"
                            ]
                        else:
                            CMD = ["samtools","faidx",f"/home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{one}/{one}.fasta",f"{chromosome}:{start}-{end}"  ]
                        result = subprocess.run(CMD, capture_output=True, text=True)
                        info_seq = result.stdout.strip()
                        seq="".join(info_seq.split('\n')[1:])
                        ###
                        f.write(f'{one}\t{query_name}\t{chromosome}\t{strand}\t{start}\t{end}\t{percent}\t{pvalue}\t{bitscore}\t{seq}\n')
                    
    with Pool(processes=thread) as pool:
        # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step, sample_list), start=1):
            # Results can be processed here, e.g., stored or printed
            progress = (i / len(sample_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush()  

if argv1=="stepall" or argv1=="step2" or argv1=="step2.2": 
    print('Summarize sequences, 2s')
    ii=0
    with open('./chr2soloLTRs/2_blastn.all.seq','w') as f3:
        with open('./chr2soloLTRs/2_blastn.all','w') as f2:
            f2.write("ID\tsample\tquery_name\tchromosome\tstrand\tstart\tend\tpercent\tpvalue\tbitscore\tseq\n")
            for one in sample_list:
                with open(f"./chr2soloLTRs/2_blastn/{one}.filter",'r') as f:
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        sample,query_name, chromosome,strand, start, end, percent,pvalue, bitscore,seq=eachline_arr
                        
                        ii+=1
                        f2.write(f"{ii}\t{eachline}\n")
                        f3.write(f">{ii}\n{seq}\n")
                    
if argv1=="stepall" or argv1=="step3" or argv1=="step3.0":    
    print('cd-hit threshold 0.8, 100s')
    if  os.path.exists('./chr2soloLTRs/3_cdhit')!=True:
        subprocess.run(["mkdir ./chr2soloLTRs/3_cdhit"], shell=True)       
    cdhid_c=0.8  ## Default 0.95, cluster at this similarity level 
    env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
    subprocess.run([f"{env_source_str} cd-hit -M 300000 -T 70 -d 100 -c {cdhid_c} -n 5 -i ./chr2soloLTRs/2_blastn.all.seq -o ./chr2soloLTRs/3_cdhit/0_cdhit"], shell=True)            
                
if argv1=="stepall" or argv1=="step3" or argv1=="step3.1":                     
    print('Parse cd-hit results 0s')            
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
    clstr_file = f"./chr2soloLTRs/3_cdhit/0_cdhit.clstr"
    dict_core_others = parse_cdhit_clstr(clstr_file)
    ######
    ii=0
    dict_seq_serial={}
    with open(f"./chr2soloLTRs/3_cdhit/0_cdhit.index",'w') as f:
        f.write(f"group1\tgroup2\trepresent_index\tseq_index\n")
        for core,others in dict_core_others.items():
            ii+=1
            f.write(f'cluster_{ii}\tcluster_{ii}_0\t{core}\t{core}\n')
            kk=0
            for one in others:
                kk+=1
                f.write(f"cluster_{ii}\tcluster_{ii}_{kk}\t{core}\t{one}\n")  
        
            
if argv1=="stepall" or argv1=="step3" or argv1=="step3.2": 
    print('Merge the initial examples and representative cluster sequences, perform cdhit together, 522s')
    subprocess.run([f"cat  ./chr2soloLTRs/soloLTR_examples.fa ./chr2soloLTRs/3_cdhit/0_cdhit >./chr2soloLTRs/3_cdhit/2_cdhit "], shell=True)
    env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
    subprocess.run([f"{env_source_str} mafft --auto ./chr2soloLTRs/3_cdhit/2_cdhit > ./chr2soloLTRs/3_cdhit/2_mafft.fasta"], shell=True)
    
if argv1=="stepall" or argv1=="step3" or argv1=="step3.3":  
    print('trimal, 6s')
    trimal_para='-automated1'
    env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate orthofinder && "
    cmd=f"{env_source_str} trimal {trimal_para} -in ./chr2soloLTRs/3_cdhit/2_mafft.fasta -out  ./chr2soloLTRs/3_cdhit/3_trimAl.fasta "
    subprocess.run([cmd], shell=True)
    env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
    subprocess.run([f'{env_source_str} fasta_formatter -i ./chr2soloLTRs/3_cdhit/3_trimAl.fasta -o ./chr2soloLTRs/3_cdhit/3_trimAl.fa'], shell=True)
    subprocess.run([f'mv ./chr2soloLTRs/3_cdhit/3_trimAl.fa ./chr2soloLTRs/3_cdhit/3_trimAl.fasta'], shell=True)


if argv1=="stepall" or argv1=="step3" or argv1=="step3.4":
    print('iqtree2')
    env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
    cmd=f"{env_source_str} iqtree2 -s ./chr2soloLTRs/3_cdhit/3_trimAl.fasta -m MFP  -bb 1000 -nt AUTO -mem 0.8 -pre ./chr2soloLTRs/3_cdhit/4_iqtree.fasta "  #-m GTR+G4
    print(cmd)
    subprocess.run([cmd], shell=True)
    if  os.path.exists(f'./chr2soloLTRs/3_cdhit/5_iTOL.Newick')!=True:
        subprocess.run([f"touch ./chr2soloLTRs/3_cdhit/5_iTOL.Newick"], shell=True)
    print('Open the tree in iTOL, click on a branch, you can use "add clade to pruned tree", export the needed branch as newick')

if argv1=="stepall" or argv1=="step3" or argv1=="step3.6":
    print('After visualizing the tree obtained from iqtree in iTOL, to get the sample order, click export to Newick format, copy into 5_iTOL.Newick')
    import re
    with open(f"./chr2soloLTRs/3_cdhit/5_iTOL.Newick2",'w') as f2:
        with open(f"./chr2soloLTRs/3_cdhit/5_iTOL.Newick",'r') as f:
            text=f.read().strip()
        cleaned_text = re.sub(r'\[.*?\]', '', text)
        f2.write(cleaned_text)
    from ete3 import Tree
    # Read Newick file
    tree = Tree(f"./chr2soloLTRs/3_cdhit/5_iTOL.Newick2")
    # Get sample order (tip labels)
    tip_order = [leaf.name for leaf in tree.get_leaves()]
    
    
    with open(f"./chr2soloLTRs/3_cdhit/6_serial",'w') as f:
        f.write('seq_index\tserial\n')
        ii=0
        for one in tip_order:
            ii+=1
            f.write(f"{one}\t{ii}\n")

if argv1=="stepall" or argv1=="step3" or argv1=="step3.7":    
    print('Now 6_serial contains only the representative sequences for Class2. Go back to extract only the Class2 part, 150s')
    seq_id_good_list=[]
    with open(f"./chr2soloLTRs/3_cdhit/6_serial",'r') as f:
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            seq_id=eachline_arr[0]
            seq_id_good_list.append(seq_id)
            
    print('Get filtered list, Tekay is about 1/5')        
    seq_id_good_list2=[]    
    with open(f"./chr2soloLTRs/3_cdhit/7_filter.index",'w') as f2:
        with open(f"./chr2soloLTRs/3_cdhit/0_cdhit.index",'r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                group1,group2,represent_index,seq_index=eachline_arr
                if represent_index in seq_id_good_list:
                    seq_id_good_list2.append(seq_index)
                    f2.write(f"{seq_index}\n")
    
    print('Extract sequences.')  
    kk=0
    with open(f"./chr2soloLTRs/3_filter.all.seq",'w') as f3:                           
        with open(f"./chr2soloLTRs/3_filter.all",'w') as f2:      
            f2.write(f"ID\tsample\tquery_name\tchromosome\tstrand\tstart\tend\tpercent\tpvalue\tbitscore\tseq\n")
            with open(f"./chr2soloLTRs/2_blastn.all",'r') as f:    
                next(f)
                for line in f:
                    kk+=1
                    print(f'{kk}',end='\r')
                    eachline=line.strip()
                    eachline_arr=line.strip().split('\t')
                    ID,sample,query_name,chromosome,strand,start,end,percent,pvalue,bitscore,seq=eachline_arr
                    if ID in seq_id_good_list2:
                        f2.write(f"{eachline}\n")
                        f3.write(f">{ID}\n{seq}\n")

#step4 continues with phylogenetic tree for ClassII
if 'step4' in    argv1 and 1==2:     
    if argv1=="stepall" or argv1=="step4" or argv1=="step4.0":    
        print('cd-hit threshold 0.9, 151s')
        if  os.path.exists('./chr2soloLTRs/4_cdhit')!=True:
            subprocess.run(["mkdir ./chr2soloLTRs/4_cdhit"], shell=True)       
        cdhid_c=0.9  ## Default 0.95, cluster at this similarity level 
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
        subprocess.run([f"{env_source_str} cd-hit -M 300000 -T 70 -d 100 -c {cdhid_c} -n 5 -i ./chr2soloLTRs/3_filter.all.seq -o ./chr2soloLTRs/4_cdhit/0_cdhit"], shell=True)            
    if argv1=="stepall" or argv1=="step4" or argv1=="step4.1":                     
        print('Parse cd-hit results 0s')            
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
        clstr_file = f"./chr2soloLTRs/4_cdhit/0_cdhit.clstr"
        dict_core_others = parse_cdhit_clstr(clstr_file)
        ######
        ii=0
        dict_seq_serial={}
        with open(f"./chr2soloLTRs/4_cdhit/0_cdhit.index",'w') as f:
            f.write(f"group1\tgroup2\trepresent_index\tseq_index\n")
            for core,others in dict_core_others.items():
                ii+=1
                f.write(f'cluster_{ii}\tcluster_{ii}_0\t{core}\t{core}\n')
                kk=0
                for one in others:
                    kk+=1
                    f.write(f"cluster_{ii}\tcluster_{ii}_{kk}\t{core}\t{one}\n")  
            
                
    if argv1=="stepall" or argv1=="step4" or argv1=="step4.2": 
        print('Do NOT merge the initial examples！！！！！！！！！, only keep representative cluster sequences, perform cdhit together, 44s')
        subprocess.run([f"cat  ./chr2soloLTRs/4_cdhit/0_cdhit >./chr2soloLTRs/4_cdhit/2_cdhit "], shell=True)   #./chr2soloLTRs/soloLTR_examples.fa 
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
        subprocess.run([f"{env_source_str} mafft --auto ./chr2soloLTRs/4_cdhit/2_cdhit > ./chr2soloLTRs/4_cdhit/2_mafft.fasta"], shell=True)
        
    if argv1=="stepall" or argv1=="step4" or argv1=="step4.3":  
        print('trimal, 16s')
        trimal_para='-automated1'
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate orthofinder && "
        cmd=f"{env_source_str} trimal {trimal_para} -in ./chr2soloLTRs/4_cdhit/2_mafft.fasta -out  ./chr2soloLTRs/4_cdhit/3_trimAl.fasta "
        subprocess.run([cmd], shell=True)
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
        subprocess.run([f'{env_source_str} fasta_formatter -i ./chr2soloLTRs/4_cdhit/3_trimAl.fasta -o ./chr2soloLTRs/4_cdhit/3_trimAl.fa'], shell=True)
        subprocess.run([f'mv ./chr2soloLTRs/4_cdhit/3_trimAl.fa ./chr2soloLTRs/4_cdhit/3_trimAl.fasta'], shell=True)
    
    
    if argv1=="stepall" or argv1=="step4" or argv1=="step4.4":
        print('iqtree2, 2.5h')
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "  ##JS is neighbor-joining, because lengths vary greatly; maximum likelihood doesn't work well
        cmd=f"{env_source_str} iqtree2 -s ./chr2soloLTRs/4_cdhit/3_trimAl.fasta -m MFP -bb 1000 -nt AUTO -mem 0.8 -pre ./chr2soloLTRs/4_cdhit/4_iqtree.fasta "  #-m C 
        print(cmd)
        subprocess.run(["rm ./chr2soloLTRs/4_cdhit/4_iqtree.fasta.ckp.gz"], shell=True)
        subprocess.run([cmd], shell=True)
        if  os.path.exists(f'./chr2soloLTRs/4_cdhit/5_iTOL.Newick')!=True:
            subprocess.run([f"touch ./chr2soloLTRs/4_cdhit/5_iTOL.Newick"], shell=True)
        print('Open the tree in iTOL, click on a branch, you can use "add clade to pruned tree", export the needed branch as newick')
    
    if argv1=="stepall" or argv1=="step4" or argv1=="step4.6":
        print('After visualizing the tree obtained from iqtree in iTOL, to get the sample order, click export to Newick format, copy into 5_iTOL.Newick')
        import re
        with open(f"./chr2soloLTRs/4_cdhit/5_iTOL.Newick2",'w') as f2:
            with open(f"./chr2soloLTRs/4_cdhit/5_iTOL.Newick",'r') as f:
                text=f.read().strip()
            cleaned_text = re.sub(r'\[.*?\]', '', text)
            f2.write(cleaned_text)
        from ete3 import Tree
        # Read Newick file
        tree = Tree(f"./chr2soloLTRs/4_cdhit/5_iTOL.Newick2")
        # Get sample order (tip labels)
        tip_order = [leaf.name for leaf in tree.get_leaves()]
        
        
        with open(f"./chr2soloLTRs/4_cdhit/6_serial",'w') as f:
            f.write('seq_index\tserial\n')
            ii=0
            for one in tip_order:
                ii+=1
                f.write(f"{one}\t{ii}\n")
    
    if argv1=="stepall" or argv1=="step4" or argv1=="step4.7":    
        print('Now 6_serial contains only the representative sequences for Class2. Go back to extract only the Class2 part, 150s')
        dict_seqindex_serial={} 
        with open(f"./chr2soloLTRs/4_cdhit/6_serial",'r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                seqindex,serial=eachline_arr
                dict_seqindex_serial[seqindex]=serial
                
        print('Get list')        
        dict_seqindex_cdhitindex={}
        with open(f"./chr2soloLTRs/4_cdhit/0_cdhit.index",'r') as f:
            next(f)
            for line in f:
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                _,_,cdhit_index,seq_index=eachline_arr
                dict_seqindex_cdhitindex[seq_index]=cdhit_index
        
        print('Extract sequences.')  
        result_list=[]        
        with open(f"./chr2soloLTRs/3_filter.all",'r') as f:
            next(f)
            for line in f:
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                seq_index=eachline_arr[0]
                cdhit_index=dict_seqindex_cdhitindex[seq_index]
                serial=dict_seqindex_serial[cdhit_index]
                result_list.append([eachline,seq_index,cdhit_index,serial])
                        
        sorted_result=sorted(result_list, key=lambda x:  int(x[3]))      
    
        with open(f"./chr2soloLTRs/4_filter.all",'w') as f:
            for one in sorted_result:
                eachline,seq_index,cdhit_index,serial=one
                f.write(f"{eachline}\t{cdhit_index}\t{serial}\n")
        
#step14 extracts 1900bp ClassII sequences and builds tree        
if argv1=="stepall" or "step14" in argv1:
    if  os.path.exists('./chr2soloLTRs/14_simple_long')!=True:
        subprocess.run(["mkdir ./chr2soloLTRs/14_simple_long"], shell=True)      
    print('Re-analyze, based on step3, now I just want to distinguish a few ancestral sequences of ClassII')
    if argv1=="stepall" or argv1=="step14" or argv1=="step14.0":  
        print("14.0 — Extract all soloLTRs longer than 1900bp, total 25290")
        with open(f"./chr2soloLTRs/14_simple_long/0_long.fasta",'w') as f3:
            with open(f"./chr2soloLTRs/14_simple_long/0_long.info",'w') as f2:
                with open(f"./chr2soloLTRs/3_filter.all",'r') as f:
                    next(f)
                    for line in f:
                        eachline=line.strip()
                        eachline_arr=eachline.split('\t')
                        ID=eachline_arr[0]
                        seq=eachline_arr[-1]
                        seq_len=len(seq)
                        if seq_len>1900:
                            f2.write(f"{eachline}\n")
                            f3.write(f">{ID}\n{seq}\n")
    if argv1=="stepall" or argv1=="step14" or argv1=="step14.1":   
        print('cdhit0.9 gives only 200+ clusters, cdhit0.95 gives 600+ clusters, 27s')
        cdhid_c=0.90  ## Default 0.95, cluster at this similarity level 
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
        subprocess.run([f"{env_source_str} cd-hit -M 300000 -T 70 -d 100 -c {cdhid_c} -n 5 -i ./chr2soloLTRs/14_simple_long/0_long.fasta -o ./chr2soloLTRs/14_simple_long/1_cdhit"], shell=True)                      
    if argv1=="stepall" or argv1=="step14" or argv1=="step14.1s":     
        # Example call
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
        clstr_file = f"./chr2soloLTRs/14_simple_long/1_cdhit.clstr"
        dict_core_others = parse_cdhit_clstr(clstr_file)
        ######
        ii=0
        dict_seq_serial={}
        with open(f"./chr2soloLTRs/14_simple_long/1_cdhit.index",'w') as f:
            f.write(f"group1\tgroup2\trepresent_index\tseq_index\n")
            for core,others in dict_core_others.items():
                ii+=1
                f.write(f'cluster_{ii}\tcluster_{ii}_0\t{core}\t{core}\n')
                kk=0
                for one in others:
                    kk+=1
                    f.write(f"cluster_{ii}\tcluster_{ii}_{kk}\t{core}\t{one}\n")  
    if argv1=="stepall" or argv1=="step14" or argv1=="step14.1s2":
        
        def parse_cdhit_clstr(clstr_file):
            dict_core_others = {}
            current_cluster = None
            core_seq = None
            other_seqs = []
        
            with open(clstr_file, "r") as f:
                for line in f:
                    line = line.strip()
                    if line.startswith(">Cluster"):
                        if current_cluster is not None:
                            dict_core_others[core_seq] = other_seqs.copy()
                            other_seqs.clear()
                        current_cluster = line
                    else:
                        parts = line.split(">")
                        seq_name = parts[1].split("...")[0]
                        if line.endswith("*"):
                            core_seq = seq_name
                        else:
                            other_seqs.append(seq_name)
                
                if current_cluster is not None:
                    dict_core_others[core_seq] = other_seqs
        
            return dict_core_others
        
        # Read sequence length information (ID -> length)
        print("Reading sequence length information...")
        seq_lengths = {}  # ID -> length
        with open("./chr2soloLTRs/14_simple_long/0_long.info", 'r') as f:
            next(f)  # Skip header
            for line in f:
                eachline_arr = line.strip().split('\t')
                ID = eachline_arr[0]
                seq = eachline_arr[-1]
                seq_lengths[ID] = len(seq)
        print(f"Read length information for {len(seq_lengths)} sequences")
        
        clstr_file = "./chr2soloLTRs/14_simple_long/1_cdhit.clstr"
        dict_core_others = parse_cdhit_clstr(clstr_file)
        
        # Statistics for each parent ID: size and length information
        print("\n" + "="*80)
        print("Cluster statistics: representative sequence -> cluster size and length information")
        print("="*80)
        
        total_clusters = 0
        total_sequences = 0
        
        # Output sorted by cluster size
        cluster_stats = []
        for core, others in dict_core_others.items():
            cluster_size = len(others) + 1  # +1 for the representative sequence itself
            
            # Get lengths of all sequences in the cluster
            all_lengths = [seq_lengths.get(core, 0)]  # Representative sequence length
            for seq_id in others:
                all_lengths.append(seq_lengths.get(seq_id, 0))
            
            # Filter out zero-length (may not be found)
            valid_lengths = [l for l in all_lengths if l > 0]
            
            if valid_lengths:
                avg_len = sum(valid_lengths) / len(valid_lengths)
                min_len = min(valid_lengths)
                max_len = max(valid_lengths)
                core_len = seq_lengths.get(core, 0)
            else:
                avg_len = min_len = max_len = core_len = 0
            
            cluster_stats.append({
                'core': core,
                'size': cluster_size,
                'core_len': core_len,
                'avg_len': avg_len,
                'min_len': min_len,
                'max_len': max_len
            })
            total_clusters += 1
            total_sequences += cluster_size
        
        # Sort by cluster size descending
        cluster_stats.sort(key=lambda x: x['size'], reverse=True)
        
        # Print header
        print(f"{'Representative':<25} {'Cluster Size':<8} {'Core Length':<12} {'Average Length':<12} {'Min Length':<12} {'Max Length':<12}")
        print("-"*80)
        
        for stat in cluster_stats:
            print(f"{stat['core']:<25} {stat['size']:<8} {stat['core_len']:<12} {stat['avg_len']:<12.1f} {stat['min_len']:<12} {stat['max_len']:<12}")
        
        print("="*80)
        print(f"Total clusters: {total_clusters}")
        print(f"Total sequences: {total_sequences}")
        print("="*80)
        
        # Save detailed statistics to file
        with open("./chr2soloLTRs/14_simple_long/1_cdhit_stats.txt", 'w') as f:
            f.write("represent_seq\tcluster_size\tcore_len\tavg_len\tmin_len\tmax_len\n")
            for stat in cluster_stats:
                f.write(f"{stat['core']}\t{stat['size']}\t{stat['core_len']}\t{stat['avg_len']:.1f}\t{stat['min_len']}\t{stat['max_len']}\n")
            f.write(f"\n# Total clusters: {total_clusters}\n")
            f.write(f"# Total sequences: {total_sequences}\n")
        
        print(f"\nStatistics saved to: ./chr2soloLTRs/14_simple_long/1_cdhit_stats.txt")
        
        # Optional: print length distribution summary
        print("\nLength distribution summary:")
        sizes = [s['size'] for s in cluster_stats]
        print(f"  Maximum cluster size: {max(sizes)}")
        print(f"  Minimum cluster size: {min(sizes)}")
        print(f"  Average cluster size: {sum(sizes)/len(sizes):.1f}")
        print(f"  Median cluster size: {sorted(sizes)[len(sizes)//2]}")
        
        # Generate index file (if needed)
        with open("./chr2soloLTRs/14_simple_long/1_cdhit.index", 'w') as f:
            f.write("group1\tgroup2\trepresent_index\tseq_index\n")
            for ii, (core, others) in enumerate(dict_core_others.items(), start=1):
                f.write(f'cluster_{ii}\tcluster_{ii}_0\t{core}\t{core}\n')
                for kk, one in enumerate(others, start=1):
                    f.write(f"cluster_{ii}\tcluster_{ii}_{kk}\t{core}\t{one}\n")
                    
    if argv1=="stepall" or argv1=="step14" or argv1=="step14.1_addtion" and 1==2:     
        print('Add outgroup')
        ##ClassIII_tekay_LTR from Chr15:12973756-12976006(-), 2250bp, LTRidentity=0.9996
        outgroup_seq=""">ClassIII_tekay_LTR
TGTAATGACCGGGTATTTTATGCCACGTTAGCAATTTTAGTTTGAAAAGTCTTAATTTATGTGATGGTTGAAAAGTCAAAAAGAAAAGCGTTTGAAGAACACGTGTCAATTTCGGATTGGTGAATTTCGTTTTATTGCGTCGGTCAATTAAACGAGTTGAAATTTTATGCCATGTCAGCCTTAGGATAGAAAATTCTTAGGATTTTAGAAATTGGAAATTTCCGGGAACGAAAAACGAAAACGAAAATGAAAAATAAATAAATTAAGATTAATTAGAATTAATTAAGAAAAATAAATTAATTAATTAAATTAAATTGATGAATAATAATATTGAGAAAATTAGTTTTTAATTGGTGTTAATAAATAAAATTTATGTGTTAAAAAAAAAATAAATTGGAAATCAAATAAAACAATTAGAATAAACAAAAGAAATAAAATAAATATATTCAATGGGCTTTGAATATGGGTTGTTAAGGTAATGGGTTGTAAAAGAAAAAAAATTTAATGGAGTTGGGCTTAAGTGGCTGGTGGCTTTTAAAAGCCTAGTGGGCTGTTAATAAAAAGAAAGAAAAAAATGGGTTGGGCTTTTGAAGTTTGGGCTAGGCTTGTGTGGCTGGATGGTTAGTAAAAAAAAAAGGGGGCAGTAGCACGCACGCAGGAGAGAGAAACAGAAAAGGGCACTGTGCATAAATTAGGGCACTCGCCGGAATTTTGATCGGACGTTTGGACGGCCAAACGATCTCCGTTTCGGCCCAAATTTGGAGGAAGTACTCTATTCGACTAGACGAACAACCCTGCGGTGGCTGATTTTGTTTTTAACGGCCAAATTTTCAGATCTGAGGTAGGGGACCCTAACCCTAAAACTCAAATTTTTATGTTTTTCCCTTGAAAAAAAATTGTAAATATTGTCATTATGAGTTAATTAGGAATTATGGATGTTGTGTGAATTTTTCTAGAATATTTGGAATTTGTTTGGAATTTTTGGAATTTCGGGAAATTGATTTTGATTGATAATGAATTGTTTTGAAAGGTTAAAAATTATTAGATGGGTTGGAAAATTCCAAAATTTGGGTTAAATTTTAAATATAATGATTTGTGAATTTTTAGGCATTGGTTCGAATATTTTTGAAAATTTTAGGGTAAATAGTTAGTTAATTTTGGATATTAGAAATTATTATAATGGTTGAAAAATTTTGTAAGTGTGGTAAAATCTCAAGATTGTTGATATGCAATTTTTTTTTTAATATTTATTTGGGATATTATCAGAAATTGTAATGTGTGAAATTGTTTGATTGGATTATTGGGGAATGTTGAGATTGTTGAATTATTGAGGAACATTGGAATTGTGGCATTCATTGCATCATGGTTATGTGAAATAAATAAATAAATAAATAAAATCATGAAGATTGTGAAAAGTGTATGAAATGGAAAAAAATGAAATAGAGATGAGAAGTGGAAGCCCGTGTGAGGCGAATTAAGAAGGGAGAGAGATCCCTAGGGTGAAAGACCCAAAAGTTTACCCCGGGAAGACTCCGAATGGGGAGTGTATTTGGGTGCGGACGCGTATCCTTCGGTGAAAACCCGAGAACAATAGTGCATTGTATGACTGTGTGTCACACGTTCATGTGCATTCATGTAATTGTTGAATATTGTGAAATTGGGTATAAAATTATATATATTAAACTATGTGGTTTGATTGACGCTGTGTGAATTGTTCCTGTGGTGTTGAAGTGTTCTTTTCGTATTCTTTGAACTTAGTAATCCCTTTAACCCCTTGGGTGTTAGGCACCCTACTGAGCAATGTGGAAATGCTCACCCCTTCCTTGTTACACCTTTTTAGATGCAGATATCTCTCCTGAGGATCCACAGGTGGGACAGGGAGGACTTCAGGATGCTCCTGGAGCGGCCTAGTAGAGCGTGGCATTTATTTGTTATTGTGTTCCTTTTTAGATTTTGGGGATTTGTTTTAGCAATTATGACTCATGGATTTGCTATGAAGCTTTATTTGATTAAGTTGTTAATTGTGAACTTTTATTTTGGACTATAATGAATATTTATTTAATCTTGTTGTTTTGAGTAGTTTTGGGATATTGTAATTCATGAGAATTCTAGTGAGATGTATAAAAAAAAAAAAAAAATCCAGAGTGTTTTGGTTGACTGCCAACATGGAATAGGAGGAACCCTTAGGGTCAAACCTTTGACCTGGGGTCACGGGTCAAATTTTGGGTCGCCCAGAATTTGGGTCGTGACA\n"""
        with open(f"./chr2soloLTRs/14_simple_long/1_outgroup",'w') as f:
            f.write(outgroup_seq)
        subprocess.run([f"cat ./chr2soloLTRs/14_simple_long/1_outgroup ./chr2soloLTRs/14_simple_long/1_cdhit > ./chr2soloLTRs/14_simple_long/1_cdhit_withoutgroup"], shell=True)       
            
    if argv1=="stepall" or argv1=="step14" or argv1=="step14.2":     
        print('200+ protein mafft 68s, 600+ protein mafft')
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
        subprocess.run([f"{env_source_str} mafft  --auto ./chr2soloLTRs/14_simple_long/1_cdhit > ./chr2soloLTRs/14_simple_long/2_mafft.fasta"], shell=True)   
    

    
    if argv1=="stepall" or argv1=="step14" or argv1=="step14.3":            
        print('trimal, 3s')
        trimal_para='-gt 0.2' #-automated1
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate orthofinder && "
        cmd=f"{env_source_str} trimal {trimal_para} -in ./chr2soloLTRs/14_simple_long/2_mafft.fasta -out  ./chr2soloLTRs/14_simple_long//3_trimAl.fasta "
        subprocess.run([cmd], shell=True)
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
        subprocess.run([f'{env_source_str} fasta_formatter -i ./chr2soloLTRs/14_simple_long//3_trimAl.fasta -o ./chr2soloLTRs/14_simple_long//3_trimAl.fa'], shell=True)
        subprocess.run([f'mv ./chr2soloLTRs/14_simple_long//3_trimAl.fa ./chr2soloLTRs/14_simple_long//3_trimAl.fasta'], shell=True)        
                
    if argv1=="stepall" or argv1=="step14" or argv1=="step14.3_trim":   
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
        subprocess.run([f"{env_source_str} seqkit subseq --region 51:-101 ./chr2soloLTRs/14_simple_long//3_trimAl.fasta -o ./chr2soloLTRs/14_simple_long//3_trimAl.fa"], shell=True)
        subprocess.run([f'mv ./chr2soloLTRs/14_simple_long//3_trimAl.fa ./chr2soloLTRs/14_simple_long//3_trimAl.fasta'], shell=True)  
        subprocess.run([f'{env_source_str} fasta_formatter -i ./chr2soloLTRs/14_simple_long//3_trimAl.fasta -o ./chr2soloLTRs/14_simple_long//3_trimAl.fa'], shell=True)
        subprocess.run([f'mv ./chr2soloLTRs/14_simple_long//3_trimAl.fa ./chr2soloLTRs/14_simple_long//3_trimAl.fasta'], shell=True) 
        
    if argv1=="stepall" or argv1=="step14" or argv1=="step14.4":              
        print('iqtree2')
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "  ##JS is neighbor-joining, because lengths vary greatly; maximum likelihood doesn't work well
        #cmd=f"{env_source_str} iqtree2 -s ./chr2soloLTRs/14_simple_long//3_trimAl.fasta -m MFP+MERGE -bb 1000 -nt AUTO --bnni -mem 0.8 -pre ./chr2soloLTRs/14_simple_long//4_iqtree.fasta -o ClassIII_tekay_LTR" #JC
        print('No outgroup set')
        cmd=f"{env_source_str} iqtree2 -s ./chr2soloLTRs/14_simple_long//3_trimAl.fasta -m MFP+MERGE -bb 1000 -nt AUTO --bnni -mem 0.8 -pre ./chr2soloLTRs/14_simple_long//4_iqtree.fasta "
        subprocess.run([cmd], shell=True)
        if  os.path.exists(f'./chr2soloLTRs/14_simple_long/5_iTOL.Newick')!=True:
            subprocess.run([f"touch ./chr2soloLTRs/14_simple_long/5_iTOL.Newick"], shell=True)
        print('Open the tree in iTOL, click on a branch, you can use "add clade to pruned tree", export the needed branch as newick')       
    if argv1=="stepall" or argv1=="step14" or argv1=="step14.5":
        print('After visualizing the tree obtained from iqtree in iTOL, to get the sample order, click export to Newick format, copy into 5_iTOL.Newick')
        import re
        with open(f"./chr2soloLTRs/14_simple_long/5_iTOL.Newick3",'w') as f2:
            with open(f"./chr2soloLTRs/14_simple_long/5_iTOL.Newick",'r') as f:
                text=f.read().strip()
            cleaned_text = re.sub(r'\[.*?\]', '', text)
            f2.write(cleaned_text)
        from ete3 import Tree
        # Read Newick file
        tree = Tree(f"./chr2soloLTRs/14_simple_long/5_iTOL.Newick3")
        # Get sample order (tip labels)
        tip_order = [leaf.name for leaf in tree.get_leaves()]
        
        
        with open(f"./chr2soloLTRs/14_simple_long/5_serial",'w') as f:
            f.write('seq_index\tserial\n')
            ii=0
            for one in tip_order:
                ii+=1
                f.write(f"{one}\t{ii}\n")
    if argv1=="stepall" or argv1=="step14" or argv1=="step14.6":
        from Bio import Phylo
        # Read tree without node names
        tree = Phylo.read("./chr2soloLTRs/14_simple_long/5_iTOL.Newick", "newick")
        # Add names to internal nodes
        for i, clade in enumerate(tree.find_clades(terminal=False)):
            clade.name = f"Node{i+1}"
        # Save tree with node names
        Phylo.write(tree, "./chr2soloLTRs/14_simple_long/5_iTOL.Newick2", "newick")

        print('treetime consensus tree')
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate tree && " 
        cmd=f"{env_source_str} treetime ancestral --aln ./chr2soloLTRs/14_simple_long/2_mafft.fasta --tree ./chr2soloLTRs/14_simple_long/5_iTOL.Newick2  --outdir ./chr2soloLTRs/14_simple_long/6_treetime_ancestral"
        subprocess.run([cmd], shell=True)
  
    if argv1=="stepall" or argv1=="step14" or argv1=="step14.7":
        print('TreeCluster to group the tree')
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate tree && " 
        cmd=f"{env_source_str} TreeCluster.py -i ./chr2soloLTRs/14_simple_long/5_iTOL.Newick2 -o ./chr2soloLTRs/14_simple_long/7_clusters.tsv -t 0.15 " # -m leaf_dist_avg -t 0.1 -v
        subprocess.run([cmd], shell=True)
        
    if argv1=="stepall" or argv1=="step14" or argv1=="step14.8" and 1==2:    
        print('step14.7 generated 17 clusters. Use iTOL to open the tree (./chr2soloLTRs/14_simple_long/5_iTOL.Newick2 with 7_clusters.tsv for coloring), manually click nodes and record')
        goodNodes=["Node2","Node8","Node10","Node12","Node16","Node18","Node21","Node23","Node60","Node68","Node71","Node87","Node141","Node178","Node184","Node190","Node200"]
        dict_node_seq={}
        with open('./chr2soloLTRs/14_simple_long/6_treetime_ancestral/ancestral_sequences.fasta','r') as f:
            for line in f:
                eachline=line.strip()
                if eachline[0]=='>':name=eachline[1:]; dict_node_seq[name]=''
                else:dict_node_seq[name]+=eachline
        #print(dict_node_seq)        
        print('Output')        
        with open ("./chr2soloLTRs/14_simple_long/8_goodnodes_2",'w') as f2:    
            with open ("./chr2soloLTRs/14_simple_long/8_goodnodes",'w') as f:    
                for one in goodNodes:
                    seq=dict_node_seq[one]
                    f.write(f'>{one}\n{seq}\n')
                    seq=seq.replace('-','')
                    f2.write(f'>{one}\n{seq}\n')
        
    if argv1=="stepall" or argv1=="step14" or argv1=="step14.9" and 1==2:  
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
        subprocess.run([f"{env_source_str} mafft  --auto ./chr2soloLTRs/14_simple_long/8_goodnodes_2 > ./chr2soloLTRs/14_simple_long/9_node_1.mafft"], shell=True)             
        ##
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate orthofinder && "
        cmd=f"{env_source_str} trimal -automated1 -in ./chr2soloLTRs/14_simple_long/9_node_1.mafft -out  ./chr2soloLTRs/14_simple_long/9_node_2.trimAl"
        subprocess.run([cmd], shell=True)
        ##
        env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "  ##JS is neighbor-joining, because lengths vary greatly; maximum likelihood doesn't work well
        cmd=f"{env_source_str} iqtree2 -s ./chr2soloLTRs/14_simple_long/9_node_2.trimAl -m JC -bb 1000 -nt AUTO -mem 0.8 -pre ./chr2soloLTRs/14_simple_long/9_node_3.iqtree"
        subprocess.run([cmd], shell=True)        
        
   
        
        
name_list=[]
name_list=[]
name_list=['PN40024','V050','V060','V008','V126']
if argv1=="part":
    
    if os.path.exists(f'./chr2soloLTRs/part/')!=True :
         subprocess.run(['mkdir ./chr2soloLTRs/part/'], shell=True)
    def run_step(one_name):
        folder=f'./chr2soloLTRs/part/{one_name}'
        if os.path.exists(f'{folder}')!=True :
             subprocess.run([f'mkdir {folder}'], shell=True)  
        if  os.path.exists(f'{folder}/0_input.fa ')!=True :   
            with open(f'{folder}/0_input.fa','w') as f3:
                with open(f'{folder}/0_input.info','w') as f2:
                    f2.write('ID\tsample\tquery_name\tchromosome\tstrand\tstart\tend\tpercent\tpvalue\tbitscore\tseq\n')
                    with open('./chr2soloLTRs/3_filter.all','r') as f:
                        next(f)
                        for line in f:
                            eachline=line.strip()
                            eachline_arr=eachline.split('\t')
                            ID,sample,query_name,chromosome,strand,start,end,percent,pvalue,bitscore,seq=eachline_arr
                            if one_name not in sample:continue
                            f2.write(f'{eachline}\n')
                            f3.write(f">{ID}\n{seq}\n")
             
        if os.path.exists(f'{folder}/1_cdhit')!=True :
            print('cd-hit threshold 0.9, 151s')
            cdhid_c=0.95  ## Default 0.95, cluster at this similarity level 
            env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
            subprocess.run([f"{env_source_str} cd-hit -M 300000 -T 70 -d 100 -c {cdhid_c} -n 5 -i {folder}/0_input.fa -o {folder}/1_cdhit"], shell=True)     
            
        if os.path.exists(f'{folder}/1_cdhit.index')!=True :    
            print('Parse cd-hit results 0s')            
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
            clstr_file = f"{folder}/1_cdhit.clstr"
            dict_core_others = parse_cdhit_clstr(clstr_file)
            ######
            ii=0
            dict_seq_serial={}
            with open(f"{folder}/1_cdhit.index",'w') as f:
                f.write(f"group1\tgroup2\trepresent_index\tseq_index\n")
                for core,others in dict_core_others.items():
                    ii+=1
                    f.write(f'cluster_{ii}\tcluster_{ii}_0\t{core}\t{core}\n')
                    kk=0
                    for one in others:
                        kk+=1
                        f.write(f"cluster_{ii}\tcluster_{ii}_{kk}\t{core}\t{one}\n")  
        if os.path.exists(f'{folder}/2_mafft')!=True :  
            print('mafft')
            env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
            subprocess.run([f"{env_source_str} mafft  --auto {folder}/1_cdhit > {folder}/2_mafft"], shell=True)
            
            
        if os.path.exists(f'{folder}/3_trimal.fasta')!=True :  
            trimal_para='-automated1'
            env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate orthofinder && "
            cmd=f"{env_source_str} trimal {trimal_para} -in {folder}/2_mafft -out  {folder}/3_trimal.fasta"
            subprocess.run([cmd], shell=True)
        
            env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
            subprocess.run([f'{env_source_str} fasta_formatter -i {folder}/3_trimal.fasta -o {folder}/3_trimal.fa'], shell=True)
            subprocess.run([f'mv {folder}/2_trimal.fa {folder}/3_trimal.fasta'], shell=True)
            
   
        if os.path.exists(f'{folder}/4_iqtree.treefile')!=True:        
            env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
            cmd=f"{env_source_str} iqtree2 -s {folder}/3_trimal.fasta -m GTR+G4 -bb 1000 -nt AUTO -mem 0.8 -pre {folder}/4_iqtree"
            print(cmd)
            subprocess.run([cmd], shell=True)
            
        if os.path.exists(f'{folder}/5_newick')!=True:                  
            subprocess.run([f'touch {folder}/5_newick'], shell=True)


        if os.path.exists(f'{folder}/6_serial')!=True:    
            import re
            with open(f"{folder}/5_newick2",'w') as f2:
                with open(f"{folder}/5_newick",'r') as f:
                    text=f.read().strip()
                cleaned_text = re.sub(r'\[.*?\]', '', text)
                f2.write(cleaned_text)
            from ete3 import Tree
            # Read Newick file
            tree = Tree(f"{folder}/5_newick2")
            # Get sample order (tip labels)
            tip_order = [leaf.name for leaf in tree.get_leaves()]
            
            
            with open(f"{folder}/6_serial",'w') as f:
                f.write('seq_index\tserial\n')
                ii=0
                for one in tip_order:
                    ii+=1
                    f.write(f"{one}\t{ii}\n")
        
        if os.path.exists(f'{folder}/7_info')!=True:  
            dict_index_represent={}
            with open (f"{folder}/1_cdhit.index",'r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    _,_,represent_index,seq_index=eachline_arr
                    dict_index_represent[seq_index]=represent_index
            dict_representindex_serial={}        
            with open (f"{folder}/6_serial",'r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    representindex,serial=eachline_arr
                    dict_representindex_serial[representindex]=serial
            
            result_list=[]     
            with open(f"{folder}/0_input.info",'r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    ID,sample,query_name,chromosome,strand,start,end,percent,pvalue,bitscore,seq=eachline_arr
                    represent_index=dict_index_represent[ID]
                    representindex_serial=int(dict_representindex_serial[represent_index])
                    end=int(end)
                    result_list.append([ID,represent_index,representindex_serial,sample,query_name,chromosome,strand,start,end,percent,pvalue,bitscore,seq])
            sorted_list = sorted(result_list, key=lambda x: (x[2],x[5],x[8]), reverse=False)     
            ii=0
            with open(f"{folder}/7_input.info",'w') as f:
                f.write('ID\trepresent_index\trepresentindex_serial\tserial_revise\tsample\thap\tquery_name\tchromosome\tstrand\tstart\tend\tpercent\tpvalue\tbitscore\tseq\n')
                for one in sorted_list:
                    ID,represent_index,representindex_serial,sample,query_name,chromosome,strand,start,end,percent,pvalue,bitscore,seq=one
                    if "hap1" in sample:hap='hap1'
                    elif "hap2" in sample:hap='hap2'
                    ii+=1
                    f.write(f"{ID}\t{represent_index}\t{representindex_serial}\t{ii}\t{sample}\t{hap}\t{query_name}\t{chromosome}\t{strand}\t{start}\t{end}\t{percent}\t{pvalue}\t{bitscore}\t{seq}\n")
            
        if os.path.exists(f'{folder}/8_plotinput')!=True:
            with open(f"{folder}/8_plotinput",'w') as f2:
                f2.write(f"hap\tchromosome\tID\tserial_revise\txmin\txmax\tymin\tymax\n")
                with open(f"{folder}/7_input.info",'r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        ID,represent_index,representindex_serial,serial_revise,sample,hap,query_name,chromosome,strand,start,end,percent,pvalue,bitscore,seq=eachline_arr
                        if "hap1" in hap:ymin=0;ymax=0.5
                        elif "hap2" in hap:ymin=-0.5;ymax=0
                        xmin=int(serial_revise)-1
                        xmax=serial_revise
                        chromosome_num=int(chromosome.split('Chr')[-1])
                        ymin=ymin-chromosome_num
                        ymax=ymax-chromosome_num
                        f2.write(f"{hap}\t{chromosome}\t{ID}\t{serial_revise}\t{xmin}\t{xmax}\t{ymin}\t{ymax}\n")
            R_txt=r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file=read.table('8_plotinput', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
        color_values=c("Chr1" = "red",  
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
                                          "Chr19" = "#FF8C00")
        # Create plot object
        p <- ggplot() +
            geom_rect(data = input_file, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = chromosome))+
            
            scale_fill_manual(values = color_values, drop = FALSE)+
           # facet_wrap( ~ chromosome) +
          
        theme_classic() +         
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
         
          
          # Save as PDF
          pdf(file = paste0('8_plot', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
          print(p)
          dev.off()
        
        }
            '''
            with open(f'{folder}/plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)  
            new_directory = folder
            os.chdir(new_directory)
            subprocess.run(['Rscript plot.R'], shell=True)    
            os.chdir('./')                            
    with Pool(processes=thread) as pool:
        # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step, name_list), start=1):
            # Results can be processed here, e.g., stored or printed
            progress = (i / len(name_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush()  


name_list=['PN40024','V050','V060','V008','V126','V031']
name_list=['V008']
if argv1=="part_nocdhit":
    if os.path.exists(f'./chr2soloLTRs/part_nocdhit/')!=True :
         subprocess.run(['mkdir ./chr2soloLTRs/part_nocdhit/'], shell=True)
    def run_step(one_name):
        folder=f'./chr2soloLTRs/part_nocdhit/{one_name}'
        if os.path.exists(f'{folder}')!=True :
             subprocess.run([f'mkdir {folder}'], shell=True)  
        if  os.path.exists(f'{folder}/0_input.fa ')!=True :   
            with open(f'{folder}/0_input.fa','w') as f3:
                with open(f'{folder}/0_input.info','w') as f2:
                    f2.write('ID\tsample\tquery_name\tchromosome\tstrand\tstart\tend\tpercent\tpvalue\tbitscore\tseq\n')
                    with open('./chr2soloLTRs/3_filter.all','r') as f:
                        next(f)
                        for line in f:
                            eachline=line.strip()
                            eachline_arr=eachline.split('\t')
                            ID,sample,query_name,chromosome,strand,start,end,percent,pvalue,bitscore,seq=eachline_arr
                            if one_name not in sample:continue
                            if "hap" not in sample:continue
                            
                            f2.write(f'{eachline}\n')
                            f3.write(f">{ID}\n{seq}\n")
             
        
        if os.path.exists(f'{folder}/2_mafft')!=True :  
            print('mafft')
            env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
            subprocess.run([f"{env_source_str} mafft --auto {folder}/0_input.fa > {folder}/2_mafft"], shell=True)
            
            
        if os.path.exists(f'{folder}/3_trimal.fasta')!=True :  
            trimal_para='-automated1'
            env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate orthofinder && "
            cmd=f"{env_source_str} trimal {trimal_para} -in {folder}/2_mafft -out  {folder}/3_trimal.fasta"
            subprocess.run([cmd], shell=True)
        
            env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
            subprocess.run([f'{env_source_str} fasta_formatter -i {folder}/3_trimal.fasta -o {folder}/3_trimal.fa'], shell=True)
            subprocess.run([f'mv {folder}/3_trimal.fa {folder}/3_trimal.fasta'], shell=True)
            
   
        if os.path.exists(f'{folder}/4_iqtree.treefile')!=True:        
            env_source_str="source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate env0 && "
            cmd=f"{env_source_str} iqtree2 -s {folder}/3_trimal.fasta -m GTR+G4 -bb 1000 -nt AUTO -mem 0.8 -pre {folder}/4_iqtree"
            print(cmd)
            subprocess.run([cmd], shell=True)
            
        if os.path.exists(f'{folder}/5_newick')!=True:                  
            subprocess.run([f'touch {folder}/5_newick'], shell=True)


        if os.path.exists(f'{folder}/5_serial')!=True:    
            import re
            with open(f"{folder}/5_newick2",'w') as f2:
                with open(f"{folder}/5_newick",'r') as f:
                    text=f.read().strip()
                cleaned_text = re.sub(r'\[.*?\]', '', text)
                f2.write(cleaned_text)
            from ete3 import Tree
            # Read Newick file
            tree = Tree(f"{folder}/5_newick2")
            # Get sample order (tip labels)
            tip_order = [leaf.name for leaf in tree.get_leaves()]
            
            
            with open(f"{folder}/5_serial",'w') as f:
                f.write('seq_index\tserial\n')
                ii=0
                for one in tip_order:
                    ii+=1
                    f.write(f"{one}\t{ii}\n")
        
        if os.path.exists(f'{folder}/6_info')!=True or 1==1:  
            # Without cd-hit, each sequence is its own representative
            dict_index_represent = {}
            with open(f"{folder}/0_input.info",'r') as f:
                next(f)
                for line in f:
                    eachline_arr = line.strip().split('\t')
                    ID = eachline_arr[0]
                    dict_index_represent[ID] = ID  # Each sequence represents itself
            
            dict_representindex_serial = {}        
            with open(f"{folder}/5_serial",'r') as f:
                next(f)
                for line in f:
                    eachline_arr = line.strip().split('\t')
                    representindex, serial = eachline_arr
                    dict_representindex_serial[representindex] = serial
            
            dict_name_mafft={}        
            with open (f"{folder}/2_mafft",'r') as f:
                for line in f:
                    eachline=line.strip()
                    if eachline[0]==">":name=eachline[1:];dict_name_mafft[name]=''
                    else:dict_name_mafft[name]+=eachline
          
          
            result_list = []     
            with open(f"{folder}/0_input.info",'r') as f:
                next(f)
                for line in f:
                    eachline_arr = line.strip().split('\t')
                    ID, sample, query_name, chromosome, strand, start, end, percent, pvalue, bitscore, seq = eachline_arr
                    represent_index = dict_index_represent[ID]
                    if ID not in dict_representindex_serial:print(f'{ID} missing from the phylogenetic tree, skipping');continue ;
                    representindex_serial = int(dict_representindex_serial[ID])
                    end = int(end)
                    result_list.append([ID, represent_index, representindex_serial, sample, query_name, chromosome, strand, start, end, percent, pvalue, bitscore, seq])
  
            
            sorted_list = sorted(result_list, key=lambda x: (x[2], x[5], x[7]), reverse=False)     
            ii = 0
            with open(f"{folder}/6_input.info",'w') as f:
                f.write('ID\trepresent_index\trepresentindex_serial\tserial_revise\tsample\thap\tquery_name\tchromosome\tstrand\tstart\tend\tpercent\tpvalue\tbitscore\tseq\tmafft\n')
                for one in sorted_list:
                    ID, represent_index, representindex_serial, sample, query_name, chromosome, strand, start, end, percent, pvalue, bitscore, seq = one
                    if "hap1" in sample: hap = 'hap1'
                    elif "hap2" in sample: hap = 'hap2'
                    ii += 1
                    mafft=dict_name_mafft[ID]
                    f.write(f"{ID}\t{represent_index}\t{representindex_serial}\t{ii}\t{sample}\t{hap}\t{query_name}\t{chromosome}\t{strand}\t{start}\t{end}\t{percent}\t{pvalue}\t{bitscore}\t{seq}\t{mafft}\n")

            subprocess.run([f'rm -rf {folder}/6_fasta_by_chr'], shell=True)    
            subprocess.run([f'mkdir {folder}/6_fasta_by_chr'], shell=True)                    
            sorted_list2 = sorted(result_list, key=lambda x: (x[3], x[5], int(x[7])), reverse=False)   
            for one in sorted_list2:
                ID, represent_index, representindex_serial, sample, query_name, chromosome, strand, start, end, percent, pvalue, bitscore, seq = one
                if "hap1" in sample: hap = 'hap1'
                elif "hap2" in sample: hap = 'hap2'
                mafft=dict_name_mafft[ID]
                with open(f"{folder}/6_fasta_by_chr/{hap}_{chromosome}",'a') as f2:
                    f2.write(f">{ID}\n{mafft}\n")
            
        if os.path.exists(f'{folder}/7_plotinput')!=True or 1==1:
            with open(f"{folder}/7_plotinput",'w') as f2:
                f2.write(f"hap\tchromosome\tID\tserial_revise\txmin\txmax\tymin\tymax\n")
                with open(f"{folder}/6_input.info",'r') as f:
                    next(f)
                    for line in f:
                        eachline_arr = line.strip().split('\t')
                        ID, represent_index, representindex_serial, serial_revise, sample, hap, query_name, chromosome, strand, start, end, percent, pvalue, bitscore, seq,mafft = eachline_arr
                        if "hap1" in hap: ymin = 0; ymax = 0.5
                        elif "hap2" in hap: ymin = -0.5; ymax = 0
                        xmin = int(serial_revise)-1
                        xmax = serial_revise
                        chromosome_num = int(chromosome.split('Chr')[-1])
                        ymin = ymin - chromosome_num
                        ymax = ymax - chromosome_num
                        f2.write(f"{hap}\t{chromosome}\t{ID}\t{serial_revise}\t{xmin}\t{xmax}\t{ymin}\t{ymax}\n")
            
            R_txt = r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file = read.table('7_plotinput', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
        color_values = c("Chr1" = "red",  
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
                          "Chr19" = "#FF8C00")
        # Create plot object
        p <- ggplot() +
            geom_rect(data = input_file, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = chromosome)) +
            scale_fill_manual(values = color_values, drop = FALSE) +
            theme_classic() +         
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
          
          # Save as PDF
          pdf(file = paste0('7_plot', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
          print(p)
          dev.off()
        }
            '''
            with open(f'{folder}/plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)  
            new_directory = folder
            os.chdir(new_directory)
            subprocess.run(['Rscript plot.R'], shell=True)    
            os.chdir('../../')                            
            
            
            
            
            
            
            
            
    with Pool(processes=thread) as pool:
        # Use imap to get results one by one
        for i, result in enumerate(pool.imap(run_step, name_list), start=1):
            # Results can be processed here, e.g., stored or printed
            progress = (i / len(name_list)) * 100
            sys.stdout.write(f"\rProgress: {progress:.2f}%")
            sys.stdout.flush()



name_list=['PN40024','V050','V060','V008','V126','V031']
if argv1=="part_pos2cluster":

    for one in name_list:
        subprocess.run([f'mkdir -p ./chr2soloLTRs/part_pos2cluster/{one}'], shell=True)
        ##
        if os.path.exists(f"./chr2soloLTRs/part_pos2cluster/{one}/1_info")!=True :
            from Bio import pairwise2
            from Bio.Seq import Seq
            results=[]
            seq_set=set()
            with open (f"./chr2soloLTRs/3_filter.all",'r') as f:
                next(f)
                for line in f: 
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    ID,sample,query_name,chromosome,strand,start,end,percent,pvalue,bitscore,seq=eachline_arr
                    if one in sample:
                        seq_set.add(seq)
                        results.append(eachline_arr)
            sorted_results = sorted(results, key=lambda x: (x[1],x[3],int(x[5])), reverse=False)   
            
            
            dict_name_seqs={}
            with open(f"./chr2soloLTRs/14_simple_long/8_goodnodes_2",'r') as f:
                for line in f :
                    eachline=line.strip()
                    if eachline[0]=='>': name =eachline[1:];dict_name_seqs[name]='';
                    else:dict_name_seqs[name]+=eachline
    
                
            def run_step(oneseq):
                best_match = None
                best_score = -1
                
                for node_name, node_seq in dict_name_seqs.items():
                    # Convert to Seq object (optional, but recommended for Biopython)
                    seq1 = Seq(oneseq)
                    seq2 = Seq(node_seq)
                    
                    # Perform global alignment (can replace with pairwise2.align.localxx for local alignment)
                    alignments = pairwise2.align.localxx(seq1, seq2)
                    
                    # Take the best alignment result
                    best_alignment = alignments[0]
                    score = best_alignment.score  # Alignment score
                    similarity = score / max(len(seq1), len(seq2))  # Normalized similarity (0~1)
                    
                    if similarity > best_score:
                        best_score = similarity
                        best_match = node_name
                
                return oneseq, best_match, best_score * 100  # Convert to percentage
            
            seq_set_list=list(seq_set)   
            dict_seq_infos={}   
            length=len(seq_set_list)
            with Pool(processes=thread) as pool:
                # Use imap to get results one by one
                for i, result in enumerate(pool.imap(run_step, seq_set_list), start=1):
                    # Results can be processed here, e.g., stored or printed
                    oneseq, best_match, best_score=result
                    dict_seq_infos[oneseq]=[best_match, best_score]
                    progress = (i / length) * 100
                    sys.stdout.write(f"\rProgress:{i}/{length}\t\t{progress:.2f}%")
                    sys.stdout.flush()  
            
            
            
            with open (f"./chr2soloLTRs/part_pos2cluster/{one}/1_info",'w') as f2:
                f2.write(f"ID\tsample\tquery_name\tchromosome\tstrand\tstart\tend\tpercent\tpvalue\tbitscore\tseq\n")
                for oneline in sorted_results:
                    ID,sample,query_name,chromosome,strand,start,end,percent,pvalue,bitscore,seq=oneline
                    if one in sample:
                        best_match,best_score=dict_seq_infos[seq]
                        f2.write(f"{ID}\t{sample}\t{query_name}\t{chromosome}\t{strand}\t{start}\t{end}\t{percent}\t{pvalue}\t{bitscore}\t{seq}\t{best_match}\t{best_score}\n")
                        
        if os.path.exists(f"./chr2soloLTRs/part_pos2cluster/{one}/2_plot.pdf")!=True or 1==1:        
            with open(f"./chr2soloLTRs/part_pos2cluster/{one}/2_plot.info",'w') as f2:
                 f2.write(f"ID\tsample\tchromosome\tchromosome_num\tstrand\tstart\tend\tbest_match\tbest_score\n")
                 with open(f"./chr2soloLTRs/part_pos2cluster/{one}/1_info",'r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        ID,sample,query_name,chromosome,strand,start,end,percent,pvalue,bitscore,seq,best_match,best_score=eachline_arr
                        chromosome_num=chromosome[3:]
                        f2.write(f"{ID}\t{sample}\t{chromosome}\t{chromosome_num}\t{strand}\t{start}\t{end}\t{best_match}\t{best_score}\n")
            R_txt = r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file = read.table('2_plot.info', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
        color_values2 = c(
        "Node2"="red",
"Node8"="red",
"Node10"="red",
"Node12"="red",
"Node16"="red",
"Node18"="red",
"Node21"="red",
"Node23"="orange",
"Node60"="orange",
"Node68"="yellow",
"Node71"="yellow",
"Node87"="green",
"Node141"="blue",
"Node178"="#cc66ff",
"Node184"="#9933ff",
"Node190"="#6600cc",
"Node200"="#666699"
)
color_values = c(
        "Node2"="red",
"Node8"="#8CC63F",
"Node10"="#CC0066",
"Node12"="#148EFF",
"Node16"="#B2DF8A",
"Node18"="#008000",
"Node21"="#FB9A99",
"Node23"="#00C08F",
"Node60"="#005FAF",
"Node68"="#FFCC00",
"Node71"="#00ccff",
"Node87"="#6A3D9A",
"Node141"="#FFFF00",
"Node178"="#A0522D",
"Node184"="#9966CC",
"Node190"="#ff4dff",
"Node200"="#ff6699"
)
        
        # Create plot object
        p <- ggplot() +
            geom_rect(data = input_file, aes(xmin = start, xmax = end, ymin = -chromosome_num-0.5, ymax = -chromosome_num, fill = best_match)) +
            scale_fill_manual(values = color_values, drop = FALSE) +
            theme_classic() +     
            facet_wrap( ~ sample) +    
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
          
          # Save as PDF
          pdf(file = paste0('2_plot', ".pdf"), width = 60/ 2.54, height = 20 / 2.54)
          print(p)
          dev.off()
        }
            '''
            folder=f"./chr2soloLTRs/part_pos2cluster/{one}/"
            with open(f'{folder}/plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)  
            new_directory = folder
            os.chdir(new_directory)
            subprocess.run(['Rscript plot.R'], shell=True)    
            subprocess.run([f'mv 2_plot.pdf 2_plot_{one}.pdf '], shell=True)    
            os.chdir('../../../')                 
                        

            
            




name_list=[
'PN40024:Chr1:15154005-15435737',
#'PN40024_hap1:Chr3:14353397-14553397',
'PN40024:Chr2:12939361-13300832',
'PN40024:Chr4:13291174-13483935',
'PN40024:Chr5:14008287-14357682',
'PN40024:Chr6:10606684-10966038',
'PN40024:Chr7:13458101-13762988', 
'PN40024:Chr8:6679384-7071344',
'PN40024:Chr9:15113362-15449399',
'PN40024:Chr13:12020240-12180415',
'PN40024:Chr14:15144034-15480228',
'PN40024:Chr17:15221832-15560758',
#'PN40024:Chr19:15995936-16195936'
'PN40024:Chr18:16390000-16620000']


 
name_list=[ 
     #'V126.hap2:Chr1:14384736-14863897',
    'V126.hap2:Chr2:14592428-15130566',
    'V126.hap2:Chr3:15258845-15877811',
    'V126.hap2:Chr4:13956925-14376580',
    'V126.hap2:Chr5:15078558-15669106',
    'V126.hap2:Chr6:12010223-12569057',
    #'V126.hap2:Chr7:14847465-15293715',
    'V126.hap2:Chr8:7790307-8223339',
    #'V126.hap2:Chr9:17432607-17687234',  
    'V126.hap2:Chr12:12175803-12899879',
    #'V126.hap2:Chr13:13781732-14322746',
    'V126.hap2:Chr14:14894034-15529422',
    #'V126.hap2:Chr17:14300851-14915357',  #It seems there really is no LIRSat
    'V126.hap2:Chr18:16444108-17145390',
    'V126.hap2:Chr19:17585851-18207892']
    

 

 
name_list=[
'PN40024:Chr1:15154005-15435737',
'PN40024:Chr13:12020240-12180415',
'PN40024:Chr7:13458101-13762988', 
'PN40024:Chr14:15144034-15480228',
'PN40024:Chr17:15221832-15560758',
'PN40024:Chr2:12939361-13300832',
'PN40024_hap1:Chr3:14353397-14553397',
'PN40024:Chr4:13291174-13483935',
'PN40024:Chr5:14008287-14357682',
'PN40024:Chr6:10606684-10966038',
'PN40024:Chr8:6679384-7071344',
'PN40024:Chr9:15113362-15449399',
'PN40024:Chr19:15995936-16195936',
'PN40024:Chr18:16390000-16620000']




name_list=[
 'WoollyGrape_hap1:Chr1:14788929-15147386',
 'WoollyGrape_hap1:Chr2:15400993-15782030',
 'WoollyGrape_hap1:Chr3:15620654-16005187',
 'WoollyGrape_hap1:Chr4:13740850-14096219',
 'WoollyGrape_hap1:Chr5:14150576-14632756',
 'WoollyGrape_hap1:Chr6:12816579-13290046',
 'WoollyGrape_hap1:Chr7:14618553-14871377',
 'WoollyGrape_hap1:Chr8:8841760-9256061',
 'WoollyGrape_hap1:Chr9:17693900-18101280',
 'WoollyGrape_hap1:Chr12:11849656-12381808',
 'WoollyGrape_hap1:Chr13:13301853-13566766',
 'WoollyGrape_hap1:Chr14:16015135-16374515',
 'WoollyGrape_hap1:Chr17:14350851-14865357',
 'WoollyGrape_hap1:Chr18:16535861-16927372',
 'WoollyGrape_hap1:Chr19:18996536-19383157']

name_list=[
'V124.hap2:Chr4:13428189-13640340',
'V124.hap2:Chr4:13606514-13849203',
'V081.hap2:Chr4:14489253-14857555',
'V079.hap1:Chr4:14158159-14697513',

'V123.hap1:Chr9:18432955-18697537',
'WoollyGrape_hap1:Chr9:17650000-17920000',
'V008.hap1:Chr9:16448098-17004837',
'V069.hap1:Chr9:15178574-15614639'  ,         
            

'V037.hap2:Chr7:14382973-14710457',
'ThompsonSeedless_hap1:Chr7:12585184-12913100',
'VHP-T2T.hap1:Chr7:14667734-15017556',
'V060.hap2:Chr7:13097624-13599746',

'V107.hap1:Chr14:16034600-16408612',
'V088.hap1:Chr14:15032479-15468677',
'V048.hap1:Chr14:13988466-14366294',
'V040.hap2:Chr14:14716115-15301017',

'V032.hap2:Chr1:15428505-15829144',
'V124.hap2:Chr1:14643236-15086673',
'V112.hap1:Chr1:14478774-14928728',
'V008.hap1:Chr1:14701815-15169418',

'V081.hap1:Chr18:16300000-16520000',
'V066.hap2:Chr18:16280000-16570000',
'V031.hap2:Chr18:16210000-16650000',
'V055.hap1:Chr18:15770000-16270000']     

name_list=[
'PN40024:Chr1:15154005-15435737',
#'PN40024_hap1:Chr3:14353397-14553397',
'PN40024:Chr2:12939361-13300832',
'PN40024:Chr4:13291174-13483935',
'PN40024:Chr5:14008287-14357682',
'PN40024:Chr6:10606684-10966038',
'PN40024:Chr7:13458101-13762988', 
'PN40024:Chr8:6679384-7071344',
'PN40024:Chr9:15113362-15449399',
'PN40024:Chr13:12020240-12180415',
'PN40024:Chr14:15144034-15480228',
'PN40024:Chr17:15221832-15560758',
#'PN40024:Chr19:15995936-16195936'
'PN40024:Chr18:16390000-16620000']

if argv1=="part_pos2cluster_plus":
    for one_part in name_list:
        one,t_chromosome,t_pos=one_part.split(':')
        t_pos1,t_pos2=t_pos.split('-')
        t_pos1=int(t_pos1)
        t_pos2=int(t_pos2)
        subprocess.run([f'mkdir -p ./chr2soloLTRs/part_pos2cluster_plus/{one_part}'], shell=True)
        ##
        if os.path.exists(f"./chr2soloLTRs/part_pos2cluster_plus/{one_part}/1_info")!=True :
            from Bio import pairwise2
            from Bio.Seq import Seq
            results=[]
            seq_set=set()
            with open (f"./chr2soloLTRs/3_filter.all",'r') as f:
                next(f)
                for line in f: 
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    ID,sample,query_name,chromosome,strand,start,end,percent,pvalue,bitscore,seq=eachline_arr
                    if one ==sample:
                        if chromosome!=t_chromosome:continue
                        if int(start)>t_pos1 and int(end)<t_pos2:
                            seq_set.add(seq)
                            results.append(eachline_arr)
            sorted_results = sorted(results, key=lambda x: (x[1],x[3],int(x[5])), reverse=False)   
            
            
            dict_name_seqs={}
            with open(f"./chr2soloLTRs/14_simple_long/8_goodnodes_2",'r') as f:
                for line in f :
                    eachline=line.strip()
                    if eachline[0]=='>': name =eachline[1:];dict_name_seqs[name]='';
                    else:dict_name_seqs[name]+=eachline
    
                
            def run_step(oneseq):
                best_match = None
                best_score = -1
                
                for node_name, node_seq in dict_name_seqs.items():
                    # Convert to Seq object (optional, but recommended for Biopython)
                    seq1 = Seq(oneseq)
                    seq2 = Seq(node_seq)
                    
                    # Perform global alignment (can replace with pairwise2.align.localxx for local alignment)
                    alignments = pairwise2.align.localxx(seq1, seq2)
                    
                    # Take the best alignment result
                    best_alignment = alignments[0]
                    score = best_alignment.score  # Alignment score
                    similarity = score / max(len(seq1), len(seq2))  # Normalized similarity (0~1)
                    
                    if similarity > best_score:
                        best_score = similarity
                        best_match = node_name
                
                return oneseq, best_match, best_score * 100  # Convert to percentage
            
            seq_set_list=list(seq_set)   
            dict_seq_infos={}   
            length=len(seq_set_list)
            with Pool(processes=thread) as pool:
                # Use imap to get results one by one
                for i, result in enumerate(pool.imap(run_step, seq_set_list), start=1):
                    # Results can be processed here, e.g., stored or printed
                    oneseq, best_match, best_score=result
                    dict_seq_infos[oneseq]=[best_match, best_score]
                    progress = (i / length) * 100
                    sys.stdout.write(f"\rProgress:{i}/{length}\t\t{progress:.2f}%")
                    sys.stdout.flush()  
            
            
            
            with open (f"./chr2soloLTRs/part_pos2cluster_plus/{one_part}/1_info",'w') as f2:
                f2.write(f"ID\tsample\tquery_name\tchromosome\tstrand\tstart\tend\tpercent\tpvalue\tbitscore\tseq\n")
                for oneline in sorted_results:
                    ID,sample,query_name,chromosome,strand,start,end,percent,pvalue,bitscore,seq=oneline
                    if one in sample:
                        best_match,best_score=dict_seq_infos[seq]
                        f2.write(f"{ID}\t{sample}\t{query_name}\t{chromosome}\t{strand}\t{start}\t{end}\t{percent}\t{pvalue}\t{bitscore}\t{seq}\t{best_match}\t{best_score}\n")
                        
        if os.path.exists(f"./chr2soloLTRs/part_pos2cluster_plus/{one_part}/2_plot.pdf")!=True :        
            with open(f"./chr2soloLTRs/part_pos2cluster_plus/{one_part}/2_plot.info",'w') as f2:
                 f2.write(f"ID\tsample\tchromosome\tchromosome_num\tstrand\tstart\tend\tbest_match\tbest_score\n")
                 with open(f"./chr2soloLTRs/part_pos2cluster_plus/{one_part}/1_info",'r') as f:
                    next(f)
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        ID,sample,query_name,chromosome,strand,start,end,percent,pvalue,bitscore,seq,best_match,best_score=eachline_arr
                        chromosome_num=chromosome[3:]
                        f2.write(f"{ID}\t{sample}\t{chromosome}\t{chromosome_num}\t{strand}\t{start}\t{end}\t{best_match}\t{best_score}\n")
            R_txt = r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file = read.table('2_plot.info', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
        color_values2 = c(
        "Node2"="red",
"Node8"="red",
"Node10"="red",
"Node12"="red",
"Node16"="red",
"Node18"="red",
"Node21"="red",
"Node23"="orange",
"Node60"="orange",
"Node68"="yellow",
"Node71"="yellow",
"Node87"="green",
"Node141"="blue",
"Node178"="#cc66ff",
"Node184"="#9933ff",
"Node190"="#6600cc",
"Node200"="#666699"
)
color_values = c(
        "Node2"="red",
"Node8"="#8CC63F",
"Node10"="#CC0066",
"Node12"="#148EFF",
"Node16"="#B2DF8A",
"Node18"="#008000",
"Node21"="#FB9A99",
"Node23"="#00C08F",
"Node60"="#005FAF",
"Node68"="#FFCC00",
"Node71"="#00ccff",
"Node87"="#6A3D9A",
"Node141"="#FFFF00",
"Node178"="#A0522D",
"Node184"="#9966CC",
"Node190"="#ff4dff",
"Node200"="#ff6699"
)
        
        # Create plot object
        p <- ggplot() +
            geom_rect(data = input_file, aes(xmin = start, xmax = end, ymin = -chromosome_num-0.5, ymax = -chromosome_num, fill = best_match)) +
            scale_fill_manual(values = color_values, drop = FALSE) +
            theme_classic() +     
            facet_wrap( ~ sample) +    
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
          
          # Save as PDF
          pdf(file = paste0('2_plot', ".pdf"), width = 60/ 2.54, height = 20 / 2.54)
          print(p)
          dev.off()
        }
            '''
            folder=f"./chr2soloLTRs/part_pos2cluster_plus/{one_part}/"
            with open(f'{folder}/plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)  
            new_directory = folder
            os.chdir(new_directory)
            subprocess.run(['Rscript plot.R'], shell=True)    
            subprocess.run([f'mv 2_plot.pdf 2_plot_{one_part}.pdf '], shell=True)    
            os.chdir('../../../')                 
    #if argv1=="part_pos2cluster_plus2":
    print('Summarize the results of part_pos2cluster_plus')
    if os.path.exists(f'./chr2soloLTRs/part_nocdhit//sum')!=True :
        subprocess.run(['mkdir -p ./chr2soloLTRs/part_pos2cluster_plus/sum'], shell=True)
    with open(f"./chr2soloLTRs/part_pos2cluster_plus/sum/all.info",'w') as f2:
        f2.write(f"serial\tID\tsample\tchromosome\tchromosome_num\tstrand\tstart\tend\tbest_match\tbest_score\tregion_start\tregion_end\tregion_len\tstart_delta\tend_delta\n")
        iii=0    
        for one_part in name_list:
            iii+=1
            one,t_chromosome,t_pos=one_part.split(':')
            t_pos1,t_pos2=t_pos.split('-')
            t_pos1=int(t_pos1)
            t_pos2=int(t_pos2)
            length=(t_pos2-t_pos1)+1
            with open(f"./chr2soloLTRs/part_pos2cluster_plus/{one_part}/2_plot.info",'r') as f:
                next(f)
                for line in f:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    
                    ID,sample,chromosome,chromosome_num,strand,start,end,best_match,best_score=eachline_arr
                    start_delta=int(start)-t_pos1
                    end_delta=int(end)-t_pos1
                    f2.write(f"{iii}\t{eachline}\t{t_pos1}\t{t_pos2}\t{length}\t{start_delta}\t{end_delta}\n")
    R_txt = r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file = read.table('all.info', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
        color_values2 = c(
        "Node2"="red",
"Node8"="red",
"Node10"="red",
"Node12"="red",
"Node16"="red",
"Node18"="red",
"Node21"="red",
"Node23"="orange",
"Node60"="orange",
"Node68"="yellow",
"Node71"="yellow",
"Node87"="green",
"Node141"="blue",
"Node178"="#cc66ff",
"Node184"="#9933ff",
"Node190"="#6600cc",
"Node200"="#666699"
)
color_values = c(
        "Node2"="red",
"Node8"="#8CC63F",
"Node10"="#CC0066",
"Node12"="#148EFF",
"Node16"="#B2DF8A",
"Node18"="#008000",
"Node21"="#FB9A99",
"Node23"="#00C08F",
"Node60"="#005FAF",
"Node68"="#FFCC00",
"Node71"="#00ccff",
"Node87"="#6A3D9A",
"Node141"="#FFFF00",
"Node178"="#A0522D",
"Node184"="#9966CC",
"Node190"="#ff4dff",
"Node200"="#ff6699"
)
        
        # Create plot object
        p <- ggplot() +
            geom_rect(data = input_file, aes(xmin = 0, xmax = region_len, ymin = -serial-0.5, ymax = -serial, fill = '#e6e6e6')) +
            geom_rect(data = input_file, aes(xmin = start_delta, xmax = end_delta, ymin = -serial-0.5, ymax = -serial, fill = best_match)) +
            scale_fill_manual(values = color_values, drop = FALSE) +
            theme_classic() +     
            #facet_wrap( ~ sample) +    
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
          
          # Save as PDF
          pdf(file = paste0('2_plot', ".pdf"), width = 60/ 2.54, height = 20 / 2.54)
          print(p)
          dev.off()
        }
            '''
    folder=f"./chr2soloLTRs/part_pos2cluster_plus/sum/"
    with open(f'{folder}/plot.R','w',encoding='utf-8') as f:
        f.write(R_txt)  
    new_directory = folder
    os.chdir(new_directory)
    subprocess.run(['Rscript plot.R'], shell=True)    
    os.chdir('../../../')  

    print('Plot forward and reverse strand diagram')
    R_txt = r'''
        library(ggplot2)
        library(dplyr)
        
        print("")
        {
          input_file = read.table('all.info', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
        color_values = c(
        "+"="red",
"-"="blue"
)
        
        # Create plot object
        p <- ggplot() +
            geom_rect(data = input_file, aes(xmin = 0, xmax = region_len, ymin = -serial-0.5, ymax = -serial, fill = '#e6e6e6')) +
            geom_rect(data = input_file, aes(xmin = start_delta, xmax = end_delta, ymin = -serial-0.5, ymax = -serial, fill = strand)) +
            scale_fill_manual(values = color_values, drop = FALSE) +
            theme_classic() +     
            #facet_wrap( ~ sample) +    
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
          
          # Save as PDF
          pdf(file = paste0('2_plot_strand', ".pdf"), width = 60/ 2.54, height = 20 / 2.54)
          print(p)
          dev.off()
        }
            '''
    folder=f"./chr2soloLTRs/part_pos2cluster_plus/sum/"
    with open(f'{folder}/plot.R','w',encoding='utf-8') as f:
        f.write(R_txt)  
    new_directory = folder
    os.chdir(new_directory)
    subprocess.run(['Rscript plot.R'], shell=True)    
    os.chdir('../../../')  



if argv1=="part_pos2cluster_plus_blast":
    
    # ========== Read color configuration file ==========
    color_config_file = "./chr2soloLTRs/14_simple_long/5_serial_color"
    seq_color_map = {}      # seq_index -> color
    seq_cluster_map = {}    # seq_index -> cluster
    seq_serial_map = {}     # seq_index -> serial
    
    if os.path.exists(color_config_file):
        print(f"Reading color configuration file: {color_config_file}")
        with open(color_config_file, 'r') as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    seq_index, serial, color, cluster = parts[0], parts[1], parts[2], parts[3]
                    seq_color_map[seq_index] = color
                    seq_cluster_map[seq_index] = cluster
                    seq_serial_map[seq_index] = serial
        print(f"  Read color annotations for {len(seq_color_map)} sequences")
    else:
        print(f"Warning: Color configuration file {color_config_file} does not exist, using default colors")
    
    # ========== Build representative sequence BLAST database ==========
    print("\nBuilding representative sequence BLAST database...")
    repr_fasta = "./chr2soloLTRs/14_simple_long/1_cdhit"
    blast_db = "./chr2soloLTRs/14_simple_long/1_cdhit_db"
    
    # Read list of representative sequence IDs (for later color lookup)
    repr_ids = []
    with open(repr_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                repr_id = line[1:].strip()
                repr_ids.append(repr_id)
    print(f"Number of representative sequences: {len(repr_ids)}")
    
    # Build BLAST database (if it doesn't exist)
    if not os.path.exists(f"{blast_db}.nhr"):
        subprocess.run([f"makeblastdb -in {repr_fasta} -dbtype nucl -out {blast_db}"], shell=True)
        print("BLAST database built")
    else:
        print("BLAST database already exists")
    
    # ========== Process each region ==========
    for one_part in name_list:
        one, t_chromosome, t_pos = one_part.split(':')
        t_pos1, t_pos2 = t_pos.split('-')
        t_pos1 = int(t_pos1)
        t_pos2 = int(t_pos2)
        
        output_dir = f"./chr2soloLTRs/part_pos2cluster_plus_blast/{one_part}"
        subprocess.run([f'mkdir -p {output_dir}'], shell=True)
        
        info_file = f"{output_dir}/1_info"
        if not os.path.exists(info_file) or 1==1:
            
            # Collect sequences for this region
            results = []
            seq_set = set()
            seq_to_record = {}  # seq -> full record
            
            print(f"\nProcessing region: {one_part}")
            with open("./chr2soloLTRs/3_filter.all", 'r') as f:
                next(f)
                for line in f:
                    eachline = line.strip()
                    eachline_arr = eachline.split('\t')
                    ID, sample, query_name, chromosome, strand, start, end, percent, pvalue, bitscore, seq = eachline_arr
                    if one == sample:
                        if chromosome != t_chromosome:
                            continue
                        if int(start) > t_pos1 and int(end) < t_pos2:
                            if seq not in seq_set:
                                seq_set.add(seq)
                                seq_to_record[seq] = eachline_arr
                            results.append(eachline_arr)
            
            if len(seq_set) == 0:
                print(f"  Warning: No sequences found in region {one_part}")
                continue
            
            print(f"  Found {len(seq_set)} unique sequences")
            
            # Write query sequence FASTA
            query_fasta = f"{output_dir}/query.fasta"
            with open(query_fasta, 'w') as f:
                for i, seq in enumerate(seq_set):
                    f.write(f">query_{i}\n{seq}\n")
            
            # Run BLAST
            blast_out = f"{output_dir}/blast_results.tsv"
            subprocess.run([
                f"blastn -query {query_fasta} -db {blast_db} "
                f"-out {blast_out} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' "
                f"-max_target_seqs 1 -num_threads {thread}"
            ], shell=True)
            print(f"  BLAST completed")
            
            # Parse BLAST results
            seq_best_match = {}  # seq -> (best_match_id, bitscore, pident)
            with open(blast_out, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 12:
                        qseqid = parts[0]
                        sseqid = parts[1]
                        bitscore = float(parts[11])
                        pident = float(parts[2])
                        # Extract original sequence
                        query_idx = int(qseqid.split('_')[1])
                        seq = list(seq_set)[query_idx]
                        if seq not in seq_best_match:
                            seq_best_match[seq] = (sseqid, bitscore, pident)
            
            # Write results file
            with open(info_file, 'w') as f2:
                f2.write("ID\tsample\tquery_name\tchromosome\tstrand\tstart\tend\tpercent\tpvalue\tbitscore\tseq\tbest_match\tbest_score\tbest_pident\tcluster\tserial\tcolor\n")
                for oneline in results:
                    ID, sample, query_name, chromosome, strand, start, end, percent, pvalue, bitscore, seq = oneline
                    if seq in seq_best_match:
                        best_match, best_bitscore, best_pident = seq_best_match[seq]
                        # Get color and cluster from color configuration file
                        color = seq_color_map.get(best_match, "#808080")  # Default gray
                        cluster = seq_cluster_map.get(best_match, "Unknown")
                        serial = seq_serial_map.get(best_match, "0")
                    else:
                        best_match, best_bitscore, best_pident = "None", 0, 0
                        color, cluster, serial = "#808080", "Unknown", "0"
                    
                    f2.write(f"{ID}\t{sample}\t{query_name}\t{chromosome}\t{strand}\t{start}\t{end}\t{percent}\t{pvalue}\t{bitscore}\t{seq}\t{best_match}\t{best_bitscore}\t{best_pident}\t{cluster}\t{serial}\t{color}\n")
            
            print(f"  Results saved to: {info_file}")
        
        # ========== Generate plot file ==========
        plot_info_file = f"{output_dir}/2_plot.info"
        if not os.path.exists(f"{output_dir}/2_plot.pdf") or 1==1:
            
            with open(plot_info_file, 'w') as f2:
                f2.write("ID\tsample\tchromosome\tchromosome_num\tstrand\tstart\tend\tbest_match\tserial\tcluster\tcolor\n")
                with open(info_file, 'r') as f:
                    next(f)
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 17:
                            ID, sample, query_name, chromosome, strand, start, end, percent, pvalue, bitscore, seq, best_match, best_score, best_pident,  cluster, serial,color = parts
                            chromosome_num = chromosome[3:]
                            f2.write(f"{ID}\t{sample}\t{chromosome}\t{chromosome_num}\t{strand}\t{start}\t{end}\t{best_match}\t{serial}\t{cluster}\t{color}\n")
            
            # Generate R plotting script (using colors read from configuration file)
            R_txt = r'''
library(ggplot2)
library(dplyr)

input_file <- read.table('2_plot.info', 
                         skip = 0, 
                         header = TRUE, 
                         stringsAsFactors = FALSE, 
                         check.names = FALSE, 
                         sep = '\t',
                         strip.white = TRUE,comment.char = "")

# Convert numeric columns
input_file$start <- as.numeric(input_file$start)
input_file$end <- as.numeric(input_file$end)
input_file$chromosome_num <- as.numeric(input_file$chromosome_num)



# Use the color column directly, not through cluster
p <- ggplot() +
    geom_rect(data = input_file, 
              aes(xmin = start, xmax = end, 
                  ymin = -chromosome_num - 0.5, ymax = -chromosome_num,
                  fill = color))+scale_fill_identity()  + 
    theme_classic()

    ##facet_wrap(~ sample, scales = "free") +    
    #theme(
    #    axis.ticks.y = element_blank(),
    #    axis.text.y = element_blank(),
    #    axis.text.x = element_blank()#,
     #   #legend.position = "none"  # No legend needed
    #) 

pdf(file = '2_plot.pdf', width = 20/2.54, height = 20/2.54)
print(p)
dev.off()
'''
            
            with open(f"{output_dir}/plot.R", 'w', encoding='utf-8') as f:
                f.write(R_txt)
            
            os.chdir(output_dir)
            #subprocess.run(['Rscript plot.R'], shell=True)   
            #subprocess.run([f'mv 2_plot.pdf 2_plot_{one_part}.pdf'], shell=True)
            os.chdir('../../../')
            print(f"  Plotting completed: {output_dir}/2_plot_{one_part}.pdf")
    
    # ========== Summarize all regions ==========
    print("\n" + "="*60)
    print("Summarizing all region results...")
    print("="*60)
    
    sum_dir = "./chr2soloLTRs/part_pos2cluster_plus_blast/sum"
    subprocess.run([f'mkdir -p {sum_dir}'], shell=True)
    
    with open(f"{sum_dir}/all.info", 'w') as f2:
        f2.write("serial\tID\tsample\tchromosome\tchromosome_num\tstrand\tstart\tend\tbest_match\tserial_in_tree\tcluster\tcolor\tregion_start\tregion_end\tregion_len\tstart_delta\tend_delta\n")
        iii = 0
        for one_part in name_list:
            iii += 1
            one, t_chromosome, t_pos = one_part.split(':')
            t_pos1, t_pos2 = t_pos.split('-')
            t_pos1 = int(t_pos1)
            t_pos2 = int(t_pos2)
            length = (t_pos2 - t_pos1) + 1
            
            plot_info_file = f"./chr2soloLTRs/part_pos2cluster_plus_blast/{one_part}/2_plot.info"
            if os.path.exists(plot_info_file) or 1==1:
                with open(plot_info_file, 'r') as f:
                    next(f)
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) == 11:
                            ID, sample, chromosome, chromosome_num, strand, start, end, best_match, serial,cluster,color = parts
                            start_delta = int(start) - t_pos1
                            end_delta = int(end) - t_pos1
                            f2.write(f"{iii}\t{ID}\t{sample}\t{chromosome}\t{chromosome_num}\t{strand}\t{start}\t{end}\t{best_match}\t{serial}\t{cluster}\t{color}\t{t_pos1}\t{t_pos2}\t{length}\t{start_delta}\t{end_delta}\n")
    
    # Generate summary plot (colored by cluster)
    R_summary = '''
library(ggplot2)
library(dplyr)

input_file <- read.table('all.info', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t', comment.char = "")

# Convert cluster to factor
input_file$cluster <- as.factor(input_file$cluster)

# Get unique clusters and their colors
unique_clusters <- unique(input_file$cluster)
cluster_colors <- setNames(
    sapply(unique_clusters, function(cl) as.character(input_file$color[input_file$cluster == cl][1])),
    unique_clusters
)

# Background color (wrap in quotes to avoid # being misinterpreted)
bg_color <- "#e6e6e6"

# Create plot object
p <- ggplot() +
    geom_rect(data = input_file, 
              aes(xmin = 0, xmax = region_len, 
                  ymin = -serial - 0.5, ymax = -serial), 
              fill = bg_color, inherit.aes = FALSE) +
    geom_rect(data = input_file, 
              aes(xmin = start_delta, xmax = end_delta, 
                  ymin = -serial - 0.5, ymax = -serial, 
                  fill = cluster)) +
    scale_fill_manual(values = cluster_colors, drop = FALSE) +
    theme_classic() +     
    theme(
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank()
    ) 

pdf(file = 'summary_plot.pdf', width = 60/2.54, height = 20/2.54)
print(p)
dev.off()

# Draw forward and reverse strand plot
p_strand <- ggplot() +
    geom_rect(data = input_file, 
              aes(xmin = 0, xmax = region_len, 
                  ymin = -serial - 0.5, ymax = -serial), 
              fill = bg_color, inherit.aes = FALSE) +
    geom_rect(data = input_file, 
              aes(xmin = start_delta, xmax = end_delta, 
                  ymin = -serial - 0.5, ymax = -serial, 
                  fill = strand)) +
    scale_fill_manual(values = c("+" = "red", "-" = "blue"), drop = FALSE) +
    theme_classic() +     
    theme(
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank()
    ) 

pdf(file = 'summary_plot_strand.pdf', width = 60/2.54, height = 20/2.54)
print(p_strand)
dev.off()
'''
    
    with open(f"{sum_dir}/plot.R", 'w', encoding='utf-8') as f:
        f.write(R_summary)
    
    os.chdir(sum_dir)
    subprocess.run(['Rscript plot.R'], shell=True)   
    os.chdir('../../../')
    
    print("\n" + "="*60)
    print("Analysis completed!")
    print(f"Results directory: {sum_dir}")
    print("="*60)







   
   
print("\n\n")        
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))      
time_end=timeit.default_timer()
print('All the running time: %.0f Seconds'%(time_end-time_start))