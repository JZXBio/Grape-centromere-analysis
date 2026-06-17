#!/usr/bin/python
# -*- coding: UTF-8 -*-


import sys
argvs=sys.argv
argv_len=len(sys.argv)
if argv_len < 2 or argvs[1] in ("h", "-h", "help", "-help") :
    print ("chr2aligntime_DIY.py-----help:")
    print ("")
    print ("Usage：")
    print ("chr2aligntime_DIY.py  inputfile(containing two fasta files)  ")
    print('Target format example:')
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

subprocess.run([f"mkdir -p chr2aligntime_DIY"], shell=True)


#python chr2aligntime_DIY.py V124.hap2:Chr4:13656514-13726514 V124.hap2:Chr4:13726514-13799203 49b
#python chr2aligntime_DIY.py V124.hap2:Chr4:13478189-13534265 V124.hap2:Chr4:13534265-13590340 49a






time_start_str = time.strftime('%Y%m%d_%H%M%S', time.localtime())



target1=sys.argv[1]            

if argv_len==2:
    output_dir=f'chr2aligntime_DIY/{time_start_str}'
else:
    output_dir=sys.argv[2]

kmer_len=41
    
if  os.path.exists(output_dir)==False:
    subprocess.run([f"mkdir {output_dir}"], shell=True)
    
def reverse_complement(seq):
    # Define complement base mapping dictionary
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
                  'N': 'N', 'n': 'n'}  # Handle unknown bases
    
    # Generate complement sequence and reverse
    rev_comp = ''.join([complement.get(base, base) for base in seq][::-1])
    return rev_comp    
    
pic='no'

    
if 'step0'=='step0':
    print('Extract main sequences')


    print(f"target1:{target1}")
    kk=0;seq1='';seq2=''
    with open(target1,'r') as f:
        for line in f:
            eachline=line.strip()
            if eachline[0]=='>':  kk+=1
            if kk==1:
                seq1+=eachline
            if kk==2:
                seq2+=eachline

    with open(f"{output_dir}/seq1.fa",'w') as f:
        f.write(f'>seq1\n{seq1}')
 
    with open(f"{output_dir}/seq2.fa",'w') as f:
        f.write(f'>seq2\n{seq2}')
        
    
        
    seq1_len=len(seq1)
    seq2=reverse_complement(seq2)
    seq2_len=len(seq2)
    raw_seq1_len=seq1_len
    raw_seq2_len=seq2_len
    print(f"seq1len:{seq1_len}")
    print(f"seq2len:{seq2_len}")   
    with open(f'{output_dir}/seqs12.fa','w') as f:
        f.write(f'>seq1\n{seq1}\n>seq2\n{seq2}')

if 'step1'=='step2':
    print('Plot using kmers 21/41')
    from collections import defaultdict
    import matplotlib.pyplot as plt
     
    def extract_kmers(seq, k):
        return [seq[i:i+k] for i in range(len(seq)-k+1)]
     
    def build_kmer_pos(seq, k):
        kmer_pos = defaultdict(list)
        for i, kmer in enumerate(extract_kmers(seq, k)):
            kmer_pos[kmer].append(i)
        return kmer_pos
     
    def mod_dotplot(ref_seq, query_seq,output, k):
        ref_kmers = build_kmer_pos(ref_seq, k)
        query_kmers = build_kmer_pos(query_seq, k)
        
        x, y = [], []
        for kmer in ref_kmers:
            if kmer in query_kmers:
                for ref_pos in ref_kmers[kmer]:
                    for query_pos in query_kmers[kmer]:
                        x.append(ref_pos)
                        y.append(query_pos)
        if 1==2:   ## Popup dot plot
            plt.scatter(x, y, s=1, alpha=0.5)
            plt.xlabel("Reference Genome Position")
            plt.ylabel("Query Genome Position")
            plt.show()
        # New code starts - Generate BED file
        # Collect all matched reference sequence positions (convert to BED format: 0-based start, 1-based end)
        bed_entries = []
        for kmer in ref_kmers:
            if kmer in query_kmers:
                for ref_pos in ref_kmers[kmer]:
                    start = ref_pos
                    end = ref_pos + k  # k-mer length
                    bed_entries.append((start, end))
        
        # Remove duplicates and sort
        bed_entries = sorted(set(bed_entries))
        
        # Write to BED file (assuming reference sequence name is "ref")
        with open(output, "w") as f:
            for ref_pos, query_pos in zip(x, y):
                f.write(f"{ref_pos}\t{query_pos}\t{k}\n")  # Two columns: ref_pos and query_pos


    # Example usage
    #ref = "ATGCGATCGATCGATCGATCG"
    #query = "CGATCGATCGATCGATCG"
    #mod_dotplot(ref, query, k=6)
    ref = seq1
    query = seq2
    mod_dotplot(ref, query,f"{output_dir}/kmer_matches_k21.bed", k=21)
    mod_dotplot(ref, query,f"{output_dir}/kmer_matches_k41.bed", k=41)

    with open(f'{output_dir}/kmer_matches.head' ,'w' ) as f:
        f.write('target1_pos\ttarget2_pos\tkmer_len\n')
    subprocess.run([f'cat {output_dir}/kmer_matches.head {output_dir}/kmer_matches_k21.bed {output_dir}/kmer_matches_k41.bed >{output_dir}/kmer_matches_kkk.bed'], shell=True)  #>n
    
    # New code ends
    R_txt=r"""
library(ggplot2)
library(dplyr)

input_file=read.table('kmer_matches_kkk.bed', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
k1=input_file %>% filter(kmer_len == 21)
k2=input_file %>% filter(kmer_len == 41)


p=ggplot()
 
p=p+    geom_point(data = k1, aes(x = target1_pos, y = target2_pos),size = 0.3, alpha = 0.3, color = "#e0e0eb",shape=16) 
p=p+    geom_point(data = k2, aes(x = target1_pos, y = target2_pos),size = 0.3, alpha = 0.3, color = "#006699",shape=16) 



p=p+theme_classic() +
  coord_equal() 

# Save as PDF
pdf(file = paste0('kmer_matches.kkk', ".pdf"), width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()        
"""
    with open(f'{output_dir}/kmer_matches.R','w') as f:
        f.write(R_txt)
    os.chdir(f'{output_dir}/')
    if 1==2:    #       pic!='no':
        subprocess.run([f'Rscript kmer_matches.R  '], shell=True)  #>n
    os.chdir(f'{current_dir}')    
if 'step2'=='step2':
    #print('step1__Use kmer to find all kmer41 points')
    #print('step2__Extend and filter the points')
    #print('step3__Sort by length, then calculate the nearest points for extension sequentially')
    #print('step4__Extract sequences separately for mafft alignment, calculate divergence time')
    #print('step5__Collect results and plot')
    if 'step2.1'=='step2.1' and 1==2:
        print('Use kmer to find all kmer points')
        from collections import defaultdict
        import matplotlib.pyplot as plt
         
        def extract_kmers(seq, k):
            return [seq[i:i+k] for i in range(len(seq)-k+1)]
         
        def build_kmer_pos(seq, k):
            kmer_pos = defaultdict(list)
            for i, kmer in enumerate(extract_kmers(seq, k)):
                kmer_pos[kmer].append(i)
            return kmer_pos
         
        def mod_dotplot(ref_seq, query_seq, k=21):
            ref_kmers = build_kmer_pos(ref_seq, k)
            query_kmers = build_kmer_pos(query_seq, k)
            
            x, y = [], []
            for kmer in ref_kmers:
                if kmer in query_kmers:
                    for ref_pos in ref_kmers[kmer]:
                        for query_pos in query_kmers[kmer]:
                            x.append(ref_pos)
                            y.append(query_pos)
            if 1==2:   ## Popup dot plot
                plt.scatter(x, y, s=1, alpha=0.5)
                plt.xlabel("Reference Genome Position")
                plt.ylabel("Query Genome Position")
                plt.show()
            # New code starts - Generate BED file
            # Collect all matched reference sequence positions (convert to BED format: 0-based start, 1-based end)
            bed_entries = []
            for kmer in ref_kmers:
                if kmer in query_kmers:
                    for ref_pos in ref_kmers[kmer]:
                        start = ref_pos
                        end = ref_pos + k  # k-mer length
                        bed_entries.append((start, end))
            
            # Remove duplicates and sort
            bed_entries = sorted(set(bed_entries))
            
            # Write to BED file (assuming reference sequence name is "ref")
            with open(f"{output_dir}/kmer_matches.bed", "w") as f:
                f.write('target1_pos\ttarget2_pos\n')
                for ref_pos, query_pos in zip(x, y):
                    f.write(f"{ref_pos}\t{query_pos}\n")  # Two columns: ref_pos and query_pos
    
            # New code ends
            R_txt=r"""
library(ggplot2)
library(dplyr)

input_file=read.table('kmer_matches.bed', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
p=ggplot(input_file, aes(x = target1_pos, y = target2_pos)) +
    geom_point(size = 0.3, alpha = 0.3, color = "#006699",shape=16) +theme_classic()

# Save as PDF
pdf(file = paste0('kmer_matches', ".pdf"), width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()        
"""
            with open(f'{output_dir}/kmer_matches.R','w') as f:
                f.write(R_txt)
            os.chdir(f'{output_dir}/')
            if 1==1:
                subprocess.run([f'Rscript kmer_matches.R  '], shell=True)  #>n
            os.chdir(f'{current_dir}')   
        # Example usage
        #ref = "ATGCGATCGATCGATCGATCG"
        #query = "CGATCGATCGATCGATCG"
        #mod_dotplot(ref, query, k=6)
        ref = seq1
        query = seq2
        mod_dotplot(ref, query, k=kmer_len)
        
    if 'step2.1'=='step2.1' :   
        print('Use kmer to find all kmer points')
        import os
        import subprocess
        from collections import defaultdict
        import numpy as np
        from multiprocessing import Pool, cpu_count
        
        def extract_kmers(seq, k):
            """Extract all k-mers"""
            return [seq[i:i + k] for i in range(len(seq) - k + 1)]
        
        def build_kmer_pos(seq, k):
            """Build k-mer position dictionary"""
            kmer_pos = defaultdict(list)
            for i, kmer in enumerate(extract_kmers(seq, k)):
                kmer_pos[kmer].append(i)
            return kmer_pos
        
        def process_kmer_chunk(args):
            """Process a chunk of k-mers for matching (for multiprocessing)"""
            ref_kmers, query_kmers, kmer_chunk = args
            x, y = [], []
            for kmer in kmer_chunk:
                if kmer in query_kmers:
                    ref_positions = ref_kmers[kmer]
                    query_positions = query_kmers[kmer]
                    # Use NumPy broadcasting to accelerate coordinate pair generation
                    x.extend(ref_positions * len(query_positions))
                    y.extend(np.repeat(query_positions, len(ref_positions)))
            return x, y
        
        def mod_dotplot_parallel(ref_seq, query_seq, k=21, num_processes=None):
            """Parallelized k-mer matching (optimized version)"""
            # Build k-mer position dictionaries
            ref_kmers = build_kmer_pos(ref_seq, k)
            query_kmers = build_kmer_pos(query_seq, k)
        
            # Find common k-mers
            common_kmers = list(set(ref_kmers.keys()) & set(query_kmers.keys()))
            if not common_kmers:
                print("No common k-mers found!")
                return
        
            # Chunk processing (to avoid a single process handling too many k-mers)
            num_processes = num_processes or cpu_count()  # Default: use all CPU cores
            chunk_size = max(1, len(common_kmers) // num_processes)
            kmer_chunks = [
                common_kmers[i:i + chunk_size]
                for i in range(0, len(common_kmers), chunk_size)
            ]
        
            # Process each k-mer chunk in parallel
            with Pool(num_processes) as pool:
                results = pool.map(
                    process_kmer_chunk,
                    [(ref_kmers, query_kmers, chunk) for chunk in kmer_chunks]
                )
        
            # Merge results from all processes
            x, y = [], []
            for chunk_x, chunk_y in results:
                x.extend(chunk_x)
                y.extend(chunk_y)
        
            # Generate BED file (only save matched coordinate pairs)
            with open(f"{output_dir}/kmer_matches.bed", "w") as f:
                f.write("ref_pos\tquery_pos\n")
                for ref_pos, query_pos in zip(x, y):
                    f.write(f"{ref_pos}\t{query_pos}\n")
        
            # R script generation and execution (unchanged)
            R_txt = r"""
        library(ggplot2)
        library(dplyr)
        
        input_file <- read.table('kmer_matches.bed', header=TRUE, sep='\t')
        p <- ggplot(input_file, aes(x=ref_pos, y=query_pos)) +
            geom_point(size=0.3, alpha=0.3, color="#006699", shape=16) + theme_classic()
        
        pdf(file='kmer_matches.pdf', width=10/2.54, height=10/2.54)
        print(p)
        dev.off()
        """
            with open(f"{output_dir}/kmer_matches.R", "w") as f:
                f.write(R_txt)
        
            # Execute R script
            os.chdir(output_dir)
            subprocess.run(["Rscript", "kmer_matches.R"], check=True)
            os.chdir(current_dir)
        
        # Example usage (seq1, seq2, kmer_len, output_dir, current_dir need to be defined in advance)
        ref = seq1
        query = seq2
        mod_dotplot_parallel(ref, query, k=kmer_len, num_processes=20)  # Use 4 processes
    if 'step2.2'=='step2.2':

        print('Manually analyze continuous points, mimicking moddotplot custom script slash')
        # Point pool 
        dict3_x_y_z={}  ## This is a numeric dictionary for xy
        with open(f"{output_dir}/kmer_matches.bed",'r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                pos1,pos2=eachline_arr
                pos1=int(pos1)
                pos2=int(pos2)
                if pos1 not in dict3_x_y_z: dict3_x_y_z[pos1]={}
                dict3_x_y_z[pos1][pos2]=1
        x_list=sorted(dict3_x_y_z.keys())    
        ####
        #
        
        x_gap_max=100         #x_ydot_neighbor also relates to this; after x moves right by +1, y can only have one point within ±5,      #Horizontal search
        y_gap_max=100                                                                                   #Vertical search
        extend_tolerance=10000
        ### Calculate diagonally upward
        with open(f"{output_dir}/kmer_matches.up_lines",'w') as fffff:
            fffff.write(f"start_x\tstart_y\tend_x\tend_y\tseq1_len\tseq2_len\tdot_num\tdotskip_num\tdot_used\tdots_string\n")  
            tmp_dot_list=[]
            used_region=[]
            for x in x_list:
                y_list=sorted(dict3_x_y_z[x].keys(), reverse=True)     ####y_list=list(dict3_x_y_z[x].keys);  y_list.sort()   
                #print(y_list)
                for y in y_list:
                    #Current point is (x,y)
                    x_start=x
                    y_start=y
                    end_mark=''
                    tolerance=0
                    serial=0;
                    while end_mark=='' and tolerance<extend_tolerance:
                        ##Check if the current starting point has been used in previous loops
                        #print(x_start,y_start)
                        if x_start not in dict3_x_y_z:              end_mark='yes';break
                        if y_start not in dict3_x_y_z[x_start]:    end_mark='yes';break 
                        #xy_distance=y_start-x_start
                        ##
                        if tmp_dot_list==[]:tmp_dot_list.append([x,y]);serial+=1
                        #Collect the next point
                        i1=0;good_mark=''
                        while i1<x_gap_max:
                            i1+=1
                            
                            if good_mark=='yes':break
                            #if serial<3:i2=1   #Slope of the first three points cannot be horizontal
                            #elif xy_distance<20:    i2=1        #The smallest duplication might be 100k, i.e., 20 points
                            #elif xy_distance<40:    i2=0
                            #elif xy_distance<90:    i2=-1
                            #else :                  i2=-2
                            i2=-2  #Strictly 45-degree upward slope, no downward points allowed in between
                            while i2<y_gap_max:
                                x_tmp=x_start+i1
                                y_tmp=y_start+i2  ###Move x one step right, y first stays (+0) then moves upward (++1)
                                i2+=1
                                if x_tmp in dict3_x_y_z:
                                    if y_tmp in dict3_x_y_z[x_tmp]:
                                        #New coordinates are (x_tmp,y_tmp)
                                        good_mark='yes';break
                        ##If the next point is collected, continue with it as the starting point            
                        if good_mark=='yes': 
                            tmp_dot_list.append([x_tmp,y_tmp])
                            x_start=x_tmp
                            y_start=y_tmp
                            serial+=1
                            tolerance=0  ##Reset tolerance
                        ##If no next point is collected, end the current line segment    
                        else:
                            tolerance+=1
                            x_start+=1
                            y_start+=1  ##
                            
                    #sys.exit()
                    ###Determine the number of points in the currently collected point sequence, slope, and goodness of fit. [[x1,y1],[x2,y2],...], calculate the slope and R² of the fitted line for all points
                    tmp_dot_list_len=len(tmp_dot_list)
                    if 1==1:#tmp_dot_list_len>30:
                        #print(tmp_dot_list)
                        #sys.exit()
                        all_start_x=tmp_dot_list[0][0];                 seq1_start=all_start_x     #seq1_start
                        all_start_y=tmp_dot_list[0][1];                 seq2_start=all_start_y     #seq2_start     #+(unit_len-1) is because the simplified start position was used
                        all_end_x=tmp_dot_list[-1][0];                  seq1_end=all_end_x      #seq1_end
                        all_end_y=tmp_dot_list[-1][1];                  seq2_end=all_end_y      #seq2_end
                        seq1_len=abs(seq1_end-seq1_start)+1
                        seq2_len=abs(seq2_end-seq2_start)+1
                        dots_string = '|'.join(f'{x},{y}' for x, y in tmp_dot_list)
                        all_x_num=tmp_dot_list[-1][0]-tmp_dot_list[0][0]+1
                        good_ratio=tmp_dot_list_len/all_x_num
                        
                        
                        used_mark=''
                        for one_region in used_region:
                            if all_start_x>=one_region[0] and all_start_y>=one_region[1] and all_end_x<=one_region[2] and all_end_y<=one_region[3]: used_mark='yes'
                        if used_mark=='':    
                            #if all_start_y>140000:print(all_start_x,all_start_y);#sys.exit()
                            fffff.write(f"{str(all_start_x)}\t{str(all_start_y)}\t{str(all_end_x)}\t{str(all_end_y)}\t{str(seq1_len)}\t{str(seq2_len)}\t{str(tmp_dot_list_len)}\t{str(all_x_num-tmp_dot_list_len)}\t{str(round(good_ratio,3))}\t{dots_string}\n")   #
                            used_region.append([all_start_x,all_start_y,all_end_x,all_end_y])
                            
                        ##Remove these points from the dictionary
                        for one_result in tmp_dot_list:
                            one_x,one_y=one_result
                            del dict3_x_y_z[one_x][one_y]
                    ##Initialize
                    tmp_dot_list=[]
                    #return False                
        R_txt=r"""
library(ggplot2)
library(dplyr)

input_file=read.table('kmer_matches.up_lines', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
dot_num_sum=sum(input_file$dot_num)
p=ggplot()+geom_segment(data=input_file,
    aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
    color = "blue",
    size = 1
  ) +theme_classic()+
  labs(
    title = paste("Total dot_num =", dot_num_sum),  # Dynamic title
    #x = "X-axis",
    #y = "Y-axis"
  ) 

# Save as PDF
pdf(file = paste0('kmer_matches.up_lines', ".pdf"), width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()        
"""
        with open(f'{output_dir}/kmer_matches.up_lines.R','w') as f:
            f.write(R_txt)
        os.chdir(f'{output_dir}/')
        if pic!='no':
            subprocess.run([f'Rscript kmer_matches.up_lines.R  '], shell=True)  #>n
        os.chdir(f'{current_dir}')              
    if 'step2.3'=='step2.3':
        print('Connect nearby line segments')
        lines=[]
        with open(f"{output_dir}/kmer_matches.up_lines",'r') as f:
            next(f)
            k=0
            for line in f:
                k+=1
                eachline_arr=line.strip().split('\t')
                start_x=    int(eachline_arr[0])	
                start_y=    int(eachline_arr[1])	
                end_x=      int(eachline_arr[2])	
                end_y=      int(eachline_arr[3])	
                seq1_len=    int(eachline_arr[4])		
                seq2_len=    int(eachline_arr[5])		
                dots_string=    eachline_arr[-1]
                dots=dots_string.split('|')
                dots_numeric = [[int(x), int(y)] for x, y in (point.split(',') for point in dots)]
                lines.append([k,start_x,start_y,end_x,end_y,seq1_len,seq2_len,dots_numeric])
        sorted_data = sorted(lines, key=lambda x: x[1], reverse=False)         
        lines=sorted_data
       
        def get_slople_R2(dot_list):
            points = np.array(dot_list)
            np_x = points[:, 0]    # Extract x and y values
            np_y = points[:, 1]
            slope, intercept = np.polyfit(np_x, np_y, 1)  # Use numpy's polyfit function for linear regression, 1 indicates first-order linear fit
            y_pred = slope * np_x + intercept  # Calculate R² goodness of fit
            ss_total = np.sum((np_y - np.mean(np_y))**2)
            ss_res = np.sum((np_y - y_pred)**2)
            if ss_total == 0:       r_squared = 0  # or other default value
            else:                   r_squared = 1 - (ss_res / ss_total)
            return slope,r_squared
        #Generate fitted line segments from list    
        def get_slope_R2_and_line(dot_list):
            points = np.array(dot_list)
            np_x = points[:, 0]
            np_y = points[:, 1]
            
            # Linear regression
            #slope, intercept = np.polyfit(np_x, np_y, 1)
            with np.errstate(invalid='ignore'):
                slope, intercept = np.polyfit(np_x, np_y, 1)

            # Calculate the start and end points of the fitted line segment
            x_start, x_end = np.min(np_x), np.max(np_x)
            y_start = slope * x_start + intercept
            y_end = slope * x_end + intercept
            
            # Calculate R²
            y_pred = slope * np_x + intercept
            ss_total = np.sum((np_y - np.mean(np_y)) ** 2)
            ss_res = np.sum((np_y - y_pred) ** 2)
            r_squared = 1 - (ss_res / ss_total) if ss_total != 0 else 0
            
            return x_start, y_start,x_end, y_end,slope,intercept,r_squared         # End point of the fitted line segment
            
        lines_nun=k
        used_k=set()
        iii=0
        with open(f"{output_dir}/kmer_matches.up_lines2",'w') as f:
            f.write(f"index\tstart_x\tstart_y\tend_x\tend_y\tdot_num\tslope\tr_squared\n")
            for one_line in lines:
                k,start_x,start_y,end_x,end_y,seq1_len,seq2_len,dots_numeric=one_line
                if k in used_k :continue
                good_dots=dots_numeric
                used_k.add(k)
                ###
                for one_line in lines:
                    one_k,one_start_x,one_start_y,one_end_x,one_end_y,one_seq1_len,one_seq2_len,one_dots_numeric=one_line
                    one_start_x,one_start_y,one_end_x,one_end_y=[one_dots_numeric[0][0],one_dots_numeric[0][1],one_dots_numeric[-1][0],one_dots_numeric[-1][1]]
                    ##
                    good_start_x,good_start_y,good_end_x,good_end_y=[good_dots[0][0],good_dots[0][1],good_dots[-1][0],good_dots[-1][1]]
                    ##
                    if one_k in used_k :continue
                    ##Filter out intervals already included
                    if good_start_x<=one_start_x and good_start_y<=one_start_y and one_end_x<=good_end_x and one_end_y<=good_end_y: used_k.add(one_k);continue
                    ##Filter out points that are too far away
                    if abs(one_start_x-good_end_x)>2000 or abs(one_start_y-good_end_y)>2000:continue
                    ###Calculate slope and intercept — new
                    tmp_newdots_list=good_dots+one_dots_numeric
                    #print(tmp_newdots_list)
                    one_slope,one_r_squared=get_slople_R2(tmp_newdots_list)
                    #print(one_slope,one_r_squared)
                    #print(one_line[:5])
                    #print(tmp_newdots_list)
                    ###
                    if one_slope<0.9 or one_slope>1.1:continue
                    if one_r_squared<0.7:continue
                    ########
                    good_dots=tmp_newdots_list
                    used_k.add(one_k)
                #######    
                #print(good_dots)    
                dot_num=len(good_dots)
                #dots_string = '|'.join(f'{x},{y}' for x, y in good_dots)
                x_start, y_start,x_end, y_end,slope,intercept,r_squared=get_slope_R2_and_line(good_dots)
                if x_end-x_start<5000:continue
                iii+=1
                f.write(f"{iii}\t{x_start}\t{int(round(y_start,0))}\t{x_end}\t{int(round(y_end,0))}\t{dot_num}\t{slope}\t{r_squared}\n")
                #if x_start==2430:break
        R_txt=r"""
library(ggplot2)
library(dplyr)

input_file=read.table('kmer_matches.up_lines2', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
dot_num_sum=sum(input_file$dot_num)
p=ggplot()+geom_segment(data=input_file,
    aes(x = start_x, y = start_y, xend = end_x, yend = end_y),
    color = "blue",
    size = 1
  ) +theme_classic()+
  labs(
    title = paste("Total dot_num =", dot_num_sum),  # Dynamic title
    #x = "X-axis",
    #y = "Y-axis"
  ) 

# Save as PDF
pdf(file = paste0('kmer_matches.up_lines2', ".pdf"), width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()        
"""
        with open(f'{output_dir}/kmer_matches.up_lines.R','w') as f:
            f.write(R_txt)
        os.chdir(f'{output_dir}/')
        if pic!='no':
            subprocess.run([f'Rscript kmer_matches.up_lines.R  '], shell=True)  #>n
        os.chdir(f'{current_dir}')  
    if 'step2.31'=='step2.31':       
        ##For segment mafft later, if there are 300k+ consecutive segments, they cannot be aligned due to memory issues. I manually split them into roughly 100k bp small segments
        def split_line(start_x, start_y, end_x, end_y, n):
            segments = []
            dx = (end_x - start_x) / n
            dy = (end_y - start_y) / n
            
            for i in range(n):
                seg_start_x = start_x + i * dx
                seg_start_y = start_y + i * dy
                seg_end_x = start_x + (i + 1) * dx
                seg_end_y = start_y + (i + 1) * dy
                segments.append((seg_start_x, seg_start_y, seg_end_x, seg_end_y))
            
            return segments
        ###
        with open(f'{output_dir}/kmer_matches.up_lines2_revise','w') as f2:
            f2.write("index\tstart_x\tstart_y\tend_x\tend_y\n")
            with open(f'{output_dir}/kmer_matches.up_lines2','r') as f:
                next(f)
                ii=1
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    if len(eachline_arr)!=8:continue
                    index,start_x,start_y,end_x,end_y,dot_num,slope,r_squared=eachline_arr
                    len1=int(end_x)-int(start_x)
                    #print(eachline_arr)
                    if len1>100000:  
                        part_num=math.floor(len1/100000)+1
                        print("Detected an overly long alignment, splitting into part_num segments to prevent mafft memory overflow",part_num)
                        segments = split_line(int(start_x), int(start_y), int(end_x), int(end_y), part_num)
                        for i, (sx, sy, ex, ey) in enumerate(segments):
                            ii+=1
                            #print(f"Segment {i+1}: start=({sx:.2f}, {sy:.2f}), end=({ex:.2f}, {ey:.2f})")
                            f2.write(f"{ii}\t{sx:.0f}\t{sy:.0f}\t{ex:.0f}\t{ey:.0f}\n")
                    else:
                        ii+=1
                        f2.write(f"{ii}\t{start_x}\t{start_y}\t{end_x}\t{end_y}\n")
                
    if 'step2.4'=='step2.4':        
        subprocess.run([f'mkdir {output_dir}/segments'], shell=True)  #>n
        index_list=[]
        with open(f'{output_dir}/kmer_matches.up_lines2_revise','r') as f:
            next(f)
            for line in f:
                eachline_arr=line.strip().split('\t')
                if len(eachline_arr)!=5:continue
                #index,start_x,start_y,end_x,end_y,dot_num,slope,r_squared=eachline_arr
                index,start_x,start_y,end_x,end_y=eachline_arr
                
                seq1_part=      seq1[max(0,int(start_x)-100):int(end_x)+100]
                seq2_part=      seq2[max(0,int(start_y)-100):int(end_y)+100]
                with open (f"{output_dir}/segments/{index}.fa",'w') as f2:
                    f2.write(f'>seq1_part\n{seq1_part}\n>seq2_part\n{seq2_part}')
                index_list.append(index)    
                
        def run_step(one_index):   # --op 4 
            os.system(f"mafft --thread 1 --genafpair   --quiet  {output_dir}/segments/{one_index}.fa > {output_dir}/segments/{one_index}.mafft")
            seq1_mafft=''
            seq2_mafft=''
            with open(f"{output_dir}/segments/{one_index}.mafft",'r') as f:
                for line in f:
                    if line.strip()=='>seq1_part':mark=1
                    elif line.strip()=='>seq2_part':mark=2            
                    else:
                        if mark==1: seq1_mafft+=line.strip()    
                        elif mark==2: seq2_mafft+=line.strip()   
            with open(f'{output_dir}/segments/{one_index}.mafft','w') as f:
                f.write(f'>seq1_part\n{seq1_mafft}\n>seq2_part\n{seq2_mafft}')    
            #Separate by indels, keep continuous segments longer than 100bp
            mafft_len=len(seq1_mafft)
            kkk=0; mid_part='';result_list=[]
            with open(f'{output_dir}/segments/{one_index}.mafft_part','w') as f:
                for z1,z2 in zip(seq1_mafft,seq2_mafft):
                    kkk+=1
                    if z1=='-'or z2=='-': 
                        if len(mid_part)>10:
                            result_list.append([kkk-len(mid_part),kkk-1,mid_part])
                            f.write(f"{kkk-len(mid_part)}\t{kkk-1}\t{mid_part}\n")
                            mid_part=''
                    else:
                        if z1==z2:zz='Y'
                        else:zz='.'
                        mid_part+=zz
                if len(mid_part)>10:
                    result_list.append([kkk-len(mid_part),kkk-1,mid_part])
                    f.write(f"{kkk-len(mid_part)}\t{kkk-1}\t{mid_part}\n")    
            ####Each sequence needs to start and end with five consecutive matches to prevent misalignment caused by INDELs
            result_list2=[]
            with open(f'{output_dir}/segments/{one_index}.mafft_part2','w') as f:
                for one_result in result_list:
                    start, end, one_seq = one_result
                    one_seq_len = len(one_seq)
                    if 'YYYYY' not in one_seq:continue
                    # Find the starting 'YYYYY' (forward direction)
                    trim_start = 0
                    for kkk in range(one_seq_len - 4):  # Need at least 5 characters
                        if one_seq[kkk:kkk+5] == 'YYYYY':
                            trim_start = kkk
                            break
                    
                    # Find the ending 'YYYYY' (reverse direction)
                    trim_end = 0
                    for kkk in range(one_seq_len - 4, -1, -1):  # Traverse from back to front
                        if one_seq[kkk:kkk+5] == 'YYYYY':
                            trim_end = one_seq_len - (kkk + 5)  # Calculate offset from the end to the end of 'YYYYY'
                            break
                    
                    # Calculate new coordinates and sequence
                    new_start = start + trim_start
                    new_end = end - trim_end
                    new_seq = one_seq[trim_start : one_seq_len - trim_end]  # Slice directly to the end
                    if len(new_seq)>10:
                        f.write(f"{new_start}\t{new_end}\t{new_seq}\n")    
                        result_list2.append([new_start,new_end,new_seq])
            ####If a 10bp segment has more than 5 SNPs, it is considered abnormal and the segment is replaced with # 
            result_list3=[]
            with open(f'{output_dir}/segments/{one_index}.mafft_part3','w') as f:
                for one_result in result_list2:
                    start, end, one_seq = one_result
                    one_seq_len=len(one_seq)  
                    #print(one_seq)
                    kkk=0
                    error_pos_set=set()
                    while kkk<one_seq_len-10:
                        seq_base10=one_seq[kkk:kkk+10]
                        if seq_base10.count('.')>5:
                            jjj=kkk
                            while jjj<kkk+10:
                                error_pos_set.add(jjj)
                                jjj+=1
                        kkk+=1
                    new_seq=one_seq    
                    for one_j in error_pos_set:
                        new_seq = new_seq[:one_j] + '#' + new_seq[one_j+1:]
                    #print(new_seq)    
                    f.write(f"{start}\t{end}\t{new_seq}\n") 
                    result_list3.append([start,end,new_seq])    
            ####Three or more consecutive mismatches within each sequence are considered abnormal and are consolidated into a single mismatch to avoid double-counting. Similar to indels, 4 dots = 4 *, 100 dots = 100 *
            result_list4=[]
            with open(f'{output_dir}/segments/{one_index}.mafft_part4','w') as f:
                for one_result in result_list3:
                    start, end, one_seq = one_result
                    one_seq = re.sub(r"\.{3,}", lambda m: '*' * len(m.group()), one_seq)
                    f.write(f"{start}\t{end}\t{one_seq}\n") 
                    result_list4.append([start,end,one_seq])
            ####Calculate final results
            with open(f'{output_dir}/segments/{one_index}.mafft_part5','w') as f:
                align_length=0
                var_num=0
                for one_result in result_list4:
                    start, end, one_seq = one_result
                    # Replace consecutive '*' with a single '*'
                    one_seq2 = re.sub(r'\*+', '*', one_seq)
                    # Replace consecutive '#' with a single '#'
                    one_seq2 = re.sub(r'#+', '#', one_seq2)
                    #
                    var_num+=one_seq2.count('#')+one_seq2.count('.')
                    align_length+=len(one_seq2)
                #############
                mutation_rate = 2.5e-9     # Mutation rate: 2.5 × 10⁻⁹ per site per generation——————————————Modified on 20250903
                if align_length!=0:
                    # Calculate mutation frequency
                    p = var_num / align_length
                    # Calculate divergence time (generations)
                    T_generations = p / (2 * mutation_rate)
                        
                        #
                    f.write(f"align_length\tvar_num\tpercent\tTime\n") 
                    f.write(f"{align_length}\t{var_num}\t{p:.4f}\t{T_generations:.0f}\n")
                    
            ####    
        with Pool(processes=20) as pool:
            # Use imap to get results one by one
            for i, result in enumerate(pool.imap(run_step, index_list), start=1):
                # Results can be processed here, e.g., stored or printed
                progress = (i / len(index_list)) * 100
                sys.stdout.write(f"\rProgress: {progress:.2f}%")
                sys.stdout.flush()    
                
        with open(f'{output_dir}/kmer_matches.up_lines3','w') as f3: 
            f3.write(f"index\tstart_x\tstart_y\tend_x\tend_y\talign_length\tvar_num\tpercent\tTime\n")
            with open(f'{output_dir}/kmer_matches.up_lines2_revise','r') as f2: 
                next(f2)
                for line in f2:
                    eachline=line.strip()
                    eachline_arr=eachline.split('\t')
                    index,start_x,start_y,end_x,end_y=eachline_arr
                    if int(start_x)<0:start_x=0
                    if int(start_y)<0:start_y=0
                    with open(f"{output_dir}/segments/{index}.mafft_part5",'r') as f:
                        next(f)
                        eachlinee= f.read().strip()
                    newline=  f"{index}\t{start_x}\t{start_y}\t{end_x}\t{end_y}\t{eachlinee}\n"
                    f3.write(newline)
    if 'step2.5'=='step2.5':    
        R_txt=r"""
library(ggplot2)
library(dplyr)

input_file=read.table('kmer_matches.up_lines3', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
dot_num_sum=sum(input_file$dot_num)

filtered_input <- input_file %>% 
mutate(
  category = case_when(
    Time<500000 ~ "<50",
    Time<1000000  ~ "50-100",
    Time<5000000  ~ "100-500",
    Time<10000000  ~ "500-1000",
    Time<20000000  ~ "1000-2000",
    Time>=20000000  ~ "2000+",
    TRUE ~ "other"
  )
) %>%filter(category!="other")

# Define color values
color_values <- c(
    "<50" = "red",
    "50-100" = "orange",
    "100-500" = "yellow",
    "500-1000" = "green",
    "1000-2000" = "blue",
    "2000+" = "purple"
)    
    
p=ggplot()+geom_segment(data=filtered_input,
    aes(x = start_x, y = start_y, xend = end_x, yend = end_y,color=category),
    size = 1
  ) +theme_classic()+
  labs(
    title = paste("Total dot_num =", dot_num_sum),  # Dynamic title
    #x = "X-axis",
    #y = "Y-axis"
  ) +
    scale_color_manual(values = color_values, drop = FALSE)+
            theme(
        #axis.ticks.y = element_blank(),
        #axis.text.y = element_blank(),
        #legend.position = "none",
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank()
        #panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Manually add border
    )

# Save as PDF
pdf(file = paste0('kmer_matches.up_lines3(sum_by_time)', ".pdf"), width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()        
"""
        with open(f'{output_dir}/kmer_matches.up_lines.R','w') as f:
            f.write(R_txt)
        os.chdir(f'{output_dir}/')
        if 1==1:
            subprocess.run([f'Rscript kmer_matches.up_lines.R  '], shell=True)  #>n
        os.chdir(f'{current_dir}')           
if 'step3'=='step3':
    # New code ends
    R_txt=r"""
library(ggplot2)
library(dplyr)
 # Get command line arguments
args <- commandArgs(TRUE)
# Input file paths
seq1_len <- args[1]
seq2_len <- args[2]
seq1_len <- as.numeric(seq1_len)
seq2_len <- as.numeric(seq2_len)
print(seq1_len)
print(seq2_len)

input_file=read.table('kmer_matches_kkk.bed', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
k1=input_file %>% filter(kmer_len == 21)
k2=input_file %>% filter(kmer_len == 41)
max_val <- max(input_file$target1_pos, input_file$target2_pos)
input_file2=read.table('kmer_matches.up_lines3', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')


filtered_input2 <- input_file2 %>% 
mutate(
  category = case_when(
    percent<0.01 ~ "<1",         #2 million years
    percent<=0.02  ~ "1-2",       #4 million years
    percent<=0.03  ~ "2-3",      #6 million years
    percent<=0.04  ~ "3-4",     #8 million years
    percent<=0.08  ~ "4-8",    #16 million years
    percent>=0.08  ~ "8+",   #16 million years+
    TRUE ~ "other"
  )
) %>%filter(category!="other")%>%filter(align_length>3000)

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
 
p=p+    geom_point(data = k1, aes(x = target1_pos, y = target2_pos),size = 0.1,  color = "#e0e0eb",shape=16) 
#p=p+    geom_point(data = k2, aes(x = target1_pos, y = target2_pos),size = 0.1,  color = "#006699",shape=16) 

p=p+    geom_segment(data=filtered_input2,
    aes(x = start_x, y = start_y, xend = end_x, yend = end_y,color=category),
    linewidth =1.5
  ) 
p=p +coord_equal() +  
    #xlim(0, max_val) +  # x-axis starts from 0, range matches y-axis
    #ylim(0, max_val) +  # y-axis starts from 0, range matches x-axis
    xlim(0, 300000) +  
    ylim(0, 300000) + 
    annotate("segment", x = 0, y = 0, xend = seq1_len, yend = 0, color = "black", linewidth = 1) +
    annotate("segment", x = 0, y = 0, xend = 0, yend = seq2_len, color = "black", linewidth = 1)+
    scale_color_manual(values = color_values, drop = FALSE)+
    # theme(aspect.ratio = 1) +  # Optional: force the graph to be square
    theme_classic()+
    theme(
        #axis.ticks.y = element_blank(),
        #axis.text.y = element_blank(),
        legend.position = "none",
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank()
        #panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Manually add border
    )


# Save as PDF
pdf(file = paste0('sum_by_percent', ".pdf"), width = 10 / 2.54, height = 10 / 2.54)
print(p)
dev.off()        
"""
    with open(f'{output_dir}/sum.R','w') as f:
        f.write(R_txt)
    os.chdir(f'{output_dir}/')
    if pic!='no':
        subprocess.run([f'Rscript sum.R {raw_seq1_len} {raw_seq2_len}'], shell=True)  #>n
    os.chdir(f'{current_dir}')         
if 'step4'=='step4':
    print('Remove overlapping small segments')
    sum_info=[]
    with open (f'{output_dir}/kmer_matches.up_lines3','r') as f: 
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            index,start_x,start_y,end_x,end_y = eachline_arr[:5]	
            sum_info.append([index,int(start_x),int(start_y),int(end_x),int(end_y),int(end_x)-int(start_x)+1,int(end_y)-int(start_y)+1])
    with open (f'{output_dir}/kmer_matches.up_lines4','w') as f2: 
        f2.write("index\tstart_x\tstart_y\tend_x\tend_y\tdot_num\tslope\tr_squared\talign_length\tvar_num\tpercent\tTime\n")
        with open (f'{output_dir}/kmer_matches.up_lines3','r') as f: 
            next(f)
            for line in f:
                eachline=line.strip()
                eachline_arr=eachline.split('\t')           
                index,start_x,start_y,end_x,end_y = eachline_arr[:5]	
                end_x=int(end_x)
                start_x=int(start_x)
                end_y=int(end_y)
                start_y=int(start_y)
                len1=end_x-start_x+1
                len2=end_y-start_y+1
                bad_mark=''
                for one in sum_info:
                    one_index,one_start_x,one_start_y,one_end_x,one_end_y,one_len1,one_len2=one
                    if one_index==index:continue
                    if len1<one_len1:
                        if start_x>one_start_x and end_x<one_end_x:
                            bad_mark='yes';break
                        elif start_x<one_start_x and end_x>one_start_x and end_x<one_end_x and (end_x-one_start_x)/len1>0.5:  
                            bad_mark='yes';break
                        elif start_x>one_start_x and start_x<one_end_x and end_x>one_end_x and (one_end_x-start_x)/len1>0.5:   
                            bad_mark='yes';break
                    if len2<one_len2:
                        if start_y>one_start_y and end_y<one_end_y:
                            bad_mark='yes';break
                        elif start_y<one_start_y and end_y>one_start_y and end_y<one_end_y and (end_y-one_start_y)/len2>0.5:  
                            bad_mark='yes';break
                        elif start_y>one_start_y and start_y<one_end_y and end_y>one_end_y and (one_end_y-start_y)/len2>0.5:   
                            bad_mark='yes';break    
                if bad_mark=='':
                    f2.write(eachline+'\n')
                    
                    
if 'step6'=='step6':
    print('Provide similarity and reference time, match distance and number of SNPs')
    result_list=[]
    with open(f"{output_dir}/kmer_matches.up_lines3",'r') as f:
        next(f)
        for line in f:
            align_length,var_num,percent=line.strip().split('\t')[-4:-1]
            align_length=int(align_length)
            var_num=int(var_num)
            percent=float(percent)
            result_list.append([align_length,var_num,percent])
    print('Strategy: keep the longest alignment line, find the maximum and minimum percent. If any other line besides the longest alignment line has percent equal to the max or min, delete it. Sum up all sequences and recalculate percent')
    result_list_sorted = sorted(result_list, key=lambda x: x[0], reverse=True)
    # Extract the third column, find max and min values
    third_col = [row[2] for row in result_list_sorted]
    
    max_val = max(third_col)
    min_val = min(third_col)
    
    #print(f"Max value of third column: {max_val}")
    #print(f"Min value of third column: {min_val}")
    kk=0;align_all=0;var_all=0
    for one in result_list_sorted:
        kk+=1
        if kk==1:
            align_all,var_all,_=   one
        else:
            if one[2]==max_val or one[2]==min_val :continue
            align_all+=one[0]
            var_all+=one[1]
    newpercent= 1-      var_all /align_all
    with open(f"{output_dir}/ALL_SUM",'w') as f:
        f.write(f'{align_all}\t{var_all}\t{newpercent:.4f}\n')
        
        
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))      
time_end=timeit.default_timer()
print('All the running time: %.0f Seconds'%(time_end-time_start))