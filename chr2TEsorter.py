#!/usr/bin/python
# -*- coding: UTF-8 -*-
import sys
argvs=sys.argv
argv_len=len(sys.argv)
if not argvs[1].startswith("part"):
    if argv_len < 2 or argvs[1] in ("h", "-h", "help", "-help") :
        print ("chr2TEsorter.py-----help:")
        print ("")
        print ("Usage：")
        print ("chr2TEsoter.py   partX   【20000-40000】")
        print ("python chr2TEsorter.py part_all_LIRinterarray_auto V069.hap1 Chr9 15178574 15614639     ####")

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
##

window_size=10000
step_size=5000
interval=50000*'N'   ## For plotting
interval_len=len(interval)

time_start=timeit.default_timer()
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  
 
## Specify whether to run TEsorter
TEsorter_step='';TEsorter_step='yes'
if TEsorter_step=='':print('Note!!!!!!! TEsorter step omitted')

## Specify whether to run moddotplot
moddotplot_step='';moddotplot_step='yes'
if moddotplot_step=='':print('Note!!!!!!! moddotplot step omitted')

blast_cen_1305='yes'


auto_i=''
target_list=[]    
if  os.path.exists('./chr2TEsorter')==False: 
    subprocess.run(["mkdir chr2TEsorter"], shell=True)

#step0
if argv1=="part1":
    target_list.append(["WoollyGrape_hap1","Chr1:13960000-14100000"])    # 0.14
    target_list.append(["WoollyGrape_hap1","Chr1:14700000-15200000"])   #0.5
    
    target_list.append(["V008.hap1","Chr1:14700000-15150000"])   #0.45
    target_list.append(["V008.hap1","Chr1:16300000-16500000"])   #0.2
    #target_list.append(["V008.hap1","Chr1:18100000-18300000"])
    
    target_list.append(["V022.hap1","Chr9:12800000-13200000"])   #0.4
    
    target_list.append(["BlackMonukka_hap2","Chr9:17500000-18300000"])  #0.8
    
    
if argv1=="part2":
    target_list.append(["V059.hap2","Chr9:16200000-18000000"]) 
   
if argv1=="part3_Chr10,11,15,16,14":
    target_list.append(['V105.hap1','Chr9:17900000-18200000'])
    target_list.append(['V030.hap1','Chr9:17850000-18150000'])
    target_list.append(['V031.hap1','Chr19:17900000-18500000'])
    target_list.append(['PN40024','Chr10:19600000-20400000'])
    target_list.append(['PN40024','Chr11:13000000-13800000'])
    
    target_list.append(['PN40024','Chr15:8900000-9400000']) 
    target_list.append(['PN40024','Chr16:12700000-13100000'])
    target_list.append(['PN40024','Chr14:14950000-15850000'])
if argv1=="part3_Chr10_11_15_16_small":    
    target_list.append(['PN40024','Chr10:20040000-20190000'])
    target_list.append(['PN40024','Chr11:13170000-13320000'])
    target_list.append(['PN40024','Chr15:8950000-9100000']) 
    target_list.append(['PN40024','Chr16:12720000-12870000'])


    
if argv1=="part4_V106.hap1_sino_Chr9":    
    target_list.append(['V106.hap1','Chr9:18800000-19500000'])
    

if argv1=="part5_VSat6_LIR1305":
    target_list.append(['PN40024','Chr18:13931960-13941270'])          ## A complete forward tekay
    target_list.append(['WoollyGrape_hap1','Chr1:14923450-14932852'])   ## A complete reverse tekay, with snapgene map
    target_list.append(['PN40024','Chr10:20040000-20190000'])    ## VSat6 region on Chr10
    target_list.append(['PN40024','Chr1:15248175-15254605'])  ### A cen1305 block on Chr1
    
    
    target_list.append(['V106.hap1','Chr9:18800000-19500000'])

if argv1=="part6_Chr7_Chr14_Chr17":
    target_list.append(['PN40024','Chr7:13450000-13800000'])
    target_list.append(['PN40024','Chr14:15100000-15500000'])
    target_list.append(['PN40024','Chr17:15200000-15600000'])
    target_list.append(['Hongmunage_hap1','Chr7:13350000-13650000'])
    target_list.append(['Hongmunage_hap1','Chr17:15100000-15400000'])

if argv1=="part7_PN40024_allchromosome":    
    target_list.append(['PN40024','Chr1:15154005-15435737'])
    target_list.append(['PN40024','Chr13:12020240-12180415'])
    target_list.append(['PN40024','Chr7:13458101-13762988']) 
    target_list.append(['PN40024','Chr14:15144034-15480228'])
    target_list.append(['PN40024','Chr17:15221832-15560758'])
    target_list.append(['PN40024','Chr2:12939361-13300832'])
    target_list.append(['PN40024_hap1','Chr3:14353397-14553397'])
    target_list.append(['PN40024','Chr4:13291174-13483935'])
    target_list.append(['PN40024','Chr5:14008287-14357682'])
    target_list.append(['PN40024','Chr6:10606684-10966038'])
    target_list.append(['PN40024','Chr8:6679384-7071344'])
    target_list.append(['PN40024','Chr9:15113362-15449399'])
    target_list.append(['PN40024','Chr19:15995936-16195936'])
    target_list.append(['PN40024','Chr18:16390000-16620000'])
    
if argv1=="V008.hap1":    
    target_list.append(['V008.hap1','Chr1:14751815-15119418'])    
    target_list.append(['V008.hap1','Chr2:13387905-13716403'])  
    target_list.append(['V008.hap1','Chr3:14579344-14970242'])  
    target_list.append(['V008.hap1','Chr4:13947089-14159125'])  
    target_list.append(['V008.hap1','Chr5:13646164-14108350'])  
    target_list.append(['V008.hap1','Chr6:10783933-11258311'])  
    target_list.append(['V008.hap1','Chr7:14425979-14425964'])  
    target_list.append(['V008.hap1','Chr8:7974898-8359064'])      
    target_list.append(['V008.hap1','Chr9:16498098-16954837'])  
    target_list.append(['V008.hap1','Chr12:12149378-12627728'])  
if argv1=="V008.hap1_Chr1":    
    target_list.append(['V008.hap1','Chr1:14751815-15119418'])  
    
    #
if argv1=="V031.hap2_Chr18":    
    target_list.append(['V031.hap2','Chr18:16000000-17000000'])  
    
if argv1=="Athila_perChrMost_V117hap2_Chr6":        
    target_list.append(['V117.hap2','Chr6:9800000-16000000'])  
    
if argv1=="V112_hap_Chr16":        
    target_list.append(['V112.hap1','Chr16:8000000-13500000']) 
    
if argv1=="V105_hap1_Chr18":        
    target_list.append(['V105.hap1','Chr18:16650000-16950000'])     
if argv1=="V124_hap1_Chr18":        
    target_list.append(['V105.hap1','Chr18:16000000-16600000']) 
if argv1=="V124_hap2_Chr4_13498189-13570340":        
    target_list.append(['V124.hap2','Chr4:13498189-13570340']) 
if argv1=="V123_hap1_Chr9_18502955-18627537":        
    target_list.append(['V123.hap1','Chr9:18502955-18627537'])     

#if argv1=="PN40024_chr2":    
   # target_list.append(['PN40024','Chr2:10000000-15000000']) 
if argv1=="PN40024_chr1":    
    target_list.append(['PN40024','Chr1:11000000-16000000']) 
if argv1=="PN40024_chr5":    
    target_list.append(['PN40024','Chr5:11500000-16500000']) 
if argv1=="PN40024_chr11":    
    target_list.append(['PN40024','Chr11:11250000-16250000'])  
  
if argv1=="PN40024_chr13":    
    target_list.append(['PN40024','Chr13:9500000-14500000']) 
if argv1=="PN40024_chr14":    
    target_list.append(['PN40024','Chr14:13000000-18000000'])     

if argv1=="PN40024_chr15":    
    target_list.append(['PN40024','Chr15:5500000-10500000'])     
#if argv1=="PN40024_chr16":    
#    target_list.append(['PN40024','Chr16:7000000-15000000'])     
if argv1=="PN40024_chr17":    
    target_list.append(['PN40024','Chr17:13000000-18000000'])     
if argv1=="PN40024_chr18":    
    target_list.append(['PN40024','Chr18:14000000-19000000'])  
    
if argv1=="PN40024_chr1_core":    
    target_list.append(['PN40024','Chr1:15050000-15550000'])     
if argv1=="PN40024_chr5_core":    
    target_list.append(['PN40024','Chr5:13900000-14400000'])     
if argv1=="PN40024_chr11_core":    
    target_list.append(['PN40024','Chr11:13000000-13500000'])
if argv1=="PN40024_chr17_core":    
    target_list.append(['PN40024','Chr17:15150000-15650000'])      
if argv1=="PN40024_chr18_core":    
    target_list.append(['PN40024','Chr18:16250000-16750000']) 
    
if argv1=="PN40024_chr17_core_part":    
    target_list.append(['PN40024','Chr17:15400000-15500000'])      
if argv1=="PN40024_chr17_core_part2":    
    target_list.append(['PN40024','Chr17:15400000-15500000'])     
if argv1=="PN40024_chr17_core_part3":    
    target_list.append(['PN40024','Chr17:15440000-15450000'])       
    
###Finding the most complete TEkay LTR, I see that V012.hap1 on chr5 has 3 consecutive LTRs, both forward and reverse, all 100%
if argv1=="NewTekayComplete":    
    target_list.append(['V012.hap1','Chr5:14000000-15000000']) 
    
####Exploring whether LIRs on Chr1, 5, 17, 18 are related to VSat6 on Chr10/11/15/16
if argv1=="compare___VSat6__LIR": 
    #target_list.append(['PN40024','Chr1:15050000-15550000']) 
    target_list.append(['PN40024','Chr5:13900000-14400000'])  
    target_list.append(['PN40024','Chr17:15150000-15650000'])      
    target_list.append(['PN40024','Chr18:16250000-16750000']) 
    #
    target_list.append(['PN40024','Chr10:20040000-20190000'])
    target_list.append(['PN40024','Chr11:13170000-13320000'])
    target_list.append(['PN40024','Chr15:8950000-9100000']) 
    target_list.append(['PN40024','Chr16:12720000-12870000'])   
    
    target_list.append(['V105.hap1','Chr9:17890000-18040000'])
    target_list.append(['V030.hap1','Chr9:17850000-18000000'])
    target_list.append(['V031.hap1','Chr19:18350000-18500000'])   
    
if argv1=="V120.hap1.Chr18_pianxie":    
    target_list.append(['V120.hap1','Chr18:16400000-16800000'])    
    
if  argv1=="double_LIR":    
    target_list.append(['WoollyGrape_hap1','Chr9:17500000-18200000']) 
if  argv1=="double_LIR3":    
    target_list.append(['WoollyGrape_hap1','Chr9:17643900-17922344'])     #18151280
    
if  argv1=="double_LIR2":    
    target_list.append(['WoollyGrape_hap1','Chr9:17700000-17900000'])        

if argv1=="bad_LIR1":  #（33b）V081.hap2 and（46a）V120.hap1
    target_list.append(['V081.hap2','Chr9:16600000-16850000']) 
if argv1=="bad_LIR2":    
    target_list.append(['V120.hap1','Chr9:18200000-18450000'])        

if argv1=="noLIR_haveVSat6":    
    target_list.append(['V048.hap1','Chr3:15300000-15500000'])  
    target_list.append(['V048.hap1','Chr5:13200000-13700000'])  

if argv1=="VSat6_LIRSat1":    
    target_list.append(['V098.hap2','Chr2:13850000-14160000'])    
    target_list.append(['PN40024','Chr10:20040000-20190000'])
    target_list.append(['PN40024','Chr11:13170000-13320000'])
if argv1=="VSat6_LIRSat1_2": 
    target_list.append(['V098.hap2','Chr2:13950000-14050000'])    
    target_list.append(['PN40024','Chr10:20040000-20190000'])

##Looking at similarities between Tekay on different chromosomes
if argv1=="Chromosomes_LIRSat": 
    #target_list.append(['V124.hap2','Chr4:13515000-13550000'])    
    #target_list.append(['V123.hap1','Chr9:184512955-18632955'])
    target_list.append(['V124.hap2','Chr4:13428189-13640340'])   
    target_list.append(['V123.hap1','Chr9:18432955-18697537'])
    target_list.append(['V037.hap2','Chr7:14382973-14710457'])
    target_list.append(['V107.hap1','Chr14:16034600-16408612'])
    target_list.append(['V032.hap2','Chr1:15428505-15829144'])
    target_list.append(['V081.hap1','Chr18:16300000-16520000'])
    
if argv1=="Chromosomes_LIRSat1_15": 
    target_list.append(['V124.hap2','Chr4:13428189-13640340'])    
    target_list.append(['V032.hap2','Chr1:15428505-15829144'])
if argv1=="Chromosomes_LIRSat1_16":    
    target_list.append(['V124.hap2','Chr4:13428189-13640340'])   
    target_list.append(['V081.hap1','Chr18:16300000-16520000'])
if argv1=="Chromosomes_LIRSat1_35":
    target_list.append(['V037.hap2','Chr7:14382973-14710457'])
    target_list.append(['V032.hap2','Chr1:15428505-15829144'])
    
##Looking at similarities between Tekay on different chromosomes
if argv1=="Chromosomes_LIRSat2": 
    target_list.append(['V105.hap1','Chr1:14900000-15100000'])   
    target_list.append(['V079.hap1','Chr19:17154840-17626446'])   
    target_list.append(['V030.hap1','Chr4:13518846-13763462'])   
    target_list.append(['V067.hap1','Chr2:13312200-13552618'])   
    target_list.append(['V123.hap1','Chr4:13592321-13830091'])   
    target_list.append(['V123.hap1','Chr9:18432955-18697537'])   
    target_list.append(['V124.hap2','Chr4:13428189-13640340']) 
    target_list.append(['V055.hap1','Chr9:17469557-17771543']) 

##Looking at similarities between Tekay on different chromosomes, LIRs on different chromosomes in woolly grape and Vitis bryoniifolia seem very similar
if argv1=="Chromosomes_LIRSat_ret":      ##Add 50kbp to both sides of the S25 plot compared to moddot
    target_list.append(['WoollyGrape_hap1','Chr1:14738929-15197386'])
    target_list.append(['WoollyGrape_hap1','Chr2:15350993-15832030'])
    target_list.append(['WoollyGrape_hap1','Chr3:15570654-16055187'])
    target_list.append(['WoollyGrape_hap1','Chr4:13690850-14146219'])
    target_list.append(['WoollyGrape_hap1','Chr5:14100576-14682756'])
    target_list.append(['WoollyGrape_hap1','Chr6:12766579-13340046'])
    target_list.append(['WoollyGrape_hap1','Chr7:14568553-14921377'])
    target_list.append(['WoollyGrape_hap1','Chr8:8791760-9306061'])
    target_list.append(['WoollyGrape_hap1','Chr9:17643900-18151280'])  
    target_list.append(['WoollyGrape_hap1','Chr12:11799656-12431808'])
    target_list.append(['WoollyGrape_hap1','Chr13:13251853-13616766'])
    target_list.append(['WoollyGrape_hap1','Chr14:16565135-16424515'])
    target_list.append(['WoollyGrape_hap1','Chr17:14300851-14915357'])
    target_list.append(['WoollyGrape_hap1','Chr18:16485861-16977372'])
    target_list.append(['WoollyGrape_hap1','Chr19:18946536-19433157'])
    
if argv1=="Chromosomes_LIRSat_bry": 
    target_list.append(['V126.hap2','Chr1:14384736-14863897'])
    target_list.append(['V126.hap2','Chr2:14592428-15130566'])
    target_list.append(['V126.hap2','Chr3:15258845-15877811'])
    target_list.append(['V126.hap2','Chr4:13956925-14376580'])
    target_list.append(['V126.hap2','Chr5:15078558-15669106'])
    target_list.append(['V126.hap2','Chr6:12010223-12569057'])
    target_list.append(['V126.hap2','Chr7:14847465-15293715'])
    target_list.append(['V126.hap2','Chr8:7790307-8223339'])
    target_list.append(['V126.hap2','Chr9:17432607-17687234'])  
    target_list.append(['V126.hap2','Chr12:12175803-12899879'])
    target_list.append(['V126.hap2','Chr13:13781732-14322746'])
    target_list.append(['V126.hap2','Chr14:14894034-15529422'])
    #target_list.append(['V126.hap2','Chr17:14300851-14915357'])  # Seems to indeed lack LIRSat
    target_list.append(['V126.hap2','Chr18:16444108-17145390'])
    target_list.append(['V126.hap2','Chr19:17585851-18207892'])
  
    
    
if argv1=="part_all_LIRinterarray":  
    with open('../samples_satellite/2_good_regions_interarray_LIR','r') as  f:
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            sample	=eachline_arr[0]
            chromosome=eachline_arr[1]
            LIR_interarrays=eachline_arr[13]
            if LIR_interarrays=='NA':continue
            LIR_interarrays_list=LIR_interarrays.split('|')
            for one in LIR_interarrays_list:
                arr=one.split('-')
                if len(arr)!=2:print(arr);print('error,238')
                start,end=arr
                #start,end=int(start),int(end)
                border1=int(start)-50000
                border2=int(end)+50000
                if os.path.exists(f"./chr2TEsorter/auto_{sample}:{chromosome}:{border1}-{border2}/zzz3——cenblast_plot_auto_{sample}:{chromosome}:{border1}-{border2}.pdf")==True:continue
                print(f"zzz2——moddotplot_plot5_auto_{sample}:{chromosome}:{border1}-{border2}.pdf")
                #sys.exit()
                subprocess.run([f"python ./chr2TEsorter.py part_all_LIRinterarray_auto {sample} {chromosome} {border1} {border2}"], shell=True)              

if argv1=="part_all_LIRinterarray_chr18":  
    with open('../samples_satellite/chr18_LIR','r') as  f:
        next(f)
        for line in f:
            eachline_arr=line.strip().split('\t')
            if len(eachline_arr)!=4:continue
            sample,start,end,one_type=eachline_arr
            if one_type=='.':continue
            start,end=float(start),float(end)
            start=10000000+start*1000000
            end=10000000+end*1000000
            chromosome='Chr18'
            if one_type=='.':continue
            border1=int(start)-50000
            border2=int(end)+50000
            #if os.path.exists(f"./chr2TEsorter/auto_{sample}:{chromosome}:{border1}-{border2}/zzz3——cenblast_plot_auto_{sample}:{chromosome}:{border1}-{border2}.pdf")==True:continue
            print(f"zzz2——moddotplot_plot5_auto_{sample}_{one_type}:{chromosome}:{border1}-{border2}.pdf")
            #sys.exit()
            subprocess.run([f"python ./chr2TEsorter.py part_all_LIRinterarray_auto {sample} {chromosome} {border1} {border2}"], shell=True)    

if argv1=="part_all_LIRinterarray_PN40024":  
    chr_num_list=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
    for one_chrnum in chr_num_list:
        #if os.path.exists(f"./chr2TEsorter/auto_{sample}:{chromosome}:{border1}-{border2}/zzz3——cenblast_plot_auto_{sample}:{chromosome}:{border1}-{border2}.pdf")==True:continue
        #print(f"zzz2——moddotplot_plot5_auto_{sample}_{one_type}:{chromosome}:{border1}-{border2}.pdf")
        #sys.exit()
        subprocess.run([f"python ./chr2TEsorter.py part_all_LIRinterarray_auto PN40024 Chr{one_chrnum} 1 99999999"], shell=True)    
            
if argv1=="part_all_LIRinterarray_auto":  
    print('Automatically executed sub-step')
    auto_i='yes'
    sample=argvs[2]
    chromosome=argvs[3]
    start=int(argvs[4])
    end=int(argvs[5])
    target_list.append([sample,f'{chromosome}:{start}-{end}'])        
    argv1=f"auto_{sample}:{chromosome}:{start}-{end}"        
print('target_list:')
print(target_list)    
print('')   
    

    
def set_ulimit_nproc(value):
    try:
        # Only takes effect for the current shell session (cannot affect parent process)
        subprocess.run(f"ulimit -u {value}", shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to call ulimit: {e}")    
print("Current nproc limit:", os.system("ulimit -u"))
set_ulimit_nproc(2061163)
print("Current nproc limit:", os.system("ulimit -u"))
    
if len(argvs)==2 or auto_i=='yes': 
    part=argv1
    subprocess.run([f"mkdir ./chr2TEsorter/{part}"], shell=True)  
    
    if os.path.exists(f'./chr2TEsorter/{part}/0_info')==True:
        print('Check if already run')
        with open(f'./chr2TEsorter/{part}/0_info','r') as f:
            old_value=f.read().strip()
            if old_value!=f"{target_list}":
                print('Same name, coordinates changed')
                subprocess.run([f"rm -r ./chr2TEsorter/{part}"], shell=True)  
                subprocess.run([f"mkdir ./chr2TEsorter/{part}"], shell=True) 
            else:
                print('Coordinates unchanged, some completed steps can be omitted')
    
    with open(f'./chr2TEsorter/{part}/0_info','w') as f:
        f.write(f"{target_list}")

    subprocess.run([f"mkdir ./chr2TEsorter/{part}/0_prepare"], shell=True)  
    subprocess.run([f"mkdir ./chr2TEsorter/{part}/1_TEsorter"], shell=True)  
    subprocess.run([f"mkdir ./chr2TEsorter/{part}/2_moddotplot"], shell=True)  
    snapgene_feature_list=[]
    with open(f"./chr2TEsorter/{part}/1_TEsorter/snapgene_feature.bed",'w') as f3:
        with open(f"./chr2TEsorter/{part}/1_TEsorter/all.info",'w') as f2:
            f2.write(f"sample\tregion\tii\tkk\tTEsorter_class3\tdomain\tTEtype_dom\tone_real_start\tone_real_end\tstrand\tshift\tinfo\trelative_start\trelative_end\n")
            ii=0
            largetarget_relative_start=1
            seq_all_list=[]
            #print(f"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            #print(target_list)
            for one in target_list:
                #print(f"{one}!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                ii+=1
                sample,pos_str = one
                chromosome=pos_str.split(':')[0]
                pos=pos_str.split(':')[1]
                largetarget_start=int(pos.split('-')[0])
                largetarget_end=int(pos.split('-')[1])
                length=largetarget_end-largetarget_start+1
                if ii>1:
                    largetarget_relative_start+=last_len+interval_len    ### Add an interval sequence of 100000 in the middle for plotting
                largetarget_relative_end= largetarget_relative_start+ length-1  
                last_len= int(length)   
                print()
                print([sample,chromosome,pos])
                subprocess.run([f"samtools faidx /home/lain/aaa_data/run0/new_work_dir/chr2EDTA/0_prepare/{sample}/{sample}.fasta {chromosome}:{pos} > ./chr2TEsorter/{part}/0_prepare/{ii}  "], shell=True)  
                #sys.exit()
                seq_one=''
                with open(f'./chr2TEsorter/{part}/0_prepare/{ii}','r') as f:
                    next(f)
                    for line in f:
                        seq_one+=line.strip()
                seq_all_list.append(seq_one)
                
                seq_one_len=len(seq_one)
                #print(seq_one)
                #Window size 10000, step size 5000, i.e., overlapping regions, last part merged into previous one
                def split_sequence(seq_one, window_size, step_size):  #
                    new_list = []
                    length = len(seq_one)
                    
                    # If sequence is empty, return empty list directly
                    if length == 0:
                        return new_list
                    
                    # Calculate starting positions
                    start = 0
                    while start < length:
                        end = start + window_size
                        chunk = seq_one[start:end]  # Current window segment
                        
                        # If it's the last block and shorter than window_size, and new_list is not empty
                        if end >= length and new_list:
                            new_list[-1] += chunk  # Merge into previous block (string concatenation)
                        else:
                            new_list.append(chunk)  # Otherwise add directly
                        
                        start += step_size  # Slide window
                    
                    return new_list
                    
                result = split_sequence(seq_one,window_size,step_size) 

                kk=0
                kk_list=[]
                for one_part in result:
                    kk+=1
                    kk_list.append(kk)
                    
                ##Whether to run TEsorter    
                if TEsorter_step=='yes':    
                    def run_kk(kk):
                        #print(kk)
                        #if  os.path.exists(f"./chr2TEsorter/{part}/1_TEsorter/{ii}_{kk}")==True:
                        #    subprocess.run([f"rm -r ./chr2TEsorter/{part}/1_TEsorter/{ii}_{kk}"], shell=True)
                        if os.path.exists(f"./chr2TEsorter/{part}/1_TEsorter/{ii}_{kk}/input.rexdb-plant.dom.gff3")==True:
                            print(f'Skipping {ii}_{kk}')    
                        else:     
                            if os.path.exists(f"./chr2TEsorter/{part}/1_TEsorter/{ii}_{kk}")==False:
                                subprocess.run([f"mkdir ./chr2TEsorter/{part}/1_TEsorter/{ii}_{kk}"], shell=True)  
                            with open(f"./chr2TEsorter/{part}/1_TEsorter/{ii}_{kk}/input",'w') as f:
                                one_part=result[kk-1]
                                f.write(f">\n{one_part}\n")
                            ########    
                            os.chdir(f"./chr2TEsorter/{part}/1_TEsorter/{ii}_{kk}")
                            
                            subprocess.run([f"source /home/lain/miniconda3/etc/profile.d/conda.sh && conda activate EDTA && TEsorter -p 1 -db rexdb-plant input > log 2>&1"], shell=True)  
                            
                            os.chdir(f"../../../../")
                    with Pool(processes=20) as pool:
                        # Use imap to get results one by one
                        for i, def_result in enumerate(pool.imap(run_kk, kk_list), start=1):
                            # Process results here, e.g., store or print
                            progress = (i / len(kk_list)) * 100
                            sys.stdout.write(f"Progress: {progress:.2f}% {str(i)} / {str(len(kk_list))}#####################################################################\n")
                            sys.stdout.flush()    
             
                ######################            
                kk=0
                for one_part in result:
                    kk+=1            
                    with open(f"./chr2TEsorter/{part}/1_TEsorter/{ii}_{kk}/input.rexdb-plant.dom.gff3",'r') as f:
                        for line in f:
                            eachline=line.strip()
                            eachline_arr=eachline.split('\t')
                            database,one_type,start,end,score,strand,shift,info=eachline_arr
                            real_step_pos_start=largetarget_start+(kk-1)*step_size
                            one_real_start=real_step_pos_start+int(start)-1
                            one_real_end=real_step_pos_start+int(end)-1
                            relative_start=largetarget_relative_start+(one_real_start-largetarget_start)
                            relative_end=largetarget_relative_start+(one_real_end-largetarget_start)
                            TEsorter_class3=info.split(';')[0].split('/')[-1].split(':')[0]
                            TEsorter_class3_dom=info.split(';')[0].split('/')[-1].split(':')[1]
                            dom=TEsorter_class3_dom.split('-')[1]
                            
                            f2.write(f"{sample}\t{pos_str}\t{ii}\t{kk}\t{TEsorter_class3}\t{TEsorter_class3_dom}\t{TEsorter_class3}_{dom}\t{one_real_start}\t{one_real_end}\t{strand}\t{shift}\t{info}\t{relative_start}\t{relative_end}\n")
                            
                            feature_one=f"{relative_start}\t{relative_end}_{TEsorter_class3}_{dom}"
                            if feature_one not in snapgene_feature_list:
                                f3.write(f"AAA\t{relative_start}\t{relative_end}\t{TEsorter_class3}_{dom}\n")
                                snapgene_feature_list.append(feature_one)
    print('len(seq_all_list)')                
    print(len(seq_all_list))
    seq_all=interval.join(seq_all_list)
    with open(f'./chr2TEsorter/{part}/2_moddotplot/seq_all.fa','w') as f:
        f.write(f">seq_all_1\n{seq_all}\n")  #>seq_all_2\n{seq_all}
    print(f'Length of synthesized sequence (including interval sequence): {len(seq_all)}')
    
    ##Plot based on info
    if 1==1:
        R_txt=f'''
            library(ggplot2)
            library(dplyr)
            
            print("")
           
              input_file=read.table('all.info', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
            
            print({largetarget_relative_end})
            # Create plot object
            p <- ggplot() 
            p <- p + geom_segment( aes(x = 1, y = 9 , xend = {largetarget_relative_end}, yend = 9),color='black',linewidth = 1)
            p <- p + geom_rect(data = input_file, aes(xmin = relative_start, xmax = relative_end, ymin = 9, ymax = 10, fill = TEsorter_class3))
            p =p+ylim(0, 20)
                #scale_color_manual(values = color_values, drop = FALSE)+coord_fixed()+
            p =p+theme_classic() +         
                theme(
                  #axis.ticks.y = element_blank(),
                  #axis.text.y = element_blank(),
                  #legend.position = "none",
                  #axis.text.x = element_blank()
                ) 
             
              
              # Save as PDF
              pdf(file = paste0('plot', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
              print(p)
              dev.off()
        '''
        with open(f'./chr2TEsorter/{part}/1_TEsorter//plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = f'./chr2TEsorter/{part}/1_TEsorter//'
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../../')         
    if 1.5==1.5:    
        R_txt=f'''
            library(ggplot2)
            library(dplyr)
            
            print("")
           
              input_file=read.table('all.info', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
            TE_Ale=   input_file   %>%  filter(TEsorter_class3== "Ale")
            TE_Ikeros=   input_file   %>%  filter(TEsorter_class3== "Ikeros")
            TE_Athila=   input_file   %>%  filter(TEsorter_class3== "Athila")
            TE_CRM=   input_file   %>%  filter(TEsorter_class3== "CRM")
            TE_Galadriel=   input_file   %>%  filter(TEsorter_class3== "Galadriel")
            TE_Ogre=   input_file   %>%  filter(TEsorter_class3== "Ogre")
            TE_Reina=   input_file   %>%  filter(TEsorter_class3== "Reina")
            TE_Retand=   input_file   %>%  filter(TEsorter_class3== "Retand")
            TE_Tekay=   input_file   %>%  filter(TEsorter_class3== "Tekay")

            print({largetarget_relative_end})
            # Create plot object
            p <- ggplot() 
            p <- p + geom_segment( aes(x = 1, y = -90 , xend = {largetarget_relative_end}, yend = -90),color='black',linewidth = 1)
            
            p <- p + geom_rect(data = TE_Ale, aes(xmin = relative_start-2000, xmax = relative_end+2000, ymin = -100, ymax =-90 ),fill ='#ffca33',alpha=1)
            p <- p + geom_rect(data = TE_Ikeros, aes(xmin = relative_start-2000, xmax = relative_end+2000, ymin = -110, ymax =-100 ),fill ='#e0d868',alpha=1)
            p <- p + geom_rect(data = TE_Athila, aes(xmin = relative_start-2000, xmax = relative_end+2000, ymin = -120, ymax =-110 ),fill ='#00aac1',alpha=1)
            p <- p + geom_rect(data = TE_CRM, aes(xmin = relative_start-2000, xmax = relative_end+2000, ymin = -130, ymax =-120 ),fill ='#891555',alpha=1)
            p <- p + geom_rect(data = TE_Galadriel, aes(xmin = relative_start-2000, xmax = relative_end+2000, ymin = -140, ymax =-130 ),fill ='#ff3300',alpha=1)
            p <- p + geom_rect(data = TE_Ogre, aes(xmin = relative_start-2000, xmax = relative_end+2000, ymin = -150, ymax =-140 ),fill ='#666699',alpha=1)
            p <- p + geom_rect(data = TE_Reina, aes(xmin = relative_start-2000, xmax = relative_end+2000, ymin = -160, ymax =-150 ),fill ='#ff9999',alpha=1)
            p <- p + geom_rect(data = TE_Retand, aes(xmin = relative_start-2000, xmax = relative_end+2000, ymin = -170, ymax =-160 ),fill ='#00cc66',alpha=1)
            p <- p + geom_rect(data = TE_Tekay, aes(xmin = relative_start-2000, xmax = relative_end+2000, ymin = -180, ymax =-170 ),fill ='#c900c4',alpha=1)             
            #p =p+ylim(-100, 20)
                #scale_color_manual(values = color_values, drop = FALSE)+coord_fixed()+
            p =p+theme_classic() +         
                theme(
                  #axis.ticks.y = element_blank(),
                  #axis.text.y = element_blank(),
                  legend.position = "none",
                  #axis.text.x = element_blank()
                ) 
             
              
              # Save as PDF
              pdf(file = paste0('plot_strip', ".pdf"), width = 20/ 2.54, height = 5 / 2.54)
              print(p)
              dev.off()
        '''
        with open(f'./chr2TEsorter/{part}/1_TEsorter//plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = f'./chr2TEsorter/{part}/1_TEsorter//'
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)  
        subprocess.run([f'cp  plot_strip.pdf ../{part}_plot_strip.pdf'], shell=True)    
        os.chdir('../../../')               
        
        
        
    ##moddotplot plotting
    if 1==1:
        
        if moddotplot_step=='yes':
            if  os.path.exists(f'./chr2TEsorter/{part}/2_moddotplot/output_dir/seq_all_1.bed')==False:
                moddotplot_software='ulimit -u 100000 && source /home/lain/software/venv/ModDotPlot/venv/bin/activate &&  /home/lain/software/venv/ModDotPlot/venv/bin/moddotplot '
                input_file=f'./chr2TEsorter/{part}/2_moddotplot/seq_all.fa'
                os.chdir(f'./chr2TEsorter/{part}/2_moddotplot/')
                subprocess.run([f"{moddotplot_software} static  --fasta seq_all.fa --kmer 21 --output-dir output_dir --identity 80  --dpi 3000 --resolution 5000 --no-plot"], shell=True)
                os.chdir(f'../../../')
            
        ###########
        if 1.1==1.1:
            with open(f'./chr2TEsorter/{part}/2_moddotplot/output_dir/input_file','w') as f2 :
                f2.write(f"x\ty\tvalue\n")
                with open(f'./chr2TEsorter/{part}/2_moddotplot/output_dir/seq_all_1.bed','r') as f :
                    next(f)
                    i=0
                    for line in f:
                        i+=1
                        eachline_arr=line.strip().split('\t')
                        query_name,query_start,query_end,reference_name,reference_start,reference_end,perID_by_events=eachline_arr
                        query_start=            float(query_start)
                        reference_start=        float(reference_start)
                        perID_by_events=        float(perID_by_events)
                        if i==1: unit=float(query_end)-query_start+1
                        if perID_by_events<97:continue
                        x=(query_start-1)/unit
                        y=(reference_start-1)/unit
                        f2.write(f"{x}\t{y}\t{perID_by_events}\n")
                        f2.write(f"{y}\t{x}\t{perID_by_events}\n")
                        
            
            R_txt=r'''
            library(ggplot2)
            library(dplyr)
            
            print("")
            {
              input_file=read.table('input_file', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
            
            # Filter data
            filtered_input <- input_file %>% 
                filter(value >95 ) %>%
              mutate(
                    category = case_when(
                    value>99~ "99+",
                    value>98~ "98+",
                    value>97~ "97+",
                    value>96~ "96+",
                    value>95~ "95+"
                    ),
                  )
            B99=   filtered_input   %>%  filter(category== "99+")
            B98=   filtered_input   %>%  filter(category== "98+")
            B97=   filtered_input   %>%  filter(category== "97+")
            #B96=   filtered_input   %>%  filter(category== "96+")
            #B95=   filtered_input   %>%  filter(category== "95+")
            
            # Define color values
            color_values <- c(
                 '99+'='red',
                 '98+'='#ffcc00',
                 '97+'='#004d66'#,
                 #'96+'='#004d66'#,
                 #'95+'='#004d66'
            )
            
        
            
             
            # Create plot object
            p <- ggplot() +
                #geom_point(data = B95,aes(x = x, y = y, color = category),shape=16,stroke=0,size = 0.5, alpha =0.5 ) +
                #geom_point(data = B96,aes(x = x, y = y, color = category),shape=16,stroke=0,size = 0.5, alpha =0.6 ) +
                geom_point(data = B97,aes(x = x, y = y, color = category),shape=16,stroke=0,size = 0.5, alpha =0.7 ) +
                geom_point(data = B98,aes(x = x, y = y, color = category),shape=16,stroke=0,size = 0.5, alpha =0.8 ) +
                geom_point(data = B99,aes(x = x, y = y, color = category),shape=16,stroke=0,size = 0.5, alpha =0.9 ) +
                scale_color_manual(values = color_values, drop = FALSE)+coord_fixed()+
            theme_classic() +         
                theme(
                  #axis.ticks.y = element_blank(),
                  #axis.text.y = element_blank(),
                  #legend.position = "none",
                  #axis.text.x = element_blank()
                ) 
             
              
              # Save as PDF
              pdf(file = paste0('plot', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
              print(p)
              dev.off()
            
            }
                '''
            with open(f'./chr2TEsorter/{part}/2_moddotplot/output_dir/plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)  
            new_directory = f'./chr2TEsorter/{part}/2_moddotplot/output_dir/'
            os.chdir(new_directory)
            if 'Rtxt'=='R':
                subprocess.run(['Rscript plot.R'], shell=True)    
            os.chdir('../../../../')     
        if 1.2==1.2:
            with open(f'./chr2TEsorter/{part}/2_moddotplot/output_dir/input_file2','w') as f2 :
                f2.write(f"x\ty\tvalue\n")
                with open(f'./chr2TEsorter/{part}/2_moddotplot/output_dir/seq_all_1.bed','r') as f :
                    next(f)
                    i=0
                    for line in f:
                        i+=1
                        eachline_arr=line.strip().split('\t')
                        query_name,query_start,query_end,reference_name,reference_start,reference_end,perID_by_events=eachline_arr
                        query_start=            float(query_start)
                        reference_start=        float(reference_start)
                        perID_by_events=        float(perID_by_events)
                        if i==1: unit=float(query_end)-query_start+1
                        if perID_by_events<90:continue
                        x=(query_start-1)/unit
                        y=(reference_start-1)/unit
                        f2.write(f"{x}\t{y}\t{perID_by_events}\n")
                        f2.write(f"{y}\t{x}\t{perID_by_events}\n")
                    
            R_txt=r'''
            library(ggplot2)
            library(dplyr)
            
            print("")
            {
              input_file=read.table('input_file2', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
            
            # Filter data
            filtered_input <- input_file %>% 
                filter(value >90 ) %>%
              mutate(
                    category = case_when(
                    value>99~ "99+",
                    value>98~ "98+",
                    value>97~ "97+",
                    value>96~ "96+",
                    value>95~ "95+",
                    value>90~ "90+"
                    ),
                  )
            B99=   filtered_input   %>%  filter(category== "99+")
            B98=   filtered_input   %>%  filter(category== "98+")
            B97=   filtered_input   %>%  filter(category== "97+")
            B96=   filtered_input   %>%  filter(category== "96+")
            B95=   filtered_input   %>%  filter(category== "95+")
            B90=   filtered_input   %>%  filter(category== "90+")
            # Define color values
            color_values <- c(
                 '99+'='#ff00ff',
                 '98+'='#ffcc00',
                 '97+'='#ccff33',
                 '96+'='#66ff99',
                 '95+'='#ccffff',
                 '90+'='black'
            )
            
        
            
             
            # Create plot object
            p <- ggplot() +
                geom_point(data = B90,aes(x = x, y = y, color = category),shape=16,stroke=0,size = 0.5, alpha =0.1 ) +
                geom_point(data = B95,aes(x = x, y = y, color = category),shape=16,stroke=0,size = 0.5, alpha =0.5 ) +
                geom_point(data = B96,aes(x = x, y = y, color = category),shape=16,stroke=0,size = 0.5, alpha =0.6 ) +
                geom_point(data = B97,aes(x = x, y = y, color = category),shape=16,stroke=0,size = 0.5, alpha =0.7 ) +
                geom_point(data = B98,aes(x = x, y = y, color = category),shape=16,stroke=0,size = 0.5, alpha =0.8 ) +
                geom_point(data = B99,aes(x = x, y = y, color = category),shape=16,stroke=0,size = 0.5, alpha =0.9 ) +
                scale_color_manual(values = color_values, drop = FALSE)+coord_fixed()+
            theme_classic() +         
                theme(
                  #axis.ticks.y = element_blank(),
                  #axis.text.y = element_blank(),
                  #legend.position = "none",
                  #axis.text.x = element_blank()
                ) 
             
              
              # Save as PDF
              pdf(file = paste0('plot2', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
              print(p)
              dev.off()
            
            }
                '''
            with open(f'./chr2TEsorter/{part}/2_moddotplot/output_dir/plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)  
            new_directory = f'./chr2TEsorter/{part}/2_moddotplot/output_dir/'
            os.chdir(new_directory)
            if 'Rtxt'=='R':
                subprocess.run(['Rscript plot.R'], shell=True)    
            os.chdir('../../../../')               
            
        ##Annotate TE positions based on bed file
        if 1.3==1.3:
            with open(f"./chr2TEsorter/{part}/2_moddotplot/output_dir/seq_all_1.bed",'r') as f:
                next(f)
                for line in f:
                    eachline_arr=line.strip().split('\t')
                    print(eachline_arr)
                    moddotplot_window=int(eachline_arr[2])-int(eachline_arr[1])+1
                    break
                
            subprocess.run([f"cp ./chr2TEsorter/{part}/1_TEsorter/all.info ./chr2TEsorter/{part}/2_moddotplot/output_dir/input_file3_TE"], shell=True)
                    
            R_txt=f'''
            library(ggplot2)
            library(dplyr)
            
            print("")
         
              input_file=read.table('input_file2', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
            
            # Filter data
            filtered_input <- input_file %>% 
                filter(value >90 ) %>%
              mutate(
                    category = case_when(
                    value>99~ "99+",
                    value>98~ "98+",
                    value>97~ "97+",
                    value>96~ "96+",
                    value>95~ "95+",
                    value>90~ "90+"
                    ),
                  )
            B99=   filtered_input   %>%  filter(category== "99+")
            B98=   filtered_input   %>%  filter(category== "98+")
            B97=   filtered_input   %>%  filter(category== "97+")
            B96=   filtered_input   %>%  filter(category== "96+")
            B95=   filtered_input   %>%  filter(category== "95+")
            B90=   filtered_input   %>%  filter(category== "90+")
            # Define color values
            color_values <- c(
                 '99+'='#ff00ff',
                 '98+'='#ffcc00',
                 '97+'='#ccff33',
                 '96+'='#66ff99',
                 '95+'='#ccffff',
                 '90+'='black'
            )
            color_values <- c(
                 '99+'='#000000',
                 '98+'='#000000',
                 '97+'='#000000',
                 '96+'='#000000',
                 '95+'='#000000',
                 '90+'='black',
                        #"Ale"="#ffca33",
                        #"Alesia"="#a57a00",
                        #"Angela"="#a57a00",
                        #"Bianca"="#a57a00",
                        #"Ikeros"="blue",
                        #"Ivana"="#a57a00",
                        #"SIRE"="#a57a00",
                        #"TAR"="#a57a00",
                        #"Athila"="#00aac1",
                        #"CRM"="#891555",
                        #"Galadriel"="#ff3300",
                        #"Ogre"="#666699",
                        #"Reina"="#ff9999",
                        #"Retand"="#00cc66",
                        "Tekay"="#c900c4",
                        
                        "Ale"="#b3d9ff",
                        "Alesia"="#b3d9ff",
                        "Angela"="#b3d9ff",
                        "Bianca"="#b3d9ff",
                        "Ikeros"="#b3d9ff",
                        "Ivana"="#b3d9ff",
                        "SIRE"="#b3d9ff",
                        "TAR"="#b3d9ff",
                        "Tork"="#b3d9ff",
                        "Athila"="#b3d9ff",
                        "CRM"="#b3d9ff",
                        "Galadriel"="#b3d9ff",
                        "Ogre"="#b3d9ff",
                        "Reina"="#b3d9ff",
                        "Retand"="#b3d9ff",
                        "LINE"="#b3d9ff",
                        "MuDR_Mutator"="#b3d9ff",
                        "non-chromo-outgroup"="#b3d9ff",
                        "pararetrovirus"="#b3d9ff"
            )            
            input_file3_TE=read.table('input_file3_TE', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
            input_file3_TE=input_file3_TE   %>%  filter(TEsorter_class3 %in% c("Ale","Alesia","Angela","Bianca","Ikeros","Ivana","SIRE","TAR","Tork","Athila","CRM","Galadriel","Ogre","Reina","Retand","Tekay"))
             
            # Create plot object
            p <- ggplot() 
                p <- p + geom_rect(data = input_file3_TE, aes(xmin = relative_start, xmax = relative_end, ymin = 0, ymax = {largetarget_relative_end}, fill = TEsorter_class3),alpha=1)
                p <- p +geom_point(data = B90,aes(x = x*{moddotplot_window}, y = y*{moddotplot_window}, color = category),shape=16,stroke=0,size = 0.5, alpha =0.1 ) +
                geom_point(data = B95,aes(x = x*{moddotplot_window}, y = y*{moddotplot_window}, color = category),shape=16,stroke=0,size = 0.5, alpha =0.5 ) +
                geom_point(data = B96,aes(x = x*{moddotplot_window}, y = y*{moddotplot_window}, color = category),shape=16,stroke=0,size = 0.5, alpha =0.6 ) +
                geom_point(data = B97,aes(x = x*{moddotplot_window}, y = y*{moddotplot_window}, color = category),shape=16,stroke=0,size = 0.5, alpha =0.7 ) +
                geom_point(data = B98,aes(x = x*{moddotplot_window}, y = y*{moddotplot_window}, color = category),shape=16,stroke=0,size = 0.5, alpha =0.8 ) +
                geom_point(data = B99,aes(x = x*{moddotplot_window}, y = y*{moddotplot_window}, color = category),shape=16,stroke=0,size = 0.5, alpha =0.9 ) +
                scale_color_manual(values = color_values, drop = TRUE)+scale_fill_manual(values = color_values, drop = TRUE)+coord_fixed()+
            theme_classic() +         
                theme(
                  #axis.ticks.y = element_blank(),
                  #axis.text.y = element_blank(),
                  #legend.position = "none",
                  #axis.text.x = element_blank()
                ) 
             
              
              # Save as PDF
              pdf(file = paste0('plot3', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
              print(p)
              dev.off()
            
    
                '''
            with open(f'./chr2TEsorter/{part}/2_moddotplot/output_dir/plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)  
            new_directory = f'./chr2TEsorter/{part}/2_moddotplot/output_dir/'
            os.chdir(new_directory)
            subprocess.run(['Rscript plot.R'], shell=True)    
            os.chdir('../../../../')               
            subprocess.run([f"cp ./chr2TEsorter/{part}/2_moddotplot/output_dir/plot3.pdf ./chr2TEsorter/{part}/zzz2——moddotplot_plot3_{part}.pdf "], shell=True) 
    
        if 1.4==1.4:
            with open(f'./chr2TEsorter/{part}/2_moddotplot/output_dir/input_file4','w') as f2 :
                f2.write(f"x\ty\tvalue\treal_relative_x_start\treal_relative_x_end\treal_relative_y_start\treal_relative_y_end\n")
                with open(f'./chr2TEsorter/{part}/2_moddotplot/output_dir/seq_all_1.bed','r') as f :
                    next(f)
                    i=0
                    for line in f:
                        i+=1
                        eachline_arr=line.strip().split('\t')
                        query_name,query_start,query_end,reference_name,reference_start,reference_end,perID_by_events=eachline_arr
                        query_start=            float(query_start)
                        reference_start=        float(reference_start)
                        perID_by_events=        float(perID_by_events)
                        if i==1: unit=float(query_end)-query_start+1
                        if perID_by_events<90:continue
                        x=(query_start-1)/unit
                        y=(reference_start-1)/unit
                        real_relative_x_start=int(round(x*moddotplot_window+1,0))
                        real_relative_x_end=int(round(real_relative_x_start+moddotplot_window-1))
                        real_relative_y_start=int(round(y*moddotplot_window+1))
                        real_relative_y_end=int(round(real_relative_y_start+moddotplot_window-1  ))                      
                        f2.write(f"{x}\t{y}\t{perID_by_events}\t{real_relative_x_start}\t{real_relative_x_end}\t{real_relative_y_start}\t{real_relative_y_end}\n")
                        f2.write(f"{y}\t{x}\t{perID_by_events}\t{real_relative_y_start}\t{real_relative_y_end}\t{real_relative_x_start}\t{real_relative_x_end}\n")

            R_txt=f'''
            library(ggplot2)
            library(dplyr)
            
            print("")
         
              input_file=read.table('input_file4', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
            
            # Filter data
            filtered_input <- input_file %>% 
                filter(value >90 ) %>%
              mutate(
                    category = case_when(
                    value>99~ "99+",
                    value>98~ "98+",
                    value>97~ "97+",
                    value>96~ "96+",
                    value>95~ "95+",
                    value>90~ "90+",
                    value>80~ "80+"
                    ),
                  )
            B99=   filtered_input   %>%  filter(category== "99+")
            B98=   filtered_input   %>%  filter(category== "98+")
            B97=   filtered_input   %>%  filter(category== "97+")
            B96=   filtered_input   %>%  filter(category== "96+")
            B95=   filtered_input   %>%  filter(category== "95+")
            B90=   filtered_input   %>%  filter(category== "90+")
            B80=   filtered_input   %>%  filter(category== "80+")
            # Define color values
            color_values <- c(
                 '99+'='#000000',
                 '98+'='#1a1a1a',
                 '97+'='#333333',
                 '96+'='#4d4d4d',
                 '95+'='#666666',
                 '90+'='#808080',
                 '80+'='#808080',
                 'Ale'='green',
                 'Athila'='red',
                 'CRM'='blue',
                 'Tekay'='purple'
            )
            
            input_file3_TE=read.table('input_file3_TE', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
            
             
            # Create plot object
            p <- ggplot() 
                p <- p + geom_rect(data = input_file3_TE, aes(xmin = relative_start, xmax = relative_end, ymin = 0, ymax = {largetarget_relative_end}, fill = TEsorter_class3),alpha=0.8)
                p <- p +
                #geom_rect(data = B80, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B90, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B95, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B96, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B97, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B98, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B99, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                scale_fill_manual(values = color_values, drop = FALSE)+
                coord_fixed()+
            theme_classic() +         
                theme(
                  #axis.ticks.y = element_blank(),
                  #axis.text.y = element_blank(),
                  #legend.position = "none",
                  #axis.text.x = element_blank()
                ) 
             
              
              # Save as PDF
              pdf(file = paste0('plot4', ".pdf"), width = 40/ 2.54, height = 40 / 2.54)
              print(p)
              dev.off()
            
    
                '''
            with open(f'./chr2TEsorter/{part}/2_moddotplot/output_dir/plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)  
            new_directory = f'./chr2TEsorter/{part}/2_moddotplot/output_dir/'
            os.chdir(new_directory)
            subprocess.run(['Rscript plot.R'], shell=True)    
            os.chdir('../../../../')      
        if 1.5==1.5:
            R_txt=f'''
            library(ggplot2)
            library(dplyr)
            
            print("")
         
              input_file=read.table('input_file4', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
            
            # Filter data
            filtered_input <- input_file %>% 
                filter(value >90 ) %>%
              mutate(
                    category = case_when(
                    value>99~ "99+",
                    value>98~ "98+",
                    value>97~ "97+",
                    value>96~ "96+",
                    value>95~ "95+",
                    value>90~ "90+",
                    value>80~ "80+"
                    ),
                  )
            B99=   filtered_input   %>%  filter(category== "99+")
            B98=   filtered_input   %>%  filter(category== "98+")
            B97=   filtered_input   %>%  filter(category== "97+")
            B96=   filtered_input   %>%  filter(category== "96+")
            B95=   filtered_input   %>%  filter(category== "95+")
            B90=   filtered_input   %>%  filter(category== "90+")
            B80=   filtered_input   %>%  filter(category== "80+")
            
            
            # Define color values
            color_values <- c(
                 '99+'='#000000',
                 '98+'='#1a1a1a',
                 '97+'='#333333',
                 '96+'='#4d4d4d',
                 '95+'='#666666',
                 '90+'='#808080',
                 '80+'='#808080',
                 'Ale_GAG'='#66ffff','Ale_PROT'='#00ffff','Ale_RT'='#33cccc','Ale_RH'='#0099cc','Ale_INT'='#3366cc',
                 'Athila_GAG'='#66ccff','Athila_PROT'='#33ccff','Athila_RT'='#00ccff','Athila_RH'='#0099cc','Athila_INT'='#336699',
                 'CRM_GAG'='#99ff99','CRM_PROT'='#66ff66','CRM_RT'='#33cc33','CRM_RH'='#009933','CRM_INT'='#006600','CRM_CHD'='#003300',
                 'Galadriel_GAG'='#ffffcc','Galadriel_PROT'='#ffff99','Galadriel_RT'='#ffff66','Galadriel_RH'='yellow','Galadriel_INT'='#cc9900','Galadriel_CHD'='#663300',
                 'Ogre_GAG'='#99ccff','Ogre_PROT'='#6699ff','Ogre_RT'='#3366ff','Ogre_RH'='#3333ff','Ogre_INT'='#000066',
                 'Reina_GAG'='#ffcccc','Reina_PROT'='#ff9999','Reina_RT'='#ff6666','Reina_RH'='#ff5050','Reina_INT'='#cc0000','Reina_CHD'='#993333',
                 'Retand_GAG'='#bfbfbf','Retand_PROT'='#8c8c8c','Retand_RT'='#595959','Retand_RH'='#333333','Retand_INT'='#0d0d0d',
                 'Tekay_GAG'='#ff99ff','Tekay_PROT'='#ff66ff','Tekay_RT'='#ff00ff','Tekay_RH'='#cc00cc','Tekay_INT'='#9900cc','Tekay_CHD'='#660066'
            )
            
            input_file3_TE=read.table('input_file3_TE', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
            TE_Ale=   input_file3_TE   %>%  filter(TEsorter_class3== "Ale")
            TE_Athila=   input_file3_TE   %>%  filter(TEsorter_class3== "Athila")
            TE_CRM=   input_file3_TE   %>%  filter(TEsorter_class3== "CRM")
            TE_Galadriel=   input_file3_TE   %>%  filter(TEsorter_class3== "Galadriel")
            TE_Ogre=   input_file3_TE   %>%  filter(TEsorter_class3== "Ogre")
            TE_Reina=   input_file3_TE   %>%  filter(TEsorter_class3== "Reina")
            TE_Retand=   input_file3_TE   %>%  filter(TEsorter_class3== "Retand")
            TE_Tekay=   input_file3_TE   %>%  filter(TEsorter_class3== "Tekay")
             
            # Create plot object
            p <- ggplot() 
                p <- p + geom_rect(data = TE_Ale, aes(xmin = relative_start, xmax = relative_end, ymin = -50000, ymax =0 , fill = TEtype_dom),alpha=1)
                p <- p + geom_rect(data = TE_Athila, aes(xmin = relative_start, xmax = relative_end, ymin = -150000, ymax =-100000 , fill = TEtype_dom),alpha=1)
                p <- p + geom_rect(data = TE_CRM, aes(xmin = relative_start, xmax = relative_end, ymin = -250000, ymax =-200000 , fill = TEtype_dom),alpha=1)
                p <- p + geom_rect(data = TE_Galadriel, aes(xmin = relative_start, xmax = relative_end, ymin = -350000, ymax =-300000 , fill = TEtype_dom),alpha=1)
                p <- p + geom_rect(data = TE_Ogre, aes(xmin = relative_start, xmax = relative_end, ymin = -450000, ymax =-400000 , fill = TEtype_dom),alpha=1)
                p <- p + geom_rect(data = TE_Reina, aes(xmin = relative_start, xmax = relative_end, ymin = -550000, ymax =-500000 , fill = TEtype_dom),alpha=1)
                p <- p + geom_rect(data = TE_Retand, aes(xmin = relative_start, xmax = relative_end, ymin = -650000, ymax =-600000 , fill = TEtype_dom),alpha=1)
                p <- p + geom_rect(data = TE_Tekay, aes(xmin = relative_start, xmax = relative_end, ymin = -750000, ymax =-700000 , fill = TEtype_dom),alpha=1)

                p <- p + geom_rect(data = TE_Ale, aes(xmin = relative_start, xmax = relative_end, ymin = -1100000, ymax =-1000000 , fill = TEtype_dom),fill ='#ffca33',alpha=1)
                p <- p + geom_rect(data = TE_Athila, aes(xmin = relative_start, xmax = relative_end, ymin = -1200000, ymax =-1100000 , fill = TEtype_dom),fill ='#00aac1',alpha=1)
                p <- p + geom_rect(data = TE_CRM, aes(xmin = relative_start, xmax = relative_end, ymin = -1300000, ymax =-1200000 , fill = TEtype_dom),fill ='#891555',alpha=1)
                p <- p + geom_rect(data = TE_Galadriel, aes(xmin = relative_start, xmax = relative_end, ymin = -1400000, ymax =-1300000 , fill = TEtype_dom),fill ='#ff3300',alpha=1)
                p <- p + geom_rect(data = TE_Ogre, aes(xmin = relative_start, xmax = relative_end, ymin = -1500000, ymax =-1400000 , fill = TEtype_dom),fill ='#666699',alpha=1)
                p <- p + geom_rect(data = TE_Reina, aes(xmin = relative_start, xmax = relative_end, ymin = -1600000, ymax =-1500000 , fill = TEtype_dom),fill ='#ff9999',alpha=1)
                p <- p + geom_rect(data = TE_Retand, aes(xmin = relative_start, xmax = relative_end, ymin = -1700000, ymax =-1600000 , fill = TEtype_dom),fill ='#00cc66',alpha=1)
                p <- p + geom_rect(data = TE_Tekay, aes(xmin = relative_start, xmax = relative_end, ymin = -1800000, ymax =-1700000 , fill = TEtype_dom),fill ='#c900c4',alpha=1)                
                p <- p +
                #geom_rect(data = B80, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B90, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B95, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B96, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B97, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B98, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B99, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                scale_fill_manual(values = color_values, drop = FALSE)+
                coord_fixed()+
            theme_classic() +         
                theme(
                  #axis.ticks.y = element_blank(),
                  #axis.text.y = element_blank(),
                  #legend.position = "none",
                  #axis.text.x = element_blank()
                ) 
             
              
              # Save as PDF
              pdf(file = paste0('plot5', ".pdf"), width = 40/ 2.54, height = 40 / 2.54)
              print(p)
              dev.off()
            
    
                '''
            with open(f'./chr2TEsorter/{part}/2_moddotplot/output_dir/plot.R','w',encoding='utf-8') as f:
                f.write(R_txt)  
            new_directory = f'./chr2TEsorter/{part}/2_moddotplot/output_dir/'
            os.chdir(new_directory)
            subprocess.run(['Rscript plot.R'], shell=True)    
            os.chdir('../../../../')          
            subprocess.run([f"cp ./chr2TEsorter/{part}/2_moddotplot/output_dir/plot5.pdf ./chr2TEsorter/{part}/zzz2——moddotplot_plot5_{part}.pdf "], shell=True) 
    
    #Display repetitive sequences within LIR
    if blast_cen_1305=='yes':  
        subprocess.run([f"mkdir ./chr2TEsorter/{part}/3_cenblast"], shell=True)  
        LTR_list=[]##tekaysWoollyGrape_hap1.LIR-LTR
        LTR_list.append(['crmLTR','TGATGGGGGAAAGCCTTTTAAAGTTTATTCAAGGAAGAAAGGGAAGACCTAAGTAAACTTTAAGTAGGATTAACATTAATTAGAATTGAATAGGGATTGGTATTTGTTGGAGTCTAATTAGGATTAGCTAGTTGAGTTATAATATAATTAGAATTTGAGTCATGATAGGTTATTAGAGTCCTAGTAGACTTTGGATTTCTTAGAGAAGCCTATAAATAGGCTAATCAATGTAAATCAAAGGAAGGAATTTTGATGAATAATATTATACTTTCTTTCATTGCAAGGTTGCAACCCTCAATGGTGAGACTCCATTGATTTTCTTCTAGGTGAGACTCCTAGAAGGCCTTAGTGAGACTCTAAGGTTTTCCATCTTTTCTTCATTGTTTCTTCTTTCTCTCTATTTCTTATCTTATAATTTTCTCTACCATAAAGTTTATTCCTTGTTCCTCTCCTTACACCCTAAAAATAAAACCCTAGCCCACCTACCCTTGAGTGTAGGCAACCATCCTAGGGTTGTCCTACATCA'])
        LTR_list.append(['tekayLTR','TGTAACACCTAGATTATTTTAGTACTTAACTTAAGCTTGCTTAATTAGATATTAAAACTAGTTTAGTGCTAATCATGTTTAATTATGAATTAAACTTTGATTAAGGTTAAGTGGGATTAGTTAATGATCTAAGCATGCTTATTAGGTTCTTAGGGACTAATTAAGGTTATGAAAGCCTAAAGGACTAATGGGCAATTATGGGTTTCATTTTTGATAACTTGAAGGACTAAAGTGCAAATTTAGAAACTTGAGTTAGTGGAATGCCACCTCCTACTTGGTGGATTGGTGGCAGCCATGTGGGCTGCCACCTCATGCTTGGGTGCTTGGAGACCCCTTTATAAGGGGGGCTGCAGCTCGTTTTGCACCATTTCACCTGCAGCAACTTGGGTTAGAGAGAGAGACTAGAGAGAGAAAGTGTAGATTTGAGGTAAGCAAGTGTTTAATTTTGTAATTTTCGTTTATTTTGGTTCCTAATGGTATGTTGAGAAGGATTAATGGTAAAAATTGTGTAAAAGTTGAGTATAATTAGCTATAAATCTGTTTTAGTTAAAAATCACAATTTATTAAAGATTAAAGTATTTTAAGTTAAGTATTTTATTAAGGGTAAGATCAGCATAAATCTTTGTTTAATTTGGTTTGGTTGGAAAATAGATTAGTAAACATGTGATGAATTGGGTTATTTGGACACTTAGGATTTGTTCATAGGTTAGGCTTTAAGAAAATCAAAAGAATTAGTGTTTGGTAATTAATTAGGGAAATTAATATGCATTGTAAATATTGTTTAATTCAATTCATATTTGATTTGAAAACATGTAGGATATAGATTAGAGAATGTGAGAATTAAAAATTTGCAACAATAAACTAAAACCTAAATTATGTTGTTTCATTTCAGTTCCTCGTAATAAGCAAAGATCGCCTACTTAGAAGAAACACAAAAAGGAAGTCGGTGTAAGGCAGGGAATTTTATGCTATTCTATGGTGTATCTGATTTCTTTCCTTGAATTATTTGTTGGAATTTTATGTTGTTTCTAGGATTTCATAAAATCTTTTAAAGTGTCATTGCATCACGGTTTATTGCATTGAAATGTTACCATGATGTTGGAACGTATATTGGTATGAAGATTCGTGGTTTTGTGTTGCCGAAATGGTTTATGGTGTGATTGCATTGCAAGAAGGATATGGAATTATCACAGATCATATTCAAAGAATTTATTGGTGGCGAAGTGCTTATGTGATGCTACAGGCCTTATCATCTGATATATTGTTTATGTGGGACCGTGGGTCCTTGGGTGAAAGTCCCTAAAGCCCTCGAGGAAGACACTCCGATGTGGGTGTATGGGGTTGGTCCTGCCCCTGGGTTTTAGTCCCTAAAGACACGGGGTTGGTCCTGCCCTTGGGTATGAGTCCCAAAAGACACGGGGTTGGTCCTGCCCTTGGGTATGAGTCCCTAAAGACACGGGGTATGACTTGCCCTTGGGTGATGTCCCAGAATAGTCATTATTATATTGATCAAGTATCTTGTGGAGCATTGATATCTGCTGTGTACATTGTCGGGAGTCAGCAATTTATCAGTTGTGAAAGGATGGATGGGGAAAAGCAATGATGAAACCAACAAACACATGCATACATGTTTGACATGATTACACATTTGTTAATGATTATGAAATATTTATCATGCATGTTATTAAACGTTTGATTCAAGGTTTTTAAGGGTATACTTAGGATGGTTATAAACTTTCCTACTGAGTTGTGAACTCACCCTATCCCTTCCACCTCTAGATGCAGGTCAGAAGTCCTATGCAGGAAAGAATGCTTGAGCAATGCTTTACTGTGATTGGATGCTGTTTGCTCATATTGAAAGTCTTTTGAGGCCATGCCTTTTGTATAAAGACTAAACCTTTTGTTGGAATGATGTAATGTACATATTTGTCAGATACACTTGTAGTATGGGCTACCCGTAGTAAACAATGGGGTTTGGTATATGTAAATTAATGTTGTACAAGTTTCTTTTGGGAAACAAAACATATTTTGTTCACTTCATAGTCACAATGTTTAATACAGGATGTACTCTTTATAACAGGTTTTCCAAATCCTCTTTTTGAACAAATTCAAGTAGGAACTCAACATTTAACTTTTGGAAAGCACTTATACTGTATTGAGTATCAATTGCTTAGTTGAAGAATAAATAATAATAATAATAAGAAAAAAAAAATAATGGGGTGTGACA'])
        LTR_list.append(['LIR-Sat1','TTGGACCATAATTCTTGGTTCTATCTTTGTTTATTCTTAGTTCCTTTAGTTACCCATCATTGGGTTGCTATCCTAGTTTTAATATGAGTTTATTCGTTTGTTCCCTCTTATTTCTCATCACTAGGATGTCATCCTAGGTTCTATATCCAATTTGAGTTAATTTTGTGTTTGATCGTCTACATACCATCCTAGGTTCTACCTAAGATTAGTTGATTGTTTGTTATGTTACTCGTCATTGTGATATCATCCTAATTTCTATCTAGGTTTGTTCTTAGTGCCTTTGGGTTCAAACATTAGGATCTAATCCTAGGTTTGGCCTATGTTCATTGTTTGTTCCTTCGATTTCCCATCATTGTGAACTATCTAAGCTTATTGATATGTTCCTTTGGTTTCTGATAATTAGGACGCCATCCTATGAGTTTTTCTTAAGCGTATGCACTATTCCTTGGGTACGTTCAGGACGTCCTCCTTCCTTGAGTTGTTTTTTGGCATTAGTTTGTTTTGGGCATTCCCCTTCATTAGGTTTATTTAGGGATCCTCTTTCCTTGAGTGTTATAAGGGTATTCCCATTCCACTACTTTTTTAGGGCATTCCTGTTCATTGGGCTATCAAACTACCTTGATTTACTTTAGTCCATCATGCTTCCTTCAGTTTCTTAGGTCATCCATCTTCATGGAGATTTTTAATTCATTCCCTTTCCTAGAGTTTGTTTACACCAACCAACTTCCTTCAATTCTTTTAGGGCATACCCTTTCCTAGGGATTACTTAGCAGACTCCCACTTCCTTTAGTGTAGTAGGACATCCCCATCTAGTAGTTTGTTTAGCGCAACTCCATTTCTTGAGGTTTTTTTAGCTCAACCGTTTTTGTTGAGATTGGGCATCCCCCTTGGTAGAGTTTGTTTGGGGCATGCCCTTTCTATGAGTCTTATTAGTGAATCCCCCTTCTTCATGTTCATTTACAGCGTTCTACTTCCATTGGCTTTTTAAGGGCATCCCACTTACTCGAATTTGCTTAGAGCATCCCCCTTCCTCGACTTTGTTTGGAGAATTCCCTTCATTAAGTTTATTAAGGGCATCCCACTTCCTTAAGTTGGTTTAGGGTATCGTTCGTCCTTCACCCCATAATCAGGATGTCATCCTATGTTCTATCTAAACATAATCATGTGTTCCCGAGATTCCTCATCATTGGAACACCATTCATAGTTCTATTCATTTTCATTCTTAGTTTATTTAGTTTCACGTTGTTGAGATGCCATCCTAGGTTCTATCTAAGCTAATTCATTTGTTCTTTTGATTGCTCATCA'])
        LTR_list.append(['LIR-Sat2',"CTTAGTTTAGTTTTTAAGCATTATTAGTTTGCTTACTTTCATGATATTATACTAAGTTTAGTTTCTTTCAAAAAGTCAGTTAATCTCATGTTTAGTTTAGTTTTGAAAAATTATTATTTTGGTTTCTTTGTTGATATTATGCTTACTTTAGTTTTTGTTTTTTGTTTTTTTTAAAAAAATCATTAGTTTGGTCATTCCTTTCACCTCATGTTTAGTTTAGTTTATAAAAATTATTAGTTTGGTTGCTTTGTTGATATTATACTTAGCTTAGTTTCTTTTAAAAATCATTAGTTATGTGGTTATTCCATTAATCTCATGTTTAGTTTAGTTTATTTTATAAAAAATTATTAATTTGGGTGCTTTGTTGATATTGTGCTTAGTTTAGT"])###TTAGTTTAGTTTCTTGTAAAAATCATTAGTAATGTGGTTATTCCATTAATCTCATGTTTAGTTTAGTTTTTAAAAAATTATGAGTTTTATTGCTTTGTTGATACTATGCTTAGTTACTTATAAAAATCATTAGTTTGGTTATTCTTTAGATCGCTTGC
        LTR_list.append(['LIR-Sat3','TCCAAGTCAATTTGTTAATAAATAAATCATTTTTCAAAATAAAAAGCTAGTTTGCCATCGAGTTGAAGCGCATGAGCCGAAATGCGGGGTCCACAAGTATTGTAATTGATATTTAAGCTAAAAATTTAATCTCCATGTACCCAAAGATTGGGTTGTGCCAAATTTATTTTTAACATGAATCCGGCCCATTGGCTTGCCCAACTTAAATTTTGGGCCTCACCTCAAATAAAAGCGAAGCTCGAACCTTAACCACTTAAAAATATTTGGGCTACGTTGGGGTAAGAAATAAATTCCTCATAATAGCCATCCAACCAAGGGAATGAAAGTACCCAATGGGCTTTCATTTTAGTTTGGCTTTTTCGGTCTTTTGTTTGTTGGAATGGAATTACCCAATGGGCTTGTATTTTTAGTTGGACCCTCCAATTATTTTTTTAACCTCATCCTAAACACCCCAAGGACCCAAAAGAAAAAAAAAATTGTTTTCTTTGAGAAGTATATCCCAACGCCTCCTAAATAAATTGTTAATATTTTACTATAAACTTTTATAAATAGAAATTTCAATATTATAATAAATAAATAAGGAGTCACACTTAGGGTTTTCTTTCTTATTTTGTTTACCCTTTAAAAATAAAACAAAAATAAGTGACGAC'])
        LTR_list.append(['cen240',"ATGAGAAACCAAAGGAACAAACAAATAAGGTTATATAGAACATAAGATGACGTCCTAATGATGAGATGAGAAATAAAAGGACCAAAGAGTAAACTTAGATAGAACCTATGATGTTGTCCCAATGAAAAGAATCCAAAGGACCAAATGATTAAGCTTAGATAGAATCGAGGATGGTGTCCCCATGATGAGAAACTAAATGAACAAAGAAGATACTTAGAACTTAGGATGATGTCCCATTGATGAGAAACTAAATGAACAAAGAAGATACTTAGAACTTAGGATGATGTCCCATTGATGAGAAACCAAAGGAACAAATGATTAAGCTTAGAGAGAACCTAGAATGGTATCCCATTGATGAGAAATGAAATGAAAAAGAATGAAGTAGGAGAAAACCTAGGGTGGTGTCCTAATGATGAGAAACTAAATGAAGAAAGAAGAAACTTAGATAGAGCCTAGGATGGTGTCCTAAAGATGAGAAACCAAAGGAACAAACAATTAAACTTAGATAGAACCTAGGACGGCATTTCAATGATAAGAAACCAAGGGAAAAAAAACGAATATTCTTAGATAGTACCTAAGAGTGTGTCCCAATGATGATAAATGAATGGAACTAAGAATAAACTTATATAGAAGTTATGATGGAGTTCCAATGATGAGAAAACAAAGGAACAAAGAATAAACTTAAATAGAACCTAGAATTGTAAGCCAATG"])
        
 

        with open(f"./chr2TEsorter/{part}/3_cenblast/0_monomer.fa",'w') as f:
            for one in LTR_list:
                LTR_name,LTR_seq=one
                f.write(f">{LTR_name}\n{LTR_seq}\n")
        subprocess.run([f"makeblastdb -dbtype nucl -in ./chr2TEsorter/{part}/2_moddotplot/seq_all.fa -out ./chr2TEsorter/{part}/3_cenblast/0_seq_all"], shell=True)          
        cmd=f"blastn -outfmt 6 -query ./chr2TEsorter/{part}/3_cenblast/0_monomer.fa -db ./chr2TEsorter/{part}/3_cenblast/0_seq_all -out ./chr2TEsorter/{part}/3_cenblast/1_blast.outfmt6"
        print(cmd)
        subprocess.run([f"{cmd}"], shell=True)     
        with open(f'./chr2TEsorter/{part}/3_cenblast/2_info.bed','w') as f3:
            with open(f'./chr2TEsorter/{part}/3_cenblast/2_info','w') as f2:
                f2.write(f"target_type\tstart\tend\tstrand\tlength\tbitscore\tevalue\tmonomer_pos1\tmonomer_pos2\n")
                with open(f'./chr2TEsorter/{part}/3_cenblast/1_blast.outfmt6','r') as f:
                    for line in f:
                        eachline_arr=line.strip().split('\t')
                        qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore=eachline_arr
                        length=float(send)-float(sstart)
                        if length>0:strand='+'
                        else:strand='-';tmp1=sstart;tmp2=send;sstart=tmp2;send=tmp1
                        length=int(abs(length))
                        f2.write(f"{qseqid}\t{sstart}\t{send}\t{strand}\t{length}\t{bitscore}\t{evalue}\t{qstart}\t{qend}\n")
                        f3.write(f'AAA\t{sstart}\t{send}\t{qseqid}|{strand}|{length}\n')
        R_txt=f'''
            library(ggplot2)
            library(dplyr)
            
            print("")
           
            input_file=read.table('2_info', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
            crmLTR   = input_file %>% filter(target_type == 'crmLTR' )
            tekayLTR   = input_file %>% filter(target_type == 'tekayLTR' )
            LIR_Sat1  = input_file %>% filter(target_type == 'LIR-Sat1' ) ###LIR-Sat1            
            LIR_Sat2  = input_file %>% filter(target_type == 'LIR-Sat2' ) ###LIR-Sat2
            LIR_Sat3  = input_file %>% filter(target_type == 'LIR-Sat3' ) ###LIR-Sat3            
            
            print({largetarget_relative_end})
            
            color_values=c(
                '+'='red','-'='blue'
                )
            # Create plot object
            p <- ggplot() 
            p <- p + geom_segment( aes(x = 1, y = 9 , xend = {largetarget_relative_end}, yend = 9),color='black',linewidth = 1)
            p <- p + geom_rect(data = crmLTR, aes(xmin = start, xmax = end, ymin = 100, ymax = 110, fill = strand))   +annotate("text", x = 10000, y = 105, label = "crmLTR", color = "#891555") 
            p <- p + geom_rect(data = tekayLTR, aes(xmin = start, xmax = end, ymin = 80, ymax = 90, fill = strand))   +annotate("text", x = 10000, y = 85, label = "tekayLTR", color = "#c900c4") 
            p <- p + geom_rect(data = LIR_Sat1, aes(xmin = start, xmax = end, ymin = 60, ymax = 70, fill = strand))   +annotate("text", x = 10000, y = 65, label = "LIR-Sat1", color = "black") 
            p <- p + geom_rect(data = LIR_Sat2, aes(xmin = start, xmax = end, ymin = 40, ymax = 50, fill = strand))   +annotate("text", x = 10000, y = 45, label = "LIR-Sat2", color = "black") 
            p <- p + geom_rect(data = LIR_Sat3, aes(xmin = start, xmax = end, ymin = 20, ymax = 30, fill = strand))  +annotate("text", x = 10000, y = 25, label = "LIR-Sat3", color = "black")             
           
            p =p+ylim(-10, 130)+scale_fill_manual(values = color_values, drop = FALSE)
            p =p+theme_classic() +         
                theme(
                  #axis.ticks.y = element_blank(),
                  #axis.text.y = element_blank(),
                  #legend.position = "none",
                  #axis.text.x = element_blank()
                ) 
             
              
              # Save as PDF
              pdf(file = paste0('plot', ".pdf"), width = 20/ 2.54, height = 20 / 2.54)
              print(p)
              dev.off()
        '''
        with open(f'./chr2TEsorter/{part}/3_cenblast/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = f'./chr2TEsorter/{part}/3_cenblast/'
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../../') 
        subprocess.run([f"cp ./chr2TEsorter/{part}/3_cenblast/plot.pdf ./chr2TEsorter/{part}/zzz3——cenblast_plot_{part}.pdf "], shell=True) 
        
    ##Plot a figure of all blocks
    if 1==1:
        subprocess.run([f"mkdir ./chr2TEsorter/{part}/4_satellite"], shell=True)  
        for one in target_list:
            sample,pos_str = one
            chromosome=pos_str.split(':')[0]
            pos=pos_str.split(':')[1]
            os.chdir('../')
            subprocess.run([f"python ./2-samples_satellite9.py step14.22_partall {sample} {chromosome} {pos}"], shell=True)  
            os.chdir('./new_work_dir')
            subprocess.run([f"mv  ../samples_satellite/14_chr_image/22_partall/{sample}:{chromosome}:{pos} ./chr2TEsorter/{part}/4_satellite/"], shell=True)  
            subprocess.run([f"cp  ./chr2TEsorter/{part}/4_satellite/{sample}:{chromosome}:{pos}/{sample}:{chromosome}:{pos}.pdf  ./chr2TEsorter/{part}/zzz999——{sample}:{chromosome}:{pos}.pdf "], shell=True)  
     
    ##Plot a figure summarizing all
    if 1==1:
        subprocess.run([f"mkdir ./chr2TEsorter/{part}/5_SUM"], shell=True)  
        subprocess.run([f"cp ./chr2TEsorter/{part}/2_moddotplot/output_dir/input_file4 ./chr2TEsorter/{part}/5_SUM/0_moddotplot"], shell=True)  
        subprocess.run([f"cp ./chr2TEsorter/{part}/2_moddotplot/output_dir/input_file3_TE ./chr2TEsorter/{part}/5_SUM/0_TEdomains"], shell=True)   
        subprocess.run([f"cp ./chr2TEsorter/{part}/3_cenblast/2_info ./chr2TEsorter/{part}/5_SUM/0_satellites_by_blast"], shell=True)   

        R_txt=f'''
            library(ggplot2)
            library(dplyr)
            
            print("")
         
              input_file=read.table('0_moddotplot', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
            
            # Filter data
            filtered_input <- input_file %>% 
                filter(value >90 ) %>%
              mutate(
                    category = case_when(
                    value>99~ "99+",
                    value>98~ "98+",
                    value>97~ "97+",
                    value>96~ "96+",
                    value>95~ "95+",
                    value>90~ "90+",
                    value>80~ "80+"
                    ),
                  )
            B99=   filtered_input   %>%  filter(category== "99+")
            B98=   filtered_input   %>%  filter(category== "98+")
            B97=   filtered_input   %>%  filter(category== "97+")
            B96=   filtered_input   %>%  filter(category== "96+")
            B95=   filtered_input   %>%  filter(category== "95+")
            B90=   filtered_input   %>%  filter(category== "90+")
            B80=   filtered_input   %>%  filter(category== "80+")
            # Define color values
            color_values <- c(
                 '99+'='#000000',
                 '98+'='#1a1a1a',
                 '97+'='#333333',
                 '96+'='#4d4d4d',
                 '95+'='#666666',
                 '90+'='#808080',
                 '80+'='#808080',
                 'Ale'='green',
                 'Athila'='red',
                 'CRM'='blue',
                 'Tekay'='purple','+'='red','-'='blue'
            )
            
            input_file3_TE=read.table('0_TEdomains', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
            ######
            satellites_by_blast=read.table('0_satellites_by_blast', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
            tekays   = satellites_by_blast %>% filter(target_type == 'WoollyGrape_hap1.LIR-LTR' )
            cen1305  = satellites_by_blast %>% filter(target_type == 'cen1305' )
            cen650  = satellites_by_blast %>% filter(target_type == 'cen650' )
            cen158  = satellites_by_blast %>% filter(target_type == 'cen158' )
            cen240  = satellites_by_blast %>% filter(target_type == 'cen240' )
            
            print({largetarget_relative_end})
            

            # Create plot object
            p <- ggplot() 
            p <- p + geom_segment( aes(x = 1, y = 9 , xend = {largetarget_relative_end}, yend = 9),color='black',linewidth = 1)
            p <- p + geom_rect(data = tekays, aes(xmin = start, xmax = end, ymin = -10000, ymax = -5000, fill = strand))
            p <- p + geom_rect(data = cen650, aes(xmin = start, xmax = end, ymin = -15000, ymax = -10000, fill = strand))
            p <- p + geom_rect(data = cen158, aes(xmin = start, xmax = end, ymin = -20000, ymax = -15000, fill = strand))
            p <- p + geom_rect(data = cen1305, aes(xmin = start, xmax = end, ymin = -25000, ymax = -20000, fill = strand))            
            p <- p + geom_rect(data = cen240, aes(xmin = start, xmax = end, ymin = -30000, ymax = -25000, fill = strand))
             
            # Create plot object

                p <- p + geom_rect(data = input_file3_TE, aes(xmin = relative_start, xmax = relative_end, ymin = 0, ymax = {largetarget_relative_end}, fill = TEsorter_class3),alpha=0.8)
                p <- p +
                #geom_rect(data = B80, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B90, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B95, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B96, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B97, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B98, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                geom_rect(data = B99, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end, fill = category))+
                scale_fill_manual(values = color_values, drop = FALSE)+
                coord_fixed()+
            theme_classic() +         
                theme(
                  #axis.ticks.y = element_blank(),
                  #axis.text.y = element_blank(),
                  #legend.position = "none",
                  #axis.text.x = element_blank()
                ) 
             
              
              # Save as PDF
              pdf(file = paste0('plot_all', ".pdf"), width = 40/ 2.54, height = 40 / 2.54)
              print(p)
              dev.off()
            
    
                '''
        with open(f'./chr2TEsorter/{part}/5_SUM/plot.R','w',encoding='utf-8') as f:
            f.write(R_txt)  
        new_directory = f'./chr2TEsorter/{part}/5_SUM/'
        os.chdir(new_directory)
        subprocess.run(['Rscript plot.R'], shell=True)    
        os.chdir('../../')            
        
        
        
        
else:
    pos_arr=argvs[2].split('-')
    start=int(pos_arr[0])
    end=int(pos_arr[1])
    part=argv1
    with open(f'./chr2TEsorter/{part}/2_moddotplot/output_dir/input_file4_filter','w') as f2 :
        f2.write(f"x\ty\tvalue\treal_relative_x_start\treal_relative_x_end\treal_relative_y_start\treal_relative_y_end\n")
        with open(f'./chr2TEsorter/{part}/2_moddotplot/output_dir/input_file4','r') as f :
            next(f)
            for line in f:
                eachline=line.strip()
                eachline_arr=eachline.split('\t')
                real_relative_x_start,real_relative_x_end,real_relative_y_start,real_relative_y_end=eachline_arr[-4:]
                mark=''
                if int(real_relative_x_start)>start and int(real_relative_x_start)<end : mark='yes'
                elif int(real_relative_x_end)>start and int(real_relative_x_end)<end : mark='yes'
                elif int(real_relative_y_start)>start and int(real_relative_y_start)<end : mark='yes'
                elif int(real_relative_y_end)>start and int(real_relative_y_end)<end : mark='yes'
                if mark=='yes':
                    f2.write(eachline+'\n')
            
    R_txt=f'''
        library(ggplot2)
        library(dplyr)
        
        print("")
     
          input_file=read.table('input_file4_filter', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\\t')
        
        # Filter data
        filtered_input <- input_file %>% 
            filter(value >80 ) %>%
          mutate(
                category = case_when(
                value>99~ "99+",
                value>98~ "98+",
                value>97~ "97+",
                value>96~ "96+",
                value>95~ "95+",
                value>90~ "90+",
                value>80~ "80+"
                ),
              )
        B99=   filtered_input   %>%  filter(category== "99+")
        B98=   filtered_input   %>%  filter(category== "98+")
        B97=   filtered_input   %>%  filter(category== "97+")
        B96=   filtered_input   %>%  filter(category== "96+")
        B95=   filtered_input   %>%  filter(category== "95+")
        B90=   filtered_input   %>%  filter(category== "90+")
        B80=   filtered_input   %>%  filter(category== "80+")
        # Define color values
        color_values <- c(
             'Ale'='green',
             'Athila'='red',
             'CRM'='blue',
             'Tekay'='purple'
        )
        
        input_file3_TE=read.table('input_file3_TE', skip = 0, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = '\t')
        
         
        # Create plot object
        p <- ggplot() 
            p <- p + geom_rect(data = input_file3_TE, aes(xmin = relative_start, xmax = relative_end, ymin = {start}, ymax = {end}, fill = TEsorter_class3),alpha=0.8)
            p <- p +
            #geom_rect(data = B80, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end),fill="#808080")+
            geom_rect(data = B90, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end),fill="#808080")+
            geom_rect(data = B95, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end),fill="#666666")+
            geom_rect(data = B96, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end),fill="#4d4d4d")+
            geom_rect(data = B97, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end),fill="#333333")+
            geom_rect(data = B98, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end),fill="#1a1a1a")+
            geom_rect(data = B99, aes(xmin = real_relative_x_start, xmax = real_relative_x_end, ymin = real_relative_y_start, ymax = real_relative_y_end),fill="#000000")+
            #scale_fill_manual(values = color_values)+
            coord_fixed()+
            scale_x_continuous(limits = c({start}, {end})) + 
            scale_y_continuous(limits = c({start}, {end}))+
        theme_classic() +         
            theme(
              #axis.ticks.y = element_blank(),
              #axis.text.y = element_blank(),
              #legend.position = "none",
              #axis.text.x = element_blank()
            ) 
         
          
          # Save as PDF
          pdf(file = paste0('plot4_{start}_{end}', ".pdf"), width = 40/ 2.54, height = 40 / 2.54)
          print(p)
          dev.off()
        

    '''
    with open(f'./chr2TEsorter/{part}/2_moddotplot/output_dir/plot.R','w',encoding='utf-8') as f:
        f.write(R_txt)  
    new_directory = f'./chr2TEsorter/{part}/2_moddotplot/output_dir/'
    os.chdir(new_directory)
    subprocess.run(['Rscript plot.R'], shell=True)    
    os.chdir('../../../../')      





    
    
    
    
    
    
    
    
    
    
    
       
print("\n\n")        
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))      
time_end=timeit.default_timer()
print('All the running time: %.0f Seconds'%(time_end-time_start))