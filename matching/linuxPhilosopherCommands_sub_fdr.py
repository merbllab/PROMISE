# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 11:38:43 2018

@author: assafk
"""
import subprocess
import os
from time import localtime, strftime, clock
from argparse import ArgumentParser
import pathlib

tic = clock()


################ recive input from the user and detemine global parameters ###################
if __name__ == "__main__":
    parser = ArgumentParser()

    parser.add_argument('batch_folder', help="working folder were the raw file exist", metavar='string')
    parser.add_argument('data', help="path to refernce fasta file", metavar='string')
    parser.add_argument('-digest', dest="digest", default = "nonspecific", choices=["stricttrypsin","Trypsin","nonspecific","chymotrypsin"], help="type of digestion", metavar='string')
    parser.add_argument('-fdr', dest="fdr", default = 0.01, help="FDR for peptides", metavar=float)
    parser.add_argument('-mods', dest="mods", default = "M:15.9949", help="sub group FDR - detemine the middle group", metavar='string')
    args = parser.parse_args()

non_specific = True if args.digest == "nonspecific" else False
batch_folder = pathlib.Path(args.batch_folder)
pepXML_flie_list = []
#
#if args.raw_file == "all.raw":
#MSFragger_folder = "MSFragger"
#create a list of raw files 
for  file in os.listdir(batch_folder):
    if file.endswith(".pepXML"):
        pepXML_flie_list.append(pathlib.Path(batch_folder / file))
            
VERSION = 1.5
"""
version update:
    0.1 - basic run on linux for peptidePhilosopher
    0.2 - add support for different FDR 
    0.3 - change commands order and parameters according to Alexey recommandation 
    1.0 - lock version with PhD proposal 
    1.5 - update philosopher version to support subgroup FDR 

"""

#philosopher = "/home/labs/yifatlab/assafk/tools/MSFragger_Philosopher/philosopher_linux_amd64"
philosopher = "/home/labs/yifatlab/assafk/PROMISE/external_software/philosopher_v2.1.2_subFDR/philosopher"




print ()
print ("      ==============================================================")
print ("      |                    Linux  philosopher                      |")
print ("      |                       Version", VERSION,"                         |")
print ("      |               ",strftime("%a, %d %b %Y %H:%M:%S", localtime()), "                  |")
print ("      ==============================================================")
print ()

print("pipline will work on the following files:")
for p in pepXML_flie_list:    
    print (p.name)
print()


#for batch_folder in  batch_folder_list:
print("------------------------------------------------------------------------------")
print(f" working on folder: {batch_folder.name}")
print("------------------------------------------------------------------------------")
   
dataFile =  pathlib.Path(args.data)

# clean
command1 =  [philosopher,"workspace","--clean"]
print(" command : " +  " ".join(command1))

subprocess.run(command1 ,cwd =batch_folder,check= True)

command2 = [philosopher ,"workspace","--init"]
print(" command :" + " ".join(command2))
subprocess.run(command2,cwd =batch_folder,check = True)

# anotate
command3 = [philosopher,"database","--annotate", str(dataFile)] # + " --prefix rev_"
print(" command : " +  " ".join(command3))
subprocess.run(command3,cwd =batch_folder,check = True)

# peptide prophet per file
for p in pepXML_flie_list:
    if non_specific:
        command4 = [philosopher,"peptideprophet","--decoy","rev_","--enzyme","nonspecific", "--decoyprobs","--ppm","--nontt",
                    "--accmass","--nonparam","--expectscore","--database", str(dataFile),str(p)] # --nonmc --nontt
    else:
        command4 = [philosopher,"peptideprophet","--decoy","rev_","--enzyme",args.digest, "--decoyprobs","--ppm",
                    "--accmass","--nonparam","--expectscore","--database", str(dataFile),str(p)] # --nonmc --nontt
    print(" command : " +  " ".join(command4))
    subprocess.run(command4,cwd =batch_folder,check = True )
    

#PTMprophet prophet per file
#for p in pepXML_flie_list:
#    command10 = [philosopher,"ptmprophet",str(batch_folder / (f"interact-{p.stem}.pep.xml"))]
#    print(" command : " +  " ".join(command10))
#    subprocess.run(command10,cwd =batch_folder,check = True )

# protein prophet 
command5 =  [philosopher,"proteinprophet","--output","interact"] + [str(batch_folder / (f"interact-{p.stem}.pep.xml")) for p in pepXML_flie_list]
#    for base_name in base_file_name_list:
#        command6 = command6 + os.path.join(batch_folder, "interact-" + base_name + ".pep.xml " )
print(" command : " + " ".join( command5))
subprocess.run(command5,cwd =batch_folder,check = True)

# FDR
if non_specific:
    if args.mods == "global":
        command6 = [philosopher,"filter","--psm",str(args.fdr), "--pep", str(args.fdr),"--pepxml",str(batch_folder)]
    else:
        #command6 = [philosopher,"filter","--psm",str(args.fdr),"--pepxml",str(batch_folder) ]
        command6 = [philosopher,"filter","--psm",str(args.fdr), "--pep", str(args.fdr),"--pepxml",str(batch_folder),"--mods",args.mods]
else:
    if args.mods == "global":
        command6 = [philosopher,"filter","--psm",str(args.fdr), "--pep", str(args.fdr),"--pepxml",str(batch_folder)]
    else:
        #command6 = [philosopher,"filter","--razor","--sequential","--psm",str(args.fdr),"--prot","0.01","--mapmods","--pepxml",str(batch_folder),"--protxml","interact.prot.xml"]
        command6 = [philosopher,"filter","--psm",str(args.fdr), "--pep", str(args.fdr),"--pepxml",str(batch_folder),"--mods",args.mods]
print(" command : " +  " ".join(command6))
subprocess.run(command6,cwd =batch_folder,check = True)

# label free quantification 
command7 =[philosopher,"freequant","--dir","./"]
print(" command : " +  " ".join(command7))
subprocess.run(command7,cwd =batch_folder,check = True)

# report
command8 = [philosopher,"report"]
print(" command : " +  " ".join(command8))
subprocess.run(command8,cwd =batch_folder,check = True)

# clean  
print(" command : " +  " ".join(command1))
subprocess.run(command1 ,cwd =batch_folder,check= True)

"""
System info:
System OS: Windows Server 2008 R2, Architecture: AMD64
Java Info: 1.8.0_151, Java HotSpot(TM) 64-Bit Server VM, Oracle Corporation

Version info:
MSFragger-GUI version 6.0
MSFragger version 20180316
Philosopher version 20180322 (build 20180322.0841)

Will execute 9 commands:

C:\MS_program\philosopher_windows_amd64.exe workspace --clean 

C:\MS_program\philosopher_windows_amd64.exe workspace --init 

C:\MS_program\philosopher_windows_amd64.exe peptideprophet --decoy rev_ --decoyprobs --ppm --accmass --nonparam --expectscore --nontt --nonmc --database D:\reference_data_experiment\parallel-test\MSFragger\group_9\group_9_decoy.fasta D:\reference_data_experiment\parallel-test\MSFragger\group_9\QEP2_YMER_190_1_170516.pepXML 

C:\MS_program\philosopher_windows_amd64.exe peptideprophet --decoy rev_ --decoyprobs --ppm --accmass --nonparam --expectscore --nontt --nonmc --database D:\reference_data_experiment\parallel-test\MSFragger\group_9\group_9_decoy.fasta D:\reference_data_experiment\parallel-test\MSFragger\group_9\QEP2_YMER_191_1_170516.pepXML 

C:\MS_program\philosopher_windows_amd64.exe proteinprophet --output interact D:\reference_data_experiment\parallel-test\MSFragger\group_9\interact-*.pep.xml 

C:\MS_program\philosopher_windows_amd64.exe database --annotate D:\reference_data_experiment\parallel-test\MSFragger\group_9\group_9_decoy.fasta --prefix rev_ 

C:\MS_program\philosopher_windows_amd64.exe filter --sequential --tag rev_ --pepxml D:\reference_data_experiment\parallel-test\MSFragger\group_9 

C:\MS_program\philosopher_windows_amd64.exe report 

"""
