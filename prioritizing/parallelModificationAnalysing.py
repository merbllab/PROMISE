# -*- coding: utf-8 -*-
"""
Created on Wed May 16 14:04:36 2018

"""

import subprocess
import os
from time import clock, sleep
import re
import pandas as pd
from argparse import ArgumentParser
import random
import string
import glob
import logging
import time
tic = clock()
import pathlib


################ recive input from the user and detemine global parameters ###################
if __name__ == "__main__":
    parser = ArgumentParser()

    parser.add_argument('MSFragger_output', help="psm file format", metavar='string')
    parser.add_argument('output_file', help="output csv file",metavar='string')
    parser.add_argument('reference_data', help="reference data file",metavar='string')
    parser.add_argument('-mode', dest="mode", default = "new",choices=['new', 'append'],help="new run or append on existing run",metavar='string')
    parser.add_argument('-batch', dest="batch_size",default = 500, help="batch size", metavar=int)
    parser.add_argument('-filter', dest="filter_string", default = "",help="ignore mass shifts", metavar='string' )
    parser.add_argument('-queue', dest="queue",default = "new-short", help="queue for bsub commands",metavar='string')
    parser.add_argument('-annotation', dest="annotation",default = "general", help="determine dictionary for annotation seperate by ,",metavar='string')
    parser.add_argument('-mem', dest="mem_size", default = 8, help="size of java virtual machine in Gigabyte", metavar=int)


    args = parser.parse_args()


# =============================================================================+==============================================
# defined path - please update the following links to the external software and the temporary folders according to your need
# ============================================================================================================================
analyze_modification = pathlib.Path(r"/home/labs/yifatlab/assafk/PROMISE/prioritizing/analysedModification.py").resolve()

working_dir = os.path.dirname(args.output_file)
print("output file " + args.output_file)
print("working directory : " + working_dir)
prefix = ''.join(random.choice(string.ascii_uppercase) for _ in range(2))
n = int(args.batch_size)  #chunk row size
MSFraggerDf = pd.read_table(args.MSFragger_output,na_filter = False)
MSFraggerDf = MSFraggerDf[MSFraggerDf['Assigned Modifications']!=""]
if args.filter_string != "":
    re_mod = re.compile('''\w\(.+?\)''')
    newColumns = ['Keep']
    MSFraggerDf = MSFraggerDf.reindex(columns = MSFraggerDf.columns.tolist() + newColumns)
    filter_mod = set(args.filter_string.split(','))
    MSFraggerDf['Keep'] = MSFraggerDf.apply(lambda row: "+" if len(set(re_mod.findall(row['Assigned Modifications'])) - filter_mod) > 0  else "-",axis=1)
    MSFraggerDf = MSFraggerDf[MSFraggerDf['Keep']=="+"]


list_df = [MSFraggerDf[i:i+n] for i in range(0,MSFraggerDf.shape[0],n)]
for i in range(len(list_df)):
    list_df[i].to_csv(args.MSFragger_output[:-4] + "_group_" + str(i) + ".csv",index = False)
# =============================================================================
# creating logging 
# =============================================================================
log_dir = os.path.join(working_dir,"log")
if not os.path.exists(log_dir):
    os.makedirs(log_dir)
logging.basicConfig(filename = os.path.join(log_dir , "log.txt"), level=logging.INFO,
                    format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p') 
header = ( "\n" + 
    "      ===============================================================" + "\n" + 
    "      |   PROMISE  PROteine Modification Integrated Search Engine   |" + "\n" + 
    "      |                       Prioritizing phase                    |" + "\n" + 
    "      |                          Version  9                         |" + "\n" + 
    "      |                   " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "                 |" + "\n" + 
    "      ===============================================================" + "\n" + 
    "\n")
print(header)
logging.info(header) 
logging.info("\noutput file " + args.output_file + "\n" + 
             "working directory : " + working_dir)  
################# run anlaysing modification for all batches  ############################################

for i in range(len(list_df)):
    command =  "bsub -q " + args.queue + " -R \"rusage[mem=" + str(args.mem_size) + "000]\" -J Ana" + prefix +  str(i) + " -o LOG-Analysing" + str(i) + "-%J.txt python3 " + str(analyze_modification) + " . " + \
    args.MSFragger_output[:-4] + "_group_" + str(i) + ".csv " + args.output_file[:-4] + "_group_" + str(i) + ".csv " + args.reference_data + " " + " -mode " + \
    args.mode + " -annotation " + args.annotation
    print(command)
    subprocess.call(command, shell = True , cwd =  working_dir)
            
################# wait again till all jobs finished  #########################################
            
finished_tasks = False
while not finished_tasks: 
    out = subprocess.check_output("bjobs -J Ana" + prefix + "*",shell = True).decode("utf-8")
    #print(out)
    if out == '':
        finished_tasks = True
    else:
        print("number of active tasks with prefix " + prefix + " : ")
        subprocess.call ("bjobs -J Ana" + prefix + "* | wc -l" , shell = True , cwd =  working_dir)
        sleep(60)

print("finished tasks")

################ check if all tasks finished successfully or more memory is needed ############
success_count = 0
for i in range(len(list_df)):
    if os.path.isfile(args.output_file[:-4] + "_group_" + str(i) + "_valid_modification.csv"):
        success_count+=1
    else:
        for file in glob.glob("./LOG-Analysing" + str(i) + "-*.txt"):
            f= open(file, mode = 'r')
            log = f.read()
            f.close()
            if log.find("job killed after reaching LSF memory usage limit") != -1:
                print( "group " + str(i) + " fail because of memory limits")

print(str(success_count) + " succeed out of " + str(len(list_df)))


################ merge output files ###########################################################
outDf = []
for i in range(len(list_df)):
    tempDf =  pd.read_csv(args.output_file[:-4] + "_group_" + str(i) + "_valid_modification.csv" ,na_filter = False)
    outDf.append(tempDf)

finalDf = pd.concat(outDf, axis=0)
finalDf = finalDf[outDf[0].columns]
finalDf.to_csv(args.output_file ,index = False)
#subprocess.call ("cat *_spectrum_log.txt >> " + args.output_file[:-4] + "_spectrum_summary_log.txt"  , shell = True , cwd =  working_dir)


############### removing temp files ############################################################

subprocess.call ("rm -f *_group_*" , shell = True , cwd =  working_dir)
subprocess.call ("rm -f *_spectrum_log.txt" , shell = True , cwd =  working_dir)
subprocess.call ("rm -f LOG-Analysing*" , shell = True , cwd =  working_dir)

