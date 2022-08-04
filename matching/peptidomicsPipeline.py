"""
PROMISE is a distributed pipeline to detect post-translation modification on Mass Spectrometry data. The pipeline is solely for academic research,
 non-commercial or educational purposes; for other uses, please contact Weizmann institution, Merbl lab.

The pipeline is implemented to work on HPC environment with LSF API. The following LSF commands should be supported: bsub, bjobs, bpeek
The user should update all the hard-coded links in the code for external software used by the pipeline.
The pipeline was implemented and tested on python 3.6.0 
For pipeline description, please refer to the publication: “Post-translational modifications reshape the antigenic landscape of MHC I-immunopeptidome in tumors” 
 in nature biotechnology  


VERSION history:
    0.1 - 7/2018
    0.2 - 2/8/2018 - add support for defining FDR and if to split the data for each peptide length. this
                    is required if there are not alot of MS\MS raw files (poor data)
    0.3 - 26/8/2018 - add two automatic options: no_split - try to run the MSFragger without the split - it is
                                                prefered interm of FDR issues
                                                combineXML - combine  MSFragger XML output before FDR
                                                combineTSV - combine philosopher outputs after FDR
    1.0 28/9/2018 - recoding all the pipeline based on MSFragger split version
    1.1 - lock version with PhD propsal and add redo capabilities
    1.15 - add suppport to control JVM size from command line
    2.00 - implement split merge batches of 10 
            support of temporary files in no-backup folder
    2.01 - add support to define the LSF queue
    2.5 - change redo option to solve conflicts and job_id duplication
        - add prefix to job_id to avoid conflict with redo mode
    2.6 - fix bug in the mergPepXml file - in case of idnetical HyperScore the hits were chosen by the smallest
          delta mass. it was changed for the smallest absulote value of the tje delta mass 
    2.7 - update MSFragger version to 20190628 , add failures reason syntax apear in the new version 
    3.0 - update Philosopher version to v1.5.0. this change the PSM file format and required fix prioritizatino code 
    3.1 - add support for subgroup FDR
    3.2 - support of copying histogram and expectation files for future add on
        - add MSFragger task summary log 
        - integrated logging 
        - change redo process to stand alone status check before redoing the MSFragger tasks
        - change append mode to merge 
"""


import itertools
import pathlib
import re
import shlex
import subprocess
import typing
import numpy as np
import pandas as pd
import random
import string
import os
import time
from argparse import ArgumentParser
from Bio import SeqIO
import summerize_MSFragger_tasks as smt
import logging

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('working_directory', help="path to folder with mzML files", metavar='string')
    parser.add_argument('modification_list', help="modification list", metavar='file')
    parser.add_argument('reference_data', help="reference fasta file", metavar='file')
    parser.add_argument('-params', dest="params", default="", help="fragger params file name, will overide default params file", metavar='file')
    parser.add_argument('-split', dest="num_parts", default = 10, help="number of fasta fragments", metavar=int)
    parser.add_argument('-merge', dest="stat_file", default ="", help="path to stat file", metavar='file')
    parser.add_argument('-redo', dest="redo_stat_file", default ="", help="path to stat file", metavar='file')
    parser.add_argument('-digest', dest="digest", default = "nonspecific", choices=["Trypsin","nonspecific"], help="type of digestion", metavar='string')
    parser.add_argument('-init', dest="init_parts", default = 1, help="number of fasta fragments for the initial phase", metavar=int)
    parser.add_argument('-fdr', dest="fdr", default = 0.01, help="FDR for peptides", metavar=float)
    parser.add_argument('-JVM', dest="JVM_size", default = 20, help="size of java virtual machine in Gigabyte", metavar=int)
    parser.add_argument('-queue', dest="queue", default = "new-short", help="MS-Fragger LSF queue", metavar=int)

    args = parser.parse_args()
    
if 0:    
    print(args)

VERSION = 3.2

# =============================================================================+==============================================
# defined path - please update the following links to the external software and the temporary folders according to your need
# ============================================================================================================================
pipeline_path = pathlib.Path(r"/home/labs/yifatlab/assafk/PROMISE")
msfragger_jar_path = pathlib.Path(pipeline_path / r"external_software/MSFragger-20190628/MSFragger-20190628.jar").resolve()
contaminate_Data = pathlib.Path(pipeline_path / r"Database/contaminants_Ab_G_with_CON_prefix.fasta")
template_for_trptic_params = pathlib.Path(r"/home/labs/yifatlab/assafk/Pipeline/params_template/tryptic-fragger.params")
template_for_nonspecific_params = pathlib.Path(pipeline_path / r"params_template/non-specific-fragger.params")
mergePepXML = pathlib.Path(pipeline_path / r"matching/mergePepXML.py")
Philosopher = pathlib.Path(pipeline_path / r"matching/linuxPhilosopherCommands_sub_fdr.py")
work_dir = pathlib.Path(args.working_directory).resolve()
tempdir =pathlib.Path( re.compile('yifatlab').sub('yifatlab/nobackup', str(work_dir))) / pathlib.Path('split_peptide_index_tempdir') # this folder will contain all the temporary files



start_time = time.time()
do_msfragger = True  if args.stat_file == "" and args.redo_stat_file == "" else False
redo_msfragger = False if args.redo_stat_file == "" else True # will be used if only merge is required 

############### check for ilegal configuration #############################################
if args.stat_file != "" and args.redo_stat_file != "":
    exit("you cannot ask for redo and append at the same time")
############################################################################################
prefix = ''.join(random.choice(string.ascii_uppercase) for _ in range(2))
#prefix ="VH"
jvm_cmd = shlex.split(f"java -Xmx{str(int(args.JVM_size) - 3)}G -jar")
bsub_cmd = shlex.split("bsub -q "+ args.queue + " -R \"rusage[mem=" + str(args.JVM_size) + "000]\" -J F" + prefix + " -o LOG-MSFragger-" + prefix +"%J.txt")
reference_data_path = pathlib.Path(args.reference_data).resolve()
stat_track_file = "task_report"
# =============================================================================
# determine the source of the fragger params file. if the pipeline didn't recive
# an MSFragger params file as an input it will choose a default one base on the 
# enzyme digest method
# =============================================================================
if args.params != "": # overide the template params file
    param_path = pathlib.Path(args.params)
    print(args.params)
else:
    if args.digest == "Trypsin" :
        param_path = template_for_trptic_params
    elif args.digest == "nonspecific":
        param_path = template_for_nonspecific_params
    else:
        exit(f"couldnt find params file {args.params}")
params_txt = param_path.read_text()
############### creating a list of input mzML files ########################################
infiles = []
for file in os.listdir(work_dir):
    if file.endswith(".mzML"):
        infiles.append(pathlib.Path(file))
# =============================================================================
# initial general parameters 
# =============================================================================
jvm_cmd, param_path, infiles
bsub_msfragger_cmd = bsub_cmd + jvm_cmd + [msfragger_jar_path]
msfragger_cmd = jvm_cmd + [msfragger_jar_path]
#tempdir = work_dir / pathlib.Path('split_peptide_index_tempdir')
tempdir.mkdir(parents=True, exist_ok=True)
nobackup_work_dir = tempdir.parent
split_directory_status_file = tempdir / pathlib.Path(args.stat_file)
output_report_topN = int(re.compile(r'^output_report_topN *= *(\d+)', re.MULTILINE).search(params_txt).group(1))
output_max_expect_mo = re.compile(r'^output_max_expect *= *(\S+)', re.MULTILINE).search(params_txt)
search_enzyme_name = re.compile(r'^search_enzyme_name *= *(\S+)', re.MULTILINE).search(params_txt).group(1)
output_max_expect = 50.0 if output_max_expect_mo is None else float(output_max_expect_mo.group(1))
num_parts = int(args.num_parts)
tempdir_parts=[]
tempdir_trak_name =[]
initial_parts = int(args.init_parts)
max_task = 500 # defualt for max tasks forwarding to the cluster
max_split = 3 # each split is for 10 fragments - so input fasta can be split up to 1000 parts
free_job_id = 1000000

fragger_task_statistic = dict()
redo = dict()

########### read modification list file ######################3 
mod_dict = dict()
with pathlib.Path(args.modification_list).open('r') as f:
    for l in f:
        mod_dict[l.split(':')[0]] = l.split(':')[1].rstrip().split(',')
        #make sure no empty modification - cause MSFragger to fail
        for m in mod_dict[l.split(':')[0]]:
            if m =="":
                exit("found empty string in modification file")
        
mod_combination = list(itertools.combinations(mod_dict.keys(), 2))
mode_array = [mod_dict[mod_combination[i][0]] + mod_dict[mod_combination[i][1]] for i in range(len(mod_combination))]
mod_dirs = [tempdir / ''.join(mod_combination[i]) for i in range(len(mod_combination))]

# =============================================================================
# creating logging 
# =============================================================================
log_dir = work_dir / pathlib.Path('log')
log_dir.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename = log_dir / "log.txt", level=logging.INFO,
                    format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

header = ( "\n" + 
    "      ===============================================================" + "\n" + 
    "      |   PROMISE  PROteine Modification Integrated Search Engine   |" + "\n" + 
    "      |               Distributed MSFragger + Philosopher           |" + "\n" + 
    "      |                      Version " + str(VERSION) +"                            |" + "\n" + 
    "      |                   " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "                 |" + "\n" + 
    "      ===============================================================" + "\n" + 
    "\n")
print(header)
logging.info(header)

defualt_params_string = "\n"
defualt_params_string += f"working directory : {work_dir}\n"
defualt_params_string += f"temporary file's directory : {tempdir}\n"
defualt_params_string += f"refernce data : {reference_data_path}\n"
defualt_params_string += f"FDR : {str(args.fdr)}\n"
print(defualt_params_string)
logging.info(defualt_params_string)

with pathlib.Path(args.modification_list).open('r') as f:
    logging.info("\nmodifications: \n" + f.read() + "\n")

output_string = "\n"
if do_msfragger:
    output_string += "pipeline activated in regular mode: split MSFragger tasks and merge tasks\n"
    output_string += "pipeline will create the following modification folders:\n"
    for d in mod_dirs:
        output_string += f"{str(d.name)},"
elif not do_msfragger  and not redo_msfragger:
    output_string += f"pipline activated in merge mode: will split merge task accorindg to file\n {args.stat_file}\n"
#elif not do_msfragger and not split_merge and not redo_msfragger:
#    output_string = f"pipline activated in merge mode\n"
elif redo_msfragger:
    output_string += "pipeline activated in redo MSFragger mode: will check which MSFragger failed in the previouse run and redo them\n"
output_string += "\npipeline will work on the following mzML files:\n"
for f in infiles:
    output_string += f"{str(f)}\n"
print(output_string)
logging.info(output_string)



############# create reference fasta file with contaminate and decoys ##########################
def decoy(input_file, output_file):
    out = open(output_file,'w')
    fasta_sequences = SeqIO.parse(open(input_file),'fasta')
    for fasta in fasta_sequences:
        SeqIO.write(fasta,out,"fasta")
        fasta.id = "rev_" + fasta.id
        fasta.seq = fasta.seq[::-1]
        SeqIO.write(fasta,out,"fasta")
    out.close()


data_with_contaminate_file_name = work_dir / (reference_data_path.stem + "-with-contaminate.fasta")
data_with_contaminate_with_decoy_name = work_dir / (reference_data_path.stem + "-with-contaminate-with-decoy.fasta")
recomp_fasta = re.compile(r'^database_name\s*=\s*(.+?)$', re.MULTILINE)
edited_params_file = work_dir / "frager.params"

if do_msfragger:
    subprocess.run(["cp " + str(reference_data_path) +" "  + str(data_with_contaminate_file_name)], shell = True)
    subprocess.run(["cat " + str(contaminate_Data) + " >> " + str(data_with_contaminate_file_name)], shell = True)
    decoy(data_with_contaminate_file_name, data_with_contaminate_with_decoy_name)
    edited_params_file.write_text(recomp_fasta.sub(f'database_name = {data_with_contaminate_with_decoy_name.name}', params_txt))


#[fasta_path_str] = recomp_fasta.findall(params_txt)
fasta_path = data_with_contaminate_with_decoy_name # pathlib.Path(fasta_path_str)
param_path = edited_params_file

def get_free_jobID():
    global free_job_id
    free_job_id+=1
    return prefix + str(free_job_id)

def rename_lob_with_new_key(path,key,new_key):
    (path / f"LOG-MSFragger-{str(key)}.txt").rename(path/ f"LOG-MSFragger-{str(new_key)}.txt")

def get_fasta_path(temp_dir : pathlib.Path):
    global fasta_path
    for root, dirs, files in os.walk(temp_dir):
        if fasta_path.name in files:
            return( temp_dir / fasta_path.name)

def get_param_path(temp_dir : pathlib.Path):
    global param_path
    for root, dirs, files in os.walk(temp_dir):
        if param_path.name in files:
            return( temp_dir / param_path.name)

# =============================================================================
# split function will be called in the initial distributed tasks and dynamically 
# in case of MSFragger task failure
# =============================================================================
    
def split(tempdir,fasta_path,param_path,track_name,num_parts):
    global tempdir_parts, tempdir_trak_name #,fragger_task_statistic
    params_txt = param_path.read_text()
    local_tempdir_parts = [tempdir / str(i) for i in range(num_parts)]
    local_trak_name = [track_name + f"{i}" for i in range(num_parts)]
    fasta_prots: typing.List[bytes] = [e.rstrip() for e in fasta_path.read_bytes()[1:].split(b'\n>')]
    fasta_part_paths: typing.List[pathlib.Path] = [tempdir / str(i) / f'{fasta_path.name}' for i in range(num_parts)]
    param_part_paths: typing.List[pathlib.Path] = [tempdir / str(i) / param_path.name for i in range(num_parts)]
    infiles_name = [e.name for e in infiles]
    infiles_symlinks_target_pairs = [(ee / e.name, e) for e in infiles for ee in local_tempdir_parts]
    cmds = [bsub_msfragger_cmd + [param_part_path.name, *infiles_name, '--partial', f'{i}']
    		for i, param_part_path in zip(range(num_parts), param_part_paths)]
    for i in range(num_parts):
        cmds[i][6] += local_trak_name[i] 



    def set_up_directories():
    	tempdir.mkdir(exist_ok=True)
    	for e in local_tempdir_parts:
    		e.mkdir()
    	for ln, target in infiles_symlinks_target_pairs:
    		ln.symlink_to(target.resolve())
    	for fasta_part, fasta_part_path in zip(np.array_split(np.array(fasta_prots, object), num_parts), fasta_part_paths):
    		with pathlib.Path(fasta_part_path).open('wb') as f:
    			f.writelines(b'>' + e + b'\n' for e in fasta_part)
    
    	for param_part_path, fasta_name in zip(param_part_paths, fasta_part_paths):
    		param_part_path.write_text(
    			recomp_fasta.sub(f'database_name = {fasta_name.name}', params_txt))
    
        
    set_up_directories()
    if tempdir in tempdir_parts: 
        tempdir_parts.remove(tempdir)
        tempdir_trak_name.remove(track_name)
    tempdir_parts += local_tempdir_parts 
    tempdir_trak_name += local_trak_name
    return(cmds,local_tempdir_parts)

# =============================================================================
# execute MSFragger    
# =============================================================================
def run_msfragger(cmds,tempdir_parts_batch):
    global fragger_task_statistic
    for cmd, cwd in zip(cmds, tempdir_parts_batch):
        launch = subprocess.check_output(cmd, cwd=cwd).decode("utf-8")
        job_id = prefix + launch.split("<")[1].split(">")[0]
        if job_id not in fragger_task_statistic.keys():
            fragger_task_statistic[job_id] = {"path" : cwd ,"NAME" : cmd[6][3:],"STAT" : "UN", "RESULT" : "UN","fail_type" : "UN","peptides" : "0" ,
                                  "modified_peptides" : "0","mem_size":"0" , "run_time":"UN"} 
        else:
            new_key = get_free_jobID()
            rename_lob_with_new_key(fragger_task_statistic[job_id]['path'],job_id,new_key) # we assume this cannot be done on a task waiting in LSF queue but old run
            fragger_task_statistic[new_key] = fragger_task_statistic.pop(job_id)
            fragger_task_statistic[job_id] = {"path" : cwd ,"NAME" : cmd[6][3:],"STAT" : "UN", "RESULT" : "UN","fail_type" : "UN","peptides" : "0" ,
                                  "modified_peptides" : "0","mem_size":"0" , "run_time":"UN"} 
           
                
        print(f"{launch}")

    
####################### create initial directories structure ####################################
tempdir.mkdir(exist_ok=True)
#print(tempdir_parts)
#for d in mod_dirs:
#    d.mkdir()
#    split(d,fasta_path,param_path,"",initial_parts)

# =============================================================================
#  edit fragger.params acording to the moficication 
# =============================================================================

if do_msfragger:
    for mod_number ,p in enumerate(mod_dirs):
        try:
            p.mkdir()
        except FileExistsError:
            print(" temp directory should be clear before running new pipline analysis")
            logging.error(" temp directory should be clear before running new pipline analysis")
            exit()
        cmds, tempdir_parts_batch = split(p,fasta_path,edited_params_file,p.name,initial_parts)
        params = [p / str(i) / edited_params_file.name for i in range(initial_parts)]
        for param in params:
            content = param.read_text()
            new_content = "I change it\n"
            mod_count = 0
            for line in content.split('\n'):
                new_line = line
                if line.find("variable_mod_0") != -1:
                    if mod_count in range(len(mode_array[mod_number])):
                        new_line = f"variable_mod_0{str(mod_count+1)} = {mode_array[mod_number][mod_count]}"
                        #print(new_line,mod_number,mod_count,param,p)
                    else:
                        new_line = f"# variable_mod_0{str(mod_count+1)} ="
                    mod_count+=1
                new_content += f"{new_line}\n"
            param.write_text (new_content)            
        run_msfragger(cmds, tempdir_parts_batch)
        #print(cmds,tempdir_parts_batch)
    
# =============================================================================
# the following utility function will gather information on the MSFragger runs
# and will verify the failure reason or success of the task
# =============================================================================

"""
    there are 3 type of MSFragger taks failures:
        1. java memory heap exception where the bsub task is stuck and should be teminate by invoke the bkill command
        2. java memory heap exception wherer the bsub finished with status EXIT
        3. MSFragger stop working because it created too many peptides , bsub command finished with status DONE
        4. java memory problem: GC overhead limit exceeded , bsub command in RUN status
        5. java: command not found
        6. too fragmented - no peptide created 
    
    in any of those failures we split the input fasta file into batches and re-run MSFragger on each batch
"""  



java_error_msg = "java.lang.OutOfMemoryError: Java heap space"
#too_many_peptide_error_msg = "We generated too many peptides"
java_error_msg2 = "java.lang.OutOfMemoryError: GC overhead limit exceeded"
java_command = "java: command not found"
no_peptide = "Reduced to 0  peptides"
#LSF_time_limit = "job killed after reaching LSF run time limit"
#LSF_memory_limit = "job killed after reaching LSF memory usage limit"

split_error_mesage_dict = {
        "Too many peptides were generated" : "MSFragger too many peptides",
        "java.lang.OutOfMemoryError: Java heap space" : "java heap log",
        "job killed after reaching LSF run time limit" : "LSF time limit",
        "job killed after reaching LSF memory usage limit" : "LSF memory limit",
        "There is insufficient memory for the Java Runtime Environment to continue" : "ENV java memory",
        "java.util.concurrent.TimeoutException" : "ENV java timeout",
        "Bus error" : "ENV Bus error",
        "job killed by owner" : "Killed by owner"}

job_cmd = ["bjobs","-J", f"F{prefix}*"]

def retrieve_peptide_count(content):
    pep = "0"
    modified_pep ="0"
    if content.find("Reduced to") != -1:
        pep = content.split("Reduced to")[1].strip().split(" ")[0]
    if content.find("Generated") != -1:
        modified_pep = content.split("Generated")[1].strip().split(" ")[0]
    elif content.find("Too many peptides were generated")!=-1:
        modified_pep = content.split("Too many peptides were generated")[1].strip().split("!")[0]
    return pep,modified_pep

def get_run_time(content):
    re_time = re.compile('''Run time :\s+?(\d+?)\ssec.''', re.MULTILINE)
    x = re_time.search(content)
    if x is None:
        return "UN"
    else:
        return x.group(1)

def get_executing_nod(content):
    re_host = re.compile('''Job was executed on host\(s\) <(.+?)>''',re.MULTILINE)
    x = re_host.search(content)
    if x is None:
        return "UN"
    else:
        return x.group(1)
    
                
def peeking(peek,jm):
    global fragger_task_statistic
    do_kill_and_split = False
    if peek.find(java_error_msg) != -1: # error type 1 (see discription above) - bsub taks has to be terminate
        print(f"{jm} stuck with java heap space exception")
        fragger_task_statistic[jm]["RESULT"] = "Fail"
        fragger_task_statistic[jm]["fail_type"] = "java heap"
        do_kill_and_split = True
        fragger_task_statistic[jm]["peptides"], fragger_task_statistic[jm]["modified_peptides"] = retrieve_peptide_count(peek)
    elif peek.find(java_error_msg2) != -1:
        print(f"{jm} stuck with java GC overhead limit exceeded")
        fragger_task_statistic[jm]["RESULT"] = "Fail"
        fragger_task_statistic[jm]["fail_type"] = "GC overhead"
        do_kill_and_split = True
        fragger_task_statistic[jm]["peptides"], fragger_task_statistic[jm]["modified_peptides"] = retrieve_peptide_count(peek)
    return do_kill_and_split

def check_log(log,jm):
    global fragger_task_statistic
    do_split = False
    known_result_reason = False
    if log.find("The output (if any) is above this job summary") != -1:# check if the log file is a complete file or it is partial - maybe task is still onoing
        # reasons that will stop the MSfragger recursive spliting 
        if log.find("Successfully completed") != -1: # success
            if set([(fragger_task_statistic[jm]["path"]  / (e.stem + '_scores_histogram.tsv')).exists() for e in infiles]):
                fragger_task_statistic[jm]["RESULT"] = "Success"
                subprocess.call( ["rm *.pepindex"], cwd = fragger_task_statistic[jm]["path"],shell = True)
            fragger_task_statistic[jm]["peptides"], fragger_task_statistic[jm]["modified_peptides"] = retrieve_peptide_count(log)
            fragger_task_statistic[jm]["run_time"] = get_run_time(log)
            fragger_task_statistic[jm]["host"] = get_executing_nod(log)
        elif log.find(java_command) != -1: # will casuse contiues failers without any stop - script should re-start 
            print(f"{fragger_task_statistic[jm]['NAME']} coudn't find java command")
            exit(1)
        elif log.find(no_peptide) != -1: # need to end without output
            print(f"{fragger_task_statistic[jm]['NAME']} too fragmented - no peptides")
            fragger_task_statistic[jm]["RESULT"] = "END"
            known_result_reason = True
            
        # reasons that will need to split fasta and re-run MSFragger 
        elif not known_result_reason: 
            for key in split_error_mesage_dict.keys():
                if log.find(key) != -1:
                    print(f"{fragger_task_statistic[jm]['NAME']} {split_error_mesage_dict[key]}")
                    fragger_task_statistic[jm]["fail_type"] = split_error_mesage_dict[key]
                    known_result_reason = True
            if not known_result_reason: # coudn't find a known reason 
                fragger_task_statistic[jm]["fail_type"] = "un known"                
            if len(fragger_task_statistic[jm]["NAME"]) < (max_split + 2): # 2 for modification letter 
                fragger_task_statistic[jm]["RESULT"] = "Fail"
                do_split = True
            else: # we split already 4 times , no more spliting 
                fragger_task_statistic[jm]["RESULT"] = "END"
            fragger_task_statistic[jm]["peptides"], fragger_task_statistic[jm]["modified_peptides"] = retrieve_peptide_count(log)
            fragger_task_statistic[jm]["host"] = get_executing_nod(log)
            fragger_task_statistic[jm]["run_time"] = get_run_time(log)
    else:
        logging.warning(f" reading partial log file {fragger_task_statistic[jm]['NAME']} , ignoring" )
    return do_split

def get_dir_size(p):
    total = 0
    for f in os.listdir(p):
        if  os.path.isfile(p / f) and not (p /f).is_symlink():
            total += os.path.getsize(p / f)
    total = total >> 20 
    return str(total)


# =============================================================================
# special function if the pipeline start in redo mode for recovery 
# =============================================================================

redo_initiate = False
if redo_msfragger:
    # create a stand alone task report   
    #tempdir =pathlib.Path( re.compile('assafk').sub('nobackup/assafk', str(work_dir))) / pathlib.Path('split_peptide_index_tempdir')
    logging.info("\nanalyzing previous pipeline run and created a recovered tasks list")
    path_list = [pathlib.Path(x[0]) for x in os.walk(tempdir)]
    for path in path_list:
        log_files =[]
        log_times = []
        key = ""
        name = "".join(str(path).split("split_peptide_index_tempdir")[1].split("/"))    
    
        for file in os.listdir(path):
            if file.startswith("LOG-MSFragger"):
                log_files.append(path / file)
        if len(log_files)>1:
            print(f"info: found more than one log at {name} , {[str(x.stem) for x in log_files]}")
            logging.debug(f"found more than one log at {name} , {[str(x.stem) for x in log_files]}")
            for log_file in log_files:
                log_times.append(os.path.getmtime(log_file))
            sorted_log_file = [x for _,x in sorted(zip(log_times,log_files),reverse = True)]  
            #print(sorted_log_file)
            key = str(sorted_log_file[0].stem).split("-")[2]
            #print(str(sorted_log_file[0].stem),key,name)
        elif len(log_files)==1:
            key = str(log_files[0].stem).split("-")[2]
            #print(str(log_files[0].stem),key,name)
        else:
            print(f"info: no log file at path {path.stem}")
            logging.debug(f"no log file at path {path.stem}")
        if key =="": # no log found , can be parent folder - check for params file for verification 
            no_params_file = True
            for file in os.listdir(path):
                if file.endswith("params"):
                    new_key = str(get_free_jobID())
                    fragger_task_statistic[new_key] = {"path" : path ,"NAME" : name,"STAT" : "UN", "RESULT" : "REDO","fail_type" : "UN","peptides" : "0" ,
                                          "modified_peptides" : "0","mem_size":"0" , "run_time":"UN"} 
                    print(f"warning: found ready to work path without log file in {name} , create new key {new_key}")
                    logging.warning(f"found ready to work path without log file in {name} , create new key {new_key}")
                    no_params_file = False
            if  no_params_file:
                print(f"info: found parent path {name} , will not enter into jobs list")
                logging.debug(f"found parent path {name} , will not enter into jobs list")
                    
        elif key not in fragger_task_statistic.keys():
            fragger_task_statistic[key] = {"path" : path ,"NAME" : name,"STAT" : "UN", "RESULT" : "UN","fail_type" : "UN","peptides" : "0" ,
                                  "modified_peptides" : "0","mem_size":"0" , "run_time":"UN"} 
        else: # duplicate key
            new_key = key[0:2] + str(get_free_jobID())
            rename_lob_with_new_key(path,key,new_key)
            fragger_task_statistic[new_key] = {"path" : path ,"NAME" : name,"STAT" : "UN", "RESULT" : "UN","fail_type" : "UN","peptides" : "0" ,
                                  "modified_peptides" : "0","mem_size":"0" , "run_time":"UN"} 
            print(f"warning: find duplicated job id between {name} and {fragger_task_statistic[key]['NAME']} alocate new key {new_key}")
            logging.warning(f"find duplicated job id between {name} and {fragger_task_statistic[key]['NAME']} alocate new key {new_key}")
                
    for jm in fragger_task_statistic.keys():
        try:
            log = (fragger_task_statistic[jm]["path"] / f"LOG-MSFragger-{jm}.txt").read_text()
            if check_log(log,jm): # found fail task - check if it already been split to jobs 
                sub_jobs_path_list = [pathlib.Path(x[0]) for x in os.walk(fragger_task_statistic[jm]["path"])]
                if len(sub_jobs_path_list) == 1: # split was not done - keep it to the redo script to find it - keep STAT as UN
                    if fragger_task_statistic[jm]["RESULT"] != "END": # make sure there wasn't any action because no more spliting aloud 
                        fragger_task_statistic[jm]["RESULT"] = "UN"
                    print(f"info: found FAILED task but no action taken in {fragger_task_statistic[jm]['NAME']} , keep its result as UN ")
                    logging.debug(f"found FAILED task but no action taken in {fragger_task_statistic[jm]['NAME']} , keep its result as UN ")
        except FileNotFoundError:
            if fragger_task_statistic[jm]["RESULT"] != "REDO": # make sure it is already known log doesn't exist
                exit(f"log doesn't exist {fragger_task_statistic[jm]['NAME']}")
     
        
    tasks_df = pd.DataFrame.from_dict(fragger_task_statistic).T  
    summary_string = (f"\nnumber of Success jobs is {len(tasks_df[tasks_df['RESULT']=='Success'])}" + "\n" +
                                                     f"number of Fail jobs is {len(tasks_df[tasks_df['RESULT']=='Fail'])}" + "\n" +
                                                     f"number of UN jobs that will be split and run is {len(tasks_df[tasks_df['RESULT']=='UN'])}" + "\n" + 
                                                     f"number of UN jobs that ended without result {len(tasks_df[tasks_df['RESULT']=='END'])}" + "\n" + 
                                                     f"number of REDO jobs that didn't start is  {len(tasks_df[tasks_df['RESULT']=='REDO'])}" + "\n" +
                                                     f"WE ESTIMATE  {len(tasks_df[tasks_df['RESULT']=='UN'])*10} JOBS ARE YET TO BE PERFORMED")               
    print(tasks_df[["NAME" ,"STAT" , "RESULT" ,"fail_type" ,"peptides" , "modified_peptides" ,"mem_size"]])
    print(summary_string) 
    logging.info(summary_string)
    
    tasks_df.to_csv(tempdir / (f"recoverd_redo_tasks.csv"))
    #statDf = pd.read_csv(args.redo_stat_file)
#    statDf = statDf.set_index('Unnamed: 0')
    #fragger_task_statistic = statDf.T.to_dict()
    infiles_name = [e.name for e in infiles]
    cmds = []
    tempdir_parts_batch = []
    for key in fragger_task_statistic.keys():
        fragger_task_statistic[key]['path'] = pathlib.Path(fragger_task_statistic[key]['path'])
        if fragger_task_statistic[key]['RESULT'] == "REDO":
            cmd = bsub_msfragger_cmd + [get_param_path(fragger_task_statistic[key]['path']).name, *infiles_name, '--partial', f"{str(fragger_task_statistic[key]['path'])[-1:]}"]
            cmd[6] = f"F{prefix}{fragger_task_statistic[key]['NAME']}"
            cmds.append(cmd)
            tempdir_parts_batch.append(fragger_task_statistic[key]['path'])
    run_msfragger(cmds, tempdir_parts_batch)
    do_msfragger = True
    if len(cmds) > 0:
        redo_initiate = True
        logging.info(f"found {str(len(cmds))} MSFragger tasks with REDO result status - redo them" )
    else:
        logging.info(f"couldn't find REDO tasks, will start directly at the UN tasks" )
        
# =============================================================================
# busy wating till all the MSFragger tasks will finsished. 
# =============================================================================

finished_MSFragger = False if do_msfragger else True
itirate = 0
snap = 1107
actual_tasks = 0
while not finished_MSFragger:
    out=(subprocess.check_output(job_cmd)).decode("utf-8")
    print(out)
    # no more bjobs - can contiue to philosopher tasks
    if out == '':
        actual_task = 0
        finished_MSFragger = True
        snapshot =  fragger_task_statistic.copy()
        for jm in snapshot.keys():
            if actual_task < max_task:
                if snapshot[jm]["RESULT"] =="UN" : # and snapshot[jm]["STAT"] !="PEND"
                    try: 
                        log = (snapshot[jm]["path"] / f"LOG-MSFragger-{jm}.txt").read_text()
                        if check_log(log,jm): # although the MSFragger jobs finished - some of the failures didn't split
                            if not redo_initiate: # will split and start handling the new jobs only at the begining of redo mode 
                                p = snapshot[jm]["path"]
                                cmds, tempdir_parts_batch = split(p,get_fasta_path(p),get_param_path(p),snapshot[jm]["NAME"],num_parts)
                                run_msfragger(cmds, tempdir_parts_batch)
                                finished_MSFragger = False
                                actual_task+=num_parts
                                redo_initiate = True
                                break
                            else:
                                print(f"although task {snapshot[jm]['NAME']}  fail , will not aplit and save for redo" )
                                redo[jm] = snapshot[jm]        
                    except FileNotFoundError:
                        print(f"log doesn't exist {snapshot[jm]['NAME']} and status is {snapshot[jm]['STAT']}, might need to redo" )
                        redo[jm] = snapshot[jm] 
                        fragger_task_statistic[jm]["RESULT"] = "REDO"
        if finished_MSFragger:
            print("finished MSFragger tasks")
    else:
        # update status 
        print("current statistics: ")
        actual_task = len(out.splitlines()[1:])
        for line in out.splitlines()[1:]: # skeep headers
            vals = line.split() # good for only first 3 fields 
            job_id = prefix + vals[0].strip()
            fragger_task_statistic[job_id]["STAT"] = vals[2].strip()
        snapshot =  fragger_task_statistic.copy()
        for jm in snapshot.keys():
            if actual_task < max_task:
                # make decisions real time peeking
                if snapshot[jm]["RESULT"] =="UN" : 
                    try:
                        peek =subprocess.check_output(["bpeek",str(jm)[2:]]).decode("utf-8")
                        if peeking(peek,jm):
                            subprocess.call(["bkill" ,jm[2:]])
                            p = snapshot[jm]["path"]
                            cmds, tempdir_parts_batch = split(p,get_fasta_path(p),get_param_path(p),snapshot[jm]["NAME"],num_parts)
                            run_msfragger(cmds, tempdir_parts_batch)
                            actual_task+=num_parts
                    except subprocess.CalledProcessError:
                        print(f"fail bpeek command {snapshot[jm]['NAME']} try log file")
                        try: 
                            log = (snapshot[jm]["path"] / f"LOG-MSFragger-{jm}.txt").read_text()
                            if check_log(log,jm):
                                p = snapshot[jm]["path"]
                                cmds, tempdir_parts_batch = split(p,get_fasta_path(p),get_param_path(p),snapshot[jm]["NAME"],num_parts)
                                run_msfragger(cmds, tempdir_parts_batch)
                                actual_task+=num_parts
                        except FileNotFoundError:
                            print(f"log doesn't exist {snapshot[jm]['NAME']}")
                            if snapshot[jm]["STAT"] !="PEND":
                                print(f"log doesn't exist {snapshot[jm]['NAME']} and status is {snapshot[jm]['STAT']}, might need to redo" )
                                redo[jm] = snapshot[jm]
                                
        snap_df = pd.DataFrame.from_dict(fragger_task_statistic).T                 
        print(snap_df[["NAME" ,"STAT" , "RESULT" ,"fail_type" ,"peptides" , "modified_peptides" ,"mem_size"]])

        if len(redo)>0:
            print("redo:")
            print(pd.DataFrame.from_dict(redo).T)
        if itirate % 10 ==0: # once every 10 min free disk space from success task and calculate disk space
            for jm in snapshot.keys():
                if snapshot[jm]["RESULT"] =="UN" :
                    try: 
                        log = (snapshot[jm]["path"] / f"LOG-MSFragger-{jm}.txt").read_text()
                        if log.find("Successfully completed") != -1:
                            if set([(snapshot[jm]["path"]  / (e.stem + '_scores_histogram.tsv')).exists() for e in infiles]):
                                fragger_task_statistic[jm]["RESULT"] = "Success"
                                subprocess.call( ["rm *.pepindex"], cwd = fragger_task_statistic[jm]["path"],shell = True)
                                fragger_task_statistic[jm]["peptides"], fragger_task_statistic[jm]["modified_peptides"] = retrieve_peptide_count(log)
                                fragger_task_statistic[jm]["run_time"] = get_run_time(log)
                    except FileNotFoundError:
                        print(f"log doesn't exist {snapshot[jm]['NAME']}")
                try:
                    fragger_task_statistic[jm]["mem_size"]= get_dir_size(fragger_task_statistic[jm]["path"])
                except FileNotFoundError: 
                    print(f"{fragger_task_statistic[jm]['NAME']} mem file wasn't updated")
                
            (pd.DataFrame.from_dict(fragger_task_statistic).T).to_csv(tempdir / (f"{stat_track_file}_{snap}.csv"))
            if len(redo)>0: 
                (pd.DataFrame.from_dict(redo).T).to_csv(tempdir / (f"redo_{stat_track_file}_{snap}.csv"))
            snap+=1
            logging.info((f"success tasks: {snap_df[snap_df.RESULT =='Success'].count()['RESULT']} , " +
                                 f"Fail tasks: {snap_df[snap_df.RESULT =='Fail'].count()['RESULT']} , " + 
                                 f"UN tasks: {snap_df[snap_df.RESULT =='UN'].count()['RESULT']} , " +
                                 f"{actual_task} tasks alive , " + 
                                 f"running time:  {round((time.time() - start_time),2)} seconds"))
        itirate+=1
        print(f"running time:  {(time.time() - start_time)} seconds ---" )
        print(f"working directory: {str(work_dir)}")
        print(f"total of success tasks: {snap_df[snap_df.RESULT =='Success'].count()['RESULT']} ")
        print(f"total of Fail tasks: {snap_df[snap_df.RESULT =='Fail'].count()['RESULT']} ")
        print(f"total of UN tasks: {snap_df[snap_df.RESULT =='UN'].count()['RESULT']} ")
        print(f"not finisihed, {actual_task} tasks alive")
        time.sleep(120)
    

print(f"start split merge tasks")

# =============================================================================
# create success / failure status report
# =============================================================================

if do_msfragger:
    statDf = pd.DataFrame.from_dict(fragger_task_statistic).T
    redoDf = pd.DataFrame.from_dict(redo).T
    statDf.to_csv(tempdir / (f"{stat_track_file}.csv"))
    split_directory_status_file = tempdir / (f"{stat_track_file}.csv")

    redoDf.to_csv(tempdir / (f"redo_{stat_track_file}.csv"))
    successful_tempdir_parts = statDf[statDf["RESULT"]=="Success"].path.tolist()
    diff = set(tempdir_parts) - (set(tempdir_parts) & set(successful_tempdir_parts))
    if len(diff)>0:
        print("the following foders fail to complete :\n ")
        for d in diff:
            print(str(d))
        print("\ncontinue with the one who succeed")
    

task_df = pd.read_csv(split_directory_status_file)
successful_df = task_df[task_df["RESULT"]=="Success"]

# =============================================================================
# start merge process by levels 
# =============================================================================

level_dict = dict()
level_one = set([name[0:2] for name in successful_df.NAME.tolist()])

print("creating level 3 merge tasks")
print(task_df)
print(successful_df)
print(level_one)
level_three_tasks = dict()
for one in level_one:
    for two in list(map(str,range(10))):
        for three in list(map(str,range(10))):
            subset_df = successful_df[successful_df['NAME'].str.startswith(one+two+three)]
            if len(subset_df)>1:
                temp_dir = pathlib.Path(subset_df.iloc[0]['path']).parent
                #print(subset_df, temp_dir)
                subset_df.to_csv(temp_dir / ("merge_report.csv"))
                level_three_tasks[one+two+three] = {'path':temp_dir,'len':len(subset_df),'STAT' : "UN"}
                successful_df.drop(subset_df.index,inplace=True)
                successful_df = successful_df.append({'NAME': one+two+three , 'RESULT' : "Merged","STAT":"","path":temp_dir}, ignore_index=True)
                
print("creating level 2 merge tasks")          
level_two_tasks = dict()
for one in level_one:
    for two in list(map(str,range(10))):
        subset_df = successful_df[successful_df['NAME'].str.startswith(one+two)]
        if len(subset_df)>1:
            temp_dir = pathlib.Path(subset_df.iloc[0]['path']).parent
            #print(subset_df, temp_dir)
            subset_df.to_csv(temp_dir / ("merge_report.csv"))
            level_two_tasks[one+two] = {'path':temp_dir, 'len':len(subset_df),'STAT' : "UN"}
            successful_df.drop(subset_df.index,inplace=True)
            successful_df = successful_df.append({'NAME': one+two , 'RESULT' : "Merged","STAT":"","path":temp_dir}, ignore_index=True)
 
print("creating level 1 merge tasks")          
level_one_tasks = dict()
for one in level_one:
    subset_df = successful_df[successful_df['NAME'].str.startswith(one)]
    if len(subset_df)>1:
        temp_dir = pathlib.Path(subset_df.iloc[0]['path']).parent
        #print(subset_df, temp_dir)
        subset_df.to_csv(temp_dir / ("merge_report.csv"))
        level_one_tasks[one] = {'path':temp_dir, 'len':len(subset_df),'STAT' : "UN"}
        successful_df.drop(subset_df.index,inplace=True)
        successful_df = successful_df.append({'NAME': one , 'RESULT' : "Merged","STAT":"","path":temp_dir}, ignore_index=True)

print("creating top level merge tasks in batches")
top_level_batches = dict()
top_level_final = dict()      
if len(successful_df)>10:
    cp = 0
    batch = 0
    while cp < len(successful_df):
        batch+=1
        ep = min(cp+10,len(successful_df))
        if ep == len(successful_df) -1: # not to leave just one task to merge
            ep = len(successful_df)
        subset_df = successful_df.iloc[list(range(cp,ep)),:]
        temp_dir = nobackup_work_dir / (f"split_peptide_index_tempdir/batch{batch}")
        temp_dir.mkdir(exist_ok=True)
        subset_df.to_csv(temp_dir / ("merge_report.csv"))
        top_level_batches[batch]= {'path':temp_dir, 'len':len(subset_df),'STAT' : "UN"}
        top_level_final[batch] = {"NAME": f"batch{str(batch)}"  , "path" : temp_dir,"RESULT" : "Merged","STAT" : ""}
        cp = ep
    

         
(pd.DataFrame.from_dict(level_one_tasks).T).to_csv(nobackup_work_dir / (f"level1_merge_tasks.csv"))        
(pd.DataFrame.from_dict(level_two_tasks).T).to_csv(nobackup_work_dir / (f"level2_merge_tasks.csv"))  
(pd.DataFrame.from_dict(level_three_tasks).T).to_csv(nobackup_work_dir / (f"level3_merge_tasks.csv")) 
(pd.DataFrame.from_dict(top_level_batches).T).to_csv(nobackup_work_dir / (f"top_level_merge_tasks.csv")) 
(pd.DataFrame.from_dict(top_level_final).T).to_csv(nobackup_work_dir / (f"merge_report.csv")) 



#successful_df.to_csv(nobackup_work_dir / ("split_peptide_index_tempdir/merge_report.csv")) 

def check_log(level,merge_stat):
    for task in merge_stat.keys():
        log = (merge_stat[task]["path"] / f"LOG-Merge-{task}-{merge_stat[task]['raw'].stem[-4:]}.txt").read_text()
        if log.find("Successfully completed") != -1: # success
            merge_stat[task]['STAT'] = "Success"
    result_df = pd.DataFrame.from_dict(merge_stat).T
    if len(result_df)>=1: # empty level
        if set(result_df['STAT'].tolist()) == {"Success"}:
            print(f"level {level} merge tasks succeed")
            result_df.to_csv(nobackup_work_dir / (f"level{level}_merge_result-Success.csv"))
        else:
            print(f"******************\nlevel {level} merge tasks has failures\n********************")
            logging.error(f"\nlevel {level} merge tasks has failures\n")
            result_df.to_csv(nobackup_work_dir / (f"level{level}_merge_result-Partial.csv"))
        
def merge_by_level (level,task_dict): 
    
    merge_task_statistic = dict()
    print(f"executing merge commands for level {level}")
    command_merge = []
    tempdir_list = []
    for task in task_dict.keys():
        tempdir = task_dict[task]['path']    
        command_merge = command_merge +  [["bsub","-q",args.queue,"-J",f"M{prefix}{task}{e.stem[-3:]}","-o", str(tempdir / f"LOG-Merge-%J-{e.stem[-4:]}.txt"),"python3",
                         mergePepXML ,str(tempdir / "merge_report.csv"),e,tempdir,tempdir,
                         "-output_report_topN",str(output_report_topN),"-output_max_expect",str(output_max_expect)] for e in infiles]
        tempdir_list = tempdir_list  + ([tempdir] * len(infiles))
    
    finished_merge = False
    run_commands = True
    cp = 0
    while not finished_merge:
        print(f"******* merge level {level} with prefix {prefix} *******")
        out=(subprocess.check_output(["bjobs","-J",f"M{prefix}*"])).decode("utf-8")
    
        if out == '' and cp < len(command_merge):
            run_commands = True
        elif out == '' and cp == len(command_merge):
            run_commands = False
            finished_merge = True
        elif len(out.splitlines()[1:]) < max_task:
            run_commands = True
        else:
            run_commands = False
            print(f"alive merge tasks = {len(out.splitlines()[1:])}")
        if run_commands:
            print(f"alive merge tasks = {len(out.splitlines()[1:])}")
            start = cp
            end = min(cp + max_task,len(command_merge))
            for cmd, cwd in zip(command_merge[start:end], tempdir_list[start:end]):
                #print(cmd,cwd)
                launch = subprocess.check_output(cmd, cwd=cwd).decode("utf-8")
                job_id = launch.split("<")[1].split(">")[0]
                merge_task_statistic[job_id] = {"path" : cwd ,"NAME" : cmd[4][3:],"STAT" : "UN", "RESULT" : "UN","fail_type" : "UN","level" : level, "raw" :  cmd[10]} 
                print (str(cmd[10])[-10:],str(cwd)[-5:],cp) 
            if cp != end:
                print(f"batch start = {start} end = {end}")
            cp= end
        time.sleep(60)    
    print(f"finished level {level} merge tasks")
    logging.info(f"finished level {level} merge tasks\n")
    return merge_task_statistic

level3_merge_stat = merge_by_level("3",level_three_tasks )
check_log("3",level3_merge_stat)

level2_merge_stat = merge_by_level("2",level_two_tasks )
check_log("2",level2_merge_stat)

level1_merge_stat = merge_by_level("1",level_one_tasks )
check_log("1",level1_merge_stat)

if len (top_level_batches) > 0:
    top_level_batches_stat = merge_by_level("top_batch",top_level_batches)
    check_log("top_batch",top_level_batches_stat)
    
    print("executing merge commands for final merge")
    tempdir = nobackup_work_dir / ("split_peptide_index_tempdir")
    command_merge_upper_level = [["bsub","-q",args.queue,"-J",f"M{prefix}{e.stem[-3:]}","-o", str(nobackup_work_dir / f"LOG-Merge-%J-{e.stem[-4:]}.txt"),"python3",
                     mergePepXML,str(nobackup_work_dir / "merge_report.csv"),e,nobackup_work_dir,work_dir,
                     "-output_report_topN",str(output_report_topN),"-output_max_expect",str(output_max_expect)] for e in infiles]
    root_merge_task_statistic = dict()
    for cmd in command_merge_upper_level:
        launch = subprocess.check_output(cmd, cwd=tempdir).decode("utf-8")
        job_id = launch.split("<")[1].split(">")[0]
        root_merge_task_statistic[job_id] = {"path" : nobackup_work_dir ,"NAME" : cmd[4][3:],"STAT" : "UN", "RESULT" : "UN","fail_type" : "UN","level" : "root", "raw" :  cmd[10]} 
    finished_merge = False
    while not finished_merge:
        print("******* merge top level *******")
        out=(subprocess.check_output(["bjobs","-J",f"M{prefix}*"])).decode("utf-8")
        if out == '':
            finished_merge = True
        else:
            actual_task = len(out.splitlines()[1:])
            print(f"Merge top level  has  {actual_task} tasks running ")
            time.sleep(60)
    check_log("root",root_merge_task_statistic)
    subprocess.call( ["cp " + str(nobackup_work_dir / "*_expectscore.tsv") + " " + str(log_dir) ], shell = True)
    subprocess.call( ["cp " + str(nobackup_work_dir / "*_scores_histogram.tsv") + " " + str(log_dir) ], shell = True)
else:
    if len(level1_merge_stat) > 0:
        subprocess.call( ["cp " + str(level1_merge_stat[list(level1_merge_stat.keys())[0]]['path'] / "*.pepXML") + " " + str(work_dir) ], shell = True)
        subprocess.call( ["cp " + str(level1_merge_stat[list(level1_merge_stat.keys())[0]]['path'] / "*_expectscore.tsv") + " " + str(log_dir) ], shell = True)
        subprocess.call( ["cp " + str(level1_merge_stat[list(level1_merge_stat.keys())[0]]['path'] / "*_scores_histogram.tsv") + " " + str(log_dir) ], shell = True)

    else:
        # no merge at all - take the path from the task_report.csv
        temp_df = pd.read_csv(tempdir / (stat_track_file +".csv"))
        if len(temp_df)>1:
            exit("the pipeline assume no merge was done - only one MSFragger run but found more than one")
        else:
            subprocess.call( ["cp " + str(pathlib.Path(temp_df.iloc[0]['path']) / "*.pepXML") + " " + str(work_dir) ], shell = True)
            subprocess.call( ["cp " + str(pathlib.Path(temp_df.iloc[0]['path']) / "*_expectscore.tsv") + " " + str(work_dir) ], shell = True)
            subprocess.call( ["cp " + str(pathlib.Path(temp_df.iloc[0]['path']) / "*_scores_histogram.tsv") + " " + str(work_dir) ], shell = True)
            
    
    
    print("finished merge tasks")

# =============================================================================
#  start philosopher program on the merged output
# =============================================================================
command_philosopher =  ["bsub","-q",args.queue,"-R","rusage[mem=10000]","-J",f"Phi{prefix}","-o","LOG-Philosopher-%J.txt","python3",
                        Philosopher,work_dir, work_dir / fasta_path,
                        "-digest",search_enzyme_name,"-fdr", str(args.fdr)]
print(" ".join([str(x) for x in command_philosopher]))
logging.info("philosopher command: \n" + " ".join([str(x) for x in command_philosopher]) + "\n")
subprocess.run(command_philosopher,cwd=work_dir,check = True)

finished_Philosopher = False
while not finished_Philosopher: 
    out=(subprocess.check_output(["bjobs","-J",f"Phi{prefix}*"])).decode("utf-8")
    if out == '':
        finished_Philosopher = True
    else:
        print("Philosopher task is running ")
        time.sleep(120)

print("finished Philosopher tasks")
logging.info("finished Philosopher tasks, summerize MSFragger tasks\n")
smt.build_MSFragger_execution_table(args.modification_list,tempdir / "task_report.csv",log_dir )
print("--- %s seconds ---" % (time.time() - start_time))
logging.info("--- %s seconds ---\n" % (time.time() - start_time))
