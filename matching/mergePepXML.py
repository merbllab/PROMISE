# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 15:35:27 2018

@author: assafk
"""
import pathlib
import shlex
from argparse import ArgumentParser
import pandas as pd
import merge as me

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('split_directory_status_file', help="path to split statistics file", metavar='string')
    parser.add_argument('infile', help="mzML file", metavar='string')
    parser.add_argument('tempdir', help="temporary working directory", metavar='string')
    parser.add_argument('output', help="output directory", metavar='string')
    parser.add_argument('-output_report_topN',dest="output_report_topN", default=5, help="output_report_topN", metavar=int)
    parser.add_argument('-output_max_expect',dest="output_max_expect", default=50.0, help="output_max_expect", metavar=float)
    
    args = parser.parse_args()

# =============================================================================+==============================================
# defined path - please update the following links to the external software and the temporary folders according to your need
# ============================================================================================================================
msfragger_jar_path = pathlib.Path(r"/home/labs/yifatlab/assafk/tools/MSFragger_Philosopher/msfragger-20180316.one-jar.jar").resolve()

    
statDf = pd.read_csv(args.split_directory_status_file)
successful_tempdir_parts_str = statDf[statDf["RESULT"].str.contains("Success|Merged")].path.tolist()
successful_tempdir_parts = [pathlib.Path(p) for p in successful_tempdir_parts_str]
output_path = pathlib.Path(args.output).resolve()
tempdir_parts = successful_tempdir_parts
infiles = [pathlib.Path(args.infile)]
tempdir = pathlib.Path(args.tempdir)
jvm_cmd = shlex.split("java -Xmx5G -jar")
msfragger_cmd = jvm_cmd + [msfragger_jar_path]
output_report_topN = int(args.output_report_topN)    
output_max_expect = float(args.output_max_expect)


me.merge(output_path,tempdir_parts,tempdir,msfragger_cmd,output_report_topN,output_max_expect,infiles)
