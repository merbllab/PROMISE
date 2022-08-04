# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 17:15:43 2018
command line example:
@author: assafk
"""
from __future__ import division
import pandas as pd
from Bio import SeqIO
from argparse import ArgumentParser
import os
import add_ranks_to_psm_file as ra
import david as d 
import PTM_window_assignment
import Philosopher_assignment_fix as pf

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('work_dir',  help="working directory", metavar='string')
    parser.add_argument('MSFragger_output',  help="psm file format", metavar='string')
    parser.add_argument('output_file',  help="output csv file", metavar='string')
    parser.add_argument('reference_data', help="reference data file", metavar='string')    
    parser.add_argument('-mode', dest="mode", default = "new", choices=['new', 'append'], help="new run or append on existing run", metavar='string')
    parser.add_argument('-annotation', dest="annotation",default = "general", help="determine dictionary for masses seperate by ,",metavar='string')

    args = parser.parse_args()


# =============================================================================+==============================================
# defined path - please update the following links to the external software and the temporary folders according to your need
# ============================================================================================================================
modification_mapping_path = r"/home/labs/yifatlab/assafk/PROMISE/prioritizing/PTM/modification_mapping.xlsx"     
modification_path = r"/home/labs/yifatlab/assafk/PROMISE/prioritizing/PTM"

###### if mode = append check what step we are ###############
step = 0

if args.mode =="append": 
    file_list = []       
    base_name = os.path.basename(args.output_file)
    for file in os.listdir(args.work_dir):
        if file.startswith(base_name[:-4] + "_"):
            file_list.append(file)
    print(file_list)
    for file in file_list:
        if file.endswith(base_name[:-4] + "_valid_modification.csv"): 
            if step < 10 : step = 10 # do nothing 
        elif file.endswith(base_name[:-4] + "_step8.csv"): 
            if step < 9 : step = 9 
        elif file.endswith(base_name[:-4] + "_collaps.csv"):
            if step < 8 : step = 8
        elif file.endswith(base_name[:-4] + "_step6.csv"):
            if step < 7 : step = 7
        elif file.endswith(base_name[:-4] + "_step5.csv"):
            if step < 6 : step = 6
        elif file.endswith(base_name[:-4] + "_step4.csv"):
            if step < 5 : step = 5
        elif file.endswith(base_name[:-4] + "_step3.csv"):
            if step < 4 : step = 4
        elif file.endswith(base_name[:-4] + "_step2.csv"):
            if step <3: step = 3
        elif file.endswith(base_name[:-4] + "_step1.csv"):
            if step <2 : step = 2
        elif file.endswith(base_name[:-4] + "_step0.csv"):
            if step <1 : step = 1
        else:
            step = 0
        #print(step)
            
    print ("script in append mode, start in step " + str(step)) 
    
if step == 0: ########## read output , filter only modified peptide and extract Uniprot Id ######### 
    '''
    this step perform:
        1. filter for only peptides with modification
        2. add lower ranks hits directly from the pepXML file if their delta score from the rank 1 is less than 1
           and identification is not decoy sequance (rev_ prefix)
    added columns: raw file , rank, delta score  
    '''
    print("step 0 : add lower hit ranks to philosopher output")
    ra.add_lower_hit_ranks(args.work_dir,args.MSFragger_output , args.output_file[:-4] + "_step0.csv")
    step+=1
    
if step == 1: ########## read output , filter only modified peptide and extract Uniprot Id ######### 
    '''
    this step fix philosopher modification assignment when total mass (aa + modification) has more than one hit
    for example D + Me = Q +Deamination
    added columns: Philosopher fixed bug
    '''
    print("step 1 : fix bug of philosopher with modification assigmnet")
    pf.fix_Pilosopher_modification_assignment(args.work_dir,args.output_file[:-4] + "_step0.csv" , args.output_file[:-4] + "_step1.csv")
    step+=1
    
if step ==2: ########### duplicate scans ID for each mapped protein ##################
    '''
    this step duplicate the rows for each uniprot ID , those extra lines will be collaps later if there are identical
    added columns: UP-ID
    '''
    print("step 2 : duplicate scans ID for each mapped protein")       
    modifiedDf = pd.read_csv(args.output_file[:-4] + "_step1.csv", na_filter = False )
    def get_UP_ID(record_header):
        try:
            if record_header.startswith("CON"):
                return record_header.split('_')[1]
            else:
                return record_header.split('|')[1]
        except IndexError:
            print("could not get the UP-ID from protein name: " + record_header )    
            return "UN"
        
    modifiedDf['UP-ID'] = modifiedDf.apply(lambda row: get_UP_ID(row['Protein']) ,axis=1)
            
    for index, row in modifiedDf.iterrows():
        if row['Mapped Proteins'] != "" : # no additonal proteins 
            additional_protein_list = row['Mapped Proteins'].split(',')
            # remove reduandants and clean list
            unique_additional_protein_set = set()
            for protein in additional_protein_list:
                p = protein.strip()
                if p != ' ' and p != '':
                    unique_additional_protein_set.add(p)
            for protein in unique_additional_protein_set:
                new_row = row
                new_row['Protein']=protein
                new_row['Mapped Proteins']=""
                try:
                    if new_row['Protein'].startswith("CON"):
                        new_row['UP-ID'] = new_row['Protein'].split('_')[1]
                    else:
                        new_row['UP-ID'] = new_row['Protein'].split('|')[1]
                    modifiedDf=modifiedDf.append(new_row)
                except IndexError:
                    print("could not get the UP-ID from protein name: " + protein )
    modifiedDf.to_csv(args.output_file[:-4] + "_step2.csv",index = False) 
    step+=1
    
modifciation_mappingDf = pd.read_excel(modification_mapping_path, headers = False)
modifciation_mappingDf = modifciation_mappingDf.fillna("")

if step == 3: ############# reteaving information from the sequance ##################
    '''
    this step perform:
        1. go over the fasta file and extract sequance information
        2. annotate the modification mass to name according to modification_mapping file
        3. translate peptide modification relative position to protein global position and store data at Modification data column
    added columns: 
        'Start position', 'End position', 'Length','Modificatyion data'
          'N cleavage window','C cleavage window', 'Amino acid before','First amino acid','Second amino acid','Second last amino acid',
          'Last amino acid','Amino acid after','Protein names'
    '''
    print("step 3 : retrieving information from the sequance") 
      
    modifiedDf = pd.read_csv(args.output_file[:-4] + "_step2.csv", na_filter = False )
    newColumns = ['Start position', 'End position', 'Length','Modificatyion data',
                  'N cleavage window','C cleavage window', 'Amino acid before','First amino acid','Second amino acid','Second last amino acid',
                  'Last amino acid','Amino acid after','Protein names']
    modifiedDf = modifiedDf.reindex(columns = modifiedDf.columns.tolist() + newColumns)
    modifiedDf = modifiedDf.fillna("")
    for record in SeqIO.parse(args.reference_data, "fasta"):
        #print(record.id)
        subsetModifiedDf = modifiedDf[modifiedDf['Protein'] == record.id]
        for index, row in subsetModifiedDf.iterrows():
            pep = row['Peptide']
            padd_seq = str(record.seq).ljust(len(str(record.seq))+8,"_")
            padd_seq = padd_seq.rjust(len(padd_seq)+8,"_")
            hit = record.seq.find(pep)
            if hit == -1:
                print("error - couldn't find the peptide " + row['Peptide'] + " in the protein " + row['Protein'])
                continue
            start = hit+1
            end = hit + len(pep)
            hit_again = record.seq.find(pep,end)
            if hit_again != -1:
                print("error - find the peptide in two location in the protein, spectrum: " + row['Spectrum'])
            padd_hit = padd_seq.find(pep)
            padd_start = padd_hit +1
            padd_end = padd_hit + len(pep)
            N_cleavage_window = padd_seq[padd_start-9:padd_start+7]
            C_cleavage_window = padd_seq[padd_end-8:padd_end+8]
            modifiedDf.loc[index,'C cleavage window' ] = C_cleavage_window
            modifiedDf.loc[index,'N cleavage window' ] = N_cleavage_window
            #UP_ID = modifiedDf.loc[modifiedDf['Spectrum'] == Spectrum,'Protein'].iloc[0].split('|')[1]
            #modifiedDf.loc[modifiedDf['Spectrum'] == Spectrum,'UP-ID'] = UP_ID
            modifiedDf.loc[index,'Start position' ] = start
            modifiedDf.loc[index,'End position' ] = end
            modifiedDf.loc[index,'Length' ] = len(pep)
            if modifiedDf.loc[index,'Assigned Modifications'] != "":
                modification_list = modifiedDf.loc[index,'Assigned Modifications'].split(',')
                #print(modification_list)
                modification_data_string = []
                for mod in modification_list:
                    #if modification in mod:
                    temp = mod.strip().split("(")
                    try:
                        amino_acid = temp[0][-1]
                    except IndexError:
                        print("fail to extract aa " + str(temp) + "," + mod)
                    if amino_acid == "n" or amino_acid =="m" : # special case where the modifcation is at the end of the protein \ peptide
                        absolute_position = hit+1
                        amino_acid = pep[0]
                    else:
                        try:
                            relative_position = int(temp[0][:-1])
                        except ValueError:
                            print("cound't convert to int " + ",".join(temp) + " mod = " + mod + ", index = " + str(index))
                            exit(1)
                        
                        absolute_position = hit + relative_position
                    mass = float(temp[1].strip(")"))
                    #print(mass)
                    try:
                        modification_name_df = modifciation_mappingDf.loc[modifciation_mappingDf['mass'] == mass]
                        if len(modification_name_df) ==1:
                            modification_name = modification_name_df['modification name'].iloc[0]
                        else:
                            modification_name = modification_name_df.loc[modifciation_mappingDf['aa'] == amino_acid,'modification name'].iloc[0]
                    except IndexError:
                        print("coudn't match modification name for mass : " + str(mass))
                        modification_name = "Unknown"
                    modification_data_string.append(str(absolute_position) + "," + amino_acid + "," + modification_name)
                modifiedDf.loc[index,'Modificatyion data' ] = ";".join(modification_data_string)
            if start > 1 :
                modifiedDf.loc[index,'Amino acid before' ] = record.seq[hit-1]
            else:
                modifiedDf.loc[index,'Amino acid before' ] = "_"
            modifiedDf.loc[index,'First amino acid' ] = record.seq[hit]
            modifiedDf.loc[index,'Second amino acid' ] =  record.seq[hit+1]
            modifiedDf.loc[index,'Second last amino acid' ] = record.seq[end-2]
            modifiedDf.loc[index,'Last amino acid' ] = record.seq[end-1]
            if end == len(record.seq):
                modifiedDf.loc[index,'Amino acid after' ] = "_"
            else:
                modifiedDf.loc[index,'Amino acid after' ] = record.seq[end]
            
            modifiedDf.loc[index,'Protein names' ] = record.description 
    
    
    modifiedDf.to_csv(args.output_file[:-4] + "_step3.csv",index = False) 
    step+=1


########## utility method ###########################################
#def calculate_tail_mass(tail): 
#    aa_list = list(tail)
#    mass = 0
#    for aa in aa_list:
#        if aa == "_":
#            mass+=0
#        else:
#            try:
#                mass+= masses.loc[aa,'mass']
#            except KeyError:
#                print("un identified aa was presented in the protein " + aa)
#    return mass


if step ==4: ######## look for decoy modification ###############
    '''
    this step search for decoy mass at after or before the peptide that equal the modification mass
    it is using modification_masses file to calculate the decoy mass
    added columns: 'modification decoy'
    '''
    print("step 4 : doing nothing - was replaced by module search for alternative modification") 
    modifiedDf = pd.read_csv(args.output_file[:-4] + "_step3.csv", na_filter = False )     
    newColumns = ['modification decoy']
    modifiedDf = modifiedDf.reindex(columns = modifiedDf.columns.tolist() + newColumns)
    modifiedDf = modifiedDf.fillna("")        
    modifiedDf.to_csv(args.output_file[:-4] + "_step4.csv",index = False) 
    step+=1
    
if step ==5: ######## validate modifciation in dpPTM ############
    '''
    this step seach for the assign modification in dbPTM file and mark + if it is found
    it usses modification_mapping file to know what is dbPTM annotation for each modification 
    added columns:
        dbPTM Modification validation','dbPTM Modification within peptide'   
    '''
    print("step 5 : validate modifciation in dpPTM")       
    modifiedDf = pd.read_csv(args.output_file[:-4] + "_step4.csv", na_filter = False )
    newColumns = ['dbPTM Modification validation','dbPTM Modification within peptide']
    modifiedDf = modifiedDf.reindex(columns = modifiedDf.columns.tolist() + newColumns)
    modifiedDf = modifiedDf.fillna("")
    UP_ID_set = set(modifiedDf['UP-ID'].unique())
    modification_file = open(os.path.join(modification_path,"dbPTM.txt"))
    for l in modification_file:
        vals = l.strip().split("\t")
        if vals[1] in UP_ID_set:
            # print ("found one " + vals[1] + "," + vals[2] + "," + vals[6])
            subset_modifiedDf = modifiedDf[modifiedDf['UP-ID']==vals[1]]
            #print(subset_modifiedDf['Spectrum'])
            for index, row in subset_modifiedDf.iterrows():
                # print("index = " , index)
                # print(row)
                if row['dbPTM Modification validation' ] != "":
                    validation_string = row['dbPTM Modification validation' ].split(";")
                else:
                    validation_string = []
                if row['dbPTM Modification within peptide'] !="":
                    modification_within_peptide_string = row['dbPTM Modification within peptide'].split(";")
                else:
                    modification_within_peptide_string =[]
                modification_data = row['Modificatyion data' ].split(";")
                start = int(row['Start position']) if row['Start position'] != '' else 0
                end = int(row['End position']) if row['End position'] != '' else  0
                #print(start,end)
                temp_validation_string = ""
                for mode in modification_data:
                    if mode != '':
                        position = mode.split(",")[0]
                        aa = mode.split(",")[1]
                        modification_name = mode.split(",")[2]
                        # convert modification name to dbPTM format
                        try:
                            modification_name_dbPTM_format = modifciation_mappingDf[modifciation_mappingDf['modification name']==modification_name]['modification name dbPTM'].iloc[0]
                        except IndexError:
                            modification_name_dbPTM_format = ""
                            print(modification_name)
                        #print(position,aa)
                        if aa == vals[6] and position == vals[2] and modification_name_dbPTM_format == vals[7] :
                            temp_validation_string += "+"
                        else:
                            temp_validation_string += "-"
                if temp_validation_string.find("+") != -1: 
                    validation_string.append(temp_validation_string)
                            
                if (int(vals[2]) <= end and int(vals[2])>=start):
                    modification_within_peptide_string.append(vals[2] + vals[6] + "," + vals[7])
                #print("validation string = " + validation_string)
                #print("mod string = " + modification_within_peptide_string)
                modifiedDf.loc[(modifiedDf['Spectrum'] == row['Spectrum']) & (modifiedDf['UP-ID'] == row['UP-ID']) , 'dbPTM Modification validation'] = ";".join(validation_string)
                modifiedDf.loc[(modifiedDf['Spectrum'] == row['Spectrum']) & (modifiedDf['UP-ID'] == row['UP-ID']) , 'dbPTM Modification within peptide'] = ";".join(modification_within_peptide_string )       
    modification_file.close()            
    modifiedDf.to_csv(args.output_file[:-4] + "_step5.csv",index = False) 
    step+=1


if step ==6: ######## validate modifciation in PhosphoSitePlus db ############
    '''
    this step seach for the assign modification in PhosphoSitePlus file and mark + if it is found
    it usses modification_mapping file to know what is PhosphoSitePlus annotation for each modification 
    added columns:
        dbPTM Modification validation','dbPTM Modification within peptide'   
    '''

    print("step 6 : validate modifciation in PhosphoSitePlus db")       
    modifiedDf = pd.read_csv(args.output_file[:-4] + "_step5.csv", na_filter = False )
    newColumns = ['PhosphoSitePlus Modification validation','PhosphoSitePlus Modification within peptide']
    modifiedDf = modifiedDf.reindex(columns = modifiedDf.columns.tolist() + newColumns)
    modifiedDf = modifiedDf.fillna("")    
    UP_ID_set = set(modifiedDf['UP-ID'].unique())
    modification_file = open(os.path.join(modification_path,"PhosphoSitePlusDb_05_2018.txt"))
    for l in modification_file:
        vals = l.strip().split("\t")
        if vals[2] in UP_ID_set:
            db_aa = vals[4][0]
            db_pos = vals[4].split("-")[0][1:]
            db_mod = vals[4].split("-")[1]
            # print ("found one " + vals[1] + "," + vals[2] + "," + vals[6])
            subset_modifiedDf = modifiedDf[modifiedDf['UP-ID']==vals[2]]
            #print(subset_modifiedDf['Spectrum'])
            for index, row in subset_modifiedDf.iterrows():
                # print("index = " , index)
                # print(row)
                if row['PhosphoSitePlus Modification validation' ] != "":
                    validation_string = row['PhosphoSitePlus Modification validation' ].split(";")
                else:
                    validation_string =[]
                if  row['PhosphoSitePlus Modification within peptide'] != "": 
                    modification_within_peptide_string = row['PhosphoSitePlus Modification within peptide'].split(";")
                else:
                    modification_within_peptide_string = []
                modification_data = row['Modificatyion data' ].split(";")
                start = int(row['Start position']) if row['Start position'] != '' else 0
                end = int(row['End position']) if row['End position'] != '' else  0
                #print(start,end)
                temp_validation_string = ["-"] * len(modification_data)
                for i in range(len(modification_data)):
                    mode = modification_data[i]
                    if mode != '':
                        position = mode.split(",")[0]
                        aa = mode.split(",")[1]
                        modification_name = mode.split(",")[2]
                        # convert modification name to PhosphoSitePlus format
                        try:
                            modification_name_PhosphoSitePlus_foramt = modifciation_mappingDf[modifciation_mappingDf['modification name']==modification_name]['modification name PhosphoSitePluse'].iloc[0]
                        except IndexError:
                            modification_name_PhosphoSitePlus_foramt = ""
#                            print(modification_name)
                        for m in modification_name_PhosphoSitePlus_foramt.split(","):
                            #print(position,aa)
                            if aa == db_aa and position == db_pos and m == db_mod :
                                temp_validation_string[i] = "+"
                if "+" in temp_validation_string: 
                    validation_string.append("".join(temp_validation_string))
                if (int(db_pos) <= end and int(db_pos)>=start):
                    modification_within_peptide_string.append(vals[4])
                #print("validation string = " + validation_string)
                #print("mod string = " + modification_within_peptide_string)
                modifiedDf.loc[(modifiedDf['Spectrum'] == row['Spectrum']) & (modifiedDf['UP-ID'] == row['UP-ID']) , 'PhosphoSitePlus Modification validation'] = ";".join(validation_string)
                modifiedDf.loc[(modifiedDf['Spectrum'] == row['Spectrum']) & (modifiedDf['UP-ID'] == row['UP-ID']) , 'PhosphoSitePlus Modification within peptide'] = ";".join(modification_within_peptide_string)       
    modification_file.close()            
    modifiedDf.to_csv(args.output_file[:-4] + "_step6.csv",index = False) 
    step+=1


def fill_row_score(row):
    score = 0
    if row['Protein'].startswith("sp"):
        score +=10
        if row['Protein'].split("|")[1].find("-") != -1:
            score-=5
    if row['dbPTM Modification validation'] != "" : score+=20
    if row['PhosphoSitePlus Modification validation'] != "" : score+=20
    if row['dbPTM Modification within peptide'] != "" : score +=5
    if row['PhosphoSitePlus Modification within peptide'] != "" : score +=5
    if row['Start position'] =="" : score-=20 # couldn't match the sequance to the protein
    return score 

if step ==7: ######## collapsing results  ############
    print("step 7 : collapsing results")       
    modifiedDf = pd.read_csv(args.output_file[:-4] + "_step6.csv", na_filter = False )
    collapsDf = pd.DataFrame(columns = modifiedDf.columns.tolist()) # create empty table to collect the collaps rows
    spectrum_rank_set = set(zip(modifiedDf['Spectrum'].tolist(),modifiedDf['rank'].tolist()))
    
    for spectrum_rank in spectrum_rank_set:
        spectrumDf = modifiedDf[(modifiedDf['Spectrum']==spectrum_rank[0]) & (modifiedDf['rank']==spectrum_rank[1]) ]
        # the final goal is to have one row per spectrum ID. the logic witch row to pick is build on scoring system according to the following logic:
        # swiss prot entry = +10 points
        # swiss prot isoform = -5 points
        # validation of modification in any db = +20 points
        # additioanl modification with in peptide = +5 points 
        # at the end choose the row with the highest score and add all the other uniprot ID as additional hits 
        spectrumDf = spectrumDf.reindex(columns = spectrumDf.columns.tolist() + ['row score'])
        spectrumDf = spectrumDf.fillna("")
        spectrumDf['row score'] = spectrumDf.apply(lambda row: fill_row_score(row),axis=1)
        max_score_spectrumDf = spectrumDf.loc[spectrumDf['row score'].idxmax()]
        # create a list of alternative proteins 
        mapped_proteins = spectrumDf['Protein'].tolist()
        mapped_proteins.remove(max_score_spectrumDf['Protein'])
        max_score_spectrumDf['Mapped Proteins'] = ", ".join(mapped_proteins)
        max_score_spectrumDf = max_score_spectrumDf.drop('row score')
        
        
        collapsDf = collapsDf.append(max_score_spectrumDf) # if there are more than one row with max value - take the first 
        # for debug:
        # print("spectrum  :" +  spectrum + " len = " + str(spectrumDf.shape))
        
    collapsDf.to_csv(args.output_file[:-4] + "_collaps.csv",index = False) 
    step+=1

if step ==8: ######## validate spectrum ions surrounding the modification  ############
    '''
    this step check spectrum quality by annotating the peaks from the mzML file and use them to validate modification position
    added columns:
      'matched ions david','y', 'b', 'ion intensity ratio','matched ions with intensity above average','unassigned ions with intensity above average',
      'matched ions in top 20 peaks','unassigned ions in top 20 peaks',
      "PTM localization",'PTM Reporter ions'
    '''
    print("step 8 : validate spectrum ions surrounding the modification") 
    d.analyzed_ptm_ions(args.output_file[:-4] + "_collaps.csv",args.output_file[:-4] + "_step8.csv",args.work_dir)  
    step+=1

if step ==9: ######## validate spectrum ions surrounding the modification  ############
    '''
    this step check define a PTM window based on the pick annotation in step 8 
    added columns:
      'PTM windows'
    '''
    print("step 9 : define PTM window localization based on pick annotation") 
    PTM_window_assignment.add_PTM_window(args.output_file[:-4] + "_step8.csv",args.output_file[:-4] + "_valid_modification.csv") 
    step+=1


