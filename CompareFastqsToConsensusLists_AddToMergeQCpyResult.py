# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 15:02:04 2021

@author: KMacDo
"""


#import packages I need
import os
import subprocess
import fnmatch
import pandas as pd
#If you need to change a directory, you can do it this way (but you don't here):
#os.chdir('C:/Users/KMacDo/Desktop/PythonTest')

#save the current working directory (cwd) to a variable to use in everything below. 
#For us, this would be the MiSeqRunID directory for each run - it changes each time we analyze a run, 
#so i want to pull this from where ever I am, 
#so I don't have to enter it each time:
cwdPath = os.getcwd()



#Define variable using last part of directory/path
#This will be used to name your files uniquely:
MiSeqRunID = os.path.basename(os.path.normpath(cwdPath))
#print(MiSeqRunID)  #works

#Run bash commands to generate the 2 files to compare (FastqList and ConsensusList):
bashCommandTest = "ls"
bashCommand1 = "ls ../../analysis_by_run/" + MiSeqRunID + "/ncov2019-artic-nf-v1.1-output/ncovIllumina_sequenceAnalysis_makeConsensus/*.primertrimmed.consensus.fa | cut -d/ -f7 | cut -d. -f1 >ConsensusList.txt" 
bashCommand2 = "ls ../../direct_fastq_symlinks_by_run/" + MiSeqRunID + "/*.fastq.gz | cut -d/ -f5 | cut -d. -f1 | cut -d_ -f1 >FastqList.txt"
#print(bashCommandTest)
#print(bashCommand1)
#print(bashCommand2)
#os.system(bashCommand1)
#os.system(bashCommand2)
#os.system is deprecated since python 2.5. Use subprocess instead:
#subprocess.run(bashCommandTest)
subprocess.run(bashCommand1, shell=True, check=True)
subprocess.run(bashCommand2, shell=True, check=True)

#df_FastqList1 = pd.read_table("FastqList.txt", header=None)


#create list of sampleIDs to remove (e.g. Undetermined) from lists:
RemoveLines = ['Undetermined', 'undetermined']
#Find those lines in the FastqList.txt file and create a new file without them:
with open('FastqList.txt') as df_FastqList1, open('FastqList2.txt', 'w') as df_FastqList2:
    for line in df_FastqList1:
        if not any(RemoveLines in line for RemoveLines in RemoveLines):
            df_FastqList2.write(line)

#Read in lists of sampleIDs:
df_FastqList2 = pd.read_table("FastqList2.txt", header=None)
#print(df_FastqList2)
#Remove duplicate sampleIDs in Fastq List:
df_FastqList3 = df_FastqList2.drop_duplicates([0], keep='first') 
#print(df_FastqList3)
#os.path.join(cwdPath, "FastqList.txt")
df_ConsensusList = pd.read_table("ConsensusList.txt", header=None)
#print(df_ConsensusList)


#Extract sampleIDs that are in FastqList but NOT in ConsensusList:
#df_FastqDiffs = df_ConsensusList[df_ConsensusList.isin(df_FastqList3)].dropna(how = 'all')
df_FastqDiffs1 = df_FastqList3.loc[~df_FastqList3[0].isin(df_ConsensusList[0])].copy()
#print(df_FastqDiffs1)
#df3 = df1.loc[~df1['ID'].isin(df2['ID'])].copy()


#Add column headers to match the other file:
#        ,sample,run_name,num_consensus_snvs,num_consensus_n,num_consensus_iupac,num_variants_snvs,num_variants_indel,num_variants_indel_triplet,mean_sequencing_depth,median_sequencing_depth,qpcr_ct,collection_date,num_weeks,scaled_variants_snvs,genome_completeness,qc_pass_x,lineage_x,watch_mutations,watchlist_id,num_observed_mutations,num_mutations_in_watchlist,proportion_watchlist_mutations_observed,note,pangoLEARN_version,pct_N_bases,pct_covered_bases,longest_no_N_run,num_aligned_reads,qc_pass_y,sample_name
df_FastqDiffs2 = df_FastqDiffs1.reindex(columns=[*df_FastqDiffs1.columns.tolist(),'run_name','num_consensus_snvs','num_consensus_n','num_consensus_iupac','num_variants_snvs','num_variants_indel','num_variants_indel_triplet','mean_sequencing_depth','median_sequencing_depth','qpcr_ct','collection_date','num_weeks','scaled_variants_snvs','genome_completeness','qc_pass_x','lineage_x','watch_mutations','watchlist_id','num_observed_mutations','num_mutations_in_watchlist','proportion_watchlist_mutations_observed','note','pangoLEARN_version','pct_N_bases','pct_covered_bases','longest_no_N_run','num_aligned_reads','qc_pass_y','sample_name'], fill_value=0)
#print(df_FastqDiffs2)
#Rename Column 0 as 'sample'
df_FastqDiffs3 = df_FastqDiffs2.rename(columns={0: "sample"})
#print(df_FastqDiffs3)

#replace 0's in qc_pass_y column with FALSE:
#df_FastqDiffs = pd.df_FastqDiffs3
df_FastqDiffs4 = df_FastqDiffs3[['qc_pass_y']].replace(0,'FALSE')
df_FastqDiffs4b = df_FastqDiffs3[['run_name']].replace(0, MiSeqRunID)
#print(df_FastqDiffs4b) #worked
#check content:
#df_FastqDiffs4.to_csv('df_FastqDiffs4.csv') #has 2 columns, and has a header ('sample')

#merge df_FastqDiffs4 with df_FastqDiffs3
df_FastqDiffs4_merge0 = [df_FastqDiffs3.iloc[:, 0:1], df_FastqDiffs4b.iloc[:, 0:1]]
df_FastqDiffs4_merge1 = pd.concat(df_FastqDiffs4_merge0, axis=1)                      
#print(df_FastqDiffs4_merge1)
df_FastqDiffs4_merge2 = [df_FastqDiffs3.iloc[:, 2:28], df_FastqDiffs4.iloc[:, 0:1]]
df_FastqDiffs4_merge3 = pd.concat(df_FastqDiffs4_merge2, axis=1)                      
#print(df_FastqDiffs4_merge3)
df_FastqDiffs4_merge4 = [df_FastqDiffs4_merge1, df_FastqDiffs4_merge3]
df_FastqDiffs4_merge5 = pd.concat(df_FastqDiffs4_merge4, axis=1)                      
#print(df_FastqDiffs4_merge5)
#df_FastqDiffs4_merge5.to_csv('df_FastqDiffs4_merge5.csv')

#replace 0's in sample_name column with the values in sample column:
df_FastqDiffs5 = df_FastqDiffs3.loc[df_FastqDiffs3["sample_name"] < 1, "sample_name"] = df_FastqDiffs3["sample"]
df_FastqDiffs5 = df_FastqDiffs5.reset_index().rename(columns={'sample': 'sample_name'})
#print(df_FastqDiffs5)
#df.loc[df['A'].isnull(), 'A'] = df['B']
#check content:
#df_FastqDiffs5.to_csv('df_FastqDiffs5.csv') #has a header ('sample'), just doesn't show in preview for whatever reason

#merge df_FastqDiffs5 with df_FastqDiffs4_merge5
df_FastqDiffs3_4_5_merge2 = pd.merge(df_FastqDiffs4_merge5.iloc[:, 0:29], df_FastqDiffs5.iloc[:, 1:2], how='left', left_on='sample', right_on='sample_name')
#print(df_FastqDiffs3_4_5_merge2) 
#check content:
#df_FastqDiffs3_4_5_merge2.to_csv('df_FastqDiffs3_4_5_merge2.csv') #has a header ('sample'), just doesn't show in preview for whatever reason


#Combine the columns I want:
# df_FastqDiffs6 = df_FastqDiffs3.iloc[:, 0:28]
# print(df_FastqDiffs6)
# df_FastqDiffs7 = df_FastqDiffs4.iloc[:, 0:1]
# print(df_FastqDiffs7)   
# df_FastqDiffs8 = df_FastqDiffs5.iloc[:, 0:1]
#print(df_FastqDiffs8)
# df_FastqDiffs9 = df_ncov_variant_lineage_artic_merge.iloc[:, 30:35]
# df_FastqDiffs10 = df_ncov_variant_lineage_artic_merge.iloc[:, 29:30]
# #print(df_FastqDiffs9)
# #print(df_FastqDiffs10)


for dir_path, dir_names, file_names in os.walk(cwdPath):
    for f in file_names:
        if fnmatch.fnmatch(f, '*_QC_lineage_VoC_OrderedFinal.csv'):
            #print(f)  # worked 
            df_ncovtoolsSummary1 = os.path.join(cwdPath, f)
            df_ncovtoolsSummary = pd.read_csv(df_ncovtoolsSummary1)
            #print(df_ncovtoolsSummary) #worked 


#MERGE the mergeQCresults.py result file with the Missing SampleIDs:
df_ncovtoolsSummary2 = df_ncovtoolsSummary.iloc[:, 1:32]
#print(df_ncovtoolsSummary2)
df_ncovtoolsSummary_plusMissing1 = [df_ncovtoolsSummary2, df_FastqDiffs3_4_5_merge2]
df_ncovtoolsSummary_plusMissing = df_ncovtoolsSummary2.append(df_FastqDiffs3_4_5_merge2, ignore_index=True)
    
# df_ncovtoolsSummary_plusMissing = pd.merge(df_ncovtoolsSummary, df_FastqDiffs, how='left', left_on='sample', right_on='0')
#print(df_ncovtoolsSummary_plusMissing)

#df_ncov_variant_merge = pd.merge(df_ncovtools, df_variantsum3, how='left', left_on='sample', right_on='sample_id')
#print(df_ncov_variant_merge) #works

#Read in final2 file from mergeQCresults.py and then
#save final3 file as csv (easier to open in Excel than tsv):
#(this has only the columns I want to link my dashboard to, for various people's purposes)

df_ncovtoolsSummary_plusMissing.to_csv(MiSeqRunID + '_' + 'MissingPlus_QC_lineage_VoC_OrderedFinal.csv')

#Run Bash Commands to Remove Unnecessary Files: 
bashCommand3 = "rm FastqList.txt; rm FastqList2.txt; rm ConsensusList.txt" 
#bashCommand4 = ""
subprocess.run(bashCommand3, shell=True, check=True)
#subprocess.run(bashCommand4, shell=True, check=True)