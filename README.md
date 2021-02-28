# MergeQCResultsWithRemovedSamples
Compares list of fastq files for the run, to the list of samples with Consensus sequences to identify what samples were removed. 
These samples are then added back to the QC summary file (appended to bottom), created using mergeQCresults.py. 
We need to track the samples that fail (for reporting to clients) as well as those that pass, 
so it's important to keep track of the files that are too small to analyze and are removed after dehosting etc. 

<b>To Run:</b>

Navigate to directory with output from mergeQCresults.py

(assuming you have pandas installed in an environment called pandas):

conda activate pandas

python3 /path/to/script/CompareFastqsToConsensusLists_AddToMergeQCpyResult.py

conda deactivate
