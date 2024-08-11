# Multiple_DNA_Sequence_Aligner
A collection of the files needed to demonstrate multiple sequence alignment.  
Usage: Download all files included in this repository to a directory on your computer. In the Command Prompt, navigate to this directory, then enter following command:  
                        python Multiple_Sequence_Aligner.py -i MSA_input.fasta -o MSA_complete.fasta -s BLOSUM50.mtx
                        
You will specify an output FASTA file that the final alignment will be written to in the same directory. Example input and scoring matrix files are included in the repository.  
                        Input file: MSA_input.fasta  
                        Scoring Matrices: BLOSUM50.mtx or BLOSUM62.mtx or nucleotide.mtx  
                        Output file: You will specify an output file name e.g. MSA_complete.fasta  
  
MSA_expected_output.fasta is the correct expected output when using the provided input FASTA file and BLOSUM50 scoring matrix. Comparing this file to the output ensures program correctness.  
  
notepad is a good way to view these files
