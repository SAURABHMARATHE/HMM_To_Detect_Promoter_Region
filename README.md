# HMM_To_Detect_Promoter_Region
HMM code in python to detect promoter regions in DNA nucelotide sequence

This program runs a HMM on a given sequence of nucleotides to find out CgP islands locations.

Few things to know while running the program:

1. The program will prompt for choice in the begining, please input 1 for HMM with pre-defined probabilities and no_of_iteraions=0 for HMM training. Please input 2 for HMM with pre-defined probabilities and no_of_iterations=50

2. Please, keep in mind that every time you run the program, the model will fit data again and both Viterbi and Baum Welch will run again. Hence result will change each time.

3. But the probabilities are carefully chosen after running the program multiple times, and getting good amount of valid CgP islands in the result sequence of nucleotides. (But this is not always guaranteed)

4. Also, please keep in mind that the start index and the last index represent the second and second-last nucleotide in the nucleotide sequence I am printing. I am printing the first and last nucleotide just to show start and end of CpG islands.

5. The model score is printed between nucleotide sequence of Veiterbi and nucleotide sequence for Baum Welch. Also, I am printing the model score at the end once again.

