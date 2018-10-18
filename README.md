## *INSTALLATION INSTRUCTIONS*

### `Linux`

Open terminal, go  to your preferred location(cd WHERE/YOU/WANT/TO/INSTALL), and run:

    git clone https://github.com/VagelisKalimeris/PARALLEL_STTC.git
    cd PARALLEL_STTC
    make

Make sure the directories:

    DATASETS
    ASTROCYTES
    RESULTS
    OBJ
have been created. If any of them is missing, create it yourself. For example:

    mkdir RESULTS

Every psm_avalanche as well as astrocytes input file should be converted from .mat to text, with the following commands in Matlab/Octave:

    load '<name>.mat'
    dlmwrite('<name>', <matrix>, 'newline', 'unix', 'delimiter', '')
where \<name\> is the name of .mat file and \<matrix\> is the 2D table (frames * cells).
If psm_avalenche has 3 dimentions, it needs to be converted to 2D (in Matlab) before the analysis, by choosing a specific one from the 3rd dimention.

The names of the final input files (dataset + astrocytes) MUST have identical names, and no extensions.

Every dataset should be placed inside DATASETS directory, and every Astrocytes file should be placed inside ASTROCYTES directory.
In case of datasets without an astrocytes file, the program runs successfully after printing a wrning message.

### `WINDOWS`

The instructions are the same as long as you have installed a terminal application (we suggest "cygwin"), and added "git", "make" and gcc/g++ support.

## *EXECUTION INSTRUCTIONS*

Run:

    ./sttc <size_of_control_group> <size_of_Dt> <name_of_dataset>
    
For example, if Dt is 3, control group has size 500, and the dataset name is M696, then run:
    
    ./sttc 500 3 M696

The resulting files will be accessible after a successfull execution, inside the RESULTS directory. 
For every psm_avalanche input file we analyze, our code produces 5 different output files (3 csv's and 2 txt's).
The first part of every output file is the same as the name of the corresponding input file.
The second part of every output file lists the analysis parameters(control group, Dt).
The third part of every output file indicates it's type.
There are five types as listed below:

1. motifs.csv: Motif analysis

2. neurons_info.txt: The total number of significant pairs and triplets as well as some general information that were collected during the analysis.

3. neurons_spikes.txt: The spikes found in each neuron.

4. pairs.csv: All the significant pairs, along with their STTC value and their percentile(position among the control group)

5. triplets.csv: All the significant triplets, along with their STTC value and their percentile(position among the control group

To delete any previously produced output files run:

    ./delete_results 
    
## **Parameters that you might want to change:**

#### After any change to the source code, you must save the changed file and run the commands "make clean" and "make".

Null distribution of the conditional STTC: If the conditional STTC of a given triplet (A->B)|C is greater than the significant threshold and the number of firing events of ‘reduced A’ is greater than 5, then this triplet is considered as significant.  
To change the significant threshold, open the file `SRC/common.cpp` with a text editor and at line 75, enter the preferred value as a real number, not integer.  
To change the value of second condition, open the file SRC/cond_null_dist.cpp with a text editor and at line 65, enter the preferred value as an integer number, not real.


### For further information read STTC_conditional.pdf
