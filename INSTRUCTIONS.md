## *INSTALLATION INSTRUCTIONS*

### Linux

Open terminal, cd to your preferred location, and type:

    git clone https://github.com/VagelisKalimeris/PARALLEL_STTC.git
    cd PARALLEL_STTC
    make

Make sure directories:

    DATASETS
    ASTROCYTES
    RESULTS
have been created. If any of them is missing, create it yourself. For example:

    mkdir RESULTS

Every psm_avalenche as well as astrocytes input file should be converted from .mat to text, with the following commands in matlab/octave:

    load '<name>.mat'
    dlmwrite('<name>', <matrix>, 'newline', 'unix', 'delimiter', '')
where \<name\> is the name of .mat file and \<matrix\> is the 2D table (frames * cells).

Every dataset should be placed inside DATASETS directory, and every Astrocytes file should be placed inside ASTROCYTES file.
The names of the final input files (dataset + astrocytes) MUST have the same name and no extension.

Run with:

    ./sttc <size_of_control_group> <size_of_Dt> <name_of_dataset>
    
For example, if Dt is 3, control group has size 500, and dataset is M696, then run:
    
    ./sttc 500 3 M696

The resulting files will be accessible after program completion inside RESULTS directory. 
For every psm_avalenche input file we analyze, our code produces 5 different output files (3 csv's and 2 txt's).
The first part of every output file is the same as the name of the corresponding input file.
The second part of every output file lists the analysis parameters(control group, Dt).
The third part of every output file indicates it's type.
There are five types as listed below:

1. motifs.csv: Motif analysis

2. neurons_info.txt: The total number of significant pairs and triplets as well as some general information that were collected during the analysis.

3. neurons_spikes.txt: The spikes found in each neuron.

4. pairs.csv(or tuplets.csv): All the significant pairs, along with their STTC value and their percentile(position among the control group)

5. triplets.csv: All the significant triplets, along with their STTC value and their percentile(position among the control group


### `WINDOWS`

The instructions are the same as long as you have installed a terminal application (we suggest "cygwin"), and added "git" and "make" support.

## **Parameters that you might want to change:**

#### After any change to the source code, you must save the changed file and run the commands "make clean" and "make".

Null distribution of the conditional STTC: If the conditional STTC of a given triplet ABC, is greater than the significant threshold and the number of firing events of ‘reduced A’ is greater than 5, then we consider this triplet as significant.
You might want to change this by opening the file SRC/cond_null_dist.cpp with a text editor, and entering your preferred number at line 65.

To change the significant threshold open the file SRC/common.cpp with a text editor and at line 75, change 3.0 to your preferred value(a real number, not integer is required).
