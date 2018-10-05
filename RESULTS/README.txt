For every psm_avalenche input file we analyze, our code produces 5 different output files(3 csv's and 2 txt's).
The first part of every output file is the same as the name of the corresponding input file.
The second part of every output file lists the analysis parameters(control group, Dt).
The third part of every output file indicates it's type.
There are five types as listed below:
	a) _motifs.csv: Motif analysis
	b) _neurons_info.txt: The total number of significant pairs and triplets as well as some general information that were collected during the analysis.
	c) _neurons_spikes.txt: The spikes found in each neuron.
	d) _pairs.csv(or _tuplets.csv): All the significant pairs, along with their STTC value and their percentile(position among the control group)
	e) _triplets.csv: All the significant triplets, along with their STTC value and their percentile(position among the control group)