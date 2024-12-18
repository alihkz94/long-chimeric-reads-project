# Data and scripts for "Overcoming chimeric sequence challenges in long-read amplicon sequencing" (Hakimzadeh et al. 2025)

### Structure
This repository contains the data and part of the analysis stack for the abovementioned paper. It is structured as follows:

[Simulated data](https://github.com/alihkz94/long-chimeric-reads-project/tree/main/Simulated%20data) holds scripts related to the simulated dataset from generating the simulated data, chimeric sequence creation, quality filtering, and chimera filtering related to the simulated dataset. Moreover, the scripts for the simulated dataset and statistical analysis were used to calculate the F1 score.

[Real data](https://github.com/alihkz94/long-chimeric-reads-project/tree/main/Real%20data) holds scripts related to real data analysis.

[BlasCh](https://github.com/alihkz94/long-chimeric-reads-project/tree/main/BlasCh) contains the BLAST scripts for alignment and specific module **BlasCh** designed for processing XML outputs to find false positive chimeras and false negative chimeras.

[Figures & tables](https://github.com/alihkz94/long-chimeric-reads-project/tree/main/Figures%20%26%20tables) contain the scripts used for generating graphs and tables.

 

The workflow we followed for the real dataset was like this:
![workflow for real dataset](workflow.jpg)
