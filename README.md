This repository contains processed data and original code for the Neuron publication Forro et al 2023: 

The file NeuronsFR.mat is a cell containing processed Firing rate data (matrices with individual cells as rows, 
columns with spatial bins) for the population of pyramidal cells, unidentified interneurons, PV+basket, 
PV+bistratified and CCK+cells. Last three columns are always reward Firing values except for the difference matrices.  

The file VRNeuronTask.m is the main script was used to run the analysis of firing rates, speed and theta coupling during 
the VR reality task. 

The file VRCCKGLM.m is the script was used extract for identified CCK+ cells the firing and theta power in defined
episodes and fit a GLM. 

VRShuffling.m is the script that was used to determine the statistical significance of different firing between 
trial types or sometimes neuron types. 

To rerun the script on our data, the data can be requested from the corresponding authors of the study. Further
instructions can be found in the scripts. 

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
