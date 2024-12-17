These Matlab codes can be used to generate the figures from 
Random Evolutionary Dynamics in Predator-Prey Systems Yields Large, Clustered Ecosystems.
Authors: Hamster, Schaap, Van Heijster and Dijksman.

The main difference between the codes for Fig. 2 and the other figures is that for Fig. 2, all timesteps are saved, which
is only feasible for short-time simulations. 

-Instruction for the files in Fig2. The script ScriptGLVexp.m is an implementation of the speciating GLV model and plots several of the observables. Just one additional function is needed, the script LV.m which is an implementation of equations (2) and (3). Note that the parameter Tend cannot be large, as in this script all data is stored. 

-Instruction for the files in Fig3-14. In comparison to the script in Fig2, in this script only the biomasses at speciation times are stored, allowing for larger values of Tend. The script ScriptGLVexpnoX.m contains the main code that can be run, followed by pieces of code that can be used to reproduce the analysis for figures 3-14. Script LV.m is an implementation of equation (2) and (3). PermMat.m is part of the GCM algorithm in GCModulMax1.m.


The files GCModulMax1.m and PermMat.m for the modularity algorithm are downloaded from https://nl.mathworks.com/matlabcentral/fileexchange/45867-community-detection-toolbox


Codes were tested with Matlab 2023b on Ubuntu 22.04 LTS
