**********************************************************************
*                        MANTA_TEC    (Ver.1.0.0)                    *
*          Scripts to calculate TEC from RINEX 2.11 files.           *
*                                                                    *
*          Institut de physique du globe de Paris (IPGP)             *
*          		Author: Fabio Manta    		          *
**********************************************************************
==== Contents of the directory ====

The structure of the working directory have to be organise as follow: 
	1) codes, contain all the script and function necessary to analyse the data.
	2) input, is the folder where to place the RINEX files to analyse, as well as the files necessary to correct the TEC for the satellite (file sp3) and receiver (file IONEX) interfrequency biases.
	3)output, it contains the .mat files resulting for the RINEX processing.
	4)plots; it contains the plots obtained by running the example provided using the script named "Test_routine.m".

===================================
Guided example:
Within the main folder, the user can find a PDF (Tutorial_Manta-TEC.pdf) file that contains a guided tutorial about MANTA-TEC with a practice example. The script "Test_routine.m" is an example that calls the function MANTA_TEC.m, to analyze some data recorded during the 2015 Mw 8.3 Illapel (Chile) earthquake (see input folder). This script returns together with the results (saved in the folder output), also a set of plots that provide the user with some example of data visualization.

===================================
Remarks:
  - Before using this toolkit, it is necessary to have preinstall some essential Matlab toolboxes:
	Image Processing Toolbox                              
	M_Map - mapping toolbox (Author: rich@eos.ubc.ca)         
	Mapping Toolbox                                       
	Optimization Toolbox                                 
	Signal Processing Toolbox       
                    
===================================



