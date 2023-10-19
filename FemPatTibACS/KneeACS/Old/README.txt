kneeACSpackage.m is the umbrella function for the automatic knee ACS algorithm.

Prior to running the kneeACSpackage MATLAB script you should have a femur *.iv 
file, a tibia *.iv file and point located on the anterior half of the tibia.

to run the MATLAB script perform the following steps.

1. make sure you set your current directory to the directory containing all of
the knee ACS algorithm scripts or set the knee ACS algorithm directory as on
of your paths.

2. type the following into the command window:

[fACS tACS] = kneeACSpackage(pathname, femurivfile, tibiaivfile, anteriortibiapt);

pathname is string variable pointing to the directory where the *.iv files are
located. this will also be the directory where the *ACS.txt files will be saved

femurivfile and tibiaivfile are string variables pointing to the 3-D femur and
tibia model file names.

anteriortibiapt is a numerical variable with the 3-D location of a point located
on the anterior half of the tibia

fACS and tACS are numerical 4x4 matrices containing the ACS axes and origin