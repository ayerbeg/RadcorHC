# RadcorHC
This code is intended to model the phase-space of a given kinematics, slice the phase-space in Es and Ep bin (for nomenclature, please refer to L.W. Mo and Y.S. Tsai, Rev. Mod. Phys 41, 205 (1969) ) The code produces and store a file with the data calculated. 
It works with the CLAS model, which this code wraps. With the Structure Functions extracted from the CLAS model, it calculates the He3 cross-section (the simple way 2p+1n). 

The code is very simple but at the same time quite versatile and easily to adapt and plug new models (but using the CLAS model options, which are handled by strucfunc.f).

For debugging, the code has the options to plot the structure functions for each nucleon, and the He3 cross-section, but only the las Es bin, which corresponds to the electron energy beam. For He3 there is the option to plot all the spectras XS vs nu for each Es bin, for a quasi-3D representation. 

All the options are handled by the options.ini file. 

# How To compile
- create a folder (this will be the SOURCE folder) -> enter into the new folder.
- clone the repository
- exit the SOURCE folder. Create another folder (this will be BUILD folder) -> enter into the new folder.
- run cmake -> cmake ../<SOURCE folder> where <SOURCE folder> is the name of the first folder you create.
- if there are not errors after the configuration, cmake should copy the files needed in the BUILD folder.
- USE the options.ini file in the BUILD folder.
- if the compilation was suscessful, run with ./a1nd2n

The options.ini is quite self-expenatory for each field, except for the CLAS model options, which you should follow the indications of the file G1 Code Overview.pdf (I will submit in a future commit)
