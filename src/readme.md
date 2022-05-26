# What this Program Does  
  The purpose of this prgram is to retireve an input file, specifically `EXPERIMENT_AME2016.dat`, and use the semi-empirical 
mass formula (SEMF) for calculating binding energies of neucleus' for all possible isotopes of known and unknown, stable and
unstable, elements. Binding energy is calculated by the SEMF method, which calculates the total binding energy of a neutron 
based on the liquid drop model. The SEMF has five terms, plus a sixth for goodness of fit later in the program. The first five 
terms relate to the volume of the liquid drop, the surface, asymmetry in distribution of protons and neutrons, coulomb neighbor
interactions, and pairing , respectively. The program uses linear algebra to solve a linear system including these terms and
their best fit parameters, solving a matrix equation in order to retrieve these five (or six) parameters. These parameters
multiply their corresponding linear terms (volume, surface, asymmetry, etc.) in order to produce the final SEMF. These 
parameters, and their errors, are printed to screen for the user.   
    After SEMF calculates all binding energies, the program writes the neutron count, proton count, calculated binding energy,
experimental binding energies (read for file), exprimental errors (read from file), and calculated uncertainty to a file named 
`results.dat`.   
    Another item of significant interest is the position of the neutron driplines and the region of most stable isotopes. In the
regions where the separation energy is positive, it costs energy to break apart a nucleus (remove neutrons or protons). Where as
when the separation energy is negative, we can observe neutrons and/or protons "dripping" off of the nucleus. This region in 
which the separation energy becomes negative is referred to as the "neutron-dripline". The region which have the highest 
separation energy is referred to as the "valley of stability" (regions where the binding energy is most negative, most stuck to 
nucleus). This program will also plot these neutron driplines and stable isotopes in `plots_analysis.ipynb` which can be 
accessed through "jupyter notebook", a python plotting package. The data for these plots can be found in `advanced_results.dat`,
which will be created and written by this program.   
    As an added attempt at decreasing error in this SEMF calculation, another term is introduced into the SEMF in order to 
increase the goodness of fit (chi-square), which is also printed in the jupyter notebook file. This extra (sixth) term is
a best guess term designed to decrease the deviation between experimental and calculated values. 

# How to Use
Please make sure all files are present in your current directory, refer to these files below in the contents section. With
the make file present, type "make" into your terminal and pres "Enter". This will compile all files and create an executable
named `nuclear_energies`. Then type "./nuclear_energies" into your terminal and press "Enter". You will then be prompted to
enter the name of your input file, this program is designed to run with `EXPERIMENT_AME2016.dat`. The program will then 
print the best fit parameters and corresponding uncertainties to screen, write binding energy results to `results.dat`,
and write the stable isotope and neutron dripline data to a file named `advanced_results.dat`. The user can now open the 
`plots_analysis.ipynb` with jupyter notebook. 


# Contents
This program contains `types.f90`, `read_write.f90`, `plots_analysis.ipynb`, `nuclear_model.f90`, `linear_algebra.f90`,
`makefile`, `main.f90`, `EXPERIMENT_AME2016.dat`.   
`read_write.f90` reads the input file, writes the binding energies to file, and writes the driplines and valley of stability
to file.  
`nuclear_model.f90` is responsible for constructing the alpha matrix and beta vector for use in finding the best fit
parameters, finding the most stable isotopes and neutron driplines, and printing to screen the best fit parameters and
corresponding uncertainties. Most importantly, this module contains the functions to calculate the SEMF and associated 
error. 
`linear_algebra.f90` performes the matrix operations and computational methods necessary to find the inverse of alpha
matrix and solve the linear system in order to retrieve the best fit parameters for use in the nuclear model module.   
`types.f90` contains the argument types, integers and reals for calculation.   
`plots_analysis.ipynb` is the file made to open in the jupyter notebook for plotting the results, click run all to see plots.
`makefile` allows the user to type "make" into terminal to compile all .f90 files and create the executable.  
`main.f90` houses the main calls to run the program, calling all write routines.  
`EXPERIMENT_AME2016.dat` houses the input data for the program to run the semi empirical methods from, hence the "semi".
