! Program: nuclear_energies
! By: Louis Andre
!-----------------------------------------------------------------------------
! Description of program: Program takes experimental data and calculates binding 
! energies of isotopes for every proton count from one to 118. It then calculates
! the stable isotope and neutron dripline for this range of isotopes. 
! More information is in the readme file. 
!-----------------------------------------------------------------------------
program nuclear_energies
use types
use read_write, only : read_exp_data, write_predictions, write_nuclear_lines!, ... your advanced subroutines would go here
use nuclear_model, only : find_best_parameters!, ... your advanced subroutines would go here
implicit none

integer, allocatable :: n_protons(:), n_neutrons(:)
real(dp), allocatable :: exp_values(:), uncertainties(:), c_parameters(:), covariance(:,:)

! Read input data.
call read_exp_data(n_protons, n_neutrons, exp_values, uncertainties)
! Find best parameters for binding energy equation.
call find_best_parameters(n_protons, n_neutrons, exp_values, uncertainties, c_parameters, covariance)
! Write predictions of binding energy and associated error. 
call write_predictions(exp_values, uncertainties, c_parameters, covariance, n_protons, n_neutrons)
! Write data corresponding to stable isotopes and neutron dripline to file. 
! The all to the DRIPLINE and STABLE ISOTOPE subroutines 
call write_nuclear_lines(c_parameters, n_protons, n_neutrons)

end program nuclear_energies