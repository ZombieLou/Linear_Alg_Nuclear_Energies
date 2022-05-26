!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This modue is responsible for reading the input file of 
!! experimental data to be used. The module also writes the calculated 
!! binding energies and corresponding errors to a file named 'results.dat'.
!! In addition to this, the module writes the neutron dripine and stable 
!! isotope data to a file named 'advanced_results.dat'.
!! 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_exp_data
!! write_predictions
!! write_nuclear_lines
!!----------------------------------------------------------------------
module read_write

use types
use nuclear_model, only : semi_empirical_mass, semi_empirical_error, n_most_stable, neutron_drip_line

implicit none

private
public :: read_exp_data, write_predictions, write_nuclear_lines

contains

!-----------------------------------------------------------------------
!! Subroutine: read_exp_data
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This subroutine reads the input file specified by user, 
!! works with a file named EXPERIMENT_AME2016.dat. Any file will work as 
!! long as it has the correct values in the ocrrect columns. What the program
!! reads in each column can be changed in this subroutine. 
!!----------------------------------------------------------------------
!! Output:
!!
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!-----------------------------------------------------------------------
subroutine read_exp_data(n_protons, n_neutrons, exp_values, uncertainties)
    implicit none
    integer, intent(out), allocatable :: n_protons(:), n_neutrons(:)
    real(dp), intent(out), allocatable :: exp_values(:), uncertainties(:)

    character(len=128) :: filename, string_trash
    logical :: file_exists
    integer :: file_unit, number_data_points, integer_trash1, i, integer_trash2
    real(dp) :: real_trash

    print *, 'This program will take information given by EXPERIMENT_AME2016.dat'
    print*, 'It takes experimental data and performs a semi-empirical calculation'
    print*, 'of the binding energies of a host of isotopes. It then uses these data'
    print*, 'to calculate the neutron dripline and the stable isotopes.'

    print *, 'Please provide the file name with the experimental data now:'
    read(*, '(a)') filename

    
    if(allocated(n_protons)) deallocate(n_protons)
    if(allocated(n_neutrons)) deallocate(n_neutrons)
    if(allocated(exp_values)) deallocate(exp_values)
    if(allocated(uncertainties)) deallocate(uncertainties)
    

    inquire(file=trim(filename), exist=file_exists)

    if (file_exists) then
        open(unit = 1, file=filename)
        ! First row, first column reads the number of rows (data points) in the input file. 
        read(1,*) number_data_points

        allocate(n_protons(1:number_data_points))
        allocate(n_neutrons(1:number_data_points))
        allocate(exp_values(1:number_data_points))
        allocate(uncertainties(1:number_data_points))
        ! Skip the two header rows. 
    do i=1,2
        read(1,*)
    enddo
    ! Loop for readig each row
        do i=1,number_data_points
    ! Data in colums we do not need is stored in a variable and never used. 
        read(1,*) integer_trash1, string_trash, integer_trash2, n_neutrons(i), n_protons(i), exp_values(i), real_trash, &
        & uncertainties(i)
        enddo 
        close(1)
    else
        ! If file does not exist, print error and abort program. 
        print*, 'The file named ',trim(filename),' could not be found.'
        stop
    endif
end subroutine read_exp_data

!-----------------------------------------------------------------------
!! Subroutine: write_predictions
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Describe what the subroutine does
!!----------------------------------------------------------------------
!! Input:
!!
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!! c_parameters     real        Array containing the parameters of the semi-empirical mass formula
!! covariance       real        Array containing the elements of the covariance matrix
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!-----------------------------------------------------------------------
subroutine write_predictions(exp_values, uncertainties, c_parameters, covariance, n_protons, n_neutrons)
    implicit none
    real(dp), intent(in) :: exp_values(:), uncertainties(:), c_parameters(:), covariance(:,:)
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    character(len=*), parameter :: file_name = 'results.dat'
    integer :: unit1, i, length
    real(dp) :: theoretical_error, theoretical_BE
    length = size(n_neutrons)

    open(newunit=unit1,file=file_name)
    write(unit1,28) 'Protons', 'Neutrons', 'Experimental BE', 'Experimental Error', 'Theory BE', 'Theory Error'
    28 format(5x,a,5x,a,5x,a,14x,a,12x,a,19x,a)
    do i=1,length
    ! Call subroutines to calculate wanted values and use do loop to print them row by row. 
        
    theoretical_BE = semi_empirical_mass(c_parameters, n_protons(i), n_neutrons(i))
    theoretical_error = semi_empirical_error(covariance, n_protons(i), n_neutrons(i))
    write(unit1,*) n_protons(i), n_neutrons(i), exp_values(i), uncertainties(i), theoretical_BE, theoretical_error

    enddo
    close(unit1)

    print *,
    print *, 'Theoretical binding energies and their errors were written in ', file_name
end subroutine write_predictions
!-----------------------------------------------------------------------
!! Subroutine: write_predictions
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Describe what the subroutine does
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters     real        Array containing the parameters of the semi-empirical mass formula
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!-----------------------------------------------------------------------
subroutine write_nuclear_lines(c_parameters, n_protons, n_neutrons)
    implicit none
    real(dp), intent(in) :: c_parameters(:)
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    character(len=*), parameter :: file_name = 'advanced_results.dat'
    integer :: unit2, i, p_max
    integer, allocatable :: n_stable(:), n_drip(:), z_values(:)
    
    allocate(n_stable(1:118))
    allocate(n_drip(1:118))
    allocate(z_values(1:118))


    ! Retrieve array of the neutron amounts associated with stable isotopes. 
    call n_most_stable(c_parameters, n_stable)
    p_max = size(n_stable)
    ! Call subroutine to calculate largest neutron count, per z value, for
    ! which the separation energy is positive. 
    call neutron_drip_line(c_parameters, p_max, n_drip)
    
    ! Do loop to fill array representing all possible proton amounts in the isotopes. 
    do i=1,p_max
        z_values(i) = i
    enddo
    ! Open file to write to. 
    open(newunit=unit2,file=file_name)
    write(unit2,2) 'Z Values', 'N Stable', 'N Drip'
    2 format(4x,a,4x,a,6x,a)

    ! Loop to write Protons, neutron of Stable Isotopes, and Neutrons of Dripline to file.
    do i=1,p_max      
    write(unit2,*) z_values(i), n_stable(i), n_drip(i)

    enddo
    close(unit2)
    print *,
    print *, 'Stable isotopes and neutron dripline was written to ', file_name
end subroutine write_nuclear_lines
    
end module read_write