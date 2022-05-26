!-----------------------------------------------------------------------
!Module: nuclear_model
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This module is responsible for the part of this program 
!! that deals with physical concepts. Here the calculation of the individual
!! linear terms within the binding energy calculation without parameters
!! (referred to in readme) is calculated. Using these linear terms, a matrix
!! is constructed consisting of the product of neighboring linear terms divided by 
!! the square of the uncertainties, then summed over all data points.
!! The beta vector is also constructed here, from the product of linear terms
!! with experimental values (one term at a time), divided by experimental
!! uncertainties, and summed over all data points in input file. 
!!
!! The matrix and vector evaluated here are used in the linear_algebra 
!! module in order to retrieve the parameters each term in the linear terms 
!! should be multiplied by in order to calculate binding energy. The 
!! linear_algebra section also calculates the inverse of this matrix
!! in order to find these parameters.
!!
!! Finally, the positions of the stable isotoped and neutron dripelines
!! are calculated in this module. These are explained in readme. 
!! 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! find_best_parameters
!! construct_alpha_beta
!! calculate_linear_termns
!! print_best_parameters
!! n_most_stable
!! neutron_drip_line
!!
!!----------------------------------------------------------------------
!! Included functions:
!!
!! volume_term
!! surface_term
!! asymmetry_term
!! coulomb_term
!! pairing_term
!! semi_empirical_mass
!! semi_empirical_error
!! delta
!! find_largest_n
!!
!!----------------------------------------------------------------------
module nuclear_model
use types
use linear_algebra, only : solve_linear_system
implicit none

private

public :: find_best_parameters, semi_empirical_mass, semi_empirical_error, n_most_stable, neutron_drip_line
contains


!-----------------------------------------------------------------------
!! Subroutine: find_best_parameters
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: Subroutine that calls for the construction of the matrix
!! alpha and the vector beta, calls for the solving of the matrix equation
!! (linear system), and finally calls for the printing of the best fit 
!! parameters (solution to the linear system).
!!----------------------------------------------------------------------
!! Input:
!!
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!-----------------------------------------------------------------------
!! Output:
!!
!! c_parameters     real        Array containing the semi-empirical mass formula parameters
!! covariance       real        Array containing the covariance matrix of the parameters
!-----------------------------------------------------------------------
subroutine find_best_parameters(n_protons, n_neutrons, exp_values, uncertainties, c_parameters, covariance)
    implicit none
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    real(dp), intent(in) :: exp_values(:), uncertainties(:)
    real(dp), intent(out), allocatable ::  c_parameters(:), covariance(:,:)

    ! Number of parameters (number of terms) in binding energy calcuation. 
    integer, parameter :: n_parameters = 6
    ! Initializing alpha matrix and beta vector with size of n_parameters
    real(dp) :: alpha(1:n_parameters,1:n_parameters), beta(1:n_parameters)
    ! If these are allocated with trash, they are cleaned out. 
    ! c_parameters and covariance were passed as arguments and need to be
    ! allocated. Deallocate them if they're allocated and then allocate them
    ! with the correct size
    if(allocated(c_parameters)) deallocate(c_parameters)
    if(allocated(covariance)) deallocate(covariance)

    ! Allocate the arrays which were just deallocated. covariance is the
    ! inverse of alpha, to be used in solve_linear_system.
    allocate(c_parameters(n_parameters))
    allocate(covariance(n_parameters, n_parameters))
    
    ! Alpha and beta are constructed.
    call construct_alpha_beta(n_protons, n_neutrons, exp_values, uncertainties, alpha, beta)

    ! The subroutine below (defined in the linear_algebra module) solves
    ! the matrix equation and returns the c_parameters and the
    ! inverse of alpha (the covariance matrix)
    call solve_linear_system(alpha,beta,c_parameters,covariance)

    ! Now just print the parameters (with it's uncertainties) to screen
    call print_best_parameters(c_parameters,covariance)
end subroutine find_best_parameters

!-----------------------------------------------------------------------
!! Subroutine: construct_alpha_beta
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This subroutine constructs the alpha matrix and beta matrix
!! discussed at the beginning of the module. 
!!----------------------------------------------------------------------
!! Input:
!!
!! n_protons        integer     Array containing the number of protons in each data point
!! n_neutrons       integer     Array containing the number of neutrons in each data point
!! exp_values       real        Array containing the binding energy in each data point
!! uncertainties    real        Array containing the statistical uncertainty in each data point
!-----------------------------------------------------------------------
!! Output:
!!
!! alpha            real        Array containing the alpha matrix
!! beta             real        Array containing the beta vector
!-----------------------------------------------------------------------
subroutine construct_alpha_beta(n_protons, n_neutrons, exp_values, uncertainties, alpha, beta)
    implicit none
    integer, intent(in) :: n_protons(:), n_neutrons(:)
    real(dp), intent(in) :: exp_values(:), uncertainties(:)
    real(dp), intent(out) :: alpha(:,:), beta(:)

    integer :: n_data, n_parameters, alpha_shape(1:2), i, j, k, beta_size
    real(dp):: linear_terms(1:size(beta))

    ! Check if the alpha array is a square matrix
    ! Also check that beta has the same number of elements as alpha has columns.
    alpha_shape = shape(alpha)
    beta_size = size(beta)

    if (alpha_shape(2) /= beta_size) then
        print*, 'Number of elements in beta do not equal number of columns in alpha'
        stop
    endif


    if (alpha_shape(1) /= alpha_shape(2)) then
        print*, 'The alpha matrix is not square'
        stop
    endif

    ! Declare variable for number of data points in file.
    n_data = size(uncertainties)
    ! Declare variable for number of terms in linear terms, also row and colums number. 
    n_parameters = alpha_shape(1)

    alpha = 0._dp
    beta = 0._dp
    ! Construct alpha with a do loop.
    do i=1,n_parameters
    ! "i" will be row identifier in matrix
        do j=1,n_parameters
    ! "j" will be column ientifier in matrix
            do k=1,n_data
    ! "k" will allow sum through all data points while one a specific "i" and "j".
                call calculate_linear_termns(n_protons(k), n_neutrons(k), linear_terms)
    ! the i'th and j'th terms of linear terms (neighboring terms) are multiplied
    ! and then divided by the experimental uncertainties for a specific k. 
    ! This is summed for one i and j while k runs through all data points. This is
    ! then repeated for the next j, until j has run to n_parameters. This is repeated for 
    ! each i, up to n_parameters. 
                alpha(i,j) = alpha(i,j) + ((linear_terms(i))*(linear_terms(j))/((uncertainties(k))**2 ))
            enddo
        enddo
    enddo
    ! Construct beta vector with a similar loop. 
        do i=1,n_parameters
    ! Find linear terms for the i'th proton and neutron number. 
        call calculate_linear_termns(n_protons(i), n_neutrons(i), linear_terms)

            do k=1,n_data
    ! Find linear terms for the k'th proton and neutron. 
                call calculate_linear_termns(n_protons(k), n_neutrons(k), linear_terms)
    ! First beta term is found by product of i'th linear term with k'th experimental
    ! value from input file, then divided by k'th uncertainty. This is summed through 
    ! all k values up to number of data points, then repeated for the next "i".
                beta(i) = beta(i) + ((linear_terms(i))*(exp_values(k)))/((uncertainties(k))**2) 
            enddo
        enddo



    
end subroutine construct_alpha_beta

!-----------------------------------------------------------------------
!! Subroutine: calculate_linear_termns
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This subroutine calls the individual functions for each 
!! term in the linear terms array. Once again, the linear terms array is
!! a list of all terms within binding energy equation, neglecting parameters. 
!!----------------------------------------------------------------------
!! Input:
!!
!! Z                integer     number of protons in an isotope
!! N                integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! linear_terms        real        Array containing the linear terms in the semi-empirical mass formula
!-----------------------------------------------------------------------
subroutine calculate_linear_termns(Z, N, linear_terms)
    implicit none
    integer, intent(in) :: Z, N
    real(dp), intent(out) :: linear_terms(:)

    ! We could write down all the formulas for each term here. However, in
    ! order to keep the code readable and easy to understand  we'll  separate
    ! them into different functions
    linear_terms(1) = volume_term(Z,N)
    linear_terms(2) = surface_term(Z,N)
    linear_terms(3) = asymmetry_term(Z,N)
    linear_terms(4) = coulomb_term(Z,N)
    linear_terms(5) = pairing_term(Z,N)
    linear_terms(6) = my_extra_term(Z,N)
end subroutine calculate_linear_termns

!-----------------------------------------------------------------------
!! function: volume_term
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: Calculated volume term of linear terms array.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        volume term
!-----------------------------------------------------------------------
real(dp) function my_extra_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N

    r = real((N - Z)**2 , kind=dp)/(real(  ( (Z+N)**(4._dp/3._dp) ),kind=dp))
    ! This term is similar to the asymettry term, but possibly dealing
    ! with asymmetry on the surface of the drop model's drop. 

end function my_extra_term

!-----------------------------------------------------------------------
!! function: volume_term
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: Calculated volume term of linear terms array.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        volume term
!-----------------------------------------------------------------------
real(dp) function volume_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N

    r = real(Z+N,kind=dp) 

end function volume_term

!-----------------------------------------------------------------------
!! function: surface_term
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: Calculates the surface term of the linear terms array. 
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        surface term
!-----------------------------------------------------------------------
real(dp) function surface_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N

    r = real(Z+N,kind=dp)**(2._dp/3._dp)

end function surface_term

!-----------------------------------------------------------------------
!! function: asymmetry_term
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: Calculates the asymmetry terms in the linear terms array.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        asymmetry term
!-----------------------------------------------------------------------
real(dp) function asymmetry_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N

    r = (real(N-Z, kind=dp)**2._dp)/real(N+Z,kind=dp)

end function asymmetry_term

!-----------------------------------------------------------------------
!! function: coulomb_term
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: Calculates the coulomb term in the linear terms array. 
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        coulomb term
!-----------------------------------------------------------------------
real(dp) function coulomb_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N

    r = real(Z*(Z - 1),kind=dp )/(real(Z+N,kind=dp)**(1._dp/3._dp))

end function coulomb_term

!-----------------------------------------------------------------------
!! function: pairing_term
!-----------------------------------------------------------------------
!! By:Louis Andre
!!
!! Description: Calculates the pairing term in the linear terms array. 
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in a nucleus
!! N            integer     number of neutrons in a nucleus
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        pairing term
!-----------------------------------------------------------------------
real(dp) function pairing_term(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N

    r = ((real(Z+N, kind=dp)**(-3._dp/4._dp)))*(delta(Z, N))
    
end function pairing_term

!-----------------------------------------------------------------------
!! Subroutine: print_best_parameters
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: Prints the best fit parameters to screen, these are the 
!! solution to the matrix equation. Each parameter is multiplied by its
!! corresponding linear term to calculate binding energy later.
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters     real        Array containing the best fit parameters
!! covariance       real        Array containing covariance matrix
!-----------------------------------------------------------------------
subroutine print_best_parameters(c_parameters, covariance)
    implicit none
    real(dp), intent(in) :: c_parameters(:), covariance(:,:)
    
    print *,
    print *, 'Best fit values:              value                 uncertainty'
    print 1, ' Volume parameter:   ', c_parameters(1),           sqrt(covariance(1,1))        
    print 1, ' Surface parameter:  ', c_parameters(2),           sqrt(covariance(2,2))
    print 1, ' Asymmetry parameter:', c_parameters(3),           sqrt(covariance(3,3))
    print 1, ' Coulomb parameter:  ', c_parameters(4),           sqrt(covariance(4,4))
    print 1, ' Pairing term:       ', c_parameters(5),           sqrt(covariance(5,5))
    print 1, ' My Extra Term:      ', c_parameters(6),           sqrt(covariance(6,6))

1 format(a,f15.8,e28.16)
end subroutine print_best_parameters



!-----------------------------------------------------------------------
!! function: semi_empirical_mass
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This subroutine uses the parameters found in the matrix 
!! equation and linear terms calculated in the functions above to calculate
!! the binding energy for a certain proton and neutron number. 
!!
!! First term in linear terms multiplies first term in parameters, second 
!! terms multiply, etc. 
!!----------------------------------------------------------------------
!! Input:
!!
!! c    real        Array containing the parameters of the semi-empirical mass formula
!! Z    integer     number of protons in an isotope
!! N    integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! r    real        Binding energy
!-----------------------------------------------------------------------
real(dp) function semi_empirical_mass(c, Z, N) result(r)
    implicit none
    real(dp), intent(in) :: c(:)
    integer, intent(in) :: Z, N
    real(dp), allocatable :: linear_terms(:)
    integer :: c_size, i
    ! Initialize variable with size equal to linear terms/parameters size. 
     c_size = size(c)
    ! Allocate array with length equal to number of terms in binding energy equation. 
    allocate(linear_terms(1:c_size))
    ! Calculate linear terms for a given proton and neutron amount. 
    call calculate_linear_termns(Z, N, linear_terms(:)) 
    ! Sum the product of corresponding terms in parameters and linear terms. 
    r = 0
    do i=1,c_size 
        r = r + c(i)*linear_terms(i)
    enddo


end function semi_empirical_mass

!-----------------------------------------------------------------------
!! function: semi_empirical_error
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This subroutine calculates the theoretical error in
!! the calculation of the binding energy in this program. 
!!----------------------------------------------------------------------
!! Input:
!!
!! covariance   real        2D array containing the parameters' covariance matrix
!! Z            integer     number of protons in an isotope
!! N            integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        statistical uncertainty in the binding energy
!-----------------------------------------------------------------------
real(dp) function semi_empirical_error(covariance, Z, N) result(r)
    implicit none
    real(dp), intent(in) :: covariance(:,:)
    integer, intent(in) :: Z, N
    real(dp), allocatable :: linear_terms(:), n_parameters(:)
    integer :: i, j, n_terms
    real(dp) :: s
    ! Allocate array holding row and column number of alpha and covariance matrix. 
    allocate(n_parameters(1:2))
    n_parameters = shape(covariance) 
    ! Rows of covariance will be equal to number of terms in binding energy. 
    n_terms = n_parameters(1)
    ! Allocate array to hold linear terms for a certain proton and neutron amount.
    allocate(linear_terms(1:n_terms))
    ! Fill array of linear terms for certain proton and neutron amount. 
    call calculate_linear_termns(Z, N, linear_terms(:))
    ! Begin do loop to calculate error for each binding energy.
    s = 0
    do i=1,n_terms
    ! First term in linear terms array will be the i'th. "i" is also row number in 
    ! covariance matrix. 
        do j=1,n_terms
    ! Second term in linear terms array will be the j'th. "j" is also column number
    ! in covariance matrix. 
            s = s + (linear_terms(i))*(linear_terms(j))*(covariance(i,j))
        enddo
    enddo
    ! The square root of the value calculated above is the error we are looking for. 
    r = sqrt(s)

end function semi_empirical_error

!-----------------------------------------------------------------------
!! function: delta
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: Function that returns -1, 0, or 1 depending on whether
!! or not the parity of the number of neutrons and protons match.
!!----------------------------------------------------------------------
!! Input:
!!
!! Z            integer     number of protons in an isotope
!! N            integer     number of neutrons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        resuls to parity delta function (both odd, both even, mis-matched)
!-----------------------------------------------------------------------
real(dp) function delta(Z, N) result(r)
    implicit none
    integer, intent(in) :: Z, N

    if ( modulo(Z,2) == 0 .and. modulo(N,2) == 0 ) then
    ! If proton and neutron numbers are both even, return a 1 for delta. 
        r = 1
    else if (modulo(Z,2) == 1 .and. modulo(N,2) == 1) then
    ! If proton and neutron numbers are both odd, return a -1 for delta. 
        r = -1
    else
    ! If parity of neutron and proton numbers dont match, return a zero for delta.
        r = 0
    endif


end function delta

!-----------------------------------------------------------------------
!! Subroutine: most_stable
!-----------------------------------------------------------------------
!! by: Louis Andre
!!
!! Description: This subroutine finds the neutron number that produces
!! the lowest binding energy per nucleons (z+n) for a given proton number. 
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters        real      1D array consisting of parameters to multiple each term in linear term to get BE.
!-----------------------------------------------------------------------
!! Output:
!!
!! n_stable       integer        1D array containing neutron number associated with lowest BE/A for a particular proton numebr
!-----------------------------------------------------------------------
subroutine n_most_stable(c_parameters, n_stable)
    implicit none
    real(dp), intent(in) :: c_parameters(:)
    integer, intent(out), allocatable :: n_stable(:)
    integer :: n, z, index, m, n_max, p_max
    real(dp) :: b, be_theory
    allocate(n_stable(1:118))
! p_max is largest number of protons seen in input file.
    p_max = size(n_stable)
z=1
index = 1
! Find what the highest amount of neutrons for an individual proton number
! will be, by calling a function in this module (see bottom of module).
n_max = find_largest_n(c_parameters, p_max  )

! Start loop to find neutron number corresponding to lowest binding energy per nucleon. 
! Outter loop will loop over proton amount, stating at one and ending at highest
! number of protons achieved. This is 118 for EXPERIMENT_AME2016.dat input file. 
do z=1,p_max
! Variable "b" to store binding energy per nucleon, starting with first binding 
! energy per nucleon value (corresponding to protons number = neutrons number =1).
    b = ( (semi_empirical_mass(c_parameters, z, z))/(real(z+z,kind=dp)) )
! Loop over all neutron amounts for z value defined by outter loop. 
    do n=1,n_max
! Calculate binding energy per nucleon (BE/A) for the outter loop's "z" and 
! inner loop's "n".
    be_theory = ( (semi_empirical_mass(c_parameters, z, n))/(real(n+z,kind=dp)) )
        if (be_theory < b) then
! If this second BE/A is smaller than the outside loop's BE/A, we store the second BE/A
! as the new "b" value. 
            b = be_theory
! We also store the current neutron number (that achieved this low BE/A) in "m".
            m=n   
        endif
! This is repeated until all neutron numbers have been checked for a single
! proton amount. 
    enddo
! Here we finished with a certain proton amount, and now will store the neutron 
! amount associated with lowest BE/A into an array "n_stable".
! The index of the array starts at one and increases by one everytime we enter a value. 
! This is done to ensure values arent overwritten by new values, on the next iteration. 
    n_stable(index)=m
    index= index +1    
enddo

end subroutine n_most_stable

!-----------------------------------------------------------------------
!! Subroutine: most_stable
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This subroutine finds the largest neutron amount for which 
!! the neutron separation energy is positive, for each proton amount. 
!! Neutron separation energy is defined here as difference between one 
!! binding energy (for a certain z and n), and the previous binding energy
!! (for the same z, but one less n). 
!! BE(Z,N-1) - BE(Z,N) > 0
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters   real        1D array containing the parameters
!! p_max          integer     largest number of protons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! n_dripline    integer        1D array containing the neutron amount associated with lowest separation energy for each proton amount
!------------------------------------------------------------------------
subroutine neutron_drip_line(c_parameters, p_max, n_dripline)
    implicit none
    integer, intent(in) :: p_max
    real(dp), intent(in) :: c_parameters(:)
    integer, intent(out), allocatable :: n_dripline(:)
    integer :: z, n, m, n_max
    real(dp) :: term_1, term_2, s
! Allocate array to house these largest neutron numbers. 
    allocate(n_dripline(1:p_max))

! Start loop to retrieve this largest neutron amount for positive separation energy.
do z=1,p_max
    n=1
! Calculate term_1 
    term_1 =  semi_empirical_mass(c_parameters, z, n)
! Calculate term_2, for an "n" one larger.  
    term_2 = semi_empirical_mass(c_parameters, z, n+1)
    do while (term_1 - term_2 > 0)
! While the difference is positive, perform these actions, then check to see if 
! they are still positive for the next pair of neutron numbers in this 
! specific isotope. 
        term_1 =  semi_empirical_mass(c_parameters, z, n)
        m=n+1
        term_2 = semi_empirical_mass(c_parameters, z, m)
        n=n+1 
    enddo
! Once this loop exits, one more than the last value of "m" we stopped on will be the largest 
! neutron number achieved for this proton amount. Now we increase "z" by one and
! repeat the process, until we reach a z=118=p_max.
    n_dripline(z) = m+1
!
enddo

end subroutine neutron_drip_line

!-----------------------------------------------------------------------
!! function: find_largest_n
!-----------------------------------------------------------------------
!! By: Louis Andre
!!
!! Description: This function takes the largest proton number in the 
!! isotopes in the input file and impliments a "do while" loop to 
!! determine the largest neutron number for which the difference
!! in BE of two neighboring neutron amounts becomes negative. This
!! number is then used to determine how large the loops over "n" should
!! be in the dripine and stable isotope subroutines, in order to not neglecting
!! any nautron amounts that may be useful. 
!!----------------------------------------------------------------------
!! Input:
!!
!! c_parameters   real        1D array containing the parameters
!! p_max          integer     largest number of protons in an isotope
!-----------------------------------------------------------------------
!! Output:
!!
!! r            real        neutron number for largest isotope of 118 proton element
!-----------------------------------------------------------------------
integer function find_largest_n(c_parameters, p_max) result(r)
    implicit none
    real(dp), intent(in) :: c_parameters(:)
    integer, intent(in) :: p_max
    integer :: p, n, m
    real(dp) :: term_2, term_1

! P is the proton number in the isotope. 
    p=p_max
    n=1
! Calculate term_1 as before
    term_1 =  semi_empirical_mass(c_parameters, p, n)
! Calculate term_2 as before, for an "n" one larger.
    m=n+1
    term_2 = semi_empirical_mass(c_parameters, p, m)
    do while (term_1 - term_2 > 0)
! While the difference is positive, perfore these actions, then check to see if 
! they are still positive for the next pair of neutron numbers in this 
! specific isotope. 
        term_1 =  semi_empirical_mass(c_parameters, p, n)
        m = n+1
        term_2 = semi_empirical_mass(c_parameters, p, m)
        n=n+1
    enddo
! Once this loop exits, the last value of "n" we stopped on will be the largest 
! nautron number achieved in this experiment. And now we know how far to loop over 
! neutron numbers for each isotope, without leaving out any data points. 
    r = m

end function find_largest_n


end module nuclear_model

