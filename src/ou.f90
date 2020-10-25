program ou
    implicit none
    ! Numerical parameters and system coefficients
    integer,  parameter :: WP   = kind(1.0D0)
    integer,  parameter :: Np   = 10000
    integer,  parameter :: Tmax = 1000
    real(WP), parameter :: tau  = 10.0_WP
    real(WP), parameter :: eta  = 0.01_WP
    real(WP), parameter :: dt   = 0.1_WP
    integer,  parameter :: Nt = int(Tmax/dt) + 1
    ! Output interval
    integer,  parameter :: Nskip = 10

    ! Array for solutions
    real(WP), dimension(Np) :: z, v
    real(WP), dimension(Nt) :: var
    ! loop variable
    integer :: i

    ! Initial values
    z(:) = 0.0_WP
    v(:) = 0.0_WP
    var(1) = variance(z)

    ! File for output
    open(unit=42, file='../output/variance_fortran.txt')

    do i = 2, Nt
        call langevin(z, v, tau, eta, dt)
        ! write time and variance every Nskip timesteps
        if (mod(i, Nskip) == 1) then
            write(42,*) (i-1)*dt, variance(z)
        end if
    end do

    ! Close output file
    close(42)

contains

    subroutine langevin(z, v, tau, eta, dt)
        ! Implements numerical solution of the SDE
        ! dz = v * dt
        ! dv = -1/tau * v * dt + eta dW
        ! by means of the Euler-Maruyama method
        ! Equivalent to Eqs. (1.9a) and (1.9b) in Gillespie (1996),
        ! Phys Rev E vol. 54 no. 2
        implicit none
        ! input variables
        real(WP), intent(in) :: tau, eta, dt
        ! Positions and velocities to be updated
        real(WP), intent(inout), dimension(:) :: z, v
        ! Number of particles
        integer  :: Np
        Np = size(z)
        ! Calculate next position and velocity
        z = z + v * dt
        v = v - (1.0_WP/tau) * v * dt + eta * dW(Np, dt)
    end subroutine

    function variance(z)
        implicit none
        ! input array
        real(WP), intent(in), dimension(:) :: z
        ! return value
        real(WP) :: variance
        ! Number of elements
        integer  :: Np
        Np = size(z)
        variance = sum((z - (sum(z)/Np))**2)/Np
    end function

    function boxmuller()
        ! Obtain two independent Gaussian random numbers with zero mean and unit variance
        ! Uses the polar form of the Box-Muller transform
        ! https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform#Polar_form
        implicit none
        ! Return variable
        real(WP), dimension(2) :: boxmuller
        ! local variables
        real(WP) :: u, v, s
        ! Draw random numbers until s is within accepted range
        s = 2
        do while ((s == 0) .or. (s >= 1))
            call random_number(u)
            call random_number(v)
            u = 2*u - 1
            v = 2*v - 1
            s = u**2 + v**2
        end do
        ! transform and return
        boxmuller = (/ u, v /) * sqrt(-2*log(s)/s)
    end function

    function dW(Np, dt)
        ! Returns an array of dimension(Np) containing
        ! independent Gaussian random numbers with zero mean and
        ! standard deviation sqrt(dt)
        implicit none
        integer,  intent(in)    :: Np
        real(WP), intent(in)    :: dt
        real(WP), dimension(Np) :: dW
        ! local variables
        real(WP), dimension(2)  :: tmp
        integer :: i
        ! Fill array with Gaussian random numbers with zero mean and unit variance
        ! (boxmuller always returns two numbers, so filling the array two at a time)
        do i = 1, int(Np/2)
            dW(2*i-1:2*i) = boxmuller()
        end do
        ! Add one last random number if array has odd length
        if (2*int(Np/2) < Np) then
            tmp = boxmuller()
            dW(Np) = tmp(1)
        endif
        ! Scale to correct standard deviation
        dW = dW * sqrt(dt)
    end function

end program
