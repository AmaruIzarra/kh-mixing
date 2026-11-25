!=============================================================================
! Kelvin-Helmholtz mixing layer, 2.5D MHD
! Setup follows Hillier & Arregui (2019, ApJ 885:101), Section 3.1:
!   - single interface at y = 0, light (rho=1) below, dense (rho=10) above
!   - tanh profiles of width 0.003 in density and x-velocity
!   - velocity jump dV = 0.1, in the "linear" zero-momentum frame
!   - uniform pressure p = 1/gamma, uniform B_z with plasma beta = 0.3
!   - v_y seeded with random noise at the level of 0.01 dV (paper's choice)
!   - passive tracer: 0 in the light fluid, 1 in the dense fluid
!
! The seed is chosen in the par file:
!   &usr_list
!     seed_type = 'random'   ! paper: noise in every cell, |v_y| <= 0.01 dV (default)
!     seed_type = 'modes'    ! our earlier test: 10 sine modes near the interface
!   /
! Compile with:  $AMRVAC_DIR/setup.pl -d=2 -v=3
!=============================================================================
module mod_usr
  use mod_mhd
  use mod_eos, only: eos

  implicit none

  ! physical parameters (code units: velocity in units of the light-gas
  ! sound speed, density in units of rho1)
  double precision, parameter :: rho1  = 1.0d0      ! density of the light gas, y < 0
  double precision, parameter :: rho2  = 10.0d0     ! density of the dense gas, y > 0
  double precision, parameter :: dv    = 0.1d0      ! velocity jump V1 - V2 (relative velocity of the two layers)
  double precision, parameter :: width = 0.003d0    ! scale length for tanh transition
  double precision, parameter :: p0    = 0.6d0      ! p0 = 1/gamma, gamma = 5/3 (Gas pressure, uniform everywhere)
  double precision, parameter :: beta  = 0.3d0      ! plasma beta = 2 p0 / B^2 (Gas pressure divided by magnetic pressure, used to compute strength of magnetic field)
  double precision, parameter :: amp   = 0.01d0*dv  ! perturbation level (amplitude of the kick in the y-direction)

  ! parameters for the 'modes' seed only other wise random noise used (like in paper)
  integer,          parameter :: nmodes   = 10      ! number of x-modes
  double precision, parameter :: pert_len = 0.05d0  ! y-localization of v_y

  ! read from &usr_list in the par file
  character(len=20) :: seed_type = 'random'

contains
  ! reads the settings from the par file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ seed_type

    do n = 1, size(files)
      open(unitpar, file=trim(files(n)), status="old")
      read(unitpar, usr_list, end=111)
111   close(unitpar)
    end do
  end subroutine usr_params_read

  ! 
  subroutine usr_init()
    integer :: iw

    usr_init_one_grid => initonegrid_usr

    call usr_params_read(par_files)
    call set_coordinate_system("Cartesian_2.5D")   ! 2D grid, 3 vector components (B_z!)
    call mhd_activate()

    if (seed_type /= 'random' .and. seed_type /= 'modes') &
         call mpistop("seed_type must be 'random' or 'modes'")

    ! print the setup once, from the first MPI process only
    if (mype == 0) then
      print *, 'Variables evolved (nwflux =', nwflux, ', nw =', nw, '):'
      do iw = 1, nw
        print *, '   ', iw, trim(cons_wnames(iw))
      end do
      print *, 'KH mixing layer, Hillier & Arregui (2019) setup'
      print *, '  --density ratio rho2/rho1 :', rho2/rho1
      print *, '  --velocity jump dV        :', dv
      print *, '  --V1 (light side)         :', rho2/(rho1+rho2)*dv
      print *, '  --B_z                     :', dsqrt(2.0d0*p0/beta)
      print *, '  --seed type               : ', trim(seed_type)
    end if

  end subroutine usr_init

  !> Uniform pseudo-random numbers in [0,1) computed from position only.
  !> Unlike Fortran's random_number, the value at a given point does not depend
  !> on which block or MPI process owns it, so a run gives the same initial
  !> noise on 1, 4 or 16 cores (needed for fair timing/scaling comparisons).
  !> (It does change with resolution, as grid-scale noise must.)
  elemental double precision function hash01(xx, yy)
    double precision, intent(in) :: xx, yy
    double precision             :: s
    s = dsin(12.9898d0*(1.0d3*xx) + 78.233d0*(1.0d3*yy)) * 43758.5453d0
    hash01 = s - floor(s)
  end function hash01

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: prof(ixG^S), vy(ixG^S)
    double precision :: v1, b0, vshift, lx, phase
    integer          :: m

    ! Light-side velocity chosen so total momentum is zero:
    ! rho1*V1 + rho2*V2 = 0 with V1 - V2 = dV  ->  V1 = rho2/(rho1+rho2) dV
    v1 = rho2/(rho1+rho2)*dv

    ! Field strength from plasma beta (AMRVAC units: p_mag = B^2/2)
    b0 = dsqrt(2.0d0*p0/beta)

    ! iprob = 2 shifts the whole flow by a constant velocity,
    ! for testing Galilean invariance (assignment item 2)
    vshift = 0.0d0
    if (iprob == 2) vshift = 1.0d0

    ! smooth step: 0 deep in the light fluid, 1 deep in the dense fluid
    prof(ixG^S) = 0.5d0*(1.0d0 + dtanh(x(ixG^S,2)/width))

    ! density and x-velocity follow the same tanh profile
    w(ixG^S,rho_)   = rho1 + (rho2 - rho1)*prof(ixG^S)
    w(ixG^S,mom(1)) = v1 - dv*prof(ixG^S) + vshift

    select case (trim(seed_type))
    case ('random')
      ! Paper (Sec. 3.1): "seeded with a random noise perturbation in v_y at
      ! the level of 0.01 dV". Uniform noise in [-amp, amp] in every cell.
      w(ixG^S,mom(2)) = amp*(2.0d0*hash01(x(ixG^S,1), x(ixG^S,2)) - 1.0d0)

    case ('modes')
      ! Earlier test seed: x-periodic sine modes with fixed phases,
      ! localized around the interface.
      lx = xprobmax1 - xprobmin1
      vy(ixG^S) = 0.0d0
      do m = 1, nmodes
        phase = 2.0d0*dpi*mod(dble(m)*0.6180339887d0, 1.0d0)
        vy(ixG^S) = vy(ixG^S) + &
             dsin(2.0d0*dpi*dble(m)*(x(ixG^S,1) - xprobmin1)/lx + phase)
      end do
      w(ixG^S,mom(2)) = amp*vy(ixG^S)/dsqrt(0.5d0*dble(nmodes)) * &
           dexp(-(x(ixG^S,2)/pert_len)**2)
    end select

    w(ixG^S,mom(3)) = 0.0d0

    ! uniform pressure and uniform out-of-plane field -> equilibrium
    w(ixG^S,p_)     = p0
    w(ixG^S,mag(1)) = 0.0d0
    w(ixG^S,mag(2)) = 0.0d0
    w(ixG^S,mag(3)) = b0

    ! passive tracer marks which fluid is where
    w(ixG^S,tracer(1)) = prof(ixG^S)

    ! everything above is primitive (rho, v, p, B); convert to
    ! conserved (rho, rho*v, total energy, B) which AMRVAC evolves
    call eos%to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr

end module mod_usr
