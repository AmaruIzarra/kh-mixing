!=============================================================================
! AGILE user file: Kelvin-Helmholtz mixing layer, full 3D MHD
! Hillier & Arregui (2019, ApJ 885:101), Section 3.1, ported from our
! MPI-AMRVAC 2.5D reference setup (paper2d/mod_usr.t).
!
!   - single interface at y = 0: light (rho=1) below, dense (rho=10) above
!   - tanh profiles of width 0.003 in density and x-velocity
!   - velocity jump dV = sqrt(0.1) (= cool-gas sound speed), zero-momentum frame
!   - uniform pressure p = 1/gamma, uniform B along z with plasma beta = 0.3
!   - v_y seeded with random noise |v_y| <= 0.01 dV in every cell (paper)
!   - passive tracer: 0 in the light fluid, 1 in the dense fluid
!
! Differences from the MPI-AMRVAC version:
!   - .fpp file, written out in 3D (no ^L / ^S shorthand); AGILE is 3D only
!   - physics chosen in agile.par (phys = 'mhd'), activated with phys_activate
!   - TRACERS ARE STORED IN CONSERVED FORM (rho * tracer) and are NOT
!     converted by phys_to_conserved, so we set rho * prof directly
!   - the random seed also depends on z, so the noise is genuinely 3D
!=============================================================================
#:mute
#:include "mod_gpu_directives.fpp"
#:endmute

module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  ! physical parameters (code units: velocity in units of the light-gas
  ! sound speed, density in units of rho1)
  double precision, parameter :: rho1  = 1.0d0         ! light fluid, y < 0
  double precision, parameter :: rho2  = 10.0d0        ! dense fluid, y > 0
  double precision, parameter :: dv    = sqrt(0.1d0)   ! velocity jump V1 - V2 (paper: sqrt(0.1))
  double precision, parameter :: width = 0.003d0       ! tanh half-width
  double precision, parameter :: p0    = 0.6d0         ! = 1/gamma, gamma = 5/3
  double precision, parameter :: beta  = 0.3d0         ! plasma beta = 2 p0 / B^2
  double precision, parameter :: amp   = 0.01d0*dv     ! perturbation level

  ! read from &usr_list in the par file (only 'random' implemented here)
  character(len=20) :: seed_type = 'random'

contains

  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /usr_list/ seed_type

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    close(unitpar)
    end do
  end subroutine usr_params_read

  subroutine usr_init()

    call usr_params_read(par_files)
    call set_coordinate_system("Cartesian_3D")

    usr_init_one_grid => initonegrid_usr

    call phys_activate()

    if (trim(seed_type) /= 'random') &
         call mpistop("this AGILE setup only implements seed_type = 'random'")

    if (mype == 0) then
      print *, 'KH mixing layer, Hillier & Arregui (2019) setup, 3D MHD (AGILE)'
      print *, '  --density ratio rho2/rho1 :', rho2/rho1
      print *, '  --velocity jump dV        :', dv
      print *, '  --V1 (light side)         :', rho2/(rho1+rho2)*dv
      print *, '  --B_z                     :', sqrt(2.0d0*p0/beta)
      print *, '  --seed type               : ', trim(seed_type)
      print *, '  --variables nw            :', nw
    end if

  end subroutine usr_init

  !> Uniform pseudo-random number in [0,1) from position only (x, y, z),
  !> so the noise is the same on any number of processes/GPUs.
  elemental double precision function hash01(xx, yy, zz)
    double precision, intent(in) :: xx, yy, zz
    double precision             :: s
    s = sin(12.9898d0*(1.0d3*xx) + 78.233d0*(1.0d3*yy) &
          + 37.719d0*(1.0d3*zz)) * 43758.5453d0
    hash01 = s - floor(s)
  end function hash01

  subroutine initonegrid_usr(&
      ixGmin1, ixGmin2, ixGmin3, ixGmax1, ixGmax2, ixGmax3,&
      ixmin1,  ixmin2,  ixmin3,  ixmax1,  ixmax2,  ixmax3,&
      w, x)
    integer, intent(in)             ::&
      ixGmin1, ixGmin2, ixGmin3, ixGmax1, ixGmax2, ixGmax3,&
      ixmin1,  ixmin2,  ixmin3,  ixmax1,  ixmax2,  ixmax3
    double precision, intent(in)    :: x(&
      ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(&
      ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:nw)

    double precision :: prof(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3)
    double precision :: v1, b0, vshift

    ! Light-side velocity for zero total momentum:
    ! rho1*V1 + rho2*V2 = 0 with V1 - V2 = dV  ->  V1 = rho2/(rho1+rho2) dV
    v1 = rho2/(rho1+rho2)*dv

    ! Field strength from plasma beta (code units: p_mag = B^2/2)
    b0 = sqrt(2.0d0*p0/beta)

    ! iprob = 2: shift the whole flow by a constant velocity (Galilean test)
    select case(iprob)
    case(1)
      vshift = 0.0d0
    case(2)
      vshift = 1.0d0
    case default
      error stop "Unknown iprob"
    end select

    associate(&
      w_ => w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,:),&
      x_ => x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,:))

      ! smooth step: 0 deep in the light fluid, 1 deep in the dense fluid
      prof = 0.5d0*(1.0d0 + tanh(x_(:,:,:,2)/width))

      ! density and x-velocity follow the same tanh profile
      w_(:,:,:,rho_)   = rho1 + (rho2 - rho1)*prof
      w_(:,:,:,mom(1)) = v1 - dv*prof + vshift

      ! paper's seed: uniform random noise in v_y, |v_y| <= amp, every cell
      w_(:,:,:,mom(2)) = amp*(2.0d0*hash01(x_(:,:,:,1), x_(:,:,:,2), x_(:,:,:,3)) - 1.0d0)
      w_(:,:,:,mom(3)) = 0.0d0

      ! uniform pressure and uniform field along z -> equilibrium
      w_(:,:,:,p_)     = p0
      w_(:,:,:,mag(1)) = 0.0d0
      w_(:,:,:,mag(2)) = 0.0d0
      w_(:,:,:,mag(3)) = b0

      ! passive tracer, already in CONSERVED form (rho * tracer) in AGILE
      w_(:,:,:,tracer(1)) = w_(:,:,:,rho_)*prof

    end associate

    ! primitive (rho, v, p, B) -> conserved (rho, rho v, e, B)
    call phys_to_conserved(&
      ixGmin1, ixGmin2, ixGmin3, ixGmax1, ixGmax2, ixGmax3,&
      ixmin1,  ixmin2,  ixmin3,  ixmax1,  ixmax2,  ixmax3,&
      w, x)

  end subroutine initonegrid_usr

end module mod_usr
