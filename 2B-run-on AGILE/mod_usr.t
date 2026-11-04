module mod_usr
  use mod_mhd
  use mod_physics

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => initonegrid_usr

    call set_coordinate_system("Cartesian_3D")
    call mhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    integer, intent(in)             :: ixG^L, ix^L
    double precision, intent(in)    :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision                :: rhodens,rholight,kx,kz,pint,dlin,vextra,sigma,delta_v
    logical                         :: first = .true.

    rholight=1.0d0
    rhodens=10.0d0
    ! kx=2 pi
    kx=4.0d0*dpi
    {^IFTHREED
    ! kz=8 pi
    kz=32d0*atan(one)
    }

    select case(iprob)
    case(1)
       vextra=0.0d0
    case(2)
       vextra=10.0d0
    case(3)
       vextra=100.0d0
    endselect

    delta_v = 0.1d0 ! Velocity jump magnitude from paper

    if (first) then
       if (mype==0) then
          print *,'3D MHD KH'
          print *,'  --density ratio:',rhodens/rholight
          print *,'  --kx:',kx
          {^IFTHREED
          print *,'  --kz:',kz
          }
          print *,'  --vextra:',vextra
       end if
       first=.false.
    end if

    ! pressure at interface
    pint=2.5d0

    ! width parameter of the smooth transition layer from image formula:
    dlin=0.003d0
    sigma=0.00125d0

    ! --- Density formula from image: rho = rho1 + 0.5*(rho2 - rho1)*(tanh(y / 0.003) + 1) ---
    w(ixG^S,rho_) = rholight + 0.5d0 * (rhodens - rholight) * &
         (dtanh(x(ixG^S,2) / dlin) + 1.0d0)

    w(ixG^S,p_)=pint

    ! Initiating the magnetic field, uniform and pointing along the z-direction:
    w(ixG^S,b1_) = 0.0d0
    w(ixG^S,b2_) = 0.0d0
    w(ixG^S,b3_) = 4.0824829d0 ! Yields beta = 0.3 with p (thermal pressure) = 2.5

    ! --- Velocity formula from image: Vx = V1 + 0.5*delta_v*(tanh(y / 0.003) + 1) ---
    w(ixG^S,mom(1)) = vextra + 0.5d0 * delta_v * &
         (dtanh(x(ixG^S,2) / dlin) + 1.0d0)

    ! Velocity perturbation (y-direction):
    w(ixG^S,mom(2)) = 0.001d0 * ( &
            dsin(kx * x(ixG^S,1) + 1.3d0 * x(ixG^S,3)) + &
            dsin(2.3d0 * kx * x(ixG^S,1) - 0.7d0 * x(ixG^S,3)) + &
            dsin(0.5d0 * kx * x(ixG^S,1) + 2.1d0 * x(ixG^S,3)) &
            ) * dexp(-0.5d0 * (x(ixG^S,2))**2 / sigma)

    call mhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr

end module mod_usr