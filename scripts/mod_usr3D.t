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
    double precision                :: rhodens,rholight,kx,kz,pint,dlin,vextra,sigma
    logical                         :: first = .true.

    rhodens=10.0d0
    rholight=1.0d0
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

    if (first) then
       if (mype==0) then
          print *,'3D MHD KH'
          print *,'  --assuming y ranging from 0-1!'
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

    ! width of the smooth transition layer:
    dlin=0.05d0
    sigma=0.00125d0

    ! --- added --- Smooth density profile using tanh:
    w(ixG^S,rho_) = rholight + 0.5d0 * (rhodens - rholight) * &
         (dtanh((x(ixG^S,2) - 0.25d0) / dlin) - dtanh((x(ixG^S,2) - 0.75d0) / dlin))

    w(ixG^S,p_)=pint

    ! Initiating the magnetic field, uniform and pointing along the z-direction:
    w(ixG^S,b1_) = 0.0d0
    w(ixG^S,b2_) = 0.0d0
    w(ixG^S,b3_) = 4.0824829d0 ! Yields beta = 0.3 with p (thermal pressure) = 2.5

    ! --- added --- Smooth velocity profile (mom(1)) using tanh:
    w(ixG^S,mom(1)) = vextra + 0.05d0 * &
         (dtanh((x(ixG^S,2) - 0.25d0) / dlin) - dtanh((x(ixG^S,2) - 0.75d0) / dlin) - 1.0d0)

    w(ixG^S,mom(2)) = 0.001d0 * ( &
            dsin(kx * x(ixG^S,1) + 1.3d0 * x(ixG^S,3)) + &
            dsin(2.3d0 * kx * x(ixG^S,1) - 0.7d0 * x(ixG^S,3)) + &
            dsin(0.5d0 * kx * x(ixG^S,1) + 2.1d0 * x(ixG^S,3)) &
            ) * (dexp(-0.5d0*(x(ixG^S,2)-0.25d0)**2/sigma) + dexp(-0.5d0*(x(ixG^S,2)-0.75d0)**2/sigma))

    call mhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr

end module mod_usr
