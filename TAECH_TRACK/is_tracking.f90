module is_tracking_mod
  use TAECH_IO
  use TAECH_LIB
  implicit none
  private

  public :: is_tracking
  public :: resonant_velocity
  public :: resonant_density
contains
    function is_tracking(fs0, M, s)
    implicit none
    integer, intent(in) :: M
    complex(kind=8), intent(in) :: fs0(-M:M)
    real(kind=8), intent(in) :: s(-M:M)
    logical :: is_tracking
    real(kind=8) :: lres, hres
    lres = v0
    hres = v1
    if (abs(resonant_velocity(fs0, M, s)-(lres+hres)/2.0D0)<0.05D0) then
       is_tracking = .True.
    else
       is_tracking = .False.
    end if
    return
  end function is_tracking

  function resonant_density(fs0, M, s)
    implicit none
    integer, intent(in) :: M
    complex(kind=8), intent(in) :: fs0(-M:M)
    real(kind=8), intent(in) :: s(-M:M)
    real(kind=8) :: resonant_density
    complex(kind=8) :: fv0(-M:M)
    real(kind=8) :: v(-M:M)
    real(kind=8) :: dv 
    real(kind=8) :: lres, hres
    integer :: i
    lres = v0
    hres = v1
    dv = pi/s(M)
    forall(i=-M:M) v(i) = i*dv 
    call IFFT1D(fv0, M, fs0, s)
    resonant_density = 0.0D0
    do i = -M, M
       if(v(i)>=lres .and. v(i)<=hres) then
          resonant_density = resonant_density + real(fv0(i))
       end if
    end do
    return
  end function resonant_density

  function resonant_velocity(fs0, M, s)
    implicit none
    integer, intent(in) :: M
    complex(kind=8), intent(in) :: fs0(-M:M)
    real(kind=8), intent(in) :: s(-M:M)
    real(kind=8) :: resonant_velocity
    complex(kind=8) :: fv0(-M:M)
    real(kind=8) :: v(-M:M)
    real(kind=8) :: dv 
    real(kind=8) :: lres, hres
    integer :: i
    lres = v0
    hres = v1
    dv = pi/s(M)
    forall(i=-M:M) v(i) = i*dv 
    call IFFT1D(fv0, M, fs0, s)
    resonant_velocity = 0.0D0
    do i = -M, M
       if(v(i)>=lres .and. v(i)<=hres) then
          resonant_velocity = resonant_velocity + real(v(i)*fv0(i))
       end if
    end do
    resonant_velocity = resonant_velocity/resonant_density(fs0, M, s)
    return
  end function resonant_velocity
end module is_tracking_mod


program test_is_tracking
  use is_tracking_mod
  use TAECH_IO
  use TAECH_LIB
  implicit none
  integer, parameter :: M=1000, N=200
  complex(kind=8) :: fsk(-M:M, 0:N)
  real(kind=8) :: s(-M:M)
  real(kind=8) :: xl, xh, vo, stdv
  real(kind=8) :: dv, vmax
  real(kind=8) :: a
  integer :: i, j, p, q
  call load_cfg('parameters.cfg')
  xl = -3.0D0
  xh = 1.0D0
  vo = -0.64D0
  stdv = 0.02D0
  vmax = 10.0D0
  dv = vmax/M
  forall(i=-M:M) s(i) = i*pi/(M*dv)
  forall(i=-M:M) fsk(i,0) = (xh-xl)*exp(imagj*s(i)*vo-s(i)**2*stdv**2/4.0D0)
  forall(i=-M:M, j=1:N) fsk(i,j)=imagj/j*(exp(-imagj*j*xh)-exp(-imagj*j*xl)) &
       * exp(imagj*s(i)*vo-s(i)**2*stdv**2/4.0D0)

  write(*,*) is_tracking(fsk(:,0), M, s)
  write(*,*) resonant_velocity(fsk(:,0), M, s)
  write(*,*) resonant_density(fsk(:,0), M, s)
  return

end program test_is_tracking
