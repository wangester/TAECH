module resonant_position_mod
  use TAECH_IO
  use TAECH_LIB
  implicit none
  private

  public :: resonant_position
  public :: resonant_density
contains
  function resonant_position(fsk, M, N, s)
    implicit none
    integer, intent(in) :: M, N
    complex(kind=8), intent(in) :: fsk(-M:M, 0:N)
    real(kind=8), intent(in) :: s(-M:M)
    real(kind=8) :: resonant_position
    complex(kind=8) :: fv(-M:M)
    complex(kind=8) :: fs(-M:M)
    real(kind=8) :: v(-M:M)
    real(kind=8) :: dv
    real(kind=8) :: lres, hres
    integer :: i
    lres = v0
    hres = v1
    dv = pi/s(M)
    forall(i=-M:M) v(i) = i*dv
    fs = (0.0D0, 0.0D0)
    do i = 1, N
       fs(-M:M) = fs(-M:M)-imagj*(-1)**i/i*(fsk(-M:M,i)-conjg(fsk(M:-M:-1,i)))
    end do
    call IFFT1D(fv, M, fs, s)
    resonant_position = 0.0D0
    do i = -M, M
       if(v(i)>=lres .and. v(i)<=hres) then
          resonant_position = resonant_position + real(fv(i))
       end if
    end do
    resonant_position = resonant_position/resonant_density(fsk(:,0), M, s)
    return
  end function resonant_position


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
          resonant_density = resonant_density+real(fv0(i))
       end if
    end do
    return
  end function resonant_density

end module resonant_position_mod


program test_resonant_position
  use resonant_position_mod
  use TAECH_IO
  use TAECH_LIB
  implicit none
  integer, parameter :: M=100, N=20
  complex(kind=8) :: fsk(-M:M, 0:N)
  real(kind=8) :: s(-M:M)
  real(kind=8) :: xl, xh, vo, stdv
  real(kind=8) :: dv, vmax
  real(kind=8) :: a
  integer :: i, j, p, q
  call load_cfg('parameters.cfg')
  xl = -3.0D0
  xh = 1.0D0
  vo = -0.68D0
  stdv = 0.02D0
  vmax = 10.0D0
  dv = vmax/M
  forall(i=-M:M) s(i) = i*pi/(M*dv)
  forall(i=-M:M) fsk(i,0) = (xh-xl)*exp(imagj*s(i)*vo-s(i)**2*stdv**2/4.0D0)
  forall(i=-M:M, j=1:N) fsk(i,j)=imagj/j*(exp(-imagj*j*xh)-exp(-imagj*j*xl)) &
       * exp(imagj*s(i)*vo-s(i)**2*stdv**2/4.0D0)

  write(*,*) resonant_position(fsk, M, N, s)
  return

end program test_resonant_position

