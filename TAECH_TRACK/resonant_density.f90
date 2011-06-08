module resonant_density_mod
  use TAECH_IO
  use TAECH_LIB
  implicit none
  private
  
  public :: resonant_density
contains
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
end module resonant_density_mod


program test_resonant_density
  use resonant_density_mod
  use TAECH_IO
  use TAECH_LIB
  implicit none
  integer, parameter :: M=1000, N=200
  complex(kind=8) :: fsk(-M:M, 0:N)
  complex(kind=8) :: fv(-M:M)
  real(kind=8) :: s(-M:M)
  real(kind=8) :: xo, vo, stdv
  real(kind=8) :: dv, vmax
  real(kind=8) :: a
  integer :: i, j, p, q
  call load_cfg('parameters.cfg')
  xo = 0.2D0
  vo = -0.68D0
  stdv = 0.02D0
  vmax = 10.0D0
  dv = vmax/M
  forall(i=-M:M) s(i) = i*pi/(M*dv)
  forall(i=-M:M, j=0:N) fsk(i,j)=exp(-imagj*j*xo) &
       * exp(imagj*s(i)*vo-s(i)**2*stdv**2/4.0D0)
  call IFFT1D(fv, M, fsk(-M:M,0), s)
  call save_output('fv.dat', fv(-M:M), 2*M+1, .False.)
  write(*,*) resonant_density(fsk(-M:M,0), M, s)
  return
end program test_resonant_density
