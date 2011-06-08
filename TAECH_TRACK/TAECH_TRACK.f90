module TAECH_TRACK
  use TAECH_IO
  use TAECH_LIB
  implicit none
  private

  public :: is_resonant
  public :: is_tracking
  public :: resonant_density
  public :: resonant_position
  public :: resonant_velocity
  public :: frame_accel

contains
  function is_resonant(subwave, twin)
    implicit none
    logical :: is_resonant
    integer, intent(in) :: twin
    complex(kind=8), intent(in) :: subwave(1:twin)
    real(kind=8) :: lres, hres
    real(kind=8) :: resonance, chirp_freq
    lres = v0
    hres = v1
    call resonant_peak(resonance, chirp_freq, subwave, twin, lres, hres)
    if (chirp_freq>lres .and. chirp_freq<hres .and. resonance>0.1D0) then
       is_resonant = .True.
    else
       is_resonant = .False.
    end if
    return
  end function is_resonant


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
    forall(i=-M:M) v(i) = i*dv
    fs = (0.0D0, 0.0D0)
    do i = 1, N
       fs(-M:M) = fs(-M:M)-imagj*((-1.0D0)**i/i*fsk(-M:M,i) &
            + (-1.0D0)**(-i)/(-i)*conjg(fsk(M:-M:-1,i)))
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


  function frame_accel(fs0, fs1, M, wave_ampl, s)
    implicit none
    integer, intent(in) :: M
    complex(kind=8), intent(in) :: fs0(-M:M), fs1(-M:M)
    complex(kind=8), intent(in) :: wave_ampl
    real(kind=8), intent(in) :: s(-M:M)
    real(kind=8) :: frame_accel
    complex(kind=8) :: accel
    complex(kind=8) :: fv1(-M:M)
    real(kind=8) :: v(-M:M)
    real(kind=8) :: dv
    real(kind=8) :: lres, hres
    integer :: i
    lres = v0
    hres = v1
    dv = pi/s(M)
    forall(i=-M:M) v(i) = i*dv
    call IFFT1D(fv1, M, fs1, s)
    accel = (0.0D0, 0.0D0)
    do i = -M, M
       if(v(i)>=lres .and. v(i)<=hres) then
          accel = accel + fv1(i)
       end if
    end do
    frame_accel = -real(wave_ampl*conjg(accel))/resonant_density(fs0, M, s)
    return
  end function frame_accel

end module TAECH_TRACK



program test_TAECH_TRACK
  use TAECH_IO
  use TAECH_LIB
  use TAECH_TRACK
  implicit none
  include "fftw3.f"
  integer, parameter :: M=1000, N=200
  complex(kind=8) :: wave_ampl
  complex(kind=8) :: fsk(-M:M, 0:N)
  complex(kind=8) :: fs0(-M:M), fs1(-M:M)
  real(kind=8) :: s(-M:M)
  real(kind=8) :: xo, vo, stdx, stdv
  real(kind=8) :: dx, dv, vmax
  real(kind=8) :: a
  integer :: i, j, p, q
  call load_cfg('parameters.cfg')
  wave_ampl = (1.0D0, -1.0D0)
  xo = 10.0D0
  vo = -0.68D0
  stdx = 0.1D0
  stdv = 0.02D0
  dx = pi/N
  vmax = 10.0D0
  dv = vmax/M
  forall(i=-M:M) s(i) = i*pi/(M*dv)
  forall(i=-M:M, j=0:N) fsk(i,j)=exp(-imagj*j*xo-j**2*stdx**2/4.0D0) &
       * exp(imagj*s(i)*vo-s(i)**2*stdv**2/4.0D0)

  write(*,*) resonant_density(fsk(-M:M,0), M, s)
  write(*,*) resonant_position(fsk, M, N, s)
  write(*,*)  resonant_velocity(fsk(:,0), M, s)
  write(*,*) frame_accel(fsk(:,0), fsk(:,1), M, wave_ampl, s)
  
  return
end program test_TAECH_TRACK
