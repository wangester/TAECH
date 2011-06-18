module TAECH_TRACK
  use TAECH_IO
  use TAECH_LIB
  implicit none
  private

  public :: is_resonant
  public :: chirp
  public :: resonant_density
  public :: resonant_position
  public :: resonant_velocity
  public :: resonant_filter

  integer, parameter :: res_size = 100
contains
  function is_resonant(subwave, twin, fs0, M, s)
    implicit none
    logical :: is_resonant
    integer, intent(in) :: twin
    complex(kind=8), intent(in) :: subwave(0:twin)
    integer, intent(in) :: M
    complex(kind=8), intent(in) :: fs0(-M:M)
    real(kind=8), intent(in) :: s(-M:M)
    real(kind=8) :: chirp_freq, chirp_ergy
    call chirp(chirp_freq, chirp_ergy, subwave, twin, resonance_bot, resonance_top)
    if ( chirp_ergy > 0.1D0 .and. &
         abs(resonant_velocity(fs0, M, s, resonance_bot, resonance_top) &
         - (resonance_bot+resonance_top)/2.0D0)<0.015D0 .and. &
         resonant_density(fs0, M, s, resonance_bot, resonance_top) &
         > 1.8D0*pi*(resonance_top-resonance_bot)**2) then
       is_resonant = .True.
    else
       is_resonant = .False.
    end if
    return
  end function is_resonant


  function resonant_density(fs0, M, s, lres, hres)
    implicit none
    integer, intent(in) :: M
    complex(kind=8), intent(in) :: fs0(-M:M)
    real(kind=8), intent(in) :: s(-M:M)
    real(kind=8), intent(in) :: lres, hres
    real(kind=8) :: resonant_density
    complex(kind=8) :: fv0(-M:M)
    complex(kind=8) :: res_fv0(0:res_size)  
    real(kind=8) :: v(-M:M)
    real(kind=8) :: dv 
    real(kind=8) :: res_v(0:res_size)
    real(kind=8) :: res_dv
    integer :: i
    dv = pi/s(M)
    forall(i=-M:M) v(i) = i*dv
    res_dv = (hres-lres)/res_size
    forall(i=0:res_size) res_v(i) = lres+i*res_dv
    call IFFT1D(fv0, M, fs0, s)
    call interpolate(res_fv0, res_v, res_size+1, fv0, v, M)
    resonant_density = 2.0D0*pi*real(sum(res_fv0)-0.5D0*(res_fv0(0)+res_fv0(res_size))) &
                     * res_dv+pi*(hres-lres)**2 
    return
  end function resonant_density
   

  function resonant_position(fsk, M, N, s, lres, hres)
    implicit none
    integer, intent(in) :: M, N
    complex(kind=8), intent(in) :: fsk(-M:M, 0:N)
    real(kind=8), intent(in) :: s(-M:M)
    real(kind=8), intent(in) :: lres, hres
    real(kind=8) :: resonant_position
    complex(kind=8) :: fv(-M:M)
    complex(kind=8) :: fs(-M:M)
    complex(kind=8) :: res_fv(0:res_size)
    real(kind=8) :: v(-M:M)
    real(kind=8) :: dv
    real(kind=8) :: res_v(0:res_size)
    real(kind=8) :: res_dv
    integer :: i
    dv = pi/s(M)
    forall(i=-M:M) v(i) = i*dv
    res_dv = (hres-lres)/res_size
    forall(i=0:res_size) res_v(i) = lres+i*res_dv
    fs = (0.0D0, 0.0D0)
    do i = 1, N
       fs(-M:M) = fs(-M:M)+imagj*(-1)**i/i*(conjg(fsk(M:-M:-1,i))-fsk(-M:M,i))
    end do
    call IFFT1D(fv, M, fs, s)
    call interpolate(res_fv, res_v, res_size+1, fv, v, M)
    resonant_position = 2.0D0*pi*real(sum(res_fv)-0.5D0*(res_fv(0)+res_fv(res_size)))*res_dv &
         / resonant_density(fsk(-M:M,0), M, s, lres, hres)
    return
  end function resonant_position


  function resonant_velocity(fs0, M, s, lres, hres)
    implicit none
    integer, intent(in) :: M
    complex(kind=8), intent(in) :: fs0(-M:M)
    real(kind=8), intent(in) :: s(-M:M)
    real(kind=8), intent(in) :: lres, hres
    real(kind=8) :: resonant_velocity
    complex(kind=8) :: fv0(-M:M)
    complex(kind=8) :: res_fv0(0:res_size)
    real(kind=8) :: v(-M:M)
    real(kind=8) :: dv 
    real(kind=8) :: res_v(0:res_size)
    real(kind=8) :: res_dv
    integer :: i
    dv = pi/s(M)
    forall(i=-M:M) v(i) = i*dv 
    res_dv = (hres-lres)/res_size
    forall(i=0:res_size) res_v(i) = lres+i*res_dv
    call IFFT1D(fv0, M, fs0, s)
    fv0 = v*fv0
    call interpolate(res_fv0, res_v, res_size+1, fv0, v, M)
    resonant_velocity = (2.0D0*pi*real(sum(res_fv0)-0.5D0*(res_fv0(0)+res_fv0(res_size)))*res_dv &
                      + 2.0D0/3.0D0*pi*(hres**3-lres**3)-pi*lres &
                      * (hres**2-lres**2))/resonant_density(fs0, M, s, lres, hres)
    return
  end function resonant_velocity


  function resonant_filter(fs0, fs1, M, s, lres, hres)
    implicit none
    integer, intent(in) :: M
    complex(kind=8), intent(in) :: fs1(-M:M)
    complex(kind=8), intent(in) :: fs0(-M:M)
    real(kind=8), intent(in) :: s(-M:M)
    real(kind=8), intent(in) :: lres, hres
    complex(kind=8) :: resonant_filter
    complex(kind=8) :: fv0(-M:M)
    complex(kind=8) :: fv1(-M:M)
    complex(kind=8) :: res_fv0(0:res_size)
    complex(kind=8) :: res_fv1(0:res_size)
    real(kind=8) :: v(-M:M)
    real(kind=8) :: dv
    real(kind=8) :: res_v(0:res_size)
    real(kind=8) :: res_dv
    real(kind=8) :: res_v0
    integer :: i
    dv = pi/s(M)
    forall(i=-M:M) v(i) = i*dv
    res_v0 = resonant_velocity(fs0, M, s, lres, hres) !(resonance_bot+resonance_top)/2.0D0   
    res_dv = (hres-lres)/res_size
    forall(i=0:res_size) res_v(i) = lres+i*res_dv
    call IFFT1D(fv0, M, fs0, s)
    call interpolate(res_fv0, res_v, res_size+1, fv0, v, M)
    call IFFT1D(fv1, M, fs1, s)
    call interpolate(res_fv1, res_v, res_size+1, fv1, v, M)
    resonant_filter = -((hres-res_v0)*res_fv1(res_size) &
                      - (lres-res_v0)*res_fv1(0) &
                      - (sum(res_fv1)-0.5D0*(res_fv1(0)+res_fv1(res_size)))*res_dv) &
                    / ( (hres-res_v0)*res_fv0(res_size) &
                      - (lres-res_v0)*res_fv0(0) &
                      - (sum(res_fv0)-0.5D0*(res_fv0(0)+res_fv0(res_size)))*res_dv)
    return
  end function resonant_filter


  subroutine chirp(chirp_freq, chirp_ergy, subwave, twin, lres, hres)
    implicit none
    include "fftw3.f"
    real(kind=8), intent(out) :: chirp_freq, chirp_ergy
    integer, intent(in) :: twin
    complex(kind=8), intent(in) :: subwave(0:twin)
    real(kind=8), intent(in) :: lres, hres
    !FFTW3 parameter
    integer(kind=8) :: plan
    complex(kind=8) :: subwave_wk(0:twin)
    complex(kind=8) :: wave_spectrum(0:twin)
    integer :: nlfreq, nhfreq
    real(kind=8) :: wave_energy(0:twin)
    real(kind=8) :: freq(0:twin)
    real(kind=8) :: dfreq
    real(kind=8) :: res_energy(0:res_size)
    real(kind=8) :: res_freq(0:res_size)
    real(kind=8) :: res_dfreq
    integer :: i, error
    nlfreq = -nint(twin/2.0)
    nhfreq = twin+nlfreq
    dfreq = 2.0*pi/(twin*dt)
    forall(i=0:twin) freq(i) = nlfreq*dfreq+i*dfreq 
    res_dfreq = (hres-lres)/res_size
    forall(i=0:res_size) res_freq(i) = lres+i*res_dfreq
    ! 1D FFT on subwave
    subwave_wk = subwave
    call dfftw_plan_dft_1d(plan, twin+1, subwave_wk, subwave_wk, FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, subwave_wk, subwave_wk)
    ! wave_spectrum fftshift
    forall(i=nlfreq:-1)  wave_spectrum(-nlfreq+i) = subwave_wk(i+twin)
    forall(i=0:nhfreq) wave_spectrum(-nlfreq+i) = subwave_wk(i)
    call dfftw_destroy_plan(plan)
    ! calculate the resonance energy(chirp_ergy) and chirp frequency(chirp_freq)
    wave_energy = 0.5D0*abs(wave_spectrum)**2
    call interpolate(res_energy, res_freq, res_size+1, wave_energy, freq, twin+1)
    chirp_ergy = (sum(res_energy)-0.5D0*(res_energy(0)+res_energy(res_size)))*res_dfreq &
         / ((sum(wave_energy)-0.5D0*(wave_energy(0)+wave_energy(twin)))*dfreq)
    chirp_freq = res_freq(maxloc(res_energy,1)-1)
    return
  end subroutine chirp

end module TAECH_TRACK


