module resonant_peak_mod
  use TAECH_IO
  implicit none
  private

  public :: resonant_peak

contains
  subroutine resonant_peak(resonance, chirp_freq, subwave, twin, lres, hres)
    implicit none
    include "fftw3.f"
    real(kind=8), intent(out) :: resonance, chirp_freq
    integer, intent(in) :: twin
    complex(kind=8), intent(in) :: subwave(1:twin)
    real(kind=8), intent(in) :: lres, hres
    !FFTW3 parameter
    integer(kind=8) :: plan
    integer :: nlfreq, nhfreq, nfreq
    real(kind=8) :: dhfreq
    complex(kind=8), allocatable :: wave_spectrum(:)
    real(kind=8), allocatable :: hfreq(:)
    integer :: i, error
    real(kind=8) :: wave_energy, peakval
    nlfreq = -int(twin/2)
    nhfreq = twin+nlfreq-1
    nfreq = twin
    dhfreq = 2.0*pi/(nfreq*dt)
    allocate(hfreq(nlfreq:nhfreq), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for array hfreq in resonant_peak"
       stop
    end if
    forall(i=nlfreq:nhfreq) hfreq(i)=i*dhfreq 
    allocate(wave_spectrum(nlfreq:nhfreq), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for array wave_spectrum in resonant_peak"
       stop
    end if
    ! 1D FFT on subwave
    call dfftw_plan_dft_1d(plan, twin, subwave, subwave, FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, subwave, subwave)
    ! wave_spectrum fftshift
    forall(i=nlfreq:-1) wave_spectrum(i) = subwave(twin+i)
    forall(i=0:nhfreq)  wave_spectrum(i) = subwave(i+1)
    call dfftw_destroy_plan(plan)
    ! calculate the resonance energy(resonance) and chirp frequency(peakpos)
    resonance = 0.0D0
    wave_energy = 0.0D0
    peakval = 0.0D0
    do i = nlfreq, nhfreq
       wave_energy = abs(wave_spectrum(i))**2
       if (hfreq(i)>lres .and. hfreq(i)<hres) then
          resonance =  resonance + wave_energy
          if (wave_energy > peakval) then
             peakval = wave_energy
             chirp_freq = hfreq(i)
          end if
       end if
    end do
    wave_energy = sum(abs(wave_spectrum)**2)
    resonance = resonance/wave_energy
    deallocate(hfreq, stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for array hfreq in resonant_peak"
       stop
    end if
    deallocate(wave_spectrum, stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for array wave_spectrum in resonant_peak"
       stop
    end if
    return


  end subroutine resonant_peak
end module resonant_peak_mod



program test_resonant_peak
  use TAECH_IO
  use resonant_peak_mod
  implicit none
  complex(kind=8) :: subwave(60001)
  integer :: twin
  real(kind=8) :: lres, hres, peakpos, Eratio
  integer :: i
  call load_cfg('parameters.cfg')
  open(unit=10, file="wave_amplitude.o", action='read', status='old')
  do i = 1, 60001
     read(10, "(ES30.20E5, 1X, ES30.20E5)") subwave(i)
  end do
  close(unit=10)
  twin = 60001
  lres = -0.8D0
  hres = -0.6D0
  call resonant_peak(Eratio, peakpos, subwave, twin, lres, hres)
  write(*,*) peakpos
  write(*,*) Eratio
  return
end program test_resonant_peak
