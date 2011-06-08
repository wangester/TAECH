module is_resonant_mod
  use TAECH_IO
  use TAECH_LIB
  implicit none
  private

  public :: is_resonant

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
  end function  is_resonant
end module is_resonant_mod


program test_is_resonant
  use TAECH_IO
  use TAECH_LIB
  use is_resonant_mod
  implicit none
  integer, parameter :: twin = 5000
  complex(kind=8) :: subwave(1:60001)
  integer :: i
  real(kind=8) :: resonance, chirp_freq
  call load_cfg('parameters.cfg')
  open(unit=10, file="wave_amplitude.o", action='read', status='old')
  do i = 1, 60001
     read(10, "(ES30.20E5, 1X, ES30.20E5)") subwave(i)
  end do
  close(unit=10)

  write(*,*) is_resonant(subwave(60001-twin+1:60001), twin)
  call  resonant_peak(resonance, chirp_freq, subwave(60001-twin+1:60001), twin, v0, v1)
  write(*,*) resonance, chirp_freq
  return
end program test_is_resonant
