program test_TAECH_TRACK
  use TAECH_IO
  use TAECH_LIB
  use TAECH_TRACK
  implicit none
  include "fftw3.f"
  complex(kind=8) :: wave_ampl
  complex(kind=8) :: fsk(-100:100, 0:20)
  complex(kind=8) :: fvx(-100:100, -20:20)
  real(kind=8) :: x(-20:20), v(-100:100)
  complex(kind=8) :: fs0(-100:100), fs1(-100:100)
  real(kind=8) :: s(-100:100)
  real(kind=8) :: dx, dv, vmax
  integer :: i, j, p, q
  integer(kind=8) :: plan
  call load_cfg('parameters.cfg')
  wave_ampl = (1.0D0, -1.0D0)
  dx = pi/20
  forall(i=-20:20) x(i) = i*dx
  vmax = 1.0D0
  dv = vmax/100
  forall(i=-100:100) v(i) = i*dv
  forall(i=-100:100) s(i) = i*pi/(100*dv)
  forall(i=-100:100, j=-20:20) fvx(i,j)=exp(-(x(j)+1.0D0)**2/0.1D0**2) &
       * exp(-(v(i)+0.68D0)**2/0.02D0**2)
  call dfftw_plan_dft_2d(plan, 201, 41, fvx, fvx, FFTW_FORWARD, FFTW_ESTIMATE)
  call dfftw_execute_dft(plan, fvx, fvx)
  call dfftw_destroy_plan(plan)
  forall(p=-100:-1, q=0:20) fsk(p,q) = fvx(p+100, q-20)
  forall(p=0:100, q=0:20) fsk(p,q) = fvx(p-100,q-20)
  
  write(*,*) resonant_density(fsk(-100:100,0), 100, s)
  !resonant_position(fsk, M, N, s)
  !a = resonant_velocity(fsk(:,0), M, s)
  !frame_accel(fsk(:,0), fsk(:,1), M, wave_ampl, s)


  return
end program test_TAECH_TRACK
