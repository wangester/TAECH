module step_fun_mod
  use TAECH_IO
  implicit none
  private

  public :: step_fun
contains
  elemental function step_fun(s, val_old, val_half, val_now)
    use TAECH_IO
    implicit none
    real(kind=8), intent(in) :: s
    complex(kind=8), intent(in) :: val_old, val_half, val_now
    complex(kind=8) :: step_fun
    if (s >= 0.0D0 .and. s < dt) then
       step_fun = (s**2/dt**2-s/(2.0*dt))*val_old+(-2.0*s**2/dt**2+2.0*s/dt)*val_half &
            + (s**2/dt**2-3.0*s/(2.0*dt)+1.0/2.0)*val_now
    else 
       step_fun = (0.0D0, 0.0D0)
    end if
    return
  end function step_fun
end module step_fun_mod


program test_step_fun
  use TAECH_IO
  use step_fun_mod
  implicit none
  real(kind=8) :: s(0:10)
  complex(kind=8) :: a(0:10)
  complex(kind=8) :: val_old, val_half, val_now 
  integer :: i
  val_old = (0.02D0, 0.02D0)
  val_half = (0.05D0, 0.05D0)
  val_now = (0.1D0, 0.1D0)
  forall(i=0:10) s(i) = i*0.01D0
  call load_cfg('parameters.cfg')
  !do i = 0, 10
  !   a(i) = step_fun(s(i), val_old, val_half, val_now)
  !end do
  !call save_output('data.dat', a, 11, .False.)
  write(*,*) step_fun(s, val_old, val_half, val_now)
  return
end program test_step_fun
