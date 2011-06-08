module conv_fgrid_mod
  use TAECH_IO
  implicit none
  private
  
  public :: conv_fgrid
contains
  subroutine conv_fgrid(fs1, s1, M1, fs0, s0, M0)
    implicit none
    integer, intent(in) :: M1
    real(kind=8), intent(in) :: s1(-M1:M1)
    complex(kind=8), intent(out) :: fs1(-M1:M1)
    integer, intent(in) :: M0
    real(kind=8), intent(in) :: s0(-M0:M0)
    complex(kind=8), intent(in) :: fs0(-M0:M0)
    real(kind=8) :: fs1_real(-M1:M1)
    real(kind=8) :: fs1_imag(-M1:M1)
    !PSPLINE parameters (r8genxpkg)
    real(kind=8) :: xpkg(-M0:M0,4)
    integer :: iper = 0
    integer :: imsg = 0
    integer :: itol = 0 
    real(kind=8) :: ztol
    integer :: ialg = -3
    integer :: ier = 1
    !(r8cspline)
    real(kind=8) :: fspl(4,-M0:M0)
    integer :: ibcxmin = 3
    integer :: ibcxmax = 3
    real(kind=8) :: bcxmin
    real(kind=8) :: bcxmax
    real(kind=8) :: wk(2*M0+1)
    integer :: ilinx
    !(r8spvec)
    integer :: ict(3)
    integer :: ivec
    integer :: ivd
    integer :: iwarn
    real(kind=8) :: xvec(-M1:M1)
    real(kind=8) :: fval(-M1:M1,1)
    ivec = 2*M1+1
    ivd = ivec
    ict = 0
    ict(1) = 1
    xvec = s1
    !interpolate fs0 onto -M1:M1 grid of fs1
    call r8genxpkg(2*M0+1, s0, xpkg, iper, imsg, itol, ztol, ialg, ier)
    if (ier.ne.0) then
       write(*,*) "error: genxpkg in conv_fgrid"
       stop
    end if
    ! real component of fs1 : fs1_real
    fspl(1,-M0:M0) = real(fs0(-M0:M0))
    call r8cspline(s0, 2*M0+1, fspl, ibcxmin, bcxmin, ibcxmax, bcxmax, wk, 2*M0+1, ilinx, ier)
    if (ier.ne.0) then
       write(*,*) "error: cspline in conv_fgrid"
       stop
    end if
    call r8spvec(ict, ivec, xvec, ivd, fval, 2*M0+1, xpkg, fspl, iwarn, ier)
    if (ier.ne.0) then
       write(*,*) "error: spvec in conv_fgrid"
       stop
    end if
    fs1_real(-M1:M1) = fval(-M1:M1,1)
    ! imaginary component of fs1 : fs1_imag
    fspl(1,-M0:M0) = aimag(fs0(-M0:M0))
    call r8cspline(s0, 2*M0+1, fspl, ibcxmin, bcxmin, ibcxmax, bcxmax, wk, 2*M0+1, ilinx, ier)
    if (ier.ne.0) then
       write(*,*) "error: cspline in conv_fgrid"
       stop
    end if
    call r8spvec(ict, ivec, xvec, ivd, fval, 2*M0+1, xpkg, fspl, iwarn, ier)
    if (ier.ne.0) then
       write(*,*) "error: spvec in conv_fgrid"
       stop
    end if
    fs1_imag(-M1:M1) = fval(-M1:M1,1)
    fs1(-M1:M1) = cmplx(fs1_real(-M1:M1), fs1_imag(-M1:M1), kind=8)
    return
  end subroutine conv_fgrid
end module conv_fgrid_mod


program test_conv_fgrid
  use conv_fgrid_mod
  use TAECH_IO
  implicit none
  complex(kind=8) :: f0(-10:10), f1(-20:20)
  real(kind=8) :: s0(-10:10), s1(-20:20)
  integer :: i
  call load_cfg('parameters.cfg')
  forall(i=-10:10) s0(i) = i*ds
  forall(i=-10:10) f0(i) = (1.0D0+imagj)*exp(-100000.0D0*s0(i)**2)
  forall(i=-20:20) s1(i) = i*ds/2.0
  call conv_fgrid(f1, s1, 20, f0, s0, 10)
  call save_output('data1.dat', f1, 41, .False.)
  call save_output('data0.dat', f0, 21, .False.)
  return
end program test_conv_fgrid
