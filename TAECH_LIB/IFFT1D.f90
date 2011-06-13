module IFFT1D_mod
  use TAECH_IO
  implicit none
  private

  public :: IFFT1D
contains
  subroutine IFFT1D(fv, M, fs, s)
    implicit none
    include "fftw3.f"
    integer, intent(in) :: M
    complex(kind=8), intent(out) :: fv(-M:M)
    complex(kind=8), intent(in) :: fs(-M:M)
    real(kind=8), intent(in) :: s(-M:M)
    real(kind=8) :: s_eq(-M:M)
    ! temporary variables
    complex(kind=8) :: work(-M:M-1)
    integer :: p
    ! FFTW3 parameters
    integer(kind=8) :: plan
    ! convert to the uniform s grid
    forall(p=-M:M) s_eq(p) = p*s(M)/M
    call conv_fgrid(fv, s_eq, M, fs, s, M)
    ! fv fftshift
    forall(p=-M:-1) work(p) = fv(-p-M)
    forall(p=0:M-1) work(p) = fv(-p+M)
    ! 1D IFFT on work matrix
    call dfftw_plan_dft_1d(plan, 2*M, work, work, FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, work, work)
    ! fv fftshift
    forall(p=-M:-1) fv(p) = work(p+M)
    forall(p=0:M) fv(p) = work(p-M)
    call dfftw_destroy_plan(plan)
    fv = fv*s(M)/M
    return
  end subroutine IFFT1D


  subroutine conv_fgrid(fs1, s1, M1, fs0, s0, M0)
    implicit none
    integer, intent(in) :: M1
    real(kind=8), intent(in) :: s1(-M1:M1)
    complex(kind=8), intent(inout) :: fs1(-M1:M1)
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

end module IFFT1D_mod


program test_IFFT1D
  use IFFT1D_mod
  use TAECH_IO
  implicit none
  complex(kind=8) :: fs(-1000:1000)
  real(kind=8) :: s(-1000:1000)
  complex(kind=8) :: fv(-1000:1000)
  integer :: i
  forall(i=-1000:1000) s(i) = i*0.01
  forall(i=-1000:1000) fs(i) = exp(-0.25D0*(s(i))**2)
  call load_cfg('parameters.cfg')
  call IFFT1D(fv, 1000, fs, s)
  call save_output('data.dat', fv, 2001, .False.)
  return
end program test_IFFT1D
