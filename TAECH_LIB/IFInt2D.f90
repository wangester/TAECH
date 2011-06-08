module IFInt2D_mod
  use TAECH_IO
  implicit none
  private

  public :: IFInt2D
contains
  subroutine IFInt2D(fvx, vdim, fsk, M, N, s)
    implicit none
    include "fftw3.f"
    integer, intent(in) :: M, N, vdim
    real(kind=8), intent(out) :: fvx(-vdim:vdim, -N-1:N+1)
    complex(kind=8), intent(in) :: fsk(-M:M, 0:N)
    real(kind=8), intent(in) :: s(-M:M)
    real(kind=8) :: s_eq(-M:M)
    ! temporary variables
    real(kind=8) :: x(-N-1:N+1), v(-vdim:vdim)
    complex(kind=8) :: fvx_work(-vdim:vdim-1, -N-1:N)
    complex(kind=8) :: fsk_work(-M:M, 0:N)
    complex(kind=8) :: f(-vdim:vdim, -N-1:N+1)
    real(kind=8) :: f_real(-M:M, 0:N)
    real(kind=8) :: f_imag(-M:M, 0:N)
    complex(kind=8) :: f1(-vdim:vdim, -N-1:N+1)
    real(kind=8) :: f1_real(-M:M, 0:N)
    real(kind=8) :: f1_imag(-M:M, 0:N)
    complex(kind=8) :: f2(-vdim:vdim, -N-1:N+1)
    real(kind=8) :: f2_real(-M:M, 0:N)
    real(kind=8) :: f2_imag(-M:M, 0:N)
    complex(kind=8) :: f3(-vdim:vdim, -N-1:N+1)
    real(kind=8) :: f3_real(-M:M, 0:N)
    real(kind=8) :: f3_imag(-M:M, 0:N)
    integer :: i, p, q
    integer :: error
    real(kind=8) :: dx, dv, ds_eq
    !FFTW3 parameters
    integer(kind=8) :: plan
    !PSPLINE parameters (r8genxpkg)
    real(kind=8) :: xpkg(-M:M,4)
    integer :: iper = 0
    integer :: imsg = 0
    integer :: itol = 0
    real(kind=8) :: ztol
    integer :: ialg = -3
    integer :: ier = 1
    !PSPLINE parameters (r8cspline)
    real(kind=8) :: fspl(4, -M:M)
    integer :: ibcxmin = 3
    integer :: ibcxmax = 3
    real(kind=8) :: bcxmin
    real(kind=8) :: bcxmax
    real(kind=8) :: wk(2*M+1)
    integer :: ilinx = 1
    real(kind=8) :: xvec(-M:M)
    ds_eq = s(M)/M
    dx = pi/(N+1)
    dv = pi/ds_eq/vdim 
    forall(i=-N-1:N+1) x(i) = i*dx
    forall(i=-vdim:vdim) v(i) = i*dv
    ! convert to the uniform s grid
    forall(q=-M:M) s_eq(q) = q*ds_eq
    do p = 0, N
       call conv_fgrid(fsk_work(:,p), s_eq, M, fsk(:,p), s, M)  
    end do
    xvec = s_eq
    call r8genxpkg(2*M+1, s_eq, xpkg, iper, imsg, itol, ztol, ialg, ier)
    if (ier .ne. 0) then
       write(*,*) "error: chirp.dist_fun.r8genxpkg"
       stop
    end if
    ! calculate spline coefficients
    do p = 0, N
       fspl(1, -M:M) = real(fsk_work(-M:M, p))
       call r8cspline(xvec, 2*M+1, fspl, ibcxmin, bcxmin, ibcxmax, bcxmax, wk, 2*M+1, ilinx, ier)
       if (ier .ne. 0) then
          write(*,*) "error: chirp.dist_fun.r8cspline"
          stop
       end if
       f1_real(-M:M, p) = fspl(2, -M:M)
       f2_real(-M:M, p) = fspl(3, -M:M)
       f3_real(-M:M, p) = fspl(4, -M:M)
       fspl(1, -M:M) = aimag(fsk_work(-M:M, p))
       call r8cspline(xvec, 2*M+1, fspl, ibcxmin, bcxmin, ibcxmax, bcxmax, wk, 2*M+1, ilinx, ier)
       if (ier .ne. 0) then
          write(*,*) "error: chirp.dist_fun.r8cspline"
          stop
       end if
       f1_imag(-M:M, p) = fspl(2, -M:M)
       f2_imag(-M:M, p) = fspl(3, -M:M)
       f3_imag(-M:M, p) = fspl(4, -M:M)
    end do
    f = (0.0D0, 0.0D0)
    f1 = (0.0D0, 0.0D0)
    f2 = (0.0D0, 0.0D0)
    f3 = (0.0D0, 0.0D0)
    f(-M:M, 0:N) = fsk_work(-M:M, 0:N)
    f(-M:M, -1:-N:-1) = conjg(fsk_work(M:-M:-1, 1:N))
    f1(-M:M, 0:N) = cmplx(f1_real(-M:M, 0:N), f1_imag(-M:M, 0:N), kind=8)
    f1(-M:M, -1:-N:-1) = -conjg(f1(M:-M:-1, 1:N))
    f2(-M:M, 0:N) = cmplx(f2_real(-M:M, 0:N), f2_imag(-M:M, 0:N), kind=8)
    f2(-M:M, -1:-N:-1) = conjg(f2(M:-M:-1, 1:N))
    f3(-M:M, 0:N) = cmplx(f3_real(-M:M, 0:N), f3_imag(-M:M,0:N), kind=8)
    f3(-M:M, -1:-N:-1) = -conjg(f3(M:-M:-1, 1:N))
    ! f fftshift 
    forall(p=-vdim:-1, q=-N-1:-1) fvx_work(p,q) = f(-p-vdim, q+N+1)
    forall(p=-vdim:-1, q=0:N)  fvx_work(p,q) = f(-p-vdim, q-N-1)
    forall(p=0:vdim-1, q=-N-1:-1)  fvx_work(p,q) = f(-p+vdim, q+N+1)
    forall(p=0:vdim-1, q=0:N)   fvx_work(p,q) = f(-p+vdim, q-N-1)
    !2D IFFT on dist matrix
    call dfftw_plan_dft_2d(plan, 2*vdim, 2*N+2, fvx_work, fvx_work, FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, fvx_work, fvx_work)
    ! f fftshift again
    forall(p=-vdim:-1, q=-N-1:-1) f(p,q) = fvx_work(p+vdim, q+N+1)
    forall(p=-vdim:-1, q=0:N+1)  f(p,q) = fvx_work(p+vdim, q-N-1)
    forall(p=0:vdim, q=-N-1:-1)  f(p,q) = fvx_work(p-vdim, q+N+1)
    forall(p=0:vdim, q=0:N+1)   f(p,q) = fvx_work(p-vdim, q-N-1)
    ! f1 fftshift 
    forall(p=-vdim:-1, q=-N-1:-1) fvx_work(p,q) = f1(-p-vdim, q+N+1)
    forall(p=-vdim:-1, q=0:N)  fvx_work(p,q) = f1(-p-vdim, q-N-1)
    forall(p=0:vdim-1, q=-N-1:-1)  fvx_work(p,q) = f1(-p+vdim, q+N+1)
    forall(p=0:vdim-1, q=0:N)   fvx_work(p,q) = f1(-p+vdim, q-N-1)
    !2D IFFT on dist matrix
    call dfftw_execute_dft(plan, fvx_work, fvx_work)
    ! f1 fftshift again
    forall(p=-vdim:-1, q=-N-1:-1) f1(p,q) = fvx_work(p+vdim, q+N+1)
    forall(p=-vdim:-1, q=0:N+1)  f1(p,q) = fvx_work(p+vdim, q-N-1)
    forall(p=0:vdim, q=-N-1:-1)  f1(p,q) = fvx_work(p-vdim, q+N+1)
    forall(p=0:vdim, q=0:N+1)   f1(p,q) = fvx_work(p-vdim, q-N-1)
    ! f2 fftshift 
    forall(p=-vdim:-1, q=-N-1:-1) fvx_work(p,q) = f2(-p-vdim, q+N+1)
    forall(p=-vdim:-1, q=0:N)  fvx_work(p,q) = f2(-p-vdim, q-N-1)
    forall(p=0:vdim-1, q=-N-1:-1)  fvx_work(p,q) = f2(-p+vdim, q+N+1)
    forall(p=0:vdim-1, q=0:N)   fvx_work(p,q) = f2(-p+vdim, q-N-1)
    !2D IFFT on dist matrix
    call dfftw_execute_dft(plan, fvx_work, fvx_work)
    ! f2 fftshift again
    forall(p=-vdim:-1, q=-N-1:-1) f2(p,q) = fvx_work(p+vdim, q+N+1)
    forall(p=-vdim:-1, q=0:N+1)  f2(p,q) = fvx_work(p+vdim, q-N-1)
    forall(p=0:vdim, q=-N-1:-1)  f2(p,q) = fvx_work(p-vdim, q+N+1)
    forall(p=0:vdim, q=0:N+1)   f2(p,q) = fvx_work(p-vdim, q-N-1)
    ! f3 fftshift 
    forall(p=-vdim:-1, q=-N-1:-1) fvx_work(p,q) = f3(-p-vdim, q+N+1)
    forall(p=-vdim:-1, q=0:N)  fvx_work(p,q) = f3(-p-vdim, q-N-1)
    forall(p=0:vdim-1, q=-N-1:-1)  fvx_work(p,q) = f3(-p+vdim, q+N+1)
    forall(p=0:vdim-1, q=0:N)   fvx_work(p,q) = f3(-p+vdim, q-N-1)
    !2D IFFT on dist matrix
    call dfftw_execute_dft(plan, fvx_work, fvx_work)
    ! f3 fftshift again
    forall(p=-vdim:-1, q=-N-1:-1) f3(p,q) = fvx_work(p+vdim, q+N+1)
    forall(p=-vdim:-1, q=0:N+1)  f3(p,q) = fvx_work(p+vdim, q-N-1)
    forall(p=0:vdim, q=-N-1:-1)  f3(p,q) = fvx_work(p-vdim, q+N+1)
    forall(p=0:vdim, q=0:N+1)   f3(p,q) = fvx_work(p-vdim, q-N-1)
    call dfftw_destroy_plan(plan)
    forall(p=-vdim:-1, q=-N-1:N+1) fvx(p,q) = real(-imagj*(1.0-exp(-imagj*v(p)*ds_eq))/v(p)*f(p,q) &
         + (exp(-imagj*v(p)*ds_eq)-1.0+imagj*v(p)*ds_eq*exp(-imagj*v(p)*ds_eq))/v(p)**2*f1(p,q) &
         + (-2.0*imagj*exp(-imagj*v(p)*ds_eq)+2.0*imagj &
         + 2.0*v(p)*ds_eq*exp(-imagj*v(p)*ds_eq)+imagj*v(p)**2*ds_eq**2*exp(-imagj*v(p)*ds_eq)) &
         / v(p)**3*f2(p,q) &
         + (-6.0*exp(-imagj*v(p)*ds_eq)+6.0-6.0*imagj*v(p)*ds_eq*exp(-imagj*v(p)*ds_eq) &
         + 3.0*v(p)**2*ds_eq**2*exp(-imagj*v(p)*ds_eq)+imagj*v(p)**3*ds_eq**3*exp(-imagj*v(p)*ds_eq)) &
         / v(p)**4*f3(p,q))
    forall(p=1:vdim, q=-N-1:N+1) fvx(p,q) = real(-imagj*(1.0-exp(-imagj*v(p)*ds_eq))/v(p)*f(p,q) &
       + (exp(-imagj*v(p)*ds_eq)-1.0+imagj*v(p)*ds_eq*exp(-imagj*v(p)*ds_eq))/v(p)**2*f1(p,q) &
       + (-2.0*imagj*exp(-imagj*v(p)*ds_eq)+2.0*imagj &
       + 2.0*v(p)*ds_eq*exp(-imagj*v(p)*ds_eq)+imagj*v(p)**2*ds_eq**2*exp(-imagj*v(p)*ds_eq)) &
       / v(p)**3*f2(p,q) &
       + (-6.0*exp(-imagj*v(p)*ds_eq)+6.0-6.0*imagj*v(p)*ds_eq*exp(-imagj*v(p)*ds_eq) &
       + 3.0*v(p)**2*ds_eq**2*exp(-imagj*v(p)*ds_eq)+imagj*v(p)**3*ds_eq**3*exp(-imagj*v(p)*ds_eq)) &
       / v(p)**4*f3(p,q))
    forall(q=-N-1:N+1) fvx(0,q) = real(ds_eq*f(0,q)+ds_eq**2/2.0*f1(0,q) &
         + ds_eq**3/3.0*f2(0,q)+ds_eq**4/4.0*f3(0,q))
    return
  end subroutine IFInt2D

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
    ilinx = 2
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

end module IFInt2D_mod


program test_IFInt2D
  use IFInt2D_mod
  use TAECH_IO
  implicit none
  real(kind=8) :: fvx(-200:200, -11:11)
  complex(kind=8) :: fsk(-10:10, 0:10)
  real(kind=8) :: s(-10:10)
  integer :: i, j
  forall(i=-10:10) s(i) = i*0.1
  forall(i=-10:10,j=0:10) fsk(i,j) = imagj/(j+1)*exp(-10.0D0*s(i)**2)
  call IFInt2D(fvx, 200, fsk, 10, 10, s)
!  write(*,*) fvx
  call save_output('data.dat', fvx, 401, 23, .False.)

  return
end program test_IFInt2D
