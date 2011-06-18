module TAECH_LIB
  use TAECH_IO
  implicit none
  private

  public :: besselj
  public :: zdotu
  public :: solve_tridiag
  public :: conv_fgrid
  public :: fspline
  public :: interpolate
  public :: IFFT1D
  public :: IFInt2D
  public :: visualize

  interface interpolate
     module procedure interpolate_Z
     module procedure interpolate_R
  end interface interpolate

contains
  elemental function besselj(X)
    !       =======================================================
    !       Purpose: Compute Bessel functions J0(x)-i*J1(x)
    !       Input :  x   --- Argument of Jn(x) 
    !       Output:  kernel = BJ0 - j*BJ1
    !       =======================================================
    implicit none
    real(kind=8), intent(in) :: X
    complex(kind=8) :: besselj
    !f2py intent(in) :: X
    !f2py intent(out) :: besselj
    real(kind=8) :: BJ0 , BJ1
    real(kind=8) :: T, T2
    real(kind=8) :: A0, P0, Q0, TA0, P1, Q1, TA1
    IF (X.EQ.0.0D0) THEN
       BJ0=1.0D0
       BJ1=0.0D0
    ELSE IF (X.LE.4.0D0) THEN
       T=X/4.0D0
       T2=T*T
       BJ0=((((((-.5014415D-3*T2+.76771853D-2)*T2 -.0709253492D0)*T2+.4443584263D0)*T2 &
            -1.7777560599D0)*T2+3.9999973021D0)*T2-3.9999998721D0)*T2+1.0D0 
       BJ1=T*(((((((-.1289769D-3*T2+.22069155D-2)*T2-.0236616773D0)*T2+.1777582922D0)*T2 &
            -.8888839649D0)*T2+2.6666660544D0)*T2-3.9999999710D0)*T2+1.9999999998D0)
    ELSE
       T=4.0D0/X
       T2=T*T
       A0=DSQRT(2.0D0/(PI*X))
       P0=((((-.9285D-5*T2+.43506D-4)*T2-.122226D-3)*T2+.434725D-3)*T2-.4394275D-2)*T2+.999999997D0 
       Q0=T*(((((.8099D-5*T2-.35614D-4)*T2+.85844D-4)*T2-.218024D-3)*T2+.1144106D-2)*T2-.031249995D0)
       TA0=X-.25D0*PI
       BJ0=A0*(P0*DCOS(TA0)-Q0*DSIN(TA0))
       P1=((((.10632D-4*T2-.50363D-4)*T2+.145575D-3)*T2-.559487D-3)*T2+.7323931D-2)*T2+1.000000004D0
       Q1=T*(((((-.9173D-5*T2+.40658D-4)*T2-.99941D-4)*T2+.266891D-3)*T2-.1601836D-2)*T2+.093749994D0)
       TA1=X-.75D0*PI
       BJ1=A0*(P1*DCOS(TA1)-Q1*DSIN(TA1))
    ENDIF
    besselj = cmplx(BJ0, -BJ1, kind=8)
    return
  END function besselj


  subroutine solve_tridiag(sub_diag, diag, sup_diag, rhs, sol, num_eqns)
    implicit none
    !      sub_diag - sub-diagonal (means it is the diagonal below the main diagonal)
    !      diag - the main diagonal
    !      sup_diag - sup-diagonal (means it is the diagonal above the main diagonal)
    !      rhs - right part
    !      sol - solution
    !      num_eqns - number of equations
    integer, intent(in) :: num_eqns
    complex(kind=8), intent(in) :: sub_diag(1:num_eqns)
    complex(kind=8), intent(in) :: diag(1:num_eqns)
    complex(kind=8), intent(in) :: sup_diag(1:num_eqns)
    complex(kind=8), intent(in) :: rhs(1:num_eqns)
    complex(kind=8), intent(out) :: sol(1:num_eqns)
    complex(kind=8) :: bp(1:num_eqns), vp(1:num_eqns)
    complex(kind=8) :: m
    integer :: i
    ! Make copies of the b and v variables so that they are unaltered by this sub
    bp(1) = diag(1)
    vp(1) = rhs(1)
    !The first pass (setting coefficients):
    firstpass: do i = 2, num_eqns
       m = sub_diag(i)/bp(i-1)
       bp(i) = diag(i) - m*sup_diag(i-1)
       vp(i) = rhs(i) - m*vp(i-1)
    end do firstpass
    sol(num_eqns) = vp(num_eqns)/bp(num_eqns)
    !The second pass (back-substition)
    backsub: do i = num_eqns-1, 1, -1
       sol(i) = (vp(i) - sup_diag(i)*sol(i+1))/bp(i)
    end do backsub
    return
  end subroutine solve_tridiag


  pure function zdotu(N,ZX,INCX,ZY,INCY)
    implicit none
    complex(kind=8) :: zdotu
    integer, intent(in) :: INCX,INCY,N
    complex(kind=8), intent(in) :: ZX(*), ZY(*)
    complex(kind=8) :: ZTEMP
    integer :: I,IX,IY
    ZTEMP = (0.0D0,0.0D0)
    zdotu = (0.0D0,0.0D0)
    IF (N.LE.0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
       DO I = 1,N
          ZTEMP = ZTEMP + ZX(I)*ZY(I)
       END DO
    ELSE
       IX = 1
       IY = 1
       IF (INCX.LT.0) IX = (-N+1)*INCX + 1
       IF (INCY.LT.0) IY = (-N+1)*INCY + 1
       DO I = 1,N
          ZTEMP = ZTEMP + ZX(IX)*ZY(IY)
          IX = IX + INCX
          IY = IY + INCY
       END DO
    END IF
    zdotu = ZTEMP
    RETURN
  END function zdotu


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
    call conv_fgrid(fv(-M:M), s_eq(-M:M), M, fs(-M:M), s(-M:M), M)
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
    fv = s(M)/M*fv
    return
  end subroutine IFFT1D
   

  subroutine IFInt2D(fvx, vlres, vhres, vdim, fsk, M, N, s)
    implicit none
    include "fftw3.f"
    integer, intent(in) :: M, N, vdim
    real(kind=8), intent(in) :: vlres, vhres
    real(kind=8), intent(out) :: fvx(0:vdim, -N-1:N+1)
    complex(kind=8), intent(in) :: fsk(-M:M, 0:N)
    real(kind=8), intent(in) :: s(-M:M)
    real(kind=8) :: s_eq(-M:M)
    ! temporary variables
    integer :: padding
    real(kind=8) :: x(-N-1:N+1)
    real(kind=8) :: v(0:vdim)
    complex(kind=8) :: fsk_wk(-M:M, 0:N)
    complex(kind=8), allocatable :: fvx_wk(:, :)
    complex(kind=8), allocatable  :: f(:, :)
    real(kind=8) :: f_real(-M:M, 0:N)
    real(kind=8) :: f_imag(-M:M, 0:N)
    complex(kind=8), allocatable :: f1(:, :)
    real(kind=8) :: f1_real(-M:M, 0:N)
    real(kind=8) :: f1_imag(-M:M, 0:N)
    complex(kind=8), allocatable :: f2(:, :)
    real(kind=8) :: f2_real(-M:M, 0:N)
    real(kind=8) :: f2_imag(-M:M, 0:N)
    complex(kind=8), allocatable :: f3(:, :)
    real(kind=8) :: f3_real(-M:M, 0:N)
    real(kind=8) :: f3_imag(-M:M, 0:N)
    integer :: i, p, q
    integer :: error = 0
    real(kind=8) :: dx, dv, ds_eq
    integer :: nvlres, nvhres
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
    integer :: ibcxmin = 0
    integer :: ibcxmax = 0
    real(kind=8) :: bcxmin 
    real(kind=8) :: bcxmax 
    real(kind=8) :: wk(2*M+1)
    integer :: ilinx = 1
    ds_eq = s(M)/M
    forall(i=-M:M) s_eq(i) = i*ds_eq
    dx = pi/(N+1)
    dv = (vhres-vlres)/vdim
    padding = nint(pi/(ds_eq*dv))
    if (padding < M) then
       write(*,*) "vdim > ", (vhres-vlres)*s(M)/pi
       stop
    end if
    forall(i=-N-1:N+1) x(i) = i*dx
    nvlres = nint(vlres/dv) 
    nvhres = vdim+nvlres
    forall(i=0:vdim) v(i) = vlres+i*dv
    allocate(fvx_wk(-padding:padding-1, -N-1:N), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for fvx_wk in TAECH_LIB.IFInt2D"
       stop
    end if
    allocate(f(-padding:padding, -N-1:N+1), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for f in TAECH_LIB.IFInt2D"
       stop
    end if
    allocate(f1(-padding:padding, -N-1:N+1), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for f1 in TAECH_LIB.IFInt2D"
       stop
    end if
    allocate(f2(-padding:padding, -N-1:N+1), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for f2 in TAECH_LIB.IFInt2D"
       stop
    end if
    allocate(f3(-padding:padding, -N-1:N+1), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for f3 in TAECH_LIB.IFInt2D"
       stop
    end if
    call r8genxpkg(2*M+1, s_eq, xpkg, iper, imsg, itol, ztol, ialg, ier)
    if (ier .ne. 0) then
       write(*,*) "error: TAECH_LIB.IFInt2D.r8genxpkg"
       stop
    end if
    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP PRIVATE(fspl, wk) &
    !$OMP FIRSTPRIVATE(ier)
    do p = 0, N
       ! nonuniform s grid converts to uniform s grid
       call conv_fgrid(fsk_wk(-M:M,p), s_eq(-M:M), M, fsk(-M:M,p), s(-M:M), M)
       ! calculate spline coefficients
       fspl(1, -M:M) = real(fsk_wk(-M:M, p))
       call r8cspline(s_eq, 2*M+1, fspl, ibcxmin, bcxmin, ibcxmax, bcxmax, wk, 2*M+1, ilinx, ier)
       if (ier .ne. 0) then
          write(*,*) "error: TAECH_LIB.IFInt2D.r8cspline"
          stop
       end if
       f1_real(-M:M, p) = fspl(2, -M:M)
       f2_real(-M:M, p) = fspl(3, -M:M)
       f3_real(-M:M, p) = fspl(4, -M:M)
       fspl(1, -M:M) = aimag(fsk_wk(-M:M, p))
       call r8cspline(s_eq, 2*M+1, fspl, ibcxmin, bcxmin, ibcxmax, bcxmax, wk, 2*M+1, ilinx, ier)
       if (ier .ne. 0) then
          write(*,*) "error: TAECH_LIB.IFInt2D.r8cspline"
          stop
       end if
       f1_imag(-M:M, p) = fspl(2, -M:M)
       f2_imag(-M:M, p) = fspl(3, -M:M)
       f3_imag(-M:M, p) = fspl(4, -M:M)
    end do
    !$OMP END PARALLEL DO
    f = (0.0D0, 0.0D0)
    f1 = (0.0D0, 0.0D0)
    f2 = (0.0D0, 0.0D0)
    f3 = (0.0D0, 0.0D0)
    f(-M:M, 0:N) = fsk_wk(-M:M, 0:N)
    f(-M:M, -1:-N:-1) = conjg(fsk_wk(M:-M:-1, 1:N))
    f1(-M:M, 0:N) = cmplx(f1_real(-M:M, 0:N), f1_imag(-M:M, 0:N), kind=8)
    f1(-M:M, -1:-N:-1) = -conjg(f1(M:-M:-1, 1:N))
    f2(-M:M, 0:N) = cmplx(f2_real(-M:M, 0:N), f2_imag(-M:M, 0:N), kind=8)
    f2(-M:M, -1:-N:-1) = conjg(f2(M:-M:-1, 1:N))
    f3(-M:M, 0:N) = cmplx(f3_real(-M:M, 0:N), f3_imag(-M:M,0:N), kind=8)
    f3(-M:M, -1:-N:-1) = -conjg(f3(M:-M:-1, 1:N))
    ! f fftshift 
    forall(p=-padding:-1, q=-N-1:-1) fvx_wk(p,q) = f(-p-padding, q+N+1)
    forall(p=-padding:-1, q=0:N)  fvx_wk(p,q) = f(-p-padding, q-N-1)
    forall(p=0:padding-1, q=-N-1:-1)  fvx_wk(p,q) = f(-p+padding, q+N+1)
    forall(p=0:padding-1, q=0:N)   fvx_wk(p,q) = f(-p+padding, q-N-1)
    !2D IFFT on dist matrix
    call dfftw_plan_dft_2d(plan, 2*padding, 2*N+2, fvx_wk, fvx_wk, FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, fvx_wk, fvx_wk)
    ! f fftshift again
    forall(p=-padding:-1, q=-N-1:-1) f(p,q) = fvx_wk(p+padding, q+N+1)
    forall(p=-padding:-1, q=0:N+1)  f(p,q) = fvx_wk(p+padding, q-N-1)
    forall(p=0:padding, q=-N-1:-1)  f(p,q) = fvx_wk(p-padding, q+N+1)
    forall(p=0:padding, q=0:N+1)   f(p,q) = fvx_wk(p-padding, q-N-1)
    ! f1 fftshift 
    forall(p=-padding:-1, q=-N-1:-1) fvx_wk(p,q) = f1(-p-padding, q+N+1)
    forall(p=-padding:-1, q=0:N)  fvx_wk(p,q) = f1(-p-padding, q-N-1)
    forall(p=0:padding-1, q=-N-1:-1)  fvx_wk(p,q) = f1(-p+padding, q+N+1)
    forall(p=0:padding-1, q=0:N)   fvx_wk(p,q) = f1(-p+padding, q-N-1)
    !2D IFFT on dist matrix
    call dfftw_execute_dft(plan, fvx_wk, fvx_wk)
    ! f1 fftshift again
    forall(p=-padding:-1, q=-N-1:-1) f1(p,q) = fvx_wk(p+padding, q+N+1)
    forall(p=-padding:-1, q=0:N+1)  f1(p,q) = fvx_wk(p+padding, q-N-1)
    forall(p=0:padding, q=-N-1:-1)  f1(p,q) = fvx_wk(p-padding, q+N+1)
    forall(p=0:padding, q=0:N+1)   f1(p,q) = fvx_wk(p-padding, q-N-1)
    ! f2 fftshift 
    forall(p=-padding:-1, q=-N-1:-1) fvx_wk(p,q) = f2(-p-padding, q+N+1)
    forall(p=-padding:-1, q=0:N)  fvx_wk(p,q) = f2(-p-padding, q-N-1)
    forall(p=0:padding-1, q=-N-1:-1)  fvx_wk(p,q) = f2(-p+padding, q+N+1)
    forall(p=0:padding-1, q=0:N)   fvx_wk(p,q) = f2(-p+padding, q-N-1)
    !2D IFFT on dist matrix
    call dfftw_execute_dft(plan, fvx_wk, fvx_wk)
    ! f2 fftshift again
    forall(p=-padding:-1, q=-N-1:-1) f2(p,q) = fvx_wk(p+padding, q+N+1)
    forall(p=-padding:-1, q=0:N+1)  f2(p,q) = fvx_wk(p+padding, q-N-1)
    forall(p=0:padding, q=-N-1:-1)  f2(p,q) = fvx_wk(p-padding, q+N+1)
    forall(p=0:padding, q=0:N+1)   f2(p,q) = fvx_wk(p-padding, q-N-1)
    ! f3 fftshift 
    forall(p=-padding:-1, q=-N-1:-1) fvx_wk(p,q) = f3(-p-padding, q+N+1)
    forall(p=-padding:-1, q=0:N)  fvx_wk(p,q) = f3(-p-padding, q-N-1)
    forall(p=0:padding-1, q=-N-1:-1)  fvx_wk(p,q) = f3(-p+padding, q+N+1)
    forall(p=0:padding-1, q=0:N)   fvx_wk(p,q) = f3(-p+padding, q-N-1)
    !2D IFFT on dist matrix
    call dfftw_execute_dft(plan, fvx_wk, fvx_wk)
    ! f3 fftshift again
    forall(p=-padding:-1, q=-N-1:-1) f3(p,q) = fvx_wk(p+padding, q+N+1)
    forall(p=-padding:-1, q=0:N+1)  f3(p,q) = fvx_wk(p+padding, q-N-1)
    forall(p=0:padding, q=-N-1:-1)  f3(p,q) = fvx_wk(p-padding, q+N+1)
    forall(p=0:padding, q=0:N+1)   f3(p,q) = fvx_wk(p-padding, q-N-1)
    call dfftw_destroy_plan(plan)
    do p = nvlres, nvhres
       if (p == 0) then
          forall(q=-N-1:N+1) fvx(p-nvlres,q) = real(ds_eq*f(0,q)+ds_eq**2/2.0*f1(0,q) &
               + ds_eq**3/3.0*f2(0,q)+ds_eq**4/4.0*f3(0,q))
       else
          forall(q=-N-1:N+1) fvx(p-nvlres,q) = real(-imagj*(1.0D0-exp(-imagj*v(p-nvlres)*ds_eq)) &
               / v(p-nvlres)*f(p,q)+(exp(-imagj*v(p-nvlres)*ds_eq)-1.0D0 &
               + imagj*v(p-nvlres)*ds_eq*exp(-imagj*v(p-nvlres)*ds_eq))/v(p-nvlres)**2*f1(p,q) &
               + (-2.0D0*imagj*exp(-imagj*v(p-nvlres)*ds_eq)+2.0D0*imagj &
               + 2.0D0*v(p-nvlres)*ds_eq*exp(-imagj*v(p-nvlres)*ds_eq) &
               + imagj*v(p-nvlres)**2*ds_eq**2*exp(-imagj*v(p-nvlres)*ds_eq)) &
               / v(p-nvlres)**3*f2(p,q) &
               + (-6.0D0*exp(-imagj*v(p-nvlres)*ds_eq)+6.0D0 &
               - 6.0D0*imagj*v(p-nvlres)*ds_eq*exp(-imagj*v(p-nvlres)*ds_eq) &
               + 3.0D0*v(p-nvlres)**2*ds_eq**2*exp(-imagj*v(p-nvlres)*ds_eq) &
               + imagj*v(p-nvlres)**3*ds_eq**3*exp(-imagj*v(p-nvlres)*ds_eq)) &
               / v(p-nvlres)**4*f3(p,q)) 
       end if
    end do
    deallocate(fvx_wk, stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for fvx_wk in TAECH_LIB.IFInt2D"
       stop
    end if
    deallocate(f, stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for f in TAECH_LIB.IFInt2D"
       stop
    end if
    deallocate(f1, stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for f1 in TAECH_LIB.IFInt2D"
       stop
    end if
    deallocate(f2, stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for f2 in TAECH_LIB.IFInt2D"
       stop
    end if
    deallocate(f3, stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for f3 in TAECH_LIB.IFInt2D"
       stop
    end if
    return
  end subroutine IFInt2D


  subroutine visualize(fvx, vlres, vhres, vdim, fbnd, fsk, M, N, nt, s) 
    implicit none
    integer, intent(in) :: M, N
    integer, intent(in) :: vdim, nt
    real(kind=8), intent(in) :: vlres, vhres
    real(kind=8), intent(in) :: s(-M:M)
    complex(kind=8), intent(in) :: fsk(-M:M, 0:N)
    complex(kind=8), intent(in) :: fbnd(1:N, 0:nt)
    real(kind=8), intent(out) :: fvx(0:vdim, -N-1:N+1)
    real(kind=8) :: fvx_wk(0:vdim,-N-1:N+1)
    complex(kind=8) :: fbnd_wk(1:N, 0:nt)
    real(kind=8) :: v(0:vdim), x(-N-1:N+1)
    real(kind=8) :: dv, dx
    integer :: i, j, p, q
    integer :: nvlres, nvhres
    dx = pi/(N+1)
    forall(i=-N-1:N+1) x(i) = i*dx
    dv = (vhres-vlres)/vdim
    nvlres = nint(vlres/dv) 
    nvhres = vdim+nvlres
    forall(i=0:vdim) v(i) = nvlres+i*dv
    call IFInt2D(fvx_wk, vlres, vhres, vdim, fsk, M, N, s)
    fvx = fvx_wk
    !do j = -N-1, N+1
    !   do i = nvlres, nvhres
    !      forall(p=1:N, q=0:nt) fbnd_wk(p,q) = exp(imagj*p*(x(j)-v(i-nvlres)*nt*dt)) &
    !           * p*exp(imagj*p*v(i-nvlres)*q*dt)*fbnd(p,q)
    !      fvx_wk(i-nvlres, j) = 2.0D0*real(exp(-imagj*s(M)*v(i-nvlres)) &
    !           * sum(dt*(sum(fbnd, dim=2)-0.5D0*(fbnd(:,0)+fbnd(:,nt))), dim=1))
    !   end do
    !end do
    !fvx = fvx + fvx_wk
    return
  end subroutine visualize

  
  subroutine interpolate_Z(finp, sgrid, nsgrid, fs, s, M)
    implicit none
    integer, intent(in) :: nsgrid, M
    complex(kind=8), intent(in) :: fs(-M:M)
    real(kind=8), intent(in) :: s(-M:M)
    real(kind=8), intent(in) :: sgrid(1:nsgrid)
    complex(kind=8), intent(out) :: finp(1:nsgrid)
    real(kind=8) :: sreal(1:nsgrid), simag(1:nsgrid)
    !internal variables (r8genxpkg)
    real(kind=8) :: xpkg(-M:M,4)
    integer :: iper = 0
    integer :: imsg = 0
    integer :: itol = 0
    real(kind=8) :: ztol
    integer :: ialg = -3
    integer :: ier = 1
    !(r8cspline)
    real(kind=8) :: fspl(4,-M:M)
    integer :: ibcxmin = 3
    integer :: ibcxmax = 3
    real(kind=8) :: bcxmin = 0.0D0
    real(kind=8) :: bcxmax = 0.0D0
    real(kind=8) :: wk(2*M+1)
    integer :: ilinx = 2
    !(r8spvec)
    integer :: ict(3)
    integer :: ivec
    integer :: ivd
    integer :: iwarn
    real(kind=8) :: xvec(1:nsgrid)
    real(kind=8) :: fval(1:nsgrid,1)
    ict = 0
    ict(1) = 1
    ivec = nsgrid
    xvec = sgrid
    ivd = ivec
    call r8genxpkg(2*M+1, s, xpkg, iper, imsg, itol, ztol, ialg, ier)
    if (ier.ne.0) then
       write(*,*) "error: genxpkg in TAECH_LIB.interpolate_Z."
       stop
    end if
    ! real component of fs: sreal
    fspl(1,-M:M) = real(fs(-M:M))
    call r8cspline(s, 2*M+1, fspl, ibcxmin, bcxmin, ibcxmax, bcxmax, wk, 2*M+1, ilinx, ier)
    if (ier.ne.0) then
       write(*,*) "error: cspline in TAECH_LIB.interpolate_Z."
       stop
    end if
    call r8spvec(ict, ivec, xvec, ivd, fval, 2*M+1, xpkg, fspl, iwarn, ier)
    if (ier.ne.0) then
       write(*,*) "error: spvec in TAECH_LIB.interpolate_Z."
       stop
    end if
    sreal(1:nsgrid) = fval(1:nsgrid, 1)
    ! imaginary component of f: simag
    fspl(1,-M:M) = aimag(fs(-M:M))
    call r8cspline(s, 2*M+1, fspl, ibcxmin, bcxmin, ibcxmax, bcxmax, wk, 2*M+1, ilinx, ier)
    if (ier.ne.0) then
       write(*,*) "error: cspline in TAECH_LIB.interpolate_Z."
       stop
    end if
    call r8spvec(ict, ivec, xvec, ivd, fval, 2*M+1, xpkg, fspl, iwarn, ier)
    if (ier.ne.0) then
       write(*,*) "error: spvec in TAECH_LIB.interpolate_Z."
       stop
    end if
    simag(1:nsgrid) = fval(1:nsgrid, 1)
    finp = cmplx(sreal, simag, kind=8)
    return
  end subroutine interpolate_Z


  subroutine interpolate_R(yinp, x1, nx1, y, x0, nx0)
    implicit none
    integer, intent(in) :: nx0, nx1
    real(kind=8), intent(in) :: y(1:nx0)
    real(kind=8), intent(in) :: x0(1:nx0)
    real(kind=8), intent(in) :: x1(1:nx1)
    real(kind=8), intent(out) :: yinp(1:nx1)
    !internal variables (r8genxpkg)
    real(kind=8) :: xpkg(1:nx0, 4)
    integer :: iper = 0
    integer :: imsg = 0
    integer :: itol = 0
    real(kind=8) :: ztol
    integer :: ialg = -3
    integer :: ier = 1
    !(r8cspline)
    real(kind=8) :: fspl(4, 1:nx0)
    integer :: ibcxmin = 0
    integer :: ibcxmax = 0
    real(kind=8) :: bcxmin 
    real(kind=8) :: bcxmax 
    real(kind=8) :: wk(nx0)
    integer :: ilinx = 1
    !(r8spvec)
    integer :: ict(3)
    integer :: ivec
    integer :: ivd
    integer :: iwarn
    real(kind=8) :: xvec(1:nx1)
    real(kind=8) :: fval(1:nx1,1)
    ict = 0
    ict(1) = 1
    ivec = nx1
    xvec = x1
    ivd = ivec
    call r8genxpkg(nx0, x0, xpkg, iper, imsg, itol, ztol, ialg, ier)
    if (ier.ne.0) then
       write(*,*) "error: genxpkg in TAECH_LIB.interpolate_R."
       stop
    end if
    fspl(1,1:nx0) = y(1:nx0)
    call r8cspline(x0, nx0, fspl, ibcxmin, bcxmin, ibcxmax, bcxmax, wk, nx0, ilinx, ier)
    if (ier.ne.0) then
       write(*,*) "error: cspline in TAECH_LIB.interpolate_R."
       stop
    end if
    call r8spvec(ict, ivec, xvec, ivd, fval, nx0, xpkg, fspl, iwarn, ier)
    if (ier.ne.0) then
       write(*,*) "error: spvec in TAECH_LIB.interpolate_R."
       stop
    end if
    yinp(1:nx1) = fval(1:nx1, 1)
    return
  end subroutine interpolate_R

    
  subroutine fspline(fs, s, M, incx)
    implicit none
    integer, intent(in) :: M
    complex(kind=8), intent(inout) :: fs(-M:M)
    real(kind=8), intent(in) :: s(-M:M)
    real(kind=8), intent(in) :: incx
    real(kind=8) :: freal(-M:M)
    real(kind=8) :: fimag(-M:M)
    !internal variables (r8genxpkg)
    real(kind=8) :: xpkg(-M:M,4)
    integer :: iper = 0
    integer :: imsg = 0
    integer :: itol = 0
    real(kind=8) :: ztol
    integer :: ialg = -3
    integer :: ier = 1
    !(r8cspline)
    real(kind=8) :: fspl(4,-M:M)
    integer :: ibcxmin = 3
    integer :: ibcxmax = 3 
    real(kind=8) :: bcxmin = 0.0D0
    real(kind=8) :: bcxmax = 0.0D0
    real(kind=8) :: wk(2*M+1)
    integer :: ilinx = 2
    !(r8spvec)
    integer :: ict(3)
    integer :: ivec
    integer :: ivd
    integer :: iwarn
    real(kind=8), allocatable :: xvec(:)
    real(kind=8), allocatable :: fval(:,:)
    integer :: i, lftbc, rgbc
    integer :: error
    !if (incx == 0.0D0) then
    !   return
    !end if
    call r8genxpkg(2*M+1, s, xpkg, iper, imsg, itol, ztol, ialg, ier)
    if (ier.ne.0) then
       write(*,*) "error: genxpkg in TAECH_LIB.fspline."
       stop
    end if
    ! real component of fs: freal
    fspl(1,-M:M) = real(fs(-M:M))
    call r8cspline(s, 2*M+1, fspl, ibcxmin, bcxmin, ibcxmax, bcxmax, wk, 2*M+1, ilinx, ier)
    if (ier.ne.0) then
       write(*,*) "error: cspline in TAECH_LIB.fspline."
       stop
    end if
    ! if the shifted location moves beyond the range of s(-M:M)
    lftbc = -M
    do i = -M, M
       if (s(i)+incx < s(-M)) then
          lftbc = lftbc + 1
       end if
    end do
    rgbc = M
    do i = -M, M
       if (s(i)+incx > s(M)) then
          rgbc = rgbc - 1
       end if
    end do
    ivec = rgbc-lftbc+1
    ivd = ivec
    allocate(xvec(lftbc:rgbc),stat=error)  
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for array xvec in TAECH_LIB.fspline"
       stop
    end if
    allocate(fval(lftbc:rgbc,1),stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for array fval in TAECH_LIB.fspline"
       stop
    end if
    ict = 0
    ict(1) = 1
    xvec(lftbc:rgbc) = s(lftbc:rgbc)+incx
    call r8spvec(ict, ivec, xvec, ivd, fval, 2*M+1, xpkg, fspl, iwarn, ier)
    if (ier.ne.0) then
       write(*,*) "error: spvec in TAECH_LIB.fspline."
       stop
    end if
    freal = 0.0D0
    freal(lftbc:rgbc) = fval(lftbc:rgbc,1)
    ! imaginary component of f: fimag
    fspl(1,-M:M) = aimag(fs(-M:M))
    call r8cspline(s, 2*M+1, fspl, ibcxmin, bcxmin, ibcxmax, bcxmax, wk, 2*M+1, ilinx, ier)
    if (ier.ne.0) then
       write(*,*) "error: cspline in TAECH_LIB.fspline."
       stop
    end if
    call r8spvec(ict, ivec, xvec, ivd, fval, 2*M+1, xpkg, fspl, iwarn, ier)
    if (ier.ne.0) then
       write(*,*) "error: spvec in TAECH_LIB.fspline."
       stop
    end if
    fimag = 0.0D0
    fimag(lftbc:rgbc) = fval(lftbc:rgbc,1)
    fs = cmplx(freal,fimag,kind=8)
    deallocate(xvec,stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for array xvec in TAECH_LIB.fspline"
       stop
    end if
    deallocate(fval,stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deabllocate memory for array fval in TAECH_LIB.fspline"
       stop
    end if
    return
  end subroutine fspline


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
    integer :: ibcxmin = 0
    integer :: ibcxmax = 0
    real(kind=8) :: bcxmin 
    real(kind=8) :: bcxmax 
    real(kind=8) :: wk(2*M0+1)
    integer :: ilinx = 2
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
       write(*,*) "error: genxpkg in  TAECH_LIB.conv_fgrid"
       stop
    end if
    ! real component of fs1 : fs1_real
    fspl(1,-M0:M0) = real(fs0(-M0:M0))
    call r8cspline(s0, 2*M0+1, fspl, ibcxmin, bcxmin, ibcxmax, bcxmax, wk, 2*M0+1, ilinx, ier)
    if (ier.ne.0) then
       write(*,*) "error: cspline in TAECH_LIB.conv_fgrid"
       stop
    end if
    call r8spvec(ict, ivec, xvec, ivd, fval, 2*M0+1, xpkg, fspl, iwarn, ier)
    if (ier.ne.0) then
       write(*,*) "error: spvec in TAECH_LIB.conv_fgrid"
       stop
    end if
    fs1_real(-M1:M1) = fval(-M1:M1,1)
    ! imaginary component of fs1 : fs1_imag
    fspl(1,-M0:M0) = aimag(fs0(-M0:M0))
    call r8cspline(s0, 2*M0+1, fspl, ibcxmin, bcxmin, ibcxmax, bcxmax, wk, 2*M0+1, ilinx, ier)
    if (ier.ne.0) then
       write(*,*) "error: cspline in TAECH_LIB.conv_fgrid"
       stop
    end if
    call r8spvec(ict, ivec, xvec, ivd, fval, 2*M0+1, xpkg, fspl, iwarn, ier)
    if (ier.ne.0) then
       write(*,*) "error: spvec in TAECH_LIB.conv_fgrid"
       stop
    end if
    fs1_imag(-M1:M1) = fval(-M1:M1,1)
    fs1(-M1:M1) = cmplx(fs1_real(-M1:M1), fs1_imag(-M1:M1), kind=8)
    return
  end subroutine conv_fgrid
end module TAECH_LIB
