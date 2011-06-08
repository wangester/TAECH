module fspline_mod
  use TAECH_IO
  implicit none
  private

  public :: fspline
contains
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
    real(kind=8) :: bcxmin 
    real(kind=8) :: bcxmax 
    real(kind=8) :: wk(2*M+1)
    integer :: ilinx
    !(r8spvec)
    integer :: ict(3)
    integer :: ivec
    integer :: ivd
    integer :: iwarn
    real(kind=8), allocatable :: xvec(:)
    real(kind=8), allocatable :: fval(:,:)
    integer :: i, lftbc, rgbc
    integer :: error
    call r8genxpkg(2*M+1, s, xpkg, iper, imsg, itol, ztol, ialg, ier)
    if (ier.ne.0) then
       write(*,*) "error: genxpkg in fspline."
       stop
    end if
    ! real component of fs: freal
    fspl(1,-M:M) = real(fs(-M:M))
    call r8cspline(s, 2*M+1, fspl, ibcxmin, bcxmin, ibcxmax, bcxmax, wk, 2*M+1, ilinx, ier)
    if (ier.ne.0) then
       write(*,*) "error: cspline in fspline."
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
       write(*,*) "error: could not allocate memory for array xvec in fspline"
       stop
    end if
    allocate(fval(lftbc:rgbc,1),stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for array fval in fspline"
       stop
    end if
    ict = 0
    ict(1) = 1
    xvec(lftbc:rgbc) = s(lftbc:rgbc)+incx
    call r8spvec(ict, ivec, xvec, ivd, fval, 2*M+1, xpkg, fspl, iwarn, ier)
    if (ier.ne.0) then
       write(*,*) "error: spvec in fspline."
       stop
    end if
    freal = 0.0D0
    freal(lftbc:rgbc) = fval(lftbc:rgbc,1)
    ! imaginary component of f: fimag
    fspl(1,-M:M) = aimag(fs(-M:M))
    call r8cspline(s, 2*M+1, fspl, ibcxmin, bcxmin, ibcxmax, bcxmax, wk, 2*M+1, ilinx, ier)
    if (ier.ne.0) then
       write(*,*) "error: cspline in fspline."
       stop
    end if
    call r8spvec(ict, ivec, xvec, ivd, fval, 2*M+1, xpkg, fspl, iwarn, ier)
    if (ier.ne.0) then
       write(*,*) "error: spvec in fspline."
       stop
    end if
    fimag = 0.0D0
    fimag(lftbc:rgbc) = fval(lftbc:rgbc,1)
    fs = cmplx(freal,fimag,kind=8)
    deallocate(xvec,stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for array xvec in fspline"
       stop
    end if
    deallocate(fval,stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deabllocate memory for array fval in fspline"
       stop
    end if
    return
  end subroutine fspline
end module fspline_mod


program test_fspline
  use fspline_mod
  use TAECH_IO
  implicit none
  real(kind=8) :: incx = -0.4D0
  real(kind=8) :: s(-100:100)
  complex(kind=8) :: fs(-100:100)
  integer :: i
  forall(i=-100:100) s(i) = (i/100.0D0)**3
  forall(i=-100:100) fs(i) = (1.0D0+imagj)*exp(-10.0*s(i)**2)
  call fspline(fs, s, 100, incx)
  call save_output('data.dat', fs, 201, .False.)
  call save_output('data1.dat', s, 201, .False.)
  return
end program test_fspline
