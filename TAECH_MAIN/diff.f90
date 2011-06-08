program diff
  use TAECH_IO
  use TAECH_LIB
  implicit none
  real(kind=8) :: v0 
  real(kind=8) :: D 
  real(kind=8) :: dt, ds
  real(kind=8) :: smax
  integer :: M, modes 
  complex(kind=8), allocatable :: fsk(:,:)
  complex(kind=8), allocatable :: spfminus(:,:), spf0(:,:), spfplus(:,:)
  real(kind=8), allocatable :: s(:)
  real(kind=8), allocatable :: fvx(:,:)
  integer :: error = 0
  integer :: nt, p, q
  real(kind=8) :: incx
  v0 = 1.0D0
  dt = 0.1D0/v0
  D = 0.025D0*v0**3
  ds = dt
  smax = 10.0D0
  M = nint(smax/ds)
  modes = 20
  allocate(fsk(-M:M, 0:modes), stat=error)
  fsk = (0.0D0, 0.0D0)
  allocate(spfminus(-M:M, -modes-1:modes+1), stat=error)
  spfminus = (0.0D0, 0.0D0)
  allocate(spf0(-M:M, -modes-1:modes+1), stat=error)
  spf0 = (0.0D0, 0.0D0)
  allocate(spfplus(-M:M, -modes-1:modes+1), stat=error)
  spfplus = (0.0D0, 0.0D0)
  allocate(s(-M:M), stat=error)
  forall(p=-M:M) s(p) = p*ds
  allocate(fvx(-M:M,-modes-1:modes+1), stat=error)
  fvx = 0.0D0
  do nt = 1, 3000
     spf0(-M:M, 0:modes) = fsk(-M:M, 0:modes)
     spf0(-M:M, -modes:-1) = conjg(fsk(M:-M:-1, modes:1:-1))
     spf0(-M:M, -modes-1) = (0.0D0, 0.0D0)
     spf0(-M:M, modes+1) = (0.0D0, 0.0D0)
     spfminus(-M:M, -modes-1:modes+1) = spf0(-M:M, -modes-1:modes+1)
     spfplus(-M:M, -modes-1:modes+1) = spf0(-M:M, -modes-1:modes+1)
     do q = -modes, modes
        incx = -q*dt
        call fspline(spfminus(-M:M, q-1), s, M, incx)
        call fspline(spf0(-M:M, q), s, M, incx)
        call fspline(spfplus(-M:M, q+1), s, M, incx)
     end do
     do p = -M, M
        do q = 1, modes
           fsk(p, q) = spf0(p,q)*exp(D/q*((s(p)-0.5D0*q*dt)**3-s(p)**3)) &
                + dt*exp(-imagj*v0*(s(p)-0.5D0*q*dt))*exp(D/q*((s(p)-0.5D0*q*dt)**3-s(p)**3))
        end do
       fsk(p,0) = spf0(p,0)*exp(-3.0D0*s(p)**2*dt*D) &
            + dt*exp(-imagj*v0*s(p))*exp(-1.5D0*s(p)**2*dt*D)
     end do
     call IFInt2D(fvx, M, fsk, M, modes, s)
     call save_output('diffusion.dat', fvx, 2*M+1, 2*modes+3, .True.)
     call save_output('fsk.dat', fsk, 2*M+1, modes+1, .True.)
     write(*,*) nt
  end do
  deallocate(fsk, stat=error)
  deallocate(spfminus, stat=error)
  deallocate(spf0, stat=error)
  deallocate(spfplus, stat=error)
  deallocate(s, stat=error)
  deallocate(fvx, stat=error)
  return
end program diff
