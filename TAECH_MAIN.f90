module TAECH_MAIN
  use TAECH_IO
  use TAECH_LIB
  use TAECH_TRACK
  use omp_lib
  implicit none
  private

  public :: initialize
  public :: finalize
  public :: diagnose
  public :: vlasov
  public :: track

  logical, public :: is_damping = .False. 
  logical, public :: is_debug = .True.
  logical, public :: is_tracking = .True.

  real(kind=8), allocatable :: s(:)
  complex(kind=8), allocatable :: fsk(:,:)
  complex(kind=8), allocatable :: wave_ampl(:)
  real(kind=8), allocatable :: xf(:), vf(:)

  complex(kind=8) :: fcorr(2) 
  logical :: is_locked 
  real(kind=8) :: elapsed
  integer :: M
  complex(kind=8), allocatable :: kernel(:)
  real(kind=8), allocatable :: t(:)

contains
  subroutine initialize()
    implicit none
    integer :: error = 0
    integer :: i
    integer :: gen
    elapsed = OMP_get_wtime()
    call load_cfg('parameters.cfg')
    allocate(t(0:tend), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for t in TAECH_MAIN.initialize"
       stop
    end if
    forall(i=0:tend) t(i) = i*dt
    allocate(kernel(0:tend), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for kernel in TAECH_MAIN.initialize"
       stop
    end if
    kernel = besselj(t)
    M = nint(dt/(0.1D0*ds)+(smax+sbnd-dt)/ds)
    allocate(s(-M:M), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for s in TAECH_MAIN.initialize"
       stop
    end if
    gen = nint(dt/(0.1D0*ds))
    forall(i=-gen:gen) s(i) = i*(0.1D0*ds)
    forall(i=-M:-gen-1) s(i) = s(-gen)+(i+gen)*ds
    forall(i=gen+1:M) s(i) = s(gen)+(i-gen)*ds
    allocate(fsk(-M:M,0:modes), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for fsk in TAECH_MAIN.initialize"
       stop
    end if
    fsk = (0.0D0, 0.0D0)
    allocate(wave_ampl(0:tend), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for wave_ampl in TAECH_MAIN.initialize"
       stop
    end if
    wave_ampl = (0.0D0, 0.0D0)
    wave_ampl(0) = deltap/(Deltam-imagj)
    allocate(xf(0:tend), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for xf in TAECH_MAIN.initialize"
       stop
    end if
    xf = 0.0D0
    allocate(vf(0:tend), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for vf in TAECH_MAIN.initialize"
       stop
    end if
    vf = 0.0D0
    if (tstart == 0) then
       call startup()
    ! resume the job when tstart/=0
    else
       call get_input('./resume_job/wave_ampl.dat', wave_ampl, tstart+1)
       call get_input('./resume_job/frame_pos.dat', xf, tstart+1)
       call get_input('./resume_job/frame_vel.dat', vf, tstart+1)
       call get_input('./resume_job/fsk.dat', fsk, 2*M+1, modes+1)
       call get_input('./resume_job/frame.dat', fcorr, 2)
       if (abs(fcorr(1)) .ne. 0.0D0) then 
          is_locked = .True.
       else
          is_locked = .False.
       end if
    end if
    return
  end subroutine initialize

  
  subroutine finalize()
    implicit none
    integer :: error
    call save_output('./results/wave_ampl.dat', wave_ampl, tend+1, .False.)
    call save_output('./results/chirp_freq.dat', vf, tend+1, .False.)
    call save_output('./resume_job/wave_ampl.dat', wave_ampl, tend+1, .False.)
    call save_output('./resume_job/frame_pos.dat', xf, tend+1, .False.)
    call save_output('./resume_job/frame_vel.dat', vf, tend+1, .False.)
    call save_output('./resume_job/fsk.dat', fsk, 2*M+1, modes+1, .False.)
    call save_output('./resume_job/frame.dat', fcorr, 2, .False.)
    deallocate(s, stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for array s in TAECH_MAIN.finalize"
       stop
    end if
    deallocate(fsk, stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for array fsk in TAECH_MAIN.finalize"
       stop
    end if
    deallocate(wave_ampl, stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for array wave_ampl in TAECH_MAIN.finalize"
       stop
    end if
    deallocate(xf, stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for array xf in TAECH_MAIN.finalize"
       stop
    end if
    deallocate(vf, stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for array vf in TAECH_MAIN.finalize"
       stop
    end if
    deallocate(kernel, stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for array kernel in TAECH_MAIN.finalize"
       stop
    end if
    deallocate(t, stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for array t in TAECH_MAIN.finalize"
       stop
    end if
    elapsed = OMP_get_wtime()-elapsed
    write(*,*) "The code elapses ", elapsed, " sec."
    return
  end subroutine finalize

  
  subroutine diagnose(nt)
    implicit none
    integer, intent(in) :: nt
    real(kind=8) :: fvx_wave(-M:M, -modes-1:modes+1)
    real(kind=8) :: fvx_lab(-M:M, -modes-1:modes+1)
    real(kind=8) :: f_fvx_lab(-M:M, -modes-1:modes+1)
    real(kind=8) :: fvx_resonance(-10*M:10*M, -modes-1:modes+1)
    complex(kind=8) :: fsk_lab(-M:M, 0:modes)
    complex(kind=8) :: fv_wave(-M:M)
    real(kind=8) :: chirp_freq, chirp_ergy
    real(kind=8) :: v(-M:M)
    real(kind=8) :: dv
    real(kind=8) :: lres, hres
    integer :: nlres, nhres
    integer :: p, q
    if (mod(nt,tsave) == 0) then    
       forall(p=-M:M, q=0:modes) fsk_lab(p,q) = fsk(p,q)*exp(-imagj*q*xf(nt)+imagj*s(p)*vf(nt))
       call IFInt2D(fvx_wave, M, fsk, M, modes, s)
       call IFInt2D(fvx_lab, M, fsk_lab, M, modes, s)
       call IFInt2D(fvx_resonance, 10*M, fsk_lab, M, modes, s)
       lres = resonance_bot+vf(nt)-0.5D0
       nlres = nint(10*s(M)*lres/pi)
       hres = resonance_top+vf(nt)+0.5D0
       nhres = nint(10*s(M)*hres/pi)
       dv = pi/s(M)
       forall(p=-M:M) v(p) = p*dv
       forall(p=-M:M) f_fvx_lab(p,-modes-1:modes+1) = fvx_lab(p,-modes-1:modes+1)+v(p)
       call save_output('./results/fvx_wave.dat', fvx_wave, 2*M+1, 2*modes+3, .True.)
       call save_output('./results/fvx_lab.dat', fvx_lab, 2*M+1, 2*modes+3, .True.)
       call save_output('./results/f_fvx_lab.dat', f_fvx_lab, 2*M+1, 2*modes+3, .True.)
       call save_output('./results/fvx_resonance.dat', & 
            fvx_resonance(nlres:nhres,:), nhres-nlres+1, 2*modes+3, .True.)
       if (is_debug) then
          call IFFT1D(fv_wave(-M:M), M, fsk(-M:M,0), s)         
          call save_output('./diagnostics/fv_wave.dat',fv_wave, 2*M+1, .True. )
          call save_output('./diagnostics/wave_ampl.dat', wave_ampl(0:nt), nt+1, .False.)
          call save_output('./diagnostics/chirp_freq.dat', vf(0:nt), nt+1, .False.)
          call save_output('./diagnostics/fsk_wave.dat', fsk, 2*M+1, modes+1, .True.)
          call save_output('./diagnostics/fsk_lab.dat', fsk_lab, 2*M+1, modes+1, .True.)
       end if
    end if
    write(*,*) nt, wave_ampl(nt) 
    if (nt>nint(100/dt) .and. is_debug) then
       call chirp(chirp_freq, chirp_ergy, wave_ampl(nt-nint(100/dt):nt), &
            nint(100/dt), resonance_bot, resonance_top)
       write(*,*) chirp_freq, chirp_ergy
       write(*,*) resonant_density(fsk(-M:M,0), M, s)
       write(*,*) resonant_position(fsk, M, modes, s)
       write(*,*) resonant_velocity(fsk(-M:M,0), M, s)
    end if
    return
  end subroutine diagnose


  subroutine set_damping(damping, N, position, is_damping)
    implicit none
    integer, intent(in) :: N
    real(kind=8), intent(out) :: damping(-N:N)
    real(kind=8), intent(in) :: position(-N:N)
    logical, intent(in) :: is_damping
    if (is_damping) then
       damping = cln*position**2
    else
       where (position<=smax .and. position>=-smax)
          damping = 0.0D0
       elsewhere (position>smax .and. position<=s(M))
          damping = dl0*((position-smax)/sbnd)**2
       elsewhere (position>=s(-M) .and. position<-smax)
          damping = dl0*((position+smax)/sbnd)**2
       elsewhere
          damping = 0.0D0
       end where
    end if
    return
  end subroutine set_damping


  subroutine set_boundary(bnd, N, position)
    implicit none
    integer, intent(in) :: N
    real(kind=8), intent(out) :: bnd(-N:N)
    real(kind=8), intent(in) :: position(-N:N)
    where (position<=smax .and. position>=-smax)  
       bnd = 1.0D0
    elsewhere (position>smax .and. position<=s(M)) 
       bnd = exp(-dl1*(position-smax)/sbnd)
    elsewhere (position>=s(-M) .and. position<-smax) 
       bnd = exp(-dl1*(-position-smax)/sbnd)
    elsewhere
       bnd = 0.0D0
    end where
    return
  end subroutine set_boundary


  subroutine startup()
    implicit none
    complex(kind=8) :: deltaf(-modes:modes,0:M)
    complex(kind=8) :: sub_diag(-modes:modes), diag(-modes:modes), sup_diag(-modes:modes)
    complex(kind=8) :: rhs(-modes:modes)
    real(kind=8) :: bnd0(-modes:modes), bnd1(-modes:modes)
    real(kind=8) :: damping0(-modes:modes), damping1(-modes:modes)
    real(kind=8) :: sft(-modes:modes)
    integer :: p, q
    is_locked = .False.
    sft = 0.0D0
    call set_damping(damping0, modes, sft, is_damping)
    wave_ampl(1) = 1.0D0/(Deltam-0.5D0*dt+imagj*(pi*eta/(1.0D0+0.5D0*dt*damping0(1))-1.0D0)) &
         * ((1.0D0+dt/(2.0D0*(Deltam-imagj)))*deltap*kernel(1)*exp(-eps*t(1)))
    do p = 0, M
       forall(q=-modes:modes) sft(q) = s(p)
       call set_boundary(bnd0, modes, sft)
       call set_damping(damping0, modes, sft, is_damping)
       forall(q=-modes:modes) sft(q) = s(p)-q*dt
       call set_boundary(bnd1, modes, sft)
       call set_damping(damping1, modes, sft, is_damping)
       forall(q=-modes:modes) sub_diag(q) = imagj*0.25D0*dt*s(p)*bnd0(q)*wave_ampl(1)
       forall(q=-modes:modes) diag(q) = 1.0D0 + 0.5D0*dt*damping0(q)
       forall(q=-modes:modes) sup_diag(q) = imagj*0.25D0*dt*s(p)*bnd0(q)*conjg(wave_ampl(1))
       forall(q=-modes:modes) rhs(q) = (0.0D0, 0.0D0)
       if (s(p) > 0.5D0*(0.1D0*ds) .and. s(p) < dt-0.5D0*(0.1D0*ds)) then
          rhs(1) = rhs(1) + 0.5D0/dt*(s(p)*wave_ampl(0)+(dt-s(p))*wave_ampl(1))
       elseif (p == 0) then
          rhs(1) = rhs(1) + 0.25D0*wave_ampl(1)
          rhs(-1) = rhs(-1) + 0.25D0*conjg(wave_ampl(1))
       elseif (abs(s(p)-dt) <= 0.5D0*(0.1D0*ds)) then 
          rhs(1) = rhs(1) + 0.25D0*wave_ampl(0)
       end if
       call solve_tridiag(sub_diag(-modes:modes), diag(-modes:modes), &
            sup_diag(-modes:modes), rhs(-modes:modes), deltaf(-modes:modes,p), 2*modes+1)
    end do
    fsk(0:M, 0:modes) = transpose(deltaf(0:modes, 0:M))+(0.0D0,0.0D0)
    fsk(-M:-1, 0:modes) = conjg(transpose(deltaf(0:-modes:-1, M:1:-1))+(0.0D0,0.0D0))
    return
  end subroutine startup
 

  subroutine track(nt)
    implicit none
    integer, intent(in) :: nt
    integer :: twin
    if (.not. is_tracking) then 
       return
    end if
    twin = nint(100/dt)
    if (.not.is_locked .and. nt>=twin .and. &
         is_resonant(wave_ampl(nt-twin:nt), twin, fsk(-M:M,0), M, s))  then
       is_locked = .True.
    end if
    if (is_locked .and. nt>=nint(4000/dt) .and. &
         (vf(nt)<-pi/dt .or. resonant_density(fsk(-M:M,0), M, s) &
         < 1.5D0*pi*(resonance_top-resonance_bot)**2))  then
       is_locked = .False.
    end if
    if (is_locked) then
       fcorr(1) = resonant_filter(fsk(-M:M,0), fsk(-M:M,1), M, s)
       fcorr(2) = conjg(fcorr(1))
       vf(nt+1) = vf(nt)+0.5D0*dt*real(wave_ampl(nt)*fcorr(2)+conjg(wave_ampl(nt))*fcorr(1))
       xf(nt+1) = xf(nt)+0.5D0*dt*(vf(nt)+vf(nt+1)) 
    else
       fcorr = (0.0D0, 0.0D0)
       vf(nt+1) = vf(nt)
       xf(nt+1) = xf(nt)+0.5D0*dt*(vf(nt)+vf(nt+1))
    end if
    return
  end subroutine track


  subroutine vlasov(nt)
    implicit none
    integer, intent(in) :: nt
    complex(kind=8) :: deltaf(-modes:modes,0:M)
    complex(kind=8) :: sub_diag(-modes:modes), diag(-modes:modes), sup_diag(-modes:modes)
    complex(kind=8) :: rhs(-modes:modes)
    real(kind=8) :: bnd0(-modes:modes), bnd1(-modes:modes)
    real(kind=8) :: damping0(-modes:modes), damping1(-modes:modes)
    real(kind=8) :: sft(-modes:modes)
    complex(kind=8) :: spfminus(-M:M, -modes-1:modes+1)
    complex(kind=8) :: spf0(-M:M, -modes-1:modes+1)
    complex(kind=8) :: spfplus(-M:M, -modes-1:modes+1)
    real(kind=8) :: incx
    integer :: p, q
    ! calculate the distribution at old time step along the characteristic line 
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
    ! boundary condition and set damping when s=0 
    sft = 0.0
    call set_boundary(bnd0, modes, sft)
    call set_damping(damping0, modes, sft, is_damping)
    forall(q=-modes:modes) sft(q) = -q*dt
    call set_boundary(bnd1, modes, sft)
    call set_damping(damping1, modes, sft, is_damping)
    do q = -modes, modes
       deltaf(q,0) = 1.0D0/(1.0D0+0.5D0*dt*damping0(q)) &
            *( -imagj*0.25D0*dt*(-q*dt)*bnd1(q)*wave_ampl(nt)*spfminus(0,q-1) &
               -0.5D0*dt*(0.5D0*imagj*(-q*dt)*bnd1(q)*(wave_ampl(nt)*fcorr(2) &
               +conjg(wave_ampl(nt))*fcorr(1))+damping1(q)+damping0(q))*spf0(0,q) &
               -imagj*0.25D0*dt*(-q*dt)*bnd1(q)*conjg(wave_ampl(nt))*spfplus(0,q+1) )
    end do
    forall(q=0:modes) fsk(0,q) = deltaf(q,0) + spf0(0,q)
    ! solve the wave equation
    wave_ampl(nt+1) = 1.0D0/(Deltam-0.5D0*dt+imagj*(pi*eta/(1.0D0+0.5D0*dt*damping0(1))-1.0D0)) &
         * (dt*exp(imagj*xf(nt+1))*zdotu(nt, wave_ampl(1:nt)*exp(-imagj*xf(1:nt)), 1, &
         kernel(1:nt)*exp(-eps*t(1:nt)), -1)+(1.0D0+dt/(2.0D0*(Deltam-imagj))) &
         * deltap*kernel(nt+1)*exp(-eps*t(nt+1)+imagj*xf(nt+1)) &
         - imagj*4.0D0*pi*eta*fsk(0,1))
    ! solve the vlasov equation 
    do p = 0, M
       forall(q=-modes:modes) sft(q) = s(p)
       call set_boundary(bnd0, modes, sft)
       call set_damping(damping0, modes, sft, is_damping)
       forall(q=-modes:modes) sft(q) = s(p)-q*dt
       call set_boundary(bnd1, modes, sft)
       call set_damping(damping1, modes, sft, is_damping)
       forall(q=-modes:modes) sub_diag(q) = imagj*0.25D0*dt*s(p)*bnd0(q)*wave_ampl(nt+1)
       forall(q=-modes:modes) diag(q) = 1.0D0+imagj*0.25D0*dt*s(p)*bnd0(q) &
            * (wave_ampl(nt+1)*fcorr(2)+conjg(wave_ampl(nt+1))*fcorr(1)) &
            + 0.5D0*dt*damping0(q)
       forall(q=-modes:modes) sup_diag(q) = imagj*0.25D0*dt*s(p)*bnd0(q)*conjg(wave_ampl(nt+1))
       forall(q=-modes:modes) rhs(q) = -imagj*0.25D0*dt &
            * (sft(q)*bnd1(q)*wave_ampl(nt)*spfminus(p,q-1)+s(p)*bnd0(q)*wave_ampl(nt+1)*spf0(p,q-1)) &
            - 0.5D0*dt*(0.5D0*imagj*sft(q)*bnd1(q) &
            * (wave_ampl(nt)*fcorr(2)+conjg(wave_ampl(nt))*fcorr(1)) &
            + damping1(q)+0.5D0*imagj*s(p)*bnd0(q)*(wave_ampl(nt+1)*fcorr(2) &
            + conjg(wave_ampl(nt+1))*fcorr(1))+damping0(q))*spf0(p,q) &
            - imagj*0.25D0*dt*(sft(q)*bnd1(q)*conjg(wave_ampl(nt))*spfplus(p,q+1) &
            + s(p)*bnd0(q)*conjg(wave_ampl(nt+1))*spf0(p,q+1))
       if (s(p) > 0.5D0*(0.1D0*ds) .and. s(p) < dt-0.5D0*(0.1D0*ds)) then
          rhs(1) = rhs(1) + 0.5D0/dt*(s(p)*wave_ampl(nt)+(dt-s(p))*wave_ampl(nt+1))
       elseif (p == 0) then
          rhs(1) = rhs(1) + 0.25D0*wave_ampl(nt+1)
          rhs(-1) = rhs(-1) + 0.25D0*conjg(wave_ampl(nt+1))
       elseif (abs(s(p)-dt) <= 0.5D0*(0.1D0*ds)) then 
          rhs(1) = rhs(1) + 0.25D0*wave_ampl(nt)
       end if
       call solve_tridiag(sub_diag(-modes:modes), diag(-modes:modes), &
            sup_diag(-modes:modes), rhs(-modes:modes), deltaf(-modes:modes,p), 2*modes+1)
    end do
    fsk(0:M, 0:modes) = transpose(deltaf(0:modes, 0:M))+spf0(0:M, 0:modes)
    fsk(-M:-1, 0:modes) = conjg(transpose(deltaf(0:-modes:-1, M:1:-1))+spf0(M:1:-1,0:-modes:-1))
    return
  end subroutine vlasov
end module TAECH_MAIN


program TAECH
  use TAECH_IO
  use TAECH_LIB
  use TAECH_TRACK
  use TAECH_MAIN
  implicit none
  integer :: nt 
  call initialize()
  do nt = max(1,tstart), tend-1
     call track(nt)
     call vlasov(nt)
     call diagnose(nt+1)
  end do
  call finalize()
  return
end program TAECH
