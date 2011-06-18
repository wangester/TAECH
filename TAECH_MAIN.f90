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
  integer :: ngen = 100

  real(kind=8), allocatable :: s(:)
  complex(kind=8), allocatable :: fsk(:,:)
  complex(kind=8), allocatable :: wave_ampl(:)
  real(kind=8), allocatable :: xf(:), vf(:)
  complex(kind=8), allocatable :: fedge(:,:)

  complex(kind=8) :: fcorr(2) 
  real(kind=8) :: res_bound(2)
  logical :: is_locked 
  real(kind=8) :: elapsed
  integer :: M, edge
  complex(kind=8), allocatable :: kernel(:)
  real(kind=8), allocatable :: t(:)

contains
  subroutine initialize()
    implicit none
    integer :: error = 0
    integer :: i
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
    M = nint(ngen+(smax+sbnd-dt)/ds)
    edge = nint(ngen+(smax-dt)/ds)
    allocate(s(-M:M), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for s in TAECH_MAIN.initialize"
       stop
    end if
    forall(i=-ngen:ngen) s(i) = i*(dt/ngen)
    forall(i=-M:-ngen-1) s(i) = s(-ngen)+(i+ngen)*ds
    forall(i=ngen+1:M) s(i) = s(ngen)+(i-ngen)*ds
    allocate(fsk(-M:M,0:modes), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for fsk in TAECH_MAIN.initialize"
       stop
    end if
    fsk = (0.0D0, 0.0D0)
    allocate(fedge(1:modes, 0:tend), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for fedge in TAECH_MAIN.initialize"
       stop
    end if
    fedge = (0.0D0, 0.0D0)
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
       call get_input('./resume_job/fedge.dat', fedge, modes, tstart+1)
       call get_input('./resume_job/frame_acel.dat', fcorr, 2)
       call get_input('./resume_job/resonance_bound.dat', res_bound, 2)
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
    integer :: error = 0
    call save_output('./results/wave_ampl.dat', wave_ampl, tend+1, .False.)
    call save_output('./results/chirp_freq.dat', vf, tend+1, .False.)
    call save_output('./resume_job/wave_ampl.dat', wave_ampl, tend+1, .False.)
    call save_output('./resume_job/frame_pos.dat', xf, tend+1, .False.)
    call save_output('./resume_job/frame_vel.dat', vf, tend+1, .False.)
    call save_output('./resume_job/fsk.dat', fsk, 2*M+1, modes+1, .False.)
    call save_output('./resume_job/fedge.dat', fedge, modes, tend+1, .False.)
    call save_output('./resume_job/frame_acel.dat', fcorr, 2, .False.)
    call save_output('./resume_job/resonance_bound.dat', res_bound, 2, .False.)
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
    deallocate(fedge, stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for array fedge in TAECH_MAIN.finalize"
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
    res_bound(1) = resonance_bot
    res_bound(2) = resonance_top
    fcorr = (0.0D0, 0.0D0)
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
       if (p > 0 .and. p < ngen) then
          rhs(1) = rhs(1) + 0.5D0/dt*(s(p)*wave_ampl(0)+(dt-s(p))*wave_ampl(1))
       elseif (p == 0) then
          rhs(1) = rhs(1) + 0.25D0*wave_ampl(1)
          rhs(-1) = rhs(-1) + 0.25D0*conjg(wave_ampl(1))
       elseif (p == ngen) then 
          rhs(1) = rhs(1) + 0.25D0*wave_ampl(0)
       end if
       call solve_tridiag(sub_diag(-modes:modes), diag(-modes:modes), &
            sup_diag(-modes:modes), rhs(-modes:modes), deltaf(-modes:modes,p), 2*modes+1)
    end do
    fsk(0:M, 0:modes) = transpose(deltaf(0:modes, 0:M))+(0.0D0,0.0D0)
    fsk(-M:-1, 0:modes) = conjg(transpose(deltaf(0:-modes:-1, M:1:-1))+(0.0D0,0.0D0))
    forall(q=1:modes) fedge(q, 1) = fsk(edge, q)
    return
  end subroutine startup


  subroutine diagnose(nt)
    implicit none
    integer, intent(in) :: nt
    integer :: error = 0
    complex(kind=8) :: fsk_lab(-edge:edge,0:modes)
    real(kind=8), allocatable :: fvx_resonance(:,:)
    real(kind=8), allocatable :: fvx_lab_resonance(:,:)
    complex(kind=8) :: fv_wave(-edge:edge)
    real(kind=8) :: chirp_freq, chirp_ergy
    real(kind=8) :: lres, hres
    integer :: nlres, nhres
    real(kind=8) :: dv
    integer :: p, q
    if (mod(nt,tsave) == 0) then    
       lres = res_bound(1)-0.1D0
       nlres = nint(10*s(M)*lres/pi)
       hres = res_bound(2)+0.5D0
       nhres = nint(10*s(M)*hres/pi)
       allocate(fvx_resonance(0:nhres-nlres, -modes-1:modes+1), stat=error)
       if (error.ne.0) then
          write(*,*) "error: could not allocate memory for fvx_resonance in TAECH_MAIN.diagnose"
          stop
       end if
       call IFInt2D(fvx_resonance(0:nhres-nlres,-modes-1:modes+1), lres, hres, nhres-nlres, &
            fsk(-edge:edge,0:modes), edge, modes, s(-edge:edge))
       !call visualize(fvx_resonance(0:nhres-nlres,-modes-1:modes+1), lres, hres, nhres-nlres, &
       !     fedge(1:modes,0:nt), fsk(-edge:edge,0:modes), edge, modes, nt, s(-edge:edge))
       lres = res_bound(1)+vf(nt)-0.1D0
       nlres = nint(10*s(M)*lres/pi)
       hres = res_bound(2)+vf(nt)+0.5D0
       nhres = nint(10*s(M)*hres/pi)
       dv = (hres-lres)/(nhres-nlres)
       allocate(fvx_lab_resonance(0:nhres-nlres, -modes-1:modes+1), stat=error)
       if (error.ne.0) then
          write(*,*) "error: could not allocate memory for fvx_lab_resonance in TAECH_MAIN.diagnose"
          stop
       end if
       forall(p=-edge:edge, q=0:modes) fsk_lab(p,q) = fsk(p,q)*exp(-imagj*q*xf(nt)+imagj*s(p)*vf(nt))
       call IFInt2D(fvx_lab_resonance(0:nhres-nlres,-modes-1:modes+1), lres, hres, &
            nhres-nlres, fsk_lab(-edge:edge,0:modes), edge, modes, s(-edge:edge))
       !call visualize(fvx_lab_resonance(0:nhres-nlres,-modes-1:modes+1), lres, hres, nhres-nlres, &
       !     fedge(1:modes,0:nt), fsk_lab(-edge:edge,0:modes), edge, modes, nt, s(-edge:edge))
       forall(p=0:nhres-nlres,q=-modes-1:modes+1) fvx_lab_resonance(p,q) = fvx_lab_resonance(p,q) &
            + lres+p*dv
       call save_output('./results/fvx_resonance.dat', fvx_resonance, nhres-nlres+1, 2*modes+3, .True.)
       call save_output('./results/fvx_lab_resonance.dat', & 
            fvx_lab_resonance, nhres-nlres+1, 2*modes+3, .True.)
       deallocate(fvx_resonance, stat=error)
       if (error.ne.0) then
          write(*,*) "error: could not deallocate memory for fvx_resonance in TAECH_MAIN.diagnose"
          stop
       end if
       deallocate(fvx_lab_resonance, stat=error)
       if (error.ne.0) then
          write(*,*) "error: could not deallocate memory for fvx_lab_resonance in TAECH_MAIN.diagnose"
          stop
       end if
       if (is_debug) then
          call IFFT1D(fv_wave(-edge:edge), edge, fsk(-edge:edge,0), s(-edge:edge))         
          call save_output('./diagnostics/fv_wave.dat',fv_wave, 2*edge+1, .True. )
          call save_output('./diagnostics/wave_ampl.dat', wave_ampl(0:nt), nt+1, .False.)
          call save_output('./diagnostics/chirp_freq.dat', vf(0:nt), nt+1, .False.)
          call save_output('./diagnostics/fsk_wave.dat', fsk, 2*M+1, modes+1, .True.)
       end if
    end if
    write(*,*) nt, wave_ampl(nt) 
    if (nt>nint(100/dt) .and. is_debug) then
       call chirp(chirp_freq, chirp_ergy, wave_ampl(nt-nint(100/dt):nt), &
            nint(100/dt), res_bound(1), res_bound(2))
       write(*,*) chirp_freq, chirp_ergy
       write(*,*) resonant_density(fsk(-edge:edge,0), edge, s(-edge:edge), &
            res_bound(1), res_bound(2))
       write(*,*) resonant_position(fsk(-edge:edge,0:modes), &
            edge, modes, s(-edge:edge), res_bound(1), res_bound(2))
       write(*,*) resonant_velocity(fsk(-edge:edge,0), edge, &
            s(-edge:edge), res_bound(1), res_bound(2))
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
       where (position<=smax .and. position>=-smax)
          damping = cln*position**2
       end where
    else
       where (position<=smax .and. position>=-smax)
          damping = 0.0D0
       end where
    end if
    where (position>smax .and. position<=s(M))
       damping = dl0*((position-smax)/sbnd)**2
    elsewhere (position>=s(-M) .and. position<-smax)
       damping = dl0*((position+smax)/sbnd)**2
    elsewhere
       damping = 0.0D0
    end where
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


  subroutine track(nt)
    implicit none
    integer, intent(in) :: nt
    integer :: twin
    real(kind=8) :: chirp_freq, chirp_ergy
    if (.not. is_tracking) then 
       return
    end if
    twin = nint(100/dt)
    if (.not.is_locked .and. nt>=twin .and. &
         is_resonant(wave_ampl(nt-twin:nt), twin, fsk(-edge:edge,0), &
         edge, s(-edge:edge)))  then
       is_locked = .True.
    end if
    if (is_locked .and. nt>=nint(4000/dt) .and. &
         (vf(nt)<-pi/dt .or. resonant_density(fsk(-edge:edge,0), edge, s(-edge:edge), &
         res_bound(1), res_bound(2))<1.5D0*pi*(resonance_top-resonance_bot)**2))  then
       is_locked = .False.
    end if
    if (is_locked) then
       call chirp(chirp_freq, chirp_ergy, wave_ampl(nt-nint(100/dt):nt), &
            nint(100/dt), res_bound(1), res_bound(2))
       do while (chirp_ergy<0.1D0)
          call chirp(chirp_freq, chirp_ergy, wave_ampl(nt-nint(100/dt):nt), &
               nint(100/dt), res_bound(1), res_bound(2))
          res_bound(2) = res_bound(2)+0.01D0
       end do
       !res_bound(1) = 0.5D0*(resonance_bot-resonance_top) &
       !     + resonant_velocity(fsk(-edge:edge,0), &
       !     edge, s(-edge:edge), res_bound(1), res_bound(2))
       !res_bound(2) = res_bound(1)+resonance_top-resonance_bot
       fcorr(1) = resonant_filter(fsk(-edge:edge,0), fsk(-edge:edge,1), edge, &
            s(-edge:edge), res_bound(1), res_bound(2))
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
    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP PRIVATE(incx)
    do q = -modes, modes
       incx = -q*dt
       call fspline(spfminus(-M:M, q-1), s, M, incx)
       call fspline(spf0(-M:M, q), s, M, incx)
       call fspline(spfplus(-M:M, q+1), s, M, incx)
    end do
    !$OMP END PARALLEL DO
    ! solve the wave equation
    ! boundary condition and set damping when s=0 
    sft = 0.0
    call set_boundary(bnd0, modes, sft)
    call set_damping(damping0, modes, sft, is_damping)
    forall(q=-modes:modes) sft(q) = -q*dt
    call set_boundary(bnd1, modes, sft)
    call set_damping(damping1, modes, sft, is_damping)
    do q = -modes, modes
       deltaf(q,0) = 1.0D0/(1.0D0+0.5D0*dt*damping0(q)) &
            * (-imagj*0.25D0*dt*(-q*dt) &
            * bnd1(q)*wave_ampl(nt)*spfminus(0,q-1) &
            - 0.5D0*dt*(0.5D0*imagj*(-q*dt) &
            * bnd1(q)*(wave_ampl(nt)*fcorr(2) &
            + conjg(wave_ampl(nt))*fcorr(1))+damping1(q)+damping0(q))*spf0(0,q) &
            - imagj*0.25D0*dt*(-q*dt)*bnd1(q) &
            * conjg(wave_ampl(nt))*spfplus(0,q+1))
    end do
    forall(q=0:modes) fsk(0,q) = deltaf(q,0)+spf0(0,q)
    wave_ampl(nt+1) = 1.0D0/(Deltam-0.5D0*dt+imagj &
         * (pi*eta/(1.0D0+0.5D0*dt*damping0(1))-1.0D0)) &
         * (dt*exp(imagj*xf(nt+1)) &
         * zdotu(nt, wave_ampl(1:nt)*exp(-imagj*xf(1:nt)), 1, &
         kernel(1:nt)*exp(-eps*t(1:nt)), -1) &
         + (1.0D0+dt/(2.0D0*(Deltam-imagj))) &
         * deltap*kernel(nt+1)*exp(-eps*t(nt+1)+imagj*xf(nt+1)) &
         - imagj*4.0D0*pi*eta*fsk(0,1))
    ! solve the vlasov equation 
    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP PRIVATE(q, sft, bnd0, bnd1, damping0, damping1, sub_diag, diag, sup_diag, rhs)
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
            * (sft(q)*bnd1(q)*wave_ampl(nt)*spfminus(p,q-1) &
            + s(p)*bnd0(q)*wave_ampl(nt+1)*spf0(p,q-1)) &
            - 0.5D0*dt*(0.5D0*imagj*sft(q)*bnd1(q) &
            * (wave_ampl(nt)*fcorr(2)+conjg(wave_ampl(nt))*fcorr(1)) &
            + damping1(q)+0.5D0*imagj*s(p)*bnd0(q)*(wave_ampl(nt+1)*fcorr(2) &
            + conjg(wave_ampl(nt+1))*fcorr(1))+damping0(q))*spf0(p,q) &
            - imagj*0.25D0*dt*(sft(q)*bnd1(q)*conjg(wave_ampl(nt))*spfplus(p,q+1) &
            + s(p)*bnd0(q)*conjg(wave_ampl(nt+1))*spf0(p,q+1))
       if ( p > 0 .and. p < ngen) then
          rhs(1) = rhs(1)+0.5D0/dt*(s(p)*wave_ampl(nt)+(dt-s(p))*wave_ampl(nt+1))
       elseif (p == 0) then
          rhs(1) = rhs(1) + 0.25D0*wave_ampl(nt+1)
          rhs(-1) = rhs(-1) + 0.25D0*conjg(wave_ampl(nt+1))
       elseif (p == ngen) then 
          rhs(1) = rhs(1) + 0.25D0*wave_ampl(nt)
       end if
       call solve_tridiag(sub_diag(-modes:modes), diag(-modes:modes), &
            sup_diag(-modes:modes), rhs(-modes:modes), deltaf(-modes:modes,p), 2*modes+1)
    end do
    !$OMP END PARALLEL DO
    fsk(0:M, 0:modes) = transpose(deltaf(0:modes, 0:M))+spf0(0:M, 0:modes)
    fsk(-M:-1, 0:modes) = conjg(transpose(deltaf(0:-modes:-1, M:1:-1))+spf0(M:1:-1,0:-modes:-1))
    forall(q=1:modes) fedge(q, nt+1) = fsk(edge, q)
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
