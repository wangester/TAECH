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
  public :: predictor
  public :: corrector

  real(kind=8), allocatable :: s(:)
  complex(kind=8), allocatable :: fsk(:,:)
  complex(kind=8), allocatable :: wave_ampl(:)
  real(kind=8), allocatable :: xf(:), vf(:), frame_boost(:)

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
    M = nint(dt/ds+(smax+sbnd)/dt-1)
    allocate(s(-M:M), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for s in TAECH_MAIN.initialize"
       stop
    end if
    gen = nint(dt/ds)
    forall(i=-gen:gen) s(i) = i*ds
    forall(i=-M:-gen-1) s(i) = s(-gen)+(i+gen)*dt
    forall(i=gen+1:M) s(i) = s(gen)+(i-gen)*dt
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
    allocate(frame_boost(0:tend), stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not allocate memory for frame_boost in TAECH_MAIN.initialize"
       stop
    end if
    frame_boost = 0.0D0
    if (t0 == 0) then
       call startup()
    ! resume the job when t0/=0
    else
       call get_input('./resume_job/wave_ampl.dat', wave_ampl, t0+1)
       call get_input('./resume_job/frame_pos.dat', xf, t0+1)
       call get_input('./resume_job/frame_vel.dat', vf, t0+1)
       call get_input('./resume_job/frame_accel.dat', frame_boost, t0+1)
       call get_input('./resume_job/fsk.dat', fsk, 2*M+1, modes+1)
    end if
    return
  end subroutine initialize

  
  subroutine finalize()
    implicit none
    integer :: error
    call save_output('./results/wave_ampl.dat', wave_ampl, tend+1, .False.)
    call save_output('./results/chirp_rate.dat', frame_boost, tend+1, .False.)
    call save_output('./results/chirp_freq.dat', vf, tend+1, .False.)
    call save_output('./resume_job/wave_ampl.dat', wave_ampl, tend+1, .False.)
    call save_output('./resume_job/frame_pos.dat', xf, tend+1, .False.)
    call save_output('./resume_job/frame_vel.dat', vf, tend+1, .False.)
    call save_output('./resume_job/frame_accel.dat', frame_boost, tend+1, .False.)
    call save_output('./resume_job/fsk.dat', fsk, 2*M+1, modes+1, .False.)
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
    deallocate(frame_boost, stat=error)
    if (error.ne.0) then
       write(*,*) "error: could not deallocate memory for array frame_boost in TAECH_MAIN.finalize"
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
    real(kind=8) :: v(-M:M)
    real(kind=8) :: dv
    real(kind=8) :: lres, hres
    integer :: nlres, nhres
    integer :: p, q
    if (mod(nt,tsave) == 0) then    
       forall(p=-M:M, q=0:modes) fsk_lab(p,q) = fsk(p,q)*exp(-imagj*q*xf(nt)+imagj*s(p)*vf(nt))
       call IFInt2D(fvx_wave(-M:M, -modes-1:modes+1), M, fsk(-M:M,0:modes), M, modes, s(-M:M))
       call IFInt2D(fvx_lab, M, fsk_lab, M, modes, s)
       call IFInt2D(fvx_resonance, 10*M, fsk_lab, M, modes, s)
       call IFFT1D(fv_wave(-M:M), M, fsk(-M:M,0), s)
       lres = v0+vf(nt)
       nlres = nint(10*s(M)*lres/pi)
       hres = v1+vf(nt)
       nhres = nint(10*s(M)*hres/pi)
       dv = pi/s(M)
       forall(p=-M:M) v(p) = p*dv
       forall(p=-M:M) f_fvx_lab(p,-modes-1:modes+1) = fvx_lab(p,-modes-1:modes+1)+v(p)
       call save_output('./diagnostics/fv_wave.dat',fv_wave, 2*M+1, .True. )
       call save_output('./diagnostics/wave_ampl.dat', wave_ampl(0:nt), nt+1, .False.)
       call save_output('./diagnostics/chirp_rate.dat', frame_boost(0:nt), nt+1, .False.)
       call save_output('./diagnostics/chirp_freq.dat', vf(0:nt), nt+1, .False.)
       call save_output('./diagnostics/fsk_wave.dat', fsk, 2*M+1, modes+1, .True.)
       call save_output('./diagnostics/fsk_lab.dat', fsk_lab, 2*M+1, modes+1, .True.)
       call save_output('./results/fvx_wave.dat', fvx_wave, 2*M+1, 2*modes+3, .True.)
       call save_output('./results/fvx_lab.dat', fvx_lab, 2*M+1, 2*modes+3, .True.)
       call save_output('./results/f_fvx_lab.dat', f_fvx_lab, 2*M+1, 2*modes+3, .True.)
       call save_output('./results/fvx_resonance.dat', & 
            fvx_resonance(nlres:nhres,:), nhres-nlres+1, 2*modes+3, .True.)
    end if
    write(*,*) nt, wave_ampl(nt) 
    !write(*,*) nt, resonant_density(fsk(-M:M,0), M, s)
    return
  end subroutine diagnose


  subroutine boundary(bnd1, bnd2, N, sgrid)
    implicit none
    integer, intent(in) :: N
    real(kind=8), intent(in) :: sgrid(-N:N)
    real(kind=8), intent(out) :: bnd1(-N:N)
    real(kind=8), intent(out) :: bnd2(-N:N)
    where (sgrid<=smax .and. sgrid>=-smax)  
       bnd1 = 1.0D0
       bnd2 = 0.0D0
    elsewhere (sgrid>smax .and. sgrid<=s(M)) 
       bnd1 = exp(-cln1*(sgrid-smax)/sbnd)
       bnd2 = 3.0D0*cln2*(sgrid-smax)**2/sbnd**2
    elsewhere (sgrid>=s(-M) .and. sgrid<-smax) 
       bnd1 = exp(-cln1*(-sgrid-smax)/sbnd)
       bnd2 = 3.0D0*cln2*(sgrid+smax)**2/sbnd**2
    elsewhere
       bnd1 = 0.0D0
       bnd2 = 0.0D0
    end where
    return
  end subroutine boundary


  subroutine startup()
    implicit none
    ! predictor
    wave_ampl(1) = 1.0D0/(Deltam-0.5D0*dt+imagj*(pi*eta-1.0D0)) &
         * ((1.0D0+dt/(2.0D0*(Deltam-imagj)))*deltap*kernel(1)*exp(-eps*t(1)))
    call vlasov(0)
    ! corrector
    wave_ampl(1) = 1.0D0/(Deltam-0.5D0*dt-imagj) &
         * ((1.0D0+dt/(2.0D0*(Deltam-imagj)))*deltap*kernel(1)*exp(-eps*t(1)) &
         - imagj*4.0D0*pi*eta*fsk(0,1))
    call vlasov(0)
    return
  end subroutine startup
 
  
  subroutine track(nt)
    implicit none
    integer, intent(in) :: nt
    logical :: locked=.False.
    integer :: twin
    twin = nint(100/dt)
    if (.not.locked) then
       if (nt>=twin .and. &
            abs(resonant_density(fsk(-M:M,0), M, s))>1.0D-2 .and. &
            is_resonant(wave_ampl(nt-twin+1:nt), twin) .and. &
            is_tracking(fsk(-M:M,0), M, s)) then
          locked = .True.
       end if
    end if
    if (locked) then   
       frame_boost(nt) = frame_accel(fsk(-M:M,1), M, wave_ampl(nt), s) &
            /resonant_density(fsk(-M:M,0), M, s)
       vf(nt+1) = vf(nt)+frame_boost(nt)*dt
       xf(nt+1) = xf(nt)+dt*(vf(nt)+vf(nt+1))/2.0D0
    else
       frame_boost(nt) = 0.0D0
       vf(nt+1) = 0.0D0
       xf(nt+1) = 0.0D0
    end if
    !IF (nt>twin) then
       !WRITE(*,*) resonant_velocity(fsk(-M:M,0), M, s)/resonant_density(fsk(-M:M,0), M, s)
       !CALL resonant_peak(resonance, chirp_freq, wave_ampl(nt-twin+1:nt), twin, v0, v1)
       !WRITE(*,*) resonance, chirp_freq
       !WRITE(*,*) is_resonant(wave_ampl(nt-twin+1:nt), twin), is_tracking(fsk(-M:M,0), M, s)
       !WRITE(*,*) resonant_density(fsk(-M:M,0), M, s)
    !END IF
    return
  end subroutine track


  subroutine vlasov(nt)
    implicit none
    integer, intent(in) :: nt
    complex(kind=8) :: deltaf(-modes:modes,0:M)
    complex(kind=8) :: spfminus(-M:M, -modes-1:modes+1)
    complex(kind=8) :: spf0(-M:M, -modes-1:modes+1)
    complex(kind=8) :: spfplus(-M:M, -modes-1:modes+1)
    complex(kind=8) :: sub_diag(-modes:modes), diag(-modes:modes), sup_diag(-modes:modes)
    complex(kind=8) :: rhs(-modes:modes)
    complex(kind=8) :: wave_mid
    real(kind=8) :: alpha
    real(kind=8) :: incx
    real(kind=8) :: bnd1(-modes:modes), bnd2(-modes:modes)
    real(kind=8) :: sft(-modes:modes)
    integer :: p, q
    wave_mid = 0.5D0*(wave_ampl(nt)+wave_ampl(nt+1))
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
    do p = 0, M
       alpha = frame_boost(nt) 
       forall(q=-modes:modes) sft(q) = s(p)-0.5D0*q*dt
       call boundary(bnd1, bnd2, modes, sft)
       forall(q=-modes:modes) sub_diag(q) = imagj*dt/4.0D0*sft(q)*bnd1(q)*wave_mid
       forall(q=-modes:modes) diag(q) = 1.0D0+imagj*dt/2.0D0 &
            * (sft(q)*bnd1(q)*alpha-imagj*bnd2(q))
       forall(q=-modes:modes) sup_diag(q) = imagj*dt/4.0D0*sft(q)*bnd1(q)*conjg(wave_mid)
       forall(q=-modes:modes) rhs(q) = -imagj*dt*sft(q)*bnd1(q) &
            * (wave_mid*(spfminus(p,q-1)+spf0(p,q-1))/4.0D0 &
            + conjg(wave_mid)*(spfplus(p,q+1)+spf0(p,q+1))/4.0D0) &
            - imagj*dt*(sft(q)*bnd1(q)*alpha-imagj*bnd2(q))*spf0(p,q)
       if (s(p) > 0.5D0*ds .and. s(p) < dt-0.5D0*ds) then
          rhs(1) = rhs(1) + 0.5D0/dt*(s(p)*wave_ampl(nt)+(dt-s(p))*wave_ampl(nt+1))
       elseif (p == 0) then
          rhs(1) = rhs(1) + 0.25D0*wave_ampl(nt+1)
          rhs(-1) = rhs(-1) + 0.25D0*conjg(wave_ampl(nt+1))
       elseif (abs(s(p)-dt) <= 0.5D0*ds) then 
          rhs(1) = rhs(1) + 0.25D0*wave_ampl(nt)
       end if
       call solve_tridiag(sub_diag(-modes:modes), diag(-modes:modes), &
            sup_diag(-modes:modes), rhs(-modes:modes), deltaf(-modes:modes,p), 2*modes+1)
    end do
    fsk(0:M, 0:modes) = transpose(deltaf(0:modes, 0:M))+spf0(0:M, 0:modes)
    fsk(-M:-1, 0:modes) = conjg(transpose(deltaf(0:-modes:-1, M:1:-1))+spf0(M:1:-1,0:-modes:-1))
    return
  end subroutine vlasov


  subroutine predictor(nt)
    implicit none
    integer, intent(in) :: nt
    complex(kind=8) :: f0, f1, f2
    !call track(nt)
    f0 = interpolate(-dt, fsk(-M:M,0), s(-M:M), 2*M+1)
    f1 = interpolate(-dt, fsk(-M:M,1), s(-M:M), 2*M+1)
    f2 = interpolate(-dt, fsk(-M:M,2), s(-M:M), 2*M+1)
    wave_ampl(nt+1) = 1.0D0/(Deltam-0.5D0*dt+imagj*(pi*eta-1.0D0)) &
         * (dt*exp(imagj*xf(nt+1))*zdotu(nt, wave_ampl(1:nt)*exp(-imagj*xf(1:nt)), 1, &
         kernel(1:nt)*exp(-eps*t(1:nt)), -1)+(1.0D0+dt/(2.0D0*(Deltam-imagj))) &
         * deltap*kernel(nt+1)*exp(-eps*t(nt+1)+imagj*xf(nt+1)) &
         + 2.0D0*pi*eta*(dt**2*frame_boost(nt)-2.0D0*imagj)*f1 &
         + pi*eta*dt**2*(wave_ampl(nt)*f0+conjg(wave_ampl(nt))*f2))
    call vlasov(nt)
    return
  end subroutine predictor


  subroutine corrector(nt)
    implicit none
    integer, intent(in) :: nt
    wave_ampl(nt+1) = 1.0D0/(Deltam-0.5D0*dt-imagj) &
         * (dt*exp(imagj*xf(nt+1))*zdotu(nt, wave_ampl(1:nt)*exp(-imagj*xf(1:nt)), 1, &
         kernel(1:nt)*exp(-eps*t(1:nt)), -1)+(1.0D0+dt/(2.0D0*(Deltam-imagj))) &
         * deltap*kernel(nt+1)*exp(-eps*t(nt+1)+imagj*xf(nt+1)) &
         - imagj*4.0D0*pi*eta*fsk(0,1))
    !call vlasov(nt)
    return
  end subroutine corrector
end module TAECH_MAIN


program TAECH
  use TAECH_IO
  use TAECH_LIB
  use TAECH_TRACK
  use TAECH_MAIN
  implicit none
  integer :: nt 
  call initialize()
  do nt = t0+1, tend-1
     call diagnose(nt)
     call predictor(nt)
     call corrector(nt)
  end do
  call finalize()
  return
end program TAECH
