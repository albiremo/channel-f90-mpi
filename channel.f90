!============================================!
!                                            !
!     Direct Numerical Simulation (DNS)      !
!       of a turbulent channel flow          !
!                                            !
!============================================!
!
! This program has been written following the
! KISS (Keep it Simple and Stupid) philosophy
!
! Author: Dr.-Ing. Davide Gatti
!

! Measure per timestep execution time
!#define chron

PROGRAM channel

  USE dnsdata
#ifdef crhon
  REAL timei,timee
#endif
  ! Stats
  INTEGER(C_SIZE_T) :: istep=0
  character(len=40) :: istattring,filename
  LOGICAL :: io
  REAL(C_DOUBLE), allocatable :: localstats(:,:), globalstats(:,:)
  COMPLEX(C_DOUBLE_COMPLEX), allocatable :: localpsd(:,:,:), globalpsd(:,:,:)
  REAL(C_DOUBLE) :: c

  ! Init MPI
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  CALL read_dnsin()
  CALL init_MPI(iproc,nproc,nx+1,nz,ny,nxd+1,nzd,nx0,nxN,nxB,nz0,nzN,nzB,block)
  CALL init_memory()

  ! Init various subroutines
  CALL init_fft(VVdz,VVdx,rVVdx,nxd,nxB,nzd,nzB)
  CALL setup_derivatives()
  CALL setup_boundary_conditions()
  CALL read_restart_file()
  IF (.NOT. time_from_restart) CALL read_dnsin()

  ! Allocate memory for stats
  ALLOCATE(localstats(-1:ny+1,1:5),globalstats(-1:ny+1,1:5))
  globalstats = 0
  ALLOCATE(localpsd(-1:ny+1,-nz:nz,2:5),globalpsd(-1:ny+1,-nz:nz,2:5))

  ! Field number (for output)
  ifield=FLOOR(time/dt_field)

IF (has_terminal) THEN
  ! Output DNS.in
  WRITE(*,*) " "
  WRITE(*,*) "!====================================================!"
  WRITE(*,*) "!                     D   N   S                      !"
  WRITE(*,*) "!====================================================!"
  WRITE(*,*) " "
  WRITE(*,"(A,I5,A,I5,A,I5)") "   nx =",nx,"   ny =",ny,"   nz =",nz
  WRITE(*,"(A,I5,A,I5)") "   nxd =",nxd,"  nzd =",nzd
  WRITE(*,"(A,F6.4,A,F6.4,A,F8.6)") "   alfa0 =",alfa0,"       beta0 =",beta0,"   ni =",ni
  WRITE(*,"(A,F6.4,A,F6.4)") "   meanpx =",meanpx,"      meanpz =",meanpz
  WRITE(*,"(A,F6.4,A,F6.4)") "   meanflowx =",meanflowx, "   meanflowz =", meanflowz
  WRITE(*,"(A,I6,A,L1,A,I6)"   ) "   nsteps =",nstep, "   time_from_restart =", time_from_restart, "nstats =", nstats
  WRITE(*,*) " "
END IF

  ! Compute CFL
  DO iy=1,ny-1
   CALL convolutions(iy,1,.TRUE.)
  END DO
  ! Time loop
  CALL outstats()
  timeloop: DO WHILE ((time<t_max-deltat/2.0) .AND. (istep<nstep))
#ifdef chron
    CALL CPU_TIME(timei)
#endif
    ! Increment number of steps
    istep=istep+1
    ! Solve
    time=time+2.0/RK1_rai(1)*deltat
    CALL buildrhs(RK1_rai,.FALSE. ); CALL linsolve(RK1_rai(1)/deltat)
    time=time+2.0/RK2_rai(1)*deltat
    CALL buildrhs(RK2_rai,.FALSE.); CALL linsolve(RK2_rai(1)/deltat)
    time=time+2.0/RK3_rai(1)*deltat
    CALL buildrhs(RK3_rai,.TRUE.); CALL linsolve(RK3_rai(1)/deltat)
    CALL outstats()
#ifdef chron
    CALL CPU_TIME(timee)
    IF (has_terminal) WRITE(*,*) timee-timei
#endif
    ! Compute statistics
    localstats = 0
    localpsd = 0;
    IF (has_terminal) localstats(:,1)=dreal(V(:,0,0,1))
    DO ix=nx0,nxN
      c = MERGE(1.0d0,2.0d0,ix==0)
      ! uu
      localpsd(:,:,2) = localpsd(:,:,2) +c*(V(:,:,ix,1)*CONJG(V(:,:,ix,1)))
      ! vv
      localpsd(:,:,3) = localpsd(:,:,3) +c*(V(:,:,ix,2)*CONJG(V(:,:,ix,2)))
      ! ww
      localpsd(:,:,4) = localpsd(:,:,4) +c*(V(:,:,ix,3)*CONJG(V(:,:,ix,3)))
      ! uv
	    localpsd(:,:,5) = localpsd(:,:,5) +c*(V(:,:,ix,1)*CONJG(V(:,:,ix,2)))
    END DO
    IF (has_terminal) THEN
      CALL MPI_Reduce(MPI_IN_PLACE,localpsd,(ny+3)*(2*nz+1)*4,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD)
    ELSE
      CALL MPI_Reduce(localpsd,localpsd,(ny+3)*(2*nz+1)*4,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD)
    END IF
    localstats(:,5) = SUM(localpsd(:,:,5),2); localstats(:,2) = SUM(localpsd(:,:,2),2)
    localstats(:,3) = SUM(localpsd(:,:,3),2); localstats(:,4) = SUM(localpsd(:,:,4),2)
    IF (has_terminal) THEN
      globalstats = globalstats+localstats; globalpsd=globalpsd+localpsd
    END IF
    CALL MPI_Barrier(MPI_COMM_WORLD)
! Compute statistics
    IF (MOD(istep,nstats)==0) THEN
      WRITE(istattring,*) istep
      IF (has_terminal) WRITE(*,*) "saving stats"
      filename="stats."//TRIM(ADJUSTL(istattring))//".dat"
      IF (has_terminal) OPEN(UNIT=108,FILE=(filename), STATUS="new",ACCESS='stream',ACTION='readwrite')
      IF (has_terminal) THEN
        WRITE(108,POS=1) istep,globalstats/nstats,globalpsd/nstats
        CLOSE(108)
      END IF
      IF (has_terminal) WRITE(*,*) "chiudo"
  !CALL MPI_Barrier(MPI_COMM_WORLD)
    globalstats = 0 !azzero le globalstats cosi da iniziare un nuovo campione
    END IF


  END DO timeloop

  IF (has_terminal) CLOSE(102)
  ! Realease memory
  CALL free_fft()
  CALL free_memory()
  CALL MPI_Finalize()


END PROGRAM channel
