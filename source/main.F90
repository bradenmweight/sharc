!
!...This code can do nonadiabatic dyanmics simulation of nulear and eletron on
!...equal footing using quasi-diabatic (QD) propagation method based on
!...SHARC. So we borrowed many subroutines from SHARC.
!...The algorithm and the details of the method are given in
!...        J. Sebastian Sandoval C. et al. JCP 2018, 149, 044115
!...We need Hamiltonian, grdients and overlap Sik=<Phi_i(t0)|Phi_k(t1)>
!...which can be obatined using SHARC.
!...Wanghuai Zhou,11/21/2018
!
program main
use definitions
use electronic
use input
use matrix
use misc
use nuclear
use qm
use restart
use output
implicit none

!> \param traj Contains all data which would be private to each trajectory in an
!ensemble
type(trajectory_type) :: traj
!> \param ctrl Contains all data which would be shared in an ensemble
type(ctrl_type) :: ctrl
!> \param i_step Loop variable for the dynamics loop
integer :: i_step
!> \param time Define the integer function time()
integer :: time
integer::istate,jstate

traj%time_start=time()
traj%time_last=traj%time_start

!read input file
call read_input(traj,ctrl)
call allocate_lapack(ctrl%nstates)

!initialize mapping variables
call init_map(ctrl%nstates,traj%state_mch,ctrl%nsteps,ctrl%windowfunc,ctrl%gammfunc)

!call get_Dmatrix(ctrl)
!write(u_den,111)0,((Dmatrix(istate,jstate),istate=1,ctrl%nstates),jstate=1,ctrl%nstates)
!call write_Dmatrix(traj,ctrl)
 
!Do initial QM calculation to get 
!Hamiltonian matrix elements: H_MCH_ss
!Gradients: grad_MCH_sad
!Nonadiabatic coupling vectors: NACdr_ssad
!Overlap matrix: overlaps_ss
!Construct Gmatrix: Gmatrix_ssqd--diagonal gradients + non-adiabatic coupling
#ifdef  Model
   write(*,*)"#We do model calculation. In model system"
   write(*,*)"#We only consider 2 or 3 states and two atoms in x direction"
   write(*,*)"#Only one atom can move and the other one is fixed at 0."
   call do_model_calc(traj,ctrl)
   !stop
#else
  call do_first_qm(traj,ctrl)  !quantum chemistry calculation
#endif

!Calculate the initial Force
call calc_force(traj,ctrl)
!write(*,*) "I performed calc_force for the first time."

!Save old variables
call save_old(traj)

call Calculate_etot(traj,ctrl) !traj%Etot=traj%Ekin+traj%Epot
call Calculate_etot_qd(traj,ctrl) !traj%Etot_qd=traj%Ekin+traj%Epot_qd
call set_time(traj)            !
call write_dat(u_dat, traj, ctrl) 
call write_list_line(u_lis,traj,ctrl)
call write_geom(u_geo, traj, ctrl)

!evolution loop
do i_step=traj%step+1,ctrl%nsteps
!do i_step=traj%step+1,1
   traj%step=i_step
   write(*,*)"Zhou,i_step=",i_step
   call write_logtimestep(u_log,i_step)

   !...Verlet step for nuclear positions
   call VelocityVerlet_xstep(traj,ctrl) 

   !...Do QM calculation to get
   !...H_MCH_ss and grad_MCH_sad
   !...In sharc, grad_MCH_sad is transformed into Gmatrix_ss and then Gmatrix_ssad
   !call do_qm_calculations(traj,ctrl)
   !Construct Gmatrix: diagonal gradients + non-adiabatic coupling
#ifdef Model
   call do_model_calc(traj,ctrl)
#else
   call do_qm_calc_qd(traj,ctrl)
#endif

   !...begin substep for mapping variable
   !...Calcualate the HQD(t) throught interpelation
   !...Propagate the forward and backward mapping variables by solving Hamilton's
   !...equation with Hamiltonian hmF and hmB
   call propagate_mapping(traj,ctrl)

   !...Transforming mapping variables into t2 QD basis
   call transform_mapping(traj,ctrl)

   !...Calculate Force of new position with new mapping variables
   call calc_force(traj,ctrl)

   !print out force for debug
   !write(33,*)i_step,traj%grad_ad(1,1),traj%geom_ad(1,1)

   !...Verlet for nuclear velocity
   call VelocityVerlet_vstep(traj,ctrl)

   call Calculate_etot(traj,ctrl)
   call Calculate_etot_qd(traj,ctrl)

   !...Save variables
   call save_old(traj)
   call set_time(traj)
   call write_list_line(u_lis,traj,ctrl)
   call write_dat(u_dat, traj, ctrl)
   call write_geom(u_geo, traj, ctrl)
   call allflushqd()
   
   !...population window
   call window(traj,ctrl)
   call write_Dmatrix(traj,ctrl)

enddo !i_step loop


call write_final(traj)

end program main
