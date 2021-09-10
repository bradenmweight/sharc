!
!...QM calculations in QD method
!
!...first time QM calculations
subroutine do_first_qm(traj,ctrl)
use definitions
use electronic
use matrix
use output
use qm
implicit none
type(trajectory_type) :: traj
type(ctrl_type) :: ctrl
integer :: i, iatom, idir

if (printlevel>1) then
   call write_logtimestep(u_log,traj%step)
endif

!...initial QM calculation
call do_qm_calc_qd(traj,ctrl)

!...now finalize the initial state, if this was not done in the read_input
!...routine
if (printlevel>1) then
  write(u_log,*)'============================================================='
  write(u_log,*) '                    First QM calculation'
  write(u_log,*)'============================================================='
endif

! check whether the initial state is active
! (active states are defined in MCH basis, initial state might be given in
! diagonal basis
if (ctrl%actstates_s(traj%state_MCH).eqv..false.) then
  write(0,*) 'Initial state is not active!'
  stop 1
endif

if (printlevel>1) then
  write(u_log,'(a,1x,i3,1x,a)') 'Initial state is ',traj%state_mch,'in the MCH basis. '
  write(u_log,*) 'Coefficients (MCH):'
  write(u_log,'(a3,1x,A12,1X,A12)') '#','Real(c)','Imag(c)'
  do i=1,ctrl%nstates
    write(u_log,'(i3,1x,F12.9,1X,F12.9)') i,traj%coeff_MCH_s(i)
  enddo
endif
!if (abs(traj%coeff_diag_s(traj%state_diag))<1.d-9) then
!  write(0,*) 'Initial state has zero population!'
!  write(0,*) 'In subroutine do_first_qm!'
!  stop 1
!endif

!save some variables

return
end subroutine do_first_qm
!
!...all the quantum chemistry calculation are done here.
!...We get Hamiltonian, gradients, overlap, non-adiabatic coupling matrix elements
!...We constrcut Gmatrix=gradients+non-adiabatic coupling
subroutine do_qm_calc_qd(traj,ctrl)
use definitions
use electronic
use matrix
use qm_out
use restart
use qm
implicit none
type(trajectory_type) :: traj
type(ctrl_type) :: ctrl
integer :: i,j,stat
integer::iatom,idir,istate,jstate
real*8::Gmatrix_ss(ctrl%nstates,ctrl%nstates)

if (printlevel>3) then
  write(u_log,*)'============================================================='
  write(u_log,*) '                       QM calculation'
  write(u_log,*)'============================================================='
  write(u_log,*) 'QMin file="QM/QM.in"'
endif

!open QM.in and write geometry + some keywords (task keywords are written by write_tasks_*)
open(u_qm_qmin,file='QM/QM.in',status='replace',action='write')
call write_infos(traj,ctrl)

!...if necessary, select the quantities for calculation
if (ctrl%calc_grad==1) call select_grad(traj,ctrl)
if (ctrl%calc_nacdr==1) call select_nacdr(traj,ctrl)
if (ctrl%calc_dipolegrad==1) call select_dipolegrad(traj,ctrl)

!...write tasks for first QM call
call write_tasks_first(traj,ctrl)
close(u_qm_qmin)

!...run QM interface
if (printlevel>3) write(u_log,*) 'Running file="QM/runQM.sh"'
call call_runqm(traj)

!...open QM.out
if (printlevel>3) write(u_log,*) 'QMout file="QM/QM.out"'
call open_qmout(u_qm_qmout, 'QM/QM.out')

!...get Hamiltonian
call get_hamiltonian(ctrl%nstates, traj%H_MCH_ss)
!...apply reference energy shift
do i=1,ctrl%nstates
  traj%H_MCH_ss(i,i)=traj%H_MCH_ss(i,i)-ctrl%ezero
enddo
if (printlevel>3) write(u_log,'(A31,A2)') 'Hamiltonian:                   ','OK'

!...take the trace of H out
traj%H0=0.d0
do i=1,ctrl%nstates
   traj%H0=traj%H0+traj%H_MCH_ss(i,i)
enddo
traj%H0=traj%H0/ctrl%nstates
do i=1,ctrl%nstates
   traj%H_MCH_ss(i,i)=traj%H_MCH_ss(i,i)-traj%H0*dble(ctrl%trace)
enddo

!...in QD, we just set U=[1]
traj%U_ss=dcmplx(0.d0,0.d0)
do i=1,ctrl%nstates
   traj%U_ss(i,i)=dcmplx(1.d0,0.d0)
enddo

!print Hamiltonian for check
!write(*,*)"Zhou,Hamiltonian="
!do i=1,ctrl%nstates
!   write(*,*)traj%H_MCH_ss(i,:)
!enddo

!...get all available Properties
call get_properties_new(ctrl, traj)

!...get gradients
call get_gradients(ctrl%nstates, ctrl%natom, traj%grad_MCH_sad)
if (printlevel>3) write(u_log,'(A31,A2)') 'Gradients:                     ','OK'
!...take the trace of G out
traj%G0=0.d0
do istate=1,ctrl%nstates
   do iatom=1,ctrl%natom
      do idir=1,3
         traj%G0(iatom,idir)=traj%G0(iatom,idir)+traj%grad_MCH_sad(istate,iatom,idir)
      enddo
   enddo
enddo
traj%G0=traj%G0/ctrl%nstates
do istate=1,ctrl%nstates
   do iatom=1,ctrl%natom
     traj%grad_MCH_sad(istate,iatom,:)=traj%grad_MCH_sad(istate,iatom,:)-&
     traj%G0(iatom,:)*dble(ctrl%trace)
   enddo
enddo

!...if this is not the initial QM calculation. Get overlap matrix.
if (traj%step>=1) then
   call get_overlap(ctrl%nstates, traj%overlaps_ss)
   if (printlevel>3) write(u_log,'(A31,A2)') 'Overlap matrix:','OK'
endif

!...get wavefunction phases
call get_phases(ctrl%nstates,traj%phases_s,stat)
if (stat==0) then
   traj%phases_found=.true.
   if(printlevel>3) write(u_log,'(A31,A2)') 'Phases:                      ','OK'
   !write(35,'(6(e23.13))')(traj%phases_s(i),i=1,ctrl%nstates)
else
   traj%phases_found=.false.
   if (printlevel>3) write(u_log,'(A31,A9)') 'Phases:                     ','NOT FOUND'
endif

!...get nonadiabatic coupling
call get_nonadiabatic_ddr(ctrl%nstates, ctrl%natom, traj%NACdr_ssad)
if (printlevel>3) write(u_log,'(A31,A2)') 'Non-adiabatic couplings (DDR):','OK'
!...write NACdr for debugging
write(81,'(3e23.13)')(real(traj%NACdr_ssad(1,2,1,i)),i=1,3)



call close_qmout
if (printlevel>3) write(u_log,*) ''

if (printlevel>4) call print_qm(u_log,traj,ctrl)

!...change the phase when we need
!...Even the phase is not changed, this step is necessary.
!...The overlap matrix is orthogonalized during the process.
if (traj%step>=1) then
  call Adjust_phases(traj,ctrl) 
endif

!...contstruc Gmatrix after we get gradients and Non-dadiabatic couplings
call construct_Gmatrix(traj,ctrl)

return
end subroutine do_qm_calc_qd
!
!...construct Gmatrix which contains gradients (force) and NACME
!...nonadiabaic coupling matrix elements
!
subroutine construct_Gmatrix(traj,ctrl)
use definitions
use matrix
implicit none
type(trajectory_type) :: traj
type(ctrl_type) :: ctrl
integer::iatom,idir,istate,jstate
real*8::Gmatrix_ss(ctrl%nstates,ctrl%nstates)

if (printlevel>3) then
  write(u_log,*)'============================================================='
  write(u_log,*) '             Constructing Gmatrix'
  write(u_log,*)'============================================================='
endif

!Gmatrix_ssqd=0.d0

do iatom=1,ctrl%natom
   do idir=1,3

      Gmatrix_ss=0.d0
      do istate=1,ctrl%nstates  !diagonal part
          Gmatrix_ss(istate,istate)=traj%grad_MCH_sad(istate,iatom,idir)
      enddo !state loop

      !off-diagonal part
      do istate=1,ctrl%nstates
         do jstate=istate+1,ctrl%nstates
            Gmatrix_ss(istate,jstate)=&
              &-(traj%H_MCH_ss(istate,istate)-traj%H_MCH_ss(jstate,jstate))&
              &*traj%NACdr_ssad(istate,jstate,iatom,idir)
              Gmatrix_ss(jstate,istate)=Gmatrix_ss(istate,jstate)
         enddo
      enddo

      if(printlevel>4) then
          write(u_log,'(A,1X,I4,1X,A,1X,I4)') 'Gmatrix calculation... iatom=',&
          iatom,'idir=',idir
          call matwrite(ctrl%nstates,Gmatrix_ss,u_log,'Gmatrix MCH','F12.9')
      endif

      !save full G matrix in traj
      traj%Gmatrix_ssqd(:,:,iatom,idir)=Gmatrix_ss

   enddo !direction loop
enddo !atom loop

!print Gmatrix
!write(*,*)"Zhou,Gmatrix="
!do istate=1,ctrl%nstates
!   do jstate=istate,ctrl%nstates
!      write(*,*)"state",istate,jstate
!      do iatom=1,ctrl%natom
!         write(*,*)traj%Gmatrix_ssqd(istate,jstate,iatom,:)
!      enddo
!   enddo
!enddo

return
end subroutine construct_Gmatrix
!
!...calculate the QD force using mapping variables
!...exactly we calculate gradients = -force
subroutine calc_force(traj,ctrl)
use definitions
use matrix
implicit none
type(trajectory_type) :: traj
type(ctrl_type) :: ctrl
integer::iatom,idir,istate,jstate
real*8::dentemp !temporary variable

traj%grad_ad=0.d0

#ifdef Model
traj%grad_ad=traj%grad_ad+traj%Gmatrix0
#else
do iatom=1,ctrl%natom   !add state-independent part
   traj%grad_ad(iatom,:)=traj%grad_ad(iatom,:)+traj%G0(iatom,:)*dble(ctrl%trace)
enddo
#endif

! Check to see if we are correcting gamma for each trajectory ~ Braden Weight
if (ctrl%gammfunc .eq. "yes") then

   do istate=1,ctrl%nstates
      do jstate=1,ctrl%nstates
         if(istate/=jstate)then
         ActionMat(istate,jstate)=0.5d0*(qF(istate)*qF(jstate)+pF(istate)*pF(jstate))
         else
         ActionMat(istate,jstate)=0.5d0*(qF(istate)*qF(jstate)+pF(istate)*pF(jstate)-2.d0*ZPE(istate))
         endif 
         do iatom=1,ctrl%natom
            do idir=1,3
               traj%grad_ad(iatom,idir)=traj%grad_ad(iatom,idir)+&
               traj%Gmatrix_ssqd(istate,jstate,iatom,idir)*ActionMat(istate,jstate)
            enddo !idir
         enddo !iatom
      enddo !jstate
   enddo !istate

else

   do istate=1,ctrl%nstates
      do jstate=1,ctrl%nstates
         if(istate/=jstate)then
         ActionMat(istate,jstate)=0.5d0*(qF(istate)*qF(jstate)+pF(istate)*pF(jstate))
         else
            ActionMat(istate,jstate)=0.5d0*(qF(istate)*qF(jstate)+pF(istate)*pF(jstate)-2.d0*gamm)
         endif 
         do iatom=1,ctrl%natom
            do idir=1,3
               traj%grad_ad(iatom,idir)=traj%grad_ad(iatom,idir)+&
               traj%Gmatrix_ssqd(istate,jstate,iatom,idir)*ActionMat(istate,jstate)
            enddo !idir
         enddo !iatom
      enddo !jstate
   enddo !istate

end if ! Check if gamma function is adjusted or not

!debug use
!traj%grad_ad(:,:)=real(traj%Gmatrix_ssad(traj%state_diag,traj%state_diag,:,:))
!Please comment the above lines when you do real simulations

!print force for check
!write(*,*)"Zhou,grad_ad="
!do iatom=1,ctrl%natom
!   write(*,*)traj%grad_ad(iatom,:)
!enddo

return
end subroutine calc_force 
!
!...save all *_old matrices when advancing to a new timestep
subroutine save_old(traj)
use definitions
implicit none
type(trajectory_type) :: traj

if (printlevel>3) then
  write(u_log,*)'============================================================='
  write(u_log,*) '                   Advancing to next step'
  write(u_log,*)'============================================================='
endif

!...initialize old variables from current ones
traj%H_MCH_old_ss=traj%H_MCH_ss  !save old Hamiltonian
traj%NACdr_old_ssad=traj%NACdr_ssad !Maybe not use

!...save transformation matrix U
!...In model system we calculate U explicitly.
!...In QM calculation we just set U = 1
traj%U_old_ss=traj%U_ss             

traj%DM_old_ssd=traj%DM_ssd
traj%phases_old_s=traj%phases_s

return
end subroutine save_old
!====================================================================
!
!...Tully model calculation
!
subroutine do_model_calc(traj,ctrl)
use definitions
use electronic
use matrix
use qm_out
use restart
use qm
implicit none
type(trajectory_type) :: traj
type(ctrl_type) :: ctrl
!integer :: i,j
integer::iatom,idir,istate,jstate
real*8::Gmatrix_ss(ctrl%nstates,ctrl%nstates),si

!...get model Hamiltonian and gradients in diabatic basis
select case (trim(ctrl%model))
  case ('tully2')
    if(ctrl%nstates/=2)then
      write(*,*)"In Tully2 model, we only consider two states and two atoms"
      stop
    endif
    call get_model_hami(ctrl%nstates,ctrl%natom,traj%geom_ad,&
                        traj%H_MCH_ss,traj%Gmatrix_ssqd)

  case default
    write(*,*)"Model we donot know!"
    stop 1
end select 

if (printlevel>3) write(u_log,'(A31,A2)') 'Hamiltonian:                   ','OK'
!write(*,*)"Zhou,Tully Hamiltonian in diabatic="
!do istate=1,ctrl%nstates
!   write(*,*)traj%H_MCH_ss(istate,:)
!enddo
!write(*,*)"Zhou,Gmatrix_ssqd="
!do istate=1,ctrl%nstates
!   write(*,*)traj%Gmatrix_ssqd(istate,istate,1,1)
!enddo

!...calculate state-independent gradient
traj%Gmatrix0=0.5d0*(traj%Gmatrix_ssqd(1,1,1,1)+traj%Gmatrix_ssqd(2,2,1,1))
!...seperate the state-independent gradient
do istate=1,2
   traj%Gmatrix_ssqd(istate,istate,1,1)=traj%Gmatrix_ssqd(istate,istate,1,1)&
                                        -traj%Gmatrix0
enddo
!write(*,*)"Zhou,Gmatrix_ssqd2="
!do istate=1,ctrl%nstates
!   write(*,*)traj%Gmatrix_ssqd(istate,istate,1,1)
!enddo

!...diagonalize the Hamiltoian
!write(84,'(9e23.13)')((real(traj%H_MCH_ss(istate,jstate)),&
!                       jstate=1,ctrl%nstates),istate=1,ctrl%nstates)
traj%H_diag_ss=traj%H_MCH_ss
call diagonalize(ctrl%nstates,traj%H_diag_ss,traj%U_ss)
traj%H_MCH_ss=traj%H_diag_ss    !now H_MCH_ss is diagnoal

!...adjust the phase of eigenvector
if (traj%step>=1) then
   do istate=1,ctrl%nstates
      call phase_match(ctrl%nstates,traj%U_old_ss,traj%U_ss(:,istate),jstate)
      si=dot_product(traj%U_old_ss(:,jstate),traj%U_ss(:,istate))
      !write(66,*)si,istate,jstate
      si=si/abs(si)
      traj%U_ss(:,istate)=traj%U_ss(:,istate)*si
   enddo
endif

!...the initial condition is in diabatic basis
!...So we need to transform the mapping variables to QD (adiabatic) basis
if(traj%step==0) call diabatic2QD(traj,ctrl)

!write(*,*)"Zhou,Tully U matrix="
!write(55,'(9e23.13)')((real(traj%U_ss(istate,jstate)),&
!                       jstate=1,ctrl%nstates),istate=1,ctrl%nstates)

!get the adiabatic force
!write(87,'(9e23.13)')((traj%Gmatrix_ssqd(istate,jstate,1,1),&
!                       jstate=1,ctrl%nstates),istate=1,ctrl%nstates)
traj%Gmatrix_ssad=traj%Gmatrix_ssqd !copy real to complex matrix
!call transform(ctrl%nstates,traj%Gmatrix_ssad(:,:,1,1),traj%U_ss,'uaut') !Maybe Wrong
call transform(ctrl%nstates,traj%Gmatrix_ssad(:,:,1,1),traj%U_ss,'utau')
traj%Gmatrix_ssqd=real(traj%Gmatrix_ssad)
!write(88,'(9e23.13)')((traj%Gmatrix_ssqd(istate,jstate,1,1),&
!                       jstate=1,ctrl%nstates),istate=1,ctrl%nstates)
!write(*,*)"Zhou,Tully Force in adiabatic="
!do istate=1,ctrl%nstates
!   write(*,*)traj%Gmatrix_ssad(istate,:,1,1)
!enddo
!write(*,*)"Zhou,Tully Force in adiabatic2="
!do istate=1,ctrl%nstates
!   write(*,*)traj%Gmatrix_ssqd(istate,:,1,1)
!enddo

!...get overlap
if (traj%step>=1) then
   call get_overlap_model(ctrl%nstates,traj%U_old_ss,traj%U_ss,traj%overlaps_ss)
   !write(81,'(9f16.8)')((real(traj%overlaps_ss(istate,jstate)),&
   !                      jstate=1,ctrl%nstates),istate=1,ctrl%nstates)
   if (printlevel>3) write(u_log,'(A31,A2)') 'Overlap matrix:','OK'
endif

return
end subroutine do_model_calc
!
!...find the match phase
subroutine phase_match(ns,U_old,U_ss,j)
implicit none
integer,intent(in)::ns
complex*16,intent(in)::U_old(ns,ns),U_ss(ns)
integer,intent(out)::j
integer::i
real*8::theta,tmp

tmp=-1.d0
!...we try to find the largest absolute value of dot_product 
do i=1,ns
   theta=dot_product(U_old(:,i),U_ss)
   if(abs(theta)>tmp)then
     j=i
     tmp=abs(theta)
   endif
enddo

end subroutine phase_match
!
!...get model Hamiltonian of Tully model II
subroutine get_model_hami(ns,natom,pos,H,Gmatrix)
implicit none
integer,intent(in)::ns !# of states and must be equal to 2.
integer,intent(in)::natom !# of atoms and must be equal to 2.
real*8,intent(in)::pos(natom,3)            !positions
complex*16,intent(out)::H(ns,ns)           !Hamiltonian
complex*16::H0 !state independent Hamiltonian
real*8,intent(out)::Gmatrix(ns,ns,natom,3) !gradients
real*8,parameter::A=0.1d0
real*8,parameter::B=0.28d0
real*8,parameter::E0=0.05d0
real*8,parameter::C=0.015d0
real*8,parameter::D=0.06d0
real*8::He,x

!begin Hamiltonian
H(1,1)=dcmplx(0.d0,0.d0)
x=pos(1,1) !x= x coordinate of first atom
He=-A*exp(-B*x*x)+E0
H(2,2)=dcmplx(He,0.d0)

He=C*exp(-D*x*x)
H(1,2)=dcmplx(He,0.d0)
H(2,1)=H(1,2)
H0=0.5d0*(H(1,1)+H(2,2))
H(1,1)=H(1,1)-H0
H(2,2)=H(2,2)-H0
!end Hamiltonian

!begin gradients. We can directedly get Gmatrix
!traj%Gmatrix_ssqd(nstates,nstates,iatom,idir)
!We only calculate the force on the first atom.
Gmatrix=0.d0

He=-2.d0*C*D*x*exp(-D*x*x)
Gmatrix(1,2,1,1)=dcmplx(He,0.d0)
Gmatrix(2,1,1,1)=dcmplx(He,0.d0)

He=2.d0*A*B*x*exp(-B*x*x)
Gmatrix(2,2,1,1)=dcmplx(He,0.d0)

return
end subroutine get_model_hami
!
!...get overlap matrix
subroutine get_overlap_model(ns,U_old,U,overlap)
use matrix
implicit none
integer,intent(in)::ns
complex*16,intent(in)::U_old(ns,ns),U(ns,ns)
complex*16,intent(out)::overlap(ns,ns)

call matmultiply(ns,U_old,U,overlap,'tn')

return
end subroutine get_overlap_model
!
!...get model Hamiltonian of Morse model
subroutine get_model_morse(ns,natom,pos,H,Gmatrix)
implicit none
integer,intent(in)::ns !# of states and must be equal to 2.
integer,intent(in)::natom !# of atoms and must be equal to 2.
real*8,intent(in)::pos(natom,3)            !positions
complex*16,intent(out)::H(ns,ns)           !Hamiltonian
complex*16::H0 !state independent Hamiltonian
real*8,intent(out)::Gmatrix(ns,ns,natom,3) !gradients
real*8::D1(3),a(3),b(3),E(3),bigA(3,3),c(3,3),d(3,3) !model parameters
real*8::He,x,temp
integer i,j

D1(1)=0.003d0; D1(2)=0.004d0; D1(3)=0.003d0
 a(1)=0.65d0;   a(2)=0.6d0;    a(3)=0.65d0
 b(1)=5.d0;     b(2)=4.d0;     b(3)=6.d0
 E(1)=0.d0;     E(2)=0.01d0;   E(3)=0.006d0
bigA(1,2)=0.002d0; bigA(1,3)=0.d0; bigA(2,3)=0.002d0 
   c(1,2)=16.d0;      c(1,3)=0.d0;    c(2,3)=16.d0
   d(1,2)=3.4d0;      d(1,3)=0.d0;    d(2,3)=4.8d0

x=pos(1,1) !x= x coordinate of first atom
!write(*,*)"x=",x

!begin Hamiltonian; diagonal part
He=D1(1)*(1.d0-exp(-a(1)*(x-b(1))))**2+E(1)
H(1,1)=dcmplx(He,0.d0)
He=D1(2)*(1.d0-exp(-a(2)*(x-b(2))))**2+E(2)
H(2,2)=dcmplx(He,0.d0)
He=D1(3)*(1.d0-exp(-a(3)*(x-b(3))))**2+E(3)
H(3,3)=dcmplx(He,0.d0)

!off-diagonal part
He=bigA(1,2)*exp(-c(1,2)*(x-d(1,2))**2)
H(1,2)=dcmplx(He,0.d0)
H(2,1)=dcmplx(He,0.d0)

He=bigA(1,3)*exp(-c(1,3)*(x-d(1,3))**2)
H(1,3)=dcmplx(He,0.d0)
H(3,1)=dcmplx(He,0.d0)

He=bigA(2,3)*exp(-c(2,3)*(x-d(2,3))**2)
H(2,3)=dcmplx(He,0.d0)
H(3,2)=dcmplx(He,0.d0)
!end Hamiltonian
!write(104,'(9e23.13)')((real(H(i,j)),j=1,ns),i=1,ns)

!begin gradients. We can directedly get Gmatrix
!traj%Gmatrix_ssqd(nstates,nstates,iatom,idir)
!We only calculate the force on the first atom.
Gmatrix=0.d0

!diagonal part
do i=1,ns
   temp=exp(-a(i)*(x-b(i)))
   He=2.d0*D1(i)*a(i)*temp*(1.d0-temp)
   Gmatrix(i,i,1,1)=dcmplx(He,0.d0)
enddo

!off-diagonal part
He=-2.d0*bigA(1,2)*c(1,2)*(x-d(1,2))*exp(-c(1,2)*(x-d(1,2))**2)
Gmatrix(1,2,1,1)=dcmplx(He,0.d0)
Gmatrix(2,1,1,1)=dcmplx(He,0.d0)
He=-2.d0*bigA(1,3)*c(1,3)*(x-d(1,3))*exp(-c(1,3)*(x-d(1,3))**2)
Gmatrix(1,3,1,1)=dcmplx(He,0.d0)
Gmatrix(3,1,1,1)=dcmplx(He,0.d0)
He=-2.d0*bigA(2,3)*c(2,3)*(x-d(2,3))*exp(-c(2,3)*(x-d(2,3))**2)
Gmatrix(2,3,1,1)=dcmplx(He,0.d0)
Gmatrix(3,2,1,1)=dcmplx(He,0.d0)
!write(105,'(9e23.13)')((real(Gmatrix(i,j,1,1)),j=1,ns),i=1,ns)

return
end subroutine get_model_morse
