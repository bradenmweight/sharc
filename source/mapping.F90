!...
!...This file contains all the subroutines relatated mapping variables
!...
!...initialize mapping variables
!...nstates=# of states
!...j=current state=We choose
subroutine init_map(nstates,j,nsteps,windowfunc,gammfunc)
use definitions
implicit none
integer,intent(in)::nstates,j,nsteps
integer::i
integer::tri
integer::dum
real*8::ran
character(len=20) windowfunc
character(len=20) gammfunc


if(j==0)stop"initial state in MCH is wrong!"
!write(*,*)"Zhou,nstates,j=",nstates,j

!QD basis
allocate(qF(nstates))  !Forward mapping variables
allocate(pF(nstates))  !Forward
allocate(Dmatrix(nstates,nstates))  !density matrix
qF=0.d0 ; pF=0.d0 ; Dmatrix=0.d0 

!diabatic basis
allocate(qFd(nstates))  !Forward mapping variables
allocate(pFd(nstates))  !Forward
allocate(ZPE(nstates))
allocate(ek(nstates))
qFd=0.d0 ; pFd=0.d0 

! Braden Weight ~ 11/15/2020
! If left at 0, window functions should become delta functions.
gamm = 0
allocate(acti(nstates),angle(nstates))
acti=0.d0 ; angle=0.d0

if (gammfunc .eq. "no") then ! If not adjusted, use original action and angle defs.
   if (windowfunc .eq. "rectangle") then
      !Square
      gamm=(sqrt(3.d0)-1.d0)*0.5d0
      do i=1,nstates
         call random_number(ran)
         acti(i)=gamm*(2.d0*ran-1.d0) ! Random number between -1 and 1
         if(i==j)acti(i)=acti(i)+1.d0  ! Active state adds one to shift square up
         call random_number(ran)
         angle(i)=ran*2.d0*pi
         !write(*,*)"acti=",acti(i)
      enddo
      do i=1,nstates
         write(*,*)"TEST:",i,acti(i)
      enddo

   elseif (windowfunc .eq. "pyramidal") then
      gamm = 1/3.d0
      do while (.TRUE.)
         call random_number(ran)
         ek(j) = ran
         call random_number(ran)
         if (1 - ek(j) >= ran) then
            exit
         endif
      enddo

      ! Unoccupied DOF
      do i=1, nstates
         call random_number(ran)
         angle(i)=ran*2.d0*pi ! Initialize random angle for all states
         if (i .ne. j) then
            call random_number(ran)
            ran = ran * (1 - ek(j))
            ek(i) = ran
         endif ! i == j?
      enddo ! nstates i

      ! Shift occupied state up
      ek(j) = ek(j) + 1

      ! Convert from ek to nk by adding window gamma
      do i=1, nstates
         acti(i) = ek(i) - gamm
      enddo

   else
      write(*,*) "'windowfunc' type not recognized or not given. Choose either 'rectangle' or 'pyramidal'."
      stop 1
   endif ! Window Type

elseif (gammfunc .eq. "yes") then
   if (windowfunc .eq. "pyramidal") then ! Code translated from Cotton/Miller Gamma-Adjusted Paper
      ! Weighted sampling of DOF for initial states
      do while (.TRUE.)
         call random_number(ran)
         ek(j) = ran
         call random_number(ran)
         if (1 - ek(j) >= ran) then
            exit
         endif
      enddo

      ! Unoccupied DOF
      do i=1, nstates
         call random_number(ran)
         angle(i)=ran*2.d0*pi ! Initialize random angle for all states
         if (i .ne. j) then
            call random_number(ran)
            ran = ran * (1 - ek(j))
            ek(i) = ran
         endif ! i == j?
      enddo ! nstates i

      ! Shift occupied state up
      ek(j) = ek(j) + 1

   elseif (windowfunc .eq. "rectangle") then
      do i=1,nstates
         call random_number(ran)
         angle(i) = ran*2.d0*pi
         call random_number(ran)
         ek(i) = 2*0.366*ran
      enddo
      ! Shift occupied state up
      ek(j) = ek(j) + 1

   else ! Adjusted pyramidal only works right now.
      write(*,*) "Adjusted pyramidal and rectangle are the only working codes right now."
      stop 1
   endif ! Pyramidal or rectangle?
else
   write(*,*) "'adjudtedGamma' choice not recognized or not given. Choose either 'yes' or 'no'."
   stop 1
endif ! Adjusted Gamma?

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! COMPELTED SAMPLING PROCEDURES !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!convert the (action,angle) into mapping viables

if (gammfunc .eq. "no") then
   do i=1,nstates
      qF(i)=sqrt(2.d0*(acti(i)+gamm))*cos(angle(i))
      pF(i)=-sqrt(2.d0*(acti(i)+gamm))*sin(angle(i))
      write(*,*) "qF,pF:",qF(i), pF(i)
   enddo

   write(*,*) "Initial Conditions: i,gamma,ek,nk"
   do i=1,nstates
      write(*,*) i, gamm, 0.5*(qF(i)**2 + pF(i)**2), 0.5*(qF(i)**2 + pF(i)**2) - gamm
      write(*,*) i, gamm, ek(i), acti(i)
   end do
   
elseif (gammfunc .eq. "yes") then
   do i=1,nstates
      qF(i) = sqrt(2.d0*ek(i)) * cos(angle(i))
      pF(i) = -sqrt(2.d0*ek(i)) * sin(angle(i)) ! Eq. 9 a,b Adj. Gamma
   enddo

   do i=1, nstates
      if (i .ne. j) then
         ZPE(i) = ek(i)
      else
         ZPE(i) = ek(i) - 1
      endif
   enddo

   write(*,*) "Initial Conditions: i,ZPE,ek,nk"
   do i=1,nstates
      write(*,*) i, ZPE(i), ek(i), 0.5*(qF(i)**2 + pF(i)**2) - ZPE(i)
   end do

endif

qFd = qF
pFd = pF



#ifndef Model
qF=qFd ; pF=pFd
#endif

allocate(histoF(nstates))
histoF=0.d0

allocate(ActionMat(nstates,nstates))

return
end subroutine init_map





!
!...propagate mapping varialbes using Verlet method
!
subroutine propagate_mapping(traj,ctrl)
use definitions
use matrix
implicit none
type(trajectory_type) :: traj
type(ctrl_type) :: ctrl
real*8,allocatable::Hnow(:,:),Hold(:,:),Hqd(:,:),Sss(:,:),Hdiff(:,:)
real*8,allocatable::vqF(:),vpF(:),accqF(:)
real*8::dt
integer::istate,jstate,i,n

!in sharc the Hamiltonian matrix is complex
!we change it to real.

!write(*,*)"Zhou,nstates=",ctrl%nstates
allocate(Hnow(ctrl%nstates,ctrl%nstates))  !current Hamiltonian in adiabatic basis t2
allocate(Hold(ctrl%nstates,ctrl%nstates))  !Old Hamiltonian in QD t1
allocate(Hqd(ctrl%nstates,ctrl%nstates))   !QD Hamiltonian
allocate(Sss(ctrl%nstates,ctrl%nstates))   !overlap matrix
allocate(Hdiff(ctrl%nstates,ctrl%nstates)) !=HQD(t2)-HQD(t1)

Hnow=real(traj%H_MCH_ss)
Hold=real(traj%H_MCH_old_ss)
Sss=real(traj%overlaps_ss)  !change complex to real

!print results for check
!write(*,*)"Zhou,Hnow="
!do istate=1,ctrl%nstates
!   write(*,*)Hnow(istate,:)
!enddo
!write(*,*)"Zhou,Hold="
!do istate=1,ctrl%nstates
!   write(*,*)Hold(istate,:)
!enddo
!write(*,*)"Zhou,overlap="
!do istate=1,ctrl%nstates
!   write(*,*)Sss(istate,:)
!enddo

!...transform the current Hamiltonian into QD basis
!...We obtain new H = S*H*S^t
call transform(ctrl%nstates,Hnow,Sss,'uaut')
!write(*,*)"Zhou,Hnow in QDt1 basis="
!do istate=1,ctrl%nstates
!   write(*,*)Hnow(istate,:)
!enddo
Hdiff=Hnow-Hold  !H(t2)-H(t1) at QD basis at t1

!write(*,*)"Zhou,ctrl%nsubsteps=",ctrl%nsubsteps
!write(*,*)"Zhou,ctrl%dtstep=",ctrl%dtstep
dt=ctrl%dtstep/ctrl%nsubsteps !in a.u. 

allocate(vqF(ctrl%nstates),vpF(ctrl%nstates)) !velocity of q and p (forward)
allocate(accqF(ctrl%nstates)) !acceleration of qF and qB 

!write(*,*)"Zhou,qF,pF=",qF,pF
do i=1,ctrl%nsubsteps

   !interpolate Hamiltonian at t in [t1,t2]
   Hnow=Hold+Hdiff*dble(i)/dble(ctrl%nsubsteps)
   
   !calculate current first derivatives
   do istate=1,ctrl%nstates
      vqF(istate)=0.d0; vpf(istate)=0.d0
      do jstate=1,ctrl%nstates
         vqF(istate)=vqF(istate)+Hnow(istate,jstate)*pF(jstate)
         vpF(istate)=vpF(istate)-Hnow(istate,jstate)*qF(jstate)
      enddo
      !write(*,*)"vq,vp=",vqF(istate),vpF(istate)
   enddo

   !calculate the current second derivatives
   do istate=1,ctrl%nstates
      accqF(istate)=0.d0 
      do jstate=1,ctrl%nstates
         accqF(istate)=accqF(istate)+Hnow(istate,jstate)*vpF(jstate)
      enddo
      !write(*,*)"accq=",accqF(istate)
   enddo

   !advance first 1/2 Verlet step
   do istate=1,ctrl%nstates
      qF(istate)=qF(istate)+vqF(istate)*dt+0.5d0*accqF(istate)*dt*dt
      pF(istate)=pF(istate)+0.5d0*vpF(istate)*dt
   enddo
   
   !calculate new first derivatives
   do istate=1,ctrl%nstates
      vpF(istate)=0.d0
      do jstate=1,ctrl%nstates
         vpF(istate)=vpF(istate)-Hnow(istate,jstate)*qF(jstate)
      enddo
   enddo
   
   !advance second 1/2 Verlet step
   do istate=1,ctrl%nstates
      pF(istate)=pF(istate)+0.5d0*vpF(istate)*dt
   enddo
   
enddo !for substeps

deallocate(Hnow,Hold,Hqd,Sss,Hdiff)
deallocate(vqF,vpF,accqF)
return
end subroutine propagate_mapping
!
!...transform mapping variables into t2 QD basis
!
subroutine transform_mapping(traj,ctrl)
use definitions
use matrix
implicit none
type(trajectory_type) :: traj
type(ctrl_type) :: ctrl
real*8,allocatable::Sss(:,:)
real*8,allocatable::qF2(:),pF2(:)
integer::istate,jstate

allocate(Sss(ctrl%nstates,ctrl%nstates))   !overlap matrix
Sss=real(traj%overlaps_ss)

allocate(qF2(ctrl%nstates),pF2(ctrl%nstates))

qF2=0.d0; pF2=0.d0

do istate=1,ctrl%nstates
   do jstate=1,ctrl%nstates
      qF2(istate)=qF2(istate)+qF(jstate)*Sss(jstate,istate)
      pF2(istate)=pF2(istate)+pF(jstate)*Sss(jstate,istate)
   enddo
enddo

qF=qF2; pF=pF2

deallocate(Sss); deallocate(qF2,pF2)
return
end subroutine transform_mapping
!
!...transform the mapping variables from diabatic to QD basis
!
subroutine diabatic2QD(traj,ctrl)
use definitions
implicit none
type(trajectory_type) :: traj
type(ctrl_type) :: ctrl
real*8,allocatable::U(:,:)
integer istate,jstate

allocate(U(ctrl%nstates,ctrl%nstates))
U=real(traj%U_ss)
!write(*,*)"U=",U
!write(*,*)"mapping in D=",qFd,pFd
!write(*,*)"mapping in D=",sum(U(:,1)*qFd(:))

do istate=1,ctrl%nstates
   qF(istate)=sum(U(:,istate)*qFd(:))
   pF(istate)=sum(U(:,istate)*pFd(:))
enddo
!write(*,*)"mapping=",qF,pF

deallocate(U)
end subroutine diabatic2QD
!
!...transform the mapping variables from QD basis to diabatic
!
subroutine QD2diabatic(traj,ctrl)
use definitions
implicit none
type(trajectory_type) :: traj
type(ctrl_type) :: ctrl
real*8,allocatable::U(:,:)
integer istate,jstate

allocate(U(ctrl%nstates,ctrl%nstates))
U=real(traj%U_ss)
!write(*,*)"U=",U
!write(*,*)"mapping in D=",qF,pF

do istate=1,ctrl%nstates
   qFd(istate)=sum(U(istate,:)*qF(:))
   pFd(istate)=sum(U(istate,:)*pF(:))
enddo
!write(*,*)"mapping=",qFd,pFd

deallocate(U)
end subroutine QD2diabatic
!
!...calculate the total energy in QD basis using mapping variables
!
subroutine Calculate_etot_qd(traj,ctrl)
  use definitions
  implicit none
  type(trajectory_type) :: traj
  type(ctrl_type) :: ctrl
  integer::is,js
  real*8::dentemp !like density matrix

  !write(*,*)"Calculate total energy in QD basis"
  !write(*,*)"Ekin=",traj%Ekin
  !write(*,*)"mapping variables=",qF,pF,qB,pB
  !write(*,*)"Ha=",traj%H_MCH_ss
  traj%Epot_qd=0.d0
  do is=1,ctrl%nstates
     do js=1,ctrl%nstates
        if(is/=js)then
           ActionMat(is,js)=0.5d0*(qF(is)*qF(js)+pF(is)*pF(js))
        else
          ActionMat(is,js)=0.5d0*(qF(is)*qF(js)+pF(is)*pF(js)-2.d0*gamm)
        endif 
        traj%Epot_qd=traj%Epot_qd+real(traj%H_MCH_ss(is,js))*ActionMat(is,js)
     enddo
  enddo
  traj%Etot_qd=traj%Ekin+traj%Epot_qd

endsubroutine Calculate_etot_qd


