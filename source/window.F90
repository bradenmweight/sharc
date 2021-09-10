!
!...calculate electronic population using window function
!
subroutine window(traj,ctrl)
use definitions
implicit none
type(trajectory_type) :: traj
type(ctrl_type) :: ctrl
integer::i,j
real*8,allocatable::xn(:) !action
real*8::val
real*8::thresh
real*8::valnew
real*8::N
real*8,external::stepF
real*8,external::stepFNew

allocate(xn(ctrl%nstates))

#ifdef Model

!...transform QD (adiabatic basis) mapping variables to diabatic mapping variables
call QD2diabatic(traj,ctrl)
write(43,'(i5,4(f16.8,1x))')traj%step,qFd,pFd

do i=1,ctrl%nstates
  xn(i)=0.5d0*(qFd(i)**2+pFd(i)**2)-gamm
enddo

do i=1,ctrl%nstates
   val=1.d0
   do j=1,ctrl%nstates
      if(j/=i)then
         val=val*stepF(gamm,0.d0,xn(j))
      else
         val=val*stepF(gamm,1.d0,xn(j))
      endif
   enddo
   histoF(i)=val
enddo

#else

!In real molecular system, we only get adiabatic electronic population.
!Because we do not get wave function explicitly.
! From Cotton, Miller 2016


if (ctrl%gammfunc .eq. "no") then ! If not adjusted, use original binning procedure
!!! For Rectangular Window
if (ctrl%windowfunc .eq. "rectangle") then
gamm=(sqrt(3.d0)-1.d0)*0.5d0

do i=1,ctrl%nstates
  xn(i)=0.5d0*(qF(i)**2+pF(i)**2)-gamm ! Eq. 3a
enddo

! 2-State Example
! W1(n1,n2) = h(gamm - |n1 - 1|) h(gamm - |n2|)
! W2(n1,n2) = h(gamm - |n1|) h(gamm - |n2| - 1)

do i=1,ctrl%nstates
   val=1.d0
   do j=1,ctrl%nstates
      if(j/=i)then
         val=val*stepF(gamm,0.d0,xn(j)) ! Eqs. 6--8
      else
         val=val*stepF(gamm,1.d0,xn(j)) ! Eqs. 6--8
      endif
   enddo
   histoF(i)=val
enddo

endif

!!! For Triangular Window
! This window only works for NSTATES = 2
if (ctrl%windowfunc .eq. "triangle") then
if (ctrl%nstates .eq. 2) then
gamm = 0.33333333

do i=1,ctrl%nstates
  xn(i)=0.5d0*(qF(i)**2+pF(i)**2)-gamm ! Eq. 3a
enddo

histoF(1) = stepFNew(xn(1) + gamm - 1) * stepFNew(xn(2) + gamm) * stepFNew(2 - 2 * gamm - xn(1) - xn(2))
histoF(2) = stepFNew(xn(2) + gamm - 1) * stepFNew(xn(1) + gamm) * stepFNew(2 - 2 * gamm - xn(1) - xn(2))
write(*,*) "TEST", xn(1), xn(2), gamm, histoF(1), histoF(2)

else
   write(*,*) "For triangle window, nstates must be 2."
   stop 1
endif ! Number of States Check
endif ! Window Type

!!! For Pyramidal Window
! This is a unit right pyramid
! This window works for any NSTATES
! Eqs. 6a-6c, Cotton and Miller, JCP 15, 104101 (2019)
if (ctrl%windowfunc .eq. "pyramidal") then
gamm = 0.33333333

! Compute propagated positive-definite action, ek
      do i=1,ctrl%nstates
         ek(i)=0.5d0*qF(i)**2 + 0.5d0*pF(i)**2 ! Eq. 8 Cotton/Miller Adjusted Gamma
      enddo

      ! Perform Widowing Procedure
      do i=1, ctrl%nstates
         histoF(i)=1
         do j=1, ctrl%nstates
            ! Check if: (j == i AND ek < 1) OR (j != i AND ek >= 1)
            if ( (j .eq. i .and. ek(j) < 1.0d0) .or. (j .ne. i .and. ek(j) >= 1.0d0) ) then
               histoF(i) = 0
            endif
         enddo !nstates j
      enddo !nstates i

endif ! Window Type

else ! If gamm function = "adjusted"
   if (ctrl%windowfunc .eq. "pyramidal") then
      ! Compute propagated positive-definite action, ek
      do i=1,ctrl%nstates
         ek(i)=0.5d0*qF(i)**2 + 0.5d0*pF(i)**2 ! Eq. 8 Cotton/Miller Adjusted Gamma
         !write(*,*) "k, qF,pF, ek:", i, qF(i), pF(i), ek(i) 
      enddo

      ! Perform Widowing Procedure
      do i=1, ctrl%nstates
         histoF(i)=1
         do j=1, ctrl%nstates
            ! Check if: (j == i AND ek < 1) OR (j != i AND ek >= 1)
            if ( (j .eq. i .and. ek(j) < 1.0d0) .or. (j .ne. i .and. ek(j) >= 1.0d0) ) then
               histoF(i) = 0
               !write(*,*)"i, j, histoF(i)", i, j, histoF(i)
            endif
         enddo !nstates j
      enddo !nstates i

   endif ! Pyramidal

   if (ctrl%windowfunc .eq. "rectangle") then
      ! Compute propagated positive-definite action, ek
      do i=1,ctrl%nstates
         ek(i)=0.5d0*qF(i)**2 + 0.5d0*pF(i)**2 ! Eq. 8 Cotton/Miller Adjusted Gamma
      enddo

      ! Perform Widowing Procedure
      do i=1, ctrl%nstates
         histoF(i)=1
         do j=1, ctrl%nstates
            ! Check if: (j == i AND 0 < ek -1 < 2*0.366) OR (j != i AND 0 < ek < 2*0.366) ! Hard-coded original gamma=0.366, not gamma_k

            if (j .eq. i) then
               if (ek(j) - 1 < 0 .or. ek(j) - 1 > 2*0.366) then
                  histoF(i) = 0
               endif
            endif
            if (j .ne. i) then
               if (ek(j) < 0 .or. ek(j) > 2*0.366) then
                  histoF(i) = 0
               endif
            endif
         enddo !nstates j
      enddo !nstates i

   endif ! Rectangle
endif ! Gamm Adjusted?



#endif

deallocate(xn)
end subroutine window
!
!
real*8 function stepF(gamm,N,acti)
implicit none
double precision,intent(in)::gamm,N,acti
real*8::val

val=gamm-abs(acti-N)
if(val>0.d0)then
   stepF=1.d0
else
   stepF=0.d0
endif
end function stepF
!
!
real*8 function stepFNew(val)
!!! BRADEN WEIGHT ~ 11/13/2020
implicit none
real*8::val

if(val>0.d0)then
   stepFNew=1.d0
else
   stepFNew=0.d0
endif
end function stepFNew


