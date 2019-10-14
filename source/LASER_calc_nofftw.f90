!******************************************
!
!    SHARC Program Suite
!
!    Copyright (c) 2019 University of Vienna
!
!    This file is part of SHARC.
!
!    SHARC is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHARC is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    inside the SHARC manual.  If not, see <http://www.gnu.org/licenses/>.
!
!******************************************

module LASER_calc

 contains

subroutine field_transform(laser_t,type_envelope,field_strength,fwhm, &
                         pulse_begin,pulse_center,pulse_center2, &
                         pulse_end,omega_0,phase,b_1, &
                         b_2,b_3,b_4,dt,t0,Nt,envelope,momentary_frequency)
  
  use LASER_definitions
  
  implicit none
      
  integer :: type_envelope
  integer :: Nt
  integer :: it
  integer :: it2
  integer :: it3
  integer :: it4
  
  real(kind=8) :: field_strength
  real(kind=8) :: fwhm
  real(kind=8) :: pulse_begin
  real(kind=8) :: pulse_center
  real(kind=8) :: pulse_center2
  real(kind=8) :: pulse_end
  real(kind=8) :: omega_0
  real(kind=8) :: phase
  real(kind=8) :: b_1
  real(kind=8) :: b_2
  real(kind=8) :: b_3
  real(kind=8) :: b_4
  real(kind=8) :: dt
  real(kind=8) :: t0
  real(kind=8) :: t
!   real(kind=8) :: E
!   real(kind=8) :: dE
  real(kind=8) :: beta
  real(kind=8) :: envelope(Nt)
  real(kind=8) :: momentary_frequency(Nt)
  real(kind=8) :: laser1
  real(kind=8) :: laser2
  real(kind=8) :: t1
  real(kind=8) :: t2
  real(kind=8) :: t_pos_intercept
  real(kind=8) :: t_neg_intercept
  
  complex(kind=8) :: laser_t(Nt)

  
  complex(kind=8), allocatable :: field(:)
  complex(kind=8), allocatable :: laser_E(:)
  !logical debug -> specifications
    
  complex(kind=8) :: field_strength_chirp  
  real(kind=8) :: a2  
  real(kind=8) :: beta_chirp  
  
  allocate (field(Nt), STAT=allocatestatus)
  if (allocatestatus /= 0) stop "*** Not enough memory ***"
  allocate (laser_E(Nt), STAT=allocatestatus)
  if (allocatestatus /= 0) stop "*** Not enough memory ***"
  
  


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! transform-limited pulse
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Choosing envelope type
  select case (type_envelope)

    !! Gaussian
    case (1)   
      beta = 4.d0*log(2.d0)/(fwhm**2)
      field_strength_chirp = field_strength * sqrt(1./(1. + (0.,1.) * b_2 * 2. * beta))
      a2 = (b_2/(1./(4.*beta**2) + b_2**2))
      beta_chirp = 1./(1./beta + 4.*beta*b_2**2)
      do it = 1, Nt
        t = t0 + (it-1) * dt
        laser_t(it) = field_strength_chirp * exp(-beta_chirp * (t - pulse_center)**2) &
                             * exp((0.D0,1.D0) * (omega_0 * (t - pulse_center) &
                             + a2 / 2.D0 * (t - pulse_center)**2 + phase))
      enddo
   
    !! Sinusoidal
    case (2)
      do it = 1,Nt
        t = t0 + (it-1) * dt
        field(it) = exp((0.d0,1.d0) * (omega_0 * (t - pulse_center) + phase))
        if (t.le.pulse_begin) then
          envelope(it) = 0.0d0
        else if (t.le.pulse_center) then
          envelope(it) = field_strength * sin(pi/2.*(t-pulse_begin)/(pulse_center-pulse_begin))**2
        else if (t.le.pulse_center2) then
          envelope(it) = field_strength
        else if (t.le.pulse_end) then
          envelope(it) = field_strength * cos(pi/2.*(t-pulse_center2)/(pulse_end-pulse_center2))**2
        else
          envelope(it)=0.0d0
        endif
        laser_t(it) = envelope(it) * field(it)
      enddo      
   
  end select


    if (b_1 /= 0. .or. b_2 /= 0. .or. b_3 /= 0. .or. b_4 /= 0.) then
      write(6,*) 
      write(6,*) '*******************************************************************************'
      write(6,*) 'Warning! Version without fftw does not support non-zero b_1, b_3, b_4.'
      write(6,*) 'Linear chirp (b_2 non-zero) only with Gaussian envelope.'
      write(6,*) '*******************************************************************************'
      write(6,*) 
    endif
     
        
       !! writing chirped field in time domain
        if (debug) then
          open(13,file='DEBUG_time_chirped.out')
          do it = 1,Nt
            t = t0 + (it-1) * dt
            write(13,*) t*au2fs, dble(laser_t(it)),aimag(laser_t(it))
          enddo
          close(13)
        endif
  
  !! Beginning of laser field, setting momentary frequency to omega_0 for first half-cycle
  it3 = 1
  do while ( (sign(1.d0,real(laser_t(it3))) == sign(1.d0,real(laser_t(it3+1)))) .and. (it3+1 < Nt) )
    momentary_frequency(it3) = omega_0
    it3 = it3 + 1
  enddo
  momentary_frequency(it3) = omega_0
  it3 = it3 + 1
  !! End of laser field, setting momentary frequency to omega_0 for last half-cycle
  it4 = 1
  momentary_frequency(Nt) = omega_0
  do while ( (sign(1.d0,real(laser_t(Nt-it4))) == sign(1.d0,real(laser_t(Nt-it4+1)))) .and. (Nt-it4 > 1))
    momentary_frequency(Nt-it4) = omega_0
    it4 = it4 + 1
  enddo
  momentary_frequency(1) = omega_0
  !! Remaining laser field, setting momentary frequency to numerical value
  do it = it3,Nt-it4
    if ( sqrt(dble(laser_t(it))**2+aimag(laser_t(it))**2) < 1.d-7 ) then
      momentary_frequency(it) = omega_0
      !print*,it
    else
      !! Find zero crossing (t-axis intercept) for times later than t(it)
      it2 = 0
      do 
        laser1 = dble(laser_t(it+it2))
        t1 = t0 + (it+it2-1) * dt
        laser2 = dble(laser_t(it+it2+1))
        t2 = t0 + (it+it2) * dt
        if ( sign(1.d0,real(laser_t(it+it2+1))) /= sign(1.d0,real(laser_t(it+it2))) ) exit
        it2 = it2 + 1
      enddo
      t_pos_intercept = t2 - laser2 * (t2-t1)/(laser2-laser1)
      !! Find zero crossing (t-axis intercept) for times earlier than t(it)
      it2 = 0
      do 
        laser2 = dble(laser_t(it-it2))
        t2 = t0 + (it-it2-1) * dt
        laser1 = dble(laser_t(it-it2-1))
        t1 = t0 + (it-it2-2) * dt
        if ( sign(1.d0,real(laser_t(it-it2))) /= sign(1.d0,real(laser_t(it-it2-1))) ) exit
        it2 = it2 + 1
      enddo
      t_neg_intercept = t2 - laser2 * (t2-t1)/(laser2-laser1)
      !print*,t_pos_intercept, t_neg_intercept
      momentary_frequency(it) = pi / (t_pos_intercept - t_neg_intercept)
    endif
  enddo
  
  
  deallocate(field)
  deallocate(laser_E)

  return
  
end subroutine field_transform

endmodule