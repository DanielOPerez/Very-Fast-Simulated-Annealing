! This program is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
! 
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with
! this program. If not, see <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------
! Created By  : Daniel O. Perez
! Created Date: between 2010 and 2015 
! email: perez.daniel.omar@gmail.com
! ---------------------------------------------------------------------------

 
 !-----------------------------------------------------------------------------
  ! En este modulo se define la funcion costo y las subrutinas y variables
  ! necesarias para su calculo
  !-----------------------------------------------------------------------------

  module fun_cost_mod

    use precision, WP => DP
    use variables, ONLY:MPAR,NPAR,MNPAR,VAR
    use count_file_lines_mod

    implicit none
    
    
    
  contains
  
    subroutine COST_FUN(VAR_IN,COST_OUT)

      implicit none
      REAL(KIND=WP), DIMENSION(:,:)::VAR_IN
      REAL(KIND=WP), INTENT(OUT)::COST_OUT

      integer::i
      real(kind=wp)::A,pi

      
      pi=4.0_wp*atan(1.0_wp)
      A=10.0_wp
      COST_OUT=0.0_wp
      
      do i=1,NPAR(1)
         COST_OUT=COST_OUT+(VAR_IN(1,i)**2-A*cos(2*pi*VAR_IN(1,i)))
      end do
      COST_OUT=COST_OUT+A*NPAR(1)
      
    end subroutine COST_FUN
    

    !==========================================================================

    subroutine COST_FUN_AUX

      
      
      
    end subroutine COST_FUN_AUX

    
    
  end module fun_cost_mod

