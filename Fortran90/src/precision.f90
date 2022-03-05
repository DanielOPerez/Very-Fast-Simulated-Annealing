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

  
  MODULE precision

    ! --------------------------------------------------------------------------
    !  SP: simple precision de la norma IEEE 754
    !  DP: doble precision de la norma IEEE 754
    !
    !  Uso: USE precision, WP => SP o USE precision, WP => DP
    ! --------------------------------------------------------------------------

  
    
    INTEGER,PARAMETER::SP=SELECTED_REAL_KIND(6,37)  
    INTEGER,PARAMETER::DP=SELECTED_REAL_KIND(15,307) 

  END MODULE precision
