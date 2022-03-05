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

 
 MODULE variables
    !---------------------------------------------------------------------------
    ! En este modulo se declaran las variables..y alguna que otra subrutina...
    !---------------------------------------------------------------------------

    USE precision, WP => DP
    IMPLICIT NONE
    

    !Parametros de entrada al VFSA------------------------------------

    INTEGER::ITMAX,ITMAX_CON,IFREQ_RE,IFREQ_TSCALE,IOP_PRT,N_PAR
    REAL(KIND=WP)::THRES,XMISFIT,QUENCH,TSCALE,TRS,TAS,ALPHA,PNORM
    CHARACTER(LEN=20)::IN1,IN2,OUT0
    
    !Variables que usa el VFSA
    
    INTEGER::MPAR,MNPAR,M_PAR,OVERLAP
    INTEGER, ALLOCATABLE :: NPAR(:)
    REAL(KIND=WP)::COST,DELTA_X
   
    TYPE :: VAR_ANNEALING
       SEQUENCE
       REAL(KIND=WP)::OLD,X0,X_TMP,X_OPT,DX,XA,XB,XMA,XMB,T,T0
       INTEGER::OLAP
    END type VAR_ANNEALING
       
    TYPE (VAR_ANNEALING), ALLOCATABLE :: VAR(:,:) 

  END MODULE variables
    
