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

