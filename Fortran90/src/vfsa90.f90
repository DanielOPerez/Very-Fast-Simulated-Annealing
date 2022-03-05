  PROGRAM VFSA90

    USE precision, WP => DP
    USE variables, ONLY:MPAR,NPAR,VAR
    USE read_parameters
    USE sub_vfsa90

  

    IMPLICIT NONE
    REAL(KIND=WP)::COST_OUT
    integer::i


    
    !Se leen los valores del archivo asa5.in y de IN1
    !Estos archivos contienen los parametros del VFSA
    !Tambien se crea el modelo inicial.

    CALL READ_PAR

    CALL COST_FUN_AUX

    CALL ASA

    
           
    DEALLOCATE(NPAR,&
         VAR)

    

    

  END PROGRAM VFSA90
