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
    
