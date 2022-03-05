!-------------------------------------------------------------------------------
! Esta subrutina lee los datos necesarios para que el VFSA funcione.
! Se lee el archivo asa5.in y el archivo asa5.i1:
!
!  FILEIN1:   input data filename (A20). 
!             Search ranges and maximum perturbation sizes.
!  FILEIN2:   initial estimate (OPTIONAL) (A20).
!             Specify 'none' (or 'NONE') if no initial
!             estimate is available.
!  FILEOUT1:  general output filename (A20). Optimization current
!             status is printed here. (see also IOP_PRT)
!  FILEOUT2:  final solution (A20). 
!  ITMAX:     1st stop condition: maximum number of iterations. 
!  ITMAX_CON: 2nd stop condition: maximum number of consecutive
!             iterations CRV<THRES. See THRES also. 
!  THRES:     threshold: if cost relative variation (CRV) between
!             consecutive iterations less than THRES, assume 
!             cost didn't change.
!             -Stop if CRV<THRES, ITMAX_CON consecutive iterations.
!  QUENCH:    quenching factor (accelerate cooling schedule):
!             -if=1, no quenching.
!             -if=N_PAR, maximum quenching.
!             -if small (say 0.01), standard exp. cooling rate.
!             Default: QUENCH=1.0 (slower convergence,
!             but maximum guarantee of finding global minimum).
!  TSCALE:    temperature scale: if too many uphill configurations
!             are accepted, reduce TSCALE (and viceversa).
!             However, its selection is not critical.
!             Default: TSCALE=1.0. 
!  IFREQ_RE:  -if>0, reanneal every IFREQ_RE generated 
!              configurations.
!             -if<0, reanneal every -IFREQ_RE accepted 
!              configurations.
!             -if=0, no reannealing is performed.
!             Defaults: IFREQ_RE=100, -1000.
!  IFREQ_TSCALE: TSCALE is adjusted every IFREQ_TSCALE iterations
!             to keep ratio of accepted configurations
!             uphill/downhill approximately constant.
!  IOP_PRT:   -if<>0, print optimization status every IOP_PRT
!              iterations on FILEOUT1. 
!             -if=0,  none is printed.   
!  ISEED:     negative integer number.           
!
!-------------------------------------------------------------------------------

  MODULE READ_PARAMETERS

    USE precision, WP => DP
    USE variables, ONLY:IN1,IN2,OUT0,ITMAX,ITMAX_CON,THRES,XMISFIT, &
          QUENCH,TSCALE,TRS,TAS,IFREQ_RE,IFREQ_TSCALE,IOP_PRT,&
          MPAR,NPAR,MNPAR,M_PAR,OVERLAP,DELTA_X,VAR
    IMPLICIT NONE

  CONTAINS

    SUBROUTINE READ_PAR

      !variables locales de la subrutina read_par
      
      LOGICAL::LEXIST
      INTEGER::EOF,CTRL_LEC,CPAR,TPAR
      INTEGER::i,j
      INTEGER, ALLOCATABLE::INIC(:)
      REAL(KIND=WP)::XAaux,XBaux,XMAaux,XMBaux
      CHARACTER(LEN=100)::AUX_CHAR1
      
      !Parametros de entrada al VFSA-----------------
      
      INQUIRE(file='op_vfsa90',exist=LEXIST)
      IF(LEXIST)THEN
         OPEN(UNIT=10, FILE='op_vfsa90', ACTION='read')
      ELSE
         WRITE(*,*)'EL ARCHIVO op_vfsa90 NO EXISTE!'
      ENDIF
      
      !read in file names 
      READ (10,*) IN1
      READ (10,*) IN2
      READ (10,*) OUT0
      
      !ASA STUFF:
      READ (10,*) ITMAX,ITMAX_CON,THRES,XMISFIT
      READ (10,*) QUENCH,TSCALE,TRS,TAS
      READ (10,*) IFREQ_RE,IFREQ_TSCALE,IOP_PRT

      
      CLOSE (UNIT=10)
      
      
      !Se leen los valores del archivo asa5.i1--------
      INQUIRE(file=IN1,exist=LEXIST)
      IF(LEXIST)THEN
         OPEN (UNIT=10, FILE=IN1, ACTION='read')
      ELSE
         WRITE(*,*)'EL ARCHIVO',IN1,' NO EXISTE!'
      ENDIF
      
      !se ve cual es el numero maximo de tipos de parametros con los que se
      !va a trabajar
      MPAR=0
      READ(10,*)AUX_CHAR1
      DO
         READ(10,*,IOSTAT=EOF)CTRL_LEC,CPAR,TPAR
         IF(EOF/=0)EXIT
         IF(TPAR > MPAR) MPAR=TPAR
      END DO
      REWIND(10)
          
      !Asigno dinamicamente la memoria de los arreglos NPAR e INIC
      ALLOCATE(NPAR(MPAR),INIC(MPAR))
      
      !Inicializo el arreglo NPAR e INIC
      NPAR=0
      INIC=1
      
      !se cuentan cuanto parametros hay de cada tipo
      READ(10,*)AUX_CHAR1
      DO
         READ(10,*,IOSTAT=EOF)CTRL_LEC,CPAR,TPAR
         IF(EOF/=0)EXIT
         IF(CTRL_LEC/=0) THEN
            NPAR(TPAR)=NPAR(TPAR)+CPAR
          END IF
      END DO

      !Se calcula el numero total de parametros que hay
      MNPAR=SUM(NPAR)
      M_PAR=MAXVAL(NPAR)

      !Se asigna dinamicamente el espacio para las variables VAR
      !son matrices de dimension TPARxMNPAR, a su vez cada variable
      !tiene la informacion de XA,XB,XMA,XMB,X0,X_OPT
      
      ALLOCATE(VAR(MPAR,M_PAR))
      
      !inicializo la variable VAR
      VAR%X0=0.0_WP
      VAR%X_TMP=0.0_WP
      VAR%X_OPT=0.0_WP
      VAR%XA=0.0_WP
      VAR%XB=0.0_WP
      VAR%XMA=0.0_WP
      VAR%XMB=0.0_WP
      
      REWIND(UNIT=10) !Rebobino el archivo param_vfsa90
      
      !Se lee el resto de la informacion y se crea el modelo inicial.
      READ(10,*)AUX_CHAR1
      DO
         READ(10,*,IOSTAT=EOF)CTRL_LEC,CPAR,TPAR,OVERLAP,DELTA_X,&
              XAaux,XBaux,XMAaux,XMBaux
         IF(EOF/=0)EXIT
         IF(CTRL_LEC/=0) THEN
            
            DO i=INIC(TPAR),CPAR+INIC(TPAR)-1
             
               VAR(TPAR,i)%OLAP=OVERLAP
               VAR(TPAR,i)%DX=DELTA_X
               VAR(TPAR,i)%XA=XAaux
               VAR(TPAR,i)%XB=XBaux
               VAR(TPAR,i)%XMA=XMAaux
               VAR(TPAR,i)%XMB=XMBaux
               VAR(TPAR,i)%X0=((XBaux-XAaux)/(CPAR+1))*(i-INIC(TPAR)+1)+XAaux
               
            ENDDO
            INIC(TPAR)=i
            
         END IF
         
      END DO

      IF(IOP_PRT < 0)THEN
         !Modelo inicial---------------------------------------------
         DO i=1,MPAR
            IF(NPAR(i)/=0)THEN
               WRITE(*,*)''
               WRITE(*,100)'MODELO INICIAL DE LOS PARAMETROS DEL TIPO',i
               DO j=1,NPAR(i)
                  WRITE(*,*)j,VAR(i,j)%X0
               ENDDO
            END IF
         ENDDO
100      FORMAT(A45,1X,I3) 
      END IF
      
      DEALLOCATE(INIC)
      

    END SUBROUTINE READ_PAR

  END MODULE READ_PARAMETERS
