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


!-------------------------------------------------------------------------------
! Este modulo contiene todas las subrutinas que hacen al VFSA
!-------------------------------------------------------------------------------

  MODULE sub_vfsa90
 
    USE precision, WP => DP
    USE variables, ONLY: MPAR,NPAR,MNPAR,M_PAR,ITMAX,ITMAX_CON,&
         THRES,XMISFIT,QUENCH,TSCALE,TRS,TAS,IOP_PRT,COST&
         ,VAR,ITMAX,OUT0
    USE fun_cost_mod

    IMPLICIT NONE

    REAL(KIND=WP)::COST_OPT,COST_0,COST_TMP,T_COST,T0_COST,T_MIN,C,&
         C_COST,QD,QD_COST,DELTA_COST,PROB,U
    INTEGER::ICONT,K
 
  CONTAINS

    SUBROUTINE ASA

      !Variables locales de la subrutina ASA, que hace todo el laburo
      INTEGER::i,j,ACC_UP,ACC_DN,REJ,ISTOP,N_SAME_CON,UNIT_VAL
            
      CALL init_random_seed() !inicializo la semilla para RANDOM_NUMBER

      ICONT=0
      DO i=1,MPAR
         DO j=1,NPAR(i)
            IF(VAR(i,j)%XA==VAR(i,j)%XB) ICONT=ICONT+1
         END DO
      END DO

      
      IF (QUENCH>=1.0_WP)THEN
         QD=QUENCH/(MNPAR-ICONT)
         C=-LOG(TRS)*EXP(-LOG(TAS)*QD)
         QD_COST=QD
         C_COST=C*TSCALE
      ELSE IF (QUENCH > 0.0_WP .AND. QUENCH < 1.0_WP) THEN
         QD=1.0_WP/(MNPAR-ICONT)
         C=-LOG(TRS)*EXP(-LOG(TAS)*QD)
         QD_COST=1.0_WP
         C_COST=QUENCH
      ELSE
         QD=1.0_WP
         C=-QUENCH
         QD_COST=QD
         C_COST=C*TSCALE
      END IF

      K=0
      ACC_UP=0
      ACC_DN=0
      REJ=0
      ISTOP=0
      
      CALL COST_FUN(VAR%X0,COST_0) ! se calcula el costo del modelo incial

      !se inicializan las temperaturas
      !--------------------------------------------------------------
      T_MIN=1.0e-150_WP !tempreatura minima
      T0_COST=COST_0    !inicializo la temperatura de costo inicial
      T_COST=T0_COST    !inicializo la temperatura del costo
      VAR%T0=1.0_WP     !inicializo la temperatura inicial
      VAR%T=VAR%T0      !inicializo la temperatura de cada variable
      !--------------------------------------------------------------

      
      VAR%X_TMP=VAR%X0 !VAR%X_TMP es la variable con la que trabajara NEWCONFIG
      VAR%X_OPT=VAR%X0 !VAR%X_OPT es el valor de X que optimiza la funcion  
      COST_OPT=COST_0  !COST_OPT es el valor optimo de la funcion costo

      !En el archivo vfsa90.out se guardan los resultados de ser neceasario
      !se guarda resultadoas cada IOP_PRT iteraciones si es >0
      UNIT_VAL=6 !salida por pantalla por 'default' 
      IF(IOP_PRT > 0) THEN
         UNIT_VAL=10
         OPEN(UNIT=UNIT_VAL,FILE=OUT0,ACTION='write')
      ENDIF
      !ENCABEZADO
      WRITE(UNIT_VAL,100)
    

      DO WHILE(ISTOP==0)
 
         !Print section
         !Solo se muestran resultado en pantalla cada IOP_PRT 
         IF(IOP_PRT/=0)THEN
            IF(MOD(K,IOP_PRT)==0) WRITE(UNIT_VAL,110)&
                 K,ACC_DN+ACC_UP,ACC_DN,ACC_UP,COST_OPT,COST_0,&
                 COST_TMP,T_COST
         END IF
         VAR%OLD=VAR%X_TMP !se guarda temporariamente el valor de 
                           !la variable sin perturbar para poder
                           !ser utilizada por el criterio de Metropolis.
        
         !Nueva configuracion, costo y temperatura
         CALL NEWCONFIG 

         !se aumenta el contador de iteraciones
         K=K+1
    

         !NewT_cost
         IF(T_COST<=T_MIN) THEN
            T_COST=T_MIN
         ELSE
            T_COST=T0_COST*EXP(-C_COST*((K-1)**QD_COST))
         END IF

         !NewT
         WHERE(VAR%T<=T_MIN)
            VAR%T=T_MIN
         ELSEWHERE
            VAR%T=VAR%T0*EXP(-C*((K-1)**QD))
         END WHERE

         CALL COST_FUN(VAR%X_TMP,COST_TMP) 

         DELTA_COST=COST_TMP-COST_0
          
         !criterio de corte 3
         IF(100.0_WP*(ABS(DELTA_COST)/COST_0) <= THRES) THEN
            N_SAME_CON=N_SAME_CON+1
         ELSE
            N_SAME_CON=0
         END IF
    

         !Metropolis
         IF(DELTA_COST<0.0_WP)THEN
            ACC_DN=ACC_DN+1
            COST_0=COST_TMP
            IF(COST_TMP<COST_OPT) THEN
               VAR%X_OPT=VAR%X_TMP
               COST_OPT=COST_TMP
            ENDIF
         ELSE
            PROB=EXP(-DELTA_COST/T_COST)
            CALL RANDOM_NUMBER(U)
            IF(U<PROB) THEN
               ACC_UP=ACC_UP+1
               COST_0=COST_TMP
               IF(COST_TMP<COST_OPT) THEN
                  VAR%X_OPT=VAR%X_TMP
                  COST_OPT=COST_TMP
               ENDIF
            ELSE
               VAR%X_TMP=VAR%OLD
               REJ=REJ+1
            ENDIF
         END IF
      
         !Criterios de corte
         IF(K > ITMAX) ISTOP=1
         IF(COST_OPT < XMISFIT) ISTOP=2
         IF(N_SAME_CON==ITMAX_CON) ISTOP=3

      END DO !FIN DEL PROCESO ITERATIVO!!!
      
      !SE GUARDAN O SE MUESTRAN LOS RESULTADOS FINALES
      WRITE(UNIT_VAL,110)K,ACC_DN+ACC_UP,ACC_DN,ACC_UP,&
           COST_OPT,COST_0,COST_TMP,T_COST
      WRITE(UNIT_VAL,*)''

      SELECT CASE(ISTOP)
      CASE(1)
         WRITE(UNIT_VAL,120) ITMAX 
      CASE(2)
         WRITE(UNIT_VAL,130) XMISFIT
      CASE(3)
         WRITE(UNIT_VAL,140) THRES,ITMAX_CON
      END SELECT

      DO i=1,MPAR
         IF(NPAR(i)/=0)THEN
            WRITE(UNIT_VAL,*)''
            WRITE(UNIT_VAL,80)"RESULTADOS VARIABLE TIPO",i
            DO j=1,NPAR(i)
               WRITE(UNIT_VAL,90)j,VAR(i,j)%X_OPT
            END DO
         END IF
      END DO


      
     
      !======================================================================
      ! FORMATOS DE SALIDA
      !======================================================================
80    FORMAT(A24,I3)
90    FORMAT(I3,F15.7)
100   FORMAT(/,38x,'SA solution',/,/,&
           3x,'gen',3x,'acc    (    dn    up )',8x,'cost_opt',11x,'cost',&
           7x,'cost_tmp',4x,'T_cost',/,2x,85('-'))
110   FORMAT(2i6,'    (',2i6,' ) ',3f15.7,e10.1)

120   FORMAT('CRITERIO DE CORTE: SE ALCANZO EL NUMERO MAXIMO DE ITERACIONES (', i6, ' )')
130   FORMAT('CRITERIO DE CORTE: SE ALCANZO EL VALOR MINIMO DEL MISFIT (', e10.1, ' )')
140   FORMAT('CRITERIO DE CORTE: EL MISFIT VARIO MENOS QUE ', f4.2, ' DURANTE',i5,' ITERACIONES')

    END SUBROUTINE ASA

    !========================================================================
    !Subrutina NEWCONFIG. Se encarga de crear los nuevos valores de las
    !variables en cada iteracion del Annealing
    !========================================================================
    SUBROUTINE NEWCONFIG
      
      !Variables locales 
      INTEGER::i,j
      REAL(KIND=WP)::xaa,xbb

      DO i=1,MPAR
         DO j=1,NPAR(i)
            SELECT CASE(VAR(i,j)%OLAP)
            CASE(0)
            !PARAMETROS SIN OVERLAP, NO SE PUEDEN SUPERPONER
               IF(j==1) xaa=VAR(i,j)%XA 
               IF(j>1) xaa=MAX(VAR(i,j-1)%X_TMP+VAR(i,j-1)%DX,VAR(i,j)%XA)
               IF(j<NPAR(i)) xbb=MIN(VAR(i,j+1)%X_TMP-VAR(i,j+1)%DX,VAR(i,j)%XB)
               IF(j==NPAR(i)) xbb=VAR(i,j)%XB
               
               VAR(i,j)%X_TMP=PERTURB(VAR(i,j)%X_TMP,VAR(i,j)%T,xaa,xbb)

            CASE DEFAULT
               !PARAMETROS CON OVERLAP, SE PUEDEN SUPERPONER
               xaa=VAR(i,j)%XA ; xbb=VAR(i,j)%XB
               
               VAR(i,j)%X_TMP=PERTURB(VAR(i,j)%X_TMP,VAR(i,j)%T,xaa,xbb)

            END SELECT
         END DO
      END DO

    CONTAINS

      !Funcion Perturb. Perturba los valores de VAR entre los limites
      !xaa_IN e xbb_IN

      FUNCTION PERTURB(X_IN,T_IN,xaa_IN,xbb_IN)
        
        IMPLICIT NONE
        REAL(KIND=WP)::xaa_IN,xbb_IN,X_IN,T_IN,PERTURB,U
        INTEGER::i1
        i1=1
        PERTURB=xaa_IN-10.0_WP !Inicializo PERTURB para ingresar al do while
        DO WHILE (PERTURB<xaa_IN.OR.PERTURB>xbb_IN)
           CALL RANDOM_NUMBER(U)
           PERTURB=X_IN+SIGN(1.0_WP,(U-0.5_WP))*T_IN*((1.0_WP+&
                1.0_WP/T_IN)**ABS(2.0_WP*U-1.0_WP)-1.0_WP)*(xbb_IN-xaa_IN)
        END DO
        RETURN
      END FUNCTION PERTURB

    END SUBROUTINE NEWCONFIG

    !========================================================================
    !La subrutina init_random_seed inicializa la semilla con la que trabajara 
    !la subrutina RAMDOM_NUMBER que genera los numeros pseudoaleatorios
    !========================================================================
    SUBROUTINE init_random_seed()
      INTEGER :: i, n, clock
      INTEGER, ALLOCATABLE :: seed(:)
      
      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))
      
      CALL SYSTEM_CLOCK(COUNT=clock)
      
      seed = clock + 1*(/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)
      
      DEALLOCATE(seed)
    END SUBROUTINE init_random_seed


  END MODULE sub_vfsa90
