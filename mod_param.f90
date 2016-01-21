module mod_constants
  USE MOD_KINDS
  implicit none
  
  private
  
  public :: zero, one, two, three, four, five, six, seven, eight, nine
  public :: ten, eleven
  public :: half, third, fourth, fifth, sixth, two_third, four_third
  public :: three_fourth, twelfth, pi, one_twentyfourth
  
  real(R8), parameter :: zero            = 0.0_R8,                 &
                        one               = 1.0_R8,                 &
                        two               = 2.0_R8,                 &
                        three             = 3.0_R8,                 &
                        four              = 4.0_R8,                 &
                        five              = 5.0_R8,                 &
                        six               = 6.0_R8,                 &
                        seven             = 7.0_R8,                 &
                        eight             = 8.0_R8,                 &
                        nine              = 9.0_R8,                 &
                        ten               = 10.0_R8,                &
                        eleven            = 11.0_R8,                &
                        half              = 0.5_R8,                 &
                        third             = 1.0_R8/ 3.0_R8,         &
                        fourth            = 1.0_R8/ 4.0_R8,         &
                        fifth             = 1.0_R8/ 5.0_R8,         &
                        sixth             = 1.0_R8/ 6.0_R8,         &
                        two_third         = 2.0_R8/ 3.0_R8,         &
                        four_third        = 4.0_R8/ 3.0_R8,         &
                        three_fourth      = 3.0_R8/ 4.0_R8,         &
                        twelfth           = 1.0_R8/12.0_R8,         &
                        one_twentyfourth  = 1.0_R8/24.0_R8
  real(R8), parameter :: pi = 3.141592653589793238_R8

end module mod_constants

MODULE MOD_GLOBAL
    USE MOD_KINDS    
    implicit none
    REAL(R8),PARAMETER                    :: REARTH=6378.1e+3_R8,pi=3.1415926_R8,SMAL=1.E-6_R8,depth_criterion=0.2_R8
    INTEGER, PARAMETER                   :: KB=3,NMAX=1100
    CHARACTER(LEN=100)                    :: ALLNAME, SUBNAME  
    REAL(R8)                               :: TPRNI,UMOL,GRAV,RAMP,TBIAS,SBIAS,PAI,SMALL,RFE,RFW,RFN,RFS
    INTEGER                               :: I,J,K       
    REAL(R8)                               :: SMOTH,HORCON  
    INTEGER, PARAMETER                   :: NUMPHYTO=1,NUMZOO=1,NUMNUTRI=2,NUMOTHER=1

    REAL(R8)                               :: TIMEB,TIMEE,TIMESTEP,TIMEALL ! IN HOUR
    INTEGER                               :: NT,NTB,NTE,IPRINT1,IPRINT3,NUMREC

    INTEGER                               :: MODE,NYRESTART
    REAL(R8)                               :: DAY0
    REAL(R8)                               :: ICOUPLE,DTI

!   interface_value
    REAL(R8)                               :: XB_SUB,YB_SUB,XE_SUB,YE_SUB  !SET AT INTERFACE
    REAL(R8)                               :: XB_ALL,YB_ALL,XE_ALL,YE_ALL
    INTEGER                               :: IEXT_ALL,IEXT_SUB,ISPLIT_SUB,ISPLIT_ALL
    INTEGER                               :: IMALL,JMALL,IMSUB,JMSUB

    REAL(R8)                               :: DTE_SUB,DTI_SUB   !INPUT AT RUNDATA
    REAL(R8)                               :: DTE_ALL,DTI_ALL   !INPUT AT RUNDATA
    INTEGER                               :: I_SUB,J_SUB,IB_ALL,JB_ALL,IE_ALL,JE_ALL  !   
    REAL(R8),ALLOCATABLE,DIMENSION(:,:)  :: ELF_I,UAF_I,VAF_I ,ELF_HALF  !ELF_I(I_SUB,J_SUB))=ELF(IB_ALL:IE_ALL,JB_ALL:JE_ALL)
    REAL(R8),ALLOCATABLE,DIMENSION(:,:,:) :: U_I,V_I,W_I,Q2_I,Q2L_I,T_I,S_I
    REAL(R8),ALLOCATABLE,DIMENSION(:,:)   :: ELF_B,UAF_B,VAF_B
    REAL(R8),ALLOCATABLE,DIMENSION(:,:,:) :: U_B,V_B,W_B,Q2_B,Q2L_B,T_B,S_B
    REAL(R8),ALLOCATABLE,DIMENSION(:,:)   :: ELF_a,UAF_a,VAF_a
    REAL(R8),ALLOCATABLE,DIMENSION(:,:,:) :: U_a,V_a,W_a,Q2_a,Q2L_a,T_a,S_a
    REAL(R8),ALLOCATABLE,DIMENSION(:,:)   :: CON1F_I,CON1F_B,CON1F_a
    REAL(R8),ALLOCATABLE,DIMENSION(:,:)   :: CON2F_I,CON2F_B,CON2F_a
    REAL(R8),ALLOCATABLE,DIMENSION(:,:)   :: CON3F_I,CON3F_B,CON3F_a

    TYPE TIME
        INTEGER            :: YEAR,MONTH,DAY,HOUR,MINUTE ! IN HOUR
    END TYPE TIME


END MODULE MOD_GLOBAL

MODULE MOD_ALLDOMAIN
    ! VAR IN VARIABLE.H
    USE MOD_KINDS
    USE MOD_GLOBAL,ONLY: NMAX,KB
    implicit none
    INTEGER, PARAMETER :: NP_BCMAX=NMAX, NUMHC=16,Nmo=1, &                  !12
                          NREC_T=Nmo,NREC_S=Nmo, &
                          NREC_U=Nmo,NREC_V=Nmo,NREC_EL=Nmo,NREC_UA=Nmo,NREC_VA=Nmo, &
                          NREC_N=Nmo,NREC_P=Nmo,NREC_Z=Nmo,NREC_O=Nmo,              &   
                          NREC_DIC=Nmo,NREC_ALK=Nmo,NREC_TIDE=80000 
    INTEGER, PARAMETER :: NREC_TAUX=Nmo,NREC_TAUY=Nmo,NREC_NETHEAT=Nmo,NREC_SSS=Nmo

    INTEGER  :: IINT,IEXT,ISPLIT,ISPADV
    REAL(R8)  :: DTE,ISPI,ISP2I,DTE2,DTI2  !1.E0/FLOAT(ISPLIT)      !,DTI

    INTEGER             :: IMM1,JMM1,KBM1
    INTEGER             :: IMM2,JMM2,KBM2
    INTEGER             :: LIJ,LIJ2,LIJK,LIJKM2
    INTEGER             :: LIJM1,LIJKM1

    REAL, DIMENSION (KB) :: Z,ZZ,DZ,DZZ
    INTEGER,ALLOCATABLE,DIMENSION(:) :: kleft(:),kright(:),ktop(:),kbottom(:)

!---------------- 2-D ARRAYS --------------------------------------
    REAL(R8),ALLOCATABLE,DIMENSION(:,:)  :: H,DX,DY,D,DT,ART,ARU,ARV,CBC,ALON,ALAT, &
      &              COR,WUSURF,WVSURF,WUBOT,WVBOT,WTSURF,WSSURF,TPS,AAM2D, &
      &      UAF,UA,UAB,VAF,VA,VAB,ELF,EL,ELB,PSI,ETF,ET,ETB,FLUXUA,FLUXVA, &
      &      EGF,EGB,Z_TIDE,LON_RHO,LAT_RHO,PM,PN
    INTEGER,ALLOCATABLE,DIMENSION(:,:)  :: DUM,DVM,FSM,DUM0,DVM0,FSM0,MBSTATE
    REAL(R8),ALLOCATABLE,DIMENSION(:,:)  :: UTB,VTB,UTF,VTF,ADVUA,ADVVA,TSURF,SSURF, &
     &                                DRX2D,DRY2D,ADX2D,ADY2D,SWRAD
    REAL(R8),ALLOCATABLE,DIMENSION(:,:)  :: LONU,LATU
    REAL(R8),ALLOCATABLE,DIMENSION(:,:)  :: LONV,LATV
!

!---------------- 3-D ARRAYS --------------------------------------
    REAL(R8),ALLOCATABLE,DIMENSION(:,:,:) :: A,C,EE,GG,UF,VF,KM,KH,KQ,L,Q2,Q2B,AAM, &
      &      Q2L,Q2LB,U,UB,W,V,VB,T,TB,S,SB,RHO,DTEF,RMEAN
    REAL(R8),ALLOCATABLE,DIMENSION(:,:,:)  :: ADVX,ADVY,DRHOX,DRHOY,TCLIM,SCLIM
     
    !REAL(R8),ALLOCATABLE,DIMENSION(:)  :: XLEV

    REAL(R8),ALLOCATABLE,DIMENSION(:,:)  :: PCO2O

    REAL(R8),ALLOCATABLE,DIMENSION(:,:,:) :: DEPTH, TEMPER, SALT, LIGHT
    REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: SURF_RAD,WIND10U,WIND10V,PCO2A,CO2EX

  !  INTEGER, PARAMETER :: IOCCB=1,IOCCE=IM,JOCCB=1,JOCCE=JM,KOCCB=1,KOCCE=KB
  !  INTEGER, PARAMETER :: IMOCC=IOCCE-IOCCB+1,JMOCC=JOCCE-JOCCB+1,KBOCC=KOCCE-KOCCB+1
    INTEGER, PARAMETER :: IOCCB=1,JOCCB=1,KOCCB=1,KOCCE=KB    
    INTEGER             :: IOCCE,JOCCE,IMOCC,JMOCC,KBOCC     
    INTEGER, DIMENSION (NMAX)        ::    IIB,JJB,IIC,JJC,IIU,JJU,IIV,JJV
    INTEGER, DIMENSION (NMAX,KB)     ::    TBC, SBC
    INTEGER, DIMENSION (NMAX,KB)     ::    UBC, VBC !UBC,VBC: 3-D VELOCITY ON BOUNDARY
    INTEGER, DIMENSION (NMAX)        ::    UBCA, VBCA, ELA !UBCA,VBCA: 2-D VELOCITY ON BOUNDARY
    INTEGER                          ::    NUMEB,numeb_c,numeb_b,numeb_l
    REAL(R8), DIMENSION (NMAX,NUMHC) ::    AMP,PHASE

    !REAL(R8),ALLOCATABLE,DIMENSION(:,:,:) :: BCDATA_T
    !REAL(R8),ALLOCATABLE,DIMENSION(:,:,:) :: BCDATA_S
    !REAL(R8),ALLOCATABLE,DIMENSION(:,:,:) :: BCDATA_U
    !REAL(R8),ALLOCATABLE,DIMENSION(:,:,:) :: BCDATA_V
    !REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: BCTEMP_T,BCTEMP_S,BCTEMP_U,BCTEMP_V
    !REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: BCDATA_EL
    !REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: BCDATA_UA
    !REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: BCDATA_VA
    !REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: BCTEMP_EL,BCTEMP_UA,BCTEMP_VA
    !REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: BCDATA_TAUX
    !REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: BCDATA_TAUY
    !REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: BCDATA_NETHEAT
    !REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: BCDATA_SSS,BCDATA_SST

    CHARACTER,ALLOCATABLE,DIMENSION(:,:) :: TYP,typ1_c,typ2_c,typ3_c
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: NTYP
    INTEGER :: NUMTYP
    INTEGER,ALLOCATABLE,DIMENSION(:) :: ITYP,JTYP

   !  REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: CON1B,CON1,CON1F,COND1,CON1_ERO,CON1_DEP
 !   REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: CON2B,CON2,CON2F,cond2,CON2_ERO,CON2_DEP
 !   REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: CON3B,CON3,CON3F,Cond3,CON3_ERO,CON3_DEP
    REAL(R8),ALLOCATABLE,DIMENSION(:,:)  :: CONB,CON,CONF,COND,CON_ERO,CON_DEP
    REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: ELR,ELFR,FSMR

    INTEGER ele,nod
    INTEGER,ALLOCATABLE,DIMENSION(:) :: IN,JN,IE,JE

END MODULE MOD_ALLDOMAIN

MODULE MOD_SUBDOMAIN
    USE MOD_KINDS
    USE MOD_GLOBAL,ONLY : NMAX,KB
    implicit NONE

    INTEGER, PARAMETER :: NP_BCMAX=NMAX, NUMHC=16,Nmo=1,                             &  !12
           &                NREC_T=Nmo,NREC_S=nmo,                                     &
           &                NREC_U=Nmo,NREC_V=Nmo,NREC_EL=Nmo,NREC_UA=Nmo,NREC_VA=Nmo, &
           &                NREC_N=Nmo,NREC_P=Nmo,NREC_Z=Nmo,NREC_O=Nmo,               &   
           &                NREC_DIC=Nmo,NREC_ALK=Nmo,NREC_TIDE=80000 
    INTEGER, PARAMETER :: NREC_TAUX=Nmo,NREC_TAUY=Nmo,NREC_NETHEAT=Nmo,NREC_SSS=Nmo
      
    INTEGER   :: IINT,IEXT,ISPLIT,ISPADV
    REAL(R8)   :: DTE,ISPI,ISP2I,DTE2,DTI2  !1.E0/FLOAT(ISPLIT)      !,DTI
    integer   :: im0,jm0
    integer   :: INDEX1,INDEX2,INDEX3   !liuy    
  
    INTEGER   :: IMM1,JMM1,KBM1
    INTEGER   :: IMM2,JMM2,KBM2
    INTEGER   :: LIJ,LIJ2,LIJK,LIJKM2
    INTEGER   :: LIJM1,LIJKM1

    REAL, DIMENSION (KB) :: Z,ZZ,DZ,DZZ
    INTEGER,ALLOCATABLE,DIMENSION(:) :: kleft(:),kright(:),ktop(:),kbottom(:)

!---------------- 2-D ARRAYS --------------------------------------
    REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: H,DX,DY,D,DT,ART,ARU,ARV,CBC,ALON,ALAT, &
    &                         COR,WUSURF,WVSURF,WUBOT,WVBOT,WTSURF,WSSURF,TPS,AAM2D, &
    &                 UAF,UA,UAB,VAF,VA,VAB,ELF,EL,ELB,PSI,ETF,ET,ETB,FLUXUA,FLUXVA, &
    &        EGF,EGB,Z_TIDE,LON_RHO,LAT_RHO,PM,PN
    REAL(R8),ALLOCATABLE,DIMENSION(:,:)  :: LONU,LATU
    REAL(R8),ALLOCATABLE,DIMENSION(:,:)  :: LONV,LATV
    INTEGER, ALLOCATABLE,DIMENSION(:,:) ::  DUM,DVM,FSM,DUM0,DVM0,FSM0,MBSTATE
    REAL(R8),ALLOCATABLE,DIMENSION(:,:)  ::  UTB,VTB,UTF,VTF,ADVUA,ADVVA,TSURF,SSURF, &
    &                           DRX2D,DRY2D,ADX2D,ADY2D,SWRAD

!---------------- 3-D ARRAYS --------------------------------------
    REAL(R8),ALLOCATABLE,DIMENSION(:,:,:) :: A,C,EE,GG,UF,VF,KM,KH,KQ,L,Q2,Q2B,AAM, &
    &        Q2L,Q2LB,U,UB,W,V,VB,T,TB,S,SB,RHO,DTEF,RMEAN
    REAL(R8),ALLOCATABLE,DIMENSION(:,:,:) :: ADVX,ADVY,DRHOX,DRHOY,TCLIM,SCLIM

    !REAL(R8),ALLOCATABLE,DIMENSION(:) :: XLEV

    REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: PCO2O
    REAL(R8),ALLOCATABLE,DIMENSION(:,:,:) :: DEPTH, TEMPER, SALT, LIGHT
    REAL(R8),ALLOCATABLE,DIMENSION(:,:) :: SURF_RAD,WIND10U,WIND10V,PCO2A,CO2EX

    INTEGER, PARAMETER :: IOCCB=1,JOCCB=1,KOCCB=1,KOCCE=KB    
    INTEGER             :: IOCCE,JOCCE,IMOCC,JMOCC,KBOCC 
    !INTEGER, PARAMETER             :: IOCCB=1,IOCCE=IM,JOCCB=1,JOCCE=JM,KOCCB=1,KOCCE=KB
    !INTEGER, PARAMETER             :: IMOCC=IOCCE-IOCCB+1,JMOCC=JOCCE-JOCCB+1,KBOCC=KOCCE-KOCCB+1
                
    INTEGER, DIMENSION (NMAX)      :: IIB,JJB,IIC,JJC,IIU,JJU,IIV,JJV
    INTEGER, DIMENSION (NMAX,NUMHC)::    AMP,PHASE
    INTEGER, DIMENSION (NMAX,KB)   :: TBC, SBC
    INTEGER, DIMENSION (NMAX,KB)   :: UBC, VBC !UBC,VBC: 3-D VELOCITY ON BOUNDARY
    INTEGER, DIMENSION (NMAX)      :: UBCA, VBCA, ELA !UBCA,VBCA: 2-D VELOCITY ON BOUNDARY
    INTEGER                        :: NUMEB,numeb_c,numeb_b,numeb_l

    REAL(R8),ALLOCATABLE,DIMENSION(:,:,:)  :: BCDATA_TAUX
    REAL(R8),ALLOCATABLE,DIMENSION(:,:,:)  :: BCDATA_TAUY
    REAL(R8),ALLOCATABLE,DIMENSION(:,:,:)  :: BCDATA_NETHEAT
    REAL(R8),ALLOCATABLE,DIMENSION(:,:,:)  :: BCDATA_SSS,BCDATA_SST

    CHARACTER,ALLOCATABLE,DIMENSION(:,:) :: TYP
    INTEGER,ALLOCATABLE,DIMENSION(:,:)  :: NTYP
    INTEGER                     :: NUMTYP
    INTEGER,ALLOCATABLE,DIMENSION(:)  :: ITYP,JTYP

    REAL(R8),ALLOCATABLE,DIMENSION(:,:)  :: CONB,CON,CONF,COND,CON_ERO,CON_DEP
    REAL(R8),ALLOCATABLE,DIMENSION(:,:)  :: ELR,ELFR,FSMR,FSMSZ

END MODULE MOD_SUBDOMAIN


