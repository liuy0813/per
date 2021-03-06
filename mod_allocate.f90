!To allocate variable used in the model
MODULE MOD_ALLOCATE
USE MOD_GLOBAL
IMPLICIT NONE
CONTAINS
SUBROUTINE MYALLOCATE_ALL
USE MOD_ALLDOMAIN
USE MOD_GLOBAL,ONLY:KB
IMPLICIT NONE
ALLOCATE(kleft  (JMALL))
ALLOCATE(kright (JMALL))
ALLOCATE(ktop   (IMALL))
ALLOCATE(kbottom(IMALL))
ALLOCATE(H(IMALL,JMALL))
ALLOCATE(DX(IMALL,JMALL))
ALLOCATE(DY(IMALL,JMALL))
ALLOCATE(D(IMALL,JMALL))
ALLOCATE(DT(IMALL,JMALL))
ALLOCATE(ART(IMALL,JMALL))
ALLOCATE(ARU(IMALL,JMALL))
ALLOCATE(ARV(IMALL,JMALL))
ALLOCATE(CBC(IMALL,JMALL))
ALLOCATE(ALON(IMALL,JMALL))
ALLOCATE(ALAT(IMALL,JMALL))
ALLOCATE(COR(IMALL,JMALL))
ALLOCATE(WUSURF(IMALL,JMALL))
ALLOCATE(WVSURF(IMALL,JMALL))
ALLOCATE(WUBOT(IMALL,JMALL))
ALLOCATE(WVBOT(IMALL,JMALL))
ALLOCATE(WTSURF(IMALL,JMALL))
ALLOCATE(WSSURF(IMALL,JMALL))
ALLOCATE(TPS(IMALL,JMALL))
ALLOCATE(AAM2D(IMALL,JMALL))
ALLOCATE(UAF(IMALL,JMALL))
ALLOCATE(UA(IMALL,JMALL))
ALLOCATE(UAB(IMALL,JMALL))
ALLOCATE(VAF(IMALL,JMALL))
ALLOCATE(VA(IMALL,JMALL))
ALLOCATE(VAB(IMALL,JMALL))
ALLOCATE(ELF(IMALL,JMALL))
ALLOCATE(EL(IMALL,JMALL))
ALLOCATE(ELB(IMALL,JMALL))
ALLOCATE(PSI(IMALL,JMALL))
ALLOCATE(ETF(IMALL,JMALL))
ALLOCATE(ET(IMALL,JMALL))
ALLOCATE(ETB(IMALL,JMALL))
ALLOCATE(FLUXUA(IMALL,JMALL))
ALLOCATE(FLUXVA(IMALL,JMALL))
ALLOCATE(EGF(IMALL,JMALL))
ALLOCATE(EGB(IMALL,JMALL))
ALLOCATE(Z_TIDE(IMALL,JMALL))
ALLOCATE(LON_RHO(IMALL,JMALL))
ALLOCATE(LAT_RHO(IMALL,JMALL))
ALLOCATE(PM(IMALL,JMALL))
ALLOCATE(PN(IMALL,JMALL))
ALLOCATE(LONU(IMALL-1,JMALL))
ALLOCATE(LATU(IMALL-1,JMALL))
ALLOCATE(LONV(IMALL,JMALL-1))
ALLOCATE(LATV(IMALL,JMALL-1))

ALLOCATE(DUM(IMALL,JMALL))
ALLOCATE(DVM(IMALL,JMALL))
ALLOCATE(FSM(IMALL,JMALL))
ALLOCATE(DUM0(IMALL,JMALL))
ALLOCATE(DVM0(IMALL,JMALL))
ALLOCATE(FSM0(IMALL,JMALL))
ALLOCATE(MBSTATE(IMALL,JMALL))
ALLOCATE(UTB(IMALL,JMALL))
ALLOCATE(VTB(IMALL,JMALL))
ALLOCATE(UTF(IMALL,JMALL))
ALLOCATE(VTF(IMALL,JMALL))
ALLOCATE(ADVUA(IMALL,JMALL))
ALLOCATE(ADVVA(IMALL,JMALL))
ALLOCATE(TSURF(IMALL,JMALL))
ALLOCATE(SSURF(IMALL,JMALL))
ALLOCATE(DRX2D(IMALL,JMALL))
ALLOCATE(DRY2D(IMALL,JMALL))
ALLOCATE(ADX2D(IMALL,JMALL))
ALLOCATE(ADY2D(IMALL,JMALL))
ALLOCATE(SWRAD(IMALL,JMALL))
!3D VAR
ALLOCATE(A(IMALL,JMALL,KB))
ALLOCATE(C(IMALL,JMALL,KB))
ALLOCATE(EE(IMALL,JMALL,KB))
ALLOCATE(GG(IMALL,JMALL,KB))
ALLOCATE(UF(IMALL,JMALL,KB))
ALLOCATE(VF(IMALL,JMALL,KB))
ALLOCATE(KM(IMALL,JMALL,KB))
ALLOCATE(KH(IMALL,JMALL,KB))
ALLOCATE(KQ(IMALL,JMALL,KB))
ALLOCATE(L(IMALL,JMALL,KB))
ALLOCATE(Q2(IMALL,JMALL,KB))
ALLOCATE(Q2B(IMALL,JMALL,KB))
ALLOCATE(AAM(IMALL,JMALL,KB))
ALLOCATE(Q2L(IMALL,JMALL,KB))
ALLOCATE(Q2LB(IMALL,JMALL,KB))
ALLOCATE(U(IMALL,JMALL,KB))
ALLOCATE(UB(IMALL,JMALL,KB))
ALLOCATE(W(IMALL,JMALL,KB))
ALLOCATE(V(IMALL,JMALL,KB))
ALLOCATE(VB(IMALL,JMALL,KB))
ALLOCATE(T(IMALL,JMALL,KB))
ALLOCATE(TB(IMALL,JMALL,KB))
ALLOCATE(S(IMALL,JMALL,KB))
ALLOCATE(SB(IMALL,JMALL,KB))
ALLOCATE(RHO(IMALL,JMALL,KB))
ALLOCATE(DTEF(IMALL,JMALL,KB))
ALLOCATE(RMEAN(IMALL,JMALL,KB))

ALLOCATE(ADVX(IMALL,JMALL,KB))
ALLOCATE(ADVY(IMALL,JMALL,KB))
ALLOCATE(DRHOX(IMALL,JMALL,KB))
ALLOCATE(DRHOY(IMALL,JMALL,KB))
ALLOCATE(TCLIM(IMALL,JMALL,KB))
ALLOCATE(SCLIM(IMALL,JMALL,KB))

ALLOCATE(PCO2O(IMALL,JMALL))
ALLOCATE(DEPTH(IMALL,JMALL,KB))
ALLOCATE(TEMPER(IMALL,JMALL,KB))
ALLOCATE(SALT(IMALL,JMALL,KB))
ALLOCATE(LIGHT(IMALL,JMALL,KB))

ALLOCATE(SURF_RAD(IMALL,JMALL))
ALLOCATE(WIND10U(IMALL,JMALL))
ALLOCATE(WIND10V(IMALL,JMALL))
ALLOCATE(PCO2A(IMALL,JMALL))
ALLOCATE(CO2EX(IMALL,JMALL))

ALLOCATE(TYP(IMALL,JMALL))
ALLOCATE(NTYP(IMALL,JMALL))
ALLOCATE(ITYP(IMALL*JMALL))
ALLOCATE(JTYP(IMALL*JMALL))

ALLOCATE(CONB(IMALL,JMALL))
ALLOCATE(CON(IMALL,JMALL))
ALLOCATE(CONF(IMALL,JMALL))
ALLOCATE(COND(IMALL,JMALL))
ALLOCATE(CON_ERO(IMALL,JMALL))
ALLOCATE(CON_DEP(IMALL,JMALL))

ALLOCATE(ELR(IMALL,JMALL))
ALLOCATE(ELFR(IMALL,JMALL))
ALLOCATE(FSMR(IMALL,JMALL))

IMM1=IMALL-1
JMM1=JMALL-1
KBM1=KB-1
IMM2=IMALL-2
JMM2=JMALL-2
KBM2=KB-2
LIJ=IMALL*JMALL
LIJ2=LIJ*2
LIJK=LIJ*KB
LIJKM2=LIJ*KBM2
LIJM1=IMALL*(JMALL-1)
LIJKM1=LIJ*KBM1
IOCCE=IMALL
JOCCE=JMALL
IMOCC=IOCCE-IOCCB+1
JMOCC=JOCCE-JOCCB+1
KBOCC=KOCCE-KOCCB+1
END SUBROUTINE MYALLOCATE_ALL

SUBROUTINE MYALLOCATE_SUB
USE MOD_SUBDOMAIN
USE MOD_GLOBAL,ONLY:KB
IMPLICIT NONE
ALLOCATE(kleft  (JMSUB))
ALLOCATE(kright (JMSUB))
ALLOCATE(ktop   (IMSUB))
ALLOCATE(kbottom(IMSUB))
ALLOCATE(H(IMSUB,JMSUB))
ALLOCATE(DX(IMSUB,JMSUB))
ALLOCATE(DY(IMSUB,JMSUB))
ALLOCATE(D(IMSUB,JMSUB))
ALLOCATE(DT(IMSUB,JMSUB))
ALLOCATE(ART(IMSUB,JMSUB))
ALLOCATE(ARU(IMSUB,JMSUB))
ALLOCATE(ARV(IMSUB,JMSUB))
ALLOCATE(CBC(IMSUB,JMSUB))
ALLOCATE(ALON(IMSUB,JMSUB))
ALLOCATE(ALAT(IMSUB,JMSUB))
ALLOCATE(COR(IMSUB,JMSUB))
ALLOCATE(WUSURF(IMSUB,JMSUB))
ALLOCATE(WVSURF(IMSUB,JMSUB))
ALLOCATE(WUBOT(IMSUB,JMSUB))
ALLOCATE(WVBOT(IMSUB,JMSUB))
ALLOCATE(WTSURF(IMSUB,JMSUB))
ALLOCATE(WSSURF(IMSUB,JMSUB))
ALLOCATE(TPS(IMSUB,JMSUB))
ALLOCATE(AAM2D(IMSUB,JMSUB))
ALLOCATE(UAF(IMSUB,JMSUB))
ALLOCATE(UA(IMSUB,JMSUB))
ALLOCATE(UAB(IMSUB,JMSUB))
ALLOCATE(VAF(IMSUB,JMSUB))
ALLOCATE(VA(IMSUB,JMSUB))
ALLOCATE(VAB(IMSUB,JMSUB))
ALLOCATE(ELF(IMSUB,JMSUB))
ALLOCATE(EL(IMSUB,JMSUB))
ALLOCATE(ELB(IMSUB,JMSUB))
ALLOCATE(PSI(IMSUB,JMSUB))
ALLOCATE(ETF(IMSUB,JMSUB))
ALLOCATE(ET(IMSUB,JMSUB))
ALLOCATE(ETB(IMSUB,JMSUB))
ALLOCATE(FLUXUA(IMSUB,JMSUB))
ALLOCATE(FLUXVA(IMSUB,JMSUB))
ALLOCATE(EGF(IMSUB,JMSUB))
ALLOCATE(EGB(IMSUB,JMSUB))
ALLOCATE(Z_TIDE(IMSUB,JMSUB))
ALLOCATE(LON_RHO(IMSUB,JMSUB))
ALLOCATE(LAT_RHO(IMSUB,JMSUB))
ALLOCATE(PM(IMSUB,JMSUB))
ALLOCATE(PN(IMSUB,JMSUB))
ALLOCATE(LONU(IMSUB-1,JMSUB))
ALLOCATE(LATU(IMSUB-1,JMSUB))
ALLOCATE(LONV(IMSUB,JMSUB-1))
ALLOCATE(LATV(IMSUB,JMSUB-1))

ALLOCATE(DUM(IMSUB,JMSUB))
ALLOCATE(DVM(IMSUB,JMSUB))
ALLOCATE(FSM(IMSUB,JMSUB))
ALLOCATE(DUM0(IMSUB,JMSUB))
ALLOCATE(DVM0(IMSUB,JMSUB))
ALLOCATE(FSM0(IMSUB,JMSUB))
ALLOCATE(MBSTATE(IMSUB,JMSUB))
ALLOCATE(UTB(IMSUB,JMSUB))
ALLOCATE(VTB(IMSUB,JMSUB))
ALLOCATE(UTF(IMSUB,JMSUB))
ALLOCATE(VTF(IMSUB,JMSUB))
ALLOCATE(ADVUA(IMSUB,JMSUB))
ALLOCATE(ADVVA(IMSUB,JMSUB))
ALLOCATE(TSURF(IMSUB,JMSUB))
ALLOCATE(SSURF(IMSUB,JMSUB))
ALLOCATE(DRX2D(IMSUB,JMSUB))
ALLOCATE(DRY2D(IMSUB,JMSUB))
ALLOCATE(ADX2D(IMSUB,JMSUB))
ALLOCATE(ADY2D(IMSUB,JMSUB))
ALLOCATE(SWRAD(IMSUB,JMSUB))
!3D VAR
ALLOCATE(A(IMSUB,JMSUB,KB))
ALLOCATE(C(IMSUB,JMSUB,KB))
ALLOCATE(EE(IMSUB,JMSUB,KB))
ALLOCATE(GG(IMSUB,JMSUB,KB))
ALLOCATE(UF(IMSUB,JMSUB,KB))
ALLOCATE(VF(IMSUB,JMSUB,KB))
ALLOCATE(KM(IMSUB,JMSUB,KB))
ALLOCATE(KH(IMSUB,JMSUB,KB))
ALLOCATE(KQ(IMSUB,JMSUB,KB))
ALLOCATE(L(IMSUB,JMSUB,KB))
ALLOCATE(Q2(IMSUB,JMSUB,KB))
ALLOCATE(Q2B(IMSUB,JMSUB,KB))
ALLOCATE(AAM(IMSUB,JMSUB,KB))
ALLOCATE(Q2L(IMSUB,JMSUB,KB))
ALLOCATE(Q2LB(IMSUB,JMSUB,KB))
ALLOCATE(U(IMSUB,JMSUB,KB))
ALLOCATE(UB(IMSUB,JMSUB,KB))
ALLOCATE(W(IMSUB,JMSUB,KB))
ALLOCATE(V(IMSUB,JMSUB,KB))
ALLOCATE(VB(IMSUB,JMSUB,KB))
ALLOCATE(T(IMSUB,JMSUB,KB))
ALLOCATE(TB(IMSUB,JMSUB,KB))
ALLOCATE(S(IMSUB,JMSUB,KB))
ALLOCATE(SB(IMSUB,JMSUB,KB))
ALLOCATE(RHO(IMSUB,JMSUB,KB))
ALLOCATE(DTEF(IMSUB,JMSUB,KB))
ALLOCATE(RMEAN(IMSUB,JMSUB,KB))

ALLOCATE(ADVX(IMSUB,JMSUB,KB))
ALLOCATE(ADVY(IMSUB,JMSUB,KB))
ALLOCATE(DRHOX(IMSUB,JMSUB,KB))
ALLOCATE(DRHOY(IMSUB,JMSUB,KB))
ALLOCATE(TCLIM(IMSUB,JMSUB,KB))
ALLOCATE(SCLIM(IMSUB,JMSUB,KB))

ALLOCATE(PCO2O(IMSUB,JMSUB))
ALLOCATE(DEPTH(IMSUB,JMSUB,KB))
ALLOCATE(TEMPER(IMSUB,JMSUB,KB))
ALLOCATE(SALT(IMSUB,JMSUB,KB))
ALLOCATE(LIGHT(IMSUB,JMSUB,KB))

ALLOCATE(SURF_RAD(IMSUB,JMSUB))
ALLOCATE(WIND10U(IMSUB,JMSUB))
ALLOCATE(WIND10V(IMSUB,JMSUB))
ALLOCATE(PCO2A(IMSUB,JMSUB))
ALLOCATE(CO2EX(IMSUB,JMSUB))

ALLOCATE(TYP(IMSUB,JMSUB))
ALLOCATE(NTYP(IMSUB,JMSUB))
ALLOCATE(ITYP(IMSUB*JMSUB))
ALLOCATE(JTYP(IMSUB*JMSUB))

ALLOCATE(CONB(IMSUB,JMSUB))
ALLOCATE(CON(IMSUB,JMSUB))
ALLOCATE(CONF(IMSUB,JMSUB))
ALLOCATE(COND(IMSUB,JMSUB))
ALLOCATE(CON_ERO(IMSUB,JMSUB))
ALLOCATE(CON_DEP(IMSUB,JMSUB))

ALLOCATE(ELR(IMSUB,JMSUB))
ALLOCATE(ELFR(IMSUB,JMSUB))
ALLOCATE(FSMR(IMSUB,JMSUB))

IMM1=IMSUB-1
JMM1=JMSUB-1
KBM1=KB-1
IMM2=IMSUB-2
JMM2=JMSUB-2
KBM2=KB-2
LIJ=IMSUB*JMSUB
LIJ2=LIJ*2
LIJK=LIJ*KB
LIJKM2=LIJ*KBM2
LIJM1=IMSUB*(JMSUB-1)
LIJKM1=LIJ*KBM1
IOCCE=IMSUB
JOCCE=JMSUB
IMOCC=IOCCE-IOCCB+1
JMOCC=JOCCE-JOCCB+1
KBOCC=KOCCE-KOCCB+1
END SUBROUTINE MYALLOCATE_SUB
END MODULE MOD_ALLOCATE