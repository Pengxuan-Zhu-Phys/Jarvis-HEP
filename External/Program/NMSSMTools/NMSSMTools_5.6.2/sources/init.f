      SUBROUTINE INITIALIZE()
      
*******************************************************************
*   This subroutine serves to
*     a) set default values for the SM parameters
*     b) set default values for the parameters (as limits on
*      sparticle masses) used for actual experimental constraints
*     c) read the tables that are used for actual experimental
*      constraints from the directory EXPCON
*******************************************************************

      IMPLICIT NONE

      CHARACTER*256 FILENAME,EXPCON_PATH,catpath
      CHARACTER dum*13

      INTEGER NBIN,I,J
      INTEGER NhZind,NhZbb,NhZll,NhZinv,NhZjj,NhZgg
      INTEGER NhA4b,NhA4tau,NhA2b2tau,NhA2tau2b
      INTEGER NAAA6b,NAAA6tau,NAAZ4b,NAAZ4tau,NAAZ2b2tau
      INTEGER Ncccc02,Ncccc04,Ncccc05,Ncccc06,Ncccc08,Ncccc1
      INTEGER Nccgg02,Nccgg04,Nccgg05,Nccgg06,Nccgg08,Nccgg1
      INTEGER Ncctt02,Ncctt04,Ncctt05,Ncctt06,Ncctt08,Ncctt1
      INTEGER Ngggg02,Ngggg04,Ngggg05,Ngggg06,Ngggg08,Ngggg1
      INTEGER Nttgg02,Nttgg04,Nttgg05,Nttgg06,Nttgg08,Nttgg1
      INTEGER Ntttt02,Ntttt04,Ntttt05,Ntttt06,Ntttt08,Ntttt1
      INTEGER Nstblsn,Nstnc,Nsbnb,Nglsq,NHGG1,NHGG2,NHGG3
      INTEGER NHAATAUS0,NHAATAUS1,NHAATAUS2,NHAATAUS3,NHAAMUS0
      INTEGER NHAAMUS1,NHAAMUS2,NHAAMUS3,NHAAMUS4,NHAAMUS5
      INTEGER NHAAMUTAU1,NHAAMUTAU2,NHAAMUB1,NHAAMUB2
      INTEGER NHAATAUB,NHAAGJ,NHAABS,NHAABS2,NHAAMUTAU3
      INTEGER NHAAMUB3,NHAAMUB4,NHAAMUS6
      INTEGER NHZH1,NHZH2,NHZH3,NHZH4,NHHH,NHHH2
      INTEGER NHZHT1,NHZHT2,NHZHT3,NHZHT4,NHHHT,NHHHT2
      INTEGER HZHT1(10000,3),HZHT2(10000,3),HZHT3(10000,3)
      INTEGER HZHT4(10000,3),HHHT(10000,3),HHHT2(10000,3)
      INTEGER IChargMAX,NCMS,IX,IY,NX,NY
      PARAMETER (NX=120,NY=80)

      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,ALEM0,MPI,MEL,MSTRANGE
      DOUBLE PRECISION GZMAX,GZINVMAX,MCHAMIN,MCMIN,SIGNEU1,SIGNEU
      DOUBLE PRECISION MSLMIN,MSTMIN,MSQMIN,MGLMIN
      DOUBLE PRECISION hZind(1000,2),hZbb(1000,2),hZll(1000,2)
      DOUBLE PRECISION hZinv(1000,2),hZjj(1000,2),hZgg(1000,2)
      DOUBLE PRECISION hA4b(10000,3),hA4tau(10000,3)
      DOUBLE PRECISION hA2b2tau(10000,3),hA2tau2b(10000,3)
      DOUBLE PRECISION AAA6b(10000,3),AAA6tau(10000,3)
      DOUBLE PRECISION AAZ4b(10000,3),AAZ4tau(10000,3)
      DOUBLE PRECISION AAZ2b2tau(10000,3)
      DOUBLE PRECISION cccc02(100,2),cccc04(100,2),cccc05(100,2)
      DOUBLE PRECISION cccc06(100,2),cccc08(100,2),cccc1(100,2)
      DOUBLE PRECISION ccgg02(100,2),ccgg04(100,2),ccgg05(100,2)
      DOUBLE PRECISION ccgg06(100,2),ccgg08(100,2),ccgg1(100,2)
      DOUBLE PRECISION cctt02(100,2),cctt04(100,2),cctt05(100,2)
      DOUBLE PRECISION cctt06(100,2),cctt08(100,2),cctt1(100,2)
      DOUBLE PRECISION gggg02(100,2),gggg04(100,2),gggg05(100,2)
      DOUBLE PRECISION gggg06(100,2),gggg08(100,2),gggg1(100,2)
      DOUBLE PRECISION ttgg02(100,2),ttgg04(100,2),ttgg05(100,2)
      DOUBLE PRECISION ttgg06(100,2),ttgg08(100,2),ttgg1(100,2)
      DOUBLE PRECISION tttt02(100,2),tttt04(100,2),tttt05(100,2)
      DOUBLE PRECISION tttt06(100,2),tttt08(100,2),tttt1(100,2)
      DOUBLE PRECISION stblsn(100,2),stnc(100,2),sbnb(100,2)
      DOUBLE PRECISION glsq(100,2),HGG1(300,2),HGG2(300,2),HGG3(300,2)
      DOUBLE PRECISION HAATAUS0(100,2),HAATAUS1(300,2),HAATAUS2(100,2)
      DOUBLE PRECISION HAATAUS3(100,2),HAAMUS0(100,2),HAAMUS1(100,2)
      DOUBLE PRECISION HAAMUS2(100,2),HAAMUS3(100,2),HAAMUS4(100,4)
      DOUBLE PRECISION HAAMUS5(100,4),HAAMUTAU1(100,2),HAAMUTAU2(100,2)
      DOUBLE PRECISION HAAMUB1(100,2),HAAMUB2(200,2),HAATAUB(100,2)
      DOUBLE PRECISION HAAGJ(100,2),HAABS(100,2),HAABS2(100,2)
      DOUBLE PRECISION HAAMUTAU3(200,2),HAAMUB3(100,2),HAAMUB4(100,2)
      DOUBLE PRECISION HAAMUS6(200,2),HZH1(10000,3),HZH2(10000,3)
      DOUBLE PRECISION HZH3(10000,3),HZH4(10000,3),HHH(10000,3)
      DOUBLE PRECISION HHH2(10000,3)
      DOUBLE PRECISION Xf,sigmaV,vcsll,vcsbb,x(100),dNdx(100),EMIN
      DOUBLE PRECISION sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      DOUBLE PRECISION OMG,OMGMIN,OMGMAX,muH2
      DOUBLE PRECISION MHmin,MHmax
      DOUBLE PRECISION LCMS,HCMS(NX,NY),XMIN,XMAX,YMIN,YMAX
      DOUBLE PRECISION DmNeutmCharg(20),mCharg(50),prosp(17,50)

      COMMON/HIGGSMS/muH2
      COMMON/HIGGSFIT/MHmin,MHmax
      COMMON/ALEM0/ALEM0
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SMEXT/MPI,MEL,MSTRANGE
      COMMON/CKM/VUS,VCB,VUB
      COMMON/LEP/GZMAX,GZINVMAX,MCHAMIN,MCMIN,SIGNEU1,SIGNEU,
     .      MSLMIN,MSTMIN,MSQMIN,MGLMIN,
     .      hZind,hZbb,hZll,hZinv,hZjj,hZgg,
     .      hA4b,hA4tau,hA2b2tau,hA2tau2b,
     .      AAA6b,AAA6tau,AAZ4b,AAZ4tau,AAZ2b2tau,
     .      cccc02,cccc04,cccc05,cccc06,cccc08,cccc1,
     .      ccgg02,ccgg04,ccgg05,ccgg06,ccgg08,ccgg1,
     .      cctt02,cctt04,cctt05,cctt06,cctt08,cctt1,
     .      gggg02,gggg04,gggg05,gggg06,gggg08,gggg1,
     .      ttgg02,ttgg04,ttgg05,ttgg06,ttgg08,ttgg1,
     .      tttt02,tttt04,tttt05,tttt06,tttt08,tttt1,
     .      stblsn,stnc,sbnb,glsq,
     .      NhZind,NhZbb,NhZll,NhZinv,NhZjj,NhZgg,
     .      NhA4b,NhA4tau,NhA2b2tau,NhA2tau2b,
     .      NAAA6b,NAAA6tau,NAAZ4b,NAAZ4tau,NAAZ2b2tau,
     .      Ncccc02,Ncccc04,Ncccc05,Ncccc06,Ncccc08,Ncccc1,
     .      Nccgg02,Nccgg04,Nccgg05,Nccgg06,Nccgg08,Nccgg1,
     .      Ncctt02,Ncctt04,Ncctt05,Ncctt06,Ncctt08,Ncctt1,
     .      Ngggg02,Ngggg04,Ngggg05,Ngggg06,Ngggg08,Ngggg1,
     .      Nttgg02,Nttgg04,Nttgg05,Nttgg06,Nttgg08,Nttgg1,
     .      Ntttt02,Ntttt04,Ntttt05,Ntttt06,Ntttt08,Ntttt1,
     .      Nstblsn,Nstnc,Nsbnb,Nglsq
      COMMON/LHCHAA/HAATAUS0,HAATAUS1,HAATAUS2,HAATAUS3,
     .      HAAMUS0,HAAMUS1,HAAMUS2,HAAMUS3,HAAMUS4,HAAMUS5,
     .      HAAMUTAU1,HAAMUTAU2,HAAMUB1,HAAMUB2,HAATAUB,HAAGJ,
     .      HAABS,HAABS2,HAAMUTAU3,HAAMUB3,HAAMUB4,HAAMUS6,
     .      NHAATAUS0,NHAATAUS1,NHAATAUS2,NHAATAUS3,NHAAMUS0,
     .      NHAAMUS1,NHAAMUS2,NHAAMUS3,NHAAMUS4,NHAAMUS5,
     .      NHAAMUTAU1,NHAAMUTAU2,NHAAMUB1,NHAAMUB2,
     .      NHAATAUB,NHAAGJ,NHAABS,NHAABS2,NHAAMUTAU3,
     .      NHAAMUB3,NHAAMUB4,NHAAMUS6
      COMMON/LHCHGG/HGG1,HGG2,HGG3,NHGG1,NHGG2,NHGG3
      COMMON/LHZH/HZH1,HZH2,HZH3,HZH4,HZHT1,HZHT2,HZHT3,HZHT4,
     .      NHZH1,NHZH2,NHZH3,NHZH4,NHZHT1,NHZHT2,NHZHT3,NHZHT4
      COMMON/LHHH/HHH,HHH2,HHHT,HHHT2,NHHH,NHHH2,NHHHT,NHHHT2
      COMMON/MICROMG/OMG,OMGMIN,OMGMAX,Xf,sigmaV,vcsll,vcsbb,
     .      x,dNdx,EMIN,NBIN
      COMMON/MICROMG2/sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      COMMON/CMSINOS/HCMS,XMIN,XMAX,YMIN,YMAX,DmNeutmCharg,
     .      mCharg,prosp,IChargMAX

* SM inputs:

*   Higgs mass squared
      muH2=(125.1d0)**2
*   +/- 3 GeV theor. errors
      MHmin=DSQRT(muH2)-3d0
      MHmax=DSQRT(muH2)+3d0
*   Alpha_EM(0)
      ALEM0= 1d0/137.036d0
*   Alpha_s at MZ:
      ALSMZ= 0.1181d0
*   Electroweak parameters:
      GF= 1.1663787d-5
*   Alpha_em at MZ (used in RGES):
      ALEMMZ= 1d0/127.92d0
*   Z, W pole masses:
      MZ= 91.187d0
      MW= 80.42d0
*   Lepton masses:            
      MTAU= 1.777d0
      MMUON= 105.6583715d-3
*   Quark pole masses:
      MS= 0.19d0
      MC= 1.4d0
      MBP= 4.54d0
      MT= 173.4d0
*   Running MS_bar bottom mass MB at the scale MB:
      MB= 4.18d0      
*   Elements of the Kobayashi-Maskawa matrix:
      VUS= 0.22d0
      VCB= 0.04d0
      VUB= 0.004d0
*   Pion/electron masses
      MPI=135d-3
      MEL=510.998928d-6
*   Strange quark running mass:
      MSTRANGE=95D-3

* Dark matter constraints
* (Used only if OMGFLAG=/=0)

      OMGMIN=0.1187d0*0.9d0   ! Planck
      OMGMAX=0.1187d0*1.1d0   ! +/-10%
      sigmaPiN=34d0!36d0(p),39d0(n),37.5d0(N)
      sigmaS=42d0!54d0
      NBIN=10
      EMIN=1d-3
      DO I=1,100
       x(I)=0d0
       dNdx(I)=0d0
      ENDDO
      DO I=1,NBIN
        x(I)=(NBIN-I)*DLOG10(Emin)/(NBIN-1d0)
      ENDDO

*   Collider constraints on sparticles:

*   Limit on the Z width for Z -> h(i) + a(j):
      GZMAX=5.78d-3
*   Limit on the inv. Z width from Z -> neutralinos:      
      GZINVMAX=0.5d-3
*   Limit on sigma(e+e- -> neutralinos (1,i)) (i>1):      
      SIGNEU1=1d-2
*   Limit on sigma(e+e- -> neutralinos (i,j)) (i,j > 1):            
      SIGNEU=1d-1
*   Lower limit on chargino masses:
      MCHAMIN=103.5d0
*   Lower limit on slepton masses:
      MSLMIN=99.9d0
*   Lower limit on stau masses:      
      MSTMIN=93.2d0
*   Lower limit on squark masses:
      MSQMIN=100d0
*   Lower limit on gluino mass:
      MGLMIN=180d0
*   Lower limit on charged Higgs mass:
      MCMIN=78.6d0

*   LEP constraints

      CALL getenv('EXPCON_PATH',EXPCON_PATH)
      if(EXPCON_PATH.eq.' ')  EXPCON_PATH='../EXPCON'

      FILENAME=catpath(EXPCON_PATH,'hZind.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 111  READ(11,*,END=112,ERR=2)(hZind(I,J),J=1,2)
      I=I+1
      GOTO 111
 112  NhZind=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hZbb.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 211  READ(11,*,END=212,ERR=2)(hZbb(I,J),J=1,2)
      I=I+1
      GOTO 211
 212  NhZbb=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hZll.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 311  READ(11,*,END=312,ERR=2)(hZll(I,J),J=1,2)
      I=I+1
      GOTO 311
 312  NhZll=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hZinv.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 411  READ(11,*,END=412,ERR=2)(hZinv(I,J),J=1,2)
      I=I+1
      GOTO 411
 412  NhZinv=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hZjj.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 511  READ(11,*,END=512,ERR=2)(hZjj(I,J),J=1,2)
      I=I+1
      GOTO 511
 512  NhZjj=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hZgg.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 611  READ(11,*,END=612,ERR=2)(hZgg(I,J),J=1,2)
      I=I+1
      GOTO 611
 612  NhZgg=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hA4b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 711  READ(11,*,END=712,ERR=2)(hA4b(I,J),J=1,3)
      I=I+1
      GOTO 711
 712  NhA4b=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hA4tau.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 811  READ(11,*,END=812,ERR=2)(hA4tau(I,J),J=1,3)
      I=I+1
      GOTO 811
 812  NhA4tau=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hA2b2tau.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 821  READ(11,*,END=822,ERR=2)(hA2b2tau(I,J),J=1,3)
      I=I+1
      GOTO 821
 822  NhA2b2tau=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'hA2tau2b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 831  READ(11,*,END=832,ERR=2)(hA2tau2b(I,J),J=1,3)
      I=I+1
      GOTO 831
 832  NhA2tau2b=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'AAA6b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 911  READ(11,*,END=912,ERR=2)(AAA6b(I,J),J=1,3)
      I=I+1
      GOTO 911
 912  NAAA6b=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'AAA6tau.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 921  READ(11,*,END=922,ERR=2)(AAA6tau(I,J),J=1,3)
      I=I+1
      GOTO 921
 922  NAAA6tau=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'AAZ4b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1101 READ(11,*,END=1102,ERR=2)(AAZ4b(I,J),J=1,3)
      I=I+1
      GOTO 1101
 1102 NAAZ4b=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'AAZ4tau.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1111 READ(11,*,END=1112,ERR=2)(AAZ4tau(I,J),J=1,3)
      I=I+1
      GOTO 1111
 1112 NAAZ4tau=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'AAZ2b2tau.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1121 READ(11,*,END=1122,ERR=2)(AAZ2b2tau(I,J),J=1,3)
      I=I+1
      GOTO 1121
 1122 NAAZ2b2tau=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cccc02.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1201 READ(11,*,END=1202)cccc02(I,1)
      READ(11,*,END=1202,ERR=2)cccc02(I,2)
      I=I+1
      GOTO 1201
 1202 Ncccc02=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cccc04.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1211 READ(11,*,END=1212)cccc04(I,1)
      READ(11,*,END=1212,ERR=2)cccc04(I,2)
      I=I+1
      GOTO 1211
 1212 Ncccc04=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cccc05.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1221 READ(11,*,END=1222)cccc05(I,1)
      READ(11,*,END=1222,ERR=2)cccc05(I,2)
      I=I+1
      GOTO 1221
 1222 Ncccc05=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cccc06.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1231 READ(11,*,END=1232)cccc06(I,1)
      READ(11,*,END=1232,ERR=2)cccc06(I,2)
      I=I+1
      GOTO 1231
 1232 Ncccc06=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cccc08.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1241 READ(11,*,END=1242)cccc08(I,1)
      READ(11,*,END=1242,ERR=2)cccc08(I,2)
      I=I+1
      GOTO 1241
 1242 Ncccc08=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cccc1.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1251 READ(11,*,END=1252)cccc1(I,1)
      READ(11,*,END=1252,ERR=2)cccc1(I,2)
      I=I+1
      GOTO 1251
 1252 Ncccc1=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ccgg02.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1301 READ(11,*,END=1302)ccgg02(I,1)
      READ(11,*,END=1302,ERR=2)ccgg02(I,2)
      I=I+1
      GOTO 1301
 1302 Nccgg02=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ccgg04.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1311 READ(11,*,END=1312)ccgg04(I,1)
      READ(11,*,END=1312,ERR=2)ccgg04(I,2)
      I=I+1
      GOTO 1311
 1312 Nccgg04=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ccgg05.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1321 READ(11,*,END=1322)ccgg05(I,1)
      READ(11,*,END=1322,ERR=2)ccgg05(I,2)
      I=I+1
      GOTO 1321
 1322 Nccgg05=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ccgg06.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1331 READ(11,*,END=1332)ccgg06(I,1)
      READ(11,*,END=1332,ERR=2)ccgg06(I,2)
      I=I+1
      GOTO 1331
 1332 Nccgg06=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ccgg08.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1341 READ(11,*,END=1342)ccgg08(I,1)
      READ(11,*,END=1342,ERR=2)ccgg08(I,2)
      I=I+1
      GOTO 1341
 1342 Nccgg08=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ccgg1.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1351 READ(11,*,END=1352)ccgg1(I,1)
      READ(11,*,END=1352,ERR=2)ccgg1(I,2)
      I=I+1
      GOTO 1351
 1352 Nccgg1=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cctt02.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1401 READ(11,*,END=1402)cctt02(I,1)
      READ(11,*,END=1402,ERR=2)cctt02(I,2)
      I=I+1
      GOTO 1401
 1402 Ncctt02=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cctt04.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1411 READ(11,*,END=1412)cctt04(I,1)
      READ(11,*,END=1412,ERR=2)cctt04(I,2)
      I=I+1
      GOTO 1411
 1412 Ncctt04=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cctt05.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1421 READ(11,*,END=1422)cctt05(I,1)
      READ(11,*,END=1422,ERR=2)cctt05(I,2)
      I=I+1
      GOTO 1421
 1422 Ncctt05=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cctt06.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1431 READ(11,*,END=1432)cctt06(I,1)
      READ(11,*,END=1432,ERR=2)cctt06(I,2)
      I=I+1
      GOTO 1431
 1432 Ncctt06=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cctt08.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1441 READ(11,*,END=1442)cctt08(I,1)
      READ(11,*,END=1442,ERR=2)cctt08(I,2)
      I=I+1
      GOTO 1441
 1442 Ncctt08=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'cctt1.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1451 READ(11,*,END=1452)cctt1(I,1)
      READ(11,*,END=1452,ERR=2)cctt1(I,2)
      I=I+1
      GOTO 1451
 1452 Ncctt1=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'gggg02.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1501 READ(11,*,END=1502)gggg02(I,1)
      READ(11,*,END=1502,ERR=2)gggg02(I,2)
      I=I+1
      GOTO 1501
 1502 Ngggg02=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'gggg04.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1511 READ(11,*,END=1512)gggg04(I,1)
      READ(11,*,END=1512,ERR=2)gggg04(I,2)
      I=I+1
      GOTO 1511
 1512 Ngggg04=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'gggg05.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1521 READ(11,*,END=1522)gggg05(I,1)
      READ(11,*,END=1522,ERR=2)gggg05(I,2)
      I=I+1
      GOTO 1521
 1522 Ngggg05=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'gggg06.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1531 READ(11,*,END=1532)gggg06(I,1)
      READ(11,*,END=1532,ERR=2)gggg06(I,2)
      I=I+1
      GOTO 1531
 1532 Ngggg06=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'gggg08.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1541 READ(11,*,END=1542)gggg08(I,1)
      READ(11,*,END=1542,ERR=2)gggg08(I,2)
      I=I+1
      GOTO 1541
 1542 Ngggg08=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'gggg1.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1551 READ(11,*,END=1552)gggg1(I,1)
      READ(11,*,END=1552,ERR=2)gggg1(I,2)
      I=I+1
      GOTO 1551
 1552 Ngggg1=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ttgg02.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1601 READ(11,*,END=1602)ttgg02(I,1)
      READ(11,*,END=1602,ERR=2)ttgg02(I,2)
      I=I+1
      GOTO 1601
 1602 Nttgg02=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ttgg04.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1611 READ(11,*,END=1612)ttgg04(I,1)
      READ(11,*,END=1612,ERR=2)ttgg04(I,2)
      I=I+1
      GOTO 1611
 1612 Nttgg04=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ttgg05.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1621 READ(11,*,END=1622)ttgg05(I,1)
      READ(11,*,END=1622,ERR=2)ttgg05(I,2)
      I=I+1
      GOTO 1621
 1622 Nttgg05=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ttgg06.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1631 READ(11,*,END=1632)ttgg06(I,1)
      READ(11,*,END=1632,ERR=2)ttgg06(I,2)
      I=I+1
      GOTO 1631
 1632 Nttgg06=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ttgg08.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1641 READ(11,*,END=1642)ttgg08(I,1)
      READ(11,*,END=1642,ERR=2)ttgg08(I,2)
      I=I+1
      GOTO 1641
 1642 Nttgg08=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'ttgg1.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1651 READ(11,*,END=1652)ttgg1(I,1)
      READ(11,*,END=1652,ERR=2)ttgg1(I,2)
      I=I+1
      GOTO 1651
 1652 Nttgg1=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'tttt02.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1701 READ(11,*,END=1702)tttt02(I,1)
      READ(11,*,END=1702,ERR=2)tttt02(I,2)
      I=I+1
      GOTO 1701
 1702 Ntttt02=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'tttt04.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1711 READ(11,*,END=1712)tttt04(I,1)
      READ(11,*,END=1712,ERR=2)tttt04(I,2)
      I=I+1
      GOTO 1711
 1712 Ntttt04=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'tttt05.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1721 READ(11,*,END=1722)tttt05(I,1)
      READ(11,*,END=1722,ERR=2)tttt05(I,2)
      I=I+1
      GOTO 1721
 1722 Ntttt05=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'tttt06.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1731 READ(11,*,END=1732)tttt06(I,1)
      READ(11,*,END=1732,ERR=2)tttt06(I,2)
      I=I+1
      GOTO 1731
 1732 Ntttt06=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'tttt08.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1741 READ(11,*,END=1742)tttt08(I,1)
      READ(11,*,END=1742,ERR=2)tttt08(I,2)
      I=I+1
      GOTO 1741
 1742 Ntttt08=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'tttt1.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1751 READ(11,*,END=1752)tttt1(I,1)
      READ(11,*,END=1752,ERR=2)tttt1(I,2)
      I=I+1
      GOTO 1751
 1752 Ntttt1=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'stblsn.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1761 READ(11,*,END=1762,ERR=2)(stblsn(I,J),J=1,2)
      I=I+1
      GOTO 1761
 1762 Nstblsn=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'stnc.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1771 READ(11,*,END=1772,ERR=2)(stnc(I,J),J=1,2)
      I=I+1
      GOTO 1771
 1772 Nstnc=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'sbnb.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1781 READ(11,*,END=1782,ERR=2)(sbnb(I,J),J=1,2)
      I=I+1
      GOTO 1781
 1782 Nsbnb=I-1
      CLOSE(11)

      FILENAME=catpath(EXPCON_PATH,'glsq.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1791 READ(11,*,END=1792,ERR=2)(glsq(I,J),J=1,2)
      I=I+1
      GOTO 1791
 1792 Nglsq=I-1
      CLOSE(11)

* Read CMS 8TeV upper limit
* HGG1(I,1): Higgs mass
* HGG1(I,2): upper on sigma[H_125]/sigma_SM*BR[H->2gamma] limit in pb
      FILENAME=catpath(EXPCON_PATH,'cms-pas-hig-17-013_5a.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1801 READ(11,*,END=1802,ERR=2)(HGG1(I,J),J=1,2)
      I=I+1
      GOTO 1801
 1802 NHGG1=I-1
      CLOSE(11)

* Read CMS 13TeV upper limit
* HGG2(I,1): Higgs mass
* HGG2(I,2): upper on sigma[H_125]/sigma_SM*BR[H->2gamma] limit in pb
      FILENAME=catpath(EXPCON_PATH,'cms-pas-hig-17-013_5c.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1803 READ(11,*,END=1804,ERR=2)(HGG2(I,J),J=1,2)
      I=I+1
      GOTO 1803
 1804 NHGG2=I-1
      CLOSE(11)

* Read ATLAS 13TeV upper limit
* HGG3(I,1): Higgs mass
* HGG3(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2gamma] in pb
      FILENAME=catpath(EXPCON_PATH,'ATLAS-CONF-2018-025_4b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1805 READ(11,*,END=1806,ERR=2)(HGG3(I,J),J=1,2)
      HGG3(I,2)=HGG3(I,2)/1d3/(.57d0-2.07d0*DEXP(-3.13d0*HGG3(I,1)/1d2))
      I=I+1
      GOTO 1805
 1806 NHGG3=I-1
      CLOSE(11)

* Read CMS upper limit
* HAAMUTAU1(I,1): CP-odd Higgs mass
* HAAMUTAU1(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2tau]*BR[A->2mu]
      FILENAME=catpath(EXPCON_PATH,'1701.02032_6f.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1807 READ(11,*,END=1808,ERR=2)(HAAMUTAU1(I,J),J=1,2)
      I=I+1
      GOTO 1807
 1808 NHAAMUTAU1=I-1
      CLOSE(11)

* Read CMS upper limit
* HAAMUTAU2(I,1): CP-odd Higgs mass
* HAAMUTAU2(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2tau]*BR[A->2mu]
      FILENAME=catpath(EXPCON_PATH,'1805.04865_4e.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1809 READ(11,*,END=1810,ERR=2)(HAAMUTAU2(I,J),J=1,2)
      I=I+1
      GOTO 1809
 1810 NHAAMUTAU2=I-1
      CLOSE(11)

* Read CMS upper limit
* HAAMUTAU3(I,1): CP-odd Higgs mass
* HAAMUTAU3(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2tau]*BR[A->2mu]
      FILENAME=catpath(EXPCON_PATH,'2005.08694_7a.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1823 READ(11,*,END=1824,ERR=2)(HAAMUTAU3(I,J),J=1,2)
      I=I+1
      GOTO 1823
 1824 NHAAMUTAU3=I-1
      CLOSE(11)

* Read CMS upper limit
* HAAMUB1(I,1): CP-odd Higgs mass
* HAAMUB1(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2mu]*BR[A->2b]
      FILENAME=catpath(EXPCON_PATH,'cms-pas-hig-14-041_4b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1811 READ(11,*,END=1812,ERR=2)(HAAMUB1(I,J),J=1,2)
      I=I+1
      GOTO 1811
 1812 NHAAMUB1=I-1
      CLOSE(11)

* Read ATLAS upper limit
* HAAMUB2(I,1): CP-odd Higgs mass
* HAAMUB2(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2mu]*BR[A->2b]
      FILENAME=catpath(EXPCON_PATH,'1807.00539_6a.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1813 READ(11,*,END=1814,ERR=2)(HAAMUB2(I,J),J=1,2)
      I=I+1
      GOTO 1813
 1814 NHAAMUB2=I-1
      CLOSE(11)

* Read ATLAS upper limit
* HAAMUB3(I,1): CP-odd Higgs mass
* HAAMUB3(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2mu]*BR[A->2b]
      FILENAME=catpath(EXPCON_PATH,'2110.00313_9.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1825 READ(11,*,END=1826,ERR=2)(HAAMUB3(I,J),J=1,2)
      I=I+1
      GOTO 1825
 1826 NHAAMUB3=I-1
      CLOSE(11)

* Read CMS upper limit
* HAAMUB4(I,1): CP-odd Higgs mass
* HAAMUB4(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2mu]*BR[A->2b]
      FILENAME=catpath(EXPCON_PATH,'1812.06359_5b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1827 READ(11,*,END=1828,ERR=2)(HAAMUB4(I,J),J=1,2)
      I=I+1
      GOTO 1827
 1828 NHAAMUB4=I-1
      CLOSE(11)

* Read CMS upper limit
* HAATAUB(I,1): CP-odd Higgs mass
* HAATAUB(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2tau]*BR[A->2b]
      FILENAME=catpath(EXPCON_PATH,'1805.10191_7d.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1815 READ(11,*,END=1816,ERR=2)(HAATAUB(I,J),J=1,2)
      I=I+1
      GOTO 1815
 1816 NHAATAUB=I-1
      CLOSE(11)

* Read ATLAS upper limit
* HAAGJ(I,1): CP-odd Higgs mass
* HAAGJ(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2j]*BR[A->2g]
      FILENAME=catpath(EXPCON_PATH,'1803.11145_2.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1817 READ(11,*,END=1818,ERR=2)(HAAGJ(I,J),J=1,2)
      I=I+1
      GOTO 1817
 1818 NHAAGJ=I-1
      CLOSE(11)

* Read ATLAS upper limit
* HAABS(I,1): CP-odd Higgs mass
* HAABS(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2b]**2
      FILENAME=catpath(EXPCON_PATH,'1806.07355_9c.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1819 READ(11,*,END=1820,ERR=2)(HAABS(I,J),J=1,2)
      HAABS(I,2)=HAABS(I,2)/2.238d0
      I=I+1
      GOTO 1819
 1820 NHAABS=I-1
      CLOSE(11)

* Read ATLAS upper limit
* HAABS2(I,1): CP-odd Higgs mass
* HAABS2(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2b]**2
      FILENAME=catpath(EXPCON_PATH,'2005.12236_10.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1821 READ(11,*,END=1822,ERR=2)(HAABS2(I,J),J=1,2)
      HAABS2(I,2)=HAABS2(I,2)/.88d0
      I=I+1
      GOTO 1821
 1822 NHAABS2=I-1
      CLOSE(11)

* Read CMS upper limit
* HAATAUS0(I,1): CP-odd Higgs mass
* HAATAUS0(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2tau]**2
      FILENAME=catpath(EXPCON_PATH,'cms-pas-hig-15-011_5f.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1899 READ(11,*,END=1900,ERR=2)(HAATAUS0(I,J),J=1,2)
      I=I+1
      GOTO 1899
 1900 NHAATAUS0=I-1
      CLOSE(11)

* Read ATLAS upper limit
* HAATAUS1(I,1): CP-odd Higgs mass
* HAATAUS1(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2tau]**2
      FILENAME=catpath(EXPCON_PATH,'1505.01609_6a.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1901 READ(11,*,END=1902,ERR=2)(HAATAUS1(I,J),J=1,2)
      I=I+1
      GOTO 1901
 1902 NHAATAUS1=I-1
      CLOSE(11)

* Read CMS upper limit
* HAATAUS2(I,1): CP-odd Higgs mass
* HAATAUS2(I,2): upper limit on sigma[H_125]*BR[H->2A]*BR[A->2tau]**2
      FILENAME=catpath(EXPCON_PATH,'1510.06534_8.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1903 READ(11,*,END=1904,ERR=2)(HAATAUS2(I,J),J=1,2)
      I=I+1
      GOTO 1903
 1904 NHAATAUS2=I-1
      CLOSE(11)

* Read CMS upper limit
* HAATAUS3(I,1): CP-odd Higgs mass
* HAATAUS3(I,2): upper limit on sigma[H_125]*BR[H->2A]*BR[A->2tau]**2
      FILENAME=catpath(EXPCON_PATH,'cms-pas-hig-14-022_7c.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1905 READ(11,*,END=1906,ERR=2)(HAATAUS3(I,J),J=1,2)
      I=I+1
      GOTO 1905
 1906 NHAATAUS3=I-1
      CLOSE(11)

* Read CMS upper limit
* HAAMUS0(I,1): CP-odd Higgs mass
* HAAMUS0(I,2): upper limit on BR[H->2A]*BR[A->2mu]**2
      FILENAME=catpath(EXPCON_PATH,'cms-pas-hig-13-010_6b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1907 READ(11,*,END=1908,ERR=2)(HAAMUS0(I,J),J=1,2)
      I=I+1
      GOTO 1907
 1908 NHAAMUS0=I-1
      CLOSE(11)

* Read CMS upper limit
* HAAMUS1(I,1): CP-odd Higgs mass
* HAAMUS1(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2mu]**2
      FILENAME=catpath(EXPCON_PATH,'1701.02032_7c.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1909 READ(11,*,END=1910,ERR=2)(HAAMUS1(I,J),J=1,2)
      I=I+1
      GOTO 1909
 1910 NHAAMUS1=I-1
      CLOSE(11)

* Read CMS upper limit
* HAAMUS2(I,1): CP-odd Higgs mass
* HAAMUS2(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2mu]**2
      FILENAME=catpath(EXPCON_PATH,'1701.02032_7b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1911 READ(11,*,END=1912,ERR=2)(HAAMUS2(I,J),J=1,2)
      I=I+1
      GOTO 1911
 1912 NHAAMUS2=I-1
      CLOSE(11)

* Read CMS upper limit
* HAAMUS3(I,1): CP-odd Higgs mass
* HAAMUS3(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2mu]**2
      FILENAME=catpath(EXPCON_PATH,'1701.02032_7a.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1913 READ(11,*,END=1914,ERR=2)(HAAMUS3(I,J),J=1,2)
      I=I+1
      GOTO 1913
 1914 NHAAMUS3=I-1
      CLOSE(11)

* Read CMS upper limit
* HAAMUS4(I,1): CP-even Higgs mass
* HAAMUS4(I,2): upper limit on sigma[H]*BR[H->2A]*BR[A->2mu]**2 in fb, M_A = 0.25GeV
* HAAMUS4(I,3): upper limit on sigma[H]*BR[H->2A]*BR[A->2mu]**2 in fb, M_A = 2GeV
* HAAMUS4(I,4): upper limit on sigma[H]*BR[H->2A]*BR[A->2mu]**2 in fb, M_A = 3.55GeV
      FILENAME=catpath(EXPCON_PATH,'1506.00424_2a.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1915 READ(11,*,END=1916,ERR=2)(HAAMUS4(I,J),J=1,4)
      I=I+1
      GOTO 1915
 1916 NHAAMUS4=I-1
      CLOSE(11)

* Read CMS upper limit
* HAAMUS5(I,1): CP-odd Higgs mass
* HAAMUS5(I,2): upper limit on sigma[H]*BR[H->2A]*BR[A->2mu]**2 in fb, M_H = 90GeV
* HAAMUS5(I,3): upper limit on sigma[H]*BR[H->2A]*BR[A->2mu]**2 in fb, M_H = 125GeV
* HAAMUS5(I,4): upper limit on sigma[H]*BR[H->2A]*BR[A->2mu]**2 in fb, M_H = 150GeV
      FILENAME=catpath(EXPCON_PATH,'cms-pas-hig-18-003_2b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1917 READ(11,*,END=1918,ERR=2)(HAAMUS5(I,J),J=1,4)
      I=I+1
      GOTO 1917
 1918 NHAAMUS5=I-1
      CLOSE(11)

* Read ATLAS upper limit
* HAAMUS6(I,1): CP-odd Higgs mass
* HAAMUS6(I,2): upper limit on sigma[H_125]/sigma_SM*BR[H->2A]*BR[A->2mu]**2
      FILENAME=catpath(EXPCON_PATH,'2110.13673_14b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1919 READ(11,*,END=1920,ERR=2)(HAAMUS6(I,J),J=1,2)
      HAAMUS6(I,2)=HAAMUS6(I,2)/48.61d0
      I=I+1
      GOTO 1919
 1920 NHAAMUS6=I-1
      CLOSE(11)

* Read CMS upper limit
* HZH1(I,1): CP-odd Higgs mass
* HZH1(I,2): CP-even Higgs mass
* HZH1(I,3): upper limit on sigma[H]*BR[H->ZA]*BR[A->2b] in pb
      FILENAME=catpath(EXPCON_PATH,'1911.03781_7b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1921 READ(11,*,END=1922,ERR=2)(HZH1(I,J),J=1,3)
      I=I+1
      GOTO 1921
 1922 NHZH1=I-1
      CLOSE(11)
      CALL TRITAB(HZH1,NHZH1,HZHT1,NHZHT1)

* Read ATLAS upper limit
* HZH2(I,1): CP-odd Higgs mass
* HZH2(I,2): CP-even Higgs mass
* HZH2(I,3): upper limit on sigma[H]*BR[H->ZA]*BR[A->2b] in pb
      FILENAME=catpath(EXPCON_PATH,'2011.05639_9b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1923 READ(11,*,END=1924,ERR=2)(HZH2(I,J),J=1,3)
      I=I+1
      GOTO 1923
 1924 NHZH2=I-1
      CLOSE(11)
      CALL TRITAB(HZH2,NHZH2,HZHT2,NHZHT2)

* Read ATLAS upper limit
* HZH3(I,1): CP-odd Higgs mass
* HZH3(I,2): CP-even Higgs mass
* HZH3(I,3): upper limit on sigma[H]*BR[H->ZA]*BR[A->2b] in pb
      FILENAME=catpath(EXPCON_PATH,'2011.05639_9d.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1925 READ(11,*,END=1926,ERR=2)(HZH3(I,J),J=1,3)
      I=I+1
      GOTO 1925
 1926 NHZH3=I-1
      CLOSE(11)
      CALL TRITAB(HZH3,NHZH3,HZHT3,NHZHT3)

* Read ATLAS upper limit
* HZH4(I,1): CP-odd Higgs mass
* HZH4(I,2): CP-even Higgs mass
* HZH4(I,3): upper limit on sigma[H]*BR[H->ZA]*BR[A->2b] in pb
      FILENAME=catpath(EXPCON_PATH,'2011.05639_13b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1927 READ(11,*,END=1928,ERR=2)(HZH4(I,J),J=1,3)
      I=I+1
      GOTO 1927
 1928 NHZH4=I-1
      CLOSE(11)
      CALL TRITAB(HZH4,NHZH4,HZHT4,NHZHT4)

* Read CMS upper limit
* HHH(I,1): Heavy Higgs mass
* HHH(I,2): Light Higgs mass
* HHH(I,3): upper limit on sigma[H]*BR[H->HSMHS]*BR[HSM->2tau]*BR[HS->2b] in pb
      FILENAME=catpath(EXPCON_PATH,'2106.10361_5.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1929 READ(11,*,END=1930,ERR=2)(HHH(I,J),J=1,3)
      I=I+1
      GOTO 1929
 1930 NHHH=I-1
      CLOSE(11)
      CALL TRITAB(HHH,NHHH,HHHT,NHHHT)

* Read CMS upper limit
* HHH(I,1): Heavy Higgs mass
* HHH(I,2): Light Higgs mass
* HHH(I,3): upper limit on sigma[H]*BR[H->HSMHS]*BR[HSM->2b]*BR[HS->2b] in pb
      FILENAME=catpath(EXPCON_PATH,'cms-pas-b2g-21-003_3b.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
      I=1
 1931 READ(11,*,END=1932,ERR=2)(HHH2(I,J),J=1,3)
      I=I+1
      GOTO 1931
 1932 NHHH2=I-1
      CLOSE(11)
      CALL TRITAB(HHH2,NHHH2,HHHT2,NHHHT2)

* CMS constraints on charginos/neutralinos
* CMS-PAS-SUS-17-004
      FILENAME=catpath(EXPCON_PATH,'higgsino-grid.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
* Read the available Delta(M_charg-M_neut) from the grid:
      READ(11,'(A13,17F13.0)',ERR=2) dum,(DmNeutmCharg(J),J=1,17)
!      WRITE(*,'(A19,17F5.0)') "DmNeutmCharg(1-17):",
!     .    (DmNeutmCharg(J),J=1,17)
* Read the available Mcharg and Prospino Xsect(1-17) from the grid:
      I=1
 2701 READ(11,'(F13.0,17F13.0)',END=2702,ERR=3) Mcharg(I),
     .    (prosp(J,I),J=1,17)
      I=I+1
      GOTO 2701
 2702 CLOSE(11)
      IchargMAX=I-1
!      WRITE(*,*) "IchargMAX, Mcharg(IchargMAX):",
!     .   IchargMAX, Mcharg(IchargMAX)

* CMS maximal cross section
      XMIN=100d0
      XMAX=700d0
      YMIN=0d0
      YMAX=400d0
      DO IX=1,NX
       DO IY=1,NY
        HCMS(IX,IY)=1d99
       ENDDO
      ENDDO
      FILENAME=catpath(EXPCON_PATH,'cms-pas-sus-17-004_7.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
 3001 READ(11,*,ERR=2,END=3002)NCMS,LCMS
      IF(HCMS(MOD(NCMS,NX+2),(NCMS-1)/(NX+2)).GT.LCMS)THEN
       HCMS(MOD(NCMS,NX+2),(NCMS-1)/(NX+2))=LCMS
      ENDIF
      GOTO 3001
 3002 CLOSE(11)
      FILENAME=catpath(EXPCON_PATH,'cms-pas-sus-17-004_8.dat')
      OPEN(11,FILE=FILENAME,STATUS='UNKNOWN',ERR=1)
 3003 READ(11,*,ERR=2,END=3004)NCMS,LCMS
      IF(HCMS(MOD(NCMS,NX+2),(NCMS-1)/(NX+2)).GT.LCMS)THEN
       HCMS(MOD(NCMS,NX+2),(NCMS-1)/(NX+2))=LCMS
      ENDIF
      GOTO 3003
 3004 CLOSE(11)

      RETURN

*   Error catch

 1    WRITE(*,*)"Cannot find the file ",FILENAME
      STOP

 2    WRITE(*,*)"Read error in the file ",FILENAME
      STOP

 3    WRITE(*,*)"Read error in the file ",FILENAME," Line: ",I
      STOP

      END


      integer function strlen(st)

*     Logical length of a string (omitting blanks)
*
*     | F | O | R | T | R | A | N |  |  |  |  |  |  |  |
*     <------------ full length returned by len() ----->
*     <---- logical length ------> <- trailing blank -->
*
      implicit none
      character      st*(*)
      strlen = len(st)
      do while (st(strlen:strlen) .eq. ' ')
       strlen = strlen - 1
      enddo
      return
      end


      character*(*) function catpath(st1,st2)

*     st1 is a string containing the path.
*     st2 is the file name or an other directory name
*     return value is the string "st1/st2"

      implicit none
      character*(*) st1
      character*(*) st2
      integer strlen
      catpath = st1(1:strlen(st1)) // '/' // st2
      return
      end


      SUBROUTINE TRITAB(TAB,NTAB,TRI,NTRI)

      IMPLICIT NONE
      INTEGER I,NTAB,NTRI,N1,N2,N3,NM,TRI(10000,3)
      DOUBLE PRECISION TAB(10000,3),X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3

      NTRI=0
      N1=NTAB
 3    X1=TAB(N1,1)
      Y1=TAB(N1,2)
      N2=N1
      X2=TAB(N2,1)
      Y2=TAB(N2,2)
      DO WHILE(N2.GT.1 .AND. DABS(X2/X1-1d0).LT.1d-2)
       N2=N2-1
       X2=TAB(N2,1)
       Y2=TAB(N2,2)
      END DO
      NM=N2

 4    IF(N1.GT.1 .AND. DABS(TAB(N1-1,1)/X1-1d0).LT.1d-2 .AND.
     .   N2.GT.1 .AND. DABS(TAB(N2-1,1)/X2-1d0).LT.1d-2)THEN
       IF(TAB(N1-1,2).GE.TAB(N2-1,2))THEN
        NTRI=NTRI+1
        N3=N1-1
        X3=TAB(N3,1)
        Y3=TAB(N3,2)
        TRI(NTRI,1)=N1
        TRI(NTRI,2)=N2
        TRI(NTRI,3)=N3
        N1=N3
        Y1=Y3
        GOTO 4
       ELSE
        NTRI=NTRI+1
        N3=N2-1
        X3=TAB(N3,1)
        Y3=TAB(N3,2)
        TRI(NTRI,1)=N1
        TRI(NTRI,2)=N2
        TRI(NTRI,3)=N3
        N2=N3
        Y2=Y3
        GOTO 4
       ENDIF
      ELSEIF(N1.GT.1 .AND. DABS(TAB(N1-1,1)/X1-1d0).LT.1d-2)THEN
       NTRI=NTRI+1
       N3=N1-1
       X3=TAB(N3,1)
       Y3=TAB(N3,2)
       TRI(NTRI,1)=N1
       TRI(NTRI,2)=N2
       TRI(NTRI,3)=N3
       N1=N3
       Y1=Y3
       GOTO 4
      ELSEIF(N2.GT.1 .AND. DABS(TAB(N2-1,1)/X2-1d0).LT.1d-2)THEN
       NTRI=NTRI+1
       N3=N2-1
       X3=TAB(N3,1)
       Y3=TAB(N3,2)
       TRI(NTRI,1)=N1
       TRI(NTRI,2)=N2
       TRI(NTRI,3)=N3
       N2=N3
       Y2=Y3
       GOTO 4
      ELSEIF(DABS(TAB(1,1)/X2-1d0).GT.1d-2)THEN
       N1=NM
       GOTO 3
      ENDIF

      END
