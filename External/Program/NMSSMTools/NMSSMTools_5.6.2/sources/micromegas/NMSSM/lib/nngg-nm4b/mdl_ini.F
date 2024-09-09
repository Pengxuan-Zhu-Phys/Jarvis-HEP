*     LanHEP output produced at Mon Apr  6 00:14:31 2015
*     Model named 'NMSSM(../spect.slha)+hgg'


      double precision function myif(a,b,c)
      double precision a,b,c
      if(a.lt.0D0) then
      myif=b
      else
      myif=c
      endif
      return 
      end

      subroutine ModelDefaults

      implicit none

      integer argc


#include "model.h"

      EE = 0.312853D0
      MZ = 91.187D0
      MW = 80.385D0
      GG = 1.117D0
      alfSMZ = 0.1172D0
      MbMb = 4.214D0
      McMc = 1.3D0
      Mtp = 173.1D0
      Ml = 1.777D0
      Mq = 0.05D0
      Q = 100D0
      Au = 0D0
      Ad = 0D0
      Mele = 0.000511D0
      Mmu = 0.1D0

      end


      integer function ModelConstIni(mode,Wh,err)
      implicit none
      integer mode, err
      real*8 Wh(5)

#include "model.h"

      double precision 
     > aPrintF,openAppend,System,slhaRead,slhaVal,
     > rDiagonal,MtRun,slhaValFormat,MassArray,MixMatrix,
     > MixMatrixU,slhaVal1,slhaVal2,slhaWidth,initQCD5,
     > MtEff,MbEff,McEff,alphaQCD,MbRun, myif, slhaValExists1
      double precision findValW	

      call ModelDefaults
      CW = MW/MZ
      SW = sqrt(1D0-CW**2)
      C2W = CW**2-SW**2
      err=0
c      LamQCD = initQCD5(alfSMZ, McMc, MbMb, Mtp)
c      rd = slhaRead('../spectr.slha', 0D0)
      if(mode.lt.2) then
      mu = slhaVal1('HMIX', MZ, 1)
      Lambda = slhaVal1('NMSSMRUN', MZ, 1)
      Kappa = slhaVal1('NMSSMRUN', MZ, 2)
      aLambda = slhaVal1('NMSSMRUN', MZ, 3)
      aKappa = slhaVal1('NMSSMRUN', MZ, 4)
      muP = slhaVal1('NMSSMRUN', MZ, 8)
      At = slhaVal2('Au', MZ, 3, 3)
      Ab = slhaVal2('Ad', MZ, 3, 3)
      Al = slhaVal2('Ae', MZ, 3, 3)

      Ml2 = slhaVal1('MSOFT', MZ, 32)
      Ml3 = slhaVal1('MSOFT', MZ, 33)
      Mr2 = slhaVal1('MSOFT', MZ, 35)
      Mr3 = slhaVal1('MSOFT', MZ, 36)
      Mq2 = slhaVal1('MSOFT', MZ, 42)
      Mq3 = slhaVal1('MSOFT', MZ, 43)
      Mu2 = slhaVal1('MSOFT', MZ, 45)
      Mu3 = slhaVal1('MSOFT', MZ, 46)
      Md2 = slhaVal1('MSOFT', MZ, 48)
      Md3 = slhaVal1('MSOFT', MZ, 49)
      if(slhaValExists1('HMIX',2).GT.0) then
          tb= slhaVal1('HMIX',Q,2) 
        else 
          tb=  slhaVal1('MINPAR',Q,3)
        endif
      Mh1 = slhaVal1('MASS', MZ, 25)
      Mh2 = slhaVal1('MASS', MZ, 35)
      Mh3 = slhaVal1('MASS', MZ, 45)
      Mha = slhaVal1('MASS', MZ, 36)
      Mhb = slhaVal1('MASS', MZ, 46)
      MHc = slhaVal1('MASS', MZ, 37)
      MNE1 = slhaVal1('MASS', MZ, 1000022)
      MNE2 = slhaVal1('MASS', MZ, 1000023)
      MNE3 = slhaVal1('MASS', MZ, 1000025)
      MNE4 = slhaVal1('MASS', MZ, 1000035)
      MNE5 = slhaVal1('MASS', MZ, 1000045)
      MC1 = slhaVal1('MASS', MZ, 1000024)
      MC2 = slhaVal1('MASS', MZ, 1000037)
      MSG = slhaVal1('MASS', MZ, 1000021)
      MSuL = slhaVal1('MASS', MZ, 1000002)
      MSuR = slhaVal1('MASS', MZ, 2000002)
      MSdL = slhaVal1('MASS', MZ, 1000001)
      MSdR = slhaVal1('MASS', MZ, 2000001)
      MScL = slhaVal1('MASS', MZ, 1000004)
      MScR = slhaVal1('MASS', MZ, 2000004)
      MSsL = slhaVal1('MASS', MZ, 1000003)
      MSsR = slhaVal1('MASS', MZ, 2000003)
      MSt1 = slhaVal1('MASS', MZ, 1000006)
      MSt2 = slhaVal1('MASS', MZ, 2000006)
      MSb1 = slhaVal1('MASS', MZ, 1000005)
      MSb2 = slhaVal1('MASS', MZ, 2000005)
      QSUSY = sqrt(MSt1*MSt2)
      MSeL = slhaVal1('MASS', MZ, 1000011)
      MSeR = slhaVal1('MASS', MZ, 2000011)
      MSmL = slhaVal1('MASS', MZ, 1000013)
      MSmR = slhaVal1('MASS', MZ, 2000013)
      MSl1 = slhaVal1('MASS', MZ, 1000015)
      MSl2 = slhaVal1('MASS', MZ, 2000015)
      MSne = slhaVal1('MASS', MZ, 1000012)
      MSnm = slhaVal1('MASS', MZ, 1000014)
      MSnl = slhaVal1('MASS', MZ, 1000016)
      Zl11 = slhaVal2('STAUMIX', QSUSY, 1, 1)
      Zl12 = slhaVal2('STAUMIX', QSUSY, 1, 2)
      Zl21 = slhaVal2('STAUMIX', QSUSY, 2, 1)
      Zl22 = slhaVal2('STAUMIX', QSUSY, 2, 2)
      Zh11 = slhaVal2('NMHMIX', QSUSY, 1, 1)
      Zh12 = slhaVal2('NMHMIX', QSUSY, 1, 2)
      Zh13 = slhaVal2('NMHMIX', QSUSY, 1, 3)
      Zh21 = slhaVal2('NMHMIX', QSUSY, 2, 1)
      Zh22 = slhaVal2('NMHMIX', QSUSY, 2, 2)
      Zh23 = slhaVal2('NMHMIX', QSUSY, 2, 3)
      Zh31 = slhaVal2('NMHMIX', QSUSY, 3, 1)
      Zh32 = slhaVal2('NMHMIX', QSUSY, 3, 2)
      Zh33 = slhaVal2('NMHMIX', QSUSY, 3, 3)
      Za11 = slhaVal2('NMAMIX', QSUSY, 1, 1)
      Za12 = slhaVal2('NMAMIX', QSUSY, 1, 2)
      Za13 = slhaVal2('NMAMIX', QSUSY, 1, 3)
      Za21 = slhaVal2('NMAMIX', QSUSY, 2, 1)
      Za22 = slhaVal2('NMAMIX', QSUSY, 2, 2)
      Za23 = slhaVal2('NMAMIX', QSUSY, 2, 3)
      Zn11 = slhaVal2('NMNMIX', QSUSY, 1, 1)
      Zn12 = slhaVal2('NMNMIX', QSUSY, 1, 2)
      Zn13 = slhaVal2('NMNMIX', QSUSY, 1, 3)
      Zn14 = slhaVal2('NMNMIX', QSUSY, 1, 4)
      Zn15 = slhaVal2('NMNMIX', QSUSY, 1, 5)
      Zn21 = slhaVal2('NMNMIX', QSUSY, 2, 1)
      Zn22 = slhaVal2('NMNMIX', QSUSY, 2, 2)
      Zn23 = slhaVal2('NMNMIX', QSUSY, 2, 3)
      Zn24 = slhaVal2('NMNMIX', QSUSY, 2, 4)
      Zn25 = slhaVal2('NMNMIX', QSUSY, 2, 5)
      Zn31 = slhaVal2('NMNMIX', QSUSY, 3, 1)
      Zn32 = slhaVal2('NMNMIX', QSUSY, 3, 2)
      Zn33 = slhaVal2('NMNMIX', QSUSY, 3, 3)
      Zn34 = slhaVal2('NMNMIX', QSUSY, 3, 4)
      Zn35 = slhaVal2('NMNMIX', QSUSY, 3, 5)
      Zn41 = slhaVal2('NMNMIX', QSUSY, 4, 1)
      Zn42 = slhaVal2('NMNMIX', QSUSY, 4, 2)
      Zn43 = slhaVal2('NMNMIX', QSUSY, 4, 3)
      Zn44 = slhaVal2('NMNMIX', QSUSY, 4, 4)
      Zn45 = slhaVal2('NMNMIX', QSUSY, 4, 5)
      Zn51 = slhaVal2('NMNMIX', QSUSY, 5, 1)
      Zn52 = slhaVal2('NMNMIX', QSUSY, 5, 2)
      Zn53 = slhaVal2('NMNMIX', QSUSY, 5, 3)
      Zn54 = slhaVal2('NMNMIX', QSUSY, 5, 4)
      Zn55 = slhaVal2('NMNMIX', QSUSY, 5, 5)
      Zu11 = slhaVal2('UMIX', QSUSY, 1, 1)
      Zu12 = slhaVal2('UMIX', QSUSY, 1, 2)
      Zu21 = slhaVal2('UMIX', QSUSY, 2, 1)
      Zu22 = slhaVal2('UMIX', QSUSY, 2, 2)
      Zv11 = slhaVal2('VMIX', QSUSY, 1, 1)
      Zv12 = slhaVal2('VMIX', QSUSY, 1, 2)
      Zv21 = slhaVal2('VMIX', QSUSY, 2, 1)
      Zv22 = slhaVal2('VMIX', QSUSY, 2, 2)
      Zt11 = slhaVal2('STOPMIX', QSUSY, 1, 1)
      Zt12 = slhaVal2('STOPMIX', QSUSY, 1, 2)
      Zt21 = slhaVal2('STOPMIX', QSUSY, 2, 1)
      Zt22 = slhaVal2('STOPMIX', QSUSY, 2, 2)
      Zb11 = slhaVal2('SBOTMIX', QSUSY, 1, 1)
      Zb12 = slhaVal2('SBOTMIX', QSUSY, 1, 2)
      Zb21 = slhaVal2('SBOTMIX', QSUSY, 2, 1)
      Zb22 = slhaVal2('SBOTMIX', QSUSY, 2, 2)
	  
	  
      la1 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'L1
     &       %lf')
      la2 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'L2
     &       %lf')
      la3 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'L3
     &       %lf')
      la4 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'L4
     &       %lf')
      la5 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'L5
     &       %lf')
      la6 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'L6
     &       %lf')
      la7 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'L7
     &       %lf')
      la1s = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'K1
     &       %lf')
      la2s = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'K2
     &       %lf')
      la3s = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'K3
     &       %lf')
      la4s = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'K4
     &       %lf')
      la5s = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'K5
     &       %lf')
      la6s = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'K6
     &       %lf')
      la7s = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'K7
     &       %lf')
      la8s = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'K8
     &       %lf')
      aa1 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'A1
     &       %lf')
      aa2 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'A2
     &       %lf')
      aa3 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'A3
     &       %lf')
      aa4 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'A4
     &       %lf')
      aa5 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'A5
     &       %lf')
      aa6 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'A6
     &       %lf')
      mB1 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'B1
     &       %lf')
      mB2 = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'B2
     &       %lf')
      X = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'X
     &       %lf')
      dMb = slhaValFormat('EFFECTIVE_COUPLINGS', QSUSY, 'DELMB
     &       %lf')
      xif = slhaVal1('NMSSMRUN', QSUSY, 6)
      xis = slhaVal1('NMSSMRUN', QSUSY, 7)
      MM3 = slhaVal1('NMSSMRUN', QSUSY, 9)
C      MSP = slhaVal1('NMSSMRUN', QSUSY, 12)
      
      else

      mu = findValW('mu')
      Lambda = findValW('Lambda')
      Kappa = findValW('Kappa')
      aLambda = findValW('aLambda')
      aKappa = findValW('aKappa')
      muP = findValW('muP')
      At = findValW('At')
      Ab = findValW('Ab')
      Al = findValW('Al')
        Ml2=findValW('Ml2')
        Ml3=findValW('Ml3')  
        Mr2=findValW('Mr2')
        Mr3=findValW('Mr3')
        Mq2=findValW('Mq2') 
        Mq3=findValW('Mq2')
        Mu2=findValW('Mu2')  
        Mu3=findValW('Mu3')      
        Md2=findValW('Md2')
        Md3=findValW('Md3')

        Mh1  = findValW('Mh1')
        Mh2  = findValW('Mh2')
        Mh3  = findValW('Mh3')

        Zh11 = findValW('Zh11')
        Zh12 = findValW('Zh12')
        Zh13 = findValW('Zh13')
        Zh21 = findValW('Zh21')
        Zh22 = findValW('Zh22')
        Zh23 = findValW('Zh23')
        Zh31 = findValW('Zh31')
        Zh32 = findValW('Zh32')
        Zh33 = findValW('Zh33')

        Mha  = findValW('Mha')
        Mhb  = findValW('Mhb')

        Za11=  findValW('Z11')
        Za12=  findValW('Z12')
        Za13=  findValW('Z13')
        Za21=  findValW('Z21')
        Za22=  findValW('Z22')
        Za23=  findValW('Z23')


        MNE1=findValW('MNE1')
        MNE2=findValW('MNE2')
        MNE3=findValW('MNE3')
        MNE4=findValW('MNE4')
        MNE5=findValW('MNE5')

        Zn11  =findValW('Zn11')
        Zn12  =findValW('Zn12')
        Zn13  =findValW('Zn13')
        Zn14  =findValW('Zn14')
        Zn15  =findValW('Zn15')
        Zn21  =findValW('Zn21')
        Zn22  =findValW('Zn22')
        Zn23  =findValW('Zn23')
        Zn24  =findValW('Zn24')
        Zn25  =findValW('Zn25')
        Zn31  =findValW('Zn31')
        Zn32  =findValW('Zn32')
        Zn33  =findValW('Zn33')
        Zn34  =findValW('Zn34')
        Zn35  =findValW('Zn35')
        Zn41  =findValW('Zn41')
        Zn42  =findValW('Zn42')
        Zn43  =findValW('Zn43')
        Zn44  =findValW('Zn44')
        Zn45  =findValW('Zn45')
        Zn51  =findValW('Zn51')
        Zn52  =findValW('Zn52')
        Zn53  =findValW('Zn53')
        Zn54  =findValW('Zn54')
        Zn55  =findValW('Zn55')


        MC1   =findValW('MC1')
        MC2   =findValW('MC2')
        Zu11  =findValW('Zu11')
        Zu21  =findValW('Zu21')
        Zu12  =findValW('Zu12')
        Zu22  =findValW('Zu22')
        Zv11  =findValW('Zv11')
        Zv21  =findValW('Zv21')
        Zv12  =findValW('Zv12')
        Zv22  =findValW('Zv22')

        MSG = findValW('MSG')

        MSeL= findValW('MSeL')
        MSmL= findValW('MSmL')
        MSl1= findValW('MSl1')
        MSuL= findValW('MSuL')
        MSdL= findValW('MSdL')
        MScL= findValW('MScL')
        MSsL= findValW('MSsL')
        MSb1= findValW('MSb1')
        MSt1= findValW('MSt1')
        MSeR= findValW('MSeR')
        MSmR= findValW('MSmR')
        MSl2= findValW('MSl2')
        MSuR= findValW('MSuR')
        MSdR= findValW('MSdR')
        MScR= findValW('MScR')
        MSsR= findValW('MSsR')
        MSb2= findValW('MSb2')
        MSt2= findValW('MSt2')
        MSne= findValW('MSne')
        MSnm= findValW('MSnm')
        MSnl= findValW('MSnl')
 
        Zl11 =findValW('Zl11')
        Zl21 =findValW('Zl21')
        Zl12 =findValW('Zl12')
        Zl22 =findValW('Zl22')

        Zb11 =findValW('Zb11')
        Zb21 =findValW('Zb21')
        Zb12 =findValW('Zb12')
        Zb22 =findValW('Zb22')

        Zt11 =findValW('Zt11')
        Zt21 =findValW('Zt21')
        Zt12 =findValW('Zt12')
        Zt22 =findValW('Zt22')

      la1 = findValW('la1')
      la2 = findValW('la2')
      la3 = findValW('la3')
      la4 = findValW('la4')
      la5 = findValW('la5')
      la6 = findValW('la6')
      la7 = findValW('la7')

      la1s = findValW('la1s')
      la2s = findValW('la2s')
      la3s = findValW('la3s')
      la4s = findValW('la4s')
      la5s = findValW('la5s')
      la6s = findValW('la6s')
      la7s = findValW('la7s')
      la8s = findValW('la8s')

      aa1 = findValW('aa1')
      aa2 = findValW('aa2')
      aa3 = findValW('aa3')
      aa4 = findValW('aa4')
      aa5 = findValW('aa5')
      aa6 = findValW('aa6')

      mB1 = findValW('mB1')
      mB2 = findValW('mB2')
      X = findValW('X')
      dMB = findValW('dMB')
      xif = findValW('xif')
      xis = findValW('xis')
      MM3 = findValW('MM3')
      MSP = findValW('MSP') 

       endif

      if(MNE1.lt.0D0) then
	     MNE1=-MNE1
		 Zn11=Zn11*(0D0,1D0)
		 Zn12=Zn12*(0D0,1D0)
		 Zn13=Zn13*(0D0,1D0)
		 Zn14=Zn14*(0D0,1D0)
		 Zn15=Zn15*(0D0,1D0)
	  endif

      hiaW(1)=Wh(4)
      hiaW(2)=Wh(5)


      sb = tB/sqrt(1D0+tB**2)
      cb = sqrt(1D0-sb**2)
      t2b = 2D0*tB/(1D0-tB**2)
      xvev = mu/Lambda
      Pa12 = Za13
      Pa22 = Za23
      Pa11 = myif(Za11*Pa22, Pa22, -Pa22)
      Pa21 = myif(Za21*Pa12, Pa12, -Pa12)
      Td3 = dMb/(1D0+dMb)
      Mt = MtEff(Q)
      Mb = MbEff(Q)
      Mc = McEff(Q)
      xH2 = MSsL**2/MSG**2
      yH2 = MSsR**2/MSG**2
      dMd = -3D0/2D0*alphaQCD(MSG)/3.1415D0*tB*mu/MSG*(xH2
     &      *log(xH2)/(1D0-xH2)-yH2*log(yH2)/(1D0-yH2))/(xH2
     &      -yH2)
      Td2 = dMd/(1D0+dMd)
      fiuu = atan(-2D0*Mq*(Au-mu/tB)/(MSuR**2-MSuL**2))/2D0
      fidd = atan(-2D0*Mq*(1D0-Td2)*(Ad-mu*tB)/(MSdR**2-MSdL
     &      **2))/2D0
      ficc = atan(-2D0*Mc*(Au-mu/tB)/(MScR**2-MScL**2))/2D0
      fiss = atan(-2D0*Mq*(1D0-Td2)*(Ad-mu*tB)/(MSsR**2-MSsL
     &      **2))/2D0
      Zuu11 = cos(fiuu)
      Zuu12 = sin(fiuu)
      Zuu21 = -sin(fiuu)
      Zuu22 = cos(fiuu)
      Zdd11 = cos(fidd)
      Zdd12 = sin(fidd)
      Zdd21 = -sin(fidd)
      Zdd22 = cos(fidd)
      Zcc11 = cos(ficc)
      Zcc12 = sin(ficc)
      Zcc21 = -sin(ficc)
      Zcc22 = cos(ficc)
      Zss11 = cos(fiss)
      Zss12 = sin(fiss)
      Zss21 = -sin(fiss)
      Zss22 = cos(fiss)
      MbMM = MbRun(sqrt(MSb1*MSb2))
      MtMM = MtRun(sqrt(MSt1*MSt2))
      NMM11 = Zn11*MNE1*Zn11+Zn21*MNE2*Zn21+Zn31*MNE3*Zn31
     &      +Zn41*MNE4*Zn41+Zn51*MNE5*Zn51
      NMM12 = Zn11*MNE1*Zn12+Zn21*MNE2*Zn22+Zn31*MNE3*Zn32
     &      +Zn41*MNE4*Zn42+Zn51*MNE5*Zn52
      NMM13 = Zn11*MNE1*Zn13+Zn21*MNE2*Zn23+Zn31*MNE3*Zn33
     &      +Zn41*MNE4*Zn43+Zn51*MNE5*Zn53
      NMM14 = Zn11*MNE1*Zn14+Zn21*MNE2*Zn24+Zn31*MNE3*Zn34
     &      +Zn41*MNE4*Zn44+Zn51*MNE5*Zn54
      NMM15 = Zn11*MNE1*Zn15+Zn21*MNE2*Zn25+Zn31*MNE3*Zn35
     &      +Zn41*MNE4*Zn45+Zn51*MNE5*Zn55
      NMM22 = Zn12*MNE1*Zn12+Zn22*MNE2*Zn22+Zn32*MNE3*Zn32
     &      +Zn42*MNE4*Zn42+Zn52*MNE5*Zn52
      NMM23 = Zn12*MNE1*Zn13+Zn22*MNE2*Zn23+Zn32*MNE3*Zn33
     &      +Zn42*MNE4*Zn43+Zn52*MNE5*Zn53
      NMM24 = Zn12*MNE1*Zn14+Zn22*MNE2*Zn24+Zn32*MNE3*Zn34
     &      +Zn42*MNE4*Zn44+Zn52*MNE5*Zn54
      NMM25 = Zn12*MNE1*Zn15+Zn22*MNE2*Zn25+Zn32*MNE3*Zn35
     &      +Zn42*MNE4*Zn45+Zn52*MNE5*Zn55
      NMM33 = Zn13*MNE1*Zn13+Zn23*MNE2*Zn23+Zn33*MNE3*Zn33
     &      +Zn43*MNE4*Zn43+Zn53*MNE5*Zn53
      NMM34 = Zn13*MNE1*Zn14+Zn23*MNE2*Zn24+Zn33*MNE3*Zn34
     &      +Zn43*MNE4*Zn44+Zn53*MNE5*Zn54
      NMM35 = Zn13*MNE1*Zn15+Zn23*MNE2*Zn25+Zn33*MNE3*Zn35
     &      +Zn43*MNE4*Zn45+Zn53*MNE5*Zn55
      NMM44 = Zn14*MNE1*Zn14+Zn24*MNE2*Zn24+Zn34*MNE3*Zn34
     &      +Zn44*MNE4*Zn44+Zn54*MNE5*Zn54
      NMM45 = Zn14*MNE1*Zn15+Zn24*MNE2*Zn25+Zn34*MNE3*Zn35
     &      +Zn44*MNE4*Zn45+Zn54*MNE5*Zn55
      NMM55 = Zn15*MNE1*Zn15+Zn25*MNE2*Zn25+Zn35*MNE3*Zn35
     &      +Zn45*MNE4*Zn45+Zn55*MNE5*Zn55
      MG1I = NMM11
      MG2I = NMM22

      call aaini01
      call mtrini
      if(2*MNE1.gt.MZ) then 
         ModelConstIni= 2
      else 
         ModelConstIni= 1
      endif 


      end

      subroutine aaini01
      implicit none
#include "model.h"

      AAABR(1) = MW*(1D0/CW-CW-SW**2/CW)
      AAABR(2) = CW**2+SW**2
      AAABR(3) = 1D0-CW**2-SW**2
      AAABR(4) = MW**2*(CW**2+2D0*SW**2+SW**4/CW**2)
      AAABR(5) = EE*MW
      AAABR(6) = EE*MW/SW
      AAABR(7) = EE*MW/SW*(Zh21*cb+Zh22*sb)
      AAABR(8) = EE*MW/SW*(Zh31*cb+Zh32*sb)
      AAABR(9) = EE*MW/SW*(Zh11*cb+Zh12*sb)
      AAABR(10) = EE*MW*(SW/CW-CW/SW)
      AAABR(11) = EE/CW*MW/SW
      AAABR(12) = EE/CW**2*MW/SW*(Zh11*cb+Zh12*sb)
      AAABR(13) = EE/CW**2*MW/SW*(Zh31*cb+Zh32*sb)
      AAABR(14) = EE/CW**2*MW/SW*(Zh21*cb+Zh22*sb)
      AAABR(15) = CW*EE/SW
      AAABR(16) = Sqrt2*Zh23*aa2+Sqrt2*Zh23*aa1*sb**2-Sqrt2
     &      *Zh23*aa2*sb**2+Sqrt2*Zh23*la2s*xvev+Sqrt2*Zh23
     &      *aa5*cb*sb+Sqrt2*Zh23*aa6*cb*sb+Sqrt2*Zh23*la1s
     &      *sb**2*xvev-Sqrt2*Zh23*la2s*sb**2*xvev+2D0*MW
     &      /EE*SW*Zh21*cb*la3+4D0*MW/EE*SW*Zh21*la6*sb-6D0
     &      *MW/EE*SW*Zh21*la6*sb**3-2D0*MW/EE*SW*Zh21*la7
     &      *sb+2D0*MW/EE*SW*Zh21*la7*sb**3-2D0*MW/EE*SW
     &      *Zh22*cb*la7+2D0*MW/EE*SW*Zh22*la2*sb-2D0*MW
     &      /EE*SW*Zh22*la2*sb**3+2D0*MW/EE*SW*Zh22*la3*sb
     &      **3-2D0*MW/EE*SW*Zh22*la4*sb+2D0*MW/EE*SW*Zh22
     &      *la4*sb**3-2D0*MW/EE*SW*Zh22*la5*sb+2D0*MW/EE
     &      *SW*Zh22*la5*sb**3+2D0*Sqrt2*Zh23*cb*la4s*sb
     &      *xvev+2D0*Sqrt2*Zh23*cb*la5s*sb*xvev+2D0*Sqrt2
     &      *Zh23*cb*la6s*sb*xvev+2D0*MW/EE*SW*Zh21*cb*la1
     &      *sb**2-2D0*MW/EE*SW*Zh21*cb*la3*sb**2-2D0*MW
     &      /EE*SW*Zh21*cb*la4*sb**2-2D0*MW/EE*SW*Zh21*cb
     &      *la5*sb**2-2D0*MW/EE*SW*Zh22*cb*la6*sb**2+6D0
     &      *MW/EE*SW*Zh22*cb*la7*sb**2
      AAABR(17) = Sqrt2*Zh33*aa2+Sqrt2*Zh33*aa1*sb**2-Sqrt2
     &      *Zh33*aa2*sb**2+Sqrt2*Zh33*la2s*xvev+Sqrt2*Zh33
     &      *aa5*cb*sb+Sqrt2*Zh33*aa6*cb*sb+Sqrt2*Zh33*la1s
     &      *sb**2*xvev-Sqrt2*Zh33*la2s*sb**2*xvev+2D0*MW
     &      /EE*SW*Zh31*cb*la3+4D0*MW/EE*SW*Zh31*la6*sb-6D0
     &      *MW/EE*SW*Zh31*la6*sb**3-2D0*MW/EE*SW*Zh31*la7
     &      *sb+2D0*MW/EE*SW*Zh31*la7*sb**3-2D0*MW/EE*SW
     &      *Zh32*cb*la7+2D0*MW/EE*SW*Zh32*la2*sb-2D0*MW
     &      /EE*SW*Zh32*la2*sb**3+2D0*MW/EE*SW*Zh32*la3*sb
     &      **3-2D0*MW/EE*SW*Zh32*la4*sb+2D0*MW/EE*SW*Zh32
     &      *la4*sb**3-2D0*MW/EE*SW*Zh32*la5*sb+2D0*MW/EE
     &      *SW*Zh32*la5*sb**3+2D0*Sqrt2*Zh33*cb*la4s*sb
     &      *xvev+2D0*Sqrt2*Zh33*cb*la5s*sb*xvev+2D0*Sqrt2
     &      *Zh33*cb*la6s*sb*xvev+2D0*MW/EE*SW*Zh31*cb*la1
     &      *sb**2-2D0*MW/EE*SW*Zh31*cb*la3*sb**2-2D0*MW
     &      /EE*SW*Zh31*cb*la4*sb**2-2D0*MW/EE*SW*Zh31*cb
     &      *la5*sb**2-2D0*MW/EE*SW*Zh32*cb*la6*sb**2+6D0
     &      *MW/EE*SW*Zh32*cb*la7*sb**2
      AAABR(18) = Sqrt2*Zh13*aa2+Sqrt2*Zh13*aa1*sb**2-Sqrt2
     &      *Zh13*aa2*sb**2+Sqrt2*Zh13*la2s*xvev+Sqrt2*Zh13
     &      *aa5*cb*sb+Sqrt2*Zh13*aa6*cb*sb+Sqrt2*Zh13*la1s
     &      *sb**2*xvev-Sqrt2*Zh13*la2s*sb**2*xvev+2D0*MW
     &      /EE*SW*Zh11*cb*la3+4D0*MW/EE*SW*Zh11*la6*sb-6D0
     &      *MW/EE*SW*Zh11*la6*sb**3-2D0*MW/EE*SW*Zh11*la7
     &      *sb+2D0*MW/EE*SW*Zh11*la7*sb**3-2D0*MW/EE*SW
     &      *Zh12*cb*la7+2D0*MW/EE*SW*Zh12*la2*sb-2D0*MW
     &      /EE*SW*Zh12*la2*sb**3+2D0*MW/EE*SW*Zh12*la3*sb
     &      **3-2D0*MW/EE*SW*Zh12*la4*sb+2D0*MW/EE*SW*Zh12
     &      *la4*sb**3-2D0*MW/EE*SW*Zh12*la5*sb+2D0*MW/EE
     &      *SW*Zh12*la5*sb**3+2D0*Sqrt2*Zh13*cb*la4s*sb
     &      *xvev+2D0*Sqrt2*Zh13*cb*la5s*sb*xvev+2D0*Sqrt2
     &      *Zh13*cb*la6s*sb*xvev+2D0*MW/EE*SW*Zh11*cb*la1
     &      *sb**2-2D0*MW/EE*SW*Zh11*cb*la3*sb**2-2D0*MW
     &      /EE*SW*Zh11*cb*la4*sb**2-2D0*MW/EE*SW*Zh11*cb
     &      *la5*sb**2-2D0*MW/EE*SW*Zh12*cb*la6*sb**2+6D0
     &      *MW/EE*SW*Zh12*cb*la7*sb**2
      AAABR(19) = EE*MW/SW*Sqrt2*cb*sb
      AAABR(20) = EE/SW*Sqrt2*(2D0*MW*Zl11/cb*sb-2D0*MW*Zl11
     &      /cb*sb**3-Lambda/MW*Ml*Zl12*xvev-Ml**2/MW*Zl11
     &      /cb*sb-Al/MW*Ml*Zl12/cb*sb)
      AAABR(21) = EE/SW*Sqrt2*(2D0*MW*Zl21/cb*sb-2D0*MW*Zl21
     &      /cb*sb**3-Lambda/MW*Ml*Zl22*xvev-Ml**2/MW*Zl21
     &      /cb*sb-Al/MW*Ml*Zl22/cb*sb)
      AAABR(22) = EE/SW*Sqrt2*(2D0*MW*Zdd11*Zuu11/cb*sb-2D0
     &      *MW*Zdd11*Zuu11/cb*sb**3-Lambda/MW*Mq*Zdd11*Zuu12
     &      *xvev-Lambda/MW*Mq*Zdd12*Zuu11*xvev-Mq**2/MW
     &      *Zdd11*Zuu11/cb/sb-Mq**2/MW*Zdd12*Zuu12/cb/sb
     &      -Ad/MW*Mq*Zdd12*Zuu11/cb*sb-Au/MW*Mq*Zdd11*Zuu12
     &      /cb/sb+Au/MW*Mq*Zdd11*Zuu12/cb*sb)
      AAABR(23) = EE/SW*Sqrt2*(2D0*MW*Zb11*Zt11/cb*sb-2D0
     &      *MW*Zb11*Zt11/cb*sb**3-Lambda/MW*MbMM*Zb12*Zt11
     &      *xvev-Lambda/MW*MtMM*Zb11*Zt12*xvev-MbMM**2/MW
     &      *Zb11*Zt11/cb*sb-MtMM**2/MW*Zb11*Zt11/cb/sb+MtMM
     &      **2/MW*Zb11*Zt11/cb*sb-Ab/MW*MbMM*Zb12*Zt11/cb
     &      *sb-At/MW*MtMM*Zb11*Zt12/cb/sb+At/MW*MtMM*Zb11
     &      *Zt12/cb*sb-MbMM/MW*MtMM*Zb12*Zt12/cb/sb)
      AAABR(24) = EE/SW*Sqrt2*(2D0*MW*Zcc11*Zss11/cb*sb-2D0
     &      *MW*Zcc11*Zss11/cb*sb**3-Lambda/MW*Mc*Zcc12*Zss11
     &      *xvev-Lambda/MW*Mq*Zcc11*Zss12*xvev-Mc**2/MW
     &      *Zcc11*Zss11/cb/sb+Mc**2/MW*Zcc11*Zss11/cb*sb
     &      -Mq**2/MW*Zcc11*Zss11/cb*sb-Ad/MW*Mq*Zcc11*Zss12
     &      /cb*sb-Au/MW*Mc*Zcc12*Zss11/cb/sb+Au/MW*Mc*Zcc12
     &      *Zss11/cb*sb-Mc/MW*Mq*Zcc12*Zss12/cb/sb)
      AAABR(25) = EE/SW*Sqrt2*(2D0*MW*Zcc11*Zss21/cb*sb-2D0
     &      *MW*Zcc11*Zss21/cb*sb**3-Lambda/MW*Mc*Zcc12*Zss21
     &      *xvev-Lambda/MW*Mq*Zcc11*Zss22*xvev-Mc**2/MW
     &      *Zcc11*Zss21/cb/sb+Mc**2/MW*Zcc11*Zss21/cb*sb
     &      -Mq**2/MW*Zcc11*Zss21/cb*sb-Ad/MW*Mq*Zcc11*Zss22
     &      /cb*sb-Au/MW*Mc*Zcc12*Zss21/cb/sb+Au/MW*Mc*Zcc12
     &      *Zss21/cb*sb-Mc/MW*Mq*Zcc12*Zss22/cb/sb)
      AAABR(26) = EE/SW*Sqrt2*(2D0*MW*Zb21*Zt11/cb*sb-2D0
     &      *MW*Zb21*Zt11/cb*sb**3-Lambda/MW*MbMM*Zb22*Zt11
     &      *xvev-Lambda/MW*MtMM*Zb21*Zt12*xvev-MbMM**2/MW
     &      *Zb21*Zt11/cb*sb-MtMM**2/MW*Zb21*Zt11/cb/sb+MtMM
     &      **2/MW*Zb21*Zt11/cb*sb-Ab/MW*MbMM*Zb22*Zt11/cb
     &      *sb-At/MW*MtMM*Zb21*Zt12/cb/sb+At/MW*MtMM*Zb21
     &      *Zt12/cb*sb-MbMM/MW*MtMM*Zb22*Zt12/cb/sb)
      AAABR(27) = EE/SW*Sqrt2*(2D0*MW*Zdd21*Zuu11/cb*sb-2D0
     &      *MW*Zdd21*Zuu11/cb*sb**3-Lambda/MW*Mq*Zdd21*Zuu12
     &      *xvev-Lambda/MW*Mq*Zdd22*Zuu11*xvev-Mq**2/MW
     &      *Zdd21*Zuu11/cb/sb-Mq**2/MW*Zdd22*Zuu12/cb/sb
     &      -Ad/MW*Mq*Zdd22*Zuu11/cb*sb-Au/MW*Mq*Zdd21*Zuu12
     &      /cb/sb+Au/MW*Mq*Zdd21*Zuu12/cb*sb)
      AAABR(28) = EE/SW*Sqrt2*(2D0*MW*Zb11*Zt21/cb*sb-2D0
     &      *MW*Zb11*Zt21/cb*sb**3-Lambda/MW*MbMM*Zb12*Zt21
     &      *xvev-Lambda/MW*MtMM*Zb11*Zt22*xvev-MbMM**2/MW
     &      *Zb11*Zt21/cb*sb-MtMM**2/MW*Zb11*Zt21/cb/sb+MtMM
     &      **2/MW*Zb11*Zt21/cb*sb-Ab/MW*MbMM*Zb12*Zt21/cb
     &      *sb-At/MW*MtMM*Zb11*Zt22/cb/sb+At/MW*MtMM*Zb11
     &      *Zt22/cb*sb-MbMM/MW*MtMM*Zb12*Zt22/cb/sb)
      AAABR(29) = EE/SW*Sqrt2*(2D0*MW*Zcc21*Zss11/cb*sb-2D0
     &      *MW*Zcc21*Zss11/cb*sb**3-Lambda/MW*Mc*Zcc22*Zss11
     &      *xvev-Lambda/MW*Mq*Zcc21*Zss12*xvev-Mc**2/MW
     &      *Zcc21*Zss11/cb/sb+Mc**2/MW*Zcc21*Zss11/cb*sb
     &      -Mq**2/MW*Zcc21*Zss11/cb*sb-Ad/MW*Mq*Zcc21*Zss12
     &      /cb*sb-Au/MW*Mc*Zcc22*Zss11/cb/sb+Au/MW*Mc*Zcc22
     &      *Zss11/cb*sb-Mc/MW*Mq*Zcc22*Zss12/cb/sb)
      AAABR(30) = EE/SW*Sqrt2*(2D0*MW*Zdd11*Zuu21/cb*sb-2D0
     &      *MW*Zdd11*Zuu21/cb*sb**3-Lambda/MW*Mq*Zdd11*Zuu22
     &      *xvev-Lambda/MW*Mq*Zdd12*Zuu21*xvev-Mq**2/MW
     &      *Zdd11*Zuu21/cb/sb-Mq**2/MW*Zdd12*Zuu22/cb/sb
     &      -Ad/MW*Mq*Zdd12*Zuu21/cb*sb-Au/MW*Mq*Zdd11*Zuu22
     &      /cb/sb+Au/MW*Mq*Zdd11*Zuu22/cb*sb)
      AAABR(31) = EE/SW*Sqrt2*(2D0*MW*Zb21*Zt21/cb*sb-2D0
     &      *MW*Zb21*Zt21/cb*sb**3-Lambda/MW*MbMM*Zb22*Zt21
     &      *xvev-Lambda/MW*MtMM*Zb21*Zt22*xvev-MbMM**2/MW
     &      *Zb21*Zt21/cb*sb-MtMM**2/MW*Zb21*Zt21/cb/sb+MtMM
     &      **2/MW*Zb21*Zt21/cb*sb-Ab/MW*MbMM*Zb22*Zt21/cb
     &      *sb-At/MW*MtMM*Zb21*Zt22/cb/sb+At/MW*MtMM*Zb21
     &      *Zt22/cb*sb-MbMM/MW*MtMM*Zb22*Zt22/cb/sb)
      AAABR(32) = EE/SW*Sqrt2*(2D0*MW*Zcc21*Zss21/cb*sb-2D0
     &      *MW*Zcc21*Zss21/cb*sb**3-Lambda/MW*Mc*Zcc22*Zss21
     &      *xvev-Lambda/MW*Mq*Zcc21*Zss22*xvev-Mc**2/MW
     &      *Zcc21*Zss21/cb/sb+Mc**2/MW*Zcc21*Zss21/cb*sb
     &      -Mq**2/MW*Zcc21*Zss21/cb*sb-Ad/MW*Mq*Zcc21*Zss22
     &      /cb*sb-Au/MW*Mc*Zcc22*Zss21/cb/sb+Au/MW*Mc*Zcc22
     &      *Zss21/cb*sb-Mc/MW*Mq*Zcc22*Zss22/cb/sb)
      AAABR(33) = EE/SW*Sqrt2*(2D0*MW*Zdd21*Zuu21/cb*sb-2D0
     &      *MW*Zdd21*Zuu21/cb*sb**3-Lambda/MW*Mq*Zdd21*Zuu22
     &      *xvev-Lambda/MW*Mq*Zdd22*Zuu21*xvev-Mq**2/MW
     &      *Zdd21*Zuu21/cb/sb-Mq**2/MW*Zdd22*Zuu22/cb/sb
     &      -Ad/MW*Mq*Zdd22*Zuu21/cb*sb-Au/MW*Mq*Zdd21*Zuu22
     &      /cb/sb+Au/MW*Mq*Zdd21*Zuu22/cb*sb)
      AAABR(34) = Pa12*Sqrt2*aa5-Pa12*Sqrt2*aa6+2D0*Pa12
     &      *Sqrt2*la5s*xvev-2D0*Pa12*Sqrt2*la6s*xvev+2D0
     &      *MW/EE*Pa11*SW*la4-2D0*MW/EE*Pa11*SW*la5
      AAABR(35) = Pa22*Sqrt2*aa5-Pa22*Sqrt2*aa6+2D0*Pa22
     &      *Sqrt2*la5s*xvev-2D0*Pa22*Sqrt2*la6s*xvev+2D0
     &      *MW/EE*Pa21*SW*la4-2D0*MW/EE*Pa21*SW*la5
      AAABR(36) = Sqrt2*Zh13*aa5+Sqrt2*Zh13*aa6-2D0*Sqrt2
     &      *Zh13*aa5*sb**2-2D0*Sqrt2*Zh13*aa6*sb**2+2D0
     &      *Sqrt2*Zh13*la4s*xvev+2D0*Sqrt2*Zh13*la5s*xvev
     &      +2D0*Sqrt2*Zh13*la6s*xvev+2D0*Sqrt2*Zh13*aa1
     &      *cb*sb-2D0*Sqrt2*Zh13*aa2*cb*sb-4D0*Sqrt2*Zh13
     &      *la4s*sb**2*xvev-4D0*Sqrt2*Zh13*la5s*sb**2*xvev
     &      -4D0*Sqrt2*Zh13*la6s*sb**2*xvev+4D0*MW/EE*SW
     &      *Zh11*cb*la6+4D0*MW/EE*SW*Zh11*la1*sb-4D0*MW
     &      /EE*SW*Zh11*la1*sb**3-4D0*MW/EE*SW*Zh11*la3*sb
     &      +4D0*MW/EE*SW*Zh11*la3*sb**3-2D0*MW/EE*SW*Zh11
     &      *la4*sb+4D0*MW/EE*SW*Zh11*la4*sb**3-2D0*MW/EE
     &      *SW*Zh11*la5*sb+4D0*MW/EE*SW*Zh11*la5*sb**3-2D0
     &      *MW/EE*SW*Zh12*cb*la4-2D0*MW/EE*SW*Zh12*cb*la5
     &      -4D0*MW/EE*SW*Zh12*la6*sb+4D0*MW/EE*SW*Zh12*la6
     &      *sb**3+8D0*MW/EE*SW*Zh12*la7*sb-12D0*MW/EE*SW
     &      *Zh12*la7*sb**3+2D0*Sqrt2*Zh13*cb*la1s*sb*xvev
     &      -2D0*Sqrt2*Zh13*cb*la2s*sb*xvev-12D0*MW/EE*SW
     &      *Zh11*cb*la6*sb**2+4D0*MW/EE*SW*Zh11*cb*la7*sb
     &      **2-4D0*MW/EE*SW*Zh12*cb*la2*sb**2+4D0*MW/EE
     &      *SW*Zh12*cb*la3*sb**2+4D0*MW/EE*SW*Zh12*cb*la4
     &      *sb**2+4D0*MW/EE*SW*Zh12*cb*la5*sb**2
      AAABR(37) = Sqrt2*Zh23*aa5+Sqrt2*Zh23*aa6-2D0*Sqrt2
     &      *Zh23*aa5*sb**2-2D0*Sqrt2*Zh23*aa6*sb**2+2D0
     &      *Sqrt2*Zh23*la4s*xvev+2D0*Sqrt2*Zh23*la5s*xvev
     &      +2D0*Sqrt2*Zh23*la6s*xvev+2D0*Sqrt2*Zh23*aa1
     &      *cb*sb-2D0*Sqrt2*Zh23*aa2*cb*sb-4D0*Sqrt2*Zh23
     &      *la4s*sb**2*xvev-4D0*Sqrt2*Zh23*la5s*sb**2*xvev
     &      -4D0*Sqrt2*Zh23*la6s*sb**2*xvev+4D0*MW/EE*SW
     &      *Zh21*cb*la6+4D0*MW/EE*SW*Zh21*la1*sb-4D0*MW
     &      /EE*SW*Zh21*la1*sb**3-4D0*MW/EE*SW*Zh21*la3*sb
     &      +4D0*MW/EE*SW*Zh21*la3*sb**3-2D0*MW/EE*SW*Zh21
     &      *la4*sb+4D0*MW/EE*SW*Zh21*la4*sb**3-2D0*MW/EE
     &      *SW*Zh21*la5*sb+4D0*MW/EE*SW*Zh21*la5*sb**3-2D0
     &      *MW/EE*SW*Zh22*cb*la4-2D0*MW/EE*SW*Zh22*cb*la5
     &      -4D0*MW/EE*SW*Zh22*la6*sb+4D0*MW/EE*SW*Zh22*la6
     &      *sb**3+8D0*MW/EE*SW*Zh22*la7*sb-12D0*MW/EE*SW
     &      *Zh22*la7*sb**3+2D0*Sqrt2*Zh23*cb*la1s*sb*xvev
     &      -2D0*Sqrt2*Zh23*cb*la2s*sb*xvev-12D0*MW/EE*SW
     &      *Zh21*cb*la6*sb**2+4D0*MW/EE*SW*Zh21*cb*la7*sb
     &      **2-4D0*MW/EE*SW*Zh22*cb*la2*sb**2+4D0*MW/EE
     &      *SW*Zh22*cb*la3*sb**2+4D0*MW/EE*SW*Zh22*cb*la4
     &      *sb**2+4D0*MW/EE*SW*Zh22*cb*la5*sb**2
      AAABR(38) = Sqrt2*Zh33*aa5+Sqrt2*Zh33*aa6-2D0*Sqrt2
     &      *Zh33*aa5*sb**2-2D0*Sqrt2*Zh33*aa6*sb**2+2D0
     &      *Sqrt2*Zh33*la4s*xvev+2D0*Sqrt2*Zh33*la5s*xvev
     &      +2D0*Sqrt2*Zh33*la6s*xvev+2D0*Sqrt2*Zh33*aa1
     &      *cb*sb-2D0*Sqrt2*Zh33*aa2*cb*sb-4D0*Sqrt2*Zh33
     &      *la4s*sb**2*xvev-4D0*Sqrt2*Zh33*la5s*sb**2*xvev
     &      -4D0*Sqrt2*Zh33*la6s*sb**2*xvev+4D0*MW/EE*SW
     &      *Zh31*cb*la6+4D0*MW/EE*SW*Zh31*la1*sb-4D0*MW
     &      /EE*SW*Zh31*la1*sb**3-4D0*MW/EE*SW*Zh31*la3*sb
     &      +4D0*MW/EE*SW*Zh31*la3*sb**3-2D0*MW/EE*SW*Zh31
     &      *la4*sb+4D0*MW/EE*SW*Zh31*la4*sb**3-2D0*MW/EE
     &      *SW*Zh31*la5*sb+4D0*MW/EE*SW*Zh31*la5*sb**3-2D0
     &      *MW/EE*SW*Zh32*cb*la4-2D0*MW/EE*SW*Zh32*cb*la5
     &      -4D0*MW/EE*SW*Zh32*la6*sb+4D0*MW/EE*SW*Zh32*la6
     &      *sb**3+8D0*MW/EE*SW*Zh32*la7*sb-12D0*MW/EE*SW
     &      *Zh32*la7*sb**3+2D0*Sqrt2*Zh33*cb*la1s*sb*xvev
     &      -2D0*Sqrt2*Zh33*cb*la2s*sb*xvev-12D0*MW/EE*SW
     &      *Zh31*cb*la6*sb**2+4D0*MW/EE*SW*Zh31*cb*la7*sb
     &      **2-4D0*MW/EE*SW*Zh32*cb*la2*sb**2+4D0*MW/EE
     &      *SW*Zh32*cb*la3*sb**2+4D0*MW/EE*SW*Zh32*cb*la4
     &      *sb**2+4D0*MW/EE*SW*Zh32*cb*la5*sb**2
      AAABR(39) = EE*MW/SW*Sqrt2*(1D0-2D0*sb**2)
      AAABR(40) = EE/SW*Sqrt2*(MW*Zl11-Ml**2/MW*Zl11-2D0
     &      *MW*Zl11*sb**2-Al/MW*Ml*Zl12+Lambda/MW*Ml*Zl12
     &      /cb*sb*xvev)
      AAABR(41) = EE/MW*Ml/SW*(Al*Zl11*Zl22-Al*Zl12*Zl21
     &      -Lambda*Zl11*Zl22/cb*sb*xvev+Lambda*Zl12*Zl21
     &      /cb*sb*xvev)
      AAABR(42) = Ml*(Lambda*Pa12*Sqrt2*Zl11*Zl22/cb*sb-Lambda
     &      *Pa12*Sqrt2*Zl12*Zl21/cb*sb+EE*Lambda/MW*Pa11
     &      /SW*Zl11*Zl22*xvev-EE*Lambda/MW*Pa11/SW*Zl12
     &      *Zl21*xvev+Al*EE/MW*Pa11/SW*Zl11*Zl22/cb*sb-Al
     &      *EE/MW*Pa11/SW*Zl12*Zl21/cb*sb)
      AAABR(43) = Ml*(Lambda*Pa22*Sqrt2*Zl11*Zl22/cb*sb-Lambda
     &      *Pa22*Sqrt2*Zl12*Zl21/cb*sb+EE*Lambda/MW*Pa21
     &      /SW*Zl11*Zl22*xvev-EE*Lambda/MW*Pa21/SW*Zl12
     &      *Zl21*xvev+Al*EE/MW*Pa21/SW*Zl11*Zl22/cb*sb-Al
     &      *EE/MW*Pa21/SW*Zl12*Zl21/cb*sb)
      AAABR(44) = EE*MW*(Zh21/SW*cb-Zh22/SW*sb-SW/CW**2*Zh21
     &      *cb+SW/CW**2*Zh22*sb)
      AAABR(45) = 2D0*EE/MW*Ml**2/SW*Zh11/cb-EE*MW/SW*Zh11
     &      *Zl11**2/cb+EE*MW/SW*Zh12*Zl11**2*sb+EE/CW**2
     &      *MW*SW*Zh11*Zl11**2/cb-2D0*EE/CW**2*MW*SW*Zh11
     &      *Zl12**2/cb-EE/CW**2*MW*SW*Zh12*Zl11**2*sb+2D0
     &      *EE/CW**2*MW*SW*Zh12*Zl12**2*sb+EE*MW/SW*Zh11
     &      *Zl11**2/cb*sb**2-EE/CW**2*MW*SW*Zh11*Zl11**2
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zh11*Zl12**2/cb
     &      *sb**2-2D0*Lambda*Ml*Sqrt2*Zh13*Zl11*Zl12/cb
     &      *sb+2D0*Al*EE/MW*Ml/SW*Zh11*Zl11*Zl12/cb-2D0
     &      *EE*Lambda/MW*Ml/SW*Zh12*Zl11*Zl12/cb*xvev
      AAABR(46) = EE*MW*(Zh11/SW*cb-Zh12/SW*sb-SW/CW**2*Zh11
     &      *cb+SW/CW**2*Zh12*sb)
      AAABR(47) = 2D0*EE/MW*Ml**2/SW*Zh21/cb-EE*MW/SW*Zh21
     &      *Zl11**2/cb+EE*MW/SW*Zh22*Zl11**2*sb+EE/CW**2
     &      *MW*SW*Zh21*Zl11**2/cb-2D0*EE/CW**2*MW*SW*Zh21
     &      *Zl12**2/cb-EE/CW**2*MW*SW*Zh22*Zl11**2*sb+2D0
     &      *EE/CW**2*MW*SW*Zh22*Zl12**2*sb+EE*MW/SW*Zh21
     &      *Zl11**2/cb*sb**2-EE/CW**2*MW*SW*Zh21*Zl11**2
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zh21*Zl12**2/cb
     &      *sb**2-2D0*Lambda*Ml*Sqrt2*Zh23*Zl11*Zl12/cb
     &      *sb+2D0*Al*EE/MW*Ml/SW*Zh21*Zl11*Zl12/cb-2D0
     &      *EE*Lambda/MW*Ml/SW*Zh22*Zl11*Zl12/cb*xvev
      AAABR(48) = EE*MW*(Zh31/SW*cb-Zh32/SW*sb-SW/CW**2*Zh31
     &      *cb+SW/CW**2*Zh32*sb)
      AAABR(49) = 2D0*EE/MW*Ml**2/SW*Zh31/cb-EE*MW/SW*Zh31
     &      *Zl11**2/cb+EE*MW/SW*Zh32*Zl11**2*sb+EE/CW**2
     &      *MW*SW*Zh31*Zl11**2/cb-2D0*EE/CW**2*MW*SW*Zh31
     &      *Zl12**2/cb-EE/CW**2*MW*SW*Zh32*Zl11**2*sb+2D0
     &      *EE/CW**2*MW*SW*Zh32*Zl12**2*sb+EE*MW/SW*Zh31
     &      *Zl11**2/cb*sb**2-EE/CW**2*MW*SW*Zh31*Zl11**2
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zh31*Zl12**2/cb
     &      *sb**2-2D0*Lambda*Ml*Sqrt2*Zh33*Zl11*Zl12/cb
     &      *sb+2D0*Al*EE/MW*Ml/SW*Zh31*Zl11*Zl12/cb-2D0
     &      *EE*Lambda/MW*Ml/SW*Zh32*Zl11*Zl12/cb*xvev
      AAABR(50) = EE*MW/SW*Zh11*Zl11*Zl21/cb-EE*MW/SW*Zh12
     &      *Zl11*Zl21*sb-EE/CW**2*MW*SW*Zh11*Zl11*Zl21/cb
     &      +2D0*EE/CW**2*MW*SW*Zh11*Zl12*Zl22/cb+EE/CW**2
     &      *MW*SW*Zh12*Zl11*Zl21*sb-2D0*EE/CW**2*MW*SW*Zh12
     &      *Zl12*Zl22*sb-EE*MW/SW*Zh11*Zl11*Zl21/cb*sb**2
     &      +Lambda*Ml*Sqrt2*Zh13*Zl11*Zl22/cb*sb+Lambda
     &      *Ml*Sqrt2*Zh13*Zl12*Zl21/cb*sb-Al*EE/MW*Ml/SW
     &      *Zh11*Zl11*Zl22/cb-Al*EE/MW*Ml/SW*Zh11*Zl12*Zl21
     &      /cb+EE/CW**2*MW*SW*Zh11*Zl11*Zl21/cb*sb**2-2D0
     &      *EE/CW**2*MW*SW*Zh11*Zl12*Zl22/cb*sb**2+EE*Lambda
     &      /MW*Ml/SW*Zh12*Zl11*Zl22/cb*xvev+EE*Lambda/MW
     &      *Ml/SW*Zh12*Zl12*Zl21/cb*xvev
      AAABR(51) = EE*MW/SW*Zh31*Zl11*Zl21/cb-EE*MW/SW*Zh32
     &      *Zl11*Zl21*sb-EE/CW**2*MW*SW*Zh31*Zl11*Zl21/cb
     &      +2D0*EE/CW**2*MW*SW*Zh31*Zl12*Zl22/cb+EE/CW**2
     &      *MW*SW*Zh32*Zl11*Zl21*sb-2D0*EE/CW**2*MW*SW*Zh32
     &      *Zl12*Zl22*sb-EE*MW/SW*Zh31*Zl11*Zl21/cb*sb**2
     &      +Lambda*Ml*Sqrt2*Zh33*Zl11*Zl22/cb*sb+Lambda
     &      *Ml*Sqrt2*Zh33*Zl12*Zl21/cb*sb-Al*EE/MW*Ml/SW
     &      *Zh31*Zl11*Zl22/cb-Al*EE/MW*Ml/SW*Zh31*Zl12*Zl21
     &      /cb+EE/CW**2*MW*SW*Zh31*Zl11*Zl21/cb*sb**2-2D0
     &      *EE/CW**2*MW*SW*Zh31*Zl12*Zl22/cb*sb**2+EE*Lambda
     &      /MW*Ml/SW*Zh32*Zl11*Zl22/cb*xvev+EE*Lambda/MW
     &      *Ml/SW*Zh32*Zl12*Zl21/cb*xvev
      AAABR(52) = EE*MW/SW*Zh21*Zl11*Zl21/cb-EE*MW/SW*Zh22
     &      *Zl11*Zl21*sb-EE/CW**2*MW*SW*Zh21*Zl11*Zl21/cb
     &      +2D0*EE/CW**2*MW*SW*Zh21*Zl12*Zl22/cb+EE/CW**2
     &      *MW*SW*Zh22*Zl11*Zl21*sb-2D0*EE/CW**2*MW*SW*Zh22
     &      *Zl12*Zl22*sb-EE*MW/SW*Zh21*Zl11*Zl21/cb*sb**2
     &      +Lambda*Ml*Sqrt2*Zh23*Zl11*Zl22/cb*sb+Lambda
     &      *Ml*Sqrt2*Zh23*Zl12*Zl21/cb*sb-Al*EE/MW*Ml/SW
     &      *Zh21*Zl11*Zl22/cb-Al*EE/MW*Ml/SW*Zh21*Zl12*Zl21
     &      /cb+EE/CW**2*MW*SW*Zh21*Zl11*Zl21/cb*sb**2-2D0
     &      *EE/CW**2*MW*SW*Zh21*Zl12*Zl22/cb*sb**2+EE*Lambda
     &      /MW*Ml/SW*Zh22*Zl11*Zl22/cb*xvev+EE*Lambda/MW
     &      *Ml/SW*Zh22*Zl12*Zl21/cb*xvev
      AAABR(53) = EE/SW*Sqrt2*(MW*Zl21-Ml**2/MW*Zl21-2D0
     &      *MW*Zl21*sb**2-Al/MW*Ml*Zl22+Lambda/MW*Ml*Zl22
     &      /cb*sb*xvev)
      AAABR(54) = EE/CW**2*MW*SW*(Zh31*cb-Zh32*sb)
      AAABR(55) = 2D0*EE/MW*Ml**2/SW*Zh21/cb-EE*MW/SW*Zh21
     &      *Zl21**2/cb+EE*MW/SW*Zh22*Zl21**2*sb+EE/CW**2
     &      *MW*SW*Zh21*Zl21**2/cb-2D0*EE/CW**2*MW*SW*Zh21
     &      *Zl22**2/cb-EE/CW**2*MW*SW*Zh22*Zl21**2*sb+2D0
     &      *EE/CW**2*MW*SW*Zh22*Zl22**2*sb+EE*MW/SW*Zh21
     &      *Zl21**2/cb*sb**2-EE/CW**2*MW*SW*Zh21*Zl21**2
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zh21*Zl22**2/cb
     &      *sb**2-2D0*Lambda*Ml*Sqrt2*Zh23*Zl21*Zl22/cb
     &      *sb+2D0*Al*EE/MW*Ml/SW*Zh21*Zl21*Zl22/cb-2D0
     &      *EE*Lambda/MW*Ml/SW*Zh22*Zl21*Zl22/cb*xvev
      AAABR(56) = EE/CW**2*MW*SW*(Zh11*cb-Zh12*sb)
      AAABR(57) = 2D0*EE/MW*Ml**2/SW*Zh31/cb-EE*MW/SW*Zh31
     &      *Zl21**2/cb+EE*MW/SW*Zh32*Zl21**2*sb+EE/CW**2
     &      *MW*SW*Zh31*Zl21**2/cb-2D0*EE/CW**2*MW*SW*Zh31
     &      *Zl22**2/cb-EE/CW**2*MW*SW*Zh32*Zl21**2*sb+2D0
     &      *EE/CW**2*MW*SW*Zh32*Zl22**2*sb+EE*MW/SW*Zh31
     &      *Zl21**2/cb*sb**2-EE/CW**2*MW*SW*Zh31*Zl21**2
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zh31*Zl22**2/cb
     &      *sb**2-2D0*Lambda*Ml*Sqrt2*Zh33*Zl21*Zl22/cb
     &      *sb+2D0*Al*EE/MW*Ml/SW*Zh31*Zl21*Zl22/cb-2D0
     &      *EE*Lambda/MW*Ml/SW*Zh32*Zl21*Zl22/cb*xvev
      AAABR(58) = EE/CW**2*MW*SW*(Zh21*cb-Zh22*sb)
      AAABR(59) = 2D0*EE/MW*Ml**2/SW*Zh11/cb-EE*MW/SW*Zh11
     &      *Zl21**2/cb+EE*MW/SW*Zh12*Zl21**2*sb+EE/CW**2
     &      *MW*SW*Zh11*Zl21**2/cb-2D0*EE/CW**2*MW*SW*Zh11
     &      *Zl22**2/cb-EE/CW**2*MW*SW*Zh12*Zl21**2*sb+2D0
     &      *EE/CW**2*MW*SW*Zh12*Zl22**2*sb+EE*MW/SW*Zh11
     &      *Zl21**2/cb*sb**2-EE/CW**2*MW*SW*Zh11*Zl21**2
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zh11*Zl22**2/cb
     &      *sb**2-2D0*Lambda*Ml*Sqrt2*Zh13*Zl21*Zl22/cb
     &      *sb+2D0*Al*EE/MW*Ml/SW*Zh11*Zl21*Zl22/cb-2D0
     &      *EE*Lambda/MW*Ml/SW*Zh12*Zl21*Zl22/cb*xvev
      AAABR(60) = EE*MW*(Zh31/SW*cb-Zh32/SW*sb+SW/CW**2*Zh31
     &      *cb-SW/CW**2*Zh32*sb)
      AAABR(61) = EE*MW*(Zh11/SW*cb-Zh12/SW*sb+SW/CW**2*Zh11
     &      *cb-SW/CW**2*Zh12*sb)
      AAABR(62) = EE*MW*(Zh21/SW*cb-Zh22/SW*sb+SW/CW**2*Zh21
     &      *cb-SW/CW**2*Zh22*sb)
      AAABR(63) = EE/SW*Sqrt2*(MW*Zcc11*Zss11+Mc**2/MW*Zcc11
     &      *Zss11-Mq**2/MW*Zcc11*Zss11-2D0*MW*Zcc11*Zss11
     &      *sb**2-Ad/MW*Mq*Zcc11*Zss12+Au/MW*Mc*Zcc12*Zss11
     &      -Lambda/MW*Mc*Zcc12*Zss11/cb/sb*xvev+Lambda/MW
     &      *Mc*Zcc12*Zss11/cb*sb*xvev+Lambda/MW*Mq*Zcc11
     &      *Zss12/cb*sb*xvev)
      AAABR(64) = EE/SW*Sqrt2*(MW*Zdd11*Zuu11-2D0*MW*Zdd11
     &      *Zuu11*sb**2-Ad/MW*Mq*Zdd12*Zuu11+Au/MW*Mq*Zdd11
     &      *Zuu12-Lambda/MW*Mq*Zdd11*Zuu12/cb/sb*xvev+Lambda
     &      /MW*Mq*Zdd11*Zuu12/cb*sb*xvev+Lambda/MW*Mq*Zdd12
     &      *Zuu11/cb*sb*xvev)
      AAABR(65) = EE/SW*Sqrt2*(MW*Zb11*Zt11-MbMM**2/MW*Zb11
     &      *Zt11+MtMM**2/MW*Zb11*Zt11-2D0*MW*Zb11*Zt11*sb
     &      **2-Ab/MW*MbMM*Zb12*Zt11+At/MW*MtMM*Zb11*Zt12
     &      +Lambda/MW*MbMM*Zb12*Zt11/cb*sb*xvev-Lambda/MW
     &      *MtMM*Zb11*Zt12/cb/sb*xvev+Lambda/MW*MtMM*Zb11
     &      *Zt12/cb*sb*xvev)
      AAABR(66) = EE/SW*Sqrt2*(MW*Zb11*Zt21-MbMM**2/MW*Zb11
     &      *Zt21+MtMM**2/MW*Zb11*Zt21-2D0*MW*Zb11*Zt21*sb
     &      **2-Ab/MW*MbMM*Zb12*Zt21+At/MW*MtMM*Zb11*Zt22
     &      +Lambda/MW*MbMM*Zb12*Zt21/cb*sb*xvev-Lambda/MW
     &      *MtMM*Zb11*Zt22/cb/sb*xvev+Lambda/MW*MtMM*Zb11
     &      *Zt22/cb*sb*xvev)
      AAABR(67) = EE/SW*Sqrt2*(MW*Zcc21*Zss11+Mc**2/MW*Zcc21
     &      *Zss11-Mq**2/MW*Zcc21*Zss11-2D0*MW*Zcc21*Zss11
     &      *sb**2-Ad/MW*Mq*Zcc21*Zss12+Au/MW*Mc*Zcc22*Zss11
     &      -Lambda/MW*Mc*Zcc22*Zss11/cb/sb*xvev+Lambda/MW
     &      *Mc*Zcc22*Zss11/cb*sb*xvev+Lambda/MW*Mq*Zcc21
     &      *Zss12/cb*sb*xvev)
      AAABR(68) = EE/SW*Sqrt2*(MW*Zdd11*Zuu21-2D0*MW*Zdd11
     &      *Zuu21*sb**2-Ad/MW*Mq*Zdd12*Zuu21+Au/MW*Mq*Zdd11
     &      *Zuu22-Lambda/MW*Mq*Zdd11*Zuu22/cb/sb*xvev+Lambda
     &      /MW*Mq*Zdd11*Zuu22/cb*sb*xvev+Lambda/MW*Mq*Zdd12
     &      *Zuu21/cb*sb*xvev)
      AAABR(69) = EE/MW*Mq/SW*(Ad*Zss11*Zss22-Ad*Zss12*Zss21
     &      -Lambda*Zss11*Zss22/cb*sb*xvev+Lambda*Zss12*Zss21
     &      /cb*sb*xvev)
      AAABR(70) = EE/MW*MbMM/SW*(Ab*Zb11*Zb22-Ab*Zb12*Zb21
     &      -Lambda*Zb11*Zb22/cb*sb*xvev+Lambda*Zb12*Zb21
     &      /cb*sb*xvev)
      AAABR(71) = EE/MW*Mq/SW*(Ad*Zdd11*Zdd22-Ad*Zdd12*Zdd21
     &      -Lambda*Zdd11*Zdd22/cb*sb*xvev+Lambda*Zdd12*Zdd21
     &      /cb*sb*xvev)
      AAABR(72) = Mq*(Lambda*Pa22*Sqrt2*Zss11*Zss22/cb*sb
     &      -Lambda*Pa22*Sqrt2*Zss12*Zss21/cb*sb+EE*Lambda
     &      /MW*Pa21/SW*Zss11*Zss22*xvev-EE*Lambda/MW*Pa21
     &      /SW*Zss12*Zss21*xvev+Ad*EE/MW*Pa21/SW*Zss11*Zss22
     &      /cb*sb-Ad*EE/MW*Pa21/SW*Zss12*Zss21/cb*sb)
      AAABR(73) = Mq*(Lambda*Pa12*Sqrt2*Zdd11*Zdd22/cb*sb
     &      -Lambda*Pa12*Sqrt2*Zdd12*Zdd21/cb*sb+EE*Lambda
     &      /MW*Pa11/SW*Zdd11*Zdd22*xvev-EE*Lambda/MW*Pa11
     &      /SW*Zdd12*Zdd21*xvev+Ad*EE/MW*Pa11/SW*Zdd11*Zdd22
     &      /cb*sb-Ad*EE/MW*Pa11/SW*Zdd12*Zdd21/cb*sb)
      AAABR(74) = MbMM*(Lambda*Pa22*Sqrt2*Zb11*Zb22/cb*sb
     &      -Lambda*Pa22*Sqrt2*Zb12*Zb21/cb*sb+EE*Lambda
     &      /MW*Pa21/SW*Zb11*Zb22*xvev-EE*Lambda/MW*Pa21
     &      /SW*Zb12*Zb21*xvev+Ab*EE/MW*Pa21/SW*Zb11*Zb22
     &      /cb*sb-Ab*EE/MW*Pa21/SW*Zb12*Zb21/cb*sb)
      AAABR(75) = Mq*(Lambda*Pa12*Sqrt2*Zss11*Zss22/cb*sb
     &      -Lambda*Pa12*Sqrt2*Zss12*Zss21/cb*sb+EE*Lambda
     &      /MW*Pa11/SW*Zss11*Zss22*xvev-EE*Lambda/MW*Pa11
     &      /SW*Zss12*Zss21*xvev+Ad*EE/MW*Pa11/SW*Zss11*Zss22
     &      /cb*sb-Ad*EE/MW*Pa11/SW*Zss12*Zss21/cb*sb)
      AAABR(76) = MbMM*(Lambda*Pa12*Sqrt2*Zb11*Zb22/cb*sb
     &      -Lambda*Pa12*Sqrt2*Zb12*Zb21/cb*sb+EE*Lambda
     &      /MW*Pa11/SW*Zb11*Zb22*xvev-EE*Lambda/MW*Pa11
     &      /SW*Zb12*Zb21*xvev+Ab*EE/MW*Pa11/SW*Zb11*Zb22
     &      /cb*sb-Ab*EE/MW*Pa11/SW*Zb12*Zb21/cb*sb)
      AAABR(77) = Mq*(Lambda*Pa22*Sqrt2*Zdd11*Zdd22/cb*sb
     &      -Lambda*Pa22*Sqrt2*Zdd12*Zdd21/cb*sb+EE*Lambda
     &      /MW*Pa21/SW*Zdd11*Zdd22*xvev-EE*Lambda/MW*Pa21
     &      /SW*Zdd12*Zdd21*xvev+Ad*EE/MW*Pa21/SW*Zdd11*Zdd22
     &      /cb*sb-Ad*EE/MW*Pa21/SW*Zdd12*Zdd21/cb*sb)
      AAABR(78) = 6D0*EE/MW*Mq**2/SW*Zh31/cb-3D0*EE*MW/SW
     &      *Zh31*Zss11**2/cb+3D0*EE*MW/SW*Zh32*Zss11**2
     &      *sb-EE/CW**2*MW*SW*Zh31*Zss11**2/cb-2D0*EE/CW
     &      **2*MW*SW*Zh31*Zss12**2/cb+EE/CW**2*MW*SW*Zh32
     &      *Zss11**2*sb+2D0*EE/CW**2*MW*SW*Zh32*Zss12**2
     &      *sb+3D0*EE*MW/SW*Zh31*Zss11**2/cb*sb**2+EE/CW
     &      **2*MW*SW*Zh31*Zss11**2/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zh31*Zss12**2/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zh33*Zss11*Zss12/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zh31
     &      *Zss11*Zss12/cb-6D0*EE*Lambda/MW*Mq/SW*Zh32*Zss11
     &      *Zss12/cb*xvev
      AAABR(79) = 6D0*EE/MW*Mq**2/SW*Zh11/cb-3D0*EE*MW/SW
     &      *Zh11*Zss11**2/cb+3D0*EE*MW/SW*Zh12*Zss11**2
     &      *sb-EE/CW**2*MW*SW*Zh11*Zss11**2/cb-2D0*EE/CW
     &      **2*MW*SW*Zh11*Zss12**2/cb+EE/CW**2*MW*SW*Zh12
     &      *Zss11**2*sb+2D0*EE/CW**2*MW*SW*Zh12*Zss12**2
     &      *sb+3D0*EE*MW/SW*Zh11*Zss11**2/cb*sb**2+EE/CW
     &      **2*MW*SW*Zh11*Zss11**2/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zh11*Zss12**2/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zh13*Zss11*Zss12/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zh11
     &      *Zss11*Zss12/cb-6D0*EE*Lambda/MW*Mq/SW*Zh12*Zss11
     &      *Zss12/cb*xvev
      AAABR(80) = 6D0*EE/MW*Mq**2/SW*Zh11/cb-3D0*EE*MW/SW
     &      *Zdd11**2*Zh11/cb+3D0*EE*MW/SW*Zdd11**2*Zh12
     &      *sb-EE/CW**2*MW*SW*Zdd11**2*Zh11/cb+EE/CW**2
     &      *MW*SW*Zdd11**2*Zh12*sb-2D0*EE/CW**2*MW*SW*Zdd12
     &      **2*Zh11/cb+2D0*EE/CW**2*MW*SW*Zdd12**2*Zh12
     &      *sb+3D0*EE*MW/SW*Zdd11**2*Zh11/cb*sb**2+EE/CW
     &      **2*MW*SW*Zdd11**2*Zh11/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zdd12**2*Zh11/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zdd11*Zdd12*Zh13/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zdd11
     &      *Zdd12*Zh11/cb-6D0*EE*Lambda/MW*Mq/SW*Zdd11*Zdd12
     &      *Zh12/cb*xvev
      AAABR(81) = 6D0*EE/MW*Mq**2/SW*Zh21/cb-3D0*EE*MW/SW
     &      *Zdd11**2*Zh21/cb+3D0*EE*MW/SW*Zdd11**2*Zh22
     &      *sb-EE/CW**2*MW*SW*Zdd11**2*Zh21/cb+EE/CW**2
     &      *MW*SW*Zdd11**2*Zh22*sb-2D0*EE/CW**2*MW*SW*Zdd12
     &      **2*Zh21/cb+2D0*EE/CW**2*MW*SW*Zdd12**2*Zh22
     &      *sb+3D0*EE*MW/SW*Zdd11**2*Zh21/cb*sb**2+EE/CW
     &      **2*MW*SW*Zdd11**2*Zh21/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zdd12**2*Zh21/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zdd11*Zdd12*Zh23/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zdd11
     &      *Zdd12*Zh21/cb-6D0*EE*Lambda/MW*Mq/SW*Zdd11*Zdd12
     &      *Zh22/cb*xvev
      AAABR(82) = 6D0*EE/MW*Mq**2/SW*Zh31/cb-3D0*EE*MW/SW
     &      *Zdd11**2*Zh31/cb+3D0*EE*MW/SW*Zdd11**2*Zh32
     &      *sb-EE/CW**2*MW*SW*Zdd11**2*Zh31/cb+EE/CW**2
     &      *MW*SW*Zdd11**2*Zh32*sb-2D0*EE/CW**2*MW*SW*Zdd12
     &      **2*Zh31/cb+2D0*EE/CW**2*MW*SW*Zdd12**2*Zh32
     &      *sb+3D0*EE*MW/SW*Zdd11**2*Zh31/cb*sb**2+EE/CW
     &      **2*MW*SW*Zdd11**2*Zh31/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zdd12**2*Zh31/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zdd11*Zdd12*Zh33/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zdd11
     &      *Zdd12*Zh31/cb-6D0*EE*Lambda/MW*Mq/SW*Zdd11*Zdd12
     &      *Zh32/cb*xvev
      AAABR(83) = 6D0*EE/MW*MbMM**2/SW*Zh11/cb-3D0*EE*MW
     &      /SW*Zb11**2*Zh11/cb+3D0*EE*MW/SW*Zb11**2*Zh12
     &      *sb-EE/CW**2*MW*SW*Zb11**2*Zh11/cb+EE/CW**2*MW
     &      *SW*Zb11**2*Zh12*sb-2D0*EE/CW**2*MW*SW*Zb12**2
     &      *Zh11/cb+2D0*EE/CW**2*MW*SW*Zb12**2*Zh12*sb+3D0
     &      *EE*MW/SW*Zb11**2*Zh11/cb*sb**2+EE/CW**2*MW*SW
     &      *Zb11**2*Zh11/cb*sb**2+2D0*EE/CW**2*MW*SW*Zb12
     &      **2*Zh11/cb*sb**2-6D0*Lambda*MbMM*Sqrt2*Zb11
     &      *Zb12*Zh13/cb*sb+6D0*Ab*EE/MW*MbMM/SW*Zb11*Zb12
     &      *Zh11/cb-6D0*EE*Lambda/MW*MbMM/SW*Zb11*Zb12*Zh12
     &      /cb*xvev
      AAABR(84) = 6D0*EE/MW*Mq**2/SW*Zh21/cb-3D0*EE*MW/SW
     &      *Zh21*Zss11**2/cb+3D0*EE*MW/SW*Zh22*Zss11**2
     &      *sb-EE/CW**2*MW*SW*Zh21*Zss11**2/cb-2D0*EE/CW
     &      **2*MW*SW*Zh21*Zss12**2/cb+EE/CW**2*MW*SW*Zh22
     &      *Zss11**2*sb+2D0*EE/CW**2*MW*SW*Zh22*Zss12**2
     &      *sb+3D0*EE*MW/SW*Zh21*Zss11**2/cb*sb**2+EE/CW
     &      **2*MW*SW*Zh21*Zss11**2/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zh21*Zss12**2/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zh23*Zss11*Zss12/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zh21
     &      *Zss11*Zss12/cb-6D0*EE*Lambda/MW*Mq/SW*Zh22*Zss11
     &      *Zss12/cb*xvev
      AAABR(85) = 6D0*EE/MW*MbMM**2/SW*Zh31/cb-3D0*EE*MW
     &      /SW*Zb11**2*Zh31/cb+3D0*EE*MW/SW*Zb11**2*Zh32
     &      *sb-EE/CW**2*MW*SW*Zb11**2*Zh31/cb+EE/CW**2*MW
     &      *SW*Zb11**2*Zh32*sb-2D0*EE/CW**2*MW*SW*Zb12**2
     &      *Zh31/cb+2D0*EE/CW**2*MW*SW*Zb12**2*Zh32*sb+3D0
     &      *EE*MW/SW*Zb11**2*Zh31/cb*sb**2+EE/CW**2*MW*SW
     &      *Zb11**2*Zh31/cb*sb**2+2D0*EE/CW**2*MW*SW*Zb12
     &      **2*Zh31/cb*sb**2-6D0*Lambda*MbMM*Sqrt2*Zb11
     &      *Zb12*Zh33/cb*sb+6D0*Ab*EE/MW*MbMM/SW*Zb11*Zb12
     &      *Zh31/cb-6D0*EE*Lambda/MW*MbMM/SW*Zb11*Zb12*Zh32
     &      /cb*xvev
      AAABR(86) = 6D0*EE/MW*MbMM**2/SW*Zh21/cb-3D0*EE*MW
     &      /SW*Zb11**2*Zh21/cb+3D0*EE*MW/SW*Zb11**2*Zh22
     &      *sb-EE/CW**2*MW*SW*Zb11**2*Zh21/cb+EE/CW**2*MW
     &      *SW*Zb11**2*Zh22*sb-2D0*EE/CW**2*MW*SW*Zb12**2
     &      *Zh21/cb+2D0*EE/CW**2*MW*SW*Zb12**2*Zh22*sb+3D0
     &      *EE*MW/SW*Zb11**2*Zh21/cb*sb**2+EE/CW**2*MW*SW
     &      *Zb11**2*Zh21/cb*sb**2+2D0*EE/CW**2*MW*SW*Zb12
     &      **2*Zh21/cb*sb**2-6D0*Lambda*MbMM*Sqrt2*Zb11
     &      *Zb12*Zh23/cb*sb+6D0*Ab*EE/MW*MbMM/SW*Zb11*Zb12
     &      *Zh21/cb-6D0*EE*Lambda/MW*MbMM/SW*Zb11*Zb12*Zh22
     &      /cb*xvev
      AAABR(87) = 3D0*EE*MW/SW*Zdd11*Zdd21*Zh21/cb-3D0*EE
     &      *MW/SW*Zdd11*Zdd21*Zh22*sb+EE/CW**2*MW*SW*Zdd11
     &      *Zdd21*Zh21/cb-EE/CW**2*MW*SW*Zdd11*Zdd21*Zh22
     &      *sb+2D0*EE/CW**2*MW*SW*Zdd12*Zdd22*Zh21/cb-2D0
     &      *EE/CW**2*MW*SW*Zdd12*Zdd22*Zh22*sb-3D0*EE*MW
     &      /SW*Zdd11*Zdd21*Zh21/cb*sb**2+3D0*Lambda*Mq*Sqrt2
     &      *Zdd11*Zdd22*Zh23/cb*sb+3D0*Lambda*Mq*Sqrt2*Zdd12
     &      *Zdd21*Zh23/cb*sb-3D0*Ad*EE/MW*Mq/SW*Zdd11*Zdd22
     &      *Zh21/cb-3D0*Ad*EE/MW*Mq/SW*Zdd12*Zdd21*Zh21
     &      /cb-EE/CW**2*MW*SW*Zdd11*Zdd21*Zh21/cb*sb**2
     &      -2D0*EE/CW**2*MW*SW*Zdd12*Zdd22*Zh21/cb*sb**2
     &      +3D0*EE*Lambda/MW*Mq/SW*Zdd11*Zdd22*Zh22/cb*xvev
     &      +3D0*EE*Lambda/MW*Mq/SW*Zdd12*Zdd21*Zh22/cb*xvev
      AAABR(88) = 3D0*EE*MW/SW*Zb11*Zb21*Zh11/cb-3D0*EE*MW
     &      /SW*Zb11*Zb21*Zh12*sb+EE/CW**2*MW*SW*Zb11*Zb21
     &      *Zh11/cb-EE/CW**2*MW*SW*Zb11*Zb21*Zh12*sb+2D0
     &      *EE/CW**2*MW*SW*Zb12*Zb22*Zh11/cb-2D0*EE/CW**2
     &      *MW*SW*Zb12*Zb22*Zh12*sb-3D0*EE*MW/SW*Zb11*Zb21
     &      *Zh11/cb*sb**2+3D0*Lambda*MbMM*Sqrt2*Zb11*Zb22
     &      *Zh13/cb*sb+3D0*Lambda*MbMM*Sqrt2*Zb12*Zb21*Zh13
     &      /cb*sb-3D0*Ab*EE/MW*MbMM/SW*Zb11*Zb22*Zh11/cb
     &      -3D0*Ab*EE/MW*MbMM/SW*Zb12*Zb21*Zh11/cb-EE/CW
     &      **2*MW*SW*Zb11*Zb21*Zh11/cb*sb**2-2D0*EE/CW**2
     &      *MW*SW*Zb12*Zb22*Zh11/cb*sb**2+3D0*EE*Lambda
     &      /MW*MbMM/SW*Zb11*Zb22*Zh12/cb*xvev+3D0*EE*Lambda
     &      /MW*MbMM/SW*Zb12*Zb21*Zh12/cb*xvev
      AAABR(89) = 3D0*EE*MW/SW*Zb11*Zb21*Zh31/cb-3D0*EE*MW
     &      /SW*Zb11*Zb21*Zh32*sb+EE/CW**2*MW*SW*Zb11*Zb21
     &      *Zh31/cb-EE/CW**2*MW*SW*Zb11*Zb21*Zh32*sb+2D0
     &      *EE/CW**2*MW*SW*Zb12*Zb22*Zh31/cb-2D0*EE/CW**2
     &      *MW*SW*Zb12*Zb22*Zh32*sb-3D0*EE*MW/SW*Zb11*Zb21
     &      *Zh31/cb*sb**2+3D0*Lambda*MbMM*Sqrt2*Zb11*Zb22
     &      *Zh33/cb*sb+3D0*Lambda*MbMM*Sqrt2*Zb12*Zb21*Zh33
     &      /cb*sb-3D0*Ab*EE/MW*MbMM/SW*Zb11*Zb22*Zh31/cb
     &      -3D0*Ab*EE/MW*MbMM/SW*Zb12*Zb21*Zh31/cb-EE/CW
     &      **2*MW*SW*Zb11*Zb21*Zh31/cb*sb**2-2D0*EE/CW**2
     &      *MW*SW*Zb12*Zb22*Zh31/cb*sb**2+3D0*EE*Lambda
     &      /MW*MbMM/SW*Zb11*Zb22*Zh32/cb*xvev+3D0*EE*Lambda
     &      /MW*MbMM/SW*Zb12*Zb21*Zh32/cb*xvev
      AAABR(90) = 3D0*EE*MW/SW*Zdd11*Zdd21*Zh11/cb-3D0*EE
     &      *MW/SW*Zdd11*Zdd21*Zh12*sb+EE/CW**2*MW*SW*Zdd11
     &      *Zdd21*Zh11/cb-EE/CW**2*MW*SW*Zdd11*Zdd21*Zh12
     &      *sb+2D0*EE/CW**2*MW*SW*Zdd12*Zdd22*Zh11/cb-2D0
     &      *EE/CW**2*MW*SW*Zdd12*Zdd22*Zh12*sb-3D0*EE*MW
     &      /SW*Zdd11*Zdd21*Zh11/cb*sb**2+3D0*Lambda*Mq*Sqrt2
     &      *Zdd11*Zdd22*Zh13/cb*sb+3D0*Lambda*Mq*Sqrt2*Zdd12
     &      *Zdd21*Zh13/cb*sb-3D0*Ad*EE/MW*Mq/SW*Zdd11*Zdd22
     &      *Zh11/cb-3D0*Ad*EE/MW*Mq/SW*Zdd12*Zdd21*Zh11
     &      /cb-EE/CW**2*MW*SW*Zdd11*Zdd21*Zh11/cb*sb**2
     &      -2D0*EE/CW**2*MW*SW*Zdd12*Zdd22*Zh11/cb*sb**2
     &      +3D0*EE*Lambda/MW*Mq/SW*Zdd11*Zdd22*Zh12/cb*xvev
     &      +3D0*EE*Lambda/MW*Mq/SW*Zdd12*Zdd21*Zh12/cb*xvev
      AAABR(91) = 3D0*EE*MW/SW*Zdd11*Zdd21*Zh31/cb-3D0*EE
     &      *MW/SW*Zdd11*Zdd21*Zh32*sb+EE/CW**2*MW*SW*Zdd11
     &      *Zdd21*Zh31/cb-EE/CW**2*MW*SW*Zdd11*Zdd21*Zh32
     &      *sb+2D0*EE/CW**2*MW*SW*Zdd12*Zdd22*Zh31/cb-2D0
     &      *EE/CW**2*MW*SW*Zdd12*Zdd22*Zh32*sb-3D0*EE*MW
     &      /SW*Zdd11*Zdd21*Zh31/cb*sb**2+3D0*Lambda*Mq*Sqrt2
     &      *Zdd11*Zdd22*Zh33/cb*sb+3D0*Lambda*Mq*Sqrt2*Zdd12
     &      *Zdd21*Zh33/cb*sb-3D0*Ad*EE/MW*Mq/SW*Zdd11*Zdd22
     &      *Zh31/cb-3D0*Ad*EE/MW*Mq/SW*Zdd12*Zdd21*Zh31
     &      /cb-EE/CW**2*MW*SW*Zdd11*Zdd21*Zh31/cb*sb**2
     &      -2D0*EE/CW**2*MW*SW*Zdd12*Zdd22*Zh31/cb*sb**2
     &      +3D0*EE*Lambda/MW*Mq/SW*Zdd11*Zdd22*Zh32/cb*xvev
     &      +3D0*EE*Lambda/MW*Mq/SW*Zdd12*Zdd21*Zh32/cb*xvev
      AAABR(92) = 3D0*EE*MW/SW*Zb11*Zb21*Zh21/cb-3D0*EE*MW
     &      /SW*Zb11*Zb21*Zh22*sb+EE/CW**2*MW*SW*Zb11*Zb21
     &      *Zh21/cb-EE/CW**2*MW*SW*Zb11*Zb21*Zh22*sb+2D0
     &      *EE/CW**2*MW*SW*Zb12*Zb22*Zh21/cb-2D0*EE/CW**2
     &      *MW*SW*Zb12*Zb22*Zh22*sb-3D0*EE*MW/SW*Zb11*Zb21
     &      *Zh21/cb*sb**2+3D0*Lambda*MbMM*Sqrt2*Zb11*Zb22
     &      *Zh23/cb*sb+3D0*Lambda*MbMM*Sqrt2*Zb12*Zb21*Zh23
     &      /cb*sb-3D0*Ab*EE/MW*MbMM/SW*Zb11*Zb22*Zh21/cb
     &      -3D0*Ab*EE/MW*MbMM/SW*Zb12*Zb21*Zh21/cb-EE/CW
     &      **2*MW*SW*Zb11*Zb21*Zh21/cb*sb**2-2D0*EE/CW**2
     &      *MW*SW*Zb12*Zb22*Zh21/cb*sb**2+3D0*EE*Lambda
     &      /MW*MbMM/SW*Zb11*Zb22*Zh22/cb*xvev+3D0*EE*Lambda
     &      /MW*MbMM/SW*Zb12*Zb21*Zh22/cb*xvev
      AAABR(93) = 3D0*EE*MW/SW*Zh31*Zss11*Zss21/cb-3D0*EE
     &      *MW/SW*Zh32*Zss11*Zss21*sb+EE/CW**2*MW*SW*Zh31
     &      *Zss11*Zss21/cb+2D0*EE/CW**2*MW*SW*Zh31*Zss12
     &      *Zss22/cb-EE/CW**2*MW*SW*Zh32*Zss11*Zss21*sb
     &      -2D0*EE/CW**2*MW*SW*Zh32*Zss12*Zss22*sb-3D0*EE
     &      *MW/SW*Zh31*Zss11*Zss21/cb*sb**2+3D0*Lambda*Mq
     &      *Sqrt2*Zh33*Zss11*Zss22/cb*sb+3D0*Lambda*Mq*Sqrt2
     &      *Zh33*Zss12*Zss21/cb*sb-3D0*Ad*EE/MW*Mq/SW*Zh31
     &      *Zss11*Zss22/cb-3D0*Ad*EE/MW*Mq/SW*Zh31*Zss12
     &      *Zss21/cb-EE/CW**2*MW*SW*Zh31*Zss11*Zss21/cb
     &      *sb**2-2D0*EE/CW**2*MW*SW*Zh31*Zss12*Zss22/cb
     &      *sb**2+3D0*EE*Lambda/MW*Mq/SW*Zh32*Zss11*Zss22
     &      /cb*xvev+3D0*EE*Lambda/MW*Mq/SW*Zh32*Zss12*Zss21
     &      /cb*xvev
      AAABR(94) = 3D0*EE*MW/SW*Zh11*Zss11*Zss21/cb-3D0*EE
     &      *MW/SW*Zh12*Zss11*Zss21*sb+EE/CW**2*MW*SW*Zh11
     &      *Zss11*Zss21/cb+2D0*EE/CW**2*MW*SW*Zh11*Zss12
     &      *Zss22/cb-EE/CW**2*MW*SW*Zh12*Zss11*Zss21*sb
     &      -2D0*EE/CW**2*MW*SW*Zh12*Zss12*Zss22*sb-3D0*EE
     &      *MW/SW*Zh11*Zss11*Zss21/cb*sb**2+3D0*Lambda*Mq
     &      *Sqrt2*Zh13*Zss11*Zss22/cb*sb+3D0*Lambda*Mq*Sqrt2
     &      *Zh13*Zss12*Zss21/cb*sb-3D0*Ad*EE/MW*Mq/SW*Zh11
     &      *Zss11*Zss22/cb-3D0*Ad*EE/MW*Mq/SW*Zh11*Zss12
     &      *Zss21/cb-EE/CW**2*MW*SW*Zh11*Zss11*Zss21/cb
     &      *sb**2-2D0*EE/CW**2*MW*SW*Zh11*Zss12*Zss22/cb
     &      *sb**2+3D0*EE*Lambda/MW*Mq/SW*Zh12*Zss11*Zss22
     &      /cb*xvev+3D0*EE*Lambda/MW*Mq/SW*Zh12*Zss12*Zss21
     &      /cb*xvev
      AAABR(95) = 3D0*EE*MW/SW*Zh21*Zss11*Zss21/cb-3D0*EE
     &      *MW/SW*Zh22*Zss11*Zss21*sb+EE/CW**2*MW*SW*Zh21
     &      *Zss11*Zss21/cb+2D0*EE/CW**2*MW*SW*Zh21*Zss12
     &      *Zss22/cb-EE/CW**2*MW*SW*Zh22*Zss11*Zss21*sb
     &      -2D0*EE/CW**2*MW*SW*Zh22*Zss12*Zss22*sb-3D0*EE
     &      *MW/SW*Zh21*Zss11*Zss21/cb*sb**2+3D0*Lambda*Mq
     &      *Sqrt2*Zh23*Zss11*Zss22/cb*sb+3D0*Lambda*Mq*Sqrt2
     &      *Zh23*Zss12*Zss21/cb*sb-3D0*Ad*EE/MW*Mq/SW*Zh21
     &      *Zss11*Zss22/cb-3D0*Ad*EE/MW*Mq/SW*Zh21*Zss12
     &      *Zss21/cb-EE/CW**2*MW*SW*Zh21*Zss11*Zss21/cb
     &      *sb**2-2D0*EE/CW**2*MW*SW*Zh21*Zss12*Zss22/cb
     &      *sb**2+3D0*EE*Lambda/MW*Mq/SW*Zh22*Zss11*Zss22
     &      /cb*xvev+3D0*EE*Lambda/MW*Mq/SW*Zh22*Zss12*Zss21
     &      /cb*xvev
      AAABR(96) = EE/SW*Sqrt2*(MW*Zcc11*Zss21+Mc**2/MW*Zcc11
     &      *Zss21-Mq**2/MW*Zcc11*Zss21-2D0*MW*Zcc11*Zss21
     &      *sb**2-Ad/MW*Mq*Zcc11*Zss22+Au/MW*Mc*Zcc12*Zss21
     &      -Lambda/MW*Mc*Zcc12*Zss21/cb/sb*xvev+Lambda/MW
     &      *Mc*Zcc12*Zss21/cb*sb*xvev+Lambda/MW*Mq*Zcc11
     &      *Zss22/cb*sb*xvev)
      AAABR(97) = EE/SW*Sqrt2*(MW*Zdd21*Zuu11-2D0*MW*Zdd21
     &      *Zuu11*sb**2-Ad/MW*Mq*Zdd22*Zuu11+Au/MW*Mq*Zdd21
     &      *Zuu12-Lambda/MW*Mq*Zdd21*Zuu12/cb/sb*xvev+Lambda
     &      /MW*Mq*Zdd21*Zuu12/cb*sb*xvev+Lambda/MW*Mq*Zdd22
     &      *Zuu11/cb*sb*xvev)
      AAABR(98) = EE/SW*Sqrt2*(MW*Zb21*Zt11-MbMM**2/MW*Zb21
     &      *Zt11+MtMM**2/MW*Zb21*Zt11-2D0*MW*Zb21*Zt11*sb
     &      **2-Ab/MW*MbMM*Zb22*Zt11+At/MW*MtMM*Zb21*Zt12
     &      +Lambda/MW*MbMM*Zb22*Zt11/cb*sb*xvev-Lambda/MW
     &      *MtMM*Zb21*Zt12/cb/sb*xvev+Lambda/MW*MtMM*Zb21
     &      *Zt12/cb*sb*xvev)
      AAABR(99) = EE/SW*Sqrt2*(MW*Zb21*Zt21-MbMM**2/MW*Zb21
     &      *Zt21+MtMM**2/MW*Zb21*Zt21-2D0*MW*Zb21*Zt21*sb
     &      **2-Ab/MW*MbMM*Zb22*Zt21+At/MW*MtMM*Zb21*Zt22
     &      +Lambda/MW*MbMM*Zb22*Zt21/cb*sb*xvev-Lambda/MW
     &      *MtMM*Zb21*Zt22/cb/sb*xvev+Lambda/MW*MtMM*Zb21
     &      *Zt22/cb*sb*xvev)
      AAABR(100) = EE/SW*Sqrt2*(MW*Zdd21*Zuu21-2D0*MW*Zdd21
     &      *Zuu21*sb**2-Ad/MW*Mq*Zdd22*Zuu21+Au/MW*Mq*Zdd21
     &      *Zuu22-Lambda/MW*Mq*Zdd21*Zuu22/cb/sb*xvev+Lambda
     &      /MW*Mq*Zdd21*Zuu22/cb*sb*xvev+Lambda/MW*Mq*Zdd22
     &      *Zuu21/cb*sb*xvev)
      AAABR(101) = EE/SW*Sqrt2*(MW*Zcc21*Zss21+Mc**2/MW*Zcc21
     &      *Zss21-Mq**2/MW*Zcc21*Zss21-2D0*MW*Zcc21*Zss21
     &      *sb**2-Ad/MW*Mq*Zcc21*Zss22+Au/MW*Mc*Zcc22*Zss21
     &      -Lambda/MW*Mc*Zcc22*Zss21/cb/sb*xvev+Lambda/MW
     &      *Mc*Zcc22*Zss21/cb*sb*xvev+Lambda/MW*Mq*Zcc21
     &      *Zss22/cb*sb*xvev)
      AAABR(102) = 6D0*EE/MW*Mq**2/SW*Zh11/cb-3D0*EE*MW/SW
     &      *Zh11*Zss21**2/cb+3D0*EE*MW/SW*Zh12*Zss21**2
     &      *sb-EE/CW**2*MW*SW*Zh11*Zss21**2/cb-2D0*EE/CW
     &      **2*MW*SW*Zh11*Zss22**2/cb+EE/CW**2*MW*SW*Zh12
     &      *Zss21**2*sb+2D0*EE/CW**2*MW*SW*Zh12*Zss22**2
     &      *sb+3D0*EE*MW/SW*Zh11*Zss21**2/cb*sb**2+EE/CW
     &      **2*MW*SW*Zh11*Zss21**2/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zh11*Zss22**2/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zh13*Zss21*Zss22/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zh11
     &      *Zss21*Zss22/cb-6D0*EE*Lambda/MW*Mq/SW*Zh12*Zss21
     &      *Zss22/cb*xvev
      AAABR(103) = 6D0*EE/MW*MbMM**2/SW*Zh31/cb-3D0*EE*MW
     &      /SW*Zb21**2*Zh31/cb+3D0*EE*MW/SW*Zb21**2*Zh32
     &      *sb-EE/CW**2*MW*SW*Zb21**2*Zh31/cb+EE/CW**2*MW
     &      *SW*Zb21**2*Zh32*sb-2D0*EE/CW**2*MW*SW*Zb22**2
     &      *Zh31/cb+2D0*EE/CW**2*MW*SW*Zb22**2*Zh32*sb+3D0
     &      *EE*MW/SW*Zb21**2*Zh31/cb*sb**2+EE/CW**2*MW*SW
     &      *Zb21**2*Zh31/cb*sb**2+2D0*EE/CW**2*MW*SW*Zb22
     &      **2*Zh31/cb*sb**2-6D0*Lambda*MbMM*Sqrt2*Zb21
     &      *Zb22*Zh33/cb*sb+6D0*Ab*EE/MW*MbMM/SW*Zb21*Zb22
     &      *Zh31/cb-6D0*EE*Lambda/MW*MbMM/SW*Zb21*Zb22*Zh32
     &      /cb*xvev
      AAABR(104) = 6D0*EE/MW*Mq**2/SW*Zh31/cb-3D0*EE*MW/SW
     &      *Zdd21**2*Zh31/cb+3D0*EE*MW/SW*Zdd21**2*Zh32
     &      *sb-EE/CW**2*MW*SW*Zdd21**2*Zh31/cb+EE/CW**2
     &      *MW*SW*Zdd21**2*Zh32*sb-2D0*EE/CW**2*MW*SW*Zdd22
     &      **2*Zh31/cb+2D0*EE/CW**2*MW*SW*Zdd22**2*Zh32
     &      *sb+3D0*EE*MW/SW*Zdd21**2*Zh31/cb*sb**2+EE/CW
     &      **2*MW*SW*Zdd21**2*Zh31/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zdd22**2*Zh31/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zdd21*Zdd22*Zh33/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zdd21
     &      *Zdd22*Zh31/cb-6D0*EE*Lambda/MW*Mq/SW*Zdd21*Zdd22
     &      *Zh32/cb*xvev
      AAABR(105) = 6D0*EE/MW*Mq**2/SW*Zh31/cb-3D0*EE*MW/SW
     &      *Zh31*Zss21**2/cb+3D0*EE*MW/SW*Zh32*Zss21**2
     &      *sb-EE/CW**2*MW*SW*Zh31*Zss21**2/cb-2D0*EE/CW
     &      **2*MW*SW*Zh31*Zss22**2/cb+EE/CW**2*MW*SW*Zh32
     &      *Zss21**2*sb+2D0*EE/CW**2*MW*SW*Zh32*Zss22**2
     &      *sb+3D0*EE*MW/SW*Zh31*Zss21**2/cb*sb**2+EE/CW
     &      **2*MW*SW*Zh31*Zss21**2/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zh31*Zss22**2/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zh33*Zss21*Zss22/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zh31
     &      *Zss21*Zss22/cb-6D0*EE*Lambda/MW*Mq/SW*Zh32*Zss21
     &      *Zss22/cb*xvev
      AAABR(106) = 6D0*EE/MW*MbMM**2/SW*Zh21/cb-3D0*EE*MW
     &      /SW*Zb21**2*Zh21/cb+3D0*EE*MW/SW*Zb21**2*Zh22
     &      *sb-EE/CW**2*MW*SW*Zb21**2*Zh21/cb+EE/CW**2*MW
     &      *SW*Zb21**2*Zh22*sb-2D0*EE/CW**2*MW*SW*Zb22**2
     &      *Zh21/cb+2D0*EE/CW**2*MW*SW*Zb22**2*Zh22*sb+3D0
     &      *EE*MW/SW*Zb21**2*Zh21/cb*sb**2+EE/CW**2*MW*SW
     &      *Zb21**2*Zh21/cb*sb**2+2D0*EE/CW**2*MW*SW*Zb22
     &      **2*Zh21/cb*sb**2-6D0*Lambda*MbMM*Sqrt2*Zb21
     &      *Zb22*Zh23/cb*sb+6D0*Ab*EE/MW*MbMM/SW*Zb21*Zb22
     &      *Zh21/cb-6D0*EE*Lambda/MW*MbMM/SW*Zb21*Zb22*Zh22
     &      /cb*xvev
      AAABR(107) = 6D0*EE/MW*MbMM**2/SW*Zh11/cb-3D0*EE*MW
     &      /SW*Zb21**2*Zh11/cb+3D0*EE*MW/SW*Zb21**2*Zh12
     &      *sb-EE/CW**2*MW*SW*Zb21**2*Zh11/cb+EE/CW**2*MW
     &      *SW*Zb21**2*Zh12*sb-2D0*EE/CW**2*MW*SW*Zb22**2
     &      *Zh11/cb+2D0*EE/CW**2*MW*SW*Zb22**2*Zh12*sb+3D0
     &      *EE*MW/SW*Zb21**2*Zh11/cb*sb**2+EE/CW**2*MW*SW
     &      *Zb21**2*Zh11/cb*sb**2+2D0*EE/CW**2*MW*SW*Zb22
     &      **2*Zh11/cb*sb**2-6D0*Lambda*MbMM*Sqrt2*Zb21
     &      *Zb22*Zh13/cb*sb+6D0*Ab*EE/MW*MbMM/SW*Zb21*Zb22
     &      *Zh11/cb-6D0*EE*Lambda/MW*MbMM/SW*Zb21*Zb22*Zh12
     &      /cb*xvev
      AAABR(108) = 6D0*EE/MW*Mq**2/SW*Zh21/cb-3D0*EE*MW/SW
     &      *Zdd21**2*Zh21/cb+3D0*EE*MW/SW*Zdd21**2*Zh22
     &      *sb-EE/CW**2*MW*SW*Zdd21**2*Zh21/cb+EE/CW**2
     &      *MW*SW*Zdd21**2*Zh22*sb-2D0*EE/CW**2*MW*SW*Zdd22
     &      **2*Zh21/cb+2D0*EE/CW**2*MW*SW*Zdd22**2*Zh22
     &      *sb+3D0*EE*MW/SW*Zdd21**2*Zh21/cb*sb**2+EE/CW
     &      **2*MW*SW*Zdd21**2*Zh21/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zdd22**2*Zh21/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zdd21*Zdd22*Zh23/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zdd21
     &      *Zdd22*Zh21/cb-6D0*EE*Lambda/MW*Mq/SW*Zdd21*Zdd22
     &      *Zh22/cb*xvev
      AAABR(109) = 6D0*EE/MW*Mq**2/SW*Zh21/cb-3D0*EE*MW/SW
     &      *Zh21*Zss21**2/cb+3D0*EE*MW/SW*Zh22*Zss21**2
     &      *sb-EE/CW**2*MW*SW*Zh21*Zss21**2/cb-2D0*EE/CW
     &      **2*MW*SW*Zh21*Zss22**2/cb+EE/CW**2*MW*SW*Zh22
     &      *Zss21**2*sb+2D0*EE/CW**2*MW*SW*Zh22*Zss22**2
     &      *sb+3D0*EE*MW/SW*Zh21*Zss21**2/cb*sb**2+EE/CW
     &      **2*MW*SW*Zh21*Zss21**2/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zh21*Zss22**2/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zh23*Zss21*Zss22/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zh21
     &      *Zss21*Zss22/cb-6D0*EE*Lambda/MW*Mq/SW*Zh22*Zss21
     &      *Zss22/cb*xvev
      AAABR(110) = 6D0*EE/MW*Mq**2/SW*Zh11/cb-3D0*EE*MW/SW
     &      *Zdd21**2*Zh11/cb+3D0*EE*MW/SW*Zdd21**2*Zh12
     &      *sb-EE/CW**2*MW*SW*Zdd21**2*Zh11/cb+EE/CW**2
     &      *MW*SW*Zdd21**2*Zh12*sb-2D0*EE/CW**2*MW*SW*Zdd22
     &      **2*Zh11/cb+2D0*EE/CW**2*MW*SW*Zdd22**2*Zh12
     &      *sb+3D0*EE*MW/SW*Zdd21**2*Zh11/cb*sb**2+EE/CW
     &      **2*MW*SW*Zdd21**2*Zh11/cb*sb**2+2D0*EE/CW**2
     &      *MW*SW*Zdd22**2*Zh11/cb*sb**2-6D0*Lambda*Mq*Sqrt2
     &      *Zdd21*Zdd22*Zh13/cb*sb+6D0*Ad*EE/MW*Mq/SW*Zdd21
     &      *Zdd22*Zh11/cb-6D0*EE*Lambda/MW*Mq/SW*Zdd21*Zdd22
     &      *Zh12/cb*xvev
      AAABR(111) = EE/MW*Mc/SW*(Au*Zcc11*Zcc22-Au*Zcc12*Zcc21
     &      -Lambda*Zcc11*Zcc22*cb/sb*xvev+Lambda*Zcc12*Zcc21
     &      *cb/sb*xvev)
      AAABR(112) = EE/MW*MtMM/SW*(At*Zt11*Zt22-At*Zt12*Zt21
     &      -Lambda*Zt11*Zt22*cb/sb*xvev+Lambda*Zt12*Zt21
     &      *cb/sb*xvev)
      AAABR(113) = EE/MW*Mq/SW*(Au*Zuu11*Zuu22-Au*Zuu12*Zuu21
     &      -Lambda*Zuu11*Zuu22*cb/sb*xvev+Lambda*Zuu12*Zuu21
     &      *cb/sb*xvev)
      AAABR(114) = Mq*(Lambda*Pa22*Sqrt2*Zuu11*Zuu22*cb/sb
     &      -Lambda*Pa22*Sqrt2*Zuu12*Zuu21*cb/sb+EE*Lambda
     &      /MW*Pa21/SW*Zuu11*Zuu22*xvev-EE*Lambda/MW*Pa21
     &      /SW*Zuu12*Zuu21*xvev+Au*EE/MW*Pa21/SW*Zuu11*Zuu22
     &      *cb/sb-Au*EE/MW*Pa21/SW*Zuu12*Zuu21*cb/sb)
      AAABR(115) = Mc*(Lambda*Pa22*Sqrt2*Zcc11*Zcc22*cb/sb
     &      -Lambda*Pa22*Sqrt2*Zcc12*Zcc21*cb/sb+EE*Lambda
     &      /MW*Pa21/SW*Zcc11*Zcc22*xvev-EE*Lambda/MW*Pa21
     &      /SW*Zcc12*Zcc21*xvev+Au*EE/MW*Pa21/SW*Zcc11*Zcc22
     &      *cb/sb-Au*EE/MW*Pa21/SW*Zcc12*Zcc21*cb/sb)
      AAABR(116) = Mq*(Lambda*Pa12*Sqrt2*Zuu11*Zuu22*cb/sb
     &      -Lambda*Pa12*Sqrt2*Zuu12*Zuu21*cb/sb+EE*Lambda
     &      /MW*Pa11/SW*Zuu11*Zuu22*xvev-EE*Lambda/MW*Pa11
     &      /SW*Zuu12*Zuu21*xvev+Au*EE/MW*Pa11/SW*Zuu11*Zuu22
     &      *cb/sb-Au*EE/MW*Pa11/SW*Zuu12*Zuu21*cb/sb)
      AAABR(117) = MtMM*(Lambda*Pa22*Sqrt2*Zt11*Zt22*cb/sb
     &      -Lambda*Pa22*Sqrt2*Zt12*Zt21*cb/sb+EE*Lambda
     &      /MW*Pa21/SW*Zt11*Zt22*xvev-EE*Lambda/MW*Pa21
     &      /SW*Zt12*Zt21*xvev+At*EE/MW*Pa21/SW*Zt11*Zt22
     &      *cb/sb-At*EE/MW*Pa21/SW*Zt12*Zt21*cb/sb)
      AAABR(118) = Mc*(Lambda*Pa12*Sqrt2*Zcc11*Zcc22*cb/sb
     &      -Lambda*Pa12*Sqrt2*Zcc12*Zcc21*cb/sb+EE*Lambda
     &      /MW*Pa11/SW*Zcc11*Zcc22*xvev-EE*Lambda/MW*Pa11
     &      /SW*Zcc12*Zcc21*xvev+Au*EE/MW*Pa11/SW*Zcc11*Zcc22
     &      *cb/sb-Au*EE/MW*Pa11/SW*Zcc12*Zcc21*cb/sb)
      AAABR(119) = MtMM*(Lambda*Pa12*Sqrt2*Zt11*Zt22*cb/sb
     &      -Lambda*Pa12*Sqrt2*Zt12*Zt21*cb/sb+EE*Lambda
     &      /MW*Pa11/SW*Zt11*Zt22*xvev-EE*Lambda/MW*Pa11
     &      /SW*Zt12*Zt21*xvev+At*EE/MW*Pa11/SW*Zt11*Zt22
     &      *cb/sb-At*EE/MW*Pa11/SW*Zt12*Zt21*cb/sb)
      AAABR(120) = 6D0*EE/MW*Mc**2/SW*Zh12/sb+3D0*EE*MW/SW
     &      *Zcc11**2*Zh11*cb-3D0*EE*MW/SW*Zcc11**2*Zh12
     &      *sb-EE/CW**2*MW*SW*Zcc11**2*Zh11*cb+EE/CW**2
     &      *MW*SW*Zcc11**2*Zh12*sb+4D0*EE/CW**2*MW*SW*Zcc12
     &      **2*Zh11*cb-4D0*EE/CW**2*MW*SW*Zcc12**2*Zh12
     &      *sb-6D0*Lambda*Mc*Sqrt2*Zcc11*Zcc12*Zh13*cb/sb
     &      +6D0*Au*EE/MW*Mc/SW*Zcc11*Zcc12*Zh12/sb-6D0*EE
     &      *Lambda/MW*Mc/SW*Zcc11*Zcc12*Zh11/sb*xvev
      AAABR(121) = 6D0*EE/MW*Mq**2/SW*Zh12/sb+3D0*EE*MW/SW
     &      *Zh11*Zuu11**2*cb-3D0*EE*MW/SW*Zh12*Zuu11**2
     &      *sb-EE/CW**2*MW*SW*Zh11*Zuu11**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh11*Zuu12**2*cb+EE/CW**2*MW*SW*Zh12
     &      *Zuu11**2*sb-4D0*EE/CW**2*MW*SW*Zh12*Zuu12**2
     &      *sb-6D0*Lambda*Mq*Sqrt2*Zh13*Zuu11*Zuu12*cb/sb
     &      +6D0*Au*EE/MW*Mq/SW*Zh12*Zuu11*Zuu12/sb-6D0*EE
     &      *Lambda/MW*Mq/SW*Zh11*Zuu11*Zuu12/sb*xvev
      AAABR(122) = 6D0*EE/MW*MtMM**2/SW*Zh12/sb+3D0*EE*MW
     &      /SW*Zh11*Zt11**2*cb-3D0*EE*MW/SW*Zh12*Zt11**2
     &      *sb-EE/CW**2*MW*SW*Zh11*Zt11**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh11*Zt12**2*cb+EE/CW**2*MW*SW*Zh12
     &      *Zt11**2*sb-4D0*EE/CW**2*MW*SW*Zh12*Zt12**2*sb
     &      -6D0*Lambda*MtMM*Sqrt2*Zh13*Zt11*Zt12*cb/sb+6D0
     &      *At*EE/MW*MtMM/SW*Zh12*Zt11*Zt12/sb-6D0*EE*Lambda
     &      /MW*MtMM/SW*Zh11*Zt11*Zt12/sb*xvev
      AAABR(123) = 6D0*EE/MW*MtMM**2/SW*Zh32/sb+3D0*EE*MW
     &      /SW*Zh31*Zt11**2*cb-3D0*EE*MW/SW*Zh32*Zt11**2
     &      *sb-EE/CW**2*MW*SW*Zh31*Zt11**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh31*Zt12**2*cb+EE/CW**2*MW*SW*Zh32
     &      *Zt11**2*sb-4D0*EE/CW**2*MW*SW*Zh32*Zt12**2*sb
     &      -6D0*Lambda*MtMM*Sqrt2*Zh33*Zt11*Zt12*cb/sb+6D0
     &      *At*EE/MW*MtMM/SW*Zh32*Zt11*Zt12/sb-6D0*EE*Lambda
     &      /MW*MtMM/SW*Zh31*Zt11*Zt12/sb*xvev
      AAABR(124) = 6D0*EE/MW*Mc**2/SW*Zh22/sb+3D0*EE*MW/SW
     &      *Zcc11**2*Zh21*cb-3D0*EE*MW/SW*Zcc11**2*Zh22
     &      *sb-EE/CW**2*MW*SW*Zcc11**2*Zh21*cb+EE/CW**2
     &      *MW*SW*Zcc11**2*Zh22*sb+4D0*EE/CW**2*MW*SW*Zcc12
     &      **2*Zh21*cb-4D0*EE/CW**2*MW*SW*Zcc12**2*Zh22
     &      *sb-6D0*Lambda*Mc*Sqrt2*Zcc11*Zcc12*Zh23*cb/sb
     &      +6D0*Au*EE/MW*Mc/SW*Zcc11*Zcc12*Zh22/sb-6D0*EE
     &      *Lambda/MW*Mc/SW*Zcc11*Zcc12*Zh21/sb*xvev
      AAABR(125) = 6D0*EE/MW*MtMM**2/SW*Zh22/sb+3D0*EE*MW
     &      /SW*Zh21*Zt11**2*cb-3D0*EE*MW/SW*Zh22*Zt11**2
     &      *sb-EE/CW**2*MW*SW*Zh21*Zt11**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh21*Zt12**2*cb+EE/CW**2*MW*SW*Zh22
     &      *Zt11**2*sb-4D0*EE/CW**2*MW*SW*Zh22*Zt12**2*sb
     &      -6D0*Lambda*MtMM*Sqrt2*Zh23*Zt11*Zt12*cb/sb+6D0
     &      *At*EE/MW*MtMM/SW*Zh22*Zt11*Zt12/sb-6D0*EE*Lambda
     &      /MW*MtMM/SW*Zh21*Zt11*Zt12/sb*xvev
      AAABR(126) = 6D0*EE/MW*Mq**2/SW*Zh32/sb+3D0*EE*MW/SW
     &      *Zh31*Zuu11**2*cb-3D0*EE*MW/SW*Zh32*Zuu11**2
     &      *sb-EE/CW**2*MW*SW*Zh31*Zuu11**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh31*Zuu12**2*cb+EE/CW**2*MW*SW*Zh32
     &      *Zuu11**2*sb-4D0*EE/CW**2*MW*SW*Zh32*Zuu12**2
     &      *sb-6D0*Lambda*Mq*Sqrt2*Zh33*Zuu11*Zuu12*cb/sb
     &      +6D0*Au*EE/MW*Mq/SW*Zh32*Zuu11*Zuu12/sb-6D0*EE
     &      *Lambda/MW*Mq/SW*Zh31*Zuu11*Zuu12/sb*xvev
      AAABR(127) = 6D0*EE/MW*Mq**2/SW*Zh22/sb+3D0*EE*MW/SW
     &      *Zh21*Zuu11**2*cb-3D0*EE*MW/SW*Zh22*Zuu11**2
     &      *sb-EE/CW**2*MW*SW*Zh21*Zuu11**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh21*Zuu12**2*cb+EE/CW**2*MW*SW*Zh22
     &      *Zuu11**2*sb-4D0*EE/CW**2*MW*SW*Zh22*Zuu12**2
     &      *sb-6D0*Lambda*Mq*Sqrt2*Zh23*Zuu11*Zuu12*cb/sb
     &      +6D0*Au*EE/MW*Mq/SW*Zh22*Zuu11*Zuu12/sb-6D0*EE
     &      *Lambda/MW*Mq/SW*Zh21*Zuu11*Zuu12/sb*xvev
      AAABR(128) = 6D0*EE/MW*Mc**2/SW*Zh32/sb+3D0*EE*MW/SW
     &      *Zcc11**2*Zh31*cb-3D0*EE*MW/SW*Zcc11**2*Zh32
     &      *sb-EE/CW**2*MW*SW*Zcc11**2*Zh31*cb+EE/CW**2
     &      *MW*SW*Zcc11**2*Zh32*sb+4D0*EE/CW**2*MW*SW*Zcc12
     &      **2*Zh31*cb-4D0*EE/CW**2*MW*SW*Zcc12**2*Zh32
     &      *sb-6D0*Lambda*Mc*Sqrt2*Zcc11*Zcc12*Zh33*cb/sb
     &      +6D0*Au*EE/MW*Mc/SW*Zcc11*Zcc12*Zh32/sb-6D0*EE
     &      *Lambda/MW*Mc/SW*Zcc11*Zcc12*Zh31/sb*xvev
      AAABR(129) = 3D0*EE*MW/SW*Zh31*Zt11*Zt21*cb-3D0*EE
     &      *MW/SW*Zh32*Zt11*Zt21*sb-EE/CW**2*MW*SW*Zh31
     &      *Zt11*Zt21*cb+4D0*EE/CW**2*MW*SW*Zh31*Zt12*Zt22
     &      *cb+EE/CW**2*MW*SW*Zh32*Zt11*Zt21*sb-4D0*EE/CW
     &      **2*MW*SW*Zh32*Zt12*Zt22*sb-3D0*Lambda*MtMM*Sqrt2
     &      *Zh33*Zt11*Zt22*cb/sb-3D0*Lambda*MtMM*Sqrt2*Zh33
     &      *Zt12*Zt21*cb/sb+3D0*At*EE/MW*MtMM/SW*Zh32*Zt11
     &      *Zt22/sb+3D0*At*EE/MW*MtMM/SW*Zh32*Zt12*Zt21
     &      /sb-3D0*EE*Lambda/MW*MtMM/SW*Zh31*Zt11*Zt22/sb
     &      *xvev-3D0*EE*Lambda/MW*MtMM/SW*Zh31*Zt12*Zt21
     &      /sb*xvev
      AAABR(130) = 3D0*EE*MW/SW*Zcc11*Zcc21*Zh11*cb-3D0*EE
     &      *MW/SW*Zcc11*Zcc21*Zh12*sb-EE/CW**2*MW*SW*Zcc11
     &      *Zcc21*Zh11*cb+EE/CW**2*MW*SW*Zcc11*Zcc21*Zh12
     &      *sb+4D0*EE/CW**2*MW*SW*Zcc12*Zcc22*Zh11*cb-4D0
     &      *EE/CW**2*MW*SW*Zcc12*Zcc22*Zh12*sb-3D0*Lambda
     &      *Mc*Sqrt2*Zcc11*Zcc22*Zh13*cb/sb-3D0*Lambda*Mc
     &      *Sqrt2*Zcc12*Zcc21*Zh13*cb/sb+3D0*Au*EE/MW*Mc
     &      /SW*Zcc11*Zcc22*Zh12/sb+3D0*Au*EE/MW*Mc/SW*Zcc12
     &      *Zcc21*Zh12/sb-3D0*EE*Lambda/MW*Mc/SW*Zcc11*Zcc22
     &      *Zh11/sb*xvev-3D0*EE*Lambda/MW*Mc/SW*Zcc12*Zcc21
     &      *Zh11/sb*xvev
      AAABR(131) = 3D0*EE*MW/SW*Zh21*Zuu11*Zuu21*cb-3D0*EE
     &      *MW/SW*Zh22*Zuu11*Zuu21*sb-EE/CW**2*MW*SW*Zh21
     &      *Zuu11*Zuu21*cb+4D0*EE/CW**2*MW*SW*Zh21*Zuu12
     &      *Zuu22*cb+EE/CW**2*MW*SW*Zh22*Zuu11*Zuu21*sb
     &      -4D0*EE/CW**2*MW*SW*Zh22*Zuu12*Zuu22*sb-3D0*Lambda
     &      *Mq*Sqrt2*Zh23*Zuu11*Zuu22*cb/sb-3D0*Lambda*Mq
     &      *Sqrt2*Zh23*Zuu12*Zuu21*cb/sb+3D0*Au*EE/MW*Mq
     &      /SW*Zh22*Zuu11*Zuu22/sb+3D0*Au*EE/MW*Mq/SW*Zh22
     &      *Zuu12*Zuu21/sb-3D0*EE*Lambda/MW*Mq/SW*Zh21*Zuu11
     &      *Zuu22/sb*xvev-3D0*EE*Lambda/MW*Mq/SW*Zh21*Zuu12
     &      *Zuu21/sb*xvev
      AAABR(132) = 3D0*EE*MW/SW*Zcc11*Zcc21*Zh21*cb-3D0*EE
     &      *MW/SW*Zcc11*Zcc21*Zh22*sb-EE/CW**2*MW*SW*Zcc11
     &      *Zcc21*Zh21*cb+EE/CW**2*MW*SW*Zcc11*Zcc21*Zh22
     &      *sb+4D0*EE/CW**2*MW*SW*Zcc12*Zcc22*Zh21*cb-4D0
     &      *EE/CW**2*MW*SW*Zcc12*Zcc22*Zh22*sb-3D0*Lambda
     &      *Mc*Sqrt2*Zcc11*Zcc22*Zh23*cb/sb-3D0*Lambda*Mc
     &      *Sqrt2*Zcc12*Zcc21*Zh23*cb/sb+3D0*Au*EE/MW*Mc
     &      /SW*Zcc11*Zcc22*Zh22/sb+3D0*Au*EE/MW*Mc/SW*Zcc12
     &      *Zcc21*Zh22/sb-3D0*EE*Lambda/MW*Mc/SW*Zcc11*Zcc22
     &      *Zh21/sb*xvev-3D0*EE*Lambda/MW*Mc/SW*Zcc12*Zcc21
     &      *Zh21/sb*xvev
      AAABR(133) = 3D0*EE*MW/SW*Zh11*Zt11*Zt21*cb-3D0*EE
     &      *MW/SW*Zh12*Zt11*Zt21*sb-EE/CW**2*MW*SW*Zh11
     &      *Zt11*Zt21*cb+4D0*EE/CW**2*MW*SW*Zh11*Zt12*Zt22
     &      *cb+EE/CW**2*MW*SW*Zh12*Zt11*Zt21*sb-4D0*EE/CW
     &      **2*MW*SW*Zh12*Zt12*Zt22*sb-3D0*Lambda*MtMM*Sqrt2
     &      *Zh13*Zt11*Zt22*cb/sb-3D0*Lambda*MtMM*Sqrt2*Zh13
     &      *Zt12*Zt21*cb/sb+3D0*At*EE/MW*MtMM/SW*Zh12*Zt11
     &      *Zt22/sb+3D0*At*EE/MW*MtMM/SW*Zh12*Zt12*Zt21
     &      /sb-3D0*EE*Lambda/MW*MtMM/SW*Zh11*Zt11*Zt22/sb
     &      *xvev-3D0*EE*Lambda/MW*MtMM/SW*Zh11*Zt12*Zt21
     &      /sb*xvev
      AAABR(134) = 3D0*EE*MW/SW*Zh21*Zt11*Zt21*cb-3D0*EE
     &      *MW/SW*Zh22*Zt11*Zt21*sb-EE/CW**2*MW*SW*Zh21
     &      *Zt11*Zt21*cb+4D0*EE/CW**2*MW*SW*Zh21*Zt12*Zt22
     &      *cb+EE/CW**2*MW*SW*Zh22*Zt11*Zt21*sb-4D0*EE/CW
     &      **2*MW*SW*Zh22*Zt12*Zt22*sb-3D0*Lambda*MtMM*Sqrt2
     &      *Zh23*Zt11*Zt22*cb/sb-3D0*Lambda*MtMM*Sqrt2*Zh23
     &      *Zt12*Zt21*cb/sb+3D0*At*EE/MW*MtMM/SW*Zh22*Zt11
     &      *Zt22/sb+3D0*At*EE/MW*MtMM/SW*Zh22*Zt12*Zt21
     &      /sb-3D0*EE*Lambda/MW*MtMM/SW*Zh21*Zt11*Zt22/sb
     &      *xvev-3D0*EE*Lambda/MW*MtMM/SW*Zh21*Zt12*Zt21
     &      /sb*xvev
      AAABR(135) = 3D0*EE*MW/SW*Zh11*Zuu11*Zuu21*cb-3D0*EE
     &      *MW/SW*Zh12*Zuu11*Zuu21*sb-EE/CW**2*MW*SW*Zh11
     &      *Zuu11*Zuu21*cb+4D0*EE/CW**2*MW*SW*Zh11*Zuu12
     &      *Zuu22*cb+EE/CW**2*MW*SW*Zh12*Zuu11*Zuu21*sb
     &      -4D0*EE/CW**2*MW*SW*Zh12*Zuu12*Zuu22*sb-3D0*Lambda
     &      *Mq*Sqrt2*Zh13*Zuu11*Zuu22*cb/sb-3D0*Lambda*Mq
     &      *Sqrt2*Zh13*Zuu12*Zuu21*cb/sb+3D0*Au*EE/MW*Mq
     &      /SW*Zh12*Zuu11*Zuu22/sb+3D0*Au*EE/MW*Mq/SW*Zh12
     &      *Zuu12*Zuu21/sb-3D0*EE*Lambda/MW*Mq/SW*Zh11*Zuu11
     &      *Zuu22/sb*xvev-3D0*EE*Lambda/MW*Mq/SW*Zh11*Zuu12
     &      *Zuu21/sb*xvev
      AAABR(136) = 3D0*EE*MW/SW*Zcc11*Zcc21*Zh31*cb-3D0*EE
     &      *MW/SW*Zcc11*Zcc21*Zh32*sb-EE/CW**2*MW*SW*Zcc11
     &      *Zcc21*Zh31*cb+EE/CW**2*MW*SW*Zcc11*Zcc21*Zh32
     &      *sb+4D0*EE/CW**2*MW*SW*Zcc12*Zcc22*Zh31*cb-4D0
     &      *EE/CW**2*MW*SW*Zcc12*Zcc22*Zh32*sb-3D0*Lambda
     &      *Mc*Sqrt2*Zcc11*Zcc22*Zh33*cb/sb-3D0*Lambda*Mc
     &      *Sqrt2*Zcc12*Zcc21*Zh33*cb/sb+3D0*Au*EE/MW*Mc
     &      /SW*Zcc11*Zcc22*Zh32/sb+3D0*Au*EE/MW*Mc/SW*Zcc12
     &      *Zcc21*Zh32/sb-3D0*EE*Lambda/MW*Mc/SW*Zcc11*Zcc22
     &      *Zh31/sb*xvev-3D0*EE*Lambda/MW*Mc/SW*Zcc12*Zcc21
     &      *Zh31/sb*xvev
      AAABR(137) = 3D0*EE*MW/SW*Zh31*Zuu11*Zuu21*cb-3D0*EE
     &      *MW/SW*Zh32*Zuu11*Zuu21*sb-EE/CW**2*MW*SW*Zh31
     &      *Zuu11*Zuu21*cb+4D0*EE/CW**2*MW*SW*Zh31*Zuu12
     &      *Zuu22*cb+EE/CW**2*MW*SW*Zh32*Zuu11*Zuu21*sb
     &      -4D0*EE/CW**2*MW*SW*Zh32*Zuu12*Zuu22*sb-3D0*Lambda
     &      *Mq*Sqrt2*Zh33*Zuu11*Zuu22*cb/sb-3D0*Lambda*Mq
     &      *Sqrt2*Zh33*Zuu12*Zuu21*cb/sb+3D0*Au*EE/MW*Mq
     &      /SW*Zh32*Zuu11*Zuu22/sb+3D0*Au*EE/MW*Mq/SW*Zh32
     &      *Zuu12*Zuu21/sb-3D0*EE*Lambda/MW*Mq/SW*Zh31*Zuu11
     &      *Zuu22/sb*xvev-3D0*EE*Lambda/MW*Mq/SW*Zh31*Zuu12
     &      *Zuu21/sb*xvev
      AAABR(138) = 6D0*EE/MW*Mq**2/SW*Zh12/sb+3D0*EE*MW/SW
     &      *Zh11*Zuu21**2*cb-3D0*EE*MW/SW*Zh12*Zuu21**2
     &      *sb-EE/CW**2*MW*SW*Zh11*Zuu21**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh11*Zuu22**2*cb+EE/CW**2*MW*SW*Zh12
     &      *Zuu21**2*sb-4D0*EE/CW**2*MW*SW*Zh12*Zuu22**2
     &      *sb-6D0*Lambda*Mq*Sqrt2*Zh13*Zuu21*Zuu22*cb/sb
     &      +6D0*Au*EE/MW*Mq/SW*Zh12*Zuu21*Zuu22/sb-6D0*EE
     &      *Lambda/MW*Mq/SW*Zh11*Zuu21*Zuu22/sb*xvev
      AAABR(139) = 6D0*EE/MW*Mq**2/SW*Zh32/sb+3D0*EE*MW/SW
     &      *Zh31*Zuu21**2*cb-3D0*EE*MW/SW*Zh32*Zuu21**2
     &      *sb-EE/CW**2*MW*SW*Zh31*Zuu21**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh31*Zuu22**2*cb+EE/CW**2*MW*SW*Zh32
     &      *Zuu21**2*sb-4D0*EE/CW**2*MW*SW*Zh32*Zuu22**2
     &      *sb-6D0*Lambda*Mq*Sqrt2*Zh33*Zuu21*Zuu22*cb/sb
     &      +6D0*Au*EE/MW*Mq/SW*Zh32*Zuu21*Zuu22/sb-6D0*EE
     &      *Lambda/MW*Mq/SW*Zh31*Zuu21*Zuu22/sb*xvev
      AAABR(140) = 6D0*EE/MW*Mc**2/SW*Zh12/sb+3D0*EE*MW/SW
     &      *Zcc21**2*Zh11*cb-3D0*EE*MW/SW*Zcc21**2*Zh12
     &      *sb-EE/CW**2*MW*SW*Zcc21**2*Zh11*cb+EE/CW**2
     &      *MW*SW*Zcc21**2*Zh12*sb+4D0*EE/CW**2*MW*SW*Zcc22
     &      **2*Zh11*cb-4D0*EE/CW**2*MW*SW*Zcc22**2*Zh12
     &      *sb-6D0*Lambda*Mc*Sqrt2*Zcc21*Zcc22*Zh13*cb/sb
     &      +6D0*Au*EE/MW*Mc/SW*Zcc21*Zcc22*Zh12/sb-6D0*EE
     &      *Lambda/MW*Mc/SW*Zcc21*Zcc22*Zh11/sb*xvev
      AAABR(141) = 6D0*EE/MW*MtMM**2/SW*Zh22/sb+3D0*EE*MW
     &      /SW*Zh21*Zt21**2*cb-3D0*EE*MW/SW*Zh22*Zt21**2
     &      *sb-EE/CW**2*MW*SW*Zh21*Zt21**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh21*Zt22**2*cb+EE/CW**2*MW*SW*Zh22
     &      *Zt21**2*sb-4D0*EE/CW**2*MW*SW*Zh22*Zt22**2*sb
     &      -6D0*Lambda*MtMM*Sqrt2*Zh23*Zt21*Zt22*cb/sb+6D0
     &      *At*EE/MW*MtMM/SW*Zh22*Zt21*Zt22/sb-6D0*EE*Lambda
     &      /MW*MtMM/SW*Zh21*Zt21*Zt22/sb*xvev
      AAABR(142) = 6D0*EE/MW*MtMM**2/SW*Zh12/sb+3D0*EE*MW
     &      /SW*Zh11*Zt21**2*cb-3D0*EE*MW/SW*Zh12*Zt21**2
     &      *sb-EE/CW**2*MW*SW*Zh11*Zt21**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh11*Zt22**2*cb+EE/CW**2*MW*SW*Zh12
     &      *Zt21**2*sb-4D0*EE/CW**2*MW*SW*Zh12*Zt22**2*sb
     &      -6D0*Lambda*MtMM*Sqrt2*Zh13*Zt21*Zt22*cb/sb+6D0
     &      *At*EE/MW*MtMM/SW*Zh12*Zt21*Zt22/sb-6D0*EE*Lambda
     &      /MW*MtMM/SW*Zh11*Zt21*Zt22/sb*xvev
      AAABR(143) = 6D0*EE/MW*Mq**2/SW*Zh22/sb+3D0*EE*MW/SW
     &      *Zh21*Zuu21**2*cb-3D0*EE*MW/SW*Zh22*Zuu21**2
     &      *sb-EE/CW**2*MW*SW*Zh21*Zuu21**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh21*Zuu22**2*cb+EE/CW**2*MW*SW*Zh22
     &      *Zuu21**2*sb-4D0*EE/CW**2*MW*SW*Zh22*Zuu22**2
     &      *sb-6D0*Lambda*Mq*Sqrt2*Zh23*Zuu21*Zuu22*cb/sb
     &      +6D0*Au*EE/MW*Mq/SW*Zh22*Zuu21*Zuu22/sb-6D0*EE
     &      *Lambda/MW*Mq/SW*Zh21*Zuu21*Zuu22/sb*xvev
      AAABR(144) = 6D0*EE/MW*Mc**2/SW*Zh32/sb+3D0*EE*MW/SW
     &      *Zcc21**2*Zh31*cb-3D0*EE*MW/SW*Zcc21**2*Zh32
     &      *sb-EE/CW**2*MW*SW*Zcc21**2*Zh31*cb+EE/CW**2
     &      *MW*SW*Zcc21**2*Zh32*sb+4D0*EE/CW**2*MW*SW*Zcc22
     &      **2*Zh31*cb-4D0*EE/CW**2*MW*SW*Zcc22**2*Zh32
     &      *sb-6D0*Lambda*Mc*Sqrt2*Zcc21*Zcc22*Zh33*cb/sb
     &      +6D0*Au*EE/MW*Mc/SW*Zcc21*Zcc22*Zh32/sb-6D0*EE
     &      *Lambda/MW*Mc/SW*Zcc21*Zcc22*Zh31/sb*xvev
      AAABR(145) = 6D0*EE/MW*Mc**2/SW*Zh22/sb+3D0*EE*MW/SW
     &      *Zcc21**2*Zh21*cb-3D0*EE*MW/SW*Zcc21**2*Zh22
     &      *sb-EE/CW**2*MW*SW*Zcc21**2*Zh21*cb+EE/CW**2
     &      *MW*SW*Zcc21**2*Zh22*sb+4D0*EE/CW**2*MW*SW*Zcc22
     &      **2*Zh21*cb-4D0*EE/CW**2*MW*SW*Zcc22**2*Zh22
     &      *sb-6D0*Lambda*Mc*Sqrt2*Zcc21*Zcc22*Zh23*cb/sb
     &      +6D0*Au*EE/MW*Mc/SW*Zcc21*Zcc22*Zh22/sb-6D0*EE
     &      *Lambda/MW*Mc/SW*Zcc21*Zcc22*Zh21/sb*xvev
      AAABR(146) = 6D0*EE/MW*MtMM**2/SW*Zh32/sb+3D0*EE*MW
     &      /SW*Zh31*Zt21**2*cb-3D0*EE*MW/SW*Zh32*Zt21**2
     &      *sb-EE/CW**2*MW*SW*Zh31*Zt21**2*cb+4D0*EE/CW
     &      **2*MW*SW*Zh31*Zt22**2*cb+EE/CW**2*MW*SW*Zh32
     &      *Zt21**2*sb-4D0*EE/CW**2*MW*SW*Zh32*Zt22**2*sb
     &      -6D0*Lambda*MtMM*Sqrt2*Zh33*Zt21*Zt22*cb/sb+6D0
     &      *At*EE/MW*MtMM/SW*Zh32*Zt21*Zt22/sb-6D0*EE*Lambda
     &      /MW*MtMM/SW*Zh31*Zt21*Zt22/sb*xvev
      AAABR(147) = Sqrt2*Zh23*aa1-Sqrt2*Zh23*aa1*sb**2+Sqrt2
     &      *Zh23*aa2*sb**2+Sqrt2*Zh23*la1s*xvev-Sqrt2*Zh23
     &      *aa5*cb*sb-Sqrt2*Zh23*aa6*cb*sb-Sqrt2*Zh23*la1s
     &      *sb**2*xvev+Sqrt2*Zh23*la2s*sb**2*xvev+2D0*MW
     &      /EE*SW*Zh21*cb*la1-6D0*MW/EE*SW*Zh21*la6*sb+6D0
     &      *MW/EE*SW*Zh21*la6*sb**3-2D0*MW/EE*SW*Zh21*la7
     &      *sb**3-2D0*MW/EE*SW*Zh22*cb*la6+2D0*MW/EE*SW
     &      *Zh22*la2*sb**3+2D0*MW/EE*SW*Zh22*la3*sb-2D0
     &      *MW/EE*SW*Zh22*la3*sb**3+2D0*MW/EE*SW*Zh22*la4
     &      *sb-2D0*MW/EE*SW*Zh22*la4*sb**3+2D0*MW/EE*SW
     &      *Zh22*la5*sb-2D0*MW/EE*SW*Zh22*la5*sb**3-2D0
     &      *Sqrt2*Zh23*cb*la4s*sb*xvev-2D0*Sqrt2*Zh23*cb
     &      *la5s*sb*xvev-2D0*Sqrt2*Zh23*cb*la6s*sb*xvev
     &      -2D0*MW/EE*SW*Zh21*cb*la1*sb**2+2D0*MW/EE*SW
     &      *Zh21*cb*la3*sb**2+2D0*MW/EE*SW*Zh21*cb*la4*sb
     &      **2+2D0*MW/EE*SW*Zh21*cb*la5*sb**2+2D0*MW/EE
     &      *SW*Zh22*cb*la6*sb**2-6D0*MW/EE*SW*Zh22*cb*la7
     &      *sb**2
      AAABR(148) = Sqrt2*Zh33*aa1-Sqrt2*Zh33*aa1*sb**2+Sqrt2
     &      *Zh33*aa2*sb**2+Sqrt2*Zh33*la1s*xvev-Sqrt2*Zh33
     &      *aa5*cb*sb-Sqrt2*Zh33*aa6*cb*sb-Sqrt2*Zh33*la1s
     &      *sb**2*xvev+Sqrt2*Zh33*la2s*sb**2*xvev+2D0*MW
     &      /EE*SW*Zh31*cb*la1-6D0*MW/EE*SW*Zh31*la6*sb+6D0
     &      *MW/EE*SW*Zh31*la6*sb**3-2D0*MW/EE*SW*Zh31*la7
     &      *sb**3-2D0*MW/EE*SW*Zh32*cb*la6+2D0*MW/EE*SW
     &      *Zh32*la2*sb**3+2D0*MW/EE*SW*Zh32*la3*sb-2D0
     &      *MW/EE*SW*Zh32*la3*sb**3+2D0*MW/EE*SW*Zh32*la4
     &      *sb-2D0*MW/EE*SW*Zh32*la4*sb**3+2D0*MW/EE*SW
     &      *Zh32*la5*sb-2D0*MW/EE*SW*Zh32*la5*sb**3-2D0
     &      *Sqrt2*Zh33*cb*la4s*sb*xvev-2D0*Sqrt2*Zh33*cb
     &      *la5s*sb*xvev-2D0*Sqrt2*Zh33*cb*la6s*sb*xvev
     &      -2D0*MW/EE*SW*Zh31*cb*la1*sb**2+2D0*MW/EE*SW
     &      *Zh31*cb*la3*sb**2+2D0*MW/EE*SW*Zh31*cb*la4*sb
     &      **2+2D0*MW/EE*SW*Zh31*cb*la5*sb**2+2D0*MW/EE
     &      *SW*Zh32*cb*la6*sb**2-6D0*MW/EE*SW*Zh32*cb*la7
     &      *sb**2
      AAABR(149) = Sqrt2*Zh13*aa1-Sqrt2*Zh13*aa1*sb**2+Sqrt2
     &      *Zh13*aa2*sb**2+Sqrt2*Zh13*la1s*xvev-Sqrt2*Zh13
     &      *aa5*cb*sb-Sqrt2*Zh13*aa6*cb*sb-Sqrt2*Zh13*la1s
     &      *sb**2*xvev+Sqrt2*Zh13*la2s*sb**2*xvev+2D0*MW
     &      /EE*SW*Zh11*cb*la1-6D0*MW/EE*SW*Zh11*la6*sb+6D0
     &      *MW/EE*SW*Zh11*la6*sb**3-2D0*MW/EE*SW*Zh11*la7
     &      *sb**3-2D0*MW/EE*SW*Zh12*cb*la6+2D0*MW/EE*SW
     &      *Zh12*la2*sb**3+2D0*MW/EE*SW*Zh12*la3*sb-2D0
     &      *MW/EE*SW*Zh12*la3*sb**3+2D0*MW/EE*SW*Zh12*la4
     &      *sb-2D0*MW/EE*SW*Zh12*la4*sb**3+2D0*MW/EE*SW
     &      *Zh12*la5*sb-2D0*MW/EE*SW*Zh12*la5*sb**3-2D0
     &      *Sqrt2*Zh13*cb*la4s*sb*xvev-2D0*Sqrt2*Zh13*cb
     &      *la5s*sb*xvev-2D0*Sqrt2*Zh13*cb*la6s*sb*xvev
     &      -2D0*MW/EE*SW*Zh11*cb*la1*sb**2+2D0*MW/EE*SW
     &      *Zh11*cb*la3*sb**2+2D0*MW/EE*SW*Zh11*cb*la4*sb
     &      **2+2D0*MW/EE*SW*Zh11*cb*la5*sb**2+2D0*MW/EE
     &      *SW*Zh12*cb*la6*sb**2-6D0*MW/EE*SW*Zh12*cb*la7
     &      *sb**2
      AAABR(150) = Pa21*Sqrt2*Zh23*aa5+Pa21*Sqrt2*Zh23*aa6
     &      -2D0*Pa21*Sqrt2*Zh23*aa5*sb**2-2D0*Pa21*Sqrt2
     &      *Zh23*aa6*sb**2+2D0*Pa21*Sqrt2*Zh23*la4s*xvev
     &      +2D0*Pa21*Sqrt2*Zh23*la5s*xvev+2D0*Pa21*Sqrt2
     &      *Zh23*la6s*xvev-Pa22*Sqrt2*Zh21*aa5*sb+Pa22*Sqrt2
     &      *Zh21*aa6*sb+Pa22*Sqrt2*Zh22*aa5*cb-Pa22*Sqrt2
     &      *Zh22*aa6*cb+2D0*Pa21*Sqrt2*Zh23*aa1*cb*sb-2D0
     &      *Pa21*Sqrt2*Zh23*aa2*cb*sb-4D0*Pa21*Sqrt2*Zh23
     &      *la4s*sb**2*xvev-4D0*Pa21*Sqrt2*Zh23*la5s*sb
     &      **2*xvev-4D0*Pa21*Sqrt2*Zh23*la6s*sb**2*xvev
     &      -2D0*Pa22*Sqrt2*Zh21*la5s*sb*xvev+2D0*Pa22*Sqrt2
     &      *Zh21*la6s*sb*xvev+2D0*Pa22*Sqrt2*Zh22*cb*la5s
     &      *xvev-2D0*Pa22*Sqrt2*Zh22*cb*la6s*xvev+4D0*MW
     &      /EE*Pa21*SW*Zh21*cb*la6+4D0*MW/EE*Pa21*SW*Zh21
     &      *la1*sb-4D0*MW/EE*Pa21*SW*Zh21*la1*sb**3-4D0
     &      *MW/EE*Pa21*SW*Zh21*la3*sb+4D0*MW/EE*Pa21*SW
     &      *Zh21*la3*sb**3-4D0*MW/EE*Pa21*SW*Zh21*la4*sb
     &      +4D0*MW/EE*Pa21*SW*Zh21*la4*sb**3+4D0*MW/EE*Pa21
     &      *SW*Zh21*la5*sb**3-4D0*MW/EE*Pa21*SW*Zh22*cb
     &      *la5-4D0*MW/EE*Pa21*SW*Zh22*la6*sb+4D0*MW/EE
     &      *Pa21*SW*Zh22*la6*sb**3+8D0*MW/EE*Pa21*SW*Zh22
     &      *la7*sb-12D0*MW/EE*Pa21*SW*Zh22*la7*sb**3+2D0
     &      *Pa21*Sqrt2*Zh23*cb*la1s*sb*xvev-2D0*Pa21*Sqrt2
     &      *Zh23*cb*la2s*sb*xvev-12D0*MW/EE*Pa21*SW*Zh21
     &      *cb*la6*sb**2+4D0*MW/EE*Pa21*SW*Zh21*cb*la7*sb
     &      **2-4D0*MW/EE*Pa21*SW*Zh22*cb*la2*sb**2+4D0*MW
     &      /EE*Pa21*SW*Zh22*cb*la3*sb**2+4D0*MW/EE*Pa21
     &      *SW*Zh22*cb*la4*sb**2+4D0*MW/EE*Pa21*SW*Zh22
     &      *cb*la5*sb**2
      AAABR(151) = Pa11*Sqrt2*Zh13*aa5+Pa11*Sqrt2*Zh13*aa6
     &      -2D0*Pa11*Sqrt2*Zh13*aa5*sb**2-2D0*Pa11*Sqrt2
     &      *Zh13*aa6*sb**2+2D0*Pa11*Sqrt2*Zh13*la4s*xvev
     &      +2D0*Pa11*Sqrt2*Zh13*la5s*xvev+2D0*Pa11*Sqrt2
     &      *Zh13*la6s*xvev-Pa12*Sqrt2*Zh11*aa5*sb+Pa12*Sqrt2
     &      *Zh11*aa6*sb+Pa12*Sqrt2*Zh12*aa5*cb-Pa12*Sqrt2
     &      *Zh12*aa6*cb+2D0*Pa11*Sqrt2*Zh13*aa1*cb*sb-2D0
     &      *Pa11*Sqrt2*Zh13*aa2*cb*sb-4D0*Pa11*Sqrt2*Zh13
     &      *la4s*sb**2*xvev-4D0*Pa11*Sqrt2*Zh13*la5s*sb
     &      **2*xvev-4D0*Pa11*Sqrt2*Zh13*la6s*sb**2*xvev
     &      -2D0*Pa12*Sqrt2*Zh11*la5s*sb*xvev+2D0*Pa12*Sqrt2
     &      *Zh11*la6s*sb*xvev+2D0*Pa12*Sqrt2*Zh12*cb*la5s
     &      *xvev-2D0*Pa12*Sqrt2*Zh12*cb*la6s*xvev+4D0*MW
     &      /EE*Pa11*SW*Zh11*cb*la6+4D0*MW/EE*Pa11*SW*Zh11
     &      *la1*sb-4D0*MW/EE*Pa11*SW*Zh11*la1*sb**3-4D0
     &      *MW/EE*Pa11*SW*Zh11*la3*sb+4D0*MW/EE*Pa11*SW
     &      *Zh11*la3*sb**3-4D0*MW/EE*Pa11*SW*Zh11*la4*sb
     &      +4D0*MW/EE*Pa11*SW*Zh11*la4*sb**3+4D0*MW/EE*Pa11
     &      *SW*Zh11*la5*sb**3-4D0*MW/EE*Pa11*SW*Zh12*cb
     &      *la5-4D0*MW/EE*Pa11*SW*Zh12*la6*sb+4D0*MW/EE
     &      *Pa11*SW*Zh12*la6*sb**3+8D0*MW/EE*Pa11*SW*Zh12
     &      *la7*sb-12D0*MW/EE*Pa11*SW*Zh12*la7*sb**3+2D0
     &      *Pa11*Sqrt2*Zh13*cb*la1s*sb*xvev-2D0*Pa11*Sqrt2
     &      *Zh13*cb*la2s*sb*xvev-12D0*MW/EE*Pa11*SW*Zh11
     &      *cb*la6*sb**2+4D0*MW/EE*Pa11*SW*Zh11*cb*la7*sb
     &      **2-4D0*MW/EE*Pa11*SW*Zh12*cb*la2*sb**2+4D0*MW
     &      /EE*Pa11*SW*Zh12*cb*la3*sb**2+4D0*MW/EE*Pa11
     &      *SW*Zh12*cb*la4*sb**2+4D0*MW/EE*Pa11*SW*Zh12
     &      *cb*la5*sb**2
      AAABR(152) = Pa21*Sqrt2*Zh33*aa5+Pa21*Sqrt2*Zh33*aa6
     &      -2D0*Pa21*Sqrt2*Zh33*aa5*sb**2-2D0*Pa21*Sqrt2
     &      *Zh33*aa6*sb**2+2D0*Pa21*Sqrt2*Zh33*la4s*xvev
     &      +2D0*Pa21*Sqrt2*Zh33*la5s*xvev+2D0*Pa21*Sqrt2
     &      *Zh33*la6s*xvev-Pa22*Sqrt2*Zh31*aa5*sb+Pa22*Sqrt2
     &      *Zh31*aa6*sb+Pa22*Sqrt2*Zh32*aa5*cb-Pa22*Sqrt2
     &      *Zh32*aa6*cb+2D0*Pa21*Sqrt2*Zh33*aa1*cb*sb-2D0
     &      *Pa21*Sqrt2*Zh33*aa2*cb*sb-4D0*Pa21*Sqrt2*Zh33
     &      *la4s*sb**2*xvev-4D0*Pa21*Sqrt2*Zh33*la5s*sb
     &      **2*xvev-4D0*Pa21*Sqrt2*Zh33*la6s*sb**2*xvev
     &      -2D0*Pa22*Sqrt2*Zh31*la5s*sb*xvev+2D0*Pa22*Sqrt2
     &      *Zh31*la6s*sb*xvev+2D0*Pa22*Sqrt2*Zh32*cb*la5s
     &      *xvev-2D0*Pa22*Sqrt2*Zh32*cb*la6s*xvev+4D0*MW
     &      /EE*Pa21*SW*Zh31*cb*la6+4D0*MW/EE*Pa21*SW*Zh31
     &      *la1*sb-4D0*MW/EE*Pa21*SW*Zh31*la1*sb**3-4D0
     &      *MW/EE*Pa21*SW*Zh31*la3*sb+4D0*MW/EE*Pa21*SW
     &      *Zh31*la3*sb**3-4D0*MW/EE*Pa21*SW*Zh31*la4*sb
     &      +4D0*MW/EE*Pa21*SW*Zh31*la4*sb**3+4D0*MW/EE*Pa21
     &      *SW*Zh31*la5*sb**3-4D0*MW/EE*Pa21*SW*Zh32*cb
     &      *la5-4D0*MW/EE*Pa21*SW*Zh32*la6*sb+4D0*MW/EE
     &      *Pa21*SW*Zh32*la6*sb**3+8D0*MW/EE*Pa21*SW*Zh32
     &      *la7*sb-12D0*MW/EE*Pa21*SW*Zh32*la7*sb**3+2D0
     &      *Pa21*Sqrt2*Zh33*cb*la1s*sb*xvev-2D0*Pa21*Sqrt2
     &      *Zh33*cb*la2s*sb*xvev-12D0*MW/EE*Pa21*SW*Zh31
     &      *cb*la6*sb**2+4D0*MW/EE*Pa21*SW*Zh31*cb*la7*sb
     &      **2-4D0*MW/EE*Pa21*SW*Zh32*cb*la2*sb**2+4D0*MW
     &      /EE*Pa21*SW*Zh32*cb*la3*sb**2+4D0*MW/EE*Pa21
     &      *SW*Zh32*cb*la4*sb**2+4D0*MW/EE*Pa21*SW*Zh32
     &      *cb*la5*sb**2
      AAABR(153) = Pa11*Sqrt2*Zh23*aa5+Pa11*Sqrt2*Zh23*aa6
     &      -2D0*Pa11*Sqrt2*Zh23*aa5*sb**2-2D0*Pa11*Sqrt2
     &      *Zh23*aa6*sb**2+2D0*Pa11*Sqrt2*Zh23*la4s*xvev
     &      +2D0*Pa11*Sqrt2*Zh23*la5s*xvev+2D0*Pa11*Sqrt2
     &      *Zh23*la6s*xvev-Pa12*Sqrt2*Zh21*aa5*sb+Pa12*Sqrt2
     &      *Zh21*aa6*sb+Pa12*Sqrt2*Zh22*aa5*cb-Pa12*Sqrt2
     &      *Zh22*aa6*cb+2D0*Pa11*Sqrt2*Zh23*aa1*cb*sb-2D0
     &      *Pa11*Sqrt2*Zh23*aa2*cb*sb-4D0*Pa11*Sqrt2*Zh23
     &      *la4s*sb**2*xvev-4D0*Pa11*Sqrt2*Zh23*la5s*sb
     &      **2*xvev-4D0*Pa11*Sqrt2*Zh23*la6s*sb**2*xvev
     &      -2D0*Pa12*Sqrt2*Zh21*la5s*sb*xvev+2D0*Pa12*Sqrt2
     &      *Zh21*la6s*sb*xvev+2D0*Pa12*Sqrt2*Zh22*cb*la5s
     &      *xvev-2D0*Pa12*Sqrt2*Zh22*cb*la6s*xvev+4D0*MW
     &      /EE*Pa11*SW*Zh21*cb*la6+4D0*MW/EE*Pa11*SW*Zh21
     &      *la1*sb-4D0*MW/EE*Pa11*SW*Zh21*la1*sb**3-4D0
     &      *MW/EE*Pa11*SW*Zh21*la3*sb+4D0*MW/EE*Pa11*SW
     &      *Zh21*la3*sb**3-4D0*MW/EE*Pa11*SW*Zh21*la4*sb
     &      +4D0*MW/EE*Pa11*SW*Zh21*la4*sb**3+4D0*MW/EE*Pa11
     &      *SW*Zh21*la5*sb**3-4D0*MW/EE*Pa11*SW*Zh22*cb
     &      *la5-4D0*MW/EE*Pa11*SW*Zh22*la6*sb+4D0*MW/EE
     &      *Pa11*SW*Zh22*la6*sb**3+8D0*MW/EE*Pa11*SW*Zh22
     &      *la7*sb-12D0*MW/EE*Pa11*SW*Zh22*la7*sb**3+2D0
     &      *Pa11*Sqrt2*Zh23*cb*la1s*sb*xvev-2D0*Pa11*Sqrt2
     &      *Zh23*cb*la2s*sb*xvev-12D0*MW/EE*Pa11*SW*Zh21
     &      *cb*la6*sb**2+4D0*MW/EE*Pa11*SW*Zh21*cb*la7*sb
     &      **2-4D0*MW/EE*Pa11*SW*Zh22*cb*la2*sb**2+4D0*MW
     &      /EE*Pa11*SW*Zh22*cb*la3*sb**2+4D0*MW/EE*Pa11
     &      *SW*Zh22*cb*la4*sb**2+4D0*MW/EE*Pa11*SW*Zh22
     &      *cb*la5*sb**2
      AAABR(154) = Pa21*Sqrt2*Zh13*aa5+Pa21*Sqrt2*Zh13*aa6
     &      -2D0*Pa21*Sqrt2*Zh13*aa5*sb**2-2D0*Pa21*Sqrt2
     &      *Zh13*aa6*sb**2+2D0*Pa21*Sqrt2*Zh13*la4s*xvev
     &      +2D0*Pa21*Sqrt2*Zh13*la5s*xvev+2D0*Pa21*Sqrt2
     &      *Zh13*la6s*xvev-Pa22*Sqrt2*Zh11*aa5*sb+Pa22*Sqrt2
     &      *Zh11*aa6*sb+Pa22*Sqrt2*Zh12*aa5*cb-Pa22*Sqrt2
     &      *Zh12*aa6*cb+2D0*Pa21*Sqrt2*Zh13*aa1*cb*sb-2D0
     &      *Pa21*Sqrt2*Zh13*aa2*cb*sb-4D0*Pa21*Sqrt2*Zh13
     &      *la4s*sb**2*xvev-4D0*Pa21*Sqrt2*Zh13*la5s*sb
     &      **2*xvev-4D0*Pa21*Sqrt2*Zh13*la6s*sb**2*xvev
     &      -2D0*Pa22*Sqrt2*Zh11*la5s*sb*xvev+2D0*Pa22*Sqrt2
     &      *Zh11*la6s*sb*xvev+2D0*Pa22*Sqrt2*Zh12*cb*la5s
     &      *xvev-2D0*Pa22*Sqrt2*Zh12*cb*la6s*xvev+4D0*MW
     &      /EE*Pa21*SW*Zh11*cb*la6+4D0*MW/EE*Pa21*SW*Zh11
     &      *la1*sb-4D0*MW/EE*Pa21*SW*Zh11*la1*sb**3-4D0
     &      *MW/EE*Pa21*SW*Zh11*la3*sb+4D0*MW/EE*Pa21*SW
     &      *Zh11*la3*sb**3-4D0*MW/EE*Pa21*SW*Zh11*la4*sb
     &      +4D0*MW/EE*Pa21*SW*Zh11*la4*sb**3+4D0*MW/EE*Pa21
     &      *SW*Zh11*la5*sb**3-4D0*MW/EE*Pa21*SW*Zh12*cb
     &      *la5-4D0*MW/EE*Pa21*SW*Zh12*la6*sb+4D0*MW/EE
     &      *Pa21*SW*Zh12*la6*sb**3+8D0*MW/EE*Pa21*SW*Zh12
     &      *la7*sb-12D0*MW/EE*Pa21*SW*Zh12*la7*sb**3+2D0
     &      *Pa21*Sqrt2*Zh13*cb*la1s*sb*xvev-2D0*Pa21*Sqrt2
     &      *Zh13*cb*la2s*sb*xvev-12D0*MW/EE*Pa21*SW*Zh11
     &      *cb*la6*sb**2+4D0*MW/EE*Pa21*SW*Zh11*cb*la7*sb
     &      **2-4D0*MW/EE*Pa21*SW*Zh12*cb*la2*sb**2+4D0*MW
     &      /EE*Pa21*SW*Zh12*cb*la3*sb**2+4D0*MW/EE*Pa21
     &      *SW*Zh12*cb*la4*sb**2+4D0*MW/EE*Pa21*SW*Zh12
     &      *cb*la5*sb**2
      AAABR(155) = Pa11*Sqrt2*Zh33*aa5+Pa11*Sqrt2*Zh33*aa6
     &      -2D0*Pa11*Sqrt2*Zh33*aa5*sb**2-2D0*Pa11*Sqrt2
     &      *Zh33*aa6*sb**2+2D0*Pa11*Sqrt2*Zh33*la4s*xvev
     &      +2D0*Pa11*Sqrt2*Zh33*la5s*xvev+2D0*Pa11*Sqrt2
     &      *Zh33*la6s*xvev-Pa12*Sqrt2*Zh31*aa5*sb+Pa12*Sqrt2
     &      *Zh31*aa6*sb+Pa12*Sqrt2*Zh32*aa5*cb-Pa12*Sqrt2
     &      *Zh32*aa6*cb+2D0*Pa11*Sqrt2*Zh33*aa1*cb*sb-2D0
     &      *Pa11*Sqrt2*Zh33*aa2*cb*sb-4D0*Pa11*Sqrt2*Zh33
     &      *la4s*sb**2*xvev-4D0*Pa11*Sqrt2*Zh33*la5s*sb
     &      **2*xvev-4D0*Pa11*Sqrt2*Zh33*la6s*sb**2*xvev
     &      -2D0*Pa12*Sqrt2*Zh31*la5s*sb*xvev+2D0*Pa12*Sqrt2
     &      *Zh31*la6s*sb*xvev+2D0*Pa12*Sqrt2*Zh32*cb*la5s
     &      *xvev-2D0*Pa12*Sqrt2*Zh32*cb*la6s*xvev+4D0*MW
     &      /EE*Pa11*SW*Zh31*cb*la6+4D0*MW/EE*Pa11*SW*Zh31
     &      *la1*sb-4D0*MW/EE*Pa11*SW*Zh31*la1*sb**3-4D0
     &      *MW/EE*Pa11*SW*Zh31*la3*sb+4D0*MW/EE*Pa11*SW
     &      *Zh31*la3*sb**3-4D0*MW/EE*Pa11*SW*Zh31*la4*sb
     &      +4D0*MW/EE*Pa11*SW*Zh31*la4*sb**3+4D0*MW/EE*Pa11
     &      *SW*Zh31*la5*sb**3-4D0*MW/EE*Pa11*SW*Zh32*cb
     &      *la5-4D0*MW/EE*Pa11*SW*Zh32*la6*sb+4D0*MW/EE
     &      *Pa11*SW*Zh32*la6*sb**3+8D0*MW/EE*Pa11*SW*Zh32
     &      *la7*sb-12D0*MW/EE*Pa11*SW*Zh32*la7*sb**3+2D0
     &      *Pa11*Sqrt2*Zh33*cb*la1s*sb*xvev-2D0*Pa11*Sqrt2
     &      *Zh33*cb*la2s*sb*xvev-12D0*MW/EE*Pa11*SW*Zh31
     &      *cb*la6*sb**2+4D0*MW/EE*Pa11*SW*Zh31*cb*la7*sb
     &      **2-4D0*MW/EE*Pa11*SW*Zh32*cb*la2*sb**2+4D0*MW
     &      /EE*Pa11*SW*Zh32*cb*la3*sb**2+4D0*MW/EE*Pa11
     &      *SW*Zh32*cb*la4*sb**2+4D0*MW/EE*Pa11*SW*Zh32
     &      *cb*la5*sb**2
      AAABR(156) = Pa11**2*Sqrt2*Zh23*aa2+Pa12**2*Sqrt2*Zh23
     &      *aa3-3D0*Pa12**2*Sqrt2*Zh23*aa4+Pa11**2*Sqrt2
     &      *Zh23*aa1*sb**2-Pa11**2*Sqrt2*Zh23*aa2*sb**2
     &      +Pa11**2*Sqrt2*Zh23*la2s*xvev+2D0*Pa12**2*Sqrt2
     &      *Zh23*la3s*xvev-12D0*Pa12**2*Sqrt2*Zh23*la7s
     &      *xvev+Pa11*Pa12*Sqrt2*Zh21*aa5*cb-Pa11*Pa12*Sqrt2
     &      *Zh21*aa6*cb+Pa11*Pa12*Sqrt2*Zh22*aa5*sb-Pa11
     &      *Pa12*Sqrt2*Zh22*aa6*sb+Pa11**2*Sqrt2*Zh23*aa5
     &      *cb*sb+Pa11**2*Sqrt2*Zh23*aa6*cb*sb+Pa11**2*Sqrt2
     &      *Zh23*la1s*sb**2*xvev-Pa11**2*Sqrt2*Zh23*la2s
     &      *sb**2*xvev+4D0*MW/EE*Pa11*Pa12*SW*Zh23*la5s
     &      -4D0*MW/EE*Pa11*Pa12*SW*Zh23*la6s+2D0*MW/EE*Pa11
     &      **2*SW*Zh21*cb*la3+2D0*MW/EE*Pa11**2*SW*Zh21
     &      *cb*la4-2D0*MW/EE*Pa11**2*SW*Zh21*cb*la5+4D0
     &      *MW/EE*Pa11**2*SW*Zh21*la6*sb-6D0*MW/EE*Pa11
     &      **2*SW*Zh21*la6*sb**3-2D0*MW/EE*Pa11**2*SW*Zh21
     &      *la7*sb+2D0*MW/EE*Pa11**2*SW*Zh21*la7*sb**3-2D0
     &      *MW/EE*Pa11**2*SW*Zh22*cb*la7+2D0*MW/EE*Pa11
     &      **2*SW*Zh22*la2*sb-2D0*MW/EE*Pa11**2*SW*Zh22
     &      *la2*sb**3+2D0*MW/EE*Pa11**2*SW*Zh22*la3*sb**3
     &      +2D0*MW/EE*Pa11**2*SW*Zh22*la4*sb**3-4D0*MW/EE
     &      *Pa11**2*SW*Zh22*la5*sb+2D0*MW/EE*Pa11**2*SW
     &      *Zh22*la5*sb**3+2D0*MW/EE*Pa12**2*SW*Zh21*cb
     &      *la1s-2D0*MW/EE*Pa12**2*SW*Zh21*la4s*sb+2D0*MW
     &      /EE*Pa12**2*SW*Zh21*la5s*sb+2D0*MW/EE*Pa12**2
     &      *SW*Zh21*la6s*sb-2D0*MW/EE*Pa12**2*SW*Zh22*cb
     &      *la4s+2D0*MW/EE*Pa12**2*SW*Zh22*cb*la5s+2D0*MW
     &      /EE*Pa12**2*SW*Zh22*cb*la6s+2D0*MW/EE*Pa12**2
     &      *SW*Zh22*la2s*sb+2D0*Pa11*Pa12*Sqrt2*Zh21*cb
     &      *la5s*xvev-2D0*Pa11*Pa12*Sqrt2*Zh21*cb*la6s*xvev
     &      +2D0*Pa11*Pa12*Sqrt2*Zh22*la5s*sb*xvev-2D0*Pa11
     &      *Pa12*Sqrt2*Zh22*la6s*sb*xvev+2D0*Pa11**2*Sqrt2
     &      *Zh23*cb*la4s*sb*xvev+2D0*Pa11**2*Sqrt2*Zh23
     &      *cb*la5s*sb*xvev+2D0*Pa11**2*Sqrt2*Zh23*cb*la6s
     &      *sb*xvev+2D0*MW/EE*Pa11**2*SW*Zh21*cb*la1*sb
     &      **2-2D0*MW/EE*Pa11**2*SW*Zh21*cb*la3*sb**2-2D0
     &      *MW/EE*Pa11**2*SW*Zh21*cb*la4*sb**2-2D0*MW/EE
     &      *Pa11**2*SW*Zh21*cb*la5*sb**2-2D0*MW/EE*Pa11
     &      **2*SW*Zh22*cb*la6*sb**2+6D0*MW/EE*Pa11**2*SW
     &      *Zh22*cb*la7*sb**2
      AAABR(157) = Pa11**2*Sqrt2*Zh33*aa2+Pa12**2*Sqrt2*Zh33
     &      *aa3-3D0*Pa12**2*Sqrt2*Zh33*aa4+Pa11**2*Sqrt2
     &      *Zh33*aa1*sb**2-Pa11**2*Sqrt2*Zh33*aa2*sb**2
     &      +Pa11**2*Sqrt2*Zh33*la2s*xvev+2D0*Pa12**2*Sqrt2
     &      *Zh33*la3s*xvev-12D0*Pa12**2*Sqrt2*Zh33*la7s
     &      *xvev+Pa11*Pa12*Sqrt2*Zh31*aa5*cb-Pa11*Pa12*Sqrt2
     &      *Zh31*aa6*cb+Pa11*Pa12*Sqrt2*Zh32*aa5*sb-Pa11
     &      *Pa12*Sqrt2*Zh32*aa6*sb+Pa11**2*Sqrt2*Zh33*aa5
     &      *cb*sb+Pa11**2*Sqrt2*Zh33*aa6*cb*sb+Pa11**2*Sqrt2
     &      *Zh33*la1s*sb**2*xvev-Pa11**2*Sqrt2*Zh33*la2s
     &      *sb**2*xvev+4D0*MW/EE*Pa11*Pa12*SW*Zh33*la5s
     &      -4D0*MW/EE*Pa11*Pa12*SW*Zh33*la6s+2D0*MW/EE*Pa11
     &      **2*SW*Zh31*cb*la3+2D0*MW/EE*Pa11**2*SW*Zh31
     &      *cb*la4-2D0*MW/EE*Pa11**2*SW*Zh31*cb*la5+4D0
     &      *MW/EE*Pa11**2*SW*Zh31*la6*sb-6D0*MW/EE*Pa11
     &      **2*SW*Zh31*la6*sb**3-2D0*MW/EE*Pa11**2*SW*Zh31
     &      *la7*sb+2D0*MW/EE*Pa11**2*SW*Zh31*la7*sb**3-2D0
     &      *MW/EE*Pa11**2*SW*Zh32*cb*la7+2D0*MW/EE*Pa11
     &      **2*SW*Zh32*la2*sb-2D0*MW/EE*Pa11**2*SW*Zh32
     &      *la2*sb**3+2D0*MW/EE*Pa11**2*SW*Zh32*la3*sb**3
     &      +2D0*MW/EE*Pa11**2*SW*Zh32*la4*sb**3-4D0*MW/EE
     &      *Pa11**2*SW*Zh32*la5*sb+2D0*MW/EE*Pa11**2*SW
     &      *Zh32*la5*sb**3+2D0*MW/EE*Pa12**2*SW*Zh31*cb
     &      *la1s-2D0*MW/EE*Pa12**2*SW*Zh31*la4s*sb+2D0*MW
     &      /EE*Pa12**2*SW*Zh31*la5s*sb+2D0*MW/EE*Pa12**2
     &      *SW*Zh31*la6s*sb-2D0*MW/EE*Pa12**2*SW*Zh32*cb
     &      *la4s+2D0*MW/EE*Pa12**2*SW*Zh32*cb*la5s+2D0*MW
     &      /EE*Pa12**2*SW*Zh32*cb*la6s+2D0*MW/EE*Pa12**2
     &      *SW*Zh32*la2s*sb+2D0*Pa11*Pa12*Sqrt2*Zh31*cb
     &      *la5s*xvev-2D0*Pa11*Pa12*Sqrt2*Zh31*cb*la6s*xvev
     &      +2D0*Pa11*Pa12*Sqrt2*Zh32*la5s*sb*xvev-2D0*Pa11
     &      *Pa12*Sqrt2*Zh32*la6s*sb*xvev+2D0*Pa11**2*Sqrt2
     &      *Zh33*cb*la4s*sb*xvev+2D0*Pa11**2*Sqrt2*Zh33
     &      *cb*la5s*sb*xvev+2D0*Pa11**2*Sqrt2*Zh33*cb*la6s
     &      *sb*xvev+2D0*MW/EE*Pa11**2*SW*Zh31*cb*la1*sb
     &      **2-2D0*MW/EE*Pa11**2*SW*Zh31*cb*la3*sb**2-2D0
     &      *MW/EE*Pa11**2*SW*Zh31*cb*la4*sb**2-2D0*MW/EE
     &      *Pa11**2*SW*Zh31*cb*la5*sb**2-2D0*MW/EE*Pa11
     &      **2*SW*Zh32*cb*la6*sb**2+6D0*MW/EE*Pa11**2*SW
     &      *Zh32*cb*la7*sb**2
      AAABR(158) = Pa21**2*Sqrt2*Zh13*aa2+Pa22**2*Sqrt2*Zh13
     &      *aa3-3D0*Pa22**2*Sqrt2*Zh13*aa4+Pa21**2*Sqrt2
     &      *Zh13*aa1*sb**2-Pa21**2*Sqrt2*Zh13*aa2*sb**2
     &      +Pa21**2*Sqrt2*Zh13*la2s*xvev+2D0*Pa22**2*Sqrt2
     &      *Zh13*la3s*xvev-12D0*Pa22**2*Sqrt2*Zh13*la7s
     &      *xvev+Pa21*Pa22*Sqrt2*Zh11*aa5*cb-Pa21*Pa22*Sqrt2
     &      *Zh11*aa6*cb+Pa21*Pa22*Sqrt2*Zh12*aa5*sb-Pa21
     &      *Pa22*Sqrt2*Zh12*aa6*sb+Pa21**2*Sqrt2*Zh13*aa5
     &      *cb*sb+Pa21**2*Sqrt2*Zh13*aa6*cb*sb+Pa21**2*Sqrt2
     &      *Zh13*la1s*sb**2*xvev-Pa21**2*Sqrt2*Zh13*la2s
     &      *sb**2*xvev+4D0*MW/EE*Pa21*Pa22*SW*Zh13*la5s
     &      -4D0*MW/EE*Pa21*Pa22*SW*Zh13*la6s+2D0*MW/EE*Pa21
     &      **2*SW*Zh11*cb*la3+2D0*MW/EE*Pa21**2*SW*Zh11
     &      *cb*la4-2D0*MW/EE*Pa21**2*SW*Zh11*cb*la5+4D0
     &      *MW/EE*Pa21**2*SW*Zh11*la6*sb-6D0*MW/EE*Pa21
     &      **2*SW*Zh11*la6*sb**3-2D0*MW/EE*Pa21**2*SW*Zh11
     &      *la7*sb+2D0*MW/EE*Pa21**2*SW*Zh11*la7*sb**3-2D0
     &      *MW/EE*Pa21**2*SW*Zh12*cb*la7+2D0*MW/EE*Pa21
     &      **2*SW*Zh12*la2*sb-2D0*MW/EE*Pa21**2*SW*Zh12
     &      *la2*sb**3+2D0*MW/EE*Pa21**2*SW*Zh12*la3*sb**3
     &      +2D0*MW/EE*Pa21**2*SW*Zh12*la4*sb**3-4D0*MW/EE
     &      *Pa21**2*SW*Zh12*la5*sb+2D0*MW/EE*Pa21**2*SW
     &      *Zh12*la5*sb**3+2D0*MW/EE*Pa22**2*SW*Zh11*cb
     &      *la1s-2D0*MW/EE*Pa22**2*SW*Zh11*la4s*sb+2D0*MW
     &      /EE*Pa22**2*SW*Zh11*la5s*sb+2D0*MW/EE*Pa22**2
     &      *SW*Zh11*la6s*sb-2D0*MW/EE*Pa22**2*SW*Zh12*cb
     &      *la4s+2D0*MW/EE*Pa22**2*SW*Zh12*cb*la5s+2D0*MW
     &      /EE*Pa22**2*SW*Zh12*cb*la6s+2D0*MW/EE*Pa22**2
     &      *SW*Zh12*la2s*sb+2D0*Pa21*Pa22*Sqrt2*Zh11*cb
     &      *la5s*xvev-2D0*Pa21*Pa22*Sqrt2*Zh11*cb*la6s*xvev
     &      +2D0*Pa21*Pa22*Sqrt2*Zh12*la5s*sb*xvev-2D0*Pa21
     &      *Pa22*Sqrt2*Zh12*la6s*sb*xvev+2D0*Pa21**2*Sqrt2
     &      *Zh13*cb*la4s*sb*xvev+2D0*Pa21**2*Sqrt2*Zh13
     &      *cb*la5s*sb*xvev+2D0*Pa21**2*Sqrt2*Zh13*cb*la6s
     &      *sb*xvev+2D0*MW/EE*Pa21**2*SW*Zh11*cb*la1*sb
     &      **2-2D0*MW/EE*Pa21**2*SW*Zh11*cb*la3*sb**2-2D0
     &      *MW/EE*Pa21**2*SW*Zh11*cb*la4*sb**2-2D0*MW/EE
     &      *Pa21**2*SW*Zh11*cb*la5*sb**2-2D0*MW/EE*Pa21
     &      **2*SW*Zh12*cb*la6*sb**2+6D0*MW/EE*Pa21**2*SW
     &      *Zh12*cb*la7*sb**2
      AAABR(159) = 2D0*Pa11*Pa21*Sqrt2*Zh23*aa2+2D0*Pa12
     &      *Pa22*Sqrt2*Zh23*aa3-6D0*Pa12*Pa22*Sqrt2*Zh23
     &      *aa4+2D0*Pa11*Pa21*Sqrt2*Zh23*aa1*sb**2-2D0*Pa11
     &      *Pa21*Sqrt2*Zh23*aa2*sb**2+2D0*Pa11*Pa21*Sqrt2
     &      *Zh23*la2s*xvev+Pa11*Pa22*Sqrt2*Zh21*aa5*cb-Pa11
     &      *Pa22*Sqrt2*Zh21*aa6*cb+Pa11*Pa22*Sqrt2*Zh22
     &      *aa5*sb-Pa11*Pa22*Sqrt2*Zh22*aa6*sb+Pa12*Pa21
     &      *Sqrt2*Zh21*aa5*cb-Pa12*Pa21*Sqrt2*Zh21*aa6*cb
     &      +Pa12*Pa21*Sqrt2*Zh22*aa5*sb-Pa12*Pa21*Sqrt2
     &      *Zh22*aa6*sb+4D0*Pa12*Pa22*Sqrt2*Zh23*la3s*xvev
     &      -24D0*Pa12*Pa22*Sqrt2*Zh23*la7s*xvev+4D0*MW/EE
     &      *Pa11*Pa22*SW*Zh23*la5s-4D0*MW/EE*Pa11*Pa22*SW
     &      *Zh23*la6s+4D0*MW/EE*Pa12*Pa21*SW*Zh23*la5s-4D0
     &      *MW/EE*Pa12*Pa21*SW*Zh23*la6s+2D0*Pa11*Pa21*Sqrt2
     &      *Zh23*aa5*cb*sb+2D0*Pa11*Pa21*Sqrt2*Zh23*aa6
     &      *cb*sb+2D0*Pa11*Pa21*Sqrt2*Zh23*la1s*sb**2*xvev
     &      -2D0*Pa11*Pa21*Sqrt2*Zh23*la2s*sb**2*xvev+2D0
     &      *Pa11*Pa22*Sqrt2*Zh21*cb*la5s*xvev-2D0*Pa11*Pa22
     &      *Sqrt2*Zh21*cb*la6s*xvev+2D0*Pa11*Pa22*Sqrt2
     &      *Zh22*la5s*sb*xvev-2D0*Pa11*Pa22*Sqrt2*Zh22*la6s
     &      *sb*xvev+2D0*Pa12*Pa21*Sqrt2*Zh21*cb*la5s*xvev
     &      -2D0*Pa12*Pa21*Sqrt2*Zh21*cb*la6s*xvev+2D0*Pa12
     &      *Pa21*Sqrt2*Zh22*la5s*sb*xvev-2D0*Pa12*Pa21*Sqrt2
     &      *Zh22*la6s*sb*xvev+4D0*MW/EE*Pa11*Pa21*SW*Zh21
     &      *cb*la3+4D0*MW/EE*Pa11*Pa21*SW*Zh21*cb*la4-4D0
     &      *MW/EE*Pa11*Pa21*SW*Zh21*cb*la5+8D0*MW/EE*Pa11
     &      *Pa21*SW*Zh21*la6*sb-12D0*MW/EE*Pa11*Pa21*SW
     &      *Zh21*la6*sb**3-4D0*MW/EE*Pa11*Pa21*SW*Zh21*la7
     &      *sb+4D0*MW/EE*Pa11*Pa21*SW*Zh21*la7*sb**3-4D0
     &      *MW/EE*Pa11*Pa21*SW*Zh22*cb*la7+4D0*MW/EE*Pa11
     &      *Pa21*SW*Zh22*la2*sb-4D0*MW/EE*Pa11*Pa21*SW*Zh22
     &      *la2*sb**3+4D0*MW/EE*Pa11*Pa21*SW*Zh22*la3*sb
     &      **3+4D0*MW/EE*Pa11*Pa21*SW*Zh22*la4*sb**3-8D0
     &      *MW/EE*Pa11*Pa21*SW*Zh22*la5*sb+4D0*MW/EE*Pa11
     &      *Pa21*SW*Zh22*la5*sb**3+4D0*MW/EE*Pa12*Pa22*SW
     &      *Zh21*cb*la1s-4D0*MW/EE*Pa12*Pa22*SW*Zh21*la4s
     &      *sb+4D0*MW/EE*Pa12*Pa22*SW*Zh21*la5s*sb+4D0*MW
     &      /EE*Pa12*Pa22*SW*Zh21*la6s*sb-4D0*MW/EE*Pa12
     &      *Pa22*SW*Zh22*cb*la4s+4D0*MW/EE*Pa12*Pa22*SW
     &      *Zh22*cb*la5s+4D0*MW/EE*Pa12*Pa22*SW*Zh22*cb
     &      *la6s+4D0*MW/EE*Pa12*Pa22*SW*Zh22*la2s*sb+4D0
     &      *Pa11*Pa21*Sqrt2*Zh23*cb*la4s*sb*xvev+4D0*Pa11
     &      *Pa21*Sqrt2*Zh23*cb*la5s*sb*xvev+4D0*Pa11*Pa21
     &      *Sqrt2*Zh23*cb*la6s*sb*xvev+4D0*MW/EE*Pa11*Pa21
     &      *SW*Zh21*cb*la1*sb**2-4D0*MW/EE*Pa11*Pa21*SW
     &      *Zh21*cb*la3*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh21
     &      *cb*la4*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh21*cb
     &      *la5*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh22*cb*la6
     &      *sb**2+12D0*MW/EE*Pa11*Pa21*SW*Zh22*cb*la7*sb**2
      AAABR(160) = Pa21**2*Sqrt2*Zh23*aa2+Pa22**2*Sqrt2*Zh23
     &      *aa3-3D0*Pa22**2*Sqrt2*Zh23*aa4+Pa21**2*Sqrt2
     &      *Zh23*aa1*sb**2-Pa21**2*Sqrt2*Zh23*aa2*sb**2
     &      +Pa21**2*Sqrt2*Zh23*la2s*xvev+2D0*Pa22**2*Sqrt2
     &      *Zh23*la3s*xvev-12D0*Pa22**2*Sqrt2*Zh23*la7s
     &      *xvev+Pa21*Pa22*Sqrt2*Zh21*aa5*cb-Pa21*Pa22*Sqrt2
     &      *Zh21*aa6*cb+Pa21*Pa22*Sqrt2*Zh22*aa5*sb-Pa21
     &      *Pa22*Sqrt2*Zh22*aa6*sb+Pa21**2*Sqrt2*Zh23*aa5
     &      *cb*sb+Pa21**2*Sqrt2*Zh23*aa6*cb*sb+Pa21**2*Sqrt2
     &      *Zh23*la1s*sb**2*xvev-Pa21**2*Sqrt2*Zh23*la2s
     &      *sb**2*xvev+4D0*MW/EE*Pa21*Pa22*SW*Zh23*la5s
     &      -4D0*MW/EE*Pa21*Pa22*SW*Zh23*la6s+2D0*MW/EE*Pa21
     &      **2*SW*Zh21*cb*la3+2D0*MW/EE*Pa21**2*SW*Zh21
     &      *cb*la4-2D0*MW/EE*Pa21**2*SW*Zh21*cb*la5+4D0
     &      *MW/EE*Pa21**2*SW*Zh21*la6*sb-6D0*MW/EE*Pa21
     &      **2*SW*Zh21*la6*sb**3-2D0*MW/EE*Pa21**2*SW*Zh21
     &      *la7*sb+2D0*MW/EE*Pa21**2*SW*Zh21*la7*sb**3-2D0
     &      *MW/EE*Pa21**2*SW*Zh22*cb*la7+2D0*MW/EE*Pa21
     &      **2*SW*Zh22*la2*sb-2D0*MW/EE*Pa21**2*SW*Zh22
     &      *la2*sb**3+2D0*MW/EE*Pa21**2*SW*Zh22*la3*sb**3
     &      +2D0*MW/EE*Pa21**2*SW*Zh22*la4*sb**3-4D0*MW/EE
     &      *Pa21**2*SW*Zh22*la5*sb+2D0*MW/EE*Pa21**2*SW
     &      *Zh22*la5*sb**3+2D0*MW/EE*Pa22**2*SW*Zh21*cb
     &      *la1s-2D0*MW/EE*Pa22**2*SW*Zh21*la4s*sb+2D0*MW
     &      /EE*Pa22**2*SW*Zh21*la5s*sb+2D0*MW/EE*Pa22**2
     &      *SW*Zh21*la6s*sb-2D0*MW/EE*Pa22**2*SW*Zh22*cb
     &      *la4s+2D0*MW/EE*Pa22**2*SW*Zh22*cb*la5s+2D0*MW
     &      /EE*Pa22**2*SW*Zh22*cb*la6s+2D0*MW/EE*Pa22**2
     &      *SW*Zh22*la2s*sb+2D0*Pa21*Pa22*Sqrt2*Zh21*cb
     &      *la5s*xvev-2D0*Pa21*Pa22*Sqrt2*Zh21*cb*la6s*xvev
     &      +2D0*Pa21*Pa22*Sqrt2*Zh22*la5s*sb*xvev-2D0*Pa21
     &      *Pa22*Sqrt2*Zh22*la6s*sb*xvev+2D0*Pa21**2*Sqrt2
     &      *Zh23*cb*la4s*sb*xvev+2D0*Pa21**2*Sqrt2*Zh23
     &      *cb*la5s*sb*xvev+2D0*Pa21**2*Sqrt2*Zh23*cb*la6s
     &      *sb*xvev+2D0*MW/EE*Pa21**2*SW*Zh21*cb*la1*sb
     &      **2-2D0*MW/EE*Pa21**2*SW*Zh21*cb*la3*sb**2-2D0
     &      *MW/EE*Pa21**2*SW*Zh21*cb*la4*sb**2-2D0*MW/EE
     &      *Pa21**2*SW*Zh21*cb*la5*sb**2-2D0*MW/EE*Pa21
     &      **2*SW*Zh22*cb*la6*sb**2+6D0*MW/EE*Pa21**2*SW
     &      *Zh22*cb*la7*sb**2
      AAABR(161) = 2D0*Pa11*Pa21*Sqrt2*Zh33*aa2+2D0*Pa12
     &      *Pa22*Sqrt2*Zh33*aa3-6D0*Pa12*Pa22*Sqrt2*Zh33
     &      *aa4+2D0*Pa11*Pa21*Sqrt2*Zh33*aa1*sb**2-2D0*Pa11
     &      *Pa21*Sqrt2*Zh33*aa2*sb**2+2D0*Pa11*Pa21*Sqrt2
     &      *Zh33*la2s*xvev+Pa11*Pa22*Sqrt2*Zh31*aa5*cb-Pa11
     &      *Pa22*Sqrt2*Zh31*aa6*cb+Pa11*Pa22*Sqrt2*Zh32
     &      *aa5*sb-Pa11*Pa22*Sqrt2*Zh32*aa6*sb+Pa12*Pa21
     &      *Sqrt2*Zh31*aa5*cb-Pa12*Pa21*Sqrt2*Zh31*aa6*cb
     &      +Pa12*Pa21*Sqrt2*Zh32*aa5*sb-Pa12*Pa21*Sqrt2
     &      *Zh32*aa6*sb+4D0*Pa12*Pa22*Sqrt2*Zh33*la3s*xvev
     &      -24D0*Pa12*Pa22*Sqrt2*Zh33*la7s*xvev+4D0*MW/EE
     &      *Pa11*Pa22*SW*Zh33*la5s-4D0*MW/EE*Pa11*Pa22*SW
     &      *Zh33*la6s+4D0*MW/EE*Pa12*Pa21*SW*Zh33*la5s-4D0
     &      *MW/EE*Pa12*Pa21*SW*Zh33*la6s+2D0*Pa11*Pa21*Sqrt2
     &      *Zh33*aa5*cb*sb+2D0*Pa11*Pa21*Sqrt2*Zh33*aa6
     &      *cb*sb+2D0*Pa11*Pa21*Sqrt2*Zh33*la1s*sb**2*xvev
     &      -2D0*Pa11*Pa21*Sqrt2*Zh33*la2s*sb**2*xvev+2D0
     &      *Pa11*Pa22*Sqrt2*Zh31*cb*la5s*xvev-2D0*Pa11*Pa22
     &      *Sqrt2*Zh31*cb*la6s*xvev+2D0*Pa11*Pa22*Sqrt2
     &      *Zh32*la5s*sb*xvev-2D0*Pa11*Pa22*Sqrt2*Zh32*la6s
     &      *sb*xvev+2D0*Pa12*Pa21*Sqrt2*Zh31*cb*la5s*xvev
     &      -2D0*Pa12*Pa21*Sqrt2*Zh31*cb*la6s*xvev+2D0*Pa12
     &      *Pa21*Sqrt2*Zh32*la5s*sb*xvev-2D0*Pa12*Pa21*Sqrt2
     &      *Zh32*la6s*sb*xvev+4D0*MW/EE*Pa11*Pa21*SW*Zh31
     &      *cb*la3+4D0*MW/EE*Pa11*Pa21*SW*Zh31*cb*la4-4D0
     &      *MW/EE*Pa11*Pa21*SW*Zh31*cb*la5+8D0*MW/EE*Pa11
     &      *Pa21*SW*Zh31*la6*sb-12D0*MW/EE*Pa11*Pa21*SW
     &      *Zh31*la6*sb**3-4D0*MW/EE*Pa11*Pa21*SW*Zh31*la7
     &      *sb+4D0*MW/EE*Pa11*Pa21*SW*Zh31*la7*sb**3-4D0
     &      *MW/EE*Pa11*Pa21*SW*Zh32*cb*la7+4D0*MW/EE*Pa11
     &      *Pa21*SW*Zh32*la2*sb-4D0*MW/EE*Pa11*Pa21*SW*Zh32
     &      *la2*sb**3+4D0*MW/EE*Pa11*Pa21*SW*Zh32*la3*sb
     &      **3+4D0*MW/EE*Pa11*Pa21*SW*Zh32*la4*sb**3-8D0
     &      *MW/EE*Pa11*Pa21*SW*Zh32*la5*sb+4D0*MW/EE*Pa11
     &      *Pa21*SW*Zh32*la5*sb**3+4D0*MW/EE*Pa12*Pa22*SW
     &      *Zh31*cb*la1s-4D0*MW/EE*Pa12*Pa22*SW*Zh31*la4s
     &      *sb+4D0*MW/EE*Pa12*Pa22*SW*Zh31*la5s*sb+4D0*MW
     &      /EE*Pa12*Pa22*SW*Zh31*la6s*sb-4D0*MW/EE*Pa12
     &      *Pa22*SW*Zh32*cb*la4s+4D0*MW/EE*Pa12*Pa22*SW
     &      *Zh32*cb*la5s+4D0*MW/EE*Pa12*Pa22*SW*Zh32*cb
     &      *la6s+4D0*MW/EE*Pa12*Pa22*SW*Zh32*la2s*sb+4D0
     &      *Pa11*Pa21*Sqrt2*Zh33*cb*la4s*sb*xvev+4D0*Pa11
     &      *Pa21*Sqrt2*Zh33*cb*la5s*sb*xvev+4D0*Pa11*Pa21
     &      *Sqrt2*Zh33*cb*la6s*sb*xvev+4D0*MW/EE*Pa11*Pa21
     &      *SW*Zh31*cb*la1*sb**2-4D0*MW/EE*Pa11*Pa21*SW
     &      *Zh31*cb*la3*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh31
     &      *cb*la4*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh31*cb
     &      *la5*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh32*cb*la6
     &      *sb**2+12D0*MW/EE*Pa11*Pa21*SW*Zh32*cb*la7*sb**2
      AAABR(162) = Pa21**2*Sqrt2*Zh33*aa2+Pa22**2*Sqrt2*Zh33
     &      *aa3-3D0*Pa22**2*Sqrt2*Zh33*aa4+Pa21**2*Sqrt2
     &      *Zh33*aa1*sb**2-Pa21**2*Sqrt2*Zh33*aa2*sb**2
     &      +Pa21**2*Sqrt2*Zh33*la2s*xvev+2D0*Pa22**2*Sqrt2
     &      *Zh33*la3s*xvev-12D0*Pa22**2*Sqrt2*Zh33*la7s
     &      *xvev+Pa21*Pa22*Sqrt2*Zh31*aa5*cb-Pa21*Pa22*Sqrt2
     &      *Zh31*aa6*cb+Pa21*Pa22*Sqrt2*Zh32*aa5*sb-Pa21
     &      *Pa22*Sqrt2*Zh32*aa6*sb+Pa21**2*Sqrt2*Zh33*aa5
     &      *cb*sb+Pa21**2*Sqrt2*Zh33*aa6*cb*sb+Pa21**2*Sqrt2
     &      *Zh33*la1s*sb**2*xvev-Pa21**2*Sqrt2*Zh33*la2s
     &      *sb**2*xvev+4D0*MW/EE*Pa21*Pa22*SW*Zh33*la5s
     &      -4D0*MW/EE*Pa21*Pa22*SW*Zh33*la6s+2D0*MW/EE*Pa21
     &      **2*SW*Zh31*cb*la3+2D0*MW/EE*Pa21**2*SW*Zh31
     &      *cb*la4-2D0*MW/EE*Pa21**2*SW*Zh31*cb*la5+4D0
     &      *MW/EE*Pa21**2*SW*Zh31*la6*sb-6D0*MW/EE*Pa21
     &      **2*SW*Zh31*la6*sb**3-2D0*MW/EE*Pa21**2*SW*Zh31
     &      *la7*sb+2D0*MW/EE*Pa21**2*SW*Zh31*la7*sb**3-2D0
     &      *MW/EE*Pa21**2*SW*Zh32*cb*la7+2D0*MW/EE*Pa21
     &      **2*SW*Zh32*la2*sb-2D0*MW/EE*Pa21**2*SW*Zh32
     &      *la2*sb**3+2D0*MW/EE*Pa21**2*SW*Zh32*la3*sb**3
     &      +2D0*MW/EE*Pa21**2*SW*Zh32*la4*sb**3-4D0*MW/EE
     &      *Pa21**2*SW*Zh32*la5*sb+2D0*MW/EE*Pa21**2*SW
     &      *Zh32*la5*sb**3+2D0*MW/EE*Pa22**2*SW*Zh31*cb
     &      *la1s-2D0*MW/EE*Pa22**2*SW*Zh31*la4s*sb+2D0*MW
     &      /EE*Pa22**2*SW*Zh31*la5s*sb+2D0*MW/EE*Pa22**2
     &      *SW*Zh31*la6s*sb-2D0*MW/EE*Pa22**2*SW*Zh32*cb
     &      *la4s+2D0*MW/EE*Pa22**2*SW*Zh32*cb*la5s+2D0*MW
     &      /EE*Pa22**2*SW*Zh32*cb*la6s+2D0*MW/EE*Pa22**2
     &      *SW*Zh32*la2s*sb+2D0*Pa21*Pa22*Sqrt2*Zh31*cb
     &      *la5s*xvev-2D0*Pa21*Pa22*Sqrt2*Zh31*cb*la6s*xvev
     &      +2D0*Pa21*Pa22*Sqrt2*Zh32*la5s*sb*xvev-2D0*Pa21
     &      *Pa22*Sqrt2*Zh32*la6s*sb*xvev+2D0*Pa21**2*Sqrt2
     &      *Zh33*cb*la4s*sb*xvev+2D0*Pa21**2*Sqrt2*Zh33
     &      *cb*la5s*sb*xvev+2D0*Pa21**2*Sqrt2*Zh33*cb*la6s
     &      *sb*xvev+2D0*MW/EE*Pa21**2*SW*Zh31*cb*la1*sb
     &      **2-2D0*MW/EE*Pa21**2*SW*Zh31*cb*la3*sb**2-2D0
     &      *MW/EE*Pa21**2*SW*Zh31*cb*la4*sb**2-2D0*MW/EE
     &      *Pa21**2*SW*Zh31*cb*la5*sb**2-2D0*MW/EE*Pa21
     &      **2*SW*Zh32*cb*la6*sb**2+6D0*MW/EE*Pa21**2*SW
     &      *Zh32*cb*la7*sb**2
      AAABR(163) = 2D0*Pa11*Pa21*Sqrt2*Zh13*aa2+2D0*Pa12
     &      *Pa22*Sqrt2*Zh13*aa3-6D0*Pa12*Pa22*Sqrt2*Zh13
     &      *aa4+2D0*Pa11*Pa21*Sqrt2*Zh13*aa1*sb**2-2D0*Pa11
     &      *Pa21*Sqrt2*Zh13*aa2*sb**2+2D0*Pa11*Pa21*Sqrt2
     &      *Zh13*la2s*xvev+Pa11*Pa22*Sqrt2*Zh11*aa5*cb-Pa11
     &      *Pa22*Sqrt2*Zh11*aa6*cb+Pa11*Pa22*Sqrt2*Zh12
     &      *aa5*sb-Pa11*Pa22*Sqrt2*Zh12*aa6*sb+Pa12*Pa21
     &      *Sqrt2*Zh11*aa5*cb-Pa12*Pa21*Sqrt2*Zh11*aa6*cb
     &      +Pa12*Pa21*Sqrt2*Zh12*aa5*sb-Pa12*Pa21*Sqrt2
     &      *Zh12*aa6*sb+4D0*Pa12*Pa22*Sqrt2*Zh13*la3s*xvev
     &      -24D0*Pa12*Pa22*Sqrt2*Zh13*la7s*xvev+4D0*MW/EE
     &      *Pa11*Pa22*SW*Zh13*la5s-4D0*MW/EE*Pa11*Pa22*SW
     &      *Zh13*la6s+4D0*MW/EE*Pa12*Pa21*SW*Zh13*la5s-4D0
     &      *MW/EE*Pa12*Pa21*SW*Zh13*la6s+2D0*Pa11*Pa21*Sqrt2
     &      *Zh13*aa5*cb*sb+2D0*Pa11*Pa21*Sqrt2*Zh13*aa6
     &      *cb*sb+2D0*Pa11*Pa21*Sqrt2*Zh13*la1s*sb**2*xvev
     &      -2D0*Pa11*Pa21*Sqrt2*Zh13*la2s*sb**2*xvev+2D0
     &      *Pa11*Pa22*Sqrt2*Zh11*cb*la5s*xvev-2D0*Pa11*Pa22
     &      *Sqrt2*Zh11*cb*la6s*xvev+2D0*Pa11*Pa22*Sqrt2
     &      *Zh12*la5s*sb*xvev-2D0*Pa11*Pa22*Sqrt2*Zh12*la6s
     &      *sb*xvev+2D0*Pa12*Pa21*Sqrt2*Zh11*cb*la5s*xvev
     &      -2D0*Pa12*Pa21*Sqrt2*Zh11*cb*la6s*xvev+2D0*Pa12
     &      *Pa21*Sqrt2*Zh12*la5s*sb*xvev-2D0*Pa12*Pa21*Sqrt2
     &      *Zh12*la6s*sb*xvev+4D0*MW/EE*Pa11*Pa21*SW*Zh11
     &      *cb*la3+4D0*MW/EE*Pa11*Pa21*SW*Zh11*cb*la4-4D0
     &      *MW/EE*Pa11*Pa21*SW*Zh11*cb*la5+8D0*MW/EE*Pa11
     &      *Pa21*SW*Zh11*la6*sb-12D0*MW/EE*Pa11*Pa21*SW
     &      *Zh11*la6*sb**3-4D0*MW/EE*Pa11*Pa21*SW*Zh11*la7
     &      *sb+4D0*MW/EE*Pa11*Pa21*SW*Zh11*la7*sb**3-4D0
     &      *MW/EE*Pa11*Pa21*SW*Zh12*cb*la7+4D0*MW/EE*Pa11
     &      *Pa21*SW*Zh12*la2*sb-4D0*MW/EE*Pa11*Pa21*SW*Zh12
     &      *la2*sb**3+4D0*MW/EE*Pa11*Pa21*SW*Zh12*la3*sb
     &      **3+4D0*MW/EE*Pa11*Pa21*SW*Zh12*la4*sb**3-8D0
     &      *MW/EE*Pa11*Pa21*SW*Zh12*la5*sb+4D0*MW/EE*Pa11
     &      *Pa21*SW*Zh12*la5*sb**3+4D0*MW/EE*Pa12*Pa22*SW
     &      *Zh11*cb*la1s-4D0*MW/EE*Pa12*Pa22*SW*Zh11*la4s
     &      *sb+4D0*MW/EE*Pa12*Pa22*SW*Zh11*la5s*sb+4D0*MW
     &      /EE*Pa12*Pa22*SW*Zh11*la6s*sb-4D0*MW/EE*Pa12
     &      *Pa22*SW*Zh12*cb*la4s+4D0*MW/EE*Pa12*Pa22*SW
     &      *Zh12*cb*la5s+4D0*MW/EE*Pa12*Pa22*SW*Zh12*cb
     &      *la6s+4D0*MW/EE*Pa12*Pa22*SW*Zh12*la2s*sb+4D0
     &      *Pa11*Pa21*Sqrt2*Zh13*cb*la4s*sb*xvev+4D0*Pa11
     &      *Pa21*Sqrt2*Zh13*cb*la5s*sb*xvev+4D0*Pa11*Pa21
     &      *Sqrt2*Zh13*cb*la6s*sb*xvev+4D0*MW/EE*Pa11*Pa21
     &      *SW*Zh11*cb*la1*sb**2-4D0*MW/EE*Pa11*Pa21*SW
     &      *Zh11*cb*la3*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh11
     &      *cb*la4*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh11*cb
     &      *la5*sb**2-4D0*MW/EE*Pa11*Pa21*SW*Zh12*cb*la6
     &      *sb**2+12D0*MW/EE*Pa11*Pa21*SW*Zh12*cb*la7*sb**2
      AAABR(164) = Pa11**2*Sqrt2*Zh13*aa2+Pa12**2*Sqrt2*Zh13
     &      *aa3-3D0*Pa12**2*Sqrt2*Zh13*aa4+Pa11**2*Sqrt2
     &      *Zh13*aa1*sb**2-Pa11**2*Sqrt2*Zh13*aa2*sb**2
     &      +Pa11**2*Sqrt2*Zh13*la2s*xvev+2D0*Pa12**2*Sqrt2
     &      *Zh13*la3s*xvev-12D0*Pa12**2*Sqrt2*Zh13*la7s
     &      *xvev+Pa11*Pa12*Sqrt2*Zh11*aa5*cb-Pa11*Pa12*Sqrt2
     &      *Zh11*aa6*cb+Pa11*Pa12*Sqrt2*Zh12*aa5*sb-Pa11
     &      *Pa12*Sqrt2*Zh12*aa6*sb+Pa11**2*Sqrt2*Zh13*aa5
     &      *cb*sb+Pa11**2*Sqrt2*Zh13*aa6*cb*sb+Pa11**2*Sqrt2
     &      *Zh13*la1s*sb**2*xvev-Pa11**2*Sqrt2*Zh13*la2s
     &      *sb**2*xvev+4D0*MW/EE*Pa11*Pa12*SW*Zh13*la5s
     &      -4D0*MW/EE*Pa11*Pa12*SW*Zh13*la6s+2D0*MW/EE*Pa11
     &      **2*SW*Zh11*cb*la3+2D0*MW/EE*Pa11**2*SW*Zh11
     &      *cb*la4-2D0*MW/EE*Pa11**2*SW*Zh11*cb*la5+4D0
     &      *MW/EE*Pa11**2*SW*Zh11*la6*sb-6D0*MW/EE*Pa11
     &      **2*SW*Zh11*la6*sb**3-2D0*MW/EE*Pa11**2*SW*Zh11
     &      *la7*sb+2D0*MW/EE*Pa11**2*SW*Zh11*la7*sb**3-2D0
     &      *MW/EE*Pa11**2*SW*Zh12*cb*la7+2D0*MW/EE*Pa11
     &      **2*SW*Zh12*la2*sb-2D0*MW/EE*Pa11**2*SW*Zh12
     &      *la2*sb**3+2D0*MW/EE*Pa11**2*SW*Zh12*la3*sb**3
     &      +2D0*MW/EE*Pa11**2*SW*Zh12*la4*sb**3-4D0*MW/EE
     &      *Pa11**2*SW*Zh12*la5*sb+2D0*MW/EE*Pa11**2*SW
     &      *Zh12*la5*sb**3+2D0*MW/EE*Pa12**2*SW*Zh11*cb
     &      *la1s-2D0*MW/EE*Pa12**2*SW*Zh11*la4s*sb+2D0*MW
     &      /EE*Pa12**2*SW*Zh11*la5s*sb+2D0*MW/EE*Pa12**2
     &      *SW*Zh11*la6s*sb-2D0*MW/EE*Pa12**2*SW*Zh12*cb
     &      *la4s+2D0*MW/EE*Pa12**2*SW*Zh12*cb*la5s+2D0*MW
     &      /EE*Pa12**2*SW*Zh12*cb*la6s+2D0*MW/EE*Pa12**2
     &      *SW*Zh12*la2s*sb+2D0*Pa11*Pa12*Sqrt2*Zh11*cb
     &      *la5s*xvev-2D0*Pa11*Pa12*Sqrt2*Zh11*cb*la6s*xvev
     &      +2D0*Pa11*Pa12*Sqrt2*Zh12*la5s*sb*xvev-2D0*Pa11
     &      *Pa12*Sqrt2*Zh12*la6s*sb*xvev+2D0*Pa11**2*Sqrt2
     &      *Zh13*cb*la4s*sb*xvev+2D0*Pa11**2*Sqrt2*Zh13
     &      *cb*la5s*sb*xvev+2D0*Pa11**2*Sqrt2*Zh13*cb*la6s
     &      *sb*xvev+2D0*MW/EE*Pa11**2*SW*Zh11*cb*la1*sb
     &      **2-2D0*MW/EE*Pa11**2*SW*Zh11*cb*la3*sb**2-2D0
     &      *MW/EE*Pa11**2*SW*Zh11*cb*la4*sb**2-2D0*MW/EE
     &      *Pa11**2*SW*Zh11*cb*la5*sb**2-2D0*MW/EE*Pa11
     &      **2*SW*Zh12*cb*la6*sb**2+6D0*MW/EE*Pa11**2*SW
     &      *Zh12*cb*la7*sb**2
      AAABR(165) = Sqrt2*Zh11**2*Zh23*aa1+Sqrt2*Zh12**2*Zh23
     &      *aa2+3D0*Sqrt2*Zh13**2*Zh23*aa3+3D0*Sqrt2*Zh13
     &      **2*Zh23*aa4-Sqrt2*Zh11*Zh12*Zh23*aa5-Sqrt2*Zh11
     &      *Zh12*Zh23*aa6+2D0*Sqrt2*Zh11*Zh13*Zh21*aa1-Sqrt2
     &      *Zh11*Zh13*Zh22*aa5-Sqrt2*Zh11*Zh13*Zh22*aa6
     &      +Sqrt2*Zh11**2*Zh23*la1s*xvev-Sqrt2*Zh12*Zh13
     &      *Zh21*aa5-Sqrt2*Zh12*Zh13*Zh21*aa6+2D0*Sqrt2
     &      *Zh12*Zh13*Zh22*aa2+Sqrt2*Zh12**2*Zh23*la2s*xvev
     &      +6D0*Sqrt2*Zh13**2*Zh23*la3s*xvev+12D0*Sqrt2
     &      *Zh13**2*Zh23*la7s*xvev+12D0*Sqrt2*Zh13**2*Zh23
     &      *la8s*xvev-2D0*Sqrt2*Zh11*Zh12*Zh23*la4s*xvev
     &      -2D0*Sqrt2*Zh11*Zh12*Zh23*la5s*xvev-2D0*Sqrt2
     &      *Zh11*Zh12*Zh23*la6s*xvev+2D0*Sqrt2*Zh11*Zh13
     &      *Zh21*la1s*xvev-2D0*Sqrt2*Zh11*Zh13*Zh22*la4s
     &      *xvev-2D0*Sqrt2*Zh11*Zh13*Zh22*la5s*xvev-2D0
     &      *Sqrt2*Zh11*Zh13*Zh22*la6s*xvev-2D0*Sqrt2*Zh12
     &      *Zh13*Zh21*la4s*xvev-2D0*Sqrt2*Zh12*Zh13*Zh21
     &      *la5s*xvev-2D0*Sqrt2*Zh12*Zh13*Zh21*la6s*xvev
     &      +2D0*Sqrt2*Zh12*Zh13*Zh22*la2s*xvev+6D0*MW/EE
     &      *SW*Zh11**2*Zh21*cb*la1-6D0*MW/EE*SW*Zh11**2
     &      *Zh21*la6*sb-6D0*MW/EE*SW*Zh11**2*Zh22*cb*la6
     &      +2D0*MW/EE*SW*Zh11**2*Zh22*la3*sb+2D0*MW/EE*SW
     &      *Zh11**2*Zh22*la4*sb+2D0*MW/EE*SW*Zh11**2*Zh22
     &      *la5*sb+2D0*MW/EE*SW*Zh12**2*Zh21*cb*la3+2D0
     &      *MW/EE*SW*Zh12**2*Zh21*cb*la4+2D0*MW/EE*SW*Zh12
     &      **2*Zh21*cb*la5-6D0*MW/EE*SW*Zh12**2*Zh21*la7
     &      *sb-6D0*MW/EE*SW*Zh12**2*Zh22*cb*la7+6D0*MW/EE
     &      *SW*Zh12**2*Zh22*la2*sb+2D0*MW/EE*SW*Zh13**2
     &      *Zh21*cb*la1s-2D0*MW/EE*SW*Zh13**2*Zh21*la4s
     &      *sb-2D0*MW/EE*SW*Zh13**2*Zh21*la5s*sb-2D0*MW
     &      /EE*SW*Zh13**2*Zh21*la6s*sb-2D0*MW/EE*SW*Zh13
     &      **2*Zh22*cb*la4s-2D0*MW/EE*SW*Zh13**2*Zh22*cb
     &      *la5s-2D0*MW/EE*SW*Zh13**2*Zh22*cb*la6s+2D0*MW
     &      /EE*SW*Zh13**2*Zh22*la2s*sb-12D0*MW/EE*SW*Zh11
     &      *Zh12*Zh21*cb*la6+4D0*MW/EE*SW*Zh11*Zh12*Zh21
     &      *la3*sb+4D0*MW/EE*SW*Zh11*Zh12*Zh21*la4*sb+4D0
     &      *MW/EE*SW*Zh11*Zh12*Zh21*la5*sb+4D0*MW/EE*SW
     &      *Zh11*Zh12*Zh22*cb*la3+4D0*MW/EE*SW*Zh11*Zh12
     &      *Zh22*cb*la4+4D0*MW/EE*SW*Zh11*Zh12*Zh22*cb*la5
     &      -12D0*MW/EE*SW*Zh11*Zh12*Zh22*la7*sb+4D0*MW/EE
     &      *SW*Zh11*Zh13*Zh23*cb*la1s-4D0*MW/EE*SW*Zh11
     &      *Zh13*Zh23*la4s*sb-4D0*MW/EE*SW*Zh11*Zh13*Zh23
     &      *la5s*sb-4D0*MW/EE*SW*Zh11*Zh13*Zh23*la6s*sb
     &      -4D0*MW/EE*SW*Zh12*Zh13*Zh23*cb*la4s-4D0*MW/EE
     &      *SW*Zh12*Zh13*Zh23*cb*la5s-4D0*MW/EE*SW*Zh12
     &      *Zh13*Zh23*cb*la6s+4D0*MW/EE*SW*Zh12*Zh13*Zh23
     &      *la2s*sb
      AAABR(166) = Sqrt2*Zh13*Zh21**2*aa1+Sqrt2*Zh13*Zh22
     &      **2*aa2+3D0*Sqrt2*Zh13*Zh23**2*aa3+3D0*Sqrt2
     &      *Zh13*Zh23**2*aa4+2D0*Sqrt2*Zh11*Zh21*Zh23*aa1
     &      -Sqrt2*Zh11*Zh22*Zh23*aa5-Sqrt2*Zh11*Zh22*Zh23
     &      *aa6-Sqrt2*Zh12*Zh21*Zh23*aa5-Sqrt2*Zh12*Zh21
     &      *Zh23*aa6+2D0*Sqrt2*Zh12*Zh22*Zh23*aa2-Sqrt2
     &      *Zh13*Zh21*Zh22*aa5-Sqrt2*Zh13*Zh21*Zh22*aa6
     &      +Sqrt2*Zh13*Zh21**2*la1s*xvev+Sqrt2*Zh13*Zh22
     &      **2*la2s*xvev+6D0*Sqrt2*Zh13*Zh23**2*la3s*xvev
     &      +12D0*Sqrt2*Zh13*Zh23**2*la7s*xvev+12D0*Sqrt2
     &      *Zh13*Zh23**2*la8s*xvev+2D0*Sqrt2*Zh11*Zh21*Zh23
     &      *la1s*xvev-2D0*Sqrt2*Zh11*Zh22*Zh23*la4s*xvev
     &      -2D0*Sqrt2*Zh11*Zh22*Zh23*la5s*xvev-2D0*Sqrt2
     &      *Zh11*Zh22*Zh23*la6s*xvev-2D0*Sqrt2*Zh12*Zh21
     &      *Zh23*la4s*xvev-2D0*Sqrt2*Zh12*Zh21*Zh23*la5s
     &      *xvev-2D0*Sqrt2*Zh12*Zh21*Zh23*la6s*xvev+2D0
     &      *Sqrt2*Zh12*Zh22*Zh23*la2s*xvev-2D0*Sqrt2*Zh13
     &      *Zh21*Zh22*la4s*xvev-2D0*Sqrt2*Zh13*Zh21*Zh22
     &      *la5s*xvev-2D0*Sqrt2*Zh13*Zh21*Zh22*la6s*xvev
     &      +6D0*MW/EE*SW*Zh11*Zh21**2*cb*la1-6D0*MW/EE*SW
     &      *Zh11*Zh21**2*la6*sb+2D0*MW/EE*SW*Zh11*Zh22**2
     &      *cb*la3+2D0*MW/EE*SW*Zh11*Zh22**2*cb*la4+2D0
     &      *MW/EE*SW*Zh11*Zh22**2*cb*la5-6D0*MW/EE*SW*Zh11
     &      *Zh22**2*la7*sb+2D0*MW/EE*SW*Zh11*Zh23**2*cb
     &      *la1s-2D0*MW/EE*SW*Zh11*Zh23**2*la4s*sb-2D0*MW
     &      /EE*SW*Zh11*Zh23**2*la5s*sb-2D0*MW/EE*SW*Zh11
     &      *Zh23**2*la6s*sb-6D0*MW/EE*SW*Zh12*Zh21**2*cb
     &      *la6+2D0*MW/EE*SW*Zh12*Zh21**2*la3*sb+2D0*MW
     &      /EE*SW*Zh12*Zh21**2*la4*sb+2D0*MW/EE*SW*Zh12
     &      *Zh21**2*la5*sb-6D0*MW/EE*SW*Zh12*Zh22**2*cb
     &      *la7+6D0*MW/EE*SW*Zh12*Zh22**2*la2*sb-2D0*MW
     &      /EE*SW*Zh12*Zh23**2*cb*la4s-2D0*MW/EE*SW*Zh12
     &      *Zh23**2*cb*la5s-2D0*MW/EE*SW*Zh12*Zh23**2*cb
     &      *la6s+2D0*MW/EE*SW*Zh12*Zh23**2*la2s*sb-12D0
     &      *MW/EE*SW*Zh11*Zh21*Zh22*cb*la6+4D0*MW/EE*SW
     &      *Zh11*Zh21*Zh22*la3*sb+4D0*MW/EE*SW*Zh11*Zh21
     &      *Zh22*la4*sb+4D0*MW/EE*SW*Zh11*Zh21*Zh22*la5
     &      *sb+4D0*MW/EE*SW*Zh12*Zh21*Zh22*cb*la3+4D0*MW
     &      /EE*SW*Zh12*Zh21*Zh22*cb*la4+4D0*MW/EE*SW*Zh12
     &      *Zh21*Zh22*cb*la5-12D0*MW/EE*SW*Zh12*Zh21*Zh22
     &      *la7*sb+4D0*MW/EE*SW*Zh13*Zh21*Zh23*cb*la1s-4D0
     &      *MW/EE*SW*Zh13*Zh21*Zh23*la4s*sb-4D0*MW/EE*SW
     &      *Zh13*Zh21*Zh23*la5s*sb-4D0*MW/EE*SW*Zh13*Zh21
     &      *Zh23*la6s*sb-4D0*MW/EE*SW*Zh13*Zh22*Zh23*cb
     &      *la4s-4D0*MW/EE*SW*Zh13*Zh22*Zh23*cb*la5s-4D0
     &      *MW/EE*SW*Zh13*Zh22*Zh23*cb*la6s+4D0*MW/EE*SW
     &      *Zh13*Zh22*Zh23*la2s*sb
      AAABR(167) = 2D0*Sqrt2*Zh11*Zh21*Zh33*aa1-Sqrt2*Zh11
     &      *Zh22*Zh33*aa5-Sqrt2*Zh11*Zh22*Zh33*aa6+2D0*Sqrt2
     &      *Zh11*Zh23*Zh31*aa1-Sqrt2*Zh11*Zh23*Zh32*aa5
     &      -Sqrt2*Zh11*Zh23*Zh32*aa6-Sqrt2*Zh12*Zh21*Zh33
     &      *aa5-Sqrt2*Zh12*Zh21*Zh33*aa6+2D0*Sqrt2*Zh12
     &      *Zh22*Zh33*aa2-Sqrt2*Zh12*Zh23*Zh31*aa5-Sqrt2
     &      *Zh12*Zh23*Zh31*aa6+2D0*Sqrt2*Zh12*Zh23*Zh32
     &      *aa2+2D0*Sqrt2*Zh13*Zh21*Zh31*aa1-Sqrt2*Zh13
     &      *Zh21*Zh32*aa5-Sqrt2*Zh13*Zh21*Zh32*aa6-Sqrt2
     &      *Zh13*Zh22*Zh31*aa5-Sqrt2*Zh13*Zh22*Zh31*aa6
     &      +2D0*Sqrt2*Zh13*Zh22*Zh32*aa2+6D0*Sqrt2*Zh13
     &      *Zh23*Zh33*aa3+6D0*Sqrt2*Zh13*Zh23*Zh33*aa4+2D0
     &      *Sqrt2*Zh11*Zh21*Zh33*la1s*xvev-2D0*Sqrt2*Zh11
     &      *Zh22*Zh33*la4s*xvev-2D0*Sqrt2*Zh11*Zh22*Zh33
     &      *la5s*xvev-2D0*Sqrt2*Zh11*Zh22*Zh33*la6s*xvev
     &      +2D0*Sqrt2*Zh11*Zh23*Zh31*la1s*xvev-2D0*Sqrt2
     &      *Zh11*Zh23*Zh32*la4s*xvev-2D0*Sqrt2*Zh11*Zh23
     &      *Zh32*la5s*xvev-2D0*Sqrt2*Zh11*Zh23*Zh32*la6s
     &      *xvev-2D0*Sqrt2*Zh12*Zh21*Zh33*la4s*xvev-2D0
     &      *Sqrt2*Zh12*Zh21*Zh33*la5s*xvev-2D0*Sqrt2*Zh12
     &      *Zh21*Zh33*la6s*xvev+2D0*Sqrt2*Zh12*Zh22*Zh33
     &      *la2s*xvev-2D0*Sqrt2*Zh12*Zh23*Zh31*la4s*xvev
     &      -2D0*Sqrt2*Zh12*Zh23*Zh31*la5s*xvev-2D0*Sqrt2
     &      *Zh12*Zh23*Zh31*la6s*xvev+2D0*Sqrt2*Zh12*Zh23
     &      *Zh32*la2s*xvev+2D0*Sqrt2*Zh13*Zh21*Zh31*la1s
     &      *xvev-2D0*Sqrt2*Zh13*Zh21*Zh32*la4s*xvev-2D0
     &      *Sqrt2*Zh13*Zh21*Zh32*la5s*xvev-2D0*Sqrt2*Zh13
     &      *Zh21*Zh32*la6s*xvev-2D0*Sqrt2*Zh13*Zh22*Zh31
     &      *la4s*xvev-2D0*Sqrt2*Zh13*Zh22*Zh31*la5s*xvev
     &      -2D0*Sqrt2*Zh13*Zh22*Zh31*la6s*xvev+2D0*Sqrt2
     &      *Zh13*Zh22*Zh32*la2s*xvev+12D0*Sqrt2*Zh13*Zh23
     &      *Zh33*la3s*xvev+24D0*Sqrt2*Zh13*Zh23*Zh33*la7s
     &      *xvev+24D0*Sqrt2*Zh13*Zh23*Zh33*la8s*xvev+12D0
     &      *MW/EE*SW*Zh11*Zh21*Zh31*cb*la1-12D0*MW/EE*SW
     &      *Zh11*Zh21*Zh31*la6*sb-12D0*MW/EE*SW*Zh11*Zh21
     &      *Zh32*cb*la6+4D0*MW/EE*SW*Zh11*Zh21*Zh32*la3
     &      *sb+4D0*MW/EE*SW*Zh11*Zh21*Zh32*la4*sb+4D0*MW
     &      /EE*SW*Zh11*Zh21*Zh32*la5*sb-12D0*MW/EE*SW*Zh11
     &      *Zh22*Zh31*cb*la6+4D0*MW/EE*SW*Zh11*Zh22*Zh31
     &      *la3*sb+4D0*MW/EE*SW*Zh11*Zh22*Zh31*la4*sb+4D0
     &      *MW/EE*SW*Zh11*Zh22*Zh31*la5*sb+4D0*MW/EE*SW
     &      *Zh11*Zh22*Zh32*cb*la3+4D0*MW/EE*SW*Zh11*Zh22
     &      *Zh32*cb*la4+4D0*MW/EE*SW*Zh11*Zh22*Zh32*cb*la5
     &      -12D0*MW/EE*SW*Zh11*Zh22*Zh32*la7*sb+4D0*MW/EE
     &      *SW*Zh11*Zh23*Zh33*cb*la1s-4D0*MW/EE*SW*Zh11
     &      *Zh23*Zh33*la4s*sb-4D0*MW/EE*SW*Zh11*Zh23*Zh33
     &      *la5s*sb-4D0*MW/EE*SW*Zh11*Zh23*Zh33*la6s*sb
     &      -12D0*MW/EE*SW*Zh12*Zh21*Zh31*cb*la6+4D0*MW/EE
     &      *SW*Zh12*Zh21*Zh31*la3*sb+4D0*MW/EE*SW*Zh12*Zh21
     &      *Zh31*la4*sb+4D0*MW/EE*SW*Zh12*Zh21*Zh31*la5
     &      *sb+4D0*MW/EE*SW*Zh12*Zh21*Zh32*cb*la3+4D0*MW
     &      /EE*SW*Zh12*Zh21*Zh32*cb*la4+4D0*MW/EE*SW*Zh12
     &      *Zh21*Zh32*cb*la5-12D0*MW/EE*SW*Zh12*Zh21*Zh32
     &      *la7*sb+4D0*MW/EE*SW*Zh12*Zh22*Zh31*cb*la3+4D0
     &      *MW/EE*SW*Zh12*Zh22*Zh31*cb*la4+4D0*MW/EE*SW
     &      *Zh12*Zh22*Zh31*cb*la5-12D0*MW/EE*SW*Zh12*Zh22
     &      *Zh31*la7*sb-12D0*MW/EE*SW*Zh12*Zh22*Zh32*cb
     &      *la7+12D0*MW/EE*SW*Zh12*Zh22*Zh32*la2*sb-4D0
     &      *MW/EE*SW*Zh12*Zh23*Zh33*cb*la4s-4D0*MW/EE*SW
     &      *Zh12*Zh23*Zh33*cb*la5s-4D0*MW/EE*SW*Zh12*Zh23
     &      *Zh33*cb*la6s+4D0*MW/EE*SW*Zh12*Zh23*Zh33*la2s
     &      *sb+4D0*MW/EE*SW*Zh13*Zh21*Zh33*cb*la1s-4D0*MW
     &      /EE*SW*Zh13*Zh21*Zh33*la4s*sb-4D0*MW/EE*SW*Zh13
     &      *Zh21*Zh33*la5s*sb-4D0*MW/EE*SW*Zh13*Zh21*Zh33
     &      *la6s*sb-4D0*MW/EE*SW*Zh13*Zh22*Zh33*cb*la4s
     &      -4D0*MW/EE*SW*Zh13*Zh22*Zh33*cb*la5s-4D0*MW/EE
     &      *SW*Zh13*Zh22*Zh33*cb*la6s+4D0*MW/EE*SW*Zh13
     &      *Zh22*Zh33*la2s*sb+4D0*MW/EE*SW*Zh13*Zh23*Zh31
     &      *cb*la1s-4D0*MW/EE*SW*Zh13*Zh23*Zh31*la4s*sb
     &      -4D0*MW/EE*SW*Zh13*Zh23*Zh31*la5s*sb-4D0*MW/EE
     &      *SW*Zh13*Zh23*Zh31*la6s*sb-4D0*MW/EE*SW*Zh13
     &      *Zh23*Zh32*cb*la4s-4D0*MW/EE*SW*Zh13*Zh23*Zh32
     &      *cb*la5s-4D0*MW/EE*SW*Zh13*Zh23*Zh32*cb*la6s
     &      +4D0*MW/EE*SW*Zh13*Zh23*Zh32*la2s*sb
      AAABR(168) = Sqrt2*Zh13*Zh31**2*aa1+Sqrt2*Zh13*Zh32
     &      **2*aa2+3D0*Sqrt2*Zh13*Zh33**2*aa3+3D0*Sqrt2
     &      *Zh13*Zh33**2*aa4+2D0*Sqrt2*Zh11*Zh31*Zh33*aa1
     &      -Sqrt2*Zh11*Zh32*Zh33*aa5-Sqrt2*Zh11*Zh32*Zh33
     &      *aa6-Sqrt2*Zh12*Zh31*Zh33*aa5-Sqrt2*Zh12*Zh31
     &      *Zh33*aa6+2D0*Sqrt2*Zh12*Zh32*Zh33*aa2-Sqrt2
     &      *Zh13*Zh31*Zh32*aa5-Sqrt2*Zh13*Zh31*Zh32*aa6
     &      +Sqrt2*Zh13*Zh31**2*la1s*xvev+Sqrt2*Zh13*Zh32
     &      **2*la2s*xvev+6D0*Sqrt2*Zh13*Zh33**2*la3s*xvev
     &      +12D0*Sqrt2*Zh13*Zh33**2*la7s*xvev+12D0*Sqrt2
     &      *Zh13*Zh33**2*la8s*xvev+2D0*Sqrt2*Zh11*Zh31*Zh33
     &      *la1s*xvev-2D0*Sqrt2*Zh11*Zh32*Zh33*la4s*xvev
     &      -2D0*Sqrt2*Zh11*Zh32*Zh33*la5s*xvev-2D0*Sqrt2
     &      *Zh11*Zh32*Zh33*la6s*xvev-2D0*Sqrt2*Zh12*Zh31
     &      *Zh33*la4s*xvev-2D0*Sqrt2*Zh12*Zh31*Zh33*la5s
     &      *xvev-2D0*Sqrt2*Zh12*Zh31*Zh33*la6s*xvev+2D0
     &      *Sqrt2*Zh12*Zh32*Zh33*la2s*xvev-2D0*Sqrt2*Zh13
     &      *Zh31*Zh32*la4s*xvev-2D0*Sqrt2*Zh13*Zh31*Zh32
     &      *la5s*xvev-2D0*Sqrt2*Zh13*Zh31*Zh32*la6s*xvev
     &      +6D0*MW/EE*SW*Zh11*Zh31**2*cb*la1-6D0*MW/EE*SW
     &      *Zh11*Zh31**2*la6*sb+2D0*MW/EE*SW*Zh11*Zh32**2
     &      *cb*la3+2D0*MW/EE*SW*Zh11*Zh32**2*cb*la4+2D0
     &      *MW/EE*SW*Zh11*Zh32**2*cb*la5-6D0*MW/EE*SW*Zh11
     &      *Zh32**2*la7*sb+2D0*MW/EE*SW*Zh11*Zh33**2*cb
     &      *la1s-2D0*MW/EE*SW*Zh11*Zh33**2*la4s*sb-2D0*MW
     &      /EE*SW*Zh11*Zh33**2*la5s*sb-2D0*MW/EE*SW*Zh11
     &      *Zh33**2*la6s*sb-6D0*MW/EE*SW*Zh12*Zh31**2*cb
     &      *la6+2D0*MW/EE*SW*Zh12*Zh31**2*la3*sb+2D0*MW
     &      /EE*SW*Zh12*Zh31**2*la4*sb+2D0*MW/EE*SW*Zh12
     &      *Zh31**2*la5*sb-6D0*MW/EE*SW*Zh12*Zh32**2*cb
     &      *la7+6D0*MW/EE*SW*Zh12*Zh32**2*la2*sb-2D0*MW
     &      /EE*SW*Zh12*Zh33**2*cb*la4s-2D0*MW/EE*SW*Zh12
     &      *Zh33**2*cb*la5s-2D0*MW/EE*SW*Zh12*Zh33**2*cb
     &      *la6s+2D0*MW/EE*SW*Zh12*Zh33**2*la2s*sb-12D0
     &      *MW/EE*SW*Zh11*Zh31*Zh32*cb*la6+4D0*MW/EE*SW
     &      *Zh11*Zh31*Zh32*la3*sb+4D0*MW/EE*SW*Zh11*Zh31
     &      *Zh32*la4*sb+4D0*MW/EE*SW*Zh11*Zh31*Zh32*la5
     &      *sb+4D0*MW/EE*SW*Zh12*Zh31*Zh32*cb*la3+4D0*MW
     &      /EE*SW*Zh12*Zh31*Zh32*cb*la4+4D0*MW/EE*SW*Zh12
     &      *Zh31*Zh32*cb*la5-12D0*MW/EE*SW*Zh12*Zh31*Zh32
     &      *la7*sb+4D0*MW/EE*SW*Zh13*Zh31*Zh33*cb*la1s-4D0
     &      *MW/EE*SW*Zh13*Zh31*Zh33*la4s*sb-4D0*MW/EE*SW
     &      *Zh13*Zh31*Zh33*la5s*sb-4D0*MW/EE*SW*Zh13*Zh31
     &      *Zh33*la6s*sb-4D0*MW/EE*SW*Zh13*Zh32*Zh33*cb
     &      *la4s-4D0*MW/EE*SW*Zh13*Zh32*Zh33*cb*la5s-4D0
     &      *MW/EE*SW*Zh13*Zh32*Zh33*cb*la6s+4D0*MW/EE*SW
     &      *Zh13*Zh32*Zh33*la2s*sb
      AAABR(169) = Sqrt2*Zh23*Zh31**2*aa1+Sqrt2*Zh23*Zh32
     &      **2*aa2+3D0*Sqrt2*Zh23*Zh33**2*aa3+3D0*Sqrt2
     &      *Zh23*Zh33**2*aa4+2D0*Sqrt2*Zh21*Zh31*Zh33*aa1
     &      -Sqrt2*Zh21*Zh32*Zh33*aa5-Sqrt2*Zh21*Zh32*Zh33
     &      *aa6-Sqrt2*Zh22*Zh31*Zh33*aa5-Sqrt2*Zh22*Zh31
     &      *Zh33*aa6+2D0*Sqrt2*Zh22*Zh32*Zh33*aa2-Sqrt2
     &      *Zh23*Zh31*Zh32*aa5-Sqrt2*Zh23*Zh31*Zh32*aa6
     &      +Sqrt2*Zh23*Zh31**2*la1s*xvev+Sqrt2*Zh23*Zh32
     &      **2*la2s*xvev+6D0*Sqrt2*Zh23*Zh33**2*la3s*xvev
     &      +12D0*Sqrt2*Zh23*Zh33**2*la7s*xvev+12D0*Sqrt2
     &      *Zh23*Zh33**2*la8s*xvev+2D0*Sqrt2*Zh21*Zh31*Zh33
     &      *la1s*xvev-2D0*Sqrt2*Zh21*Zh32*Zh33*la4s*xvev
     &      -2D0*Sqrt2*Zh21*Zh32*Zh33*la5s*xvev-2D0*Sqrt2
     &      *Zh21*Zh32*Zh33*la6s*xvev-2D0*Sqrt2*Zh22*Zh31
     &      *Zh33*la4s*xvev-2D0*Sqrt2*Zh22*Zh31*Zh33*la5s
     &      *xvev-2D0*Sqrt2*Zh22*Zh31*Zh33*la6s*xvev+2D0
     &      *Sqrt2*Zh22*Zh32*Zh33*la2s*xvev-2D0*Sqrt2*Zh23
     &      *Zh31*Zh32*la4s*xvev-2D0*Sqrt2*Zh23*Zh31*Zh32
     &      *la5s*xvev-2D0*Sqrt2*Zh23*Zh31*Zh32*la6s*xvev
     &      +6D0*MW/EE*SW*Zh21*Zh31**2*cb*la1-6D0*MW/EE*SW
     &      *Zh21*Zh31**2*la6*sb+2D0*MW/EE*SW*Zh21*Zh32**2
     &      *cb*la3+2D0*MW/EE*SW*Zh21*Zh32**2*cb*la4+2D0
     &      *MW/EE*SW*Zh21*Zh32**2*cb*la5-6D0*MW/EE*SW*Zh21
     &      *Zh32**2*la7*sb+2D0*MW/EE*SW*Zh21*Zh33**2*cb
     &      *la1s-2D0*MW/EE*SW*Zh21*Zh33**2*la4s*sb-2D0*MW
     &      /EE*SW*Zh21*Zh33**2*la5s*sb-2D0*MW/EE*SW*Zh21
     &      *Zh33**2*la6s*sb-6D0*MW/EE*SW*Zh22*Zh31**2*cb
     &      *la6+2D0*MW/EE*SW*Zh22*Zh31**2*la3*sb+2D0*MW
     &      /EE*SW*Zh22*Zh31**2*la4*sb+2D0*MW/EE*SW*Zh22
     &      *Zh31**2*la5*sb-6D0*MW/EE*SW*Zh22*Zh32**2*cb
     &      *la7+6D0*MW/EE*SW*Zh22*Zh32**2*la2*sb-2D0*MW
     &      /EE*SW*Zh22*Zh33**2*cb*la4s-2D0*MW/EE*SW*Zh22
     &      *Zh33**2*cb*la5s-2D0*MW/EE*SW*Zh22*Zh33**2*cb
     &      *la6s+2D0*MW/EE*SW*Zh22*Zh33**2*la2s*sb-12D0
     &      *MW/EE*SW*Zh21*Zh31*Zh32*cb*la6+4D0*MW/EE*SW
     &      *Zh21*Zh31*Zh32*la3*sb+4D0*MW/EE*SW*Zh21*Zh31
     &      *Zh32*la4*sb+4D0*MW/EE*SW*Zh21*Zh31*Zh32*la5
     &      *sb+4D0*MW/EE*SW*Zh22*Zh31*Zh32*cb*la3+4D0*MW
     &      /EE*SW*Zh22*Zh31*Zh32*cb*la4+4D0*MW/EE*SW*Zh22
     &      *Zh31*Zh32*cb*la5-12D0*MW/EE*SW*Zh22*Zh31*Zh32
     &      *la7*sb+4D0*MW/EE*SW*Zh23*Zh31*Zh33*cb*la1s-4D0
     &      *MW/EE*SW*Zh23*Zh31*Zh33*la4s*sb-4D0*MW/EE*SW
     &      *Zh23*Zh31*Zh33*la5s*sb-4D0*MW/EE*SW*Zh23*Zh31
     &      *Zh33*la6s*sb-4D0*MW/EE*SW*Zh23*Zh32*Zh33*cb
     &      *la4s-4D0*MW/EE*SW*Zh23*Zh32*Zh33*cb*la5s-4D0
     &      *MW/EE*SW*Zh23*Zh32*Zh33*cb*la6s+4D0*MW/EE*SW
     &      *Zh23*Zh32*Zh33*la2s*sb
      AAABR(170) = Sqrt2*Zh33**3*aa3+Sqrt2*Zh33**3*aa4+Sqrt2
     &      *Zh31**2*Zh33*aa1+Sqrt2*Zh32**2*Zh33*aa2+2D0
     &      *Sqrt2*Zh33**3*la3s*xvev+4D0*Sqrt2*Zh33**3*la7s
     &      *xvev+4D0*Sqrt2*Zh33**3*la8s*xvev-Sqrt2*Zh31
     &      *Zh32*Zh33*aa5-Sqrt2*Zh31*Zh32*Zh33*aa6+Sqrt2
     &      *Zh31**2*Zh33*la1s*xvev+Sqrt2*Zh32**2*Zh33*la2s
     &      *xvev+2D0*MW/EE*SW*Zh31**3*cb*la1-2D0*MW/EE*SW
     &      *Zh31**3*la6*sb-2D0*MW/EE*SW*Zh32**3*cb*la7+2D0
     &      *MW/EE*SW*Zh32**3*la2*sb-2D0*Sqrt2*Zh31*Zh32
     &      *Zh33*la4s*xvev-2D0*Sqrt2*Zh31*Zh32*Zh33*la5s
     &      *xvev-2D0*Sqrt2*Zh31*Zh32*Zh33*la6s*xvev+2D0
     &      *MW/EE*SW*Zh31*Zh32**2*cb*la3+2D0*MW/EE*SW*Zh31
     &      *Zh32**2*cb*la4+2D0*MW/EE*SW*Zh31*Zh32**2*cb
     &      *la5-6D0*MW/EE*SW*Zh31*Zh32**2*la7*sb+2D0*MW
     &      /EE*SW*Zh31*Zh33**2*cb*la1s-2D0*MW/EE*SW*Zh31
     &      *Zh33**2*la4s*sb-2D0*MW/EE*SW*Zh31*Zh33**2*la5s
     &      *sb-2D0*MW/EE*SW*Zh31*Zh33**2*la6s*sb-6D0*MW
     &      /EE*SW*Zh31**2*Zh32*cb*la6+2D0*MW/EE*SW*Zh31
     &      **2*Zh32*la3*sb+2D0*MW/EE*SW*Zh31**2*Zh32*la4
     &      *sb+2D0*MW/EE*SW*Zh31**2*Zh32*la5*sb-2D0*MW/EE
     &      *SW*Zh32*Zh33**2*cb*la4s-2D0*MW/EE*SW*Zh32*Zh33
     &      **2*cb*la5s-2D0*MW/EE*SW*Zh32*Zh33**2*cb*la6s
     &      +2D0*MW/EE*SW*Zh32*Zh33**2*la2s*sb
      AAABR(171) = Sqrt2*Zh13**3*aa3+Sqrt2*Zh13**3*aa4+Sqrt2
     &      *Zh11**2*Zh13*aa1+Sqrt2*Zh12**2*Zh13*aa2+2D0
     &      *Sqrt2*Zh13**3*la3s*xvev+4D0*Sqrt2*Zh13**3*la7s
     &      *xvev+4D0*Sqrt2*Zh13**3*la8s*xvev-Sqrt2*Zh11
     &      *Zh12*Zh13*aa5-Sqrt2*Zh11*Zh12*Zh13*aa6+Sqrt2
     &      *Zh11**2*Zh13*la1s*xvev+Sqrt2*Zh12**2*Zh13*la2s
     &      *xvev+2D0*MW/EE*SW*Zh11**3*cb*la1-2D0*MW/EE*SW
     &      *Zh11**3*la6*sb-2D0*MW/EE*SW*Zh12**3*cb*la7+2D0
     &      *MW/EE*SW*Zh12**3*la2*sb-2D0*Sqrt2*Zh11*Zh12
     &      *Zh13*la4s*xvev-2D0*Sqrt2*Zh11*Zh12*Zh13*la5s
     &      *xvev-2D0*Sqrt2*Zh11*Zh12*Zh13*la6s*xvev+2D0
     &      *MW/EE*SW*Zh11*Zh12**2*cb*la3+2D0*MW/EE*SW*Zh11
     &      *Zh12**2*cb*la4+2D0*MW/EE*SW*Zh11*Zh12**2*cb
     &      *la5-6D0*MW/EE*SW*Zh11*Zh12**2*la7*sb+2D0*MW
     &      /EE*SW*Zh11*Zh13**2*cb*la1s-2D0*MW/EE*SW*Zh11
     &      *Zh13**2*la4s*sb-2D0*MW/EE*SW*Zh11*Zh13**2*la5s
     &      *sb-2D0*MW/EE*SW*Zh11*Zh13**2*la6s*sb-6D0*MW
     &      /EE*SW*Zh11**2*Zh12*cb*la6+2D0*MW/EE*SW*Zh11
     &      **2*Zh12*la3*sb+2D0*MW/EE*SW*Zh11**2*Zh12*la4
     &      *sb+2D0*MW/EE*SW*Zh11**2*Zh12*la5*sb-2D0*MW/EE
     &      *SW*Zh12*Zh13**2*cb*la4s-2D0*MW/EE*SW*Zh12*Zh13
     &      **2*cb*la5s-2D0*MW/EE*SW*Zh12*Zh13**2*cb*la6s
     &      +2D0*MW/EE*SW*Zh12*Zh13**2*la2s*sb
      AAABR(172) = Sqrt2*Zh21**2*Zh33*aa1+Sqrt2*Zh22**2*Zh33
     &      *aa2+3D0*Sqrt2*Zh23**2*Zh33*aa3+3D0*Sqrt2*Zh23
     &      **2*Zh33*aa4-Sqrt2*Zh21*Zh22*Zh33*aa5-Sqrt2*Zh21
     &      *Zh22*Zh33*aa6+2D0*Sqrt2*Zh21*Zh23*Zh31*aa1-Sqrt2
     &      *Zh21*Zh23*Zh32*aa5-Sqrt2*Zh21*Zh23*Zh32*aa6
     &      +Sqrt2*Zh21**2*Zh33*la1s*xvev-Sqrt2*Zh22*Zh23
     &      *Zh31*aa5-Sqrt2*Zh22*Zh23*Zh31*aa6+2D0*Sqrt2
     &      *Zh22*Zh23*Zh32*aa2+Sqrt2*Zh22**2*Zh33*la2s*xvev
     &      +6D0*Sqrt2*Zh23**2*Zh33*la3s*xvev+12D0*Sqrt2
     &      *Zh23**2*Zh33*la7s*xvev+12D0*Sqrt2*Zh23**2*Zh33
     &      *la8s*xvev-2D0*Sqrt2*Zh21*Zh22*Zh33*la4s*xvev
     &      -2D0*Sqrt2*Zh21*Zh22*Zh33*la5s*xvev-2D0*Sqrt2
     &      *Zh21*Zh22*Zh33*la6s*xvev+2D0*Sqrt2*Zh21*Zh23
     &      *Zh31*la1s*xvev-2D0*Sqrt2*Zh21*Zh23*Zh32*la4s
     &      *xvev-2D0*Sqrt2*Zh21*Zh23*Zh32*la5s*xvev-2D0
     &      *Sqrt2*Zh21*Zh23*Zh32*la6s*xvev-2D0*Sqrt2*Zh22
     &      *Zh23*Zh31*la4s*xvev-2D0*Sqrt2*Zh22*Zh23*Zh31
     &      *la5s*xvev-2D0*Sqrt2*Zh22*Zh23*Zh31*la6s*xvev
     &      +2D0*Sqrt2*Zh22*Zh23*Zh32*la2s*xvev+6D0*MW/EE
     &      *SW*Zh21**2*Zh31*cb*la1-6D0*MW/EE*SW*Zh21**2
     &      *Zh31*la6*sb-6D0*MW/EE*SW*Zh21**2*Zh32*cb*la6
     &      +2D0*MW/EE*SW*Zh21**2*Zh32*la3*sb+2D0*MW/EE*SW
     &      *Zh21**2*Zh32*la4*sb+2D0*MW/EE*SW*Zh21**2*Zh32
     &      *la5*sb+2D0*MW/EE*SW*Zh22**2*Zh31*cb*la3+2D0
     &      *MW/EE*SW*Zh22**2*Zh31*cb*la4+2D0*MW/EE*SW*Zh22
     &      **2*Zh31*cb*la5-6D0*MW/EE*SW*Zh22**2*Zh31*la7
     &      *sb-6D0*MW/EE*SW*Zh22**2*Zh32*cb*la7+6D0*MW/EE
     &      *SW*Zh22**2*Zh32*la2*sb+2D0*MW/EE*SW*Zh23**2
     &      *Zh31*cb*la1s-2D0*MW/EE*SW*Zh23**2*Zh31*la4s
     &      *sb-2D0*MW/EE*SW*Zh23**2*Zh31*la5s*sb-2D0*MW
     &      /EE*SW*Zh23**2*Zh31*la6s*sb-2D0*MW/EE*SW*Zh23
     &      **2*Zh32*cb*la4s-2D0*MW/EE*SW*Zh23**2*Zh32*cb
     &      *la5s-2D0*MW/EE*SW*Zh23**2*Zh32*cb*la6s+2D0*MW
     &      /EE*SW*Zh23**2*Zh32*la2s*sb-12D0*MW/EE*SW*Zh21
     &      *Zh22*Zh31*cb*la6+4D0*MW/EE*SW*Zh21*Zh22*Zh31
     &      *la3*sb+4D0*MW/EE*SW*Zh21*Zh22*Zh31*la4*sb+4D0
     &      *MW/EE*SW*Zh21*Zh22*Zh31*la5*sb+4D0*MW/EE*SW
     &      *Zh21*Zh22*Zh32*cb*la3+4D0*MW/EE*SW*Zh21*Zh22
     &      *Zh32*cb*la4+4D0*MW/EE*SW*Zh21*Zh22*Zh32*cb*la5
     &      -12D0*MW/EE*SW*Zh21*Zh22*Zh32*la7*sb+4D0*MW/EE
     &      *SW*Zh21*Zh23*Zh33*cb*la1s-4D0*MW/EE*SW*Zh21
     &      *Zh23*Zh33*la4s*sb-4D0*MW/EE*SW*Zh21*Zh23*Zh33
     &      *la5s*sb-4D0*MW/EE*SW*Zh21*Zh23*Zh33*la6s*sb
     &      -4D0*MW/EE*SW*Zh22*Zh23*Zh33*cb*la4s-4D0*MW/EE
     &      *SW*Zh22*Zh23*Zh33*cb*la5s-4D0*MW/EE*SW*Zh22
     &      *Zh23*Zh33*cb*la6s+4D0*MW/EE*SW*Zh22*Zh23*Zh33
     &      *la2s*sb
      AAABR(173) = Sqrt2*Zh11**2*Zh33*aa1+Sqrt2*Zh12**2*Zh33
     &      *aa2+3D0*Sqrt2*Zh13**2*Zh33*aa3+3D0*Sqrt2*Zh13
     &      **2*Zh33*aa4-Sqrt2*Zh11*Zh12*Zh33*aa5-Sqrt2*Zh11
     &      *Zh12*Zh33*aa6+2D0*Sqrt2*Zh11*Zh13*Zh31*aa1-Sqrt2
     &      *Zh11*Zh13*Zh32*aa5-Sqrt2*Zh11*Zh13*Zh32*aa6
     &      +Sqrt2*Zh11**2*Zh33*la1s*xvev-Sqrt2*Zh12*Zh13
     &      *Zh31*aa5-Sqrt2*Zh12*Zh13*Zh31*aa6+2D0*Sqrt2
     &      *Zh12*Zh13*Zh32*aa2+Sqrt2*Zh12**2*Zh33*la2s*xvev
     &      +6D0*Sqrt2*Zh13**2*Zh33*la3s*xvev+12D0*Sqrt2
     &      *Zh13**2*Zh33*la7s*xvev+12D0*Sqrt2*Zh13**2*Zh33
     &      *la8s*xvev-2D0*Sqrt2*Zh11*Zh12*Zh33*la4s*xvev
     &      -2D0*Sqrt2*Zh11*Zh12*Zh33*la5s*xvev-2D0*Sqrt2
     &      *Zh11*Zh12*Zh33*la6s*xvev+2D0*Sqrt2*Zh11*Zh13
     &      *Zh31*la1s*xvev-2D0*Sqrt2*Zh11*Zh13*Zh32*la4s
     &      *xvev-2D0*Sqrt2*Zh11*Zh13*Zh32*la5s*xvev-2D0
     &      *Sqrt2*Zh11*Zh13*Zh32*la6s*xvev-2D0*Sqrt2*Zh12
     &      *Zh13*Zh31*la4s*xvev-2D0*Sqrt2*Zh12*Zh13*Zh31
     &      *la5s*xvev-2D0*Sqrt2*Zh12*Zh13*Zh31*la6s*xvev
     &      +2D0*Sqrt2*Zh12*Zh13*Zh32*la2s*xvev+6D0*MW/EE
     &      *SW*Zh11**2*Zh31*cb*la1-6D0*MW/EE*SW*Zh11**2
     &      *Zh31*la6*sb-6D0*MW/EE*SW*Zh11**2*Zh32*cb*la6
     &      +2D0*MW/EE*SW*Zh11**2*Zh32*la3*sb+2D0*MW/EE*SW
     &      *Zh11**2*Zh32*la4*sb+2D0*MW/EE*SW*Zh11**2*Zh32
     &      *la5*sb+2D0*MW/EE*SW*Zh12**2*Zh31*cb*la3+2D0
     &      *MW/EE*SW*Zh12**2*Zh31*cb*la4+2D0*MW/EE*SW*Zh12
     &      **2*Zh31*cb*la5-6D0*MW/EE*SW*Zh12**2*Zh31*la7
     &      *sb-6D0*MW/EE*SW*Zh12**2*Zh32*cb*la7+6D0*MW/EE
     &      *SW*Zh12**2*Zh32*la2*sb+2D0*MW/EE*SW*Zh13**2
     &      *Zh31*cb*la1s-2D0*MW/EE*SW*Zh13**2*Zh31*la4s
     &      *sb-2D0*MW/EE*SW*Zh13**2*Zh31*la5s*sb-2D0*MW
     &      /EE*SW*Zh13**2*Zh31*la6s*sb-2D0*MW/EE*SW*Zh13
     &      **2*Zh32*cb*la4s-2D0*MW/EE*SW*Zh13**2*Zh32*cb
     &      *la5s-2D0*MW/EE*SW*Zh13**2*Zh32*cb*la6s+2D0*MW
     &      /EE*SW*Zh13**2*Zh32*la2s*sb-12D0*MW/EE*SW*Zh11
     &      *Zh12*Zh31*cb*la6+4D0*MW/EE*SW*Zh11*Zh12*Zh31
     &      *la3*sb+4D0*MW/EE*SW*Zh11*Zh12*Zh31*la4*sb+4D0
     &      *MW/EE*SW*Zh11*Zh12*Zh31*la5*sb+4D0*MW/EE*SW
     &      *Zh11*Zh12*Zh32*cb*la3+4D0*MW/EE*SW*Zh11*Zh12
     &      *Zh32*cb*la4+4D0*MW/EE*SW*Zh11*Zh12*Zh32*cb*la5
     &      -12D0*MW/EE*SW*Zh11*Zh12*Zh32*la7*sb+4D0*MW/EE
     &      *SW*Zh11*Zh13*Zh33*cb*la1s-4D0*MW/EE*SW*Zh11
     &      *Zh13*Zh33*la4s*sb-4D0*MW/EE*SW*Zh11*Zh13*Zh33
     &      *la5s*sb-4D0*MW/EE*SW*Zh11*Zh13*Zh33*la6s*sb
     &      -4D0*MW/EE*SW*Zh12*Zh13*Zh33*cb*la4s-4D0*MW/EE
     &      *SW*Zh12*Zh13*Zh33*cb*la5s-4D0*MW/EE*SW*Zh12
     &      *Zh13*Zh33*cb*la6s+4D0*MW/EE*SW*Zh12*Zh13*Zh33
     &      *la2s*sb
      AAABR(174) = Sqrt2*Zh23**3*aa3+Sqrt2*Zh23**3*aa4+Sqrt2
     &      *Zh21**2*Zh23*aa1+Sqrt2*Zh22**2*Zh23*aa2+2D0
     &      *Sqrt2*Zh23**3*la3s*xvev+4D0*Sqrt2*Zh23**3*la7s
     &      *xvev+4D0*Sqrt2*Zh23**3*la8s*xvev-Sqrt2*Zh21
     &      *Zh22*Zh23*aa5-Sqrt2*Zh21*Zh22*Zh23*aa6+Sqrt2
     &      *Zh21**2*Zh23*la1s*xvev+Sqrt2*Zh22**2*Zh23*la2s
     &      *xvev+2D0*MW/EE*SW*Zh21**3*cb*la1-2D0*MW/EE*SW
     &      *Zh21**3*la6*sb-2D0*MW/EE*SW*Zh22**3*cb*la7+2D0
     &      *MW/EE*SW*Zh22**3*la2*sb-2D0*Sqrt2*Zh21*Zh22
     &      *Zh23*la4s*xvev-2D0*Sqrt2*Zh21*Zh22*Zh23*la5s
     &      *xvev-2D0*Sqrt2*Zh21*Zh22*Zh23*la6s*xvev+2D0
     &      *MW/EE*SW*Zh21*Zh22**2*cb*la3+2D0*MW/EE*SW*Zh21
     &      *Zh22**2*cb*la4+2D0*MW/EE*SW*Zh21*Zh22**2*cb
     &      *la5-6D0*MW/EE*SW*Zh21*Zh22**2*la7*sb+2D0*MW
     &      /EE*SW*Zh21*Zh23**2*cb*la1s-2D0*MW/EE*SW*Zh21
     &      *Zh23**2*la4s*sb-2D0*MW/EE*SW*Zh21*Zh23**2*la5s
     &      *sb-2D0*MW/EE*SW*Zh21*Zh23**2*la6s*sb-6D0*MW
     &      /EE*SW*Zh21**2*Zh22*cb*la6+2D0*MW/EE*SW*Zh21
     &      **2*Zh22*la3*sb+2D0*MW/EE*SW*Zh21**2*Zh22*la4
     &      *sb+2D0*MW/EE*SW*Zh21**2*Zh22*la5*sb-2D0*MW/EE
     &      *SW*Zh22*Zh23**2*cb*la4s-2D0*MW/EE*SW*Zh22*Zh23
     &      **2*cb*la5s-2D0*MW/EE*SW*Zh22*Zh23**2*cb*la6s
     &      +2D0*MW/EE*SW*Zh22*Zh23**2*la2s*sb
      AAABR(175) = EE*(SW/CW-CW/SW)
      AAABR(176) = EE*Pa21/SW
      AAABR(177) = EE*Pa11/SW
      AAABR(178) = EE/SW*(Zh31*sb-Zh32*cb)
      AAABR(179) = EE/SW*(Zh21*sb-Zh22*cb)
      AAABR(180) = EE/SW*(Zh11*sb-Zh12*cb)
      AAABR(181) = EE*(SW/CW*Zl11**2+2D0*SW/CW*Zl12**2-CW
     &      /SW*Zl11**2)
      AAABR(182) = EE*(SW/CW*Zl11*Zl21+2D0*SW/CW*Zl12*Zl22
     &      -CW/SW*Zl11*Zl21)
      AAABR(183) = EE/SW*Sqrt2*Zl11
      AAABR(184) = EE/SW*Sqrt2
      AAABR(185) = EE/CW*SW
      AAABR(186) = EE*(SW/CW*Zl21**2+2D0*SW/CW*Zl22**2-CW
     &      /SW*Zl21**2)
      AAABR(187) = EE/SW*Sqrt2*Zl21
      AAABR(188) = EE*(SW/CW+CW/SW)
      AAABR(189) = EE*(SW/CW*Zdd11**2-2D0*SW/CW*Zdd12**2
     &      +3D0*CW/SW*Zdd11**2)
      AAABR(190) = EE*(SW/CW*Zb11**2-2D0*SW/CW*Zb12**2+3D0
     &      *CW/SW*Zb11**2)
      AAABR(191) = EE*(SW/CW*Zss11**2-2D0*SW/CW*Zss12**2
     &      +3D0*CW/SW*Zss11**2)
      AAABR(192) = EE*(SW/CW*Zdd11*Zdd21-2D0*SW/CW*Zdd12
     &      *Zdd22+3D0*CW/SW*Zdd11*Zdd21)
      AAABR(193) = EE*(SW/CW*Zss11*Zss21-2D0*SW/CW*Zss12
     &      *Zss22+3D0*CW/SW*Zss11*Zss21)
      AAABR(194) = EE*(SW/CW*Zb11*Zb21-2D0*SW/CW*Zb12*Zb22
     &      +3D0*CW/SW*Zb11*Zb21)
      AAABR(195) = EE/SW*Sqrt2*Zb11*Zt11
      AAABR(196) = EE/SW*Sqrt2*Zdd11*Zuu11
      AAABR(197) = EE/SW*Sqrt2*Zcc11*Zss11
      AAABR(198) = EE/SW*Sqrt2*Zdd11*Zuu21
      AAABR(199) = EE/SW*Sqrt2*Zcc21*Zss11
      AAABR(200) = EE/SW*Sqrt2*Zb11*Zt21
      AAABR(201) = EE*(SW/CW*Zss21**2-2D0*SW/CW*Zss22**2
     &      +3D0*CW/SW*Zss21**2)
      AAABR(202) = EE*(SW/CW*Zb21**2-2D0*SW/CW*Zb22**2+3D0
     &      *CW/SW*Zb21**2)
      AAABR(203) = EE*(SW/CW*Zdd21**2-2D0*SW/CW*Zdd22**2
     &      +3D0*CW/SW*Zdd21**2)
      AAABR(204) = EE/SW*Sqrt2*Zcc11*Zss21
      AAABR(205) = EE/SW*Sqrt2*Zdd21*Zuu11
      AAABR(206) = EE/SW*Sqrt2*Zb21*Zt11
      AAABR(207) = EE/SW*Sqrt2*Zcc21*Zss21
      AAABR(208) = EE/SW*Sqrt2*Zb21*Zt21
      AAABR(209) = EE/SW*Sqrt2*Zdd21*Zuu21
      AAABR(210) = EE*(SW/CW*Zcc11**2+4D0*SW/CW*Zcc12**2
     &      -3D0*CW/SW*Zcc11**2)
      AAABR(211) = EE*(SW/CW*Zuu11**2+4D0*SW/CW*Zuu12**2
     &      -3D0*CW/SW*Zuu11**2)
      AAABR(212) = EE*(SW/CW*Zt11**2+4D0*SW/CW*Zt12**2-3D0
     &      *CW/SW*Zt11**2)
      AAABR(213) = EE*(SW/CW*Zt11*Zt21+4D0*SW/CW*Zt12*Zt22
     &      -3D0*CW/SW*Zt11*Zt21)
      AAABR(214) = EE*(SW/CW*Zcc11*Zcc21+4D0*SW/CW*Zcc12
     &      *Zcc22-3D0*CW/SW*Zcc11*Zcc21)
      AAABR(215) = EE*(SW/CW*Zuu11*Zuu21+4D0*SW/CW*Zuu12
     &      *Zuu22-3D0*CW/SW*Zuu11*Zuu21)
      AAABR(216) = EE*(SW/CW*Zt21**2+4D0*SW/CW*Zt22**2-3D0
     &      *CW/SW*Zt21**2)
      AAABR(217) = EE*(SW/CW*Zcc21**2+4D0*SW/CW*Zcc22**2
     &      -3D0*CW/SW*Zcc21**2)
      AAABR(218) = EE*(SW/CW*Zuu21**2+4D0*SW/CW*Zuu22**2
     &      -3D0*CW/SW*Zuu21**2)
      AAABR(219) = EE/SW
      AAABR(220) = EE/SW*(Zh11*cb+Zh12*sb)
      AAABR(221) = EE/SW*(Zh21*cb+Zh22*sb)
      AAABR(222) = EE/SW*(Zh31*cb+Zh32*sb)
      AAABR(223) = EE*(SW/CW*Zh21*cb+SW/CW*Zh22*sb+CW/SW
     &      *Zh21*cb+CW/SW*Zh22*sb)
      AAABR(224) = EE*(SW/CW*Zh31*cb+SW/CW*Zh32*sb+CW/SW
     &      *Zh31*cb+CW/SW*Zh32*sb)
      AAABR(225) = EE*(SW/CW*Zh11*cb+SW/CW*Zh12*sb+CW/SW
     &      *Zh11*cb+CW/SW*Zh12*sb)
      AAABR(226) = EE*Pa11*(SW/CW*Zh11*sb-SW/CW*Zh12*cb+CW
     &      /SW*Zh11*sb-CW/SW*Zh12*cb)
      AAABR(227) = EE*Pa11*(SW/CW*Zh31*sb-SW/CW*Zh32*cb+CW
     &      /SW*Zh31*sb-CW/SW*Zh32*cb)
      AAABR(228) = EE*Pa21*(SW/CW*Zh21*sb-SW/CW*Zh22*cb+CW
     &      /SW*Zh21*sb-CW/SW*Zh22*cb)
      AAABR(229) = EE*Pa21*(SW/CW*Zh31*sb-SW/CW*Zh32*cb+CW
     &      /SW*Zh31*sb-CW/SW*Zh32*cb)
      AAABR(230) = EE*Pa21*(SW/CW*Zh11*sb-SW/CW*Zh12*cb+CW
     &      /SW*Zh11*sb-CW/SW*Zh12*cb)
      AAABR(231) = EE*Pa11*(SW/CW*Zh21*sb-SW/CW*Zh22*cb+CW
     &      /SW*Zh21*sb-CW/SW*Zh22*cb)
      AAABR(232) = EE/CW*MW*SW
      AAABR(233) = EE*MW*(2D0*SW*Zh11*cb+2D0*SW*Zh12*sb+SW
     &      **3/CW**2*Zh11*cb+SW**3/CW**2*Zh12*sb+CW**2/SW
     &      *Zh11*cb+CW**2/SW*Zh12*sb)
      AAABR(234) = EE*MW*(2D0*SW*Zh21*cb+2D0*SW*Zh22*sb+SW
     &      **3/CW**2*Zh21*cb+SW**3/CW**2*Zh22*sb+CW**2/SW
     &      *Zh21*cb+CW**2/SW*Zh22*sb)
      AAABR(235) = EE*MW*(2D0*SW*Zh31*cb+2D0*SW*Zh32*sb+SW
     &      **3/CW**2*Zh31*cb+SW**3/CW**2*Zh32*sb+CW**2/SW
     &      *Zh31*cb+CW**2/SW*Zh32*sb)
      AAABC(1) = EE/MW/SW*Sqrt2*(NMM23*Zu22*Zv21+NMM24*Zu21
     &      *Zv22)
      AAABC(2) = EE/MW/SW*Sqrt2*(NMM23*Zu12*Zv21+NMM24*Zu11
     &      *Zv22)
      AAABC(3) = EE/MW/SW*Sqrt2*(dconjg(NMM23)*Zu22*Zv11
     &      +dconjg(NMM24)*Zu21*Zv12)
      AAABC(4) = EE/MW/SW*Sqrt2*(NMM23*Zu12*Zv11+NMM24*Zu11
     &      *Zv12)
      AAABC(5) = Sqrt2*(NMM34*Pa22*Zu12*Zv12/xvev+EE/MW*NMM23
     &      *Pa21/SW*Zu12*Zv11/cb*sb-EE/MW*NMM24*Pa21/SW
     &      *Zu11*Zv12/cb/sb+EE/MW*NMM24*Pa21/SW*Zu11*Zv12
     &      /cb*sb)
      AAABC(6) = Sqrt2*(NMM34*Pa22*Zu22*Zv12/xvev+EE/MW*NMM23
     &      *Pa21/SW*Zu22*Zv11/cb*sb-EE/MW*NMM24*Pa21/SW
     &      *Zu21*Zv12/cb/sb+EE/MW*NMM24*Pa21/SW*Zu21*Zv12
     &      /cb*sb)
      AAABC(7) = Sqrt2*(dconjg(NMM34)*Pa22*Zu12*Zv22/xvev
     &      +EE/MW*dconjg(NMM23)*Pa21/SW*Zu12*Zv21/cb*sb
     &      -EE/MW*dconjg(NMM24)*Pa21/SW*Zu11*Zv22/cb/sb
     &      +EE/MW*dconjg(NMM24)*Pa21/SW*Zu11*Zv22/cb*sb)
      AAABC(8) = Sqrt2*(NMM34*Pa12*Zu22*Zv22/xvev+EE/MW*NMM23
     &      *Pa11/SW*Zu22*Zv21/cb*sb-EE/MW*NMM24*Pa11/SW
     &      *Zu21*Zv22/cb/sb+EE/MW*NMM24*Pa11/SW*Zu21*Zv22
     &      /cb*sb)
      AAABC(9) = Sqrt2*(NMM34*Pa22*Zu22*Zv22/xvev+EE/MW*NMM23
     &      *Pa21/SW*Zu22*Zv21/cb*sb-EE/MW*NMM24*Pa21/SW
     &      *Zu21*Zv22/cb/sb+EE/MW*NMM24*Pa21/SW*Zu21*Zv22
     &      /cb*sb)
      AAABC(10) = Sqrt2*(NMM34*Pa12*Zu12*Zv22/xvev+EE/MW
     &      *NMM23*Pa11/SW*Zu12*Zv21/cb*sb-EE/MW*NMM24*Pa11
     &      /SW*Zu11*Zv22/cb/sb+EE/MW*NMM24*Pa11/SW*Zu11
     &      *Zv22/cb*sb)
      AAABC(11) = Sqrt2*(dconjg(NMM34)*Pa12*Zu22*Zv12/xvev
     &      +EE/MW*dconjg(NMM23)*Pa11/SW*Zu22*Zv11/cb*sb
     &      -EE/MW*dconjg(NMM24)*Pa11/SW*Zu21*Zv12/cb/sb
     &      +EE/MW*dconjg(NMM24)*Pa11/SW*Zu21*Zv12/cb*sb)
      AAABC(12) = Sqrt2*(NMM34*Pa12*Zu12*Zv12/xvev+EE/MW
     &      *NMM23*Pa11/SW*Zu12*Zv11/cb*sb-EE/MW*NMM24*Pa11
     &      /SW*Zu11*Zv12/cb/sb+EE/MW*NMM24*Pa11/SW*Zu11
     &      *Zv12/cb*sb)
      AAABC(13) = Sqrt2*(NMM34*Zh13*Zu22*Zv22/xvev-EE/MW
     &      *NMM23/SW*Zh11*Zu22*Zv21/cb+EE/MW*NMM24/SW*Zh12
     &      *Zu21*Zv22/sb)
      AAABC(14) = Sqrt2*(NMM34*Zh23*Zu22*Zv22/xvev-EE/MW
     &      *NMM23/SW*Zh21*Zu22*Zv21/cb+EE/MW*NMM24/SW*Zh22
     &      *Zu21*Zv22/sb)
      AAABC(15) = Sqrt2*(NMM34*Zh13*Zu12*Zv12/xvev-EE/MW
     &      *NMM23/SW*Zh11*Zu12*Zv11/cb+EE/MW*NMM24/SW*Zh12
     &      *Zu11*Zv12/sb)
      AAABC(16) = Sqrt2*(NMM34*Zh23*Zu22*Zv12/xvev-EE/MW
     &      *NMM23/SW*Zh21*Zu22*Zv11/cb+EE/MW*NMM24/SW*Zh22
     &      *Zu21*Zv12/sb)
      AAABC(17) = Sqrt2*(dconjg(NMM34)*Zh23*Zu12*Zv22/xvev
     &      -EE/MW*dconjg(NMM23)/SW*Zh21*Zu12*Zv21/cb+EE
     &      /MW*dconjg(NMM24)/SW*Zh22*Zu11*Zv22/sb)
      AAABC(18) = Sqrt2*(NMM34*Zh33*Zu22*Zv22/xvev-EE/MW
     &      *NMM23/SW*Zh31*Zu22*Zv21/cb+EE/MW*NMM24/SW*Zh32
     &      *Zu21*Zv22/sb)
      AAABC(19) = Sqrt2*(NMM34*Zh23*Zu12*Zv12/xvev-EE/MW
     &      *NMM23/SW*Zh21*Zu12*Zv11/cb+EE/MW*NMM24/SW*Zh22
     &      *Zu11*Zv12/sb)
      AAABC(20) = Sqrt2*(NMM34*Zh13*Zu22*Zv12/xvev-EE/MW
     &      *NMM23/SW*Zh11*Zu22*Zv11/cb+EE/MW*NMM24/SW*Zh12
     &      *Zu21*Zv12/sb)
      AAABC(21) = Sqrt2*(dconjg(NMM34)*Zh13*Zu12*Zv22/xvev
     &      -EE/MW*dconjg(NMM23)/SW*Zh11*Zu12*Zv21/cb+EE
     &      /MW*dconjg(NMM24)/SW*Zh12*Zu11*Zv22/sb)
      AAABC(22) = Sqrt2*(NMM34*Zh33*Zu22*Zv12/xvev-EE/MW
     &      *NMM23/SW*Zh31*Zu22*Zv11/cb+EE/MW*NMM24/SW*Zh32
     &      *Zu21*Zv12/sb)
      AAABC(23) = Sqrt2*(dconjg(NMM34)*Zh33*Zu12*Zv22/xvev
     &      -EE/MW*dconjg(NMM23)/SW*Zh31*Zu12*Zv21/cb+EE
     &      /MW*dconjg(NMM24)/SW*Zh32*Zu11*Zv22/sb)
      AAABC(24) = Sqrt2*(NMM34*Zh33*Zu12*Zv12/xvev-EE/MW
     &      *NMM23/SW*Zh31*Zu12*Zv11/cb+EE/MW*NMM24/SW*Zh32
     &      *Zu11*Zv12/sb)
      AAABR(236) = EE/SW*Zu11
      AAABR(237) = EE/SW*Zu21
      AAABR(238) = EE/SW*(2D0*Zl11*Zu11-Ml/MW*Sqrt2*Zl12
     &      *Zu12/cb)
      AAABR(239) = EE/SW*(2D0*Zl11*Zu21-Ml/MW*Sqrt2*Zl12
     &      *Zu22/cb)
      AAABR(240) = EE/SW*(2D0*Zl21*Zu11-Ml/MW*Sqrt2*Zl22
     &      *Zu12/cb)
      AAABR(241) = EE/SW*(2D0*Zl21*Zu21-Ml/MW*Sqrt2*Zl22
     &      *Zu22/cb)
      AAABC(25) = EE/MW/SW/cb*(2D0*dconjg(NMM24)*dconjg(Zn24
     &      )*Zv11/sb-2D0*dconjg(NMM24)*dconjg(Zn24)*Zv11
     &      *sb-dconjg(NMM14)*Sqrt2*dconjg(Zn21)*Zv12/sb
     &      +dconjg(NMM14)*Sqrt2*dconjg(Zn21)*Zv12*sb+dconjg(NMM24
     &      )*Sqrt2*dconjg(Zn22)*Zv12/sb-dconjg(NMM24)*Sqrt2
     &      *dconjg(Zn22)*Zv12*sb+dconjg(NMM45)*Sqrt2*dconjg(Zn25
     &      )*Zv12*sb)
      AAABC(26) = EE/MW/SW/cb*(2D0*NMM23*Zn23*Zu11*sb+NMM13
     &      *Sqrt2*Zn21*Zu12*sb-NMM23*Sqrt2*Zn22*Zu12*sb
     &      -NMM35*Sqrt2*Zn25*Zu12/sb+NMM35*Sqrt2*Zn25*Zu12
     &      *sb)
      AAABC(27) = EE/MW/SW/cb*(2D0*dconjg(NMM24)*dconjg(Zn54
     &      )*Zv11/sb-2D0*dconjg(NMM24)*dconjg(Zn54)*Zv11
     &      *sb-dconjg(NMM14)*Sqrt2*dconjg(Zn51)*Zv12/sb
     &      +dconjg(NMM14)*Sqrt2*dconjg(Zn51)*Zv12*sb+dconjg(NMM24
     &      )*Sqrt2*dconjg(Zn52)*Zv12/sb-dconjg(NMM24)*Sqrt2
     &      *dconjg(Zn52)*Zv12*sb+dconjg(NMM45)*Sqrt2*dconjg(Zn55
     &      )*Zv12*sb)
      AAABC(28) = EE/MW/SW/cb*(2D0*NMM23*Zn53*Zu11*sb+NMM13
     &      *Sqrt2*Zn51*Zu12*sb-NMM23*Sqrt2*Zn52*Zu12*sb
     &      -NMM35*Sqrt2*Zn55*Zu12/sb+NMM35*Sqrt2*Zn55*Zu12
     &      *sb)
      AAABC(29) = EE/MW/SW/cb*(2D0*dconjg(NMM24)*dconjg(Zn44
     &      )*Zv11/sb-2D0*dconjg(NMM24)*dconjg(Zn44)*Zv11
     &      *sb-dconjg(NMM14)*Sqrt2*dconjg(Zn41)*Zv12/sb
     &      +dconjg(NMM14)*Sqrt2*dconjg(Zn41)*Zv12*sb+dconjg(NMM24
     &      )*Sqrt2*dconjg(Zn42)*Zv12/sb-dconjg(NMM24)*Sqrt2
     &      *dconjg(Zn42)*Zv12*sb+dconjg(NMM45)*Sqrt2*dconjg(Zn45
     &      )*Zv12*sb)
      AAABC(30) = EE/MW/SW/cb*(2D0*NMM23*Zn43*Zu11*sb+NMM13
     &      *Sqrt2*Zn41*Zu12*sb-NMM23*Sqrt2*Zn42*Zu12*sb
     &      -NMM35*Sqrt2*Zn45*Zu12/sb+NMM35*Sqrt2*Zn45*Zu12
     &      *sb)
      AAABC(31) = EE/MW/SW/cb*(2D0*dconjg(NMM24)*dconjg(Zn34
     &      )*Zv11/sb-2D0*dconjg(NMM24)*dconjg(Zn34)*Zv11
     &      *sb-dconjg(NMM14)*Sqrt2*dconjg(Zn31)*Zv12/sb
     &      +dconjg(NMM14)*Sqrt2*dconjg(Zn31)*Zv12*sb+dconjg(NMM24
     &      )*Sqrt2*dconjg(Zn32)*Zv12/sb-dconjg(NMM24)*Sqrt2
     &      *dconjg(Zn32)*Zv12*sb+dconjg(NMM45)*Sqrt2*dconjg(Zn35
     &      )*Zv12*sb)
      AAABC(32) = EE/MW/SW/cb*(2D0*NMM23*Zn33*Zu11*sb+NMM13
     &      *Sqrt2*Zn31*Zu12*sb-NMM23*Sqrt2*Zn32*Zu12*sb
     &      -NMM35*Sqrt2*Zn35*Zu12/sb+NMM35*Sqrt2*Zn35*Zu12
     &      *sb)
      AAABC(33) = EE/MW/SW/cb*(2D0*dconjg(NMM24)*dconjg(Zn14
     &      )*Zv11/sb-2D0*dconjg(NMM24)*dconjg(Zn14)*Zv11
     &      *sb-dconjg(NMM14)*Sqrt2*dconjg(Zn11)*Zv12/sb
     &      +dconjg(NMM14)*Sqrt2*dconjg(Zn11)*Zv12*sb+dconjg(NMM24
     &      )*Sqrt2*dconjg(Zn12)*Zv12/sb-dconjg(NMM24)*Sqrt2
     &      *dconjg(Zn12)*Zv12*sb+dconjg(NMM45)*Sqrt2*dconjg(Zn15
     &      )*Zv12*sb)
      AAABC(34) = EE/MW/SW/cb*(2D0*NMM23*Zn13*Zu11*sb+NMM13
     &      *Sqrt2*Zn11*Zu12*sb-NMM23*Sqrt2*Zn12*Zu12*sb
     &      -NMM35*Sqrt2*Zn15*Zu12/sb+NMM35*Sqrt2*Zn15*Zu12
     &      *sb)
      AAABC(35) = EE/MW/SW/cb*(2D0*dconjg(NMM24)*dconjg(Zn34
     &      )*Zv21/sb-2D0*dconjg(NMM24)*dconjg(Zn34)*Zv21
     &      *sb-dconjg(NMM14)*Sqrt2*dconjg(Zn31)*Zv22/sb
     &      +dconjg(NMM14)*Sqrt2*dconjg(Zn31)*Zv22*sb+dconjg(NMM24
     &      )*Sqrt2*dconjg(Zn32)*Zv22/sb-dconjg(NMM24)*Sqrt2
     &      *dconjg(Zn32)*Zv22*sb+dconjg(NMM45)*Sqrt2*dconjg(Zn35
     &      )*Zv22*sb)
      AAABC(36) = EE/MW/SW/cb*(2D0*NMM23*Zn33*Zu21*sb+NMM13
     &      *Sqrt2*Zn31*Zu22*sb-NMM23*Sqrt2*Zn32*Zu22*sb
     &      -NMM35*Sqrt2*Zn35*Zu22/sb+NMM35*Sqrt2*Zn35*Zu22
     &      *sb)
      AAABC(37) = EE/MW/SW/cb*(2D0*dconjg(NMM24)*dconjg(Zn24
     &      )*Zv21/sb-2D0*dconjg(NMM24)*dconjg(Zn24)*Zv21
     &      *sb-dconjg(NMM14)*Sqrt2*dconjg(Zn21)*Zv22/sb
     &      +dconjg(NMM14)*Sqrt2*dconjg(Zn21)*Zv22*sb+dconjg(NMM24
     &      )*Sqrt2*dconjg(Zn22)*Zv22/sb-dconjg(NMM24)*Sqrt2
     &      *dconjg(Zn22)*Zv22*sb+dconjg(NMM45)*Sqrt2*dconjg(Zn25
     &      )*Zv22*sb)
      AAABC(38) = EE/MW/SW/cb*(2D0*NMM23*Zn23*Zu21*sb+NMM13
     &      *Sqrt2*Zn21*Zu22*sb-NMM23*Sqrt2*Zn22*Zu22*sb
     &      -NMM35*Sqrt2*Zn25*Zu22/sb+NMM35*Sqrt2*Zn25*Zu22
     &      *sb)
      AAABC(39) = EE/MW/SW/cb*(2D0*dconjg(NMM24)*dconjg(Zn14
     &      )*Zv21/sb-2D0*dconjg(NMM24)*dconjg(Zn14)*Zv21
     &      *sb-dconjg(NMM14)*Sqrt2*dconjg(Zn11)*Zv22/sb
     &      +dconjg(NMM14)*Sqrt2*dconjg(Zn11)*Zv22*sb+dconjg(NMM24
     &      )*Sqrt2*dconjg(Zn12)*Zv22/sb-dconjg(NMM24)*Sqrt2
     &      *dconjg(Zn12)*Zv22*sb+dconjg(NMM45)*Sqrt2*dconjg(Zn15
     &      )*Zv22*sb)
      AAABC(40) = EE/MW/SW/cb*(2D0*NMM23*Zn13*Zu21*sb+NMM13
     &      *Sqrt2*Zn11*Zu22*sb-NMM23*Sqrt2*Zn12*Zu22*sb
     &      -NMM35*Sqrt2*Zn15*Zu22/sb+NMM35*Sqrt2*Zn15*Zu22
     &      *sb)
      AAABC(41) = EE/MW/SW/cb*(2D0*dconjg(NMM24)*dconjg(Zn44
     &      )*Zv21/sb-2D0*dconjg(NMM24)*dconjg(Zn44)*Zv21
     &      *sb-dconjg(NMM14)*Sqrt2*dconjg(Zn41)*Zv22/sb
     &      +dconjg(NMM14)*Sqrt2*dconjg(Zn41)*Zv22*sb+dconjg(NMM24
     &      )*Sqrt2*dconjg(Zn42)*Zv22/sb-dconjg(NMM24)*Sqrt2
     &      *dconjg(Zn42)*Zv22*sb+dconjg(NMM45)*Sqrt2*dconjg(Zn45
     &      )*Zv22*sb)
      AAABC(42) = EE/MW/SW/cb*(2D0*NMM23*Zn43*Zu21*sb+NMM13
     &      *Sqrt2*Zn41*Zu22*sb-NMM23*Sqrt2*Zn42*Zu22*sb
     &      -NMM35*Sqrt2*Zn45*Zu22/sb+NMM35*Sqrt2*Zn45*Zu22
     &      *sb)
      AAABC(43) = EE/MW/SW/cb*(2D0*dconjg(NMM24)*dconjg(Zn54
     &      )*Zv21/sb-2D0*dconjg(NMM24)*dconjg(Zn54)*Zv21
     &      *sb-dconjg(NMM14)*Sqrt2*dconjg(Zn51)*Zv22/sb
     &      +dconjg(NMM14)*Sqrt2*dconjg(Zn51)*Zv22*sb+dconjg(NMM24
     &      )*Sqrt2*dconjg(Zn52)*Zv22/sb-dconjg(NMM24)*Sqrt2
     &      *dconjg(Zn52)*Zv22*sb+dconjg(NMM45)*Sqrt2*dconjg(Zn55
     &      )*Zv22*sb)
      AAABC(44) = EE/MW/SW/cb*(2D0*NMM23*Zn53*Zu21*sb+NMM13
     &      *Sqrt2*Zn51*Zu22*sb-NMM23*Sqrt2*Zn52*Zu22*sb
     &      -NMM35*Sqrt2*Zn55*Zu22/sb+NMM35*Sqrt2*Zn55*Zu22
     &      *sb)
      AAABC(45) = EE/MW/SW*(2D0*dconjg(NMM24)*dconjg(Zn54
     &      )*Zv21-dconjg(NMM14)*Sqrt2*dconjg(Zn51)*Zv22
     &      +dconjg(NMM24)*Sqrt2*dconjg(Zn52)*Zv22-dconjg(NMM45
     &      )*Sqrt2*dconjg(Zn55)*Zv22)
      AAABC(46) = EE/MW/SW*(2D0*NMM23*Zn53*Zu21+NMM13*Sqrt2
     &      *Zn51*Zu22-NMM23*Sqrt2*Zn52*Zu22+NMM35*Sqrt2
     &      *Zn55*Zu22)
      AAABC(47) = EE/MW/SW*(2D0*dconjg(NMM24)*dconjg(Zn34
     &      )*Zv21-dconjg(NMM14)*Sqrt2*dconjg(Zn31)*Zv22
     &      +dconjg(NMM24)*Sqrt2*dconjg(Zn32)*Zv22-dconjg(NMM45
     &      )*Sqrt2*dconjg(Zn35)*Zv22)
      AAABC(48) = EE/MW/SW*(2D0*NMM23*Zn33*Zu21+NMM13*Sqrt2
     &      *Zn31*Zu22-NMM23*Sqrt2*Zn32*Zu22+NMM35*Sqrt2
     &      *Zn35*Zu22)
      AAABC(49) = EE/MW/SW*(2D0*dconjg(NMM24)*dconjg(Zn24
     &      )*Zv21-dconjg(NMM14)*Sqrt2*dconjg(Zn21)*Zv22
     &      +dconjg(NMM24)*Sqrt2*dconjg(Zn22)*Zv22-dconjg(NMM45
     &      )*Sqrt2*dconjg(Zn25)*Zv22)
      AAABC(50) = EE/MW/SW*(2D0*NMM23*Zn23*Zu21+NMM13*Sqrt2
     &      *Zn21*Zu22-NMM23*Sqrt2*Zn22*Zu22+NMM35*Sqrt2
     &      *Zn25*Zu22)
      AAABC(51) = EE/MW/SW*(2D0*dconjg(NMM24)*dconjg(Zn44
     &      )*Zv11-dconjg(NMM14)*Sqrt2*dconjg(Zn41)*Zv12
     &      +dconjg(NMM24)*Sqrt2*dconjg(Zn42)*Zv12-dconjg(NMM45
     &      )*Sqrt2*dconjg(Zn45)*Zv12)
      AAABC(52) = EE/MW/SW*(2D0*NMM23*Zn43*Zu11+NMM13*Sqrt2
     &      *Zn41*Zu12-NMM23*Sqrt2*Zn42*Zu12+NMM35*Sqrt2
     &      *Zn45*Zu12)
      AAABC(53) = EE/MW/SW*(2D0*dconjg(NMM24)*dconjg(Zn34
     &      )*Zv11-dconjg(NMM14)*Sqrt2*dconjg(Zn31)*Zv12
     &      +dconjg(NMM24)*Sqrt2*dconjg(Zn32)*Zv12-dconjg(NMM45
     &      )*Sqrt2*dconjg(Zn35)*Zv12)
      AAABC(54) = EE/MW/SW*(2D0*NMM23*Zn33*Zu11+NMM13*Sqrt2
     &      *Zn31*Zu12-NMM23*Sqrt2*Zn32*Zu12+NMM35*Sqrt2
     &      *Zn35*Zu12)
      AAABC(55) = EE/MW/SW*(2D0*dconjg(NMM24)*dconjg(Zn14
     &      )*Zv21-dconjg(NMM14)*Sqrt2*dconjg(Zn11)*Zv22
     &      +dconjg(NMM24)*Sqrt2*dconjg(Zn12)*Zv22-dconjg(NMM45
     &      )*Sqrt2*dconjg(Zn15)*Zv22)
      AAABC(56) = EE/MW/SW*(2D0*NMM23*Zn13*Zu21+NMM13*Sqrt2
     &      *Zn11*Zu22-NMM23*Sqrt2*Zn12*Zu22+NMM35*Sqrt2
     &      *Zn15*Zu22)
      AAABC(57) = EE/MW/SW*(2D0*dconjg(NMM24)*dconjg(Zn14
     &      )*Zv11-dconjg(NMM14)*Sqrt2*dconjg(Zn11)*Zv12
     &      +dconjg(NMM24)*Sqrt2*dconjg(Zn12)*Zv12-dconjg(NMM45
     &      )*Sqrt2*dconjg(Zn15)*Zv12)
      AAABC(58) = EE/MW/SW*(2D0*NMM23*Zn13*Zu11+NMM13*Sqrt2
     &      *Zn11*Zu12-NMM23*Sqrt2*Zn12*Zu12+NMM35*Sqrt2
     &      *Zn15*Zu12)
      AAABC(59) = EE/MW/SW*(2D0*dconjg(NMM24)*dconjg(Zn54
     &      )*Zv11-dconjg(NMM14)*Sqrt2*dconjg(Zn51)*Zv12
     &      +dconjg(NMM24)*Sqrt2*dconjg(Zn52)*Zv12-dconjg(NMM45
     &      )*Sqrt2*dconjg(Zn55)*Zv12)
      AAABC(60) = EE/MW/SW*(2D0*NMM23*Zn53*Zu11+NMM13*Sqrt2
     &      *Zn51*Zu12-NMM23*Sqrt2*Zn52*Zu12+NMM35*Sqrt2
     &      *Zn55*Zu12)
      AAABC(61) = EE/MW/SW*(2D0*dconjg(NMM24)*dconjg(Zn24
     &      )*Zv11-dconjg(NMM14)*Sqrt2*dconjg(Zn21)*Zv12
     &      +dconjg(NMM24)*Sqrt2*dconjg(Zn22)*Zv12-dconjg(NMM45
     &      )*Sqrt2*dconjg(Zn25)*Zv12)
      AAABC(62) = EE/MW/SW*(2D0*NMM23*Zn23*Zu11+NMM13*Sqrt2
     &      *Zn21*Zu12-NMM23*Sqrt2*Zn22*Zu12+NMM35*Sqrt2
     &      *Zn25*Zu12)
      AAABC(63) = EE/MW/SW*(2D0*dconjg(NMM24)*dconjg(Zn44
     &      )*Zv21-dconjg(NMM14)*Sqrt2*dconjg(Zn41)*Zv22
     &      +dconjg(NMM24)*Sqrt2*dconjg(Zn42)*Zv22-dconjg(NMM45
     &      )*Sqrt2*dconjg(Zn45)*Zv22)
      AAABC(64) = EE/MW/SW*(2D0*NMM23*Zn43*Zu21+NMM13*Sqrt2
     &      *Zn41*Zu22-NMM23*Sqrt2*Zn42*Zu22+NMM35*Sqrt2
     &      *Zn45*Zu22)
      AAABR(242) = EE/SW*(2D0*Zb11*Zu11-MbMM/MW*Sqrt2*Zb12
     &      *Zu12/cb)
      AAABR(243) = EE/MW*MtMM/SW*Sqrt2*Zb11*Zv12/sb
      AAABR(244) = EE/SW*(2D0*Zss11*Zu11-Mq/MW*Sqrt2*Zss12
     &      *Zu12/cb)
      AAABR(245) = EE/MW*Mc/SW*Sqrt2*Zss11*Zv12/sb
      AAABR(246) = EE/SW*(2D0*Zss11*Zu21-Mq/MW*Sqrt2*Zss12
     &      *Zu22/cb)
      AAABR(247) = EE/MW*Mc/SW*Sqrt2*Zss11*Zv22/sb
      AAABR(248) = EE/SW*(2D0*Zb11*Zu21-MbMM/MW*Sqrt2*Zb12
     &      *Zu22/cb)
      AAABR(249) = EE/MW*MtMM/SW*Sqrt2*Zb11*Zv22/sb
      AAABR(250) = EE/SW*(2D0*Zdd11*Zu21-Mq/MW*Sqrt2*Zdd12
     &      *Zu22/cb)
      AAABR(251) = EE/MW*Mq/SW*Sqrt2*Zdd11*Zv22/sb
      AAABR(252) = EE/SW*(2D0*Zdd11*Zu11-Mq/MW*Sqrt2*Zdd12
     &      *Zu12/cb)
      AAABR(253) = EE/MW*Mq/SW*Sqrt2*Zdd11*Zv12/sb
      AAABR(254) = EE/SW*(2D0*Zb21*Zu21-MbMM/MW*Sqrt2*Zb22
     &      *Zu22/cb)
      AAABR(255) = EE/MW*MtMM/SW*Sqrt2*Zb21*Zv22/sb
      AAABR(256) = EE/SW*(2D0*Zss21*Zu11-Mq/MW*Sqrt2*Zss22
     &      *Zu12/cb)
      AAABR(257) = EE/MW*Mc/SW*Sqrt2*Zss21*Zv12/sb
      AAABR(258) = EE/SW*(2D0*Zb21*Zu11-MbMM/MW*Sqrt2*Zb22
     &      *Zu12/cb)
      AAABR(259) = EE/MW*MtMM/SW*Sqrt2*Zb21*Zv12/sb
      AAABR(260) = EE/SW*(2D0*Zdd21*Zu21-Mq/MW*Sqrt2*Zdd22
     &      *Zu22/cb)
      AAABR(261) = EE/MW*Mq/SW*Sqrt2*Zdd21*Zv22/sb
      AAABR(262) = EE/SW*(2D0*Zdd21*Zu11-Mq/MW*Sqrt2*Zdd22
     &      *Zu12/cb)
      AAABR(263) = EE/MW*Mq/SW*Sqrt2*Zdd21*Zv12/sb
      AAABR(264) = EE/SW*(2D0*Zss21*Zu21-Mq/MW*Sqrt2*Zss22
     &      *Zu22/cb)
      AAABR(265) = EE/MW*Mc/SW*Sqrt2*Zss21*Zv22/sb
      AAABR(266) = EE/MW*Ml/SW
      AAABR(267) = EE/MW*Ml*Pa11/SW*tB
      AAABR(268) = EE/MW*Ml*Pa21/SW*tB
      AAABR(269) = EE/MW*Ml/SW*Zh31/cb
      AAABR(270) = EE/MW*Ml/SW*Zh11/cb
      AAABR(271) = EE/MW*Ml/SW*Zh21/cb
      AAABR(272) = EE/MW*Ml/SW*Sqrt2*tB
      AAABR(273) = EE/MW*Ml/SW*Sqrt2
      AAABC(65) = EE*Sqrt2*(dconjg(Zn21)/CW+dconjg(Zn22)/SW)
      AAABC(66) = EE*Sqrt2*(dconjg(Zn11)/CW+dconjg(Zn12)/SW)
      AAABC(67) = EE*Sqrt2*(dconjg(Zn41)/CW+dconjg(Zn42)/SW)
      AAABC(68) = EE*Sqrt2*(Zl11/CW*dconjg(Zn41)+Zl11/SW
     &      *dconjg(Zn42)-Ml/MW/SW*Zl12*dconjg(Zn43)/cb)
      AAABC(69) = EE*Sqrt2*(2D0*Zl12/CW*Zn41+Ml/MW/SW*Zl11
     &      *Zn43/cb)
      AAABC(70) = EE*Sqrt2*(Zl11/CW*dconjg(Zn11)+Zl11/SW
     &      *dconjg(Zn12)-Ml/MW/SW*Zl12*dconjg(Zn13)/cb)
      AAABC(71) = EE*Sqrt2*(2D0*Zl12/CW*Zn11+Ml/MW/SW*Zl11
     &      *Zn13/cb)
      AAABC(72) = EE*Sqrt2*(dconjg(Zn51)/CW+dconjg(Zn52)/SW)
      AAABC(73) = EE*Sqrt2*(dconjg(Zn31)/CW+dconjg(Zn32)/SW)
      AAABC(74) = EE*Sqrt2*(Zl11/CW*dconjg(Zn21)+Zl11/SW
     &      *dconjg(Zn22)-Ml/MW/SW*Zl12*dconjg(Zn23)/cb)
      AAABC(75) = EE*Sqrt2*(2D0*Zl12/CW*Zn21+Ml/MW/SW*Zl11
     &      *Zn23/cb)
      AAABC(76) = EE*Sqrt2*(Zl11/CW*dconjg(Zn31)+Zl11/SW
     &      *dconjg(Zn32)-Ml/MW/SW*Zl12*dconjg(Zn33)/cb)
      AAABC(77) = EE*Sqrt2*(2D0*Zl12/CW*Zn31+Ml/MW/SW*Zl11
     &      *Zn33/cb)
      AAABC(78) = EE*Sqrt2*(Zl11/CW*dconjg(Zn51)+Zl11/SW
     &      *dconjg(Zn52)-Ml/MW/SW*Zl12*dconjg(Zn53)/cb)
      AAABC(79) = EE*Sqrt2*(2D0*Zl12/CW*Zn51+Ml/MW/SW*Zl11
     &      *Zn53/cb)
      AAABC(80) = EE/CW*Sqrt2*Zn11
      AAABC(81) = EE*Sqrt2*(Zl21/CW*dconjg(Zn41)+Zl21/SW
     &      *dconjg(Zn42)-Ml/MW/SW*Zl22*dconjg(Zn43)/cb)
      AAABC(82) = EE*Sqrt2*(2D0*Zl22/CW*Zn41+Ml/MW/SW*Zl21
     &      *Zn43/cb)
      AAABC(83) = EE/CW*Sqrt2*Zn31
      AAABC(84) = EE/CW*Sqrt2*Zn41
      AAABC(85) = EE*Sqrt2*(Zl21/CW*dconjg(Zn31)+Zl21/SW
     &      *dconjg(Zn32)-Ml/MW/SW*Zl22*dconjg(Zn33)/cb)
      AAABC(86) = EE*Sqrt2*(2D0*Zl22/CW*Zn31+Ml/MW/SW*Zl21
     &      *Zn33/cb)
      AAABC(87) = EE/CW*Sqrt2*Zn51
      AAABC(88) = EE*Sqrt2*(Zl21/CW*dconjg(Zn21)+Zl21/SW
     &      *dconjg(Zn22)-Ml/MW/SW*Zl22*dconjg(Zn23)/cb)
      AAABC(89) = EE*Sqrt2*(2D0*Zl22/CW*Zn21+Ml/MW/SW*Zl21
     &      *Zn23/cb)
      AAABC(90) = EE/CW*Sqrt2*Zn21
      AAABC(91) = EE*Sqrt2*(Zl21/CW*dconjg(Zn51)+Zl21/SW
     &      *dconjg(Zn52)-Ml/MW/SW*Zl22*dconjg(Zn53)/cb)
      AAABC(92) = EE*Sqrt2*(2D0*Zl22/CW*Zn51+Ml/MW/SW*Zl21
     &      *Zn53/cb)
      AAABC(93) = EE*Sqrt2*(Zl21/CW*dconjg(Zn11)+Zl21/SW
     &      *dconjg(Zn12)-Ml/MW/SW*Zl22*dconjg(Zn13)/cb)
      AAABC(94) = EE*Sqrt2*(2D0*Zl22/CW*Zn11+Ml/MW/SW*Zl21
     &      *Zn13/cb)
      AAABR(274) = EE/SW*Zv11
      AAABR(275) = EE/MW*Ml/SW*Sqrt2*Zu12/cb
      AAABR(276) = EE/SW*Zv21
      AAABR(277) = EE/MW*Ml/SW*Sqrt2*Zu22/cb
      AAABC(95) = EE*Sqrt2*(dconjg(Zn51)/CW-dconjg(Zn52)/SW)
      AAABC(96) = EE*Sqrt2*(dconjg(Zn41)/CW-dconjg(Zn42)/SW)
      AAABC(97) = EE*Sqrt2*(dconjg(Zn31)/CW-dconjg(Zn32)/SW)
      AAABC(98) = EE*Sqrt2*(dconjg(Zn21)/CW-dconjg(Zn22)/SW)
      AAABC(99) = EE*Sqrt2*(dconjg(Zn11)/CW-dconjg(Zn12)/SW)
      AAABC(100) = EE*Sqrt2*(dconjg(Zn21)/CW*Zss11-3D0*dconjg(Zn22
     &      )/SW*Zss11+3D0*Mq/MW/SW*dconjg(Zn23)*Zss12/cb)
      AAABC(101) = EE*Sqrt2*(2D0*Zn21/CW*Zss12+3D0*Mq/MW
     &      /SW*Zn23*Zss11/cb)
      AAABC(102) = EE*Sqrt2*(Zdd11/CW*dconjg(Zn41)-3D0*Zdd11
     &      /SW*dconjg(Zn42)+3D0*Mq/MW/SW*Zdd12*dconjg(Zn43
     &      )/cb)
      AAABC(103) = EE*Sqrt2*(2D0*Zdd12/CW*Zn41+3D0*Mq/MW
     &      /SW*Zdd11*Zn43/cb)
      AAABC(104) = EE*Sqrt2*(Zb11/CW*dconjg(Zn11)-3D0*Zb11
     &      /SW*dconjg(Zn12)+3D0*MbMM/MW/SW*Zb12*dconjg(Zn13
     &      )/cb)
      AAABC(105) = EE*Sqrt2*(2D0*Zb12/CW*Zn11+3D0*MbMM/MW
     &      /SW*Zb11*Zn13/cb)
      AAABC(106) = EE*Sqrt2*(Zb11/CW*dconjg(Zn41)-3D0*Zb11
     &      /SW*dconjg(Zn42)+3D0*MbMM/MW/SW*Zb12*dconjg(Zn43
     &      )/cb)
      AAABC(107) = EE*Sqrt2*(2D0*Zb12/CW*Zn41+3D0*MbMM/MW
     &      /SW*Zb11*Zn43/cb)
      AAABC(108) = EE*Sqrt2*(dconjg(Zn11)/CW*Zss11-3D0*dconjg(Zn12
     &      )/SW*Zss11+3D0*Mq/MW/SW*dconjg(Zn13)*Zss12/cb)
      AAABC(109) = EE*Sqrt2*(2D0*Zn11/CW*Zss12+3D0*Mq/MW
     &      /SW*Zn13*Zss11/cb)
      AAABC(110) = EE*Sqrt2*(Zdd11/CW*dconjg(Zn21)-3D0*Zdd11
     &      /SW*dconjg(Zn22)+3D0*Mq/MW/SW*Zdd12*dconjg(Zn23
     &      )/cb)
      AAABC(111) = EE*Sqrt2*(2D0*Zdd12/CW*Zn21+3D0*Mq/MW
     &      /SW*Zdd11*Zn23/cb)
      AAABC(112) = EE*Sqrt2*(Zb11/CW*dconjg(Zn31)-3D0*Zb11
     &      /SW*dconjg(Zn32)+3D0*MbMM/MW/SW*Zb12*dconjg(Zn33
     &      )/cb)
      AAABC(113) = EE*Sqrt2*(2D0*Zb12/CW*Zn31+3D0*MbMM/MW
     &      /SW*Zb11*Zn33/cb)
      AAABC(114) = EE*Sqrt2*(Zdd11/CW*dconjg(Zn51)-3D0*Zdd11
     &      /SW*dconjg(Zn52)+3D0*Mq/MW/SW*Zdd12*dconjg(Zn53
     &      )/cb)
      AAABC(115) = EE*Sqrt2*(2D0*Zdd12/CW*Zn51+3D0*Mq/MW
     &      /SW*Zdd11*Zn53/cb)
      AAABC(116) = EE*Sqrt2*(dconjg(Zn41)/CW*Zss11-3D0*dconjg(Zn42
     &      )/SW*Zss11+3D0*Mq/MW/SW*dconjg(Zn43)*Zss12/cb)
      AAABC(117) = EE*Sqrt2*(2D0*Zn41/CW*Zss12+3D0*Mq/MW
     &      /SW*Zn43*Zss11/cb)
      AAABC(118) = EE*Sqrt2*(dconjg(Zn31)/CW*Zss11-3D0*dconjg(Zn32
     &      )/SW*Zss11+3D0*Mq/MW/SW*dconjg(Zn33)*Zss12/cb)
      AAABC(119) = EE*Sqrt2*(2D0*Zn31/CW*Zss12+3D0*Mq/MW
     &      /SW*Zn33*Zss11/cb)
      AAABC(120) = EE*Sqrt2*(Zdd11/CW*dconjg(Zn11)-3D0*Zdd11
     &      /SW*dconjg(Zn12)+3D0*Mq/MW/SW*Zdd12*dconjg(Zn13
     &      )/cb)
      AAABC(121) = EE*Sqrt2*(2D0*Zdd12/CW*Zn11+3D0*Mq/MW
     &      /SW*Zdd11*Zn13/cb)
      AAABC(122) = EE*Sqrt2*(Zdd11/CW*dconjg(Zn31)-3D0*Zdd11
     &      /SW*dconjg(Zn32)+3D0*Mq/MW/SW*Zdd12*dconjg(Zn33
     &      )/cb)
      AAABC(123) = EE*Sqrt2*(2D0*Zdd12/CW*Zn31+3D0*Mq/MW
     &      /SW*Zdd11*Zn33/cb)
      AAABC(124) = EE*Sqrt2*(Zb11/CW*dconjg(Zn21)-3D0*Zb11
     &      /SW*dconjg(Zn22)+3D0*MbMM/MW/SW*Zb12*dconjg(Zn23
     &      )/cb)
      AAABC(125) = EE*Sqrt2*(2D0*Zb12/CW*Zn21+3D0*MbMM/MW
     &      /SW*Zb11*Zn23/cb)
      AAABC(126) = EE*Sqrt2*(Zb11/CW*dconjg(Zn51)-3D0*Zb11
     &      /SW*dconjg(Zn52)+3D0*MbMM/MW/SW*Zb12*dconjg(Zn53
     &      )/cb)
      AAABC(127) = EE*Sqrt2*(2D0*Zb12/CW*Zn51+3D0*MbMM/MW
     &      /SW*Zb11*Zn53/cb)
      AAABC(128) = EE*Sqrt2*(dconjg(Zn51)/CW*Zss11-3D0*dconjg(Zn52
     &      )/SW*Zss11+3D0*Mq/MW/SW*dconjg(Zn53)*Zss12/cb)
      AAABC(129) = EE*Sqrt2*(2D0*Zn51/CW*Zss12+3D0*Mq/MW
     &      /SW*Zn53*Zss11/cb)
      AAABC(130) = EE*Sqrt2*(Zb21/CW*dconjg(Zn51)-3D0*Zb21
     &      /SW*dconjg(Zn52)+3D0*MbMM/MW/SW*Zb22*dconjg(Zn53
     &      )/cb)
      AAABC(131) = EE*Sqrt2*(2D0*Zb22/CW*Zn51+3D0*MbMM/MW
     &      /SW*Zb21*Zn53/cb)
      AAABC(132) = EE*Sqrt2*(dconjg(Zn41)/CW*Zss21-3D0*dconjg(Zn42
     &      )/SW*Zss21+3D0*Mq/MW/SW*dconjg(Zn43)*Zss22/cb)
      AAABC(133) = EE*Sqrt2*(2D0*Zn41/CW*Zss22+3D0*Mq/MW
     &      /SW*Zn43*Zss21/cb)
      AAABC(134) = EE*Sqrt2*(Zb21/CW*dconjg(Zn41)-3D0*Zb21
     &      /SW*dconjg(Zn42)+3D0*MbMM/MW/SW*Zb22*dconjg(Zn43
     &      )/cb)
      AAABC(135) = EE*Sqrt2*(2D0*Zb22/CW*Zn41+3D0*MbMM/MW
     &      /SW*Zb21*Zn43/cb)
      AAABC(136) = EE*Sqrt2*(Zdd21/CW*dconjg(Zn31)-3D0*Zdd21
     &      /SW*dconjg(Zn32)+3D0*Mq/MW/SW*Zdd22*dconjg(Zn33
     &      )/cb)
      AAABC(137) = EE*Sqrt2*(2D0*Zdd22/CW*Zn31+3D0*Mq/MW
     &      /SW*Zdd21*Zn33/cb)
      AAABC(138) = EE*Sqrt2*(Zdd21/CW*dconjg(Zn21)-3D0*Zdd21
     &      /SW*dconjg(Zn22)+3D0*Mq/MW/SW*Zdd22*dconjg(Zn23
     &      )/cb)
      AAABC(139) = EE*Sqrt2*(2D0*Zdd22/CW*Zn21+3D0*Mq/MW
     &      /SW*Zdd21*Zn23/cb)
      AAABC(140) = EE*Sqrt2*(dconjg(Zn51)/CW*Zss21-3D0*dconjg(Zn52
     &      )/SW*Zss21+3D0*Mq/MW/SW*dconjg(Zn53)*Zss22/cb)
      AAABC(141) = EE*Sqrt2*(2D0*Zn51/CW*Zss22+3D0*Mq/MW
     &      /SW*Zn53*Zss21/cb)
      AAABC(142) = EE*Sqrt2*(Zb21/CW*dconjg(Zn11)-3D0*Zb21
     &      /SW*dconjg(Zn12)+3D0*MbMM/MW/SW*Zb22*dconjg(Zn13
     &      )/cb)
      AAABC(143) = EE*Sqrt2*(2D0*Zb22/CW*Zn11+3D0*MbMM/MW
     &      /SW*Zb21*Zn13/cb)
      AAABC(144) = EE*Sqrt2*(Zb21/CW*dconjg(Zn21)-3D0*Zb21
     &      /SW*dconjg(Zn22)+3D0*MbMM/MW/SW*Zb22*dconjg(Zn23
     &      )/cb)
      AAABC(145) = EE*Sqrt2*(2D0*Zb22/CW*Zn21+3D0*MbMM/MW
     &      /SW*Zb21*Zn23/cb)
      AAABC(146) = EE*Sqrt2*(dconjg(Zn31)/CW*Zss21-3D0*dconjg(Zn32
     &      )/SW*Zss21+3D0*Mq/MW/SW*dconjg(Zn33)*Zss22/cb)
      AAABC(147) = EE*Sqrt2*(2D0*Zn31/CW*Zss22+3D0*Mq/MW
     &      /SW*Zn33*Zss21/cb)
      AAABC(148) = EE*Sqrt2*(Zb21/CW*dconjg(Zn31)-3D0*Zb21
     &      /SW*dconjg(Zn32)+3D0*MbMM/MW/SW*Zb22*dconjg(Zn33
     &      )/cb)
      AAABC(149) = EE*Sqrt2*(2D0*Zb22/CW*Zn31+3D0*MbMM/MW
     &      /SW*Zb21*Zn33/cb)
      AAABC(150) = EE*Sqrt2*(dconjg(Zn21)/CW*Zss21-3D0*dconjg(Zn22
     &      )/SW*Zss21+3D0*Mq/MW/SW*dconjg(Zn23)*Zss22/cb)
      AAABC(151) = EE*Sqrt2*(2D0*Zn21/CW*Zss22+3D0*Mq/MW
     &      /SW*Zn23*Zss21/cb)
      AAABC(152) = EE*Sqrt2*(Zdd21/CW*dconjg(Zn11)-3D0*Zdd21
     &      /SW*dconjg(Zn12)+3D0*Mq/MW/SW*Zdd22*dconjg(Zn13
     &      )/cb)
      AAABC(153) = EE*Sqrt2*(2D0*Zdd22/CW*Zn11+3D0*Mq/MW
     &      /SW*Zdd21*Zn13/cb)
      AAABC(154) = EE*Sqrt2*(dconjg(Zn11)/CW*Zss21-3D0*dconjg(Zn12
     &      )/SW*Zss21+3D0*Mq/MW/SW*dconjg(Zn13)*Zss22/cb)
      AAABC(155) = EE*Sqrt2*(2D0*Zn11/CW*Zss22+3D0*Mq/MW
     &      /SW*Zn13*Zss21/cb)
      AAABC(156) = EE*Sqrt2*(Zdd21/CW*dconjg(Zn41)-3D0*Zdd21
     &      /SW*dconjg(Zn42)+3D0*Mq/MW/SW*Zdd22*dconjg(Zn43
     &      )/cb)
      AAABC(157) = EE*Sqrt2*(2D0*Zdd22/CW*Zn41+3D0*Mq/MW
     &      /SW*Zdd21*Zn43/cb)
      AAABC(158) = EE*Sqrt2*(Zdd21/CW*dconjg(Zn51)-3D0*Zdd21
     &      /SW*dconjg(Zn52)+3D0*Mq/MW/SW*Zdd22*dconjg(Zn53
     &      )/cb)
      AAABC(159) = EE*Sqrt2*(2D0*Zdd22/CW*Zn51+3D0*Mq/MW
     &      /SW*Zdd21*Zn53/cb)
      AAABR(278) = EE/MW*Mq/SW
      AAABR(279) = EE/MW*Mb/SW
      AAABR(280) = EE/MW*Mq*Pa21/SW/cb*(sb-Td2/sb)
      AAABR(281) = EE/MW*Mq*Pa11/SW/cb*(sb-Td2/sb)
      AAABR(282) = EE/MW*Mb*Pa21/SW/cb*(sb-Td3/sb)
      AAABR(283) = EE/MW*Mb*Pa11/SW/cb*(sb-Td3/sb)
      AAABR(284) = EE/MW*Mq/SW*(Zh21/cb-Td2*Zh21/cb+Td2*Zh22
     &      /sb)
      AAABR(285) = EE/MW*Mq/SW*(Zh11/cb-Td2*Zh11/cb+Td2*Zh12
     &      /sb)
      AAABR(286) = EE/MW*Mq/SW*(Zh31/cb-Td2*Zh31/cb+Td2*Zh32
     &      /sb)
      AAABR(287) = EE/MW*Mb/SW*(Zh11/cb-Td3*Zh11/cb+Td3*Zh12
     &      /sb)
      AAABR(288) = EE/MW*Mb/SW*(Zh31/cb-Td3*Zh31/cb+Td3*Zh32
     &      /sb)
      AAABR(289) = EE/MW*Mb/SW*(Zh21/cb-Td3*Zh21/cb+Td3*Zh22
     &      /sb)
      AAABR(290) = EE/MW*Mq/SW*Sqrt2/cb*(1D0/sb-sb)
      AAABR(291) = EE/MW*Mq/SW*Sqrt2/cb*(sb-Td2/sb)
      AAABR(292) = EE/MW*Mc/SW*Sqrt2/cb*(1D0/sb-sb)
      AAABR(293) = EE/MW*Mb/SW*Sqrt2/cb*(sb-Td3/sb)
      AAABR(294) = EE/MW*Mt/SW*Sqrt2/cb*(1D0/sb-sb)
      AAABR(295) = EE/MW*Mc/SW*Sqrt2
      AAABR(296) = EE/MW*Mq/SW*Sqrt2
      AAABR(297) = EE/MW*Mb/SW*Sqrt2
      AAABR(298) = EE/MW*Mt/SW*Sqrt2
      AAABR(299) = EE/SW*(2D0*Zt11*Zv11-MtMM/MW*Sqrt2*Zt12
     &      *Zv12/sb)
      AAABR(300) = EE/MW*MbMM/SW*Sqrt2*Zt11*Zu12/cb
      AAABR(301) = EE/SW*(2D0*Zuu11*Zv11-Mq/MW*Sqrt2*Zuu12
     &      *Zv12/sb)
      AAABR(302) = EE/MW*Mq/SW*Sqrt2*Zu12*Zuu11/cb
      AAABR(303) = EE/SW*(2D0*Zcc11*Zv11-Mc/MW*Sqrt2*Zcc12
     &      *Zv12/sb)
      AAABR(304) = EE/MW*Mq/SW*Sqrt2*Zcc11*Zu12/cb
      AAABR(305) = EE/SW*(2D0*Zcc21*Zv11-Mc/MW*Sqrt2*Zcc22
     &      *Zv12/sb)
      AAABR(306) = EE/MW*Mq/SW*Sqrt2*Zcc21*Zu12/cb
      AAABR(307) = EE/SW*(2D0*Zuu21*Zv11-Mq/MW*Sqrt2*Zuu22
     &      *Zv12/sb)
      AAABR(308) = EE/MW*Mq/SW*Sqrt2*Zu12*Zuu21/cb
      AAABR(309) = EE/SW*(2D0*Zt21*Zv11-MtMM/MW*Sqrt2*Zt22
     &      *Zv12/sb)
      AAABR(310) = EE/MW*MbMM/SW*Sqrt2*Zt21*Zu12/cb
      AAABR(311) = EE/SW*(2D0*Zt11*Zv21-MtMM/MW*Sqrt2*Zt12
     &      *Zv22/sb)
      AAABR(312) = EE/MW*MbMM/SW*Sqrt2*Zt11*Zu22/cb
      AAABR(313) = EE/SW*(2D0*Zcc11*Zv21-Mc/MW*Sqrt2*Zcc12
     &      *Zv22/sb)
      AAABR(314) = EE/MW*Mq/SW*Sqrt2*Zcc11*Zu22/cb
      AAABR(315) = EE/SW*(2D0*Zuu11*Zv21-Mq/MW*Sqrt2*Zuu12
     &      *Zv22/sb)
      AAABR(316) = EE/MW*Mq/SW*Sqrt2*Zu22*Zuu11/cb
      AAABR(317) = EE/SW*(2D0*Zt21*Zv21-MtMM/MW*Sqrt2*Zt22
     &      *Zv22/sb)
      AAABR(318) = EE/MW*MbMM/SW*Sqrt2*Zt21*Zu22/cb
      AAABR(319) = EE/SW*(2D0*Zuu21*Zv21-Mq/MW*Sqrt2*Zuu22
     &      *Zv22/sb)
      AAABR(320) = EE/MW*Mq/SW*Sqrt2*Zu22*Zuu21/cb
      AAABR(321) = EE/SW*(2D0*Zcc21*Zv21-Mc/MW*Sqrt2*Zcc22
     &      *Zv22/sb)
      AAABR(322) = EE/MW*Mq/SW*Sqrt2*Zcc21*Zu22/cb
      AAABR(323) = GG*Sqrt2*Zdd11
      AAABR(324) = GG*Sqrt2*Zdd12
      AAABR(325) = GG*Sqrt2*Zss11
      AAABR(326) = GG*Sqrt2*Zss12
      AAABR(327) = GG*Sqrt2*Zb11
      AAABR(328) = GG*Sqrt2*Zb12
      AAABR(329) = GG*Sqrt2*Zss21
      AAABR(330) = GG*Sqrt2*Zss22
      AAABR(331) = GG*Sqrt2*Zb21
      AAABR(332) = GG*Sqrt2*Zb22
      AAABR(333) = GG*Sqrt2*Zdd21
      AAABR(334) = GG*Sqrt2*Zdd22
      AAABC(160) = EE*Sqrt2*(dconjg(Zn11)/CW*Zt11+3D0*dconjg(Zn12
     &      )/SW*Zt11+3D0*MtMM/MW/SW*dconjg(Zn14)*Zt12/sb)
      AAABC(161) = EE*Sqrt2*(4D0*Zn11/CW*Zt12-3D0*MtMM/MW
     &      /SW*Zn14*Zt11/sb)
      AAABC(162) = EE*Sqrt2*(dconjg(Zn11)/CW*Zuu11+3D0*dconjg(Zn12
     &      )/SW*Zuu11+3D0*Mq/MW/SW*dconjg(Zn14)*Zuu12/sb)
      AAABC(163) = EE*Sqrt2*(4D0*Zn11/CW*Zuu12-3D0*Mq/MW
     &      /SW*Zn14*Zuu11/sb)
      AAABC(164) = EE*Sqrt2*(Zcc11/CW*dconjg(Zn21)+3D0*Zcc11
     &      /SW*dconjg(Zn22)+3D0*Mc/MW/SW*Zcc12*dconjg(Zn24
     &      )/sb)
      AAABC(165) = EE*Sqrt2*(4D0*Zcc12/CW*Zn21-3D0*Mc/MW
     &      /SW*Zcc11*Zn24/sb)
      AAABC(166) = EE*Sqrt2*(dconjg(Zn51)/CW*Zt11+3D0*dconjg(Zn52
     &      )/SW*Zt11+3D0*MtMM/MW/SW*dconjg(Zn54)*Zt12/sb)
      AAABC(167) = EE*Sqrt2*(4D0*Zn51/CW*Zt12-3D0*MtMM/MW
     &      /SW*Zn54*Zt11/sb)
      AAABC(168) = EE*Sqrt2*(dconjg(Zn21)/CW*Zt11+3D0*dconjg(Zn22
     &      )/SW*Zt11+3D0*MtMM/MW/SW*dconjg(Zn24)*Zt12/sb)
      AAABC(169) = EE*Sqrt2*(4D0*Zn21/CW*Zt12-3D0*MtMM/MW
     &      /SW*Zn24*Zt11/sb)
      AAABC(170) = EE*Sqrt2*(Zcc11/CW*dconjg(Zn41)+3D0*Zcc11
     &      /SW*dconjg(Zn42)+3D0*Mc/MW/SW*Zcc12*dconjg(Zn44
     &      )/sb)
      AAABC(171) = EE*Sqrt2*(4D0*Zcc12/CW*Zn41-3D0*Mc/MW
     &      /SW*Zcc11*Zn44/sb)
      AAABC(172) = EE*Sqrt2*(Zcc11/CW*dconjg(Zn51)+3D0*Zcc11
     &      /SW*dconjg(Zn52)+3D0*Mc/MW/SW*Zcc12*dconjg(Zn54
     &      )/sb)
      AAABC(173) = EE*Sqrt2*(4D0*Zcc12/CW*Zn51-3D0*Mc/MW
     &      /SW*Zcc11*Zn54/sb)
      AAABC(174) = EE*Sqrt2*(Zcc11/CW*dconjg(Zn11)+3D0*Zcc11
     &      /SW*dconjg(Zn12)+3D0*Mc/MW/SW*Zcc12*dconjg(Zn14
     &      )/sb)
      AAABC(175) = EE*Sqrt2*(4D0*Zcc12/CW*Zn11-3D0*Mc/MW
     &      /SW*Zcc11*Zn14/sb)
      AAABC(176) = EE*Sqrt2*(dconjg(Zn31)/CW*Zt11+3D0*dconjg(Zn32
     &      )/SW*Zt11+3D0*MtMM/MW/SW*dconjg(Zn34)*Zt12/sb)
      AAABC(177) = EE*Sqrt2*(4D0*Zn31/CW*Zt12-3D0*MtMM/MW
     &      /SW*Zn34*Zt11/sb)
      AAABC(178) = EE*Sqrt2*(dconjg(Zn41)/CW*Zt11+3D0*dconjg(Zn42
     &      )/SW*Zt11+3D0*MtMM/MW/SW*dconjg(Zn44)*Zt12/sb)
      AAABC(179) = EE*Sqrt2*(4D0*Zn41/CW*Zt12-3D0*MtMM/MW
     &      /SW*Zn44*Zt11/sb)
      AAABC(180) = EE*Sqrt2*(dconjg(Zn51)/CW*Zuu11+3D0*dconjg(Zn52
     &      )/SW*Zuu11+3D0*Mq/MW/SW*dconjg(Zn54)*Zuu12/sb)
      AAABC(181) = EE*Sqrt2*(4D0*Zn51/CW*Zuu12-3D0*Mq/MW
     &      /SW*Zn54*Zuu11/sb)
      AAABC(182) = EE*Sqrt2*(Zcc11/CW*dconjg(Zn31)+3D0*Zcc11
     &      /SW*dconjg(Zn32)+3D0*Mc/MW/SW*Zcc12*dconjg(Zn34
     &      )/sb)
      AAABC(183) = EE*Sqrt2*(4D0*Zcc12/CW*Zn31-3D0*Mc/MW
     &      /SW*Zcc11*Zn34/sb)
      AAABC(184) = EE*Sqrt2*(dconjg(Zn41)/CW*Zuu11+3D0*dconjg(Zn42
     &      )/SW*Zuu11+3D0*Mq/MW/SW*dconjg(Zn44)*Zuu12/sb)
      AAABC(185) = EE*Sqrt2*(4D0*Zn41/CW*Zuu12-3D0*Mq/MW
     &      /SW*Zn44*Zuu11/sb)
      AAABC(186) = EE*Sqrt2*(dconjg(Zn31)/CW*Zuu11+3D0*dconjg(Zn32
     &      )/SW*Zuu11+3D0*Mq/MW/SW*dconjg(Zn34)*Zuu12/sb)
      AAABC(187) = EE*Sqrt2*(4D0*Zn31/CW*Zuu12-3D0*Mq/MW
     &      /SW*Zn34*Zuu11/sb)
      AAABC(188) = EE*Sqrt2*(dconjg(Zn21)/CW*Zuu11+3D0*dconjg(Zn22
     &      )/SW*Zuu11+3D0*Mq/MW/SW*dconjg(Zn24)*Zuu12/sb)
      AAABC(189) = EE*Sqrt2*(4D0*Zn21/CW*Zuu12-3D0*Mq/MW
     &      /SW*Zn24*Zuu11/sb)
      AAABC(190) = EE*Sqrt2*(dconjg(Zn11)/CW*Zuu21+3D0*dconjg(Zn12
     &      )/SW*Zuu21+3D0*Mq/MW/SW*dconjg(Zn14)*Zuu22/sb)
      AAABC(191) = EE*Sqrt2*(4D0*Zn11/CW*Zuu22-3D0*Mq/MW
     &      /SW*Zn14*Zuu21/sb)
      AAABC(192) = EE*Sqrt2*(Zcc21/CW*dconjg(Zn21)+3D0*Zcc21
     &      /SW*dconjg(Zn22)+3D0*Mc/MW/SW*Zcc22*dconjg(Zn24
     &      )/sb)
      AAABC(193) = EE*Sqrt2*(4D0*Zcc22/CW*Zn21-3D0*Mc/MW
     &      /SW*Zcc21*Zn24/sb)
      AAABC(194) = EE*Sqrt2*(Zcc21/CW*dconjg(Zn51)+3D0*Zcc21
     &      /SW*dconjg(Zn52)+3D0*Mc/MW/SW*Zcc22*dconjg(Zn54
     &      )/sb)
      AAABC(195) = EE*Sqrt2*(4D0*Zcc22/CW*Zn51-3D0*Mc/MW
     &      /SW*Zcc21*Zn54/sb)
      AAABC(196) = EE*Sqrt2*(dconjg(Zn51)/CW*Zuu21+3D0*dconjg(Zn52
     &      )/SW*Zuu21+3D0*Mq/MW/SW*dconjg(Zn54)*Zuu22/sb)
      AAABC(197) = EE*Sqrt2*(4D0*Zn51/CW*Zuu22-3D0*Mq/MW
     &      /SW*Zn54*Zuu21/sb)
      AAABC(198) = EE*Sqrt2*(dconjg(Zn21)/CW*Zuu21+3D0*dconjg(Zn22
     &      )/SW*Zuu21+3D0*Mq/MW/SW*dconjg(Zn24)*Zuu22/sb)
      AAABC(199) = EE*Sqrt2*(4D0*Zn21/CW*Zuu22-3D0*Mq/MW
     &      /SW*Zn24*Zuu21/sb)
      AAABC(200) = EE*Sqrt2*(Zcc21/CW*dconjg(Zn41)+3D0*Zcc21
     &      /SW*dconjg(Zn42)+3D0*Mc/MW/SW*Zcc22*dconjg(Zn44
     &      )/sb)
      AAABC(201) = EE*Sqrt2*(4D0*Zcc22/CW*Zn41-3D0*Mc/MW
     &      /SW*Zcc21*Zn44/sb)
      AAABC(202) = EE*Sqrt2*(Zcc21/CW*dconjg(Zn31)+3D0*Zcc21
     &      /SW*dconjg(Zn32)+3D0*Mc/MW/SW*Zcc22*dconjg(Zn34
     &      )/sb)
      AAABC(203) = EE*Sqrt2*(4D0*Zcc22/CW*Zn31-3D0*Mc/MW
     &      /SW*Zcc21*Zn34/sb)
      AAABC(204) = EE*Sqrt2*(dconjg(Zn31)/CW*Zt21+3D0*dconjg(Zn32
     &      )/SW*Zt21+3D0*MtMM/MW/SW*dconjg(Zn34)*Zt22/sb)
      AAABC(205) = EE*Sqrt2*(4D0*Zn31/CW*Zt22-3D0*MtMM/MW
     &      /SW*Zn34*Zt21/sb)
      AAABC(206) = EE*Sqrt2*(dconjg(Zn51)/CW*Zt21+3D0*dconjg(Zn52
     &      )/SW*Zt21+3D0*MtMM/MW/SW*dconjg(Zn54)*Zt22/sb)
      AAABC(207) = EE*Sqrt2*(4D0*Zn51/CW*Zt22-3D0*MtMM/MW
     &      /SW*Zn54*Zt21/sb)
      AAABC(208) = EE*Sqrt2*(Zcc21/CW*dconjg(Zn11)+3D0*Zcc21
     &      /SW*dconjg(Zn12)+3D0*Mc/MW/SW*Zcc22*dconjg(Zn14
     &      )/sb)
      AAABC(209) = EE*Sqrt2*(4D0*Zcc22/CW*Zn11-3D0*Mc/MW
     &      /SW*Zcc21*Zn14/sb)
      AAABC(210) = EE*Sqrt2*(dconjg(Zn11)/CW*Zt21+3D0*dconjg(Zn12
     &      )/SW*Zt21+3D0*MtMM/MW/SW*dconjg(Zn14)*Zt22/sb)
      AAABC(211) = EE*Sqrt2*(4D0*Zn11/CW*Zt22-3D0*MtMM/MW
     &      /SW*Zn14*Zt21/sb)
      AAABC(212) = EE*Sqrt2*(dconjg(Zn41)/CW*Zuu21+3D0*dconjg(Zn42
     &      )/SW*Zuu21+3D0*Mq/MW/SW*dconjg(Zn44)*Zuu22/sb)
      AAABC(213) = EE*Sqrt2*(4D0*Zn41/CW*Zuu22-3D0*Mq/MW
     &      /SW*Zn44*Zuu21/sb)
      AAABC(214) = EE*Sqrt2*(dconjg(Zn21)/CW*Zt21+3D0*dconjg(Zn22
     &      )/SW*Zt21+3D0*MtMM/MW/SW*dconjg(Zn24)*Zt22/sb)
      AAABC(215) = EE*Sqrt2*(4D0*Zn21/CW*Zt22-3D0*MtMM/MW
     &      /SW*Zn24*Zt21/sb)
      AAABC(216) = EE*Sqrt2*(dconjg(Zn31)/CW*Zuu21+3D0*dconjg(Zn32
     &      )/SW*Zuu21+3D0*Mq/MW/SW*dconjg(Zn34)*Zuu22/sb)
      AAABC(217) = EE*Sqrt2*(4D0*Zn31/CW*Zuu22-3D0*Mq/MW
     &      /SW*Zn34*Zuu21/sb)
      AAABC(218) = EE*Sqrt2*(dconjg(Zn41)/CW*Zt21+3D0*dconjg(Zn42
     &      )/SW*Zt21+3D0*MtMM/MW/SW*dconjg(Zn44)*Zt22/sb)
      AAABC(219) = EE*Sqrt2*(4D0*Zn41/CW*Zt22-3D0*MtMM/MW
     &      /SW*Zn44*Zt21/sb)
      AAABR(335) = EE/MW*Mt/SW
      AAABR(336) = EE/MW*Mc/SW
      AAABR(337) = EE/MW*Mt*Pa11/SW/tB
      AAABR(338) = EE/MW*Mt*Pa21/SW/tB
      AAABR(339) = EE/MW*Mq*Pa21/SW/tB
      AAABR(340) = EE/MW*Mc*Pa11/SW/tB
      AAABR(341) = EE/MW*Mc*Pa21/SW/tB
      AAABR(342) = EE/MW*Mq*Pa11/SW/tB
      AAABR(343) = EE/MW*Mt/SW*Zh32/sb
      AAABR(344) = EE/MW*Mq/SW*Zh32/sb
      AAABR(345) = EE/MW*Mq/SW*Zh12/sb
      AAABR(346) = EE/MW*Mc/SW*Zh12/sb
      AAABR(347) = EE/MW*Mc/SW*Zh22/sb
      AAABR(348) = EE/MW*Mt/SW*Zh22/sb
      AAABR(349) = EE/MW*Mc/SW*Zh32/sb
      AAABR(350) = EE/MW*Mq/SW*Zh22/sb
      AAABR(351) = EE/MW*Mt/SW*Zh12/sb
      AAABR(352) = GG*Sqrt2*Zt11
      AAABR(353) = GG*Sqrt2*Zt12
      AAABR(354) = GG*Sqrt2*Zuu11
      AAABR(355) = GG*Sqrt2*Zuu12
      AAABR(356) = GG*Sqrt2*Zcc11
      AAABR(357) = GG*Sqrt2*Zcc12
      AAABR(358) = GG*Sqrt2*Zt21
      AAABR(359) = GG*Sqrt2*Zt22
      AAABR(360) = GG*Sqrt2*Zuu21
      AAABR(361) = GG*Sqrt2*Zuu22
      AAABR(362) = GG*Sqrt2*Zcc21
      AAABR(363) = GG*Sqrt2*Zcc22
      AAABC(220) = EE/MW/SW*(NMM13*Zn11*Zn13-NMM14*Zn11*Zn14
     &      +NMM23*Zn12*Zn13-NMM24*Zn12*Zn14+NMM35*Zn13*Zn15
     &      -NMM45*Zn14*Zn15)
      AAABC(221) = EE/MW/SW*(NMM13*Zn11*Zn43+NMM13*Zn13*Zn41
     &      -NMM14*Zn11*Zn44-NMM14*Zn14*Zn41+NMM23*Zn12*Zn43
     &      +NMM23*Zn13*Zn42-NMM24*Zn12*Zn44-NMM24*Zn14*Zn42
     &      +NMM35*Zn13*Zn45+NMM35*Zn15*Zn43-NMM45*Zn14*Zn45
     &      -NMM45*Zn15*Zn44)
      AAABC(222) = EE/MW/SW*(NMM13*Zn11*Zn53+NMM13*Zn13*Zn51
     &      -NMM14*Zn11*Zn54-NMM14*Zn14*Zn51+NMM23*Zn12*Zn53
     &      +NMM23*Zn13*Zn52-NMM24*Zn12*Zn54-NMM24*Zn14*Zn52
     &      +NMM35*Zn13*Zn55+NMM35*Zn15*Zn53-NMM45*Zn14*Zn55
     &      -NMM45*Zn15*Zn54)
      AAABC(223) = EE/MW/SW*(NMM13*Zn11*Zn23+NMM13*Zn13*Zn21
     &      -NMM14*Zn11*Zn24-NMM14*Zn14*Zn21+NMM23*Zn12*Zn23
     &      +NMM23*Zn13*Zn22-NMM24*Zn12*Zn24-NMM24*Zn14*Zn22
     &      +NMM35*Zn13*Zn25+NMM35*Zn15*Zn23-NMM45*Zn14*Zn25
     &      -NMM45*Zn15*Zn24)
      AAABC(224) = EE/MW/SW*(NMM13*Zn11*Zn33+NMM13*Zn13*Zn31
     &      -NMM14*Zn11*Zn34-NMM14*Zn14*Zn31+NMM23*Zn12*Zn33
     &      +NMM23*Zn13*Zn32-NMM24*Zn12*Zn34-NMM24*Zn14*Zn32
     &      +NMM35*Zn13*Zn35+NMM35*Zn15*Zn33-NMM45*Zn14*Zn35
     &      -NMM45*Zn15*Zn34)
      AAABC(225) = NMM34*Pa12*Sqrt2*Zn13*Zn54/xvev+NMM34
     &      *Pa12*Sqrt2*Zn14*Zn53/xvev+NMM55*Pa12*Sqrt2*Zn15
     &      *Zn55/xvev-Pa12*Sqrt2*Zn15*Zn55*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn11*Zn53/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn13*Zn51/cb*sb-EE/MW*NMM14*Pa11/SW*Zn11
     &      *Zn54/cb/sb+EE/MW*NMM14*Pa11/SW*Zn11*Zn54/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn14*Zn51/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn14*Zn51/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn12*Zn53/cb*sb-EE/MW*NMM23*Pa11/SW*Zn13
     &      *Zn52/cb*sb-EE/MW*NMM24*Pa11/SW*Zn12*Zn54/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn12*Zn54/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn14*Zn52/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn14*Zn52/cb*sb+EE/MW*NMM35*Pa11/SW*Zn13
     &      *Zn55/cb/sb-EE/MW*NMM35*Pa11/SW*Zn13*Zn55/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn15*Zn53/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn15*Zn53/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn14*Zn55/cb*sb+EE/MW*NMM45*Pa11/SW*Zn15
     &      *Zn54/cb*sb
      AAABC(226) = NMM34*Pa22*Sqrt2*Zn13*Zn24/xvev+NMM34
     &      *Pa22*Sqrt2*Zn14*Zn23/xvev+NMM55*Pa22*Sqrt2*Zn15
     &      *Zn25/xvev-Pa22*Sqrt2*Zn15*Zn25*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn11*Zn23/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn13*Zn21/cb*sb-EE/MW*NMM14*Pa21/SW*Zn11
     &      *Zn24/cb/sb+EE/MW*NMM14*Pa21/SW*Zn11*Zn24/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn14*Zn21/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn14*Zn21/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn12*Zn23/cb*sb-EE/MW*NMM23*Pa21/SW*Zn13
     &      *Zn22/cb*sb-EE/MW*NMM24*Pa21/SW*Zn12*Zn24/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn12*Zn24/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn14*Zn22/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn14*Zn22/cb*sb+EE/MW*NMM35*Pa21/SW*Zn13
     &      *Zn25/cb/sb-EE/MW*NMM35*Pa21/SW*Zn13*Zn25/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn15*Zn23/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn15*Zn23/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn14*Zn25/cb*sb+EE/MW*NMM45*Pa21/SW*Zn15
     &      *Zn24/cb*sb
      AAABC(227) = NMM34*Pa22*Sqrt2*Zn13*Zn54/xvev+NMM34
     &      *Pa22*Sqrt2*Zn14*Zn53/xvev+NMM55*Pa22*Sqrt2*Zn15
     &      *Zn55/xvev-Pa22*Sqrt2*Zn15*Zn55*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn11*Zn53/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn13*Zn51/cb*sb-EE/MW*NMM14*Pa21/SW*Zn11
     &      *Zn54/cb/sb+EE/MW*NMM14*Pa21/SW*Zn11*Zn54/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn14*Zn51/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn14*Zn51/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn12*Zn53/cb*sb-EE/MW*NMM23*Pa21/SW*Zn13
     &      *Zn52/cb*sb-EE/MW*NMM24*Pa21/SW*Zn12*Zn54/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn12*Zn54/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn14*Zn52/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn14*Zn52/cb*sb+EE/MW*NMM35*Pa21/SW*Zn13
     &      *Zn55/cb/sb-EE/MW*NMM35*Pa21/SW*Zn13*Zn55/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn15*Zn53/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn15*Zn53/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn14*Zn55/cb*sb+EE/MW*NMM45*Pa21/SW*Zn15
     &      *Zn54/cb*sb
      AAABC(228) = NMM34*Pa12*Sqrt2*Zn13*Zn34/xvev+NMM34
     &      *Pa12*Sqrt2*Zn14*Zn33/xvev+NMM55*Pa12*Sqrt2*Zn15
     &      *Zn35/xvev-Pa12*Sqrt2*Zn15*Zn35*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn11*Zn33/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn13*Zn31/cb*sb-EE/MW*NMM14*Pa11/SW*Zn11
     &      *Zn34/cb/sb+EE/MW*NMM14*Pa11/SW*Zn11*Zn34/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn14*Zn31/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn14*Zn31/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn12*Zn33/cb*sb-EE/MW*NMM23*Pa11/SW*Zn13
     &      *Zn32/cb*sb-EE/MW*NMM24*Pa11/SW*Zn12*Zn34/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn12*Zn34/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn14*Zn32/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn14*Zn32/cb*sb+EE/MW*NMM35*Pa11/SW*Zn13
     &      *Zn35/cb/sb-EE/MW*NMM35*Pa11/SW*Zn13*Zn35/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn15*Zn33/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn15*Zn33/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn14*Zn35/cb*sb+EE/MW*NMM45*Pa11/SW*Zn15
     &      *Zn34/cb*sb
      AAABC(229) = NMM55*Pa12*Sqrt2*Zn15**2/xvev-Pa12*Sqrt2
     &      *Zn15**2*muP/xvev+2D0*NMM34*Pa12*Sqrt2*Zn13*Zn14
     &      /xvev-2D0*EE/MW*NMM13*Pa11/SW*Zn11*Zn13/cb*sb
     &      -2D0*EE/MW*NMM14*Pa11/SW*Zn11*Zn14/cb/sb+2D0
     &      *EE/MW*NMM14*Pa11/SW*Zn11*Zn14/cb*sb-2D0*EE/MW
     &      *NMM23*Pa11/SW*Zn12*Zn13/cb*sb-2D0*EE/MW*NMM24
     &      *Pa11/SW*Zn12*Zn14/cb/sb+2D0*EE/MW*NMM24*Pa11
     &      /SW*Zn12*Zn14/cb*sb+2D0*EE/MW*NMM35*Pa11/SW*Zn13
     &      *Zn15/cb/sb-2D0*EE/MW*NMM35*Pa11/SW*Zn13*Zn15
     &      /cb*sb+2D0*EE/MW*NMM45*Pa11/SW*Zn14*Zn15/cb*sb
      AAABC(230) = NMM34*Pa12*Sqrt2*Zn13*Zn44/xvev+NMM34
     &      *Pa12*Sqrt2*Zn14*Zn43/xvev+NMM55*Pa12*Sqrt2*Zn15
     &      *Zn45/xvev-Pa12*Sqrt2*Zn15*Zn45*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn11*Zn43/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn13*Zn41/cb*sb-EE/MW*NMM14*Pa11/SW*Zn11
     &      *Zn44/cb/sb+EE/MW*NMM14*Pa11/SW*Zn11*Zn44/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn14*Zn41/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn14*Zn41/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn12*Zn43/cb*sb-EE/MW*NMM23*Pa11/SW*Zn13
     &      *Zn42/cb*sb-EE/MW*NMM24*Pa11/SW*Zn12*Zn44/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn12*Zn44/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn14*Zn42/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn14*Zn42/cb*sb+EE/MW*NMM35*Pa11/SW*Zn13
     &      *Zn45/cb/sb-EE/MW*NMM35*Pa11/SW*Zn13*Zn45/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn15*Zn43/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn15*Zn43/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn14*Zn45/cb*sb+EE/MW*NMM45*Pa11/SW*Zn15
     &      *Zn44/cb*sb
      AAABC(231) = NMM55*Pa22*Sqrt2*Zn15**2/xvev-Pa22*Sqrt2
     &      *Zn15**2*muP/xvev+2D0*NMM34*Pa22*Sqrt2*Zn13*Zn14
     &      /xvev-2D0*EE/MW*NMM13*Pa21/SW*Zn11*Zn13/cb*sb
     &      -2D0*EE/MW*NMM14*Pa21/SW*Zn11*Zn14/cb/sb+2D0
     &      *EE/MW*NMM14*Pa21/SW*Zn11*Zn14/cb*sb-2D0*EE/MW
     &      *NMM23*Pa21/SW*Zn12*Zn13/cb*sb-2D0*EE/MW*NMM24
     &      *Pa21/SW*Zn12*Zn14/cb/sb+2D0*EE/MW*NMM24*Pa21
     &      /SW*Zn12*Zn14/cb*sb+2D0*EE/MW*NMM35*Pa21/SW*Zn13
     &      *Zn15/cb/sb-2D0*EE/MW*NMM35*Pa21/SW*Zn13*Zn15
     &      /cb*sb+2D0*EE/MW*NMM45*Pa21/SW*Zn14*Zn15/cb*sb
      AAABC(232) = NMM34*Pa22*Sqrt2*Zn13*Zn34/xvev+NMM34
     &      *Pa22*Sqrt2*Zn14*Zn33/xvev+NMM55*Pa22*Sqrt2*Zn15
     &      *Zn35/xvev-Pa22*Sqrt2*Zn15*Zn35*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn11*Zn33/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn13*Zn31/cb*sb-EE/MW*NMM14*Pa21/SW*Zn11
     &      *Zn34/cb/sb+EE/MW*NMM14*Pa21/SW*Zn11*Zn34/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn14*Zn31/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn14*Zn31/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn12*Zn33/cb*sb-EE/MW*NMM23*Pa21/SW*Zn13
     &      *Zn32/cb*sb-EE/MW*NMM24*Pa21/SW*Zn12*Zn34/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn12*Zn34/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn14*Zn32/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn14*Zn32/cb*sb+EE/MW*NMM35*Pa21/SW*Zn13
     &      *Zn35/cb/sb-EE/MW*NMM35*Pa21/SW*Zn13*Zn35/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn15*Zn33/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn15*Zn33/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn14*Zn35/cb*sb+EE/MW*NMM45*Pa21/SW*Zn15
     &      *Zn34/cb*sb
      AAABC(233) = NMM34*Pa22*Sqrt2*Zn13*Zn44/xvev+NMM34
     &      *Pa22*Sqrt2*Zn14*Zn43/xvev+NMM55*Pa22*Sqrt2*Zn15
     &      *Zn45/xvev-Pa22*Sqrt2*Zn15*Zn45*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn11*Zn43/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn13*Zn41/cb*sb-EE/MW*NMM14*Pa21/SW*Zn11
     &      *Zn44/cb/sb+EE/MW*NMM14*Pa21/SW*Zn11*Zn44/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn14*Zn41/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn14*Zn41/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn12*Zn43/cb*sb-EE/MW*NMM23*Pa21/SW*Zn13
     &      *Zn42/cb*sb-EE/MW*NMM24*Pa21/SW*Zn12*Zn44/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn12*Zn44/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn14*Zn42/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn14*Zn42/cb*sb+EE/MW*NMM35*Pa21/SW*Zn13
     &      *Zn45/cb/sb-EE/MW*NMM35*Pa21/SW*Zn13*Zn45/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn15*Zn43/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn15*Zn43/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn14*Zn45/cb*sb+EE/MW*NMM45*Pa21/SW*Zn15
     &      *Zn44/cb*sb
      AAABC(234) = NMM34*Pa12*Sqrt2*Zn13*Zn24/xvev+NMM34
     &      *Pa12*Sqrt2*Zn14*Zn23/xvev+NMM55*Pa12*Sqrt2*Zn15
     &      *Zn25/xvev-Pa12*Sqrt2*Zn15*Zn25*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn11*Zn23/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn13*Zn21/cb*sb-EE/MW*NMM14*Pa11/SW*Zn11
     &      *Zn24/cb/sb+EE/MW*NMM14*Pa11/SW*Zn11*Zn24/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn14*Zn21/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn14*Zn21/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn12*Zn23/cb*sb-EE/MW*NMM23*Pa11/SW*Zn13
     &      *Zn22/cb*sb-EE/MW*NMM24*Pa11/SW*Zn12*Zn24/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn12*Zn24/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn14*Zn22/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn14*Zn22/cb*sb+EE/MW*NMM35*Pa11/SW*Zn13
     &      *Zn25/cb/sb-EE/MW*NMM35*Pa11/SW*Zn13*Zn25/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn15*Zn23/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn15*Zn23/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn14*Zn25/cb*sb+EE/MW*NMM45*Pa11/SW*Zn15
     &      *Zn24/cb*sb
      AAABC(235) = NMM34*Sqrt2*Zh33*Zn13*Zn54/xvev+NMM34
     &      *Sqrt2*Zh33*Zn14*Zn53/xvev+NMM55*Sqrt2*Zh33*Zn15
     &      *Zn55/xvev-Sqrt2*Zh33*Zn15*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn11*Zn53/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn13*Zn51/cb+EE/MW*NMM14/SW*Zh32*Zn11*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn14*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn12*Zn53/cb+EE/MW*NMM23/SW*Zh31*Zn13
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh32*Zn12*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn14*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn13*Zn55/sb+EE/MW*NMM35/SW*Zh32*Zn15*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn14*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn15*Zn54/cb
      AAABC(236) = NMM34*Sqrt2*Zh33*Zn13*Zn34/xvev+NMM34
     &      *Sqrt2*Zh33*Zn14*Zn33/xvev+NMM55*Sqrt2*Zh33*Zn15
     &      *Zn35/xvev-Sqrt2*Zh33*Zn15*Zn35*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn11*Zn33/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn13*Zn31/cb+EE/MW*NMM14/SW*Zh32*Zn11*Zn34/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn14*Zn31/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn12*Zn33/cb+EE/MW*NMM23/SW*Zh31*Zn13
     &      *Zn32/cb+EE/MW*NMM24/SW*Zh32*Zn12*Zn34/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn14*Zn32/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn13*Zn35/sb+EE/MW*NMM35/SW*Zh32*Zn15*Zn33
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn14*Zn35/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn15*Zn34/cb
      AAABC(237) = NMM34*Sqrt2*Zh13*Zn13*Zn34/xvev+NMM34
     &      *Sqrt2*Zh13*Zn14*Zn33/xvev+NMM55*Sqrt2*Zh13*Zn15
     &      *Zn35/xvev-Sqrt2*Zh13*Zn15*Zn35*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn11*Zn33/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn13*Zn31/cb+EE/MW*NMM14/SW*Zh12*Zn11*Zn34/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn14*Zn31/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn12*Zn33/cb+EE/MW*NMM23/SW*Zh11*Zn13
     &      *Zn32/cb+EE/MW*NMM24/SW*Zh12*Zn12*Zn34/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn14*Zn32/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn13*Zn35/sb+EE/MW*NMM35/SW*Zh12*Zn15*Zn33
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn14*Zn35/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn15*Zn34/cb
      AAABC(238) = NMM34*Sqrt2*Zh23*Zn13*Zn54/xvev+NMM34
     &      *Sqrt2*Zh23*Zn14*Zn53/xvev+NMM55*Sqrt2*Zh23*Zn15
     &      *Zn55/xvev-Sqrt2*Zh23*Zn15*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn11*Zn53/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn13*Zn51/cb+EE/MW*NMM14/SW*Zh22*Zn11*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn14*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn12*Zn53/cb+EE/MW*NMM23/SW*Zh21*Zn13
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh22*Zn12*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn14*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn13*Zn55/sb+EE/MW*NMM35/SW*Zh22*Zn15*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn14*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn15*Zn54/cb
      AAABC(239) = NMM34*Sqrt2*Zh23*Zn13*Zn44/xvev+NMM34
     &      *Sqrt2*Zh23*Zn14*Zn43/xvev+NMM55*Sqrt2*Zh23*Zn15
     &      *Zn45/xvev-Sqrt2*Zh23*Zn15*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn11*Zn43/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn13*Zn41/cb+EE/MW*NMM14/SW*Zh22*Zn11*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn14*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn12*Zn43/cb+EE/MW*NMM23/SW*Zh21*Zn13
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh22*Zn12*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn14*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn13*Zn45/sb+EE/MW*NMM35/SW*Zh22*Zn15*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn14*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn15*Zn44/cb
      AAABC(240) = NMM34*Sqrt2*Zh23*Zn13*Zn24/xvev+NMM34
     &      *Sqrt2*Zh23*Zn14*Zn23/xvev+NMM55*Sqrt2*Zh23*Zn15
     &      *Zn25/xvev-Sqrt2*Zh23*Zn15*Zn25*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn11*Zn23/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn13*Zn21/cb+EE/MW*NMM14/SW*Zh22*Zn11*Zn24/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn14*Zn21/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn12*Zn23/cb+EE/MW*NMM23/SW*Zh21*Zn13
     &      *Zn22/cb+EE/MW*NMM24/SW*Zh22*Zn12*Zn24/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn14*Zn22/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn13*Zn25/sb+EE/MW*NMM35/SW*Zh22*Zn15*Zn23
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn14*Zn25/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn15*Zn24/cb
      AAABC(241) = NMM34*Sqrt2*Zh33*Zn13*Zn44/xvev+NMM34
     &      *Sqrt2*Zh33*Zn14*Zn43/xvev+NMM55*Sqrt2*Zh33*Zn15
     &      *Zn45/xvev-Sqrt2*Zh33*Zn15*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn11*Zn43/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn13*Zn41/cb+EE/MW*NMM14/SW*Zh32*Zn11*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn14*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn12*Zn43/cb+EE/MW*NMM23/SW*Zh31*Zn13
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh32*Zn12*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn14*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn13*Zn45/sb+EE/MW*NMM35/SW*Zh32*Zn15*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn14*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn15*Zn44/cb
      AAABC(242) = NMM34*Sqrt2*Zh13*Zn13*Zn24/xvev+NMM34
     &      *Sqrt2*Zh13*Zn14*Zn23/xvev+NMM55*Sqrt2*Zh13*Zn15
     &      *Zn25/xvev-Sqrt2*Zh13*Zn15*Zn25*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn11*Zn23/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn13*Zn21/cb+EE/MW*NMM14/SW*Zh12*Zn11*Zn24/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn14*Zn21/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn12*Zn23/cb+EE/MW*NMM23/SW*Zh11*Zn13
     &      *Zn22/cb+EE/MW*NMM24/SW*Zh12*Zn12*Zn24/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn14*Zn22/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn13*Zn25/sb+EE/MW*NMM35/SW*Zh12*Zn15*Zn23
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn14*Zn25/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn15*Zn24/cb
      AAABC(243) = NMM34*Sqrt2*Zh13*Zn13*Zn54/xvev+NMM34
     &      *Sqrt2*Zh13*Zn14*Zn53/xvev+NMM55*Sqrt2*Zh13*Zn15
     &      *Zn55/xvev-Sqrt2*Zh13*Zn15*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn11*Zn53/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn13*Zn51/cb+EE/MW*NMM14/SW*Zh12*Zn11*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn14*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn12*Zn53/cb+EE/MW*NMM23/SW*Zh11*Zn13
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh12*Zn12*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn14*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn13*Zn55/sb+EE/MW*NMM35/SW*Zh12*Zn15*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn14*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn15*Zn54/cb
      AAABC(244) = NMM55*Sqrt2*Zh33*Zn15**2/xvev-Sqrt2*Zh33
     &      *Zn15**2*muP/xvev+2D0*NMM34*Sqrt2*Zh33*Zn13*Zn14
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh31*Zn11*Zn13/cb+2D0
     &      *EE/MW*NMM14/SW*Zh32*Zn11*Zn14/sb+2D0*EE/MW*NMM23
     &      /SW*Zh31*Zn12*Zn13/cb+2D0*EE/MW*NMM24/SW*Zh32
     &      *Zn12*Zn14/sb+2D0*EE/MW*NMM35/SW*Zh32*Zn13*Zn15
     &      /sb+2D0*EE/MW*NMM45/SW*Zh31*Zn14*Zn15/cb
      AAABC(245) = NMM34*Sqrt2*Zh33*Zn13*Zn24/xvev+NMM34
     &      *Sqrt2*Zh33*Zn14*Zn23/xvev+NMM55*Sqrt2*Zh33*Zn15
     &      *Zn25/xvev-Sqrt2*Zh33*Zn15*Zn25*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn11*Zn23/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn13*Zn21/cb+EE/MW*NMM14/SW*Zh32*Zn11*Zn24/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn14*Zn21/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn12*Zn23/cb+EE/MW*NMM23/SW*Zh31*Zn13
     &      *Zn22/cb+EE/MW*NMM24/SW*Zh32*Zn12*Zn24/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn14*Zn22/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn13*Zn25/sb+EE/MW*NMM35/SW*Zh32*Zn15*Zn23
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn14*Zn25/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn15*Zn24/cb
      AAABC(246) = NMM34*Sqrt2*Zh23*Zn13*Zn34/xvev+NMM34
     &      *Sqrt2*Zh23*Zn14*Zn33/xvev+NMM55*Sqrt2*Zh23*Zn15
     &      *Zn35/xvev-Sqrt2*Zh23*Zn15*Zn35*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn11*Zn33/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn13*Zn31/cb+EE/MW*NMM14/SW*Zh22*Zn11*Zn34/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn14*Zn31/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn12*Zn33/cb+EE/MW*NMM23/SW*Zh21*Zn13
     &      *Zn32/cb+EE/MW*NMM24/SW*Zh22*Zn12*Zn34/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn14*Zn32/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn13*Zn35/sb+EE/MW*NMM35/SW*Zh22*Zn15*Zn33
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn14*Zn35/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn15*Zn34/cb
      AAABC(247) = NMM55*Sqrt2*Zh13*Zn15**2/xvev-Sqrt2*Zh13
     &      *Zn15**2*muP/xvev+2D0*NMM34*Sqrt2*Zh13*Zn13*Zn14
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh11*Zn11*Zn13/cb+2D0
     &      *EE/MW*NMM14/SW*Zh12*Zn11*Zn14/sb+2D0*EE/MW*NMM23
     &      /SW*Zh11*Zn12*Zn13/cb+2D0*EE/MW*NMM24/SW*Zh12
     &      *Zn12*Zn14/sb+2D0*EE/MW*NMM35/SW*Zh12*Zn13*Zn15
     &      /sb+2D0*EE/MW*NMM45/SW*Zh11*Zn14*Zn15/cb
      AAABC(248) = NMM34*Sqrt2*Zh13*Zn13*Zn44/xvev+NMM34
     &      *Sqrt2*Zh13*Zn14*Zn43/xvev+NMM55*Sqrt2*Zh13*Zn15
     &      *Zn45/xvev-Sqrt2*Zh13*Zn15*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn11*Zn43/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn13*Zn41/cb+EE/MW*NMM14/SW*Zh12*Zn11*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn14*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn12*Zn43/cb+EE/MW*NMM23/SW*Zh11*Zn13
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh12*Zn12*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn14*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn13*Zn45/sb+EE/MW*NMM35/SW*Zh12*Zn15*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn14*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn15*Zn44/cb
      AAABC(249) = NMM55*Sqrt2*Zh23*Zn15**2/xvev-Sqrt2*Zh23
     &      *Zn15**2*muP/xvev+2D0*NMM34*Sqrt2*Zh23*Zn13*Zn14
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh21*Zn11*Zn13/cb+2D0
     &      *EE/MW*NMM14/SW*Zh22*Zn11*Zn14/sb+2D0*EE/MW*NMM23
     &      /SW*Zh21*Zn12*Zn13/cb+2D0*EE/MW*NMM24/SW*Zh22
     &      *Zn12*Zn14/sb+2D0*EE/MW*NMM35/SW*Zh22*Zn13*Zn15
     &      /sb+2D0*EE/MW*NMM45/SW*Zh21*Zn14*Zn15/cb
      AAABC(250) = EE/MW/SW*(NMM13*Zn21*Zn33+NMM13*Zn23*Zn31
     &      -NMM14*Zn21*Zn34-NMM14*Zn24*Zn31+NMM23*Zn22*Zn33
     &      +NMM23*Zn23*Zn32-NMM24*Zn22*Zn34-NMM24*Zn24*Zn32
     &      +NMM35*Zn23*Zn35+NMM35*Zn25*Zn33-NMM45*Zn24*Zn35
     &      -NMM45*Zn25*Zn34)
      AAABC(251) = EE/MW/SW*(NMM13*Zn21*Zn43+NMM13*Zn23*Zn41
     &      -NMM14*Zn21*Zn44-NMM14*Zn24*Zn41+NMM23*Zn22*Zn43
     &      +NMM23*Zn23*Zn42-NMM24*Zn22*Zn44-NMM24*Zn24*Zn42
     &      +NMM35*Zn23*Zn45+NMM35*Zn25*Zn43-NMM45*Zn24*Zn45
     &      -NMM45*Zn25*Zn44)
      AAABC(252) = EE/MW/SW*(NMM13*Zn21*Zn23-NMM14*Zn21*Zn24
     &      +NMM23*Zn22*Zn23-NMM24*Zn22*Zn24+NMM35*Zn23*Zn25
     &      -NMM45*Zn24*Zn25)
      AAABC(253) = EE/MW/SW*(NMM13*Zn21*Zn53+NMM13*Zn23*Zn51
     &      -NMM14*Zn21*Zn54-NMM14*Zn24*Zn51+NMM23*Zn22*Zn53
     &      +NMM23*Zn23*Zn52-NMM24*Zn22*Zn54-NMM24*Zn24*Zn52
     &      +NMM35*Zn23*Zn55+NMM35*Zn25*Zn53-NMM45*Zn24*Zn55
     &      -NMM45*Zn25*Zn54)
      AAABC(254) = NMM34*Pa12*Sqrt2*Zn23*Zn34/xvev+NMM34
     &      *Pa12*Sqrt2*Zn24*Zn33/xvev+NMM55*Pa12*Sqrt2*Zn25
     &      *Zn35/xvev-Pa12*Sqrt2*Zn25*Zn35*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn21*Zn33/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn23*Zn31/cb*sb-EE/MW*NMM14*Pa11/SW*Zn21
     &      *Zn34/cb/sb+EE/MW*NMM14*Pa11/SW*Zn21*Zn34/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn24*Zn31/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn24*Zn31/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn22*Zn33/cb*sb-EE/MW*NMM23*Pa11/SW*Zn23
     &      *Zn32/cb*sb-EE/MW*NMM24*Pa11/SW*Zn22*Zn34/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn22*Zn34/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn24*Zn32/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn24*Zn32/cb*sb+EE/MW*NMM35*Pa11/SW*Zn23
     &      *Zn35/cb/sb-EE/MW*NMM35*Pa11/SW*Zn23*Zn35/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn25*Zn33/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn25*Zn33/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn24*Zn35/cb*sb+EE/MW*NMM45*Pa11/SW*Zn25
     &      *Zn34/cb*sb
      AAABC(255) = NMM34*Pa22*Sqrt2*Zn23*Zn44/xvev+NMM34
     &      *Pa22*Sqrt2*Zn24*Zn43/xvev+NMM55*Pa22*Sqrt2*Zn25
     &      *Zn45/xvev-Pa22*Sqrt2*Zn25*Zn45*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn21*Zn43/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn23*Zn41/cb*sb-EE/MW*NMM14*Pa21/SW*Zn21
     &      *Zn44/cb/sb+EE/MW*NMM14*Pa21/SW*Zn21*Zn44/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn24*Zn41/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn24*Zn41/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn22*Zn43/cb*sb-EE/MW*NMM23*Pa21/SW*Zn23
     &      *Zn42/cb*sb-EE/MW*NMM24*Pa21/SW*Zn22*Zn44/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn22*Zn44/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn24*Zn42/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn24*Zn42/cb*sb+EE/MW*NMM35*Pa21/SW*Zn23
     &      *Zn45/cb/sb-EE/MW*NMM35*Pa21/SW*Zn23*Zn45/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn25*Zn43/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn25*Zn43/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn24*Zn45/cb*sb+EE/MW*NMM45*Pa21/SW*Zn25
     &      *Zn44/cb*sb
      AAABC(256) = NMM34*Pa12*Sqrt2*Zn23*Zn54/xvev+NMM34
     &      *Pa12*Sqrt2*Zn24*Zn53/xvev+NMM55*Pa12*Sqrt2*Zn25
     &      *Zn55/xvev-Pa12*Sqrt2*Zn25*Zn55*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn21*Zn53/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn23*Zn51/cb*sb-EE/MW*NMM14*Pa11/SW*Zn21
     &      *Zn54/cb/sb+EE/MW*NMM14*Pa11/SW*Zn21*Zn54/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn24*Zn51/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn24*Zn51/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn22*Zn53/cb*sb-EE/MW*NMM23*Pa11/SW*Zn23
     &      *Zn52/cb*sb-EE/MW*NMM24*Pa11/SW*Zn22*Zn54/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn22*Zn54/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn24*Zn52/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn24*Zn52/cb*sb+EE/MW*NMM35*Pa11/SW*Zn23
     &      *Zn55/cb/sb-EE/MW*NMM35*Pa11/SW*Zn23*Zn55/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn25*Zn53/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn25*Zn53/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn24*Zn55/cb*sb+EE/MW*NMM45*Pa11/SW*Zn25
     &      *Zn54/cb*sb
      AAABC(257) = NMM55*Pa22*Sqrt2*Zn25**2/xvev-Pa22*Sqrt2
     &      *Zn25**2*muP/xvev+2D0*NMM34*Pa22*Sqrt2*Zn23*Zn24
     &      /xvev-2D0*EE/MW*NMM13*Pa21/SW*Zn21*Zn23/cb*sb
     &      -2D0*EE/MW*NMM14*Pa21/SW*Zn21*Zn24/cb/sb+2D0
     &      *EE/MW*NMM14*Pa21/SW*Zn21*Zn24/cb*sb-2D0*EE/MW
     &      *NMM23*Pa21/SW*Zn22*Zn23/cb*sb-2D0*EE/MW*NMM24
     &      *Pa21/SW*Zn22*Zn24/cb/sb+2D0*EE/MW*NMM24*Pa21
     &      /SW*Zn22*Zn24/cb*sb+2D0*EE/MW*NMM35*Pa21/SW*Zn23
     &      *Zn25/cb/sb-2D0*EE/MW*NMM35*Pa21/SW*Zn23*Zn25
     &      /cb*sb+2D0*EE/MW*NMM45*Pa21/SW*Zn24*Zn25/cb*sb
      AAABC(258) = NMM55*Pa12*Sqrt2*Zn25**2/xvev-Pa12*Sqrt2
     &      *Zn25**2*muP/xvev+2D0*NMM34*Pa12*Sqrt2*Zn23*Zn24
     &      /xvev-2D0*EE/MW*NMM13*Pa11/SW*Zn21*Zn23/cb*sb
     &      -2D0*EE/MW*NMM14*Pa11/SW*Zn21*Zn24/cb/sb+2D0
     &      *EE/MW*NMM14*Pa11/SW*Zn21*Zn24/cb*sb-2D0*EE/MW
     &      *NMM23*Pa11/SW*Zn22*Zn23/cb*sb-2D0*EE/MW*NMM24
     &      *Pa11/SW*Zn22*Zn24/cb/sb+2D0*EE/MW*NMM24*Pa11
     &      /SW*Zn22*Zn24/cb*sb+2D0*EE/MW*NMM35*Pa11/SW*Zn23
     &      *Zn25/cb/sb-2D0*EE/MW*NMM35*Pa11/SW*Zn23*Zn25
     &      /cb*sb+2D0*EE/MW*NMM45*Pa11/SW*Zn24*Zn25/cb*sb
      AAABC(259) = NMM34*Pa22*Sqrt2*Zn23*Zn34/xvev+NMM34
     &      *Pa22*Sqrt2*Zn24*Zn33/xvev+NMM55*Pa22*Sqrt2*Zn25
     &      *Zn35/xvev-Pa22*Sqrt2*Zn25*Zn35*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn21*Zn33/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn23*Zn31/cb*sb-EE/MW*NMM14*Pa21/SW*Zn21
     &      *Zn34/cb/sb+EE/MW*NMM14*Pa21/SW*Zn21*Zn34/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn24*Zn31/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn24*Zn31/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn22*Zn33/cb*sb-EE/MW*NMM23*Pa21/SW*Zn23
     &      *Zn32/cb*sb-EE/MW*NMM24*Pa21/SW*Zn22*Zn34/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn22*Zn34/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn24*Zn32/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn24*Zn32/cb*sb+EE/MW*NMM35*Pa21/SW*Zn23
     &      *Zn35/cb/sb-EE/MW*NMM35*Pa21/SW*Zn23*Zn35/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn25*Zn33/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn25*Zn33/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn24*Zn35/cb*sb+EE/MW*NMM45*Pa21/SW*Zn25
     &      *Zn34/cb*sb
      AAABC(260) = NMM34*Pa22*Sqrt2*Zn23*Zn54/xvev+NMM34
     &      *Pa22*Sqrt2*Zn24*Zn53/xvev+NMM55*Pa22*Sqrt2*Zn25
     &      *Zn55/xvev-Pa22*Sqrt2*Zn25*Zn55*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn21*Zn53/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn23*Zn51/cb*sb-EE/MW*NMM14*Pa21/SW*Zn21
     &      *Zn54/cb/sb+EE/MW*NMM14*Pa21/SW*Zn21*Zn54/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn24*Zn51/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn24*Zn51/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn22*Zn53/cb*sb-EE/MW*NMM23*Pa21/SW*Zn23
     &      *Zn52/cb*sb-EE/MW*NMM24*Pa21/SW*Zn22*Zn54/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn22*Zn54/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn24*Zn52/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn24*Zn52/cb*sb+EE/MW*NMM35*Pa21/SW*Zn23
     &      *Zn55/cb/sb-EE/MW*NMM35*Pa21/SW*Zn23*Zn55/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn25*Zn53/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn25*Zn53/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn24*Zn55/cb*sb+EE/MW*NMM45*Pa21/SW*Zn25
     &      *Zn54/cb*sb
      AAABC(261) = NMM34*Pa12*Sqrt2*Zn23*Zn44/xvev+NMM34
     &      *Pa12*Sqrt2*Zn24*Zn43/xvev+NMM55*Pa12*Sqrt2*Zn25
     &      *Zn45/xvev-Pa12*Sqrt2*Zn25*Zn45*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn21*Zn43/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn23*Zn41/cb*sb-EE/MW*NMM14*Pa11/SW*Zn21
     &      *Zn44/cb/sb+EE/MW*NMM14*Pa11/SW*Zn21*Zn44/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn24*Zn41/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn24*Zn41/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn22*Zn43/cb*sb-EE/MW*NMM23*Pa11/SW*Zn23
     &      *Zn42/cb*sb-EE/MW*NMM24*Pa11/SW*Zn22*Zn44/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn22*Zn44/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn24*Zn42/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn24*Zn42/cb*sb+EE/MW*NMM35*Pa11/SW*Zn23
     &      *Zn45/cb/sb-EE/MW*NMM35*Pa11/SW*Zn23*Zn45/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn25*Zn43/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn25*Zn43/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn24*Zn45/cb*sb+EE/MW*NMM45*Pa11/SW*Zn25
     &      *Zn44/cb*sb
      AAABC(262) = NMM34*Sqrt2*Zh33*Zn23*Zn44/xvev+NMM34
     &      *Sqrt2*Zh33*Zn24*Zn43/xvev+NMM55*Sqrt2*Zh33*Zn25
     &      *Zn45/xvev-Sqrt2*Zh33*Zn25*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn21*Zn43/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn23*Zn41/cb+EE/MW*NMM14/SW*Zh32*Zn21*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn24*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn22*Zn43/cb+EE/MW*NMM23/SW*Zh31*Zn23
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh32*Zn22*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn24*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn23*Zn45/sb+EE/MW*NMM35/SW*Zh32*Zn25*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn24*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn25*Zn44/cb
      AAABC(263) = NMM34*Sqrt2*Zh23*Zn23*Zn34/xvev+NMM34
     &      *Sqrt2*Zh23*Zn24*Zn33/xvev+NMM55*Sqrt2*Zh23*Zn25
     &      *Zn35/xvev-Sqrt2*Zh23*Zn25*Zn35*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn21*Zn33/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn23*Zn31/cb+EE/MW*NMM14/SW*Zh22*Zn21*Zn34/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn24*Zn31/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn22*Zn33/cb+EE/MW*NMM23/SW*Zh21*Zn23
     &      *Zn32/cb+EE/MW*NMM24/SW*Zh22*Zn22*Zn34/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn24*Zn32/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn23*Zn35/sb+EE/MW*NMM35/SW*Zh22*Zn25*Zn33
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn24*Zn35/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn25*Zn34/cb
      AAABC(264) = NMM34*Sqrt2*Zh33*Zn23*Zn54/xvev+NMM34
     &      *Sqrt2*Zh33*Zn24*Zn53/xvev+NMM55*Sqrt2*Zh33*Zn25
     &      *Zn55/xvev-Sqrt2*Zh33*Zn25*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn21*Zn53/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn23*Zn51/cb+EE/MW*NMM14/SW*Zh32*Zn21*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn24*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn22*Zn53/cb+EE/MW*NMM23/SW*Zh31*Zn23
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh32*Zn22*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn24*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn23*Zn55/sb+EE/MW*NMM35/SW*Zh32*Zn25*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn24*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn25*Zn54/cb
      AAABC(265) = NMM34*Sqrt2*Zh23*Zn23*Zn54/xvev+NMM34
     &      *Sqrt2*Zh23*Zn24*Zn53/xvev+NMM55*Sqrt2*Zh23*Zn25
     &      *Zn55/xvev-Sqrt2*Zh23*Zn25*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn21*Zn53/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn23*Zn51/cb+EE/MW*NMM14/SW*Zh22*Zn21*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn24*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn22*Zn53/cb+EE/MW*NMM23/SW*Zh21*Zn23
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh22*Zn22*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn24*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn23*Zn55/sb+EE/MW*NMM35/SW*Zh22*Zn25*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn24*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn25*Zn54/cb
      AAABC(266) = NMM34*Sqrt2*Zh13*Zn23*Zn44/xvev+NMM34
     &      *Sqrt2*Zh13*Zn24*Zn43/xvev+NMM55*Sqrt2*Zh13*Zn25
     &      *Zn45/xvev-Sqrt2*Zh13*Zn25*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn21*Zn43/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn23*Zn41/cb+EE/MW*NMM14/SW*Zh12*Zn21*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn24*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn22*Zn43/cb+EE/MW*NMM23/SW*Zh11*Zn23
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh12*Zn22*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn24*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn23*Zn45/sb+EE/MW*NMM35/SW*Zh12*Zn25*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn24*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn25*Zn44/cb
      AAABC(267) = NMM34*Sqrt2*Zh23*Zn23*Zn44/xvev+NMM34
     &      *Sqrt2*Zh23*Zn24*Zn43/xvev+NMM55*Sqrt2*Zh23*Zn25
     &      *Zn45/xvev-Sqrt2*Zh23*Zn25*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn21*Zn43/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn23*Zn41/cb+EE/MW*NMM14/SW*Zh22*Zn21*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn24*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn22*Zn43/cb+EE/MW*NMM23/SW*Zh21*Zn23
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh22*Zn22*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn24*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn23*Zn45/sb+EE/MW*NMM35/SW*Zh22*Zn25*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn24*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn25*Zn44/cb
      AAABC(268) = NMM34*Sqrt2*Zh13*Zn23*Zn34/xvev+NMM34
     &      *Sqrt2*Zh13*Zn24*Zn33/xvev+NMM55*Sqrt2*Zh13*Zn25
     &      *Zn35/xvev-Sqrt2*Zh13*Zn25*Zn35*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn21*Zn33/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn23*Zn31/cb+EE/MW*NMM14/SW*Zh12*Zn21*Zn34/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn24*Zn31/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn22*Zn33/cb+EE/MW*NMM23/SW*Zh11*Zn23
     &      *Zn32/cb+EE/MW*NMM24/SW*Zh12*Zn22*Zn34/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn24*Zn32/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn23*Zn35/sb+EE/MW*NMM35/SW*Zh12*Zn25*Zn33
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn24*Zn35/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn25*Zn34/cb
      AAABC(269) = NMM55*Sqrt2*Zh33*Zn25**2/xvev-Sqrt2*Zh33
     &      *Zn25**2*muP/xvev+2D0*NMM34*Sqrt2*Zh33*Zn23*Zn24
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh31*Zn21*Zn23/cb+2D0
     &      *EE/MW*NMM14/SW*Zh32*Zn21*Zn24/sb+2D0*EE/MW*NMM23
     &      /SW*Zh31*Zn22*Zn23/cb+2D0*EE/MW*NMM24/SW*Zh32
     &      *Zn22*Zn24/sb+2D0*EE/MW*NMM35/SW*Zh32*Zn23*Zn25
     &      /sb+2D0*EE/MW*NMM45/SW*Zh31*Zn24*Zn25/cb
      AAABC(270) = NMM55*Sqrt2*Zh23*Zn25**2/xvev-Sqrt2*Zh23
     &      *Zn25**2*muP/xvev+2D0*NMM34*Sqrt2*Zh23*Zn23*Zn24
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh21*Zn21*Zn23/cb+2D0
     &      *EE/MW*NMM14/SW*Zh22*Zn21*Zn24/sb+2D0*EE/MW*NMM23
     &      /SW*Zh21*Zn22*Zn23/cb+2D0*EE/MW*NMM24/SW*Zh22
     &      *Zn22*Zn24/sb+2D0*EE/MW*NMM35/SW*Zh22*Zn23*Zn25
     &      /sb+2D0*EE/MW*NMM45/SW*Zh21*Zn24*Zn25/cb
      AAABC(271) = NMM55*Sqrt2*Zh13*Zn25**2/xvev-Sqrt2*Zh13
     &      *Zn25**2*muP/xvev+2D0*NMM34*Sqrt2*Zh13*Zn23*Zn24
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh11*Zn21*Zn23/cb+2D0
     &      *EE/MW*NMM14/SW*Zh12*Zn21*Zn24/sb+2D0*EE/MW*NMM23
     &      /SW*Zh11*Zn22*Zn23/cb+2D0*EE/MW*NMM24/SW*Zh12
     &      *Zn22*Zn24/sb+2D0*EE/MW*NMM35/SW*Zh12*Zn23*Zn25
     &      /sb+2D0*EE/MW*NMM45/SW*Zh11*Zn24*Zn25/cb
      AAABC(272) = NMM34*Sqrt2*Zh33*Zn23*Zn34/xvev+NMM34
     &      *Sqrt2*Zh33*Zn24*Zn33/xvev+NMM55*Sqrt2*Zh33*Zn25
     &      *Zn35/xvev-Sqrt2*Zh33*Zn25*Zn35*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn21*Zn33/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn23*Zn31/cb+EE/MW*NMM14/SW*Zh32*Zn21*Zn34/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn24*Zn31/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn22*Zn33/cb+EE/MW*NMM23/SW*Zh31*Zn23
     &      *Zn32/cb+EE/MW*NMM24/SW*Zh32*Zn22*Zn34/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn24*Zn32/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn23*Zn35/sb+EE/MW*NMM35/SW*Zh32*Zn25*Zn33
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn24*Zn35/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn25*Zn34/cb
      AAABC(273) = NMM34*Sqrt2*Zh13*Zn23*Zn54/xvev+NMM34
     &      *Sqrt2*Zh13*Zn24*Zn53/xvev+NMM55*Sqrt2*Zh13*Zn25
     &      *Zn55/xvev-Sqrt2*Zh13*Zn25*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn21*Zn53/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn23*Zn51/cb+EE/MW*NMM14/SW*Zh12*Zn21*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn24*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn22*Zn53/cb+EE/MW*NMM23/SW*Zh11*Zn23
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh12*Zn22*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn24*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn23*Zn55/sb+EE/MW*NMM35/SW*Zh12*Zn25*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn24*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn25*Zn54/cb
      AAABC(274) = EE/MW/SW*(NMM13*Zn31*Zn33-NMM14*Zn31*Zn34
     &      +NMM23*Zn32*Zn33-NMM24*Zn32*Zn34+NMM35*Zn33*Zn35
     &      -NMM45*Zn34*Zn35)
      AAABC(275) = EE/MW/SW*(NMM13*Zn31*Zn43+NMM13*Zn33*Zn41
     &      -NMM14*Zn31*Zn44-NMM14*Zn34*Zn41+NMM23*Zn32*Zn43
     &      +NMM23*Zn33*Zn42-NMM24*Zn32*Zn44-NMM24*Zn34*Zn42
     &      +NMM35*Zn33*Zn45+NMM35*Zn35*Zn43-NMM45*Zn34*Zn45
     &      -NMM45*Zn35*Zn44)
      AAABC(276) = EE/MW/SW*(NMM13*Zn31*Zn53+NMM13*Zn33*Zn51
     &      -NMM14*Zn31*Zn54-NMM14*Zn34*Zn51+NMM23*Zn32*Zn53
     &      +NMM23*Zn33*Zn52-NMM24*Zn32*Zn54-NMM24*Zn34*Zn52
     &      +NMM35*Zn33*Zn55+NMM35*Zn35*Zn53-NMM45*Zn34*Zn55
     &      -NMM45*Zn35*Zn54)
      AAABC(277) = NMM34*Pa12*Sqrt2*Zn33*Zn54/xvev+NMM34
     &      *Pa12*Sqrt2*Zn34*Zn53/xvev+NMM55*Pa12*Sqrt2*Zn35
     &      *Zn55/xvev-Pa12*Sqrt2*Zn35*Zn55*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn31*Zn53/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn33*Zn51/cb*sb-EE/MW*NMM14*Pa11/SW*Zn31
     &      *Zn54/cb/sb+EE/MW*NMM14*Pa11/SW*Zn31*Zn54/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn34*Zn51/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn34*Zn51/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn32*Zn53/cb*sb-EE/MW*NMM23*Pa11/SW*Zn33
     &      *Zn52/cb*sb-EE/MW*NMM24*Pa11/SW*Zn32*Zn54/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn32*Zn54/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn34*Zn52/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn34*Zn52/cb*sb+EE/MW*NMM35*Pa11/SW*Zn33
     &      *Zn55/cb/sb-EE/MW*NMM35*Pa11/SW*Zn33*Zn55/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn35*Zn53/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn35*Zn53/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn34*Zn55/cb*sb+EE/MW*NMM45*Pa11/SW*Zn35
     &      *Zn54/cb*sb
      AAABC(278) = NMM55*Pa12*Sqrt2*Zn35**2/xvev-Pa12*Sqrt2
     &      *Zn35**2*muP/xvev+2D0*NMM34*Pa12*Sqrt2*Zn33*Zn34
     &      /xvev-2D0*EE/MW*NMM13*Pa11/SW*Zn31*Zn33/cb*sb
     &      -2D0*EE/MW*NMM14*Pa11/SW*Zn31*Zn34/cb/sb+2D0
     &      *EE/MW*NMM14*Pa11/SW*Zn31*Zn34/cb*sb-2D0*EE/MW
     &      *NMM23*Pa11/SW*Zn32*Zn33/cb*sb-2D0*EE/MW*NMM24
     &      *Pa11/SW*Zn32*Zn34/cb/sb+2D0*EE/MW*NMM24*Pa11
     &      /SW*Zn32*Zn34/cb*sb+2D0*EE/MW*NMM35*Pa11/SW*Zn33
     &      *Zn35/cb/sb-2D0*EE/MW*NMM35*Pa11/SW*Zn33*Zn35
     &      /cb*sb+2D0*EE/MW*NMM45*Pa11/SW*Zn34*Zn35/cb*sb
      AAABC(279) = NMM34*Pa12*Sqrt2*Zn33*Zn44/xvev+NMM34
     &      *Pa12*Sqrt2*Zn34*Zn43/xvev+NMM55*Pa12*Sqrt2*Zn35
     &      *Zn45/xvev-Pa12*Sqrt2*Zn35*Zn45*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn31*Zn43/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn33*Zn41/cb*sb-EE/MW*NMM14*Pa11/SW*Zn31
     &      *Zn44/cb/sb+EE/MW*NMM14*Pa11/SW*Zn31*Zn44/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn34*Zn41/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn34*Zn41/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn32*Zn43/cb*sb-EE/MW*NMM23*Pa11/SW*Zn33
     &      *Zn42/cb*sb-EE/MW*NMM24*Pa11/SW*Zn32*Zn44/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn32*Zn44/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn34*Zn42/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn34*Zn42/cb*sb+EE/MW*NMM35*Pa11/SW*Zn33
     &      *Zn45/cb/sb-EE/MW*NMM35*Pa11/SW*Zn33*Zn45/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn35*Zn43/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn35*Zn43/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn34*Zn45/cb*sb+EE/MW*NMM45*Pa11/SW*Zn35
     &      *Zn44/cb*sb
      AAABC(280) = NMM55*Pa22*Sqrt2*Zn35**2/xvev-Pa22*Sqrt2
     &      *Zn35**2*muP/xvev+2D0*NMM34*Pa22*Sqrt2*Zn33*Zn34
     &      /xvev-2D0*EE/MW*NMM13*Pa21/SW*Zn31*Zn33/cb*sb
     &      -2D0*EE/MW*NMM14*Pa21/SW*Zn31*Zn34/cb/sb+2D0
     &      *EE/MW*NMM14*Pa21/SW*Zn31*Zn34/cb*sb-2D0*EE/MW
     &      *NMM23*Pa21/SW*Zn32*Zn33/cb*sb-2D0*EE/MW*NMM24
     &      *Pa21/SW*Zn32*Zn34/cb/sb+2D0*EE/MW*NMM24*Pa21
     &      /SW*Zn32*Zn34/cb*sb+2D0*EE/MW*NMM35*Pa21/SW*Zn33
     &      *Zn35/cb/sb-2D0*EE/MW*NMM35*Pa21/SW*Zn33*Zn35
     &      /cb*sb+2D0*EE/MW*NMM45*Pa21/SW*Zn34*Zn35/cb*sb
      AAABC(281) = NMM34*Pa22*Sqrt2*Zn33*Zn54/xvev+NMM34
     &      *Pa22*Sqrt2*Zn34*Zn53/xvev+NMM55*Pa22*Sqrt2*Zn35
     &      *Zn55/xvev-Pa22*Sqrt2*Zn35*Zn55*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn31*Zn53/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn33*Zn51/cb*sb-EE/MW*NMM14*Pa21/SW*Zn31
     &      *Zn54/cb/sb+EE/MW*NMM14*Pa21/SW*Zn31*Zn54/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn34*Zn51/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn34*Zn51/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn32*Zn53/cb*sb-EE/MW*NMM23*Pa21/SW*Zn33
     &      *Zn52/cb*sb-EE/MW*NMM24*Pa21/SW*Zn32*Zn54/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn32*Zn54/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn34*Zn52/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn34*Zn52/cb*sb+EE/MW*NMM35*Pa21/SW*Zn33
     &      *Zn55/cb/sb-EE/MW*NMM35*Pa21/SW*Zn33*Zn55/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn35*Zn53/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn35*Zn53/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn34*Zn55/cb*sb+EE/MW*NMM45*Pa21/SW*Zn35
     &      *Zn54/cb*sb
      AAABC(282) = NMM34*Pa22*Sqrt2*Zn33*Zn44/xvev+NMM34
     &      *Pa22*Sqrt2*Zn34*Zn43/xvev+NMM55*Pa22*Sqrt2*Zn35
     &      *Zn45/xvev-Pa22*Sqrt2*Zn35*Zn45*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn31*Zn43/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn33*Zn41/cb*sb-EE/MW*NMM14*Pa21/SW*Zn31
     &      *Zn44/cb/sb+EE/MW*NMM14*Pa21/SW*Zn31*Zn44/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn34*Zn41/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn34*Zn41/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn32*Zn43/cb*sb-EE/MW*NMM23*Pa21/SW*Zn33
     &      *Zn42/cb*sb-EE/MW*NMM24*Pa21/SW*Zn32*Zn44/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn32*Zn44/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn34*Zn42/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn34*Zn42/cb*sb+EE/MW*NMM35*Pa21/SW*Zn33
     &      *Zn45/cb/sb-EE/MW*NMM35*Pa21/SW*Zn33*Zn45/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn35*Zn43/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn35*Zn43/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn34*Zn45/cb*sb+EE/MW*NMM45*Pa21/SW*Zn35
     &      *Zn44/cb*sb
      AAABC(283) = NMM55*Sqrt2*Zh23*Zn35**2/xvev-Sqrt2*Zh23
     &      *Zn35**2*muP/xvev+2D0*NMM34*Sqrt2*Zh23*Zn33*Zn34
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh21*Zn31*Zn33/cb+2D0
     &      *EE/MW*NMM14/SW*Zh22*Zn31*Zn34/sb+2D0*EE/MW*NMM23
     &      /SW*Zh21*Zn32*Zn33/cb+2D0*EE/MW*NMM24/SW*Zh22
     &      *Zn32*Zn34/sb+2D0*EE/MW*NMM35/SW*Zh22*Zn33*Zn35
     &      /sb+2D0*EE/MW*NMM45/SW*Zh21*Zn34*Zn35/cb
      AAABC(284) = NMM34*Sqrt2*Zh23*Zn33*Zn54/xvev+NMM34
     &      *Sqrt2*Zh23*Zn34*Zn53/xvev+NMM55*Sqrt2*Zh23*Zn35
     &      *Zn55/xvev-Sqrt2*Zh23*Zn35*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn31*Zn53/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn33*Zn51/cb+EE/MW*NMM14/SW*Zh22*Zn31*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn34*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn32*Zn53/cb+EE/MW*NMM23/SW*Zh21*Zn33
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh22*Zn32*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn34*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn33*Zn55/sb+EE/MW*NMM35/SW*Zh22*Zn35*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn34*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn35*Zn54/cb
      AAABC(285) = NMM55*Sqrt2*Zh13*Zn35**2/xvev-Sqrt2*Zh13
     &      *Zn35**2*muP/xvev+2D0*NMM34*Sqrt2*Zh13*Zn33*Zn34
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh11*Zn31*Zn33/cb+2D0
     &      *EE/MW*NMM14/SW*Zh12*Zn31*Zn34/sb+2D0*EE/MW*NMM23
     &      /SW*Zh11*Zn32*Zn33/cb+2D0*EE/MW*NMM24/SW*Zh12
     &      *Zn32*Zn34/sb+2D0*EE/MW*NMM35/SW*Zh12*Zn33*Zn35
     &      /sb+2D0*EE/MW*NMM45/SW*Zh11*Zn34*Zn35/cb
      AAABC(286) = NMM34*Sqrt2*Zh13*Zn33*Zn54/xvev+NMM34
     &      *Sqrt2*Zh13*Zn34*Zn53/xvev+NMM55*Sqrt2*Zh13*Zn35
     &      *Zn55/xvev-Sqrt2*Zh13*Zn35*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn31*Zn53/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn33*Zn51/cb+EE/MW*NMM14/SW*Zh12*Zn31*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn34*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn32*Zn53/cb+EE/MW*NMM23/SW*Zh11*Zn33
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh12*Zn32*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn34*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn33*Zn55/sb+EE/MW*NMM35/SW*Zh12*Zn35*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn34*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn35*Zn54/cb
      AAABC(287) = NMM34*Sqrt2*Zh33*Zn33*Zn44/xvev+NMM34
     &      *Sqrt2*Zh33*Zn34*Zn43/xvev+NMM55*Sqrt2*Zh33*Zn35
     &      *Zn45/xvev-Sqrt2*Zh33*Zn35*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn31*Zn43/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn33*Zn41/cb+EE/MW*NMM14/SW*Zh32*Zn31*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn34*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn32*Zn43/cb+EE/MW*NMM23/SW*Zh31*Zn33
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh32*Zn32*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn34*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn33*Zn45/sb+EE/MW*NMM35/SW*Zh32*Zn35*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn34*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn35*Zn44/cb
      AAABC(288) = NMM34*Sqrt2*Zh33*Zn33*Zn54/xvev+NMM34
     &      *Sqrt2*Zh33*Zn34*Zn53/xvev+NMM55*Sqrt2*Zh33*Zn35
     &      *Zn55/xvev-Sqrt2*Zh33*Zn35*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn31*Zn53/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn33*Zn51/cb+EE/MW*NMM14/SW*Zh32*Zn31*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn34*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn32*Zn53/cb+EE/MW*NMM23/SW*Zh31*Zn33
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh32*Zn32*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn34*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn33*Zn55/sb+EE/MW*NMM35/SW*Zh32*Zn35*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn34*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn35*Zn54/cb
      AAABC(289) = NMM34*Sqrt2*Zh13*Zn33*Zn44/xvev+NMM34
     &      *Sqrt2*Zh13*Zn34*Zn43/xvev+NMM55*Sqrt2*Zh13*Zn35
     &      *Zn45/xvev-Sqrt2*Zh13*Zn35*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn31*Zn43/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn33*Zn41/cb+EE/MW*NMM14/SW*Zh12*Zn31*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn34*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn32*Zn43/cb+EE/MW*NMM23/SW*Zh11*Zn33
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh12*Zn32*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn34*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn33*Zn45/sb+EE/MW*NMM35/SW*Zh12*Zn35*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn34*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn35*Zn44/cb
      AAABC(290) = NMM55*Sqrt2*Zh33*Zn35**2/xvev-Sqrt2*Zh33
     &      *Zn35**2*muP/xvev+2D0*NMM34*Sqrt2*Zh33*Zn33*Zn34
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh31*Zn31*Zn33/cb+2D0
     &      *EE/MW*NMM14/SW*Zh32*Zn31*Zn34/sb+2D0*EE/MW*NMM23
     &      /SW*Zh31*Zn32*Zn33/cb+2D0*EE/MW*NMM24/SW*Zh32
     &      *Zn32*Zn34/sb+2D0*EE/MW*NMM35/SW*Zh32*Zn33*Zn35
     &      /sb+2D0*EE/MW*NMM45/SW*Zh31*Zn34*Zn35/cb
      AAABC(291) = NMM34*Sqrt2*Zh23*Zn33*Zn44/xvev+NMM34
     &      *Sqrt2*Zh23*Zn34*Zn43/xvev+NMM55*Sqrt2*Zh23*Zn35
     &      *Zn45/xvev-Sqrt2*Zh23*Zn35*Zn45*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn31*Zn43/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn33*Zn41/cb+EE/MW*NMM14/SW*Zh22*Zn31*Zn44/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn34*Zn41/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn32*Zn43/cb+EE/MW*NMM23/SW*Zh21*Zn33
     &      *Zn42/cb+EE/MW*NMM24/SW*Zh22*Zn32*Zn44/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn34*Zn42/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn33*Zn45/sb+EE/MW*NMM35/SW*Zh22*Zn35*Zn43
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn34*Zn45/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn35*Zn44/cb
      AAABC(292) = EE/MW/SW*(NMM13*Zn41*Zn53+NMM13*Zn43*Zn51
     &      -NMM14*Zn41*Zn54-NMM14*Zn44*Zn51+NMM23*Zn42*Zn53
     &      +NMM23*Zn43*Zn52-NMM24*Zn42*Zn54-NMM24*Zn44*Zn52
     &      +NMM35*Zn43*Zn55+NMM35*Zn45*Zn53-NMM45*Zn44*Zn55
     &      -NMM45*Zn45*Zn54)
      AAABC(293) = EE/MW/SW*(NMM13*Zn41*Zn43-NMM14*Zn41*Zn44
     &      +NMM23*Zn42*Zn43-NMM24*Zn42*Zn44+NMM35*Zn43*Zn45
     &      -NMM45*Zn44*Zn45)
      AAABC(294) = NMM34*Pa22*Sqrt2*Zn43*Zn54/xvev+NMM34
     &      *Pa22*Sqrt2*Zn44*Zn53/xvev+NMM55*Pa22*Sqrt2*Zn45
     &      *Zn55/xvev-Pa22*Sqrt2*Zn45*Zn55*muP/xvev-EE/MW
     &      *NMM13*Pa21/SW*Zn41*Zn53/cb*sb-EE/MW*NMM13*Pa21
     &      /SW*Zn43*Zn51/cb*sb-EE/MW*NMM14*Pa21/SW*Zn41
     &      *Zn54/cb/sb+EE/MW*NMM14*Pa21/SW*Zn41*Zn54/cb
     &      *sb-EE/MW*NMM14*Pa21/SW*Zn44*Zn51/cb/sb+EE/MW
     &      *NMM14*Pa21/SW*Zn44*Zn51/cb*sb-EE/MW*NMM23*Pa21
     &      /SW*Zn42*Zn53/cb*sb-EE/MW*NMM23*Pa21/SW*Zn43
     &      *Zn52/cb*sb-EE/MW*NMM24*Pa21/SW*Zn42*Zn54/cb
     &      /sb+EE/MW*NMM24*Pa21/SW*Zn42*Zn54/cb*sb-EE/MW
     &      *NMM24*Pa21/SW*Zn44*Zn52/cb/sb+EE/MW*NMM24*Pa21
     &      /SW*Zn44*Zn52/cb*sb+EE/MW*NMM35*Pa21/SW*Zn43
     &      *Zn55/cb/sb-EE/MW*NMM35*Pa21/SW*Zn43*Zn55/cb
     &      *sb+EE/MW*NMM35*Pa21/SW*Zn45*Zn53/cb/sb-EE/MW
     &      *NMM35*Pa21/SW*Zn45*Zn53/cb*sb+EE/MW*NMM45*Pa21
     &      /SW*Zn44*Zn55/cb*sb+EE/MW*NMM45*Pa21/SW*Zn45
     &      *Zn54/cb*sb
      AAABC(295) = NMM34*Pa12*Sqrt2*Zn43*Zn54/xvev+NMM34
     &      *Pa12*Sqrt2*Zn44*Zn53/xvev+NMM55*Pa12*Sqrt2*Zn45
     &      *Zn55/xvev-Pa12*Sqrt2*Zn45*Zn55*muP/xvev-EE/MW
     &      *NMM13*Pa11/SW*Zn41*Zn53/cb*sb-EE/MW*NMM13*Pa11
     &      /SW*Zn43*Zn51/cb*sb-EE/MW*NMM14*Pa11/SW*Zn41
     &      *Zn54/cb/sb+EE/MW*NMM14*Pa11/SW*Zn41*Zn54/cb
     &      *sb-EE/MW*NMM14*Pa11/SW*Zn44*Zn51/cb/sb+EE/MW
     &      *NMM14*Pa11/SW*Zn44*Zn51/cb*sb-EE/MW*NMM23*Pa11
     &      /SW*Zn42*Zn53/cb*sb-EE/MW*NMM23*Pa11/SW*Zn43
     &      *Zn52/cb*sb-EE/MW*NMM24*Pa11/SW*Zn42*Zn54/cb
     &      /sb+EE/MW*NMM24*Pa11/SW*Zn42*Zn54/cb*sb-EE/MW
     &      *NMM24*Pa11/SW*Zn44*Zn52/cb/sb+EE/MW*NMM24*Pa11
     &      /SW*Zn44*Zn52/cb*sb+EE/MW*NMM35*Pa11/SW*Zn43
     &      *Zn55/cb/sb-EE/MW*NMM35*Pa11/SW*Zn43*Zn55/cb
     &      *sb+EE/MW*NMM35*Pa11/SW*Zn45*Zn53/cb/sb-EE/MW
     &      *NMM35*Pa11/SW*Zn45*Zn53/cb*sb+EE/MW*NMM45*Pa11
     &      /SW*Zn44*Zn55/cb*sb+EE/MW*NMM45*Pa11/SW*Zn45
     &      *Zn54/cb*sb
      AAABC(296) = NMM55*Pa22*Sqrt2*Zn45**2/xvev-Pa22*Sqrt2
     &      *Zn45**2*muP/xvev+2D0*NMM34*Pa22*Sqrt2*Zn43*Zn44
     &      /xvev-2D0*EE/MW*NMM13*Pa21/SW*Zn41*Zn43/cb*sb
     &      -2D0*EE/MW*NMM14*Pa21/SW*Zn41*Zn44/cb/sb+2D0
     &      *EE/MW*NMM14*Pa21/SW*Zn41*Zn44/cb*sb-2D0*EE/MW
     &      *NMM23*Pa21/SW*Zn42*Zn43/cb*sb-2D0*EE/MW*NMM24
     &      *Pa21/SW*Zn42*Zn44/cb/sb+2D0*EE/MW*NMM24*Pa21
     &      /SW*Zn42*Zn44/cb*sb+2D0*EE/MW*NMM35*Pa21/SW*Zn43
     &      *Zn45/cb/sb-2D0*EE/MW*NMM35*Pa21/SW*Zn43*Zn45
     &      /cb*sb+2D0*EE/MW*NMM45*Pa21/SW*Zn44*Zn45/cb*sb
      AAABC(297) = NMM55*Pa12*Sqrt2*Zn45**2/xvev-Pa12*Sqrt2
     &      *Zn45**2*muP/xvev+2D0*NMM34*Pa12*Sqrt2*Zn43*Zn44
     &      /xvev-2D0*EE/MW*NMM13*Pa11/SW*Zn41*Zn43/cb*sb
     &      -2D0*EE/MW*NMM14*Pa11/SW*Zn41*Zn44/cb/sb+2D0
     &      *EE/MW*NMM14*Pa11/SW*Zn41*Zn44/cb*sb-2D0*EE/MW
     &      *NMM23*Pa11/SW*Zn42*Zn43/cb*sb-2D0*EE/MW*NMM24
     &      *Pa11/SW*Zn42*Zn44/cb/sb+2D0*EE/MW*NMM24*Pa11
     &      /SW*Zn42*Zn44/cb*sb+2D0*EE/MW*NMM35*Pa11/SW*Zn43
     &      *Zn45/cb/sb-2D0*EE/MW*NMM35*Pa11/SW*Zn43*Zn45
     &      /cb*sb+2D0*EE/MW*NMM45*Pa11/SW*Zn44*Zn45/cb*sb
      AAABC(298) = NMM55*Sqrt2*Zh23*Zn45**2/xvev-Sqrt2*Zh23
     &      *Zn45**2*muP/xvev+2D0*NMM34*Sqrt2*Zh23*Zn43*Zn44
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh21*Zn41*Zn43/cb+2D0
     &      *EE/MW*NMM14/SW*Zh22*Zn41*Zn44/sb+2D0*EE/MW*NMM23
     &      /SW*Zh21*Zn42*Zn43/cb+2D0*EE/MW*NMM24/SW*Zh22
     &      *Zn42*Zn44/sb+2D0*EE/MW*NMM35/SW*Zh22*Zn43*Zn45
     &      /sb+2D0*EE/MW*NMM45/SW*Zh21*Zn44*Zn45/cb
      AAABC(299) = NMM34*Sqrt2*Zh33*Zn43*Zn54/xvev+NMM34
     &      *Sqrt2*Zh33*Zn44*Zn53/xvev+NMM55*Sqrt2*Zh33*Zn45
     &      *Zn55/xvev-Sqrt2*Zh33*Zn45*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh31*Zn41*Zn53/cb+EE/MW*NMM13/SW*Zh31
     &      *Zn43*Zn51/cb+EE/MW*NMM14/SW*Zh32*Zn41*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh32*Zn44*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh31*Zn42*Zn53/cb+EE/MW*NMM23/SW*Zh31*Zn43
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh32*Zn42*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh32*Zn44*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh32*Zn43*Zn55/sb+EE/MW*NMM35/SW*Zh32*Zn45*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh31*Zn44*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh31*Zn45*Zn54/cb
      AAABC(300) = NMM34*Sqrt2*Zh23*Zn43*Zn54/xvev+NMM34
     &      *Sqrt2*Zh23*Zn44*Zn53/xvev+NMM55*Sqrt2*Zh23*Zn45
     &      *Zn55/xvev-Sqrt2*Zh23*Zn45*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh21*Zn41*Zn53/cb+EE/MW*NMM13/SW*Zh21
     &      *Zn43*Zn51/cb+EE/MW*NMM14/SW*Zh22*Zn41*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh22*Zn44*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh21*Zn42*Zn53/cb+EE/MW*NMM23/SW*Zh21*Zn43
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh22*Zn42*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh22*Zn44*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh22*Zn43*Zn55/sb+EE/MW*NMM35/SW*Zh22*Zn45*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh21*Zn44*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh21*Zn45*Zn54/cb
      AAABC(301) = NMM34*Sqrt2*Zh13*Zn43*Zn54/xvev+NMM34
     &      *Sqrt2*Zh13*Zn44*Zn53/xvev+NMM55*Sqrt2*Zh13*Zn45
     &      *Zn55/xvev-Sqrt2*Zh13*Zn45*Zn55*muP/xvev+EE/MW
     &      *NMM13/SW*Zh11*Zn41*Zn53/cb+EE/MW*NMM13/SW*Zh11
     &      *Zn43*Zn51/cb+EE/MW*NMM14/SW*Zh12*Zn41*Zn54/sb
     &      +EE/MW*NMM14/SW*Zh12*Zn44*Zn51/sb+EE/MW*NMM23
     &      /SW*Zh11*Zn42*Zn53/cb+EE/MW*NMM23/SW*Zh11*Zn43
     &      *Zn52/cb+EE/MW*NMM24/SW*Zh12*Zn42*Zn54/sb+EE
     &      /MW*NMM24/SW*Zh12*Zn44*Zn52/sb+EE/MW*NMM35/SW
     &      *Zh12*Zn43*Zn55/sb+EE/MW*NMM35/SW*Zh12*Zn45*Zn53
     &      /sb+EE/MW*NMM45/SW*Zh11*Zn44*Zn55/cb+EE/MW*NMM45
     &      /SW*Zh11*Zn45*Zn54/cb
      AAABC(302) = NMM55*Sqrt2*Zh13*Zn45**2/xvev-Sqrt2*Zh13
     &      *Zn45**2*muP/xvev+2D0*NMM34*Sqrt2*Zh13*Zn43*Zn44
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh11*Zn41*Zn43/cb+2D0
     &      *EE/MW*NMM14/SW*Zh12*Zn41*Zn44/sb+2D0*EE/MW*NMM23
     &      /SW*Zh11*Zn42*Zn43/cb+2D0*EE/MW*NMM24/SW*Zh12
     &      *Zn42*Zn44/sb+2D0*EE/MW*NMM35/SW*Zh12*Zn43*Zn45
     &      /sb+2D0*EE/MW*NMM45/SW*Zh11*Zn44*Zn45/cb
      AAABC(303) = NMM55*Sqrt2*Zh33*Zn45**2/xvev-Sqrt2*Zh33
     &      *Zn45**2*muP/xvev+2D0*NMM34*Sqrt2*Zh33*Zn43*Zn44
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh31*Zn41*Zn43/cb+2D0
     &      *EE/MW*NMM14/SW*Zh32*Zn41*Zn44/sb+2D0*EE/MW*NMM23
     &      /SW*Zh31*Zn42*Zn43/cb+2D0*EE/MW*NMM24/SW*Zh32
     &      *Zn42*Zn44/sb+2D0*EE/MW*NMM35/SW*Zh32*Zn43*Zn45
     &      /sb+2D0*EE/MW*NMM45/SW*Zh31*Zn44*Zn45/cb
      AAABC(304) = EE/MW/SW*(NMM13*Zn51*Zn53-NMM14*Zn51*Zn54
     &      +NMM23*Zn52*Zn53-NMM24*Zn52*Zn54+NMM35*Zn53*Zn55
     &      -NMM45*Zn54*Zn55)
      AAABC(305) = NMM55*Pa12*Sqrt2*Zn55**2/xvev-Pa12*Sqrt2
     &      *Zn55**2*muP/xvev+2D0*NMM34*Pa12*Sqrt2*Zn53*Zn54
     &      /xvev-2D0*EE/MW*NMM13*Pa11/SW*Zn51*Zn53/cb*sb
     &      -2D0*EE/MW*NMM14*Pa11/SW*Zn51*Zn54/cb/sb+2D0
     &      *EE/MW*NMM14*Pa11/SW*Zn51*Zn54/cb*sb-2D0*EE/MW
     &      *NMM23*Pa11/SW*Zn52*Zn53/cb*sb-2D0*EE/MW*NMM24
     &      *Pa11/SW*Zn52*Zn54/cb/sb+2D0*EE/MW*NMM24*Pa11
     &      /SW*Zn52*Zn54/cb*sb+2D0*EE/MW*NMM35*Pa11/SW*Zn53
     &      *Zn55/cb/sb-2D0*EE/MW*NMM35*Pa11/SW*Zn53*Zn55
     &      /cb*sb+2D0*EE/MW*NMM45*Pa11/SW*Zn54*Zn55/cb*sb
      AAABC(306) = NMM55*Pa22*Sqrt2*Zn55**2/xvev-Pa22*Sqrt2
     &      *Zn55**2*muP/xvev+2D0*NMM34*Pa22*Sqrt2*Zn53*Zn54
     &      /xvev-2D0*EE/MW*NMM13*Pa21/SW*Zn51*Zn53/cb*sb
     &      -2D0*EE/MW*NMM14*Pa21/SW*Zn51*Zn54/cb/sb+2D0
     &      *EE/MW*NMM14*Pa21/SW*Zn51*Zn54/cb*sb-2D0*EE/MW
     &      *NMM23*Pa21/SW*Zn52*Zn53/cb*sb-2D0*EE/MW*NMM24
     &      *Pa21/SW*Zn52*Zn54/cb/sb+2D0*EE/MW*NMM24*Pa21
     &      /SW*Zn52*Zn54/cb*sb+2D0*EE/MW*NMM35*Pa21/SW*Zn53
     &      *Zn55/cb/sb-2D0*EE/MW*NMM35*Pa21/SW*Zn53*Zn55
     &      /cb*sb+2D0*EE/MW*NMM45*Pa21/SW*Zn54*Zn55/cb*sb
      AAABC(307) = NMM55*Sqrt2*Zh33*Zn55**2/xvev-Sqrt2*Zh33
     &      *Zn55**2*muP/xvev+2D0*NMM34*Sqrt2*Zh33*Zn53*Zn54
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh31*Zn51*Zn53/cb+2D0
     &      *EE/MW*NMM14/SW*Zh32*Zn51*Zn54/sb+2D0*EE/MW*NMM23
     &      /SW*Zh31*Zn52*Zn53/cb+2D0*EE/MW*NMM24/SW*Zh32
     &      *Zn52*Zn54/sb+2D0*EE/MW*NMM35/SW*Zh32*Zn53*Zn55
     &      /sb+2D0*EE/MW*NMM45/SW*Zh31*Zn54*Zn55/cb
      AAABC(308) = NMM55*Sqrt2*Zh13*Zn55**2/xvev-Sqrt2*Zh13
     &      *Zn55**2*muP/xvev+2D0*NMM34*Sqrt2*Zh13*Zn53*Zn54
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh11*Zn51*Zn53/cb+2D0
     &      *EE/MW*NMM14/SW*Zh12*Zn51*Zn54/sb+2D0*EE/MW*NMM23
     &      /SW*Zh11*Zn52*Zn53/cb+2D0*EE/MW*NMM24/SW*Zh12
     &      *Zn52*Zn54/sb+2D0*EE/MW*NMM35/SW*Zh12*Zn53*Zn55
     &      /sb+2D0*EE/MW*NMM45/SW*Zh11*Zn54*Zn55/cb
      AAABC(309) = NMM55*Sqrt2*Zh23*Zn55**2/xvev-Sqrt2*Zh23
     &      *Zn55**2*muP/xvev+2D0*NMM34*Sqrt2*Zh23*Zn53*Zn54
     &      /xvev+2D0*EE/MW*NMM13/SW*Zh21*Zn51*Zn53/cb+2D0
     &      *EE/MW*NMM14/SW*Zh22*Zn51*Zn54/sb+2D0*EE/MW*NMM23
     &      /SW*Zh21*Zn52*Zn53/cb+2D0*EE/MW*NMM24/SW*Zh22
     &      *Zn52*Zn54/sb+2D0*EE/MW*NMM35/SW*Zh22*Zn53*Zn55
     &      /sb+2D0*EE/MW*NMM45/SW*Zh21*Zn54*Zn55/cb
      AAABR(364) = EE*(SW/CW*Zu12*Zu22-2D0*CW/SW*Zu11*Zu21
     &      -CW/SW*Zu12*Zu22)
      AAABR(365) = EE*(SW/CW*Zv12*Zv22-2D0*CW/SW*Zv11*Zv21
     &      -CW/SW*Zv12*Zv22)
      AAABR(366) = EE*(SW/CW*Zu12**2-2D0*CW/SW*Zu11**2-CW
     &      /SW*Zu12**2)
      AAABR(367) = EE*(SW/CW*Zv12**2-2D0*CW/SW*Zv11**2-CW
     &      /SW*Zv12**2)
      AAABR(368) = EE*(SW/CW*Zu22**2-2D0*CW/SW*Zu21**2-CW
     &      /SW*Zu22**2)
      AAABR(369) = EE*(SW/CW*Zv22**2-2D0*CW/SW*Zv21**2-CW
     &      /SW*Zv22**2)
      AAABC(310) = EE/SW*(2D0*dconjg(Zn42)*Zu21+Sqrt2*dconjg(Zn43
     &      )*Zu22)
      AAABC(311) = EE/SW*(2D0*Zn42*Zv21-Sqrt2*Zn44*Zv22)
      AAABC(312) = EE/SW*(2D0*dconjg(Zn52)*Zu11+Sqrt2*dconjg(Zn53
     &      )*Zu12)
      AAABC(313) = EE/SW*(2D0*Zn52*Zv11-Sqrt2*Zn54*Zv12)
      AAABC(314) = EE/SW*(2D0*dconjg(Zn32)*Zu21+Sqrt2*dconjg(Zn33
     &      )*Zu22)
      AAABC(315) = EE/SW*(2D0*Zn32*Zv21-Sqrt2*Zn34*Zv22)
      AAABC(316) = EE/SW*(2D0*dconjg(Zn12)*Zu21+Sqrt2*dconjg(Zn13
     &      )*Zu22)
      AAABC(317) = EE/SW*(2D0*Zn12*Zv21-Sqrt2*Zn14*Zv22)
      AAABC(318) = EE/SW*(2D0*dconjg(Zn22)*Zu11+Sqrt2*dconjg(Zn23
     &      )*Zu12)
      AAABC(319) = EE/SW*(2D0*Zn22*Zv11-Sqrt2*Zn24*Zv12)
      AAABC(320) = EE/SW*(2D0*dconjg(Zn22)*Zu21+Sqrt2*dconjg(Zn23
     &      )*Zu22)
      AAABC(321) = EE/SW*(2D0*Zn22*Zv21-Sqrt2*Zn24*Zv22)
      AAABC(322) = EE/SW*(2D0*dconjg(Zn42)*Zu11+Sqrt2*dconjg(Zn43
     &      )*Zu12)
      AAABC(323) = EE/SW*(2D0*Zn42*Zv11-Sqrt2*Zn44*Zv12)
      AAABC(324) = EE/SW*(2D0*dconjg(Zn52)*Zu21+Sqrt2*dconjg(Zn53
     &      )*Zu22)
      AAABC(325) = EE/SW*(2D0*Zn52*Zv21-Sqrt2*Zn54*Zv22)
      AAABC(326) = EE/SW*(2D0*dconjg(Zn32)*Zu11+Sqrt2*dconjg(Zn33
     &      )*Zu12)
      AAABC(327) = EE/SW*(2D0*Zn32*Zv11-Sqrt2*Zn34*Zv12)
      AAABC(328) = EE/SW*(2D0*dconjg(Zn12)*Zu11+Sqrt2*dconjg(Zn13
     &      )*Zu12)
      AAABC(329) = EE/SW*(2D0*Zn12*Zv11-Sqrt2*Zn14*Zv12)
      AAABR(370) = EE*(SW/CW+3D0*CW/SW)
      AAABR(371) = EE*(SW/CW-3D0*CW/SW)
      AAABC(330) = EE*(SW/CW*dconjg(Zn13)*Zn23-SW/CW*dconjg(Zn14
     &      )*Zn24+CW/SW*dconjg(Zn13)*Zn23-CW/SW*dconjg(Zn14
     &      )*Zn24)
      AAABC(331) = EE*(SW/CW*dconjg(Zn13)*Zn33-SW/CW*dconjg(Zn14
     &      )*Zn34+CW/SW*dconjg(Zn13)*Zn33-CW/SW*dconjg(Zn14
     &      )*Zn34)
      AAABC(332) = EE*(SW/CW*dconjg(Zn13)*Zn53-SW/CW*dconjg(Zn14
     &      )*Zn54+CW/SW*dconjg(Zn13)*Zn53-CW/SW*dconjg(Zn14
     &      )*Zn54)
      AAABC(333) = EE*(SW/CW*dconjg(Zn13)*Zn43-SW/CW*dconjg(Zn14
     &      )*Zn44+CW/SW*dconjg(Zn13)*Zn43-CW/SW*dconjg(Zn14
     &      )*Zn44)
      AAABC(334) = EE*(SW/CW*Zn13*dconjg(Zn13)-SW/CW*Zn14
     &      *dconjg(Zn14)+CW/SW*Zn13*dconjg(Zn13)-CW/SW*Zn14
     &      *dconjg(Zn14))
      AAABC(335) = EE*(SW/CW*dconjg(Zn23)*Zn53-SW/CW*dconjg(Zn24
     &      )*Zn54+CW/SW*dconjg(Zn23)*Zn53-CW/SW*dconjg(Zn24
     &      )*Zn54)
      AAABC(336) = EE*(SW/CW*Zn23*dconjg(Zn23)-SW/CW*Zn24
     &      *dconjg(Zn24)+CW/SW*Zn23*dconjg(Zn23)-CW/SW*Zn24
     &      *dconjg(Zn24))
      AAABC(337) = EE*(SW/CW*dconjg(Zn23)*Zn43-SW/CW*dconjg(Zn24
     &      )*Zn44+CW/SW*dconjg(Zn23)*Zn43-CW/SW*dconjg(Zn24
     &      )*Zn44)
      AAABC(338) = EE*(SW/CW*dconjg(Zn23)*Zn33-SW/CW*dconjg(Zn24
     &      )*Zn34+CW/SW*dconjg(Zn23)*Zn33-CW/SW*dconjg(Zn24
     &      )*Zn34)
      AAABC(339) = EE*(SW/CW*dconjg(Zn33)*Zn43-SW/CW*dconjg(Zn34
     &      )*Zn44+CW/SW*dconjg(Zn33)*Zn43-CW/SW*dconjg(Zn34
     &      )*Zn44)
      AAABC(340) = EE*(SW/CW*dconjg(Zn33)*Zn53-SW/CW*dconjg(Zn34
     &      )*Zn54+CW/SW*dconjg(Zn33)*Zn53-CW/SW*dconjg(Zn34
     &      )*Zn54)
      AAABC(341) = EE*(SW/CW*Zn33*dconjg(Zn33)-SW/CW*Zn34
     &      *dconjg(Zn34)+CW/SW*Zn33*dconjg(Zn33)-CW/SW*Zn34
     &      *dconjg(Zn34))
      AAABC(342) = EE*(SW/CW*Zn43*dconjg(Zn43)-SW/CW*Zn44
     &      *dconjg(Zn44)+CW/SW*Zn43*dconjg(Zn43)-CW/SW*Zn44
     &      *dconjg(Zn44))
      AAABC(343) = EE*(SW/CW*dconjg(Zn43)*Zn53-SW/CW*dconjg(Zn44
     &      )*Zn54+CW/SW*dconjg(Zn43)*Zn53-CW/SW*dconjg(Zn44
     &      )*Zn54)
      AAABC(344) = EE*(SW/CW*Zn53*dconjg(Zn53)-SW/CW*Zn54
     &      *dconjg(Zn54)+CW/SW*Zn53*dconjg(Zn53)-CW/SW*Zn54
     &      *dconjg(Zn54))
      AAABR(372) = EE**2*(SW/CW-CW/SW)
      AAABR(373) = EE**2
      AAABR(374) = EE**2*(2D0-SW**2/CW**2-CW**2/SW**2)
      AAABR(375) = EE**2/SW**2
      AAABR(376) = EE**2/CW*Pa21
      AAABR(377) = EE**2*Pa21/SW
      AAABR(378) = EE**2/CW*Pa11
      AAABR(379) = EE**2*Pa11/SW
      AAABR(380) = EE**2/CW*(Zh21*sb-Zh22*cb)
      AAABR(381) = EE**2/SW*(Zh21*sb-Zh22*cb)
      AAABR(382) = EE**2/CW*(Zh11*sb-Zh12*cb)
      AAABR(383) = EE**2/SW*(Zh11*sb-Zh12*cb)
      AAABR(384) = EE**2/CW*(Zh31*sb-Zh32*cb)
      AAABR(385) = EE**2/SW*(Zh31*sb-Zh32*cb)
      AAABR(386) = EE**2*(SW/CW*Zl11**2+2D0*SW/CW*Zl12**2
     &      -CW/SW*Zl11**2)
      AAABR(387) = EE**2*(2D0*Zl11**2-SW**2/CW**2*Zl11**2
     &      -4D0*SW**2/CW**2*Zl12**2-CW**2/SW**2*Zl11**2)
      AAABR(388) = EE**2/SW**2*Zl11**2
      AAABR(389) = EE**2*(SW/CW*Zl11*Zl21+2D0*SW/CW*Zl12
     &      *Zl22-CW/SW*Zl11*Zl21)
      AAABR(390) = EE**2*(2D0*Zl11*Zl21-SW**2/CW**2*Zl11
     &      *Zl21-4D0*SW**2/CW**2*Zl12*Zl22-CW**2/SW**2*Zl11
     &      *Zl21)
      AAABR(391) = EE**2/SW**2*Zl11*Zl21
      AAABR(392) = EE**2/CW*Sqrt2*Zl11
      AAABR(393) = EE**2/SW*Sqrt2*Zl11
      AAABR(394) = EE**2/CW*Sqrt2
      AAABR(395) = EE**2/SW*Sqrt2
      AAABR(396) = EE**2/CW*SW
      AAABR(397) = EE**2*(SW/CW*Zl21**2+2D0*SW/CW*Zl22**2
     &      -CW/SW*Zl21**2)
      AAABR(398) = EE**2/CW**2*SW**2
      AAABR(399) = EE**2*(2D0*Zl21**2-SW**2/CW**2*Zl21**2
     &      -4D0*SW**2/CW**2*Zl22**2-CW**2/SW**2*Zl21**2)
      AAABR(400) = EE**2/SW**2*Zl21**2
      AAABR(401) = EE**2/CW*Sqrt2*Zl21
      AAABR(402) = EE**2/SW*Sqrt2*Zl21
      AAABR(403) = EE**2*(2D0+SW**2/CW**2+CW**2/SW**2)
      AAABR(404) = EE**2*(SW/CW*Zss11**2-2D0*SW/CW*Zss12
     &      **2+3D0*CW/SW*Zss11**2)
      AAABR(405) = EE**2*(SW/CW*Zb11**2-2D0*SW/CW*Zb12**2
     &      +3D0*CW/SW*Zb11**2)
      AAABR(406) = EE**2*(SW/CW*Zdd11**2-2D0*SW/CW*Zdd12
     &      **2+3D0*CW/SW*Zdd11**2)
      AAABR(407) = EE*GG*(SW/CW*Zb11**2-2D0*SW/CW*Zb12**2
     &      +3D0*CW/SW*Zb11**2)
      AAABR(408) = EE*GG
      AAABR(409) = EE*GG*(SW/CW*Zss11**2-2D0*SW/CW*Zss12
     &      **2+3D0*CW/SW*Zss11**2)
      AAABR(410) = EE*GG*(SW/CW*Zdd11**2-2D0*SW/CW*Zdd12
     &      **2+3D0*CW/SW*Zdd11**2)
      AAABR(411) = EE**2*(6D0*Zdd11**2+SW**2/CW**2*Zdd11
     &      **2+4D0*SW**2/CW**2*Zdd12**2+9D0*CW**2/SW**2
     &      *Zdd11**2)
      AAABR(412) = EE**2*(6D0*Zss11**2+SW**2/CW**2*Zss11
     &      **2+4D0*SW**2/CW**2*Zss12**2+9D0*CW**2/SW**2
     &      *Zss11**2)
      AAABR(413) = EE**2*(6D0*Zb11**2+SW**2/CW**2*Zb11**2
     &      +4D0*SW**2/CW**2*Zb12**2+9D0*CW**2/SW**2*Zb11
     &      **2)
      AAABR(414) = GG**2
      AAABR(415) = EE**2/SW**2*Zdd11**2
      AAABR(416) = EE**2/SW**2*Zb11**2
      AAABR(417) = EE**2/SW**2*Zss11**2
      AAABR(418) = EE**2*(SW/CW*Zss11*Zss21-2D0*SW/CW*Zss12
     &      *Zss22+3D0*CW/SW*Zss11*Zss21)
      AAABR(419) = EE**2*(SW/CW*Zb11*Zb21-2D0*SW/CW*Zb12
     &      *Zb22+3D0*CW/SW*Zb11*Zb21)
      AAABR(420) = EE**2*(SW/CW*Zdd11*Zdd21-2D0*SW/CW*Zdd12
     &      *Zdd22+3D0*CW/SW*Zdd11*Zdd21)
      AAABR(421) = EE*GG*(SW/CW*Zdd11*Zdd21-2D0*SW/CW*Zdd12
     &      *Zdd22+3D0*CW/SW*Zdd11*Zdd21)
      AAABR(422) = EE*GG*(SW/CW*Zb11*Zb21-2D0*SW/CW*Zb12
     &      *Zb22+3D0*CW/SW*Zb11*Zb21)
      AAABR(423) = EE*GG*(SW/CW*Zss11*Zss21-2D0*SW/CW*Zss12
     &      *Zss22+3D0*CW/SW*Zss11*Zss21)
      AAABR(424) = EE**2*(6D0*Zb11*Zb21+SW**2/CW**2*Zb11
     &      *Zb21+4D0*SW**2/CW**2*Zb12*Zb22+9D0*CW**2/SW
     &      **2*Zb11*Zb21)
      AAABR(425) = EE**2*(6D0*Zdd11*Zdd21+SW**2/CW**2*Zdd11
     &      *Zdd21+4D0*SW**2/CW**2*Zdd12*Zdd22+9D0*CW**2
     &      /SW**2*Zdd11*Zdd21)
      AAABR(426) = EE**2*(6D0*Zss11*Zss21+SW**2/CW**2*Zss11
     &      *Zss21+4D0*SW**2/CW**2*Zss12*Zss22+9D0*CW**2
     &      /SW**2*Zss11*Zss21)
      AAABR(427) = EE**2/SW**2*Zb11*Zb21
      AAABR(428) = EE**2/SW**2*Zdd11*Zdd21
      AAABR(429) = EE**2/SW**2*Zss11*Zss21
      AAABR(430) = EE**2/CW*Sqrt2*Zcc11*Zss11
      AAABR(431) = EE**2/SW*Sqrt2*Zcc11*Zss11
      AAABR(432) = EE**2/CW*Sqrt2*Zdd11*Zuu11
      AAABR(433) = EE**2/SW*Sqrt2*Zdd11*Zuu11
      AAABR(434) = EE**2/CW*Sqrt2*Zb11*Zt11
      AAABR(435) = EE**2/SW*Sqrt2*Zb11*Zt11
      AAABR(436) = EE*GG/SW*Sqrt2*Zcc11*Zss11
      AAABR(437) = EE*GG/SW*Sqrt2*Zdd11*Zuu11
      AAABR(438) = EE*GG/SW*Sqrt2*Zb11*Zt11
      AAABR(439) = EE**2/CW*Sqrt2*Zdd11*Zuu21
      AAABR(440) = EE**2/SW*Sqrt2*Zdd11*Zuu21
      AAABR(441) = EE**2/CW*Sqrt2*Zcc21*Zss11
      AAABR(442) = EE**2/SW*Sqrt2*Zcc21*Zss11
      AAABR(443) = EE**2/CW*Sqrt2*Zb11*Zt21
      AAABR(444) = EE**2/SW*Sqrt2*Zb11*Zt21
      AAABR(445) = EE*GG/SW*Sqrt2*Zdd11*Zuu21
      AAABR(446) = EE*GG/SW*Sqrt2*Zb11*Zt21
      AAABR(447) = EE*GG/SW*Sqrt2*Zcc21*Zss11
      AAABR(448) = EE**2*(SW/CW*Zdd21**2-2D0*SW/CW*Zdd22
     &      **2+3D0*CW/SW*Zdd21**2)
      AAABR(449) = EE**2*(SW/CW*Zb21**2-2D0*SW/CW*Zb22**2
     &      +3D0*CW/SW*Zb21**2)
      AAABR(450) = EE**2*(SW/CW*Zss21**2-2D0*SW/CW*Zss22
     &      **2+3D0*CW/SW*Zss21**2)
      AAABR(451) = EE*GG*(SW/CW*Zss21**2-2D0*SW/CW*Zss22
     &      **2+3D0*CW/SW*Zss21**2)
      AAABR(452) = EE*GG*(SW/CW*Zdd21**2-2D0*SW/CW*Zdd22
     &      **2+3D0*CW/SW*Zdd21**2)
      AAABR(453) = EE*GG*(SW/CW*Zb21**2-2D0*SW/CW*Zb22**2
     &      +3D0*CW/SW*Zb21**2)
      AAABR(454) = EE**2*(6D0*Zdd21**2+SW**2/CW**2*Zdd21
     &      **2+4D0*SW**2/CW**2*Zdd22**2+9D0*CW**2/SW**2
     &      *Zdd21**2)
      AAABR(455) = EE**2*(6D0*Zb21**2+SW**2/CW**2*Zb21**2
     &      +4D0*SW**2/CW**2*Zb22**2+9D0*CW**2/SW**2*Zb21
     &      **2)
      AAABR(456) = EE**2*(6D0*Zss21**2+SW**2/CW**2*Zss21
     &      **2+4D0*SW**2/CW**2*Zss22**2+9D0*CW**2/SW**2
     &      *Zss21**2)
      AAABR(457) = EE**2/SW**2*Zdd21**2
      AAABR(458) = EE**2/SW**2*Zb21**2
      AAABR(459) = EE**2/SW**2*Zss21**2
      AAABR(460) = EE**2/CW*Sqrt2*Zcc11*Zss21
      AAABR(461) = EE**2/SW*Sqrt2*Zcc11*Zss21
      AAABR(462) = EE**2/CW*Sqrt2*Zb21*Zt11
      AAABR(463) = EE**2/SW*Sqrt2*Zb21*Zt11
      AAABR(464) = EE**2/CW*Sqrt2*Zdd21*Zuu11
      AAABR(465) = EE**2/SW*Sqrt2*Zdd21*Zuu11
      AAABR(466) = EE*GG/SW*Sqrt2*Zdd21*Zuu11
      AAABR(467) = EE*GG/SW*Sqrt2*Zb21*Zt11
      AAABR(468) = EE*GG/SW*Sqrt2*Zcc11*Zss21
      AAABR(469) = EE**2/CW*Sqrt2*Zdd21*Zuu21
      AAABR(470) = EE**2/SW*Sqrt2*Zdd21*Zuu21
      AAABR(471) = EE**2/CW*Sqrt2*Zcc21*Zss21
      AAABR(472) = EE**2/SW*Sqrt2*Zcc21*Zss21
      AAABR(473) = EE**2/CW*Sqrt2*Zb21*Zt21
      AAABR(474) = EE**2/SW*Sqrt2*Zb21*Zt21
      AAABR(475) = EE*GG/SW*Sqrt2*Zcc21*Zss21
      AAABR(476) = EE*GG/SW*Sqrt2*Zb21*Zt21
      AAABR(477) = EE*GG/SW*Sqrt2*Zdd21*Zuu21
      AAABR(478) = EE**2*(SW/CW*Zt11**2+4D0*SW/CW*Zt12**2
     &      -3D0*CW/SW*Zt11**2)
      AAABR(479) = EE**2*(SW/CW*Zcc11**2+4D0*SW/CW*Zcc12
     &      **2-3D0*CW/SW*Zcc11**2)
      AAABR(480) = EE**2*(SW/CW*Zuu11**2+4D0*SW/CW*Zuu12
     &      **2-3D0*CW/SW*Zuu11**2)
      AAABR(481) = EE*GG*(SW/CW*Zuu11**2+4D0*SW/CW*Zuu12
     &      **2-3D0*CW/SW*Zuu11**2)
      AAABR(482) = EE*GG*(SW/CW*Zt11**2+4D0*SW/CW*Zt12**2
     &      -3D0*CW/SW*Zt11**2)
      AAABR(483) = EE*GG*(SW/CW*Zcc11**2+4D0*SW/CW*Zcc12
     &      **2-3D0*CW/SW*Zcc11**2)
      AAABR(484) = EE**2*(6D0*Zcc11**2-SW**2/CW**2*Zcc11
     &      **2-16D0*SW**2/CW**2*Zcc12**2-9D0*CW**2/SW**2
     &      *Zcc11**2)
      AAABR(485) = EE**2*(6D0*Zt11**2-SW**2/CW**2*Zt11**2
     &      -16D0*SW**2/CW**2*Zt12**2-9D0*CW**2/SW**2*Zt11
     &      **2)
      AAABR(486) = EE**2*(6D0*Zuu11**2-SW**2/CW**2*Zuu11
     &      **2-16D0*SW**2/CW**2*Zuu12**2-9D0*CW**2/SW**2
     &      *Zuu11**2)
      AAABR(487) = EE**2/SW**2*Zt11**2
      AAABR(488) = EE**2/SW**2*Zuu11**2
      AAABR(489) = EE**2/SW**2*Zcc11**2
      AAABR(490) = EE**2*(SW/CW*Zt11*Zt21+4D0*SW/CW*Zt12
     &      *Zt22-3D0*CW/SW*Zt11*Zt21)
      AAABR(491) = EE**2*(SW/CW*Zuu11*Zuu21+4D0*SW/CW*Zuu12
     &      *Zuu22-3D0*CW/SW*Zuu11*Zuu21)
      AAABR(492) = EE**2*(SW/CW*Zcc11*Zcc21+4D0*SW/CW*Zcc12
     &      *Zcc22-3D0*CW/SW*Zcc11*Zcc21)
      AAABR(493) = EE*GG*(SW/CW*Zcc11*Zcc21+4D0*SW/CW*Zcc12
     &      *Zcc22-3D0*CW/SW*Zcc11*Zcc21)
      AAABR(494) = EE*GG*(SW/CW*Zt11*Zt21+4D0*SW/CW*Zt12
     &      *Zt22-3D0*CW/SW*Zt11*Zt21)
      AAABR(495) = EE*GG*(SW/CW*Zuu11*Zuu21+4D0*SW/CW*Zuu12
     &      *Zuu22-3D0*CW/SW*Zuu11*Zuu21)
      AAABR(496) = EE**2*(6D0*Zcc11*Zcc21-SW**2/CW**2*Zcc11
     &      *Zcc21-16D0*SW**2/CW**2*Zcc12*Zcc22-9D0*CW**2
     &      /SW**2*Zcc11*Zcc21)
      AAABR(497) = EE**2*(6D0*Zuu11*Zuu21-SW**2/CW**2*Zuu11
     &      *Zuu21-16D0*SW**2/CW**2*Zuu12*Zuu22-9D0*CW**2
     &      /SW**2*Zuu11*Zuu21)
      AAABR(498) = EE**2*(6D0*Zt11*Zt21-SW**2/CW**2*Zt11
     &      *Zt21-16D0*SW**2/CW**2*Zt12*Zt22-9D0*CW**2/SW
     &      **2*Zt11*Zt21)
      AAABR(499) = EE**2/SW**2*Zcc11*Zcc21
      AAABR(500) = EE**2/SW**2*Zt11*Zt21
      AAABR(501) = EE**2/SW**2*Zuu11*Zuu21
      AAABR(502) = EE**2*(SW/CW*Zcc21**2+4D0*SW/CW*Zcc22
     &      **2-3D0*CW/SW*Zcc21**2)
      AAABR(503) = EE**2*(SW/CW*Zt21**2+4D0*SW/CW*Zt22**2
     &      -3D0*CW/SW*Zt21**2)
      AAABR(504) = EE**2*(SW/CW*Zuu21**2+4D0*SW/CW*Zuu22
     &      **2-3D0*CW/SW*Zuu21**2)
      AAABR(505) = EE*GG*(SW/CW*Zcc21**2+4D0*SW/CW*Zcc22
     &      **2-3D0*CW/SW*Zcc21**2)
      AAABR(506) = EE*GG*(SW/CW*Zt21**2+4D0*SW/CW*Zt22**2
     &      -3D0*CW/SW*Zt21**2)
      AAABR(507) = EE*GG*(SW/CW*Zuu21**2+4D0*SW/CW*Zuu22
     &      **2-3D0*CW/SW*Zuu21**2)
      AAABR(508) = EE**2*(6D0*Zt21**2-SW**2/CW**2*Zt21**2
     &      -16D0*SW**2/CW**2*Zt22**2-9D0*CW**2/SW**2*Zt21
     &      **2)
      AAABR(509) = EE**2*(6D0*Zuu21**2-SW**2/CW**2*Zuu21
     &      **2-16D0*SW**2/CW**2*Zuu22**2-9D0*CW**2/SW**2
     &      *Zuu21**2)
      AAABR(510) = EE**2*(6D0*Zcc21**2-SW**2/CW**2*Zcc21
     &      **2-16D0*SW**2/CW**2*Zcc22**2-9D0*CW**2/SW**2
     &      *Zcc21**2)
      AAABR(511) = EE**2/SW**2*Zt21**2
      AAABR(512) = EE**2/SW**2*Zuu21**2
      AAABR(513) = EE**2/SW**2*Zcc21**2
      AAABR(514) = EE**2/CW
      AAABR(515) = EE**2/SW
      AAABR(516) = EE**2/CW*(Zh31*cb+Zh32*sb)
      AAABR(517) = EE**2/SW*(Zh31*cb+Zh32*sb)
      AAABR(518) = EE**2/CW*(Zh11*cb+Zh12*sb)
      AAABR(519) = EE**2/SW*(Zh11*cb+Zh12*sb)
      AAABR(520) = EE**2/CW*(Zh21*cb+Zh22*sb)
      AAABR(521) = EE**2/SW*(Zh21*cb+Zh22*sb)
      AAABR(522) = EE**2*Pa11**2*(2D0+SW**2/CW**2+CW**2/SW
     &      **2)
      AAABR(523) = EE**2*Pa11*Pa21*(2D0+SW**2/CW**2+CW**2
     &      /SW**2)
      AAABR(524) = EE**2*Pa21**2*(2D0+SW**2/CW**2+CW**2/SW
     &      **2)
      AAABR(525) = EE**2*Pa11**2/SW**2
      AAABR(526) = EE**2*Pa11*Pa21/SW**2
      AAABR(527) = EE**2*Pa21**2/SW**2
      AAABR(528) = EE**2*(2D0*Zh21*Zh31+2D0*Zh22*Zh32+SW
     &      **2/CW**2*Zh21*Zh31+SW**2/CW**2*Zh22*Zh32+CW
     &      **2/SW**2*Zh21*Zh31+CW**2/SW**2*Zh22*Zh32)
      AAABR(529) = EE**2*(2D0*Zh11*Zh21+2D0*Zh12*Zh22+SW
     &      **2/CW**2*Zh11*Zh21+SW**2/CW**2*Zh12*Zh22+CW
     &      **2/SW**2*Zh11*Zh21+CW**2/SW**2*Zh12*Zh22)
      AAABR(530) = EE**2*(2D0*Zh31**2+2D0*Zh32**2+SW**2/CW
     &      **2*Zh31**2+SW**2/CW**2*Zh32**2+CW**2/SW**2*Zh31
     &      **2+CW**2/SW**2*Zh32**2)
      AAABR(531) = EE**2*(2D0*Zh11*Zh31+2D0*Zh12*Zh32+SW
     &      **2/CW**2*Zh11*Zh31+SW**2/CW**2*Zh12*Zh32+CW
     &      **2/SW**2*Zh11*Zh31+CW**2/SW**2*Zh12*Zh32)
      AAABR(532) = EE**2*(2D0*Zh21**2+2D0*Zh22**2+SW**2/CW
     &      **2*Zh21**2+SW**2/CW**2*Zh22**2+CW**2/SW**2*Zh21
     &      **2+CW**2/SW**2*Zh22**2)
      AAABR(533) = EE**2*(2D0*Zh11**2+2D0*Zh12**2+SW**2/CW
     &      **2*Zh11**2+SW**2/CW**2*Zh12**2+CW**2/SW**2*Zh11
     &      **2+CW**2/SW**2*Zh12**2)
      AAABR(534) = EE**2/SW**2*(Zh21**2+Zh22**2)
      AAABR(535) = EE**2/SW**2*(Zh11*Zh31+Zh12*Zh32)
      AAABR(536) = EE**2/SW**2*(Zh11*Zh21+Zh12*Zh22)
      AAABR(537) = EE**2/SW**2*(Zh11**2+Zh12**2)
      AAABR(538) = EE**2/SW**2*(Zh31**2+Zh32**2)
      AAABR(539) = EE**2/SW**2*(Zh21*Zh31+Zh22*Zh32)
      AAABR(540) = CW*EE**2/SW
      AAABR(541) = CW**2*EE**2/SW**2
      end

      subroutine mtrini
      implicit none
#include "model.h"

      integer m1,m2,m3,m4

      do m1 = 1,3
      MTR005(m1)=0D0
      MTR006(m1)=0D0
      MTR007(m1)=0D0
      MTR008(m1)=0D0
      MTR009(m1)=0D0
      MTR011(m1)=0D0
      MTR012(m1)=0D0
      MTR018(m1)=0D0
      MTR019(m1)=0D0
      MTR020(m1)=0D0
      MTR024(m1)=0D0
      MTR025(m1)=0D0
      MTR027(m1)=0D0
      MTR032(m1)=0D0
      MTR037(m1)=0D0
      MTR038(m1)=0D0
      MTR039(m1)=0D0
      MTR040(m1)=0D0
      MTR041(m1)=0D0
      MTR042(m1)=0D0
      MTR043(m1)=0D0
      MTR044(m1)=0D0
      MTR045(m1)=0D0
      MTR046(m1)=0D0
      MTR047(m1)=0D0
      MTR048(m1)=0D0
      MTR049(m1)=0D0
      MTR050(m1)=0D0
      MTR051(m1)=0D0
      MTR052(m1)=0D0
      MTR054(m1)=0D0
      MTR084(m1)=0D0
      MTR087(m1)=0D0
      MTR088(m1)=0D0
      MTR089(m1)=0D0
      MTR090(m1)=0D0
      MTR095(m1)=0D0
      MTR096(m1)=0D0
      MTR097(m1)=0D0
      MTR098(m1)=0D0
      MTR103(m1)=0D0
      MTR106(m1)=0D0
      MTR107(m1)=0D0
      MTR108(m1)=0D0
      MTR109(m1)=0D0
      MTR143(m1)=0D0
      MTR144(m1)=0D0
      MTR145(m1)=0D0
      MTR146(m1)=0D0
      MTR147(m1)=0D0
      MTR148(m1)=0D0
      MTR149(m1)=0D0
      MTR150(m1)=0D0
      MTR151(m1)=0D0
      MTR152(m1)=0D0
      MTR153(m1)=0D0
      MTR154(m1)=0D0
      MTR155(m1)=0D0
      MTR156(m1)=0D0
      MTR157(m1)=0D0
      MTR158(m1)=0D0
      MTR159(m1)=0D0
      MTR160(m1)=0D0
      MTR161(m1)=0D0
      MTR162(m1)=0D0
      MTR163(m1)=0D0
      MTR164(m1)=0D0
      MTR165(m1)=0D0
      MTR166(m1)=0D0
      MTR167(m1)=0D0
      MTR168(m1)=0D0
      MTR169(m1)=0D0
      MTR170(m1)=0D0
      MTR171(m1)=0D0
      MTR172(m1)=0D0
      MTR173(m1)=0D0
      MTR174(m1)=0D0
      MTR175(m1)=0D0
      MTR176(m1)=0D0
      MTR177(m1)=0D0
      MTR178(m1)=0D0
      MTR179(m1)=0D0
      MTR180(m1)=0D0
      MTR181(m1)=0D0
      MTR182(m1)=0D0
      MTR183(m1)=0D0
      MTR184(m1)=0D0
      MTR185(m1)=0D0
      MTR186(m1)=0D0
      MTR187(m1)=0D0
      MTR188(m1)=0D0
      MTR189(m1)=0D0
      enddo
      do m1 = 1,2
      do m2 = 1,3
      MTR013(m1,m2)=0D0
      MTR021(m1,m2)=0D0
      MTR028(m1,m2)=0D0
      MTR033(m1,m2)=0D0
      MTR053(m1,m2)=0D0
      MTR061(m1,m2)=0D0
      MTR062(m1,m2)=0D0
      MTR067(m1,m2)=0D0
      MTR068(m1,m2)=0D0
      MTR069(m1,m2)=0D0
      MTR070(m1,m2)=0D0
      MTR077(m1,m2)=0D0
      MTR078(m1,m2)=0D0
      MTR091(m1,m2)=0D0
      MTR092(m1,m2)=0D0
      MTR093(m1,m2)=0D0
      MTR094(m1,m2)=0D0
      enddo
      enddo
      do m1 = 1,3
      do m2 = 1,3
      MTR015(m1,m2)=0D0
      MTR016(m1,m2)=0D0
      MTR017(m1,m2)=0D0
      MTR022(m1,m2)=0D0
      MTR023(m1,m2)=0D0
      MTR026(m1,m2)=0D0
      MTR029(m1,m2)=0D0
      MTR030(m1,m2)=0D0
      MTR031(m1,m2)=0D0
      MTR072(m1,m2)=0D0
      MTR086(m1,m2)=0D0
      MTR105(m1,m2)=0D0
      MTR192(m1,m2)=0D0
      MTR193(m1,m2)=0D0
      enddo
      enddo
      do m1 = 1,3
      do m2 = 1,2
      MTR071(m1,m2)=0D0
      MTR085(m1,m2)=0D0
      MTR104(m1,m2)=0D0
      enddo
      enddo
      do m1 = 1,5
      do m2 = 1,3
      MTR074(m1,m2)=0D0
      MTR075(m1,m2)=0D0
      MTR076(m1,m2)=0D0
      MTR079(m1,m2)=0D0
      MTR080(m1,m2)=0D0
      MTR081(m1,m2)=0D0
      MTR082(m1,m2)=0D0
      MTR083(m1,m2)=0D0
      MTR099(m1,m2)=0D0
      MTR100(m1,m2)=0D0
      MTR101(m1,m2)=0D0
      MTR102(m1,m2)=0D0
      MTR110(m1,m2)=0D0
      MTR111(m1,m2)=0D0
      MTR112(m1,m2)=0D0
      MTR113(m1,m2)=0D0
      MTR114(m1,m2)=0D0
      MTR115(m1,m2)=0D0
      MTR116(m1,m2)=0D0
      MTR117(m1,m2)=0D0
      MTR118(m1,m2)=0D0
      MTR119(m1,m2)=0D0
      MTR120(m1,m2)=0D0
      MTR121(m1,m2)=0D0
      MTR122(m1,m2)=0D0
      enddo
      enddo
      quuMass(1)=Mq
      quuMass(2)=Mc
      quuMass(3)=Mt
      qudMass(1)=Mq
      qudMass(2)=Mq
      qudMass(3)=Mb
      lpdMass(1)=Mele
      lpdMass(2)=Mmu
      lpdMass(3)=Ml
      neuMass(1)=MNE1
      neuMass(2)=MNE2
      neuMass(3)=MNE3
      neuMass(4)=MNE4
      neuMass(5)=MNE5
      chaMass(1)=MC1
      chaMass(2)=MC2
      sluMass(1)=MSne
      sluMass(2)=MSnm
      sluMass(3)=MSnl
      sldMass(1)=MSeL
      sldMass(2)=MSmL
      sldMass(3)=MSl1
      sleMass(1)=MSeR
      sleMass(2)=MSmR
      sleMass(3)=MSl2
      squMass(1)=MSuL
      squMass(2)=MScL
      squMass(3)=MSt1
      sqvMass(1)=MSuR
      sqvMass(2)=MScR
      sqvMass(3)=MSt2
      sqdMass(1)=MSdL
      sqdMass(2)=MSsL
      sqdMass(3)=MSb1
      sqeMass(1)=MSdR
      sqeMass(2)=MSsR
      sqeMass(3)=MSb2
      hisMass(1)=Mh1
      hisMass(2)=Mh2
      hisMass(3)=Mh3
      hiaMass(1)=Mha
      hiaMass(2)=Mhb
      MTR001(1)=AAABR(9)
      MTR001(2)=AAABR(7)
      MTR001(3)=AAABR(8)
      MTR002(1)=AAABR(12)
      MTR002(2)=AAABR(14)
      MTR002(3)=AAABR(13)
      MTR003(1)=AAABR(18)
      MTR003(2)=AAABR(16)
      MTR003(3)=AAABR(17)
      MTR004(1)=2D0*AAABR(19)
      MTR004(2)=2D0*AAABR(19)
      MTR004(3)=AAABR(20)
      MTR005(3)=1D0
      MTR006(1)=AAABR(22)
      MTR006(2)=AAABR(24)
      MTR006(3)=AAABR(23)
      MTR007(1)=AAABR(27)
      MTR007(2)=AAABR(25)
      MTR007(3)=AAABR(26)
      MTR008(1)=AAABR(30)
      MTR008(2)=AAABR(29)
      MTR008(3)=AAABR(28)
      MTR009(1)=AAABR(33)
      MTR009(2)=AAABR(32)
      MTR009(3)=AAABR(31)
      MTR010(1)=AAABR(34)
      MTR010(2)=AAABR(35)
      MTR011(1)=AAABR(36)
      MTR011(2)=AAABR(37)
      MTR011(3)=AAABR(38)
      MTR012(1)=AAABR(39)
      MTR012(2)=AAABR(39)
      MTR012(3)=AAABR(40)
      MTR013(1,3)=AAABR(42)
      MTR013(2,3)=AAABR(43)
      MTR014(1,1)=AAABR(46)
      MTR014(1,2)=AAABR(46)
      MTR014(1,3)=-AAABR(45)
      MTR014(2,1)=AAABR(44)
      MTR014(2,2)=AAABR(44)
      MTR014(2,3)=-AAABR(47)
      MTR014(3,1)=AAABR(48)
      MTR014(3,2)=AAABR(48)
      MTR014(3,3)=-AAABR(49)
      MTR015(1,3)=AAABR(50)
      MTR015(2,3)=AAABR(52)
      MTR015(3,3)=AAABR(51)
      MTR016(1,1)=2D0*AAABR(56)
      MTR016(1,2)=2D0*AAABR(56)
      MTR016(1,3)=-AAABR(59)
      MTR016(2,1)=2D0*AAABR(58)
      MTR016(2,2)=2D0*AAABR(58)
      MTR016(2,3)=-AAABR(55)
      MTR016(3,1)=2D0*AAABR(54)
      MTR016(3,2)=2D0*AAABR(54)
      MTR016(3,3)=-AAABR(57)
      MTR017(1,1)=AAABR(61)
      MTR017(1,2)=AAABR(61)
      MTR017(1,3)=AAABR(61)
      MTR017(2,1)=AAABR(62)
      MTR017(2,2)=AAABR(62)
      MTR017(2,3)=AAABR(62)
      MTR017(3,1)=AAABR(60)
      MTR017(3,2)=AAABR(60)
      MTR017(3,3)=AAABR(60)
      MTR018(1)=AAABR(64)
      MTR018(2)=AAABR(63)
      MTR018(3)=AAABR(65)
      MTR019(1)=AAABR(68)
      MTR019(2)=AAABR(67)
      MTR019(3)=AAABR(66)
      MTR020(1)=AAABR(71)
      MTR020(2)=AAABR(69)
      MTR020(3)=AAABR(70)
      MTR021(1,1)=AAABR(73)
      MTR021(1,2)=AAABR(75)
      MTR021(1,3)=AAABR(76)
      MTR021(2,1)=AAABR(77)
      MTR021(2,2)=AAABR(72)
      MTR021(2,3)=AAABR(74)
      MTR022(1,1)=AAABR(80)
      MTR022(1,2)=AAABR(79)
      MTR022(1,3)=AAABR(83)
      MTR022(2,1)=AAABR(81)
      MTR022(2,2)=AAABR(84)
      MTR022(2,3)=AAABR(86)
      MTR022(3,1)=AAABR(82)
      MTR022(3,2)=AAABR(78)
      MTR022(3,3)=AAABR(85)
      MTR023(1,1)=AAABR(90)
      MTR023(1,2)=AAABR(94)
      MTR023(1,3)=AAABR(88)
      MTR023(2,1)=AAABR(87)
      MTR023(2,2)=AAABR(95)
      MTR023(2,3)=AAABR(92)
      MTR023(3,1)=AAABR(91)
      MTR023(3,2)=AAABR(93)
      MTR023(3,3)=AAABR(89)
      MTR024(1)=AAABR(97)
      MTR024(2)=AAABR(96)
      MTR024(3)=AAABR(98)
      MTR025(1)=AAABR(100)
      MTR025(2)=AAABR(101)
      MTR025(3)=AAABR(99)
      MTR026(1,1)=AAABR(110)
      MTR026(1,2)=AAABR(102)
      MTR026(1,3)=AAABR(107)
      MTR026(2,1)=AAABR(108)
      MTR026(2,2)=AAABR(109)
      MTR026(2,3)=AAABR(106)
      MTR026(3,1)=AAABR(104)
      MTR026(3,2)=AAABR(105)
      MTR026(3,3)=AAABR(103)
      MTR027(1)=AAABR(113)
      MTR027(2)=AAABR(111)
      MTR027(3)=AAABR(112)
      MTR028(1,1)=AAABR(116)
      MTR028(1,2)=AAABR(118)
      MTR028(1,3)=AAABR(119)
      MTR028(2,1)=AAABR(114)
      MTR028(2,2)=AAABR(115)
      MTR028(2,3)=AAABR(117)
      MTR029(1,1)=AAABR(121)
      MTR029(1,2)=AAABR(120)
      MTR029(1,3)=AAABR(122)
      MTR029(2,1)=AAABR(127)
      MTR029(2,2)=AAABR(124)
      MTR029(2,3)=AAABR(125)
      MTR029(3,1)=AAABR(126)
      MTR029(3,2)=AAABR(128)
      MTR029(3,3)=AAABR(123)
      MTR030(1,1)=AAABR(135)
      MTR030(1,2)=AAABR(130)
      MTR030(1,3)=AAABR(133)
      MTR030(2,1)=AAABR(131)
      MTR030(2,2)=AAABR(132)
      MTR030(2,3)=AAABR(134)
      MTR030(3,1)=AAABR(137)
      MTR030(3,2)=AAABR(136)
      MTR030(3,3)=AAABR(129)
      MTR031(1,1)=AAABR(138)
      MTR031(1,2)=AAABR(140)
      MTR031(1,3)=AAABR(142)
      MTR031(2,1)=AAABR(143)
      MTR031(2,2)=AAABR(145)
      MTR031(2,3)=AAABR(141)
      MTR031(3,1)=AAABR(139)
      MTR031(3,2)=AAABR(144)
      MTR031(3,3)=AAABR(146)
      MTR032(1)=AAABR(149)
      MTR032(2)=AAABR(147)
      MTR032(3)=AAABR(148)
      MTR033(1,1)=AAABR(151)
      MTR033(1,2)=AAABR(153)
      MTR033(1,3)=AAABR(155)
      MTR033(2,1)=AAABR(154)
      MTR033(2,2)=AAABR(150)
      MTR033(2,3)=AAABR(152)
      MTR034(1,1,1)=2D0*AAABR(164)
      MTR034(1,1,2)=2D0*AAABR(156)
      MTR034(1,1,3)=2D0*AAABR(157)
      MTR034(1,2,1)=AAABR(163)
      MTR034(1,2,2)=AAABR(159)
      MTR034(1,2,3)=AAABR(161)
      MTR034(2,1,1)=MTR034(1,2,1)
      MTR034(2,1,2)=MTR034(1,2,2)
      MTR034(2,1,3)=MTR034(1,2,3)
      MTR034(2,2,1)=2D0*AAABR(158)
      MTR034(2,2,2)=2D0*AAABR(160)
      MTR034(2,2,3)=2D0*AAABR(162)
      MTR035(1,1,1)=6D0*AAABR(171)
      MTR035(1,1,2)=2D0*AAABR(165)
      MTR035(1,1,3)=2D0*AAABR(173)
      MTR035(1,2,1)=MTR035(1,1,2)
      MTR035(1,2,2)=2D0*AAABR(166)
      MTR035(1,2,3)=AAABR(167)
      MTR035(1,3,1)=MTR035(1,1,3)
      MTR035(1,3,2)=MTR035(1,2,3)
      MTR035(1,3,3)=2D0*AAABR(168)
      MTR035(2,1,1)=MTR035(1,1,2)
      MTR035(2,1,2)=MTR035(1,2,2)
      MTR035(2,1,3)=MTR035(1,2,3)
      MTR035(2,2,1)=MTR035(1,2,2)
      MTR035(2,2,2)=6D0*AAABR(174)
      MTR035(2,2,3)=2D0*AAABR(172)
      MTR035(2,3,1)=MTR035(1,2,3)
      MTR035(2,3,2)=MTR035(2,2,3)
      MTR035(2,3,3)=2D0*AAABR(169)
      MTR035(3,1,1)=MTR035(1,1,3)
      MTR035(3,1,2)=MTR035(1,2,3)
      MTR035(3,1,3)=MTR035(1,3,3)
      MTR035(3,2,1)=MTR035(1,2,3)
      MTR035(3,2,2)=MTR035(2,2,3)
      MTR035(3,2,3)=MTR035(2,3,3)
      MTR035(3,3,1)=MTR035(1,3,3)
      MTR035(3,3,2)=MTR035(2,3,3)
      MTR035(3,3,3)=6D0*AAABR(170)
      MTR036(1)=AAABR(177)
      MTR036(2)=AAABR(176)
      MTR037(1)=AAABR(180)
      MTR037(2)=AAABR(179)
      MTR037(3)=AAABR(178)
      MTR038(1)=AAABR(175)
      MTR038(2)=AAABR(175)
      MTR038(3)=AAABR(181)
      MTR039(1)=AAABR(184)
      MTR039(2)=AAABR(184)
      MTR039(3)=AAABR(183)
      MTR040(1)=2D0*AAABR(185)
      MTR040(2)=2D0*AAABR(185)
      MTR040(3)=AAABR(186)
      MTR041(1)=AAABR(189)
      MTR041(2)=AAABR(191)
      MTR041(3)=AAABR(190)
      MTR042(1)=AAABR(192)
      MTR042(2)=AAABR(193)
      MTR042(3)=AAABR(194)
      MTR043(1)=AAABR(196)
      MTR043(2)=AAABR(197)
      MTR043(3)=AAABR(195)
      MTR044(1)=AAABR(198)
      MTR044(2)=AAABR(199)
      MTR044(3)=AAABR(200)
      MTR045(1)=AAABR(203)
      MTR045(2)=AAABR(201)
      MTR045(3)=AAABR(202)
      MTR046(1)=AAABR(205)
      MTR046(2)=AAABR(204)
      MTR046(3)=AAABR(206)
      MTR047(1)=AAABR(209)
      MTR047(2)=AAABR(207)
      MTR047(3)=AAABR(208)
      MTR048(1)=AAABR(211)
      MTR048(2)=AAABR(210)
      MTR048(3)=AAABR(212)
      MTR049(1)=AAABR(215)
      MTR049(2)=AAABR(214)
      MTR049(3)=AAABR(213)
      MTR050(1)=AAABR(218)
      MTR050(2)=AAABR(217)
      MTR050(3)=AAABR(216)
      MTR051(1)=AAABR(220)
      MTR051(2)=AAABR(221)
      MTR051(3)=AAABR(222)
      MTR052(1)=AAABR(225)
      MTR052(2)=AAABR(223)
      MTR052(3)=AAABR(224)
      MTR053(1,1)=AAABR(226)
      MTR053(1,2)=AAABR(231)
      MTR053(1,3)=AAABR(227)
      MTR053(2,1)=AAABR(230)
      MTR053(2,2)=AAABR(228)
      MTR053(2,3)=AAABR(229)
      MTR054(1)=AAABR(233)
      MTR054(2)=AAABR(234)
      MTR054(3)=AAABR(235)
      MTR055(1,1)=AAABC(4)
      MTR055(1,2)=AAABC(2)
      MTR055(2,1)=dconjg(AAABC(3))
      MTR055(2,2)=AAABC(1)
      MTR056(1,1)=dconjg(AAABC(4))
      MTR056(1,2)=AAABC(3)
      MTR056(2,1)=dconjg(AAABC(2))
      MTR056(2,2)=dconjg(AAABC(1))
      MTR057(1,1,1)=AAABC(12)
      MTR057(1,1,2)=AAABC(5)
      MTR057(1,2,1)=AAABC(10)
      MTR057(1,2,2)=dconjg(AAABC(7))
      MTR057(2,1,1)=dconjg(AAABC(11))
      MTR057(2,1,2)=AAABC(6)
      MTR057(2,2,1)=AAABC(8)
      MTR057(2,2,2)=AAABC(9)
      MTR058(1,1,1)=dconjg(AAABC(12))
      MTR058(1,1,2)=dconjg(AAABC(5))
      MTR058(1,2,1)=AAABC(11)
      MTR058(1,2,2)=dconjg(AAABC(6))
      MTR058(2,1,1)=dconjg(AAABC(10))
      MTR058(2,1,2)=AAABC(7)
      MTR058(2,2,1)=dconjg(AAABC(8))
      MTR058(2,2,2)=dconjg(AAABC(9))
      MTR059(1,1,1)=AAABC(15)
      MTR059(1,1,2)=AAABC(19)
      MTR059(1,1,3)=AAABC(24)
      MTR059(1,2,1)=dconjg(AAABC(21))
      MTR059(1,2,2)=dconjg(AAABC(17))
      MTR059(1,2,3)=dconjg(AAABC(23))
      MTR059(2,1,1)=AAABC(20)
      MTR059(2,1,2)=AAABC(16)
      MTR059(2,1,3)=AAABC(22)
      MTR059(2,2,1)=AAABC(13)
      MTR059(2,2,2)=AAABC(14)
      MTR059(2,2,3)=AAABC(18)
      MTR060(1,1,1)=dconjg(AAABC(15))
      MTR060(1,1,2)=dconjg(AAABC(19))
      MTR060(1,1,3)=dconjg(AAABC(24))
      MTR060(1,2,1)=dconjg(AAABC(20))
      MTR060(1,2,2)=dconjg(AAABC(16))
      MTR060(1,2,3)=dconjg(AAABC(22))
      MTR060(2,1,1)=AAABC(21)
      MTR060(2,1,2)=AAABC(17)
      MTR060(2,1,3)=AAABC(23)
      MTR060(2,2,1)=dconjg(AAABC(13))
      MTR060(2,2,2)=dconjg(AAABC(14))
      MTR060(2,2,3)=dconjg(AAABC(18))
      MTR061(1,1)=2D0*AAABR(236)
      MTR061(1,2)=2D0*AAABR(236)
      MTR061(1,3)=AAABR(238)
      MTR061(2,1)=2D0*AAABR(237)
      MTR061(2,2)=2D0*AAABR(237)
      MTR061(2,3)=AAABR(239)
      MTR062(1,3)=AAABR(240)
      MTR062(2,3)=AAABR(241)
      MTR063(1,1)=AAABC(33)
      MTR063(1,2)=AAABC(25)
      MTR063(1,3)=AAABC(31)
      MTR063(1,4)=AAABC(29)
      MTR063(1,5)=AAABC(27)
      MTR063(2,1)=AAABC(39)
      MTR063(2,2)=AAABC(37)
      MTR063(2,3)=AAABC(35)
      MTR063(2,4)=AAABC(41)
      MTR063(2,5)=AAABC(43)
      MTR064(1,1)=AAABC(34)
      MTR064(1,2)=AAABC(26)
      MTR064(1,3)=AAABC(32)
      MTR064(1,4)=AAABC(30)
      MTR064(1,5)=AAABC(28)
      MTR064(2,1)=AAABC(40)
      MTR064(2,2)=AAABC(38)
      MTR064(2,3)=AAABC(36)
      MTR064(2,4)=AAABC(42)
      MTR064(2,5)=AAABC(44)
      MTR065(1,1)=AAABC(57)
      MTR065(1,2)=AAABC(61)
      MTR065(1,3)=AAABC(53)
      MTR065(1,4)=AAABC(51)
      MTR065(1,5)=AAABC(59)
      MTR065(2,1)=AAABC(55)
      MTR065(2,2)=AAABC(49)
      MTR065(2,3)=AAABC(47)
      MTR065(2,4)=AAABC(63)
      MTR065(2,5)=AAABC(45)
      MTR066(1,1)=AAABC(58)
      MTR066(1,2)=AAABC(62)
      MTR066(1,3)=AAABC(54)
      MTR066(1,4)=AAABC(52)
      MTR066(1,5)=AAABC(60)
      MTR066(2,1)=AAABC(56)
      MTR066(2,2)=AAABC(50)
      MTR066(2,3)=AAABC(48)
      MTR066(2,4)=AAABC(64)
      MTR066(2,5)=AAABC(46)
      MTR067(1,1)=AAABR(252)
      MTR067(1,2)=AAABR(244)
      MTR067(1,3)=AAABR(242)
      MTR067(2,1)=AAABR(250)
      MTR067(2,2)=AAABR(246)
      MTR067(2,3)=AAABR(248)
      MTR068(1,1)=AAABR(253)
      MTR068(1,2)=AAABR(245)
      MTR068(1,3)=AAABR(243)
      MTR068(2,1)=AAABR(251)
      MTR068(2,2)=AAABR(247)
      MTR068(2,3)=AAABR(249)
      MTR069(1,1)=AAABR(262)
      MTR069(1,2)=AAABR(256)
      MTR069(1,3)=AAABR(258)
      MTR069(2,1)=AAABR(260)
      MTR069(2,2)=AAABR(264)
      MTR069(2,3)=AAABR(254)
      MTR070(1,1)=AAABR(263)
      MTR070(1,2)=AAABR(257)
      MTR070(1,3)=AAABR(259)
      MTR070(2,1)=AAABR(261)
      MTR070(2,2)=AAABR(265)
      MTR070(2,3)=AAABR(255)
      MTR071(3,1)=AAABR(267)
      MTR071(3,2)=AAABR(268)
      MTR072(3,1)=AAABR(270)
      MTR072(3,2)=AAABR(271)
      MTR072(3,3)=AAABR(269)
      MTR073(1,1)=AAABC(66)
      MTR073(1,2)=AAABC(66)
      MTR073(1,3)=AAABC(70)
      MTR073(2,1)=AAABC(65)
      MTR073(2,2)=AAABC(65)
      MTR073(2,3)=AAABC(74)
      MTR073(3,1)=AAABC(73)
      MTR073(3,2)=AAABC(73)
      MTR073(3,3)=AAABC(76)
      MTR073(4,1)=AAABC(67)
      MTR073(4,2)=AAABC(67)
      MTR073(4,3)=AAABC(68)
      MTR073(5,1)=AAABC(72)
      MTR073(5,2)=AAABC(72)
      MTR073(5,3)=AAABC(78)
      MTR074(1,3)=AAABC(71)
      MTR074(2,3)=AAABC(75)
      MTR074(3,3)=AAABC(77)
      MTR074(4,3)=AAABC(69)
      MTR074(5,3)=AAABC(79)
      MTR075(1,1)=2D0*AAABC(80)
      MTR075(1,2)=2D0*AAABC(80)
      MTR075(1,3)=AAABC(94)
      MTR075(2,1)=2D0*AAABC(90)
      MTR075(2,2)=2D0*AAABC(90)
      MTR075(2,3)=AAABC(89)
      MTR075(3,1)=2D0*AAABC(83)
      MTR075(3,2)=2D0*AAABC(83)
      MTR075(3,3)=AAABC(86)
      MTR075(4,1)=2D0*AAABC(84)
      MTR075(4,2)=2D0*AAABC(84)
      MTR075(4,3)=AAABC(82)
      MTR075(5,1)=2D0*AAABC(87)
      MTR075(5,2)=2D0*AAABC(87)
      MTR075(5,3)=AAABC(92)
      MTR076(1,3)=AAABC(93)
      MTR076(2,3)=AAABC(88)
      MTR076(3,3)=AAABC(85)
      MTR076(4,3)=AAABC(81)
      MTR076(5,3)=AAABC(91)
      MTR077(1,1)=AAABR(274)
      MTR077(1,2)=AAABR(274)
      MTR077(1,3)=AAABR(274)
      MTR077(2,1)=AAABR(276)
      MTR077(2,2)=AAABR(276)
      MTR077(2,3)=AAABR(276)
      MTR078(1,3)=AAABR(275)
      MTR078(2,3)=AAABR(277)
      MTR079(1,1)=AAABC(99)
      MTR079(1,2)=AAABC(99)
      MTR079(1,3)=AAABC(99)
      MTR079(2,1)=AAABC(98)
      MTR079(2,2)=AAABC(98)
      MTR079(2,3)=AAABC(98)
      MTR079(3,1)=AAABC(97)
      MTR079(3,2)=AAABC(97)
      MTR079(3,3)=AAABC(97)
      MTR079(4,1)=AAABC(96)
      MTR079(4,2)=AAABC(96)
      MTR079(4,3)=AAABC(96)
      MTR079(5,1)=AAABC(95)
      MTR079(5,2)=AAABC(95)
      MTR079(5,3)=AAABC(95)
      MTR080(1,1)=AAABC(120)
      MTR080(1,2)=AAABC(108)
      MTR080(1,3)=AAABC(104)
      MTR080(2,1)=AAABC(110)
      MTR080(2,2)=AAABC(100)
      MTR080(2,3)=AAABC(124)
      MTR080(3,1)=AAABC(122)
      MTR080(3,2)=AAABC(118)
      MTR080(3,3)=AAABC(112)
      MTR080(4,1)=AAABC(102)
      MTR080(4,2)=AAABC(116)
      MTR080(4,3)=AAABC(106)
      MTR080(5,1)=AAABC(114)
      MTR080(5,2)=AAABC(128)
      MTR080(5,3)=AAABC(126)
      MTR081(1,1)=AAABC(121)
      MTR081(1,2)=AAABC(109)
      MTR081(1,3)=AAABC(105)
      MTR081(2,1)=AAABC(111)
      MTR081(2,2)=AAABC(101)
      MTR081(2,3)=AAABC(125)
      MTR081(3,1)=AAABC(123)
      MTR081(3,2)=AAABC(119)
      MTR081(3,3)=AAABC(113)
      MTR081(4,1)=AAABC(103)
      MTR081(4,2)=AAABC(117)
      MTR081(4,3)=AAABC(107)
      MTR081(5,1)=AAABC(115)
      MTR081(5,2)=AAABC(129)
      MTR081(5,3)=AAABC(127)
      MTR082(1,1)=AAABC(152)
      MTR082(1,2)=AAABC(154)
      MTR082(1,3)=AAABC(142)
      MTR082(2,1)=AAABC(138)
      MTR082(2,2)=AAABC(150)
      MTR082(2,3)=AAABC(144)
      MTR082(3,1)=AAABC(136)
      MTR082(3,2)=AAABC(146)
      MTR082(3,3)=AAABC(148)
      MTR082(4,1)=AAABC(156)
      MTR082(4,2)=AAABC(132)
      MTR082(4,3)=AAABC(134)
      MTR082(5,1)=AAABC(158)
      MTR082(5,2)=AAABC(140)
      MTR082(5,3)=AAABC(130)
      MTR083(1,1)=AAABC(153)
      MTR083(1,2)=AAABC(155)
      MTR083(1,3)=AAABC(143)
      MTR083(2,1)=AAABC(139)
      MTR083(2,2)=AAABC(151)
      MTR083(2,3)=AAABC(145)
      MTR083(3,1)=AAABC(137)
      MTR083(3,2)=AAABC(147)
      MTR083(3,3)=AAABC(149)
      MTR083(4,1)=AAABC(157)
      MTR083(4,2)=AAABC(133)
      MTR083(4,3)=AAABC(135)
      MTR083(5,1)=AAABC(159)
      MTR083(5,2)=AAABC(141)
      MTR083(5,3)=AAABC(131)
      MTR084(1)=AAABR(278)
      MTR084(2)=AAABR(278)
      MTR084(3)=AAABR(279)
      MTR085(1,1)=AAABR(281)
      MTR085(1,2)=AAABR(280)
      MTR085(2,1)=AAABR(281)
      MTR085(2,2)=AAABR(280)
      MTR085(3,1)=AAABR(283)
      MTR085(3,2)=AAABR(282)
      MTR086(1,1)=AAABR(285)
      MTR086(1,2)=AAABR(284)
      MTR086(1,3)=AAABR(286)
      MTR086(2,1)=AAABR(285)
      MTR086(2,2)=AAABR(284)
      MTR086(2,3)=AAABR(286)
      MTR086(3,1)=AAABR(287)
      MTR086(3,2)=AAABR(289)
      MTR086(3,3)=AAABR(288)
      MTR087(1)=AAABR(290)
      MTR087(2)=AAABR(292)
      MTR087(3)=AAABR(294)
      MTR088(1)=AAABR(291)
      MTR088(2)=AAABR(291)
      MTR088(3)=AAABR(293)
      MTR089(1)=AAABR(296)
      MTR089(2)=AAABR(295)
      MTR089(3)=AAABR(298)
      MTR090(1)=AAABR(296)
      MTR090(2)=AAABR(296)
      MTR090(3)=AAABR(297)
      MTR091(1,1)=AAABR(301)
      MTR091(1,2)=AAABR(303)
      MTR091(1,3)=AAABR(299)
      MTR091(2,1)=AAABR(315)
      MTR091(2,2)=AAABR(313)
      MTR091(2,3)=AAABR(311)
      MTR092(1,1)=AAABR(302)
      MTR092(1,2)=AAABR(304)
      MTR092(1,3)=AAABR(300)
      MTR092(2,1)=AAABR(316)
      MTR092(2,2)=AAABR(314)
      MTR092(2,3)=AAABR(312)
      MTR093(1,1)=AAABR(307)
      MTR093(1,2)=AAABR(305)
      MTR093(1,3)=AAABR(309)
      MTR093(2,1)=AAABR(319)
      MTR093(2,2)=AAABR(321)
      MTR093(2,3)=AAABR(317)
      MTR094(1,1)=AAABR(308)
      MTR094(1,2)=AAABR(306)
      MTR094(1,3)=AAABR(310)
      MTR094(2,1)=AAABR(320)
      MTR094(2,2)=AAABR(322)
      MTR094(2,3)=AAABR(318)
      MTR095(1)=AAABR(323)
      MTR095(2)=AAABR(325)
      MTR095(3)=AAABR(327)
      MTR096(1)=AAABR(324)
      MTR096(2)=AAABR(326)
      MTR096(3)=AAABR(328)
      MTR097(1)=AAABR(333)
      MTR097(2)=AAABR(329)
      MTR097(3)=AAABR(331)
      MTR098(1)=AAABR(334)
      MTR098(2)=AAABR(330)
      MTR098(3)=AAABR(332)
      MTR099(1,1)=AAABC(162)
      MTR099(1,2)=AAABC(174)
      MTR099(1,3)=AAABC(160)
      MTR099(2,1)=AAABC(188)
      MTR099(2,2)=AAABC(164)
      MTR099(2,3)=AAABC(168)
      MTR099(3,1)=AAABC(186)
      MTR099(3,2)=AAABC(182)
      MTR099(3,3)=AAABC(176)
      MTR099(4,1)=AAABC(184)
      MTR099(4,2)=AAABC(170)
      MTR099(4,3)=AAABC(178)
      MTR099(5,1)=AAABC(180)
      MTR099(5,2)=AAABC(172)
      MTR099(5,3)=AAABC(166)
      MTR100(1,1)=AAABC(163)
      MTR100(1,2)=AAABC(175)
      MTR100(1,3)=AAABC(161)
      MTR100(2,1)=AAABC(189)
      MTR100(2,2)=AAABC(165)
      MTR100(2,3)=AAABC(169)
      MTR100(3,1)=AAABC(187)
      MTR100(3,2)=AAABC(183)
      MTR100(3,3)=AAABC(177)
      MTR100(4,1)=AAABC(185)
      MTR100(4,2)=AAABC(171)
      MTR100(4,3)=AAABC(179)
      MTR100(5,1)=AAABC(181)
      MTR100(5,2)=AAABC(173)
      MTR100(5,3)=AAABC(167)
      MTR101(1,1)=AAABC(190)
      MTR101(1,2)=AAABC(208)
      MTR101(1,3)=AAABC(210)
      MTR101(2,1)=AAABC(198)
      MTR101(2,2)=AAABC(192)
      MTR101(2,3)=AAABC(214)
      MTR101(3,1)=AAABC(216)
      MTR101(3,2)=AAABC(202)
      MTR101(3,3)=AAABC(204)
      MTR101(4,1)=AAABC(212)
      MTR101(4,2)=AAABC(200)
      MTR101(4,3)=AAABC(218)
      MTR101(5,1)=AAABC(196)
      MTR101(5,2)=AAABC(194)
      MTR101(5,3)=AAABC(206)
      MTR102(1,1)=AAABC(191)
      MTR102(1,2)=AAABC(209)
      MTR102(1,3)=AAABC(211)
      MTR102(2,1)=AAABC(199)
      MTR102(2,2)=AAABC(193)
      MTR102(2,3)=AAABC(215)
      MTR102(3,1)=AAABC(217)
      MTR102(3,2)=AAABC(203)
      MTR102(3,3)=AAABC(205)
      MTR102(4,1)=AAABC(213)
      MTR102(4,2)=AAABC(201)
      MTR102(4,3)=AAABC(219)
      MTR102(5,1)=AAABC(197)
      MTR102(5,2)=AAABC(195)
      MTR102(5,3)=AAABC(207)
      MTR103(1)=AAABR(278)
      MTR103(2)=AAABR(336)
      MTR103(3)=AAABR(335)
      MTR104(1,1)=AAABR(342)
      MTR104(1,2)=AAABR(339)
      MTR104(2,1)=AAABR(340)
      MTR104(2,2)=AAABR(341)
      MTR104(3,1)=AAABR(337)
      MTR104(3,2)=AAABR(338)
      MTR105(1,1)=AAABR(345)
      MTR105(1,2)=AAABR(350)
      MTR105(1,3)=AAABR(344)
      MTR105(2,1)=AAABR(346)
      MTR105(2,2)=AAABR(347)
      MTR105(2,3)=AAABR(349)
      MTR105(3,1)=AAABR(351)
      MTR105(3,2)=AAABR(348)
      MTR105(3,3)=AAABR(343)
      MTR106(1)=AAABR(354)
      MTR106(2)=AAABR(356)
      MTR106(3)=AAABR(352)
      MTR107(1)=AAABR(355)
      MTR107(2)=AAABR(357)
      MTR107(3)=AAABR(353)
      MTR108(1)=AAABR(360)
      MTR108(2)=AAABR(362)
      MTR108(3)=AAABR(358)
      MTR109(1)=AAABR(361)
      MTR109(2)=AAABR(363)
      MTR109(3)=AAABR(359)
      MTR110(1,1)=dconjg(AAABC(120))
      MTR110(1,2)=dconjg(AAABC(108))
      MTR110(1,3)=dconjg(AAABC(104))
      MTR110(2,1)=dconjg(AAABC(110))
      MTR110(2,2)=dconjg(AAABC(100))
      MTR110(2,3)=dconjg(AAABC(124))
      MTR110(3,1)=dconjg(AAABC(122))
      MTR110(3,2)=dconjg(AAABC(118))
      MTR110(3,3)=dconjg(AAABC(112))
      MTR110(4,1)=dconjg(AAABC(102))
      MTR110(4,2)=dconjg(AAABC(116))
      MTR110(4,3)=dconjg(AAABC(106))
      MTR110(5,1)=dconjg(AAABC(114))
      MTR110(5,2)=dconjg(AAABC(128))
      MTR110(5,3)=dconjg(AAABC(126))
      MTR111(1,1)=dconjg(AAABC(121))
      MTR111(1,2)=dconjg(AAABC(109))
      MTR111(1,3)=dconjg(AAABC(105))
      MTR111(2,1)=dconjg(AAABC(111))
      MTR111(2,2)=dconjg(AAABC(101))
      MTR111(2,3)=dconjg(AAABC(125))
      MTR111(3,1)=dconjg(AAABC(123))
      MTR111(3,2)=dconjg(AAABC(119))
      MTR111(3,3)=dconjg(AAABC(113))
      MTR111(4,1)=dconjg(AAABC(103))
      MTR111(4,2)=dconjg(AAABC(117))
      MTR111(4,3)=dconjg(AAABC(107))
      MTR111(5,1)=dconjg(AAABC(115))
      MTR111(5,2)=dconjg(AAABC(129))
      MTR111(5,3)=dconjg(AAABC(127))
      MTR112(1,1)=dconjg(AAABC(152))
      MTR112(1,2)=dconjg(AAABC(154))
      MTR112(1,3)=dconjg(AAABC(142))
      MTR112(2,1)=dconjg(AAABC(138))
      MTR112(2,2)=dconjg(AAABC(150))
      MTR112(2,3)=dconjg(AAABC(144))
      MTR112(3,1)=dconjg(AAABC(136))
      MTR112(3,2)=dconjg(AAABC(146))
      MTR112(3,3)=dconjg(AAABC(148))
      MTR112(4,1)=dconjg(AAABC(156))
      MTR112(4,2)=dconjg(AAABC(132))
      MTR112(4,3)=dconjg(AAABC(134))
      MTR112(5,1)=dconjg(AAABC(158))
      MTR112(5,2)=dconjg(AAABC(140))
      MTR112(5,3)=dconjg(AAABC(130))
      MTR113(1,1)=dconjg(AAABC(153))
      MTR113(1,2)=dconjg(AAABC(155))
      MTR113(1,3)=dconjg(AAABC(143))
      MTR113(2,1)=dconjg(AAABC(139))
      MTR113(2,2)=dconjg(AAABC(151))
      MTR113(2,3)=dconjg(AAABC(145))
      MTR113(3,1)=dconjg(AAABC(137))
      MTR113(3,2)=dconjg(AAABC(147))
      MTR113(3,3)=dconjg(AAABC(149))
      MTR113(4,1)=dconjg(AAABC(157))
      MTR113(4,2)=dconjg(AAABC(133))
      MTR113(4,3)=dconjg(AAABC(135))
      MTR113(5,1)=dconjg(AAABC(159))
      MTR113(5,2)=dconjg(AAABC(141))
      MTR113(5,3)=dconjg(AAABC(131))
      MTR114(1,1)=dconjg(AAABC(162))
      MTR114(1,2)=dconjg(AAABC(174))
      MTR114(1,3)=dconjg(AAABC(160))
      MTR114(2,1)=dconjg(AAABC(188))
      MTR114(2,2)=dconjg(AAABC(164))
      MTR114(2,3)=dconjg(AAABC(168))
      MTR114(3,1)=dconjg(AAABC(186))
      MTR114(3,2)=dconjg(AAABC(182))
      MTR114(3,3)=dconjg(AAABC(176))
      MTR114(4,1)=dconjg(AAABC(184))
      MTR114(4,2)=dconjg(AAABC(170))
      MTR114(4,3)=dconjg(AAABC(178))
      MTR114(5,1)=dconjg(AAABC(180))
      MTR114(5,2)=dconjg(AAABC(172))
      MTR114(5,3)=dconjg(AAABC(166))
      MTR115(1,1)=dconjg(AAABC(163))
      MTR115(1,2)=dconjg(AAABC(175))
      MTR115(1,3)=dconjg(AAABC(161))
      MTR115(2,1)=dconjg(AAABC(189))
      MTR115(2,2)=dconjg(AAABC(165))
      MTR115(2,3)=dconjg(AAABC(169))
      MTR115(3,1)=dconjg(AAABC(187))
      MTR115(3,2)=dconjg(AAABC(183))
      MTR115(3,3)=dconjg(AAABC(177))
      MTR115(4,1)=dconjg(AAABC(185))
      MTR115(4,2)=dconjg(AAABC(171))
      MTR115(4,3)=dconjg(AAABC(179))
      MTR115(5,1)=dconjg(AAABC(181))
      MTR115(5,2)=dconjg(AAABC(173))
      MTR115(5,3)=dconjg(AAABC(167))
      MTR116(1,1)=dconjg(AAABC(190))
      MTR116(1,2)=dconjg(AAABC(208))
      MTR116(1,3)=dconjg(AAABC(210))
      MTR116(2,1)=dconjg(AAABC(198))
      MTR116(2,2)=dconjg(AAABC(192))
      MTR116(2,3)=dconjg(AAABC(214))
      MTR116(3,1)=dconjg(AAABC(216))
      MTR116(3,2)=dconjg(AAABC(202))
      MTR116(3,3)=dconjg(AAABC(204))
      MTR116(4,1)=dconjg(AAABC(212))
      MTR116(4,2)=dconjg(AAABC(200))
      MTR116(4,3)=dconjg(AAABC(218))
      MTR116(5,1)=dconjg(AAABC(196))
      MTR116(5,2)=dconjg(AAABC(194))
      MTR116(5,3)=dconjg(AAABC(206))
      MTR117(1,1)=dconjg(AAABC(191))
      MTR117(1,2)=dconjg(AAABC(209))
      MTR117(1,3)=dconjg(AAABC(211))
      MTR117(2,1)=dconjg(AAABC(199))
      MTR117(2,2)=dconjg(AAABC(193))
      MTR117(2,3)=dconjg(AAABC(215))
      MTR117(3,1)=dconjg(AAABC(217))
      MTR117(3,2)=dconjg(AAABC(203))
      MTR117(3,3)=dconjg(AAABC(205))
      MTR117(4,1)=dconjg(AAABC(213))
      MTR117(4,2)=dconjg(AAABC(201))
      MTR117(4,3)=dconjg(AAABC(219))
      MTR117(5,1)=dconjg(AAABC(197))
      MTR117(5,2)=dconjg(AAABC(195))
      MTR117(5,3)=dconjg(AAABC(207))
      MTR118(1,1)=dconjg(AAABC(66))
      MTR118(1,2)=dconjg(AAABC(66))
      MTR118(1,3)=dconjg(AAABC(70))
      MTR118(2,1)=dconjg(AAABC(65))
      MTR118(2,2)=dconjg(AAABC(65))
      MTR118(2,3)=dconjg(AAABC(74))
      MTR118(3,1)=dconjg(AAABC(73))
      MTR118(3,2)=dconjg(AAABC(73))
      MTR118(3,3)=dconjg(AAABC(76))
      MTR118(4,1)=dconjg(AAABC(67))
      MTR118(4,2)=dconjg(AAABC(67))
      MTR118(4,3)=dconjg(AAABC(68))
      MTR118(5,1)=dconjg(AAABC(72))
      MTR118(5,2)=dconjg(AAABC(72))
      MTR118(5,3)=dconjg(AAABC(78))
      MTR119(1,3)=dconjg(AAABC(71))
      MTR119(2,3)=dconjg(AAABC(75))
      MTR119(3,3)=dconjg(AAABC(77))
      MTR119(4,3)=dconjg(AAABC(69))
      MTR119(5,3)=dconjg(AAABC(79))
      MTR120(1,1)=2D0*dconjg(AAABC(80))
      MTR120(1,2)=2D0*dconjg(AAABC(80))
      MTR120(1,3)=dconjg(AAABC(94))
      MTR120(2,1)=2D0*dconjg(AAABC(90))
      MTR120(2,2)=2D0*dconjg(AAABC(90))
      MTR120(2,3)=dconjg(AAABC(89))
      MTR120(3,1)=2D0*dconjg(AAABC(83))
      MTR120(3,2)=2D0*dconjg(AAABC(83))
      MTR120(3,3)=dconjg(AAABC(86))
      MTR120(4,1)=2D0*dconjg(AAABC(84))
      MTR120(4,2)=2D0*dconjg(AAABC(84))
      MTR120(4,3)=dconjg(AAABC(82))
      MTR120(5,1)=2D0*dconjg(AAABC(87))
      MTR120(5,2)=2D0*dconjg(AAABC(87))
      MTR120(5,3)=dconjg(AAABC(92))
      MTR121(1,3)=dconjg(AAABC(93))
      MTR121(2,3)=dconjg(AAABC(88))
      MTR121(3,3)=dconjg(AAABC(85))
      MTR121(4,3)=dconjg(AAABC(81))
      MTR121(5,3)=dconjg(AAABC(91))
      MTR122(1,1)=dconjg(AAABC(99))
      MTR122(1,2)=dconjg(AAABC(99))
      MTR122(1,3)=dconjg(AAABC(99))
      MTR122(2,1)=dconjg(AAABC(98))
      MTR122(2,2)=dconjg(AAABC(98))
      MTR122(2,3)=dconjg(AAABC(98))
      MTR122(3,1)=dconjg(AAABC(97))
      MTR122(3,2)=dconjg(AAABC(97))
      MTR122(3,3)=dconjg(AAABC(97))
      MTR122(4,1)=dconjg(AAABC(96))
      MTR122(4,2)=dconjg(AAABC(96))
      MTR122(4,3)=dconjg(AAABC(96))
      MTR122(5,1)=dconjg(AAABC(95))
      MTR122(5,2)=dconjg(AAABC(95))
      MTR122(5,3)=dconjg(AAABC(95))
      MTR123(1,1)=dconjg(AAABC(33))
      MTR123(1,2)=dconjg(AAABC(25))
      MTR123(1,3)=dconjg(AAABC(31))
      MTR123(1,4)=dconjg(AAABC(29))
      MTR123(1,5)=dconjg(AAABC(27))
      MTR123(2,1)=dconjg(AAABC(39))
      MTR123(2,2)=dconjg(AAABC(37))
      MTR123(2,3)=dconjg(AAABC(35))
      MTR123(2,4)=dconjg(AAABC(41))
      MTR123(2,5)=dconjg(AAABC(43))
      MTR124(1,1)=dconjg(AAABC(34))
      MTR124(1,2)=dconjg(AAABC(26))
      MTR124(1,3)=dconjg(AAABC(32))
      MTR124(1,4)=dconjg(AAABC(30))
      MTR124(1,5)=dconjg(AAABC(28))
      MTR124(2,1)=dconjg(AAABC(40))
      MTR124(2,2)=dconjg(AAABC(38))
      MTR124(2,3)=dconjg(AAABC(36))
      MTR124(2,4)=dconjg(AAABC(42))
      MTR124(2,5)=dconjg(AAABC(44))
      MTR125(1,1)=dconjg(AAABC(57))
      MTR125(1,2)=dconjg(AAABC(61))
      MTR125(1,3)=dconjg(AAABC(53))
      MTR125(1,4)=dconjg(AAABC(51))
      MTR125(1,5)=dconjg(AAABC(59))
      MTR125(2,1)=dconjg(AAABC(55))
      MTR125(2,2)=dconjg(AAABC(49))
      MTR125(2,3)=dconjg(AAABC(47))
      MTR125(2,4)=dconjg(AAABC(63))
      MTR125(2,5)=dconjg(AAABC(45))
      MTR126(1,1)=dconjg(AAABC(58))
      MTR126(1,2)=dconjg(AAABC(62))
      MTR126(1,3)=dconjg(AAABC(54))
      MTR126(1,4)=dconjg(AAABC(52))
      MTR126(1,5)=dconjg(AAABC(60))
      MTR126(2,1)=dconjg(AAABC(56))
      MTR126(2,2)=dconjg(AAABC(50))
      MTR126(2,3)=dconjg(AAABC(48))
      MTR126(2,4)=dconjg(AAABC(64))
      MTR126(2,5)=dconjg(AAABC(46))
      MTR127(1,1)=2D0*AAABC(220)
      MTR127(1,2)=AAABC(223)
      MTR127(1,3)=AAABC(224)
      MTR127(1,4)=AAABC(221)
      MTR127(1,5)=AAABC(222)
      MTR127(2,1)=AAABC(223)
      MTR127(2,2)=2D0*AAABC(252)
      MTR127(2,3)=AAABC(250)
      MTR127(2,4)=AAABC(251)
      MTR127(2,5)=AAABC(253)
      MTR127(3,1)=AAABC(224)
      MTR127(3,2)=AAABC(250)
      MTR127(3,3)=2D0*AAABC(274)
      MTR127(3,4)=AAABC(275)
      MTR127(3,5)=AAABC(276)
      MTR127(4,1)=AAABC(221)
      MTR127(4,2)=AAABC(251)
      MTR127(4,3)=AAABC(275)
      MTR127(4,4)=2D0*AAABC(293)
      MTR127(4,5)=AAABC(292)
      MTR127(5,1)=AAABC(222)
      MTR127(5,2)=AAABC(253)
      MTR127(5,3)=AAABC(276)
      MTR127(5,4)=AAABC(292)
      MTR127(5,5)=2D0*AAABC(304)
      MTR128(1,1)=2D0*dconjg(AAABC(220))
      MTR128(1,2)=dconjg(AAABC(223))
      MTR128(1,3)=dconjg(AAABC(224))
      MTR128(1,4)=dconjg(AAABC(221))
      MTR128(1,5)=dconjg(AAABC(222))
      MTR128(2,1)=dconjg(AAABC(223))
      MTR128(2,2)=2D0*dconjg(AAABC(252))
      MTR128(2,3)=dconjg(AAABC(250))
      MTR128(2,4)=dconjg(AAABC(251))
      MTR128(2,5)=dconjg(AAABC(253))
      MTR128(3,1)=dconjg(AAABC(224))
      MTR128(3,2)=dconjg(AAABC(250))
      MTR128(3,3)=2D0*dconjg(AAABC(274))
      MTR128(3,4)=dconjg(AAABC(275))
      MTR128(3,5)=dconjg(AAABC(276))
      MTR128(4,1)=dconjg(AAABC(221))
      MTR128(4,2)=dconjg(AAABC(251))
      MTR128(4,3)=dconjg(AAABC(275))
      MTR128(4,4)=2D0*dconjg(AAABC(293))
      MTR128(4,5)=dconjg(AAABC(292))
      MTR128(5,1)=dconjg(AAABC(222))
      MTR128(5,2)=dconjg(AAABC(253))
      MTR128(5,3)=dconjg(AAABC(276))
      MTR128(5,4)=dconjg(AAABC(292))
      MTR128(5,5)=2D0*dconjg(AAABC(304))
      MTR129(1,1,1)=AAABC(229)
      MTR129(1,1,2)=AAABC(231)
      MTR129(1,2,1)=AAABC(234)
      MTR129(1,2,2)=AAABC(226)
      MTR129(1,3,1)=AAABC(228)
      MTR129(1,3,2)=AAABC(232)
      MTR129(1,4,1)=AAABC(230)
      MTR129(1,4,2)=AAABC(233)
      MTR129(1,5,1)=AAABC(225)
      MTR129(1,5,2)=AAABC(227)
      MTR129(2,1,1)=AAABC(234)
      MTR129(2,1,2)=AAABC(226)
      MTR129(2,2,1)=AAABC(258)
      MTR129(2,2,2)=AAABC(257)
      MTR129(2,3,1)=AAABC(254)
      MTR129(2,3,2)=AAABC(259)
      MTR129(2,4,1)=AAABC(261)
      MTR129(2,4,2)=AAABC(255)
      MTR129(2,5,1)=AAABC(256)
      MTR129(2,5,2)=AAABC(260)
      MTR129(3,1,1)=AAABC(228)
      MTR129(3,1,2)=AAABC(232)
      MTR129(3,2,1)=AAABC(254)
      MTR129(3,2,2)=AAABC(259)
      MTR129(3,3,1)=AAABC(278)
      MTR129(3,3,2)=AAABC(280)
      MTR129(3,4,1)=AAABC(279)
      MTR129(3,4,2)=AAABC(282)
      MTR129(3,5,1)=AAABC(277)
      MTR129(3,5,2)=AAABC(281)
      MTR129(4,1,1)=AAABC(230)
      MTR129(4,1,2)=AAABC(233)
      MTR129(4,2,1)=AAABC(261)
      MTR129(4,2,2)=AAABC(255)
      MTR129(4,3,1)=AAABC(279)
      MTR129(4,3,2)=AAABC(282)
      MTR129(4,4,1)=AAABC(297)
      MTR129(4,4,2)=AAABC(296)
      MTR129(4,5,1)=AAABC(295)
      MTR129(4,5,2)=AAABC(294)
      MTR129(5,1,1)=AAABC(225)
      MTR129(5,1,2)=AAABC(227)
      MTR129(5,2,1)=AAABC(256)
      MTR129(5,2,2)=AAABC(260)
      MTR129(5,3,1)=AAABC(277)
      MTR129(5,3,2)=AAABC(281)
      MTR129(5,4,1)=AAABC(295)
      MTR129(5,4,2)=AAABC(294)
      MTR129(5,5,1)=AAABC(305)
      MTR129(5,5,2)=AAABC(306)
      MTR130(1,1,1)=dconjg(AAABC(229))
      MTR130(1,1,2)=dconjg(AAABC(231))
      MTR130(1,2,1)=dconjg(AAABC(234))
      MTR130(1,2,2)=dconjg(AAABC(226))
      MTR130(1,3,1)=dconjg(AAABC(228))
      MTR130(1,3,2)=dconjg(AAABC(232))
      MTR130(1,4,1)=dconjg(AAABC(230))
      MTR130(1,4,2)=dconjg(AAABC(233))
      MTR130(1,5,1)=dconjg(AAABC(225))
      MTR130(1,5,2)=dconjg(AAABC(227))
      MTR130(2,1,1)=dconjg(AAABC(234))
      MTR130(2,1,2)=dconjg(AAABC(226))
      MTR130(2,2,1)=dconjg(AAABC(258))
      MTR130(2,2,2)=dconjg(AAABC(257))
      MTR130(2,3,1)=dconjg(AAABC(254))
      MTR130(2,3,2)=dconjg(AAABC(259))
      MTR130(2,4,1)=dconjg(AAABC(261))
      MTR130(2,4,2)=dconjg(AAABC(255))
      MTR130(2,5,1)=dconjg(AAABC(256))
      MTR130(2,5,2)=dconjg(AAABC(260))
      MTR130(3,1,1)=dconjg(AAABC(228))
      MTR130(3,1,2)=dconjg(AAABC(232))
      MTR130(3,2,1)=dconjg(AAABC(254))
      MTR130(3,2,2)=dconjg(AAABC(259))
      MTR130(3,3,1)=dconjg(AAABC(278))
      MTR130(3,3,2)=dconjg(AAABC(280))
      MTR130(3,4,1)=dconjg(AAABC(279))
      MTR130(3,4,2)=dconjg(AAABC(282))
      MTR130(3,5,1)=dconjg(AAABC(277))
      MTR130(3,5,2)=dconjg(AAABC(281))
      MTR130(4,1,1)=dconjg(AAABC(230))
      MTR130(4,1,2)=dconjg(AAABC(233))
      MTR130(4,2,1)=dconjg(AAABC(261))
      MTR130(4,2,2)=dconjg(AAABC(255))
      MTR130(4,3,1)=dconjg(AAABC(279))
      MTR130(4,3,2)=dconjg(AAABC(282))
      MTR130(4,4,1)=dconjg(AAABC(297))
      MTR130(4,4,2)=dconjg(AAABC(296))
      MTR130(4,5,1)=dconjg(AAABC(295))
      MTR130(4,5,2)=dconjg(AAABC(294))
      MTR130(5,1,1)=dconjg(AAABC(225))
      MTR130(5,1,2)=dconjg(AAABC(227))
      MTR130(5,2,1)=dconjg(AAABC(256))
      MTR130(5,2,2)=dconjg(AAABC(260))
      MTR130(5,3,1)=dconjg(AAABC(277))
      MTR130(5,3,2)=dconjg(AAABC(281))
      MTR130(5,4,1)=dconjg(AAABC(295))
      MTR130(5,4,2)=dconjg(AAABC(294))
      MTR130(5,5,1)=dconjg(AAABC(305))
      MTR130(5,5,2)=dconjg(AAABC(306))
      MTR131(1,1,1)=AAABC(247)
      MTR131(1,1,2)=AAABC(249)
      MTR131(1,1,3)=AAABC(244)
      MTR131(1,2,1)=AAABC(242)
      MTR131(1,2,2)=AAABC(240)
      MTR131(1,2,3)=AAABC(245)
      MTR131(1,3,1)=AAABC(237)
      MTR131(1,3,2)=AAABC(246)
      MTR131(1,3,3)=AAABC(236)
      MTR131(1,4,1)=AAABC(248)
      MTR131(1,4,2)=AAABC(239)
      MTR131(1,4,3)=AAABC(241)
      MTR131(1,5,1)=AAABC(243)
      MTR131(1,5,2)=AAABC(238)
      MTR131(1,5,3)=AAABC(235)
      MTR131(2,1,1)=AAABC(242)
      MTR131(2,1,2)=AAABC(240)
      MTR131(2,1,3)=AAABC(245)
      MTR131(2,2,1)=AAABC(271)
      MTR131(2,2,2)=AAABC(270)
      MTR131(2,2,3)=AAABC(269)
      MTR131(2,3,1)=AAABC(268)
      MTR131(2,3,2)=AAABC(263)
      MTR131(2,3,3)=AAABC(272)
      MTR131(2,4,1)=AAABC(266)
      MTR131(2,4,2)=AAABC(267)
      MTR131(2,4,3)=AAABC(262)
      MTR131(2,5,1)=AAABC(273)
      MTR131(2,5,2)=AAABC(265)
      MTR131(2,5,3)=AAABC(264)
      MTR131(3,1,1)=AAABC(237)
      MTR131(3,1,2)=AAABC(246)
      MTR131(3,1,3)=AAABC(236)
      MTR131(3,2,1)=AAABC(268)
      MTR131(3,2,2)=AAABC(263)
      MTR131(3,2,3)=AAABC(272)
      MTR131(3,3,1)=AAABC(285)
      MTR131(3,3,2)=AAABC(283)
      MTR131(3,3,3)=AAABC(290)
      MTR131(3,4,1)=AAABC(289)
      MTR131(3,4,2)=AAABC(291)
      MTR131(3,4,3)=AAABC(287)
      MTR131(3,5,1)=AAABC(286)
      MTR131(3,5,2)=AAABC(284)
      MTR131(3,5,3)=AAABC(288)
      MTR131(4,1,1)=AAABC(248)
      MTR131(4,1,2)=AAABC(239)
      MTR131(4,1,3)=AAABC(241)
      MTR131(4,2,1)=AAABC(266)
      MTR131(4,2,2)=AAABC(267)
      MTR131(4,2,3)=AAABC(262)
      MTR131(4,3,1)=AAABC(289)
      MTR131(4,3,2)=AAABC(291)
      MTR131(4,3,3)=AAABC(287)
      MTR131(4,4,1)=AAABC(302)
      MTR131(4,4,2)=AAABC(298)
      MTR131(4,4,3)=AAABC(303)
      MTR131(4,5,1)=AAABC(301)
      MTR131(4,5,2)=AAABC(300)
      MTR131(4,5,3)=AAABC(299)
      MTR131(5,1,1)=AAABC(243)
      MTR131(5,1,2)=AAABC(238)
      MTR131(5,1,3)=AAABC(235)
      MTR131(5,2,1)=AAABC(273)
      MTR131(5,2,2)=AAABC(265)
      MTR131(5,2,3)=AAABC(264)
      MTR131(5,3,1)=AAABC(286)
      MTR131(5,3,2)=AAABC(284)
      MTR131(5,3,3)=AAABC(288)
      MTR131(5,4,1)=AAABC(301)
      MTR131(5,4,2)=AAABC(300)
      MTR131(5,4,3)=AAABC(299)
      MTR131(5,5,1)=AAABC(308)
      MTR131(5,5,2)=AAABC(309)
      MTR131(5,5,3)=AAABC(307)
      MTR132(1,1,1)=dconjg(AAABC(247))
      MTR132(1,1,2)=dconjg(AAABC(249))
      MTR132(1,1,3)=dconjg(AAABC(244))
      MTR132(1,2,1)=dconjg(AAABC(242))
      MTR132(1,2,2)=dconjg(AAABC(240))
      MTR132(1,2,3)=dconjg(AAABC(245))
      MTR132(1,3,1)=dconjg(AAABC(237))
      MTR132(1,3,2)=dconjg(AAABC(246))
      MTR132(1,3,3)=dconjg(AAABC(236))
      MTR132(1,4,1)=dconjg(AAABC(248))
      MTR132(1,4,2)=dconjg(AAABC(239))
      MTR132(1,4,3)=dconjg(AAABC(241))
      MTR132(1,5,1)=dconjg(AAABC(243))
      MTR132(1,5,2)=dconjg(AAABC(238))
      MTR132(1,5,3)=dconjg(AAABC(235))
      MTR132(2,1,1)=dconjg(AAABC(242))
      MTR132(2,1,2)=dconjg(AAABC(240))
      MTR132(2,1,3)=dconjg(AAABC(245))
      MTR132(2,2,1)=dconjg(AAABC(271))
      MTR132(2,2,2)=dconjg(AAABC(270))
      MTR132(2,2,3)=dconjg(AAABC(269))
      MTR132(2,3,1)=dconjg(AAABC(268))
      MTR132(2,3,2)=dconjg(AAABC(263))
      MTR132(2,3,3)=dconjg(AAABC(272))
      MTR132(2,4,1)=dconjg(AAABC(266))
      MTR132(2,4,2)=dconjg(AAABC(267))
      MTR132(2,4,3)=dconjg(AAABC(262))
      MTR132(2,5,1)=dconjg(AAABC(273))
      MTR132(2,5,2)=dconjg(AAABC(265))
      MTR132(2,5,3)=dconjg(AAABC(264))
      MTR132(3,1,1)=dconjg(AAABC(237))
      MTR132(3,1,2)=dconjg(AAABC(246))
      MTR132(3,1,3)=dconjg(AAABC(236))
      MTR132(3,2,1)=dconjg(AAABC(268))
      MTR132(3,2,2)=dconjg(AAABC(263))
      MTR132(3,2,3)=dconjg(AAABC(272))
      MTR132(3,3,1)=dconjg(AAABC(285))
      MTR132(3,3,2)=dconjg(AAABC(283))
      MTR132(3,3,3)=dconjg(AAABC(290))
      MTR132(3,4,1)=dconjg(AAABC(289))
      MTR132(3,4,2)=dconjg(AAABC(291))
      MTR132(3,4,3)=dconjg(AAABC(287))
      MTR132(3,5,1)=dconjg(AAABC(286))
      MTR132(3,5,2)=dconjg(AAABC(284))
      MTR132(3,5,3)=dconjg(AAABC(288))
      MTR132(4,1,1)=dconjg(AAABC(248))
      MTR132(4,1,2)=dconjg(AAABC(239))
      MTR132(4,1,3)=dconjg(AAABC(241))
      MTR132(4,2,1)=dconjg(AAABC(266))
      MTR132(4,2,2)=dconjg(AAABC(267))
      MTR132(4,2,3)=dconjg(AAABC(262))
      MTR132(4,3,1)=dconjg(AAABC(289))
      MTR132(4,3,2)=dconjg(AAABC(291))
      MTR132(4,3,3)=dconjg(AAABC(287))
      MTR132(4,4,1)=dconjg(AAABC(302))
      MTR132(4,4,2)=dconjg(AAABC(298))
      MTR132(4,4,3)=dconjg(AAABC(303))
      MTR132(4,5,1)=dconjg(AAABC(301))
      MTR132(4,5,2)=dconjg(AAABC(300))
      MTR132(4,5,3)=dconjg(AAABC(299))
      MTR132(5,1,1)=dconjg(AAABC(243))
      MTR132(5,1,2)=dconjg(AAABC(238))
      MTR132(5,1,3)=dconjg(AAABC(235))
      MTR132(5,2,1)=dconjg(AAABC(273))
      MTR132(5,2,2)=dconjg(AAABC(265))
      MTR132(5,2,3)=dconjg(AAABC(264))
      MTR132(5,3,1)=dconjg(AAABC(286))
      MTR132(5,3,2)=dconjg(AAABC(284))
      MTR132(5,3,3)=dconjg(AAABC(288))
      MTR132(5,4,1)=dconjg(AAABC(301))
      MTR132(5,4,2)=dconjg(AAABC(300))
      MTR132(5,4,3)=dconjg(AAABC(299))
      MTR132(5,5,1)=dconjg(AAABC(308))
      MTR132(5,5,2)=dconjg(AAABC(309))
      MTR132(5,5,3)=dconjg(AAABC(307))
      MTR133(1,1)=AAABR(366)
      MTR133(1,2)=AAABR(364)
      MTR133(2,1)=AAABR(364)
      MTR133(2,2)=AAABR(368)
      MTR134(1,1)=AAABR(367)
      MTR134(1,2)=AAABR(365)
      MTR134(2,1)=AAABR(365)
      MTR134(2,2)=AAABR(369)
      MTR135(1,1)=AAABC(328)
      MTR135(1,2)=AAABC(318)
      MTR135(1,3)=AAABC(326)
      MTR135(1,4)=AAABC(322)
      MTR135(1,5)=AAABC(312)
      MTR135(2,1)=AAABC(316)
      MTR135(2,2)=AAABC(320)
      MTR135(2,3)=AAABC(314)
      MTR135(2,4)=AAABC(310)
      MTR135(2,5)=AAABC(324)
      MTR136(1,1)=AAABC(329)
      MTR136(1,2)=AAABC(319)
      MTR136(1,3)=AAABC(327)
      MTR136(1,4)=AAABC(323)
      MTR136(1,5)=AAABC(313)
      MTR136(2,1)=AAABC(317)
      MTR136(2,2)=AAABC(321)
      MTR136(2,3)=AAABC(315)
      MTR136(2,4)=AAABC(311)
      MTR136(2,5)=AAABC(325)
      MTR137(1,1)=dconjg(AAABC(328))
      MTR137(1,2)=dconjg(AAABC(318))
      MTR137(1,3)=dconjg(AAABC(326))
      MTR137(1,4)=dconjg(AAABC(322))
      MTR137(1,5)=dconjg(AAABC(312))
      MTR137(2,1)=dconjg(AAABC(316))
      MTR137(2,2)=dconjg(AAABC(320))
      MTR137(2,3)=dconjg(AAABC(314))
      MTR137(2,4)=dconjg(AAABC(310))
      MTR137(2,5)=dconjg(AAABC(324))
      MTR138(1,1)=dconjg(AAABC(329))
      MTR138(1,2)=dconjg(AAABC(319))
      MTR138(1,3)=dconjg(AAABC(327))
      MTR138(1,4)=dconjg(AAABC(323))
      MTR138(1,5)=dconjg(AAABC(313))
      MTR138(2,1)=dconjg(AAABC(317))
      MTR138(2,2)=dconjg(AAABC(321))
      MTR138(2,3)=dconjg(AAABC(315))
      MTR138(2,4)=dconjg(AAABC(311))
      MTR138(2,5)=dconjg(AAABC(325))
      MTR139(1,1)=AAABC(334)
      MTR139(1,2)=dconjg(AAABC(330))
      MTR139(1,3)=dconjg(AAABC(331))
      MTR139(1,4)=dconjg(AAABC(333))
      MTR139(1,5)=dconjg(AAABC(332))
      MTR139(2,1)=AAABC(330)
      MTR139(2,2)=AAABC(336)
      MTR139(2,3)=dconjg(AAABC(338))
      MTR139(2,4)=dconjg(AAABC(337))
      MTR139(2,5)=dconjg(AAABC(335))
      MTR139(3,1)=AAABC(331)
      MTR139(3,2)=AAABC(338)
      MTR139(3,3)=AAABC(341)
      MTR139(3,4)=dconjg(AAABC(339))
      MTR139(3,5)=dconjg(AAABC(340))
      MTR139(4,1)=AAABC(333)
      MTR139(4,2)=AAABC(337)
      MTR139(4,3)=AAABC(339)
      MTR139(4,4)=AAABC(342)
      MTR139(4,5)=dconjg(AAABC(343))
      MTR139(5,1)=AAABC(332)
      MTR139(5,2)=AAABC(335)
      MTR139(5,3)=AAABC(340)
      MTR139(5,4)=AAABC(343)
      MTR139(5,5)=AAABC(344)
      MTR140(1,1)=AAABC(334)
      MTR140(1,2)=AAABC(330)
      MTR140(1,3)=AAABC(331)
      MTR140(1,4)=AAABC(333)
      MTR140(1,5)=AAABC(332)
      MTR140(2,1)=dconjg(AAABC(330))
      MTR140(2,2)=AAABC(336)
      MTR140(2,3)=AAABC(338)
      MTR140(2,4)=AAABC(337)
      MTR140(2,5)=AAABC(335)
      MTR140(3,1)=dconjg(AAABC(331))
      MTR140(3,2)=dconjg(AAABC(338))
      MTR140(3,3)=AAABC(341)
      MTR140(3,4)=AAABC(339)
      MTR140(3,5)=AAABC(340)
      MTR140(4,1)=dconjg(AAABC(333))
      MTR140(4,2)=dconjg(AAABC(337))
      MTR140(4,3)=dconjg(AAABC(339))
      MTR140(4,4)=AAABC(342)
      MTR140(4,5)=AAABC(343)
      MTR140(5,1)=dconjg(AAABC(332))
      MTR140(5,2)=dconjg(AAABC(335))
      MTR140(5,3)=dconjg(AAABC(340))
      MTR140(5,4)=dconjg(AAABC(343))
      MTR140(5,5)=AAABC(344)
      MTR141(1)=AAABR(379)
      MTR141(2)=AAABR(377)
      MTR142(1)=AAABR(378)
      MTR142(2)=AAABR(376)
      MTR143(1)=AAABR(383)
      MTR143(2)=AAABR(381)
      MTR143(3)=AAABR(385)
      MTR144(1)=AAABR(382)
      MTR144(2)=AAABR(380)
      MTR144(3)=AAABR(384)
      MTR145(1)=AAABR(372)
      MTR145(2)=AAABR(372)
      MTR145(3)=AAABR(386)
      MTR146(1)=AAABR(375)
      MTR146(2)=AAABR(375)
      MTR146(3)=AAABR(388)
      MTR147(1)=AAABR(374)
      MTR147(2)=AAABR(374)
      MTR147(3)=AAABR(387)
      MTR148(1)=AAABR(395)
      MTR148(2)=AAABR(395)
      MTR148(3)=AAABR(393)
      MTR149(1)=AAABR(394)
      MTR149(2)=AAABR(394)
      MTR149(3)=AAABR(392)
      MTR150(1)=2D0*AAABR(396)
      MTR150(2)=2D0*AAABR(396)
      MTR150(3)=AAABR(397)
      MTR151(1)=4D0*AAABR(398)
      MTR151(2)=4D0*AAABR(398)
      MTR151(3)=-AAABR(399)
      MTR152(1)=AAABR(406)
      MTR152(2)=AAABR(404)
      MTR152(3)=AAABR(405)
      MTR153(1)=AAABR(410)
      MTR153(2)=AAABR(409)
      MTR153(3)=AAABR(407)
      MTR154(1)=AAABR(415)
      MTR154(2)=AAABR(417)
      MTR154(3)=AAABR(416)
      MTR155(1)=AAABR(411)
      MTR155(2)=AAABR(412)
      MTR155(3)=AAABR(413)
      MTR156(1)=AAABR(420)
      MTR156(2)=AAABR(418)
      MTR156(3)=AAABR(419)
      MTR157(1)=AAABR(421)
      MTR157(2)=AAABR(423)
      MTR157(3)=AAABR(422)
      MTR158(1)=AAABR(428)
      MTR158(2)=AAABR(429)
      MTR158(3)=AAABR(427)
      MTR159(1)=AAABR(425)
      MTR159(2)=AAABR(426)
      MTR159(3)=AAABR(424)
      MTR160(1)=AAABR(433)
      MTR160(2)=AAABR(431)
      MTR160(3)=AAABR(435)
      MTR161(1)=AAABR(437)
      MTR161(2)=AAABR(436)
      MTR161(3)=AAABR(438)
      MTR162(1)=AAABR(432)
      MTR162(2)=AAABR(430)
      MTR162(3)=AAABR(434)
      MTR163(1)=AAABR(440)
      MTR163(2)=AAABR(442)
      MTR163(3)=AAABR(444)
      MTR164(1)=AAABR(445)
      MTR164(2)=AAABR(447)
      MTR164(3)=AAABR(446)
      MTR165(1)=AAABR(439)
      MTR165(2)=AAABR(441)
      MTR165(3)=AAABR(443)
      MTR166(1)=AAABR(448)
      MTR166(2)=AAABR(450)
      MTR166(3)=AAABR(449)
      MTR167(1)=AAABR(452)
      MTR167(2)=AAABR(451)
      MTR167(3)=AAABR(453)
      MTR168(1)=AAABR(457)
      MTR168(2)=AAABR(459)
      MTR168(3)=AAABR(458)
      MTR169(1)=AAABR(454)
      MTR169(2)=AAABR(456)
      MTR169(3)=AAABR(455)
      MTR170(1)=AAABR(465)
      MTR170(2)=AAABR(461)
      MTR170(3)=AAABR(463)
      MTR171(1)=AAABR(466)
      MTR171(2)=AAABR(468)
      MTR171(3)=AAABR(467)
      MTR172(1)=AAABR(464)
      MTR172(2)=AAABR(460)
      MTR172(3)=AAABR(462)
      MTR173(1)=AAABR(470)
      MTR173(2)=AAABR(472)
      MTR173(3)=AAABR(474)
      MTR174(1)=AAABR(477)
      MTR174(2)=AAABR(475)
      MTR174(3)=AAABR(476)
      MTR175(1)=AAABR(469)
      MTR175(2)=AAABR(471)
      MTR175(3)=AAABR(473)
      MTR176(1)=AAABR(480)
      MTR176(2)=AAABR(479)
      MTR176(3)=AAABR(478)
      MTR177(1)=AAABR(481)
      MTR177(2)=AAABR(483)
      MTR177(3)=AAABR(482)
      MTR178(1)=AAABR(488)
      MTR178(2)=AAABR(489)
      MTR178(3)=AAABR(487)
      MTR179(1)=AAABR(486)
      MTR179(2)=AAABR(484)
      MTR179(3)=AAABR(485)
      MTR180(1)=AAABR(491)
      MTR180(2)=AAABR(492)
      MTR180(3)=AAABR(490)
      MTR181(1)=AAABR(495)
      MTR181(2)=AAABR(493)
      MTR181(3)=AAABR(494)
      MTR182(1)=AAABR(501)
      MTR182(2)=AAABR(499)
      MTR182(3)=AAABR(500)
      MTR183(1)=AAABR(497)
      MTR183(2)=AAABR(496)
      MTR183(3)=AAABR(498)
      MTR184(1)=AAABR(504)
      MTR184(2)=AAABR(502)
      MTR184(3)=AAABR(503)
      MTR185(1)=AAABR(507)
      MTR185(2)=AAABR(505)
      MTR185(3)=AAABR(506)
      MTR186(1)=AAABR(512)
      MTR186(2)=AAABR(513)
      MTR186(3)=AAABR(511)
      MTR187(1)=AAABR(509)
      MTR187(2)=AAABR(510)
      MTR187(3)=AAABR(508)
      MTR188(1)=AAABR(519)
      MTR188(2)=AAABR(521)
      MTR188(3)=AAABR(517)
      MTR189(1)=AAABR(518)
      MTR189(2)=AAABR(520)
      MTR189(3)=AAABR(516)
      MTR190(1,1)=AAABR(525)
      MTR190(1,2)=AAABR(526)
      MTR190(2,1)=MTR190(1,2)
      MTR190(2,2)=AAABR(527)
      MTR191(1,1)=AAABR(522)
      MTR191(1,2)=AAABR(523)
      MTR191(2,1)=MTR191(1,2)
      MTR191(2,2)=AAABR(524)
      MTR192(1,1)=AAABR(537)
      MTR192(1,2)=AAABR(536)
      MTR192(1,3)=AAABR(535)
      MTR192(2,1)=MTR192(1,2)
      MTR192(2,2)=AAABR(534)
      MTR192(2,3)=AAABR(539)
      MTR192(3,1)=MTR192(1,3)
      MTR192(3,2)=MTR192(2,3)
      MTR192(3,3)=AAABR(538)
      MTR193(1,1)=AAABR(533)
      MTR193(1,2)=AAABR(529)
      MTR193(1,3)=AAABR(531)
      MTR193(2,1)=MTR193(1,2)
      MTR193(2,2)=AAABR(532)
      MTR193(2,3)=AAABR(528)
      MTR193(3,1)=MTR193(1,3)
      MTR193(3,2)=MTR193(2,3)
      MTR193(3,3)=AAABR(530)

      end

***********************************************

      subroutine ModelVarIni(fail, sqrtS)
      implicit none
      double precision sqrtS
      integer fail
      double precision Alfas

#include "model.h"

c      double precision ALPHAS2
c      external ALPHAS2

c      Alfas = ALPHAS2(sqrtS)
c      GG = sqrt(4*pi*Alfas)
      fail=0
      end

************************************************

      subroutine ModelDigest
      implicit none

#include "model.h"

      write(16,*) 'EE=',EE
      write(16,*) 'MZ=',MZ
      write(16,*) 'MW=',MW
      write(16,*) 'GG=',GG
      write(16,*) 'alfSMZ=',alfSMZ
      write(16,*) 'MbMb=',MbMb
      write(16,*) 'McMc=',McMc
      write(16,*) 'Mtp=',Mtp
      write(16,*) 'Ml=',Ml
      write(16,*) 'Mq=',Mq
      write(16,*) 'CW=',CW
      write(16,*) 'SW=',SW
      write(16,*) 'C2W=',C2W
      write(16,*) 'LamQCD=',LamQCD
      write(16,*) 'rd=',rd
      write(16,*) 'mu=',mu
      write(16,*) 'Lambda=',Lambda
      write(16,*) 'Kappa=',Kappa
      write(16,*) 'aLambda=',aLambda
      write(16,*) 'aKappa=',aKappa
      write(16,*) 'muP=',muP
      write(16,*) 'At=',At
      write(16,*) 'Ab=',Ab
      write(16,*) 'Al=',Al
      write(16,*) 'Ml2=',Ml2
      write(16,*) 'Ml3=',Ml3
      write(16,*) 'Mr2=',Mr2
      write(16,*) 'Mr3=',Mr3
      write(16,*) 'Mq2=',Mq2
      write(16,*) 'Mq3=',Mq3
      write(16,*) 'Mu2=',Mu2
      write(16,*) 'Mu3=',Mu3
      write(16,*) 'Md2=',Md2
      write(16,*) 'Md3=',Md3
      write(16,*) 'tB=',tB
      write(16,*) 'Mh1=',Mh1
      write(16,*) 'Mh2=',Mh2
      write(16,*) 'Mh3=',Mh3
      write(16,*) 'Mha=',Mha
      write(16,*) 'Mhb=',Mhb
      write(16,*) 'MHc=',MHc
      write(16,*) 'MNE1=',MNE1
      write(16,*) 'MNE2=',MNE2
      write(16,*) 'MNE3=',MNE3
      write(16,*) 'MNE4=',MNE4
      write(16,*) 'MNE5=',MNE5
      write(16,*) 'MC1=',MC1
      write(16,*) 'MC2=',MC2
      write(16,*) 'MSG=',MSG
      write(16,*) 'MSuL=',MSuL
      write(16,*) 'MSuR=',MSuR
      write(16,*) 'MSdL=',MSdL
      write(16,*) 'MSdR=',MSdR
      write(16,*) 'MScL=',MScL
      write(16,*) 'MScR=',MScR
      write(16,*) 'MSsL=',MSsL
      write(16,*) 'MSsR=',MSsR
      write(16,*) 'MSt1=',MSt1
      write(16,*) 'MSt2=',MSt2
      write(16,*) 'MSb1=',MSb1
      write(16,*) 'MSb2=',MSb2
      write(16,*) 'QSUSY=',QSUSY
      write(16,*) 'MSeL=',MSeL
      write(16,*) 'MSeR=',MSeR
      write(16,*) 'MSmL=',MSmL
      write(16,*) 'MSmR=',MSmR
      write(16,*) 'MSl1=',MSl1
      write(16,*) 'MSl2=',MSl2
      write(16,*) 'MSne=',MSne
      write(16,*) 'MSnm=',MSnm
      write(16,*) 'MSnl=',MSnl
      write(16,*) 'Zl11=',Zl11
      write(16,*) 'Zl12=',Zl12
      write(16,*) 'Zl21=',Zl21
      write(16,*) 'Zl22=',Zl22
      write(16,*) 'Zh11=',Zh11
      write(16,*) 'Zh12=',Zh12
      write(16,*) 'Zh13=',Zh13
      write(16,*) 'Zh21=',Zh21
      write(16,*) 'Zh22=',Zh22
      write(16,*) 'Zh23=',Zh23
      write(16,*) 'Zh31=',Zh31
      write(16,*) 'Zh32=',Zh32
      write(16,*) 'Zh33=',Zh33
      write(16,*) 'Za11=',Za11
      write(16,*) 'Za12=',Za12
      write(16,*) 'Za13=',Za13
      write(16,*) 'Za21=',Za21
      write(16,*) 'Za22=',Za22
      write(16,*) 'Za23=',Za23
      write(16,*) 'Zn11=',Zn11
      write(16,*) 'Zn12=',Zn12
      write(16,*) 'Zn13=',Zn13
      write(16,*) 'Zn14=',Zn14
      write(16,*) 'Zn15=',Zn15
      write(16,*) 'Zn21=',Zn21
      write(16,*) 'Zn22=',Zn22
      write(16,*) 'Zn23=',Zn23
      write(16,*) 'Zn24=',Zn24
      write(16,*) 'Zn25=',Zn25
      write(16,*) 'Zn31=',Zn31
      write(16,*) 'Zn32=',Zn32
      write(16,*) 'Zn33=',Zn33
      write(16,*) 'Zn34=',Zn34
      write(16,*) 'Zn35=',Zn35
      write(16,*) 'Zn41=',Zn41
      write(16,*) 'Zn42=',Zn42
      write(16,*) 'Zn43=',Zn43
      write(16,*) 'Zn44=',Zn44
      write(16,*) 'Zn45=',Zn45
      write(16,*) 'Zn51=',Zn51
      write(16,*) 'Zn52=',Zn52
      write(16,*) 'Zn53=',Zn53
      write(16,*) 'Zn54=',Zn54
      write(16,*) 'Zn55=',Zn55
      write(16,*) 'Zu11=',Zu11
      write(16,*) 'Zu12=',Zu12
      write(16,*) 'Zu21=',Zu21
      write(16,*) 'Zu22=',Zu22
      write(16,*) 'Zv11=',Zv11
      write(16,*) 'Zv12=',Zv12
      write(16,*) 'Zv21=',Zv21
      write(16,*) 'Zv22=',Zv22
      write(16,*) 'Zt11=',Zt11
      write(16,*) 'Zt12=',Zt12
      write(16,*) 'Zt21=',Zt21
      write(16,*) 'Zt22=',Zt22
      write(16,*) 'Zb11=',Zb11
      write(16,*) 'Zb12=',Zb12
      write(16,*) 'Zb21=',Zb21
      write(16,*) 'Zb22=',Zb22
      write(16,*) 'la1=',la1
      write(16,*) 'la2=',la2
      write(16,*) 'la3=',la3
      write(16,*) 'la4=',la4
      write(16,*) 'la5=',la5
      write(16,*) 'la6=',la6
      write(16,*) 'la7=',la7
      write(16,*) 'la1s=',la1s
      write(16,*) 'la2s=',la2s
      write(16,*) 'la3s=',la3s
      write(16,*) 'la4s=',la4s
      write(16,*) 'la5s=',la5s
      write(16,*) 'la6s=',la6s
      write(16,*) 'la7s=',la7s
      write(16,*) 'la8s=',la8s
      write(16,*) 'aa1=',aa1
      write(16,*) 'aa2=',aa2
      write(16,*) 'aa3=',aa3
      write(16,*) 'aa4=',aa4
      write(16,*) 'aa5=',aa5
      write(16,*) 'aa6=',aa6
      write(16,*) 'B1=',mB1
      write(16,*) 'B2=',mB2
      write(16,*) 'X=',X
      write(16,*) 'dMb=',dMb
      write(16,*) 'x=',xif
      write(16,*) 'xis=',xis
      write(16,*) 'MM3=',MM3
      write(16,*) 'MSP=',MSP
      write(16,*) 'sb=',sb
      write(16,*) 'cb=',cb
      write(16,*) 't2b=',t2b
      write(16,*) 'xvev=',xvev
      write(16,*) 'Pa12=',Pa12
      write(16,*) 'Pa22=',Pa22
      write(16,*) 'Pa11=',Pa11
      write(16,*) 'Pa21=',Pa21
      write(16,*) 'Td3=',Td3
      write(16,*) 'Q=',Q
      write(16,*) 'Mt=',Mt
      write(16,*) 'Mb=',Mb
      write(16,*) 'Mc=',Mc
      write(16,*) 'xH2=',xH2
      write(16,*) 'yH2=',yH2
      write(16,*) 'dMd=',dMd
      write(16,*) 'Td2=',Td2
      write(16,*) 'Au=',Au
      write(16,*) 'Ad=',Ad
      write(16,*) 'fiuu=',fiuu
      write(16,*) 'fidd=',fidd
      write(16,*) 'ficc=',ficc
      write(16,*) 'fiss=',fiss
      write(16,*) 'Zuu11=',Zuu11
      write(16,*) 'Zuu12=',Zuu12
      write(16,*) 'Zuu21=',Zuu21
      write(16,*) 'Zuu22=',Zuu22
      write(16,*) 'Zdd11=',Zdd11
      write(16,*) 'Zdd12=',Zdd12
      write(16,*) 'Zdd21=',Zdd21
      write(16,*) 'Zdd22=',Zdd22
      write(16,*) 'Zcc11=',Zcc11
      write(16,*) 'Zcc12=',Zcc12
      write(16,*) 'Zcc21=',Zcc21
      write(16,*) 'Zcc22=',Zcc22
      write(16,*) 'Zss11=',Zss11
      write(16,*) 'Zss12=',Zss12
      write(16,*) 'Zss21=',Zss21
      write(16,*) 'Zss22=',Zss22
      write(16,*) 'Mele=',Mele
      write(16,*) 'Mmu=',Mmu
      write(16,*) 'MbMM=',MbMM
      write(16,*) 'MtMM=',MtMM
      write(16,*) 'NMM11=',NMM11
      write(16,*) 'NMM12=',NMM12
      write(16,*) 'NMM13=',NMM13
      write(16,*) 'NMM14=',NMM14
      write(16,*) 'NMM15=',NMM15
      write(16,*) 'NMM22=',NMM22
      write(16,*) 'NMM23=',NMM23
      write(16,*) 'NMM24=',NMM24
      write(16,*) 'NMM25=',NMM25
      write(16,*) 'NMM33=',NMM33
      write(16,*) 'NMM34=',NMM34
      write(16,*) 'NMM35=',NMM35
      write(16,*) 'NMM44=',NMM44
      write(16,*) 'NMM45=',NMM45
      write(16,*) 'NMM55=',NMM55
      write(16,*) 'MG1I=',MG1I
      write(16,*) 'MG2I=',MG2I

      end

