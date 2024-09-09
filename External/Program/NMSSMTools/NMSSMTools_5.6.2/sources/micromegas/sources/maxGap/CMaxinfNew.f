      Real Function CMaxinf(CL,Ifmin,mu)
C Uses table CMaxf.in, which was made by CMaxinfgen, to compute
C the CL confidence level value of CMaxinf for minimum f of the
C Ifmin value of fmin and for total expected number of events, mu.
C fmins has the list of values of Nfmin values of fmin.
C As of this writing, correspondence between Ifminin and fmin:
C Ifmin  fmin
C   1    0.00
C   2    0.01
C   3    0.02
C   4    0.05
C   5    0.10
C   6    0.20
C   7    0.50 = fmins(Nfmin)
C
C As of this writing, CLs consist of NCLs=40 values from .8 to .995.
C CL is not allowed to be less than CLs(1) or greater than CLs(NCls).
C If mu is too low, we can expect the result to be unreliable, but it's
C permitted.  If mu > mus(Nmus), give the same result as mu=mus(Nmus)
C or, if extrapolation is available, use A + B/sqrt(mu).  If mu<mus(Nmus),
C linearly interpolate in log(mus),CLs.
C
C Feb 14, 2007: Also read in table CMaxfLow.in, also made with CMaxinfgen.
C It tabulates mu in uniform steps of .1 starting at mu=1.6 and ending at
C mu=54.6.  Use it for mu<=mucut=?54.5, and interpolate between table entries
C N and N+1, where N=10.*(mu-1.5) for N=1 to 530.
C If there is not CMaxfLow.in, set UseTable2 to .false.
C
C Jan 2010: Also read in table CMaxfLowNew.in, made with CMaxinfgenNew.  It
C has NCLs=80 values from .8 to .9975, mu in uniform steps of .0125 starting at
C mu=1.6 and ending at mu=12.  Use it for mu < 12, ifmin=1, interpolate between
C table entries N and N+1, where N=80.*(mu-1.5875) for N=1 to 832.
C If there is not CMaxfLowNew.in, set UseTable3 to .false.
      Implicit None
      Real CL,fmin,mu
      Integer Ifmin,Nmus,NmuMax,NCLs,NCLmax,Nfmin,NfminMax,I,J,
     1 Imu,ICL,Imu24,Imu14,Imu34,Imustart,Nmus2,NmuMax2,NCLmax2,
     2 NCL3s,NmuMax3,Nmus3
      Parameter (NCLmax=40)
      Parameter (NCL3s=80) ! For CMaxfLowNew.in
      Parameter (NfminMax=7)
      Parameter (NmuMax=50)
      Parameter (NmuMax2=730) ! 723 should be needed
      Parameter (NmuMax3=840) ! 833 should be needed
      Real CLs(NCLmax),mus(NmuMax),Table(NCLmax,NfminMax,NmuMax),
     1 fmins(NfminMax),logmu(NmuMax),lmu,xCL,xmu,sqrtmu,
     2 A(NClmax,NfminMax),B(NClmax,NfminMax),mus2(NmuMax2),
     3 Table2(NCLmax,NFminMax,NmuMax2),mucut/54.5/,CL3s(NCL3s),
     4 mucut3/12./,mus3(NmuMax3),Table3(NCL3s,NmuMax3),UseNewCMax
      Logical first/.true./,UseTable2/.true./,UseTable3/.true./,
     1 use
      Common/CMaxinfcom/NCLs,Nfmin,CLs,fmins,mus,Table,Nmus,Imu24,
     1 Imu14,Imu34,logmu,A,B,CL3s,nmus3,mus3,table3

      character*1000 PATH, FILENAME,catpath
      COMMON/maxGapPath/PATH,FILENAME    
      external catpath
           
      If(first) Then
         first=.false.
         FILENAME=catpath(PATH,'CMaxf.in')
         Open(20,file=FILENAME,status='OLD',form='UNFORMATTED')
         Read(20) NCLs,Nfmin
         Read(20) (CLs(I),I=1,NCLs),(fmins(I),I=1,Nfmin)
         Nmus=0
         A(1,1)=0. ! A(1,1)=0 implies extrapolation data is unavailable
 10      Continue
         Nmus=Nmus+1
         read(20,End=100) mus(Nmus)
         If(mus(Nmus).eq.0.) Then
C mus(Nmus)=0 ==> no more Table entries but maybe extrapolation data.
            Read(20,End=100) ((A(ICL,J),B(ICL,J),ICL=1,NCLs),
     1       J=1,Nfmin)
            Go to 100
         EndIf
         logmu(Nmus)=log(mus(Nmus))
         Do J=1,Nfmin
           Read(20) (Table(I,J,Nmus),I=1,NCLs)
         EndDo
         Go to 10
 100     close(20)
         Nmus=Nmus-1
         Imu24=Nmus/2
         Imu14=Imu24/2
         Imu34=3*Nmus/4
C If it exists, read mus2 and Table2 in, assuming they have the same
C CLs and fmins as went with Table.  Table2 is from CMaxfLow.in. 
         If(UseTable2) Then
            FILENAME=catpath(PATH,'CMaxfLow.in')
            Open(20,file=FILENAME,status='OLD',
     1        form='UNFORMATTED')
            Read(20)
            Read(20)
            Nmus2=0
 104        Continue
            Nmus2=Nmus2+1
            read(20,end=107) mus2(Nmus2)
            Do J=1,Nfmin
               Read(20) (Table2(I,J,Nmus2),I=1,NCLs)
            EndDo
            Go to 104
 107        Close(20)
            Nmus2=Nmus2-1
         EndIf
         If(UseTable3) Then
            FILENAME=catpath(PATH,'CMaxflowNew.in')
            Open(20,file=FILENAME,status='OLD',
     1        form='UNFORMATTED')
            Read(20) CL3s
            Nmus3=0
 108        Continue
            Nmus3=Nmus3+1
            Read(20,end=109) mus3(Nmus3),(Table3(I,Nmus3),I=1,NCL3s)
            Go to 108
 109        Close(20)
            Nmus3=Nmus3-1
            Do I=2,112
               Table3(1,I)=0.8
            EndDo
         EndIf
      EndIf
      If(Ifmin.gt.Nfmin .or. Ifmin.lt.1) Then
          Write(6,115) Ifmin
 115      Format(I9,' is outside the allowed range of Ifmin values')
          Stop
      EndIf
      If(UseTable3 .and. mu.le.mucut3) Then
         ICL=(CL-.7975)*400.
         If(ICL.eq.0) ICL=1
         If(ICL.lt.1 .or. ICL.ge.NCL3s) Then
            Write(6,120) CL
            Stop
         EndIf
         xCL=(CL-CL3s(ICL))/(CL3s(ICL+1)-CL3s(ICL)) ! fraction of way to ICL+1
         xmu=80.*mu - 127. ! 80*1.6-127.=1.  Imu=xmu=1 <--> mu=1.6
         Imu=xmu
         xmu=xmu-float(Imu) ! Now xmu is the fraction of the way to Imu+1
C Interpolate Table with CL between CL3s(ICL) and CL3s(ICL+1) and with
C mu between mus3(Imu) and mus3(Imu+1).
          CMaxinf=(1.-xCL)*((1.-xmu)*Table3(ICL,Imu) +
     1   xmu*Table3(ICL,Imu+1)) +
     2   xCL*((1.-xmu)*Table3(ICL+1,Imu) +
     3   xmu*Table3(ICL+1,Imu+1))
         Goto 300
      Else
C Find between what CL's we should interpolate
C Here we assume 40 CL's from .8 to .995.
        ICL=(CL-.795)*200.
C Allow CL as low as .795 with extrapolation.
        If(ICL.eq.0) ICL=1
        If(ICL.lt.1 .or. ICL.ge.NCLs) Then
          Write(6,120) CL
 120      Format(F9.4,' is outside the allowed range of CL values')
          Stop
        EndIF
      EndIf
C Find between what mu values we should interpolate
      If(UseTable2.and.mu.le.mucut) Then
C Initially xmu will be approximately the index of mu in the table
C We will want Imu to be the index of the table entry whose mus is
C just below mu.  We will interpolate between Imu and Imu+1.  Imu will
C be the highest integer below xmu.
         If(mu.lt.6.0) Then
            xmu=40.*mu-63. ! 40*1.6 - 63. = 1.  Imu=xmu=1 <--> mu=1.6
C For mu < 6, xmu is < 40*6-63=177.  Imu=xmu=177 is where mu=6.  So
C this part of the calculation only gets Imu up to 176, and up to Imu+1=177
C is the last of the series where spacing is .025 = 1/40.
            If(xmu.lt.1.) Then
              CMaxinf=0. ! mu is below the lowest entry in the table
              Return
            EndIf
         ElseIf(mu.lt.12) Then
            xmu=20.*mu + 57. ! 20*6+57=177=40*6-63
C For mu < 12, xmu is < 20*12+57=297.  This part of the calculation only
C gets Imu up to 296, and up to Imu+1=297 the spacing is .05 = 1/20. 
         Else
            xmu=10.*mu + 177. ! 10*12+177=297=20*12+57
C For mu < 54.6, xmu is < 10*54.6+177 = 723 = Nmus2, the number of mu entries
C in the entire Table2.   From the above comments, we should have:
C mus2(1)=1.6, mus2(177)=6, mus2(297)=12, mus2(723)=54.6.
         EndIf
         Imu=xmu
C Now xmu is between Imu and Imu+1.
         xmu=xmu-float(Imu) ! Now xmu is the fraction of the way to Imu+1
         GoTo 200
      EndIf
      lmu=log(mu)
C To speed this up, first find what quarter of the range should
C be searched.  Then for simplicity over a complete binary search,
C complete the search in a dumb way.
      If(mus(Imu24).ge.mu) Then
         If(mus(Imu14).ge.mu) Then
            Imustart=2
         Else
            Imustart=Imu14
         EndIf
      Elseif(mus(Imu34).ge.mu) Then
         Imustart=Imu24
      Else
         Imustart=Imu34
      EndIf
      Do I=Imustart,Nmus
          Imu=I
          If(mus(I).ge.mu) Then
             Go to 200
          EndIf
      EndDo
C mu is greater than all the tabulated values; so use the top one or
C use A and B
      lmu=logmu(Nmus)
 200  Continue
      xCL=(CL-CLS(ICL))/(CLS(ICL+1)-CLS(ICL))
      If(mu.gt.mus(Nmus) .and. A(1,1).ne.0.) Then
C Interpolate the computed extrapolated value between CLs(ICL) and
C CLs(ICL+1).
       sqrtmu=sqrt(mu)
       CMaxinf=(1.-xCL)*(A(ICL,Ifmin) + B(ICL,Ifmin)/sqrtmu) +
     1   xCL*(A(ICL+1,Ifmin) + B(ICL+1,Ifmin)/sqrtmu)
      Elseif(UseTable2 .and. mu.le.mucut) Then
C Interpolate Table with CL between CLs(ICL) and CLs(ICL+1) and with
C mu between mus2(Imu) and mus2(Imu+1).
         CMaxinf=(1.-xCL)*((1.-xmu)*Table2(ICL,Ifmin,Imu) +
     1   xmu*Table2(ICL,Ifmin,Imu+1)) +
     2   xCL*((1.-xmu)*Table2(ICL+1,Ifmin,Imu) +
     3   xmu*Table2(ICL+1,Ifmin,Imu+1))
      Else
        Imu=Imu-1
C Interpolate Table with CL between CLs(ICL) and CLs(ICL+1) and with
C lmu between logmu(Imu) and logmu(Imu+1).
        xmu=(lmu-logmu(Imu))/(logmu(Imu+1)-logmu(Imu))
         CMaxinf=(1.-xCL)*((1.-xmu)*Table(ICL,Ifmin,Imu) +
     1   xmu*Table(ICL,Ifmin,Imu+1)) +
     2   xCL*((1.-xmu)*Table(ICL+1,Ifmin,Imu) +
     3   xmu*Table(ICL+1,Ifmin,Imu+1))
      EndIf
 300  CMaxinf=max(CMaxinf,CL)
C The above statement may be needed if mu is too low for the table.
C If the table is made with steps in mu of, say 0.025, then between
C the threshold of where CMaxinf should be defined and where the
C first non-zero entry of the table is there should be a legitimate
C return of CL.
      return
      Entry  UseNewCMax(use)
C Turn on or off use of new, more accurate table.
      UseTable3=use
      first=.true.
      UseNewCMax=0.
      Return
      End
