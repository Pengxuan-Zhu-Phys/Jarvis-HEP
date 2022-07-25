      Real Function UpperLim(CL,If,N,FC,muB,FB,Iflag)
C Calls y_vs_CLf, which calls DGAUSN, and ConfLev, which calls GAMDIS,
C both of CERNLIB (libmathlib and libkernlib).  But in Jan 2010 changed
C to use ConfLevNew, which calls ConfLev if necessary.
C
C Suppose you have a set of N events distributed in some 1-d variable
C and want to know the CL confidence level upper limit on the mean of
C the expected number of events.  Assume there's expected distribution
C characterized by some cumulative probability function, on top of
C which is an unknown background.  UpperLim is the optimum interval
C upper limit, taking into account deviation between the observed
C distribution and the predicted one.
C
C CL is the confidence level desired
C If says which minimum fraction of the cumulative probability
C   is allowed for seeking the optimum interval.  If=1,2,3,4,5,6,7
C   corresponds to minimum cumulative probability interval =
C   .00, .01, .02, .05, .10, .20, .50.
C N is the number of events
C FC: Given the foreground distribution whose shape is known, but whose
C    normalization is to have its upper limit total expected number of
C    events determined, FC(0) to FC(N+1), with FC(0)=0, FC(N+1)=1, and with 
C    FC(i) the increasing ordered set of cumulative probabilities for the
C    foreground distribution for event i, i=1 to N.
C muB is the total expected number of events from known background.
C FB is like FC but assuming the distribution shape from known background.
C Iflag is the return code flag.:
C   0: Normal return
C  The flag bits correspond to:
C   1: More than 5 iterations required.  Unusual, but should be ok.
C   2: More than 10 iterations needed, but not obtained.  This may be serious.
C    Other bits in Iflag refer to what happened in the last iteration
C   4: y_vs_CLf returned status 1 at some time (extrapolated from f0=0.01)
C   8: y_vs_CLf returned status 2 at some time (extrapolated from f0=1)
C  16: The optimum interval had status 1.
C  32: The optimum interval had status 2.
C  64: Failure to solve CMax = CMaxbar.
C 128: Couldn't solve CMax=CMaxbar because upperlim wants to be <0.
C 256: More than one solution.  The other two are in Exclude_low, which gives
C      the excluded range below the absolute upper limit.
C If something goes wrong which prevents return of correct results, the
C program prints an error message and stops.
C  There is also a Common/UpperLimCom/EndPoints(2),Exclude_low(2).
C  EndPoints(1) is an integer giving the I of FC(I) at which the optimum
C  interval started and EndPoints(2) is the I at which the optimum interval
C  ended.  When Iflag has bit 256 on, everything above UpperLim is excluded,
C  as is everything between Exclude_low(1) and Exclude_low(2).  But between
C  Exclude_low(2) and UpperLim is allowed.
      Implicit None
      Integer N,If,Iflag,NMax,I,m,Niter,Istat,IflagOpt,MaxF,
     1 NMax1,N1,If1,NCalls,I1,I2,mmax,EndPoints
      Parameter (NMax=150000)
      Real MeanMax
C The number of iteration, Niter, can be 5 and it only occasionally needs
C more for low mu and with lots of background.  Even then, 5 is enough to
C almost always get almost exactly the same answer.
      Parameter (Niter=10)
      Real CL,FC(0:N+1),y_vs_CLf,CMaxinf,mu,mu0,f(0:NMax),y,y2,x,
     1 CMax,fmin(7)/.00,.01,.02,.05,.1,.2,.5/,eps/.0001/,Cdbg,
     2 fdebug,ydebug,mutmp,R,Topmu,Botmu,CL1,muB,FB(0:N+1),FC1,FB1,muB1,
     3 FUp,fmin1(7),mucand(3),TopNow,BotNow,FUPbot,FUPnow,Exclude_low,
     4 fhigh,fdiff
      Common/Fupcom/f,N1,CL1,If1,MeanMax,NCalls,Istat,FC1(0:NMax),
     1 FB1(0:NMax),muB1,fmin1,Cdbg(0:10),mmax
      Common/UpperLimcom/EndPoints(2),Exclude_low(2)
      Logical debug/.false./
      External FUp
      If(N.ge.NMax) Then
         Write(6,*) N,", the number of events, is above ",NMax-1
         Stop
      EndIf
      Do I=1,3
         mucand(I)=-1.
      EndDo
      MeanMax=54.5
      FC(0)=0.
      FC(N+1)=1.
      muB1=muB
      Do I=1,7
         fmin1(I)=fmin(I)
      EndDo
      If(muB.ne.0.) Then
         FB(0)=0.
         FB(N+1)=1.
         Do I=0,N+1
            FC1(I)=FC(I)
            FB1(I)=FB(I)
         EndDo
      EndIf
C For each m=0 to N find f(m), the maximum over all I from 0 to N-m of
C FC(I+m+1)-FC(I).  Start out with mu0=float(N) and evaluate
C CMax=CMaxinf(CL,If,mu).  For each m, y=y_vs_CLf(CMax,f(m)), and x from
C y=(m-x)/sqrt(x).  Find the smallest value of x/f(m) and call
C it the new mu.  Iterate until the fractional change of mu < eps, at which
C time take UpperLim=mu.
      If(CL.lt.0.8 .or. CL.gt. 0.995) Then
         Write(6,*) CL,
     1    " is out of the permissible confidence level range."
         Stop
      EndIf
      If(If.lt.1 .or. If.gt.7) Then
         Write(6,*) If,
     1    " is out of the permissible value selecting fmin."
         Stop
      EndIf
      Iflag=0
      If(N.eq.0) Then
         UpperLim=log(1./(1.-CL))
         EndPoints(1)=0
         EndPoints(2)=0
         Return
      EndIf
      If(muB.eq.0.) Then
       Do m=0,N
         f(m)=0.
         Do I=0,N-m
            f(m)=Max(f(m),FC(I+m+1)-FC(I))
         EndDo
       EndDo
      EndIf
C For some reason, the quick method of convergence sometimes fails
C with muB>0.
      I1=1
      If(muB.ne.0.) Go to 50
      mu0=Float(N)
      Do I=1,Niter
         If(mu0.lt.MeanMax) GoTo 50
         CMax=CMaxinf(CL,If,mu0)
         mu=1.E10
         Do m=0,N
            If(muB.ne.0.) Then
               f(m)=0.
               Do I2=0,N-m
                  f(m)=Max(f(m),(1.-(muB/mu0))*(FC(I2+m+1)-FC(I2))+
     1                (muB/mu0)*(FB(I2+m+1)-FB(I2)))
               EndDo
            EndIf
            If(f(m).gt.fmin(If)) Then
               y=y_vs_CLf(CMax,f(m),Istat)
               If(Istat.gt.2) Then
                  Write(6,*) "y_vs_CLf returned with status",Istat
                  Go to 50
               EndIf
               Iflag=Or(Iflag,4*Istat)
               y2=y*y
               x=float(m)+ .5*(y2+sqrt(y2*(4.*float(m)+y2)))
               mutmp=x/f(m)
               If(mutmp.lt.mu) Then
                  mu=mutmp
C IflagOpt will have the status of the optimum interval
                  IflagOpt=Istat
                  mmax=m
                  fdebug=f(m)
                  ydebug=y
               EndIf
             EndIf
         EndDo
         If(abs(mu0-mu)/mu .lt. eps .and. I.gt.1) Go to 100
         mu0=max(muB,mu)
         If(mu0.lt.muB) Go to 50
         If(I.eq.5) Then
            Iflag=Iflag+1
            GoTo 50 ! It looks like it won't converge (often with muB>0)
         EndIf
      EndDo
      Iflag=Iflag+2
      Write(6,*) "UpperLim did the maximum number of iterations,",
     1 Niter
 50   Continue
C Come here if it's starting to look like mu<54.5, or if convergence
C fails for mu>54.5.
      IflagOpt=0
      MAXF=500
      N1=N
      CL1=CL
      If1=If
      Topmu=float(N)+4.*sqrt(Float(N))+5.
      Botmu=max(muB,log(1./(1.-CL)))
      If( muB.gt.0. .and. FUp(Botmu,1) .gt.0.) Then
C It looks like UpperLim wants to be negative
         UpperLim=0.
         Iflag=IFlag+192
         Return
      EndIf
      Call RZERO(Botmu,Topmu,mu,R,EPS,MAXF,FUp)
      If(R.lt.0. .or. Istat.gt.4) Iflag=Iflag+64
      If(mu.gt.2.5 .and. mu.lt.Max(12.,15.*CL-1.5)) Then
C For small mu look for multiple solutions.  CL=.8 has its first jump at
C mu=3 and its 6th around 10.6; CL=.9 has its first jump at 3.9 and its
C 6th at 12.01; the table of CMax ends at 12 and covers CL=.8+ to near 1.
         BotNow=Min(2.5,mu-1.5)
         TopNow=BotNow
         FUpbot=FUp(BotNow,1)
         Do I=1,300
C         Do I=1,1000
            TopNow=TopNow+.01
C            TopNow=TopNow+.003
            If(TopNow.gt.12.) Go to 100
            FUpnow=FUp(TopNow,1)
            If(FUpbot*FUpnow.lt.0.) Then
               Call RZERO(BotNow,TopNow,mucand(I1),R,EPS,MAXF,FUp)
               mucand(I1)=mucand(I1)-muB
               BotNow=TopNow
               FUpbot=FUpnow
               I1=I1+1
               If(I1.gt.3) Go to 100
            EndIf
         EndDo
      EndIf
 100  UpperLim=mu-muB
      If(I1.gt.2) Then
         UpperLim=mucand(3)-muB
         Exclude_low(1)=mucand(1)-muB
         Exclude_low(2)=mucand(2)-muB
         Iflag=Or(Iflag,256)
      EndIf
      Iflag=Or(Iflag,16*IflagOpt)
C mmax is now the number of events in the optimum interval.  Find the
C endpoints of the optimum interval.
      fhigh=-1.
      Do I=0,N-mmax
         fdiff=FC(I+mmax+1)-FC(I)
         if(fdiff.gt.fhigh) Then
            fhigh=fdiff
            EndPoints(1)=I
            EndPoints(2)=I+mmax+1
         EndIf
      EndDo
      If(debug) Then
          Write(6,200) N,mmax,mu,fdebug,ydebug,fhigh,EndPoints
 200     Format('N, m, mu, f(m), y_vs_CLf(CMax,f(m)), fhigh, Endpoints',
     1 2I5,F12.2,3F9.5,2I6)
C        Write(6,210) mmax,(FC(I),I=1,N)
C 210     Format('m, FC: ',I4,(10F8.5))
C         Write(6,215) (f(I),Cdbg(I),I=0,N)
C 215     Format('f,C: ',10F8.5)
      EndIf
      Return
      End
      Real Function Fup(x,I)
      Implicit None
      Integer N1,NMax,I,If1,Istat,Icode,m,NCalls,I1,mmax
      Parameter (NMax=150000)
C      Real f(0:NMax),x,C,y,CMax,CL1,MeanMax,Cinf,CMaxinf,ConfLev,
      Real f(0:NMax),x,C,y,CMax,CL1,MeanMax,Cinf,CMaxinf,ConfLevNew,
     1 CMxinf,FC1,FB1,muB1,fmin1(7),Cdbg
      Common/Fupcom/f,N1,CL1,If1,MeanMax,NCalls,Istat,FC1(0:NMax),
     1 FB1(0:NMax),muB1,fmin1,Cdbg(0:10),mmax
      Logical Debug2/.false./
      Real FUpsave
      save FUpsave
      If(I.eq.1) NCalls=0
      CMax=0.
      CMxinf=CMaxinf(CL1,If1,x)
      Do m=0,N1
         If(muB1.ne.0.) Then
            f(m)=0.
            Do I1=0,N1-m
               f(m)=Max(f(m),(1.-(muB1/x))*(FC1(I1+m+1)-FC1(I1))+
     1             (muB1/x)*(FB1(I1+m+1)-FB1(I1)))
            EndDo
         EndIf
         If(f(m).gt.fmin1(If1)) Then
          If(x.le.MeanMax) Then
C            C=ConfLev(m,x*f(m),x,istat)
            C=ConfLevNew(m,x*f(m),x,istat)
            If(Istat.eq.2) Then
C ConfLev(New) failed.
               Istat=5
               Fup=0.
               Return
            EndIf
          Else
            y=(float(m)-x*f(m))/sqrt(x*f(m))
            C=Cinf(y,f(m),Istat)
          EndIf
          If(m.lt.10) Cdbg(m)=C
          If(C.gt.CMax) Then
             CMax=C
             mmax=m
          EndIf
          If(f(m).ge.1.) Go to 10
         EndIf
      EndDo
 10   FUp=CMax - CMxinf
      If(Debug2) Then
         If(I.eq.3) Then
            Write(6,*) "ncalls, mmax, f(mmax), CMax:",
     1       ncalls,mmax,f(mmax),CMax
         EndIf
         If(NCalls.eq.0) FUpsave=FUp
         If(NCalls.eq.1 .and. FUp*FUPsave.gt.0.) Then
           Write(6,*) "RZERO will fail: ",x,muB1,CMax,mmax
         EndIf
      EndIf
      NCalls=NCalls+1
      Return
      End

