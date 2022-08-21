      Real Function ConfLevNew(m,x,mu,icode)
C Evaluate the probability that the maximum expected contents of an interval
C with at most m points is less than x given that the total range is mu.
C The return code, icode, is 0 if all is ok, 1 if the true value is guaranteed
C to not be above the returned value, 2 if the program can't supply any useful
C information.  ConfLev is supposed to be used only for finding confidence
C levels above 50%; if the program finds it below this, it may not bother
C finding a more accurate result.
C
C ConfLevNew changed from ConfLev.f in Jan 2010 to use CLtableNew.in
C when that table made with more statistics and finer binning is useable.
C Otherwise it uses CLtable.in via the old ConfLev.f
      Implicit None
      Integer m,icode
      Real x,mu
      Integer N,Nmeans,NCLs
C      Parameter(N=50)
C      Parameter(Nmeans=41)
C      Parameter(NCLs=22)
      Parameter(N=5) ! Jan 2010
      Parameter(Nmeans=51) ! Jan 2010 (From CtableNew where L is 0-50)
      Parameter(NCLs=44) ! Jan 2010
C      Common/CLtab/CL,Mean,xcl,Nxcl,meanlog ! For aid in debugging
      Common/CLtabNew/CL,Mean,xcl,Nxcl,meanlog ! For aid in debugging, save
      Real CL(NCLs),xcl(0:N,NCLs,Nmeans),mean(Nmeans),ConfLev
      Real meanlog(Nmeans),mulog,xcl2(NCLs+1),clxmu,GAMDIS,
     x CL2(NCLs+1),UseNewConfL
      Logical first/.true./,UseNew/.true./,use
      Integer I,J,JJ,K,L,LL,Nxcl(0:N,Nmeans),Nxcl2
      character*1000 PATH, FILENAME,catpath
      COMMON/maxGapPath/PATH,FILENAME    
      external catpath
      
      If(first) Then
C Input the table.  xcl(m,J,I) is the maximum interval containing m points
C corresponding to mu=mean(I) and confidence level CL(J).  If Nxcl(m,I)=0,
C then the m is so large compared with mu that one hardly ever gets m
C events; so if x<mu, the confidence level is very small.  I assume that if
C xcl(m,J,I) is meaningful then so is xcl(m,J,I+1) and xcl(m,J-1,I)
C (unless J=1).
       First=.false.
       If(UseNew) Then
        FILENAME=catpath(PATH,'CLtableNew.in')
        Open(21,file=FILENAME,status='OLD',form='UNFORMATTED') ! Jan 2010
        Read(21) CL ! There are places where I assume CL(2) .le. .5
        Do I=1,Nmeans
          Read(21) mean(I),(Nxcl(K,I),K=0,N)
          Do K=0,N
             If(Nxcl(K,I).gt.0) Read(21) (xcl(K,J,I),J=1,Nxcl(K,I))
          EndDo
          meanlog(I)=log(mean(I))
        EndDo
        Close(21)
       EndIf
      EndIf
      If(mu.gt.12. .or. m.gt.5 .or. .not. UseNew) Then
         ConfLevNew = ConfLev(m,x,mu,icode)
         Return
      EndIf
C Make some simple checks on reasonableness of the input parameters.  If
C they are unreasonable, either return a reasonable output anyway, or give up.
      If(mu.lt.0. .or. m.lt.0) Then
         Go to 1000
      ElseIf(x.gt.1.00001*mu) Then
         ConfLevNew=1.
      ElseIf(x.le.0.) Then
         ConfLevNew=0.
      ElseIf(mu.lt.1. .or. mu.ge.mean(Nmeans) .or. m.gt.N) Then
C The table hasn't the needed information, but at least we know that the
C answer is smaller than the confidence level for x=mu-.
         ConfLevNew=GAMDIS(mu,float(m+1))
         Go to 2000
      Else ! The table might include this m,mu
C Find which mean(I) is just below mu.  Since mu is less than mean(Nmeans),
C the I must be less than Nmeans; so there is information for I+1.
         mulog=log(mu) ! Since mu.ge.1, mulog is at least 0.
C         I=10.*mulog + 1. ! mu = exp(.1*L) L=0,1,...  I=L+1
         I=20.*mulog + 1. ! mu = exp(.05*L) Jan 2010
         Nxcl2=min(Nxcl(m,I),Nxcl(m,I+1))
         If(Nxcl2.lt.2) Then
C The m must be too large for the mu; CL is very low.
            ConfLevNew=CL(2)
            Go to 2000
         EndIf
C Nxcl(m,I) is at least 2 after this point.
C Find clxmu, the confidence level for x=mu-.
         clxmu=GAMDIS(mu,float(m+1))
C Return clxmu if x is equal to mu to within one part in 100,000.
         If(x .gt. .99999*mu) Then
            ConfLevNew=clxmu
            Go to 500
         EndIf
C Interpolate table entries to apply to this value of mu
         Do J=1,Nxcl2
            xcl2(J)=xcl(m,J,I) + (mulog-meanlog(I))*
     x      (xcl(m,J,I+1)-xcl(m,J,I))/(meanlog(I+1)-meanlog(I))
            CL2(J)=CL(J)
         EndDo
         xcl2(Nxcl2+1)=mu
         CL2(Nxcl2+1)=clxmu
C Interpolate xcl2 to find ConfLev at x.  K = 2 or 3 xcl2's are used,
C beginning with xcl2(J).  xcl2(1) and xcl2(2) must exist to get here.
C If x is before xcl2(2), then xcl2(1), xcl2(2), and xcl2(3) are
C used.  Analogously, if x is after xcl2(Nxcl(m,I)), use that xcl2, the one
C above it, and the one below it.  Otherwise use the two xcl2's 
C between which x lies, and the nearest other one to x.
         K=3
         If(x.lt.xcl2(1)) Then
C Only linearly extrapolate before xcl2(1); otherwise the parabola could
C result in smaller x getting larger confidence level
            J=1
            K=2
         ElseIf(x.le.xcl2(2)) Then
            J=1
C         ElseIf(Nxcl(m,I).eq.NCLs .and. x.gt.xcl2(NCLs)) Then
C Only linearly extrapolate beyond xcl2(NCLs); otherwise the parabola could
C result in larger x getting smaller confidence level.
C            J=NCLs
C            K=2
         ElseIf(x.ge.xcl2(Nxcl2)) Then
            J=Nxcl2-1
         Else ! x is between xcl2(2) and xcl2(Nxcl2); Nxcl2>3.
            Do L=3,Nxcl2 ! Find the first xcl2 that's above x
               If(x.le.xcl2(L)) Then
                  If(xcl2(L+1)-x .lt. x-xcl2(L-2)) Then
                     J=L-1
                  Else
                     J=L-2
                  EndIf
                  Go to 100
               EndIf
            EndDo
 100        Continue
         EndIf
C Maybe at this point I should improve the xcl2 that are used: quadratically
C interpolate the xcl over three neighboring I.  I tried it; no improvement.
C Find JJ, the first of the triplet used for quadratic interpolation
C in improving the xcl2's.  Normally use I, I+1, and I+2 for the interpolation,
C for if the xcl exist for I and I+1, then they do for I+2 unless I+2>Nmeans.
C         If(I+2.gt.Nmeans) Then
C            JJ=I-1
C         Else
C            JJ=I
C         EndIf
C         Do L=1,K
C            LL=J+L-1
C            If(LL.LE.Nxcl(m,I)) Then
C             xcl2(LL)=xcl(m,LL,JJ) + ((mulog-meanlog(JJ))/(meanlog(JJ+1)
C     x        - meanlog(JJ+2))) * ( (xcl(m,LL,JJ+1)-xcl(m,LL,JJ))*
C     x        (mulog-meanlog(JJ+2))/(meanlog(JJ+1)-meanlog(JJ)) +
C     x        (xcl(m,LL,JJ+2)-xcl(m,LL,JJ))*(meanlog(JJ+1)-mulog)/
C     x        (meanlog(JJ+2)-meanlog(JJ)) )
C            EndIf
C         EndDo
         If(K.eq.2) Then ! Linearly interpolate
            ConfLevNew=CL2(J) + (CL2(J+1)-CL2(J))*(x - xcl2(J))/
     x       (xcl2(J+1) - xcl2(J))
         Else ! Quadratically interpolate
            ConfLevNew=CL2(J) + ((x-xcl2(J))/(xcl2(J+1)-xcl2(J+2))) * (
     x       (CL2(J+1)-CL2(J))*(x-xcl2(J+2))/(xcl2(J+1)-xcl2(J)) +
     x       (CL2(J+2)-CL2(J))*(xcl2(J+1)-x)/(xcl2(J+2)-xcl2(J)) )
         EndIf
         If(ConfLevNew.lt.0.) ConfLevNew=0.
         If(ConfLevNew.gt.clxmu) ConfLevNew=clxmu
      EndIf ! End of the case in which the table can include the m and mu.
 500  icode=0
      return
 1000 ConfLevNew=-1.
      icode=2
      return
 2000 icode=1
      return
      Entry UseNewConfL(use)
C Turn on or off use of new, more accurate table.
      UseNew=use
      First=.true.
      UseNewConfL=0.
      Return
      end
