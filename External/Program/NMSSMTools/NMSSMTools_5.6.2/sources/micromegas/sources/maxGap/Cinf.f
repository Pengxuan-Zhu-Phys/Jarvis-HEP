      Real Function Cinf(y,fin,Istat)
C Calculates C_\infty from a table.
C Input: y = deviation, fin = fraction of the range
C Output: Cinf=C_\infty(y,fin), Istat = return status code
C Istat  Meaning
C   0    Result interpolated from table ymintable.in
C   1    Result extrapolated from f0=.01 (lowest f in the table)
C   2    y especially low.  Returns Cinf=1.
C   3    y especially high.  Returns Cinf=0.
C   4    Result estimated using derfc (not implemented)
C   5    failure: fin > 1
C   6    failure: fin < 1.E-10
      Implicit none
      Integer Nf,Nmin,Ntrials,Ntable,If,Ibin,Istat
      Real y,ytemp,f,ytable,ylow,yhigh,flog,Table(2000,100),FNf,
     1 fbin,ybin,dfbin,dIbin,dy,fin,f0/.01/
      Logical first/.true./
      Common/Cinfcom/Nf,Ntable,ylow,yhigh,dy,FNf,Table

      character*1000 PATH, FILENAME,catpath
      COMMON/maxGapPath/PATH,FILENAME    
      external catpath
      
      If(first) Then
       first=.false.
       FILENAME=catpath(PATH,'ymintable.in')
       Open(50,file=FILENAME,status='OLD',form='UNFORMATTED')
       Read(50) ylow,yhigh,Nmin,Nf,NTable,Ntrials
       Do If=1,Nf
          Read(50) (Table(Ibin,If),Ibin=1,Ntable)
       EndDo
       dy=(yhigh-ylow)/float(Ntable)
       FNf=float(Nf)
       Close(50)
      EndIf
      Istat=0 ! Default for success
      Cinf=1. ! Default for failure
      f=fin
      If(f.gt.1.) Then
         Istat=4
         Return
      EndIf
      If(f.lt.1.E-10) Then
         Istat=5
         Return
      EndIf
      If(f.lt.0.01) Then
         Istat=1
         f=f0
      EndIf
      flog=log(f)
      ytemp=y*(1. - 0.3*flog) - 1.7*flog
      If(ytemp.lt.ylow) Then
         Istat=2
         Cinf=1.
         Return
      Elseif(ytemp.gt.yhigh) Then
         Istat=3
         Cinf=0.
         Return
      EndIf
      fbin=f*FNf
      If=fbin
      If(If.lt.1) If=1
C      If(If.gt.Nf-2) If=Nf-2
      If(If.gt.Nf-1) If=Nf-1
      dfbin=fbin-Float(If)
      ybin=(ytemp-ylow)/dy
      Ibin=ybin
      If(Ibin.lt.2) Ibin=2
      If(Ibin.gt.Ntable-1) Ibin=Ntable-1
      dIbin=ybin-Float(Ibin)
C To check the following interpolation formula, verify that it gives
C the correct 4 table entries when dfbin is close to 0 and 1 and when
C dIbin is close to 0 and 1.  Table(Ibin,If)=Cinf evaluated at
C ytemp=ylow+dy*Ibin, f=If/FNf.
      Cinf=(1.-dfbin)*(dIbin*Table(Ibin+1,If)+(1.-dIbin)*Table(Ibin,If))
     1    +dfbin*(dIbin*Table(Ibin+1,If+1)+(1.-dIbin)*Table(Ibin,If+1))
      If(Istat.eq.1) Cinf=Cinf**((1./fin - 0.94)/(1./f0 - 0.94))
      Return
      End

