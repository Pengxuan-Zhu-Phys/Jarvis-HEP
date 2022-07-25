      Real Function y_vs_CLf(CLin,fin,Istat)
C Computes y such that CLin=Cinf(y,fin).  .8 < CLin < .9999, .005 < f < 1
C Input: CLin=Confidence level, fin=fraction of the range, sometimes called f.
C Output: y=y_vs_CLf=y_{infinity}(CLin,fin), Istat = indicator of quality.
C         If Istat .ge. 3, the program failed, and returns y_vs_CLf=-10.
C Istat  Meaning
C   0    y_vs_CLf.in table used to interpolate
C   1    y_vs_CLf.in table used to extrapolate from f0=0.01 to low fin
C   2    Normal distribution used to extrapolate fron f0=1 to very low fin
C   3    failure: CLin > 0.9999
C   4    failure: fin > 1
C   5    failure: fin > 0.01 but CLin < 0.8
C   6    failure: unable to extrapolate for very low fin.
      Implicit None
      Real CL,f,ylow,yhigh,table(425,100),xf,xCL,FNf,fbin,dfbin,
     1 CLs(425),dCL,flog,ytemp,CLin,fin,f0/.01/,f0a/.1/,Ctop/.9999/,
     2 Cbot/0.8/
C As the program was written in 2004, the lowest value of f in the
C y_vs_CLf.in table was 0.01. The table goes from CL=.8 to CL=.999.
C Sometimes f0 is too small, in which case use f0a.
C      Real*8 PPND16,C,df
      Real*8 DGAUSN,C,df
      Integer Istat,Nmin,Nf,Ntable,Ntrials,I,J,If,ICL
      Common/y_vs_CLfcom/ylow,yhigh,Nmin,Nf,NTable,Ntrials,FNf,table,
     1 CLs
     
      character*1000 PATH, FILENAME,catpath
      COMMON/maxGapPath/PATH,FILENAME    
      external catpath
                        
      Logical first/.true./
      If(first) Then
         first=.false.
         FILENAME=catpath(PATH,'y_vs_CLf.in')
         open(20,file=FILENAME,status='OLD',form='UNFORMATTED')
         read(20) ylow,yhigh,Nmin,Nf,NTable,Ntrials
         FNf=float(Nf)
         Do If=1,Nf
            Read(20) (table(J,If),J=1,425)
         EndDo
         Close(20)
         Do ICL=1,175
            CLs(ICL)=.799+.001*Float(ICL)
         EndDo
         Do ICL=176,425
            CLs(ICL)=.9574+.0001*Float(ICL)
         EndDo
      EndIf
      y_vs_CLf=-10. ! Default for failure
      Istat=0 ! Default for success (but 1 and 2 aren't too bad)
      If(CLin .gt. Ctop) Then
         Istat=3
         Return
      ElseIf(fin .gt. 1.) Then
         Istat=4
         Return
      ElseIf(fin.ge.0.01 .and. CLin .lt. Cbot) Then
         Istat=5
         Return
      EndIf
      If(fin.lt.0.01) Then
C If f<0.01 and CL<0.9999, extrapolate from f0 or f0a.
         Istat=1
         CL=CLin**(((1./f0)-.94)/((1./fin)-.94))
         f=f0
         If(CL.lt.Cbot) Then
            CL=CLin**(((1./f0a)-.94)/((1./fin)-.94))
            f=f0a
         EndIf
         If(CL.gt.Ctop .or. CL.lt.Cbot) Then
C For very small f, CL can get too close to 1 for the table to handle.
C Even greater than CL=.999 involves extrapolation of the table.
C But for CL>0.9999, extrapolate with the inverse of the normal distribution.
C Give a rough estimate of y_vs_CLf based on extrapolation from f=1.
C This part depends on a routine from http://lib.stat.cmu.edu/apstat/241,
C published in Applied Statistics, vol. 37, pp. 477-484, 1988.  If the
C optimum interval method is used with fewer than, say, 1000 events, I
C expect DGAUSN never to be called, so the block calling it can be removed.
           Istat=2
           df=fin
           C=CLin
C The constants in the next executable line were tuned so that the
C extrapolated results would be pretty good compared with the table
C (for df=.01) and with the extrapolation from f0=.01 for df=0.0005.
           C=C**(0.051D0/((1.D0/df)-0.946D0))
           y_vs_CLf=-DGAUSN(C,I)
C DGAUSN shouldn't fail if f>1.E-12 and CLin<0.999.  This formula
C is essentially y_vs_CLf(C,1).
           If(I.eq.1) Then
              y_vs_CLf=-10.
              Istat=6
           EndIf
           Return
         EndIf
      Else
         CL=CLin
         f=fin
      EndIf
      If(f.gt. 1.) Then
C Shouldn't be able to get here
         Write(6,*) f," is an out of range value of f in y_vs_CLf"
         stop
      EndIf
      fbin=f*FNf
      If=fbin
      If(If.lt.1) If=1
      If(If.gt.Nf-1) If=Nf-1
      xf=fbin-float(If)
      If(CL.lt.Cbot .or. CL.gt.Ctop) Then
C Shouldn't be able to get here
         Write(6,*) CL," is an out of range value of CL in y_vs_CLf"
         stop
      ElseIf(CL.lt..975) Then
         ICL=1000.*(CL-.799) ! ICL is from 1 to 175
         dCL=.001
      Else
         ICL=10000.*(CL-.9574) ! ICL is from 176 to 424
         dCL=.0001
      EndIf
      xCL=(CL-CLs(ICL))/dCL
      ytemp=(1.-xCL)*((1.-xf)*table(ICL,If)+xf*table(ICL,If+1))
     1 + xCL*((1.-xf)*table(ICL+1,If)+xf*table(ICL+1,If+1))
      flog=log(f)
      y_vs_CLf=(ytemp+1.7*flog)/(1.-0.3*flog)
      Return
      End
