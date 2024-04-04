      SUBROUTINE NS_GRAVITINO

************************************************************************
*
*     This subroutine computes the gravitino decays
*
************************************************************************

      IMPLICIT NONE

      INTEGER i,j,GRFLAG

      DOUBLE PRECISION gravul,gravur,gravdl,gravdr,gravt1,gravt2,gravb1,
     .         gravb2,gravel,graver,gravm1,gravm2,gravl1,gravl2,gravne,
     .         gravnm,gravnl,gravch(2),gravcw(2),gravng(5),gravnz(5),
     .         gravnh(5,3),gravna(5,2),gravgg
      DOUBLE PRECISION brgravul,brgravur,brgravdl,brgravdr,brgravt1,
     .         brgravt2,brgravb1,brgravb2,brgravel,brgraver,brgravm1,
     .         brgravm2,brgravl1,brgravl2,brgravne,brgravnm,brgravnl,
     .         brgravch(2),brgravcw(2),brgravng(5),brgravnz(5),
     .         brgravnh(5,3),brgravna(5,2),brgravgg,gravtot
      DOUBLE PRECISION KNG(5),KNZ(5),KNH(5,3),KNA(5,2),KCW(2),KCH(2)
      DOUBLE PRECISION M32,CGR,MPL
      DOUBLE PRECISION mh(3),S(3,3),ma(2),P2(2,2),mhc
      DOUBLE PRECISION ms,mc,mb,amb,amt,amtau,ammuon,amz,amw
      DOUBLE PRECISION mgluino,amneut(5),xmneut(5),amchar(2),xmchar(2)
      DOUBLE PRECISION asupr,asupl,asdownr,asdownl,aser,asel,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

      COMMON/GRAVITINO_BR/brgravul,brgravur,brgravdl,brgravdr,brgravt1,
     .         brgravt2,brgravb1,brgravb2,brgravel,brgraver,brgravm1,
     .         brgravm2,brgravl1,brgravl2,brgravne,brgravnm,brgravnl,
     .         brgravch,brgravcw,brgravng,brgravnz,brgravnh,brgravna,
     .         brgravgg
      COMMON/GRAVITINO_WIDTH/gravtot
      COMMON/GRAVCOUP/KNG,KNZ,KNH,KNA,KCW,KCH
      COMMON/M32/M32,CGR,MPL,GRFLAG
      COMMON/HIGGSPEC/mh,S,ma,P2,mhc
      COMMON/SMSPEC/ms,mc,mb,amb,amt,amtau,ammuon,amz,amw
      COMMON/NS_massgino/mgluino,amneut,xmneut,amchar,xmchar
      COMMON/SFSPEC/asupr,asupl,asdownr,asdownl,aser,asel,asne1,
     .         ast1,ast2,asb1,asb2,astau1,astau2,asntau1,
     .         CST,CSB,CSL,asmu1,asmu2,asnmu1,csmu

c ---Initialization
      gravul=0d0
      gravur=0d0
      gravdl=0d0
      gravdr=0d0
      gravt1=0d0
      gravt2=0d0
      gravb1=0d0
      gravb2=0d0
      gravel=0d0
      graver=0d0
      gravm1=0d0
      gravm2=0d0
      gravl1=0d0
      gravl2=0d0
      gravne=0d0
      gravnm=0d0
      gravnl=0d0
      do i=1,2
        gravch(i)=0d0
      enddo
      do i=1,2
        gravcw(i)=0d0
      enddo
      do i=1,5
        gravng(i)=0d0
      enddo
      do i=1,5
        gravnz(i)=0d0
      enddo
      do i=1,5
        do j=1,3
          gravnh(i,j)=0d0
        enddo
      enddo
      do i=1,5
        do j=1,2
          gravna(i,j)=0d0
        enddo
      enddo
      gravgg=0d0

c --- gravitino --> squark quark
      if(M32.GT.asupl)then
        gravul=0d0
      endif
      if(M32.GT.asupr)then
        gravur=0d0
      endif
      if(M32.GT.asdownl)then
        gravdl=0d0
      endif
      if(M32.GT.asdownr)then
        gravdr=0d0
      endif
      if(M32.GT.(ast1+amt))then
        gravt1=0d0
      endif
      if(M32.GT.(ast2+amt))then
        gravt2=0d0
      endif
      if(M32.GT.(asb1+amb))then
        gravb1=0d0
      endif
      if(M32.GT.(asb2+amb))then
        gravb2=0d0
      endif

c --- gravitino --> slepton lepton
      if(M32.GT.asel)then
        gravel=0d0
      endif
      if(M32.GT.aser)then
        graver=0d0
      endif
      if(M32.GT.asmu1)then
        gravm1=0d0
      endif
      if(M32.GT.asmu2)then
        gravm2=0d0
      endif
      if(M32.GT.astau1+amtau)then
        gravl1=0d0
      endif
      if(M32.GT.astau2+amtau)then
        gravl2=0d0
      endif
      if(M32.GT.asne1)then
        gravne=0d0
      endif
      if(M32.GT.asnmu1)then
        gravnm=0d0
      endif
      if(M32.GT.asntau1)then
        gravnl=0d0
      endif

c --- gravitino --> chargino H+
      do i=1,2
        if(M32.GT.(amchar(i)+mhc))then
          gravch(i)=0d0
        endif
      enddo

c --- gravitino --> chargino W+
      do i=1,2
        if(M32.GT.(amchar(i)+amw))then
          gravcw(i)=0d0
        endif
      enddo

c --- gravitino --> neutralino gamma
      do i=1,5
        if(M32.GT.amneut(i))then
          gravng(i)=0d0
        endif
      enddo

c --- gravitino --> neutralino Z
      do i=1,5
        if(M32.GT.(amneut(i)+amz))then
          gravnz(i)=0d0
        endif
      enddo

c --- gravitino --> neutralino H
      do i=1,5
        do j=1,3
          if(M32.GT.(amneut(i)+mh(j)))then
            gravnh(i,j)=0d0
          endif
        enddo
      enddo

c --- gravitino --> neutralino A
      do i=1,5
        do j=1,2
          if(M32.GT.(amneut(i)+ma(j)))then
            gravna(i,j)=0d0
          endif
        enddo
      enddo

c --- gravitino --> gluino gluon
      if(M32.GT.mgluino)then
        gravgg=0d0
      endif

c --- total width
      gravtot=4d0*(gravul+gravur+gravdl+gravdr)+
     .        2d0*(gravt1+gravt2+gravb1+gravb2+gravel+graver+
     .        gravm1+gravm2+gravl1+gravl2+gravne+gravnm+gravnl)

      do i=1,2
        gravtot=gravtot+2d0*gravch(i)
      enddo

      do i=1,2
        gravtot=gravtot+2d0*gravcw(i)
      enddo

      do i=1,5
        gravtot=gravtot+gravng(i)
      enddo

      do i=1,5
        gravtot=gravtot+gravnz(i)
      enddo

      do i=1,5
        do j=1,3
          gravtot=gravtot+gravnh(i,j)
        enddo
      enddo

      do i=1,5
        do j=1,2
          gravtot=gravtot+gravna(i,j)
        enddo
      enddo

      gravtot=gravtot+gravgg

c --- branching ratios

      if(gravtot.ne.0d0)then

      brgravul=gravul/gravtot
      brgravur=gravur/gravtot
      brgravdl=gravdl/gravtot
      brgravdr=gravdr/gravtot
      brgravt1=gravt1/gravtot
      brgravt2=gravt2/gravtot
      brgravb1=gravb1/gravtot
      brgravb2=gravb2/gravtot

      brgravel=gravel/gravtot
      brgraver=graver/gravtot
      brgravm1=gravm1/gravtot
      brgravm2=gravm2/gravtot
      brgravl1=gravl1/gravtot
      brgravl2=gravl2/gravtot
      brgravne=gravne/gravtot
      brgravnm=gravnm/gravtot
      brgravnl=gravnl/gravtot

      do i=1,2
        brgravch(i)=gravch(i)/gravtot
      enddo

      do i=1,2
        brgravcw(i)=gravcw(i)/gravtot
      enddo

      do i=1,5
        brgravng(i)=gravng(i)/gravtot
      enddo

      do i=1,5
        brgravnz(i)=gravnz(i)/gravtot
      enddo

      do i=1,5
        do j=1,3
          brgravnh(i,j)=gravnh(i,j)/gravtot
        enddo
      enddo

      do i=1,5
        do j=1,2
          brgravna(i,j)=gravna(i,j)/gravtot
        enddo
      enddo

      brgravgg=gravgg/gravtot

      else

      brgravul= 0d0
      brgravur= 0d0
      brgravdl= 0d0
      brgravdr= 0d0
      brgravt1= 0d0
      brgravt2= 0d0
      brgravb1= 0d0
      brgravb2= 0d0

      brgravel= 0d0
      brgraver= 0d0
      brgravm1= 0d0
      brgravm2= 0d0
      brgravl1= 0d0
      brgravl2= 0d0
      brgravne= 0d0
      brgravnm= 0d0
      brgravnl= 0d0

      do i=1,2
        brgravch(i)= 0d0
      enddo

      do i=1,2
        brgravcw(i)= 0d0
      enddo

      do i=1,5
        brgravng(i)= 0d0
      enddo

      do i=1,5
        brgravnz(i)= 0d0
      enddo

      do i=1,5
        do j=1,3
          brgravnh(i,j)= 0d0
        enddo
      enddo

      do i=1,5
        do j=1,2
          brgravna(i,j)= 0d0
        enddo
      enddo

      brgravgg= 0d0

      endif

      end
