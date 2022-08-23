      function hard_dijet() result(hard)
      use nrtype
      use comvar
      use pdf
      use ff
      use qed
      use qcd
      use sudakov
      implicit none
      real(sp) :: hard
      real(sp), dimension(-6:6) :: PFa,PFb,FFc,FFd
      real(sp) :: qq,qg,gg
      real(sp) :: SIG,SIGt,SIGu
      real(sp) :: DIS,DISt,DISu
      real(sp) :: SUD,SUDa,SUDb      
      real(sp) :: channel1,channel2,channel3,channel4
      real(sp) :: channel5,channel6,channel7,channel8
      integer :: i,j

      ! CTEQ18 PDF
      ! +1: proton
      ! -1: anti-proton
      call getpdf(+1,xa,mufac,PFa)
      call getpdf(+1,xb,mufac,PFb)
      ! akk fragmentation
      ! 6: pi^0
      ! 9: charged hadron
      if(nproc/10 .eq. 2) then
        call getff(6,zc,mufrg,FFc)
      else
        FFc = 1.0
      endif
      if(mod(nproc,10) .eq. 2) then
        call getff(6,zd,mufrg,FFd)
      else
        FFd = 1.0
      endif

      qq = 0.0; qg = 0.0; gg = 0.0

      if(qqch) then ! quark+quark initiated channel

      ! q q' -> q q'
      SIGt = 4.0/9.0 * (mans**2+manu**2)/mant**2
      SIGu = 4.0/9.0 * (mans**2+mant**2)/manu**2
      DISt = 0.0; DISu = 0.0
      do i = 1,Nf
      do j = 1,Nf
        if(i.ne.j) then
        DISt = DISt 
     &       + PFa(+i) * PFb(+j) * FFc(+i) * FFd(+j)
     &       + PFa(-i) * PFb(+j) * FFc(-i) * FFd(+j)
     &       + PFa(+i) * PFb(-j) * FFc(+i) * FFd(-j)
     &       + PFa(-i) * PFb(-j) * FFc(-i) * FFd(-j)
        DISu = DISu 
     &       + PFa(+i) * PFb(+j) * FFc(+j) * FFd(+i)
     &       + PFa(-i) * PFb(+j) * FFc(+j) * FFd(-i)
     &       + PFa(+i) * PFb(-j) * FFc(-j) * FFd(+i)
     &       + PFa(-i) * PFb(-j) * FFc(-j) * FFd(-i)
        endif
      enddo
      enddo
      SUD = sudfac(1,1,1,1)
      channel1 = SIGt * DISt * SUD
     &         + SIGu * DISu * SUD

      ! q qb -> q' qb'
      SIG = 4.0/9.0 * (mant**2+manu**2)/mans**2
      DISt = 0.0; DISu = 0.0
      do i = 1,Nf
      do j = 1,Nf
        if(i.ne.j) then
        DISt = DISt
     &       + PFa(+i) * PFb(-i) * FFc(+j) * FFd(-j)
     &       + PFa(-i) * PFb(+i) * FFc(-j) * FFd(+j)
        DISu = DISu
     &       + PFa(+i) * PFb(-i) * FFc(-j) * FFd(+j)
     &       + PFa(-i) * PFb(+i) * FFc(+j) * FFd(-j)
        endif
      enddo
      enddo
      channel2 = SIG * DISt * SUD
     &         + SIG * DISu * SUD

      ! q qb -> q qb
      SIGt = 4.0/9.0 * ( (mans**2+manu**2)/mant**2
     &                  +(mant**2+manu**2)/mans**2
     &                  -2.0/3.0 * manu**2/mans/mant )
      SIGu = 4.0/9.0 * ( (mans**2+mant**2)/manu**2
     &                  +(manu**2+mant**2)/mans**2
     &                  -2.0/3.0 * mant**2/mans/manu )
      DISt = 0.0; DISu = 0.0
      do i = 1,Nf
        DISt = DISt
     &       + PFa(+i) * PFb(-i) * FFc(+i) * FFd(-i)
     &       + PFa(-i) * PFb(+i) * FFc(-i) * FFd(+i)
        DISu = DISu
     &       + PFa(+i) * PFb(-i) * FFc(-i) * FFd(+i)
     &       + PFa(-i) * PFb(+i) * FFc(+i) * FFd(-i)
      enddo
      channel3 = SIGt * DISt * SUD
     &         + SIGu * DISu * SUD

      ! q q -> q q
      SIG = 4.0/9.0*( (mans**2+manu**2)/mant**2
     &               +(mans**2+mant**2)/manu**2
     &               -2.0/3.0*mans**2/mant/manu )
      DIS = 0.0
      do i = 1,Nf
        DIS = DIS
     &      + PFa(+i) * PFb(+i) * FFc(+i) * FFd(+i)
     &      + PFa(-i) * PFb(-i) * FFc(-i) * FFd(-i)
      enddo
      channel4 = SIG * DIS * SUD

      ! q qb -> g g
      SIG = 8.0/3.0 * (4.0/9.0*(mant**2+manu**2)/mant/manu
     &                       -(mant**2+manu**2)/mans**2 )
      DIS = 0.0
      do i = 1,Nf
        DIS = DIS
     &      + PFa(+i) * PFb(-i) * FFc(0) * FFd(0)
     &      + PFa(-i) * PFb(+i) * FFc(0) * FFd(0)
      enddo
      SUD = sudfac(1,1,0,0)
      channel5 = SIG * DIS * SUD

      qq = channel1
     &   + channel2
     &   + channel3
     &   + channel4
     &   + channel5

      endif

      if(qgch) then ! quark+gluon initiated channel

      ! q g -> q g
      SIGt = (mans**2+manu**2)/mant**2 
     &       -4.0/9.0*(mans**2+manu**2)/mans/manu
      SIGu = (mans**2+mant**2)/manu**2 
     &       -4.0/9.0*(mans**2+mant**2)/mans/mant
      SUDa = sudfac(1,0,1,0)
      SUDb = sudfac(1,0,0,1)
      DISt = 0.0; DISu = 0.0
      do i = 1,Nf
        DISt = DISt 
     &       + PFa(+i) * PFb( 0) * FFc(+i) * FFd( 0) * SUDa
     &       + PFa(-i) * PFb( 0) * FFc(-i) * FFd( 0) * SUDa
     &       + PFa( 0) * PFb(+i) * FFc( 0) * FFd(+i) * SUDb
     &       + PFa( 0) * PFb(-i) * FFc( 0) * FFd(-i) * SUDb
        DISu = DISu
     &       + PFa(+i) * PFb( 0) * FFc( 0) * FFd(+i) * SUDb
     &       + PFa(-i) * PFb( 0) * FFc( 0) * FFd(-i) * SUDb
     &       + PFa( 0) * PFb(+i) * FFc(+i) * FFd( 0) * SUDa
     &       + PFa( 0) * PFb(-i) * FFc(-i) * FFd( 0) * SUDa
      enddo
      channel6 = SIGt * DISt
     &         + SIGu * DISu

      qg = channel6
      
      endif

      if(ggch) then ! gluon+gluon initiated channel

      ! g g -> q qb
      SIG = 3.0/8.0 * (4.0/9.0*(mant**2+manu**2)/mant/manu
     &                   -(mant**2+manu**2)/mans**2 )
      DISt = 0.0; DISu = 0.0
      do i = 1,Nf
        DISt = DISt
     &       + PFa(0) * PFb(0) * FFc(+i) * FFd(-i)
        DISu = DISu
     &       + PFa(0) * PFb(0) * FFc(-i) * FFd(+i)
      enddo
      SUD = sudfac(0,0,1,1)
      channel7 = SIG * DISt * SUD
     &         + SIG * DISu * SUD

      ! g g -> g g
      SIG = 9.0/2.0 * ( 3.0 - manu*mant/mans**2
     &                      - manu*mans/mant**2
     &                      - mans*mant/manu**2 )
      DIS = PFa(0) * PFb(0) * FFc(0) * FFd(0)
      SUD = sudfac(0,0,0,0)
      channel8 = SIG * DIS * SUD

      gg = channel7
     &   + channel8

      endif

      ! sum all channels and overall normalization
      hard = qq + qg + gg
      hard = hard * alphas(muren)**2 / mans**2 * PI
      hard = hard * xa*xb * 2.0*ptj

      ! debugging options (nan) (before filling histogram)
      ! compile without -ffast-math option
      if(isnan(hard)) then
        write(*,*) 'nan error:',hard
        write(*,*) 'momfrac:',xa,xb,zc,zd
        write(*,*) 'mandelstam:',mans,mant,manu
        write(*,*) 'scales:',muren,mufac,mures,mufrg
        write(*,*) 'channels:'
        write(*,*) channel1,channel2,channel3,channel4
        write(*,*) channel5,channel6,channel7,channel8
        write(*,*) qq,qg,gg
        write(*,*) alphas(muren)
      endif

      return

      end function hard_dijet
