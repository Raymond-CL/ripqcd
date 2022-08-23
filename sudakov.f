      module sudakov

      use nrtype
      use comvar
      use qcd, only : Ca,Cf,Tr,Nf,qcdloop,b0,b1,Lqcd2

      implicit none
      real(dp), parameter :: C1 = 2d0*exp(-Euler)
      real(dp), parameter :: A1q = Cf/twopi
      real(dp), parameter :: A1g = Ca/twopi
      real(dp), parameter :: Dq =  Cf/twopi
      real(dp), parameter :: Dg =  Ca/twopi
      real(dp), parameter :: bmax = 1.5d0
      real(dp) :: beta,K,A2q,A2g,B1q,B1g
      public :: setsud,sudfac,bmax

      contains

      subroutine setsud()
      beta = (33d0-2d0*dble(Nf))/36d0
      K = (67d0/18d0 - pi2/6d0)*Ca - 10d0/9d0*dble(Nf)*Tr
      B1q = -3d0/2d0*Cf/twopi
      B1g = -2d0*beta*Ca/twopi
      A2q = K * Cf/twopi2
      A2g = K * Ca/twopi2
      end subroutine setsud

      function sudfac(ia,ib,oc,od)
      real(dp) :: sudfac
      integer(i4b), intent(in) :: ia,ib,oc,od
      real(sp) :: A1a,A2a,B1a,Gfa
      real(sp) :: A1b,A2b,B1b,Gfb
      real(sp) :: A1c,A2c,B1c,Gfc,Gqc,Djc
      real(sp) :: A1d,A2d,B1d,Gfd,Gqd,Djd
      real(sp) :: A1s,A2s,B1s,            Gfs,    Djs
      real(sp) :: Lqr,Lqf,Lql,Lrf,Lrl,Lfl, LLrl,LLfl
      real(sp) :: SP,SNP,Smed

      ! incoming parton: a
      if(ia.eq.0) then
        A1a = A1g;      A2a = A2g;        B1a = B1g
        Gfa = CA/CF
      else
        A1a = A1q;      A2a = A2q;        B1a = B1q
        Gfa = 1.0
      endif
      ! incoming parton: b
      if(ib.eq.0) then
        A1b = A1g;      A2b = A2g;        B1b = B1g
        Gfb = CA/CF
      else
        A1b = A1q;      A2b = A2q;        B1b = B1q
        Gfb = 1.0
      endif
      ! outgoing parton: c
      if(oc.eq.0) then
        A1c = A1g/2.0;  A2c = A2g/2.0;    B1c = B1g/2.0
        Gfc = CA/CF;    Gqc = CA/CF
        Djc = Dg;
      else
        A1c = A1q/2.0;  A2c = A2q/2.0;    B1c = B1q/2.0
        Gfc = 1.0;      Gqc = 1.0
        Djc = Dq
      endif
      ! outgoing parton: d
      if(od.eq.0) then
        A1d = A1g/2.0;  A2d = A2g/2.0;    B1d = B1g/2.0
        Gfd = CA/CF;    Gqd = CA/CF
        Djd = Dg
      else
        A1d = A1q/2.0;  A2d = A2q/2.0;    B1d = B1q/2.0
        Gfd = 1.0;      Gqd = 1.0
        Djd = Dq
      endif

      ! setting: A1s,A2s,B1s,Djs
      A1s = A1a + A1b;  A2s = A2a + A2b;  B1s = B1a + B1b
      Gfs = Gfa + Gfb
      Djs = 0.0
      ! for hadronic final state
      if(nproc.ge.20 .and. nproc.lt.30) then
        A1s = A1s + A1c;  A2s = A2s + A2c;  B1s = B1s + B1c
        Gfs = Gfs + Gfc
      endif
      if(mod(nproc,10) .eq. 2) then
        A1s = A1s + A1d;  A2s = A2s + A2d;  B1s = B1s + B1d
        Gfs = Gfs + Gfd
      endif
      ! for jet final state
      if(nproc.ge.10 .and. nproc.lt.20) then
        Djs = Djs + Djc*log(muren2/ptj/ptj/Rcone**2)
      endif
      if(mod(nproc,10) .eq. 1) then
        Djs = Djs + Djd*log(muren2/ptj/ptj/Rcone**2)
      endif

      ! define logarithm short-hands
      ! Lqcd2 << mufac2 <= mures2 <= Q2
      Lqr = log(muren2/mures2)
      Lqf = log(muren2/mufac2)
      Lql = log(muren2/Lqcd2)
      Lrf = log(mures2/mufac2)
      Lrl = log(mures2/Lqcd2)
      Lfl = log(mufac2/Lqcd2)
      if(Lrl.le.0.0 .or. Lfl.le.0.0) then  
        write(*,*) 'Lqcd too large!!!'
        write(*,*) ptj,ptc,ptd,qt
      endif
      LLrl = log(Lrl)
      LLfl = log(Lfl)

      ! perturbative Sudakov factor (LL:A1, NLL:A1,A2,B1)
      ! A1-term @ 1-loop a_s
      SP = A1s/b0 * (Lql*log(Lrl/Lfl)-Lrf)
      ! A1-term @ 2-loop a_s
     &   + A1s*b1/b0**3 *( !A1@2L
     &     + Lqr/Lrl*(1d0+LLrl) - Lqf/Lfl*(1d0+LLfl) 
     &     + log(Lrl/Lfl) + (LLrl**2 - LLfl**2)/2.0 )
      ! A2-term @ 1-loop a_s
     &   + A2s/b0**2 * (Lqf/Lfl-Lqr/Lrl-log(Lrl/Lfl))
      ! A2-term @ 2-loop a_s
c     &   + A2s*b1/b0**4/2d0 *(
c     &      + Lqr*(1d0+2d0*LLrl)/Lrl**2 - Lqf*(1d0+2d0*LLfl)/Lfl**2
c     &      - (3d0+2d0*LLrl)/Lrl + (3d0+2d0*LLfl)/Lfl )
c     &   + A2s*b1**2/b0**6/108d0 *(
c     &      + (19d0+30d0*LLrl+18d0*LLrl**2)/Lrl**2
c     &      - (19d0+30d0*LLfl+18d0*LLfl**2)/Lfl**2
c     &      - 4d0*Lqr*(2d0+6d0*LLrl+9d0*LLrl**2)/Lrl**3
c     &      + 4d0*Lqf*(2d0+6d0*LLfl+9d0*LLfl**2)/Lfl**3 )
      ! B1-term @ 1-loop a_s
     &   + (B1s+Djs)/b0 * log(Lrl/Lfl)
      ! B1-term @ 2-loop a_s
     &   + (B1s+Djs)*b1/b0**3 * ((1.0+LLrl)/Lrl-(1.0+LLfl)/Lfl)

      ! non-perturbative Sudakov factor
      SNP = 0.212*bp**2 + 0.84/2.0*log(muren2/2.4)*log(bp/bstar)
      SNP = SNP * Gfs/2.0

c      ! Non-perturbative Sudakov terms
c      if(NPschm.eq.1) then              ! DWS-g
c        SNP = bp**2 * ( g1 + g2*log(mures/2d0/Q0) )
c      elseif(NPschm.eq.2) then          ! LY-g
c        SNP = bp**2 * ( g1 + g2*log(mures/2d0/Q0) )
c     v       +bp * g1*g3*log(100d0*xa*xb)
c      elseif(NPschm.eq.3) then          ! BLNY-g
c        SNP = bp**2 * ( g1 + g2*log(mures/2d0/Q0) 
c     v       +g1*g3*log(100d0*xa*xb) )
c      elseif(NPschm.eq.4) then          ! SIYY-DY
c        SNP = g1*bp**2 + g2*log(bp/bstar)*log(mures/Q0)
c      else
c        write(*,*) "error in NP scheme configuration"
c        SNP = 0d0
c      endif

      ! medium induced broadening
c      Smed = 0.0
c      if(inmedium) then
        Smed = bp**2/4.0 * (Gqc+Gqd)*qhat
c      endif

      sudfac = exp(-SP -SNP -Smed)

      return 
      end function sudfac

      end module sudakov
