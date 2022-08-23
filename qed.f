      module qed

      ! QED module (with weak)
      ! includes chrg(q) for quark electric charge
      ! and tau3(q) for weak isospin (3rd component)

      use nrtype
      implicit none
      public
      ! alpha_e at Q_0 or m_Z
      real(sp), parameter :: AeQ0 = 1d0/137.036d0
      real(sp), parameter :: AemT = 1d0/133.472d0
      real(sp), parameter :: AemZ = 1d0/127.952d0
      ! weak boson masses
      real(sp), parameter :: mZ = 91.1876d0
      real(sp), parameter :: mW = 80.379d0
      real(sp), parameter :: mH = 125.10d0
      ! Fermi coupling constant
      real(sp), parameter :: Gf = 1.1663787e-5
      ! Weinberg angle constants
      real(sp), parameter :: cos2w = mW**2/mZ**2
      real(sp), parameter :: sin2w = 1d0 - cos2w
      real(sp), parameter :: cosw = sqrt(cos2w)
      real(sp), parameter :: sinw = sqrt(sin2w)
      real(sp), parameter :: thetaw = acos(cosw)
      real(sp), parameter :: gw = sqrt(AemZ*4d0*pi/sin2w)
      ! CKM matrix elements
      real(sp), parameter :: Vud = 0.97401d0
      real(sp), parameter :: Vus = 0.22650d0
      real(sp), parameter :: Vub = 0.00361d0
      real(sp), parameter :: Vcd = 0.22636d0
      real(sp), parameter :: Vcs = 0.97320d0
      real(sp), parameter :: Vcb = 0.04053d0
      real(sp), parameter :: Vtd = 0.00854d0
      real(sp), parameter :: Vts = 0.03978d0
      real(sp), parameter :: Vtb = 0.99917d0
      ! Cabibbo angle constants
      real(sp), parameter :: tan2c = Vus*Vus/Vud/Vud
      real(sp), parameter :: cos2c = 1d0/(1d0+tan2c)
      real(sp), parameter :: sin2c = 1d0/(1d0+1d0/tan2c)
      real(sp), parameter :: thetac = atan(sqrt(tan2c))
      ! Z/W to dilepton decay ratio
      real(sp), parameter :: Z2ee = 0.033632d0
      real(sp), parameter :: Z2uu = 0.033662d0
      real(sp), parameter :: Z2tt = 0.033696d0
      real(sp), parameter :: Z2ll = 0.033658d0
      real(sp), parameter :: W2ev = 0.1071d0
      real(sp), parameter :: W2uv = 0.1063d0
      real(sp), parameter :: W2tv = 0.1138d0
      real(sp), parameter :: W2lv = 0.1086d0

      contains

      function chrg(q)
        real(sp) :: chrg
        integer(i4b) :: q
        chrg = 0d0
        if(q.eq.0) return
        if(abs(q).eq.1) chrg = +2d0/3d0
        if(abs(q).eq.2) chrg = -1d0/3d0
        if(abs(q).eq.3) chrg = -1d0/3d0
        if(abs(q).eq.4) chrg = +2d0/3d0
        if(abs(q).eq.5) chrg = -1d0/3d0
        if(abs(q).eq.6) chrg = +2d0/3d0
        chrg = chrg * (q/abs(q))
        return
      end function chrg

      function tau3(q)
        real(sp) :: tau3
        integer(i4b) :: q
        tau3 = 0d0
        if(q.eq.0) return
        if(abs(q).eq.1) tau3 = +1d0/2d0
        if(abs(q).eq.2) tau3 = -1d0/2d0
        if(abs(q).eq.3) tau3 = -1d0/2d0
        if(abs(q).eq.4) tau3 = +1d0/2d0
        if(abs(q).eq.5) tau3 = -1d0/2d0
        if(abs(q).eq.6) tau3 = +1d0/2d0
        tau3 = tau3 * (q/abs(q))
        return
      end function tau3

      function coupqz(q)
        real(sp) :: coupqz
        integer(i4b) :: q
        real(sp) :: t3q,Vq,Aq,temp
        t3q = tau3(q)
        Vq = t3q - 2d0*chrg(q)*sin2w
        Aq = t3q
        temp = gw/2d0/cosw
        coupqz = temp**2 * (Vq**2 + Aq**2)
      end function coupqz

      function coupqp(q)
        real(sp) :: coupqp
        integer(i4b) :: q
        coupqp = AemZ * chrg(q) * chrg(q)
      end function coupqp

      end module qed
