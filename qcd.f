      module qcd

      ! QCD module
      ! initialize by setting Nf and qcdloop, then call setqcd() 
      ! includes alpha_s 

      use nrtype
      implicit none
      public
      integer(i4b), parameter :: Nc = 3
      real(sp), parameter :: Ca = dble(Nc)
      real(sp), parameter :: Cf = (dble(Nc)**2-1d0)/2d0/dble(Nc)
      real(sp), parameter :: Tr = 1d0/2d0
      real(sp), parameter :: AsmZ = 0.1181d0
      integer(i4b), public :: Nf,qcdloop
      real(sp), protected :: b0,b1,b2
      real(sp), public :: Lqcd,Lqcd2

      contains

      subroutine setqcd(x,y)
      integer(i4b), intent(in) :: x,y
      Nf = x
      qcdloop = y
      b0 = (33d0-2d0*dble(Nf))/(12d0*PI)
      b1 = (153d0-19d0*dble(Nf))/(24d0*PI2)
      b2 = (2857d0-5033d0/9d0*dble(Nf)+325d0/27d0*dble(Nf)**2)
     >  /(128d0*PI3)
      Lqcd = 0d0
      if(qcdloop .eq. 1) then
        if(Nf.eq.3) Lqcd = 0.247d0
        if(Nf.eq.4) Lqcd = 0.155d0
        if(Nf.eq.5) Lqcd = 0.0884d0
        if(Nf.eq.6) Lqcd = 0.0456d0
      elseif(qcdloop .eq. 2) then
        if(Nf.eq.3) Lqcd = 0.741d0
        if(Nf.eq.4) Lqcd = 0.438d0
        if(Nf.eq.5) Lqcd = 0.228d0
        if(Nf.eq.6) Lqcd = 0.0987d0
      elseif(qcdloop .eq. 3) then
        if(Nf.eq.3) Lqcd = 0.641d0
        if(Nf.eq.4) Lqcd = 0.389d0
        if(Nf.eq.5) Lqcd = 0.210d0
        if(Nf.eq.6) Lqcd = 0.0953d0
      endif
      Lqcd2 = Lqcd**2
      end subroutine setqcd

      function alphas(q)
      real(sp) :: alphas
      real(sp), intent(in) :: q
      real(sp) :: as,q2,l2,t,l
      q2 = q**2
      l2 = lqcd**2
      t = log(q2/l2)
      l = log(t)
      as = 1d0/(b0*t)
      if(qcdloop.ge.2) then
        as = as * -b1*l/(b0**2*t)
      endif
      if(qcdloop.ge.3) then
        as = as * +(b1**2*(l**2-l-1d0)+b0*b2)/(b0**4*t**2)
      endif
      alphas = as
      return
      end function alphas

      end module qcd
