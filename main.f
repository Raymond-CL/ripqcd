*234567890**********1234567890**********1234567890**********1234567890**

      program main
      use, intrinsic :: iso_fortran_env, only: 
     &  compiler_version, compiler_options, stdout=>output_unit
      use nrtype
      use nr
      use comvar
      use sysio
      use bookkeep
      use glauber
      use qed
      use qcd
      use pdf
      use sudakov
      implicit none
      real :: tstart,tstop
      integer :: u
      interface
        function fxn(pt,wgt)
        use nrtype
        implicit none
        real(sp), dimension(:), intent(in) :: pt
        real(sp), intent(in) :: wgt
        real(sp) :: fxn
        end function fxn
      end interface
      call cpu_time(tstart)

      ! initialize
      call readin()
      if(info) call printinfo(stdout)
      call setqcd(nfin,asloopin)
      call setsud()
      call setpdf()
      if(inmedium) call setgeo(208)
      if(tofile) then 
        u = getu()
        open(u,file=ofile,status='replace')
        close(u)
      endif
      tempcounter = 0 ! debug counter

      ! set histogram
      call initbook()
      call newbook(1, .false., 'ds/dphi:ptc=[ 75,100]', 64, PI/2.0, PI)
      call newbook(2, .false., 'ds/dphi:ptc=[100,130]', 64, PI/2.0, PI)
      call newbook(3, .false., 'ds/dphi:ptc=[130,180]', 64, PI/2.0, PI)
      call newbook(4, .false., 'ds/dphi:ptc=[180,980]', 64, PI/2.0, PI)

      ! define vegas integration dimensions
      ndim = 3 ! for LO: only yc,yd,ptj
      if(doresum) ndim = ndim + 3 ! for resum: add qt,phiq,bp
      if(nproc/10 .eq. 2) ndim = ndim + 1 ! for c hadronic: add zc
      if(mod(nproc,10) .eq. 2) ndim = ndim + 1 ! for d hadronic: add zd
      if(inmedium) ndim = ndim + 3! for geometry: add geox,geoy,phi
      ! for e-loss: add de
      ! for smear: add gpeak,gwidth
      ! for Drell-Yan: maybe one less dimension?

      ! define vegas integration UL-limits
      dind = 1
      region(dind) = ycmin;  region(ndim+dind) = ycmax;  dind=dind+1
      region(dind) = ydmin;  region(ndim+dind) = ydmax;  dind=dind+1
      region(dind) = ptjmin; region(ndim+dind) = ptjmax; dind=dind+1
      if(doresum) then
        region(dind) = qtmin;   region(ndim+dind) = qtmax;   dind=dind+1
        region(dind) = phiqmin; region(ndim+dind) = phiqmax; dind=dind+1
        region(dind) = bpmin;   region(ndim+dind) = bpmax;   dind=dind+1
      endif
      if(nproc/10 .eq. 2) then
        region(dind) = zcmin; region(ndim+dind) = zcmax; dind=dind+1
      endif
      if(mod(nproc,10) .eq. 2) then
        region(dind) = zdmin; region(ndim+dind) = zdmax; dind=dind+1
      endif
      if(inmedium) then
        region(dind) = geoxmin; region(ndim+dind) = geoxmax; dind=dind+1
        region(dind) = geoymin; region(ndim+dind) = geoymax; dind=dind+1
        region(dind) = phimin;  region(ndim+dind) = phimax;  dind=dind+1
      endif
      if(dind-1 .ne. ndim) call exit(32)

      ! warm-up vegas grid, don't fill histogram
      fill_hist = .false.
      init = -1
      call vegas(region(1:2*ndim),fxn,init,ncall1,
     >  itmax1,nprn,avgi,sd,chi2a)
      if(tocons) then
        u = stdout
        write(u,*) 'total cross-section (warm-up):'
        write(u,*) avgi,'+-',sd,new_line('a')
      endif
      if(tofile) then
        u = getu()
        open(unit=u,file=ofile,position='append')
        write(u,*) 'total cross-section (warm-up):'
        write(u,*) avgi,'+-',sd,new_line('a')
        close(u)
      endif

      ! end-run vegas, fill histogram
      fill_hist = .true.
      init = +1
      call vegas(region(1:2*ndim),fxn,init,ncall2,
     >  itmax2,nprn,avgi,sd,chi2a)
      if(tocons) then
        u = stdout
        write(u,*) 'total cross-section (end-run):'
        write(u,*) avgi,'+-',sd,new_line('a')
      endif
      if(tofile) then
        u = getu()
        open(unit=u,file=ofile,position='append')
        write(u,*) 'total cross-section (end-run):'
        write(u,*) avgi,'+-',sd,new_line('a')
        close(u)
      endif

      ! print histogram
      if(tocons) then
        u = stdout
        call printbook(u)
      endif
      if(tofile) then
        u = getu()
        open(unit=u,file=ofile,access='append')
        call printbook(u)
        close(u)
      endif

      ! temporary counter for debugging
      !write(*,*) 'temp counter=',tempcounter

      call cpu_time(tstop)
      write(*,*) 'time elapsed: ',tstop-tstart,'seconds.'
      end program main

*234567890**********1234567890**********1234567890**********1234567890**

      function fxn(pt,wgt)
      use nrtype
      use comvar
      use bookkeep
      use sysio
      use glauber
      use pdf
      use ff
      use qed
      use qcd
      use sudakov
      implicit none
      real(sp), dimension(:), intent(in) :: pt
      real(sp), intent(in) :: wgt
      real(sp) :: fxn
      integer :: ptind
      real(sp) :: hard_dijet
      !real(sp) :: Etrig

***********************************
* setting kinematic variables
* try not to touch this part
***********************************

      ! initialize some default values
      fxn = 0.0
      qt = 0.0; phiq=0.0; bp=0.0
      zc = 1.0; zd = 1.0
      geox=0.0; geoy=0.0; phi=0.0

      ! vegas input from random vector
      ptind = 1
      yc = pt(ptind);  ptind=ptind+1
      yd = pt(ptind);  ptind=ptind+1
      ptj = pt(ptind); ptind=ptind+1
      if(doresum) then
        qt = pt(ptind);   ptind=ptind+1
        phiq = pt(ptind); ptind=ptind+1
        bp = pt(ptind);   ptind=ptind+1
      endif
      if(nproc/10 .eq. 2) then
        zc = pt(ptind); ptind=ptind+1
      endif
      if(mod(nproc,10) .eq. 2) then
        zd = pt(ptind); ptind=ptind+1
      endif
      if(inmedium) then
        geox = pt(ptind); ptind=ptind+1
        geoy = pt(ptind); ptind=ptind+1
        phi = pt(ptind);  ptind=ptind+1
      endif

      call set_kinematics()

      ! some useful partonic scale
      mc = 0.0; md = 0.0 ! massless unless for other processes
      if(nproc/10 .eq. 5) mc = mZ
      if(nproc/10 .eq. 6) mc = mW
      if(nproc/10 .eq. 7) mc = mH
      muc2 = ptc*ptc + mc*mc
      mud2 = ptd*ptd + md*md
      muc = sqrt(muc2)
      mud = sqrt(mud2)
      mujet = max(ptc,ptd) !ptj
      mujet2 = mujet*mujet

      ! set momentum fraction
      !xa = (muc*exp(+yc) + mud*exp(+yd)) / CME
      !xb = (muc*exp(-yc) + mud*exp(-yd)) / CME
      xa = mujet * (exp(+yc) + exp(+yd)) / CME
      xb = mujet * (exp(-yc) + exp(-yd)) / CME

      ! set Mandelstam variable
      mans = +xa*xb * CME2
      !mant = -xa * CME * muc * exp(-yc)
      !manu = -xa * CME * mud * exp(-yd)
      mant = -xa * CME * mujet * exp(-yc)
      manu = -xa * CME * mujet * exp(-yd)

      ! setting some important scales below
      ! renormalization scale
      muren2 = mans
      muren = sqrt(muren2)
      ! factorization scale
      bstar = bp/sqrt(1.0 + bp*bp/bmax/bmax)
      if(bstar.le.0.0) return
      mufac = 2.0*exp(-EULER)/bstar
      if(mufac.ge.1e+5) return
      mufac2 = mufac*mufac
      ! resummation scale (same as muren, could be different)
      mures = muren
      mures2 = mures*mures
      ! fragmentation scale (same as mufac, or pth)
      mufrg = merge(1.0,mufac,mufac.lt.1.0)
      mufrg2 = mufrg*mufrg

***********************************
* end of setting kinematics
* now we apply cuts suitable to experiment
***********************************

      ! apply kinematic and experimental cuts below
      ! tune to maximize MC integration efficiency
      ! momentum fraction cuts
      if(xa.le.xamin .or. xb.le.xbmin) return !does nothing
      if(xa.gt.xamax .or. xb.gt.xbmax) return

      ! large log(ptj/qt) to preserve applicability of resummation
      !if(qt.ge.ptj) return

      ! kinematic cut (trig and asso, parton/hadron)
      if(nproc/10 .eq. 2) then
        !Etrig = pthc*cosh(yhc)
        !if(Etrig.lt.ptTmin .or. Etrig.gt.ptTmax) return
        if(pthc.lt.ptTmin .or. pthc.gt.ptTmax) return
        if(yhc.lt.yTmin .or. yhc.gt.yTmax) return
      else
        if(ptc.lt.ptTmin .or. ptc.gt.ptTmax) return
        if(yc.lt.yTmin .or. yc.gt.yTmax) return
      endif
      if(mod(nproc,10) .eq. 2) then
        if(pthd.le.ptAmin .or. pthd.gt.ptAmax) return
        if(yhd.le.yAmin .or. yhd.gt.yAmax) return
      else
        if(ptd.le.ptAmin .or. ptd.gt.ptAmax) return
        if(yd.le.yAmin .or. yd.gt.yAmax) return
      endif

***********************************
* end of cuts
* calling hard kernel that calculates matrix elements
***********************************

      if(nproc.lt.40) then
        fxn = hard_dijet()
      elseif(nproc.gt.40 .and. mod(nproc,10).ne.0) then
        !fxn = hard_bosjet()
      elseif(nproc.gt.40 .and. mod(nproc,10).eq.0) then
        !fxn = hard_DY()
      else
        !fxn = hard_other()
      endif

***********************************
* apply overall normalization
***********************************

      if(doresum) fxn = fxn * qt*bp * bessel_jn(0,qt*bp) / twoPI
      if(inmedium) fxn = fxn * TAA(geox,geoy,0.0) / taa0 / twoPI

      ! unit conversion (38938573)
      if(unitset.eq.1) unitfac = 1.0
      if(unitset.eq.2) unitfac = 0.389379372E-3 * 1E3
      if(unitset.eq.3) unitfac = 0.389379372E-3 * 1E9
      if(unitset.eq.4) unitfac = 0.389379372E-3 * 1E15
      fxn = fxn * unitfac
     
      ! debugging options (nan) (before filling histogram)
      ! compile without -ffast-math option
      if(isnan(fxn)) then
        write(*,*) 'nan error:',fxn
        write(*,*) 'input:',pt(1),pt(2),pt(3)
        write(*,*) 'res:',pt(4),pt(5),pt(6)
        write(*,*) 'c parton:',ptc,phic
        write(*,*) 'd parton:',ptd,phid
        write(*,*) 'momfrac:',xa,xb,zc,zd
        write(*,*) 'mandelstam:',mans,mant,manu
        write(*,*) 'scales:',muren,mufac,mures,mufrg
      endif

***********************************
* define observable and fill histogram
***********************************

      ! define observable
      dphi = abs(phid-phic)
      dphi = merge(twoPI-dphi,dphi,dphi.gt.PI) ! folded dphi
      if( abs(dphi-PI) .le. epsilon(dphi)) dphi = PI-epsilon(dphi)
      if( abs(dphi-PIo2) .le. epsilon(dphi)) dphi = PIo2+epsilon(dphi)

      ! fill histogram
      if(fill_hist) then

        if(ptc.ge.75.0 .and. ptc.lt.100.0) then
          call fillbook(1, dphi, fxn*wgt/itmax2)
        elseif(ptc.ge.100.0 .and. ptc.lt.130.0) then
          call fillbook(2, dphi, fxn*wgt/itmax2)
        elseif(ptc.ge.130.0 .and. ptc.lt.180.0) then
          call fillbook(3, dphi, fxn*wgt/itmax2)
        elseif(ptc.ge.180.0) then
          call fillbook(4, dphi, fxn*wgt/itmax2)
        endif

      endif

      ! return differential cross-section
      return
      end function fxn 

*234567890**********1234567890**********1234567890**********1234567890**
