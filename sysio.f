      module sysio

      use, intrinsic :: iso_fortran_env, only: 
     &  compiler_version, compiler_options
      use comvar
      implicit none
      public

      contains

      function getu()
      integer getu,u
      logical ex
      do u=10,100
        inquire(unit=u,opened=ex)
        if(.not. ex) then
          getu = u
          return
        endif
      enddo
      end function getu

      subroutine readin
      integer :: u,stat
      u = getu()
      ! read data from input file
      open(u,iostat=stat,file=ifile,status='old')
      if(stat.ne.0) call exit(11)
      read(u,*) nproc
      read(u,*) CME
      read(u,*) ptTmin,ptTmax
      read(u,*) yTmin,yTmax
      read(u,*) ptAmin,ptAmax
      read(u,*) yAmin,yAmax
      read(u,*) qhat
      read(u,*) Rcone
      read(u,*) nfin,asloopin
      read(u,*) ncall1,itmax1
      read(u,*) ncall2,itmax2
      read(u,*) unitset
      read(u,*) qqch
      read(u,*) qgch
      read(u,*) ggch
      read(u,*) info
      read(u,*) tocons
      read(u,*) tofile
      read(u,*) doresum
      read(u,*) inmedium
      close(u)
      ! verify data integrity
      !if(CME.le.0d0) call exit(21)
      !if(ncall1.lt.100 .or. ncall2.lt.100) call exit(21)
      !if(CME.ge.huge(CME)) call exit(22)
      !if(ncall1.ge.huge(1) .or. ncall2.ge.huge(1)) call exit(22)
      !if(nfin.lt.3 .or. nfin.gt.6) call exit(23)
      !if(asloopin.lt.1 .or. asloopin.gt.3) call exit(23)
      ! set associate values
      CME2 = CME*CME
      ycmin = yTmin;                ycmax = yTmax
      ydmin = yAmin;                ydmax = yAmax
      ptjmin = 1d0;                 ptjmax = CME/2d0
      qtmin = 0d0;                  qtmax = ptjmax
      phiqmin = 0d0;                phiqmax = twoPI
      bpmin = 0d0;                  bpmax = 10d0
      xamin = 0d0;                  xamax = 1d0
      xbmin = 0d0;                  xbmax = 1d0
      zcmin = 0d0;                  zcmax = 1d0
      zdmin = 0d0;                  zdmax = 1d0
      geoxmin = -20d0;              geoxmax = +20d0
      geoymin = -20d0;              geoymax = +20d0
      phimin = 0d0;                 phimax = twoPI
      nprn = -1
      avgi = 0d0
      sd = 0d0
      chi2a = 0d0
      end subroutine readin

      subroutine printinfo(u)
      integer, intent(in) :: u
      write(*,*) 'program is compiled by:',compiler_version()
      write(*,*) 'compiler options are:',compiler_options()
      write(u,*) '*******************'
      write(u,*) 'resummation program'
      write(u,*) '  Raymond CL 2022  '
      write(u,*) '*******************'
      write(u,*) '===kinematics(GeV)==='
      write(u,*) 'center-of-mass energy:',CME
      write(u,*) 'trigger pt:',ptTmin,'<ptT<',ptTmax
      write(u,*) 'trigger rap:',yTmin,'<yT<',yTmax
      write(u,*) 'associate pt:',ptAmin,'<ptA<',ptAmax
      write(u,*) 'associate rap:',yAmin,'<yA<',yAmax
      write(u,*) 'qhat:',qhat
      write(u,*) 'Rcone:',Rcone
      write(u,*) nfin,'number of quark flavours'
      write(u,*) asloopin,'loop alpha_s expression'
      write(u,*) '===Vegas integration==='
      write(u,*) 'warm-up:',ncall1,'calls',itmax1,'iterations'
      write(u,*) 'end-run:',ncall2,'calls',itmax2,'iterations'
      write(u,*) '===flag settings==='
      if(unitset.eq.1) write(u,*) 'result in natural unit'
      if(unitset.eq.2) write(u,*) 'result in mb'
      if(unitset.eq.3) write(u,*) 'result in nb'
      if(unitset.eq.4) write(u,*) 'result in fb'
      write(u,*) 'qq-channel ',merge('opened','closed',qqch)
      write(u,*) 'qg-channel ',merge('opened','closed',qgch)
      write(u,*) 'gg-channel ',merge('opened','closed',ggch)
      if(tocons) write(u,*) 'results written to console'
      if(tofile) write(u,*) 'results written to file'
      if(doresum) write(u,*) 'performing resummation'
      if(inmedium) write(u,*) 'process in AA collision'
      write(u,*) new_line('a')
      end subroutine printinfo

      subroutine set_kinematics

      ! set trigger parton info (pt,phi)
      if(nproc/10 .ge. 4) then ! boson jet, qt contribute to ptd only
        ptc = ptj
        phic = phi
      else ! dijet, qt/2 contribute to both ptc+ptd
        ptc = ptj**2 + qt**2/4d0 + ptj*qt*cos(phiq) 
        ptc = merge(0.0,sqrt(ptc),ptc.le.0.0)
        phic = qt/ptc/2d0*sin(phiq)
        if(phic.le.-1.0) phic = -1.0
        if(phic.ge.+1.0) phic = +1.0
        phic = phi + asin(phic)
      endif

      ! set associate parton info (pt,phi)
      if(nproc/10 .ge. 4) then
        ptd = ptj**2 + qt**2 - 2d0*ptj*qt*cos(phiq)
        ptd = merge(0.0,sqrt(ptd),ptd.le.0.0)
        phid = qt/ptd*sin(phiq)
      else
        ptd = ptj**2 + qt**2/4d0 - ptj*qt*cos(phiq)
        ptd = merge(0.0,sqrt(ptd),ptd.le.0.0)
        phid = qt/ptd/2d0*sin(phiq)
      endif
      if(phid.le.-1.0) phid = -1.0
      if(phid.ge.+1.0) phid = +1.0
      phid = phi + PI - asin(phid)

      ! set hadron info (pt,rap)
      if(nproc/10 .eq. 2)      pthc = ptc*zc;  yhc = yc
      if(mod(nproc,10) .eq. 2) pthd = ptd*zd;  yhd = yd

      end subroutine set_kinematics

      end module sysio
