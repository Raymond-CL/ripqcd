      module comvar
      ! common (global) variables
      use nrtype
      implicit none
      public

      integer :: tempcounter,nancounter
      integer :: nproc

      ! basic kinematic variables:
      ! process: a + b -> c + d
      ! a: + moving (plus direction, +z momentum)
      ! b: - moving (minus direction, -z momentum)
      ! c: trigger particle
      ! d: associate particle

      ! center of mass energy
      real(sp) :: CME,CME2

      ! trigger and associate kinematics
      real(sp) :: ptTmin,ptTmax,yTmin,yTmax
      real(sp) :: ptAmin,ptAmax,yAmin,yAmax

      ! particle kinematics parton level
      real(sp) :: ptjmin,ptj,ptjmax,yjmin,yj,yjmax,phij,mj
      real(sp) :: ptcmin,ptc,ptcmax,ycmin,yc,ycmax,phic,mc
      real(sp) :: ptdmin,ptd,ptdmax,ydmin,yd,ydmax,phid,md

      ! hadron kinematics
      real(sp) :: pthcmin,pthc,pthcmax,yhcmin,yhc,yhcmax,phihc
      real(sp) :: pthdmin,pthd,pthdmax,yhdmin,yhd,yhdmax,phihd

      ! momentum fraction
      real(sp) :: xamin,xa,xamax
      real(sp) :: xbmin,xb,xbmax
      real(sp) :: zcmin,zc,zcmax
      real(sp) :: zdmin,zd,zdmax

      ! angles
      real(sp) :: dphimin,dphi,dphimax !azim-ang betwn c&d
      real(sp) :: thetamin,theta,thetamax

      ! imbalance
      real(sp) :: qtmin,qt,qtmax
      real(sp) :: phiqmin,phiq,phiqmax !azim-ang betwn qtnet&xaxis
      real(sp) :: bpmin,bp,bpmax
      real(sp) :: bstar

      ! Mandelstam variables
      real(sp) :: mans,mant,manu

      ! scales
      real(sp) :: muc,muc2,mud,mud2
      real(sp) :: mujet,mujet2,HT
      real(sp) :: mufac,mures,muren,mufrg
      real(sp) :: mufac2,mures2,muren2,mufrg2

      ! geometry (in AA collisions)
      real(sp) :: geoxmin,geox,geoxmax
      real(sp) :: geoymin,geoy,geoymax
      real(sp) :: phimin,phi,phimax !azim-ang betwn c&xaxis

      ! medium
      real(sp) :: qhat

      ! jet cone radius size
      real(sp) :: Rcone
      
      ! qcd set
      integer(i4b) :: nfin,asloopin

      ! vegas
      integer(i4b) :: ncall1,itmax1
      integer(i4b) :: ncall2,itmax2
      integer(i4b) :: init,ndim,nprn,dind
      real(sp) :: avgi,chi2a,sd
      real(sp), dimension(30) :: region

      ! unit conversion
      integer(i4b) :: unitset
      real(sp) :: unitfac
     
      ! io-file name 
      character(len=*), parameter :: ifile = 'input.dat'
      character(len=*), parameter :: ofile = 'output.dat'

      ! flags
      logical :: qqch         ! quark-quark channel
      logical :: qgch         ! quark-gluon channel
      logical :: ggch         ! gluon-gluon channel
      logical :: fill_hist    ! fill histogram
      logical :: info         ! print program info
      logical :: tocons       ! print to console
      logical :: tofile       ! print to file
      logical :: doresum      ! do resummation
      logical :: inmedium     ! include medium effect

      end module comvar
