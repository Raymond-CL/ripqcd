      module pdf

      ! PDF module to return corresponding pdf matrix
      ! call setpdf() first that read pdf to memory buffer
      ! getpdf(): is the normal proton pdf
      ! getnpdf_s is the pdf with nuclear effect
      ! getnpdf_g is the npdf with geometry correction
      ! both npdf are polymorphic using getnpdf() interface

      use nrtype
      use glauber
      implicit none
      public
      character(len=40) pdfgrid

      interface getnpdf
        module procedure getnpdf_s,getnpdf_g
      end interface getnpdf

      contains

      subroutine setpdf()
        !pdfgrid = './data-grids/i2Tn3.00.pds'
        pdfgrid = './data-grids/CT14nlo_NF4.pds'
        !pdfgrid = './data-grids/CT14LN.pds'
        call setCT18(pdfgrid)
      end subroutine setpdf

      subroutine getpdf(sig,x,Q,g)
        integer(i4b), intent(in) :: sig
        real(sp), intent(in) :: x,Q
        real(sp), intent(out), dimension(-6:6) :: g
        real(dp) :: CT18Pdf
        integer(i4b) :: i
        do i=-6,6
          g(i) = 0d0
        enddo
        if(x.le.1e-9 .or. x.ge.1d0) return
        !if(Q.le.1.3d0 .or. Q.ge.1e5) return
        do i = -5,5 
          g(i) = CT18Pdf(i*sig,dble(x),dble(Q))
        enddo
      end subroutine getpdf

      subroutine getnpdf_s(sig,x,Q,ng)
      integer(i4b), intent(in) :: sig
      real(sp), intent(in) :: x,Q
      real(sp), intent(out), dimension(-6:6) :: ng
      real(sp), dimension(-6:6) :: g
      integer(i4b) :: ord,set,A,Z
      real(sp) :: pfac,nfac,ruv,rdv,ru,rd,rs,rc,rb,rg
      if(x.le.1e-7 .or. x.ge.1d0) return
      if(Q.le.1.3d0 .or. Q.ge.1e4) return
      call getpdf(sig,x,Q,g)
      A = 208
      Z = 82
      pfac = dble(Z)/A
      nfac = dble(A-Z)/A
      ord = 1
      set = 1
      call EPPS16(ord,set,A,x,Q,ruv,rdv,ru,rd,rs,rc,rb,rg)
      ng(-5) = rb*g(-5)
      ng(-4) = rc*g(-4)
      ng(-3) = rs*g(-3)
      ng(-2) = pfac * rd*g(-2) + nfac * ru*g(-1)
      ng(-1) = pfac * ru*g(-1) + nfac * rd*g(-2)
      ng( 0) = rg*g( 0)
      ng(+1) = pfac * (ruv*(g(+1)-g(-1)) + ru*g(-1))
     >        +nfac * (rdv*(g(+2)-g(-2)) + rd*g(-2))
      ng(+2) = pfac * (rdv*(g(+2)-g(-2)) + rd*g(-2))
     >        +nfac * (ruv*(g(+1)-g(-1)) + ru*g(-1))
      ng(+3) = rs*g(+3)
      ng(+4) = rc*g(+4)
      ng(+5) = rb*g(+5)
      end subroutine getnpdf_s

      subroutine getnpdf_g(sig,x,Q,geox,geoy,ng)
      integer(i4b), intent(in) :: sig
      real(sp), intent(in) :: x,Q
      real(sp), intent(in) :: geox,geoy
      real(sp), intent(out), dimension(-6:6) :: ng
      real(sp), dimension(-6:6) :: g
      integer(i4b) :: ord,set,A,Z
      real(sp) :: pfac,nfac,ruv,rdv,ru,rd,rs,rc,rb,rg
      real(sp) :: geofac
      if(x.le.1e-7 .or. x.ge.1d0) return
      if(Q.le.1.3d0 .or. Q.ge.1e4) return
      call getpdf(sig,x,Q,g)
      A = 208
      Z = 82
      pfac = dble(Z)/A
      nfac = dble(A-Z)/A
      ord = 1
      set = 1
      call EPPS16(ord,set,A,x,Q,ruv,rdv,ru,rd,rs,rc,rb,rg)
      geofac = dble(A)*TA(geox,geoy)/taa0
      ruv = 1d0+(ruv-1d0)*geofac
      rdv = 1d0+(rdv-1d0)*geofac
      ru  = 1d0+(ru -1d0)*geofac
      rd  = 1d0+(rd -1d0)*geofac
      rs  = 1d0+(rs -1d0)*geofac
      rc  = 1d0+(rc -1d0)*geofac
      rb  = 1d0+(rb -1d0)*geofac
      rg  = 1d0+(rg -1d0)*geofac     
      ng(-5) = rb*g(-5)
      ng(-4) = rc*g(-4)
      ng(-3) = rs*g(-3)
      ng(-2) = pfac * rd*g(-2) + nfac * ru*g(-1)
      ng(-1) = pfac * ru*g(-1) + nfac * rd*g(-2)
      ng( 0) = rg*g( 0)
      ng(+1) = pfac * (ruv*(g(+1)-g(-1)) + ru*g(-1))
     >        +nfac * (rdv*(g(+2)-g(-2)) + rd*g(-2))
      ng(+2) = pfac * (rdv*(g(+2)-g(-2)) + rd*g(-2))
     >        +nfac * (ruv*(g(+1)-g(-1)) + ru*g(-1))
      ng(+3) = rs*g(+3)
      ng(+4) = rc*g(+4)
      ng(+5) = rb*g(+5)
      end subroutine getnpdf_g
 
      end module pdf
