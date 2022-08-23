      module ff

      ! sig<0 for neutral hadrons
      ! sig=0 for all hadrons
      ! sig>0 for charged hadrons

      use nrtype
      implicit none
      private
      public :: getff

      contains

      subroutine getff(sig,z,Q,d)
        integer(i4b), intent(in) :: sig
        real(sp), intent(inout) :: z,Q
        real(sp), intent(out), dimension(-6:6) :: d
        real(sp), dimension(0:10) :: dt
        integer(i4b) :: j
        !do i=0,10
        !  d(i) = 0.0
        !enddo
        if(z.le.0.05 .or. z.ge.1.0) z = 0.05
        !if(z.le.0.05 .or. z.ge.1.0) return
        if(Q.le.1.00 .or. Q.ge.1e3) return
c        if(sig.le.0) then
c          do i=5,8
c            call akk(i,dble(z),dble(Q),dt)
c            do j=0,10
c              d(j) = d(j) + dt(j)
c            enddo
c          enddo
c        endif
c        if(sig.ge.0) then
          call akk(sig,dble(z),dble(Q),dt)
          d(0) = dt(0)
          do j=1,5
            d(+j) = dt(2*j-1)
            d(-j) = dt( 2*j )
          enddo
c        endif
      end subroutine getff

      end module ff
