!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CAO -> SAO (Cartesian -> real spherical harmonic) transformation of a
! shell-pair integral block. Supports s,p,d,f shells (li,lj = 0,1,2,3).
!
! Cartesian ordering (matches intpack.f/bas.f90):
!   d (6): xx,yy,zz,xy,xz,yz
!   f (10): xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz
! (f Cartesian components must already carry the sqrt(5)/sqrt(15)
!  reference-normalization factors applied in bas.f90.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine dtrf2(s,li,lj)
      implicit none
      real*8, intent(inout) :: s(10,10)
      integer,intent(in)    :: li,lj

      integer :: cartd(0:3), sphd(0:3)
      data cartd /1,3,6,10/
      data sphd  /1,3,5,7/

c --- d shell: CAO(xx,yy,zz,xy,xz,yz) -> SAO(dx2-y2,dz2,dxy,dxz,dyz)
      real*8 :: Td(6,5)
      parameter (Td = reshape((/
     . 0.5d0*sqrt(3.0d0),-0.5d0*sqrt(3.0d0),0.0d0,0.0d0,0.0d0,0.0d0,
     .-0.5d0,-0.5d0,1.0d0,0.0d0,0.0d0,0.0d0,
     . 0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,
     . 0.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,
     . 0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0 /), (/6,5/)))

c --- f shell: CAO(xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz) ->
c     SAO(fz3, fxz2, fyz2, fz(x2-y2), fxyz, fx(x2-3y2), fy(3x2-y2))
      real*8 :: Tf(10,7)
      parameter (Tf = reshape((/
     . 0.0d0,0.0d0,2.0d0,0.0d0,-3.0d0,0.0d0,-3.0d0,0.0d0,0.0d0,0.0d0,
     .-1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-1.0d0,0.0d0,4.0d0,0.0d0,0.0d0,
     . 0.0d0,-1.0d0,0.0d0,-1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,4.0d0,0.0d0,
     . 0.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,-1.0d0,0.0d0,0.0d0,0.0d0,
     . 0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,2.0d0,
     . 1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-3.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     . 0.0d0,-1.0d0,0.0d0,3.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0
     . /), (/10,7/)))

      real*8 :: tmp(10,10),out(10,10)
      integer :: nci,ncj,nsi,nsj

!     transformation not needed for pure s/p overlap -> do nothing
      if (li.lt.2.and.lj.lt.2) return

      nci = cartd(li)
      ncj = cartd(lj)
      nsi = sphd(li)
      nsj = sphd(lj)

c     transform "i" (column/ket) index: cartesian -> spherical
      tmp = 0.0d0
      if (li.eq.2) then
         tmp(1:ncj,1:nsi) = matmul(s(1:ncj,1:nci),Td(1:nci,1:nsi))
      else if (li.eq.3) then
         tmp(1:ncj,1:nsi) = matmul(s(1:ncj,1:nci),Tf(1:nci,1:nsi))
      else
         tmp(1:ncj,1:nsi) = s(1:ncj,1:nci)
      endif

c     transform "j" (row/bra) index: cartesian -> spherical
      out = 0.0d0
      if (lj.eq.2) then
         out(1:nsj,1:nsi) = matmul(transpose(Td(1:ncj,1:nsj)),
     .                              tmp(1:ncj,1:nsi))
      else if (lj.eq.3) then
         out(1:nsj,1:nsi) = matmul(transpose(Tf(1:ncj,1:nsj)),
     .                              tmp(1:ncj,1:nsi))
      else
         out(1:nsj,1:nsi) = tmp(1:ncj,1:nsi)
      endif

      s = 0.0d0
      s(1:nsj,1:nsi) = out(1:nsj,1:nsi)

      end subroutine dtrf2
