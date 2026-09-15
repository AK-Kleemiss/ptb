!! ------------------------------------------------------------------------
! Writes the AO density matrix, overlap matrix, and basis-set metadata
! needed to reconstruct rho(r) on an arbitrary real-space grid (e.g. to
! compare against a reference DFT density cube), in the PTB_DENMAT text
! format consumed by tools/lanthanide_fit/ptb_lnf/denmat.py.
!
! sint (ints.f90) already rescales S in-place by norm(i)*norm(j) before
! returning it (its final loop, easy to miss), so S/P as received here are
! already expressed w.r.t. the unit-normalized SAOs (S_ii = 1 exactly) --
! confirmed empirically (trace(P*S) = nel, and S's diagonal is exactly 1).
! xnorm(i) is still needed on the Python side: denmat.py's evaluate_density
! builds each AO from raw (unnormalized) cartesian primitive combinations
! and multiplies by ao.norm to match this same normalized-SAO convention.
! P and S themselves need no further rescaling here.
!! ------------------------------------------------------------------------
subroutine wrdenmat(fname,n,at,xyz,ndim,P,S,xnorm)
   use iso_fortran_env, only : wp => real64
   use bascom
   implicit none
   character(len=*), intent(in) :: fname
   integer, intent(in) :: n, at(n), ndim
   real(wp), intent(in) :: xyz(3,n)
   real(wp), intent(in) :: P(ndim*(ndim+1)/2)
   real(wp), intent(in) :: S(ndim*(ndim+1)/2)
   real(wp), intent(in) :: xnorm(ndim)

   integer :: lladr(0:3)
   data lladr  /1,3,6,10/

   ! Cartesian power (px,py,pz) per component, by shell type. Ordering
   ! matches bas.f90/dtrf2.f: s=(0,0,0); p=x,y,z; d=xx,yy,zz,xy,xz,yz;
   ! f=xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz.
   integer :: powS(3,1), powP(3,3), powD(3,6), powF(3,10)
   data powS /0,0,0/
   data powP /1,0,0, 0,1,0, 0,0,1/
   data powD /2,0,0, 0,2,0, 0,0,2, 1,1,0, 1,0,1, 0,1,1/
   data powF /3,0,0, 0,3,0, 0,0,3, 2,1,0, 2,0,1, 1,2,0, 0,2,1, 1,0,2, 0,1,2, 1,1,1/

   ! CAO->SAO transform coefficients, copied verbatim from dtrf2.f so the
   ! exported basis matches exactly what PTB itself used to build S and P.
   real(wp) :: Td(6,5), Tf(10,7)
   parameter (Td = reshape((/ &
    0.5d0*sqrt(3.0d0),-0.5d0*sqrt(3.0d0),0.0d0,0.0d0,0.0d0,0.0d0, &
   -0.5d0,-0.5d0,1.0d0,0.0d0,0.0d0,0.0d0, &
    0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0, &
    0.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0, &
    0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0 /), (/6,5/)))
   parameter (Tf = reshape((/ &
    0.0d0,0.0d0,2.0d0,0.0d0,-3.0d0,0.0d0,-3.0d0,0.0d0,0.0d0,0.0d0, &
   -1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-1.0d0,0.0d0,4.0d0,0.0d0,0.0d0, &
    0.0d0,-1.0d0,0.0d0,-1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,4.0d0,0.0d0, &
    0.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,-1.0d0,0.0d0,0.0d0,0.0d0, &
    0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,2.0d0, &
    1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,-3.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
    0.0d0,-1.0d0,0.0d0,3.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0 &
    /), (/10,7/)))

   integer :: iu, i, j, ia, iat, ish, l, m, c, nterms, np, pr, icao, iprim, lin
   integer :: sao, ijidx, ccomp
   real(wp) :: transform

   iu = 87
   open(unit=iu,file=fname)
   write(iu,'(a)') '# PTB_DENMAT 1'
   write(iu,'(a,1x,i0)') 'NATOMS ', n
   do i=1,n
      write(iu,'(2i6,i6,3f20.12)') i, at(i), at(i), xyz(1,i), xyz(2,i), xyz(3,i)
   enddo
   write(iu,'(a,1x,i0)') 'NDIM ', ndim

   do sao=1,ndim
      ia  = aoat(sao)
      iat = at(ia)
      ish = shell2ao(sao)
      l   = bas_lsh(ish,iat)
      m   = sao - aoshell(ish,ia)
      icao = caoshell(ish,ia)
      np = bas_npr(ish,iat)
      ! s/p: cartesian and spherical bases coincide 1:1 (no CAO->SAO mixing,
      ! dtrf2 itself is a no-op for li,lj < 2), so there's exactly one
      ! cartesian component -- the m'th one -- not lladr(l) of them.
      if (l.lt.2) then
         nterms = 1
      else
         nterms = lladr(l)
      endif
      write(iu,'(a,6i8,f20.12)') 'AO ', sao, ia, ish, l, m, nterms, xnorm(sao)
      do c=1,nterms
         if (l.lt.2) then
            ccomp = m
            transform = 1.0_wp
         else if (l.eq.2) then
            ccomp = c
            transform = Td(c,m)
         else
            ccomp = c
            transform = Tf(c,m)
         endif
         if (l.eq.0) then
            write(iu,'(a,3i4,i6,f20.12)') 'TERM ', powS(1,ccomp),powS(2,ccomp),powS(3,ccomp), np, transform
         else if (l.eq.1) then
            write(iu,'(a,3i4,i6,f20.12)') 'TERM ', powP(1,ccomp),powP(2,ccomp),powP(3,ccomp), np, transform
         else if (l.eq.2) then
            write(iu,'(a,3i4,i6,f20.12)') 'TERM ', powD(1,ccomp),powD(2,ccomp),powD(3,ccomp), np, transform
         else
            write(iu,'(a,3i4,i6,f20.12)') 'TERM ', powF(1,ccomp),powF(2,ccomp),powF(3,ccomp), np, transform
         endif
         do pr=1,np
            iprim = pr + prim_count(icao+ccomp)
            write(iu,'(2e24.14)') prim_exp(iprim), prim_cnt(iprim)
         enddo
      enddo
   enddo

   write(iu,'(a,1x,i0)') 'DENSITY_PACKED ', ndim*(ndim+1)/2
   do i=1,ndim
      do j=1,i
         ijidx = lin(i,j)
         write(iu,'(e24.14)') P(ijidx)
      enddo
   enddo

   write(iu,'(a,1x,i0)') 'OVERLAP_PACKED ', ndim*(ndim+1)/2
   do i=1,ndim
      do j=1,i
         ijidx = lin(i,j)
         write(iu,'(e24.14)') S(ijidx)
      enddo
   enddo

   close(iu)
end subroutine wrdenmat
