module density_export
      use iso_fortran_env, only : wp => real64
      use bascom
      implicit none
      private
      public :: write_density_export

contains

subroutine write_density_export(fname,n,ndim,at,xyz,P,S,xnorm)
      implicit none
      character(len=*), intent(in) :: fname
      integer, intent(in) :: n,ndim
      integer, intent(in) :: at(n)
      real(wp), intent(in) :: xyz(3,n)
      real(wp), intent(in) :: P(ndim*(ndim+1)/2)
      real(wp), intent(in) :: S(ndim*(ndim+1)/2)
      real(wp), intent(in) :: xnorm(ndim)

      integer :: io,iat,ish,iao,isao,icao,l,msph
      integer :: sph_count(0:3)

      data sph_count /1,3,5,7/

      io = 271
      open(unit=io,file=trim(fname),status='replace',action='write')

      write(io,'(a)') '# PTB_DENMAT 1'
      write(io,'(a,1x,i0)') 'NATOMS', n
      do iat=1,n
         write(io,'(a,1x,i0,1x,i0,3(1x,es24.16))') 'ATOM',iat,at(iat),xyz(1:3,iat)
      enddo

      write(io,'(a,1x,i0)') 'NDIM', ndim
      isao = 0
      do iat=1,n
         do ish=1,bas_nsh(at(iat))
            l = bas_lsh(ish,at(iat))
            icao = caoshell(ish,iat)
            do msph=1,sph_count(l)
               isao = isao + 1
               call write_ao(io,isao,iat,ish,l,msph,icao,xnorm(isao))
            enddo
         enddo
      enddo

      write(io,'(a,1x,i0)') 'DENSITY_PACKED', ndim*(ndim+1)/2
      call write_packed(io,ndim,P)
      write(io,'(a,1x,i0)') 'OVERLAP_PACKED', ndim*(ndim+1)/2
      call write_packed(io,ndim,S)
      write(io,'(a)') 'END'
      close(io)
end subroutine write_density_export

subroutine write_ao(io,isao,iat,ish,l,msph,icao,norm)
      implicit none
      integer, intent(in) :: io,isao,iat,ish,l,msph,icao
      real(wp), intent(in) :: norm

      integer :: t,nterms,nemit,cart_idx,px,py,pz,ipr,iprim
      real(wp) :: coeff
      real(wp) :: ccoef(10)

      call spherical_terms(l,msph,ccoef,nterms)
      nemit = 0
      do t=1,nterms
         if(abs(ccoef(t)).ge.1.0d-14) nemit = nemit + 1
      enddo
      write(io,'(a,6(1x,i0),1x,es24.16)') 'AO',isao,iat,ish,l,msph,nemit,norm

      do t=1,nterms
         cart_idx = t
         coeff = ccoef(t)
         if(abs(coeff).lt.1.0d-14) cycle
         call cart_powers(l,cart_idx,px,py,pz)
         write(io,'(a,4(1x,i0),1x,es24.16)') 'TERM',px,py,pz,prim_npr(icao+cart_idx),coeff
         do ipr=1,prim_npr(icao+cart_idx)
            iprim = prim_count(icao+cart_idx) + ipr
            write(io,'(2(1x,es24.16))') prim_exp(iprim), prim_cnt(iprim)
         enddo
      enddo
end subroutine write_ao

subroutine write_packed(io,ndim,mat)
      implicit none
      integer, intent(in) :: io,ndim
      real(wp), intent(in) :: mat(ndim*(ndim+1)/2)
      integer :: i

      do i=1,ndim*(ndim+1)/2
         write(io,'(es24.16)') mat(i)
      enddo
end subroutine write_packed

subroutine spherical_terms(l,msph,ccoef,nterms)
      implicit none
      integer, intent(in) :: l,msph
      real(wp), intent(out) :: ccoef(10)
      integer, intent(out) :: nterms

      real(wp), parameter :: sq3 = sqrt(3.0_wp)

      ccoef = 0.0_wp
      select case(l)
      case(0)
         nterms = 1
         ccoef(1) = 1.0_wp
      case(1)
         nterms = 3
         ccoef(msph) = 1.0_wp
      case(2)
         nterms = 6
         select case(msph)
         case(1) ! dx2-y2
            ccoef(1) = 0.5_wp*sq3
            ccoef(2) =-0.5_wp*sq3
         case(2) ! dz2
            ccoef(1) =-0.5_wp
            ccoef(2) =-0.5_wp
            ccoef(3) = 1.0_wp
         case(3) ! dxy
            ccoef(4) = 1.0_wp
         case(4) ! dxz
            ccoef(5) = 1.0_wp
         case(5) ! dyz
            ccoef(6) = 1.0_wp
         end select
      case(3)
         nterms = 10
         select case(msph)
         case(1) ! fz3
            ccoef(3) = 2.0_wp
            ccoef(5) =-3.0_wp
            ccoef(7) =-3.0_wp
         case(2) ! fxz2
            ccoef(1) =-1.0_wp
            ccoef(6) =-1.0_wp
            ccoef(8) = 4.0_wp
         case(3) ! fyz2
            ccoef(2) =-1.0_wp
            ccoef(4) =-1.0_wp
            ccoef(9) = 4.0_wp
         case(4) ! fz(x2-y2)
            ccoef(5) = 1.0_wp
            ccoef(7) =-1.0_wp
         case(5) ! fxyz
            ccoef(10)= 2.0_wp
         case(6) ! fx(x2-3y2)
            ccoef(1) = 1.0_wp
            ccoef(6) =-3.0_wp
         case(7) ! fy(3x2-y2)
            ccoef(2) =-1.0_wp
            ccoef(4) = 3.0_wp
         end select
      end select
end subroutine spherical_terms

subroutine cart_powers(l,idx,px,py,pz)
      implicit none
      integer, intent(in) :: l,idx
      integer, intent(out) :: px,py,pz

      px = 0
      py = 0
      pz = 0
      select case(l)
      case(0)
         return
      case(1)
         if(idx.eq.1) px = 1
         if(idx.eq.2) py = 1
         if(idx.eq.3) pz = 1
      case(2)
         if(idx.eq.1) px = 2
         if(idx.eq.2) py = 2
         if(idx.eq.3) pz = 2
         if(idx.eq.4) then; px = 1; py = 1; endif
         if(idx.eq.5) then; px = 1; pz = 1; endif
         if(idx.eq.6) then; py = 1; pz = 1; endif
      case(3)
         if(idx.eq.1) px = 3
         if(idx.eq.2) py = 3
         if(idx.eq.3) pz = 3
         if(idx.eq.4) then; px = 2; py = 1; endif
         if(idx.eq.5) then; px = 2; pz = 1; endif
         if(idx.eq.6) then; px = 1; py = 2; endif
         if(idx.eq.7) then; py = 2; pz = 1; endif
         if(idx.eq.8) then; px = 1; pz = 2; endif
         if(idx.eq.9) then; py = 1; pz = 2; endif
         if(idx.eq.10) then; px = 1; py = 1; pz = 1; endif
      end select
end subroutine cart_powers

end module density_export
