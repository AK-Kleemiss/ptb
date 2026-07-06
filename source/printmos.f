ccccccccccccccccccccccccccccccccccccccccccc
!    write out unformatted sTDA input     c
ccccccccccccccccccccccccccccccccccccccccccc           
! ncent  : # atoms
! nmo    : # MOs
! nbf    : # AOs
! nprims : # primitives (in total)
! xyz(4,ncent) : Cartesian coordinates & nuclear charge
! cont(nprims) : contraction coefficients of primitives
! alp(nprims) : exponents of primitives
! cmo(nbf,nmo) : LCAO-MO coefficients 
! eval(nmo)    : orbital eigenvalues
! occ(nmo)     : occupation # of MO
! ipty(nprims) : angular momentum of primitive function
! ipao(nbf)    : # primitives in contracted AO 
! ibf(ncent)   : # of contracted AOs on atom

      subroutine printmos(nc,at,xyz,nmo,homo,norm_,mowrcut,eval,occ,tmp)

      use bascom
      implicit none

      integer, intent ( in ) :: nmo,nc,at(nc),homo
      real*8,  intent ( in ) :: xyz(3,nc)
      real*8,  intent ( in ) :: norm_(nmo)
      real*8,  intent ( in ) :: mowrcut
      real*8,  intent ( in ) :: eval(nmo)
      real*8,  intent ( in ) :: occ(nmo)
      real*8,  intent ( inout ) :: tmp(nmo,nmo)

      ! temporary variables
      integer nbf
      integer i,j,k,nprims,nmomax,iat,iwfn,iao
      real*8 dum
      character*2 atyp
      integer lao(ncao),nprim(ncao),aoatcart(ncao)
      integer lladr(0:3),ll(0:3)
      data lladr  /1,3,6,10/
      data ll     /0,1,4,10/
      real*8,allocatable :: cmo(:,:)

      allocate(cmo(ncao,nmo))

      !Check for the sizes of matrizes
      if (size(tmp, 2) /= size(norm_, 1)) then
         error stop "Error: Dimensions do not match for multiplication."
      endif
      if (size(tmp, 1) /= nmo) then
         error stop "Error: Dimensions do not match for multiplication."
      endif

      do i=1,nmo
         tmp(i,:)=tmp(i,:)*norm_(i)
      enddo
      call sao2cao(nmo,tmp,cmo,nc,at)

      nbf=ncao
      nprim=prim_npr
      nprims=npr

      iwfn=29
      open(unit=iwfn,file='wfn.xtb',form='unformatted',
     .     status='replace')

! only print out virtuals below cutoff
      nmomax=nmo

                    !***********
                    ! RHF case *
                    !***********
! write dimensions
      write(iwfn)1
      write(iwfn)nc,nbf,nmomax,nprims
! now write coordinates & atom symbol
      do i = 1,nc
         call aasym(at(i),atyp)
         write(iwfn) atyp
      enddo

      do i = 1,nc
         do j=1,3
            dum=xyz(j,i)
            write(iwfn) dum
         enddo
         write(iwfn) at(i)
      enddo
! Now print basis set data

      k=0
      do i=1, nc
         iat=at(i)
         do j=1, bas_nsh(iat)
            do iao=1,lladr(bas_lsh(j,iat))
               k=k+1
               lao(k)=ll(bas_lsh(j,iat))+iao
               aoatcart(k)=i
            enddo
         enddo
      enddo

! print ipty
      do i=1,nbf
         k = lao(i)
         do j=1,nprim(i)
            write(iwfn) k
         enddo
      enddo
! iaoat
      do i=1,nbf
         k=aoatcart(i)
         do j=1,nprim(i)
            write(iwfn) k
         enddo
      enddo
! ipao
      do i=1,nbf
         k=i
         do j=1,nprim(i)
            write(iwfn) k
         enddo
      enddo

! exponents and coefficients
      write(iwfn) prim_exp(1:nprims)
      write(iwfn) prim_cnt(1:nprims)

! now the mo data
      write(iwfn) occ(1:nmomax)
      write(iwfn) eval(1:nmomax)
      write(iwfn) cmo(1:nbf,1:nmomax)
      close(iwfn)

      return
      end

C     *****************************************************************         

      subroutine AASYM(I,asy)
      integer, intent ( in ) :: i
      CHARACTER*2 ASY
      CHARACTER*2 ELEMNT(107), AS
      DATA ELEMNT/'h ','he',
     1 'li','be','b ','c ','n ','o ','f ','ne',
     2 'na','mg','al','si','p ','s ','cl','ar',
     3 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     4 'zn','ga','ge','as','se','br','kr',
     5 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     6 'cd','in','sn','sb','te','i ','xe',
     7 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     8 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
     9 'au','hg','tl','pb','bi','po','at','rn',
     1 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx',
     2 'fm','md','cb','xx','xx','xx','xx','xx'/
      AS=ELEMNT(I)
      CALL UPPER(AS)
      ASY=AS
      if(i.eq.103) asy='XX'
      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c transforms sao(5d) integrals to cao(6d) basis 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sao2cao(nbf,s,x,ncent,at)
      use bascom
      implicit none

      integer nbf,new,ncent,at(ncent)
      real*8  s(nbf,nbf),x(ncao,nbf)
      real*8  xcart
      integer lll(20),firstd(nbf),blockl(nbf),idprev
      integer i,j,k,jj,mm,m
      data lll/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/

      real*8 trafo(5,6)
      real*8 trafof(7,10)

      integer lao(ncao),iat,iao
      integer lladr(0:3),ll(0:3)
      data lladr  /1,3,6,10/
      data ll     /0,1,4,10/

      k=0
      do i=1, ncent
         iat=at(i)
         do j=1, bas_nsh(iat)
            do iao=1,lladr(bas_lsh(j,iat))
               k=k+1
               lao(k)=ll(bas_lsh(j,iat))+iao
            enddo
         enddo
      enddo



! Sign of 'trafo(2,:)' changed with respect to xTB
      trafo = 0.0d0
! x2
      trafo(1,1)=1./dsqrt(2.d0)*dsqrt(3.d0/2.d0)
      trafo(2,1)=-0.50d0
! y2
      trafo(1,2)=-1./dsqrt(2.d0)*dsqrt(3.d0/2.d0)
      trafo(2,2)=-0.50d0
! z2
      trafo(1,3)=0.0d0
      trafo(2,3)=1.0d0
c rest
      trafo(3,4)=1.0d0
      trafo(4,5)=1.0d0
      trafo(5,6)=1.0d0

! f-shell CAO(xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz) <- SAO(fz3,fxz2,fyz2,
! fz(x2-y2),fxyz,fx(x2-3y2),fy(3x2-y2)); transpose of the trafo used in dtrf2.f
      trafof = 0.0d0
      trafof(1,3)= 2.0d0
      trafof(1,5)=-3.0d0
      trafof(1,7)=-3.0d0
      trafof(2,1)=-1.0d0
      trafof(2,6)=-1.0d0
      trafof(2,8)= 4.0d0
      trafof(3,2)=-1.0d0
      trafof(3,4)=-1.0d0
      trafof(3,9)= 4.0d0
      trafof(4,5)= 1.0d0
      trafof(4,7)=-1.0d0
      trafof(5,10)=2.0d0
      trafof(6,1)= 1.0d0
      trafof(6,6)=-3.0d0
      trafof(7,2)=-1.0d0
      trafof(7,4)= 3.0d0

      new=ncao-nbf

      if(new.eq.0) then
         x=s
         return
      endif

      firstd = 0
      blockl = 0
      i=1
      j=0
      ! lao is still in old dimensions (i.e., ncao) while s comes with nsao
 42   if(lao(i).gt.4.and.lao(i).le.10)then
         firstd(i-j:i-j+4)=i-j
         blockl(i-j)=2
         j=j+1
         i=i+5
      else if(lao(i).gt.10)then
         firstd(i-j:i-j+6)=i-j
         blockl(i-j)=3
         j=j+3
         i=i+9
      endif
      i=i+1
      if(i.lt.ncao)goto 42
      ! sanity check
      if(new.ne.j) stop 'error in sao2cao trafo'

      x=0.0d0

      do i=1,nbf ! go through eigenvectors
         k = 0
         idprev=0
         do j=1,nbf ! go through LCAO-MO coefficients
            if(idprev.gt.0.and.firstd(j).eq.idprev) cycle
            if(firstd(j).gt.idprev)then ! a new d- or f-shell block starts here
               if(blockl(firstd(j)).eq.3)then ! f-shell: 7 spherical -> 10 cartesian
                  do jj=1,10
                    k=k+1
                    xcart=0.0d0
                    do m=1,7
                       mm=firstd(j)-1+m
                       xcart=xcart+trafof(m,jj)*s(mm,i)
                    enddo
                    x(k,i)=xcart
                  enddo
               else                          ! d-shell: 5 spherical -> 6 cartesian
                  do jj=1,6
                    xcart=0.0d0
                    k=k+1
                    do m=1,5
                       mm=firstd(j)-1+m
                       xcart=xcart+trafo(m,jj)*s(mm,i)
                    enddo
                    x(k,i)=xcart
                  enddo
               endif
                idprev=firstd(j) ! guarantees the rest of this block's spherical fns are skipped
                cycle
            endif
            k=k+1
            x(k,i)= s(j,i)
         enddo
         if (k.ne.ncao) stop 'error in eigenvector dimension'
      enddo

      end
