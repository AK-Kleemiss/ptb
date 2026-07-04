
!! ------------------------------------------------------------------------
!! ------------------------------------------------------------------------
subroutine setupbas0(n,at,ndim)   
      use bascom
      implicit none
      integer n,at(n),ndim

      integer i,iat,ish,iao,pr,lladr(0:3),lladr2(0:3)
      data lladr  /1,3,6,10/
      data lladr2 /1,3,5, 7/

      ncao=0
      nsao=0
      npr =0
      nsh =0
      do i=1,n   
         iat = at(i)
         do ish=1,bas_nsh(iat)
            nsh  = nsh + 1
            nsao = nsao + lladr2(bas_lsh(ish,iat))
            do iao=1,lladr(bas_lsh(ish,iat))
               ncao = ncao + 1
               do pr=1,bas_npr(ish, iat)
                  npr = npr + 1
               enddo
            enddo
         enddo
      enddo

      if(ndim.ne.0.and.nsao.ne.ndim) then
         write(*,*) nsao,ndim
         stop 'ndim <> nsao'
      endif
      ndim = nsao

      allocate(prim_npr  (ncao))
      allocate(prim_count(ncao))
      allocate(prim_exp  (npr) )
      allocate(prim_cnt  (npr) )
      allocate(caoshell(10,n)  )  ! 10 shells max
      allocate( aoshell(10,n)  )  
      allocate( shmap  (10,n)  )  
      allocate(shell2ao(nsao)  )  
      allocate(aoat    (nsao)  )  
end

!! ------------------------------------------------------------------------
!! ------------------------------------------------------------------------
subroutine setupbas(n,at,ndim)
      use bascom
      implicit none
      integer n,at(n),ndim

      integer i,j,k,l,iat,ish,iish,iao,nao,lshell,pr,lladr(0:3),lladr2(0:3)
      data lladr  /1,3,6,10/
      data lladr2 /1,3,5, 7/
      real*8 DEX(-1:96),XNORM,DEX2,PI,EXPO,IAM,alp,tmp

      PI=4.D0*ATAN(1.D0)
      DO I=-1,10
         DEX(I)=DEX2(I)
      ENDDO   

      k  = 0
      j  = 0
      iish=0
      npr= 0
      ncao=0
      nao =0
      do i=1,n   
         iat = at(i)
         do ish=1,bas_nsh(iat)
            iish = iish + 1
            shmap(ish,i)   =iish  ! needed for shell ES
            caoshell(ish,i)=k
             aoshell(ish,i)=j
            lshell=bas_lsh(ish,iat)
            k=k+lladr (lshell)          
            j=j+lladr2(lshell)               
            do iao=1,lladr2(lshell)              
               nao = nao + 1
               aoat(nao) = i
            enddo
            do iao=1,lladr(lshell)              
               ncao = ncao + 1
               prim_npr(ncao)=bas_npr(ish,iat)
               alp = 0
               tmp = 0
               do pr=1,bas_npr(ish,iat)
                  npr= npr+ 1
                  expo = bas_ec(1,pr,ish,iat)
                  iam = lshell
                  XNORM=(2.D0*EXPO/PI)**0.75D0*(4.D0*EXPO)**(IAM/2.D0)/SQRT(DEX(2*IAM-1))
                  prim_exp(npr) = bas_ec(1,pr,ish,iat)
                  prim_cnt(npr) = bas_ec(2,pr,ish,iat)*xnorm
                  if(lshell.eq.2.and.iao.gt.3) prim_cnt(npr)=prim_cnt(npr)*sqrt(3.0d0)
                  if(lshell.eq.3.and.iao.ge.4.and.iao.le.9) prim_cnt(npr)=prim_cnt(npr)*sqrt(5.0d0)
                  if(lshell.eq.3.and.iao.eq.10) prim_cnt(npr)=prim_cnt(npr)*sqrt(15.0d0)
!                 alp=alp+prim_exp(npr)*abs(prim_cnt(npr))
!                 tmp=tmp+              abs(prim_cnt(npr))
               enddo
!              alp=alp/tmp
            enddo
!           bas_aexp(ish,iat) = alp
         enddo
      enddo

      if(iish.ne.nsh) stop 'weird iish'

      nsao=0
      do i=1,n   
         iat = at(i)
         do ish=1,bas_nsh(iat)
            lshell=bas_lsh(ish,iat)
            do iao=1,lladr2(lshell)              
               nsao = nsao + 1
               shell2ao(nsao) = ish
!              write(*,*) 'at,ao,sh',i,nsao,ish
            enddo
         enddo
      enddo

      k=0  
      do i=1,ncao
         prim_count(i)=k
         k=k+prim_npr(i)
      enddo   

      write(*,*) 'ncao ',ncao
      write(*,*) 'nsao ',nsao
      write(*,*) 'npr  ',npr 

end

!! ------------------------------------------------------------------------
!! ------------------------------------------------------------------------
subroutine modbas(n,at,itv)
      use bascom
      use parcom
      implicit none
      integer n,at(n)
      integer itv    

      integer i,j,k,l,iat,ish,iao,lshell,pr,lladr(0:3),lladr2(0:3)
      data lladr  /1,3,6,10/
      data lladr2 /1,3,5, 7/
      real*8 DEX(-1:96),XNORM,DEX2,PI,EXPO,IAM

      PI=4.D0*ATAN(1.D0)
      DO I=-1,10
         DEX(I)=DEX2(I)
      ENDDO   

      k  = 0
      j  = 0
      npr= 0
      do i=1,n   
         iat = at(i)
         do ish=1,bas_nsh(iat)
            lshell=bas_lsh(ish,iat)
            k=k+lladr (lshell)          
            j=j+lladr2(lshell)               
            do iao=1,lladr(lshell)              
               do pr=1,bas_npr(ish,iat)
                  npr= npr+ 1
                  expo = bas_ec(1,pr,ish,iat)*expscal(itv,ish,iat)
                  iam = lshell
                  XNORM=(2.D0*EXPO/PI)**0.75D0*(4.D0*EXPO)**(IAM/2.D0)/SQRT(DEX(2*IAM-1))
                  prim_exp(npr) = bas_ec(1,pr,ish,iat)*expscal(itv,ish,iat)
                  prim_cnt(npr) = bas_ec(2,pr,ish,iat)*xnorm
                  if(lshell.eq.2.and.iao.gt.3) prim_cnt(npr)=prim_cnt(npr)*sqrt(3.0d0)
                  if(lshell.eq.3.and.iao.ge.4.and.iao.le.9) prim_cnt(npr)=prim_cnt(npr)*sqrt(5.0d0)
                  if(lshell.eq.3.and.iao.eq.10) prim_cnt(npr)=prim_cnt(npr)*sqrt(15.0d0)
               enddo
            enddo
         enddo
      enddo
end
!! ------------------------------------------------------------------------
! scale AO exponents individually for each atom
!! ------------------------------------------------------------------------

subroutine modbasd(n,at,scal)
      use bascom
      use parcom
      implicit none
      integer n,at(n)
      integer itv    
      real*8 scal(10,n)

      integer i,j,k,l,iat,ish,iao,lshell,pr,lladr(0:3),lladr2(0:3)
      data lladr  /1,3,6,10/
      data lladr2 /1,3,5, 7/
      real*8 DEX(-1:96),XNORM,DEX2,PI,EXPO,IAM

      PI=4.D0*ATAN(1.D0)
      DO I=-1,10
         DEX(I)=DEX2(I)
      ENDDO   

      k  = 0
      j  = 0
      npr= 0
      do i=1,n   
         iat = at(i)
         do ish=1,bas_nsh(iat)
            lshell=bas_lsh(ish,iat)
            k=k+lladr (lshell)          
            j=j+lladr2(lshell)               
            do iao=1,lladr(lshell)              
               do pr=1,bas_npr(ish,iat)
                  npr= npr+ 1
                  expo = bas_ec(1,pr,ish,iat)*scal(ish,i)
                  iam = lshell
                  XNORM=(2.D0*EXPO/PI)**0.75D0*(4.D0*EXPO)**(IAM/2.D0)/SQRT(DEX(2*IAM-1))
                  prim_exp(npr) = bas_ec(1,pr,ish,iat)*scal(ish,i)
                  prim_cnt(npr) = bas_ec(2,pr,ish,iat)*xnorm
                  if(lshell.eq.2.and.iao.gt.3) prim_cnt(npr)=prim_cnt(npr)*sqrt(3.0d0)
                  if(lshell.eq.3.and.iao.ge.4.and.iao.le.9) prim_cnt(npr)=prim_cnt(npr)*sqrt(5.0d0)
                  if(lshell.eq.3.and.iao.eq.10) prim_cnt(npr)=prim_cnt(npr)*sqrt(15.0d0)
               enddo
            enddo
         enddo
      enddo
end

!! ------------------------------------------------------------------------
!! ------------------------------------------------------------------------

subroutine rdbas(fname)
      use bascom
      use default_files, only: default_basis
      implicit none

      character(len=*), intent(in) :: fname
      character*80 atmp
      integer nn,i,iat,np,l
      real*8 xx(10)
      logical :: ex, use_default, at_end
      integer :: line_idx

      bas_nsh = 0
      bas_lsh = 0
      bas_npr = 0
      bas_ec  = 0

      inquire(file=fname, exist=ex)
      use_default = .not. ex
      if (use_default) then
         print '(a)', "Basis file '"//trim(fname)//"' not found, using built-in defaults."
      else
         open(unit=44,file=fname)
      end if
      line_idx = 0

 10   call get_line(atmp,at_end)
      if (at_end) goto 20
      if(index(atmp,'*').ne.0) then
         call get_line(atmp,at_end)
         if (at_end) goto 20
         read(atmp,*) iat
 12      call get_line(atmp,at_end)
         if (at_end) goto 20
         if(index(atmp,'*').ne.0) goto 10
         call readl(atmp,xx,nn)
         if(nn.eq.1)then
            np = idint(xx(1))
            if(np.gt.5) stop 'contraction > 5 (recompile -> bascom.f90)'
            if(index(atmp,'s').ne.0) l=0
            if(index(atmp,'p').ne.0) l=1
            if(index(atmp,'d').ne.0) l=2
            if(index(atmp,'f').ne.0) l=3
            if(index(atmp,'g').ne.0) stop 'lmax = f'
            bas_nsh(iat)=bas_nsh(iat)+1
            bas_lsh(bas_nsh(iat),iat)=l
            bas_npr(bas_nsh(iat),iat)=np
            do i=1,np
               call get_line(atmp,at_end)
               read(atmp,*) bas_ec(1,i,bas_nsh(iat),iat),bas_ec(2,i,bas_nsh(iat),iat)
            enddo
         endif
         goto 12
      endif
      goto 10
 20   if (.not.use_default) close(44)
      return

contains

      subroutine get_line(line,at_end)
      character(len=*), intent(out) :: line
      logical, intent(out) :: at_end
      integer :: ios
      at_end = .false.
      if (use_default) then
         line_idx = line_idx + 1
         if (line_idx.gt.size(default_basis)) then
            at_end = .true.
         else
            line = default_basis(line_idx)
         end if
      else
         read(44,'(a)',iostat=ios) line
         at_end = (ios.ne.0)
      end if
      end subroutine get_line

end

      DOUBLE PRECISION FUNCTION DEX2(M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF(M .LT. 2) THEN
         DEX2=1
      ELSE
         DEX2=1
         DO 10 I=1,M,2
   10    DEX2=DEX2*I
      ENDIF
      RETURN
      END

