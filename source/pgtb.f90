!! ------------------------------------------------------------------------
! model Hamiltonian to provide atomic shell populations and density matrix
! as well as dipole response
! fitted to shell populations/q/BO/dipole/IR/alpha/Raman/Ekin/matched orbital energies
! SG, 2020-2022
!
! last change on method or HCNO/global parameters: Sat Feb  5 09:38:39 CET 2022
!
!! ------------------------------------------------------------------------
! CPU  profiling yields: 50 % for 2 Eigensolves, 
!                        25 % for ML-matrix Eigensolve and setup
!                        25 % for 5 gemm (Dmat, Vpauli, pops)
!                        non-linalg is < 5 %
!! ------------------------------------------------------------------------
!        read(1,*) shell_xi  (1:10,j)    ! 71-80  
!        read(1,*) shell_cnf1(1:10,j)    ! 81-90    
!        read(1,*) shell_cnf2(1:10,j)    ! 91-100   
!        read(1,*) shell_cnf3(1:10,j)    ! 101-110  
!        read(1,*) expscal  (3,1:10,j)   ! 111-120  
!        read(1,*) shell_cnf4(1:10,j)    ! 121-130  
!        read(1,*) shell_resp(1:10,j,1)  ! 131-140  
!        read(1,*) empty                 ! 141-150


subroutine pgtb(pr,prop,n,ndim,nel,nopen,homo,at,chrg,xyz,z,rab, & 
&               pnt,norm,S,T,D,efield,S1,S2,psh,pa,&
&               P,H,eps,eel,ecoul,wbo,dip,alp)
   use iso_fortran_env, only : wp => real64
   use parcom
   use bascom
   use com
   use mocom ! fit only
   implicit none

!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
   logical, intent(in)    :: pr                 ! print
   integer, intent(in)    :: prop               ! type of property calc (0: p,q,P,WBO  1/-1: +dip  2/-2: + alpha  
                                                !                        3: beta       4   : stda write (if < 0, no WBO)
   integer, intent(in)    :: n                  ! number of atoms 
   integer, intent(in)    :: ndim               ! number of AOs       
   integer, intent(in)    :: nel                ! number of electrons 
   integer, intent(in)    :: nopen              ! number of open shells
   integer, intent(in)    :: homo               ! as the name says...
   integer, intent(in)    :: at(n)              ! ordinal number of atoms
   real(wp),intent(in)    :: chrg               ! system charge           
   real(wp),intent(in)    :: xyz(3,n)           ! coordinates (not used)
   real(wp),intent(in)    :: z(n)               ! nuclear charges          
   real(wp),intent(in)    :: rab(n*(n+1)/2)     ! distances  
   real(wp),intent(in)    :: pnt(3)             ! property reference point
   real(wp),intent(in)    :: norm(ndim)         ! SAO normalization factors
   real(wp),intent(in)    :: S(ndim*(ndim+1)/2) ! exact overlap maxtrix in SAO
   real(wp),intent(in)    :: T(ndim*(ndim+1)/2) ! exact T maxtrix in SAO
   real(wp),intent(in)    :: D(ndim*(ndim+1)/2,3)!dipole integrals
   real(wp),intent(in)    :: efield(3)          ! electric field
!! ------------------------------------------------------------------------
!  Output
!! ------------------------------------------------------------------------
   real*4  ,intent(out)   :: S1(ndim,ndim)      ! Mull-Loew trafo   
   real*4  ,intent(out)   :: S2(ndim,ndim)      !   "   "    "
   real(wp),intent(out)   :: psh(10,n)          ! shell populations 
   real(wp),intent(inout) :: pa(n)              ! atom      "
   real(wp),intent(inout) :: P(ndim*(ndim+1)/2) ! density matrix
   real(wp),intent(inout) :: H(ndim*(ndim+1)/2) ! unperturbed Hamilton Matrix
   real(wp),intent(out)   :: eps(ndim)          ! eigenvalues
   real(wp),intent(out)   :: eel                ! electronic energy
   real(wp),intent(out)   :: ecoul              ! Coulomb energy
   real(wp),intent(out)   :: wbo(n,n)           ! WBOs                 
   real(wp),intent(out)   :: dip(3)             ! dipole moment
   real(wp),intent(out)   :: alp(6)             ! dipole polarizability tensor

!! ------------------------------------------------------------------------
!  local
!! ------------------------------------------------------------------------
   real(wp),allocatable :: focc(:)              ! occupations
   real(wp),allocatable :: Hdiag(:)             ! diagonal
   real(wp),allocatable :: SS(:)                ! scaled overlap/perturbed H
   real(wp),allocatable :: Vecp(:)              ! ECP ints
   real(wp),allocatable :: Htmp(:)              ! modified H                    
   real(wp),allocatable :: P1  (:)              ! perturbed P     
   real(wp)             :: alpha(3,3)         
   real(wp)             :: dip1(3),dip2(3)

   real(wp),parameter   :: eT     = 300.0_wp   ! electronic temp. 
   real(wp),parameter   :: ffs    = 0.001_wp   ! finite field step
   real(wp),parameter   :: erfs   =-2.000_wp   ! erf expo for CN (about 1/5 orig. value)

   integer  :: i,j,k,ish,jsh,ati,ia,ii,lin
   integer  :: ns, iter, ibeta
   integer  :: llao(4)
   data llao /1,3,5,7 /
   real(wp) :: r,r2,rcovij,arg,tmp,hi,hj,hij,pol,t8,t9,xk,w0,w1,t0,t1
   real(wp) :: scfpar(8)
   real(wp) :: scal(10,n)
   real(wp),allocatable :: cn(:), cns(:), cnorg(:), xnrm(:), qeeq(:)
   real(wp),allocatable :: ves0(:), ves1(:), gab(:,:)
   real(wp),allocatable :: patmp(:), pshtmp(:,:)
   real(wp),parameter   :: au2ev = 27.2113957_wp
   logical fail, highsym

!! ------------------------------------------------------------------------
!  initizialization
!! ------------------------------------------------------------------------
   if(pr)then
   write(*,*) '----------------'
   write(*,*) '| P-gTB  model |'
   write(*,*) '| SG 2020-2021 |'
   write(*,*) '----------------'
   write(*,*) 'prop :',prop
   endif

!  call prmat(6,S,ndim,0,'S')

   totmatch=0
   scfpar = 0
   ibeta  = 1
   if(prop.eq.3.or.prop.eq.102) ibeta = 2  ! different correction (parameters) for beta

   allocate( SS(ndim*(ndim+1)/2), focc(ndim), Hdiag(ndim), &
  &          ves0(nsh), ves1(nsh), gab(nsh,nsh), &
  &          Htmp(ndim*(ndim+1)/2),xnrm(ndim), &
  &          Vecp(ndim*(ndim+1)/2), &
  &          qeeq(n), cn(n), cns(n), cnorg(n), source = 0.0_wp )

   if(prop.gt.100) goto 999  ! just property evaluation for given H,P,pa...
                             ! this saves twoscf step in beta evaluation

!! ------------------------------------------------------------------------
!  calculation of P,q,... starts here
!! ------------------------------------------------------------------------

   if(pr)write(*,'('' PA setup with Mull-Loew partitioning '',f10.6,'' ...'')') mull_loew14
!  the ML can be done in real*4 without loss in accuracy
   call mlpop14(ndim,S,S1,S2) ! ML precalc, x=1/4
   if(pr)write(*,*) 'done.'

!  CNs
   call ncoord_erf(n,at,rab,erfs,cns) ! the org. EEQ was parameterized with erfs=-7.5, same expo. here for both new CN
                                      ! for this erfs, the avcn values must be generated, see avcn.f and setavcn!
   call ncoord_erf(n,at,rab,-7.5_wp,cnorg) 
!  special CN with fitted radii
   cn = 0_wp
   do i = 2, n
      do j = 1, i-1 
         r = rab(lin(i,j))
         rcovij=abs(shell_cnf4(3,at(i)))+abs(shell_cnf4(3,at(j))) ! to avoid num. problems in fit if radius -> 0
         arg = (r-rcovij)/rcovij
         tmp = 0.5_wp * (1_wp + erf(erfs*arg)) 
         cn(i) = cn(i) + tmp      
         cn(j) = cn(j) + tmp 
      enddo
   enddo

   call eeq(n,at,rab,chrg,cnorg,.false.,shell_cnf1(9,:),shell_cnf1(8,:),qeeq)    !  slightly modified EEQ charges qeeq for first iter Ves
!  call eeq(n,at,rab,chrg,cnorg,.true. ,shell_cnf1(9,:),shell_cnf1(8,:),qeeq_org)!  orig EEQ

   if(pr)then
   write(*,'(''EEQ done. sum q : '',f8.3,i5)') sum(qeeq)
   write(*,'(''    atom   Zeff  qEEQ  mod   srCN   noHCN   CN(EEQ)'')')
   do i=1,n
       write(*,'(2i5,f5.1,2x,6f8.3)') i,at(i),z(i),qeeq(i),cns(i),cn(i),cnorg(i)
   enddo
   endif

!  special overlap matrix for H0
   call modbas(n,at,3) 
   call sint(n,ndim,at,xyz,rab,SS,eps) ! scaled S 
   call modbas(n,at,4) 

   ! atomic H0 
   jsh = 0
   ii  = 0
   do i = 1, n
      ati = at(i)
      t9  =         cns(i)*shell_cnf4( 9,ati) ! shift
      t8  = cn(i) + cns(i)*shell_cnf4(10,ati) ! shell-wise
      do ish=1,bas_nsh(ati)
          jsh = jsh + 1
          tmp = shell_xi(ish,ati) + shell_cnf1(ish,ati) * t8 + t9 
          do j=1,llao(bas_lsh(ish,ati)+1) ! AO loop
             ii = ii + 1
             Hdiag(ii) = tmp
          enddo
      enddo
   enddo

!  simple ECP
   if(pr)write(*,*) 'computing Vecp ...'
   call calcvecp (n,ndim,at,xyz,rab,norm,Vecp)

!! ------------------------------------------------------------------------
!  set up the H matrix twice
!! ------------------------------------------------------------------------

!    0.3876991369   0.9416938677   0.4680146099   0.0421282612   0.0000000000   1.0883439565  -0.3320393628
   pa = qeeq ! initialization of atomic charges, true pa on output
   scfpar(1) =  glob_par(11)  ! gpol 
   scfpar(2) =  glob_par(14)  ! Wolfsberg l dep.
   scfpar(4) =  glob_par(16)  ! iter1 off-diag
   scfpar(5) =  0.020_wp      ! glob_par(14)  ! gamscal 24
   scfpar(6) =  glob_par(13)  ! iter1 two-center
   scfpar(7) =  0.00_wp       ! gamscal in onescf
   scfpar(8) =  glob_par(12)  ! gpol in onescf
   call twoscf(pr,prop,n,ndim,nel,nopen,homo,at,xyz,z,rab,cns,S,T,SS,Vecp,Hdiag,focc,&
               norm,eT,scfpar,S1,S2,psh,pa,P,H,ves0,gab,eps,eel,ecoul)  

!! ------------------------------------------------------------------------
!  done
!! ------------------------------------------------------------------------

!! ------------------------------------------------------------------------
!  properties   
!! ------------------------------------------------------------------------

999 continue  ! entry point for beta calcs when only a field is added

   Htmp = H 

! MO match for fit
   if(cmo_ref(1,1).gt.-98.999d0.and.prop.ge.0.and.prop.lt.4) then
      call getsymmetry(pr,n,at,xyz,0.01d0,50,highsym) ! get PG to check for MO degen.
      call momatch(pr,highsym,ndim,homo,0,S)          ! which modifies match routine
   endif

! stda
   if(prop.eq.4) call printmos(n,at,xyz,ndim,homo,norm,2d0) ! cut virt. > 2 Eh because very high lying gTB MOs are crap
   if(prop.eq.5) call wr_tm_mos(ndim,n,nel,at,nopen,ndim)   ! write for TM all mos (incl. virts!)
   close(42,status='delete')

! just with field
   if(sum(abs(efield)).gt.1.d-6) then
      do j=1,3
         Htmp(:)=Htmp(:)-efield(j)*D(:,j)  ! the field perturbation on unperturbed H
      enddo
      if(prop.ne.102) then  ! if this is not a beta calc, p,q,mu... are computed with field P
         call solve2(2,n,ndim,nel,nopen,homo,at,eT,focc,Htmp,S,P,eps,eel,fail) ! solve with efield
      endif
   endif

! WBO
  if(prop.ge.0) call wiberg(n,ndim,at,rab,P,S,wbo)

! dipole moment 
   if(abs(prop).gt.0)then
     call dipmom2(n,ndim,xyz,z,norm,P,D,pnt,dip) ! get dipole moment at point pnt
   endif

! polarizability by simple perturbative treatment 
! this is only done in alpha,beta cases
   if(abs(prop).eq.2.or.prop.eq.102)then
! six perturbed dipole moment calcs H = H_final + H_resp + field1 + field2
     allocate(P1(ndim*(ndim+1)/2),patmp(n),pshtmp(10,n))
     do k=1,3
        call addsym(ndim, ffs,Htmp,D(1,k),H)                                           ! perturb field free H with field 
        call solve2(1,n,ndim,nel,nopen,homo,at,eT,focc,H,S,P1,eps,tmp,fail)            ! solve 
        call mlpop2(n,ndim,P1,S1,S2,patmp,pshtmp)                                      ! pops
        patmp = z - patmp
        H = Vecp 
        call adddsym(ndim, ffs,D(1,k),H)                                               ! perturb H with field only
        call onescf(n,ndim,nel,nopen,homo,at,rab,cns,&                                 ! and add 2nd iter part
     &              S,SS,H,Hdiag,focc,eT,scfpar,ves0,pshtmp,patmp,P1)  
        call dipmom2(n,ndim,xyz,z,norm,P1,D,pnt,dip1)                                  ! get dipole moment

        call addsym(ndim,-ffs,Htmp,D(1,k),H)                                           ! other direction
        call solve2(1,n,ndim,nel,nopen,homo,at,eT,focc,H,S,P1,eps,tmp,fail)   
        call mlpop2(n,ndim,P1,S1,S2,patmp,pshtmp)                         
        patmp = z - patmp
        H = Vecp 
        call adddsym(ndim,-ffs,D(1,k),H)                                
        call onescf(n,ndim,nel,nopen,homo,at,rab,cns,&
     &              S,SS,H,Hdiag,focc,eT,scfpar,ves0,pshtmp,patmp,P1)  
        call dipmom2(n,ndim,xyz,z,norm,P1,D,pnt,dip2)                    

        alpha(k,1:3)=-(dip1(1:3)-dip2(1:3))/(2_wp*ffs)                               ! numerical diff. dmu/dfield
     enddo
     alp(1)=alpha(1,1)
     alp(2)=0.5*(alpha(2,1)+alpha(1,2))
     alp(3)=alpha(2,2)
     alp(4)=0.5*(alpha(3,1)+alpha(1,3))
     alp(5)=0.5*(alpha(3,2)+alpha(2,3))
     alp(6)=alpha(3,3)
   endif

end

!! ------------------------------------------------------------------------
!  set up and diag the H matrix twice
!! ------------------------------------------------------------------------

subroutine twoscf(pr,prop,n,ndim,nel,nopen,homo,at,xyz,z,rab,cn,S,T,SS,Vecp,Hdiag,focc,&
                  norm,eT,scfpar,S1,S2,psh,pa,P,Hmat,ves,gab,eps,eel,ecoul)
   use iso_fortran_env, only : wp => real64
   use bascom
   use parcom
   use com
   implicit none 
!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
   logical, intent(in)    :: pr                    ! print
   integer, intent(in)    :: prop                  ! calc type
   integer, intent(in)    :: n                     ! number of atoms 
   integer, intent(in)    :: ndim                  ! number of AOs       
   integer, intent(in)    :: nel                   ! number of electrons 
   integer, intent(in)    :: nopen                 ! number of open shells
   integer, intent(in)    :: homo                  ! as the name says...
   integer, intent(in)    :: at(n)                 ! ordinal number of atoms
   real(wp),intent(in)    :: xyz(3,n)              ! coordinates
   real(wp),intent(in)    :: z(n)                  ! nuclear charges          
   real(wp),intent(in)    :: rab(n*(n+1)/2)        ! distances  
   real(wp),intent(in)    :: cn(n)                 ! CN           
   real(wp),intent(in)    :: S(ndim*(ndim+1)/2)    ! exact overlap maxtrix in SAO
   real(wp),intent(in)    :: T(ndim*(ndim+1)/2)    ! exact T maxtrix in SAO
   real(wp),intent(in)    :: SS(ndim*(ndim+1)/2)   ! scaled overlap maxtrix in SAO
   real(wp),intent(in)    :: Vecp(ndim*(ndim+1)/2) ! ECP ints
   real(wp),intent(in)    :: Hdiag(ndim)           ! diagonal of H0
   real(wp),intent(in)    :: focc (ndim)           ! fractional occ.
   real(wp),intent(in)    :: norm (ndim)           ! SAO normalization factors
   real(wp),intent(in)    :: eT                    ! el. temp.
   real(wp),intent(in)    :: scfpar(8)             ! parameters
   real*4  ,intent(in)    :: S1(ndim,ndim)         ! ML trafo
   real*4  ,intent(in)    :: S2(ndim,ndim)         ! "   "          

!! ------------------------------------------------------------------------
!  Output
!! ------------------------------------------------------------------------
   real(wp),intent(inout)   :: psh(10,n)             ! shell populations 
   real(wp),intent(inout)   :: pa(n)                 ! atom      "
   real(wp),intent(inout)   :: P   (ndim*(ndim+1)/2) ! density matrix
   real(wp),intent(inout)   :: Hmat(ndim*(ndim+1)/2) ! Hamiltonian matrix
   real(wp),intent(inout)   :: ves(nsh)              ! 
   real(wp),intent(inout)   :: gab(nsh,nsh)          ! 
   real(wp),intent(out)     :: eps (ndim)            ! orbital energies
   real(wp),intent(out)     :: eel                   ! electronic energy = sum eps
   real(wp),intent(out)     :: ecoul                 ! Coulomb energy 

!  local
   logical  :: fail
   integer  :: i,j,k,l,ish,ati,atj,ia,ib,jsh,ii,jj,lin,ij,li,lj,iter,iish,jjsh,mode
!  real(wp),parameter :: cok   = 0.95_wp  ! OK mixing in 1. iter, not as important as in second iter
!  real(wp),parameter :: cmn   = 1_wp-cok     
   real(wp),parameter :: au2ev = 27.2113957_wp
   real(wp) :: r,tmp,pol,hi,hj,hij,xk,t8,t9,qa,qb,keav,eh1,dmp,tmp2
   real(wp) :: xiter(2),ziter(2),ssh
   real(wp) :: t0,t1,w0,w1
   real(wp) :: gq(n), geff(n), xab(n*(n+1)/2)


   ziter(1)=scfpar(4)
   ziter(2)=1_wp
   xiter(1)=scfpar(6)
   xiter(2)=1_wp

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do iter=1, 2         ! two "iterations": in the first, q (=pa) = q(EEQ) and NO P (=+U)

   if(iter.eq.1) then   ! atom-wise ES, no shell pop available
      do i=1, n
         geff(i) = 1_wp / gam(at(i))
      enddo
      do i=1, n
         eh1 = 0_wp
         do l=1,n
            k=lin(l,i)
            xk =2_wp /(geff(i) + geff(l))                   ! harm. av.
!           tmp= cok/sqrt(rab(k)**2+1_wp/xk**2) + cmn/(rab(k)+1_wp/xk)  
            tmp= 1d0/sqrt(rab(k)**2+1_wp/xk**2)
            eh1=eh1+pa(l)*tmp                               ! contraction with charge (DFTB2 atomic term)
         enddo
         gq(i) = eh1 * 0.5_wp 
      enddo
   else                  ! shell-wise xTB like
      call setgab  (n,at,rab,pa,scfpar(5),gab)  ! the gab contain q as higher order effect on Ves
      call setespot(n,at,psh,gab,ves) 
      ves = ves * 0.5_wp
   endif    

! H0 +  third-order (atomic charge exists in 1. AND 2. iter)
   do i=1, n
      geff(i) = pa(i)**2*shell_cnf4(1,at(i)) ! geff is temp.
   enddo

   Hmat = 0_wp
   ij = 0
   do i=1,ndim
      ia = aoat(i)
      ati= at(ia)
      ish= shell2ao(i)
      li = bas_lsh(ish,ati)
      hi = Hdiag(i)
      do j=1,i  
         ij = ij + 1
         ib = aoat(j)
         r  = rab(lin(ia,ib))
         if(r.gt.50d0) cycle
         hj = Hdiag(j)
         hij= hi+hj
         ssh= hij * SS(ij)
         atj= at(ib)
         jsh= shell2ao(j)
         if(ia.ne.ib) then            ! different atoms
            lj  = bas_lsh(jsh,atj)
            xk  = (shell_cnf4(2,ati)+shell_cnf4(2,atj)) * xiter(iter)
            pol = ((hi-hj)/hij)**2
            keav= 0.5_wp*(shell_cnf2(9,ati)+dble(li)*scfpar(2) + shell_cnf2(9,atj)+dble(lj)*scfpar(2))
            tmp = ssh * keav * (1_wp-pol*scfpar(1)) * (1_wp+xk/r) ! fit yields same values for iter1,2
         else                         ! same atoms
            if(ish.ne.jsh) then       ! s-s', p-p', d-d' off-diagonal, li=lj because S=0 otherwise
               tmp2= shell_cnf4(4+li,ati) * ziter(iter)
               tmp = ssh * tmp2 + shell_cnf3(9,ati)* tmp2 * hij * SS(ij)**2 ! second term only for more than 2 shells of same l
            else
               tmp = ssh
            endif
         endif
!                                               third order diagonal
         Hmat(ij) = Hmat(ij) + tmp + Vecp(ij) - S(ij)*(geff(ia)+geff(ib))
      enddo
   enddo

! H1
   if(iter.eq.1)then ! atom-wise, no +U (no P exists)
    call calcpauli1(n,ndim,at,z,pa      ,S,Hdiag,Hmat) ! add valence X correction based on three-index ECP formula using qeeq
    ij = 0
    do i=1,ndim
      ia = aoat(i)
      do j=1,i  
         ib = aoat(j)
         ij = ij + 1
         Hmat(ij) = Hmat(ij) - (gq(ia)+gq(ib)) * S(ij)
      enddo
    enddo
   else ! shell-wise ES + U terms with density from 1. iter
!   for +U LR damping
    k = 0
    do i=1,n
      gq(i) = 1_wp-(shell_xi(9,at(i))*pa(i)+shell_xi(10,at(i))*pa(i)**2) ! gq is temp., important charge scaling
      hi = shell_cnf3(10,at(i)) + (cn(i)-avcn(at(i)))*shell_resp(10,at(i),1)
      do j=1,i
         k = k + 1
         r = hi + shell_cnf3(10,at(j)) + (cn(j)-avcn(at(j)))*shell_resp(10,at(j),1)
         t8= (rab(k)-r)/r
         xab(k) = 0.5_wp*(1_wp+erf(-1.8_wp*t8)) ! paramter not very important
      enddo
    enddo
    call calcpauli2(n,ndim,at,psh,S,Hdiag,Hmat) ! add valence X correction based on three-index ECP formula, now using psh
    k = 0
    do i=1,ndim
      ia = aoat(i)
      ati= at(ia)
      ish= shell2ao(i)
      iish=shmap(ish,ia)
      hi = shell_cnf2(ish,ati)*gq(ia)
      do j=1,i
         k  = k + 1
         ib = aoat(j)
         atj= at(ib)
         jsh= shell2ao(j)
         jjsh=shmap(jsh,ib)
         hj = shell_cnf2(jsh,atj)*gq(ib)
!                            this part is INDO two-c like         shell ES
         Hmat(k) = Hmat(k) + P(k) * (hi + hj) * xab(lin(ib,ia)) - S(k)*(ves(iish)+ves(jjsh)) 
      enddo
    enddo
   endif

   if(pr)write(*,*) 'gTB H matrix iteration ',iter, ' done. Now diag ...'

!  solve 
   mode = iter
   if(iter.eq.2.and.prop.eq.4) mode = 3     ! stda write
   if(iter.eq.2.and.prop.eq.5) mode = 4     ! TM write
   if(              prop.lt.0) mode = -iter ! IR/Raman  
   call solve2 (mode,n,ndim,nel,nopen,homo,at,eT,focc,Hmat,S,P,eps,eel,fail) 
   if(fail) stop 'diag error'

   if(iter.eq.2.and.pr)then
     ii=max(homo-9,1)
     jj=min(homo+2,ndim)
     xk=eps(homo+1)-eps(homo)
     write(*,'('' frontier MO occupations   : '',14f8.4)') focc(ii:jj)
     write(*,'('' (shifted) level energies  : '',14f8.4)')  eps(ii:jj)
     write(*,'('' gap (eV)                  : '',f9.3)') xk*au2ev    
     write(*,'('' virtual MO correction (eV): '',f9.3)') (glob_par(4) + glob_par(5)*xk)*au2ev     
     if(xk.lt.0.05) then
        write(*,*) 'WARNING WARNING WARNING WARNING'
        write(*,*) ':::::::   small HL gap  :::::::'
        write(*,*) 'WARNING WARNING WARNING WARNING'
     endif
   endif

!  pop
   call mlpop2(n, ndim, P, S1, S2, pa, psh)
   if(iter.eq.1) pa = z - pa  ! output are populations but here q is used for convenience

   enddo
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  call calces(n,at,scfpar(5),rab,z,psh,gab,ecoul)

end

!! ------------------------------------------------------------------------
!  add Pauli Term to H in first iter
!! ------------------------------------------------------------------------

subroutine calcpauli1(n,nao,at,z,q,S,Hdiag,Hmat)
      use  bascom
      use  parcom
      use gtb_la, only : la_symm
      implicit none          
      integer, intent(in)   :: nao,n,at(n)
      real*8,  intent(in)   :: z(n),q(n)
      real*8,  intent(in)   :: S(nao*(nao+1)/2)    
      real*8,  intent(in)   :: Hdiag(nao)    
      real*8,  intent(inout):: Hmat(nao*(nao+1)/2)    

      integer i,j,k,l,nl,m,atk,jsh,llao2(0:3)
      data llao2/1,3,5,7 /
      real*8 atocc(10), f1,f2
      real*4,allocatable :: stmp(:,:), sdum(:,:), xtmp(:,:)

      allocate(stmp(nao,nao),sdum(nao,nao),xtmp(nao,nao))

      call blowsym84(nao,S,sdum)

!     N^2 step
      do i=1,nao                          
         m=0 
         do k=1,n       
            atk=at(k)
            f1=(z(k)-q(k))/z(k)             ! Nel factor
            f1=f1*shell_cnf2(10,atk)        ! element scaling
            call shellocc_ref(atk,atocc)    ! ref. atomic pop.
            do jsh=1,bas_nsh(atk)           ! shells of atom nn
               nl=llao2(bas_lsh(jsh,atk))
               f2=f1/dble(nl)               ! AO degen
               do l=1,nl                    ! AOs of shell jsh
                  m = m + 1
                  stmp(m,i)= Hdiag(m) * sdum(m,i) * atocc(jsh) * f2
               enddo
            enddo
         enddo
      enddo

!     N^3 step
!     call la_gemm('N','N',nao,nao,nao,1.0e0,sdum,nao,stmp,nao,0.0e0,xtmp,nao)
      call la_symm('L','L',nao,nao,1e0,sdum,nao,stmp,nao,0e0,xtmp,nao)   

      k = 0 
      do i=1, nao 
         do j=1, i
            k = k + 1 
            Hmat(k) = Hmat(k) + xtmp(j,i)                         
         enddo
      enddo

      end

!! ------------------------------------------------------------------------
!  add Pauli Term to H in second iter
!! ------------------------------------------------------------------------

subroutine calcpauli2(n,nao,at,psh,S,Hdiag,Hmat)
      use  bascom
      use  parcom
      use gtb_la, only : la_gemm, la_symm
      implicit none          
      integer, intent(in)   :: nao,n,at(n)
      real*8,  intent(in)   :: psh(10,n)
      real*8,  intent(in)   :: S(nao*(nao+1)/2)    
      real*8,  intent(in)   :: Hdiag(nao)    
      real*8,  intent(inout):: Hmat(nao*(nao+1)/2)    

      integer i,j,k,l,m,nl,atk,jsh,llao2(0:3)
      data llao2/1,3,5,7 /
      real*8 f1
      real*4,allocatable :: stmp(:,:), sdum(:,:), xtmp(:,:)

      allocate(stmp(nao,nao),sdum(nao,nao),xtmp(nao,nao))

      call blowsym84(nao,S,sdum)

!     N^2 step
      do i=1,nao                          
         m=0 
         do k=1,n       
            atk=at(k)
            do jsh=1,bas_nsh(atk)                ! shells of atom nn
               l =bas_lsh(jsh,atk)
               nl=llao2(l)
!              f1=shell_resp(l+1,atk,1)/dble(nl) ! element,l scaling
               f1=shell_resp(jsh,atk,1)/dble(nl) ! element shell wise scaling
               do l=1,nl                         ! AOs of shell jsh
                  m = m + 1
                  stmp(m,i)= Hdiag(m) * sdum(m,i) * psh(jsh,k) * f1
               enddo
            enddo
         enddo
      enddo

!     N^3 step
!     call la_gemm('N','N',nao,nao,nao,1.0e0,sdum,nao,stmp,nao,0.0e0,xtmp,nao)
      call la_symm('L','L',nao,nao,1.0e0,sdum,nao,stmp,nao,0.0e0,xtmp,nao)   

      k = 0 
      do i=1, nao 
         do j=1, i
            k = k + 1 
            Hmat(k) = Hmat(k) + xtmp(j,i)                         
         enddo
      enddo

      end

!! ------------------------------------------------------------------------
!  set up Coulomb potential due to 2nd order fluctuation shell-wise
!! ------------------------------------------------------------------------

subroutine setgab(n,at,rab,q,gamscal,gab)
   use bascom
   use parcom
   use com
      implicit none 
      integer, intent(in)  :: n
      integer, intent(in)  :: at(n)
      real*8,  intent(in)  :: rab(n*(n+1)/2)  
      real*8,  intent(in)  :: q(n)
      real*8,  intent(in)  :: gamscal
      real*8,  intent(out) :: gab(nsh,nsh)

      integer i,j,k,ati,atj,ish,jsh,ii,jj,lin
      real*8 gish,gjsh,xk,r2,geff(n),ff
      real*8,parameter :: cok   = 0.9d0
      real*8,parameter :: cmn   = 1d0-cok

      do i=1,n
         geff(i) = (1d0 - gamscal*q(i))*gam(at(i))
      enddo

!     DFTB second order term J matrix
      ii = 0
      do i=1, n
      ati = at(i)
      do ish=1, bas_nsh(ati)
         ii = ii + 1
         gish = shell_cnf3(ish,ati) * geff(i) ! important higher-order effect
         jj = 0
         do j=1,n
            k = lin(j,i)
            r2= rab(k)**2
            atj = at(j)
            do jsh=1, bas_nsh(atj)
               jj = jj + 1
               if (jj.gt.ii) cycle
               gjsh = shell_cnf3(jsh,atj) * geff(j)
               xk   = 2d0 /(1d0/gish + 1d0/gjsh)    ! harm. av.
               gab(jj,ii)= cok/sqrt(r2+1d0/xk**2) + cmn/(rab(k)+1d0/xk) !Ohno-Klopman-Mataga average
               gab(ii,jj)= gab(jj,ii)
            enddo
         enddo
      enddo
      enddo

end

subroutine setespot(n,at,qsh,gab,ves)
   use iso_fortran_env, only : wp => real64
   use bascom
      implicit none 
      integer, intent(in)  :: n
      integer, intent(in)  :: at(n)
      real*8,  intent(in)  :: qsh(10,n)
      real*8,  intent(in)  :: gab(nsh,nsh)
      real*8,  intent(out) :: ves(nsh) 

      integer i,j,ati,ish,iish
      real*8  vesi,qshi
      real*8  atocc(10)
      real*8  qshtmp(nsh)

      iish = 0
      do i=1,n
         ati = at(i)
         do ish=1,bas_nsh(ati)
            iish = iish + 1
            call shellocc_ref(ati,atocc) ! ref. atomic pop.
            qshtmp(iish)=atocc(ish)-qsh(ish,i)
         enddo
      enddo

!     taken from xtb GFN1 part
      ves = 0.0_wp
      do i=1,nsh
         qshi=qshtmp(i)
         vesi=0.0_wp
         do j=1,i-1
            ves(j)=ves(j)+qshi*gab(j,i) 
            vesi=vesi+qshtmp(j)*gab(j,i)
         enddo
         vesi=vesi+qshi*gab(i,i)
         ves(i)=ves(i)+vesi
      enddo

end  

subroutine calces(n,at,scfpar,rab,z,qsh,gab,ecoul)
   use iso_fortran_env, only : wp => real64
   use bascom
      implicit none 
      integer, intent(in)  :: n
      integer, intent(in)  :: at(n)
      real*8,  intent(in)  :: scfpar          
      real*8,  intent(in)  :: rab(n*(n+1)/2)  
      real*8,  intent(in)  :: z(n)
      real*8,  intent(in)  :: qsh(10,n)
      real*8,  intent(in)  :: gab(nsh,nsh)
      real*8,  intent(out) :: ecoul    

      integer i,j,ati,ish,iish
      real*8  atocc(10)
      real*8  qshtmp(nsh)
      real*8  q(n)

      iish = 0
      do i=1,n
         ati = at(i)
         q(i)= 0
         do ish=1,bas_nsh(ati)
            iish = iish + 1
            call shellocc_ref(ati,atocc) ! ref. atomic pop.
            qshtmp(iish)=atocc(ish)-qsh(ish,i)
            q(i)=q(i)+qsh(ish,i)
         enddo
         q(i)=z(i)-q(i)
      enddo

      call setgab  (n,at,rab,q,scfpar,gab)  ! the gab contain q as higher order effect on Ves

      ecoul=0
      do i=1,iish-1
         do j=i+1,iish
            ecoul =ecoul + qshtmp(i)*qshtmp(j)*gab(j,i)
         enddo
      enddo
      ecoul=ecoul*2.0d0
      do i=1,iish
         ecoul =ecoul + qshtmp(i)*qshtmp(i)*gab(i,i)
      enddo
      ecoul=ecoul*0.5d0

end  

!! ------------------------------------------------------------------------
!  set the average (common) CN of elements
!  use code avcn.f
!! ------------------------------------------------------------------------

subroutine setavcn ! same but with erfs=-2
      use com
      avcn = 4d0
      avcn( 1)= 0.8571
      avcn( 6)= 3.1576
      avcn( 7)= 2.6221
      avcn( 8)= 1.5218
      avcn( 2)= 0.5576
      avcn(10)= 0.6325
      avcn(18)= 0.9971
      avcn(36)= 1.2780
      avcn(54)= 1.4142
      avcn(86)= 1.4813
      avcn( 9)= 1.1580
      avcn(17)= 1.6193
      avcn(35)= 1.7517
      avcn(53)= 1.9027
      avcn(85)= 1.6944
      avcn(16)= 2.2996
      avcn(34)= 2.2906
      avcn(52)= 2.3380
      avcn(84)= 2.3199
      avcn(15)= 3.4848
      avcn(33)= 3.1793
      avcn(51)= 3.1065
      avcn(83)= 3.1497
      avcn(14)= 3.1467
      avcn(32)= 3.1015
      avcn(50)= 3.3329
      avcn(82)= 3.0994
      avcn( 5)= 2.7596
      avcn(13)= 3.4564
      avcn(31)= 3.3721
      avcn(49)= 3.4365
      avcn(81)= 3.2478
      avcn( 4)= 2.2666
      avcn(12)= 2.8134
      avcn(20)= 4.0839
      avcn(38)= 2.8011
      avcn(56)= 2.5869
      avcn( 3)= 2.5636
      avcn(11)= 2.4115
      avcn(19)= 2.7819
      avcn(37)= 3.0927
      avcn(55)= 3.2766
      avcn(30)= 2.4055
      avcn(29)= 2.3970
      avcn(28)= 3.3283
      avcn(27)= 4.6666
      avcn(26)= 5.4032
      avcn(25)= 5.1648
      avcn(24)= 4.8825
      avcn(23)= 5.2103
      avcn(22)= 4.1372
      avcn(21)= 4.6882
      avcn(48)= 2.5331
      avcn(47)= 2.5553
      avcn(46)= 3.3691
      avcn(45)= 5.1645
      avcn(44)= 5.5080
      avcn(43)= 5.1423
      avcn(42)= 5.7108
      avcn(41)= 5.2953
      avcn(40)= 4.4574
      avcn(39)= 5.3701
      avcn(80)= 2.5980
      avcn(79)= 2.4969
      avcn(78)= 3.5058
      avcn(77)= 4.9224
      avcn(76)= 5.5049
      avcn(75)= 5.1728
      avcn(74)= 5.4683
      avcn(73)= 5.1258
      avcn(72)= 4.3281
      avcn(57)= 5.7731

end

!! ------------------------------------------------------------------------
!  normalize shell occ. to q_atom = 0, just to make shell_occ_ML 
!  exclude semi-core for group 1,2 and TMs
!! ------------------------------------------------------------------------

subroutine qshnorm(at,z,nsh,qshref)
      implicit none
      integer at,nsh
      real*8  z,qshref(10)

      integer i, j
      real*8  tot, pol, vale, core
      real*8  qshnew(10)     

      integer corel(10,86)

      corel=0
! this changes if the basis is changed!
! Li,Be 1s is semi-core
      corel(1, 3: 4)=1
! for gr.1/2 and TMs, shells 1-2 (2s) and 4 (2p) are (semi)core. 
      corel(1,11:12)=1
      corel(2,11:12)=1
      corel(4,11:12)=1
      corel(1,19:30)=1
      corel(2,19:30)=1
      corel(4,19:30)=1
      corel(1,37:48)=1
      corel(2,37:48)=1
      corel(4,37:48)=1
      corel(1,55:80)=1
      corel(2,55:80)=1
      corel(4,55:80)=1
      
      if(sum(corel(:,at)).eq.0) then
! no core shells, as before
         tot  = sum(qshref(1:nsh))
         qshnew(1:nsh) = qshref(1:nsh) * z / tot 
      else
! take core as is i.e. do not charge it and only re-normalize the valence
         core = 0
         tot  = sum(qshref(1:nsh))
         do i=1,nsh        
            core = core + qshref(i)*corel(i,at)
         enddo
         vale = tot - core
         do i=1,nsh
            if(corel(i,at).eq.0) then
               qshnew(i) = qshref(i) * ( z - core ) / vale
            else
               qshnew(i) = qshref(i)
            endif
         enddo
!        write(*,*) 'vale,core,all', vale,core,sum(qshnew)
      endif
      write(*,*) 'SUM',sum(qshnew(1:nsh))

      do i=1,nsh
         write(*,'(6x,''socc('',i2,'')='',F16.12)') i, qshnew(i)
      enddo
      do i=1,nsh
         write(*,'(2i4,F20.14,'' ##'')') at, i, qshnew(i)
      enddo

end

!! ------------------------------------------------------------------------
! returns row in PSE
!! ------------------------------------------------------------------------
INTEGER FUNCTION iTabRow6(i)
      implicit none
      INTEGER i
 
      iTabRow6=0
      If (i.gt. 0 .and. i.le. 2) Then
         iTabRow6=1
      Else If (i.gt. 2 .and. i.le.10) Then
         iTabRow6=2
      Else If (i.gt.10 .and. i.le.18) Then
         iTabRow6=3
      Else If (i.gt.18 .and. i.le.36) Then
         iTabRow6=4
      Else If (i.gt.36 .and. i.le.54) Then
         iTabRow6=5
      Else If (i.gt.54) Then
         iTabRow6=6
      End If
 
End

!! ------------------------------------------------------------------------
!  solve eigenvalue problem
!! ------------------------------------------------------------------------
subroutine solve2(mode,n,ndim,nel,nopen,homo,at,et,focc,H,S,P,e,eel,fail)
      use parcom
      use bascom, only: nsh
      use gtb_la, only : la_sygvx, la_sygvd
      implicit none
      integer mode,n,ndim,nel,nopen,homo
      integer at(n)
      real*8 H(ndim*(ndim+1)/2)
      real*8 S(ndim*(ndim+1)/2)
      real*8 P(ndim*(ndim+1)/2)
      real*8 e(ndim)
      real*8 focc(ndim)
      real*8 eel       
      real*8 et        
      real*8 gappar    
      logical fail

      integer i,j,info,lwork,liwork,ij,iu
      integer ihomoa,ihomob
      real*8 nfoda,nfodb,ga,gb,efa,efb,gap,w1,w0,t1,t0
      integer,allocatable ::iwork(:),ifail(:)
      real*8 ,allocatable ::D(:,:),hdum(:,:),sdum(:,:),work(:)
      real*8 ,allocatable ::focca(:), foccb(:)

      fail =.false.
      allocate (D(ndim,ndim),hdum(ndim,ndim),sdum(ndim,ndim),focca(ndim),foccb(ndim))  

      call blowsym(ndim,H,hdum)
      call blowsym(ndim,S,sdum)

      iu=min(homo+4,ndim)                 ! normal case 
      if(mode.eq.3) iu = min(ndim,5*homo) ! sTDA write, guess for highest relevant virt
      if(mode.eq.4) iu = ndim             ! all virts MUST be given to TM (otherwise virts are strange orthogonalized and results
                                          ! are worse! (i.e. iu = ndim)
! diag case branch
      if(iu.eq.ndim) then                  
! full diag (faster if all eigenvalues are taken)
       allocate (work(1),iwork(1))
       call la_sygvd(1,'V','U',ndim,hdum,ndim,sdum,ndim,e,work,-1,IWORK,LIWORK,INFO)
       lwork=int(work(1))
       liwork=iwork(1)
       deallocate(work,iwork)
       allocate (work(lwork),iwork(liwork)) 
       call la_sygvd(1,'V','U',ndim,hdum,ndim,sdum,ndim,e,work,LWORK,IWORK,LIWORK,INFO)
       D = hdum
      else
! for a large basis, taking only the occ.+few virt eigenvalues is faster than a full diag
       allocate(iwork(5*ndim),ifail(ndim),work(1))
       call la_sygvx(1,'V','I','U',ndim, hdum, ndim, sdum, ndim, ga, gb, &
     &               1, IU, 1d-7, ij, e, D, ndim, WORK, -1, IWORK, &
     &              IFAIL, INFO )
       lwork=idint(work(1))
       deallocate(work)
       allocate(work(lwork))          
       call la_sygvx(1,'V','I','U',ndim, hdum, ndim, sdum, ndim, ga, gb, &
     &               1, IU, 1d-7, ij, e, D, ndim, WORK, LWORK, IWORK, &
     &               IFAIL, INFO )
       do i=iu+1,ndim     
         e(i)=10d0
         D(1:ndim,i)=0d0
         D(i,i)     =1d0
       enddo
! end of diag case branch      
      endif

      if(info.ne.0) fail=.true.
      deallocate(hdum)

      if(mode.ge.2)then ! only if SP calc and in 2. iter
!        global shift to match DFT absolute orbital energies (this has no effect on anything in gTB)
         e = e + glob_par(15) + glob_par(17)*e
         open(unit=42,file='gtb_tmpmos',form='unformatted') 
         write(42) D
         write(42) e
      endif

      ga=0
      gb=0
! Fermi smearing                                          
!     convert restricted occ first to alpha/beta             
      if(nel.gt.0) then
         call occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
      else
         focca=0.0d0
         foccb=0.0d0
         ihomoa=0
         ihomob=0
      endif
      if(ihomoa+1.le.ndim) then 
         call FERMISMEAR(.false.,ndim,ihomoa,et,e,focca,nfoda,efa,ga)
      endif
      if(ihomob+1.le.ndim.and.nel.gt.1) then
         call FERMISMEAR(.false.,ndim,ihomob,et,e,foccb,nfodb,efb,gb)
      endif
      focc = focca + foccb

      call dmat(ndim,focc,D,sdum)
      call packsym(ndim,sdum,P)

      eel = sum(focc*e)

      end

!! ------------------------------------------------------------------------
!  response analogue of the twoscf routine
!! ------------------------------------------------------------------------

subroutine onescf(n,ndim,nel,nopen,homo,at,rab,cn,S,SS,Hmat,Hdiag,focc,&
                  eT,scfpar,ves0,psh,pa,P)                      
   use iso_fortran_env, only : wp => real64
   use bascom
   use parcom
   use com
   use gtb_la, only : la_gemm, la_symm
   implicit none 
!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
   integer, intent(in)    :: n                     ! number of atoms 
   integer, intent(in)    :: ndim                  ! number of AOs       
   integer, intent(in)    :: nel                   ! number of electrons 
   integer, intent(in)    :: nopen                 ! number of open shells
   integer, intent(in)    :: homo                  ! as the name says...
   integer, intent(in)    :: at(n)                 ! ordinal number of atoms
   real(wp),intent(in)    :: rab(n*(n+1)/2)        ! distances  
   real(wp),intent(in)    :: cn(n)                 ! CN           
   real(wp),intent(in)    :: S(ndim*(ndim+1)/2)    ! exact overlap maxtrix in SAO
   real(wp),intent(in)    :: SS(ndim*(ndim+1)/2)   ! scaled overlap maxtrix in SAO
   real(wp),intent(inout) :: Hmat(ndim*(ndim+1)/2) ! Vecp + field initialized
   real(wp),intent(in)    :: Hdiag(ndim)           ! diagonal of H0
   real(wp),intent(in)    :: focc (ndim)           ! fractional occ.
   real(wp),intent(in)    :: eT                    ! el. temp.
   real(wp),intent(in)    :: scfpar(8)             ! parameters
   real(wp),intent(in)    :: ves0(nsh)             ! ES potential field free      
   real(wp),intent(in)    :: psh(10,n)             ! shell populations with field
   real(wp),intent(in)    :: pa(n)                 ! atom      "         "   "

!! ------------------------------------------------------------------------
!  Output
!! ------------------------------------------------------------------------
   real(wp),intent(inout)   :: P (ndim*(ndim+1)/2) ! density matrix as a result of the field perturbation

!  local
   logical  :: fail
   integer  :: i,j,k,l,ish,ati,atj,ia,ib,jsh,ii,jj,lin,ij,li,iish,jjsh,mode
   real(wp) :: r,tmp,pol,hi,hj,hij,xk,t8,t9,qa,qb,keav,eh1,dmp,tmp2,ssh,eel
   real(wp) :: vi,vj         
   real(wp) :: gq(n),geff(n),ves(nsh),eps(ndim),xab(nsh,nsh)

   call setgab  (n,at,rab,pa,scfpar(7),xab)  ! the gab contain q as higher order effect on Ves, parameter different from twoscf
   call setespot(n,at,psh,xab,ves) 
   ves = ves * 0.5_wp

! H0 +  third-order (atomic charge exists in 1. AND 2. iter)
   do i=1, n
      geff(i) = pa(i)**2*shell_cnf4(1,at(i)) ! geff is temp.
   enddo

   ij = 0
   do i=1,ndim
      ia = aoat(i)
      ati= at(ia)
      ish= shell2ao(i)
      li = bas_lsh(ish,ati)+1
      hi = Hdiag(i)
      do j=1,i  
         ij = ij + 1
         ib = aoat(j)
         r  = rab(lin(ia,ib))
         if(r.gt.50d0) cycle
         hj = Hdiag(j)
         hij= hi+hj
         ssh= hij * SS(ij)
         atj= at(ib)
         if(ia.ne.ib) then            ! different atoms
            xk  = (shell_cnf4(2,ati)+shell_cnf4(2,atj)) 
            pol = ((hi-hj)/hij)**2
            keav= 0.5_wp*(shell_cnf2(9,ati) + shell_cnf2(9,atj))
            tmp = ssh * keav * (1_wp-pol*scfpar(8)) * (1_wp+xk/r) ! fit yields same values for iter1,2, parameter different from twoscf
         else                         ! same atoms
            jsh = shell2ao(j)
            if(ish.ne.jsh) then       ! s-s', p-p', d-d' off-diagonal, li=lj because S=0 otherwise
               tmp2= shell_cnf4(3+li,ati) 
               tmp = ssh * tmp2 + shell_cnf3(9,ati)* tmp2 * hij * SS(ij)**2
            else
               tmp = ssh 
            endif
         endif
!                                               third order diagonal
         Hmat(ij) = Hmat(ij) + tmp - S(ij)*(geff(ia)+geff(ib))
      enddo
   enddo

! H1
    k = 0
    do i=1,n
      gq(i) = 1_wp-(shell_xi(9,at(i))*pa(i)+shell_xi(10,at(i))*pa(i)**2) 
      hi = shell_cnf3(10,at(i)) + (cn(i)-avcn(at(i)))*shell_resp(10,at(i),1)
      do j=1,i
         k = k + 1
         r = hi + shell_cnf3(10,at(j)) + (cn(j)-avcn(at(j)))*shell_resp(10,at(j),1)
         t8= (rab(k)-r)/r
         xab(j,i) = 0.5_wp*(1_wp+erf(-2_wp*t8)) 
         xab(i,j) = xab(j,i)
      enddo
    enddo
    call calcpauli2(n,ndim,at,psh,S,Hdiag,Hmat) 

    k = 0
    do i=1,ndim
      ia = aoat(i)
      ati= at(ia)
      ish= shell2ao(i)
      iish=shmap(ish,ia)
      hi = shell_cnf2(ish,ati)*gq(ia)*expscal(3,10,ati)                  ! adapted +U scaling
      vi = ves(iish)*expscal(3,9,ati)+ves0(iish)*(1_wp-expscal(3,9,ati)) ! mixing of field perturbed and field free ES potential 
      do j=1,i                                                           ! with element wise parameter
         k  = k + 1
         ib = aoat(j)
         atj= at(ib)
         jsh= shell2ao(j)
         jjsh=shmap(jsh,ib)
         hj = shell_cnf2(jsh,atj)*gq(ib)*expscal(3,10,atj)
         vj = ves(jjsh)*expscal(3,9,atj)+ves0(jjsh)*(1_wp-expscal(3,9,atj))
!                            this part is INDO two-c like         shell ES
         Hmat(k) = Hmat(k) + P(k) * (hi + hj) * xab(ib,ia) - S(k)*(vi+vj) 
      enddo
    enddo

   call solve2 (1,n,ndim,nel,nopen,homo,at,eT,focc,Hmat,S,P,eps,eel,fail) 
   if(fail) stop 'diag error onescf'

end

!! ------------------------------------------------------------------------
!  MO match score for occupied space and scaled LUMO contribution
!! ------------------------------------------------------------------------

subroutine momatch(pr,ex,ndim,nocc,nvirt,S)
      use mocom  
      use gtb_la, only : la_symm, la_gemm
      implicit none          
      logical, intent(in)  :: pr,ex
      integer, intent(in)  :: ndim,nocc,nvirt
      real*8 , intent(in)  :: S(ndim*(ndim+1)/2)
      
      integer i,j,k,l,match
      real*8,allocatable :: SS(:,:), C2(:,:), C(:,:), e(:), srt(:)
      real*8 norm, wei, ss2, ss4, smax, gap_ref, homo_ref, ff, lumoweight
      real*8,parameter :: mothr      = 0.3d0    ! lower weight for MOs .lt. HOMO of this vaue
      real*8,parameter :: au2ev = 27.2113957d0 

      allocate(SS(ndim,ndim),C2(ndim,ndim),C(ndim,ndim),e(ndim),srt(nocc))

      lumoweight = 0.40d0   ! LUMO weight for fit
      ff = 1.d0
      if(increase_eps_weight) then
         ff = 3.d0 ! increase eps weight for metals to improve LS procedure
         lumoweight = 0.20d0
      endif

      if(pr) then
         write(*,*) 'running MO match with DFT reference', ff
         if(ex) write(*,*) 'excluding spread penalty because HIGHSYM was found'
         write(*,*) ' MO  DFT #   eps (DFT)    eps      penalty       Smax(1,2)'
      endif

      rewind 42   
      read  (42) C  ! gTB MOs, DFT on cmo_ref        (see mocom)
      read  (42) e  !  "  eigenvalues, DFT on epsref   "   2

      call blowsym(ndim,S,SS)
      CALL la_symm('L','L',ndim,ndim,1.D0,SS,ndim,C,ndim,0.D0,C2,ndim)  
      call la_gemm('T','N',ndim,ndim,ndim,1.0d0,cmo_ref,ndim,C2,ndim,0.0d0,SS,ndim)

      homo_ref=epsref(nocc)
      totmatch = 0
      do i=1, nocc
         norm = 0
         smax = 0
         ss4  = 0
         match= i
         do j=1, nocc
            ss2 = SS(j,i)**2
            srt(j) = -ss2
            if(homo_ref-epsref(j).gt.mothr) then
               wei = 0.4 * ff
            else
               wei = 1.2 * ff
            endif
            norm = norm + ss2 * abs(epsref(j) - e(i)) * wei
            if(ss2.gt.smax)then ! get best matching DFT MO number for printout
               smax = ss2
               match = j
            endif
            ss4 = ss4 + ss2**2
         enddo
         if(.not.ex)then
          if(homo_ref-epsref(i).gt.mothr) then
             norm = norm + (1d0-ss4)*0.02d0 ! add spread of overlap i.e. 10000... is best
          else
             norm = norm + (1d0-ss4)*0.06d0 ! 
          endif
         endif
         call qqsort(srt,1,nocc)
         if(pr)write(*,'(2i4,6f11.4)') i,match,epsref(match),e(i),norm,-srt(1),-srt(2)
         totmatch = totmatch + norm 
!        write(133,'(6f11.4)') epsref(match),e(i)
      enddo
!     write(133,'(6f11.4)') epsref(nocc+1),e(nocc+1)

      if(pr)write(*,*) 'MO match score occupied MOs   :', totmatch  

      gap_ref = (epsref(nocc+1) - epsref(nocc))*au2ev
      if(gap_ref.lt.13d0) then ! only reasonable gaps included in fit
         do j=nocc+1,nocc+2
          totmatch = totmatch + abs(epsref(j) - e(j)) * dble(nocc) * lumoweight  ! add LUMO deviation without assignment
         enddo                                                                   ! ie assume its right (because it should be the EA)
      endif

      if(pr) then
         write(*,*) 'total MO match score with gap :', totmatch  ! should be zero for perfect MOs
         write(*,'('' gap (eV)  DFT vs. gTB         : '',2f9.5)') gap_ref,(e(nocc+1) - e(nocc))*au2ev
      endif
      return
         
!     add some virt. but overlap only in penalty
      deallocate(srt)
      allocate(srt(ndim-nocc+1))
      do i=nocc + 1, nocc + nvirt
         k = 0
         ss4 = 0
         smax = 0
         do j=nocc + 1, ndim        
            ss2 = SS(j,i)**2
            k = k + 1
            srt(k) = -ss2
            ss4 = ss4 + ss2**2
            if(ss2.gt.smax)then 
               smax = ss2
               match = j
            endif
         enddo
         call qqsort(srt,1,ndim-nocc+1)
         if(abs(e(i)-epsref(match)).lt.0.2d0)then
             norm = (1d0-ss4)*0.05d0 
         else
             norm = 0d0
         endif
         totmatch = totmatch + norm
         if(pr)write(*,'(2i4,6f11.4)') i,match,epsref(match),e(i),norm,-srt(1),-srt(2)
      enddo 

      if(pr)then
         write(*,*) 'total MO match score incl virt:', totmatch  
      endif

      end

