program gTB

   use bascom   ! AO basis
   use cbascom  ! core AO basis
   use parcom   ! TB method parameters
   use com      ! general stuff
   use mocom    ! ref MOs for fit and momatch value
   use dftd4
   use pgtb_
   use density_export, only: write_density_export
   use iso_fortran_env, only : wp => real64
   use purification_settings, only: tPurificationSet
   use default_files, only: default_atompara
   implicit none

   real(wp),allocatable :: xyz(:,:),rab(:),z(:), wbo(:,:), cn(:)
   real(wp),allocatable :: psh(:,:),q(:), psh_ref(:,:), q_ref(:), wbo_ref(:,:), qd4(:)
   real(wp),allocatable :: S(:),P(:),F(:),D3(:,:)
   real(wp),allocatable :: eps(:),focc(:),xnorm(:)
   real*4  ,allocatable :: ML1(:,:),ML2(:,:)

   integer, allocatable :: at(:)

   integer :: n
   integer ndim
   integer nopen
   integer na,nb,nel,ihomo
   integer prop
   integer i,j,ns,nf,nl
   integer :: par_idx
   real(wp) chrg ! could be fractional for model systems

   real(wp) t1,w1,t00,w00
   real(wp) pnt(3),dip(3)
   real(wp) alp(6)
   real(wp) efield(3)
   real(wp) floats(10)
   character*80 str(10)
   logical :: ex
   logical :: stda
   logical :: write_denmat
   logical :: logicals(10)
   logical :: used_default_par

   character(len=256)      :: fname,pname,bname,atmp,arg1,denmat_name
   character(len=220)      :: pline

   !> Handle purification
   type(tPurificationSet), allocatable :: pur

   call timing(t00,w00)

   stda    =.false.
   write_denmat = .false.
   prop = 1
   pnt  = 0
   chrg = 0
   nopen= 0
   pname='~/.atompara'
   bname="~/.basis_vDZP"
   denmat_name = 'ptb.denmat'

   expscal=0
   expscal(4,1:10,1:86)=1.0d0  ! back to standard exp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! input options
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do i = 1, command_argument_count()
      call getarg(i, arg1)
      select case(trim(arg1))
       case('-help')
         call help
         stop
       case('-version')
         call head
         stop
      end select
   end do

   if (command_argument_count() == 0) then
      call help
      error stop
   end if

   call head
   call getarg(1,fname)
   do i=2,command_argument_count()
      call getarg(i,arg1)
      if(index(arg1,'-stda')   .ne.0)stda=.true.  ! stda write
      if(index(arg1,'-purify').ne.0) then ! purification modus
         allocate(pur)
      endif
      if(index(arg1,'-par').ne.0)then
         call getarg(i+1,pname)
      endif
      if(index(arg1,'-bas').ne.0)then
         call getarg(i+1,bname)
      endif
      if(index(arg1,'-denmat').ne.0)then
         write_denmat = .true.
         call getarg(i+1,denmat_name)
      endif
      if(index(arg1,'-chrg').ne.0)then
         call getarg(i+1,atmp)
         call readline(atmp, floats, str, logicals, ns, nf, nl)
         chrg=floats(1)
      endif
      if(index(arg1,'-uhf').ne.0)then
         call getarg(i+1,atmp)
         call readline(atmp, floats, str, logicals, ns, nf, nl)
         nopen=floats(1)
      endif
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read parameter file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   inquire(file=pname,exist=ex)
   used_default_par = .not.ex
   if (used_default_par) then
      write(*,'(a)') "Parameter file '"//trim(pname)//"' not found, using built-in defaults."
   else
      write(*,*) pname
      open(unit=1,file=pname)
   end if
   par_idx = 0
   pline = parline(); read(pline,*) glob_par (1:10)
   pline = parline(); read(pline,*) glob_par(11:20)
   do i=1,86 ! change here for new elements
      pline = parline(); read(pline,*) j
      pline = parline(); read(pline,*) ener_par1 (1:10,j)    ! 1 -10
      pline = parline(); read(pline,*) ener_par2 (1:10,j)    ! 11-20
      pline = parline(); read(pline,*) expscal  (1,1:10,j)   ! 21-30
      pline = parline(); read(pline,*) ener_par6 (1:10,j)    ! 31-40
      pline = parline(); read(pline,*) ener_par4 (1:10,j)    ! 41-50
      pline = parline(); read(pline,*) ener_par5 (1:10,j)    ! 51-60
      pline = parline(); read(pline,*) expscal  (2,1:10,j)   ! 61-70  PTB
      pline = parline(); read(pline,*) shell_xi  (1:13,j)    ! slots 1-9 real per-shell, 10-13 special  "
      pline = parline(); read(pline,*) shell_cnf1(1:13,j)    ! "
      pline = parline(); read(pline,*) shell_cnf2(1:13,j)    ! "
      pline = parline(); read(pline,*) shell_cnf3(1:13,j)    ! "
      pline = parline(); read(pline,*) expscal  (3,1:10,j)   ! 111-120  "
      pline = parline(); read(pline,*) shell_cnf4(1:10,j)    ! 121-130  "
      pline = parline(); read(pline,*) shell_resp(1:13,j,1)  ! "
      pline = parline(); read(pline,*) shell_resp(1:13,j,2)  ! "
   enddo
   if (.not.used_default_par) close(1)

! mol. charge
   inquire(file='.CHRG',exist=ex)
   if(ex)then
      open(unit=1,file='.CHRG')
      read(1,'(a)')atmp
      close(1)
      call readline(atmp, floats, str, logicals, ns, nf, nl)
      chrg=floats(1)
   endif
! electric field
   inquire(file='.EFIELD',exist=ex)
   efield=0
   if(ex)then
      open(unit=1,file='.EFIELD')
      read(1,'(a)')atmp
      close(1)
      call readline(atmp, floats, str, logicals, ns, nf, nl)
      if(nf.lt.3) stop '.EFIELD read error'
      efield(1:3)=floats(1:3)
      write(*,'(''.EFIELD :'',3f12.6)') efield
   endif
   inquire(file='.UHF',exist=ex)
   if(ex)then
      open(unit=1,file='.UHF')
      read(1,'(a)')atmp
      close(1)
      call readline(atmp, floats, str, logicals, ns, nf, nl)
      nopen=int(floats(1))
   endif

! how many atoms?
   call rd0(fname,n)

   allocate(at(n),xyz(3,n),z(n),q(n),cn(n),rab(n*(n+1)/2),        &
   &         psh(10,n),psh_ref(10,n),q_ref(n),qd4(n),wbo(n,n),     &
   &         wbo_ref(n,n))

! read coordinates
   call rd(.true.,fname,n,xyz,at)
   call calcrab(n,at,xyz,rab)

! lanthanides currently use experimental PTB parameters. Warn loudly so
! results are never mistaken for validated production PTB output.
   do i=1,n
      if (at(i).ge.57.and.at(i).le.71) then
         print '(a,i0,a,i0,a)', "WARNING: atom ",i," is Z=",at(i), &
         & " (La-Lu). PTB lanthanide parameters are experimental; treat "// &
         & "density results as a fitting baseline until this element range "// &
         & "has been revalidated."
      end if
   enddo

   call setavcn   ! av. el. CNs with erfs=-2.0

   increase_eps_weight = .false.
   do i=1,n
      if (metal(at(i)).ne.0) increase_eps_weight = .true.  ! increase orbital energy weight in fit
      z(i)=valel(at(i))
   enddo
   nel=int(sum(z))-int(chrg)

   ndim=0
   call rdbas(bname)                      ! file: ~/.basis_vDZP
   write(*,*) 'basis read done.'
   call setupbas0(n,at,ndim)

   allocate(S(ndim*(ndim+1)/2),P(ndim*(ndim+1)/2),F(ndim*(ndim+1)/2), &
   &         D3(ndim*(ndim+1)/2,3),ML1(ndim,ndim),ML2(ndim,ndim),xnorm(ndim),focc(ndim),eps(ndim))

! valence basis
   call setupbas (n,at,ndim)
   write(*,*) 'basis setup done.'
   write(*,*) 'Ndim',ndim

! core basis
   call setupcbas0(n,at)
   call setupcbas (n,at)

! determine occupations
   call occ(ndim,nel,nopen,ihomo,na,nb,focc)
   write(*,*)
   write(*,*) 'nalpha ',na
   write(*,*) 'nbeta  ',nb
   write(*,*) 'ntotal ',na+nb
   if(nopen.eq.0.and.na.ne.nb) nopen = na - nb ! case .UHF does not exist i.e. radical
   write(*,*) 'nopen  ',nopen

! exact S
   call sint (n,ndim,at,xyz,rab,S,xnorm)
   call timing(t1,w1)
   call prtime(6,t1-t00,w1-w00,'startup and initial S')

   if(n.eq.1) call prmat(6,S,ndim,0,'overlap matrix')

   if(stda) prop = 4

! SINGLE POINT PTB
   if(prop.gt.0) call dipint(n,ndim,at,xyz,rab,xnorm,pnt,D3)! dipole integrals
   call pgtb(.true.,prop,n,ndim,nel,nopen,ihomo,at,chrg,xyz,z,rab,pnt,xnorm,S,D3,&
   &          efield,ML1,ML2,psh,q,P,F,eps,wbo,dip,alp, pur)
   if(write_denmat) call write_density_export(trim(denmat_name),n,ndim,at,xyz,P,S,xnorm)

   call timing(t1,w1)
   call prtime(6,t1-t00,w1-w00,'all')

contains

   function parline() result(line)
      character(len=220) :: line
      if (used_default_par) then
         par_idx = par_idx + 1
         line = default_atompara(par_idx)
      else
         read(1,'(a)') line
      end if
   end function parline

end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine help
   use, intrinsic :: iso_fortran_env, only : output_unit
   implicit none
   write(output_unit, '(a)') &
      "Usage: ptb <input> [options]...", &
      "", &
      "Geometry input must be given with the first argument.", &
      "Accepted formats are Turbomole coord and xyz format.", &
      "", &
      "Options:", &
      "", &
      "-chrg <int>        specify systems total charge", &
      "-uhf <int>         specify systems #open shells", &
      "-stda              output stda/TM compatible format", &
      "-par <file>        read parameters from provided file", &
      "-bas <file>        read basis set from provided file", &
      "-denmat <file>     write PTB density matrix and AO metadata", &
      "-purify            use density matrix purification instead of diagonalization", &
      "-version           print version header and exit", &
      "-help              show this help message", &
      ""
end subroutine help


subroutine head
   implicit none
   character(len=40),parameter:: date='08. Mar 2024'
   character(len=10),parameter:: version='3.7'

   write(*,*)
   write(*,'(7x,''=============================================='')')
   write(*,'(7x,''|                 P T B                      |'')')
   write(*,'(7x,''|                S.Grimme                    |'')')
   write(*,'(7x,''|          Universitaet Bonn, MCTC           |'')')
   write(*,'(7x,''=============================================='')')
   write(*,'(7x,''Version '',a,'', '',a)')trim(version),trim(date)
   write(*,*)
   write(*,'(7x,''Cite work conducted with this code as'')')
   write(*,'(7x,''S. Grimme, M. Mueller, A. Hansen, J. Chem. Phys., 2023, 158, 124111.'')')
   write(*,*) '                      DOI: 10.1063/5.0137838'
   write(*,*)
   write(*,*) '     This Version was made for use in NoSpherA2 and has been'
   write(*,*) '     modified to match requirements of IO, without altering the'
   write(*,*) '     original PTB model. The original PTB model is available at:'
   write(*,*) '                https://github.com/grimme-lab/ptb'
   write(*,*) '     The modified repository is available at:'
   write(*,*) '                https://github.com/AK-Kleemiss/ptb'
   write(*,*)

end
