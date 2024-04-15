! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_sanity
  use, intrinsic :: iso_c_binding, only: C_PTR
  use mpi
  use decomp_2d
  use mod_bound     , only: boundp,bounduvw,updt_rhs_b
  use mod_chkdiv    , only: chkdiv
  use mod_common_mpi, only: myid,ierr,ipencil_axis
  use mod_correc    , only: correc
  use mod_debug     , only: chk_helmholtz
  use mod_fft       , only: fftend
  use mod_fillps    , only: fillps
  use mod_initflow  , only: add_noise
  use mod_initmpi   , only: initmpi
  use mod_initsolver, only: initsolver
  use mod_param     , only: small
#if !defined(_OPENACC)
  use mod_solver    , only: solver
#else
  use mod_solver_gpu, only: solver => solver_gpu
#endif
  use mod_wmodel    , only: cmpt_bcuvw,cmpt_bcp
  use mod_typedef   , only: cond_bound
  use mod_precision
  implicit none
  private
  public test_sanity_input
  ! public test_sanity_input,test_sanity_solver
  contains
  subroutine test_sanity_input(ng,dims,sgstype,stop_type,cbcvel,cbcpre,cbcsgs,bcvel,bcpre,bcsgs,n,is_bound,lwm,l,zc,dl,h,is_forced)
    !
    ! performs some a priori checks of the input files before the calculation starts
    !
    implicit none
    integer , intent(in), dimension(3) :: ng
    integer , intent(in), dimension(2) :: dims
    character(len=*), intent(in) :: sgstype
    logical , intent(in), dimension(3) :: stop_type
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3)   :: cbcpre,cbcsgs
    real(rp), intent(in), dimension(0:1,3,3) :: bcvel
    real(rp), intent(in), dimension(0:1,3) :: bcpre,bcsgs
    integer, intent(in), dimension(3) :: n
    logical, intent(in), dimension(0:1,3) :: is_bound
    integer, intent(in), dimension(0:1,3) :: lwm
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in) :: h
    logical, intent(in), dimension(3) :: is_forced
    logical :: passed
    !
    call chk_dims(ng,dims,cbcvel,sgstype,passed); if(.not.passed) call abortit
    call chk_stop_type(stop_type,passed); if(.not.passed) call abortit
    call chk_bc(cbcvel,cbcpre,cbcsgs,bcvel,bcpre,bcsgs,n,is_bound,lwm,l,zc,dl,h,passed); if(.not.passed) call abortit
    call chk_forcing(cbcpre,is_forced,passed); if(.not.passed) call abortit
#if defined(_IMPDIFF_1D) && !defined(_IMPDIFF)
    if(myid == 0)  print*, 'ERROR: `_IMPDIFF_1D` cpp macro requires building with `_IMPDIFF` too.'; call abortit
#endif
#if defined(_IMPDIFF_1D) && !defined(_DECOMP_Z)
    if(myid == 0)  print*, 'WARNING: a run with implicit Z diffusion (`_IMPDIFF_1D`) is much more efficient &
                                   & when combined with a Z-pencils parallelization (`_DECOMP_Z`).'
#endif
  end subroutine test_sanity_input
  !
  subroutine chk_stop_type(stop_type,passed)
    implicit none
    logical, intent(in), dimension(3) :: stop_type
    logical, intent(out) :: passed
    passed = .true.
    if(.not.any(stop_type(:))) then
      if(myid == 0) print*, 'ERROR: stopping criterion not chosen.'
      passed = .false.
    end if
  end subroutine chk_stop_type
  !
  subroutine chk_dims(ng,dims,cbcvel,sgstype,passed)
    implicit none
    integer, intent(in), dimension(3) :: ng
    integer, intent(in), dimension(2) :: dims
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=*), intent(in) :: sgstype
    logical, intent(out) :: passed
    integer, dimension(2) :: ii
    logical :: passed_loc
    character(len=2) :: bc01v
    integer :: idir,ivel,i
    passed = .true.
    ii = pack([1,2,3],[1,2,3] /= ipencil_axis)
    passed_loc = all(dims(:)<=ng(ii)).and.all(dims(:)>=1)
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: 1 <= dims(:) <= [itot,jtot], or [itot,ktot], or [jtot ktot] depending on the decomposition.'
    passed = passed.and.passed_loc
    !
    if(sgstype=='smag') then
      passed_loc = .true.
      do i = 1,2
        idir = ii(i)
        ivel = ii(i)
        bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
        if(bc01v == 'DD') then
          passed_loc = passed_loc.and.(dims(i)<=2)
        end if
      end do
      if(myid == 0.and.(.not.passed_loc)) &
        print*, 'ERROR: more than two subdomains between two opposite walls.'
      passed = passed.and.passed_loc
    end if
    !
  end subroutine chk_dims
  !
  subroutine chk_bc(cbcvel,cbcpre,cbcsgs,bcvel,bcpre,bcsgs,n,is_bound,lwm,l,zc,dl,h,passed)
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3) :: cbcpre,cbcsgs
    real(rp)        , intent(in), dimension(0:1,3,3) :: bcvel
    real(rp)        , intent(in), dimension(0:1,3) :: bcpre,bcsgs
    integer, intent(in), dimension(3) :: n
    logical, intent(in), dimension(0:1,3) :: is_bound
    integer, intent(in), dimension(0:1,3) :: lwm
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in) :: h
    logical, intent(out) :: passed
    character(len=2) :: bc01v,bc01p,bc01s
    integer :: i,ivel,idir
    logical :: passed_loc
    passed = .true.
    !
    ! check validity of pressure and velocity BCs
    !
    passed_loc = .true.
    do ivel = 1,3
      do idir=1,3
        bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
        passed_loc = passed_loc.and.( (bc01v == 'PP').or. &
                                      (bc01v == 'ND').or. &
                                      (bc01v == 'DN').or. &
                                      (bc01v == 'NN').or. &
                                      (bc01v == 'DD') )
      end do
    end do
    if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: velocity BCs not valid.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.( (bc01p == 'PP').or. &
                                    (bc01p == 'ND').or. &
                                    (bc01p == 'DN').or. &
                                    (bc01p == 'NN').or. &
                                    (bc01p == 'DD') )
    end do
    if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: pressure BCs not valid.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      ivel = idir
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.( (bc01v == 'PP'.and.bc01p == 'PP').or. &
                                    (bc01v == 'ND'.and.bc01p == 'DN').or. &
                                    (bc01v == 'DN'.and.bc01p == 'ND').or. &
                                    (bc01v == 'DD'.and.bc01p == 'NN').or. &
                                    (bc01v == 'NN'.and.bc01p == 'DD') )
    end do
    if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: velocity and pressure BCs not compatible.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      bc01s = cbcsgs(0,idir)//cbcsgs(1,idir)
      passed_loc = passed_loc.and.( (bc01s == 'PP').or. &
                                    (bc01s == 'ND').or. &
                                    (bc01s == 'DN').or. &
                                    (bc01s == 'DD').or. &
                                    (bc01s == 'NN') )
    end do
    if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: sgs BCs not valid.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      ivel = idir
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      bc01s = cbcsgs(0,idir)//cbcsgs(1,idir)
      passed_loc = passed_loc.and.( (bc01v == 'PP'.and.bc01s == 'PP').or. &
                                    (bc01v == 'ND'.and.bc01s == 'DD').or. &
                                    (bc01v == 'DN'.and.bc01s == 'DD').or. &
                                    (bc01v == 'DD'.and.bc01s == 'DD').or. &
                                    (bc01v == 'NN'.and.bc01s == 'DD') )
    end do
    if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: velocity and sgs BCs not compatible.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,2
      passed_loc = passed_loc.and.((bcpre(0,idir) == 0.).and.(bcpre(1,idir) == 0.))
    end do
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: pressure BCs in directions x and y must be homogeneous (value = 0.).'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir = 1,3
      do i = 0,1
        if(lwm(i,idir)/=0) then
          do ivel = 1,3
            if(ivel==idir) then
              passed_loc = passed_loc.and.cbcvel(i,idir,ivel)=='D'
            else
              passed_loc = passed_loc.and.cbcvel(i,idir,ivel)=='N'
            end if
          end do
        end if
      end do
    end do
    if(myid == 0.and.(.not.passed_loc)) &
    print*, 'ERROR: wall model BCs must be Neumann (wall-parallel) and Dirichlet (wall-normal).'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    if(is_bound(0,1).and.lwm(0,1)/=0) passed_loc = passed_loc.and.(h>0.5_rp*dl(1) .and.h<(n(1)-0.5_rp)*dl(1))
    if(is_bound(1,1).and.lwm(1,1)/=0) passed_loc = passed_loc.and.(h>0.5_rp*dl(1) .and.h<(n(1)-0.5_rp)*dl(1))
    if(is_bound(0,2).and.lwm(0,2)/=0) passed_loc = passed_loc.and.(h>0.5_rp*dl(2) .and.h<(n(2)-0.5_rp)*dl(2))
    if(is_bound(1,2).and.lwm(1,2)/=0) passed_loc = passed_loc.and.(h>0.5_rp*dl(2) .and.h<(n(2)-0.5_rp)*dl(2))
    if(is_bound(0,3).and.lwm(0,3)/=0) passed_loc = passed_loc.and.(h>zc(1)        .and.h<zc(n(3))  )
    if(is_bound(1,3).and.lwm(1,3)/=0) passed_loc = passed_loc.and.(h>l(3)-zc(n(3)).and.h<l(3)-zc(1))
    if(.not.passed_loc) print*, 'ERROR: invalid wall model height.'
    passed = passed.and.passed_loc
    !
#if defined(_IMPDIFF) && !defined(_IMPDIFF_1D)
    !implicit: bcs in x and y must be 'D0' or 'P',
    !so implicit mode can not be used for WMLES of square ducts
    passed_loc = .true.
    do ivel = 1,3
      do idir=1,2
        bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
        passed_loc = passed_loc.and.(bc01v /= 'NN')
      end do
    end do
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: Neumann-Neumann velocity BCs with implicit diffusion currently not supported in x and y; only in z.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do ivel = 1,3
      do idir=1,2
        passed_loc = passed_loc.and.((bcvel(0,idir,ivel) == 0.).and.(bcvel(1,idir,ivel) == 0.))
      end do
    end do
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: velocity BCs with implicit diffusion in directions x and y must be homogeneous (value = 0.).'
    passed = passed.and.passed_loc
#endif
#if defined(_OPENACC)
    do idir=1,2
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and..not.( (bc01p == 'DN').or. &
                                         (bc01p == 'ND') )
    end do
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: pressure BCs "ND" or "DN" along x or y not implemented on GPUs yet.'
#endif
  end subroutine chk_bc
  !
  subroutine chk_forcing(cbcpre,is_forced,passed)
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbcpre
    logical         , intent(in), dimension(3) :: is_forced
    logical         , intent(out) :: passed
    integer :: idir
    passed = .true.
    !
    ! 1) check for compatibility between pressure BCs and flow forcing
    !
    do idir=1,3
      if(is_forced(idir)) then
        passed = passed.and.(cbcpre(0,idir)//cbcpre(1,idir) == 'PP')
      end if
    end do
    if(myid == 0.and.(.not.passed)) &
    print*, 'ERROR: Flow cannot be forced in a non-periodic direction; check the BCs and is_forced in `input.nml`.'
  end subroutine chk_forcing
  !
!   subroutine test_sanity_solver(ng,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z,dli,dzc,dzf,dzci,dzfi,dzci_g,dzfi_g, &
!                                 nb,is_bound,cbcvel,cbcpre,bcvel,bcpre)
! #if defined(_OPENACC)
!     use mod_workspaces, only: set_cufft_wspace
!     use openacc,        only: acc_get_cuda_stream
! #endif
!     implicit none
!     integer , intent(in), dimension(3) :: ng,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z
!     real(rp), intent(in), dimension(3) :: dli
!     real(rp), intent(in), dimension(0:) :: dzc,dzf,dzci,dzfi,dzci_g,dzfi_g
!     integer , intent(in), dimension(0:1,3) :: nb
!     logical , intent(in), dimension(0:1,3) :: is_bound
!     character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
!     character(len=1), intent(in), dimension(0:1,3)   :: cbcpre
!     real(rp), intent(in), dimension(0:1,3,3)         :: bcvel
!     real(rp), intent(in), dimension(0:1,3)           :: bcpre
!     real(rp), allocatable, dimension(:,:,:) :: u,v,w,p
!     type(cond_bound) :: bcu,bcv,bcw,bcp
! #if !defined(_OPENACC)
!     type(C_PTR), dimension(2,2) :: arrplan
! #else
!     integer    , dimension(2,2) :: arrplan
! #endif
!     real(rp) :: normfft
!     real(rp), allocatable, dimension(:,:) ::  lambdaxy
!     real(rp), allocatable, dimension(:) :: a,b,c,bb
!     real(rp), allocatable, dimension(:,:,:) :: rhsbx,rhsby,rhsbz
!     real(rp), dimension(3) :: dl
!     real(rp) :: dt,dti,alpha
!     real(rp) :: divtot,divmax,resmax
!     logical :: passed,passed_loc
!     passed = .true.
!     !$acc wait
!     allocate(u(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
!              v(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
!              w(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
!              p(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
!              lambdaxy(n_z(1),n_z(2)), &
!              a(n_z(3)),b(n_z(3)),c(n_z(3)), &
!              rhsbx(n(2),n(3),0:1), &
!              rhsby(n(1),n(3),0:1), &
!              rhsbz(n(1),n(2),0:1))
!     allocate(bcu%x(0:n(2)+1,0:n(3)+1,0:1), & 
!              bcv%x(0:n(2)+1,0:n(3)+1,0:1), &
!              bcw%x(0:n(2)+1,0:n(3)+1,0:1), &
!              bcu%y(0:n(1)+1,0:n(3)+1,0:1), &
!              bcv%y(0:n(1)+1,0:n(3)+1,0:1), &
!              bcw%y(0:n(1)+1,0:n(3)+1,0:1), &
!              bcu%z(0:n(1)+1,0:n(2)+1,0:1), &
!              bcv%z(0:n(1)+1,0:n(2)+1,0:1), &
!              bcw%z(0:n(1)+1,0:n(2)+1,0:1))
!     allocate(bcp%x(0:n(2)+1,0:n(3)+1,0:1), &
!              bcp%y(0:n(1)+1,0:n(3)+1,0:1), &
!              bcp%z(0:n(1)+1,0:n(2)+1,0:1))
!     !
!     ! initialize velocity below with some random noise
!     !
!     u(:,:,:) = 0.
!     v(:,:,:) = 0.
!     w(:,:,:) = 0.
!     p(:,:,:) = 0.
!     call add_noise(ng,lo,123,.5_rp,u(1:n(1),1:n(2),1:n(3)))
!     call add_noise(ng,lo,456,.5_rp,v(1:n(1),1:n(2),1:n(3)))
!     call add_noise(ng,lo,789,.5_rp,w(1:n(1),1:n(2),1:n(3)))
!     !$acc enter data copyin(u,v,w,p)
!     !
!     ! test pressure correction
!     !
!     call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcpre,bcpre(:,:), &
!                     lambdaxy,['c','c','c'],a,b,c,arrplan,normfft,rhsbx,rhsby,rhsbz)
!     !$acc enter data copyin(lambdaxy,a,b,c,rhsbx,rhsby,rhsbz)
!     !@acc call set_cufft_wspace(pack(arrplan,.true.),acc_get_cuda_stream(1))
!     dl  = dli**(-1)
!     dt  = acos(-1.) ! value is irrelevant
!     dti = dt**(-1)
!     call cmpt_bcuvw(cbcvel,n,bcvel,is_bound,u,v,w,bcu,bcv,bcw) !may call several times in this subroutine 
!     call cmpt_bcp  (cbcpre,n,bcpre,is_bound,p    ,bcp) !may call several times in this subroutine
!     call bounduvw(cbcvel,n,bcu,bcv,bcw,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
!     call fillps(n,dli,dzfi,dti,u,v,w,p)
!     call updt_rhs_b(['c','c','c'],cbcpre,n,is_bound,rhsbx,rhsby,rhsbz,p)
!     call solver(n,ng,arrplan,normfft,lambdaxy,a,b,c,cbcpre,['c','c','c'],p)
!     call boundp(cbcpre,n,bcp,nb,is_bound,dl,dzc,p)
!     call correc(n,dli,dzci,dt,p,u,v,w)
!     call bounduvw(cbcvel,n,bcu,bcv,bcw,nb,is_bound,.true.,dl,dzc,dzf,u,v,w)
!     call chkdiv(lo,hi,dli,dzfi,u,v,w,divtot,divmax)
!     passed_loc = divmax < small
!     if(myid == 0.and.(.not.passed_loc)) &
!     print*, 'ERROR: Pressure correction: Divergence is too large, with maximum = ', divmax
!     passed = passed.and.passed_loc
!     call fftend(arrplan)
! #if defined(_IMPDIFF) && !defined(_IMPDIFF_1D)
!     allocate(bb(n_z(3)))
!     alpha = acos(-1.) ! irrelevant
!     !$acc kernels default(present)
!     u(:,:,:) = 0.
!     v(:,:,:) = 0.
!     w(:,:,:) = 0.
!     !$acc end kernels
!     !$acc update self(u,v,w)
!     call add_noise(ng,lo,123,.5_rp,u(1:n(1),1:n(2),1:n(3)))
!     call add_noise(ng,lo,456,.5_rp,v(1:n(1),1:n(2),1:n(3)))
!     call add_noise(ng,lo,789,.5_rp,w(1:n(1),1:n(2),1:n(3)))
!     !$acc update device(u,v,w)
!     !$acc enter data create(bb)
!     call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,1),bcvel(:,:,1), &
!                     lambdaxy,['f','c','c'],a,b,c,arrplan,normfft,rhsbx,rhsby,rhsbz)
!     !$acc update device(lambdaxy,a,b,c,rhsbx,rhsby,rhsbz)
!     !@acc call set_cufft_wspace(pack(arrplan,.true.),acc_get_cuda_stream(1))
!     call bounduvw(cbcvel,n,bcu,bcv,bcw,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
!     !$acc kernels default(present)
!     u(:,:,:) = u(:,:,:)*alpha
!     p(:,:,:) = u(:,:,:)
!     bb(:) = b(:) + alpha
!     !$acc end kernels
!     call updt_rhs_b(['f','c','c'],cbcvel(:,:,1),n,is_bound,rhsbx,rhsby,rhsbz,u)
!     call solver(n,ng,arrplan,normfft,lambdaxy,a,bb,c,cbcvel(:,:,1),['f','c','c'],u)
!     call fftend(arrplan)
!     call bounduvw(cbcvel,n,bcu,bcv,bcw,nb,is_bound,.false.,dl,dzc,dzf,u,v,w) ! actually we are only interested in boundary condition in u
!     call chk_helmholtz(lo,hi,dli,dzci,dzfi,alpha,p,u,cbcvel(:,:,1),is_bound,['f','c','c'],resmax)
!     passed_loc = resmax < small
!     if(myid == 0.and.(.not.passed_loc)) &
!     print*, 'ERROR: wrong solution of Helmholtz equation in x direction.'
!     passed = passed.and.passed_loc
!     !
!     call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,2),bcvel(:,:,2), &
!                     lambdaxy,['c','f','c'],a,b,c,arrplan,normfft,rhsbx,rhsby,rhsbz)
!     !$acc update device(lambdaxy,a,b,c,rhsbx,rhsby,rhsbz)
!     !@acc call set_cufft_wspace(pack(arrplan,.true.),acc_get_cuda_stream(1))
!     call bounduvw(cbcvel,n,bcu,bcv,bcw,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
!     !$acc kernels default(present)
!     v(:,:,:) = v(:,:,:)*alpha
!     p(:,:,:) = v(:,:,:)
!     bb(:) = b(:) + alpha
!     !$acc end kernels
!     call updt_rhs_b(['c','f','c'],cbcvel(:,:,2),n,is_bound,rhsbx,rhsby,rhsbz,v)
!     call solver(n,ng,arrplan,normfft,lambdaxy,a,bb,c,cbcvel(:,:,2),['c','f','c'],v)
!     call fftend(arrplan)
!     call bounduvw(cbcvel,n,bcu,bcv,bcw,nb,is_bound,.false.,dl,dzc,dzf,u,v,w) ! actually we are only interested in boundary condition in v
!     call chk_helmholtz(lo,hi,dli,dzci,dzfi,alpha,p,v,cbcvel(:,:,2),is_bound,['c','f','c'],resmax)
!     passed_loc = resmax < small
!     if(myid == 0.and.(.not.passed_loc)) &
!     print*, 'ERROR: wrong solution of Helmholtz equation in y direction.'
!     passed = passed.and.passed_loc
!     !
!     call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,3),bcvel(:,:,3), &
!                     lambdaxy,['c','c','f'],a,b,c,arrplan,normfft,rhsbx,rhsby,rhsbz)
!     !$acc update device(lambdaxy,a,b,c,rhsbx,rhsby,rhsbz)
!     !@acc call set_cufft_wspace(pack(arrplan,.true.),acc_get_cuda_stream(1))
!     call bounduvw(cbcvel,n,bcu,bcv,bcw,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
!     !$acc kernels default(present)
!     w(:,:,:) = w(:,:,:)*alpha
!     p(:,:,:) = w(:,:,:)
!     bb(:) = b(:) + alpha
!     !$acc end kernels
!     call updt_rhs_b(['c','c','f'],cbcvel(:,:,3),n,is_bound,rhsbx,rhsby,rhsbz,w)
!     call solver(n,ng,arrplan,normfft,lambdaxy,a,bb,c,cbcvel(:,:,3),['c','c','f'],w)
!     call fftend(arrplan)
!     call bounduvw(cbcvel,n,bcu,bcv,bcw,nb,is_bound,.false.,dl,dzc,dzf,u,v,w) ! actually we are only interested in boundary condition in w
!     call chk_helmholtz(lo,hi,dli,dzci,dzfi,alpha,p,w,cbcvel(:,:,3),is_bound,['c','c','f'],resmax)
!     passed_loc = resmax < small
!     if(myid == 0.and.(.not.passed_loc)) &
!     print*, 'ERROR: wrong solution of Helmholtz equation in z direction.'
!     passed = passed.and.passed_loc
!     !$acc exit data delete(bb)
! #endif
!     !$acc exit data delete(u,v,w,p,lambdaxy,a,b,c,rhsbx,rhsby,rhsbz)
!     if(.not.passed) then
!       call decomp_2d_finalize
!       call MPI_FINALIZE(ierr)
!       error stop
!     end if
!   end subroutine test_sanity_solver
  !
  subroutine abortit
    implicit none
    if(myid == 0) print*, ''
    if(myid == 0) print*, '*** Simulation aborted due to errors in the input file ***'
    if(myid == 0) print*, '    check `input.nml`.'
    call decomp_2d_finalize
    call MPI_FINALIZE(ierr)
    error stop
  end subroutine abortit
end module mod_sanity
