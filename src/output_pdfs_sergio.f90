module mod_output_pdfs
! unused at the moment
#if defined(_USEPOST)
  use mpi
  use decomp_2d_io
  use cudecomp
  use mod_param     , only: dx,dy,dz,itot,jtot,ktot
  use mod_common_mpi, only:ierr,myid,offsets,i_pencil
  use mod_precision
  implicit none
  private
  public pdfs_sergio
  contains
  subroutine pdfs_sergio(fname,lo,hi,ng,z_g,u,v,w,p)
    use cudafor
    use mod_param, only: visc,lz
    implicit none
    integer, parameter :: nvars=5,npdf=200,nplanes=4
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3)  :: lo,hi,ng
    real(rp), intent(in), dimension(0:) :: z_g
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w,p
    real(rp), dimension(nvars) :: pdfmin,pdfmax,dupdf
    real(rp), dimension(nvars) :: uu
    integer :: ll1,ll2,ll3,ll4,ll5,llv(nvars)
    integer , dimension(nplanes*2) :: kplanes
    integer :: npts
    integer :: i,j,k,l,m,ll,mm,kk,kkk,ib
    real(rp), allocatable, save :: pdf(:,:,:),jpdf(:,:,:,:,:)
    real(rp) :: var1,var2
    character(len=1) :: navar
    integer :: rlen
    real(rp) :: rey
    !
    !  pdf of flow variables (ut, ur, uz, T, p)
    !
    pdfmin(1)=-0.25 - 1._rp ! convective reference frame
    pdfmin(2)=-0.25
    pdfmin(3)=-0.1
    pdfmin(4)=-0.1
    pdfmin(5)=-0.03
    pdfmax(1)= 0.25 - 1._rp ! convective reference frame
    pdfmax(2)= 0.25
    pdfmax(3)= 0.8
    pdfmax(4)= 0.8
    pdfmax(5)= 0.03
    !
    dupdf(:)=(pdfmax(:)-pdfmin(:))/npdf ! bin width
    !$acc wait
    !$acc enter data copyin(pdfmin,dupdf)
    !
    !  pdfs 
    !
    if(.not.allocated(pdf)) then
      allocate(pdf(npdf,ng(3),nvars))
      !$acc enter data create(pdf)
    end if
    !$acc enter data copyin(ng)
    !
    !$acc kernels default(present)
    pdf(:,:,:) = 0._rp
    !$acc end kernels
!   do k=lo(3),hi(3)
!     do j=lo(2),hi(2)
!       do i=lo(1),hi(1)
!         uu(1)=u(i,j,k)
!         uu(2)=v(i,j,k)
!         uu(3)=w(i,j,k)
!          uu(4)=dens(ic,jc,kc)
!         uu(4) = 0._rp
!         uu(5)=p(i,j,k)
!         do l=1,nvars
!           ll = int((uu(l)-pdfmin(l))/dupdf(l)) + 1
!           if (ll.gt.0.and.ll.le.npdf) pdf(ll,k,l)=pdf(ll,k,l) + 1._rp
!         end do
!       end do
!     end do
!   end do
    !$acc parallel loop collapse(3) default(present) private(llv)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          llv(1) = int((u(i,j,k)-pdfmin(1))/dupdf(1)) + 1
          llv(2) = int((v(i,j,k)-pdfmin(2))/dupdf(2)) + 1
          llv(3) = int((w(i,j,k)-pdfmin(3))/dupdf(3)) + 1
          llv(4) = int((   0._rp-pdfmin(4))/dupdf(4)) + 1
          llv(5) = int((p(i,j,k)-pdfmin(5))/dupdf(5)) + 1
          !$acc loop seq
          do l=1,nvars
            !$acc atomic update
            if (llv(l).gt.0.and.llv(l).le.npdf) pdf(llv(l),k,l)=pdf(llv(l),k,l) + 1._rp
            !$acc end atomic
          end do
        end do
      end do
    end do
    !
    !$acc kernels default(present)
    pdf(:,:,:)=pdf(:,:,:)/(1._rp*ng(1)*ng(2))
    !$acc end kernels
    !$acc update self(pdf)
    npts = nvars*ng(3)*npdf
    if(myid == 0) then
     call mpi_reduce(MPI_IN_PLACE,pdf(1,1,1),npts,MPI_REAL_RP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    else
     call mpi_reduce(  pdf(1,1,1),pdf(1,1,1),npts,MPI_REAL_RP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    end if
    !call mpi_allreduce(MPI_IN_PLACE,pdf(1,1,1),nvars*ng(3)*npdf,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    !
    if (myid.eq.0) then
      inquire(iolength=rlen) 1._rp
      open(10,file=trim(fname)//'_pdf.bin',access='direct',recl=npdf*ng(3)*nvars*rlen)
      write(10) pdf
      close(10)
      open(10,file=fname//'_pdf_meta.out',form='formatted')
        write(10,100) 1._rp*npdf,(pdfmin(i),pdfmax(i),i=1,5)
      close(10)
      !
      !do l=1,nvars
      !  write(navar,'(I1.1)') l
      !  open(10,file=trim(fname)//'_pdf_var'//navar//'.out',form='formatted')
      !  write(10,*) 'zone i=',npdf,', j=',ng(3)
      !  do k=1,ng(3)
      !    do i=1,npdf
      !      uu(l) = pdfmin(l)+dupdf(l)*(i-0.5)
      !      write(10,100) uu(l),z_g(k),pdf(i,k,l)
      !    end do
      !  end do
      !  close(10)
      !end do
    end if
    !
    ! joint pdfs at control planes
    !
    !
    ! define radial positions for theta-z planes
    !
    rey = 1.*lz/visc
    kplanes(1) = 1  ! first off-wall cell
    do k=1,ng(3)
      if (z_g(k)/(lz/2._rp)*(0.09*rey**0.88)>15) exit ! z^+ ~= 15 ! assuming lz = 2h
    end do
    kplanes(2) = k
    do k=1,ng(3)
      if (z_g(k).gt.0.2) exit ! z = 0.2
    end do
    kplanes(3) = k
    do k=1,ng(3)
      if (z_g(k).gt.0.5) exit ! z = 0.5
    end do
    kplanes(4) = k
    do kk=1,nplanes ! symmetric planes
      kplanes(nplanes+kk) = ng(3)-kplanes(nplanes-kk+1)+1
    end do
    !
    if (myid.eq.0) then
      open(10,file=fname//'_jpdf_meta.out',form='formatted')
      do l=1,nplanes
        write(10,100) 1._rp*npdf,1._rp*kplanes(l),z_g(kplanes(l)),(pdfmin(i),pdfmax(i),i=1,5)
      end do
      close(10)
    end if
    if(.not.allocated(jpdf)) then
      allocate(jpdf(npdf,npdf,nplanes,nvars,nvars))
      !$acc enter data create(jpdf)
    end if
    ! 
    !$acc kernels default(present)
    jpdf(:,:,:,:,:) = 0._rp
    !$acc end kernels
    do kk = 1,nplanes*2
      k=kplanes(kk)
      kkk = min(kk,nplanes*2-kk+1)
      if(k < lo(3) .or. k > hi(3)) cycle
      !$acc parallel loop collapse(2) default(present) private(llv,l,m)
      do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          llv(1)=int((0.5_rp*(u(i,j,k)+u(i-1,j,k))-pdfmin(1))/dupdf(1)) + 1
          llv(2)=int((0.5_rp*(v(i,j,k)+v(i,j-1,k))-pdfmin(2))/dupdf(2)) + 1
          llv(3)=int((0.5_rp*(w(i,j,k)+w(i,j,k-1))-pdfmin(3))/dupdf(3)) + 1
          llv(4)=int((0._rp                       -pdfmin(4))/dupdf(4)) + 1
          llv(5)=int((p(i,j,k)                    -pdfmin(5))/dupdf(5)) + 1
          !$acc loop seq
          do l=1,nvars
            !$acc loop seq
            do m=l+1,nvars
              if ((llv(l) >= 1.and.llv(l) <= npdf) .and. (llv(m) >= 1.and.llv(m) <= npdf)) then
                !$acc atomic update
                jpdf(llv(l),llv(m),kkk,l,m)=jpdf(llv(l),llv(m),kkk,l,m)+1._rp
                !$acc end atomic
              end if
            end do
          end do
        end do
      end do
    end do
    !$acc kernels default(present)
    jpdf(:,:,:,:,:)=jpdf(:,:,:,:,:)/(1._rp*ng(1)*ng(2)*2)
    !$acc end kernels
    !
    !$acc update self(jpdf)
    npts = nvars*nvars*nplanes*npdf*npdf
    if(myid == 0) then
       call mpi_reduce(MPI_IN_PLACE,jpdf(1,1,1,1,1),npts,MPI_REAL_RP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      else
       call mpi_reduce(jpdf(1,1,1,1,1),jpdf(1,1,1,1,1),npts,MPI_REAL_RP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      end if
    !call mpi_allreduce(MPI_IN_PLACE,jpdf(1,1,1,1,1),npts,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    !
    if (myid.eq.0) then
      inquire(iolength=rlen) 1._rp
      open(10,file=trim(fname)//'_jpdf.bin',access='direct',recl=npdf*npdf*nplanes*nvars*nvars*rlen)
      write(10) jpdf
      close(10)
      !!  print *,'done with jpdf'
      !open(10,file=trim(fname)//'_jpdfrz.out',form='formatted')
      !write(10,*) 'zone i=',npdf,', j=',npdf
      !k=3
      !do m=1,npdf
      ! var2 = pdfmin(3)+dupdf(3)*(m-0.5)
      !  do l=1,npdf
      !    var1 = pdfmin(2)+dupdf(2)*(l-0.5)
      !    write(10,100) var1, var2, jpdf(l,m,k,2,3)
      !  end do
      !end do
      !close(10)
    end if
    !$acc exit data delete(pdfmin,dupdf)
    !$acc exit data delete(ng)
    !$acc wait
100 format(40E20.10)
  end subroutine pdfs_sergio
#endif
end module mod_output_pdfs
