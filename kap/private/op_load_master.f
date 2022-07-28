      module op_load_master

      use const_def, only: dp


      implicit none
      logical :: loaded_op_master = .false.
      contains



      subroutine load_op_master(emesh_data_for_op_mono_path, iz,ite,jne,epatom,amamu,sig,eumesh,ierr)

      character (len=*), intent(in) :: emesh_data_for_op_mono_path

      integer, intent(inout) :: ierr
      integer , pointer, intent(out) :: iz(:),ite(:),jne(:)
      real(dp), pointer, intent(out) :: sig(:,:,:)
      real(dp), pointer, intent(out):: epatom(:,:),amamu(:),eumesh(:,:,:)

      integer :: n, m, ke, ik
      CHARACTER(LEN=72) :: FMT
      integer :: nel, nptot, np
      parameter(nel = 17, nptot = 10000, np=1648) !number of elements and number of u-mesh points.
      real(dp), allocatable :: amamu_f(:,:)
      integer, allocatable  :: iz_f(:,:)


      if (loaded_op_master) return

      allocate(iz_f(nel,np),iz(nel),ite(np),jne(np),stat=ierr)
      allocate(sig(nel,np,nptot),stat=ierr)
      allocate(epatom(nel,np),amamu_f(nel,np),amamu(nel),eumesh(nel,np,nptot),stat=ierr)

      FMT = '(i2,1x,i3,1x,i3,1x,F14.10,1x,F14.10,10000(1x,E12.6E3),10000(1x,E13.6E3))'

      write(*,*) 'Opening file...'
      open(1,file = emesh_data_for_op_mono_path,form = 'formatted', action ='read')
      write(*,*) 'Loading OP mono data...'

      do ke =1, nel
        do n =1,np
        read(1,FMT)iz_f(ke,n),ite(n),jne(n),epatom(ke,n),amamu_f(ke,n),(sig(ke,n,m), m=1,nptot),(eumesh(ke,n,m), m=1,nptot)
        enddo
      enddo

      close(1)

      do ke=1,nel
        amamu(ke) = amamu_f(ke,1)
        iz(ke)    = iz_f(ke,1)
      enddo

      write(*,*) 'OP mono data loaded.'
      ierr = 0
      loaded_op_master = .true.
      end subroutine load_op_master
      end module op_load_master
