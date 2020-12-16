
      
      subroutine read_inlist_convective_bdy(ierr)
         integer, intent(out) :: ierr
         character (len=256) :: filename, message
         integer :: unit         
         filename = 'inlist_convective_bdy'
         
         write(*,*) 'read_inlist_convective_bdy'
         
         ! set defaults
         convective_bdy_weight = 0
         convective_bdy_dq_limit = 1d-4
         convective_bdy_min_dt_yrs = 1d-3

         open(newunit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open control namelist file ', trim(filename)
         else
            read(unit, nml=convective_bdy, iostat=ierr)  
            close(unit)
            if (ierr /= 0) then
               write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
               write(*, '(a)') &
                  'The following runtime error message might help you find the problem'
               write(*, *) 
               open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=convective_bdy)
               close(unit)
            end if  
         end if

      end subroutine read_inlist_convective_bdy


      
      ! here is an example that adds a mesh function for log(opacity)
      subroutine how_many_other_mesh_fcns(id, n)
         integer, intent(in) :: id
         integer, intent(out) :: n
         n = 2
      end subroutine how_many_other_mesh_fcns
      
      
      subroutine other_mesh_fcn_data( &
            id, nfcns, names, gval_is_xa_function, vals1, ierr)
         integer, intent(in) :: id
         integer, intent(in) :: nfcns
         character (len=*) :: names(:)
         logical, intent(out) :: gval_is_xa_function(:) ! (nfcns)
         real(dp), pointer :: vals1(:) ! =(nz, nfcns)
         integer, intent(out) :: ierr
         integer :: nz, k
         real(dp), pointer :: vals(:,:)
         real(dp), parameter :: weight = 20
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (nfcns /= 2) stop 'bad nfcns for arg to other_mesh_fcn_data'
         gval_is_xa_function(1:nfcns) = .false.
         nz = s% nz
         vals(1:nz,1:nfcns) => vals1(1:nz*nfcns)
         names(1) = 'kap_function'
         do k=1,nz
            vals(k,1) = weight*log10(s% opacity(k))
         end do
         names(2) = 'conv_and_newly_nonconv'
         call do_conv_and_newly_nonconv(s, 2, vals, ierr)
      end subroutine other_mesh_fcn_data


      subroutine do_conv_and_newly_nonconv(s, i, vals, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: i ! mesh function index
         real(dp), pointer :: vals(:,:)
         integer, intent(out) :: ierr
         integer :: k, ktop, j, nz
         logical :: in_region, convective_below, convective_above

         logical, parameter :: dbg = .false.

         include 'formats'
         
         ierr = 0
         
         nz = s% nz
         vals(1:nz,i) = 0d0
         
         if (convective_bdy_weight == 0d0) return
         if (s% dt < convective_bdy_min_dt_yrs*secyer) return
         
         in_region = s% newly_nonconvective(1)
         ktop = 1
         convective_above = .false.

         if (dbg) then ! debug
            write(*,2) 'initial check', i
            do k=1,nz
               if (vals(k,i) /= 0) write(*,3) 'vals(k,i)', k, i, vals(k,i)
            end do
            do k=1,nz
               if (s% newly_nonconvective(k)) write(*,2) 's% newly_nonconvective(k)', k
            end do
            write(*,*)
         end if

         do k=2,nz
            if (in_region) then
               if (.not. s% newly_nonconvective(k)) then ! end of region
                  convective_below = (s% mixing_type(k) == convective_mixing)
                  if (dbg) write(*,3) 'call check_newly_nonconv_region', ktop, k-1
                  call check_newly_nonconv_region( &
                     s, vals, i, ktop, k-1, convective_below, convective_above)
                  in_region = .false.
               end if
            else if (s% newly_nonconvective(k)) then ! start of region
               if (dbg) write(*,2) 'start of newly_nonconvective region', k
               convective_above = (s% mixing_type(k-1) == convective_mixing)
               in_region = .true.
               ktop = max(1,k-1)
            end if
         end do

         if (dbg) then ! debug
            write(*,2) 'middle check', i
            do k=1,nz
               if (vals(k,i) /= 0) write(*,3) 'vals(k,i)', k, i, vals(k,i)
            end do
            write(*,*)
            stop 'do_conv_and_newly_nonconv'
         end if

         if (in_region) then
            convective_below = .false.
            if (dbg) write(*,3) 'call check_newly_nonconv_region', ktop, nz
            call check_newly_nonconv_region( &
               s, vals, i, ktop, nz, convective_below, convective_above)
         end if

         if (dbg) write(*,2) 's% num_conv_boundaries', s% num_conv_boundaries
         do j = 1, s% num_conv_boundaries
            k = s% conv_bdy_loc(j)
            if (dbg) write(*,3) 'convective boundaries', j, k
            vals(k,i) = convective_bdy_weight
            if (k < nz .and. s% top_conv_bdy(i)) then
               vals(k+1,i) = convective_bdy_weight
            else if (k > 1 .and. .not. s% top_conv_bdy(i)) then
               vals(k-1,i) = convective_bdy_weight
            end if
         end do

         if (dbg) then ! debug
            write(*,2) 'final check', i
            do k=1,nz
               if (vals(k,i) /= 0) write(*,3) 'vals(k,i)', k, i, vals(k,i)
            end do
            write(*,*)
            stop 'do_conv_and_newly_nonconv'
         end if

         do k=2,nz
            vals(k,i) = vals(k,i) + vals(k-1,i)
            !if (vals(k,i) > 0d0) write(*,2) 'conv mesh fcn', k, vals(k,i)
         end do

      end subroutine do_conv_and_newly_nonconv


      subroutine check_newly_nonconv_region( &
            s, vals, i, ktop, kbot, convective_below, convective_above)
         type (star_info), pointer :: s
         real(dp), pointer :: vals(:,:)
         integer, intent(in) :: i, ktop, kbot
         logical, intent(in) :: convective_below, convective_above
         integer :: k
         real(dp) :: qtop, qbot, dq, q_start, q_end
         include 'formats'
         qtop = s% q(ktop)
         qbot = s% q(kbot)
         dq = sum(s% dq(ktop:kbot-1))
         q_start = qtop
         q_end = qbot
         q_end = max(0d0, q_end - dq)
         q_start = min(1d0, q_start + dq)
         do k=1,s% nz
            if (s% q(k) <= q_start .and. s% q(k) >= q_end) then
               if (s% dq(k) >= convective_bdy_dq_limit) then
                  vals(k,i) = convective_bdy_weight
               else if (convective_bdy_dq_limit > 0) then
                  vals(k,i) = convective_bdy_weight * &
                     s% dq(k)/convective_bdy_dq_limit
               end if
            end if
         end do
      end subroutine check_newly_nonconv_region
