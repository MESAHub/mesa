
      
      subroutine read_inlist_xtrans_mesh_factor(ierr)
         integer, intent(out) :: ierr
         character (len=256) :: filename, message
         integer :: unit         
         filename = 'inlist_xtrans_mesh_factor'
         
         write(*,*) 'read_inlist_xtrans_mesh_factor'
         
         ! set defaults
         xtra_coef_above_xtrans = 1d0
         xtra_coef_below_xtrans = 1d0
         xtra_dist_above_xtrans = 0d0
         xtra_dist_below_xtrans = 0d0

         open(newunit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open control namelist file ', trim(filename)
         else
            read(unit, nml=xtrans_mesh_factor, iostat=ierr)  
            close(unit)
            if (ierr /= 0) then
               write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
               write(*, '(a)') &
                  'The following runtime error message might help you find the problem'
               write(*, *) 
               open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=xtrans_mesh_factor)
               close(unit)
            end if  
         end if

      end subroutine read_inlist_xtrans_mesh_factor

      
      subroutine other_mesh_delta_coeff_factor(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, nz, prev_jmax, jmax
         real(dp) :: xtra_coef, xtra_dist, Hp, r_extra
         logical, parameter :: dbg = .false.
         include 'formats'
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         !write(*,2) 'other_mesh_delta_coeff_factor s% model_number', s% model_number
         if (s% model_number < 20 .or. s% model_number > 40) return
         s% mesh_delta_coeff_factor(:) = 0.75d0
         !write(*,*) 'set other_mesh_delta_coeff_factor'

         if (xtra_coef_above_xtrans >= 1d0 .and. xtra_coef_below_xtrans >= 1d0) return
         
         if (dbg) write (*,*) 'do_xtrans_coef'

         nz = s% nz

         xtra_coef = xtra_coef_below_xtrans
         xtra_dist = xtra_dist_below_xtrans
         if (xtra_coef /= 1) then ! first go from surface to center doing below transitions
            prev_jmax = maxloc(s% xa(:,1),dim=1)
            k = 2
            do while (k <= nz)
               jmax = maxloc(s% xa(:,k),dim=1)
               if (jmax /= prev_jmax) then ! at a change in max abundance species
                  prev_jmax = jmax
                  Hp = s% Peos(k)/(s% rho(k)*s% grav(k))
                  r_extra = max(0d0, s% r(k) - xtra_dist*Hp)
                  do ! toward center from a transition
                     if (s% r(k) < r_extra) exit
                     if (xtra_coef < s% mesh_delta_coeff_factor(k)) &
                        s% mesh_delta_coeff_factor(k) = xtra_coef
                     if (k == nz) exit
                     k = k+1
                  end do
               end if
               k = k+1
            end do
         end if

         xtra_coef = xtra_coef_above_xtrans
         xtra_dist = xtra_dist_above_xtrans
         if (xtra_coef /= 1) then ! now go from center to surface doing above transitions
            prev_jmax = maxloc(s% xa(:,nz),dim=1)
            k = nz-1
            do while (k >= 1)
               jmax = maxloc(s% xa(:,k),dim=1)
               if (jmax /= prev_jmax) then ! at a change in max abundance species
                  prev_jmax = jmax
                  Hp = s% Peos(k)/(s% rho(k)*s% grav(k))
                  r_extra = max(0d0, s% r(k) + xtra_dist*Hp)
                  do ! toward surface from a transition
                     if (s% r(k) > r_extra) exit
                     if (xtra_coef < s% mesh_delta_coeff_factor(k)) &
                        s% mesh_delta_coeff_factor(k) = xtra_coef
                     if (k == 1) exit
                     k = k-1
                  end do
               end if
               k = k-1
            end do
         end if
         
      end subroutine other_mesh_delta_coeff_factor
