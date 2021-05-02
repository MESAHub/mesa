
      
      subroutine read_inlist_overshoot_dbl_exp(ierr)
         integer, intent(out) :: ierr
         character (len=256) :: filename, message
         integer :: unit         
         filename = 'inlist_overshoot_dbl_exp'
         
         write(*,*) 'read_inlist_overshoot_dbl_exp'
         
         ! set defaults
         overshoot_f2 = 0d0
         overshoot_D2 = 0d0

         open(newunit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open control namelist file ', trim(filename)
         else
            read(unit, nml=overshoot_dbl_exp, iostat=ierr)  
            close(unit)
            if (ierr /= 0) then
               write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
               write(*, '(a)') &
                  'The following runtime error message might help you find the problem'
               write(*, *) 
               open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=overshoot_dbl_exp)
               close(unit)
            end if  
         end if

      end subroutine read_inlist_overshoot_dbl_exp


      subroutine eval_overshoot_dbl_exp (id, i, j, k_a, k_b, D, vc, ierr)

       integer, intent(in)      :: id
       integer, intent(in)      :: i
       integer, intent(in)      :: j
       integer, intent(out)     :: k_a
       integer, intent(out)     :: k_b
       real(dp), intent(out)    :: D(:)
       real(dp), intent(out)    :: vc(:)
       integer, intent(out)     :: ierr

       type(star_info), pointer :: s
       real(dp) :: f
       real(dp) :: f0
       real(dp) :: f2
       real(dp) :: D0
       real(dp) :: Delta0
       real(dp) :: D2
       real(dp) :: w
       real(dp) :: factor
       real(dp) :: Hp_cb
       integer  :: k_ob
       real(dp) :: r_ob
       real(dp) :: D_ob
       real(dp) :: vc_ob
       real(dp) :: dr2
       logical  :: outward
       integer  :: dk
       integer  :: k
       real(dp) :: r
       real(dp) :: dr

       ! Evaluate the overshoot diffusion coefficient D(k_a:k_b) and
       ! mixing velocity vc(k_a:k_b) at the i'th convective boundary,
       ! using the j'th set of overshoot parameters. The overshoot
       ! follows the exponential scheme description by Falk Herwig, "The
       ! evolution of AGB stars with convective overshoot", A&A, 360,
       ! 952-968 (2000), as modified with the second exponential tail
       ! described by by Battino et al., "Application of a Theory and
       ! Simulation-based Convective Boundary Mixing Model for AGB Star
       ! Evolution and Nucleosynthesis", ApJ, 827, 30 (2016)

       ierr = 0
       call star_ptr(id, s, ierr)
       if (ierr /= 0) return

       ! Extract parameters

       f = s%overshoot_f(j)
       f0 = s%overshoot_f0(j)
       f2 = overshoot_f2(j)

       D0 = s%overshoot_D0(j)
       Delta0 = s%overshoot_Delta0(j)
       D2 = overshoot_D2(j)

       if (f <= 0._dp .OR. f0 <= 0._dp .OR. f2 <= 0._dp) then
          write(*,*) 'ERROR: for double exponential overshooting, must set f, f0 and f2 > 0'
          write(*,*) 'see description of overshooting in star/defaults/control.defaults'
          ierr = -1
          return
       end if

       if (D0 <= 0._dp .AND. Delta0 <= 0._dp) then
          write(*,*) 'ERROR: for double exponential overshooting, must set D0 or Delta0 > 0'
          write(*,*) 'see description of overshooting in star/defaults/control.defaults'
          ierr = -1
          return
       end if

       if (D2 <= 0._dp) then
          write(*,*) 'ERROR: for double exponential overshooting, must set D2 > 0'
          write(*,*) 'see description of overshooting in star/defaults/control.defaults'
          ierr = -1
          return
       endif

       ! Apply mass limits

       if (s%mstar/msun < s%overshoot_mass_full_on(j)) then
          if (s%mstar/msun > s%overshoot_mass_full_off(j)) then
             w = (s%mstar/msun - s%overshoot_mass_full_off(j)) / &
                 (s%overshoot_mass_full_on(j) - s%overshoot_mass_full_off(j))
             factor = 0.5_dp*(1._dp - cospi(w))
             f = f*factor
             f0 = f0*factor
             f2 = f2*factor
          else
             f = 0._dp
             f0 = 0._dp
             f2 = 0._dp
          endif
       endif

       ! Evaluate convective boundary (_cb) parameters

       call star_eval_conv_bdy_Hp(s, i, Hp_cb, ierr)
       if (ierr /= 0) return

       ! Evaluate overshoot boundary (_ob) parameters

       call star_eval_over_bdy_params(s, i, f0, k_ob, r_ob, D_ob, vc_ob, ierr)
       if (ierr /= 0) return

       ! Evaluate additional parameters

       if (f2 > 0._dp .AND. D2 > 0._dp .AND. D2 < (D0 + Delta0*D_ob)) then
          dr2 = f*Hp_cb*log((D0 + Delta0*D_ob)/D2)/2._dp
       else
          dr2 = HUGE(0._dp)
       endif

       ! Loop over cell faces, adding overshoot until D <= overshoot_D_min

       outward = s%top_conv_bdy(i)

       if (outward) then
          k_a = k_ob
          k_b = 1
          dk = -1
       else
          k_a = k_ob+1
          k_b = s%nz
          dk = 1
       endif

       face_loop : do k = k_a, k_b, dk

          ! Evaluate the exponential factor

          r = s%r(k)

          if (outward) then
             dr = r - r_ob
          else
             dr = r_ob - r
          endif

          if (f > 0._dp .AND. f2 > 0._dp) then
             if (dr >= dr2) then
                factor = exp(-2._dp*(dr2/(f*Hp_cb) + (dr-dr2)/(f2*Hp_cb)))
             else
                factor = exp(-2._dp*dr/(f*Hp_cb))
             endif
          else
             factor = 0._dp
          endif

          ! Store the diffusion coefficient and velocity

          D(k) = (D0 + Delta0*D_ob)*factor
          if(D_ob /= 0d0) then
             vc(k) = (D0/D_ob + Delta0)*vc_ob*factor
          else
             vc(k) = 0d0
          end if

          ! Check for early overshoot completion

          if (D(k) < s%overshoot_D_min) then
             k_b = k
             exit face_loop
          endif
       
       end do face_loop

       ! Finish

       ierr = 0

       return

      end subroutine eval_overshoot_dbl_exp
