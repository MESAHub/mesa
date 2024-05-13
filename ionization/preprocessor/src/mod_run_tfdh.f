!
! run_tfdh
! francois hebert, dec 22 2010
!
! this is the function you call!
! turns ne, T, composition ==> znet
!
! optionally passes back profile information
!

      module mod_run_tfdh

      use def_args
      use def_type
      use def_const
      use lib_util
      use mod_tfdh
      use root_chi
      use root_vex

      implicit none

      contains


      subroutine run_tfdh(ztr, tmp, ne, xs, zs, as, profile_flag, args)

         real (dp), intent(in) :: ztr, tmp, ne
         real (dp), intent(in), dimension(:) :: xs, zs, as
         logical, intent(in) :: profile_flag
         type(datastruct), intent(inout) :: args

         real (dp) :: kt, mu, tau, chi, rws
         integer :: ns


         !
         ! first:
         ! read and process input to prepare for calling tfdh codes
         !
         kt = kb * tmp
         ns = size(xs)

         ! check consistency of input
         if (sum(xs) /= 1.0_dp) call alert(1, '(run_tfdh) mass fractions xs dont add to 1')
         if (ns /= size(as))  call alert(1, '(run_tfdh) composition vectors xs,as of different length')
         if (ns /= size(zs))  call alert(1, '(run_tfdh) composition vectors xs,zs of different length')

         ! find mu_electron, tau, chi
         mu = 1.0_dp/sum( zs(:)*xs(:)/as(:) )
         tau = kt / (me * clight**2)
         chi = rootfind_chi(ne, tau)
         rws = (3*ztr/(4*pi*ne))**(1.0_dp/3.0_dp)

         ! allocate things
         args % mu    = mu
         args % neinf = ne
         args % kt    = kt
         args % ztr   = ztr
         args % tau   = tau
         args % chi   = chi
         args % rws   = rws

         args % species = ns
         call reallocate_dp(args%as, ns)
         call reallocate_dp(args%zs, ns)
         call reallocate_dp(args%xs, ns)
         call reallocate_dp(args%niinf, ns)
         call reallocate_dp(args%qiloc, ns)
         call reallocate_dp(args%rex,   ns)

         args%as(:) = as(:)
         args%zs(:) = zs(:)
         args%xs(:) = xs(:)
         args%niinf(:) = mu*ne * xs(:)/as(:)
         args%qiloc(:) = 0
         args%rex(:) = 0

         !
         ! second:
         ! use rootfind_vex to do a rootfind for the value of dv0
         ! (this calls the fast_tfdh module)
         !
         call rootfind_vex(args)

         !
         ! finally:
         ! call the primary integration code!
         !
         call tfdh(args, profile_flag)


      end subroutine run_tfdh


      end module mod_run_tfdh



   
   