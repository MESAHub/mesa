      
      ! use the Tolman–Oppenheimer–Volkoff (TOV) equation. 
      ! See first equation in https://en.wikipedia.org/wiki/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation.  
      ! want to replace -G*m/r^2 by -G*m/r^2*(1 + P/(rho c^2))(1 + 4 pi r^3 P /(m c^2))/(1 - 2 G m/(r c^2))
      subroutine my_other_cgrav(id, ierr)
         use star_def
         use utils_lib, only: is_bad
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type(star_info), pointer :: s
         integer :: k
         real(dp) :: r, G, m, rho, P, f1, f2, f3
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         do k=1,s% nz
            r = s% r(k)
            G = standard_cgrav
            m = s% m(k) ! m_grav hasn't been set when other_cgrav is called
            if (k > 1) then ! get approximate P and rho at face
               rho = 0.5d0*(s% rho(k-1) + s% rho(k))
               P = 0.5d0*(s% P(k-1) + s% P(k))
               f1 = 1d0 + P/(rho*clight**2)
               f2 = 1d0 + pi4*pow3(r)*P/(m*clight**2)
            else ! k == 1
               f1 = 1d0
               f2 = 1d0
            end if
            f3 = 1d0 - 2d0*G*m/(r*clight**2)
            s% cgrav(k) = G*f1*f2/f3
            if (is_bad(s% cgrav(k))) then
               write(*,2) 's% cgrav(k)', k, s% cgrav(k)
               write(*,2) 'f1', k, f1
               write(*,2) 'f2', k, f2
               write(*,2) 'f3', k, f3
               write(*,2) 'r', k, r
               write(*,2) 'G', k, G
               write(*,2) 'm', k, m
               if (k > 1) then
                  write(*,2) 'P', k, P
                  write(*,2) 'rho', k, rho
               end if
               stop 'my_other_cgrav'
            end if
         end do
         !write(*,*) 'done my_other_cgrav', s% model_number
      end subroutine my_other_cgrav
