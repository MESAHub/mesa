subroutine  precompute_some_zs()
    integer :: int_i
    real(dp) :: i
 
      do int_i=1,max_precomp_ints
         i = real(int_i,kind=dp)

         pre_z(int_i)% z2 = i*i
         pre_z(int_i)% z3 = pre_z(int_i)% z2 * i
         pre_z(int_i)% z4 = pre_z(int_i)% z3 * i
         pre_z(int_i)% z5 = pre_z(int_i)% z4 * i
         pre_z(int_i)% z6 = pre_z(int_i)% z5 * i
         pre_z(int_i)% z7 = pre_z(int_i)% z6 * i
         pre_z(int_i)% z8 = pre_z(int_i)% z7 * i

         pre_z(int_i)% z1_3 = pow(i,1d0/3d0)
         pre_z(int_i)% z2_3 = pow(i,2d0/3d0)
         pre_z(int_i)% z4_3 = pow(i,4d0/3d0)
         pre_z(int_i)% z5_3 = pow(i,5d0/3d0)
         pre_z(int_i)% z7_3 = pow(i,7d0/3d0)

         pre_z(int_i)% zm1_3 = pow(i,-1d0/3d0)
         pre_z(int_i)% zm2_3 = pow(i,-2d0/3d0)
         pre_z(int_i)% zm4_3 = pow(i,-4d0/3d0)
         pre_z(int_i)% zm5_3 = pow(i,-5d0/3d0)
         pre_z(int_i)% zm7_3 = pow(i,-7d0/3d0)

         pre_z(int_i)% z7_6 = pow(i,7d0/6d0)

         pre_z(int_i)% z1_5 = pow(i,1d0/5d0)
         pre_z(int_i)% z2_5 = pow(i,2d0/5d0)
         pre_z(int_i)% z3_5 = pow(i,3d0/5d0)
         pre_z(int_i)% z4_5 = pow(i,4d0/5d0)

         pre_z(int_i)% logz = log(i)
         pre_z(int_i)% sqlogz = pow2(pre_z(int_i)% logz)
         pre_z(int_i)% logz_3_2 = pow(pre_z(int_i)% logz, 3d0/2d0)

         pre_z(int_i)% z1_2 = pow(i,0.5d0)
         pre_z(int_i)% z3_2 = pow(i,1.5d0)

         pre_z(int_i)% zm1_2 = pow(i,-0.5d0)
         pre_z(int_i)% zm1_4 = pow(i,-0.25d0)
         pre_z(int_i)% zm3_2 = pow(i,-3d0/2d0)

         pre_z(int_i)% z0p475 = pow(i,0.475d0)
         pre_z(int_i)% zp1_3_2 = pow(1d0 + i, 3d0/2d0)
         pre_z(int_i)% zm0p267 = pow(i, -0.267d0)

      end do
 
 end subroutine  precompute_some_zs
