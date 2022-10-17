module accurate_sum_auto_diff_star_order1
   
   use const_def
   use auto_diff
   
   implicit none
   
   private
   public :: accurate_auto_diff_real_star_order1, neumaier_sum, operator(+), operator(-), assignment(=), &
      operator(<), operator(>), operator(*), operator(/)
   
   ! Type for easily using Neumaier's summation algorithm in place of normal addition.
   type accurate_auto_diff_real_star_order1
      type(auto_diff_real_star_order1) :: sum, compensator
   
   contains
      
      procedure :: value
   
   end type accurate_auto_diff_real_star_order1
   
   interface operator(+)
      procedure add_acc_adr
      procedure add_adr_acc
      procedure add_acc_acc
   end interface operator(+)
   
   interface operator(-)
      procedure sub_acc_adr
      procedure sub_adr_acc
      procedure sub_acc_acc
   end interface operator(-)
   
   interface operator(*)
      procedure mult_acc_adr
      procedure mult_adr_acc
      procedure mult_acc_acc
      procedure mult_acc_rdp
      procedure mult_rdp_acc
   end interface operator(*)
   
   interface operator(/)
      procedure div_acc_adr
      procedure div_adr_acc
      procedure div_acc_acc
      procedure div_acc_rdp
      procedure div_rdp_acc
   end interface operator(/)
   
   
   interface assignment(=)
      procedure set_acc_adr
      procedure set_adr_acc
      procedure set_acc_acc
      end interface assignment(=)
      
      interface operator(<)
      procedure acc_less_than_adr
      procedure num_less_than_acc
      end interface operator(<)
      
      interface operator(>)
   procedure acc_greater_than_adr
   procedure num_greater_than_acc
   end interface operator(>)
      
      
      contains
      
      ! Helper method for evaluating an accurate_auto_diff_real_star_order1.
      function value(this) result(res)
      ! Inputs
      class(accurate_auto_diff_real_star_order1), intent(in) :: this
   type(auto_diff_real_star_order1) :: res
   
   res = this % sum + this % compensator
   end function value
   
   ! Performs one step of Neumaier's summation algorithm
   ! (see doi:10.1002/zamm.19740540106)
   ! NOTE: Do not use --ffast-math or the like with this routine.
      ! The compiler will optimize away the trick which preserves
      ! accuracy.
      subroutine neumaier_sum(sum, compensator, summand)
      ! Inputs
      type(auto_diff_real_star_order1) sum, compensator, summand
      
      ! Intermediates
      type(auto_diff_real_star_order1) provisional
   
   provisional = sum + summand
   if (abs(sum) >= abs(summand)) then
   compensator = compensator + (sum - provisional) + summand
   else
   compensator = compensator + (summand - provisional) + sum
   end if
   sum = provisional
      
      end subroutine neumaier_sum
      
      ! The remaining are helper methods which overload +,-,*,/,=,<,>
      ! to work with accurate_auto_diff_real_star_order1 and type(auto_diff_real_star_order1) numbers interchangeably.
      
      ! *************
      type(accurate_auto_diff_real_star_order1) function mult_acc_acc(op1, op2) result (ret)
      type(accurate_auto_diff_real_star_order1), intent(in) :: op1
      type(accurate_auto_diff_real_star_order1), intent(in) :: op2
   
   ret%sum = op1%sum*op2%sum
   ret%compensator = op1%compensator * op2%sum
   
   call neumaier_sum(ret%sum, ret%compensator, op1%sum*op2%compensator)
      call neumaier_sum(ret%sum, ret%compensator, op1%compensator*op2%compensator)
      
      end function mult_acc_acc
      
      type(accurate_auto_diff_real_star_order1) function mult_acc_adr(op1, op2) result (ret)
      type(accurate_auto_diff_real_star_order1), intent(in) :: op1
      type(auto_diff_real_star_order1), intent(in) :: op2
   
   ret%sum = op1%sum * op2
   ret%compensator = op1%compensator * op2
   
   end function mult_acc_adr
   
   type(accurate_auto_diff_real_star_order1) function mult_adr_acc(op1, op2) result (ret)
      type(auto_diff_real_star_order1), intent(in) :: op1
   type(accurate_auto_diff_real_star_order1), intent(in) :: op2
   
   ret%sum = op2%sum * op1
   ret%compensator = op2%compensator * op1
   
   end function mult_adr_acc
   
   type(accurate_auto_diff_real_star_order1) function mult_acc_rdp(op1, op2) result (ret)
      type(accurate_auto_diff_real_star_order1), intent(in) :: op1
      real(dp), intent(in) :: op2
   
   ret%sum = op1%sum * op1
   ret%compensator = op1%compensator * op1
   
   end function mult_acc_rdp
   
   type(accurate_auto_diff_real_star_order1) function mult_rdp_acc(op1, op2) result (ret)
      real(dp), intent(in) :: op1
      type(accurate_auto_diff_real_star_order1), intent(in) :: op2
      
      ret%sum = op2%sum * op1
      ret%compensator = op2%compensator * op1
      
      end function mult_rdp_acc
      
      ! ///////////
      type(accurate_auto_diff_real_star_order1) function div_acc_acc(op1, op2) result (ret)
         type(accurate_auto_diff_real_star_order1), intent(in) :: op1
   type(accurate_auto_diff_real_star_order1), intent(in) :: op2
   
   ret % sum = op1 % value() / op2 % value()
      ret % compensator = 0
      
      end function div_acc_acc
      
      type(accurate_auto_diff_real_star_order1) function div_acc_adr(op1, op2) result (ret)
      type(accurate_auto_diff_real_star_order1), intent(in) :: op1
      type(auto_diff_real_star_order1), intent(in) :: op2
   
   ret%sum = op1%sum / op2
   ret%compensator = op1%compensator / op2
   
   end function div_acc_adr
   
   type(accurate_auto_diff_real_star_order1) function div_adr_acc(op1, op2) result (ret)
   type(auto_diff_real_star_order1), intent(in) :: op1
   type(accurate_auto_diff_real_star_order1), intent(in) :: op2
   
   ret%sum = op1 / op2%value()
   ret%compensator = 0
   
   end function div_adr_acc
   
   type(accurate_auto_diff_real_star_order1) function div_acc_rdp(op1, op2) result (ret)
      type(accurate_auto_diff_real_star_order1), intent(in) :: op1
   real(dp), intent(in) :: op2
   
   ret%sum = op1%sum / op2
   ret%compensator = op1%compensator / op2
   
   end function div_acc_rdp
   
   type(accurate_auto_diff_real_star_order1) function div_rdp_acc(op1, op2) result (ret)
      real(dp), intent(in) :: op1
      type(accurate_auto_diff_real_star_order1), intent(in) :: op2
   
   ret%sum = op1 / op2%value()
      ret%compensator = 0
      
      end function div_rdp_acc
      
      ! ============
      subroutine set_acc_adr(this, new)
      type(accurate_auto_diff_real_star_order1), intent(out) :: this
   type(auto_diff_real_star_order1), intent(in) :: new
   this % sum = new
   this % compensator = 0
   end subroutine set_acc_adr
   
   subroutine set_adr_acc(this, new)
      type(auto_diff_real_star_order1), intent(out) :: this
   type(accurate_auto_diff_real_star_order1), intent(in) :: new
   this = new % value()
      end subroutine set_adr_acc
      
      subroutine set_acc_acc(this, new)
      type(accurate_auto_diff_real_star_order1), intent(out) :: this
      type(accurate_auto_diff_real_star_order1), intent(in) :: new
      this % sum = new % sum
         this % compensator = new % compensator
         end subroutine set_acc_acc
   
   ! <<<<<<<<<<<
   logical function acc_less_than_adr(acc, num) result (ret)
      type(accurate_auto_diff_real_star_order1), intent(in) :: acc
      type(auto_diff_real_star_order1), intent(in) :: num
   if (acc % value() < num) then
   ret = .true.
   else
   ret = .false.
   end if
   end function acc_less_than_adr
   
   
   logical function num_less_than_acc(num, acc) result (ret)
      type(accurate_auto_diff_real_star_order1), intent(in) :: acc
      type(auto_diff_real_star_order1), intent(in) :: num
      if (num < acc % value()) then
   ret = .true.
   else
   ret = .false.
   end if
   end function num_less_than_acc
   
   ! >>>>>>>>>>>
      logical function acc_greater_than_adr(acc, num) result (ret)
      type(accurate_auto_diff_real_star_order1), intent(in) :: acc
   type(auto_diff_real_star_order1), intent(in) :: num
   if (acc % value() > num) then
   ret = .true.
   else
   ret = .false.
   end if
   end function acc_greater_than_adr
   
   logical function num_greater_than_acc(num, acc) result (ret)
      type(accurate_auto_diff_real_star_order1), intent(in) :: acc
      type(auto_diff_real_star_order1), intent(in) :: num
   if (num > acc % value()) then
   ret = .true.
   else
   ret = .false.
   end if
   end function num_greater_than_acc
   
   ! +++++++++++
      type(accurate_auto_diff_real_star_order1) function add_acc_acc(op1, op2) result (ret)
      type(accurate_auto_diff_real_star_order1), intent(in) :: op1
      type(accurate_auto_diff_real_star_order1), intent(in) :: op2
      
      ret%sum = op1%sum
         ret%compensator = op1%compensator
         call neumaier_sum(ret%sum, ret%compensator, op2%sum)
         call neumaier_sum(ret%sum, ret%compensator, op2%compensator)
         
         end function add_acc_acc
         
         type(accurate_auto_diff_real_star_order1) function add_acc_adr(op1, op2) result (ret)
      type(accurate_auto_diff_real_star_order1), intent(in) :: op1
      type(auto_diff_real_star_order1), intent(in) :: op2
      
      ret%sum = op1%sum
         ret%compensator = op1%compensator
         call neumaier_sum(ret%sum, ret%compensator, op2)
         
         end function add_acc_adr
         
         type(accurate_auto_diff_real_star_order1) function add_adr_acc(op1, op2) result (ret)
         type(auto_diff_real_star_order1), intent(in) :: op1
      type(accurate_auto_diff_real_star_order1), intent(in) :: op2
      
      ret%sum = op2%sum
         ret%compensator = op2%compensator
         call neumaier_sum(ret%sum, ret%compensator, op1)
      
      end function add_adr_acc
      
      ! -----------
      type(accurate_auto_diff_real_star_order1) function sub_acc_acc(op1, op2) result (ret)
         type(accurate_auto_diff_real_star_order1), intent(in) :: op1
   type(accurate_auto_diff_real_star_order1), intent(in) :: op2
   
   ret%sum = op1%sum
      ret%compensator = op1%compensator
      call neumaier_sum(ret%sum, ret%compensator, -op2%sum)
      call neumaier_sum(ret%sum, ret%compensator, -op2%compensator)
   
   end function sub_acc_acc
   
   type(accurate_auto_diff_real_star_order1) function sub_acc_adr(op1, op2) result (ret)
      type(accurate_auto_diff_real_star_order1), intent(in) :: op1
   type(auto_diff_real_star_order1), intent(in) :: op2
   
   ret%sum = op1%sum
      ret%compensator = op1%compensator
      call neumaier_sum(ret%sum, ret%compensator, -op2)
      
      end function sub_acc_adr
      
      type(accurate_auto_diff_real_star_order1) function sub_adr_acc(op1, op2) result (ret)
      type(auto_diff_real_star_order1), intent(in) :: op1
      type(accurate_auto_diff_real_star_order1), intent(in) :: op2
   
   ret%sum = op1
   ret%compensator = -op2%compensator
   call neumaier_sum(ret%sum, ret%compensator, -op2%sum)

end function sub_adr_acc

end module accurate_sum_auto_diff_star_order1
