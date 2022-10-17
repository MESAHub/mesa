module accurate_sum
   
   use const_def
   
   implicit none
   
   private
   public :: accurate_real, neumaier_sum, operator(+), operator(-), assignment(=), &
      operator(<), operator(>), operator(*), operator(/)
   
   ! Type for easily using Neumaier's summation algorithm in place of normal addition.
   type accurate_real
      real(qp) :: sum
      real(qp) :: compensator
   
   contains
      
      procedure :: value
   
   end type accurate_real
   
   interface operator(+)
      procedure add_accurate_real
      procedure add_accurate_real_num
      procedure add_accurate_real_num_rev
   end interface operator(+)
   
   interface operator(-)
      procedure sub_accurate_real
      procedure sub_accurate_real_num
      procedure sub_accurate_real_num_rev
   end interface operator(-)
   
   interface operator(*)
      procedure mult_acc_real
      procedure mult_real_acc
      procedure mult_acc_acc
   end interface operator(*)
   
   interface operator
      (/)
   procedure div_acc_real
   procedure div_real_acc
   procedure div_acc_acc
end interface operator(/)
   
   
   interface assignment(=)
      procedure set_accurate_real_num
      procedure set_num_accurate_real
      procedure set_acc_acc
      end interface assignment(=)
      
      interface operator(<)
      procedure acc_less_than_num
      procedure num_less_than_acc
      end interface operator(<)
      
      interface operator(>)
   procedure acc_greater_than_num
   procedure num_greater_than_acc
   end interface operator(>)
      
      
      contains
      
      ! Helper method for evaluating an accurate_real.
      function value(this) result(res)
      ! Inputs
      class(accurate_real), intent(in) :: this
   real(qp) res
   
   res = this % sum + this % compensator
   end function value
   
   ! Performs one step of Neumaier's summation algorithm
   ! (see doi:10.1002/zamm.19740540106)
   ! NOTE: Do not use --ffast-math or the like with this routine.
   ! The compiler will optimize away the trick which preserves
      ! accuracy.
      subroutine neumaier_sum(sum, compensator, summand)
      ! Inputs
      real(qp) sum, compensator, summand
      
      ! Intermediates
      real(qp) provisional
   
   provisional = sum + summand
   if (abs(sum) >= abs(summand)) then
   compensator = compensator + (sum - provisional) + summand
   else
   compensator = compensator + (summand - provisional) + sum
   end if
   sum = provisional
      
      end subroutine neumaier_sum
      
      ! The remaining are helper methods which overload +,-,*,/,=,<,>
      ! to work with accurate_real and real(qp) numbers interchangeably.
      
      ! *************
      type(accurate_real) function mult_acc_acc(op1, op2) result (ret)
      type(accurate_real), intent(in) :: op1
      type(accurate_real), intent(in) :: op2
   
   ret%sum = op1%sum*op2%sum
   ret%compensator = op1%compensator * op2%sum
   
   call neumaier_sum(ret%sum, ret%compensator, op1%sum*op2%compensator)
      call neumaier_sum(ret%sum, ret%compensator, op1%compensator*op2%compensator)
      
      end function mult_acc_acc
      
      type(accurate_real) function mult_acc_real(op1, op2) result (ret)
      type(accurate_real), intent(in) :: op1
      real(qp), intent(in) :: op2
   
   ret%sum = op1%sum * op2
   ret%compensator = op1%compensator * op2
   
   end function mult_acc_real
   
   type(accurate_real) function mult_real_acc(op1, op2) result (ret)
      real(qp), intent(in) :: op1
   type(accurate_real), intent(in) :: op2
   
   ret%sum = op2%sum * op1
   ret%compensator = op2%compensator * op1
   
   end function mult_real_acc
   
   ! ///////////
   type(accurate_real) function div_acc_acc(op1, op2) result (ret)
      type(accurate_real), intent(in) :: op1
      type(accurate_real), intent(in) :: op2
   
   ret % sum = op1 % value() / op2 % value()
      ret % compensator = 0
      
      end function div_acc_acc
      
      type(accurate_real) function div_acc_real(op1, op2) result (ret)
      type(accurate_real), intent(in) :: op1
      real(qp), intent(in) :: op2
   
   ret%sum = op1%sum / op2
   ret%compensator = op1%compensator / op2
   
   end function div_acc_real
   
   type(accurate_real) function div_real_acc(op1, op2) result (ret)
      real(qp), intent(in) :: op1
   type(accurate_real), intent(in) :: op2
   
   ret%sum = op1 / op2%value()
      ret%compensator = 0
      
      end function div_real_acc
      
      ! ============
      subroutine set_accurate_real_num(this, new)
         type(accurate_real), intent(out) :: this
      real(qp), intent(in) :: new
      this % sum = new
      this % compensator = 0
      end subroutine set_accurate_real_num
      
      subroutine set_num_accurate_real(this, new)
         real(qp), intent(out) :: this
         type(accurate_real), intent(in) :: new
      this = new % value()
         end subroutine set_num_accurate_real
         
         subroutine set_acc_acc(this, new)
         type(accurate_real), intent(out) :: this
      type(accurate_real), intent(in) :: new
      this % sum = new % sum
         this % compensator = new % compensator
         end subroutine set_acc_acc
         
         ! <<<<<<<<<<<
         logical function acc_less_than_num(acc, num) result (ret)
      type(accurate_real), intent(in) :: acc
      real(qp), intent(in) :: num
      if (acc % value() < num) then
      ret = .true.
      else
      ret = .false.
      end if
      end function acc_less_than_num
      
      
      logical function num_less_than_acc(num, acc) result (ret)
         type(accurate_real), intent(in) :: acc
   real(qp), intent(in) :: num
   if (num < acc % value()) then
   ret = .true.
   else
   ret = .false.
   end if
   end function num_less_than_acc
   
   ! >>>>>>>>>>>
      logical function acc_greater_than_num(acc, num) result (ret)
      type(accurate_real), intent(in) :: acc
      real(qp), intent(in) :: num
   if (acc % value() > num) then
   ret = .true.
   else
   ret = .false.
   end if
   end function acc_greater_than_num
   
   logical function num_greater_than_acc(num, acc) result (ret)
      type(accurate_real), intent(in) :: acc
   real(qp), intent(in) :: num
   if (num > acc % value()) then
   ret = .true.
   else
   ret = .false.
   end if
   end function num_greater_than_acc
   
   ! +++++++++++
      type(accurate_real) function add_accurate_real(op1, op2) result (ret)
         type(accurate_real), intent(in) :: op1
      type(accurate_real), intent(in) :: op2
      
      ret%sum = op1%sum
         ret%compensator = op1%compensator
         call neumaier_sum(ret%sum, ret%compensator, op2%sum)
         call neumaier_sum(ret%sum, ret%compensator, op2%compensator)
      
      end function add_accurate_real
      
      type(accurate_real) function add_accurate_real_num(op1, op2) result (ret)
         type(accurate_real), intent(in) :: op1
      real(qp), intent(in) :: op2
      
      ret%sum = op1%sum
         ret%compensator = op1%compensator
         call neumaier_sum(ret%sum, ret%compensator, op2)
         
         end function add_accurate_real_num
         
         type(accurate_real) function add_accurate_real_num_rev(op1, op2) result (ret)
         real(qp), intent(in) :: op1
         type(accurate_real), intent(in) :: op2
      
      ret%sum = op2%sum
         ret%compensator = op2%compensator
         call neumaier_sum(ret%sum, ret%compensator, op1)
         
         end function add_accurate_real_num_rev
         
         ! -----------
         type(accurate_real) function sub_accurate_real(op1, op2) result (ret)
      type(accurate_real), intent(in) :: op1
      type(accurate_real), intent(in) :: op2
   
   ret%sum = op1%sum
      ret%compensator = op1%compensator
      call neumaier_sum(ret%sum, ret%compensator, -op2%sum)
      call neumaier_sum(ret%sum, ret%compensator, -op2%compensator)
      
      end function sub_accurate_real
      
      type(accurate_real) function sub_accurate_real_num(op1, op2) result (ret)
      type(accurate_real), intent(in) :: op1
      real(qp), intent(in) :: op2
   
   ret%sum = op1%sum
      ret%compensator = op1%compensator
      call neumaier_sum(ret%sum, ret%compensator, -op2)
      
      end function sub_accurate_real_num
      
      type(accurate_real) function sub_accurate_real_num_rev(op1, op2) result (ret)
   real(qp), intent(in) :: op1
   type(accurate_real), intent(in) :: op2
   
   ret%sum = op1
   ret%compensator = -op2%compensator
   call neumaier_sum(ret%sum, ret%compensator, -op2%sum)

end function sub_accurate_real_num_rev

end module accurate_sum
