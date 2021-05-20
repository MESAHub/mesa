! ***********************************************************************
!
!   Copyright (C) 2012  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module mod_random
      use const_def, only: dp
      

      contains


      subroutine get_seed ( seed )

      !*****************************************************************************80
      !
      !! GET_SEED returns a seed for the random number generator.
      !
      !  Discussion:
      !
      !    The seed depends on the current time, and ought to be (slightly)
      !    different every millisecond.  Once the seed is obtained, a random
      !    number generator should be called a few times to further process
      !    the seed.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    02 August 2004
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Output, integer ( kind = 4 ) SEED, a pseudorandom seed value.
      !
        implicit none

        integer ( kind = 4 ) seed
        real ( kind = 8 ) temp
        character ( len = 10 ) time
        character ( len = 8 ) today
        integer ( kind = 4 ) values(8)
        character ( len = 5 ) zone

        call date_and_time ( today, time, zone, values )

        temp = 0.0D+00

        temp = temp + real ( values(2) - 1, kind = 8 ) /  11.0D+00
        temp = temp + real ( values(3) - 1, kind = 8 ) /  30.0D+00
        temp = temp + real ( values(5),     kind = 8 ) /  23.0D+00
        temp = temp + real ( values(6),     kind = 8 ) /  59.0D+00
        temp = temp + real ( values(7),     kind = 8 ) /  59.0D+00
        temp = temp + real ( values(8),     kind = 8 ) / 999.0D+00
        temp = temp                                    /   6.0D+00

        do while ( temp <= 0.0D+00 )
          temp = temp + 1.0D+00
        end do

        do while ( 1.0D+00 < temp )
          temp = temp - 1.0D+00
        end do

        seed = int ( real ( huge ( 1 ), kind = 8 ) * temp )
      !
      !  Never use a seed of 0 or maximum integer ( kind = 4 ).
      !
        if ( seed == 0 ) then
          seed = 1
        end if

        if ( seed == huge ( 1 ) ) then
          seed = seed - 1
        end if

        return
      end subroutine get_seed


      function i4_uniform ( a, b, seed )

      !*****************************************************************************80
      !
      !! I4_UNIFORM returns a scaled pseudorandom I4.
      !
      !  Discussion:
      !
      !    An I4 is an integer ( kind = 4 ) value.
      !
      !    The pseudorandom number will be scaled to be uniformly distributed
      !    between A and B.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    31 May 2007
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Reference:
      !
      !    Paul Bratley, Bennett Fox, Linus Schrage,
      !    A Guide to Simulation,
      !    Second Edition,
      !    Springer, 1987,
      !    ISBN: 0387964673,
      !    LC: QA76.9.C65.B73.
      !
      !    Bennett Fox,
      !    Algorithm 647:
      !    Implementation and Relative Efficiency of Quasirandom
      !    Sequence Generators,
      !    ACM Transactions on Mathematical Software,
      !    Volume 12, Number 4, December 1986, pages 362-376.
      !
      !    Pierre L'Ecuyer,
      !    Random Number Generation,
      !    in Handbook of Simulation,
      !    edited by Jerry Banks,
      !    Wiley, 1998,
      !    ISBN: 0471134031,
      !    LC: T57.62.H37.
      !
      !    Peter Lewis, Allen Goodman, James Miller,
      !    A Pseudo-Random Number Generator for the System/360,
      !    IBM Systems Journal,
      !    Volume 8, Number 2, 1969, pages 136-143.
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) A, B, the limits of the interval.
      !
      !    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
      !    should NOT be 0.  On output, SEED has been updated.
      !
      !    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
      !
        implicit none

        integer ( kind = 4 ) a
        integer ( kind = 4 ) b
        integer ( kind = 4 ), parameter :: i4_huge = 2147483647
        integer ( kind = 4 ) i4_uniform
        integer ( kind = 4 ) k
        real ( kind = 4 ) r
        integer ( kind = 4 ) seed
        integer ( kind = 4 ) value

        if ( seed == 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
          write ( *, '(a)' ) '  Input value of SEED = 0.'
          stop
        end if

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed < 0 ) then
          seed = seed + i4_huge
        end if

        r = real ( seed, kind = 4 ) * 4.656612875E-10
      !
      !  Scale R to lie between A-0.5 and B+0.5.
      !
        r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) &
          +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
      !
      !  Use rounding to convert R to an integer between A and B.
      !
        value = nint ( r, kind = 4 )

        value = max ( value, min ( a, b ) )
        value = min ( value, max ( a, b ) )

        i4_uniform = value

        return
      end function i4_uniform


      subroutine perm_uniform ( n, base, seed, p )

      !*****************************************************************************80
      !
      !! PERM_UNIFORM selects a random permutation of N objects.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    18 November 2008
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Reference:
      !
      !    Albert Nijenhuis, Herbert Wilf,
      !    Combinatorial Algorithms,
      !    Academic Press, 1978, second edition,
      !    ISBN 0-12-519260-6.
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
      !
      !    Input, integer ( kind = 4 ) BASE, is 0 for a 0-based permutation and 1 for
      !    a 1-based permutation.
      !
      !    Input/output, integer ( kind = 4 ) SEED, a seed for the random
      !    number generator.
      !
      !    Output, integer ( kind = 4 ) P(N), the permutation.  P(I) is the "new"
      !    location of the object originally at I.
      !
        implicit none

        integer ( kind = 4 ) n

        integer ( kind = 4 ) base
        integer ( kind = 4 ) i
        integer ( kind = 4 ) j
        integer ( kind = 4 ) k
        integer ( kind = 4 ) p(n)
        integer ( kind = 4 ) seed

        do i = 1, n
          p(i) = ( i - 1 ) + base
        end do

        do i = 1, n
          j = i4_uniform ( i, n, seed )
          k    = p(i)
          p(i) = p(j)
          p(j) = k
        end do

        return
      end subroutine perm_uniform



      function r8_uniform_01 ( seed )

      !*****************************************************************************80
      !
      !! R8_UNIFORM_01 returns a unit pseudorandom R8.
      !
      !  Discussion:
      !
      !    An R8 is a real ( kind = 8 ) value.
      !
      !    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
      !
      !    This routine implements the recursion
      !
      !      seed = 16807 * seed mod ( 2**31 - 1 )
      !      r8_uniform_01 = seed / ( 2**31 - 1 )
      !
      !    The integer ( kind = 4 ) arithmetic never requires more than 32 bits,
      !    including a sign bit.
      !
      !    If the initial seed is 12345, then the first three computations are
      !
      !      Input     Output      R8_UNIFORM_01
      !      SEED      SEED
      !
      !         12345   207482415  0.096616
      !     207482415  1790989824  0.833995
      !    1790989824  2035175616  0.947702
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    05 July 2006
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Reference:
      !
      !    Paul Bratley, Bennett Fox, Linus Schrage,
      !    A Guide to Simulation,
      !    Springer Verlag, pages 201-202, 1983.
      !
      !    Pierre L'Ecuyer,
      !    Random Number Generation,
      !    in Handbook of Simulation,
      !    edited by Jerry Banks,
      !    Wiley Interscience, page 95, 1998.
      !
      !    Bennett Fox,
      !    Algorithm 647:
      !    Implementation and Relative Efficiency of Quasirandom
      !    Sequence Generators,
      !    ACM Transactions on Mathematical Software,
      !    Volume 12, Number 4, pages 362-376, 1986.
      !
      !    Peter Lewis, Allen Goodman, James Miller
      !    A Pseudo-Random Number Generator for the System/360,
      !    IBM Systems Journal,
      !    Volume 8, pages 136-143, 1969.
      !
      !  Parameters:
      !
      !    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
      !    NOT be 0. On output, SEED has been updated.
      !
      !    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
      !    strictly between 0 and 1.
      !
        implicit none

        integer ( kind = 4 ) k
        real ( kind = 8 ) r8_uniform_01
        integer ( kind = 4 ) seed

        if ( seed == 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
          write ( *, '(a)' ) '  Input value of SEED = 0.'
          stop
        end if

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed < 0 ) then
          seed = seed + 2147483647
        end if
      !
      !  Although SEED can be represented exactly as a 32 bit integer ( kind = 4 ),
      !  it generally cannot be represented exactly as a 32 bit real number!
      !
        r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

        return
      end function r8_uniform_01


      end module mod_random



























