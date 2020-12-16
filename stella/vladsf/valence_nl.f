c This routine determines the principal quantum number and angular momentum
c quantum number of an element of nuclear charge z and in ionization i,
c assuming that shells are filled sequentially.
      subroutine valence_nl(z, i, n_princ, l_ang)
      implicit integer (a-z)

      integer z, i, n_princ, l_ang

c Number of electrons:
      ne = z - i + 1
c Assume that shells are filled in order of increasing
c principal quantum number and angular momentum state.
c This is not always the case however...
c Determine n,l of the outermost partially filled shell.

      do n = 1, 99
c This is the number electrons up to and including a full
c n shell.
         ne_full_shell = (n * (n + 1) * (2 * n + 1)) / 3

         if (ne .le. ne_full_shell) then
            do l = 0, n - 1
c This is the number of electrons which are in all lower filled
c subshells.
               ne_lim = ((n - 1) * n * (2 * n - 1)) / 3
     .              + 2 * l**2
c max_shell is the maximum number of electrons which would fit in this shell.
               max_shell = 2 * (2 * l + 1)
c excess is the number of electrons available for filling this shell.
               excess = ne - ne_lim

               if (excess .le. max_shell) then
c If the number of available electrons is less than or equal to what
c will fit in this shell, then it is the valence shell.
c nholes is the number of holes.
                  n_princ = n
                  l_ang = l
                  return
               end if
            end do
         end if
      end do
      
      write(0,*) " couldn't find n,l for z = ", z,
     .     " i = ", i

      stop

      end
