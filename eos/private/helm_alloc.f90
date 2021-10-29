! ***********************************************************************
!
!   Copyright (C) 2006, 2007-2019  Frank Timmes & The MESA Team
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module helm_alloc
      use const_def, only: dp, use_mesa_temp_cache
      use math_lib
      
      implicit none


      contains
            
      
      subroutine alloc_helm_table(h, imax, jmax, ierr)
         ! This routine allocates a Helm_Table and places pointer to it in h.
         ! It also allocates the arrays in the Helm_Table record.
         
         use eos_def
         
         type (Helm_Table), pointer :: h
         integer, intent(in) :: imax, jmax
         integer, intent(out) :: ierr ! 0 means AOK.
         
         ierr = 0
         
         allocate(h,stat=ierr)
         if (ierr /= 0) return
         
         h% imax = imax
         h% jmax = jmax 
         h% with_coulomb_corrections = .true.
         
         call alloc_1d_array(h% d, imax)
         call alloc_1d_array(h% t, jmax)
         
         !..for the helmholtz free energy tables
         call alloc_2d_array(h% f, imax, jmax)
         call alloc_2d_array(h% fd, imax, jmax)
         call alloc_2d_array(h% ft, imax, jmax)
         call alloc_2d_array(h% fdd, imax, jmax)
         call alloc_2d_array(h% ftt, imax, jmax)
         call alloc_2d_array(h% fdt, imax, jmax)
         call alloc_2d_array(h% fddt, imax, jmax)
         call alloc_2d_array(h% fdtt, imax, jmax)
         call alloc_2d_array(h% fddtt, imax, jmax)
         
         !..for the pressure derivative with density tables
         call alloc_2d_array(h% dpdf, imax, jmax)
         call alloc_2d_array(h% dpdfd, imax, jmax)
         call alloc_2d_array(h% dpdft, imax, jmax)
         call alloc_2d_array(h% dpdfdt, imax, jmax)

         !..for chemical potential tables
         call alloc_2d_array(h% ef, imax, jmax)
         call alloc_2d_array(h% efd, imax, jmax)
         call alloc_2d_array(h% eft, imax, jmax)
         call alloc_2d_array(h% efdt, imax, jmax)

         !..for the number density tables
         call alloc_2d_array(h% xf, imax, jmax)
         call alloc_2d_array(h% xfd, imax, jmax)
         call alloc_2d_array(h% xft, imax, jmax)
         call alloc_2d_array(h% xfdt, imax, jmax)
            
         !..for storing the differences
         call alloc_1d_array(h% dt_sav, jmax)
         call alloc_1d_array(h% dt2_sav, jmax)
         call alloc_1d_array(h% dti_sav, jmax)
         call alloc_1d_array(h% dt2i_sav, jmax)
         call alloc_1d_array(h% dt3i_sav, jmax)
         
         call alloc_1d_array(h% dd_sav, imax)
         call alloc_1d_array(h% dd2_sav, imax)
         call alloc_1d_array(h% ddi_sav, imax)
         call alloc_1d_array(h% dd2i_sav, imax)
         call alloc_1d_array(h% dd3i_sav, imax)

         contains
         
         subroutine alloc_1d_array(ptr,sz)
            real(dp), dimension(:), pointer :: ptr
            integer, intent(in) :: sz        
            allocate(ptr(sz),stat=ierr)
         end subroutine alloc_1d_array
         
         subroutine alloc_2d_array(ptr,sz1,sz2)
            real(dp), dimension(:,:), pointer :: ptr
            integer, intent(in) :: sz1,sz2         
            allocate(ptr(sz1,sz2),stat=ierr)
         end subroutine alloc_2d_array
      
      
      end subroutine alloc_helm_table
      
      
      subroutine setup_td_deltas(h, imax, jmax)
         use eos_def
         type (Helm_Table), pointer :: h
         integer, intent(in) :: imax, jmax
         integer :: i, j
         real(dp) dth,dt2,dti,dt2i,dt3i,dd,dd2,ddi,dd2i,dd3i
         !..construct the temperature and density deltas and their inverses 
         do j=1,jmax-1
            dth         = h% t(j+1) - h% t(j)
            dt2         = dth * dth
            dti         = 1.0d0/dth
            dt2i        = 1.0d0/dt2
            dt3i        = dt2i*dti
            h% dt_sav(j)   = dth
            h% dt2_sav(j)  = dt2
            h% dti_sav(j)  = dti
            h% dt2i_sav(j) = dt2i
            h% dt3i_sav(j) = dt3i
         end do
         do i=1,imax-1
            dd          = h% d(i+1) - h% d(i)
            dd2         = dd * dd
            ddi         = 1.0d0/dd
            dd2i        = 1.0d0/dd2
            dd3i        = dd2i*ddi
            h% dd_sav(i)   = dd
            h% dd2_sav(i)  = dd2
            h% ddi_sav(i)  = ddi
            h% dd2i_sav(i) = dd2i
            h% dd3i_sav(i) = dd3i
         enddo
      end subroutine setup_td_deltas


      subroutine read_helm_table(h, data_dir, cache_dir, temp_cache_dir, use_cache, ierr)
      use eos_def
      use utils_lib, only: mv, switch_str
      
      implicit none

      type (Helm_Table), pointer :: h
      character(*), intent(IN) :: data_dir, cache_dir, temp_cache_dir
      logical, intent(IN) :: use_cache
      integer, intent(out) :: ierr

!..this routine reads the helmholtz eos file, and 
!..must be called once before the helmeos routine is invoked.

!..declare local variables
      character (len=256) :: filename, message, temp_filename
      character (len=500) :: buf
      character (len=26) :: s26
      real(dp), target :: vec_ary(20)
      real(dp), pointer :: vec(:)
      integer          i,j,k,ios,imax,jmax,n
      real(dp) tsav,dsav
      logical, parameter :: dmp = .false.
      
       ierr = 0
       vec => vec_ary
       
!..read the normal helmholtz free energy table
       h% logtlo   = 3.0d0
       h% logthi   = 13.0d0
       h% logdlo   = -12.0d0
       h% logdhi   = 15.0d0
       
       h% templo = exp10(h% logtlo)
       h% temphi = exp10(h% logthi)
       h% denlo = exp10(h% logdlo)
       h% denhi = exp10(h% logdhi)

       imax = h% imax
       jmax = h% jmax
       h% logtstp  = (h% logthi - h% logtlo)/real(jmax-1,kind=dp)
       h% logtstpi = 1.0d0/h% logtstp
       h% logdstp  = (h% logdhi - h% logdlo)/real(imax-1,kind=dp)
       h% logdstpi = 1.0d0/h% logdstp
       
       ios = -1
       if (use_cache) then
         write(filename,'(2a)') trim(cache_dir), '/helm_table.bin'
         open(unit=19,file=trim(filename), &
               action='read',status='old',iostat=ios,form='unformatted')
       end if
       
       if (ios .eq. 0) then
         
          read(19) imax
          read(19) jmax
         
         if (imax /= h% imax .or. jmax /= h% jmax) then
            ios = 1 ! wrong cached info
         else
             read(19) h% f(1:imax,1:jmax)
             read(19) h% fd(1:imax,1:jmax)
             read(19) h% ft(1:imax,1:jmax)
             read(19) h% fdd(1:imax,1:jmax)
             read(19) h% ftt(1:imax,1:jmax)
             read(19) h% fdt(1:imax,1:jmax)
             read(19) h% fddt(1:imax,1:jmax)
             read(19) h% fdtt(1:imax,1:jmax)
             read(19) h% fddtt(1:imax,1:jmax)
             read(19) h% dpdf(1:imax,1:jmax)
             read(19) h% dpdfd(1:imax,1:jmax)
             read(19) h% dpdft(1:imax,1:jmax)
             read(19) h% dpdfdt(1:imax,1:jmax)
             read(19) h% ef(1:imax,1:jmax)
             read(19) h% efd(1:imax,1:jmax)
             read(19) h% eft(1:imax,1:jmax)
             read(19) h% efdt(1:imax,1:jmax)
             read(19) h% xf(1:imax,1:jmax)
             read(19) h% xfd(1:imax,1:jmax)
             read(19) h% xft(1:imax,1:jmax)
             read(19) h% xfdt(1:imax,1:jmax)
         
            do j=1,jmax
               tsav = h% logtlo + (j-1)*h% logtstp
               h% t(j) = exp10(tsav)
            enddo
            do i=1,imax
               dsav = h% logdlo + (i-1)*h% logdstp
               h% d(i) = exp10(dsav)
            enddo
         end if
         
         close(unit=19)
       end if

       if (ios .ne. 0) then
         
          write(filename,'(2a)') trim(data_dir), '/helm_table.dat'
          write(*,*) 'read  ', trim(filename) 
         
          ios = 0
          open(unit=19,file=trim(filename),action='read',status='old',iostat=ios)
          if (ios .ne. 0) then 
            write(*,'(3a,i6)') 'failed to open ', trim(filename), ' : ios ', ios
            ierr = -1
            return
          end if
         
          do j=1,jmax
            tsav = h% logtlo + (j-1)*h% logtstp
            h% t(j) = exp10(tsav)
            do i=1,imax
               dsav = h% logdlo + (i-1)*h% logdstp
               h% d(i) = exp10(dsav)
               read(19,'(a)',iostat=ierr) buf
               if (ierr == 0) call str_to_vector(buf, vec, n, ierr)
               if (ierr /= 0 .or. n /= 9) then
                  write(*,'(a)') 'failed while reading ' // trim(filename)
                  close(19)
                  return
               end if
               h% f(i,j) = vec(1)
               h% fd(i,j) = vec(2)
               h% ft(i,j) = vec(3)
               h% fdd(i,j) = vec(4)
               h% ftt(i,j) = vec(5)
               h% fdt(i,j) = vec(6)
               h% fddt(i,j) = vec(7)
               h% fdtt(i,j) = vec(8)
               h% fddtt(i,j) = vec(9)
               if (dmp) then
                  do k=1,9
                     write(*,'(1pd24.16)',advance='no') vec(k)
                  end do
                  write(*,*)
               end if
            enddo
          enddo

         !..read the pressure derivative with density table
          do j=1,jmax
           do i=1,imax
            read(19,'(a)',iostat=ierr) buf
            if (ierr == 0) call str_to_vector(buf, vec, n, ierr)
            if (ierr /= 0 .or. n /= 4) then
               write(*,'(a)') 'failed while reading ' // trim(filename)
               close(19)
               return
            end if
            h% dpdf(i,j) = vec(1)
            h% dpdfd(i,j) = vec(2)
            h% dpdft(i,j) = vec(3)
            h% dpdfdt(i,j) = vec(4)
            if (dmp) then
               do k=1,4
                  write(*,'(1pd24.16)',advance='no') vec(k)
               end do
               write(*,*)
            end if
           enddo
          enddo

         !..read the electron chemical potential table
          do j=1,jmax
           do i=1,imax
            read(19,'(a)',iostat=ierr) buf
            if (ierr == 0) call str_to_vector(buf, vec, n, ierr)
            if (ierr /= 0 .or. n /= 4) then
               write(*,'(a)') 'failed while reading ' // trim(filename)
               close(19)
               return
            end if
            h% ef(i,j) = vec(1)
            h% efd(i,j) = vec(2)
            h% eft(i,j) = vec(3)
            h% efdt(i,j) = vec(4)
            if (dmp) then
               do k=1,4
                  write(*,'(1pd24.16)',advance='no') vec(k)
               end do
               write(*,*)
            end if
           enddo
          enddo

         !..read the number density table
          do j=1,jmax
           do i=1,imax
            read(19,'(a)',iostat=ierr) buf
            if (ierr == 0) call str_to_vector(buf, vec, n, ierr)
            if (ierr /= 0 .or. n /= 4) then
               write(*,'(a)') 'failed while reading ' // trim(filename)
               close(19)
               return
            end if
            h% xf(i,j) = vec(1)
            h% xfd(i,j) = vec(2)
            h% xft(i,j) = vec(3)
            h% xfdt(i,j) = vec(4)
            if (dmp) then
               do k=1,4
                  write(*,'(1pd24.16)',advance='no') vec(k)
               end do
               write(*,*)
            end if
           enddo
          enddo

          close(unit=19)
          !..write cachefile
          
          if (dmp) call mesa_error(__FILE__,__LINE__,'helm_alloc')
          
          ios = -1
          if (use_cache) then
            write(filename,'(2a)') trim(cache_dir), '/helm_table.bin'
            write(temp_filename,'(2a)') trim(temp_cache_dir), '/helm_table.bin'
            write(*,*) 'write ', trim(filename) 
            open(unit=19,file=trim(switch_str(temp_filename, filename, use_mesa_temp_cache)),  &
             status='replace', iostat=ios,action='write',form='unformatted')
          end if
          
          if (ios == 0) then      
             write(19) imax
             write(19) jmax
             write(19) h% f(1:imax,1:jmax)
             write(19) h% fd(1:imax,1:jmax)
             write(19) h% ft(1:imax,1:jmax)
             write(19) h% fdd(1:imax,1:jmax)
             write(19) h% ftt(1:imax,1:jmax)
             write(19) h% fdt(1:imax,1:jmax)
             write(19) h% fddt(1:imax,1:jmax)
             write(19) h% fdtt(1:imax,1:jmax)
             write(19) h% fddtt(1:imax,1:jmax)
             write(19) h% dpdf(1:imax,1:jmax)
             write(19) h% dpdfd(1:imax,1:jmax)
             write(19) h% dpdft(1:imax,1:jmax)
             write(19) h% dpdfdt(1:imax,1:jmax)
             write(19) h% ef(1:imax,1:jmax)
             write(19) h% efd(1:imax,1:jmax)
             write(19) h% eft(1:imax,1:jmax)
             write(19) h% efdt(1:imax,1:jmax)
             write(19) h% xf(1:imax,1:jmax)
             write(19) h% xfd(1:imax,1:jmax)
             write(19) h% xft(1:imax,1:jmax)
             write(19) h% xfdt(1:imax,1:jmax)
             close(unit=19)
             
             if (use_mesa_temp_cache) call mv(temp_filename,filename,.true.)
             
          end if            
       
       end if 

       call setup_td_deltas(h, imax, jmax)

      end subroutine read_helm_table




      subroutine free_helm_table(h)
         use eos_def
         
         type (Helm_Table), pointer :: h
         
         call do_free(h% d)
         call do_free(h% t)
         
         call do_free2(h% f)
         call do_free2(h% fd)
         call do_free2(h% ft)
         call do_free2(h% fdd)
         call do_free2(h% ftt)
         call do_free2(h% fdt)
         call do_free2(h% fddt)
         call do_free2(h% fdtt)
         call do_free2(h% fddtt)

         !..for the pressure derivative with density tables
         call do_free2(h% dpdf)
         call do_free2(h% dpdfd)
         call do_free2(h% dpdft)
         call do_free2(h% dpdfdt)

         !..for chemical potential tables
         call do_free2(h% ef)
         call do_free2(h% efd)
         call do_free2(h% eft)
         call do_free2(h% efdt)

         !..for the number density tables
         call do_free2(h% xf)
         call do_free2(h% xfd)
         call do_free2(h% xft)
         call do_free2(h% xfdt)

         !..for storing the differences
         call do_free(h% dt_sav)
         call do_free(h% dt2_sav)
         call do_free(h% dti_sav)
         call do_free(h% dt2i_sav)
         call do_free(h% dt3i_sav)
         call do_free(h% dd_sav)
         call do_free(h% dd2_sav)
         call do_free(h% ddi_sav)
         call do_free(h% dd2i_sav)
         call do_free(h% dd3i_sav)
         
         deallocate(h)
         nullify(h)
         
         contains
         
         subroutine do_free(array_ptr)
            real(dp), pointer :: array_ptr(:)
            if (associated(array_ptr)) deallocate(array_ptr)
         end subroutine do_free
         
         subroutine do_free2(array_ptr)
            real(dp), pointer :: array_ptr(:,:)
            if (associated(array_ptr)) deallocate(array_ptr)
         end subroutine do_free2


      end subroutine free_helm_table
      
      

      end module helm_alloc
