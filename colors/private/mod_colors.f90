! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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
!
! ***********************************************************************

   module mod_colors
      use colors_def
      use math_lib
      use const_def
      use utils_lib
      
      implicit none
      private
      public :: do_colors_init, free_colors_all, Eval_Colors

      contains
      
      subroutine do_colors_init(num_files,fnames,num_colors,ierr)
         integer,intent(in) :: num_files
         integer,dimension(:),intent(in) :: num_colors
         character(len=*),dimension(:),intent(in) :: fnames 
         character(len=strlen) :: fname
         type (lgt_list), pointer :: thead =>null()
         
         integer, intent(out) :: ierr
         integer :: i
         
         ierr = 0
         num_thead=num_files
         
         if(num_thead<1)THEN
            color_is_initialized=.false.
            !Not a failure just dont have any files to read
            ierr=0
            return
         END IF
         
         ALLOCATE(thead_all(1:num_thead),STAT=ierr)
         if (ierr /= 0) then 
            write(*,*) 'colors can not allocate memory in do_colors_init'
            ierr=-1
            return
         end if
         
         bc_total_num_colors=0
         do i=1,num_thead
            allocate(thead_all(i)%thead)
            
            thead=>thead_all(i)%thead
            
            fname=fnames(i)
            thead_all(i)%n_colors=num_colors(i)
            if(len(fname)==0)THEN
               exit
            end if
            
            if(thead_all(i)%n_colors<1)THEN
               write(*,*) "num_colors must be > 0"
               ierr=-1
               return
            end if
            
            call init_colors(fname,thead,thead_all(i)%color_names,thead_all(i)%n_colors,ierr)

            thead_all(i)%thead=>thead
            bc_total_num_colors=bc_total_num_colors+thead_all(i)%n_colors
            nullify(thead)
         end do

         color_is_initialized=.true.
         
      end subroutine do_colors_init


      subroutine init_colors(fname, thead, col_names, n_colors, ierr)
         integer, intent(out) :: ierr
         character(len=*),intent(in) :: fname
         type (lgt_list), pointer :: thead 
         character(len=*),dimension(:) :: col_names
         integer, intent(in) :: n_colors
         call Read_Colors_Data(fname, thead, col_names, n_colors, ierr) 
      end subroutine init_colors
      
      subroutine free_colors_all
         integer :: i
         type (lgt_list), pointer :: thead => null()
         
         do i=1,num_thead
            if(associated(thead_all(i)%thead))then
               thead=>thead_all(i)%thead
               call free_colors(thead)
            end if
         end do
         deallocate(thead_all)
         
      end subroutine free_colors_all

      subroutine free_colors(thead)
         type (lgt_list), pointer :: thead
         type(lgt_list),pointer :: tlist => null(), tnxt => null()

         tlist => thead
         do while (associated(tlist))
            tnxt => tlist% nxt
            call free_glist(tlist% glist)
            deallocate(tlist)
            tlist => tnxt
         end do
         nullify(thead)
         
         contains
         
         subroutine free_glist(gptr)
            type (lgg_list), pointer :: gptr          
            type (lgg_list), pointer :: glist   => null()         
            type (lgg_list), pointer :: gnxt => null()   
            glist => gptr     
            do while (associated(glist))
               gnxt => glist% nxt
               call free_zlist(glist% zlist)
               deallocate(glist)
               glist => gnxt
            end do         
         end subroutine free_glist
         
         subroutine free_zlist(zptr)
            type (lgz_list), pointer :: zptr
            type (lgz_list), pointer :: zlist => null()
            type (lgz_list), pointer :: znxt => null()
            zlist => zptr        
            do while (associated(zlist))
               znxt => zlist% nxt
               deallocate(zlist)
               zlist => znxt
            end do
         end subroutine free_zlist
         
      end subroutine free_colors


      subroutine show_tree(thead)
         type (lgt_list), pointer :: thead
         type (lgt_list), pointer :: tlist => null(), tnxt => null()
         
         tlist => thead
         do while (associated(tlist))
            write(*,*) 'tlist lgt', exp10(tlist% lgt)
            tnxt => tlist% nxt
            call show_glist(tlist% glist)
            tlist => tnxt
         end do
         
         contains
         
         subroutine show_glist(gptr)
            type (lgg_list), pointer :: gptr
            type (lgg_list), pointer :: glist => null(), gnxt => null()
            
            glist => gptr
            do while (associated(glist))
               write(*,*) 'glist% lgg', glist% lgg
               gnxt => glist% nxt
               call show_zlist(glist% zlist)
               glist => gnxt
            end do         
         end subroutine show_glist
         
         subroutine show_zlist(zptr)
            type (lgz_list), pointer :: zptr          
            type (lgz_list), pointer :: zlist => null()           
            type (lgz_list), pointer :: znxt => null()   
            zlist => zptr     
            do while (associated(zlist))
               write(*,*) 'zlist% lgz', zlist% lgz
               znxt => zlist% nxt
               zlist => znxt
            end do         
         end subroutine show_zlist
         
      end subroutine show_tree
      
      
      subroutine Read_One_Colors_Data(fname, thead, n_colors, col_names, ierr)
         integer, intent(out) :: ierr ! 0 means ok
         type (lgt_list), pointer,intent(inout) :: thead

         ! read file and build lists
         integer :: ios, cnt
         integer,intent(in) :: n_colors
         character (len=*) :: fname
         real(dp) :: lgt, lgg, lgz
         type (lgg_list), pointer :: glist => null()
         type (lgt_list), pointer :: tlist => null()
         type (lgz_list), pointer :: zlist => null()
         real(dp), dimension(max_num_bcs_per_file) :: colors
         character(len=*),dimension(:),intent(out) :: col_names
         character(len=256) :: tmp_cols(3+n_colors)

         character(len=4096) :: tmp
         integer :: num_entries, num_made, IO_UBV
         
         include 'formats'
         
              ! Try local folder first
         open(NEWUNIT=IO_UBV, FILE=trim(fname), ACTION='READ', STATUS='OLD', IOSTAT=ios)
         if(ios/=0) THEN
            ! Try colors data dir
            open(NEWUNIT=IO_UBV, FILE=trim(mesa_data_dir)//'/colors_data/'//trim(fname), ACTION='READ', STATUS='OLD', IOSTAT=ios)
            if (ios /= 0) then ! failed to open lcb98cor.dat
               write(*,*) 'colors_init: failed to open ' // trim(fname)
               write(*,*) 'or ',trim(mesa_data_dir)//'/colors_data/'//trim(fname)
               write(*,*)
               write(*,*) 'Please check that the file exists in your local work directory'
               write(*,*) 'or in ',trim(mesa_data_dir)//'/colors_data/'
               ierr = 1; return
            endif
         end if
         
         ierr = 0
         num_entries = 0
         cnt = 0
         tmp = ''
         !First line should be a header and containing the name of the colours
         !#teff logg m_div_h col_1 col_2 etc
         read(IO_UBV,'(a)') tmp

         call split_line(tmp,3+n_colors,tmp_cols)
      
         col_names(1:n_colors) = tmp_cols(4:n_colors+3)
         
         do while (.true.)
            read(IO_UBV,fmt=*,iostat=ios) lgt, lgg, lgz, colors(1:n_colors)
            if (ios /= 0) exit
            cnt = cnt + 1
            lgt = log10(lgt)  

            if(cnt==1) thead%lgt = lgt
            
            call get_tlist(thead, lgt, tlist, ierr)
            if (ierr /= 0) exit            
            
            call get_glist(tlist% glist, lgg, glist, ierr)
            if (ierr /= 0) exit
            
            call get_zlist(glist% zlist, lgz, zlist, num_entries, ierr)
            if (ierr /= 0) exit
            
            if(zlist% colors(1) > -1d98) then
               write(*,*) "Warning found duplicated color data for (T, log g, M/H)=", 10**lgT, lgg, lgz
            end if

            zlist% colors = colors

         end do
         
         close(IO_UBV)
         if (ierr /= 0) return
         
         num_made = 0
         
         tlist => thead
         lgt = 1d99
         do while (associated(tlist))
            if (tlist% lgt >= lgt) then ! bad tlist order
               ierr = -1; return
            end if
            lgt = tlist% lgt
            glist => tlist% glist
            lgg = 1d99
            do while (associated(glist))
               if (glist% lgg >= lgg) then ! bad glist order
                  ierr = -2; return
               end if
               lgg = glist% lgg
               zlist => glist% zlist
               lgz = 1d99
               do while (associated(zlist))
                  if (zlist% lgz >= lgz) then ! bad zlist order
                     ierr = -3; return
                  end if
                  lgz = zlist% lgz
                  num_made = num_made + 1
                  zlist => zlist% nxt
               end do
               glist => glist% nxt
            end do
            tlist => tlist% nxt
         end do
         
         if(num_entries /= num_made)then
            write(0,*) "Error fond less colors than expected ",num_entries,num_made  
            stop
         end if
         
         
      end subroutine Read_One_Colors_Data
      
      
      subroutine Read_Colors_Data(fname, thead, col_names, n_colors, ierr)
         use const_def, only: mesa_data_dir
         integer, intent(out) :: ierr ! 0 means ok
         type (lgt_list), pointer,intent(inout) :: thead
         character (len=*),intent(in) :: fname
         character(len=*),dimension(:),intent(out) :: col_names
         integer, intent(in) :: n_colors
         
         Call Read_One_Colors_Data(fname, thead, n_colors, col_names, ierr)
         if (ierr /= 0) THEN
            write(*,*) "Read_Colors_One_Data error"
            stop
         end if

      end subroutine Read_Colors_Data

      subroutine get_tlist(head, lgt, tlist, ierr)
         type (lgt_list), pointer :: head
         real(dp), intent(in) :: lgt
         type (lgt_list), pointer :: tlist
         integer, intent(out) :: ierr ! 0 means ok
         
         type (lgt_list), pointer :: t1=>null(), t2=>null()
         
         ierr = 0
         
         if (.not. associated(head)) then ! first time
            if (.not. alloc_tlist()) return
            head => tlist
            return
         end if
      
         if (head% lgt == lgt) then ! matches head of list
            tlist => head
            return
         end if
         
         if (head% lgt < lgt) then ! becomes new head of list
            if (.not. alloc_tlist()) return
            tlist% nxt => head
            head => tlist
            return
         end if
         
         ! check list
         t1 => head
         do while (associated(t1% nxt))
            t2 => t1% nxt
            if (t2% lgt == lgt) then
               tlist => t2; return
            end if
            if (t2% lgt < lgt) then ! insert new one before t2 
               if (.not. alloc_tlist()) return
               tlist% nxt => t2
               t1% nxt => tlist
               return
            end if
            t1 => t2
         end do
         ! add to end of list after t1
         if (.not. alloc_tlist()) return
         t1% nxt => tlist
         
         contains
         
         logical function alloc_tlist()
            integer :: istat
            allocate(tlist,stat=istat)
            if (istat /= 0) then
               alloc_tlist = .false.; ierr = -1; return
            end if
            nullify(tlist% glist)
            nullify(tlist% nxt)
            tlist% lgt = lgt        
            alloc_tlist = .true.
         end function alloc_tlist

      end subroutine get_tlist


      subroutine get_glist(head, lgg, glist, ierr)
         type (lgg_list), pointer :: head
         real(dp), intent(in) :: lgg
         type (lgg_list), pointer :: glist 
         integer, intent(out) :: ierr
         
         type (lgg_list), pointer :: g1 => null(), g2 => null()
         
         ierr = 0
         
         if (.not. associated(head)) then ! first time
            if (.not. alloc_glist()) return
            head => glist
            return
         end if
         
         if (head% lgg == lgg) then ! matches head of list
            glist => head
            return
         end if
         
         if (head% lgg < lgg) then ! becomes new head of list
            if (.not. alloc_glist()) return
            glist% nxt => head
            head => glist
            return
         end if
         
         ! check list
         g1 => head
         do while (associated(g1% nxt))
            g2 => g1% nxt
            if (g2% lgg == lgg) then
               glist => g2; return
            end if
            if (g2% lgg < lgg) then ! insert new one before g2 
               if (.not. alloc_glist()) return
               glist% nxt => g2
               g1% nxt => glist
               return
            end if
            g1 => g2
         end do
         ! add to end of list after g1
         if (.not. alloc_glist()) return
         g1% nxt => glist
         
         contains
         
         logical function alloc_glist()
            integer :: istat
            allocate(glist,stat=istat)
            if (istat /= 0) then ! allocate failed in alloc_glist
               alloc_glist = .false.; ierr = -1; return
            end if
            nullify(glist% zlist)
            nullify(glist% nxt)
            glist% lgg = lgg  
            alloc_glist = .true.    
         end function alloc_glist
         
      end subroutine get_glist                     

      subroutine get_zlist(head, lgz, zlist, num_entries, ierr)
         type (lgz_list), pointer :: head
         real(dp), intent(in) :: lgz
         type (lgz_list), pointer :: zlist
         integer, intent(out) :: ierr ! 0 means ok
         integer,intent(inout) :: num_entries
         
         type (lgz_list), pointer :: z1=>null(), z2=>null()
         
         ierr = 0
         
         if (.not. associated(head)) then ! first time
            if (.not. alloc_zlist()) return
            head => zlist
            return
         end if
         
         if (head% lgz == lgz) then ! matches head of list
            zlist => head
            return
         end if
         
         if (head% lgz < lgz) then ! becomes new head of list
            if (.not. alloc_zlist()) return
            zlist% nxt => head
            head => zlist
            return
         end if
         
         ! check list
         z1 => head
         do while (associated(z1% nxt))
            z2 => z1% nxt
            if (z2% lgz == lgz) then
               zlist => z2; return
            end if
            if (z2% lgz < lgz) then ! insert new one before z2 
               if (.not. alloc_zlist()) return
               zlist% nxt => z2
               z1% nxt => zlist
               return
            end if
            z1 => z2
         end do
         ! add to end of list after z1
         if (.not. alloc_zlist()) return
         z1% nxt => zlist
         
         contains
         
         logical function alloc_zlist()
            integer :: istat
            allocate(zlist,stat=istat)
            if (istat /= 0) then
               alloc_zlist = .false.; ierr = -1; return
            end if
            nullify(zlist% nxt)
            num_entries=num_entries+1
            zlist% lgz = lgz     
            alloc_zlist = .true.
         end function alloc_zlist
         
      end subroutine get_zlist
                     

      subroutine Eval_Colors(log_Teff,log_g, M_div_h_in, results, thead, n_colors, ierr)
         real(dp), intent(in)  :: log_Teff ! log10 of surface temp
         real(dp), intent(in)  :: M_div_h_in ! [M/H]
         real(dp),dimension(:), intent(out) :: results
         real(dp), intent(in) :: log_g
         integer, intent(in) :: n_colors
         integer, intent(out) :: ierr
         type (lgt_list), pointer,intent(inout) :: thead
         
         !real(dp), parameter :: Zsol = 0.02d0, colors_bol_sun = 4.746d0
         real(dp) :: colors_bol
         
         real(dp) :: lgg, lgz, lgt, alfa, beta
         real(dp),dimension(max_num_bcs_per_file) :: results1, results2
         type (lgt_list), pointer :: tlist => null(), tnxt => null()
         
         lgg=log_g
         lgz=M_div_h_in
         lgt=log_Teff
         
         if (.not. associated(thead)) then
            ierr = -1; return
         end if
                  
         ierr = 0
         
         if (lgt >= thead% lgt) then ! use the largest lgt
            call get_glist_results( &
               thead% glist, lgg, lgz, results1, thead%n_colors, ierr)
            if (ierr /= 0) return
            results = results1
         else
            
            tlist => thead
            do while (associated(tlist% nxt))
               tnxt => tlist% nxt
               if (lgt == tnxt% lgt) then ! use tnxt
                  call get_glist_results( &
                     tnxt% glist, lgg, lgz, results1, n_colors, ierr)
                     if (ierr /= 0) return
                  results = results1
                  return
               end if
               if (lgt >= tnxt% lgt) then ! interpolate between tlist and tnxt
                  call get_glist_results(tlist% glist, lgg, lgz, results1, n_colors, ierr)
                  if (ierr /= 0) return
                  call get_glist_results(tnxt% glist, lgg, lgz, results2, n_colors, ierr)
                  if (ierr /= 0) return
                  alfa = (lgt - tnxt% lgt) / (tlist% lgt - tnxt% lgt)
                  beta = 1.d0 - alfa
                  results(1:n_colors) = alfa * results1(1:n_colors) + beta * results2(1:n_colors)
                  return
               end if
               tlist => tnxt
            end do
      
            if (.not. (associated(tlist% nxt))) then
               ! use the smallest lgt
               call get_glist_results( &
                  tlist% glist, lgg, lgz, results1, n_colors,ierr)
               if (ierr /= 0) return
               results = results1
            end if
         
         end if
         
         !results(bol) = colors_bol
                  
      end subroutine Eval_Colors
      
      
      subroutine get_glist_results(gptr, lgg, lgz, results, n_colors, ierr)
         type (lgg_list), pointer :: gptr
         real(dp), intent(in) :: lgg, lgz
         real(dp),dimension(:), intent(out) :: results
         integer,intent(in) :: n_colors
         integer,intent(out) :: ierr
         
         type (lgg_list), pointer :: glist => null(), gnxt => null()
         real(dp),dimension(max_num_bcs_per_file) :: results1, results2
         real(dp) :: alfa, beta
               
         
         glist => gptr
         ierr = -1
         if (.not. associated(glist)) then
            write(*,*) 'bad glist for get_glist_results'
            ierr = -1
            return
         end if
         
         if (lgg >= glist% lgg) then ! use the largest lgg
            call get_zlist_results( &
               glist% zlist, lgz, results1, n_colors, ierr)
            if (ierr /= 0) return
            results = results1
            return
         end if

         do while (associated(glist% nxt))
            gnxt => glist% nxt
            if (lgg == gnxt% lgg) then ! use gnxt
               call get_zlist_results(gnxt% zlist, lgz, results1, n_colors, ierr)
               if (ierr /= 0) return
               results=results1
               return
            end if
            if (lgg >= gnxt% lgg) then ! interpolate between lgg's
               call get_zlist_results(glist% zlist, lgz, results1, n_colors, ierr)
               if (ierr /= 0) return
               call get_zlist_results(gnxt% zlist, lgz, results2, n_colors, ierr)
               if (ierr /= 0) return               
               alfa = (lgg - gnxt% lgg) / (glist% lgg - gnxt% lgg)
               beta = 1.d0 - alfa
               results(1:n_colors) = alfa * results1(1:n_colors) + beta * results2(1:n_colors)
               return
            end if
            glist => gnxt
         end do
      
         ! use the smallest lgg
         call get_zlist_results(glist% zlist, lgz, results1, n_colors, ierr)
         if (ierr /= 0) return
         results = results1
               
      end subroutine get_glist_results
      
      
      subroutine get_zlist_results(zptr, lgz, results, n_colors, ierr)
         type (lgz_list), pointer :: zptr
         real(dp), intent(in) :: lgz
         real(dp),dimension(:), intent(out) :: results
         integer, intent(in) :: n_colors
         integer, intent(out) :: ierr
         
         type (lgz_list), pointer :: zlist => null(), znxt => null()
         real(dp) :: alfa, beta
         
         zlist => zptr
         
         ierr = 0
         if (.not. associated(zlist)) then
            write(*,*) 'bad zlist for get_zlist_results'
            ierr = -1
            return
         end if
         
         if (lgz >= zlist% lgz) then ! use the largest lgz
            results = zlist% colors
            return
         end if

         do while (associated(zlist% nxt))
            znxt => zlist% nxt
            if (lgz == znxt% lgz) then ! use znxt
               results(1:n_colors) = znxt% colors(1:n_colors)
               return
            end if
            if (lgz >= znxt% lgz) then ! interpolate between zlist and znxt
               alfa = (lgz - znxt% lgz) / (zlist% lgz - znxt% lgz)
               beta = 1 - alfa
               results(1:n_colors) = alfa * zlist% colors(1:n_colors) + beta * znxt% colors(1:n_colors)
               return
            end if
            zlist => znxt
         end do
      
         ! use the smallest lgz
         results = zlist% colors
      
      end subroutine get_zlist_results


      end module mod_colors

