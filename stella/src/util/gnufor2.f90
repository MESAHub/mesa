! Author: Alexey Kuznetsov
! Modified: 28/12/2008
! this Fortran90 module contains a collection of subroutines for plotting data,
! including 2D, 3D plots, surfaces, polar coordinates, histograms
! it is a modification of the GNUFOR interface written by John Burkardt:
! http://orion.math.iastate.edu/burkardt/g_src/gnufor/gnufor.html 
!***********************************************************************************
	module gnufor2
	implicit none
!***********************************************************************************
! these are default parameters which control linewidth, colors and terminal
!***********************************************************************************
	character(len=3), parameter	:: default_linewidth='1'
	character(len=100), parameter	:: default_color1='blue'
	character(len=100), parameter	:: default_color2='dark-green'
	character(len=100), parameter	:: default_color3='orange-red'
	character(len=100), parameter	:: default_color4='dark-salmon'
	character(len=100), parameter	:: default_terminal='wxt'
	character(len=100), parameter	:: default_palette='CMY'
!***********************************************************************************
	interface plot
		module procedure plot_1
		module procedure plot_2
		module procedure plot_3
		module procedure plot_4
	end interface
!***********************************************************************************
	interface surf
		module procedure surf_1
		module procedure surf_2
		module procedure surf_3
	end interface
!***********************************************************************************
	interface image
		module procedure image_1
		module procedure image_2
		module procedure image_3
		module procedure image_4
	end interface
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	contains
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	function my_date_and_time() result(f_result)
!***********************************************************************************
! this function creates a string with current date and time
! it is a default method to name output files
!***********************************************************************************
	implicit none
        character(len=8)  :: date
        character(len=10) :: time
        character(len=33) :: f_result
!***********************************************************************************
	call date_and_time(date,time)
	f_result= 'date_'//date(7:8)//'-'//date(5:6)//'-'//date(1:4)//'_time_'//time(1:2)//':'//time(3:4)//':'//time(5:10)
!***********************************************************************************
	end function my_date_and_time
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	function output_terminal(terminal) result(f_result)
!***********************************************************************************
	implicit none
	character(len=*),intent(in)	:: terminal
	integer, parameter		:: Nc=35
	character(len=Nc)		:: f_result
!***********************************************************************************
	select case(terminal)
		case('ps')
			f_result='postscript landscape color'
		case default
			f_result=terminal
	end select
!***********************************************************************************
	end function output_terminal
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine image_4(x,y,rgb,pause,terminal,filename,persist,input)
!***********************************************************************************
! this is the most general subroutine for generating 3D plots.
! The data is contained in a 3D array z(:,:,:)
!***********************************************************************************
	implicit none
	real(kind=4), intent(in)	:: x(:), y(:)
	integer,intent(in)		:: rgb(:,:,:)
	real(kind=4), optional		:: pause
	character(len=*),optional	:: terminal, filename, persist, input
	integer				:: nx, ny
	integer 			:: i, j, ierror, ios, file_unit
	character(len=100)		:: data_file_name, command_file_name, my_pause, my_persist,&
		& xrange1,xrange2,yrange1,yrange2
!***********************************************************************************
	nx=size(rgb(1,:,1))
	ny=size(rgb(1,1,:))
	if ((size(x).ne.nx).or.(size(y).ne.ny)) then
		ierror=1
		print *,'image2 ERROR: sizes of x(:),y(:) and gray(:,:) are not compatible'
		stop
	end if
	write (xrange1,'(e15.7)') minval(x)
	write (xrange2,'(e15.7)') maxval(x)
	write (yrange1,'(e15.7)') minval(y)
	write (yrange2,'(e15.7)') maxval(y)
	do j=1,ny
	do i=1,nx
		if ((maxval(rgb(:,i,j))>255).or.(minval(rgb(:,i,j))<0)) then
			print *,'image ERROR: a value of rgb(:,:,:) array is outside [0,255]'
			stop
		end if
	end do
	end do
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************	
	do j=1,ny
		do i=1,nx
			write (file_unit,'(2E12.4,3I5)') x(i),y(j),rgb(1,i,j),rgb(2,i,j),rgb(3,i,j)
		end do
		write (file_unit,'(a)')
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'set nokey'
	write ( file_unit, '(a)' ) 'set xrange ['// trim(xrange1) // ':'// trim(xrange2) //']'
	write ( file_unit, '(a)' ) 'set yrange ['// trim(yrange1) // ':'// trim(yrange2) //']'
	write ( file_unit, '(a)' ) 'unset colorbox'
!***********************************************************************************	
	write ( file_unit, '(a)' ) 'plot "' // trim ( data_file_name ) // &
	& '" with rgbimage'
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine image_4
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine image_3(rgb,pause,terminal,filename,persist,input)
!***********************************************************************************
! this is the most general subroutine for generating 3D plots.
! The data is contained in a 3D array z(:,:,:)
!***********************************************************************************
	implicit none
	integer, intent(in)		:: rgb(:,:,:)
	real(kind=4), optional		:: pause
	character(len=*),optional	:: terminal, filename, persist, input
	integer				:: nx, ny
	integer 			:: i, j, ierror, ios, file_unit
	character(len=100)		:: data_file_name, command_file_name, my_pause, my_persist
!***********************************************************************************
	nx=size(rgb(1,:,1))
	ny=size(rgb(1,1,:))
	do j=1,ny
	do i=1,nx
		if ((maxval(rgb(:,i,j))>255).or.(minval(rgb(:,i,j))<0)) then
			print *,'image ERROR: a value of rgb(:,:,:) array is outside [0,255]'
			stop
		end if
	end do
	end do
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************	
	do j=1,ny
		do i=1,nx
			write (file_unit,'(5I5)') i,j,rgb(1,i,j),rgb(2,i,j),rgb(3,i,j)
		end do
		write (file_unit,'(a)')
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'set nokey'
	write ( file_unit, '(a)' ) 'unset border'
	write ( file_unit, '(a)' ) 'unset xtics'
	write ( file_unit, '(a)' ) 'unset ytics'
	write ( file_unit, '(a)' ) 'unset colorbox'
!***********************************************************************************	
	write ( file_unit, '(a)' ) 'plot "' // trim ( data_file_name ) // &
	& '" with rgbimage'
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine image_3
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine image_2(x,y,gray,pause,palette,terminal,filename,persist,input)
!***********************************************************************************
! this is the most general subroutine for generating 3D plots.
! The data is contained in a 3D array z(:,:,:)
!***********************************************************************************
	implicit none
	real(kind=4), intent(in)	:: x(:), y(:), gray(:,:)
	real(kind=4), optional		:: pause
	character(len=*),optional	:: palette, terminal, filename, persist, input
	integer				:: nx, ny
	integer 			:: i, j, ierror, ios, file_unit
	character(len=100)		:: data_file_name, command_file_name, my_pause, my_persist,&
		& xrange1,xrange2,yrange1,yrange2
!***********************************************************************************
	nx=size(gray(:,1))
	ny=size(gray(1,:))
	if ((size(x).ne.nx).or.(size(y).ne.ny)) then
		ierror=1
		print *,'image2 ERROR: sizes of x(:),y(:) and gray(:,:) are not compatible'
		stop
	end if
	write (xrange1,'(e15.7)') minval(x)
	write (xrange2,'(e15.7)') maxval(x)
	write (yrange1,'(e15.7)') minval(y)
	write (yrange2,'(e15.7)') maxval(y)
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************	
	do j=1,ny
		do i=1,nx
			write (file_unit,'(3E12.4)') x(i), y(j), gray(i,j)
		end do
		write (file_unit,'(a)')
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'set nokey'
	write ( file_unit, '(a)' ) 'set xrange ['// trim(xrange1) // ':'// trim(xrange2) //']'
	write ( file_unit, '(a)' ) 'set yrange ['// trim(yrange1) // ':'// trim(yrange2) //']'
	write ( file_unit, '(a)' ) 'unset colorbox'
	if (present(palette)) then
	if ((trim(palette).ne.'RGB').and.(trim(palette).ne.'HSV').and.(trim(palette).ne.'CMY').and.&
		& (trim(palette).ne.'YIQ').and.(trim(palette).ne.'XYZ')) then
		write ( file_unit, '(a)' ) 'set palette '// trim(palette)
	else
		write ( file_unit, '(a)' ) 'set palette model '// trim(palette)
	end if
	else
		write ( file_unit, '(a)' ) 'set palette model '// trim(default_palette)
	end if
!***********************************************************************************	
	write ( file_unit, '(a)' ) 'plot "' // trim ( data_file_name ) // &
	& '" with image'
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine image_2
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine image_1(gray,pause,palette,terminal,filename,persist,input)
!***********************************************************************************
! this is the most general subroutine for generating 3D plots.
! The data is contained in a 3D array z(:,:,:)
!***********************************************************************************
	implicit none
	real(kind=4), intent(in)	:: gray(:,:)
	real(kind=4), optional		:: pause
	character(len=*),optional	:: palette, terminal, filename, persist, input
	integer				:: nx, ny
	integer 			:: i, j, ierror, ios, file_unit
	character(len=100)		:: data_file_name, command_file_name, my_pause, my_persist
!***********************************************************************************
	nx=size(gray(:,1))
	ny=size(gray(1,:))
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************	
	do j=1,ny
		do i=1,nx
			write (file_unit,'(I5,I5,E15.7)') i,j,gray(i,j)
		end do
		write (file_unit,'(a)')
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'set nokey'
	write ( file_unit, '(a)' ) 'unset border'
	write ( file_unit, '(a)' ) 'unset xtics'
	write ( file_unit, '(a)' ) 'unset ytics'
	write ( file_unit, '(a)' ) 'unset colorbox'
	if (present(palette)) then
	if ((trim(palette).ne.'RGB').and.(trim(palette).ne.'HSV').and.(trim(palette).ne.'CMY').and.&
		& (trim(palette).ne.'YIQ').and.(trim(palette).ne.'XYZ')) then
		write ( file_unit, '(a)' ) 'set palette '// trim(palette)
	else
		write ( file_unit, '(a)' ) 'set palette model '// trim(palette)
	end if
	else
		write ( file_unit, '(a)' ) 'set palette model '// trim(default_palette)
	end if
!***********************************************************************************	
	write ( file_unit, '(a)' ) 'plot "' // trim ( data_file_name ) // &
	& '" with image'
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine image_1
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine plot3d(x,y,z,pause,color,terminal,filename,persist,input,linewidth)
!***********************************************************************************
! this subroutine plots 3D curve, given by three arrays x,y,z
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: x(:),y(:),z(:)
	real(kind=4), optional		:: pause, linewidth
	character(len=*),optional	:: color, terminal, filename, persist, input
	integer 			:: i, j, ierror, ios, file_unit, nx
	character(len=100)		:: data_file_name, command_file_name, my_color, my_pause, my_persist, my_linewidth
!***********************************************************************************
! prepare the data
	nx=size(x)
	if ((size(x).ne.size(y)).or.(size(x).ne.size(z))) then
		print *,'subroutine plot3d ERROR: incompatible sizes of x(:),y(:) and z(:)'
		stop
	end if
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************	
	do i=1,nx
		write (file_unit,'(3E15.7)') x(i), y(i), z(i)
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_persist='persist'
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
		&// trim(my_persist) // ' title "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'set nokey'
!***********************************************************************************
	if (present(linewidth)) then
		write (	my_linewidth,'(e9.3)') linewidth
	else
		my_linewidth=trim(default_linewidth)
	end if	
	if (present(color)) then
		my_color='"'//trim(color)//'"'
	else
		my_color='"'//trim(default_color1)//'"'
	end if	
	write ( file_unit, '(a)' ) 'splot "' // trim ( data_file_name ) // &
	'" using 1:2:3 with lines linecolor rgb' // my_color //' linewidth ' // my_linewidth
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine plot3d
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine hist(x,n,pause,color,terminal,filename,persist,input)
!***********************************************************************************
! this subroutine plots the histogram of data contained in array x, using n bins
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: x(:) !the data to plot
	integer, intent(in)		:: n !the number of intervals
	real(kind=4), optional		:: pause
	character(len=*),optional	:: color, terminal, filename, persist, input
	integer 			:: i, j, ierror, ios, file_unit, nx
	character(len=100)		:: data_file_name, command_file_name, yrange, xrange1, xrange2, my_color, &
				& xtic_start, dxtic, xtic_end, my_pause, my_persist
	real(kind=8)			:: xmin, xmax, xhist(0:n), yhist(n+1), dx
!***********************************************************************************
! prepare the data
	nx=size(x)
	xmin=minval(x)
	xmax=maxval(x)
	dx=(xmax-xmin)/n
	do i=0,n
		xhist(i)=xmin+i*dx
	end do
	yhist=0.0d0
	do i=1,nx
		j=floor((x(i)-xmin)/dx)+1
		yhist(j)=yhist(j)+1
	end do
!***********************************************************************************
	write (dxtic,'(e15.7)') dx
	write (yrange,'(e15.7)') maxval(yhist)*1.2
	write (xrange1,'(e15.7)') xmin-(n/10.0)*dx
	write (xrange2,'(e15.7)') xmax+(n/10.0)*dx
	xtic_start=xrange1
	xtic_end=xrange2
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************	
	do i=1,n
		write (file_unit,'(2E15.7)') (xhist(i-1)+0.5*dx), yhist(i)
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_persist='persist'
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'set nokey'
	write ( file_unit, '(a)' ) 'set yrange [0.0:'// trim(yrange) //']'
	write ( file_unit, '(a)' ) 'set xrange ['// trim(xrange1) // ':'// trim(xrange2) //']'
	write ( file_unit, '(a)' ) 'set xtic nomirror rotate by -45 '
	write ( file_unit, '(a)' ) 'set xtics '// trim(xtic_start) // ','// trim(dxtic) // ','// trim(xtic_end)
	write ( file_unit, '(a)' ) 'set style data histograms'
	write ( file_unit, '(a)' ) 'set style fill solid border -1'
!***********************************************************************************
	if (present(color)) then
		my_color='"'//color//'"'
	else
		my_color='"'//trim(default_color1)//'"'
	end if	
	write ( file_unit, '(a)' ) 'plot "' // trim ( data_file_name ) // &
	'" using 1:2 with boxes linecolor rgb' // trim(my_color)
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine hist
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine surf_3(x,y,z,pause,palette,terminal,filename,pm3d,contour,persist,input)
!***********************************************************************************
! this subroutine plots a surface. x and y are arrays needed to generate the x-y grid
! z(:,:) is a 2D array
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: x(:),y(:),z(:,:)
	real(kind=4), optional		:: pause
	real(kind=8)			:: xyz(3,size(z(:,1)),size(z(1,:)))
	character(len=*),optional	:: palette, terminal, filename, pm3d, contour, persist, input
	integer				:: nx, ny
	integer 			:: i, j
!***********************************************************************************
	nx=size(z(:,1))
	ny=size(z(1,:))
	if ((size(x).ne.nx).or.(size(y).ne.ny)) then
		print *,'subroutine surf_3 ERROR: sizes of x(:),y(:), and z(:,:) are incompatible'
		stop
	end if 
!***********************************************************************************
	do i=1,nx
		do j=1,ny
			xyz(1,i,j)=x(i)
			xyz(2,i,j)=y(j)
			xyz(3,i,j)=z(i,j)
		end do
	end do
	call surf_1(xyz,pause,palette,terminal,filename,pm3d,contour,persist,input)
!***********************************************************************************	
	end subroutine surf_3
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine surf_2(z,pause,palette,terminal,filename,pm3d,contour,persist,input)
!***********************************************************************************
! this subroutine plots a surface. The only input is a 2D array z(:,:), the x-y grid 
! is generated automatically
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: z(:,:)
	real(kind=4), optional		:: pause
	real(kind=8)			:: xyz(3,size(z(:,1)),size(z(1,:)))
	character(len=*),optional	:: palette, terminal, filename, pm3d, contour, persist, input
	integer				:: nx, ny
	integer 			:: i, j
!***********************************************************************************
	nx=size(z(:,1))
	ny=size(z(1,:))
	do i=1,nx
		do j=1,ny
			xyz(1,i,j)=dble(i)
			xyz(2,i,j)=dble(j)
			xyz(3,i,j)=z(i,j)
		end do
	end do
	call surf_1(xyz,pause,palette,terminal,filename,pm3d,contour,persist,input)
!***********************************************************************************	
	end subroutine surf_2
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine surf_1(xyz,pause,palette,terminal,filename,pm3d,contour,persist,input)
!***********************************************************************************
! this is the most general subroutine for generating 3D plots.
! The data is contained in a 3D array z(:,:,:)
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: xyz(:,:,:)
	real(kind=4), optional		:: pause
	character(len=*),optional	:: palette, terminal, filename, pm3d, contour, persist, input
	integer				:: nx, ny, nrow
	integer 			:: i, j, ierror, ios, file_unit
	character(len=100)		:: data_file_name, command_file_name, my_pause, my_persist
!***********************************************************************************
	nx=size(xyz(1,:,1))
	ny=size(xyz(1,1,:))
	nrow=nx*ny
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************	
	do j=1,ny
		do i=1,nx
			write (file_unit,'(3E15.7)') xyz(1:3,i,j)
		end do
		write (file_unit,'(a)')
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'set nokey'
	if (present(palette)) then
	if ((trim(palette).ne.'RGB').and.(trim(palette).ne.'HSV').and.(trim(palette).ne.'CMY').and.&
		& (trim(palette).ne.'YIQ').and.(trim(palette).ne.'XYZ')) then
		write ( file_unit, '(a)' ) 'set palette '// trim(palette)
	else
		write ( file_unit, '(a)' ) 'set palette model '// trim(palette)
	end if
	else
		write ( file_unit, '(a)' ) 'set palette model '// trim(default_palette)
	end if
!***********************************************************************************
	if (present(pm3d)) then
		write ( file_unit, '(a)' ) 'set '// pm3d	
	else
		write ( file_unit, '(a)' ) 'set surface'
		if (present(contour)) then
			if (contour=='surface') then 
				write ( file_unit, '(a)' ) 'set contour surface'
			elseif (contour=='both') then
				write ( file_unit, '(a)' ) 'set contour both'
			else 
				write ( file_unit, '(a)' ) 'set contour'
			end if
		end if
	end if
	write ( file_unit, '(a)' ) 'set hidden3d'
	write ( file_unit, '(a)' ) 'set parametric'
!***********************************************************************************	
	write ( file_unit, '(a)' ) 'splot "' // trim ( data_file_name ) // &
	& '" using 1:2:3 with lines palette'
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine surf_1
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine plot_4(x1,y1,x2,y2,x3,y3,x4,y4,style,pause,color1,color2,color3,color4,&
						& terminal,filename,polar,persist,input,linewidth)
!***********************************************************************************
! this subroutine plots 4 two-dimensional graphs in the same coordinate system
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: x1(:), y1(:), x2(:), y2(:), x3(:), y3(:), x4(:), y4(:)
	real(kind=4), optional		:: pause,linewidth
	character(len=*),optional	:: style, color1, color2, color3, color4, terminal, filename, polar,&
						& persist, input
	integer 			:: i, ierror, ios, file_unit, Nx1, Nx2, Nx3, Nx4, Nmax
	character(len=100)		:: data_file_name, command_file_name, my_linewidth
	integer, parameter		:: Nc=20
	character(len=Nc)		:: my_line_type1, my_line_type2, my_line_type3, my_line_type4, &
					& my_color1, my_color2, my_color3,  my_color4, my_range, my_pause, my_persist
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	Nx1=size(x1)
	Nx2=size(x2)
	Nx3=size(x3)
	Nx4=size(x4)
	if ((size(x1).ne.size(y1)).or.(size(x2).ne.size(y2)).or.(size(x3).ne.size(y3)).or.(size(x4).ne.size(y4))) then
		print *,'subroutine plot ERROR: size(x) is not equal to size(y)'
		stop
	end if
	if (present(style).and.(len(style).ne.12)) then
		print *,'subroutine plot ERROR: argument "style" has wrong number of characters'
		stop
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************
	Nmax=max(Nx1,Nx2,Nx3,Nx4)	
	do i=1,Nmax
		write (file_unit,'(8E15.7)') x1(min(i,Nx1)), y1(min(i,Nx1)), x2(min(i,Nx2)), y2(min(i,Nx2)), &
		& x3(min(i,Nx3)), y3(min(i,Nx3)), x4(min(i,Nx4)), y4(min(i,Nx4))
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_line_type1='lines'
	if (present(style)) then
	if ((style(3:3)=='-')) then
		my_line_type1='linespoints'
	else
		my_line_type1='points'
	end if
	end if
	my_line_type2='lines'
	if (present(style)) then
	if ((style(6:6)=='-')) then
		my_line_type2='linespoints'
	else
		my_line_type2='points'
	end if
	end if
	my_line_type3='lines'
	if (present(style)) then
	if ((style(9:9)=='-')) then
		my_line_type3='linespoints'
	else
		my_line_type3='points'
	end if
	end if
	my_line_type4='lines'
	if (present(style)) then
	if ((style(12:12)=='-')) then
		my_line_type4='linespoints'
	else
		my_line_type4='points'
	end if
	end if
	if (present(linewidth)) then
		write (	my_linewidth,'(e9.3)') linewidth
	else
		my_linewidth=trim(default_linewidth)
	end if	
	if (present(color1)) then
		my_color1='"'//trim(color1)//'"'
	else
		my_color1='"'//trim(default_color1)//'"'
	end if
	if (present(color2)) then
		my_color2='"'//trim(color2)//'"'
	else
		my_color2='"'//trim(default_color2)//'"'
	end if
	if (present(color3)) then
		my_color3='"'//trim(color3)//'"'
	else
		my_color3='"'//trim(default_color3)//'"'
	end if
	if (present(color4)) then
		my_color4='"'//trim(color4)//'"'
	else
		my_color4='"'//trim(default_color4)//'"'
	end if
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'unset key'
	if (present(polar).and.(polar=='yes')) then
		write (my_range,'(e15.7)') max(maxval(abs(y1)),maxval(abs(y2)),maxval(abs(y3)),maxval(abs(y4)))
		write ( file_unit, '(a)' ) 'set xrange [-'//trim(my_range)//':'//trim(my_range)//']'
		write ( file_unit, '(a)' ) 'set yrange [-'//trim(my_range)//':'//trim(my_range)//']'
		write ( file_unit, '(a)' ) 'set size square'
		write ( file_unit, '(a)' ) 'set polar'
		write ( file_unit, '(a)' ) 'set grid polar'
	else
		write ( file_unit, '(a)' ) 'set grid'
	end if
!***********************************************************************************
	if (present(style)) then
		write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
		&//'" using 1:2 with ' // trim(my_line_type1) // ' pointtype ' // &
		& style(1:2) // ' linecolor rgb ' // trim(my_color1) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		&//'" using 3:4 with ' // trim(my_line_type2) // ' pointtype ' &
		&// style(4:5) // ' linecolor rgb ' // trim(my_color2) // ' linewidth '// trim(my_linewidth) //',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		&//'" using 5:6 with ' // trim(my_line_type3) // ' pointtype ' &
		&// style(7:8) // ' linecolor rgb ' // trim(my_color3) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		&//'" using 7:8 with ' // trim(my_line_type4) // ' pointtype ' &
		&// style(10:11) // ' linecolor rgb '// trim(my_color4)// ' linewidth '// trim(my_linewidth) 
	else 
		write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
		& //'" using 1:2 with ' // trim(my_line_type1)  // ' linecolor rgb '&
		& // trim(my_color1) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		& //'" using 3:4 with ' // trim(my_line_type2)  // ' linecolor rgb '&
		& // trim(my_color2) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		& //'" using 5:6 with ' // trim(my_line_type3)  // ' linecolor rgb '&
		& // trim(my_color3) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		& //'" using 7:8 with ' // trim(my_line_type4)  // ' linecolor rgb '&
		& // trim(my_color4) // ' linewidth '// trim(my_linewidth) 
	end if
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine plot_4
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine plot_3(x1,y1,x2,y2,x3,y3,style,pause,color1,color2,color3,terminal,filename,polar,persist,input,linewidth)
!***********************************************************************************
! this subroutine plots 3 two-dimensional graphs in the same coordinate system
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: x1(:), y1(:), x2(:), y2(:), x3(:), y3(:)
	real(kind=4), optional		:: pause,linewidth
	character(len=*),optional	:: style, color1, color2, color3, terminal, filename, polar, persist, input
	integer 			:: i, ierror, ios, file_unit, Nx1, Nx2, Nx3, Nmax
	character(len=100)		:: data_file_name, command_file_name, my_linewidth
	integer, parameter		:: Nc=20
	character(len=Nc)		:: my_line_type1, my_line_type2, my_line_type3, my_color1, my_color2,&
					& my_color3, my_range, my_pause, my_persist
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	Nx1=size(x1)
	Nx2=size(x2)
	Nx3=size(x3)
	if ((size(x1).ne.size(y1)).or.(size(x2).ne.size(y2)).or.(size(x3).ne.size(y3))) then
		print *,'subroutine plot ERROR: size(x) is not equal to size(y)'
		stop
	end if
	if (present(style).and.(len(style).ne.9)) then
		print *,'subroutine plot ERROR: argument "style" has wrong number of characters'
		stop
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************
	Nmax=max(Nx1,Nx2,Nx3)	
	do i=1,Nmax
		write (file_unit,'(6E15.7)') x1(min(i,Nx1)), y1(min(i,Nx1)), x2(min(i,Nx2)), y2(min(i,Nx2)), &
		& x3(min(i,Nx3)), y3(min(i,Nx3))
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_line_type1='lines'
	if (present(style)) then
	if ((style(3:3)=='-')) then
		my_line_type1='linespoints'
	else
		my_line_type1='points'
	end if
	end if
	my_line_type2='lines'
	if (present(style)) then
	if ((style(6:6)=='-')) then
		my_line_type2='linespoints'
	else
		my_line_type2='points'
	end if
	end if
	my_line_type3='lines'
	if (present(style)) then
	if ((style(9:9)=='-')) then
		my_line_type3='linespoints'
	else
		my_line_type3='points'
	end if
	end if
	if (present(linewidth)) then
		write (	my_linewidth,'(e9.3)') linewidth
	else
		my_linewidth=trim(default_linewidth)
	end if	
	if (present(color1)) then
		my_color1='"'//trim(color1)//'"'
	else
		my_color1='"'//trim(default_color1)//'"'
	end if
	if (present(color2)) then
		my_color2='"'//trim(color2)//'"'
	else
		my_color2='"'//trim(default_color2)//'"'
	end if
	if (present(color3)) then
		my_color3='"'//trim(color3)//'"'
	else
		my_color3='"'//trim(default_color3)//'"'
	end if
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' & 
			&//trim(my_persist) // ' title  "Gnuplot"' 
	end if

!***********************************************************************************
	write ( file_unit, '(a)' ) 'unset key'
	if (present(polar).and.(polar=='yes')) then
		write (my_range,'(e15.7)') max(maxval(abs(y1)),maxval(abs(y2)),maxval(abs(y3)))
		write ( file_unit, '(a)' ) 'set xrange [-'//trim(my_range)//':'//trim(my_range)//']'
		write ( file_unit, '(a)' ) 'set yrange [-'//trim(my_range)//':'//trim(my_range)//']'
		write ( file_unit, '(a)' ) 'set size square'
		write ( file_unit, '(a)' ) 'set polar'
		write ( file_unit, '(a)' ) 'set grid polar'
	else
		write ( file_unit, '(a)' ) 'set grid'
	end if	
!***********************************************************************************
	if (present(style)) then
		write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
		&//'" using 1:2 with ' // trim(my_line_type1) // ' pointtype ' // &
		& style(1:2) // ' linecolor rgb ' // trim(my_color1) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		&//'" using 3:4 with ' // trim(my_line_type2) // ' pointtype ' &
		&// style(4:5) // ' linecolor rgb ' // trim(my_color2) // ' linewidth '// trim(my_linewidth) //',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		&//'" using 5:6 with ' // trim(my_line_type3) // ' pointtype ' &
		&// style(7:8) // ' linecolor rgb ' // trim(my_color3) // ' linewidth '// trim(my_linewidth) 
	else 
		write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
		& //'" using 1:2 with ' // trim(my_line_type1)  // ' linecolor rgb '&
		& // trim(my_color1) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		& //'" using 3:4 with ' // trim(my_line_type2)  // ' linecolor rgb '&
		& // trim(my_color2) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		& //'" using 5:6 with ' // trim(my_line_type3)  // ' linecolor rgb '&
		& // trim(my_color3) // ' linewidth '// trim(my_linewidth) 
	end if
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine plot_3
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine plot_2(x1,y1,x2,y2,style,pause,color1,color2,terminal,filename,polar,persist,input,linewidth)
!***********************************************************************************
! this subroutine plots 2 two-dimensional graphs in the same coordinate system
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: x1(:), y1(:), x2(:), y2(:)
	real(kind=4), optional		:: pause,linewidth
	character(len=*),optional	:: style, color1, color2, terminal, filename, polar, persist, input
	integer 			:: i, ierror, ios, file_unit, Nx1, Nx2, Nmax
	character(len=100)		:: data_file_name, command_file_name, my_linewidth
	integer, parameter		:: Nc=20
	character(len=Nc)		:: my_line_type1, my_line_type2, my_color1, my_color2, my_range, my_pause, my_persist
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	Nx1=size(x1)
	Nx2=size(x2)
	if ((size(x1).ne.size(y1)).or.(size(x2).ne.size(y2))) then
		print *,'subroutine plot ERROR: size(x) is not equal to size(y)'
		stop
	end if
	if (present(style).and.(len(style).ne.6)) then
		print *,'subroutine plot ERROR: argument "style" has wrong number of characters'
		stop
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************
	Nmax=max(Nx1,Nx2)	
	do i=1,Nmax
		write (file_unit,'(4E15.7)') x1(min(i,Nx1)), y1(min(i,Nx1)), x2(min(i,Nx2)), y2(min(i,Nx2))
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_line_type1='lines'
	if (present(style)) then
	if ((style(3:3)=='-')) then
		my_line_type1='linespoints'
	else
		my_line_type1='points'
	end if
	end if
	my_line_type2='lines'
	if (present(style)) then
	if ((style(6:6)=='-')) then
		my_line_type2='linespoints'
	else
		my_line_type2='points'
	end if
	end if
	if (present(linewidth)) then
		write (	my_linewidth,'(e9.3)') linewidth
	else
		my_linewidth=trim(default_linewidth)
	end if	
	if (present(color1)) then
		my_color1='"'//trim(color1)//'"'
	else
		my_color1='"'//trim(default_color1)//'"'
	end if
	if (present(color2)) then
		my_color2='"'//trim(color2)//'"'
	else
		my_color2='"'//trim(default_color2)//'"'
	end if
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"' 
	end if

!***********************************************************************************
	write ( file_unit, '(a)' ) 'unset key'
	write ( file_unit, '(a)' ) 'unset key'
	if (present(polar).and.(polar=='yes')) then
		write (my_range,'(e15.7)') max(maxval(abs(y1)),maxval(abs(y2)))
		write ( file_unit, '(a)' ) 'set xrange [-'//trim(my_range)//':'//trim(my_range)//']'
		write ( file_unit, '(a)' ) 'set yrange [-'//trim(my_range)//':'//trim(my_range)//']'
		write ( file_unit, '(a)' ) 'set size square'
		write ( file_unit, '(a)' ) 'set polar'
		write ( file_unit, '(a)' ) 'set grid polar'
	else
		write ( file_unit, '(a)' ) 'set grid'
	end if		
!***********************************************************************************
	if (present(style)) then
		write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
		&//'" using 1:2 with ' // trim(my_line_type1) // ' pointtype ' // &
		& style(1:2) // ' linecolor rgb ' // trim(my_color1) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		&//'" using 3:4 with ' // trim(my_line_type2) // ' pointtype ' &
		&// style(4:5) // ' linecolor rgb ' // trim(my_color2) // ' linewidth '// trim(my_linewidth) 
	else 
		write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
		& //'" using 1:2 with ' // trim(my_line_type1)  // ' linecolor rgb '&
		& // trim(my_color1) // ' linewidth '// trim(my_linewidth) // ',\' 
		write ( file_unit, '(a,i2,a)' ) '     "'// trim (data_file_name) &
		& //'" using 3:4 with ' // trim(my_line_type2)  // ' linecolor rgb '&
		& // trim(my_color2) // ' linewidth '// trim(my_linewidth) 
	end if
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name)
!***********************************************************************************
	end subroutine plot_2
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine plot_1(x1,y1,style,pause,color1,terminal,filename,polar,persist,input,linewidth)
!***********************************************************************************
! this subroutine plots a two-dimensional graph
!***********************************************************************************
	implicit none
	real(kind=8), intent(in)	:: x1(:), y1(:)
	real(kind=4), optional		:: pause,linewidth
	character(len=*),optional	:: style, color1, terminal, filename, polar, persist, input
	integer 			:: i, ierror, ios, file_unit, Nx1
	character(len=100)		:: data_file_name, command_file_name, my_linewidth
	integer, parameter		:: Nc=20
	character(len=Nc)		:: my_line_type1, my_color1, my_range, my_pause, my_persist
!***********************************************************************************
	if (present(input)) then
		data_file_name='data_file_'//input//'.txt'
		command_file_name='command_file_'//input//'.txt'		
	else
		data_file_name='data_file.txt'
		command_file_name='command_file.txt'
	end if
!***********************************************************************************
	Nx1=size(x1)
	if ((size(x1).ne.size(y1))) then
		print *,'subroutine plot ERROR: size(x) is not equal to size(y)'
		stop
	end if
	if (present(style).and.(len(style).ne.3)) then
		print *,'subroutine plot ERROR: argument "style" has wrong number of characters'
		stop
	end if
!***********************************************************************************
	ierror=0	
	call get_unit(file_unit)	
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=data_file_name, status='replace', iostat=ios)	
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal data file.'
		stop
	end if
!***********************************************************************************	
! here we write the date to the data_file - the gnuplot will read this data later
!***********************************************************************************	
	do i=1,Nx1
		write (file_unit,'(2E15.7)') x1(i), y1(i)
	end do
!***********************************************************************************	
	close (unit=file_unit)
!***********************************************************************************
	ierror = 0
	call get_unit(file_unit)
	if (file_unit==0) then
		ierror=1
		print *,'write_vector_data - fatal error! Could not get a free FORTRAN unit.'
		stop
	end if
	open (unit=file_unit, file=command_file_name, status='replace',	iostat=ios)
	if (ios/=0) then
		ierror=2
		print *,'write_vector_data - fatal error! Could not open the terminal command file.'
		stop
	end if	
!***********************************************************************************
! here we write the commands to the commands file which gnuplot will execute
!***********************************************************************************
	my_line_type1='lines'
	if (present(style)) then
	if ((style(3:3)=='-')) then
		my_line_type1='linespoints'
	else
		my_line_type1='points'
	end if
	end if
	if (present(linewidth)) then
		write (	my_linewidth,'(e9.3)') linewidth
	else
		my_linewidth=trim(default_linewidth)
	end if	
	if (present(color1)) then
		my_color1='"'//trim(color1)//'"'
	else
		my_color1='"'//trim(default_color1)//'"'
	end if
!***********************************************************************************
	my_persist='persist '
	if (present(persist).and.(persist=='no')) my_persist=' '
	if (present(terminal)) then
		write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal)) 
	if (present(filename)) then
		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"' 
	else
		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"' 
	end if
	else
		write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) //' title  "Gnuplot"' 
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'unset key'
	if (present(polar).and.(polar=='yes')) then
		write (my_range,'(e15.7)') maxval(abs(y1))
		write ( file_unit, '(a)' ) 'set xrange [-'//trim(my_range)//':'//trim(my_range)//']'
		write ( file_unit, '(a)' ) 'set yrange [-'//trim(my_range)//':'//trim(my_range)//']'
		write ( file_unit, '(a)' ) 'set size square'
		write ( file_unit, '(a)' ) 'set polar'
		write ( file_unit, '(a)' ) 'set grid polar'
	else
		write ( file_unit, '(a)' ) 'set grid'
	end if	
!***********************************************************************************
	if (present(style)) then
		write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
		&//'" using 1:2 with ' // trim(my_line_type1) // ' pointtype ' // &
		& style(1:2) // ' linecolor rgb ' // trim(my_color1) // ' linewidth '// trim(my_linewidth) 
	else 
		write ( file_unit, '(a,i2,a)' ) 'plot "' // trim (data_file_name) &
		& //'" using 1:2 with ' // trim(my_line_type1)  // ' linecolor rgb '&
		& // trim(my_color1) // ' linewidth '// trim(my_linewidth)  
	end if
!***********************************************************************************
	if (present(pause)) then
		if (pause<0.0) then
			write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
		else 
			write (	my_pause,'(e9.3)') pause
			write ( file_unit, '(a)' ) 'pause ' // trim(my_pause) 
		end if
	else
		write ( file_unit, '(a)' ) 'pause 0'
	end if
!***********************************************************************************
	write ( file_unit, '(a)' ) 'q'
	close ( unit = file_unit )
!***********************************************************************************
	call run_gnuplot (command_file_name) 
!***********************************************************************************
	end subroutine plot_1
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine run_gnuplot(command_file_name)
!***********************************************************************************
	implicit none
	character (len = 100) command
	character (len = *) command_file_name
	integer status
	integer system
!***********************************************************************************
!  Issue a command to the system that will startup GNUPLOT, using
!  the file we just wrote as input.
!***********************************************************************************
	write (command, *) 'gnuplot ' // trim (command_file_name)		
	status=system(trim(command))	
	if (status.ne.0) then
		print *,'RUN_GNUPLOT - Fatal error!'
		stop
	end if	
	return
!***********************************************************************************
	end subroutine run_gnuplot
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	subroutine get_unit(iunit)
!***********************************************************************************
	implicit none
	integer i
	integer ios
	integer iunit
	logical lopen
!***********************************************************************************	
	iunit=0
	do i=1,99
		if (i/= 5 .and. i/=6) then	
			inquire (unit=i, opened=lopen, iostat=ios)
			if (ios==0) then
				if (.not.lopen) then
					iunit=i
					return
				end if
			end if
		
		end if
	end do	
	return
	end subroutine get_unit
!***********************************************************************************
	end module gnufor2