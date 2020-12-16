subroutine r8_pardump (num, value)
  integer,intent(in) :: num ! channel to write
  real(r8),intent(in) :: value (lbx:ubx, lby:uby, lbz:ubz)

  integer n
  real(r8) tval(nx,ny,nz) 

  if (iomaster) then ! iomaster -- process that writes to disk
    tval = value(1:nx, 1:ny, 1:nz)
    write (num) tval
    do n = me_1d+1, me_1d+chunksize-1
      call mpi_recv &
        (tval, nx*ny*nz, mpi_real8, rank(n), 0, mycomm, status, dum)
          ! rank(n) -- vector of  
      write (num) tval
    end do
  else ! process that does not write to disk
    tval = value(1:nx, 1:ny, 1:nz)
    call mpi_send (tval, nx*ny*nz, mpi_real8, io_rank, 0, mycomm, dum)
     ! they send data to the writing process
  endif
end subroutine