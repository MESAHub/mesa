      function timer ( )

c*********************************************************************72
c
cc TIMER computes the elapsed execution time.
c
c  Modified:
c
c    10 December 2007
c
      real result
      real tarray(2)
      real timer

      call etime(tarray, result )

      timer = result

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end

      subroutine y12mae(n,z,a,snr,nn,rnr,nn1,pivot,ha,iha,aflag,iflag,
     &  b,ifail)

c*********************************************************************72
c
cc Y12MAE solves a single sparse linear system with one right hand side.
c
      integer iha
      integer n
      integer nn
      integer nn1

      real a(nn)
      real aflag(8)
      real b(n)
      integer ha(iha,11)
      integer ifail
      integer iflag(10)
      real pivot(n)
      integer rnr(nn1)
      integer snr(nn)
      integer z

      aflag(1)=16.0
      aflag(2)=1.0e-12 ! THR threshold
      aflag(3)=1.0e+16
      aflag(4)=1.0e-12

      iflag(2)=2
      iflag(3)=1
      iflag(4)=0
      iflag(5)=1

      call y12mbe(n,z,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)

      if(ifail.ne.0) then
        write(*,*)' '
        write(*,*)'Y12MAE - Fatal error!'
        write(*,*)'  Nonzero value of IFAIL returned by Y12MBE.'
        return
      end if

      call y12mce(n,z,a,snr,nn,rnr,nn1,pivot,b,ha,iha,aflag,iflag,
     &  ifail)

      if(ifail.ne.0) then
        write(*,*)' '
        write(*,*)'Y12MAE - Fatal error!'
        write(*,*)'  Nonzero value of IFAIL returned by Y12MCE.'
        return
      end if

      call y12mde(n,a,nn,b,pivot,snr,ha,iha,iflag,ifail)

      if(ifail.ne.0) then
        write(*,*)' '
        write(*,*)'Y12MAE - Fatal error!'
        write(*,*)'  Nonzero value of IFAIL returned by Y12MDE.'
        return
      end if

      return
      end

      subroutine y12maf(n,z,a,snr,nn,rnr,nn1,pivot,ha,iha,aflag,iflag,
     &  b,ifail)

c*********************************************************************72
c
cc Y12MAF solves a single sparse linear system with one right hand side.
c
      integer iha
      integer n
      integer nn
      integer nn1

      double precision a(nn)
      double precision aflag(8)
      double precision b(n)
      integer ha(iha,11)
      integer ifail
      integer iflag(10)
      double precision pivot(n)
      integer rnr(nn1)
      integer snr(nn)
      integer z

      aflag(1)=16.0d0
      aflag(2)=1.0d-12  ! THR threshold
      aflag(3)=1.0d+16
      aflag(4)=1.0d-12
      iflag(2)=2
      iflag(3)=1
      iflag(4)=0
      iflag(5)=1

      call y12mbf(n,z,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)

      if(ifail.ne.0) then
        write(*,*)' '
        write(*,*)'Y12MAF - Fatal error!'
        write(*,*)'  Nonzero value of IFAIL returned by Y12MBF.'
        return
      end if

      call y12mcf(n,z,a,snr,nn,rnr,nn1,pivot,b,ha,iha,aflag,iflag,
     &  ifail)

      if(ifail.ne.0) then
        write(*,*)' '
        write(*,*)'Y12MAF - Fatal error!'
        write(*,*)'  Nonzero value of IFAIL returned by Y12MCF.'
        return
      end if

      call y12mdf(n,a,nn,b,pivot,snr,ha,iha,iflag,ifail)

      if(ifail.ne.0) then
        write(*,*)' '
        write(*,*)'Y12MAF - Fatal error!'
        write(*,*)'  Nonzero value of IFAIL returned by Y12MDF.'
        return
      end if

      return
      end

      subroutine y12mbe(n,z,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)

c*********************************************************************72
c
cc Y12MBE prepares a sparse linear system to be factored and solved.
c
c  Discussion:
c
c    the non-zero elements of a sparse matrix a are prepared  in order to
c    solve the system ax=b by use of sparse matrix techniques.
c
      integer iha
      integer nn
      integer nn1

      real a(nn)
      real aflag(8)
      real gt1
      integer ha(iha,11)
      integer i
      integer ifail
      integer iflag(10)
      integer index
      integer j
      integer l1
      integer l2
      integer l3
      integer l4
      integer l5
      integer mode
      integer n
      integer r
      integer rnr(nn1)
      integer snr(nn)
      real t
      integer z

      mode=iflag(4)
      ifail=0

      if(n.lt.2) then
        ifail=12
        return
      end if

      if(z.le.0) then
        ifail=13
        return
      end if

      if(nn.lt.2*z) then
        ifail=5
        write(*,*)' '
        write(*,*)'Y12MBE - Fatal error!'
        write(*,*)'  NN, the declared dimension of A and SNR'
        write(*,*)'  is too small.'
        write(*,*)'  NN must be at least 2*Z.'
        write(*,*)'  NN = ',nn
        write(*,*)'  Z = ',z
        return
      end if

      if(nn1.lt.z) then
        ifail=6
        return
      end if

      if(ifail.eq.0.and.n.gt.z) then
        ifail=14
        return
      end if

      if(iha.lt.n) then
        ifail=15
        return
      end if

      if(mode.lt.0) then
        ifail=16
        return
      end if

      if(mode.gt.2) then
        ifail=16
        return
      end if

      gt1=0.0

      do i=1,n
        ha(i,2)=0
        ha(i,3)=0
        ha(i,6)=0
      end do
c
c  Find the number of the non-zero elements in each row and column;
c  Move the non-zero elements in the end of the arrays a and snr;find the
c  largest non-zero element in a(in absolute value).
c
      do i=1,z

        t=abs(a(i))
        l3=rnr(i)
        l4=snr(i)

        if(l4.gt.n.or.l4.lt.1) then
          ifail=24
          write(*,*)' '
          write(*,*)'Y12MBE - Fatal error!'
          write(*,*)'  IFAIL = ',ifail
          write(*,*)'  Column number I=',i,' is illegal.'
          write(*,*)'  SNR(I)=',snr(i)
          return
        end if

        if(l3.gt.n.or.l3.lt.1) then
          ifail=25
          write(*,*)' '
          write(*,*)'Y12MBE - Fatal error!'
          write(*,*)'  IFAIL = ',ifail
          write(*,*)'  Row number I=',i,' is illegal.'
          write(*,*)'  RNR(I)=',rnr(i)
          return
        end if


        ha(l3,3)=ha(l3,3)+1
        ha(l4,6)=ha(l4,6)+1
        if(t.gt.gt1)gt1=t
        a(z+i)=a(i)
        snr(z+i)=snr(i)

      end do

      if(ifail.gt.0)return
c
c  Store the information of the row starts in HA(I,1) and of the column
c  starts in HA(I,4).
c
      l1=1
      l2=1

      do i=1,n

        l3=ha(i,3)
        l4=ha(i,6)

        if(l3.le.0) then
          ifail=17
          return
        end if

        if(l4.le.0) then
          ifail=18
          return
        end if

        if(mode.ne.2) then
          ha(i,9)=l3
          ha(i,10)=l4
          ha(i,11)=0
          ha(l3,2)=ha(l3,2)+1
         ha(i,5)=l3
        end if

        ha(i,1)=l1
        ha(i,4)=l2
        l1=l1+l3
        l2=l2+l4
        ha(i,3)=0
        ha(i,6)=0

      end do
c
c  Store the non-zero elements of matrix A (ordered in rows) in the
c  first Z locations of the array A.
c
c  Do the same for their column numbers
c
      do i=1,z
        l1=z+i
        l3=rnr(i)
        l2=ha(l3,1)+ha(l3,3)
        a(l2)=a(l1)
        snr(l2)=snr(l1)
        ha(l3,3)=ha(l3,3)+1
      end do
c
c  Store the row numbers of the non-zero elements ordered by columns in
c  the first z locations of the array rnr. store information about row
c  ends(in ha(i,3)).
c
      l4=1
      do i=1,n

        if(mode.ne.2) then
          if(ha(i,2).ne.0) then
            ha(i,11)=l4
            l4=l4+ha(i,2)
            ha(i,2)=ha(i,11)
          end if
        end if

        ha(i,3)=ha(i,1)+ha(i,3)-1
        l1=ha(i,1)
        l2=ha(i,3)

        do j=l1,l2

          l3=snr(j)
          r=ha(l3,6)
          index=ha(l3,4)+r
          rnr(index)=i

          if(r.ne.0) then
            if(j.ne.l1) then
              if(rnr(index-1).eq.i) then
                ifail=11
                return
              end if
            end if
          end if

          ha(l3,6)=r+1

        end do

      end do

      do i=1,n

        if(mode.ne.2) then
          l3=ha(i,5)
          l5=ha(l3,2)
          ha(l5,8)=i
          ha(i,7)=l5
          ha(l3,2)=ha(l3,2)+1
        end if

        ha(i,6)=ha(i,4)+ha(i,6)-1

      end do

      aflag(6)=gt1
      iflag(6)=0
      iflag(7)=0
      iflag(8)=z
      iflag(1)=-1

      return
      end

      subroutine y12mbf(n,z,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)

c*********************************************************************72
c
cc Y12MBF prepares a sparse linear system to be factored and solved.
c
c  Discussion:
c
c    the non-zero elements of a sparse matrix a are prepared  in order to
c    solve the system ax=b by use of sparse matrix techniques.
c
      integer iha
      integer nn
      integer nn1

      double precision a(nn)
      double precision aflag(8)
      double precision gt1
      integer ha(iha,11)
      integer i
      integer ifail
      integer iflag(10)
      integer index
      integer j
      integer l1
      integer l2
      integer l3
      integer l4
      integer l5
      integer mode
      integer n
      integer r
      integer rnr(nn1)
      integer snr(nn)
      double precision t
      integer z

      mode=iflag(4)
      ifail=0

      if(n.lt.2) then
        ifail=12
        return
      end if

      if(z.le.0) then
        ifail=13
        return
      end if

      if(nn.lt.2*z) then
        ifail=5
        write(*,*)' '
        write(*,*)'Y12MBF - Fatal error!'
        write(*,*)'  NN, the declared dimension of A and SNR'
        write(*,*)'  is too small.'
        write(*,*)'  NN must be at least 2*Z.'
        write(*,*)'  NN = ',nn
        write(*,*)'  Z = ',z
        return
      end if

      if(nn1.lt.z) then
        ifail=6
        return
      end if

      if(ifail.eq.0.and.n.gt.z) then
        ifail=14
        return
      end if

      if(iha.lt.n) then
        ifail=15
        return
      end if

      if(mode.lt.0) then
        ifail=16
        return
      end if

      if(mode.gt.2) then
        ifail=16
        return
      end if

      gt1=0.0

      do i=1,n
        ha(i,2)=0
        ha(i,3)=0
        ha(i,6)=0
      end do
c
c  find the number of the non-zero elements in each row and column;move
c  the non-zero elements in the end of the arrays a and snr;find the
c  largest non-zero element in a(in absolute value).
c
      do i=1,z
        t=abs(a(i))
        l3=rnr(i)
        l4=snr(i)

        if(l4.gt.n.or.l4.lt.1) then
          ifail=24
          write(*,*)' '
          write(*,*)'Y12MBF - Fatal error!'
          write(*,*)'  IFAIL = ',ifail
          write(*,*)'  Column number I=',i,' is illegal.'
          write(*,*)'  SNR(I)=',snr(i)
          return
        end if

        if(l3.gt.n.or.l3.lt.1) then
          ifail=25
          write(*,*)' '
          write(*,*)'Y12MBF - Fatal error!'
          write(*,*)'  IFAIL = ',ifail
          write(*,*)'  Row number I=',i,' is illegal.'
          write(*,*)'  RNR(I)=',rnr(i)
          return
        end if

        ha(l3,3)=ha(l3,3)+1
        ha(l4,6)=ha(l4,6)+1
        if(t.gt.gt1)gt1=t
        a(z+i)=a(i)
        snr(z+i)=snr(i)
      end do

      if(ifail.gt.0)return
c
c  store the information of the row starts(in ha(i,1))and of the column
c  starts(in ha(i,4)).
c
      l1=1
      l2=1

      do i=1,n

        l3=ha(i,3)
        l4=ha(i,6)

        if(l3.le.0) then
          ifail=17
          return
        end if

        if(l4.le.0) then
          ifail=18
          return
        end if

        if(mode.ne.2) then
          ha(i,9)=l3
          ha(i,10)=l4
          ha(i,11)=0
          ha(l3,2)=ha(l3,2)+1
         ha(i,5)=l3
        end if

        ha(i,1)=l1
        ha(i,4)=l2
        l1=l1+l3
        l2=l2+l4
        ha(i,3)=0
        ha(i,6)=0

      end do
c
c  Store the non-zero elements of matrix a(ordered in rows) in the
c  first z locations of the array a.
c
c  Do the same for their column numbers
c
      do i=1,z
        l1=z+i
        l3=rnr(i)
        l2=ha(l3,1)+ha(l3,3)
        a(l2)=a(l1)
        snr(l2)=snr(l1)
        ha(l3,3)=ha(l3,3)+1
      end do
c
c  Store the row numbers of the non-zero elements ordered by columns in
c  the first z locations of the array rnr. store information about row
c  ends(in ha(i,3)).
c
      l4=1
      do i=1,n

        if(mode.ne.2) then
          if(ha(i,2).ne.0) then
            ha(i,11)=l4
            l4=l4+ha(i,2)
            ha(i,2)=ha(i,11)
          end if
        end if

        ha(i,3)=ha(i,1)+ha(i,3)-1
        l1=ha(i,1)
        l2=ha(i,3)

        do j=l1,l2

          l3=snr(j)
          r=ha(l3,6)
          index=ha(l3,4)+r
          rnr(index)=i

          if(r.ne.0) then
            if(j.ne.l1) then
              if(rnr(index-1).eq.i) then
                ifail=11
                return
              end if
            end if
          end if

          ha(l3,6)=r+1

        end do

      end do

      do i=1,n

        if(mode.ne.2) then
          l3=ha(i,5)
          l5=ha(l3,2)
          ha(l5,8)=i
          ha(i,7)=l5
          ha(l3,2)=ha(l3,2)+1
        end if

        ha(i,6)=ha(i,4)+ha(i,6)-1

      end do

      aflag(6)=gt1
      iflag(6)=0
      iflag(7)=0
      iflag(8)=z
      iflag(1)=-1

      return
      end

      subroutine y12mce(n,z,a,snr,nn,rnr,nn1,pivot,b,ha,iha,aflag,
     &  iflag,ifail)

c*********************************************************************72
c
cc Y12MCE finds LU factors of a sparse matrix preprocessed by Y12MBE.
c
      integer iha
      integer n
      integer nn
      integer nn1

      real a(nn)
      real aflag(8)
      real b(n)
      integer c1
      integer c2
      integer cr1
      integer cr2
      integer cr3
      integer cr4
      real grmin
      integer ha(iha,11)
      integer i
      integer i1
      integer ifail
      integer iflag(10)
      integer index
      integer j
      integer jj
      integer k
      integer kk
      integer l
      integer l1
      integer l2
      integer l3
      integer l4
      integer l5
      integer l6
      integer l7
      integer lfc
      integer lfr
      integer ll
      integer n7
      integer n8
      integer nr
      real pivot(n)
      integer r
      integer r1
      integer r10
      integer r2
      integer r3
      integer r4
      integer r5
      integer r6
      integer r7
      integer r8
      integer r9
      integer rcoll
      integer rpivot
      integer rr
      integer rr1
      integer rr2
      integer rr3
      integer rr4
      integer rrow
      integer rnr(nn1)
      integer slut
      integer snr(nn)
      real t
      real td
      real td1
      real tol1
      real tol2
      real tol3
      real u
      real v
      integer z
      integer zz
c
c  Information which is necessary to begin the elimination is stored.
c
      ifail=0

      if(iflag(1).ne.-1)ifail=2

      if(aflag(1).lt.1.0) then
        aflag(1)=1.0005
      end if

      if(aflag(3).lt.1.0e+5) then
        aflag(3)=1.0e+5
      end if

      if(aflag(4).lt.0.0) then
        aflag(4)=-aflag(4)
      end if

      if(iflag(2).lt.1)ifail=19

      if(iflag(3).lt.0.or.iflag(3).gt.2)ifail=20

      if(iflag(5).lt.1.or.iflag(5).gt.3) then
        ifail=21
        write(*,*)' '
        write(*,*)'Y12MCE - Fatal error!'
        write(*,*)'  Input value of IFLAG(5)=',iflag(5)
        write(*,*)'  Legal values are 1, 2, or 3.'
        return
      end if

      if(iflag(5).eq.3)ifail=22

      if(ifail.gt.0)return

      snr(z+1)=0
      rnr(z+1)=0
      n8=n+1
      n7=n-1
      u=aflag(1)
      grmin=aflag(4)*aflag(6)
      zz=z
c
c  Use the information about fill-ins if it is possible.
c
      nr=n*n
      if(iflag(4).ne.2)go to 100

      if(iflag(10).le.nn) then

        l1=iflag(10)
        l5=l1+1
        if(l5.le.nn)snr(l5)=0

        do i=1,n

          l=n8-i
          l2=ha(l,3)+1
          l3=l2-ha(l,1)

          do j=1,l3
            snr(l5-j)=snr(l2-j)
            a(l5-j)=a(l2-j)
          end do

          ha(l,3)=l1
          ha(l,1)=l5-l3
          l6=l1-l3
          l5=l5-ha(l,9)

          do j=l5,l6
            snr(j)=0
          end do

          l1=l5-1

        end do

      end if

      if(iflag(9).le.nn1) then

        l2=iflag(9)
        l5=l2+1
        if(l5.le.nn1)rnr(l5)=0

        do i=1,n
          l=n8-i
          l1=ha(l,6)+1

          l4=l1-ha(l,4)
          do j=1,l4
            rnr(l5-j)=rnr(l1-j)
          end do

          ha(l,4)=l5-l4
          ha(l,6)=l2

          l6=l2-l4
          l5=l5-ha(l,10)
          do j=l5,l6
            rnr(j)=0
          end do

          l2=l5-1
        end do

      end if

  100 continue

      r4=ha(n,3)
      r5=ha(n,6)
      aflag(7)=aflag(6)
      aflag(8)=aflag(6)

      do i=1,n
        pivot(i)=0.0
        ha(i,2)=ha(i,1)
        ha(i,5)=ha(i,4)
      end do

      index=ha(n,8)
      slut=ha(index,3)-ha(index,2)+1
c
c  Start of gaussian elimination.
c
      do 950 i=1,n7

      rr3=ha(i,2)
      rr4=ha(i,3)
      c1=ha(i,4)
      cr4=ha(i,6)
      if(iflag(3).eq.0)go to 350

      if(iflag(4).eq.2) then
        rrow=ha(i,7)
        rcoll=ha(i,8)
        go to 220
      end if

      l4=ha(i,8)

      if(iflag(3).ne.1) then
        rrow=l4
        rcoll=rrow
        rpivot=i
        go to 170
      end if

      r=nr
      v=0.0
      index=iflag(2)

      do kk=1,index

        l1=i-1+kk
        if(l1.gt.n)go to 170
        j=ha(l1,8)
        r7=ha(j,2)
        r8=ha(j,3)
        r9=r8-r7

        t=0.0
        do k=r7,r8
          td=abs(a(k))
          if(t.lt.td)t=td
        end do

        t=t/u

        do k=r7,r8

          td=abs(a(k))
          if(td.lt.t)go to 150
          r6=snr(k)
          r3=r9*(ha(r6,6)-ha(r6,5))
          if(r3.gt.r)go to 150
          if(r3.lt.r)go to 151
          if(v.ge.td)go to 150
  151     v=td
          rrow=j
          rcoll=r6
          r=r3
          rpivot=l1

  150     continue

        end do

      end do

  170 r3=ha(rcoll,10)
      ha(rcoll,10)=ha(i,10)
      ha(i,10)=r3
      r3=ha(rrow,9)
      ha(rrow,9)=ha(i,9)
      ha(i,9)=r3
c
c  remove the pivot row of the list where the rows are ordered by
c  increasing numbers of non-zero elements.
c
      l1=0
      l=i
      l2=ha(l4,3)-ha(l4,2)+1
  180 l=l+1
      if(l2.gt.l1)ha(l2,11)=l
      if(l.gt.n)go to 190
      l5=ha(l,8)
      l3=ha(l5,3)-ha(l5,2)+1

      if(rpivot.ge.l) then
        ha(l4,7)=l
        ha(l,8)=l4
        l4=l5
        l1=l2
        l2=l3
        l3=n8
        go to 180
      end if

  190 if(l2.eq.l1)go to 200
      if(l3.eq.l2)go to 200
      ha(l2,11)=0
  200 l5=ha(i,7)
      if(rrow.eq.i)go to 210
      ha(l5,8)=rrow
      ha(rrow,7)=l5
  210 ha(i,7)=rrow
      ha(i,8)=rcoll
c
c  Row interchanges.
c
  220 continue

      if(rrow.eq.i)go to 290
      t=b(rrow)
      b(rrow)=b(i)
      b(i)=t

      do j=rr3,rr4
        l1=snr(j)
        r=ha(l1,5)-1
        r10=ha(l1,6)

  240   continue
        r=r+1
        if(rnr(r).ne.i)go to 240

        rnr(r)=rnr(r10)
        rnr(r10)=rrow
      end do

      rr3=ha(rrow,2)
      rr4=ha(rrow,3)

      do j=rr3,rr4
        l1=snr(j)
        r=ha(l1,5)-1

  260   continue
        r=r+1
        if(rnr(r).ne.rrow)go to 260

        rnr(r)=i
      end do

      do j=1,3
        r3=ha(rrow,j)
        ha(rrow,j)=ha(i,j)
        ha(i,j)=r3
      end do
c
c  column interchanges.
c
  290 if(rcoll.eq.i)go to 350

      do j=c1,cr4
        l1=rnr(j)
        r=ha(l1,2)-1
        r10=ha(l1,3)
  300   r=r+1
        if(snr(r).ne.i)go to 300
        t=a(r10)
        a(r10)=a(r)
        a(r)=t
        snr(r)=snr(r10)
        snr(r10)=rcoll
      end do

      c1=ha(rcoll,4)
      cr4=ha(rcoll,6)

      do j=c1,cr4
        l1=rnr(j)
        r=ha(l1,2)-1

  320   continue
        r=r+1
        if(snr(r).ne.rcoll)go to 320

        snr(r)=i
      end do

      do j=4,6
        r3=ha(rcoll,j)
        ha(rcoll,j)=ha(i,j)
        ha(i,j)=r3
      end do
c
c end of the interchanges.
c the row ordered list and the column ordered list are prepared to
c begin step i of the elimination.
c
  350 r9=rr4-rr3

      do rr=rr3,rr4
        if(snr(rr).eq.i)go to 370
      end do

      ifail=9
      go to 1110

  370 v=a(rr)
      pivot(i)=v
      td=abs(v)
      if(td.lt.aflag(8))aflag(8)=td

      if(td.lt.grmin) then
        ifail=3
        go to 1110
      end if

      r2=ha(i,1)
      a(rr)=a(rr3)
      snr(rr)=snr(rr3)
      a(rr3)=a(r2)
      snr(rr3)=snr(r2)
      snr(r2)=0
      z=z-1
      rr3=rr3+1
      ha(i,2)=rr3
      ha(i,1)=r2+1
      cr3=ha(i,5)
      if(r9.le.0)go to 431

      do j=rr3,rr4
        index=snr(j)
        pivot(index)=a(j)
      end do

  431 r7=cr4-cr3+1
      do 880 k=1,r7
      r1=rnr(cr3-1+k)
      if(r1.eq.i)go to 870
      i1=ha(r1,1)
      rr1=ha(r1,2)
      rr2=ha(r1,3)
      l2=rr2-rr1+1
      l=rr1-1
  390 l=l+1
      if(snr(l).ne.i)go to 390
      t=a(l)/v
      if(iflag(5).eq.2)go to 400
      a(l)=a(i1)
      snr(l)=snr(i1)
      snr(i1)=0
      i1=i1+1
      ha(r1,1)=i1
      z=z-1
      go to 410
  400 a(l)=a(rr1)
      a(rr1)=t
      r3=snr(rr1)
      snr(rr1)=snr(l)
      snr(l)=r3
  410 rr1=rr1+1
      ha(r1,2)=rr1
      b(r1)=b(r1)-b(i)*t
      if(r9.le.0)go to 669
      r=rr1

      do l=r,rr2
        l1=snr(l)
        td=pivot(l1)
        if(td.eq.0.0)go to 450
        pivot(l1)=0.0
        td=a(l)-td*t
        a(l)=td
        td1=abs(td)
        if(td1.gt.aflag(7))aflag(7)=td1
        if(td1.gt.aflag(2))go to 450   ! THR threshold
c
c  Too small element is created.  Remove it from the lists.
c
        z=z-1
        a(l)=a(rr1)
        snr(l)=snr(rr1)
        a(rr1)=a(i1)
        snr(rr1)=snr(i1)
        snr(i1)=0
        rr1=rr1+1
        i1=i1+1
        ha(r1,2)=rr1
        ha(r1,1)=i1
        r3=ha(l1,5)
        r2=r3-1
        l4=ha(l1,4)
        l5=rnr(l4)
        l6=rnr(r3)
  440   r2=r2+1
        if(rnr(r2).ne.r1)go to 440
        rnr(r2)=l6
        rnr(r3)=l5
        rnr(l4)=0
        ha(l1,5)=r3+1
        ha(l1,4)=l4+1
  450   continue
      end do

      do 750 j=1,r9
      r=rr3-1+j
      r2=snr(r)
      tol2=pivot(r2)
      pivot(r2)=a(r)
      if(tol2.eq.0.0)go to 740
      tol3=-tol2*t
      tol1=abs(tol3)
      if(tol1.lt.aflag(2))go to 740  ! THR threshold
      c2=ha(r2,4)
      cr2=ha(r2,6)
      cr1=ha(r2,5)
      lfr=rr2-i1+2
      lfc=cr2-c2+2

      if(iflag(4).eq.1) then
        if(lfr.gt.ha(r1,9))ha(r1,9)=lfr
        if(lfc.gt.ha(r2,10))ha(r2,10)=lfc
      end if

      if(i1.eq.1)go to 490
      if(snr(i1-1).eq.0)go to 600
  490 if(rr2.eq.nn)go to 500
      if(snr(rr2+1).eq.0)go to 580
  500 r10=nn-lfr
c
c  collection in row ordered list.
c
      if(r10.ge.r4)go to 560

      iflag(6)=iflag(6)+1

      do jj=1,n
        l1=ha(jj,3)

        if(l1.ge.ha(jj,1)) then
          ha(jj,3)=snr(l1)
          snr(l1)=-jj
        end if

      end do

      l3=0
      l4=1

      do jj=1,r4

        if(snr(jj).ne.0) then

          l3=l3+1

          if(snr(jj).le.0) then
            l5=-snr(jj)
            snr(jj)=ha(l5,3)
            ha(l5,3)=l3
            l6=l4+ha(l5,2)-ha(l5,1)
            ha(l5,2)=l6
            ha(l5,1)=l4
            l4=l3+1
          end if

          a(l3)=a(jj)
          snr(l3)=snr(jj)

        end if

      end do

      r4=l3
      snr(l3+1)=0
      rr3=ha(i,2)
      rr4=ha(i,3)
      i1=ha(r1,1)
      rr1=ha(r1,2)
      r=rr3-1+j

      if(r10.lt.r4) then
        ifail=5
        write(*,*)' '
        write(*,*)'Y12MCE - Fatal error!'
        write(*,*)'  NN, the declared dimension of A and SNR,'
        write(*,*)'  is too small.'
        write(*,*)'  Current value is NN = ',nn
        return
      end if
c
c  Fill-in takes place in the row ordered list.
c
560   continue

      r8=lfr-1
      rr2=r4+lfr

      l3=i1-1

      do ll=1,r8
        l4=r4+ll
        l5=l3+ll
        a(l4)=a(l5)
        snr(l4)=snr(l5)
        snr(l5)=0
      end do

      rr1=r4+rr1-i1+1
      ha(r1,3)=rr2
      ha(r1,2)=rr1
      i1=r4+1
      ha(r1,1)=i1
      l1=rr2
      go to 590

  580 rr2=rr2+1
      ha(r1,3)=rr2
      l1=rr2
      if(rr2.le.r4)go to 610
  590 r4=rr2
      if(r4.lt.nn)snr(r4+1)=0
      go to 610
  600 rr1=rr1-1
      i1=i1-1
      ha(r1,1)=i1
      ha(r1,2)=rr1
      l1=rr1
      snr(i1)=snr(l1)
      a(i1)=a(l1)
  610 a(l1)=tol3
      snr(l1)=snr(r)
      td=abs(a(l1))
      if(td.gt.aflag(7))aflag(7)=td
      z=z+1
      if(iflag(8).lt.z)      iflag(8)=z
      if(c2.eq.1)go to 620
      if(rnr(c2-1).eq.0)go to 720
  620 if(cr2.eq.nn1)go to 630
      if(rnr(cr2+1).eq.0)go to 700
  630 r10=nn1-lfc
c
c  collection in column ordered list.
c
      if(r10.ge.r5)go to 680
      iflag(7)=iflag(7)+1

      do jj=i,n
        l1=ha(jj,6)
        ha(jj,6)=rnr(l1)
        rnr(l1)=-jj
      end do

      l3=0
      l4=1

      do jj=1,r5

        if(rnr(jj).ne.0) then

          l3=l3+1

          if(rnr(jj).le.0) then
            l5=-rnr(jj)
            rnr(jj)=ha(l5,6)
            ha(l5,6)=l3
            l6=l4+ha(l5,5)-ha(l5,4)
            ha(l5,5)=l6
            ha(l5,4)=l4
            l4=l3+1
          end if

          rnr(l3)=rnr(jj)

        end if

      end do

      r5=l3
      rnr(r5+1)=0
      c2=ha(r2,4)
      cr3=ha(i,5)
      cr4=ha(i,6)
      cr1=ha(r2,5)

      if(r10.lt.r5) then
        ifail=6
        z=zz
        return
      end if
c
c  Fill-in takes place in the column ordered list.
c
680   continue

      r8=lfc-1
      cr2=r5+lfc

      l3=c2-1

      do l=1,r8
        l4=r5+l
        l5=l3+l
        rnr(l4)=rnr(l5)
        rnr(l5)=0
      end do

      cr1=r5+cr1-c2+1
      c2=r5+1
      ha(r2,6)=cr2
      ha(r2,4)=c2
      ha(r2,5)=cr1
      r=cr2
      go to 710

  700 cr2=cr2+1
      ha(r2,6)=cr2
      r=cr2
      if(cr2.le.r5)go to 730
  710 r5=cr2
      if(r5.lt.nn1)rnr(r5+1)=0
      go to 730
  720 cr1=cr1-1
      c2=c2-1
      ha(r2,4)=c2
      ha(r2,5)=cr1
      r=cr1
      rnr(c2)=rnr(r)
  730 rnr(r)=r1
  740 continue
  750 continue
  669 if(rr1.le.rr2)go to 760
      ifail=7
      go to 1110
c
c  update the information in the list where the rows are ordered by
c  increasing numbers of the non-zero elements.
c
  760 if(iflag(4).eq.2)go to 870
      if(iflag(3).eq.0)go to 870
      l1=rr2-rr1+1
      if(l1.eq.l2)go to 870
      l6=ha(r1,7)
      l4=ha(l2,11)
      if(l1.gt.l2)go to 820
      if(l6.gt.l4)go to 780
      if(l4.eq.n)go to 770
      l=ha(l4+1,8)
      l5=ha(l,3)-ha(l,2)+1
      if(l5.eq.l2)go to 790
  770 ha(l2,11)=0
      go to 800
  780 l5=ha(l4,8)
      l3=ha(l6,8)
      ha(l4,8)=l3
      ha(l6,8)=l5
      ha(l5,7)=l6
      ha(l3,7)=l4
      l6=l4
  790 ha(l2,11)=l4+1
  800 if(l4.eq.i+1)go to 810
      l=ha(l6-1,8)
      l2=ha(l,3)-ha(l,2)+1
      l4=ha(l2,11)
      if(l1.lt.l2)go to 780
  810 if(l1.ne.l2)ha(l1,11)=l6
      go to 870
  820 if(l6.gt.l4)go to 840
      if(l4.eq.n)go to 830
      l=ha(l4+1,8)
      l5=ha(l,3)-ha(l,2)+1
      if(l5.eq.l2)go to 840
  830 ha(l2,11)=0
  840 l2=l2+1
      if(l2.le.slut)go to 850
      l3=n
      slut=l1
      l2=l1
      go to 860
  850 l3=ha(l2,11)-1
      if(l3.eq.-1)go to 840
      if(l2.gt.l1)l2=l1
  860 ha(l2,11)=l3
      l4=ha(l3,8)
      l7=ha(l6,8)
      ha(l3,8)=l7
      ha(l6,8)=l4
      ha(l7,7)=l3
      ha(l4,7)=l6
      l6=l3
      if(l2.lt.l1)go to 840
  870 continue
  880 continue

      if(r9.gt.0) then
        do j=rr3,rr4
          index=snr(j)
          pivot(index)=0.0
        end do
      end if

      cr3=ha(i,4)

      do j=cr3,cr4
        rnr(j)=0
      end do

      l2=ha(i,2)-1

      do ll=1,r9

        r=snr(l2+ll)
        r1=ha(r,5)
        r2=ha(r,6)

        if(r2.le.r1) then
          ifail=8
          z=zz
          return
        end if

        ha(r,5)=r1+1
        r3=r1-1

  910   r3=r3+1
        if(rnr(r3).ne.i)go to 910

        rnr(r3)=rnr(r1)
        rnr(r1)=i

      end do

      aflag(5)=aflag(7)/aflag(6)

      if(aflag(5).ge.aflag(3)) then
        ifail=4
        z=zz
        return
      end if

  950 continue
c
c  preparation to begin the back substitution.
c
      index=ha(n,2)
      pivot(n)=a(index)
      a(index)=0.0
      td=abs(pivot(n))
      if(td.gt.aflag(7))aflag(7)=td
      if(td.lt.aflag(8))aflag(8)=td

      if(td.le.grmin) then
        ifail=3
        z=zz
        return
      end if

      if(iflag(4).eq.1) then

        iflag(10)=ha(n,9)
        iflag(9)=ha(n,10)

        do i=1,n7

          r1=n-i
          iflag(10)=iflag(10)+ha(r1,9)
          iflag(9)=iflag(9)+ha(r1,10)

          if(iflag(3).ne.0) then

            do j=9,10
              r2=ha(r1,j-2)
              r6=ha(r2,j)
              ha(r2,j)=ha(r1,j)
              ha(r1,j)=r6
            end do

          end if

        end do

      end if

      aflag(5)=aflag(7)/aflag(6)
      iflag(1)=-2
 1110 z=zz
      return
      end

      subroutine y12mcf(n,z,a,snr,nn,rnr,nn1,pivot,b,ha,iha,aflag,
     &  iflag,ifail)

c*********************************************************************72
c
cc Y12MCF finds LU factors of a sparse matrix preprocessed by Y12MBF.
c
      integer iha
      integer n
      integer nn
      integer nn1

      double precision a(nn)
      double precision aflag(8)
      double precision b(n)
      integer c1
      integer c2
      integer cr1
      integer cr2
      integer cr3
      integer cr4
      double precision grmin
      integer ha(iha,11)
      integer i
      integer i1
      integer ifail
      integer iflag(10)
      integer index
      integer j
      integer jj
      integer k
      integer kk
      integer l
      integer l1
      integer l2
      integer l3
      integer l4
      integer l5
      integer l6
      integer l7
      integer lfc
      integer lfr
      integer ll
      integer n7
      integer n8
      integer nr
      double precision pivot(n)
      integer r
      integer r1
      integer r10
      integer r2
      integer r3
      integer r4
      integer r5
      integer r6
      integer r7
      integer r8
      integer r9
      integer rcoll
      integer rpivot
      integer rr
      integer rr1
      integer rr2
      integer rr3
      integer rr4
      integer rrow
      integer rnr(nn1)
      integer slut
      integer snr(nn)
      double precision t
      double precision td
      double precision td1
      double precision tol1
      double precision tol2
      double precision tol3
      double precision u
      double precision v
      integer z
      integer zz
c
c  Information which is necessary to begin the elimination is stored.
c
      ifail=0

      if(iflag(1).ne.-1)ifail=2

      if(aflag(1).lt.1.0) then
        aflag(1)=1.0005
      end if

      if(aflag(3).lt.1.0e+5) then
        aflag(3)=1.0e+5
      end if

      if(aflag(4).lt.0.0) then
        aflag(4)=-aflag(4)
      end if

      if(iflag(2).lt.1)ifail=19

      if(iflag(3).lt.0.or.iflag(3).gt.2)ifail=20

      if(iflag(5).lt.1.or.iflag(5).gt.3) then
        ifail=21
        write(*,*)' '
        write(*,*)'Y12MCF - Fatal error!'
        write(*,*)'  Input value of IFLAG(5)=',iflag(5)
        write(*,*)'  Legal values are 1, 2, or 3.'
        return
      end if

      if(iflag(5).eq.3)ifail=22

      if(ifail.gt.0)return

      snr(z+1)=0
      rnr(z+1)=0
      n8=n+1
      n7=n-1
      u=aflag(1)
      grmin=aflag(4)*aflag(6)
      zz=z
c
c  Use the information about fill-ins if it is possible.
c
      nr=n*n
      if(iflag(4).ne.2)go to 100

      if(iflag(10).le.nn) then

        l1=iflag(10)
        l5=l1+1
        if(l5.le.nn)snr(l5)=0

        do i=1,n

          l=n8-i
          l2=ha(l,3)+1
          l3=l2-ha(l,1)

          do j=1,l3
            snr(l5-j)=snr(l2-j)
            a(l5-j)=a(l2-j)
          end do

          ha(l,3)=l1
          ha(l,1)=l5-l3
          l6=l1-l3
          l5=l5-ha(l,9)

          do j=l5,l6
            snr(j)=0
          end do

          l1=l5-1

        end do

      end if

      if(iflag(9).le.nn1) then

        l2=iflag(9)
        l5=l2+1
        if(l5.le.nn1)rnr(l5)=0

        do i=1,n
          l=n8-i
          l1=ha(l,6)+1

          l4=l1-ha(l,4)
          do j=1,l4
            rnr(l5-j)=rnr(l1-j)
          end do

          ha(l,4)=l5-l4
          ha(l,6)=l2

          l6=l2-l4
          l5=l5-ha(l,10)
          do j=l5,l6
            rnr(j)=0
          end do

          l2=l5-1
        end do

      end if

  100 continue

      r4=ha(n,3)
      r5=ha(n,6)
      aflag(7)=aflag(6)
      aflag(8)=aflag(6)

      do i=1,n
        pivot(i)=0.0
        ha(i,2)=ha(i,1)
        ha(i,5)=ha(i,4)
      end do

      index=ha(n,8)
      slut=ha(index,3)-ha(index,2)+1
c
c  Start of gaussian elimination.
c
      do 950 i=1,n7

      rr3=ha(i,2)
      rr4=ha(i,3)
      c1=ha(i,4)
      cr4=ha(i,6)
      if(iflag(3).eq.0)go to 350

      if(iflag(4).eq.2) then
        rrow=ha(i,7)
        rcoll=ha(i,8)
        go to 220
      end if

      l4=ha(i,8)

      if(iflag(3).ne.1) then
        rrow=l4
        rcoll=rrow
        rpivot=i
        go to 170
      end if

      r=nr
      v=0.0
      index=iflag(2)

      do kk=1,index

        l1=i-1+kk
        if(l1.gt.n)go to 170
        j=ha(l1,8)
        r7=ha(j,2)
        r8=ha(j,3)
        r9=r8-r7

        t=0.0
        do k=r7,r8
          td=abs(a(k))
          if(t.lt.td)t=td
        end do

        t=t/u

        do k=r7,r8

          td=abs(a(k))
          if(td.lt.t)go to 150
          r6=snr(k)
          r3=r9*(ha(r6,6)-ha(r6,5))
          if(r3.gt.r)go to 150
          if(r3.lt.r)go to 151
          if(v.ge.td)go to 150
  151     v=td
          rrow=j
          rcoll=r6
          r=r3
          rpivot=l1

  150     continue

        end do

      end do

  170 r3=ha(rcoll,10)
      ha(rcoll,10)=ha(i,10)
      ha(i,10)=r3
      r3=ha(rrow,9)
      ha(rrow,9)=ha(i,9)
      ha(i,9)=r3
c
c  remove the pivot row of the list where the rows are ordered by
c  increasing numbers of non-zero elements.
c
      l1=0
      l=i
      l2=ha(l4,3)-ha(l4,2)+1
  180 l=l+1
      if(l2.gt.l1)ha(l2,11)=l
      if(l.gt.n)go to 190
      l5=ha(l,8)
      l3=ha(l5,3)-ha(l5,2)+1
      if(rpivot.lt.l)go to 190
      ha(l4,7)=l
      ha(l,8)=l4
      l4=l5
      l1=l2
      l2=l3
      l3=n8
      go to 180
  190 if(l2.eq.l1)go to 200
      if(l3.eq.l2)go to 200
      ha(l2,11)=0
  200 l5=ha(i,7)
      if(rrow.eq.i)go to 210
      ha(l5,8)=rrow
      ha(rrow,7)=l5
  210 ha(i,7)=rrow
      ha(i,8)=rcoll
c
c  Row interchanges.
c
  220 continue

      if(rrow.eq.i)go to 290
      t=b(rrow)
      b(rrow)=b(i)
      b(i)=t

      do j=rr3,rr4
        l1=snr(j)
        r=ha(l1,5)-1
        r10=ha(l1,6)
  240   r=r+1
        if(rnr(r).ne.i)go to 240
        rnr(r)=rnr(r10)
        rnr(r10)=rrow
      end do

      rr3=ha(rrow,2)
      rr4=ha(rrow,3)

      do j=rr3,rr4
        l1=snr(j)
        r=ha(l1,5)-1
  260   r=r+1
        if(rnr(r).ne.rrow)go to 260
        rnr(r)=i
      end do

      do j=1,3
        r3=ha(rrow,j)
        ha(rrow,j)=ha(i,j)
        ha(i,j)=r3
      end do
c
c  column interchanges.
c
  290 if(rcoll.eq.i)go to 350

      do j=c1,cr4
        l1=rnr(j)
        r=ha(l1,2)-1
        r10=ha(l1,3)
  300   r=r+1
        if(snr(r).ne.i)go to 300
        t=a(r10)
        a(r10)=a(r)
        a(r)=t
        snr(r)=snr(r10)
        snr(r10)=rcoll
      end do

      c1=ha(rcoll,4)
      cr4=ha(rcoll,6)

      do j=c1,cr4
        l1=rnr(j)
        r=ha(l1,2)-1
  320   r=r+1
        if(snr(r).ne.rcoll)go to 320
        snr(r)=i
      end do

      do j=4,6
        r3=ha(rcoll,j)
        ha(rcoll,j)=ha(i,j)
        ha(i,j)=r3
      end do
c
c end of the interchanges.
c the row ordered list and the column ordered list are prepared to
c begin step i of the elimination.
c
  350 r9=rr4-rr3

      do rr=rr3,rr4
        if(snr(rr).eq.i)go to 370
      end do

      ifail=9
      go to 1110
  370 v=a(rr)
      pivot(i)=v
      td=abs(v)
      if(td.lt.aflag(8))aflag(8)=td
      if(td.ge.grmin)go to 380
      ifail=3
      go to 1110
  380 r2=ha(i,1)
      a(rr)=a(rr3)
      snr(rr)=snr(rr3)
      a(rr3)=a(r2)
      snr(rr3)=snr(r2)
      snr(r2)=0
      z=z-1
      rr3=rr3+1
      ha(i,2)=rr3
      ha(i,1)=r2+1
      cr3=ha(i,5)
      if(r9.le.0)go to 431

      do j=rr3,rr4
        index=snr(j)
        pivot(index)=a(j)
      end do

  431 r7=cr4-cr3+1
      do 880 k=1,r7
      r1=rnr(cr3-1+k)
      if(r1.eq.i)go to 870
      i1=ha(r1,1)
      rr1=ha(r1,2)
      rr2=ha(r1,3)
      l2=rr2-rr1+1
      l=rr1-1
  390 l=l+1
      if(snr(l).ne.i)go to 390
      t=a(l)/v
      if(iflag(5).eq.2)go to 400
      a(l)=a(i1)
      snr(l)=snr(i1)
      snr(i1)=0
      i1=i1+1
      ha(r1,1)=i1
      z=z-1
      go to 410
  400 a(l)=a(rr1)
      a(rr1)=t
      r3=snr(rr1)
      snr(rr1)=snr(l)
      snr(l)=r3
  410 rr1=rr1+1
      ha(r1,2)=rr1
      b(r1)=b(r1)-b(i)*t
      if(r9.le.0)go to 669
      r=rr1

      do l=r,rr2
        l1=snr(l)
        td=pivot(l1)
        if(td.eq.0.0)go to 450
        pivot(l1)=0.0
        td=a(l)-td*t
        a(l)=td
        td1=abs(td)
        if(td1.gt.aflag(7))aflag(7)=td1
        if(td1.gt.aflag(2))go to 450  ! THR threshold
c
c  Too small element is created.  Remove it from the lists.
c
        z=z-1
        a(l)=a(rr1)
        snr(l)=snr(rr1)
        a(rr1)=a(i1)
        snr(rr1)=snr(i1)
        snr(i1)=0
        rr1=rr1+1
        i1=i1+1
        ha(r1,2)=rr1
        ha(r1,1)=i1
        r3=ha(l1,5)
        r2=r3-1
        l4=ha(l1,4)
        l5=rnr(l4)
        l6=rnr(r3)
  440   r2=r2+1
        if(rnr(r2).ne.r1)go to 440
        rnr(r2)=l6
        rnr(r3)=l5
        rnr(l4)=0
        ha(l1,5)=r3+1
        ha(l1,4)=l4+1
  450   continue
      end do

      do 750 j=1,r9
      r=rr3-1+j
      r2=snr(r)
      tol2=pivot(r2)
      pivot(r2)=a(r)
      if(tol2.eq.0.0)go to 740
      tol3=-tol2*t
      tol1=abs(tol3)
      if(tol1.lt.aflag(2))go to 740  ! THR threshold
      c2=ha(r2,4)
      cr2=ha(r2,6)
      cr1=ha(r2,5)
      lfr=rr2-i1+2
      lfc=cr2-c2+2

      if(iflag(4).eq.1) then
        if(lfr.gt.ha(r1,9))ha(r1,9)=lfr
        if(lfc.gt.ha(r2,10))ha(r2,10)=lfc
      end if

      if(i1.eq.1)go to 490
      if(snr(i1-1).eq.0)go to 600
  490 if(rr2.eq.nn)go to 500
      if(snr(rr2+1).eq.0)go to 580
  500 r10=nn-lfr
c
c  collection in row ordered list.
c
      if(r10.ge.r4)go to 560
      iflag(6)=iflag(6)+1

      do jj=1,n
        l1=ha(jj,3)

        if(l1.ge.ha(jj,1)) then
          ha(jj,3)=snr(l1)
          snr(l1)=-jj
        end if

      end do

      l3=0
      l4=1

      do jj=1,r4

        if(snr(jj).ne.0) then

          l3=l3+1

          if(snr(jj).le.0) then
            l5=-snr(jj)
            snr(jj)=ha(l5,3)
            ha(l5,3)=l3
            l6=l4+ha(l5,2)-ha(l5,1)
            ha(l5,2)=l6
            ha(l5,1)=l4
            l4=l3+1
          end if

          a(l3)=a(jj)
          snr(l3)=snr(jj)

        end if

      end do

      r4=l3
      snr(l3+1)=0
      rr3=ha(i,2)
      rr4=ha(i,3)
      i1=ha(r1,1)
      rr1=ha(r1,2)
      r=rr3-1+j

      if(r10.lt.r4) then
        ifail=5
        write(*,*)' '
        write(*,*)'Y12MCF - Fatal error!'
        write(*,*)'  NN, the declared dimension of A and SNR,'
        write(*,*)'  is too small.'
        write(*,*)'  Current value is NN = ',nn
        return
      end if
c
c  Fill-in takes place in the row ordered list.
c
560   continue

      r8=lfr-1
      rr2=r4+lfr

      l3=i1-1

      do ll=1,r8
        l4=r4+ll
        l5=l3+ll
        a(l4)=a(l5)
        snr(l4)=snr(l5)
        snr(l5)=0
      end do

      rr1=r4+rr1-i1+1
      ha(r1,3)=rr2
      ha(r1,2)=rr1
      i1=r4+1
      ha(r1,1)=i1
      l1=rr2
      go to 590

  580 rr2=rr2+1
      ha(r1,3)=rr2
      l1=rr2
      if(rr2.le.r4)go to 610
  590 r4=rr2
      if(r4.lt.nn)snr(r4+1)=0
      go to 610
  600 rr1=rr1-1
      i1=i1-1
      ha(r1,1)=i1
      ha(r1,2)=rr1
      l1=rr1
      snr(i1)=snr(l1)
      a(i1)=a(l1)
  610 a(l1)=tol3
      snr(l1)=snr(r)
      td=abs(a(l1))
      if(td.gt.aflag(7))aflag(7)=td
      z=z+1
      if(iflag(8).lt.z)      iflag(8)=z
      if(c2.eq.1)go to 620
      if(rnr(c2-1).eq.0)go to 720
  620 if(cr2.eq.nn1)go to 630
      if(rnr(cr2+1).eq.0)go to 700
  630 r10=nn1-lfc
c
c  collection in column ordered list.
c
      if(r10.ge.r5)go to 680
      iflag(7)=iflag(7)+1

      do jj=i,n
        l1=ha(jj,6)
        ha(jj,6)=rnr(l1)
        rnr(l1)=-jj
      end do

      l3=0
      l4=1

      do jj=1,r5

        if(rnr(jj).ne.0) then

          l3=l3+1

          if(rnr(jj).le.0) then
            l5=-rnr(jj)
            rnr(jj)=ha(l5,6)
            ha(l5,6)=l3
            l6=l4+ha(l5,5)-ha(l5,4)
            ha(l5,5)=l6
            ha(l5,4)=l4
            l4=l3+1
          end if

          rnr(l3)=rnr(jj)

        end if

      end do

      r5=l3
      rnr(r5+1)=0
      c2=ha(r2,4)
      cr3=ha(i,5)
      cr4=ha(i,6)
      cr1=ha(r2,5)

      if(r10.lt.r5) then
        ifail=6
        z=zz
        return
      end if
c
c  Fill-in takes place in the column ordered list.
c
680   continue

      r8=lfc-1
      cr2=r5+lfc

      l3=c2-1

      do l=1,r8
        l4=r5+l
        l5=l3+l
        rnr(l4)=rnr(l5)
        rnr(l5)=0
      end do

      cr1=r5+cr1-c2+1
      c2=r5+1
      ha(r2,6)=cr2
      ha(r2,4)=c2
      ha(r2,5)=cr1
      r=cr2
      go to 710

  700 cr2=cr2+1
      ha(r2,6)=cr2
      r=cr2
      if(cr2.le.r5)go to 730
  710 r5=cr2
      if(r5.lt.nn1)rnr(r5+1)=0
      go to 730
  720 cr1=cr1-1
      c2=c2-1
      ha(r2,4)=c2
      ha(r2,5)=cr1
      r=cr1
      rnr(c2)=rnr(r)
  730 rnr(r)=r1
  740 continue
  750 continue
  669 if(rr1.le.rr2)go to 760
      ifail=7
      go to 1110
c
c  update the information in the list where the rows are ordered by
c  increasing numbers of the non-zero elements.
c
  760 if(iflag(4).eq.2)go to 870
      if(iflag(3).eq.0)go to 870
      l1=rr2-rr1+1
      if(l1.eq.l2)go to 870
      l6=ha(r1,7)
      l4=ha(l2,11)
      if(l1.gt.l2)go to 820
      if(l6.gt.l4)go to 780
      if(l4.eq.n)go to 770
      l=ha(l4+1,8)
      l5=ha(l,3)-ha(l,2)+1
      if(l5.eq.l2)go to 790
  770 ha(l2,11)=0
      go to 800
  780 l5=ha(l4,8)
      l3=ha(l6,8)
      ha(l4,8)=l3
      ha(l6,8)=l5
      ha(l5,7)=l6
      ha(l3,7)=l4
      l6=l4
  790 ha(l2,11)=l4+1
  800 if(l4.eq.i+1)go to 810
      l=ha(l6-1,8)
      l2=ha(l,3)-ha(l,2)+1
      l4=ha(l2,11)
      if(l1.lt.l2)go to 780
  810 if(l1.ne.l2)ha(l1,11)=l6
      go to 870
  820 if(l6.gt.l4)go to 840
      if(l4.eq.n)go to 830
      l=ha(l4+1,8)
      l5=ha(l,3)-ha(l,2)+1
      if(l5.eq.l2)go to 840
  830 ha(l2,11)=0
  840 l2=l2+1
      if(l2.le.slut)go to 850
      l3=n
      slut=l1
      l2=l1
      go to 860
  850 l3=ha(l2,11)-1
      if(l3.eq.-1)go to 840
      if(l2.gt.l1)l2=l1
  860 ha(l2,11)=l3
      l4=ha(l3,8)
      l7=ha(l6,8)
      ha(l3,8)=l7
      ha(l6,8)=l4
      ha(l7,7)=l3
      ha(l4,7)=l6
      l6=l3
      if(l2.lt.l1)go to 840
  870 continue
  880 continue

      if(r9.gt.0) then
        do j=rr3,rr4
          index=snr(j)
          pivot(index)=0.0
        end do
      end if

      cr3=ha(i,4)

      do j=cr3,cr4
        rnr(j)=0
      end do

      l2=ha(i,2)-1

      do ll=1,r9

        r=snr(l2+ll)
        r1=ha(r,5)
        r2=ha(r,6)

        if(r2.le.r1) then
          ifail=8
          z=zz
          return
        end if

        ha(r,5)=r1+1
        r3=r1-1

  910   r3=r3+1
        if(rnr(r3).ne.i)go to 910

        rnr(r3)=rnr(r1)
        rnr(r1)=i

      end do

      aflag(5)=aflag(7)/aflag(6)

      if(aflag(5).ge.aflag(3)) then
        ifail=4
        z=zz
        return
      end if

  950 continue
c
c  preparation to begin the back substitution.
c
      index=ha(n,2)
      pivot(n)=a(index)
      a(index)=0.0
      td=abs(pivot(n))
      if(td.gt.aflag(7))aflag(7)=td
      if(td.lt.aflag(8))aflag(8)=td

      if(td.le.grmin) then
        ifail=3
        z=zz
        return
      end if

      if(iflag(4).eq.1) then

        iflag(10)=ha(n,9)
        iflag(9)=ha(n,10)

        do i=1,n7

          r1=n-i
          iflag(10)=iflag(10)+ha(r1,9)
          iflag(9)=iflag(9)+ha(r1,10)

          if(iflag(3).ne.0) then

            do j=9,10
              r2=ha(r1,j-2)
              r6=ha(r2,j)
              ha(r2,j)=ha(r1,j)
              ha(r1,j)=r6
            end do

          end if

        end do

      end if

      aflag(5)=aflag(7)/aflag(6)
      iflag(1)=-2
 1110 z=zz
      return
      end
      subroutine y12mde(n,a,nn,b,pivot,snr,ha,iha,iflag,ifail)

c*********************************************************************72
c
cc Y12MDE solves a linear system which has been factored by Y12MCE.
c
      integer iha
      integer n
      integer nn

      real a(nn)
      real b(n)
      integer ha(iha,11)
      integer i
      integer ifail
      integer iflag(10)
      integer ipiv
      integer j
      integer l1
      integer n7
      integer n8
      real pivot(n)
      integer r1
      integer r2
      integer rr1
      integer rr2
      integer snr(nn)
      integer state
      real t
c
      if(iflag(1).ne.-2) then
        ifail=1
        return
      end if

      ipiv=iflag(3)
      n7=n-1
      state=iflag(5)
c
c  solve the system with lower triangular matrix  l  (if the
c  lu-factorization is available).
c
      if(state.eq.3) then

        if(ipiv.ne.0) then

          do i=1,n7
            l1=ha(i,7)
            t=b(l1)
            b(l1)=b(i)
            b(i)=t
         end do

        end if

        do i=1,n

          rr1=ha(i,1)
          rr2=ha(i,2)-1

          do j=rr1,rr2
            l1=snr(j)
            b(i)=b(i)-a(j)*b(l1)
          end do

        end do

      end if
c
c  Solve the system with upper triagular matrix.
c
      n8=n+1

      do i=1,n

        r1=n8-i
        rr1=ha(r1,2)
        rr2=ha(r1,3)

        do j=rr1,rr2
          r2=snr(j)
          b(r1)=b(r1)-a(j)*b(r2)
        end do

        b(r1)=b(r1)/pivot(r1)

      end do
c
c  If interchanges were used during the elimination then a reordering in
c  the solution vector is made.
c
      if(ipiv.ne.0) then

        do i=1,n7
          r1=n-i
          r2=ha(r1,8)
          t=b(r2)
          b(r2)=b(r1)
          b(r1)=t
        end do

      end if

      return
      end

      subroutine y12mdf(n,a,nn,b,pivot,snr,ha,iha,iflag,ifail)

c*********************************************************************72
c
cc Y12MDF solves a linear system which has been factored by Y12MCF.
c
      integer iha
      integer n
      integer nn

      double precision a(nn)
      double precision b(n)
      integer ha(iha,11)
      integer i
      integer ifail
      integer iflag(10)
      integer ipiv
      integer j
      integer l1
      integer n7
      integer n8
      double precision pivot(n)
      integer r1
      integer r2
      integer rr1
      integer rr2
      integer snr(nn)
      integer state
      double precision t

      if(iflag(1).ne.-2) then
        ifail=1
        return
      end if

      ipiv=iflag(3)
      n7=n-1
      state=iflag(5)
c
c  solve the system with lower triangular matrix  l  (if the
c  lu-factorization is available).
c
      if(state.eq.3) then

        if(ipiv.ne.0) then

          do i=1,n7
            l1=ha(i,7)
            t=b(l1)
            b(l1)=b(i)
            b(i)=t
         end do

        end if

        do i=1,n

          rr1=ha(i,1)
          rr2=ha(i,2)-1

          do j=rr1,rr2
            l1=snr(j)
            b(i)=b(i)-a(j)*b(l1)
          end do

        end do

      end if
c
c  Solve the system with upper triagular matrix.
c
      n8=n+1

      do i=1,n

        r1=n8-i
        rr1=ha(r1,2)
        rr2=ha(r1,3)

        do j=rr1,rr2
          r2=snr(j)
          b(r1)=b(r1)-a(j)*b(r2)
        end do

        b(r1)=b(r1)/pivot(r1)

      end do
c
c  If interchanges were used during the elimination then a reordering in
c  the solution vector is made.
c
      if(ipiv.ne.0) then

        do i=1,n7
          r1=n-i
          r2=ha(r1,8)
          t=b(r2)
          b(r2)=b(r1)
          b(r1)=t
        end do

      end if

      return
      end

      subroutine y12mfe(n,a,snr,nn,rnr,nn1,a1,sn,nz,ha,iha,b,b1,x,y,
     &  aflag,iflag,ifail)

c*********************************************************************72
c
cc Y12MFE solves a single linear system; includes iterative refinement.
c
      integer iha
      integer n
      integer nn
      integer nn1
      integer nz

      real a(nn)
      real a1(nz)
      real aflag(11)
      real b(n)
      real b1(n)
      real d
      real dd
      real dres
      double precision er
      double precision er1
      double precision er2
      real gt1
      real gt2
      integer ha(iha,13)
      integer i
      integer ifail
      integer iflag(12)
      integer it
      integer j
      integer kit
      integer l1
      integer l2
      integer l3
      integer nres
      integer rnr(nn1)
      integer sn(nz)
      integer snr(nn)
      integer state
      real x(n)
      real xm
      real y(n)
c
c  Store the non-zero elements, their column numbers, information about
c  row starts, information about row ends and the right-hand side.
c
      ifail=0
      nres=0
      dres=0.0
      state=iflag(5)
      kit=1
      it=iflag(11)
      if(state.eq.1)ifail=10
      if(it.lt.2)ifail=23
      if(ifail.ne.0)go to 160

      do i=1,n
        b1(i)=b(i)
      end do

      if(state.eq.3)go to 70
      call y12mbe(n,nz,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)
      if(ifail.ne.0)go to 160

      do i=1,nz
        sn(i)=snr(i)
        a1(i)=a(i)
      end do

      do i=1,n
        ha(i,12)=ha(i,1)
        ha(i,13)=ha(i,3)
      end do

      if(aflag(2).ge.0.0)go to 60

      gt1=aflag(6)

      do i=1,n

        l1=ha(i,1)
        l2=ha(i,3)
        gt2=0.0
        do j=l1,l2
          gt2=max(gt2,abs(a(j)))
        end do

        gt1=min(gt1,gt2)

      end do

      aflag(2)=-gt1*aflag(2)
c
c  Find the first solution.
c
   60 call y12mce(n,nz,a,snr,nn,rnr,nn1,y,b,ha,iha,aflag,iflag,ifail)
      if(ifail.ne.0)go to 160

   70 call y12mde(n,a,nn,b,y,snr,ha,iha,iflag,ifail)
      if(ifail.ne.0)go to 160
c
c  prepare the data in order to begin the iterations.
c
      dd=0.0
      do i=1,n
        x(i)=b(i)
        dd=max(dd,abs(b(i)))
      end do

      xm=dd
      if(dd.eq.0.)go to 160
c
c  begin to iterate.
c
   90 d=dd
      dres=0.0

      do i=1,n

        er=b1(i)
        l1=ha(i,12)
        l2=ha(i,13)

        do j=l1,l2
          er1=a1(j)
          l3=sn(j)
          er2=x(l3)
          er=er-er1*er2
        end do
c
c  Store residuals rounded to single precision.
c
        b(i)=er
        dres=max(dres,abs(b(i)))

      end do

      if(dres.eq.0.)go to 160
      if(nres.eq.1) go to 150
      if(dres.gt.1.0e+4*xm)go to 150
      kit=kit+1
      iflag(5)=3

      call y12mde(n,a,nn,b,y,snr,ha,iha,iflag,ifail)

      if(ifail.ne.0)go to 160
c
c  Compute the uniform norm of the current solution vector.
c
      dd=0.0
      do i=1,n
        dd=max(dd,abs(b(i)))
      end do
      if(dd.eq.0.0)go to 160
c
c  Check the convergence criterion.
c
      if(dd.gt.d.and.kit.gt.2)go to 160
c
c  Calculate an improved solution.
c
      xm=0.0
      do i=1,n
        x(i)=x(i)+b(i)
        xm=max(xm,abs(x(i)))
      end do
c
c  check the stopping criteria.
c
      if(10.0+dd/xm.eq.10.0) go to 140
      if(kit.lt.it) go to 90
c
c  end of the iterations.
c
  140 nres=1
      go to 90

  150 dd=abs(dd)

  160 iflag(5)=state
      iflag(12)=kit
      aflag(9)=dd
      aflag(10)=dres
      aflag(11)=xm
      return
      end

      subroutine y12mff(n,a,snr,nn,rnr,nn1,a1,sn,nz,ha,iha,b,b1,x,y,
     &  aflag,iflag,ifail)

c*********************************************************************72
c
cc Y12MFF solves a single linear system; includes iterative refinement.
c
      implicit none
      integer iha
      integer n
      integer nn
      integer nn1
      integer nz

      double precision a(nn)
      double precision a1(nz)
      double precision aflag(11)
      double precision b(n)
      double precision b1(n)
      double precision d
      double precision dd
      double precision dres
      real(16) :: er ! quadruple precision
      real(16) :: er1
      real(16) :: er2
      double precision gt1
      double precision gt2
      integer ha(iha,13)
      integer i
      integer ifail
      integer iflag(12)
      integer it
      integer j
      integer kit
      integer l1
      integer l2
      integer l3
      integer nres
      integer rnr(nn1)
      integer sn(nz)
      integer snr(nn)
      integer state
      double precision x(n)
      double precision xm
      double precision y(n)

      aflag(1)=16.0d0
!       aflag(2)=1.0d-12  ! THR threshold
      aflag(2)=1.0d-16  ! THR threshold
      aflag(3)=1.0d+16
!       aflag(4)=1.0d-12
      aflag(4)=1.0d-16
      ifail=0
!       iflag(2)=2
      iflag(3)=1
      iflag(4)=0
!       iflag(5)=2

C             The value of IFLAG(2)  should  be  a  positive
C                  integer (IFLAG(2) = 3 is recommended).
      iflag(2)=3

C             IFLAG(5) - If   the  LU  factorization  of  the
C                        coefficient matrix is not available,
C                        then  IFLAG(5)  must  be set to 2 on
C                        entry. If the  LU  factorization  of
C                        the coefficient matrix is available,
C                        then IFLAG(5) must be set  to  3  on
C                        entry.  Unchanged on exit.
C
      iflag(5)=2

c
c  Store the non-zero elements, their column numbers, information about
c  row starts, information about row ends and the right-hand side.
c
      nres=0
      dres=0.0
      state=iflag(5)
      kit=1
      it=iflag(11)
      if(state.eq.1)ifail=10
      if(it.lt.2)ifail=23
      if(ifail.ne.0)go to 160

      do i=1,n
        b1(i)=b(i)
      end do

      if(state.eq.3)go to 70
      call y12mbf(n,nz,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)
!       write(*,*)'1 iflag(5)=',iflag(5)
      if(ifail.ne.0)go to 160

!       write(*,*)' nz=',nz,' nn=',nn
      do i=1,nz
        sn(i)=snr(i) ! dim sn = nz, dim snr = nn
        a1(i)=a(i)   ! dim a1 = nz, dim a = nn
!         write(*,*)'i, iflag(5)=',i,iflag(5)
      end do

!       write(*,*)'2 iflag(5)=',iflag(5)


      do i=1,n
        ha(i,12)=ha(i,1)
        ha(i,13)=ha(i,3)
      end do
!       write(*,*)'3 iflag(5)=',iflag(5)

      if(aflag(2).ge.0.0)go to 60

      gt1=aflag(6)

      do i=1,n

        l1=ha(i,1)
        l2=ha(i,3)
        gt2=0.0
        do j=l1,l2
          gt2=max(gt2,abs(a(j)))
        end do

        gt1=min(gt1,gt2)

      end do

      aflag(2)=-gt1*aflag(2)
c
c  Find the first solution.
c
   60 continue
!     write(*,*)'4 iflag(5)=',iflag(5)
!       pause
      call y12mcf(n,nz,a,snr,nn,rnr,nn1,y,b,ha,iha,aflag,iflag,ifail)
      if(ifail.ne.0)go to 160

   70 call y12mdf(n,a,nn,b,y,snr,ha,iha,iflag,ifail)
      if(ifail.ne.0)go to 160
c
c  prepare the data in order to begin the iterations.
c
      dd=0.0
      do i=1,n
        x(i)=b(i)
        dd=max(dd,abs(b(i)))
      end do

      xm=dd
      if(dd.eq.0.)go to 160
c
c  begin to iterate.
c
   90 d=dd
      dres=0.0

      do i=1,n

        er=b1(i)
        l1=ha(i,12)
        l2=ha(i,13)

        do j=l1,l2
          er1=a1(j)
          l3=sn(j)
          er2=x(l3)
          er=er-er1*er2
        end do
c
c  Store residuals rounded to double precision.
c
        b(i)=er
        dres=max(dres,abs(b(i)))

      end do

      if(dres.eq.0.)go to 160
      if(nres.eq.1) go to 150
      if(dres.gt.1.0e+4*xm)go to 150
      kit=kit+1
      iflag(5)=3

      call y12mdf(n,a,nn,b,y,snr,ha,iha,iflag,ifail)

      if(ifail.ne.0)go to 160
c
c  Compute the uniform norm of the current solution vector.
c
      dd=0.0
      do i=1,n
        dd=max(dd,abs(b(i)))
      end do
      if(dd.eq.0.0)go to 160
c
c  Check the convergence criterion.
c
      if(dd.gt.d.and.kit.gt.2)go to 160
c
c  Calculate an improved solution.
c
      xm=0.0
      do i=1,n
        x(i)=x(i)+b(i)
        xm=max(xm,abs(x(i)))
      end do
c
c  check the stopping criteria.
c
      if(10.0+dd/xm.eq.10.0) go to 140
      if(kit.lt.it) go to 90
c
c  end of the iterations.
c
  140 nres=1
      go to 90

  150 dd=abs(dd)

  160 iflag(5)=state
      iflag(12)=kit
      aflag(9)=dd
      aflag(10)=dres
      aflag(11)=xm
      return
      end


      subroutine y12mge(n,nn,a,snr,w,pivot,anorm,rcond,iha,ha,
     &  iflag,ifail)

c*********************************************************************72
c
cc Y12MGE estimates the condition number of a sparse matrix.
c
c  Discussion:
c
c  Y12MGE computes a number called RCOND by Dongarra et al (1979).
c  This number is the reciprocal of the estimated condition number of
c  the matrix A.
c
c  If Y12MGE is to be called, then the user should do two things
c  first:
c
c    1) Call Y12MHE to compute the one-norm of the matrix A.
c    This must be done before A is factored.
c
c    2) Call Y12MCE to compute the LU decomposition of A.
c
c  Since other routines may alter or destroy the LU decomposition
c  computed by Y12MCE, the call to Y12MGE should come immediately
c  after the call to Y12MCE.
c
c
c  Reference:
c
c  Dongarra, Bunch, Moler and Stewart,
c  LINPACK User's Guide", SIAM, Philadelphia, 1979.
c
c
c  N      Input, INTEGER N.  N contains the number of equations in the
c         system Ax=b.
c
c  NN     Input, INTEGER NN.  NN must contain the length of arrays A and SNR.
c         Restriction: NN must be at least 2*Z.
c         Recommended value: 2*Z <= NN <= 3*Z.
c
c  A      Input, REAL A(NN)
c
c         A contains the LU decomposition of the original matrix,
c         as computed by Y12MCE.
c
c  SNR    Input, INTEGER SNR(NN).
c
c         SNR contains the column numbers of the non-zero elements
c         of the upper triangular matrix U (without the column numbers
c         of the diagonal elements of matrix U).
c
c  W      Workspace, REAL W(N).
c
c  PIVOT  Input, REAL PIVOT(N), contains the pivotal elements (the diagonal
c         elements of matrix U).
c
c  ANORM  Input, REAL ANORM, the one-norm of the matrix A, as computed
c         by Y12MHE.
c
c  RCOND  Output, REAL RCOND, the estimate of the reciprocal of the
c         condition number of A.
c
c  IHA    Input, INTEGER IHA.  The first dimension of array HA.
c
c  HA     Input, INTEGER HA(IHA,11), contains information computed
c         by Y12MCE.
c
c  IFLAG  Output, INTEGER IFLAG(10).  The components of this array can be
c         described as follows.
c
c  IFAIL  Output, INTEGER IFAIL, Error diagnostic parameter.
c
c         IFAIL = 0 if the subroutine has not detected any error.
c
c         Positive values of IFAIL on exit show that some error has been
c         detected by the subroutine.  Many of the error diagnostics are common
c         for all subroutines in the package.  Therefore the error diagnostics
c         are listed in a separate section, Section 7, of this book.  We advise
c         the user to check the value of this parameter on exit.
c
      integer n
      integer nn

      real a(nn)
      real anorm
      integer ha(iha,3)
      integer i
      integer ifail
      integer iflag(5)
      integer iha
      integer j
      integer l
      integer l1
      integer l2
      integer l3
      integer n7
      integer n8
      real pivot(n)
      real rcond
      integer snr(nn)
      real t
      real w(n)
      real ynorm
      real znorm
c
c  Check whether the entry is correct or not.
c
      if(ifail.ne.0) then
        rcond=-1.0
        return
      end if

      if(iflag(5).eq.1) then
        ifail=26
        rcond=-1.0
        return
      end if

      n8=n+1
      n7=n-1
c
c  Solve a system of the form u1*w=e where u1 is the
c  transpose of matrix u in the lu-factorization of matrix
c  a and e is a vector whose components are equal to +1
c  or -1.
c
      w(1)=1.0/pivot(1)
      do i=2,n
        w(i)=0.0
      end do

      do i=2,n

        l1=ha(i,2)
        l2=ha(i,3)

        t=w(i-1)
        do j=l1,l2
          l=snr(j)
          w(l)=w(l)+t*a(j)
        end do

        if(w(i).gt.0.0) then
          w(i)=w(i)+1.0
        else
          w(i)=w(i)-1.0
        end if

        w(i)=-w(i)/pivot(i)

      end do
c
c  Solve a system of the form l1*y=w where l1 is the
c  transpose of matrix l in the lu-factorization of
c  matrix a.   the components of vector y are stored
c  array w (thus, the contents of array w are overwritten
c  by the components of vector y).
c
      do i=1,n7

        l=n-i
        l1=ha(l,1)
        l2=ha(l,2)-1

        t=w(l+1)
        do j=l1,l2
          l3=snr(j)
          w(l3)=w(l3)-t*a(j)
        end do

      end do
c
c  Calculate the one-norm of Y.
c
      ynorm=0.0
      do i=1,n
        ynorm=ynorm+abs(w(i))
      end do
c
c   Compute the solution of (lu)z=y.  This means that
c   two systems with triangular matrices are solved using the
c   same ideas as above.  The components of the calculated solution
c   are stored in array W.
c
      do i=1,n

        l1=ha(i,1)
        l2=ha(i,2)-1

      do j=l1,l2
        l=snr(j)
        w(i)=w(i)-a(j)*w(l)
      end do

      end do

      do i=1,n

        l3=n8-i
        l1=ha(l3,2)
        l2=ha(l3,3)

        do j=l1,l2
          l=snr(j)
          w(l3)=w(l3)-a(j)*w(l)
        end do

        w(l3)=w(l3)/pivot(l3)

      end do
c
c  Compute the one-norm of Z, the vector calculated above and stored in W.
c
      znorm=0.0
      do i=1,n
        znorm=znorm+abs(w(i))
      end do
c
c  Find the value of the required estimate for the reciprocal
c  of the condition number of matrix A.
c
      rcond=(ynorm/anorm)/znorm

      return
      end
      subroutine y12mgf(n,nn,a,snr,w,pivot,anorm,rcond,iha,ha,
     &  iflag,ifail)

c*********************************************************************72
c
cc Y12MGF estimates the condition number of a sparse matrix.
c
c  Discussion:
c
c  Y12MGF computes a number called RCOND by Dongarra et al (1979).
c  This number is the reciprocal of the estimated condition number of
c  the matrix A.
c
c  If Y12MGF is to be called, then the user should do two things
c  first:
c
c    1) Call Y12MHF to compute the one-norm of the matrix A.
c    This must be done before A is factored.
c
c    2) Call Y12MCF to compute the LU decomposition of A.
c
c  Since other routines may alter or destroy the LU decomposition
c  computed by Y12MCF, the call to Y12MGF should come immediately
c  after the call to Y12MCF.
c
c
c  Reference:
c
c  Dongarra, Bunch, Moler and Stewart,
c  LINPACK User's Guide", SIAM, Philadelphia, 1979.
c
c
c  N      Input, INTEGER N.  N contains the number of equations in the
c         system Ax=b.
c
c  NN     Input, INTEGER NN.  NN must contain the length of arrays A and SNR.
c         Restriction: NN must be at least 2*Z.
c         Recommended value: 2*Z <= NN <= 3*Z.
c
c  A      Input, REAL A(NN)
c
c         A contains the LU decomposition of the original matrix,
c         as computed by Y12MCF.
c
c  SNR    Input, INTEGER SNR(NN).
c
c         SNR contains the column numbers of the non-zero elements
c         of the upper triangular matrix U (without the column numbers
c         of the diagonal elements of matrix U).
c
c  W      Workspace, REAL W(N).
c
c  PIVOT  Input, REAL PIVOT(N), contains the pivotal elements (the diagonal
c         elements of matrix U).
c
c  ANORM  Input, REAL ANORM, the one-norm of the matrix A, as computed
c         by Y12MHF.
c
c  RCOND  Output, REAL RCOND, the estimate of the reciprocal of the
c         condition number of A.
c
c  IHA    Input, INTEGER IHA.  The first dimension of array HA.
c
c  HA     Input, INTEGER HA(IHA,11), contains information computed
c         by Y12MCF.
c
c  IFLAG  Output, INTEGER IFLAG(10).  The components of this array can be
c         described as follows.
c
c  IFAIL  Output, INTEGER IFAIL, Error diagnostic parameter.
c
c         IFAIL = 0 if the subroutine has not detected any error.
c
c         Positive values of IFAIL on exit show that some error has been
c         detected by the subroutine.  Many of the error diagnostics are common
c         for all subroutines in the package.  Therefore the error diagnostics
c         are listed in a separate section, Section 7, of this book.  We advise
c         the user to check the value of this parameter on exit.
c
      integer n
      integer nn

      double precision a(nn)
      double precision anorm
      integer ha(iha,3)
      integer i
      integer ifail
      integer iflag(5)
      integer iha
      integer j
      integer l
      integer l1
      integer l2
      integer l3
      integer n7
      integer n8
      double precision pivot(n)
      double precision rcond
      integer snr(nn)
      double precision t
      double precision w(n)
      double precision ynorm
      double precision znorm
c
c  Check whether the entry is correct or not.
c
      if(ifail.ne.0) then
        rcond=-1.0
        return
      end if

      if(iflag(5).eq.1) then
        ifail=26
        rcond=-1.0
        return
      end if

      n8=n+1
      n7=n-1
c
c  Solve a system of the form u1*w=e where u1 is the
c  transpose of matrix u in the lu-factorization of matrix
c  a and e is a vector whose components are equal to +1
c  or -1.
c
      w(1)=1.0/pivot(1)
      do i=2,n
        w(i)=0.0
      end do

      do i=2,n

        l1=ha(i,2)
        l2=ha(i,3)

        t=w(i-1)
        do j=l1,l2
          l=snr(j)
          w(l)=w(l)+t*a(j)
        end do

        if(w(i).gt.0.0) then
          w(i)=w(i)+1.0
        else
          w(i)=w(i)-1.0
        end if

        w(i)=-w(i)/pivot(i)

      end do
c
c  Solve a system of the form l1*y=w where l1 is the
c  transpose of matrix l in the lu-factorization of
c  matrix a.   the components of vector y are stored
c  array w (thus, the contents of array w are overwritten
c  by the components of vector y).
c
      do i=1,n7

        l=n-i
        l1=ha(l,1)
        l2=ha(l,2)-1

        t=w(l+1)
        do j=l1,l2
          l3=snr(j)
          w(l3)=w(l3)-t*a(j)
        end do

      end do
c
c  Calculate the one-norm of Y.
c
      ynorm=0.0
      do i=1,n
        ynorm=ynorm+abs(w(i))
      end do
c
c   Compute the solution of (lu)z=y.  This means that
c   two systems with triangular matrices are solved using the
c   same ideas as above.  The components of the calculated solution
c   are stored in array W.
c
      do i=1,n

        l1=ha(i,1)
        l2=ha(i,2)-1

      do j=l1,l2
        l=snr(j)
        w(i)=w(i)-a(j)*w(l)
      end do

      end do

      do i=1,n

        l3=n8-i
        l1=ha(l3,2)
        l2=ha(l3,3)

        do j=l1,l2
          l=snr(j)
          w(l3)=w(l3)-a(j)*w(l)
        end do

        w(l3)=w(l3)/pivot(l3)

      end do
c
c  Compute the one-norm of Z, the vector calculated above and stored in W.
c
      znorm=0.0
      do i=1,n
        znorm=znorm+abs(w(i))
      end do
c
c  Find the value of the required estimate for the reciprocal
c  of the condition number of matrix A.
c
      rcond=(ynorm/anorm)/znorm

      return
      end
      subroutine y12mhe(n,nz,a,snr,work,anorm)

c*********************************************************************72
c
cc Y12MHE computes the one-norm of a sparse matrix.
c
c  Discussion:
c
c  Let R(I) be the sum of the absolute values of the elements of row I
c  of the matrix.  Then ANORM, the one-norm of the matrix, is defined
c  to be the maximum of R(I) for I from 1 to N.
c
c
c  N      Input, INTEGER N.  N contains the number of equations in the
c         system Ax=b.
c
c  Z      Input, INTEGER Z.  Z contains the number of non-zero elements in the
c         coefficient matrix A of the system Ax = b.
c
c  A      Input, REAL A(NN)
c
c         On entry, the first Z locations of array A must contain the non-zero
c         elements of the coefficient matrix A of the system Ax = b. The
c         order of the non-zero elements may be completely arbitrary.
c
c  SNR    Input, INTEGER SNR(NN).
c
c         On entry SNR(j), j = 1 to Z, must contain the column number of the
c         non-zero element stored in A(j).
c
c  WORK   Workspace, REAL WORK(N).
c
c  ANORM  Output, REAL ANORM, the one-norm of the matrix A.
c
      integer n
      integer nz

      real a(nz)
      real anorm
      integer i
      integer l
      integer snr(nz)
      real work(n)

      do i=1,n
        work(i)=0.0
      end do
c
c  Calculate the sums of the absolute values of the non-zero
c  elements in each row of the matrix.
c
      do i=1,nz
        l=snr(i)
        work(l)=work(l)+abs(a(i))
      end do
c
c  Calculate the one-norm of the matrix.
c
      anorm=0.0
      do i=1,n
        anorm=max(anorm,work(i))
      end do

      return
      end
      subroutine y12mhf(n,nz,a,snr,work,anorm)

c*********************************************************************72
c
cc Y12MHF computes the one-norm of a sparse matrix.
c
c  Discussion:
c
c  Let R(I) be the sum of the absolute values of the elements of row I
c  of the matrix.  Then ANORM, the one-norm of the matrix, is defined
c  to be the maximum of R(I) for I from 1 to N.
c
c
c  N      Input, INTEGER N.  N contains the number of equations in the
c         system Ax=b.
c
c  Z      Input, INTEGER Z.  Z contains the number of non-zero elements in the
c         coefficient matrix A of the system Ax = b.
c
c  A      Input, REAL A(NN)
c
c         On entry, the first Z locations of array A must contain the non-zero
c         elements of the coefficient matrix A of the system Ax = b. The
c         order of the non-zero elements may be completely arbitrary.
c
c  SNR    Input, INTEGER SNR(NN).
c
c         On entry SNR(j), j = 1 to Z, must contain the column number of the
c         non-zero element stored in A(j).
c
c  WORK   Workspace, REAL WORK(N).
c
c  ANORM  Output, REAL ANORM, the one-norm of the matrix A.
c
      integer n
      integer nz

      double precision a(nz)
      double precision anorm
      integer i
      integer l
      integer snr(nz)
      double precision work(n)
c
      do i=1,n
        work(i)=0.0
      end do
c
c  Calculate the sums of the absolute values of the non-zero
c  elements in each row of the matrix.
c
      do i=1,nz
        l=snr(i)
        work(l)=work(l)+abs(a(i))
      end do
c
c  Calculate the one-norm of the matrix.
c
      anorm=0.0
      do i=1,n
        anorm=max(anorm,work(i))
      end do

      return
      end
      subroutine y12mve(n,z,a,snr,nn,rnr,nn1,pivot,ha,iha,aflag,
     &  iflag,ifail)

c*********************************************************************72
c
cc Y12MVE finds LU factors of a matrix preprocessed by Y12MBE.
c
c  Discussion:
c
c    This routine is the same as Y12MCE, except that there is no vector B.
c
      integer iha
      integer n
      integer nn
      integer nn1

      real a(nn)
      real aflag(8)
      integer c1
      integer c2
      integer cr1
      integer cr2
      integer cr3
      integer cr4
      real grmin
      integer ha(iha,11)
      integer i
      integer i1
      integer ifail
      integer iflag(10)
      integer index
      integer j
      integer jj
      integer k
      integer kk
      integer l
      integer l1
      integer l2
      integer l3
      integer l4
      integer l5
      integer l6
      integer l7
      integer lfc
      integer lfr
      integer ll
      integer n7
      integer n8
      integer nr
      real pivot(n)
      integer r
      integer r1
      integer r10
      integer r2
      integer r3
      integer r4
      integer r5
      integer r6
      integer r7
      integer r8
      integer r9
      integer rcoll
      integer rpivot
      integer rr
      integer rr1
      integer rr2
      integer rr3
      integer rr4
      integer rrow
      integer rnr(nn1)
      integer slut
      integer snr(nn)
      real t
      real td
      real td1
      real tol1
      real tol2
      real tol3
      real u
      real v
      integer z
      integer zz
c
c  Information which is necessary to begin the elimination is stored.
c
      ifail=0

      if(iflag(1).ne.-1)ifail=2

      if(aflag(1).lt.1.0) then
        aflag(1)=1.0005
      end if

      if(aflag(3).lt.1.0e+5) then
        aflag(3)=1.0e+5
      end if

      if(aflag(4).lt.0.0) then
        aflag(4)=-aflag(4)
      end if

      if(iflag(2).lt.1)ifail=19

      if(iflag(3).lt.0.or.iflag(3).gt.2)ifail=20

      if(iflag(5).lt.1.or.iflag(5).gt.3) then
        ifail=21
        write(*,*)' '
        write(*,*)'Y12MCE - Fatal error!'
        write(*,*)'  Input value of IFLAG(5)=',iflag(5)
        write(*,*)'  Legal values are 1, 2, or 3.'
        return
      end if

      if(iflag(5).eq.3)ifail=22
      if(ifail.gt.0)return

      snr(z+1)=0
      rnr(z+1)=0
      n8=n+1
      n7=n-1
      u=aflag(1)
      grmin=aflag(4)*aflag(6)
      zz=z
c
c  Use the information about fill-ins if it is possible.
c
      nr=n*n
      if(iflag(4).ne.2)go to 100

      if(iflag(10).le.nn) then

        l1=iflag(10)
        l5=l1+1
        if(l5.le.nn)snr(l5)=0

        do i=1,n

          l=n8-i
          l2=ha(l,3)+1
          l3=l2-ha(l,1)

          do j=1,l3
            snr(l5-j)=snr(l2-j)
            a(l5-j)=a(l2-j)
          end do

          ha(l,3)=l1
          ha(l,1)=l5-l3
          l6=l1-l3
          l5=l5-ha(l,9)

          do j=l5,l6
            snr(j)=0
          end do

          l1=l5-1

        end do

      end if

      if(iflag(9).le.nn1) then

        l2=iflag(9)
        l5=l2+1
        if(l5.le.nn1)rnr(l5)=0

        do i=1,n
          l=n8-i
          l1=ha(l,6)+1

          l4=l1-ha(l,4)
          do j=1,l4
            rnr(l5-j)=rnr(l1-j)
          end do

          ha(l,4)=l5-l4
          ha(l,6)=l2

          l6=l2-l4
          l5=l5-ha(l,10)
          do j=l5,l6
            rnr(j)=0
          end do

          l2=l5-1
        end do

      end if

  100 continue

      r4=ha(n,3)
      r5=ha(n,6)
      aflag(7)=aflag(6)
      aflag(8)=aflag(6)

      do i=1,n
        pivot(i)=0.0
        ha(i,2)=ha(i,1)
        ha(i,5)=ha(i,4)
      end do

      index=ha(n,8)
      slut=ha(index,3)-ha(index,2)+1
c
c  Start of gaussian elimination.
c
      do 950 i=1,n7

      rr3=ha(i,2)
      rr4=ha(i,3)
      c1=ha(i,4)
      cr4=ha(i,6)
      if(iflag(3).eq.0)go to 350

      if(iflag(4).eq.2) then
        rrow=ha(i,7)
        rcoll=ha(i,8)
        go to 220
      end if

      l4=ha(i,8)

      if(iflag(3).ne.1) then
        rrow=l4
        rcoll=rrow
        rpivot=i
        go to 170
      end if

      r=nr
      v=0.0
      index=iflag(2)

      do kk=1,index

        l1=i-1+kk
        if(l1.gt.n)go to 170
        j=ha(l1,8)
        r7=ha(j,2)
        r8=ha(j,3)
        r9=r8-r7

        t=0.0
        do k=r7,r8
          td=abs(a(k))
          if(t.lt.td)t=td
        end do

        t=t/u

        do k=r7,r8

          td=abs(a(k))
          if(td.lt.t)go to 150
          r6=snr(k)
          r3=r9*(ha(r6,6)-ha(r6,5))
          if(r3.gt.r)go to 150
          if(r3.lt.r)go to 151
          if(v.ge.td)go to 150
  151     v=td
          rrow=j
          rcoll=r6
          r=r3
          rpivot=l1

  150     continue

        end do

      end do

  170 r3=ha(rcoll,10)
      ha(rcoll,10)=ha(i,10)
      ha(i,10)=r3
      r3=ha(rrow,9)
      ha(rrow,9)=ha(i,9)
      ha(i,9)=r3
c
c  remove the pivot row of the list where the rows are ordered by
c  increasing numbers of non-zero elements.
c
      l1=0
      l=i
      l2=ha(l4,3)-ha(l4,2)+1
  180 l=l+1
      if(l2.gt.l1)ha(l2,11)=l
      if(l.gt.n)go to 190
      l5=ha(l,8)
      l3=ha(l5,3)-ha(l5,2)+1
      if(rpivot.lt.l)go to 190
      ha(l4,7)=l
      ha(l,8)=l4
      l4=l5
      l1=l2
      l2=l3
      l3=n8
      go to 180
  190 if(l2.eq.l1)go to 200
      if(l3.eq.l2)go to 200
      ha(l2,11)=0
  200 l5=ha(i,7)
      if(rrow.eq.i)go to 210
      ha(l5,8)=rrow
      ha(rrow,7)=l5
  210 ha(i,7)=rrow
      ha(i,8)=rcoll
c
c  Row interchanges.
c
  220 continue

      if(rrow.eq.i)go to 290

      do j=rr3,rr4
        l1=snr(j)
        r=ha(l1,5)-1
        r10=ha(l1,6)
  240   r=r+1
        if(rnr(r).ne.i)go to 240
        rnr(r)=rnr(r10)
        rnr(r10)=rrow
      end do

      rr3=ha(rrow,2)
      rr4=ha(rrow,3)

      do j=rr3,rr4
        l1=snr(j)
        r=ha(l1,5)-1
  260   r=r+1
        if(rnr(r).ne.rrow)go to 260
        rnr(r)=i
      end do

      do j=1,3
        r3=ha(rrow,j)
        ha(rrow,j)=ha(i,j)
        ha(i,j)=r3
      end do
c
c  column interchanges.
c
  290 if(rcoll.eq.i)go to 350

      do j=c1,cr4
        l1=rnr(j)
        r=ha(l1,2)-1
        r10=ha(l1,3)
  300   r=r+1
        if(snr(r).ne.i)go to 300
        t=a(r10)
        a(r10)=a(r)
        a(r)=t
        snr(r)=snr(r10)
        snr(r10)=rcoll
      end do

      c1=ha(rcoll,4)
      cr4=ha(rcoll,6)

      do j=c1,cr4
        l1=rnr(j)
        r=ha(l1,2)-1
  320   r=r+1
        if(snr(r).ne.rcoll)go to 320
        snr(r)=i
      end do

      do j=4,6
        r3=ha(rcoll,j)
        ha(rcoll,j)=ha(i,j)
        ha(i,j)=r3
      end do
c
c end of the interchanges.
c the row ordered list and the column ordered list are prepared to
c begin step i of the elimination.
c
  350 r9=rr4-rr3

      do rr=rr3,rr4
        if(snr(rr).eq.i)go to 370
      end do

      ifail=9
      go to 1110
  370 v=a(rr)
      pivot(i)=v
      td=abs(v)
      if(td.lt.aflag(8))aflag(8)=td
      if(td.ge.grmin)go to 380
      ifail=3
      go to 1110
  380 r2=ha(i,1)
      a(rr)=a(rr3)
      snr(rr)=snr(rr3)
      a(rr3)=a(r2)
      snr(rr3)=snr(r2)
      snr(r2)=0
      z=z-1
      rr3=rr3+1
      ha(i,2)=rr3
      ha(i,1)=r2+1
      cr3=ha(i,5)
      if(r9.le.0)go to 431

      do j=rr3,rr4
        index=snr(j)
        pivot(index)=a(j)
      end do

  431 r7=cr4-cr3+1
      do 880 k=1,r7
      r1=rnr(cr3-1+k)
      if(r1.eq.i)go to 870
      i1=ha(r1,1)
      rr1=ha(r1,2)
      rr2=ha(r1,3)
      l2=rr2-rr1+1
      l=rr1-1
  390 l=l+1
      if(snr(l).ne.i)go to 390
      t=a(l)/v
      if(iflag(5).eq.2)go to 400
      a(l)=a(i1)
      snr(l)=snr(i1)
      snr(i1)=0
      i1=i1+1
      ha(r1,1)=i1
      z=z-1
      go to 410
  400 a(l)=a(rr1)
      a(rr1)=t
      r3=snr(rr1)
      snr(rr1)=snr(l)
      snr(l)=r3
  410 rr1=rr1+1
      ha(r1,2)=rr1
      if(r9.le.0)go to 669
      r=rr1

      do l=r,rr2
        l1=snr(l)
        td=pivot(l1)
        if(td.eq.0.0)go to 450
        pivot(l1)=0.0
        td=a(l)-td*t
        a(l)=td
        td1=abs(td)
        if(td1.gt.aflag(7))aflag(7)=td1
        if(td1.gt.aflag(2))go to 450  ! THR threshold
c
c  Too small element is created.  Remove it from the lists.
c
        z=z-1
        a(l)=a(rr1)
        snr(l)=snr(rr1)
        a(rr1)=a(i1)
        snr(rr1)=snr(i1)
        snr(i1)=0
        rr1=rr1+1
        i1=i1+1
        ha(r1,2)=rr1
        ha(r1,1)=i1
        r3=ha(l1,5)
        r2=r3-1
        l4=ha(l1,4)
        l5=rnr(l4)
        l6=rnr(r3)
  440   r2=r2+1
        if(rnr(r2).ne.r1)go to 440
        rnr(r2)=l6
        rnr(r3)=l5
        rnr(l4)=0
        ha(l1,5)=r3+1
        ha(l1,4)=l4+1
  450   continue
      end do

      do 750 j=1,r9
      r=rr3-1+j
      r2=snr(r)
      tol2=pivot(r2)
      pivot(r2)=a(r)
      if(tol2.eq.0.0)go to 740
      tol3=-tol2*t
      tol1=abs(tol3)
      if(tol1.lt.aflag(2))go to 740  ! THR threshold
      c2=ha(r2,4)
      cr2=ha(r2,6)
      cr1=ha(r2,5)
      lfr=rr2-i1+2
      lfc=cr2-c2+2

      if(iflag(4).eq.1) then
        if(lfr.gt.ha(r1,9))ha(r1,9)=lfr
        if(lfc.gt.ha(r2,10))ha(r2,10)=lfc
      end if

      if(i1.eq.1)go to 490
      if(snr(i1-1).eq.0)go to 600
  490 if(rr2.eq.nn)go to 500
      if(snr(rr2+1).eq.0)go to 580
  500 r10=nn-lfr
c
c  collection in row ordered list.
c
      if(r10.ge.r4)go to 560
      iflag(6)=iflag(6)+1

      do jj=1,n
        l1=ha(jj,3)

        if(l1.ge.ha(jj,1)) then
          ha(jj,3)=snr(l1)
          snr(l1)=-jj
        end if

      end do

      l3=0
      l4=1

      do jj=1,r4

        if(snr(jj).ne.0) then

          l3=l3+1

          if(snr(jj).le.0) then
            l5=-snr(jj)
            snr(jj)=ha(l5,3)
            ha(l5,3)=l3
            l6=l4+ha(l5,2)-ha(l5,1)
            ha(l5,2)=l6
            ha(l5,1)=l4
            l4=l3+1
          end if

          a(l3)=a(jj)
          snr(l3)=snr(jj)

        end if

      end do

      r4=l3
      snr(l3+1)=0
      rr3=ha(i,2)
      rr4=ha(i,3)
      i1=ha(r1,1)
      rr1=ha(r1,2)
      r=rr3-1+j

      if(r10.lt.r4) then
        ifail=5
        write(*,*)' '
        write(*,*)'Y12MVE - Fatal error!'
        write(*,*)'  NN, the declared dimension of A and SNR,'
        write(*,*)'  is too small.'
        write(*,*)'  Current value is NN = ',nn
        return
      end if
c
c  Fill-in takes place in the row ordered list.
c
560   continue

      r8=lfr-1
      rr2=r4+lfr

      l3=i1-1

      do ll=1,r8
        l4=r4+ll
        l5=l3+ll
        a(l4)=a(l5)
        snr(l4)=snr(l5)
        snr(l5)=0
      end do

      rr1=r4+rr1-i1+1
      ha(r1,3)=rr2
      ha(r1,2)=rr1
      i1=r4+1
      ha(r1,1)=i1
      l1=rr2
      go to 590

  580 rr2=rr2+1
      ha(r1,3)=rr2
      l1=rr2
      if(rr2.le.r4)go to 610
  590 r4=rr2
      if(r4.lt.nn)snr(r4+1)=0
      go to 610
  600 rr1=rr1-1
      i1=i1-1
      ha(r1,1)=i1
      ha(r1,2)=rr1
      l1=rr1
      snr(i1)=snr(l1)
      a(i1)=a(l1)
  610 a(l1)=tol3
      snr(l1)=snr(r)
      td=abs(a(l1))
      if(td.gt.aflag(7))aflag(7)=td
      z=z+1
      if(iflag(8).lt.z)      iflag(8)=z
      if(c2.eq.1)go to 620
      if(rnr(c2-1).eq.0)go to 720
  620 if(cr2.eq.nn1)go to 630
      if(rnr(cr2+1).eq.0)go to 700
  630 r10=nn1-lfc
c
c  collection in column ordered list.
c
      if(r10.ge.r5)go to 680
      iflag(7)=iflag(7)+1

      do jj=i,n
        l1=ha(jj,6)
        ha(jj,6)=rnr(l1)
        rnr(l1)=-jj
      end do

      l3=0
      l4=1

      do jj=1,r5

        if(rnr(jj).ne.0) then

          l3=l3+1

          if(rnr(jj).le.0) then
            l5=-rnr(jj)
            rnr(jj)=ha(l5,6)
            ha(l5,6)=l3
            l6=l4+ha(l5,5)-ha(l5,4)
            ha(l5,5)=l6
            ha(l5,4)=l4
            l4=l3+1
          end if

          rnr(l3)=rnr(jj)

        end if

      end do

      r5=l3
      rnr(r5+1)=0
      c2=ha(r2,4)
      cr3=ha(i,5)
      cr4=ha(i,6)
      cr1=ha(r2,5)

      if(r10.lt.r5) then
        ifail=6
        z=zz
        return
      end if
c
c  Fill-in takes place in the column ordered list.
c
680   continue

      r8=lfc-1
      cr2=r5+lfc

      l3=c2-1

      do l=1,r8
        l4=r5+l
        l5=l3+l
        rnr(l4)=rnr(l5)
        rnr(l5)=0
      end do

      cr1=r5+cr1-c2+1
      c2=r5+1
      ha(r2,6)=cr2
      ha(r2,4)=c2
      ha(r2,5)=cr1
      r=cr2
      go to 710

  700 cr2=cr2+1
      ha(r2,6)=cr2
      r=cr2
      if(cr2.le.r5)go to 730
  710 r5=cr2
      if(r5.lt.nn1)rnr(r5+1)=0
      go to 730
  720 cr1=cr1-1
      c2=c2-1
      ha(r2,4)=c2
      ha(r2,5)=cr1
      r=cr1
      rnr(c2)=rnr(r)
  730 rnr(r)=r1
  740 continue
  750 continue
  669 if(rr1.le.rr2)go to 760
      ifail=7
      go to 1110
c
c  update the information in the list where the rows are ordered by
c  increasing numbers of the non-zero elements.
c
  760 if(iflag(4).eq.2)go to 870
      if(iflag(3).eq.0)go to 870
      l1=rr2-rr1+1
      if(l1.eq.l2)go to 870
      l6=ha(r1,7)
      l4=ha(l2,11)
      if(l1.gt.l2)go to 820
      if(l6.gt.l4)go to 780
      if(l4.eq.n)go to 770
      l=ha(l4+1,8)
      l5=ha(l,3)-ha(l,2)+1
      if(l5.eq.l2)go to 790
  770 ha(l2,11)=0
      go to 800
  780 l5=ha(l4,8)
      l3=ha(l6,8)
      ha(l4,8)=l3
      ha(l6,8)=l5
      ha(l5,7)=l6
      ha(l3,7)=l4
      l6=l4
  790 ha(l2,11)=l4+1
  800 if(l4.eq.i+1)go to 810
      l=ha(l6-1,8)
      l2=ha(l,3)-ha(l,2)+1
      l4=ha(l2,11)
      if(l1.lt.l2)go to 780
  810 if(l1.ne.l2)ha(l1,11)=l6
      go to 870
  820 if(l6.gt.l4)go to 840
      if(l4.eq.n)go to 830
      l=ha(l4+1,8)
      l5=ha(l,3)-ha(l,2)+1
      if(l5.eq.l2)go to 840
  830 ha(l2,11)=0
  840 l2=l2+1
      if(l2.le.slut)go to 850
      l3=n
      slut=l1
      l2=l1
      go to 860
  850 l3=ha(l2,11)-1
      if(l3.eq.-1)go to 840
      if(l2.gt.l1)l2=l1
  860 ha(l2,11)=l3
      l4=ha(l3,8)
      l7=ha(l6,8)
      ha(l3,8)=l7
      ha(l6,8)=l4
      ha(l7,7)=l3
      ha(l4,7)=l6
      l6=l3
      if(l2.lt.l1)go to 840
  870 continue
  880 continue

      if(r9.gt.0) then
        do j=rr3,rr4
          index=snr(j)
          pivot(index)=0.0
        end do
      end if

      cr3=ha(i,4)

      do j=cr3,cr4
        rnr(j)=0
      end do

      l2=ha(i,2)-1

      do ll=1,r9

        r=snr(l2+ll)
        r1=ha(r,5)
        r2=ha(r,6)

        if(r2.le.r1) then
          ifail=8
          z=zz
          return
        end if

        ha(r,5)=r1+1
        r3=r1-1

  910   r3=r3+1
        if(rnr(r3).ne.i)go to 910

        rnr(r3)=rnr(r1)
        rnr(r1)=i

      end do

      aflag(5)=aflag(7)/aflag(6)

      if(aflag(5).ge.aflag(3)) then
        ifail=4
        z=zz
        return
      end if

  950 continue
c
c  preparation to begin the back substitution.
c
      index=ha(n,2)
      pivot(n)=a(index)
      a(index)=0.0
      td=abs(pivot(n))
      if(td.gt.aflag(7))aflag(7)=td
      if(td.lt.aflag(8))aflag(8)=td

      if(td.le.grmin) then
        ifail=3
        z=zz
        return
      end if

      if(iflag(4).eq.1) then

        iflag(10)=ha(n,9)
        iflag(9)=ha(n,10)

        do i=1,n7

          r1=n-i
          iflag(10)=iflag(10)+ha(r1,9)
          iflag(9)=iflag(9)+ha(r1,10)

          if(iflag(3).ne.0) then

            do j=9,10
              r2=ha(r1,j-2)
              r6=ha(r2,j)
              ha(r2,j)=ha(r1,j)
              ha(r1,j)=r6
            end do

          end if

        end do

      end if

      aflag(5)=aflag(7)/aflag(6)
      iflag(1)=-2
 1110 z=zz
      return
      end
      subroutine y12mvf(n,z,a,snr,nn,rnr,nn1,pivot,ha,iha,aflag,
     &  iflag,ifail)

c*********************************************************************72
c
cc Y12MVF finds LU factors of a matrix preprocessed by Y12MBF.
c
c  Discussion:
c
c    This routine is the same as Y12MCF, except that there is no vector B.
c
      integer iha
      integer n
      integer nn
      integer nn1

      double precision a(nn)
      double precision aflag(8)
      integer c1
      integer c2
      integer cr1
      integer cr2
      integer cr3
      integer cr4
      double precision grmin
      integer ha(iha,11)
      integer i
      integer i1
      integer ifail
      integer iflag(10)
      integer index
      integer j
      integer jj
      integer k
      integer kk
      integer l
      integer l1
      integer l2
      integer l3
      integer l4
      integer l5
      integer l6
      integer l7
      integer lfc
      integer lfr
      integer ll
      integer n7
      integer n8
      integer nr
      double precision pivot(n)
      integer r
      integer r1
      integer r10
      integer r2
      integer r3
      integer r4
      integer r5
      integer r6
      integer r7
      integer r8
      integer r9
      integer rcoll
      integer rpivot
      integer rr
      integer rr1
      integer rr2
      integer rr3
      integer rr4
      integer rrow
      integer rnr(nn1)
      integer slut
      integer snr(nn)
      double precision t
      double precision td
      double precision td1
      double precision tol1
      double precision tol2
      double precision tol3
      double precision u
      double precision v
      integer z
      integer zz
c
c  Information which is necessary to begin the elimination is stored.
c
      ifail=0

      if(iflag(1).ne.-1)ifail=2

      if(aflag(1).lt.1.0) then
        aflag(1)=1.0005
      end if

      if(aflag(3).lt.1.0e+5) then
        aflag(3)=1.0e+5
      end if

      if(aflag(4).lt.0.0) then
        aflag(4)=-aflag(4)
      end if

      if(iflag(2).lt.1)ifail=19

      if(iflag(3).lt.0.or.iflag(3).gt.2)ifail=20

      if(iflag(5).lt.1.or.iflag(5).gt.3) then
        ifail=21
        write(*,*)' '
        write(*,*)'Y12MCF - Fatal error!'
        write(*,*)'  Input value of IFLAG(5)=',iflag(5)
        write(*,*)'  Legal values are 1, 2, or 3.'
        return
      end if

      if(iflag(5).eq.3)ifail=22
      if(ifail.gt.0)return

      snr(z+1)=0
      rnr(z+1)=0
      n8=n+1
      n7=n-1
      u=aflag(1)
      grmin=aflag(4)*aflag(6)
      zz=z
c
c  Use the information about fill-ins if it is possible.
c
      nr=n*n
      if(iflag(4).ne.2)go to 100

      if(iflag(10).le.nn) then

        l1=iflag(10)
        l5=l1+1
        if(l5.le.nn)snr(l5)=0

        do i=1,n

          l=n8-i
          l2=ha(l,3)+1
          l3=l2-ha(l,1)

          do j=1,l3
            snr(l5-j)=snr(l2-j)
            a(l5-j)=a(l2-j)
          end do

          ha(l,3)=l1
          ha(l,1)=l5-l3
          l6=l1-l3
          l5=l5-ha(l,9)

          do j=l5,l6
            snr(j)=0
          end do

          l1=l5-1

        end do

      end if

      if(iflag(9).le.nn1) then

        l2=iflag(9)
        l5=l2+1
        if(l5.le.nn1)rnr(l5)=0

        do i=1,n
          l=n8-i
          l1=ha(l,6)+1

          l4=l1-ha(l,4)
          do j=1,l4
            rnr(l5-j)=rnr(l1-j)
          end do

          ha(l,4)=l5-l4
          ha(l,6)=l2

          l6=l2-l4
          l5=l5-ha(l,10)
          do j=l5,l6
            rnr(j)=0
          end do

          l2=l5-1
        end do

      end if

  100 continue

      r4=ha(n,3)
      r5=ha(n,6)
      aflag(7)=aflag(6)
      aflag(8)=aflag(6)

      do i=1,n
        pivot(i)=0.0
        ha(i,2)=ha(i,1)
        ha(i,5)=ha(i,4)
      end do

      index=ha(n,8)
      slut=ha(index,3)-ha(index,2)+1
c
c  Start of gaussian elimination.
c
      do 950 i=1,n7

      rr3=ha(i,2)
      rr4=ha(i,3)
      c1=ha(i,4)
      cr4=ha(i,6)
      if(iflag(3).eq.0)go to 350

      if(iflag(4).eq.2) then
        rrow=ha(i,7)
        rcoll=ha(i,8)
        go to 220
      end if

      l4=ha(i,8)

      if(iflag(3).ne.1) then
        rrow=l4
        rcoll=rrow
        rpivot=i
        go to 170
      end if

      r=nr
      v=0.0
      index=iflag(2)

      do kk=1,index

        l1=i-1+kk
        if(l1.gt.n)go to 170
        j=ha(l1,8)
        r7=ha(j,2)
        r8=ha(j,3)
        r9=r8-r7

        t=0.0
        do k=r7,r8
          td=abs(a(k))
          if(t.lt.td)t=td
        end do

        t=t/u

        do k=r7,r8

          td=abs(a(k))
          if(td.lt.t)go to 150
          r6=snr(k)
          r3=r9*(ha(r6,6)-ha(r6,5))
          if(r3.gt.r)go to 150
          if(r3.lt.r)go to 151
          if(v.ge.td)go to 150
  151     v=td
          rrow=j
          rcoll=r6
          r=r3
          rpivot=l1

  150     continue

        end do

      end do

  170 r3=ha(rcoll,10)
      ha(rcoll,10)=ha(i,10)
      ha(i,10)=r3
      r3=ha(rrow,9)
      ha(rrow,9)=ha(i,9)
      ha(i,9)=r3
c
c  remove the pivot row of the list where the rows are ordered by
c  increasing numbers of non-zero elements.
c
      l1=0
      l=i
      l2=ha(l4,3)-ha(l4,2)+1
  180 l=l+1
      if(l2.gt.l1)ha(l2,11)=l
      if(l.gt.n)go to 190
      l5=ha(l,8)
      l3=ha(l5,3)-ha(l5,2)+1
      if(rpivot.lt.l)go to 190
      ha(l4,7)=l
      ha(l,8)=l4
      l4=l5
      l1=l2
      l2=l3
      l3=n8
      go to 180
  190 if(l2.eq.l1)go to 200
      if(l3.eq.l2)go to 200
      ha(l2,11)=0
  200 l5=ha(i,7)
      if(rrow.eq.i)go to 210
      ha(l5,8)=rrow
      ha(rrow,7)=l5
  210 ha(i,7)=rrow
      ha(i,8)=rcoll
c
c  Row interchanges.
c
  220 continue

      if(rrow.eq.i)go to 290

      do j=rr3,rr4
        l1=snr(j)
        r=ha(l1,5)-1
        r10=ha(l1,6)
  240   r=r+1
        if(rnr(r).ne.i)go to 240
        rnr(r)=rnr(r10)
        rnr(r10)=rrow
      end do

      rr3=ha(rrow,2)
      rr4=ha(rrow,3)

      do j=rr3,rr4
        l1=snr(j)
        r=ha(l1,5)-1
  260   r=r+1
        if(rnr(r).ne.rrow)go to 260
        rnr(r)=i
      end do

      do j=1,3
        r3=ha(rrow,j)
        ha(rrow,j)=ha(i,j)
        ha(i,j)=r3
      end do
c
c  column interchanges.
c
  290 if(rcoll.eq.i)go to 350

      do j=c1,cr4
        l1=rnr(j)
        r=ha(l1,2)-1
        r10=ha(l1,3)
  300   r=r+1
        if(snr(r).ne.i)go to 300
        t=a(r10)
        a(r10)=a(r)
        a(r)=t
        snr(r)=snr(r10)
        snr(r10)=rcoll
      end do

      c1=ha(rcoll,4)
      cr4=ha(rcoll,6)
      do j=c1,cr4
        l1=rnr(j)
        r=ha(l1,2)-1
  320   r=r+1
        if(snr(r).ne.rcoll)go to 320
        snr(r)=i
      end do

      do j=4,6
        r3=ha(rcoll,j)
        ha(rcoll,j)=ha(i,j)
        ha(i,j)=r3
      end do
c
c end of the interchanges.
c the row ordered list and the column ordered list are prepared to
c begin step i of the elimination.
c
  350 r9=rr4-rr3

      do rr=rr3,rr4
        if(snr(rr).eq.i)go to 370
      end do

      ifail=9
      go to 1110
  370 v=a(rr)
      pivot(i)=v
      td=abs(v)
      if(td.lt.aflag(8))aflag(8)=td
      if(td.ge.grmin)go to 380
      ifail=3
      go to 1110
  380 r2=ha(i,1)
      a(rr)=a(rr3)
      snr(rr)=snr(rr3)
      a(rr3)=a(r2)
      snr(rr3)=snr(r2)
      snr(r2)=0
      z=z-1
      rr3=rr3+1
      ha(i,2)=rr3
      ha(i,1)=r2+1
      cr3=ha(i,5)
      if(r9.le.0)go to 431

      do j=rr3,rr4
        index=snr(j)
        pivot(index)=a(j)
      end do

  431 r7=cr4-cr3+1
      do 880 k=1,r7
      r1=rnr(cr3-1+k)
      if(r1.eq.i)go to 870
      i1=ha(r1,1)
      rr1=ha(r1,2)
      rr2=ha(r1,3)
      l2=rr2-rr1+1
      l=rr1-1
  390 l=l+1
      if(snr(l).ne.i)go to 390
      t=a(l)/v
      if(iflag(5).eq.2)go to 400
      a(l)=a(i1)
      snr(l)=snr(i1)
      snr(i1)=0
      i1=i1+1
      ha(r1,1)=i1
      z=z-1
      go to 410
  400 a(l)=a(rr1)
      a(rr1)=t
      r3=snr(rr1)
      snr(rr1)=snr(l)
      snr(l)=r3
  410 rr1=rr1+1
      ha(r1,2)=rr1
      if(r9.le.0)go to 669
      r=rr1

      do l=r,rr2
        l1=snr(l)
        td=pivot(l1)
        if(td.eq.0.0)go to 450
        pivot(l1)=0.0
        td=a(l)-td*t
        a(l)=td
        td1=abs(td)
        if(td1.gt.aflag(7))aflag(7)=td1
        if(td1.gt.aflag(2))go to 450  ! THR threshold
c
c  Too small element is created.  Remove it from the lists.
c
        z=z-1
        a(l)=a(rr1)
        snr(l)=snr(rr1)
        a(rr1)=a(i1)
        snr(rr1)=snr(i1)
        snr(i1)=0
        rr1=rr1+1
        i1=i1+1
        ha(r1,2)=rr1
        ha(r1,1)=i1
        r3=ha(l1,5)
        r2=r3-1
        l4=ha(l1,4)
        l5=rnr(l4)
        l6=rnr(r3)
  440   r2=r2+1
        if(rnr(r2).ne.r1)go to 440
        rnr(r2)=l6
        rnr(r3)=l5
        rnr(l4)=0
        ha(l1,5)=r3+1
        ha(l1,4)=l4+1
  450   continue
      end do

      do 750 j=1,r9
      r=rr3-1+j
      r2=snr(r)
      tol2=pivot(r2)
      pivot(r2)=a(r)
      if(tol2.eq.0.0)go to 740
      tol3=-tol2*t
      tol1=abs(tol3)
      if(tol1.lt.aflag(2))go to 740  ! THR threshold
      c2=ha(r2,4)
      cr2=ha(r2,6)
      cr1=ha(r2,5)
      lfr=rr2-i1+2
      lfc=cr2-c2+2

      if(iflag(4).eq.1) then
        if(lfr.gt.ha(r1,9))ha(r1,9)=lfr
        if(lfc.gt.ha(r2,10))ha(r2,10)=lfc
      end if

      if(i1.eq.1)go to 490
      if(snr(i1-1).eq.0)go to 600
  490 if(rr2.eq.nn)go to 500
      if(snr(rr2+1).eq.0)go to 580
  500 r10=nn-lfr
c
c  collection in row ordered list.
c
      if(r10.ge.r4)go to 560
      iflag(6)=iflag(6)+1

      do jj=1,n
        l1=ha(jj,3)

        if(l1.ge.ha(jj,1)) then
          ha(jj,3)=snr(l1)
          snr(l1)=-jj
        end if

      end do

      l3=0
      l4=1

      do jj=1,r4

        if(snr(jj).ne.0) then

          l3=l3+1

          if(snr(jj).le.0) then
            l5=-snr(jj)
            snr(jj)=ha(l5,3)
            ha(l5,3)=l3
            l6=l4+ha(l5,2)-ha(l5,1)
            ha(l5,2)=l6
            ha(l5,1)=l4
            l4=l3+1
          end if

          a(l3)=a(jj)
          snr(l3)=snr(jj)

        end if

      end do

      r4=l3
      snr(l3+1)=0
      rr3=ha(i,2)
      rr4=ha(i,3)
      i1=ha(r1,1)
      rr1=ha(r1,2)
      r=rr3-1+j

      if(r10.lt.r4) then
        ifail=5
        write(*,*)' '
        write(*,*)'Y12MVF - Fatal error!'
        write(*,*)'  NN, the declared dimension of A and SNR,'
        write(*,*)'  is too small.'
        write(*,*)'  Current value is NN = ',nn
        return
      end if
c
c  Fill-in takes place in the row ordered list.
c
560   continue

      r8=lfr-1
      rr2=r4+lfr

      l3=i1-1

      do ll=1,r8
        l4=r4+ll
        l5=l3+ll
        a(l4)=a(l5)
        snr(l4)=snr(l5)
        snr(l5)=0
      end do

      rr1=r4+rr1-i1+1
      ha(r1,3)=rr2
      ha(r1,2)=rr1
      i1=r4+1
      ha(r1,1)=i1
      l1=rr2
      go to 590

  580 rr2=rr2+1
      ha(r1,3)=rr2
      l1=rr2
      if(rr2.le.r4)go to 610
  590 r4=rr2
      if(r4.lt.nn)snr(r4+1)=0
      go to 610
  600 rr1=rr1-1
      i1=i1-1
      ha(r1,1)=i1
      ha(r1,2)=rr1
      l1=rr1
      snr(i1)=snr(l1)
      a(i1)=a(l1)
  610 a(l1)=tol3
      snr(l1)=snr(r)
      td=abs(a(l1))
      if(td.gt.aflag(7))aflag(7)=td
      z=z+1
      if(iflag(8).lt.z)iflag(8)=z
      if(c2.eq.1)go to 620
      if(rnr(c2-1).eq.0)go to 720
  620 if(cr2.eq.nn1)go to 630
      if(rnr(cr2+1).eq.0)go to 700
  630 r10=nn1-lfc
c
c  collection in column ordered list.
c
      if(r10.ge.r5)go to 680
      iflag(7)=iflag(7)+1

      do jj=i,n
        l1=ha(jj,6)
        ha(jj,6)=rnr(l1)
        rnr(l1)=-jj
      end do

      l3=0
      l4=1

      do jj=1,r5

        if(rnr(jj).ne.0) then

          l3=l3+1

          if(rnr(jj).le.0) then
            l5=-rnr(jj)
            rnr(jj)=ha(l5,6)
            ha(l5,6)=l3
            l6=l4+ha(l5,5)-ha(l5,4)
            ha(l5,5)=l6
            ha(l5,4)=l4
            l4=l3+1
          end if

          rnr(l3)=rnr(jj)

        end if

      end do

      r5=l3
      rnr(r5+1)=0
      c2=ha(r2,4)
      cr3=ha(i,5)
      cr4=ha(i,6)
      cr1=ha(r2,5)

      if(r10.lt.r5) then
        ifail=6
        z=zz
        return
      end if
c
c  Fill-in takes place in the column ordered list.
c
680   continue

      r8=lfc-1
      cr2=r5+lfc

      l3=c2-1

      do l=1,r8
        l4=r5+l
        l5=l3+l
        rnr(l4)=rnr(l5)
        rnr(l5)=0
      end do

      cr1=r5+cr1-c2+1
      c2=r5+1
      ha(r2,6)=cr2
      ha(r2,4)=c2
      ha(r2,5)=cr1
      r=cr2
      go to 710

  700 cr2=cr2+1
      ha(r2,6)=cr2
      r=cr2
      if(cr2.le.r5)go to 730
  710 r5=cr2
      if(r5.lt.nn1)rnr(r5+1)=0
      go to 730
  720 cr1=cr1-1
      c2=c2-1
      ha(r2,4)=c2
      ha(r2,5)=cr1
      r=cr1
      rnr(c2)=rnr(r)
  730 rnr(r)=r1
  740 continue
  750 continue
  669 if(rr1.le.rr2)go to 760
      ifail=7
      go to 1110
c
c  update the information in the list where the rows are ordered by
c  increasing numbers of the non-zero elements.
c
  760 if(iflag(4).eq.2)go to 870
      if(iflag(3).eq.0)go to 870
      l1=rr2-rr1+1
      if(l1.eq.l2)go to 870
      l6=ha(r1,7)
      l4=ha(l2,11)
      if(l1.gt.l2)go to 820
      if(l6.gt.l4)go to 780
      if(l4.eq.n)go to 770
      l=ha(l4+1,8)
      l5=ha(l,3)-ha(l,2)+1
      if(l5.eq.l2)go to 790
  770 ha(l2,11)=0
      go to 800
  780 l5=ha(l4,8)
      l3=ha(l6,8)
      ha(l4,8)=l3
      ha(l6,8)=l5
      ha(l5,7)=l6
      ha(l3,7)=l4
      l6=l4
  790 ha(l2,11)=l4+1
  800 if(l4.eq.i+1)go to 810
      l=ha(l6-1,8)
      l2=ha(l,3)-ha(l,2)+1
      l4=ha(l2,11)
      if(l1.lt.l2)go to 780
  810 if(l1.ne.l2)ha(l1,11)=l6
      go to 870
  820 if(l6.gt.l4)go to 840
      if(l4.eq.n)go to 830
      l=ha(l4+1,8)
      l5=ha(l,3)-ha(l,2)+1
      if(l5.eq.l2)go to 840
  830 ha(l2,11)=0
  840 l2=l2+1
      if(l2.le.slut)go to 850
      l3=n
      slut=l1
      l2=l1
      go to 860
  850 l3=ha(l2,11)-1
      if(l3.eq.-1)go to 840
      if(l2.gt.l1)l2=l1
  860 ha(l2,11)=l3
      l4=ha(l3,8)
      l7=ha(l6,8)
      ha(l3,8)=l7
      ha(l6,8)=l4
      ha(l7,7)=l3
      ha(l4,7)=l6
      l6=l3
      if(l2.lt.l1)go to 840
  870 continue
  880 continue

      if(r9.gt.0) then
        do j=rr3,rr4
          index=snr(j)
          pivot(index)=0.0
        end do
      end if

      cr3=ha(i,4)

      do j=cr3,cr4
        rnr(j)=0
      end do

      l2=ha(i,2)-1

      do ll=1,r9

        r=snr(l2+ll)
        r1=ha(r,5)
        r2=ha(r,6)

        if(r2.le.r1) then
          ifail=8
          z=zz
          return
        end if

        ha(r,5)=r1+1
        r3=r1-1

  910   r3=r3+1
        if(rnr(r3).ne.i)go to 910

        rnr(r3)=rnr(r1)
        rnr(r1)=i

      end do

      aflag(5)=aflag(7)/aflag(6)

      if(aflag(5).ge.aflag(3)) then
        ifail=4
        z=zz
        return
      end if

  950 continue
c
c  preparation to begin the back substitution.
c
      index=ha(n,2)
      pivot(n)=a(index)
      a(index)=0.0
      td=abs(pivot(n))
      if(td.gt.aflag(7))aflag(7)=td
      if(td.lt.aflag(8))aflag(8)=td

      if(td.le.grmin) then
        ifail=3
        z=zz
        return
      end if

      if(iflag(4).eq.1) then

        iflag(10)=ha(n,9)
        iflag(9)=ha(n,10)

        do i=1,n7

          r1=n-i
          iflag(10)=iflag(10)+ha(r1,9)
          iflag(9)=iflag(9)+ha(r1,10)

          if(iflag(3).ne.0) then

            do j=9,10
              r2=ha(r1,j-2)
              r6=ha(r2,j)
              ha(r2,j)=ha(r1,j)
              ha(r1,j)=r6
            end do

          end if

        end do

      end if

      aflag(5)=aflag(7)/aflag(6)
      iflag(1)=-2
 1110 z=zz
      return
      end

      subroutine y12mwe(n,z,a,snr,nn,rnr,nn1,pivot,ha,iha,aflag,iflag,
     &  ifail)

c*********************************************************************72
c
cc Y12MWE reorders and factors a sparse matrix, no right hand side.
c
      integer iha
      integer n
      integer nn
      integer nn1

      real a(nn)
      real aflag(8)
      integer ha(iha,11)
      integer ifail
      integer iflag(10)
      real pivot(n)
      integer rnr(nn1)
      integer snr(nn)
      integer z
c
      aflag(1)=16.0
      aflag(2)=1.0e-12
      aflag(3)=1.0e+16
      aflag(4)=1.0e-12
      ifail=0
      iflag(2)=2
      iflag(3)=1
      iflag(4)=0
      iflag(5)=2

      call y12mbe(n,z,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)

      if(ifail.ne.0) then
        write(*,*)' '
        write(*,*)'Y12MWE - Fatal error!'
        write(*,*)'  Nonzero value of IFAIL returned by Y12MBE.'
        return
      end if

      call y12mve(n,z,a,snr,nn,rnr,nn1,pivot,ha,iha,aflag,iflag,
     &  ifail)

      if(ifail.ne.0) then
        write(*,*)' '
        write(*,*)'Y12MWE - Fatal error!'
        write(*,*)'  Nonzero value of IFAIL returned by Y12MVE.'
        return
      end if

      return
      end

      subroutine y12mwf(n,z,a,snr,nn,rnr,nn1,pivot,ha,iha,aflag,iflag,
     &  ifail)

c*********************************************************************72
c
cc Y12MWF reorders and factors a sparse matrix, no right hand side.
c
      integer iha
      integer n
      integer nn
      integer nn1

      double precision a(nn)
      double precision aflag(8)
      integer ha(iha,11)
      integer ifail
      integer iflag(10)
      double precision pivot(n)
      integer rnr(nn1)
      integer snr(nn)
      integer z

      aflag(1)=16.0d0
      aflag(2)=1.0d-12  ! THR threshold
      aflag(3)=1.0d+16
      aflag(4)=1.0d-12
      ifail=0
      iflag(2)=2
      iflag(3)=1
      iflag(4)=0
      iflag(5)=2

      call y12mbf(n,z,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)

      if(ifail.ne.0) then
        write(*,*)' '
        write(*,*)'Y12MWF - Fatal error!'
        write(*,*)'  Nonzero value of IFAIL returned by Y12MBF.'
        return
      end if

      call y12mvf(n,z,a,snr,nn,rnr,nn1,pivot,ha,iha,aflag,iflag,
     &  ifail)

      if(ifail.ne.0) then
        write(*,*)' '
        write(*,*)'Y12MWF - Fatal error!'
        write(*,*)'  Nonzero value of IFAIL returned by Y12MVF.'
        return
      end if

      return
      end
      subroutine y12mxe(n,z,a,snr,nn,rnr,nn1,pivot,ha,iha,aflag,iflag,
     &  b,ifail)

c*********************************************************************72
c
cc Y12MXE factors a sparse matrix, solves multiple right hand sides.
c
      integer iha
      integer n
      integer nn
      integer nn1

      real a(nn)
      real aflag(8)
      real b(n)
      integer ha(iha,11)
      integer ifail
      integer iflag(10)
      real pivot(n)
      integer rnr(nn1)
      integer snr(nn)
      integer z

      aflag(1)=16.0
      aflag(2)=1.0e-12  ! THR threshold
      aflag(3)=1.0e+16
      aflag(4)=1.0e-12
      ifail=0
      iflag(2)=2
      iflag(3)=1
      iflag(4)=0
      iflag(5)=2

      call y12mbe(n,z,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)

      if(ifail.ne.0) then
        write(*,*)' '
        write(*,*)'Y12MXE - Fatal error!'
        write(*,*)'  Nonzero value of IFAIL returned by Y12MBE.'
        return
      end if

      call y12mce(n,z,a,snr,nn,rnr,nn1,pivot,b,ha,iha,aflag,iflag,
     &  ifail)

      if(ifail.ne.0) then
        write(*,*)' '
        write(*,*)'Y12MXE - Fatal error!'
        write(*,*)'  Nonzero value of IFAIL returned by Y12MCE.'
        return
      end if

      call y12mde(n,a,nn,b,pivot,snr,ha,iha,iflag,ifail)

      if(ifail.ne.0) then
        write(*,*)' '
        write(*,*)'Y12MXE - Fatal error!'
        write(*,*)'  Nonzero value of IFAIL returned by Y12MDE.'
        return
      end if

      return
      end

      subroutine y12mxf(n,z,a,snr,nn,rnr,nn1,pivot,ha,iha,aflag,iflag,
     &  b,ifail)

c*********************************************************************72
c
cc Y12MXF factors a sparse matrix, solves multiple right hand sides.
c
      integer iha
      integer n
      integer nn
      integer nn1

      double precision a(nn)
      double precision aflag(8)
      double precision b(n)
      integer ha(iha,11)
      integer ifail
      integer iflag(10)
      double precision pivot(n)
      integer rnr(nn1)
      integer snr(nn)
      integer z

      aflag(1)=16.0d0
      aflag(2)=1.0d-12  ! THR threshold
      aflag(3)=1.0d+16
      aflag(4)=1.0d-12
      ifail=0
      iflag(2)=2
      iflag(3)=1
      iflag(4)=0
      iflag(5)=2

      call y12mbf(n,z,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)

      if(ifail.ne.0) then
        write(*,*)' '
        write(*,*)'Y12MXF - Fatal error!'
        write(*,*)'  Nonzero value of IFAIL returned by Y12MBF.'
        return
      end if

      call y12mcf(n,z,a,snr,nn,rnr,nn1,pivot,b,ha,iha,aflag,iflag,
     &  ifail)

      if(ifail.ne.0) then
        write(*,*)' '
        write(*,*)'Y12MXF - Fatal error!'
        write(*,*)'  Nonzero value of IFAIL returned by Y12MCF.'
        return
      end if

      call y12mdf(n,a,nn,b,pivot,snr,ha,iha,iflag,ifail)

      if(ifail.ne.0) then
        write(*,*)' '
        write(*,*)'Y12MXF - Fatal error!'
        write(*,*)'  Nonzero value of IFAIL returned by Y12MDF.'
        return
      end if

      return
      end
