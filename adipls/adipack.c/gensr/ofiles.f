      subroutine ofiles
c
c  reads unit numbers and file names from standard input, in
c  format
c
c  <unit number>   <file name>
c
c  if file name is given as 0, /dev/null is used for the file.
c  input ends with EOF or a line containing -1.
c  returns number of files in nfiles, unit numbers in ids(.),
c  and file names in file(.).
c
c  19/8/87: modified for HP9000 by taking out action option
c     in s/r openf.
c
c  21/8/87: modified for HP9000 by replacing multiple occurences
c     of /dev/null by scratch files, since the HP9000, unfortunately,
c     does not allow several unit numbers to be associated with
c     /dev/null
c
c  4/5/95:  Add s/r openfs, which opens file and stores information
c     in common /cofile/
c
c  9/9/96: Add array iopen as flag for files being actually open.
c     iopen(n) =  1 or 2 for normal file (iopen(n) = 2 is used to
c     flag for newly opened file).
c     iopen(n) = -1 for scratch file
c
c  14/10/04: Add option to add trailer to file name (mainly for
c     iteration with evolution code, etc.). The trailer must be
c     set in common/trl_param/ and is added in s/r openf if
c     status is `nt', `ut' or `ot'
c
c  21/10/04: Add s/r newfil to test whether the file has been newly
c     opened, by testing the value of the relevant iopen(k). AFter
c     call of newfil iopen is reset to flag for old file.
c
c  17/3/05: Include istdin for standard input
c  ..............................................
c
      character*280 file, filein, filess
      common/cofile/ nfiles, ids(20), file(20), iopen(20), filess(20)
c
c  nfiles and iopen initialized in block data blopen below
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  the following line required on CR32
c
c..      save nfiles, ids, file
c
      write(istdou,*) 'Input format: <unit number>   <file name>'
      write(istdou,*) 'input ends with EOF or a line containing -1.'
      write(istdou,*)
     *  'if file name is given as 0, /dev/null is used for the file.'
c
   10 read(istdin,*, end=30) idsin, filein
      if(idsin.lt.0) go to 30
c
c  test for /dev/null
c
      if(filein.eq.'0') then
        filein='/dev/null'
      end if
c
      nfiles=nfiles+1
      ids(nfiles)=idsin
      file(nfiles)=filein
      go to 10
c
   30 continue
c
c  output file information
c
      write(istdou,100)
      do 40 n=1,nfiles
   40 write(istdou,110) ids(n),file(n)
      write(istdou,120)
      return
  100 format(/' files set in s/r ofiles:'/)
  110 format(i3,2x,'''',a,'''      @')
  120 format('-1      ''''          @')
      end
      subroutine stfile(idsst, nfst)
c
c  find number of file nfst corresponding to unit number idsst.
c  list of unit numbers and file names in ids and file must have been
c  set up in common/cofile/ by call of ofiles.
c
      character*280 file, filein, filess
      common/cofile/ nfiles, ids(20), file(20), iopen(20), filess(20)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  the following line required on CR32
c
c..      save nfiles, ids, file
c
      do 10 i=1,nfiles
      if(idsst.eq.ids(i)) then
        nfst=i
        go to 20
      end if
   10 continue
c
c  idsst not found, print diagnostics
c
      write(istdou,'(/i4,'' not found''//
     *  '' List of files available:'')') idsst
      do 15 i=1,nfiles
      lfile=length(file(i))
   15 write(istdou,'(i4,3x,a)') ids(i),file(i)(1:lfile)
c
      nfst=-1
c
   20 continue
      return
      end
      subroutine addfil(idsnew, new_file)
c
c  Adds file file_name to list of files set by ofiles.
c  If unit number ids already exists, replaces file name,
c  otherwise adds to list of unit numbers.
c 
c  Original version: 11/2/06
c
      character*(*) new_file
      character*280 file, filein, filess
      common/cofile/ nfiles, ids(20), file(20), iopen(20), filess(20)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      do 10 i=1,nfiles
      if(idsnew.eq.ids(i)) then
        nfst=i
        go to 20
      end if
   10 continue
c
c  ids not found, add number to list
c
      if(istdpr.gt.0) write(istdpr,'(/'' File number '', i3,
     *  '' added to list''/)') idsnew
c
      nfiles=nfiles+1
      nfst=nfiles
      ids(nfst)=idsnew
c
   20 continue
c
      file(nfst)=new_file
      if(istdpr.gt.0) write(istdpr,'(/'' New unit, file:'',
     *  i3,2x,a/)') ids(nfst), file(nfst)
      iopen(nfst)=0
      return
      end
      logical function opnfil(idsst)
c
c  Returns .true. if file has been opened, false otherwise.
c  Diagnostics are returned depending on the hard-coded idiag
c
c  Original version: 3/8/05
c
      character*280 file, filein, filess
      common/cofile/ nfiles, ids(20), file(20), iopen(20), filess(20)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  the following line required on CR32
c
c..      save nfiles, ids, file
c
      idiag=0
c
      nfst=-1
      do 10 i=1,nfiles
      if(idsst.eq.ids(i)) then
        nfst=i
      end if
   10 continue
c
c  idsst not found, print diagnostics
c
      if(idiag.gt.0.and.istdpr.gt.0) then
        write(istdpr,*) idsst,' not found'
        write(istdpr,*) 'List of files available:'
        do 15 i=1,nfiles
   15   write(istdpr,*) ids(i),'  ',file(i)
      end if
c
      opnfil=nfst.ne.-1
c
      return
      end
      logical function nscfil(idsst)
c
c  returns true if file corresponding to unit number idsst
c  is not a scratch file, false otherwise, including if the
c  unit number has not been set or is negative.
c
      character*280 file, filein, filess
      common/cofile/ nfiles, ids(20), file(20), iopen(20), filess(20)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(idsst.le.0) then
	nscfil=.false.
      else
        call stfile(idsst, nfst)
        if(istdpr.gt.0) 
     *    write(istdpr,*) 'In nscfil, idsst, nfst, iopen(nfst) =',
     *    idsst, nfst, iopen(nfst)
        if(nfst.eq.-1) then
          nscfil=.false.
        else
          nscfil = iopen(nfst).gt.0
        end if
      end if
      return
      end
      logical function newfil(idsst)
c
c  returns true if file corresponding to unit number idsst
c  is newly opened and not a scratch file and resets relevant iopen
c  to flag file as old.
c
      character*280 file, filein, filess
      common/cofile/ nfiles, ids(20), file(20), iopen(20), filess(20)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(idsst.le.0) then
	newfil=.false.
      else
        call stfile(idsst, nfst)
        if(istdpr.gt.0) 
     *    write(istdpr,*) 'In newfil, idsst, nfst, iopen(nfst) =',
     *    idsst, nfst, iopen(nfst)
        if(nfst.eq.-1) then
          newfil=.false.
        else 
  	  if(iopen(nfst).lt.0) then
	    newfil=.false.
	  else
            newfil = iopen(nfst).gt.1
	    iopen(nfst) = 1
          end if
        end if
      end if
      return
      end
      subroutine openf(id,status,form)
c
c  open file with unit number id, status as in string status,
c  and format as in string form.
c  status and form may be abbreviated to a single character,
c  as i.e.'n' for 'new', 'u' for 'unformatted'.
c
c  s/r ofiles must have been called previously to set up
c  nfiles, ids and file in common /cofile/.
c
c  original version 30/9/86
c
c            ....................................
c
      character*(*) status, form
      character*280 stat1, form1, file, files, trailer_par, filess,
     *  strcompr
      character*24 fdate
      character ss1*1, ss2*2, stime*10
      common/cofile/ nfiles, ids(20), file(20), iopen(20), filess(20)
      common/trl_param/ trailer_par
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      !external time
      data trailer_par /''/
c
      save 
      data nnul /0/
c
c  find file name
c
      call stfile(id,nfin)
      if(nfin.lt.0) go to 90
c
c  set full status
c
      if(status(1:1).eq.'o') then
        stat1='old'
      else if(status(1:1).eq.'n') then
        stat1='new'
      else if(status(1:1).eq.'s') then
        stat1='scratch'
      else 
        stat1='unknown'
      end if
c
c  for /dev/null, set status to old
c
      if(file(nfin).eq.'/dev/null') then
        nnul = min(9,nnul+1)
	if(istdpr.gt.0) write(istdpr,*) 'Now nnul =',nnul
        if(nnul.eq.1) then
          stat1='old'
        else
          stat1='scratch'
        end if
      end if
c
c  set full format
c
      if(form(1:1).eq.'u') then
        form1='unformatted'
      else
        form1='formatted'
      end if 
c
c  open file, if not already open
c
      if(stat1.eq.'scratch') then
	write(ss1,'(i1)') nnul
	     call system_clock(itime)
        !itime = time()
	write(stime,'(i10)') itime
c..        files='scratch/ttt.'//fdate()//'.'//ss1
c..        files=strcompr(files)
        files='scratch/ttt.'//stime//'.'//ss1
	if(iopen(nfin).eq.-1) then
          write(istdou,'(2a)') 'Scratch file already open: ',
     *      filess(nfin)
        else
          write(istdou,'(2a)') 'Scratch file: ',files
          open(id,status='unknown',file=files, form=form1)
	  iopen(nfin)=-1
	  filess(nfin)=files
        end if
c
c  diagnostic output
c
        if(istdpr.gt.0) write(istdpr,110) id,stat1,form1
c
      else
	if(status(1:2).eq.'ut'.or.status(1:2).eq.'nt'.or.
     *    status(1:2).eq.'ot') then
	  files=file(nfin)
	  if(files.ne.'/dev/null') then
            files=
     *        files(1:length(files))//trailer_par(1:length(trailer_par))
	  end if
c..	  write(6,*) ' Set files to ',files
        else
	  files=file(nfin)
c..	  write(6,*) ' Set files to ',files
	end if
        open(id,file=files,status=stat1,form=form1)
	iopen(nfin)=2
c
	filess(nfin)=files
c
c  diagnostic output
c
        if(istdpr.gt.0) write(istdpr,100) id,files,stat1,form1
c
      end if
c
      return
c
c  error in locating file name. exit
c
   90 stop
  100 format(' open(',i3,',file=',a30,',status=',a10,',form=',a10,')')
  110 format(' open(',i3,',status=',a,',form=',a,')')
      end
      subroutine openfc(id,idp,status,form)
c
c  open file with unit number id, status as in string status,
c  format as in string form. For details, see s/r openf.
c
c  open only takes place if id .ne. idp or trailer_par has been
c  changed. If idp .gt. 0, also closes unit idp. 
c  idp is returned as id.
c
      logical newfile
      character*(*) status, form
      character*280 stat1, form1, file, files, filess, trailer_par
      common/cofile/ nfiles, ids(20), file(20), iopen(20), filess(20)
      common/trl_param/ trailer_par
c
      save
c
      newfile=.false.
c
c  find file name
c
      call stfile(id,nfin)
      if(nfin.lt.0) go to 90
      if(status(1:2).eq.'ut'.or.status(1:2).eq.'nt'.or.
     *  status(1:2).eq.'ot') then
        files=file(nfin)
	if(files.ne.'/dev/null') files=
     *      files(1:length(files))//trailer_par(1:length(trailer_par))
c..	write(6,*) files
c..	write(6,*) filess(nfin)
	newfile=files.ne.filess(nfin)
      end if
c
      newfile=newfile.or.(id.ne.idp)
      if(newfile) then
        if(idp.gt.0) close(idp)
        call openf(id,status,form)
        idp=id
      end if
c
      return
c
   90 stop 'Error in openfc'
      end
      subroutine openfs(id,fn,status,form)
c
c  open file with unit number id, filename fn,
c  status as in string status,
c  format as in string form. For details, see s/r openf.
c
c  Also stores file name in list in common/cofile/, for later
c  access by, say, stfile.
c
      character*(*) fn,status, form
      character*280 file, filess
      common/cofile/ nfiles, ids(20), file(20), iopen(20), filess(20)
c
      nfiles=nfiles+1
      ids(nfiles)=id
      file(nfiles)=fn
c
      call openf(id,status,form)
c
      return
      end
      subroutine dmpfil(idsdmp)
c
c  Output list of files, on unit idsdmp,
c  established by s/r ofiles in common /cofile/ and opened.
c
c  Original version: 9/9/96
c  ..............................................
c
      character*280 file, filout, filess
      common/cofile/ nfiles, ids(20), file(20), iopen(20), filess(20)
c
      write(idsdmp,100)
      do 20 n=1,nfiles
      if(iopen(n).ne.0) then
        ll=length(file(n))
        filout=file(n)
        write(idsdmp,110) ids(n),filout(1:ll)
      end if
   20 continue
      return
  100 format(' Assignment of unit numbers to files:'/)
  110 format(i3,':',2x,a)
      end
      block data blopen
c
c  initialize array iopen
c
c  New version: 7/10/97
c
      character*280 file, filess
      common/cofile/ nfiles, ids(20), file(20), iopen(20), filess(20)
      data nfiles /0/
      data iopen /20*0/
      data filess /20*''/
      end
      subroutine close_file(idsst)
c
c  closes file corresponding to unit ids and resets flag for open file
c
      character*280 file, filein, filess
      common/cofile/ nfiles, ids(20), file(20), iopen(20), filess(20)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      call stfile(idsst, nfs)
      iopen(nfs)=0
      close(idsst)
      return
      end
