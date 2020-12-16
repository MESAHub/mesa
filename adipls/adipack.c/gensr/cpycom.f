      subroutine cpycom(ids, idsout)
c
c  copies initial comment lines on unit ids to unit idsout.
c  a comment is flagged bu a "#" in column no 1.
c  idsout is assumed to be open, before call of cpycom.
c
c  on exit file is positioned at first non-comment line
c
      character i*256,ic*1
      data ic /'#'/
c
   10 read(ids,100,end=90,err=90) i
      if(i(1:1).eq.ic) then
        write(idsout,100) i(1:length(i))
        go to 10
      end if
c
      backspace ids
      return
c
c  diagnostics for errors
c
   90 write(6,*) 'error or EOF detected in s/r cpycom on unit',
     *  ids
      return
  100 format(a)
      end
      subroutine cpycmn(ids, idsout, nline)
c
c  copies nline first comment lines on unit ids to unit idsout.
c  and skips the rest.
c  a comment is flagged by a "#" in column no 1.
c  idsout is assumed to be open, before call of cpycom.
c
c  on exit file is positioned at first non-comment line
c
      character i*256,ic*1
      data ic /'#'/
c
      nl=0
c
   10 read(ids,100,end=90,err=90) i
      nl=nl+1
      if(i(1:1).eq.ic) then
        if(nl.le.nline) write(idsout,100) i(1:length(i))
        go to 10
      end if
c
      backspace ids
      return
c
c  diagnostics for errors
c
   90 write(6,*) 'error or EOF detected in s/r cpycom on unit',
     *  ids
      return
  100 format(a)
      end
