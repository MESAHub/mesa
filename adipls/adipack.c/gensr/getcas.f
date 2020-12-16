      integer function getcas(string)
c
c  extracts and returns case number, defined as number coming after 
c  last ":" in string. If no :, or no number after :, returns 0
c
c  Original version: 27/3/96
c
      character*(*) string
c
      l=len(string)
c
c  locate :
c
      do 10 i=1,l
      k=l+1-i
      if(string(k:k).eq.':') go to 15
   10 continue
c
c  no colon found. Return 0
c
      getcas=0
      return
c
c Try to read integer from last part of string
c
   15 continue
      read(string(k+1:l),*,err=20) i
c
c  if successful, return integer
c
      getcas=i
      return
c
c  for error, return 0
c
   20 getcas=0
      return
      end
