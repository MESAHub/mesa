	function length(t)
c
c  returns length of non-blank portion of string t
c
	character*(*) t
c
	l=len(t)
	do 1 i=1,l
	k = l-i+1
	if (t(k:k).ne.' ') then 
	  go to 2
	endif
    1	continue
	k = 0
    2	length = k
	return
	end
