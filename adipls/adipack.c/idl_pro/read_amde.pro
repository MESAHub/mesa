pro read_amde,lun,x,y,n,cs,file=kfile,double=kdouble,diag=kdiag,mcase=mcase, $
   xtest=xtest

;  reads binary eigenfunctions from file.
;  if file is set, opens file and reads first eigenfunction
;  if the unit number is et in lun, this is used,
;  otherwise lun is set to 1.
;  if mcase is set and is equal to 2, read reduced set
;  Otherwise full set is assumed

;  if xtest is set and mcase = 2, test for each new record whether
;  it is a record containing the mesh or a mdoe record.

;  if file is not set, read is from unit lun, which is assumed
;  to be open

;  if double is set, the data are in double precision,
;  otherwise in single precision

common c_read_amde,init,nn,xx

if keyword_set(mcase) then icase = mcase else icase = 1
if keyword_set(xtest) then xt = 1 else xt = 0

if keyword_set(kfile) then begin
	if keyword_set(lun) eq 0 then lun=1
	openr,lun,kfile,/f77_unformatted
        init=1
endif else begin
	init=0
endelse

if EOF(lun) then begin
	print,' ****** EOF on unit',lun,'^G'
	n=-1l
	return
endif

n=1l
if keyword_set(kdouble) then cs=dblarr(50) else cs=fltarr(50)

point_lun,-lun,pos

if icase eq 2 then begin
	if init eq 1 then begin
	  readu,lun,n
	  point_lun,lun,pos
          if keyword_set(kdouble) then x=dblarr(n) else x=fltarr(n)
	  readu,lun,n,x
	  nn=n
	  xx=x
	endif else if xt then begin
	  nt=1l
	  xxt=1.
	  readu,lun,nt
	  point_lun,lun,pos
	  readu,lun,xxt
	  point_lun,lun,pos
	  print,nt,xxt
	  if nt gt 0 and nt lt 1000000 then begin
	    print,' New mesh record'
	    n=nt
            if keyword_set(kdouble) then x=dblarr(n) else x=fltarr(n)
	    readu,lun,n,x
	    nn=n
	    xx=x
	  endif else begin
	    n=nn
          endelse
	endif else begin
	  n=nn
        endelse
        if keyword_set(kdouble) then y=dblarr(2,nn) else y=fltarr(2,nn)
	readu,lun,cs,y
	x=xx
endif else begin
	readu,lun,cs,n
	if keyword_set(kdouble) then a=dblarr(7,n) else a=fltarr(7,n)
	point_lun,lun,pos
	readu,lun,cs,n,a
	x=a(0,*)
	y=a(1:6,*)
endelse

if keyword_set(kdiag) then $
  print,' In read_amde l, n, freq  =', cs(17),cs(18),cs(26)
return
end
