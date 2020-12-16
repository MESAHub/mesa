pro read_amdes,lun,x,y,nst,cs,file=kfile,double=kdouble,diag=kdiag, $
  mcase=mcase,nskip=nskip,nrdmax=nrdmax

;  reads several binary eigenfunctions from file.
;  if file is set, opens file and reads first eigenfunction
;  if the unit number is set in lun, this is used,
;  otherwise lun is set to 1.
;  if mcase is set and is equal to 2, read reduced set
;  Otherwise full set is assumed
;  If nskip is set, return kernels at every nskip-th point
;  If keyword nrdmax is set, read at most nrdmax eigenfunctions.

;  if file is not set, read is from unit lun, which is assumed
;  to be open

;  if double is set, the data are in double precision,
;  otherwise in single precision

;  Note: it is assumed that x is the same for all eigenfunctions.
;  x(0:nn-1) returns the radial mesh (r/R)
;  nrd+1 returns the number of modes
;  The eigenfunctions are returned in y(0:ivar-1,0:nn-1,0:nrd)
;  cs(0:49,0:nrd) return the grand summaries for the modes.

;  x
;  
common c_read_amde,init,nn,xx

numamdemax=1000

if keyword_set(mcase) then icase = mcase else icase = 1

if keyword_set(nrdmax) then mrdmax = nrdmax else mrdmax = numamdemax

if keyword_set(kfile) then begin
	if keyword_set(lun) eq 0 then lun=1
	openr,lun,kfile,/f77_unformatted
        init=1
endif else begin
	init=0
endelse

if EOF(lun) then begin
	print,' ****** EOF on unit',lun,'^G'
	n=-1
	return
endif

n=1l
if keyword_set(kdouble) then begin
  csr=dblarr(50) 
  cs =dblarr(50,mrdmax+1)
endif else begin
  csr=fltarr(50)
  cs =fltarr(50,mrdmax+1)
endelse

point_lun,-lun,pos

if icase eq 2 then begin
	if init eq 1 then begin
	  readu,lun,n
	  point_lun,lun,pos
          if keyword_set(kdouble) then x=dblarr(n) else x=fltarr(n)
	  readu,lun,n,x
	  nn=n
	  xx=x
	endif
        if keyword_set(kdouble) then yr=dblarr(2,nn) else yr=fltarr(2,nn)
	iy=2
	readu,lun,csr,yr
	x=xx
endif else begin
	readu,lun,csr,n
	if keyword_set(kdouble) then a=dblarr(7,n) else a=fltarr(7,n)
	iy=6
	point_lun,lun,pos
	readu,lun,csr,n,a
	x=a(0,*)
	yr=a(1:6,*)
endelse

if keyword_set(nskip) then nsk=nskip else nsk = 1

if nsk gt 1 then begin
	nst=1+(n-1)/nsk
	sel=nsk*indgen(nst)
	x=x(sel)
endif else begin
	nst=n
	sel=indgen(n)
endelse

;  set storage y array and set first mode

if keyword_set(double) then y=dblarr(iy,nst,mrdmax+1) else $
			    y=fltarr(iy,nst,mrdmax+1) 

cs(*,0)=csr
y(*,*,0)=yr(*,sel)
if keyword_set(kdiag) then $
  print,' In read_amde l, n, freq  =', csr(17),csr(18),csr(26)

imd = 0

while (NOT EOF(lun) and imd lt mrdmax) do begin
  if icase eq 2 then begin
	readu,lun,csr,yr
  endif else begin
	readu,lun,csr,n,a
	yr=a(1:6,*)
  endelse
  imd=imd+1
  cs(*,imd)=csr
  y(*,*,imd)=yr(*,sel)

  if keyword_set(kdiag) then $
    print,' In read_amde l, n, freq  =', csr(17),csr(18),csr(26)
endwhile

cs=cs(*,0:imd)
y=y(*,*,0:imd)

return
end
