pro read_amdl,lun,nmodel,nn,data,aa,file=kfile, number=knum,ivarr=ivarr
;  reads unformatted adiabatic oscillation model from unit lun,
;  if file is set, opens file and reads model
;  if the unit number is set in lun, this is used,
;  otherwise lun is set to 1.
;  If ivarr is set, assume that the model contains ivarr variables
;
;  nmodel returns a (largely irrelevant) model number
;  nn returns number of mesh points.
;  data(0:7) return D_1 - D_8
;  aa(0,.) returns x = r/R
;  aa(1:ivar,.) return A_1 - A_ivar
;  where ivar depends on D_8

common numamdl, nread

if keyword_set(kfile) then begin
	if keyword_set(lun) eq 0 then lun=1
	openr,lun,kfile,/f77_unformatted
	nread=0
endif

if EOF(lun) then begin
	print,' ****** EOF on unit',lun,''
	nn=-1
	return
endif

point_lun,-lun,pos
data=dblarr(8)
nmodel=1l
nn=1l
readu,lun,nmodel,nn,data

;  test for number of variables
if keyword_set(ivarr) then ivar = ivarr else $
if data(7) lt 10 then ivar = 6 else if data(7) lt 100 then ivar = 7 $
else ivar = 9

aa=dblarr(ivar,nn)
point_lun,lun,pos
readu,lun,nmodel,nn,data,aa
nread=nread+1

if keyword_set(knum) then begin
  while nread lt knum do begin
	readu,lun,nmodel,nn,data,aa
	nread=nread+1
  endwhile
endif

print,'nread, nn =',nread, nn
return
end
