pro read_gsm,lun,cs,file=kfile

;  reads binary grand summary from file.
;  if file is set, opens file and reads first grand summary
;  if the unit number is et in lun, this is used,
;  otherwise lun is set to 1.

;  if file is not set, read is from unit lun, which is assumed
;  to be open

; On EOF, cs(0) and cs(1) are returned as -1.

;  data assumed to be in double precision


if keyword_set(kfile) then begin
	if keyword_set(lun) eq 0 then lun=1
	openr,lun,kfile,/f77_unformatted
endif

cs=dblarr(50) 

if EOF(lun) then begin
	print,' ****** EOF on unit',lun
	cs(0:1)=-1
	return
endif

point_lun,-lun,pos

readu,lun,cs

return
end
