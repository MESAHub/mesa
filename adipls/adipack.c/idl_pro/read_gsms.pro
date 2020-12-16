pro read_gsms,a,file=file,modpar=modpar, gsms=gsms
;  reads all grand summaries on file file
;  returns:
;  a(0,*): radius
;  a(1,*): l
;  a(2,*): order
;  a(3,*): sigma^2
;  a(4,*): nu (microHz)
;  a(5,*): E

;  If keyword modpar is set, in addition returns model parameters in
;  modpar:

;  modpar(0,*): mass
;  modpar(1,*): radius
;  modpar(2,*): central pressure
;  modpar(3,*): central density
;  modpar(4,*): scaled central second derivative of pressure
;  modpar(5,*): scaled central second derivative of density
;  modpar(6,*): polytropic index

;  If keyword gsms is set, in addition stores entire
;  grand summary in gsms

a=dblarr(6,10000)
if keyword_set(modpar) then modpar=dblarr(7,10000)
if keyword_set(gsms) then gsms=dblarr(50,10000)

read_gsm,1,cs,file=file
i=-1

while cs(1) gt 0 do begin
  i=i+1
  a(0,i)=cs(2)
  a(1,i)=cs(17)
  a(2,i)=cs(18)
  a(3,i)=cs(20)
  a(4,i)=1000.*cs(36)
  a(5,i)=cs(23)
  if keyword_set(modpar) then modpar(0:6,i)=cs(1:7)
  if keyword_set(gsms) then gsms(*,i)=cs
  read_gsm,1,cs
endwhile

a=a(*,0:i)
if keyword_set(modpar) then modpar=modpar(*,0:i)
if keyword_set(gsms) then gsms=gsms(*,0:i)
close,1
return
end
