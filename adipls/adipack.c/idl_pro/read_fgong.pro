pro read_fgong,head,glob,var,iconst,ivar,nn,xr,q,file=file,format=format
;  Read formatted GONG file on standard form
;
;  On input, file must be set to the file name to be read from
;
;  format is an optional string parameter, defining the input format
;  for the real variables.
;    default: '(5e16.9)'
;    if set to 'free', input is without format.
;    Otherwise must be set to the appropriate format.
;
;  head returns the character header
;  glob returns the global varameters
;  var returns the variables given at each mesh point
;  iconst returns the number of global parameters
;  ivar returns the number of variables at each mesh point
;  nn returns the number of mesh points
;  xr returns r/R
;  q returns m/M

openr,1,file
head=strarr(4)
ss=''
for i=0,3 do begin $
  readf,1,ss & head(i)=ss 
endfor

for i = 0,3 do print,head(i)

if keyword_set(format) eq 0 then $
	cformat='(5e16.9)' $
else cformat=format

nn=0l
iconst=0l
ivar=0l
ivers=0l
readf,1,nn,iconst,ivar,ivers
glob=dblarr(iconst)
var=dblarr(ivar,nn)

if cformat eq 'free' then begin
  readf,1,glob 
  readf,1,var 
endif else begin 
  readf,1,glob,format=cformat 
  readf,1,var,format=cformat 
endelse

close,1
xr=var(0,*)/glob(1)
q=exp(var(1,*))
return
end
