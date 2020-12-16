import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, interpolate
import astropy.io.fits as fits

c=299792458.0 #m/s
teffsun=5777.0 #k
rsun=6.95700*10**8 #m
h=6.6260693*10**(-34) # j s
k=1.380658*10**(-23) # J/K
pc=3.0857*10**16 # m
stefan=(2.0*(np.pi**5)*(k**4))/(15*(c**2)*(h**3)) # W k^-4 m^-2
lsun=4.0*np.pi*stefan*rsun**2*teffsun**4 # j/s

#Vega defines the colours to be 0
magZero=0.03
mbolsun=4.77
Rdsq=6.247*10**-17

vega=fits.open('vega.fits')

vega_wave=vega[1].data.WAVELENGTH*10**-10
vega_flux=vega[1].data.FLUX

#Convert to si 1 FLAM is 10^-7 W m^-2 m^-1
vega_flux=vega_flux*10**7

vega_interp=interpolate.interp1d(vega_wave,vega_flux,bounds_error=False,fill_value=0.0)


def do_one(f_in,vega_interp,teff):
	j=np.genfromtxt(f_in,names=["wave_nm",'flux'])
      
	j_m=j['wave_nm']*10**-10 # meters
	j_f=j['flux']
	j_f_interp=interpolate.interp1d(j_m,j_f,bounds_error=False,fill_value=0.0)

	j_norm=integrate.romberg(j_f_interp,j_m[0],j_m[-1],divmax=50,vec_func=True,tol=1e-10,rtol=1e-10)


	def bb_wave(wave,teff):
		return ((2*np.pi*h*c**2)/(wave**5))*(1.0/(np.exp((h*c)/(wave*k*teff))-1.0))


	def const(t):
		return mbolsun-2.5*np.log10((4.0*np.pi*(10.0*pc)**2*stefan*t**4)/lsun)-magZero

	def topfunc(wave,t):
		return wave*(bb_wave(wave,t))*(j_f_interp(wave)/j_norm)
	def botfunc(wave):
		return wave*(vega_interp(wave))*(j_f_interp(wave)/j_norm)	


	def bc(t):
		constants=const(t)
		top=integrate.romberg(topfunc,j_m[0],j_m[-1],args=(t,),divmax=50,vec_func=True,tol=1e-10,rtol=1e-10)
		bot=integrate.romberg(botfunc,j_m[0],j_m[-1],divmax=50,vec_func=True,tol=1e-10,rtol=1e-10)
		print(t,top,bot)
		return constants+2.5*np.log10(top/bot)

	y=[]
	for i in teff:
		y.append(bc(i))

	return np.array(y)

teff=np.arange(100.0,50001.0,100.0)

base='bess'
u=do_one(base+'-u.dat',vega_interp,teff)
b=do_one(base+'-b.dat',vega_interp,teff)
v=do_one(base+'-v.dat',vega_interp,teff)
r=do_one(base+'-r.dat',vega_interp,teff)
i=do_one(base+'-i.dat',vega_interp,teff)


with open('blackbody_johnson.dat','w') as f:
	print("#Teff logg mdivh bb_U bb_B bb_V bb_R bb_I",file=f)
	for i in zip(teff,u,b,v,r,i):
		print(i[0],0.0,0.0,i[1],i[2],i[3],i[4],i[5],file=f)


		
