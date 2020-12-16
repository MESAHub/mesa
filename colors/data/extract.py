import numpy as np

#Download the photon weigthed values from http://vizier.cfa.harvard.edu/viz-bin/VizieR-3?-source=J/A%2bAS/130/65/lcb98cor

file_in="lcb98cor_download.data"
file_out="lcb98cor.dat"


x=np.genfromtxt(file_in,names=True)

# Diff in Mag Say B-V implies -BCb+BCv

V=x['BCVmagp']

B=V-x['BVp']
U=B-x['UBp']

R=x['VRp']+V
I=x['VIp']+V
K=x['VKp']+V

H=K-x['HKp']
J=H-x['JHp']
L=x['JLp']+J
Lprime=x['JLp_1']+J
M=x['KMp']+K

with open(file_out,'w') as f:
    print("#Teff logg M_div_H   U  B  V  R  I  J  H  K  L  Lprime  M",file=f)
    #for y in zip(x['Teff'],x['logg'],x['FeH'],U,B,V,R,I,K,H,J,L,Lprime,M):
    for y in zip(x['Teff'],x['logg'],x['FeH'],U,B,V,R,I,J,K,H,L,Lprime,M):
        print(" ".join([str(z) for z in y]),file=f)




