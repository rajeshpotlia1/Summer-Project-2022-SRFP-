# Simulation Codes in Python
### Code to: Generate 100 realization full skymaps from an input $C_\ell$, Generate a mask (section:6.1)

```markdown
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
plt.ioff()
#plt.ion()


dir_maps = "/home/rajesh/"
ellgiven,dlgiven = np.loadtxt(dir_maps+'clgiven.txt',usecols=(0,1),unpack=True,skiprows=1)
clgiven2 = [0,0]
clgiven1= dlgiven*2*np.pi/ellgiven/(ellgiven+1)

clgiven = np.append(clgiven2,clgiven1)
print(clgiven)
print(np.shape(clgiven))
ell=np.arange(2509)

p1 =plt.plot(ell, ell*(ell+1)*clgiven/(2.*np.pi), linestyle='solid', color='b', label = r'Input $C_\ell$')
plt.xlabel(r'Multipole, $\ell$')
plt.ylabel(r'$\ell (\ell+1) C_{\ell}/2\pi$ [in $\mu K^2$]')
plt.legend()
plt.savefig("givenclgraph.png", format="png", bbox_inches="tight")
plt.show()
# Plotting the given Cl with l

cl=np.zeros((100,768))
# Generate 100  temperature maps using the given Cl by hp.synfast
# Then generate Cl from each of such map using hp.anafast
for i in range(100):
    map = hp.synfast(clgiven,256, fwhm = np.radians(15./60.), pixwin = True)
    cl[i] = hp.anafast(map, lmax=767)

print(np.shape(map))
avgsum=np.mean(cl,axis=0) # Find average Cl value from these 100 generated Cl

beam=hp.gauss_beam(np.radians(15./60.),lmax=767)
ell=np.arange(768)

p2 =plt.plot(ell, ell*(ell+1)*avgsum/((2.*np.pi)*(beam)**2), linestyle='solid', color='b', label = r'Recovered $C_\ell$') # Plot the avg Cl vs l
plt.xlim(xmin=3,xmax=800.)
plt.ylim(ymin=0,ymax=7000)
plt.xlabel(r'Multipole, $\ell$')
plt.ylabel(r'$\ell (\ell+1) C_{\ell}/2\pi$ [in $\mu K^2$]')
plt.legend()
plt.savefig("avgcl.png", format="png")
plt.show()


clgiven_new = clgiven[:768]
beam=hp.gauss_beam(np.radians(15./60.),lmax=767)
ell=np.arange(768)

p3 =plt.plot(ell, ell*(ell+1)*clgiven_new/(2.*np.pi), linestyle='solid', color='b', label = r'Input $C_\ell$')
p3 =plt.plot(ell, ell*(ell+1)*avgsum/((2.*np.pi)*(beam)**2), linestyle='solid', color='r', label = r'Recovered $C_\ell$') # Plot the avg Cl vs l
plt.xlim(xmin=0,xmax=800.)
plt.ylim(ymin=0,ymax=7000)
plt.xlabel(r'Multipole, $\ell$')
plt.ylabel(r'$\ell (\ell+1) C_{\ell}/2\pi$ [in $\mu K^2$]')
plt.legend()
plt.savefig("givenvsobtainedCl.png", format="png", bbox_inches="tight")
plt.show()


p_l = hp.pixwin(nside = 256 , lmax = 767)
print(np.shape(p_l))
#print(p_l)

clgiven_new = clgiven[:768]
beam=hp.gauss_beam(np.radians(15./60.),lmax=767)
ell=np.arange(768)
p4 =plt.plot(ell, ell*(ell+1)*clgiven_new/(2.*np.pi), linestyle='solid', color='b', label = r'Input $C_\ell$')
p4 =plt.plot(ell, ell*(ell+1)*avgsum/(((2.*np.pi)*(beam)**2)*((p_l)**2)), linestyle='solid', color='r', label = r'Recovered $C_\ell$') # Plot the avg Cl vs l
plt.xlim(xmin=0,xmax=800.)
plt.ylim(ymin=0,ymax=11000)
plt.xlabel(r'Multipole, $\ell$')
plt.ylabel(r'$\ell (\ell+1) C_{\ell}/2\pi$ [in $\mu K^2$]')
plt.legend()

plt.savefig("givenvsobtainedClafterPleffect.png", format="png")

plt.show()

#To find variance from these 100 Cl and comparing it with anlytical formula.
std_cl = np.std(cl,axis=0) # To find standard deviation
#print(std_cl)

# Calculation of std by analytical formula
import math
expstd_cl = []
for h in range(768):
  expstd_cl.append( avgsum[h] * math.sqrt(2/((2*h)+1)))

#print(expstd_cl)
  
# Comparing the both standard deviations by plotting them
ell=np.arange(768)
p5 =plt.plot(ell, std_cl, linestyle='solid', color='b', label = 'Std calculated from 100 realiations')
p5 =plt.plot(ell, expstd_cl, linestyle='solid', color='r', label = 'Expected std by analytical formula')
plt.xlim(xmin=0,xmax=50.)
plt.ylim(ymin=0,ymax=1300)
plt.xlabel(r'Multipole, $\ell$')
plt.ylabel(r'Standard deviation')
plt.legend()
plt.show()
plt.savefig("givenvsobtainedCl_variance.png", format="png")

#Masking a region of the sky and then obtaining a power spetrum
cl=np.zeros((100,768))
maps= []
# Generate 100  temperature maps using the given Cl by hp.synfast
# Then generate Cl from each of such map using hp.anafast
for i in range(100):
    maps.append(hp.pixelfunc.reorder(hp.synfast(clgiven,256, fwhm = np.radians(15./60.), pixwin = True), inp=None, out=None, r2n=True, n2r=None))

print(np.shape(maps))
print(maps[0])
hp.mollview(maps[0], cbar=True,nest=True,norm='hist')
'''
#MASK1

mask1 = np.ones(786432)

for j in range(786432):

    theta,phi = hp.pix2ang(nside = 256 , ipix = j , nest=True , lonlat = True)
    if theta < 30 :
      mask1[j] = 0
    else:
      mask1[j] = 1
      
n_zeros = np.count_nonzero(mask1==0)
# display the count of zeros
print(n_zeros)

# Masking all temp maps obtained from given C_l
maskedtempmaps1 = []
for i in range(100):
  letproduct = np.multiply(maps[i],mask1)
  maskedtempmaps1.append(letproduct)
print(np.shape(maskedtempmaps1))

#Now, we have 100 temp maps masked with mask1


#MASK2

mask2 = np.ones(786432)

for j in range(786432):

    theta,phi = hp.pix2ang(nside = 256 , ipix = j , nest=True , lonlat = True)
    if phi > 45  :
      mask2[j] = 0
    else:
      mask2[j] = 1


n_zeros = np.count_nonzero(mask2==0)
# display the count of zeros
print(n_zeros)


# Masking all temp maps obtained from given C_l
maskedtempmaps2 = []
for i in range(100):
  letproduct = np.multiply(maps[i],mask2)
  maskedtempmaps2.append(letproduct)
print(np.shape(maskedtempmaps2))

#Now, we have 100 temp maps masked with mask2
'''
#MASK3
# Let us create another masking of temp maps
mask3 = np.ones(786432)

for j in range(786432):

    theta,phi = hp.pix2ang(nside = 256 , ipix = j , nest=True , lonlat = True)
    if -15 < phi < 15  :
      mask3[j] = 0
    else:
      mask3[j] = 1

n_zeros = np.count_nonzero(mask3==0)
# display the count of zeros
print(n_zeros)



'''
#MASK4
# Let us create another masking of temp maps
mask4 = np.ones(786432)
for j in range(786432):

    theta,phi = hp.pix2ang(nside = 256 , ipix = j , nest=True , lonlat = True)
    if phi < 0  :
      mask4[j] = 0
    else:
      mask4[j] = 1


#n_zeros = np.count_nonzero(mask4==0)
# display the count of zeros
#print(n_zeros)

# Masking all temp maps obtained from given C_l
maskedtempmaps4 = []
for i in range(100):
  maskedtempmaps4.append(np.multiply(maps[i],mask4))
  #maskedtempmaps4.append(letproduct)
#print(np.shape(maskedtempmaps4))

#Now, we have 100 temp maps masked with mask4
'''
fitsfile = []
for i in range(100):
  ok = "mapinfits3M"+str(i)
  notok = ok + str(".fits")
  fitsfile.append(notok)


#Creating a FITS file from an array map
#fitsfile = ["mapinfits1.fits","mapinfits2.fits","mapinfits3.fits","mapinfits4.fits","mapinfits5.fits","mapinfits6.fits","mapinfits7.fits","mapinfits8.fits","mapinfits9.fits","mapinfits10.fits"]

for i in range(100):
  hp.fitsfunc.write_map(fitsfile[i], maps[i], nest=True, dtype=None, fits_IDL=True, coord=None, partial=False, column_names=None, column_units=None, extra_header=(), overwrite=False)
hp.fitsfunc.write_map("mask3infits.fits", mask3, nest=True, dtype=None, fits_IDL=True, coord=None, partial=False, column_names=None, column_units=None, extra_header=(), overwrite=False)

```
### Code to apodize a mask (section:6.2)


```markdown
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
m100=hp.read_map("/home/rajesh/mask3infits.fits", nest = True)
sm100 = hp.smoothing(m100, fwhm=np.radians(np.sqrt(5.**2-(5./60.)**2)), nest = True)
hp.fitsfunc.write_map("apodmask3infits.fits", sm100, nest=True, dtype=None, fits_IDL=True, coord=None, partial=False, column_names=None, column_units=None, extra_header=(), overwrite=False)
```

### PolSpice code (section:6.2.1)
```markdown

import sys
sys.path.append('/home/rajesh/polspice/PolSpice_v03-07-02/src/')
import ispice
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

fitsfile = []
for i in range(100):
  ok = "/home/rajesh/apodmapsinfitspw/mapinfits3M"+str(i)
  notok = ok + str(".fits")
  fitsfile.append(notok)


clcutsky100 = []

for i in fitsfile:
    ispice.ispice(i , 'cls_tmp0.fits', nlmax=767, weightfile1='/home/rajesh/apodmapsinfitspw/apodmask3infits.fits',beam1=np.radians(15./60.),label='spice', polarization='NO', thetamax = 120.,apodizesigma=120.,apodizetype=1)
    tmp=hp.read_cl('cls_tmp0.fits')
    clcutsky100.append(tmp)

avgcl = np.mean(clcutsky100,axis=0) # Find average Cl value from these 100 generated Cl


fitsfile1 = []
for i in range(100):
  ok = "/home/rajesh/masked3maps100/mapinfits3M"+str(i)
  notok = ok + str(".fits")
  fitsfile1.append(notok)


clcutsky = []

for i in fitsfile:
    ispice.ispice(i , 'cls_tmp0.fits', nlmax=767, maskfile1='/home/rajesh/masked3maps100/mask3infits.fits',beam1=np.radians(15./60.),label='spice', polarization='NO', thetamax = 120.)
    tmp=hp.read_cl('cls_tmp0.fits')
    clcutsky.append(tmp)

avgcl1 = np.mean(clcutsky,axis=0) # Find average Cl value from these 100 generated Cl



ell=np.arange(768)
beam=hp.gauss_beam(np.radians(15./60.),lmax=767)


dl= ell*(ell+1)*avgcl/(((2.*np.pi)*(beam)**2))
dl1= ell*(ell+1)*avgcl1/(((2.*np.pi)*(beam)**2))



dir_maps = "/home/rajesh/"
ellgiven,dlgiven = np.loadtxt(dir_maps+'clgiven.txt',usecols=(0,1),unpack=True,skiprows=1)
clgiven2 = [0,0]
clgiven1= dlgiven*2*np.pi/ellgiven/(ellgiven+1)
clgiven = np.append(clgiven2,clgiven1)
clgiven_new = clgiven[:768]



p1 =plt.plot(ell,dl,linestyle='solid',color='r',label='Recovered $C_\ell$ without apodization')
p2 =plt.plot(ell,dl1,linestyle='solid',color='b',label='Recovered $C_\ell$ with apodization') 
p4 =plt.plot(ell, ell*(ell+1)*clgiven_new/(2.*np.pi), linestyle='solid', color='k', label = 'Input $C_\ell$')





plt.xlim(xmin=3,xmax=600.)
plt.ylim(ymin=200,ymax=10000)
plt.xlabel(r'Multipole, $\ell$')
plt.ylabel(r'$\ell (\ell+1) C_{\ell}/2\pi$ [in $\mu K^2$]')
plt.legend()
plt.show()







```


### XPol code (section:6.2.2)

```markdown
import sys
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import healpy as hp

path_bin   = '/home/tghosh/softwares/Xpol/Xpol-master-3b64027b3d53aa03ff8ca3773f8ba3a69b6d506b/'
path_simu  = 'outputs/'

fitsfile = []
for i in range(100):
  ok = "mapinfits3M"+str(i)
  notok = ok + str(".fits")
  fitsfile.append(notok)

dlcutsky100 = []

for i in fitsfile:

            nside    = 256
            nstokes  = 1
            nprocs   = 64
            lmax     = 767
            nmaps    = 1

            path_maps   = '/home/rajesh/apodmapsinfitspw/'
            file_mask = "/home/rajesh/apodmapsinfitspw/apodmask3infits.fits"


            f = open("%sxpol_%s.par" % (path_simu,'cl'), "w")
            f.write( "nside    = %d\n" % nside)
            f.write( "lmax     = %d\n" % lmax)
            f.write( "nstokes  = %d\n" % nstokes)
            f.write( "nmaps    = %d\n" % nmaps)
            f.write( f"mapfile1 = %s{i}\n" % path_maps)
            f.write( "weightI1 = %s\n" % file_mask)
            f.write( "cross    = %s%s\n"  % (path_simu,'cl'))
            f.close()
            os.system( "mpirun -n %d %sxpol %sxpol_cl.par" % (nprocs,path_bin,path_simu))


            ##outputcl file
            filename='%scl_0_0.fits' %(path_simu)
            f=fits.open(filename)
            datTT=f[1].data
            ell  = datTT.field(0)
            dlTT = datTT.field(1)
            print(np.max(ell))
            dlcutsky100.append(dlTT)

dlcutsky = []
for i in dlcutsky100:
    let = [0,0]
    dlcutsky.append(np.append(let,i))




bindl=[]

for i in dlcutsky:
  newarr = np.array_split(i, 64)
  newavgarr = []
  for j in newarr:
    newavgarr.append(np.mean(j))
  bindl.append(newavgarr)






avgdl = np.mean(bindl,axis=0) # Find average Cl value from these 100 generated Cl
std_dl = np.std(bindl,axis=0)


ell=np.arange(768)
newarr1 = np.array_split(ell, 64)
binell2 = []
for i in newarr1:
  avg1 = np.mean(i)
  binell2.append(avg1)
binell = np.array(binell2)


beam=hp.gauss_beam(np.radians(15./60.),lmax=767)
newarr2 = np.array_split(beam, 64)
binbeam2 = []
for i in newarr2:
  avg2 = np.mean(i)
  binbeam2.append(avg2)
binbeam = np.array(binbeam2)

dl=avgdl/((binbeam)**2)

dir_maps = "/home/rajesh/"
ellgiven,dlgiven = np.loadtxt(dir_maps+'clgiven.txt',usecols=(0,1),unpack=True,skiprows=1)
clgiven2 = [0,0]
clgiven1= dlgiven*2*np.pi/ellgiven/(ellgiven+1)
clgiven = np.append(clgiven2,clgiven1)
clgiven_new = clgiven[:768]



p1 =plt.plot(binell,dl,linestyle='solid',color='r',label='Recovered $C_\ell$')
p2 =plt.plot(ell, ell*(ell+1)*clgiven_new/(2.*np.pi), linestyle='solid', color='k', label = 'Input $C_\ell$')
p3 =plt.errorbar(binell,dl,yerr = std_dl , ecolor = 'g',elinewidth = 3)





plt.xlim(xmin=0.,xmax=600.)
plt.ylim(ymin=0.,ymax=8000.)
plt.xlabel(r'Multipole, $\ell$')
plt.ylabel(r'$\ell (\ell+1) C_{\ell}/2\pi$ [in $\mu K^2$]')
plt.legend()
plt.show()
plt.savefig("givenvsobtainedCl.png", format="png")
plt.savefig("givenvsobtainedCl.pdf", format="pdf")





```
