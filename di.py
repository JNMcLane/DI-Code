"""
di.py

Written by Jacob N. McLane

Code for doing basic analysis of direct imaging data.
Written for Planetary Astrophysics course.

No inputs, just a few booleans to toggle.
Program assumes you have a you list of
fits files sorted by target and stored as
'target'.txt

Booleans

basic:
program runs basic shifting + stack and
shift + rotate + stack routine. Produces
median averaged and summed fits files.

azm :
Program preforms median averaging using
concentric anulii of 1 pixel, and subtracts
from the image before shifting + stacking
and shifting + rotating + stacking routine.
Produces median averaged and summed fits
files.

median:
Creates a Master Median file for the target, which
is then subtracted from each observation before
rotating and stacking. Produces median averaged
and summed fits files

gFWHM:
Generates FWHM map of files for use in FWHMd.

FWHMd:
Subtracts the file with the most similar psf from
each file before shifting, stacking, and rotating.
Produces median averaged fits files.

Program utilizes cntrd.py, a python recreation of
a IDL program. cntrd.py was translated by
D. Jones.

Code was written and tested with Python 3.4.1.
Compatibility not gauranteed  for any other
release
"""
###Import Modules###
import numpy as np
from math import *
from cntrd import *
import scipy.ndimage.interpolation as inter
import os
from scipy.interpolate import UnivariateSpline
import time
from astropy.io import fits

start_time = time.time() #Timer to time program

#Set Target
#tar='ROXs12'
tar='ROXs42b'

#Choose type of reduction
basic = False
azm = False
median = False
gFWHM = False
FWHMd = True

#Load fits file
fstar=np.loadtxt(tar+'.txt',dtype=str)


for i in range(len(fstar)):
    fstar[i]=fstar[i].strip("b''")

log=open(tar+'_log.txt','w')


PAR=[]
ROT=[]
EL=[]
INS=[]
#Use centroid to find center of star and write data to file
for i in range(len(fstar)):
    hdulists=fits.open(fstar[i])
    sdat=hdulists[0].data
    a=cntrd(sdat,612,472,65)
    #Pull Position angle parameters from header and write to table
    PAR.append(hdulists[0].header['PARANG'])
    ROT.append(hdulists[0].header['ROTPPOSN'])
    EL.append(hdulists[0].header['EL'])
    INS.append(hdulists[0].header['INSTANGL'])
    print(fstar[i],a[0],a[1], file=log)

log.close()

#Load fits file w/ star center parameters
a=np.loadtxt(tar+'_log.txt', dtype={'names':('file','xcen','ycen'),'formats' : ('S26','f16','f16')})

image=a['file']
xcen=a['xcen']
ycen=a['ycen']

b=[]

for i in range(len(image)):
    b.append(str(image[i]))
    b[i]=b[i].strip("b''")

#Shifting + Shift and Rotate Image processing
if basic == True:
    
    data=[]
    rot=[]

    #open files
    for i in range(len(b)):
        img=fits.open(b[i])
        dat=img[0].data
        xshift=511.5-xcen[i]
        yshift=511.5-ycen[i]
        PA=-(PAR[i]+ROT[i]-EL[i]-INS[i])
        ndat=inter.shift(dat,(yshift,xshift))
        rdat=inter.rotate(ndat, PA, reshape=False)
        data.append(ndat)
        rot.append(rdat)

    #Find sum and median
    dsum=np.sum(data, axis=0)
    dmed=np.median(data, axis=0)

    #Write to .fits files
    #Shifted and stacked images
    hdu1 = fits.PrimaryHDU(dsum)
    hdulist = fits.HDUList([hdu1])
    hdulist.writeto(tar+'_sum.fits' ,clobber=True)

    hdu2 = fits.PrimaryHDU(dmed)
    hdulist = fits.HDUList([hdu2])
    hdulist.writeto(tar+'_median.fits' ,clobber=True)


    dsumr=np.sum(rot, axis=0)
    dmedr=np.median(rot, axis=0)

    #Shifted, Rotated, and Stacked images
    hdu3 = fits.PrimaryHDU(dsumr)
    hdulistr = fits.HDUList([hdu3])
    hdulistr.writeto(tar+'_sum_rot.fits' ,clobber=True)

    hdu4 = fits.PrimaryHDU(dmedr)
    hdulistr = fits.HDUList([hdu4])
    hdulistr.writeto(tar+'_median_rot.fits' ,clobber=True)

#Do median subtraction over annulus
if azm == True:
    rlog=open(tar+'radial_prof.txt','w')
    data=[]
    rot=[]

    for i in range(len(b)):
        img=fits.open(b[i])
        dat=img[0].data
        xshift=511.5-xcen[i]
        yshift=511.5-ycen[i]
        PA=-(PAR[i]+ROT[i]-EL[i]-INS[i])

#Actual median subtraction over anulii done here
        np.indices(dat.shape)
        x,y=np.indices(dat.shape)
        r = np.hypot(x - int(xcen[i]), y - int(ycen[i]))
        nbins = (np.round(r.max())+1)
        bins = np.linspace(0,nbins,nbins+1)
        bin_centers = (bins[1:]+bins[:-1])/2.0
        whichbin = np.digitize(r.flat,bins)
        nr = np.bincount(whichbin)[1:]
        r_proc=[]
        for k in range(1,int(nbins+1)):
            r_proc.append(np.median(np.array((dat).flat[whichbin==k])))
        for m in range(1024):
            for n in range(1024):
                ch=int(r[m,n])
                dat[m,n]=dat[m,n]-r_proc[ch]
            if i ==0:
                for h in range(len(bin_centers)):
                    print(int(bin_centers[h]),r_proc[h], file=rlog)
###########################################

        print('starting shift/rotation')
        ndat=inter.shift(dat,(yshift,xshift))
        rdat=inter.rotate(ndat, PA, reshape=False)
        data.append(ndat)
        rot.append(rdat)
        print('finished shift/rotation')
        print('Finished loop '+str(i+1)+' of '+str(len(b))+ '.')
        
    rlog.close()
    dsum=np.sum(data, axis=0)
    dmed=np.median(data, axis=0)

    #Produce image files
    hdu1 = fits.PrimaryHDU(dsum)
    hdulist = fits.HDUList([hdu1])
    hdulist.writeto(tar+'_rad_sum.fits' ,clobber=True)

    hdu2 = fits.PrimaryHDU(dmed)
    hdulist = fits.HDUList([hdu2])
    hdulist.writeto(tar+'_rad_median.fits' ,clobber=True)


    dsumr=np.sum(rot, axis=0)
    dmedr=np.median(rot, axis=0)


    hdu3 = fits.PrimaryHDU(dsumr)
    hdulistr = fits.HDUList([hdu3])
    hdulistr.writeto(tar+'_rad_sum_rot.fits' ,clobber=True)

    hdu4 = fits.PrimaryHDU(dmedr)
    hdulistr = fits.HDUList([hdu4])
    hdulistr.writeto(tar+'_rad_median_rot.fits' ,clobber=True)

#Do subtraction with Master median
if median == True:
    print('Starting Median Subtraction')
    
    data=[]

    #Generate master median

    for i in range(len(b)):
        img=fits.open(b[i])
        dat=img[0].data
        xshift=511.5-xcen[i]
        yshift=511.5-ycen[i]
        ndat=inter.shift(dat,(yshift,xshift))
        data.append(ndat)

    dmed=np.median(data, axis=0)

    print('Median Image Generated')
    ##########################
    data=[]
    rot=[]

    print('Beginning Loop')
    
    for i in range(len(b)):
        img=fits.open(b[i])
        dat=img[0].data
        xshift=511.5-xcen[i]
        yshift=511.5-ycen[i]
        PA=-(PAR[i]+ROT[i]-EL[i]-INS[i])
        ndat=inter.shift(dat,(yshift,xshift))
        ndat=ndat-dmed
        rdat=inter.rotate(ndat, PA, reshape=False)
        data.append(ndat)
        rot.append(rdat)
        print('Finished  Loop '+str(i+1)+' of '+str(len(b))+ '.')

        #Make .fits files
    print('Making .fits files')
    dsum=np.sum(data, axis=0)
    dmed=np.median(data, axis=0)

    hdu1 = fits.PrimaryHDU(dsum)
    hdulist = fits.HDUList([hdu1])
    hdulist.writeto(tar+'_med_sum.fits' ,clobber=True)

    hdu2 = fits.PrimaryHDU(dmed)
    hdulist = fits.HDUList([hdu2])
    hdulist.writeto(tar+'_med_median.fits' ,clobber=True)


    dsumr=np.sum(rot, axis=0)
    dmedr=np.median(rot, axis=0)


    hdu3 = fits.PrimaryHDU(dsumr)
    hdulistr = fits.HDUList([hdu3])
    hdulistr.writeto(tar+'_med_sum_rot.fits' ,clobber=True)

    hdu4 = fits.PrimaryHDU(dmedr)
    hdulistr = fits.HDUList([hdu4])
    hdulistr.writeto(tar+'_med_median_rot.fits' ,clobber=True)

    print('Finished')

    #Generate table of FWHM values

if gFWHM == True:

    logf=open(tar+'_fwhm.txt','w')
    p=np.loadtxt(tar+'_planet.txt', dtype={'names':('file','xcor','ycor'),'formats' : ('S26','f16','f16')})

    fit=p['file']
    x=p['xcor']
    y=p['ycor']

    f=[]

    for i in range(len(fit)):
        f.append(str(fit[i]))
        f[i]=f[i].strip("b''")
        
    mfwhm=[]
    for j in range(len(p)):
        img=fits.open(f[j])
        dat=img[0].data
        
        x=[]
        for i in range(1024):
            x.append(i)
    
        row=int(y[j])
        x1=dat[row]
        x2=dat[row+1]
        x3=dat[row-1]
        vec=[]
    
        for i in range(1024):
            vec.append(x1+x2+x3)
            
        nvec=np.median(vec,axis=0)   	
        spline = UnivariateSpline(x, nvec-np.max(nvec)/2, s=0)
        if len(spline.roots()) == 2:
            r1, r2 = spline.roots()
            fwhm=r2-r1

        elif len(spline.roots()) == 4:
            r1, r2, r3, r4 = spline.roots()
            fwhm=r3-r2

        mfwhm.append(fwhm)
    
    for j in range(len(p)):
        mat=mfwhm
        val=mat[j]
        mat[j]=0
        for h in range(len(mat)):
            mat[h]=abs(mat[h]-val)
        comp=np.argmin(mat)
        p1=str(f[j])
        p2=str(f[comp])
        print(p1,p2, file=logf)

    logf.close()

    #Use table of FWHM values to find most similar image and subtract
if FWHMd == True:

    p=np.loadtxt(tar+'_fwhm.txt', dtype={'names':('file1','file2'),'formats' : ('S26','S26')})

    cal=p['file2']
    f=[]

    for i in range(len(cal)):
        f.append(str(cal[i]))
        f[i]=f[i].strip("b''")
    
    data=[]
    rot=[]

    for i in range(len(b)):
        img=fits.open(b[i])
        dat=img[0].data
        cali=fits.open(f[i])
        cald=cali[0].data
        xshift=511.5-xcen[i]
        yshift=511.5-ycen[i]
        PA=-(PAR[i]+ROT[i]-EL[i]-INS[i])
        dat=dat-cald
        ndat=inter.shift(dat,(yshift,xshift))
        rdat=inter.rotate(ndat, PA, reshape=False)
        data.append(ndat)
        rot.append(rdat)

    dsumr=np.sum(rot, axis=0)
    dmedr=np.median(rot, axis=0)

#Produce files

    hdu3 = fits.PrimaryHDU(dsumr)
    hdulist = fits.HDUList([hdu3])
    hdulist.writeto(tar+'_cal_sum_rot.fits' ,clobber=True)

    hdu4 = fits.PrimaryHDU(dmedr)
    hdulistr = fits.HDUList([hdu4])
    hdulistr.writeto(tar+'_cal_median_rot.fits' ,clobber=True)

 #Print time to run code       
print("--- %s seconds ---" % (time.time() - start_time))

