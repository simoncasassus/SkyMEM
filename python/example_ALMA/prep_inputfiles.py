import sys
import os
import numpy as np
from casatasks import tclean, exportfits, imhead, importfits, immath, imstat
from casatools import image
from casatools import quanta
from astropy.io import fits


def getdirty(
        input_ms="",
        imsize=512,
        cellsize='0.05arcsec',
        stokes="I",
        workdir="./skymem/",
        rmsnoise=9.801996e-06,  # Jy/beam thermal noise
        weighting="robust",
        robust=2.0):
    ia = image()
    qa = quanta()
    os.system("rm -rf  " + workdir)
    os.system("mkdir " + workdir)

    basename_dirty_image = workdir + os.path.basename(input_ms)
    print("basename_dirty_image", basename_dirty_image)
    psf = workdir + os.path.basename(input_ms) + ".psf"
    PB = workdir + os.path.basename(input_ms) + ".flux"

    os.system("rm -rf *.log *.last " + basename_dirty_image + ".*  ")

    #importfits(imagename="model_out", fitsimage=model_fits, overwrite=True)
    #cdelt2 = imhead(imagename="model_out", mode="get", hdkey="cdelt2")["value"]
    #cdelt_string = str(qa.convert(v=cdelt2,
    #outunit="arcsec")["value"]) + "arcsec"
    #size = imhead(imagename="model_out", mode="get", hdkey="shape")

    tclean(vis=input_ms,
           imagename=basename_dirty_image,
           specmode='mfs',
           deconvolver='hogbom',
           niter=0,
           stokes=stokes,
           nterms=1,
           weighting=weighting,
           robust=robust,
           imsize=imsize,
           cell=cellsize,
           datacolumn='data')

    exportfits(imagename=basename_dirty_image + ".image",
               fitsimage=basename_dirty_image + ".dirty.fits",
               overwrite=True,
               history=False)

    exportfits(imagename=basename_dirty_image + ".psf",
               fitsimage=basename_dirty_image + ".psf.fits",
               overwrite=True,
               history=False)

    exportfits(imagename=basename_dirty_image + ".pb",
               fitsimage=basename_dirty_image + ".pb.fits",
               overwrite=True,
               history=False)

    hdu0 = fits.open(basename_dirty_image + ".pb.fits")
    hdr0 = hdu0[0].header
    im0 = hdu0[0].data
    sensitivity = rmsnoise * (1. / im0)
    hdr0['BUNIT'] = '1/Jy/beam'
    hdu0[0].data = sensitivity
    hdu0[0].header = hdr0
    hdu0.writeto(basename_dirty_image + ".sens.fits", overwrite=True)

    #ia.open(infile=dirty_casa_image)
    #record_beam = ia.restoringbeam()
    #ia.done()
    #ia.close()
    #
    #image_name_list = ["convolved_model_out", dirty_casa_image + ".fits"]
    #
    #immath(imagename=image_name_list,
    #       expr=" (IM0   + IM1) ",
    #       outfile=restored_image,
    #       imagemd=dirty_casa_image + ".fits")
    #
    #exportfits(imagename=restored_image,
    #           fitsimage=restored_image + ".fits",
    #           overwrite=True,
    #           history=False)
    #
    #peak = imstat(imagename=restored_image)["max"][0]
    #rms = imstat(imagename=dirty_casa_image)["rms"][0]
    #psnr = peak / rms
    return basename_dirty_image


def scaleprior(file_prior, workdir="./", basename_dirty_image=""):
    hdu0 = fits.open(basename_dirty_image + ".dirty.fits")
    hdr0 = hdu0[0].header
    omegabeam = (np.pi/(4.*np.log(2)))*hdr0['BMAJ']*hdr0['BMIN']
    hdu = fits.open(file_prior)
    prior = hdu[0].data
    hdr= hdu[0].header
    omegapixel=hdr['CDELT2']**2
    prior *= omegabeam/omegapixel
    hdr['BMAJ']=hdr0['BMAJ']
    hdr['BMIN']=hdr0['BMIN']
    hdr['BPA']=hdr0['BPA']
    hdr['BUNIT']=hdr0['Jy/beam']
    hdu[0].data=prior
    hdu[0].header=hdr
    file_prior_scaled = workdir+'prior.fits'
    hdu.writeto(file_prior_scaled, overwrite=True)


input_ms = '/home/simon/debris_SM/data/corrected/TYC9340/TYC9340.12m.cal.bary.continuum.fav.tav.SMGsub.corrected.ms'
workdir = "./skymem/"
basename_dirty_image = getdirty(input_ms=input_ms,
                                imsize=512,
                                cellsize='0.05arcsec',
                                stokes="I",
                                workdir=workdir,
                                weighting="briggs",
                                robust=2.0)

prior = '/home/simon/debris_SM/guvmem_runs/ACA_corrected/mem_lS0.0_lL0.0_nogrid_nc1.2_niter3_finepix/mod_out.fits'

scaleprior(prior, workdir=workdir, basename_dirty_image=basename_dirty_image)

#print("PSNR: {0:0.3f}".format(psnr))
#print("Peak: {0:0.3f} mJy/beam".format(peak * 1000.0))
#print("RMS: {0:0.3f} mJy/beam".format(rms * 1000.0))
