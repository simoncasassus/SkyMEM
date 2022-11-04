import sys
import numpy as np
import os
import os.path
from scipy import ndimage
from astropy.io import fits
import re
from scipy.signal import convolve as scipy_convolve
import scipy.optimize as op
import time
from copy import deepcopy

include_path = '/home/simon/common/python/include/'
sys.path.append(include_path)
from ImUtils.Resamp import gridding
from ImUtils.Cube2Im import slice0
import Gausssmooth
import PyVtools.Vtools as Vtools

#from astropy.wcs import WCS
#from pylab import *


def pass_setup(ZMpass):
    global ZM
    ZM = ZMpass


def exec_prep_files(ZM):

    os.system("rm -rf " + ZM.workdir)
    os.system("mkdir " + ZM.workdir)

    datapass = ZM.datapass
    reffield = ZM.datafieldref  # datapass[ZM.ifieldref]

    hdu0_fieldref = fits.open(reffield['image'])
    hdu1_fieldref = slice0(hdu0_fieldref, ReturnHDUList=True)
    ZM.canvas_shape = hdu1_fieldref[0].data.shape
    hcanvas = hdu1_fieldref[0].header

    if (reffield['headfile']):
        hdu0 = fits.open(reffield['headfile'])
    else:
        hdu0 = fits.open(reffield['image'])

    hrstor = hdu0[0].header
    reffield['rstorheader'] = hrstor
    ZM.rstor_head = reffield['rstorheader']
    beam_deg2 = (np.pi / (4. * np.log(2.))) * (hrstor['BMAJ'] * hrstor['BMIN'])
    ZM.beam_deg2 = beam_deg2
    ZM.pixscale = hrstor['CDELT2']

    sensmos = np.zeros(ZM.canvas_shape)

    for afield in datapass:
        hdu0 = slice0(fits.open(afield['image']), ReturnHDUList=True)
        hdu1 = gridding(hdu0, hcanvas, ReturnHDUList=True)
        afield['data_image'] = np.double(hdu1[0].data)

        hdu0 = slice0(fits.open(afield['sens']), ReturnHDUList=True)
        sens0 = np.double(hdu0[0].data)

        sens0[np.isnan(sens0)] = -1
        sens0[(sens0 <= 0.)] = -1
        afield['sigma'] = np.min(sens0[sens0 > 0.])

        if (ZM.pbcutoff > 0.):
            sens0[(sens0 > (afield['sigma'] * ZM.pbcutoff))] = -1

        hdu0[0].data = sens0
        hdu1 = gridding(hdu0, hcanvas, ReturnHDUList=True)
        sens = hdu1[0].data
        sens[(sens <= 0.)] = -1

        PB = np.double(1. / sens)
        PB[(PB < 0.)] = 0.

        #hdu0[0].data=PB0
        #hdu1=gridding(hdu0,hcanvas,ReturnHDUList=True)
        #PB=hdu1[0].data
        PB /= np.max(PB)

        if ZM.ViewInit:
            print("original PB")
            Vtools.View(PB)

        if ZM.GaussPB:
            hdr_sens = hdu0[0].header
            (nx, ny) = sens0.shape
            x = np.arange(0, nx)
            y = np.arange(0, ny)
            XX, YY = np.meshgrid(x, y)
            rrs = np.sqrt((XX - hdr_sens['CRPIX1'] + 1)**2 +
                          (YY - hdr_sens['CRPIX2'] + 1)**2)
            stdev = (0.0193 / (np.sqrt(2. * np.log(2)))) / hdr_sens['CDELT2']
            GaussPB0 = np.exp(-0.5 * (rrs**2 / stdev**2))
            hdu0[0].data = GaussPB0
            hdu1 = gridding(hdu0, hcanvas, ReturnHDUList=True)
            if ZM.PBextend:
                PB1 = hdu1[0].data
                if (ZM.pbcutoff <= 0.):
                    sys.exit("need a pbcutoff to extend the PB")
                continuityvalue = np.min(PB[(PB > 1. / ZM.pbcutoff)])
                continuityvalue1 = np.min(PB1[(PB > 1. / ZM.pbcutoff)])
                #continuityvalue=np.max(PB[(PB<=1./ZM.pbcutoff)])
                #continuityvalue1=np.max(PB1[(PB<=1./ZM.pbcutoff)])
                #print("continuityvalue",continuityvalue)
                #print("continuityvalue1",continuityvalue1)
                PB[(PB <= 1. / ZM.pbcutoff)] = continuityvalue * PB1[
                    (PB <= 1. / ZM.pbcutoff)] / continuityvalue1
            else:
                PB = hdu1[0].data
            if ZM.ViewInit:
                print("Gauss PB")
                Vtools.View(PB)

        fits.writeto(ZM.workdir + "PB_" + os.path.basename(afield['image']),
                     PB,
                     hrstor,
                     overwrite=True)

        afield['data_sens'] = np.double(sens)
        afield['data_PB'] = np.double(PB)

        dum1 = 1. / sens**2
        dum1[(dum1 <= 0.)] = 0.
        sensmos += dum1

        hdu0 = slice0(fits.open(afield['psf']), ReturnHDUList=True)
        psf = hdu0[0].data
        afield['data_psf'] = np.double(psf)

        knorm = ZM.beam_deg2 / (ZM.pixscale**2)  # inherited from multimem.pl
        kernel = (psf / knorm)

        ### IMPORTANT KERNEL NORMALIZATION IS NOT knorm=np.sum(psf) BECAUSE THIS IS NOT SMOOTHING. THE CORRECT NORMALIZATION DERIVES FROM THE DIRTY MAP AND DIRTY BEAM FORMULAE.

        if ZM.ViewInit:
            print("kernel")
            Vtools.View(kernel)

        afield['kernel'] = np.double(kernel)

    noise_mosaic = np.sqrt(1. / sensmos)

    noise_mosaic[(noise_mosaic > 1E-1)] = -1  # inherited from multimem.pl

    PBmosaic = 1. / noise_mosaic
    PBmosaic[PBmosaic < 0.] = 0.
    PBmosaic /= np.max(PBmosaic)

    ZM.PBmosaic = PBmosaic
    ZM.noise_mosaic = noise_mosaic
    ZM.noise = np.min(noise_mosaic[(noise_mosaic > 0)])
    ZM.minpix = ZM.minpix_noise * ZM.noise
    print("minpix is ", ZM.minpix)

    fits.writeto(ZM.workdir + "noise_mosaic.fits",
                 noise_mosaic,
                 hrstor,
                 overwrite=True)
    fits.writeto(ZM.workdir + "attenuation.fits",
                 PBmosaic,
                 hrstor,
                 overwrite=True)

    if isinstance(ZM.prior, str):
        hdu0 = fits.open(ZM.prior)
        hdu1 = gridding(hdu0, hcanvas, ReturnHDUList=True)
        im_prior = hdu1[0].data
        if ZM.doclip:
            print("clipping prior to > ", ZM.minpix)
            im_prior[(im_prior < ZM.minpix)] = ZM.minpix
        ZM.data_prior = im_prior / np.e


def restorematrix(xvec):
    im = xvec.reshape(ZM.canvas_shape)
    return im


#def proc_afield(modelimage,afield):
#    bim = modelimage.copy()
#    PB = afield['data_PB']
#    bim *= PB
#    kernel=afield['data_psf']
#    kernel/=np.sum(kernel)
#    dirtymodel=scipy_convolve(bim, kernel, mode='same')
#    residual=afield['data_image']-dirtymodel
#    achi2=np.sum(residual**2/(afield['sigma']**2))
#    return achi2


def neglnlike(xfree, Restoring=False):

    if (ZM.doclip):
        xfree[xfree < ZM.minpix] = ZM.minpix

    modelimage = restorematrix(xfree)
    chi2 = 0.
    for afield in ZM.datapass:
        bim = modelimage.copy()
        PB = afield['data_PB']
        bim *= PB
        kernel = afield['kernel']
        #kernel/=np.sum(kernel)
        dirtymodel = scipy_convolve(bim, kernel,
                                    mode='same')  #,boundary='fill'
        residual = afield['data_image'] - dirtymodel
        chi2 += np.sum(residual**2 / (afield['sigma']**2))

        if Restoring:
            afield['residuals'] = residual
    if (ZM.lambdaS > 0.) & (not Restoring):
        #print("np.min(xfree)",np.min(xfree),"np.min(ZM.prior_vec)",np.min(ZM.prior_vec),"np.min(xfree/ZM.prior_vec)",np.min(xfree/ZM.prior_vec)) # DEV
        entropy = np.sum(xfree * np.log(xfree / ZM.prior_vec))
    else:
        entropy = 0.

    L = chi2 - ZM.lambdaS * entropy
    # print("chi2 %e  entropy %e  L %e = " % ( chi2, entropy, L))
    return L


def dneglnlike(xfree):
    if (ZM.doclip):
        xfree[xfree < ZM.minpix] = ZM.minpix

    # my $I = &$aux($x);
    modelimage = restorematrix(xfree)

    if ZM.DumpAllIters:
        fits.writeto(ZM.workdir + "MEM_" + str(ZM.iterdf) + ".fits",
                     modelimage,
                     ZM.rstor_head,
                     overwrite=True)
    #else:
    #    fits.writeto(ZM.workdir+"MEM_last.fits",modelimage,ZM.rstor_head,overwrite=True)

    View = False
    if ((ZM.iterdf % 100) == 0):
        print("ZM.iterdf", ZM.iterdf)
        View = ZM.View

    chi2 = 0.
    deriv = np.zeros(xfree.shape)
    deriv2Dstack = np.zeros(modelimage.shape)
    for afield in ZM.datapass:
        bim = modelimage.copy()
        if View:
            print("inputmodel")
            Vtools.View(bim)

        PB = afield['data_PB']
        bim *= PB
        kernel = afield['kernel']
        #kernel/=np.sum(kernel)

        if View:
            Vtools.View(PB)
            Vtools.View(bim)

        dirtymodel = scipy_convolve(bim, kernel,
                                    mode='same')  #,boundary='fill'
        if View:
            print("dirtymodel")
            Vtools.View(dirtymodel)

        residual = afield['data_image'] - dirtymodel
        chi2 += np.sum(residual**2 / (afield['sigma']**2))

        A = 2 * (dirtymodel - afield['data_image']) / afield['sigma']**2
        deriv2D = scipy_convolve(A, kernel, mode='same')
        deriv2D *= PB
        if View:
            print("deriv2D")
            Vtools.View(deriv2D)

        deriv2Dstack += deriv2D

        View = False

    deriv = deriv2Dstack.flatten()

    View = False
    if (((ZM.iterdf % 100) == 0) and View):
        View = True
    if View:
        print("modelimage")
        Vtools.View(modelimage)
        print("total deriv image")
        Vtools.View(deriv2Dstack)

    if ZM.lambdaS > 0.:
        xfree[(xfree < ZM.minpix)] = ZM.minpix
        memderiv = np.log(xfree / ZM.prior_vec) + 1.
    else:
        memderiv = 0.

    #retvect = deriv - lambdaS * memderiv
    #retvect = 0.5*deriv
    retvect = deriv - ZM.lambdaS * memderiv
    #print("norm deriv ", np.sum(deriv**2), " norm memderiv ",
    #      np.sum(memderiv**2), " ZM.lambdaS", ZM.lambdaS)
    retvect[(ZM.PBmosaic.flatten() <= 0.)] = 0.
    ZM.iterdf += 1
    #print("ZM.iterdf",ZM.iterdf)
    return retvect


def restore_skymem(ZM):

    neglnlike(ZM.modout, Restoring=True)
    dum = np.zeros(ZM.canvas_shape)
    norm = np.zeros(ZM.canvas_shape)
    for afield in ZM.datapass:
        dum += afield['residuals'] * afield['data_PB'] / (afield['sigma']**2)
        norm += afield['data_PB']**2 / (afield['sigma']**2)
    mos_residuals = dum / norm
    fits.writeto(ZM.workdir + "residual_mosaic.fits",
                 mos_residuals,
                 ZM.rstor_head,
                 overwrite=True)

    modout = ZM.modout
    pixscale = ZM.rstor_head['CDELT2']
    bmaj = ZM.rstor_head['BMAJ']
    bmin = ZM.rstor_head['BMIN']
    bpa = ZM.rstor_head['BPA']
    stdev_x = (bmaj / (2. * np.sqrt(2. * np.log(2.)))) / pixscale
    stdev_y = (bmin / (2. * np.sqrt(2. * np.log(2.)))) / pixscale

    smodout = Gausssmooth.Gauss_filter(modout, stdev_x, stdev_y, bpa)
    restored = smodout + mos_residuals

    fits.writeto(ZM.workdir + "rstor_out.fits",
                 restored,
                 ZM.rstor_head,
                 overwrite=True)

    restored_PBmult = restored * ZM.PBmosaic

    fits.writeto(ZM.workdir + "rstor_out_PBmult.fits",
                 restored_PBmult,
                 ZM.rstor_head,
                 overwrite=True)


def run_scipy_optimize_minimize(ZM):

    print("starting op.minimize")
    start_time = time.time()
    #ftol=1E-8 # 1e-10 too small leads to abnormal termination
    xfree = ZM.xfree
    #eps=0.1*np.ones(len(xfree))

    #options={'gtol': 1e-05, 'norm': inf, 'eps': 1.4901161193847656e-08, 'maxiter': None, 'disp': False, 'return_all': False, 'finite_diff_rel_step': None}
    #method='Newton-CG'
    #method='BFGS'
    #method='L-BFGS-B'
    method = ZM.scipy_optimize_method
    if ZM.VerboseInit:
        op.show_options(solver="minimize", method=method)
    bounds = None
    if (method == 'CG'):
        options = {
            'gtol': 1e-05,
            'norm': np.inf,
            'eps': 1.4901161193847656e-08,
            'maxiter': 1000,
            'disp': True,
            'return_all': False,
            'finite_diff_rel_step': None
        }
    elif (method == 'Newton-CG'):
        options = {
            'disp': True,
            'xtol': 6e-7,
            'maxiter': 100,
            'eps': 1e-8,
            'return_all': False
        }
    elif (method == 'BFGS'):  ## broken in scipy
        options = {
            'disp': True,
            'maxiter': 100,
            'gtol': 1e-05,
            'norm': 1.,
            'eps': 1e-08,
            'return_all': False,
            'finite_diff_rel_step': None
        }
    elif (method == 'L-BFGS-B'):
        options = {
            'disp': True,
            'maxcor': 10,
            'ftol': 2.220446049250313e-09,
            'gtol': 1e-05,
            'eps': 1e-08,
            'maxfun': 15000,
            'maxiter': 1000,
            'iprint': -1,
            'maxls': 20,
            'finite_diff_rel_step': None
        }
        bounds = np.zeros((len(xfree), 2))
        bounds[:, 0] = ZM.minpix
        bounds[:, 1] = ZM.noise * 1E6
    else:
        sys.exit("not a recognised optim method", scipy_optimize_method)

    result = op.minimize(neglnlike,
                         xfree,
                         method=method,
                         jac=dneglnlike,
                         tol=None,
                         callback=None,
                         bounds=bounds,
                         options=options)  # , args=ZM

    result_ml = result["x"]
    print("Optim done in (elapsed time):", time.time() - start_time)
    return result_ml


def my_f(v, params):
    f = neglnlike(v)
    return f


def my_df(v, params):
    df = dneglnlike(v)
    return df


def my_fdf(v, params):
    #print( "\t\t-- fdf -- ")
    f = my_f(v, params)
    df = my_df(v, params)
    #print( "\t\t-- fdf -- ")
    return f, df


def run_pygsl_minimize(ZM):

    from pygsl import _numobj as numx
    from pygsl import multiminimize, errno

    print("starting pygsl.multiminimize")
    start_time = time.time()
    #ftol=1E-8 # 1e-10 too small leads to abnormal termination
    xfree = ZM.xfree
    size = len(xfree)

    params = numx.array((1., 2.), )  # dummy
    sys = multiminimize.gsl_multimin_function_fdf(my_f, my_df, my_fdf, params,
                                                  size)
    solver = multiminimize.conjugate_fr(sys, size)

    #start = numx.array( (5., 7.), )
    start = numx.array(xfree)
    #solver.set(start, 0.01, 1e-4)
    solver.set(start, 1., 1e-1)
    print("Using solver ", solver.name())
    print("%5s %9s %9s  %9s %9s %9s" % ("iter", "x", "y", "f", "dx", "dy"))
    fprev = 0.  # SIMON MODIF
    for iter in range(200):
        status = solver.iterate()
        gradient = solver.gradient()
        x = solver.getx()
        f = solver.getf()

        #SIMON MODIF START
        normgradient = np.sqrt(np.sum(gradient**2))
        if (fprev == 0.):
            if (f == 0.):
                status = errno.GSL_SUCCESS
                break
        else:
            reldiff = np.fabs((f - fprev) / fprev)
            print("normgradient", normgradient, "nelem", len(x), "f", f,
                  "reldiff", reldiff)  #SIMON MODIF
            if (reldiff < 1E-8):
                status = errno.GSL_SUCCESS
                print(
                    "CONVERGENCE BY VALUE - not by gradient, as the GSL default ! \n"
                )
                break


### SIMON MODIF END

        status = multiminimize.test_gradient(gradient, 1e-2)
        if status == errno.GSL_SUCCESS:
            print("Converged ")
        print("%5d % .7f % .7f  % .7f % .7f % .7f" %
              (iter, x[0], x[1], f, gradient[0], gradient[1]))
        if status == errno.GSL_SUCCESS:
            break
        fprev = f
    else:
        raise ValueError("Number of Iterations exceeded!")

    print("Optim done in (elapsed time):", time.time() - start_time)
    result_ml = solver.getx()

    return result_ml


def exec_ConjGrad(ZM):

    pass_setup(ZM)

    ZM.prep_files()
    print("Init ConjGrad:")

    if ('numpy' in str(type(ZM.data_prior))):
        ZM.xfree = ZM.data_prior.flatten()
        ZM.prior_vec = ZM.data_prior.flatten()
    else:
        ZM.xfree = ZM.minpix * np.ones(ZM.canvas_shape).flatten()
        ZM.prior_vec = ZM.minpix * np.ones(ZM.canvas_shape).flatten()

    if ZM.DoGSL:
        (result_ml) = run_pygsl_minimize(ZM)
    else:
        (result_ml) = run_scipy_optimize_minimize(ZM)

    if (ZM.doclip):
        result_ml[result_ml < ZM.minpix] = ZM.minpix

    modout = restorematrix(result_ml)
    ZM.modout = modout
    fits.writeto(ZM.workdir + "mod_out.fits",
                 modout,
                 ZM.rstor_head,
                 overwrite=True)

    restore_skymem(ZM)

    return result_ml


class Setup():

    def __init__(
            self,
            datapass=[],
            prior=False,  # pass string
            datafieldref=False,  # pass a dictionary
            pbcutoff=20.,
            GaussPB=True,  # replace input PB with Gaussian
            PBextend=True,  # extend input PB with Gaussian or replace with Gaussian
            minpix=0.,  # absolute in Jy/beam
            minpix_noise=0.,  # relative to noise
            lambdaS=0.,
            VerboseInit=True,
            workdir='ouput/',
            doclip=True,
            scipy_optimize_method='CG',
            DoGSL=False,  ## broken in pygsl
            View=False,  # view intermediate results
            ViewInit=False,  # view mosaic setup
            DumpAllIters=False):

        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            if VerboseInit:
                print("SkyMEM setting ", a_attribute, " to ",
                      initlocals[a_attribute])
            setattr(self, a_attribute, initlocals[a_attribute])

        if (scipy_optimize_method == 'L-BFGS-B'):
            self.doclip = False

        if (not datafieldref):
            self.datafieldref = datapass[0]

        self.data_prior = False
        self.prior_vec = False
        self.canvas_shape = [0, 0]
        self.rstor_head = False
        self.beam_deg2 = 1.
        self.pixscale = 1.
        self.noise_mosaic = False
        self.PBmosaic = False
        self.xfree = False
        self.iterdf = 0
        self.modout = False

    def prep_files(self):
        return exec_prep_files(self)

    def run(self):
        return exec_ConjGrad(self)
