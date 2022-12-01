import sys
import numpy as np
from copy import deepcopy

include_path = '/home/simon/gitcommon/SkyMEM/python/'
sys.path.append(include_path)

import SkyMEM


sourcedir = './skymem/'
basemsname = 'TYC9340.12m.cal.bary.continuum.fav.tav.SMGsub.corrected.ms'
datapass = [{
    'image': sourcedir + basemsname + '.dirty.fits',
    'psf': sourcedir + basemsname + '.psf.fits',
    'sens': sourcedir + basemsname + '.sens.fits',
    'headfile': sourcedir + basemsname + '.dirty.fits'
}]

prior = sourcedir+'prior.fits'
### prior = '/home/simon/gitcommon/SkyMEM/example_data/mod_in_prior.fits'

fieldref = datapass[0]

GaussPB = False
pbcutoff = -1
workdir = 'mem_l1_CG_noprior'
PBextend = False
if GaussPB:
    PBextend = True
    pbcutoff = 20.
    workdir += '_GaussPB'
workdir += '/'

M = SkyMEM.Setup(
    datapass=datapass,
    prior=False, # prior, #False, # prior,
    datafieldref=fieldref,
    pbcutoff=pbcutoff,
    GaussPB=GaussPB,
    PBextend=PBextend,
    workdir=workdir,
    #scipy_optimize_method='Newton-CG',
    #scipy_optimize_method='L-BFGS-B',
    minpix_noise=1E-3,
    scipy_optimize_method='CG',
    DoGSL=False,  ## broken in pygsl
    DumpAllIters=True,
    View=False,
    ViewInit=False,
    VerboseInit=True,
    lambdaS=4000.) #5000

M.run()
