import sys
import numpy as np
from copy import deepcopy

include_path='/home/simon/common/python/include/'
sys.path.append(include_path)

import SkyMEM


freq='20160'
#freq='39157'

sourcedirPSF='/home/simon/common/ROPH2/ATCA/mem/MemSky/data_128/'
sourcedir='/home/simon/common/ROPH2/ATCA/mem/MemSky/data_256/'
labels = ('1_1',  '1_2', '1_3',  '2_1', '2_2', '2_3')
datapass = map(lambda alabel: {'image':sourcedir+'oph_'+alabel+'.'+freq+'.fits',
                              'psf':sourcedirPSF+'oph_'+alabel+'.'+freq+'.beam.fits',
                              'sens':sourcedir+'oph_'+alabel+'.'+freq+'.sensitivity.fits',
                              'headfile':sourcedir+'oph_'+alabel+'.'+freq+'.rstor.fits'}, labels)


datapass=list(datapass)

prior='/home/simon/common/ROPH2/ATCA/mem/MemSky/out_mem7_'+freq+'_l0_fix3/mod_in_prior.fits'

reflabel = '1_2'
fieldref = {'image':sourcedirPSF+"/oph_"+reflabel+"."+freq+".fits",
            'psf':sourcedirPSF+"/oph_"+reflabel+"."+freq+".beam.fits",
            'sens':sourcedirPSF+"/oph_"+reflabel+"."+freq+".sensitivity.fits",
            'headfile':sourcedirPSF+"/oph_"+reflabel+"."+freq+".rstor.fits"}


GaussPB=False
pbcutoff=-1
workdir='mem_l0_'+freq+'_ATCAtrial'
#workdir='mem_l0_'+freq
PBextend=False
if GaussPB:
    PBextend=True
    pbcutoff=20.
    workdir+='_GaussPB'
workdir+='/'
    

M=SkyMEM.Setup(
    datapass=datapass,
    prior=prior,
    datafieldref=fieldref,
    pbcutoff=pbcutoff,
    GaussPB=GaussPB,
    PBextend=PBextend,
    workdir=workdir, 
    #scipy_optimize_method='Newton-CG',
    #scipy_optimize_method='L-BFGS-B',
    scipy_optimize_method='CG',
    View=False,
    ViewInit=True,
    VerboseInit=True,
    lambdaS=0)


M.run()
