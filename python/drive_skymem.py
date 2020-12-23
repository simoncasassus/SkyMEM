import sys
import numpy as np
from copy import deepcopy

include_path='/home/simon/common/python/include/'
sys.path.append(include_path)

import SkyMEM


freq='20160'
sourcedir='/home/simon/common/ROPH2/ATCA/mem/MemSky/data_128/'
labels = ('1_1',  '1_2', '1_3',  '2_1', '2_2', '2_3')
datapass = map(lambda alabel: {'image':sourcedir+'oph_'+alabel+'.'+freq+'.fits',
                              'psf':sourcedir+'oph_'+alabel+'.'+freq+'.beam.fits',
                              'sens':sourcedir+'oph_'+alabel+'.'+freq+'.sensitivity.fits',
                              'headfile':sourcedir+'oph_'+alabel+'.'+freq+'.rstor.fits'}, labels)


datapass=list(datapass)

prior='/home/simon/common/ROPH2/ATCA/mem/MemSky/out_mem7_20160_l0_fix3/mod_in_prior.fits'

#my $reflabel = '1_2';
#fieldref=datapass[1]


M=SkyMEM.Setup(
    datapass=datapass,
    #prior=prior,
    ifieldref=1,
    #pbcutoff=20.,
    pbcutoff=20.,
    workdir='mem_l0/',
    #scipy_optimize_method='Newton-CG',
    #scipy_optimize_method='L-BFGS-B',
    scipy_optimize_method='CG',
    #scipy_optimize_method='BFGS', <<broken in scipy
    #scipy_optimize_method='CG', 
    DoGSL=False,
    View=False,
    lambdaS=0)


M.run()
