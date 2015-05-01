from cosmoslik import *
from numpy import interp, identity, exp, inf, arange, hstack, loadtxt, zeros, ones, log, invert
import sys
import math        
param = param_shortcut('start','scale')

import os.path as osp
        
param = param_shortcut('start','scale')
   
class main(SlikPlugin):

    def __init__(self):
        super(SlikPlugin,self).__init__()

        d = {}
        d['P_k_ini type']='external_Pk'
        d['modes'] = 's'
        d['output']='tCl, lCl, pCl'
        d['command'] = 'python ../inflarec/pkinflarec.py'
        d['ombh2'] = param(0.0221, scale = 0.01)
        d['omch2'] = param(0.12, scale = 0.05)
        d['tau'] = param(0.09, scale = 0.03, range=(0.05,0.15))
        #d['theta'] = param(0.010413)
        d['custom1'] = param(-2, scale = 2.4, range = (-3, 3)) #xi0
        d['custom2'] = param(3.2, scale = 2.4, range = (0, 6)) #k1
        d['custom3'] = param(0.4, scale = 2.4, range = (-3, 3)) #xi2
        d['custom4'] = param(2.8, scale = 2.4, range = (0, 6)) #k3
        d['custom5'] = param(2, scale = 2.4, range = (-3, 3)) #xi4
        d['custom6'] = param(2.4, scale = 2.4, range = (0, 6)) #k5
        d['custom7'] = param(0.99, scale = 2.4, range = (-3, 3)) #xi6
        d['custom8'] = param(2, scale = 2.4, range = (0, 6)) #k7
        d['custom9'] = param(0.98, scale = 2.4, range = (-3, 3)) #xi8
        d['custom10'] = param(2.5, scale = 2.4, range = (0, 6)) #A
        d['massless_neutrinos']=3.046
        d['l_max_scalar']=3000
        d['l_max_tensor']=3000
        d['pivot_scalar']=0.05
        d['lensing'] = 'yes'
        d['Tcmb']=2.7255
        d['omk']=0
        #d['w']=-1.0 #This is used fot the theta to hubble conversion but not in CLASS
        d['Omega_fld'] = 0 #This and the line below specify a cosmological constant (in CLASS), consistent with w=-1 above
	d['Omega_scf'] = 0
        
        self.parameters = get_plugin('models.CLASSparams')(**d)

        self.cosmo = get_plugin('models.cosmology')(
            theta = param(0.010413),
            logA = None,
            As = None,
            ns = None,
            k_c = None,
            alpha_exp = None,
            w=-1.0,
            r=None,
            nrun=None,
            Yp=None,
            **d
        )
        
	#print self.cosmo
	#print self.parameters

	#print 'loading likelihoods'
        self.camspec = get_plugin('likelihoods.clik')(
            clik_file='/software/mint15/cosmomc/likelihoods/clik_0313/data/CAMspec_v6.2TN_2013_02_26_dist.clik',
            A_ps_100=param(150,min=0),
            A_ps_143=param(60,min=0),
            A_ps_217=param(60,min=0),
            A_cib_143=param(10,min=0),
            A_cib_217=param(40,min=0),
            A_sz=param(5,scale=1,range=(0,20)),
            r_ps=param(0.7,range=(0,1)),
            r_cib=param(0.7,range=(0,1)),
            n_Dl_cib=param(0.8,scale=0.2,gaussian_prior=(0.8,0.2)),
            cal_100=param(1,scale=0.001),
            cal_217=param(1,scale=0.001),
            xi_sz_cib=param(0.5,scale=0.2,range=(-1,1)),
            A_ksz=param(1,range=(0,5)),
            Bm_1_1=param(0,scale=1,gaussian_prior=(0,1))
        )

        self.lowl = get_plugin('likelihoods.clik')(
          clik_file='/software/mint15/cosmomc/likelihoods/clik_0313/data/commander_v4.1_lm49.clik'
        )
        
        self.pol = get_plugin('likelihoods.clik')(
          clik_file='/software/mint15/cosmomc/likelihoods/clik_0313/data/lowlike_v222.clik'
        )

    	#print 'loading cosmology'

        self.get_cmb = get_plugin('models.classy')()

	#print 'loading derivers'
        self.bbn = get_plugin('models.bbn_consistency')()
        self.hubble_theta = get_plugin('models.hubble_theta')()
        self.priors = get_plugin('likelihoods.priors')(self)

	#print 'loading sampler'
        self.sampler = get_plugin('samplers.metropolis_hastings')(
             self,
             num_samples=400,
             output_file='chains/infrecMC.chain',
             proposal_cov='r2cov.covmat',
             proposal_scale=2.4,
             #print_level=0,
             output_extra_params=['parameters.Yp','parameters.H0']
	)
        
    def __call__(self):
        self.parameters.Yp = self.bbn(**self.cosmo)
        self.parameters.H0 = self.hubble_theta.theta_to_hubble(**self.cosmo)
	    #print 'getting cmb'
        self.cmb_result = self.get_cmb(**self.parameters)
        
        return lsum(lambda: self.priors(self),
                    lambda: self.camspec(self.cmb_result),
                    lambda: self.lowl(self.cmb_result),
                    lambda: self.pol(self.cmb_result)
                    )

if __name__=='__main__':
     #run the chain
     for _ in Slik(main()).sample(): pass
