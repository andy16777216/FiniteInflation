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

        name_mapping = {'As':'A_s',
                    'ns':'n_s',
                    'r':'r',
                    'xi0':'custom1',
                    'k1':'custom2',
                    'xi2':'custom3',
                    'k3':'custom4',
                    'xi4':'custom5',
                    'k5':'custom6',
                    'xi6':'custom7',
                    'k7':'custom8',
                    'xi8':'custom9',
                    'A':'custom10',
                    'nt':'n_t',
                    'ombh2':'omega_b',
                    'omch2':'omega_cdm',
                    'omnuh2':'omega_ncdm',
                    'tau':'tau_reio',
                    'H0':'H0',
                    'massive_neutrinos':'N_ncdm',
                    'massless_neutrinos':'N_ur',
                    'Yp':'YHe',
                    'pivot_scalar':'k_pivot',
                    }
                  
                  
        #d={name_mapping[k]:v for k,v in locals().items() 
        #if k in name_mapping and v is not None}
        d['P_k_ini type']='external_Pk'
        d['modes'] = 's'
        d['output']='tCl, lCl, pCl'
        d['lensing'] ='yes'
        l_max_scalar=3000
        d['l_max_scalars']=l_max_scalar
        d['command'] = 'python ../inflarec/pkinflarec.py'

	ombh2 = param(0.0221)
        omch2 = param(0.12)
        tau = param(0.09, range=(0.05,0.15))
        theta = param(0.010413)
        xi0 = param(-2, scale = 2.4, range = (-3, 3))
        k1 = param(3.2, scale = 2.4, range = (0, 6))
        xi2 = param(0.4, scale = 2.4, range = (-3, 3))
        k3 = param(2.8, scale = 2.4, range = (0, 6))
        xi4 = param(2, scale = 2.4, range = (-3, 3))
        k5 = param(2.4, scale = 2.4, range = (0, 6))
        xi6 = param(0.99, scale = 2.4, range = (-3, 3))
        k7 = param(2, scale = 2.4, range = (0, 6))
        xi8 = param(0.98, scale = 2.4, range = (-3, 3))
        A = param(2.5, scale = 2.4, range = (0, 6))
        logA = None
        As = None
        ns = None
        k_c = None
        alpha_exp = None
        massless_neutrinos=3.046 #param(3,.2)
        l_max_tensor=3000
        pivot_scalar=0.05
        w=-1.0
        r=None
        nrun=None
        omk=0
        Yp=None
        Tcmb=2.7255

        self.cosmo = get_plugin('models.cosmology')(
            {name_mapping[k]:v for k,v in locals().items() 
        if k in name_mapping and v is not None},
            **d
        )
        


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
             output_extra_params=['cosmo.Yp','cosmo.H0']
	)
        
    def __call__(self):
        self.cosmo.Yp = self.bbn(**self.cosmo)
        self.cosmo.H0 = self.hubble_theta.theta_to_hubble(**self.cosmo)
	    #print 'getting cmb'
        self.cmb_result = self.get_cmb(force = True, outputs=['cl_TT','cl_TE','cl_EE','cl_BB','cl_PP','cl_TP'],**self.cosmo)
        
        return lsum(lambda: self.priors(self),
                    lambda: self.camspec(self.cmb_result),
                    lambda: self.lowl(self.cmb_result),
                    lambda: self.pol(self.cmb_result)
                    )

if __name__=='__main__':
     #run the chain
     for _ in Slik(main()).sample(): pass
