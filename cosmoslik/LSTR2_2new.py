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
        d['modes'] = 's,t'
        d['output']='tCl, lCl, pCl'
        d['command'] = '../LSTfiniteR2/pk'
        d['ombh2'] = param(0.0221)
        d['omch2'] = param(0.12)
        d['tau'] = param(0.09, range=(0.05,0.15))
        #d['theta'] = param(0.010413)
        d['custom1'] = param(9.8, scale = 2.4, range = (9.2, 10.0)) #phi0
        d['custom2'] = param(3.5, scale = 2.4, range = (2.5, 4.5)) #L
        d['custom3'] = param(8.4, scale = 2.4, range = (6.6, 10.4)) #logkc
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
        
        self.parameters = d

        self.cosmo = get_plugin('models.cosmology')(
            #ombh2 = param(0.0221),
            #omch2 = param(0.12),
            #tau = param(0.09, range=(0.05,0.15)),
            theta = param(0.010413),
            #phi0 = param(9.8, scale = 2.4, range = (9.2, 10.0)),
            #L = param(3.5, scale = 2.4, range = (2.5, 4.5)),
            #logkc = param(8.4, scale = 2.4, range = (6.6, 10.4)),
            logA = None,
            As = None,
            ns = None,
            k_c = None,
            alpha_exp = None,
            #massless_neutrinos=3.046, #param(3,.2)
            #l_max_scalar=3000,  #These variables are not set here, but in classy.py, must be edited there!!
            #l_max_tensor=3000,
            #pivot_scalar=0.05,
            w=-1.0,
            r=None,
            nrun=None,
            #omk=0,
            Yp=None,
            #Tcmb=2.7255,
            #lensing = 'yes',
            **d
        )
        
	print self.cosmo

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
             num_samples=4000,
             output_file='chains/LSTR2_2.chain',
             proposal_cov='r2cov.covmat',
             proposal_scale=2.4,
             #print_level=0,
             output_extra_params=['cosmo.Yp','cosmo.H0']
	)
        
    def __call__(self):
        self.cosmo.Yp = self.bbn(**self.cosmo)
        self.cosmo.H0 = self.hubble_theta.theta_to_hubble(**self.cosmo)
	    #print 'getting cmb'
        self.cmb_result = self.get_cmb(self.cosmo.Yp,self.cosmo.H0,self.parameters)
        
        return lsum(lambda: self.priors(self),
                    lambda: self.camspec(self.cmb_result),
                    lambda: self.lowl(self.cmb_result),
                    lambda: self.pol(self.cmb_result)
                    )

if __name__=='__main__':
     #run the chain
     for _ in Slik(main()).sample(): pass
