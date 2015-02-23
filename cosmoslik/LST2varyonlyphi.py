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

        self.cosmo = get_plugin('models.cosmology')(
            phi0 = param(20.5, scale = 0.1, range = (19.7,21.3)),
            ombh2 = 0.022153290145926357,
            omch2 = 0.11804578115504939,
            tau = 0.088064215067749391,
            theta = 0.010417970388328802,
            m6 = 5.962621883596591,
            logA = None,
            As = None,
            ns = None,
            k_c = None,
            alpha_exp = None,
            massless_neutrinos=3.046, #param(3,.2)
            l_max_scalar=3000,  #These variables are not set here, but in classy.py, must be edited there!!
            l_max_tensor=3000,
            pivot_scalar=0.05,
            w=-1.0,
            r=None,
            nrun=None,
            omk=0,
            Yp=None,
            Tcmb=2.7255,
            lensing = 'yes'#,
        )

	#print 'loading likelihoods'
        self.camspec = get_plugin('likelihoods.clik')(
            clik_file='/software/mint15/cosmomc/likelihoods/clik_0313/data/CAMspec_v6.2TN_2013_02_26_dist.clik',
            A_ps_100=149.86133092590069,
            A_ps_143=42.921337942177814,
            A_ps_217=98.131910433229166,
            A_cib_143=8.8925518782832018,
            A_cib_217=35.03414675434886,
            A_sz=98.131910433229166,
            r_ps=0.86985724867605718,
            r_cib=0.45674724749543938,
            n_Dl_cib=0.61502611065131785,
            cal_100=1.0005917033465019,
            cal_217=0.99631759854151036,
            xi_sz_cib=-0.065567371172318839,
            A_ksz=2.2992427075566164,
            Bm_1_1=0.39688376738950149
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
             num_samples=1000,
             output_file='chains/LST2vop.chain',
             proposal_cov='../data/proposal.covmat',
             proposal_scale=1,
             #print_level=0,
             output_extra_params=['cosmo.Yp','cosmo.H0','cl_TT2','cl_TT3','cl_TT4','cl_TT5','cl_TT6','cl_TT7','cl_TT8','cl_TT20','cl_TT40','cl_TT80','cl_TT120','cl_TT200','cl_TT500']
	)


        
    def __call__(self):
        self.cosmo.Yp = self.bbn(**self.cosmo)
        self.cosmo.H0 = self.hubble_theta.theta_to_hubble(**self.cosmo)
	    #print 'getting cmb'
        self.cmb_result = self.get_cmb(force = True, outputs=['cl_TT','cl_TE','cl_EE','cl_BB','cl_PP','cl_TP'],**self.cosmo)
        self.cl_TT = self.cmb_result['cl_TT']
        self.cl_TT2 = self.cl_TT[2]
        self.cl_TT3 = self.cl_TT[3]
        self.cl_TT4 = self.cl_TT[4]
        self.cl_TT5 = self.cl_TT[5]
        self.cl_TT6 = self.cl_TT[6]
        self.cl_TT7 = self.cl_TT[7]
        self.cl_TT8 = self.cl_TT[8]
        self.cl_TT20 = self.cl_TT[20]
        self.cl_TT40 = self.cl_TT[40]
        self.cl_TT80 = self.cl_TT[80]
        self.cl_TT120 = self.cl_TT[120]
        self.cl_TT200 = self.cl_TT[200]
        self.cl_TT500 = self.cl_TT[500]
        
        return lsum(lambda: self.priors(self),
                    lambda: self.camspec(self.cmb_result),
                    lambda: self.lowl(self.cmb_result),
                    lambda: self.pol(self.cmb_result)
                    )


if __name__=='__main__':
     #run the chain
     for _ in Slik(main()).sample(): pass
