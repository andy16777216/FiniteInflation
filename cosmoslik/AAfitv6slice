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
            logA = 3.0828303593955471,
            ns = 0.95995802437054512,
            k_c = param(-9, scale = 1, range = (-13,-7)),
            alpha_exp = param(0.8, scale = 0.1, range=(0.4,0.9)),
            ombh2 = 0.022062077811123734,
            omch2 = 0.11974442183029269,
            tau = 0.086558937968886043,
            theta = 0.010415977212164643,
            #omnuh2 = 0,	#0.000645,
            #massive_neutrinos=0,#param( 3, .2),
            massless_neutrinos=3.046, #param(3,.2)
            l_max_scalar=3000,
            l_max_tensor=3000,
            pivot_scalar=0.05,
            w=-1.0,
            r=None,
            nrun=None,
            omk=0,
            Yp=None,
            Tcmb=2.7255,
            #P_k_ini type = analytic_Pk,
            lensing = 'yes'#,
            #P_k_max_hinvMpc = 1.
        )

	#print 'setting phase template'
        #self.phase_template = SlikDict()
        #self.phase_template.alpha = 5.5
        #self.phase_template.A = 8.43
        #self.phase_template.B = 6.5e-4
        #if 'neff' in model: self.cosmo.massless_neutrinos = param(3,.2)
        #if 'yp' in model: self.cosmo.Yp = param(.24,0.1)
        #if 'mnu' in model: self.cosmo.omnuh2 = param(0,0.001,range=(0,1))

	#print 'loading likelihoods'
        self.camspec = get_plugin('likelihoods.clik')(
            clik_file='/software/mint15/cosmomc/likelihoods/clik_0313/data/CAMspec_v6.2TN_2013_02_26_dist.clik',
            A_ps_100=155.74974564586893,
            A_ps_143=42.990837159987066,
            A_ps_217=97.695142281342967,
            A_cib_143=9.0610630022180896,
            A_cib_217=35.756240161532823,
            A_sz=7.6487534282036265,
            r_ps=0.84792479189221115,
            r_cib=0.48718759077123552,
            n_Dl_cib=0.61145370143786038,
            cal_100=1.0005742128940385,
            cal_217=0.99635794254317223,
            xi_sz_cib=-0.11312275434269781,
            A_ksz=2.4495876748338397,
            Bm_1_1=0.48688296219336374
        )

        self.lowl = get_plugin('likelihoods.clik')(
          clik_file='/software/mint15/cosmomc/likelihoods/clik_0313/data/commander_v4.1_lm49.clik'
        )

        # self.s12 = get_plugin('likelihoods.spt_lowl')(
        #          which='s12',
        #          cal = param(1,0.02),
        #          fgs = get_plugin('models.clust_poisson_egfs')
        #                  (
        #                              Aps = param(10,10,min=0),
        #                              Acib = param(10,10,min=0),
        #                              ncib = 0.8
        #                  )
        #  )
    	#print 'loading cosmology'

        self.get_cmb = get_plugin('models.classy')()

	#print 'loading derivers'
        self.bbn = get_plugin('models.bbn_consistency')()
        self.hubble_theta = get_plugin('models.hubble_theta')()
        self.priors = get_plugin('likelihoods.priors')(self)

	#print 'loading sampler'
        self.sampler = get_plugin('samplers.metropolis_hastings')(
             self,
             num_samples=1000000,
             output_file='chains/andyAAfitv6.chain',
             proposal_cov='../data/proposal.covmat',
             proposal_scale=1,
             #print_level=0,
             output_extra_params=['cosmo.Yp','cosmo.H0','cosmo.kcactual','cosmo.alphaactual','cl_TT2','cl_TT3','cl_TT4','cl_TT5','cl_TT6','cl_TT7','cl_TT8','cl_TT20','cl_TT40','cl_TT80','cl_TT120','cl_TT200','cl_TT500','cl_TT1000','cl_TT2000','cl_TT2999']
	)


        
    def __call__(self):
        self.cosmo.As = exp(self.cosmo.logA)*1e-10
        self.cosmo.Yp = self.bbn(**self.cosmo)
        self.cosmo.H0 = self.hubble_theta.theta_to_hubble(**self.cosmo)
        self.cosmo.kcactual = exp(self.cosmo.k_c)
        self.cosmo.alphaactual = 1./(1.-self.cosmo.alpha_exp)-1.
        #self.cosmo.neff_phase = self.amp_to_neff()
        #self.cosmo.leq = 125
	    #print 'getting cmb'
        self.cmb_result = self.get_cmb(force = True, outputs=['cl_TT'],**self.cosmo)
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
        self.cl_TT1000 = self.cl_TT[1000]
        self.cl_TT2000 = self.cl_TT[2000]
        self.cl_TT2999 = self.cl_TT[2999]
        
        return lsum(lambda: self.priors(self),
                    lambda: self.camspec(self.cmb_result),
                    lambda: self.lowl(self.cmb_result))


if __name__=='__main__':
     #run the chain
     for _ in Slik(main()).sample(): pass
