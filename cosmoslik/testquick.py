from cosmoslik import *

import os.path as osp
        
param = param_shortcut('start','scale')
   
class planck(SlikPlugin):

    def __init__(self, camspec_clik_file, model='lcdm'):
        super(SlikPlugin,self).__init__(**all_kw(locals()))

        self.cosmo = get_plugin('models.cosmology')(
            logA = param(3.2),
            ns = param(0.96),
            ombh2 = param(0.0221),
            omch2 = param(0.12),
            tau = param(0.09,min=0,gaussian_prior=(0.085,0.015)),
            theta = param(0.010413),
            omnuh2 = 0.000645,
            massive_neutrinos=1,
            massless_neutrinos=2.046,
        )
        if 'neff' in model: self.cosmo.massless_neutrinos = param(3,.2)
        if 'yp' in model: self.cosmo.Yp = param(.24,0.1)
        if 'mnu' in model: self.cosmo.omnuh2 = param(0,0.001,range=(0,1))

        self.camspec = get_plugin('likelihoods.clik_like')(
            clik_file=camspec_clik_file,
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
            xi_sz_cib=param(0.5,range=(-1,1),scale=0.2),
            A_ksz=param(1,range=(0,5)),
            Bm_1_1=param(0,gaussian_prior=(0,1),scale=1)
        )

        self.get_cmb = get_plugin('models.camb')()
        self.bbn = get_plugin('models.bbn_consistency')()
        self.hubble_theta = get_plugin('models.hubble_theta')()
        self.priors = get_plugin('likelihoods.priors')(self)

        self.sampler = get_plugin('samplers.metropolis_hastings')(
             self,
             num_samples=1000000,
             output_file='chains/chain_%s_%s.chain'%(model,osp.basename(camspec_clik_file)),
             proposal_cov='/software/mint15/cosmomc/11.13/cosmomc/covmats/planck.covmat',
             proposal_scale=1,
             print_level=1,
             output_extra_params=['cosmo.Yp','cosmo.H0']
        )

    def __call__(self):
        self.cosmo.As = exp(self.cosmo.logA)*1e-10
        if 'yp' not in self.model: self.cosmo.Yp = self.bbn(**self.cosmo)
        self.cosmo.H0 = self.hubble_theta.theta_to_hubble(**self.cosmo)

        self.cmb_result = self.get_cmb(outputs=['cl_TT'],**self.cosmo)

        return lsum(lambda: self.priors(self),
                    lambda: self.camspec(self.cmb_result))
