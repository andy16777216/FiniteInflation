from cosmoslik import *
from numpy import interp, identity, exp, inf, arange, hstack, loadtxt, zeros, ones, log, invert, isinf
import sys
import math
import csv
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
        d['command'] = '../LSTR2_2/pk'
        d['ombh2'] = param(0.02207863858642484, scale = 0.00025176398454803333)
        d['omch2'] = param(0.11957079477757622, scale = 0.00207309880071336)
        d['tau'] = param(0.087773083101280264, scale = 0.012187648687260483, range=(0.05,0.15))
        #d['theta'] = param(0.010413)
        d['custom1'] = param(9.6442079425117182, scale = 0.13033213062684648, range = (9.2, 10.0)) #phi0
        d['custom2'] = param(3.4954418402782967, scale = 0.20282354381662321, range = (2.5, 4.5)) #L
        d['custom3'] = param(8.876773202948355, scale = 0.68803250879040911, range = (7.0, 10.4)) #logkc
        d['massless_neutrinos']=3.046
        d['l_max_scalar']=3000
        d['l_max_tensor']=3000
        d['pivot_scalar']=0.05
        d['lensing'] = 'yes'
        d['Tcmb']=2.7255
        d['omk']=0
        #d['k_per_decade_primordial'] = 100.
        #d['w']=-1.0 #This is used fot the theta to hubble conversion but not in CLASS
        #d['Omega_fld'] = 0 #This and the line below specify a cosmological constant (in CLASS), consistent with w=-1 above
	#d['Omega_scf'] = 0
        
        d['classparamlist']=d.keys()
        d['classparamlist'].extend(['cosmo.Yp','cosmo.H0'])
        
        self.cosmo = get_plugin('models.cosmology')(
            theta = param(0.01041647232923798, scale = 5.7229094877058284e-06),
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
            clik_file='/software/cosmomc/likelihoods/clik_0313/data/CAMspec_v6.2TN_2013_02_26_dist.clik',
            A_ps_100=param(151.84993105553934, scale = 64.753365940060135,min=0),
            A_ps_143=param(41.82997089071911, scale = 16.954706637330407,min=0),
            A_ps_217=param(95.768959993076891, scale = 19.418585842152446,min=0),
            A_cib_143=param(9.420343804906306, scale = 8.1999445289557755,min=0),
            A_cib_217=param(37.085359206553534, scale = 9.2834146629028602,min=0),
            A_sz=param(7.855364936949421, scale = 4.8787513109523815,range=(0,20)),
            r_ps=param(0.83710552570838892, scale = 0.14025187983106405,range=(0,1)),
            r_cib=param(0.49469482845803275, scale = 0.2459406113725135,range=(0,1)),
            n_Dl_cib=param(0.6260361205728302, scale = 0.10416066352042916,gaussian_prior=(0.8,0.2)),
            cal_100=param(1.0005741662931826, scale = 0.00038589325907093972),
            cal_217=param(0.99641097150505309, scale = 0.0013591071204870334),
            xi_sz_cib=param(-0.18021245310622472, scale = 0.56241393555434716,range=(-1,1)),
            A_ksz=param(2.4371002035474274, scale = 1.4508314235334643,range=(0,5)),
            Bm_1_1=param(0.4310444057902954, scale = 0.51288786510446649,gaussian_prior=(0,1))
        )

        self.lowl = get_plugin('likelihoods.clik')(
          clik_file='/software/cosmomc/likelihoods/clik_0313/data/commander_v4.1_lm49.clik'
        )
        
        self.pol = get_plugin('likelihoods.clik')(
          clik_file='/software/cosmomc/likelihoods/clik_0313/data/lowlike_v222.clik'
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
             num_samples=1000000,
             output_file='chains/LSTR2_2new2_5.chain',
             #proposal_cov='../data/proposal.covmat',
             proposal_cov='r2cov.covmat',
             proposal_scale=1,
             print_level=2,
             proposal_update_start=1000,
             mpi_comm_freq=10,
             debug_output=True,
             output_extra_params=['cosmo.Yp','cosmo.H0']
	)
        
    def __call__(self):
        self.cosmo.Yp = self.bbn(**self.cosmo)
        self.cosmo.H0 = self.hubble_theta.theta_to_hubble(**self.cosmo)
	    #print 'getting cmb'
	#print (self.cosmo.custom1, self.cosmo.custom2, self.cosmo.custom3)
	
        #print self.cosmo.classparamlist	
	
	#print self.priors(self)
	
	if isinf(self.priors(self)):
		return inf
	else:
        	self.cmb_result = self.get_cmb(**self.cosmo)
        
        	self.cosmo.loglike = lsum(lambda: self.priors(self),
                    lambda: self.camspec(self.cmb_result),
                    lambda: self.lowl(self.cmb_result),
                    lambda: self.pol(self.cmb_result)
                    )
        
	#with open('LSTR2_2new.csv', 'ab') as csvfile:
    		#spamwriter = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    		#spamwriter.writerow([self.parameters.custom1, self.parameters.custom2, self.parameters.custom3,self.parameters.loglike])
    		#for item in self.parameters:
    		#	spamwriter.writerow(item)
    	
        return self.cosmo.loglike

if __name__=='__main__':
     #run the chain
     for _ in Slik(main()).sample(): pass
