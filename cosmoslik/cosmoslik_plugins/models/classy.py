from __future__ import absolute_import
from cosmoslik import SlikPlugin
from numpy import arange, pi

from cosmoslik import *

class classy(SlikPlugin):
    """
    Plugin for CLASS.
    Credit: Brent Follin, Teresa Hamill, Andy Scacco
    """

    #{cosmoslik name : class name} - This needs to be done even for variables with the same name (because of for loop in self.model.set)!
    name_mapping = {#'As':'A_s',
                    #'ns':'n_s',
                    #'r':'r',
                    #'nt':'n_t',
                    'ombh2':'omega_b',
                    'omch2':'omega_cdm',
                    'omnuh2':'omega_ncdm',
                    'tau':'tau_reio',
                    'H0':'H0',
                    'massive_neutrinos':'N_ncdm',
                    'massless_neutrinos':'N_ur',
                    'Yp':'YHe',
                    'pivot_scalar':'k_pivot',
                    #'Tcmb':'T_cmb',
                    #'P_k_max_hinvMpc':'P_k_max_h/Mpc'
                    'w':'w0_fld',
                    'nrun':'alpha_s',
                    'omk':'Omega_k',
                    #'l_max_scalar':'l_max_scalars',
                    #'l_max_tensor':'l_max_tensors',
                    'phi0':'custom1',
                    'm6':'custom2'#,
                    #'modes':'modes',
                    #'lensing':'lensing'
                    }


    def __init__(self):
        super(classy,self).__init__()

        try:
            from classy import Class
        except ImportError:
            raise Exception("Failed to import CLASS python wrapper 'Classy'.")

        self.model = Class()


    def __call__(self,
                 ombh2,
                 omch2,
                 H0,
                 tau,
                 phi0,
                 m6,
                 As = None,
                 ns = None,
                 k_c = None,
                 alpha_exp = None,
                 #omnuh2=0, #0.006  #None means that Class will take the default for this, maybe?
                 w=None,
                 r=None,
                 nrun=None,
                 omk=0,
                 Yp=None,
                 Tcmb=2.7255,
                 #massive_neutrinos=0,
                 massless_neutrinos=3.046,
                 l_max_scalar=3000,
                 l_max_tensor=500,
                 pivot_scalar=0.05,
                 outputs=[],
                 **kwargs):


        d={self.name_mapping[k]:v for k,v in locals().items() 
                          if k in self.name_mapping and v is not None}
        d['P_k_ini type']='external_Pk'
        d['tensor method'] = 'massless'
        self.model.set(output='tCl, lCl, pCl',
                       lensing='yes',
                       #modes = 's,t',
                       l_max_scalars=l_max_scalar,
                       l_max_tensors=l_max_tensor,
                       command = '../LSODA/pk',
                       #custom1 = phi0
                       #custom2 = m6,
                       **d)
        self.model.compute()

        ell = arange(l_max_scalar+1)
        self.cmb_result = {'cl_%s'%x:(self.model.lensed_cl(l_max_scalar)[x.lower()])*Tcmb**2*1e12*ell*(ell+1)/2/pi
                           for x in ['TT','TE','EE','BB','PP','TP']}

        self.model.struct_cleanup()
        self.model.empty()
        
        return self.cmb_result

    def get_bao_observables(self, z):
        return {'H':self.model.Hubble(z),
                'D_A':self.model.angular_distance(z),
                'c':1.0,
                'r_d':(self.model.get_current_derived_parameters(['rs_rec']))['rs_rec']}
