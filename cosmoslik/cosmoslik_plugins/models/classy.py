from __future__ import absolute_import
from cosmoslik import SlikPlugin
from numpy import arange, pi

class classy(SlikPlugin):
    """
    Plugin for CLASS.
    Credit: Brent Follin, Teresa Hamill, Andy Scacco
    """

    #{cosmoslik name : class name} - This needs to be done even for variables with the same name (because of for loop in self.model.set)!
    name_mapping = {#'As':'A_s',
                    #'ns':'n_s',
                    #'r':'r',
                    'custom1':'custom1',
                    'custom2':'custom2',
                    'custom3':'custom3',
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
                    'omk':'Omega_k',
                    'l_max_scalar':'l_max_scalars',
                    'l_max_tensor':'l_max_tensors',
                    'Tcmb':'T_cmb'
                    }


    def __init__(self):
        super(classy,self).__init__()

        try:
            from classy import Class
        except ImportError:
            raise Exception("Failed to import CLASS python wrapper 'Classy'.")

        self.model = Class()

    #def __call__(self,
    #             **kwargs):
    
    #    d={}
     #   for k, v in kwargs.iteritems():
      #      if k in self.name_mapping and v is not None:
       #         d[self.name_mapping[k]]=v
        #    else:
         #       d[k]=v
    
    #def __call__(self,
                 #ombh2,
                 #omch2,
                 #H0,
                 #As,
                 #ns,
                 #custom1,
                 #custom2,
                 #custom3,
                 #tau,
                 #w=None,
                 #r=None,
                 #nrun=None,
                 #omk=0,
                 #Yp=None,
                 #Tcmb=2.7255,
                 #massless_neutrinos=3.046,
                 #l_max_scalar=3000,
                 #l_max_tensor=3000,
                 #pivot_scalar=0.05,
                 #outputs=[],
                 #**kwargs):

        #print kwargs
        
    def __call__(self,**kwargs):
        #print kwargs
        #print kwargs['classparamlist']
        #print kwargs['d']
        
        d={}
        for k,v in kwargs.iteritems():
            if k in kwargs['classparamlist']:
                if k in self.name_mapping and v is not None:
                    d[self.name_mapping[k]]=v
                else:
                    d[k]=v
            
        
        #d['P_k_ini type']='external_Pk'
        #d['modes'] = 's,t'
        self.model.set(**d)
                       
        l_max = d['l_max_scalars']
        Tcmb =  d['T_cmb']
        
        #print l_max

        #print d
        
        self.model.compute()

        ell = arange(l_max+1)
        self.cmb_result = {'cl_%s'%x:(self.model.lensed_cl(l_max)[x.lower()])*Tcmb**2*1e12*ell*(ell+1)/2/pi
                           for x in ['TT','TE','EE','BB','PP','TP']}

        self.model.struct_cleanup()
        self.model.empty()
        
        return self.cmb_result

    def get_bao_observables(self, z):
        return {'H':self.model.Hubble(z),
                'D_A':self.model.angular_distance(z),
                'c':1.0,
                'r_d':(self.model.get_current_derived_parameters(['rs_rec']))['rs_rec']}
