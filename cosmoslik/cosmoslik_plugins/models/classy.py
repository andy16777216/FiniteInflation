from __future__ import absolute_import
from cosmoslik import SlikPlugin
from numpy import arange, pi

class classy(SlikPlugin):
    """
    Plugin for CLASS.
    Credit: Brent Follin, Teresa Hamill, Andy Scacco
    """
    
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
                 As,
                 ns,
                 tau,
                 w=None,
                 r=None,
                 nrun=None,
                 omk=0,
                 Yp=None,
                 Tcmb=2.7255,
                 massless_neutrinos=3.046,
                 l_max_scalar=3000,
                 l_max_tensor=3000,
                 pivot_scalar=0.05,
                 outputs=[],
                 **kwargs):


        #self.model.set(**d)
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
