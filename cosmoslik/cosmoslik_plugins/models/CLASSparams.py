from cosmoslik import SlikPlugin, all_kw, param_shortcut, SlikFunction
param = param_shortcut('start','scale')

class CLASSparams(SlikPlugin):    
    
    def __init__(self,**kwargs):
        super(cosmology,self).__init__(**dict(all_kw(locals()),**kwargs))
