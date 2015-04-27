#!/usr/bin/python
import sys
import math
import cmath
#from scipy import special
import scipy.special

# README:
#
# This is an example python script for the external_Pk mode of Class.
# It generates the primordial spectrum of LambdaCDM.
# It can be edited and used directly, though keeping a copy of it is recommended.
#
# Two (maybe three) things need to be edited:
#
# 1. The name of the parameters needed for the calculation of Pk.
#    "sys.argv[1]" corresponds to "custom1" in Class, an so on

try :
    xi0           = float(sys.argv[1])
    k1             = float(sys.argv[2])
    xi2           = float(sys.argv[3])
    k3             = float(sys.argv[4])
    xi4           = float(sys.argv[5])
    k5             = float(sys.argv[6])
    xi6           = float(sys.argv[7])
    k7             = float(sys.argv[8])
    xi8           = float(sys.argv[9])

    parlist = [xi0,k1,xi2,k3,xi4,k5,xi6,k7,xi8]
    
    D1 = -math.sqrt(math.pi/4)*1j
    D2 = -math.sqrt(math.pi/4)*1j

#    parlist = float(sys.argv[1:len(sys.argv)])

#print parlist[5]
#for eachArg in sys.argv[1:max[sys.argv]]:   
        

# Error control, no need to touch
except IndexError :
    raise IndexError("It seems you are calling this script with too few arguments.")
except ValueError :
    raise ValueError("It seems some of the arguments are not correctly formatted. "+
                     "Remember that they must be floating point numbers.")

# 2. The function giving P(k), including the necessary import statements.
#    Inside this function, you can use the parameters named in the previous step.


#scipy.special.yv(1.5, -3+0j)



def P(k) :
    global D1,D2    
    i = 1
    while i < len(parlist)-6:
        xi = parlist[i-1]
        r = k*pow(10,parlist[i])
        xi1 = parlist[i+1]
        D1new = (D1*math.pi*(r*scipy.special.jv(1.5 + 1./xi, r/xi+0j)*scipy.special.yv(0.5 + 1./xi1, r/xi1+0j) - 
        scipy.special.jv(0.5 + 1./xi, r/xi+0j)*((xi - xi1)*scipy.special.yv(0.5 + 1./xi1, r/xi1+0j) + 
        r*scipy.special.yv(1.5 + 1./xi1, r/xi1+0j)))/(2*cmath.sqrt(xi*xi1)) - 
        1j*D2*math.pi*(-r*scipy.special.yv(1.5 + 1./xi, r/xi+0j)*scipy.special.yv(0.5 + 1./xi1, r/xi1+0j) + 
        scipy.special.yv(0.5 + 1./xi, r/xi+0j)*((xi - xi1)*scipy.special.yv(0.5 + 1./xi1, r/xi1+0j) + 
        r*scipy.special.yv(1.5 + 1./xi1, r/xi1+0j)))/(2*cmath.sqrt(xi*xi1)))
        D2new = (-1j*D1*math.pi*(-r*scipy.special.jv(1.5 + 1./xi, r/xi+0j)*scipy.special.jv(0.5 + 1./xi1, r/xi1+0j) + 
        scipy.special.jv(0.5 + 1./xi, r/xi+0j)*((xi - xi1)*scipy.special.jv(0.5 + 1./xi1, r/xi1+0j) + 
        r*scipy.special.jv(1.5 + 1./xi1, r/xi1+0j)))/(2*cmath.sqrt(xi*xi1)) + 
        D2*math.pi*(r*scipy.special.jv(1.5 + 1./xi1, r/xi1+0j)*scipy.special.yv(0.5 + 1./xi, r/xi+0j) + 
        scipy.special.jv(0.5 + 1./xi1, r/xi1+0j)*((xi - xi1)*scipy.special.yv(0.5 + 1./xi, r/xi+0j) - 
        r*scipy.special.yv(1.5 + 1./xi, r/xi+0j)))/(2*cmath.sqrt(xi*xi1)))
        D1=D1new
        D2=D2new
        i+=2
        print D1,D2


#    return A * (k/k_0)**(n_s-1.)
    return pow(k,(2. - 2./xi1))*pow(abs(D2),2)




# 3. Limits for k and precision:
#    Check that the boundaries are correct for your case.
#    It is safer to set k_per_decade primordial slightly bigger than that of Class.

k_min  = 1.e-6
k_max  = 10.
k_per_decade_primordial = 2.


#for eachArg in parlist:
#    print eachArg

#
# And nothing should need to be edited from here on.
#

# Filling the array of k's
ks = [float(k_min)]
while ks[-1] <= float(k_max) :
    ks.append(ks[-1]*10.**(1./float(k_per_decade_primordial)))

# Filling the array of Pk's
for k in ks :
    P_k = P(k)
    print "%.18g %.18g" % (k, P_k)
