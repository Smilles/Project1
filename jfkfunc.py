from jasp import *
from ase.utils.eos import EquationOfState
import numpy as np
from scipy.special import erf
from scipy.optimize import fsolve


def eosanal(volumes,energies,CI):
    eos=EquationOfState(volumes,energies) #Performs the EOS calcs
    v0,e0,B0=eos.fit() #Gives us the values at minimum energy
    def func(x): #sets up the function to be solved
        return CI-erf(.701707*x)
    #Proposes the function dictating the normal random variable
    z=fsolve(func,0) #solves for the normal random variable
    v=np.power(volumes,-.3333333) #sets up the variable as the one in the EOS
    fit3=np.polyder(((np.poly1d(np.polyfit(v,energies,3)))),2)
    #Creates an equation that solves for the bulk modulus in the same manner as the EOS
    BM=fit3(v)/9*v**5
    #solves for the bulk modulus corresponding to each volume value
    var=[np.std(volumes),np.std(energies),np.std(BM)]
    #solves for the standard deviation of each set of data of interest
    SS=[np.count_nonzero(volumes),np.count_nonzero(energies),np.count_nonzero(BM)]
    #determines the sample size of each set of data of interest.
    R=[z*var[0]*SS[0]**-.5,z*var[1]*SS[1]**-.5,z*var[2]*SS[2]**-.5]
    #Determines the radius of the confidence interval
    Limits=[v0-R[0],v0+R[0],e0-R[1],e0+R[1],B0-R[2],B0+R[2]]
    #sets up each confidence interval
    print 'The {0} confidence interval around the volume at minimum energy of the selected structure is between {1} and {2} A^3.The {0} confidence interval around the minimum energy of the selected structure is between {3} and {4} eV.The {0} confidence interval around the bulk modulus of the selected structure is between {5} and {6} eV/A^3  '.format(CI,Limits[0],Limits[1],Limits[2],Limits[3],Limits[4],Limits[5]) #displays the results
    return;
print eosanal([6.7500000000000009, 7.0931562499999989, 7.447750000000001, 7.813968749999999, 8.1920000000000002, 8.5820312500000018, 8.9842499999999959, 9.3988437499999975, 9.825999999999997, 10.265906249999999, 10.71875, 11.184718749999998, 11.663999999999998, 12.156781250000002, 12.663250000000003, 13.183593749999996, 13.717999999999996, 14.266656249999997, 14.829749999999995, 15.407468749999998, 15.999999999999998, 16.607531250000001, 17.230249999999995, 17.868343749999998, 18.521999999999998, 19.19140625, 19.876749999999998, 20.578218749999994, 21.296000000000003, 22.030281250000002, 22.781250000000004, 23.549093749999994, 24.333999999999985, 25.136156249999999, 25.955750000000002, 26.79296875, 27.647999999999996, 28.521031249999993, 29.41225, 30.32184375000001, 31.250000000000007, 32.196906249999984, 33.162749999999996, 34.147718750000003, 35.152000000000008, 36.175781249999993, 37.219250000000002, 38.282593749999982, 39.366000000000021, 40.46965625, 41.59375, 42.738468750000003, 43.903999999999982, 45.090531250000026, 46.298250000000017, 47.527343750000007, 48.777999999999992, 50.050406249999995, 51.344749999999998, 52.661218750000003, 54.000000000000007],[6.897374, 5.974071, 5.143004, 4.396682, 3.727342, 3.12759, 2.591493, 2.113506, 1.688708, 1.312134, 0.979491, 0.686502, 0.429396, 0.20462, 0.008959, -0.160481, -0.306446, -0.431322, -0.537025, -0.625937, -0.699817, -0.760306, -0.808954, -0.84708, -0.875895, -0.896482, -0.90987, -0.916977, -0.918637, -0.915588, -0.908672, -0.898053, -0.884525, -0.868575, -0.850635, -0.831079, -0.810219, -0.78833, -0.765678, -0.742508, -0.719026, -0.695402, -0.671729, -0.648126, -0.624743, -0.601677, -0.579013, -0.556806, -0.5351, -0.513934, -0.493357, -0.473414, -0.454138, -0.43556, -0.417677, -0.400595, -0.384044, -0.368179, -0.352994, -0.338464, -0.324563],.95)
