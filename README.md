Project1
========

The documentation required by JFK for project 1.

I assumed a large part of your opinion of the code was the reduction of lines, so I tried to make as few as possible.

The goal of this project is to implement an uncertainty analysis on the estimated parameters in an equation of state, 
e.g. what is the 95% confidence interval on the minimum volume and bulk modulus.

draft.py contains a lot of my pre work and attempts that lead to my final project. 
SMI-Project-1.py is what was going to be my final project without the drafting work.
jfkfunc.py is my final project. I made it into an executable function where all that needs to be put in are the volume and
energy vectors as well as the desired confidence interval. At the end of the file, I run the function as it pertains to the 
FCC Cl system I made in SMI-Project-1.py to show the use.

I chose to base my calculations around a simple structure, FCC Cl;however, it should be applicable to any structure.
Up until line 36 of the my initial code, the code is mostly to set up the problem. Changes to the structure should be made in here. The actual 
unique code as it pertains to this project begins on line 41 and ends at line 54 of that file. Outside of this interval, as long 
as the appropriate names are given to the volumes and energies vectors, it should work. I represent that by making an 
executable function in my final write up. I begin by using the equation of state solver from ase. This fits the volumes 
and energies to the following polynomial.

E= C1 + C2*t + C3*t^2 + C4*t^3
        Where
              E=Energy in eV
              t=(1/Volume)^(1/3) in A^-1
              C denotes constants
This returns the minimum energy (e0), the volume at minimum energy (v0), and the bulk modulus at minimum energy (B0).
Next I define the variable CI as the confidence interval. It is currently set to .95, but can be changed to any value 
between 0-1.
The next few lines solve for the normal random variable corresponding to the denoted confidence interval. I'm actually 
quite proud of this part. Initially my code made the user look up the normal random variable, but this allows for much
easier use. The drawback is that the solver is not quite as accurate as many literature sources or other solvers. At .95,
the value was off by ~.7%, which is acceptable for most calculations. At .99 the solver was off by ~.6%, which should 
still be acceptable. Also, this is assuming a standard normal distribution. If a different distribution is required, 
this equation and possibly the solver must be changed.
Next I define the variable that is used in the EOS fit for future use.
The next few lines are used to calculate the bulk modulus (labelled BM *snicker) using the same methods as the EOS.
The difference here is that I now have a vector of the bulk modulus that corresponds to each value of the volume and 
the energy. The code now contains all of the data needed to calculate the confidence intervals.
The next 2 lines calculate the standard deviation and sample size, respectively, for each data set for use in the 
CI equation.
The next line uses the results above and the data to calculate the radius of the confidence interval for each data set.
The next line takes the results of the previous line and applies the radius to each of the values determined earlier by 
the EOS fit.
The remainder of the code simply tells the results of the application in a clear format.

Below is the example code I ran, which uses the results of SMI-Project-1.py pertaining to the FCC Cl system to solve.

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


The result of this code is as follows

The 0.95 confidence interval around the volume at minimum energy of the selected structure is between [ 18.003] and 
[ 25.021] A^3.The 0.95 confidence interval around the minimum energy of the selected structure is between [-1.366] and 
[-0.474] eV.The 0.95 confidence interval around the bulk modulus of the selected structure is between [-0.111] and 
[ 0.487] eV/A^3  




