from jasp import *
from ase import Atom, Atoms
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf
from scipy.optimize import fsolve

# the above imports the nescessary components to solve the problem.

LC=[3,3.05,3.1,3.15,3.2,3.25,3.3,3.35,3.4,3.45,3.5,3.55,3.6,3.65,3.7,3.75,3.8,3.85,3.9,3.95,4,4.05,4.1,4.15,4.2,4.25,4.3,4.35,4.4,4.45,4.5,4.55,4.6,4.65,4.7,4.75,4.8,4.85,4.9,4.95,5,5.05,5.1,5.15,5.2,5.25,5.3,5.35,5.4,5.45,5.5,5.55,5.6,5.65,5.7,5.75,5.8,5.85,5.9,5.95,6]
ones=[]
energies=[]
volumes=[]
ready=True
for a in LC:
    atoms=Atoms([Atom('Cl',(0,0,0))],
                cell=.5*a*np.array([[1,1,0],
                                    [0,1,1],
                                    [1,0,1]]))
    with jasp('bulk/Cl-{0}'.format(a),
              xc='PBE',
              encut=350,
              kpts=(8,8,8),
              atoms=atoms) as calc:
        try:
                e=atoms.get_potential_energy()
                energies.append(e)
                ones.append(1)
        except (VaspSubmitted,VaspQueued):
                ready=False
if not ready:
    import sys; sys.exit()
for m in LC:
    with jasp('bulk/Cl-{0}'.format(m)) as calc:
        atoms=calc.get_atoms()
        volumes.append(atoms.get_volume())

# all of the above sets up the problem for the code below. It simply creates an FCC Cl. Changing the above can easily make a different structure.
print volumes
print energies
eos=EquationOfState(volumes,energies) #Performs the EOS calcs
v0,e0,B0=eos.fit() #Gives us the values at minimum energy
CI=.95 #Determines the level of Confidence Interval
def func(x): #sets up the function to be solved
    return CI-erf(.701707*x) #Proposes the function dictating the normal random variable
z=fsolve(func,0) #solves for the normal random variable
v=np.power(volumes,-.3333333) #sets up the variable as the one in the EOS
fit3=np.polyder(((np.poly1d(np.polyfit(v,energies,3)))),2) #Creates an equation that solves for the bulk modulus in the same manner as the EOS
BM=fit3(v)/9*v**5 #solves for the bulk modulus corresponding to each volume value
var=[np.std(volumes),np.std(energies),np.std(BM)] #solves for the standard deviation of each set of data of interest
SS=[np.count_nonzero(volumes),np.count_nonzero(energies),np.count_nonzero(BM)] #determines the sample size of each set of data of interest.
R=[z*var[0]*SS[0]**-.5,z*var[1]*SS[1]**-.5,z*var[2]*SS[2]**-.5] #Determines the radius of the confidence interval
Limits=[v0-R[0],v0+R[0],e0-R[1],e0+R[1],B0-R[2],B0+R[2]] #sets up each confidence interval
print 'The {0} confidence interval around the volume at minimum energy of the selected structure is between {1} and {2} A^3.The {0} confidence interval around the minimum energy of the selected structure is between {3} and {4} eV.The {0} confidence interval around the bulk modulus of the selected structure is between {5} and {6} eV/A^3  '.format(CI,Limits[0],Limits[1],Limits[2],Limits[3],Limits[4],Limits[5]) #displays the results#solves for the standard deviation of each set of data of interest
