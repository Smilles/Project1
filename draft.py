from jasp import *
from ase import Atom, Atoms
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf
from scipy.optimize import fsolve
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


eos=EquationOfState(volumes,energies)
v0,e0,B0=eos.fit()
CI=.95
def func(x):
    return CI-erf(.701707*x)
z=fsolve(func,0)
v=np.power(volumes,-.3333333)
fit3=np.polyder(((np.poly1d(np.polyfit(v,energies,3)))),2)
BM=fit3(v)/9*v**5
var=[np.std(volumes),np.std(energies),np.std(BM)]
SS=[np.count_nonzero(volumes),np.count_nonzero(energies),np.count_nonzero(BM)]
R=[z*var[0]*SS[0]**-.5,z*var[1]*SS[1]**-.5,z*var[2]*SS[2]**-.5]
Limits=[v0-R[0],v0+R[0],e0-R[1],e0+R[1],B0-R[2],B0+R[2]]
print 'The {0} confidence interval around the volume at minimum energy of the selected structure is between {1} and {2} A^3.The {0} confidence interval around the minimum energy of the selected structure is between {3} and {4} eV.The {0} confidence interval around the bulk modulus of the selected structure is between {5} and {6} eV/A^3  '.format(CI,Limits[0],Limits[1],Limits[2],Limits[3],Limits[4],Limits[5])


#assumes standard Z
#print '''
#v0={0} A^3
#e0={1} eV
#B={2}eV/A^3'''.format(v0,e0,B)
#t1=np.mean(volumes)#t2=np.std(volumes)
#t3=np.count_nonzero(volumes)
#t5=np.sqrt(t3)
#lb=v0-1.96*t2/t5
#ub=v0+1.96*t2/t5
#print 1.96*t2/t5
#print lb
#print ub
#y9=np.divide(energies,volumes)
#m9=[energies,volumes,y9]
#print '| Lattice Constant | Total Energy (eV) | Volumes (A^3) | BM? |'
#for w1, w2, w3, w4 in zip(LC,energies,volumes,y9):
#    print '| {0} | {1} | {2} | {3}|'.format(w1, w2,w3,w4)
#eos.plot()
#plt.show()
#h=[ones,volumes,energies]
#h2=zip(*h)
#print np.dot(h,h2)
#vol1=np.power(volumes,-.33333)
#vol2=np.power(volumes,-.66666)
#vol3=np.power(volumes,-1)
#Xprime=[ones,vol1,vol2,vol3]
#print Xprime
#X=zip(*Xprime)
#XX=np.dot(Xprime,X)
#Xy=np.dot(Xprime,energies)
#XXinv=np.linalg.inv(XX)
#B=np.dot(XXinv,Xy)
#x0=[1,v0**-.33333,v0**-.66666,v0**-1]
#u=np.dot(x0,B)
#v=np.var(energies)
#n4=[1,8,275]
#n5=[[.214653,-.007491,-.000340],
#    [-.007491,.001671,-.000019],
#    [-.000340,-.000019,.0000015]]
#n6=[[1,0,0],
#    [8,0,0],
#    [275,0,0]]
#n7=np.cross(n6,n5)
#n8=np.dot(n5,n4)
#print n8
