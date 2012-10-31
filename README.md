Project1
========

The documentation required by JFK for project 1.

I assumed a large part of your opinion of the code was the reduction of lines, so I tried to make as few as possible.

The goal of this project is to implement an uncertainty analysis on the estimated parameters in an equation of state, 
e.g. what is the 95% confidence interval on the minimum volume and bulk modulus.

I chose to base my calculations around a simple structure, FCC Cl;however, it should be applicable to any structure.
Up until line 36, the code is mostly to set up the problem. Changes to the structure should be made in here. The actual 
unique code as it pertains to this project begins on line 41 and ends at line 54. Outside of this interval, as long 
as the appropriate names are given to the volumes and energies vectors, it should work. I begin by using the equation 
of state solver from ase. This fits the volumes and energies to the following polynomial.

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
easier use. The drawback is that the solver is not as accurate as many literature sources or other solvers. This error
increases the closer you get to the extremes of the function's range (0 and 1). At .95, the value was off by ~.7%, which
is acceptable for most calculations, but at .99 the solver was off by ~11.4%, which is significant for many systems.
I hope to try some of the other solvers before I turn htis in to optimize. Also, this is assuming a standard normal
distribution. If a different distribution is required, this equation and possibly the solver must be changed.
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

