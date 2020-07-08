# thesis
In my Master's thesis, the main stochastic models in population biology that I have studied are the competition and predator-prey systems in both the homogeneous and heterogeneous cases.

I have focused on the individual behavior of the system that can be simulated with Monte Carlo algorithms. It is well known that when the sample size is large, the computational cost of the Monte Carlo algorithm is high and so it is convenient to describe the dynamics with the mean-field equations. I have shown how to derive this equations, that are either ordinary differential equations or partial differential equations and how to solve them.

In the last part, I have studied the oscillatory behavior of the solutions in a predator-prey model focusing on the different nature of oscillation between the simulations and the numerical solutions and on the regularity of the period.


Folder description:
# homogeneous
Simulated solutions for both the one and two species model in the homogeneous case, obtained with the Monte Carlo algorithm.
# nonHomogeneous1Species
# simula
Simulated solutions for the one species model in the heterogeneous case, obtained with the Monte Carlo algorithm.
# numeric
Numerical solutions obtained combining a finite difference method in space and implicit or explicit methods in time. 

# nonHomogeneous2Species
# simula
Simulated solutions for the one species model in the heterogeneous case, obtained with the Monte Carlo algorithm.
# numeric
Numerical solutions obtained combining a finite difference method in space and implicit or explicit methods in time. 

# predatorPrey
# simula
Simulated solutions for the predator prey model, obtained with the Monte Carlo algorithm.
# twoD
Two dimensional predator prey model
# regularity
Regularity of the periods in the simulation as the parameters vary. 
The regularity values are related to the choice of the competition parameters.  
# oscillations
Signals with and without noise.
The signal with noise behaves as the simulations while the one without noise behaves as the numerical solution of a predator prey model.  
