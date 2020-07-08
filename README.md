# thesis
Competition and predator-prey models in both the homogeneous and heterogeneous cases.

The individual behavior of the system can be simulated with Monte Carlo algorithms. Here you can find the codes that implement the algorithm. 
When the sample size is large, the dynamics can be descibed by the mean field equation. Here I show how to solve the mean field equation numerically in the heterogeneous case. 

The solutions of the predator-prey model present an oscillatory behavior. First, I focous on the regularity of the period and then, studying Langevin equation, I recover the behaviour of the simulations and the numerical solutions. 


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
