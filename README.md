# Optimal Operation of LNG Refrigeration Cycles

This repository corresponds to the MATLAB simulation codes for modelling and optimizing the Mini - LNG System proposed by Nekså, 2010.
A previous mdoel of this system was developed by Leguizamon, 2016.
This works expands and improves the previous formulation.

## Work Scope

The model is first created using a pseudo modular approach, which then is used as an initial point for an equations oriented (EO) solution of the plantwide simulation.

* Thermodynamics: the solution is done by using an SRK implementation.

### Pseudo modular approach
The system is modelled by a series of flash calculations. These are performed by using the EO approach proposed by Kamath, 2010. For this solution, the initial points are given by a priori knowledge of the process or by previous results for other streams. <br>
Each flash calculations is solved using MATLAB's fmincon with default options (unless noted)

### EO Optimization
A plantwide EO solution can be achieved once there is complete information for each stream. This solution is then optimized for a range of disturbances.
**This part has not been implemented yet.**
#### Challenges
This problem is very stiff and the objective function is not convex. The following challenges arise:
* Find convergence for each point
* Make sure the optimum found is global

## Analysis
The results from the optimization are then split into sections depending on the active constraints obtained for each point.

#### Ideas for analysis
* One challenge in the areas of active constraints is to define with pression when the boundary lies. For this reason, a clustering algorithm could be implemented to predict the region in which a new point will be.
* Due to the size and stiffness of the optimization problem, global convergence for each point is not guaranteed. An SVR machine could be set up to approximate the optimum for new points.

## References
* Nekså, P., Brendeng, E., Drescher, M., & Norberg, B. (2010). Development and analysis of a natural gas reliquefaction plant for small gas carriers. Journal of Natural Gas Science and Engineering, 2(2–3), 143–149. <url>http://doi.org/10.1016/j.jngse.2010.05.001</url>
* Kamath, R. S., Biegler, L. T., & Grossmann, I. E. (2010). An equation-oriented approach for handling thermodynamics based on cubic equation of state in process optimization. Computers & Chemical Engineering, 34(12), 2085–2096. <url>http://doi.org/10.1016/j.compchemeng.2010.07.028</url>
* Leguizamón, R. A., (2016) Optimal operation of LNG refriegeration cycles. Master Thesis. NTNU. Trondheim, Norway.
<url> http://hdl.handle.net/11250/2413532</url>
