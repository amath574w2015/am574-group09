# am574-group09

Title:
Augmented approximate Riemann solvers for the shallow water equations with variable bathymetry.

Authors:
Xin Chen
Jacob Ortega-Gingrich

Abstract:

We describe an augmented Riemann solver for the one-dimensional shallow wa- ter equations with variable topography suggested by David George which addresses a number of the needs of applications involving flooding and small perturbations of delicate steady states. Fluid flows over varying topography, for example add a source term which the numerical method must balance with jumps in momentum and depth in order to preserve delicate steady states, such as an ocean at rest. Furthermore, applications involving flooding, such as the modeling of tsunami inundation, require a Riemann solver that can handle dry cells and preserve depth non-negativity. The augmented solver herein discussed satisfies these properties by amalgamating various aspects of existing solvers such as the HLLE solvers and f-wave approaches and the addition of a stationary steady state wave to account for the source term.  Additionally, we discuss some specific techniques which may be use to handle various challenging situations which may arise in a model such as the handling of steep shorelines. Finally, we discuss an implementation of this augmented Riemann solver for the one-dimensional shallow water equations and present a few numerical demonstrations.