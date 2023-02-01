# DC_dFBA
Supporting codes of the paper: Nonlinear programming reformulation of dynamic flux balance analysis models.

The paper describes a nonlinear programming formulation for solving dynamic Flux Balance Analysis models. The
methodology consists of incorporating a KKT reformulation of a parsimonious FBA problem and applying a collocation technique to the ODE system. Please refer to the paper for further information.

## Case studies:
The files here are organized by the case studies solved on the paper with different metabolic networks models: 
1. *E. coli* core
2. *E. coli* iJR904
3. *E. coli* iJO1366
4. *S. cerevisiae* iND750
5. Yeast 8.3

## Scripts:
- main.jl: main file that defines all the parameters and initial conditions for the problem.
- pFBA_KKT_flux.jl: solves the DC_dFBA optimization problem in JuMP.
- plotar.jl: creates figures in PDF from the output.

## Julia and packages versions:
- julia 1.5.3
- JuMP v0.21.4
- Ipopt v0.7.0
- HSL MA27