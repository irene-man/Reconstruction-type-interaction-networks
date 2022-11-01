# Reconstruction-type-interaction-networks

### Brief description
Python code for simulating the equilibrium of infection disease dynamics of multi-type pathogens and, subsequently, sampling random cross-sectional datasets from the equilibrium. The code is used to compute the performance of statistical network inference models in reconstructing pathogen type interaction networks from co-occurrence data proposed in "Reconstructing heterogeneous pathogen interactions from co-occurrence data via statistical network inference", *bioRXiv, 2021 (doi: 10.1101/2021.11.15.468692)*.

### Prerequisites
- numpy
- scipy
- pandas
- matplotlib

### Use
`parameter_sets_generation.py`: This script generates random parameter sets used by `equilibrium_simulation.py`, `equilibrium_simulation_assortative.py`, and `equilibrium_simulation_twice_connected.py`. 

`equilibrium_simulation.py`, `equilibrium_simulation_assortative.py`, and `equilibrium_simulation_twice_connected.py`: These scripts simulate the equilibrium of multi-type dynamical systems corresponding to the inputted generated parameter set, using `infection_state.py`, and either `ode_model.py` or `ode_model_assortative.py`. 

### Authors
Irene Man<sup>1,2</sup>, Elisa Beninca<sup>1</sup>, Mirjam E. Kretzschmar<sup>2</sup>, Johannes A. Bogaards<sup>1,3</sup>
<table>
  <tr>
    <td>1</th>
    <td>Centre for Infectious Disease Control, National Institute for Public Health and the Environment, Bilthoven</th>
    <td>The Netherlands</td>
  </tr>
  <tr>
    <td>2</td>
    <td>Julius Centre, UMC Utrecht, Utrecht University, Utrecht, The Netherlands</td>
    <td>The Netherlands</td>
  </tr>
  <tr>
    <td>3</th>
    <td>Department of Epidemiology & Data Science, Amsterdam UMC, Amsterdam, The Netherlands</th>
    <td>The Netherlands</td>
  </tr>
</table>
