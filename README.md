# Reconstruction-type-interaction-networks

### Brief description
Python code for simulating multi-type pathogen dynamics and for sampling random cross-sectional datasets from the ensuing endemic equilibrium. The underlying type interaction networks are subsequently reconstructed from co-occurrence data via statistical network inference in R or SAS. The corresponding codes were used to evaluate the performance of network inference models as proposed in "Reconstructing heterogeneous pathogen interactions from co-occurrence data via statistical network inference", *bioRXiv, 2021 (doi: 10.1101/2021.11.15.468692)*.

### Python prerequisites 
- numpy
- scipy
- pandas
- matplotlib

### R prerequisites
- IsingFit
- IsingSampler
- boot

### Use of files under subdirectory `scripts`
`parameter_sets_generation.py`: This script generates random parameter sets used by `equilibrium_simulation.py`, `equilibrium_simulation_assortative.py`, and `equilibrium_simulation_twice_connected.py`. 

`equilibrium_simulation.py`, `equilibrium_simulation_assortative.py`, and `equilibrium_simulation_twice_connected.py`: These scripts simulate the equilibrium of multi-type dynamical systems corresponding to the inputted generated parameter set, using `infection_state.py`, and either `ode_model.py` or `ode_model_assortative.py`.

`validation_inference_performance.R`: This script validates the performance of statistical network inference.

`inference_GEE.sas`: This script performs inference with Generalized Estimating Equations models (GEE).

### Authors
Irene Man<sup>*,1,2</sup>, Elisa Benincà<sup>1</sup>, Mirjam E. Kretzschmar<sup>2</sup>, Johannes A. Bogaards<sup>1,3</sup>
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
* Corresponding author
