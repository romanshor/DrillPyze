# DrillPyze - The Drillstring Dynamics Library

This is an object oriented library that is able to simulation coupled dynamics of a drillstring in real-time and provides basic estimation capabilities for a specified parameter.

## References
Torsional drillstring dynamics is based on the model proposed by Aarsnes and Shor:
- Aarsnes, U.J.F. and Shor, R.J., 2018. Torsional vibrations with bit off bottom: Modeling, characterization and field data validation. Journal of Petroleum Science and Engineering, 163, pp.712-721. [link](https://www.sciencedirect.com/science/article/pii/S0920410517309075)

Estimation is based on the technique proposed by Aarsnes et al.:
- Aarsnes, U.J.F., Auriol, J., Di Meglio, F. and Shor, R.J., 2019. Estimating friction factors while drilling. Journal of Petroleum Science and Engineering, 179, pp.80-91. [link](https://www.sciencedirect.com/science/article/pii/S092041051930350X)

An initial bit-rock interaction model can be found in Auriol et al.: 
- Auriol, J., Aarsnes, U.J.F., and Shor, R.J., "Self-Tuning Torsional Drilling Model for Real-Time Applications," 2020 American Control Conference (ACC), Denver, CO, USA, 2020, pp. 3091-3096.  [link](https://ieeexplore.ieee.org/abstract/document/9147242)

A summary of surface control systems may be found in Aarsnes et al.: 
- Aarsnes, U.J.F., di Meglio, F., and Shor, R.J. Benchmarking of Industrial Stick-Slip Mitigation Controllers, IFAC-PapersOnLine, Volume 51, Issue 8, 2018, Pages 233-238.  [link](https://www.sciencedirect.com/science/article/pii/S2405896318307122)

## Directory Structure

```

data
|   bha.xlsx - example bottom hole assembly
|   survey.xlsx - example wellbore survey
|   data.csv - 30 seconds of example data based on the included BHA and survey
|   data_long.csv - 240 seconds of example data
docs
├───build/html/ - on running ./docs_generator.sh, local copy of documentation
├───src/ - source files for sphinx to generate documentation
src
├───DrillPyze - the DrillPyze library
|   Create_Data.py - create an example dataset for use with the estimator
|   Run_Estimator_Parallel.py - example usage to run the estimator with a ThreadingPool
|   Run_Estimator.py - example usage to run the estimator as a single thread
|   Run_Sim.py - example useage to run the simulator without estimation

```


## Dependencies
This package requires the following libraries:
1) [well_profile](https://github.com/pro-well-plan/well_profile) - an open source library from Pro Well Plan for well path calculations
2) [numba](https://numba.pydata.org/) - an open source JIT compiler for Python


## Output
Utilizng the provided plotting utilities, this package can create an animation of drillstring dynamics behavior in the following visual.  
![Sample Output](example.png?raw=true "An example output of multiple estimation")

## Function Documentation
The documentation of the `DrillPyze` package can be found on [github pages](https://romanshor.github.io/DrillPyze).  Example usage is available as `Run_Sim.py` and `Run_Estimator.py`.  A python notebook with example usage is coming soon.

## Using the library
The library has the following components:
1) Steady-state / Static Torque and Drag
2) Steady-state Hydraulics
3) Rotational Dynamics

### Using the Torque and Drag Module
First, the BHA and wellpath need to initialized.  An example BHA and wellbore survey are available in the data directory.
```python
well1 = well.wellbore({'xlsx':'../data/survey.xlsx'}, casing_depth, bit_size)
bha1 = drillstring.BHA('BHA 1', '../data/bha.xlsx')
```
Next, the BHA must be discretized within the wellpath and then the wellpath must be interpolated so that each drillstring component has an inclinationa and azimuth.  The bit depth should be less than the maximum depth of the wellpath.

```python
bha1.discretize(10, bit_depth)
well1.get_interpolated(bha1.x)
```
Finally, static torque and drag calcualtion can be carried out.
```python
(Ft, Ft_pu, Ft_so, Fn, M) = statics.normal_forces(well1, bha1, mu=mu_f)
```

Torque and Drag calculations are based on:
- Johancsik, C.A., Friesen, D.B., and Rapier Dawson. "Torque and Drag in Directional Wells-Prediction and Measurement." J Pet Technol 36 (1984): 987–992
- Sheppard, M. C., Wick, C., and T. Burgess. "Designing Well Paths to Reduce Drag and Torque." SPE Drill Eng 2 (1987): 344–350.


### Using the Hydraulics Module
Steady state hydraulic calcuations for *Newtonian Fluids* may be performed using the hydraulics module.  Frictional pressure drop may be calculated via:
```python
hydraulics.frictional_pressure_drop(rho, Q, mu, L, Do, Di)
```
Utilities are also included to compute flow velocity and Reynolds Number.

Other fluid rheologies may be implemented in a future release.

### Using the Rotational D?ynamics Module
The rotational dynamics module allows for transient simulation of drillstring rotation.  Similar to the torque & drag module, first, the BHA and wellpath need to be imported:
```python
well1 = well.wellbore({'xlsx':'../data/survey.xlsx'}, casing_depth, bit_size)
bha1 = drillstring.BHA('BHA 1', '../data/bha.xlsx')
```

Next, the simulation needs to be set up and initialized.
```python
sim = dynamics.Simulation('Test Simulation', well1, bha1)
sim.initialize_simulation(bit_depth, simulation_length)
```

If plotting is required, set up the figure
```python
(fig, axes) = util.setup_plotting(live=PLOT_LIVE)
```

Next, a loop should be set up to step the simulaton forward.  The simulation may be stepped forward by an arbitrary deltaT (milliseconds to seconds).  During the step, no inputs or setpoints may be changed, however, the simulation does include a top drive control system (operating in stiff configuration by default; Soft Torque and Ztorque will be added in the future).  Inside the lop, the step function should be called, and if plotting is desired, a call to the plot funciton.
```python
sim.step(1/plot_frequency, TD_RPM_SP=TD_Setpoint, bit_constant=rock_constant)
util.plot(fig, axes, sim, bha1, M)
```


## Useful commands
The following are notes about useful commands for generating and posting documentation and for creating a distributable package.

### Documentation
The source files are set up so that *sphinx* can be used to automatically create user documentation.  To create the documentation, run
```
./docs_generator.sh
```
in the main directory.  This will create documentation in html format which can be found in `docs/build/`.  Open `index.html` in your browser to view.

This package also uses `github actions` to autogenerate updated documentation on any `push` to the main branch of this repository. 

The workflow was adapted from the write-up [here](https://betterprogramming.pub/quickly-generate-documentation-with-sphinx-cli-tools-99db0cb5994c) and [here](https://betterprogramming.pub/three-ways-to-publish-sphinx-documentation-c276ffb22e5b)

### Github actions
Documentation is automatically recompiled and updated on the [DrillPyze](https://romanshor.github.io/DrillPyze) page using the github action found in `.github/workflows`.  This workflow will generate the documentation, push the `docs/build` director to the `gh-pages` branch and will then publish the sphinx documentation.

### Distributable Library
To create a distributable library, first run
``` 
python setup.py check
```
To check if there are any issues.  Then run
```
python setup.py sdist
```
To create a distributable package (`DrillPyze-<VERSION>.tar.gz`) in the `dist/` directory.

Further information on how to deploy the package to pypi can be found [here](https://betterscientificsoftware.github.io/python-for-hpc/tutorials/python-pypi-packaging/)