# DrillPyze
Drillstring Dynamics Library (**DrillPyze**)

This is an object oriented library that is able to simulation coupled dynamics of a drillstring in real-time and provides basic estimation capabilities for a specified parameter.

## References
Torsional drillstring dynamics is based on the model proposed by Aarsnes and Shor:
- Aarsnes, U.J.F. and Shor, R.J., 2018. Torsional vibrations with bit off bottom: Modeling, characterization and field data validation. Journal of Petroleum Science and Engineering, 163, pp.712-721. [link](https://www.sciencedirect.com/science/article/pii/S0920410517309075)

Estimation is based on the technique proposed by Aarsnes et al.:
- Aarsnes, U.J.F., Auriol, J., Di Meglio, F. and Shor, R.J., 2019. Estimating friction factors while drilling. Journal of Petroleum Science and Engineering, 179, pp.80-91. [link](https://www.sciencedirect.com/science/article/pii/S092041051930350X)


## Directory Structure

```

data
|   bha.xlsx - example bottom hole assembly
|   survey.xlsx - example wellbore survey
docs
├───_build/html/index.html - documentation
src
├───DrillPyze - the DrillPyze library
|   Run_Estimator_Parallel.py - example usage to run the estimator with a ThreadingPool
|   Run_Estimator.py - example usage to run the estimator as a single thread
|   Run_Sim.py - example useage to run the simulator without estimation

```


## Dependencies
Requires:
1) [well_profile](https://github.com/pro-well-plan/well_profile) - an open source library from Pro Well Plan for well path calculations
2) [numba](https://numba.pydata.org/) - an open source JIT compiler for Python


## Output

![Sample Output](example.png?raw=true "An example output of multiple estimation")



## Useful commands

### Documentation
The source files are set up so that *sphinx* can be used to automatically create user documentation.  To create the documentation, run
```
make html
```
in the `docs` directory.  This will create documentation in html format which can be found in `docs/_build/html/`.  Open `index.html` in your browser to view.


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

Further informatin on how to deploy the package to pypi can be found [here](https://betterscientificsoftware.github.io/python-for-hpc/tutorials/python-pypi-packaging/)