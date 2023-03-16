# DrillPyze
Drillstring Dynamics Library (**DrillPyze**)

This is an object oriented library that is able to simulation coupled dynamics of a drillstring in real-time and provides basic estimation capabilities for a specified parameter.

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
1) well_profile 


## Output

![Sample Output](example.png?raw=true "An example output of multiple estimation")