flm
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.com/REPLACE_WITH_OWNER_ACCOUNT/flm.svg?branch=master)](https://travis-ci.com/REPLACE_WITH_OWNER_ACCOUNT/flm)
[![Documentation Status](https://readthedocs.org/projects/flm/badge/?version=latest)](https://flm.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/flm/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/flm/branch/master)

Documentation and testing are under development!

Simulate fractional Levy motion

Check out the basic functionality with the following commands:

```
import flm

# Initialize fractional levy motion
fractional_levy_motion = flm.FLM(1000, 0.4, 1.75)  # trajectories of length 1000, hurst parameter = 0.4, alpha = 1.75
fractional_levy_motion.generate_realizations(10)  # generate 10 independent trajectories
fractional_levy_motion.plot_marginal(show=False, bounds=(-5, 5))  # plot the marginal distribution of hops
fractional_levy_motion.plot_autocorrelation()  # calculate and plot the autocorrelation function of the increments
fractional_levy_motion.plot_trajectory(0)  # plot one of the realizations
fractional_levy_motion.plot_msd(show=True)  # calculate and plot the mean squared displacement of the trajectories
```

### Copyright

Copyright (c) 2020, Benjamin Coscia


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
