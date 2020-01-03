flm
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/bencoscia/flm.svg?branch=master)](https://travis-ci.org/bencoscia/flm)
[![Documentation Status](https://readthedocs.org/projects/flm/badge/?version=latest)](https://flm.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/bencoscia/flm/branch/master/graph/badge.svg)](https://codecov.io/gh/bencoscia/flm/branch/master)

Documentation and testing are under development!

Simulate fractional Levy motion

Check out the basic functionality with the following commands:

```
import flm

# Initialize fractional levy motion with trajectories of length 1000, hurst parameter = 0.4, alpha = 1.75
fractional_levy_motion = flm.FLM(1000, 0.4, 1.75)

# generate 10 independent trajectories
fractional_levy_motion.generate_realizations(10)

# plot the marginal distribution of hops
fractional_levy_motion.plot_marginal(show=False, bounds=(-5, 5))

# calculate and plot the autocorrelation function of the increments
fractional_levy_motion.plot_autocorrelation()

# plot the first 2 realizations
fractional_levy_motion.plot_trajectory([0, 1])

# calculate and plot the mean squared displacement of the trajectories
fractional_levy_motion.plot_msd(show=True)
```

### Copyright

Copyright (c) 2020, Benjamin Coscia


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
