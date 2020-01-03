.. _flm:

flm
===

Usage
-----

The FLM class of flm.py is the primary way one should generate realizations of fractional Levy motion. Here is a simple example which shows how to generate realizations and how to perform some time series analysis on them.

.. code-block:: python

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
 
.. autoclass:: flm.FLM
   :members: __init__, generate_realizations, plot_marginal, autocorrelation, plot_autocorrelation, plot_trajectory, plot_msd
   

