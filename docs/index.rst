.. flm documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Simulate Fractional Levy Motion
================================

Welcome to the documentation page for flm, a package designed to approximately
simulate fractional L\'evy motion (FLM) for negatively correlated increments. 
Simulate positively correlated increments at your own risk. 

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   api

.. _theory:

What is fractional Levy motion?
-------------------------------

It is easiest to understand FLM by first considering Brownian and fractional
Brownian motion. During Brownian motion, consecutive increments (the distance a
particle displaces) are drawn indepenently from a Gaussian distribution. The
big assumption here is that each draw from the distribution is completely
independent of previous draws. However, this is frequently not the case in
complex systems. While the distribution of increments might appear Gaussian,
it's possible that draws from the distribution are correlated. Correlated draws
from a Gaussian distribution are well described by fractional Brownian motion (FBM).
In fact, the autocorrelation function has the analytical form:

.. math::

        \gamma(k) = \dfrac{\sigma^2}{2}\bigg[|k-1|^{2H} - 2|k|^{2H} + |k+1|^{2H}\bigg]

where :math:`\sigma` is the width of the Gaussian distribution and H is the
Hurst parameter, a quantity which describes the degree of anti-correlation
between draws from the underlying Gaussian distribution. When H = 0.5, we
recover Brownian motion. When H < 0.5, we observe negatively correlated draws
from the Gaussian distribution. When H > 0.5, draws are positively correlated.
An excellent python package already exists for exactly simulating fractional
Brownian motion: see the `fbm <https://github.com/crflynn/fbm>`_ package

While FBM is quite useful, it is not always the case that the marginal distribution
of increments is Gaussian. But luckily, there is a more general class of 
distributions that might apply, called L\'evy stable distributions. Somewhat
unfortunately, the probability density function for the general L\'evy stable distribution
is only described through its characteristic function in Fourier space:

.. math::

    p_{\alpha_h, \beta}(k;\mu,\sigma) =\exp\left[i\mu k - \sigma^{\alpha_h}|k|^{\alpha_h}\left(1 - i\beta\frac{k}{|k|}\omega(k, \alpha_h)\right)\right]

where

.. math::
  \omega(k, \alpha_h) = \begin{cases}
        \tan{\frac{\pi \alpha_h}{2}} & \text{if}~\alpha_h \neq 1, 0 < \alpha_h < 2, \\
        -\frac{2}{\pi}\ln |k| & \text{if}~\alpha_h = 1
         \end{cases}
  
:math:`\alpha_h` is the index of stability or L\'evy index, :math:`\beta` is
the skewness parameter, :math:`\mu` is the shift parameter and :math:`\sigma`
is a scale parameter. 

More to come ...

.. Note thae while this package will accept values of H that are greater than 0.5,
  it has not been explicitly tested with positively correlated increments (yet). 

.. Note that to my knowledge, there are no exact
  simulation techniques for generating linear

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
