#!/usr/bin/env python

""" Generate Linear Fractional Stable Noise
"""

import numpy as np
from scipy.stats import levy_stable
import matplotlib.pyplot as plt
from flm import timeseries, stats
import tqdm


class FLM:

    def __init__(self, n, hurst, alpha, scale=1, correct_hurst=False, truncate=None, correct_truncation=True, m=256,
                 M=6000, C=1):
        """ Generate realizations of fractional levy motion, also know as linear fractional stable motion

        :param n: length of each trajectory
        :param hurst: Hurst parameter. Also known as the self-similarity parameter
        :param alpha: the tail-exponent of the stable distribution (between 0 and 2). Lower alpha = heavier tails
        :param m: 1/m is the mesh size
        :param M: kernel cut-off parameter
        :param C: normalization parameter
        :param scale: scale parameter of Levy distribution
        :param correct_hurst: Correct the input Hurst parameter so the output correlation structure is consistent with \
        the analytical autocorrelation.
        :param truncate: the largest value of emissions
        :param correct_truncation: Correct the truncation parameter so the max value of emissions is close to truncate

        :type n: int
        :type hurst: float
        :type alpha: float
        :type m: int
        :type M: int
        :type C: float
        :type scale: float
        :type correct_hurst: bool
        :type truncate: NoneType or float
        :type correct_truncation: bool
        """

        self.truncate = truncate

        if self.truncate is not None:

            raise NotImplementedError('Truncation of the levy distribution has not yet been implemented! Check back'
                                      'later!')

        # TODO: need to automate database creation
        # if truncate is not None and correct_truncation:

            # trunc = TruncateLevy()
            # # H values recorded in database are not the corrected values. They were corrected in order to get the output
            # # truncation value, but you must read the database using the uncorrected H value
            # self.truncate = trunc.interpolate(H, alpha, truncate, scale)
            #
            # if not correct_hurst:
            #     raise Warning('You are correcting the truncation parameter, but not correcting the Hurst parameter. If '
            #                   'you corrected the Hurst parameter before passing it to FLM, then you can ignore this '
            #                   'message. If you did not already correct the Hurst parameter, the truncation correction '
            #                   'will not be accurate.')

        # else:
        #     self.truncate = truncate

        # TODO: automate database creation
        if correct_hurst:

            raise NotImplementedError('The Hurst correction yet been implemented! Check back later!')

            # # Interpolate a database of input and output H parameters
            # interpolator = HurstCorrection()
            # H = interpolator.interpolate(H, alpha)

        nM = n + M
        pow = np.log(nM) / np.log(2)
        newM = int(2 ** np.ceil(pow) - n)

        if newM != M:
            print('I adjusted M from %d to %d so that m * (M + n) is a power of 2' % (M, newM))
            M = newM

        self.H = hurst
        self.alpha = alpha
        self.m = m
        self.M = M
        self.N = n
        self.Na = m * (M + n)

        if alpha < 0 or alpha > 2:
            raise Exception('Alpha must be greater than 0 and less than or equal to 2!')

        mh = 1 / m
        d = self.H - 1 / self.alpha
        t0 = np.linspace(mh, 1, m) ** d
        t1 = np.linspace(1 + mh, M, int((M - (1 + mh)) / mh) + 1)
        t1 = t1 ** d - (t1 - 1) ** d
        self.A = mh ** (1 / alpha) * np.concatenate((t0, t1))
        self.C = C * (np.abs(self.A) ** alpha).sum() ** (-1 / alpha)
        self.A *= self.C
        self.A = np.fft.fft(self.A, n=self.Na)

        self.realizations = None
        self.noise = None
        self.scale = scale
        self.acf = None
        self.autocov = None

    def generate_realizations(self, n, progress=True, truncated_distribution=None):
        """ Generate realizations of fractional levy motion

        :param n: Number of realizations to generate
        :param progress: show progress bar
        :param truncated_distribution: TruncatedLevyDistribution object from flm.sampling (not implemented)

        :type n: int
        :type progress: bool
        """
        np.random.seed(1)
        self.noise = np.zeros([n, self.N])

        for i in tqdm.tqdm(range(n), disable=(not progress)):

                self.noise[i, :] = self._realization(truncated_distribution)  # generate fractional levy noise

        self.realizations = np.cumsum(self.noise, axis=1)  # cumulative sum to get fraction levy motion trajectories

    def _realization(self, truncated_distribution):

        if self.alpha == 2:

            z = np.random.normal(0, scale=self.scale, size=self.Na)

        else:

            if self.truncate is not None:

                if truncated_distribution is not None:

                    z = truncated_distribution.sample(self.Na)

                else:

                    z = sampling.truncated_levy_distribution(self.truncate, self.alpha, self.scale, self.Na)

            else:

                z = levy_stable.rvs(self.alpha, 0, loc=0, scale=self.scale, size=self.Na)

        z = np.fft.fft(z, self.Na)

        w = np.fft.ifft(z * self.A, self.Na).real

        return w[:self.N * self.m:self.m]

    def plot_marginal(self, bounds=(-.5, .5), bins=50, show=False):
        """ Plot a histogram of the marginal distribution of increments, with the expected PDF overlayed on top of it

        :param bounds: largest increments to be included in histogram
        :param bins: number of bins in histogram of empirical distribution
        :param show: show the plot when done

        :type bounds: tuple
        :type bins: int
        :type show: bool
        """

        x = np.linspace(bounds[0], bounds[1], 1000)

        hist, bin_edges = np.histogram(self.noise.flatten(), bins=bins, range=bounds, density=True)

        # account for part of PDF that is chopped off. Using density=True makes hist sum to 1
        area_covered = levy_stable.cdf(bounds[1], self.alpha, 0, loc=0, scale=self.scale) - \
                       levy_stable.cdf(bounds[0], self.alpha, 0, loc=0, scale=self.scale)

        hist *= area_covered

        # plot bars. Can't use plt.hist since I needed to modify the bin heights
        bin_width = bin_edges[1] - bin_edges[0]
        bin_centers = [i + bin_width / 2 for i in bin_edges[:-1]]
        plt.figure()
        plt.bar(bin_centers, hist, width=bin_width)
        plt.plot(x, levy_stable.pdf(x, self.alpha, 0, loc=0, scale=self.scale), '--', color='black', lw=2)

        # formatting
        plt.xlabel('Step Size', fontsize=14)
        plt.ylabel('Frequency', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()

        if show:
            plt.show()

    def autocorrelation(self):
        """ Calculate the autocorrelation function of realizations generated by generate_realizations(). This will
        populate self.acf which is a numpy array with shape (nT, n) where nT is the number of solute trajectories and n
        is the length of each trajectory.
        """

        # calculate acf of each trajectory
        ntraj = self.noise.shape[0]
        try:
            self.acf = np.zeros([ntraj, self.N - 1])
            for i in range(ntraj):
                self.acf[i, :] = timeseries.acf(self.noise[i, :])
        except ValueError:
            length = timeseries.acf(self.noise[0, :]).size
            self.acf = np.zeros([ntraj, length])
            for i in range(ntraj):
                self.acf[i, :] = timeseries.acf(self.noise[i, :])

    def autocovariance(self):

        ntraj = self.noise.shape[0]
        try:
            self.autocov = np.zeros([ntraj, self.N - 1])
            for i in range(ntraj):
                self.autocov[i, :] = timeseries.autocovariance(self.noise[i, :])
        except ValueError:
            length = timeseries.autocovariance(self.noise[0, :]).size
            self.autocov = np.zeros([ntraj, length])
            for i in range(ntraj):
                self.autocov[i, :] = timeseries.autocovariance(self.noise[i, :])

    def plot_autocorrelation(self, max_k=25, nboot=200, confidence=68.27, show=False, overlay=False, label=None,
                             fontsize=14):
        """ Calculate the mean autocorrelation function, if it hasn't been done already, generate confidence intervals
        via statistical bootstrapping, then plot.

        :param max_k: maximum lag time to plot
        :param nboot: number of bootstrap trials for generating confidence intervals
        :param confidence: confidence interval of shaded error region (percent)
        :param show: show the plot when done
        :param overlay: If True, a new figure will not be created and it will be plotted on the active figure.
        :param label: legend label
        :param fontsize: size of font an axis ticks and labels

        :type max_k: int
        :type nboot: int
        :type confidence: float
        :type show: bool
        :type overlay: bool
        :type label: str
        :type fontsize: int
        """

        if self.acf is None:
            self.autocorrelation()

        ntraj = self.acf.shape[0]

        # bootstrap
        boot = np.zeros([nboot, self.acf.shape[1]])
        for i in range(nboot):
            ndx = np.random.randint(ntraj, size=ntraj)
            boot[i, :] = self.acf[ndx, :].mean(axis=0)

        errorbars = stats.confidence_interval(boot, confidence)

        if not overlay:
            plt.figure()
        avg = boot.mean(axis=0)
        plt.plot(np.arange(self.acf.shape[1]), avg, lw=2, label=label)
        plt.fill_between(np.arange(self.acf.shape[1]), avg + errorbars[0, :], avg - errorbars[1, :], alpha=0.25)

        # formatting
        plt.xlim(-0.5, max_k)
        plt.ylim(-0.6, 1)
        plt.xlabel('Lag Time (steps)', fontsize=fontsize)
        plt.ylabel('Correlation', fontsize=fontsize)
        plt.gcf().get_axes()[0].tick_params(labelsize=fontsize)
        plt.tight_layout()

        if show:
            plt.show()

    def plot_trajectory(self, traj_no, show=False):
        """ Plot position versus time of FLM realization(s)

        :param traj_no: trajectory number or list of trajetory numbers to plot
        :param show: show the plot when done

        :type traj_no: int or list of int
        :type show: bool
        """

        if type(traj_no) is int:
            traj_no = [traj_no]

        plt.figure()
        for i in traj_no:
            plt.plot(self.realizations[i, :], lw=2)

        plt.xlabel('Time', fontsize=14)
        plt.ylabel('Position', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()

        if show:
            plt.show()

    def plot_msd(self, frac=0.4, nboot=200, confidence=68, show=False):
        """ Calculate and plot the mean squared displacement of the FLM realizations

        :param frac: fraction of MSD plot to show
        :param nboot: number of bootstrap trials
        :param confidence: percent confidence interval
        :param show: show the plot when done

        :type frac: float
        :type nboot: int
        :type confidence: float
        :type show: bool
        """

        plt.figure()

        msds = timeseries.msd(self.realizations.T[..., np.newaxis], 0)

        errorbars = timeseries.bootstrap_msd(msds, nboot, confidence=confidence)
        end = int(frac*self.N)
        avg = msds.mean(axis=1)

        plt.fill_between(np.arange(end), avg[:end] + errorbars[0, :end], avg[:end] - errorbars[1, :end], alpha=0.25)
        plt.plot(np.arange(end), avg[:end])
        plt.xlabel('Time', fontsize=14)
        plt.ylabel('MSD (nm$^2$)', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()

        if show:
            plt.show()


