#!/usr/bin/env python

import numpy as np
import tqdm


def acf(t, largest_prime=500, autocov=False):
    """ Quickly calculated the autocorrelation function of a time series, t. This gives the same results as acf_slow()
    but uses FFTs. This method is faster than numpy.correlate.

    :param t: time series array (npoints, nseries)
    :param largest_prime : the largest prime factor of array length allowed. The smaller the faster. 1.6M points takes
    about 5 seconds with largest_prime=1000. Just be aware that you are losing data by truncating. But 5-6 data points
    isn't a big deal for large arrays.
    :param autocov: return autocovariance function instead (which is just the unnormalized autocorrelation)

    :type t: numpy.ndarray
    :type largest_prime: int
    :type autocov: bool

    """

    T = np.array(t)

    # Don't allow a prime factor larger than 'largest_prime'. Truncate data until that condition is met
    l = 2 * T.shape[0] - 1

    while largest_prime_factor(l) >= largest_prime or l % 2 == 0:
        l -= 1

    T = T[:(l + 1) // 2, ...]  # '...' allows for no second dimension if only a single time series is analysed
    length = T.shape[0] * 2 - 1

    T -= np.mean(T, axis=0)

    fftx = np.fft.fft(T, n=length, axis=0)
    ret = np.fft.ifft(fftx * np.conjugate(fftx), axis=0)
    ret = np.fft.fftshift(ret, axes=(0,))

    autocorr_fxn = ret[length // 2:].real

    if len(autocorr_fxn.shape) > 1:
        autocorr_fxn /= np.arange(T.shape[0], 0, -1)[:, None]
    else:
        autocorr_fxn /= np.arange(T.shape[0], 0, -1)

    if not autocov:
            autocorr_fxn /= np.var(T, axis=0)

    return autocorr_fxn  # normalized


def autocovariance(t, largest_prime=500):

    return acf(t, largest_prime=largest_prime, autocov=True)


def autocorrFFT(x):
    """ Function used for fast calculation of mean squared displacement

    :param x:
    :return:
    """

    N = len(x)
    F = np.fft.fft(x, n=2*N)  # 2*N because of zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res = (res[:N]).real   # now we have the autocorrelation in convention B
    n = N*np.ones(N) - np.arange(0, N)  # divide res(m) by (N-m)

    return res / n  # this is the autocorrelation in convention A


def bootstrap_msd(msds, nboot, confidence=68):
    """ Estimate error at each point in the MSD curve using bootstrapping

    :param msds: mean squared discplacements to sample
    :param N: number of bootstrap trials
    :param confidence: percentile for error calculation

    :type msds: np.ndarray
    :type N: int
    :type confidence: float
    """

    nT, nparticles = msds.shape

    msd_average = msds.mean(axis=1)

    eMSDs = np.zeros([nT, nboot], dtype=float)  # create n bootstrapped trajectories

    print('Bootstrapping MSD curves...')
    for b in tqdm.tqdm(range(nboot)):
        indices = np.random.randint(0, nparticles, nparticles)  # randomly choose particles with replacement
        for n in range(nparticles):
            eMSDs[:, b] += msds[:, indices[n]]  # add the MSDs of a randomly selected particle
        eMSDs[:, b] /= nparticles  # average the MSDs

    lower_confidence = (100 - confidence) / 2
    upper_confidence = 100 - lower_confidence

    limits = np.zeros([2, nT], dtype=float)  # upper and lower bounds at each point along MSD curve
    # determine error bound for each tau (out of n MSD's, use that for the error bars)
    for t in range(nT):
        limits[0, t] = np.abs(np.percentile(eMSDs[t, :], lower_confidence) - msd_average[t])
        limits[1, t] = np.abs(np.percentile(eMSDs[t, :], upper_confidence) - msd_average[t])

    return limits


def ensemble_msd(x0, x, size):

    if size == 1:

        return (x - x0) ** 2

    else:

        return np.linalg.norm(x0 - x, axis=1) ** 2


def largest_prime_factor(n):
    i = 2
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
    return n


def msd(x, axis, ensemble=False):
    """ Calculate mean square displacement based on particle positions

    :param x: particle positions
    :param axis: axis along which you want MSD (0, 1, 2, [0, 1], [0, 2], [1, 2], [0, 1, 2])
    :param ensemble: if True, calculate the ensemble MSD instead of the time-averaged MSD

    :type x: ndarray (n_frames, n_particles, 3)
    :type axis: int or list of ints
    :type ensemble: bool

    :return: MSD of each particle
    """

    frames = x.shape[0]  # number of trajectory frames
    ntraj = x.shape[1]  # number of trajectories
    MSD = np.zeros([frames, ntraj], dtype=float)  # a set of MSDs per particle

    size = len(x[0, :, axis].shape)  # number of axes in array where MSDs will be calculate

    if ensemble:

        for n in range(ntraj):  # start at 1 since all row 0 will be all zeros
            MSD[:, n] = ensemble_msd(x[0, n, axis], x[:, n, axis], size)

    else:

        for n in tqdm.tqdm(range(ntraj)):
            MSD[:, n] = msd_fft((x[:, n, :], axis))

    return MSD


def msd_fft(args):
    """ Calculate msd using a fast fourier transform algorithm

    :param x: trajectory of particle positions, equispaced in time
    :param axis: axis along which to calculate msd ({x:0, y:1, z:2})

    :type x: np.ndarray
    :type axis: int

    :return: msd as a function of time
    """

    x, axis = args

    r = np.copy(x)
    r = r[:, axis]

    if len(r.shape) == 1:
        r = r[:, np.newaxis]

    N = len(r)
    D = np.square(r).sum(axis=1)
    D = np.append(D, 0)
    S2 = sum([autocorrFFT(r[:, i]) for i in range(r.shape[1])])
    Q = 2 * D.sum()
    S1 = np.zeros(N)
    for m in range(N):
      Q = Q - D[m - 1] - D[N - m]
      S1[m] = Q / (N - m)

    return S1 - 2 * S2
