"""
Unit and regression test for the flm package.
"""

# Import package, test suite, and other packages as needed
import flm
import unittest
import sys
import pytest
import numpy as np


def test_flm_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "flm" in sys.modules


# can write different sets of parameters named as different functions
def parameters1():

    params = dict()

    inputparams = dict()
    outputparams = dict()

    # __init__()
    inputparams['alpha'] = 1.5
    inputparams['hurst'] = 0.4
    inputparams['nsteps'] = 50
    inputparams['scale'] = 1
    inputparams['correct_hurst'] = False
    inputparams['truncate'] = None
    inputparams['correct_truncation'] = True
    inputparams['m'] = 2
    inputparams['M'] = 60
    inputparams['C'] = 1

    # generate_realizations()
    inputparams['ntraj'] = 2
    inputparams['progress'] = True  # use progress bar
    inputparams['truncated_distribution'] = None

    outputparams['M'] = 78
    outputparams['A'] = 0.3199853 + 0j
    outputparams['realizations'] = -18.3440877

    params['input_params'] = inputparams
    params['output_params'] = outputparams

    return params


def parameters2():

    params = dict()

    inputparams = dict()
    outputparams = dict()

    # __init__()
    inputparams['alpha'] = 1.5
    inputparams['hurst'] = 0.4
    inputparams['nsteps'] = 50
    inputparams['scale'] = 1
    inputparams['correct_hurst'] = False
    inputparams['truncate'] = None
    inputparams['correct_truncation'] = True
    inputparams['m'] = 2
    inputparams['M'] = 60
    inputparams['C'] = 1

    # generate_realizations()
    inputparams['ntraj'] = 2
    inputparams['progress'] = True  # use progress bar
    inputparams['truncated_distribution'] = None

    outputparams['M'] = 78
    outputparams['A'] = 0.3199853 + 0j

    params['input_params'] = inputparams
    params['output_params'] = outputparams
    outputparams['realizations'] = -18.3440877

    return params


test_params = [parameters1()]  # all of the sets of parameters to test. Just add to list


def initialize(inputparams):

    return flm.FLM(inputparams['nsteps'], inputparams['hurst'], inputparams['alpha'], scale=inputparams['scale'],
                   correct_hurst=inputparams['correct_hurst'], truncate=inputparams['truncate'],
                   correct_truncation=inputparams['correct_truncation'], m=inputparams['m'], M=inputparams['M'],
                   C=inputparams['C'])


@pytest.mark.parametrize("params", test_params)
def test_flm_init(params):

    FLM = initialize(params['input_params'])
    output_params = params['output_params']

    np.testing.assert_equal(output_params['M'], FLM.M)
    np.testing.assert_almost_equal(FLM.A[0], output_params['A'], 6)


@pytest.mark.parametrize("params", test_params)
def test_generate_realizations(params):

    np.random.seed(1)

    input_params = params['input_params']
    output_params = params['output_params']

    FLM = initialize(input_params)

    FLM.generate_realizations(input_params['ntraj'], progress=input_params['progress'],
                              truncated_distribution=input_params['truncated_distribution'])

    np.testing.assert_almost_equal(FLM.realizations[:, -1].sum(), output_params['realizations'], 6)


if __name__ == "__main__":

    unittest.main()
