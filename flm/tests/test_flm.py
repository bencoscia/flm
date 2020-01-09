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

    input_params = dict()
    output_params = dict()

    # __init__()
    input_params['alpha'] = 1.5
    input_params['hurst'] = 0.4
    input_params['n_steps'] = 50
    input_params['scale'] = 1
    input_params['correct_hurst'] = False
    input_params['truncate'] = None
    input_params['correct_truncation'] = True
    input_params['m'] = 2
    input_params['M'] = 60
    input_params['C'] = 1

    # generate_realizations()
    input_params['n_traj'] = 2
    input_params['progress'] = True  # use progress bar
    input_params['truncated_distribution'] = None

    output_params['M'] = 78
    output_params['A'] = 0.3199853 + 0j
    output_params['realizations'] = -18.3440877

    params['input_params'] = input_params
    params['output_params'] = output_params

    return params


def parameters2():

    params = dict()

    input_params = dict()
    output_params = dict()

    # __init__()
    input_params['alpha'] = 1.5
    input_params['hurst'] = 0.4
    input_params['n_steps'] = 50
    input_params['scale'] = 1
    input_params['correct_hurst'] = False
    input_params['truncate'] = None
    input_params['correct_truncation'] = True
    input_params['m'] = 2
    input_params['M'] = 60
    input_params['C'] = 1

    # generate_realizations()
    input_params['n_traj'] = 2
    input_params['progress'] = True  # use progress bar
    input_params['truncated_distribution'] = None

    output_params['M'] = 78
    output_params['A'] = 0.3199853 + 0j

    params['input_params'] = input_params
    params['output_params'] = output_params
    output_params['realizations'] = -18.3440877

    return params


test_params = [parameters1()]  # all of the sets of parameters to test. Just add to list


def initialize(input_params):

    return flm.FLM(input_params['n_steps'], input_params['hurst'], input_params['alpha'], scale=input_params['scale'],
                   correct_hurst=input_params['correct_hurst'], truncate=input_params['truncate'],
                   correct_truncation=input_params['correct_truncation'], m=input_params['m'], M=input_params['M'],
                   C=input_params['C'])


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

    FLM.generate_realizations(input_params['n_traj'], progress=input_params['progress'],
                              truncated_distribution=input_params['truncated_distribution'])

    np.testing.assert_almost_equal(FLM.realizations[:, -1].sum(), output_params['realizations'], 6)


if __name__ == "__main__":

    unittest.main()
