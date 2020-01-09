"""
Unit and regression test for the flm package.
"""

# Import package, test suite, and other packages as needed
import flm
import unittest
import sys
import pytest
import numpy as np
import matplotlib.pyplot as plt


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

    output_params['M'] = 78
    output_params['A'] = 0.3199853 + 0j

    # generate_realizations()
    input_params['n_traj'] = 2
    input_params['progress'] = True  # use progress bar
    input_params['truncated_distribution'] = None

    output_params['realizations'] = -18.3440877

    # plot_marginal()
    input_params['marginal_bounds'] = (-5, 5)
    input_params['marginal_bins'] = 50

    output_params['marginal_x'] = 2502.5025025
    output_params['marginal_y'] = 95.77742125

    # plot_autocorrelation()
    input_params['k_acf'] = 3
    input_params['nboot_acf'] = 2
    input_params['confidence_acf'] = 68.27

    output_params['acf_x'] = 1225.0
    output_params['acf_y'] = 0.356708

    # plot_trajectory()
    input_params['traj_no'] = [0, [0, 1]]

    output_params['traj_x'] = [[1225.0], [1225.0, 1225.0]]
    output_params['traj_y'] = [[-833.20001517], [-833.20001517, 31.533102]]

    # plot_msd()
    input_params['msd_frac'] = 0.4
    input_params['msd_nboot'] = 2
    input_params['msd_confidence'] = 68

    output_params['msd_x'] = 190.0
    output_params['msd_y'] = 1494.172745

    params['input_params'] = input_params
    params['output_params'] = output_params

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

    return FLM


@pytest.mark.parametrize("params", test_params)
def test_plot_marginal(params):

    input_params = params['input_params']
    output_params = params['output_params']

    FLM = test_generate_realizations(params)

    fig, ax = plt.subplots()
    FLM.plot_marginal(bounds=input_params['marginal_bounds'], bins=input_params['marginal_bins'])
    x_plot, y_plot = ax.lines[0].get_xydata().T

    np.testing.assert_almost_equal(np.abs(x_plot).sum(), output_params['marginal_x'], 6)
    np.testing.assert_almost_equal(y_plot.sum(), output_params['marginal_y'], 6)


@pytest.mark.parametrize("params", test_params)
def test_plot_autocorrelation(params):

    input_params = params['input_params']
    output_params = params['output_params']

    FLM = test_generate_realizations(params)

    fig, ax = plt.subplots()
    # This also tests the autocorrelation function calculation
    FLM.plot_autocorrelation(max_k=input_params['k_acf'], nboot=input_params['nboot_acf'],
                             confidence=input_params['confidence_acf'], overlay=True)

    x_plot, y_plot = ax.lines[0].get_xydata().T

    np.testing.assert_almost_equal(x_plot.sum(), output_params['acf_x'], 6)
    np.testing.assert_almost_equal(y_plot.sum(), output_params['acf_y'], 6)


@pytest.mark.parametrize("params", test_params)
def test_plot_trajectory(params):

    input_params = params['input_params']
    output_params = params['output_params']

    FLM = test_generate_realizations(params)

    for i, t in enumerate(input_params['traj_no']):  # test integer and list inputs

        fig, ax = plt.subplots()

        FLM.plot_trajectory(t, overlay=True)

        for j, y in enumerate(ax.lines):

            x_plot, y_plot = y.get_xydata().T

            np.testing.assert_almost_equal(x_plot.sum(), output_params['traj_x'][i][j], 6)
            np.testing.assert_almost_equal(y_plot.sum(), output_params['traj_y'][i][j], 6)


@pytest.mark.parametrize("params", test_params)
def test_plot_msd(params):

    input_params = params['input_params']
    output_params = params['output_params']

    FLM = test_generate_realizations(params)

    fig, ax = plt.subplots()
    # This also tests the autocorrelation function calculation
    FLM.plot_msd(frac=input_params['msd_frac'], nboot=input_params['msd_nboot'],
                 confidence=input_params['msd_confidence'], overlay=True)

    x_plot, y_plot = ax.lines[0].get_xydata().T

    np.testing.assert_almost_equal(x_plot.sum(), output_params['msd_x'], 6)
    np.testing.assert_almost_equal(y_plot.sum(), output_params['msd_y'], 6)


if __name__ == "__main__":

    unittest.main()
