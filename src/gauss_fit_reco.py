import statusAlert, recoPreparation, color_schemes
from os import mkdir, path
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import leastsq


def fit_function_caller(pmt_position_class, snippet_class, muon_points, out_path):
    statusAlert.processStatus("search gauss patterns in snippet")

    '''photon data pre-processed and ordered into time snippets'''
    snippet_array = np.asarray(snippet_class.time_snippets)

    '''iterate snippets'''
    fit_results = []
    # for snippet in range(len(snippet_array)):
    for snippet in range(30):
        statusAlert.processStatus("processing snippet: " + str(snippet))
        gauss_fit_data_horizontal = combine_pmts_horizontal(pmt_position_class, snippet_class, snippet)
        if not path.isdir(out_path + "test/horizontal/" + str(snippet)):
            mkdir(out_path + "test/horizontal/" + str(snippet))
        gauss_fit(gauss_fit_data_horizontal, snippet, out_path + "test/horizontal/" + str(snippet) + "/")


def gauss_fit(fit_data, snippet, out_path):
    '''get data into numpy arrays'''
    horizontal_fit_results_per_snippet = []
    for horizontal_sector in range(len(fit_data)):
        print "Fit horizontal sector: " + str(horizontal_sector)
        pmts_per_sector = fit_data[horizontal_sector]

        horizontal_fit_results_per_snippet.append(fit_double_gauss_in_sector(pmts_per_sector,
                                                                             snippet,
                                                                             horizontal_sector,
                                                                             out_path))

    return horizontal_fit_results_per_snippet



def fit_double_gauss_in_sector(sector_pmts, snipp, sector, out_path):
    # draw 1 actual data
    hist_list = []
    for pmt in sector_pmts:
        for hits in range(pmt.hits):
            hist_list.append(pmt.phi)
    n, bins, patches = plt.hist(hist_list, bins=60, range=(0., math.pi))

    bin_centers = bins[:-1] + 0.5 * (bins[1:] - bins[:-1])
    plt.xlim([-4, 4])

    single_gauss_params = [c, mu, sigma] = [0, 1, 1]
    double_gauss_params = [c1, mu1, sigma1, c2, mu2, sigma2] = [0, 1, 1, 1, 1, 1] # Initial guesses for leastsq

    plsq = leastsq(res_double_gauss_dist, double_gauss_params, args=(n, bin_centers))
    # print plsq[0]
    print plsq[1]

    y_est = double_gaussian_dist(bin_centers, plsq[0])

    plt.plot(bin_centers, y_est, c='r')

    plt.savefig(out_path + str(sector) + ".png")
    plt.close()


def res_double_gauss_dist(p, y, x):
    (c1, mu1, sigma1, c2, dmu, sigma2) = p
    y_fit = double_gaussian_dist(x, p)
    err = y - y_fit
    return err


def res_single_gauss(p, x, y):
    (c, mu, sigma) = p
    y_fit = single_gaussian(x, p)
    err = y - y_fit
    return err


def single_gaussian(x, params):
    (c, mu, sigma) = params
    res = c * np.exp(-(x - mu) ** 2.0 / (2.0 * sigma ** 2.0))
    return res


def double_gaussian_dist(x, params):
    (c1, mu1, sigma1, c2, dmu, sigma2) = params
    res = c1 * np.exp(-(x - mu1) ** 2.0 / (2.0 * sigma1 ** 2.0)) + c2 * np.exp(
        -(x - mu1+dmu) ** 2.0 / (2.0 * sigma2 ** 2.0))
    return res


class pmt_data_for_fit_class:
    def __init__(self):  # this method creates the class object.
        self.x_sector = 0
        self.y_sector = 0
        self.hits = 0
        self.phi = 0
        self.theta = 0
        self.hist_array = []


def combine_pmts_vertical(pmt_position_class, snippet_class, snippet):
    return_array = [[] for _ in range(pmt_position_class.x_sectors)]
    for pmt_id in range(len(pmt_position_class.id)):
        pmt_data = pmt_data_for_fit_class()
        pmt_data.x_sector = pmt_position_class.is_x_sector[pmt_id]
        pmt_data.y_sector = pmt_position_class.is_y_sector[pmt_id]
        pmt_data.hits = snippet_class.time_snippets[snippet][pmt_id]
        pmt_data.phi = pmt_position_class.phi_position[pmt_id]
        pmt_data.theta = pmt_position_class.theta_position[pmt_id]

        return_array[pmt_position_class.is_x_sector[pmt_id]].append(pmt_data)
    return return_array


def combine_pmts_horizontal(pmt_position_class, snippet_class, snippet):
    return_array = [[] for _ in range(pmt_position_class.y_sectors)]
    for pmt_id in range(len(pmt_position_class.id)):
        pmt_data = pmt_data_for_fit_class()
        pmt_data.x_sector = pmt_position_class.is_x_sector[pmt_id]
        pmt_data.y_sector = pmt_position_class.is_y_sector[pmt_id]
        pmt_data.hits = snippet_class.time_snippets[snippet][pmt_id]
        pmt_data.phi = pmt_position_class.phi_position[pmt_id]
        pmt_data.theta = pmt_position_class.theta_position[pmt_id]

        return_array[pmt_position_class.is_y_sector[pmt_id]].append(pmt_data)
    return return_array
