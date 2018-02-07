import statusAlert, recoPreparation, color_schemes
from os import mkdir, path
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy.signal import argrelmax


def fit_function_caller(pmt_position_class, snippet_class, muon_points, out_path):
    statusAlert.processStatus("search gauss patterns in snippet")

    '''photon data pre-processed and ordered into time snippets'''
    snippet_array = np.asarray(snippet_class.time_snippets)

    '''iterate snippets'''
    fit_results = []
    # for snippet in range(len(snippet_array)):
    for snippet in range(9, 10):
        statusAlert.processStatus("processing snippet: " + str(snippet))

        statusAlert.processStatus("     combining pmt data to sectors")
        gauss_fit_data_horizontal = combine_pmts_horizontal(pmt_position_class, snippet_class, snippet)
        gauss_fit_data_vertical = combine_pmts_vertical(pmt_position_class, snippet_class, snippet)

        if not path.isdir(out_path + "test/horizontal/" + str(snippet)):
            mkdir(out_path + "test/horizontal/" + str(snippet))
        if not path.isdir(out_path + "test/vertical/" + str(snippet)):
            mkdir(out_path + "test/vertical/" + str(snippet))

        statusAlert.processStatus("     performing fit")

        statusAlert.processStatus("         horizontal")
        # gauss_fit(gauss_fit_data_horizontal, snippet, out_path + "test/horizontal/" + str(snippet) + "/")
        multifit(gauss_fit_data_horizontal, snippet, out_path + "test/horizontal/" + str(snippet) + "/")

        statusAlert.processStatus("         vertical")
        # gauss_fit(gauss_fit_data_vertical, snippet, out_path + "test/vertical/" + str(snippet) + "/")
        multifit(gauss_fit_data_vertical, snippet, out_path + "test/vertical/" + str(snippet) + "/")


class SectorHistoData:
    def __init__(self):  # this method creates the class object.
        self.entries = []
        self.bins = []
        self.patches = []


def get_sector_histo(sector_pmts):
    return_data = SectorHistoData
    hist_list = []
    for pmt in sector_pmts:
        for hits in range(pmt.hits):
            hist_list.append(pmt.phi)
    n, bins, patches = plt.hist(hist_list, bins=60, range=(0, math.pi))

    return_data.entries = n
    return_data.bins = bins
    return_data.patches = patches

    plt.clf()

    return return_data


def gauss_fit(fit_data, snippet, out_path):
    '''get data into numpy arrays'''
    for sector in range(len(fit_data)):
        print "Fit sector: " + str(sector)
        pmts_per_sector = fit_data[sector]
        sector_data = get_sector_histo(pmts_per_sector)

        single_fit = fit_single_gauss_in_sector(sector_data)
        double_fit = fit_double_gauss_in_sector(sector_data)

        picture_drawer(sector_data, sector, single_fit, double_fit, out_path)


def multifit(fit_data, snippet, out_path):
    '''get data into numpy arrays'''
    for sector in range(len(fit_data)):
        print "Fit sector: " + str(sector)
        pmts_per_sector = fit_data[sector]
        sector_pmts = get_sector_histo(pmts_per_sector)

        bin_centers = 0.5 * (sector_pmts.bins[1:] - sector_pmts.bins[:-1])
        plt.xlim([0, math.pi])

        possible_gauss_positions = get_fit_estimates(sector_pmts.entries)

        y_values = plt.errorbar(sector_pmts.bins[:-1], sector_pmts.entries, yerr=np.sqrt(sector_pmts.entries))
        plt.clf()

        print(possible_gauss_positions)

        fit_results = []

        for peak in range(len(possible_gauss_positions)):
            print(bin_centers[[possible_gauss_positions[peak]]])
            print(sector_pmts.entries[[possible_gauss_positions[peak]]])

            # single_gauss_params = [c, mu, sigma] = [sector_pmts.entries[possible_gauss_positions[peak]],
            #                                         possible_gauss_positions[peak]/float(len(sector_pmts.bins))*3.1415,
            #                                         0.05]
            # plsq = leastsq(res_single_gauss, single_gauss_params, args=(bin_centers, sector_pmts.entries))
            #
            # print(single_gauss_params)
            # print(plsq)
            # fit_results.append(plsq)

            fit_parameters = (sector_pmts.entries[possible_gauss_positions[peak]], (possible_gauss_positions[peak]/60.0*3.1415), 0.1)
            try:
                fit, err = curve_fit(gauss,
                                     bin_centers,
                                     sector_pmts.entries,
                                     # bin_centers[sector_pmts.entries[possible_gauss_positions[peak]]-5:sector_pmts.entries[possible_gauss_positions[peak]]+5],
                                     # sector_pmts.entries[sector_pmts.entries[possible_gauss_positions[peak]]-5:sector_pmts.entries[possible_gauss_positions[peak]]+5],
                                     p0=fit_parameters,
                                     # bounds=(
                                     #     [0, (possible_gauss_positions[peak]/60.0*3.1415), 0],
                                     #     [sector_pmts.entries[possible_gauss_positions[peak]], (possible_gauss_positions[peak]/60.0*3.1415)+0.2, 0.2])
                                     )

                # print(1, possible_gauss_positions[peak] / 60.0 * 3.1415, 1)
                print(fit)
                print(fit_parameters)
            except:
                fit = [1, 1, 1]
                print("Fit no work")
            fit_results.append(fit)

        picture_drawer_2(sector_pmts, sector, fit_results, out_path)


def fit_single_gauss_in_sector(sector_pmts):
    bin_centers = sector_pmts.bins[:-1] + 0.5 * (sector_pmts.bins[1:] - sector_pmts.bins[:-1])
    plt.xlim([-4, 4])

    single_gauss_params = [c, mu, sigma] = [1, -3, 1]
    plsq = leastsq(res_single_gauss, single_gauss_params, args=(sector_pmts.entries, bin_centers))

    # print("Fit parameters: " + str(plsq[0]))
    # print("Fit quality parameter: " + str(plsq[1]))

    plt.clf()

    return plsq


def fit_double_gauss_in_sector(sector_pmts):

    bin_centers = sector_pmts.bins[:-1] + 0.5 * (sector_pmts.bins[1:] - sector_pmts.bins[:-1])
    plt.xlim([-4, 4])

    double_gauss_params = [c1, mu1, sigma1, c2, dmu, sigma2] = [1, 0, 1, 1, 2, 1] # Initial guesses for leastsq

    plsq = leastsq(res_double_gauss_dist, double_gauss_params, args=(sector_pmts.entries, bin_centers))

    # print("Fit parameters: " + str(plsq[0]))
    # print("Fit quality parameter: " + str(plsq[1])

    plt.clf()

    return plsq


def get_fit_estimates(hist_entries):
    max_array = argrelmax(hist_entries)

    if len(max_array[0]):
        return max_array[0]
    else:
        return [0]


def picture_drawer(sector_pmts, sector, single_fit, double_fit, out_path):
    # draw 1 actual data
    plt.hist(sector_pmts.bins[:-1], len(sector_pmts.bins)-1, weights=sector_pmts.entries)

    bin_centers = sector_pmts.bins[:-1] + 0.5 * (sector_pmts.bins[1:] - sector_pmts.bins[:-1])
    # plt.xlim([-4, 4])

    single_gauss_fit_graph = single_gaussian(bin_centers, single_fit[0])
    double_gauss_fit_graph = double_gaussian_dist(bin_centers, double_fit[0])

    # plt.plot(bin_centers, single_gauss_fit_graph, c='k')
    plt.plot(sector_pmts.bins[:-1], double_gauss_fit_graph, c='r')

    plt.ylabel("theta (deg)")
    plt.xlabel("phi (deg)")

    plt.savefig(out_path + str(sector) + ".png")
    plt.close()


def picture_drawer_2(sector_pmts, sector, single_fit, out_path):
    # draw 1 actual data
    plt.xlim(0, math.pi)
    bin_centers = 0.5 * (sector_pmts.bins[1:] + sector_pmts.bins[:-1])
    # plt.hist(sector_pmts.bins, len(sector_pmts.bins)-1, weights=sector_pmts.entries, color='b')
    plt.bar(bin_centers, sector_pmts.entries, width=1.0*math.pi/(len(sector_pmts.bins)-1))

    plt.errorbar(bin_centers, sector_pmts.entries, yerr=np.sqrt(sector_pmts.entries),
                 fmt='b', linestyle=''
                 )
    for param in range(len(single_fit)):
        single_gauss_fit_graph = gauss(bin_centers,
                                       single_fit[param][0],
                                       single_fit[param][1],
                                       single_fit[param][2])
        plt.plot(bin_centers, single_gauss_fit_graph, 'r--')

        plt.ylabel("Number of Photons")
        plt.xlabel("theta (deg)")

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
    res = (c/(math.sqrt(2.0 * math.pi * sigma ** 2.0))) * np.exp(-(x - mu) ** 2.0 / (2.0 * sigma ** 2.0))
    return res


def gauss(x, c, mu, sigma):
    # res = (c / (math.sqrt(2.0 * math.pi * sigma ** 2.0))) * np.exp(-(x - mu) ** 2.0 / (2.0 * sigma ** 2.0))
    res = c * np.exp(-(x - mu) ** 2.0 / (sigma ** 2.0))
    return res


def lorentz(x, A, x0, gamma):
    return A * gamma**2 / ((x - x0) ** 2 + gamma ** 2)


def double_gaussian_dist(x, params):
    (c1, mu1, sigma1, c2, dmu, sigma2) = params
    res = (c1/(math.sqrt(2.0 * math.pi * sigma1 ** 2.0))) * np.exp(-(x - mu1) ** 2.0 / (2.0 * sigma1 ** 2.0))\
          + (c2/(math.sqrt(2.0 * math.pi * sigma2 ** 2.0))) * np.exp(-(x - (mu1+dmu)) ** 2.0 / (2.0 * sigma2 ** 2.0))
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
