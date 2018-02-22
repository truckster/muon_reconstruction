import statusAlert, recoPreparation, color_schemes
from os import mkdir, path
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy.signal import argrelmax


def fit_function_caller(pmt_position_class, snippet_class, muon_points, out_path, result_file):
    statusAlert.processStatus("search gauss patterns in snippet")

    '''photon data pre-processed and ordered into time snippets'''
    snippet_array = np.asarray(snippet_class.time_snippets)

    '''iterate snippets'''
    fit_results = []
    # for snippet in range(len(snippet_array)):
    for snippet in range(35):
        statusAlert.processStatus("processing snippet: " + str(snippet))

        statusAlert.processStatus("     combining pmt data to sectors")
        gauss_fit_data_horizontal = combine_pmts_horizontal(pmt_position_class, snippet_class, snippet)
        gauss_fit_data_vertical = combine_pmts_vertical(pmt_position_class, snippet_class, snippet)

        if not path.isdir(out_path + "horizontal/"):
            mkdir(out_path + "horizontal/")
        if not path.isdir(out_path + "vertical/"):
            mkdir(out_path + "vertical/")

        if not path.isdir(out_path + "horizontal/" + str(snippet)):
            mkdir(out_path + "horizontal/" + str(snippet))
        if not path.isdir(out_path + "vertical/" + str(snippet)):
            mkdir(out_path + "vertical/" + str(snippet))

        statusAlert.processStatus("     performing fit")

        statusAlert.processStatus("         horizontal")
        fit, x_range_lower = \
            multifit(gauss_fit_data_horizontal, -math.pi, snippet, out_path + "horizontal/" + str(snippet) + "/", result_file)

        statusAlert.processStatus("         vertical")
        fit, x_range_lower = \
            multifit(gauss_fit_data_vertical, 0, snippet, out_path + "vertical/" + str(snippet) + "/", result_file)


class SectorHistoData:
    def __init__(self):  # this method creates the class object.
        self.entries = []
        self.bins = []
        self.patches = []


def get_sector_histo(sector_pmts, x_lower_limit):
    return_data = SectorHistoData
    hist_list = []
    for pmt in sector_pmts:
        for hits in range(pmt.hits):
            hist_list.append(pmt.phi)
    n, bins, patches = plt.hist(hist_list, bins=60, range=(x_lower_limit, math.pi))

    return_data.entries = n
    return_data.bins = bins
    return_data.patches = patches

    plt.clf()

    return return_data


def multifit(fit_data, x_range_lower, snippet, out_path, result_file):
    '''get data into numpy arrays'''
    for sector in range(len(fit_data)):
        print("----------------------------------------------------------------------------------------")
        print "Fit sector: " + str(sector)
        pmts_per_sector = fit_data[sector]
        sector_pmts = get_sector_histo(pmts_per_sector, x_range_lower)

        bin_centers = 0.5 * (sector_pmts.bins[1:] + sector_pmts.bins[:-1])
        plt.xlim([x_range_lower, math.pi])

        possible_gauss_positions = get_fit_estimates(sector_pmts.entries)

        y_values = plt.errorbar(sector_pmts.bins[:-1], sector_pmts.entries, yerr=np.sqrt(sector_pmts.entries))
        plt.clf()

        print(possible_gauss_positions)

        fit_results = []

        for peak in range(len(possible_gauss_positions)):
            bin_range = 3
            interesting_bins = bin_centers[possible_gauss_positions[peak]-bin_range:
                                           possible_gauss_positions[peak]+bin_range]
            interesting_entries = sector_pmts.entries[possible_gauss_positions[peak]-bin_range:
                                                      possible_gauss_positions[peak]+bin_range]

            # single_gauss_params = [c, mu, sigma] = [sector_pmts.entries[possible_gauss_positions[peak]],
            #                                         possible_gauss_positions[peak]/float(len(sector_pmts.bins))*3.1415,
            #                                         0.05]
            # plsq = leastsq(res_single_gauss, single_gauss_params, args=(bin_centers, sector_pmts.entries))
            #
            # print(single_gauss_params)
            # print(plsq)
            # fit_results.append(plsq)

            fit_parameters = (sector_pmts.entries[possible_gauss_positions[peak]],
                              (possible_gauss_positions[peak]/60.0*3.1415),
                              0.1)

            def gauss_in(x, c, mu, sigma):
                res = c * np.exp(-(x - mu) ** 2.0 / (sigma ** 2.0))
                return res

            try:
                fit, err = curve_fit(gauss_in,
                                     interesting_bins,
                                     interesting_entries,
                                     # bin_centers[sector_pmts.entries[possible_gauss_positions[peak]]-5:sector_pmts.entries[possible_gauss_positions[peak]]+5],
                                     # sector_pmts.entries[sector_pmts.entries[possible_gauss_positions[peak]]-5:sector_pmts.entries[possible_gauss_positions[peak]]+5],
                                     # p0=fit_parameters,
                                     # bounds=(
                                     #     [sector_pmts.entries[possible_gauss_positions[peak]]-50, 0, 0],
                                     #     [sector_pmts.entries[possible_gauss_positions[peak]]+50, np.inf, 5])
                                     )

                # print(1, possible_gauss_positions[peak] / 60.0 * 3.1415, 1)
                print(fit[0]/fit[2])
                print(peak/60.0*math.pi)
            except:
                fit = [1, 1, 1]
                # fit = [1, 1, 1, 1]
                print("Fit no work")
            fit_results.append(fit)
            if abs(fit[0] / fit[2] / fit[2]) > 1000000 and abs(fit[2]) < 0.2:
                result_output_writer(result_file, fit, x_range_lower, snippet)

        picture_drawer_2(sector_pmts, sector, fit_results, x_range_lower, out_path)

    return fit, x_range_lower


def get_fit_estimates(hist_entries):
    max_array = argrelmax(hist_entries)

    if len(max_array[0]):
        return max_array[0]
    else:
        return [0]


def picture_drawer_2(sector_pmts, sector, single_fit, x_range_lower, out_path):
    # draw 1 actual data
    plt.xlim(x_range_lower, math.pi)
    bin_centers = 0.5 * (sector_pmts.bins[1:] + sector_pmts.bins[:-1])
    # plt.hist(sector_pmts.bins, len(sector_pmts.bins)-1, weights=sector_pmts.entries, color='b')
    plt.bar(bin_centers, sector_pmts.entries, width=(math.pi-x_range_lower)/math.pi*math.pi/(len(sector_pmts.bins)-1))

    plt.errorbar(bin_centers, sector_pmts.entries, yerr=np.sqrt(sector_pmts.entries),
                 fmt='b', linestyle=''
                 )
    for param in range(len(single_fit)):
        single_gauss_fit_graph = gauss(bin_centers,
                                       single_fit[param][0],
                                       single_fit[param][1],
                                       single_fit[param][2])

        # single_gauss_fit_graph = gauss_d(bin_centers,
        #                                single_fit[param][0],
        #                                single_fit[param][1],
        #                                single_fit[param][2],
        #                                single_fit[param][3])

        plt.plot(bin_centers, single_gauss_fit_graph, 'r--')

        plt.ylabel("Number of Photons")
        if x_range_lower < 0:
            plt.xlabel("theta (deg)")
        else:
            plt.xlabel("phi (deg)")
        # plt.figtext(0.15, 0.75-0.12*param, ("Height/Width: %.1f \nWidth: %.5f \nPosition: %.3f"
        #                         % (single_fit[param][0]/single_fit[param][2],
        #                            single_fit[param][2],
        #                            single_fit[param][1])), fontsize=10)
        plt.figtext(0.15, 0.75 - 0.12 * param, ("Height/Width/Width: %.1f"
                                                % (single_fit[param][0] / single_fit[param][2]/single_fit[param][2])),
                    fontsize=10)

    plt.savefig(out_path + str(sector) + ".png")
    plt.close()


def gauss(x, c, mu, sigma):
    # res = (c / (math.sqrt(2.0 * math.pi * sigma ** 2.0))) * np.exp(-(x - mu) ** 2.0 / (2.0 * sigma ** 2.0))
    res = c * np.exp(-(x - mu) ** 2.0 / (sigma ** 2.0))
    return res


def gauss_d(x, c, mu, sigma, d):
    # res = (c / (math.sqrt(2.0 * math.pi * sigma ** 2.0))) * np.exp(-(x - mu) ** 2.0 / (2.0 * sigma ** 2.0))
    res = d + c * np.exp(-(x - mu) ** 2.0 / (sigma ** 2.0))
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


def calc_muon_points_sphere(muon_points_raw):
    muon_points = []
    for i in range(len(muon_points_raw)):
        for j in range(len(muon_points_raw[i])):
            muon_points.append("Phi")
            muon_points.append(recoPreparation.calcPMTPolarPhi(muon_points_raw[i][j]))
            muon_points.append("Theta")
            muon_points.append(recoPreparation.calcPMTPolarTheta(muon_points_raw[i][j]))
    return muon_points


def result_output_writer(write_file, reconstructed_muon, orientation, snippet):
    write_file.write("Snippet: " + str(snippet) + '\n')
    if orientation < 0:
        write_file.write("Phi: " + str(reconstructed_muon[1]) + '\n')
    if orientation is 0:
        write_file.write("Theta: " + str(reconstructed_muon[1]) + '\n')
