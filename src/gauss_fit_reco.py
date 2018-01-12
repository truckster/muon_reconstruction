import statusAlert, recoPreparation, color_schemes
from os import mkdir
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import optimize


def fit_function_caller(pmt_position_class, snippet_class, muon_points, out_path):
    statusAlert.processStatus("search gauss patterns in snippet")

    '''photon data pre-processed and ordered into time snippets'''
    snippet_array = np.asarray(snippet_class.time_snippets)

    '''iterate snippets'''
    fit_results = []
    # for snippet in range(len(snippet_array)):
    for snippet in range(10):
        statusAlert.processStatus("processing snippet: " + str(snippet))
        gauss_fit_data_horizontal = combine_pmts_horizontal(pmt_position_class, snippet_class, snippet)
        horizontal_gauss_fit(gauss_fit_data_horizontal, snippet, out_path + "test/horizontal/" + str(snippet) + "/")


def horizontal_gauss_fit(fit_data, snippet, out_path):
    '''get data into numpy arrays'''
    horizontal_fit_results_per_snippet = []
    for horizontal_sector in range(len(fit_data)):
        print "Fit horizontal sector: " + str(horizontal_sector)
        pmts_per_sector = fit_data[horizontal_sector]

        horizontal_fit_results_per_snippet.append(fit_gauss_in_sector_horizontal(pmts_per_sector,
                                                                                 snippet,
                                                                                 horizontal_sector,
                                                                                 out_path))

    return horizontal_fit_results_per_snippet



def fit_gauss_in_sector_horizontal(sector_pmts, snipp, sector, out_path):

    # draw actual data
    hist_list = []
    for pmt in sector_pmts:
        for hits in range(pmt.hits):
            hist_list.append(pmt.theta)
    n, bins, patches = plt.hist(hist_list, bins=60, range=(0., math.pi))

    # get estimate gauss values for data: NO FIT!!
    moms = moments(n)

    # draw gauss from moments
    bin_centers = bins[:-1] + 0.5 * (bins[1:] - bins[:-1])
    plt.xlim([-4, 4])
    # plt.plot(bin_centers, gauss(bin_centers, moms[0], moms[1] / n.size * math.pi, moms[2] / n.size * math.pi), 'r--')

    # print gauss parameters to image
    # plt.figtext(0.65, 0.7, ("Height: " + str(moms[0])
    #                         + "\nCenter: " + str(moms[1] / n.size * math.pi)
    #                         + "\nWidth: " + str(moms[2] / n.size * math.pi)
    #                         + "\n\nHeight/Width: " + str(moms[0] / (moms[2] / n.size * math.pi))))

    plt.savefig(out_path + str(sector) + ".png")
    plt.close()


# Function to be fitted
def gauss(x, height, center, width):
    return height * np.exp(-(((center - x) / width) ** 2) / 2)


def gaussian(height, center, width):
    width = float(width)
    return lambda x: height * np.exp(-(((center - x) / width) ** 2) / 2)


def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    X = np.arange(data.size)
    x = np.sum(X * data) / np.sum(data)
    width = np.sqrt(np.abs(np.sum((X - x) ** 2 * data) / np.sum(data)))
    height = data.max()

    # print(x/data.size * math.pi)
    # print(width/data.size * math.pi)
    # print(height)

    return height, x, width


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
