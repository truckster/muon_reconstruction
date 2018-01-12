import statusAlert, recoPreparation, color_schemes
from os import mkdir
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.colors as mco
from scipy.optimize import curve_fit
from scipy import optimize
from scipy.stats import gamma


def reco_by_fittig_gauss(pmt_position_class, snippet_class, muon_points, out_path):

    statusAlert.processStatus("search gauss patterns in snippet")

    '''photon data pre-processed and ordered into time snippets'''
    snippet_array = np.asarray(snippet_class.time_snippets)

    '''iterate snippets'''
    fit_results = []
    # for snippet in range(len(snippet_array)):
    for snippet in range(10):
        statusAlert.processStatus("processing snippet: " + str(snippet))

        # gauss_fit_data_vertical = combine_pmts_vertical(pmt_position_class, snippet_class, snippet)
        # mkdir(out_path + "test/vertical2/" + str(snippet))
        # fit_results.append(vertical_gauss_fit(gauss_fit_data_vertical, snippet, out_path + "test/vertical/" + str(snippet) + "/"))

        gauss_fit_data_horizontal = combine_pmts_horizontal(pmt_position_class, snippet_class, snippet)
        # mkdir(out_path + "test/horizontal2/" + str(snippet))
        fit_results.append(horizontal_gauss_fit(gauss_fit_data_horizontal, snippet, out_path + "test/horizontal/" + str(snippet) + "/"))

    print len(fit_results)


def vertical_gauss_fit(fit_data, snippet, out_path):
    '''get data into numpy arrays'''

    vertical_fit_results_per_snippet = []

    for vertical_sector in range(len(fit_data)):
        print "Fit vertical sector: " + str(vertical_sector)
        pmts_per_sector = fit_data[vertical_sector]
        # vertical_fit_results_per_snippet.append(fit_gauss_in_sector_vertical(pmts_per_sector, snippet, vertical_sector, out_path))
        find_all_gauss_in_data(pmts_per_sector, snippet, vertical_sector)


    return vertical_fit_results_per_snippet


def horizontal_gauss_fit(fit_data, snippet, out_path):
    '''get data into numpy arrays'''
    horizontal_fit_results_per_snippet = []
    for horizontal_sector in range(len(fit_data)):
        # print "Fit horizontal vector: " + str(horizontal_sector)
        pmts_per_sector = fit_data[horizontal_sector]

        horizontal_fit_results_per_snippet.append(fit_gauss_in_sector_horizontal(pmts_per_sector, snippet, horizontal_sector, out_path))
        find_all_gauss_in_data(pmts_per_sector, snippet, horizontal_sector)

    return horizontal_fit_results_per_snippet


class GaussFitResults:
    def __init__(self):
        self.height = 0
        self.location = 0
        self.width = 0


def find_all_gauss_in_data(grouped_sector_pmt_data, snippet, sector):
    histogramm_data = []
    for pmt in grouped_sector_pmt_data:
        for hits in range(pmt.hits):
            histogramm_data.append(pmt.theta)
    n, bins, patches = plt.hist(histogramm_data, bins=60, range=(0., math.pi))

    all_gauss_found = False
    peak = 0
    while not all_gauss_found:
        gauss_parameters = fitgaussian(n)

        if gauss_parameters[0] > 10 and gauss_parameters[1] > 0:
            print n
            single_gauss_fit_parameters = GaussFitResults()
            single_gauss_fit_parameters.height = gauss_parameters[0]
            single_gauss_fit_parameters.location = gauss_parameters[1]/n.size*math.pi
            single_gauss_fit_parameters.width = gauss_parameters[2]/n.size*math.pi

            n = n[int(gauss_parameters[1]+5):]

            peak += 1
            plt.savefig("/home/gpu/Analysis/muReconstruction/Output/test/horizontal2/" + str(snippet) + "/" + str(sector) + "+" + str(peak) + ".png")
            plt.close()
        else:
            all_gauss_found = True


def fit_gauss_in_sector_vertical(sector_pmts, snippet, vertical_sector, out_path):
    hist_list = []
    for pmt in sector_pmts:
        for hits in range(pmt.hits):
            hist_list.append(pmt.theta)
    n, bins, patches = plt.hist(hist_list, bins=60, range=(0., math.pi))

    p = fitgaussian(n)
    bin_centers = bins[:-1] + 0.5 * (bins[1:] - bins[:-1])
    plt.plot(bin_centers, gauss(bin_centers, p[0], p[1] / n.size * math.pi, p[2] / n.size * math.pi), 'r--')
    plt.figtext(0.65, 0.7, ("Height: " + str(p[0])
                            + "\nCenter: " + str(p[1] / n.size * math.pi)
                            + "\nWidth: " + str(p[2] / n.size * math.pi)
                            + "\n\nHeight/Width: " + str(p[0] / (p[2] / n.size * math.pi))))

    if 12000 > p[0]/(p[2] / n.size * math.pi) > 3000 and 0.03 > (p[2] / n.size * math.pi) > 0 and n.max() > 25:
        print("Muon entry or exit found in vertical sector: " + str(vertical_sector) + " at theta: " + str(p[1] / n.size * math.pi))

    plt.savefig(out_path + str(vertical_sector) + ".png")
    plt.close()


def fit_gauss_in_sector_horizontal(sector_pmts, snippet, horizontal_sector, out_path):
    hist_list = []
    for pmt in sector_pmts:
        for hits in range(pmt.hits):
            hist_list.append(pmt.phi)
    n, bins, patches = plt.hist(hist_list, bins=60, range=(0., math.pi))


    p = fitgaussian(n)
    bin_centers = bins[:-1] + 0.5 * (bins[1:] - bins[:-1])
    plt.xlim([-4, 4])
    plt.plot(bin_centers, gauss(bin_centers, p[0], p[1] / n.size * math.pi, p[2] / n.size * math.pi), 'r--')
    plt.figtext(0.65, 0.7, ("Height:  %f \nCenter:  %f \nWidth:  %f\n\nHeight/Width: %f"
                            % (p[0], p[1] / n.size * math.pi, p[2] / n.size * math.pi,
                               p[0] / (p[2] / n.size * math.pi))))

    if 12000 > p[0]/(p[2] / n.size * math.pi) > 3000 and 0.03 > (p[2] / n.size * math.pi) > 0 and n.max() > 25:
        # print("n.max = " + str(n.max()))
        # print p[0] / (p[2] / n.size * math.pi)
        print("Muon entry or exit found in horizontal sector: " + str(horizontal_sector) + " at phi: " + str(p[1] / n.size * math.pi))

    plt.savefig(out_path + str(horizontal_sector) + ".png")
    plt.close()


# Function to be fitted
def gauss(x, height, center, width):
    return height*np.exp(-(((center-x)/width)**2)/2)


def gaussian(height, center, width):
    width = float(width)
    return lambda x: height*np.exp(-(((center-x)/width)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    X = np.arange(data.size)
    x = np.sum(X * data) / np.sum(data)
    width = np.sqrt(np.abs(np.sum((X - x) ** 2 * data) / np.sum(data)))
    height = data.max()
    return height, x, width


def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) - data)
    p, success = optimize.leastsq(errorfunction, params)
    return p


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
