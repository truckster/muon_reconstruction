import statusAlert, recoPreparation, color_schemes
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os


class PatternPosition:
    def __init__(self):  # this method creates the class object.
        self.hits = 0
        self.phi = 0
        self.theta = 0


def pattern_detector(pmt_position_class, snippet_class, muon_points, out_path):
    statusAlert.processStatus("Iterate snippets and draw")

    '''photon data pre-processed and ordered into time snippets'''
    snippet_array = np.asarray(snippet_class.time_snippets)

    '''iterate snippets'''
    # for snippet in range(len(snippet_array)):
    for snippet in range(50):
        statusAlert.processStatus("processing snippet: " + str(snippet))

        '''Draw the detector picture for the certain time snippet'''
        draw_snippet_picture(pmt_position_class, snippet_class.time_snippets[snippet], muon_points, snippet,
                             out_path, "absolute")


def pattern_detector_difference(pmt_position_class, snippet_class, muon_points, out_path):
    statusAlert.processStatus("Iterate snippets and draw")

    '''photon data pre-processed and ordered into time snippets'''
    snippet_array = np.asarray(snippet_class.time_snippets)

    '''iterate snippets'''
    # for snippet in range(len(snippet_array)):
    for snippet in range(50):
        if snippet > 0:
            snippet_diff = np.asarray(snippet_class.time_snippets[snippet]) \
                           - np.asarray(snippet_class.time_snippets[snippet-1])
            statusAlert.processStatus("processing snippet: " + str(snippet))

            '''Draw the detector picture for the certain time snippet'''
            draw_snippet_picture(pmt_position_class, snippet_diff, muon_points, snippet, out_path, "differential")


def draw_snippet_picture(pmt_position_class, snippet_array, muon_points, snippet, out_path, mode=None):
    statusAlert.processStatus("Creating graphics")
    fig = plt.figure(num=None, figsize=(20, 10))

    '''Analysis design'''
    ax1 = fig.add_subplot(111, axisbg='gainsboro')
    # scatterHits = ax1.scatter(pmt_position_class.phi_position, pmt_position_class.theta_position, marker='o', c='k')
    scatterHits = ax1.scatter(pmt_position_class.phi_position, pmt_position_class.theta_position,
                              c=snippet_array, cmap=color_schemes.analysis_design(),
                              edgecolor='k', label='hit')
    cb = fig.colorbar(scatterHits)

    '''Nice design'''
    # ax1 = fig.add_subplot(111)
    # scatter_hits_new = ax1.scatter(pmt_position_class.phi_hammer_aitov, pmt_position_class.theta_hammer_aitov,
    #                                c=snippet_class.time_snippets[snippet], cmap=color_analysis, label='hit')
    # cb = fig.colorbar(scatter_hits_new)

    draw_sector_lines(pmt_position_class)
    draw_muon_points(muon_points)
    # draw_pattern_results(sector_pattern_array, template_radius, ax1)
    # draw_reco_muons(muon_finder_psnippet)

    plt.ylabel("theta (deg)")
    plt.xlabel("phi (deg)")

    if mode is "absolute":
        cb.set_label("Number of photons")

        # plt.savefig(out_path + str(snippet) + ".pdf")
        try:
            os.chdir(out_path + "absolute hits/")
        except:
            os.makedirs(out_path + "absolute hits/")
        plt.savefig(out_path + "absolute hits/" + str(snippet) + ".png", bbox_inches='tight')

        plt.close()

    elif mode is "differential":
        cb.set_label("Number of photon difference")

        # plt.savefig(out_path + str(snippet) + ".pdf")
        try:
            os.chdir(out_path + "differential hits/")
        except:
            os.makedirs(out_path + "differential hits/")
        plt.savefig(out_path + "differential hits/" + str(snippet) + ".png", bbox_inches='tight')

        plt.close()

    else:
        print("Nothing to print")


def draw_muon_points(muon_points):
    phi_muon_entry = []
    theta_muon_entry = []
    phi_muon_exit = []
    theta_muon_exit = []

    for muon_event in muon_points:
        if muon_event.enters is True:
            phi_muon_entry.append(muon_event.phi)
            theta_muon_entry.append(muon_event.theta)
        if muon_event.leaves is True:
            phi_muon_exit.append(muon_event.phi)
            theta_muon_exit.append(muon_event.theta)

    scatter_muon_entry = plt.scatter(phi_muon_entry, theta_muon_entry,
                                     facecolors='none', edgecolors='white', marker='v', s=120)
    scatter_muon_exit = plt.scatter(phi_muon_exit, theta_muon_exit,
                                   facecolors='none', edgecolors='white', marker='^', s=120)



def draw_sector_lines(pmt_position_class):
    for i in range(pmt_position_class.x_sectors + 1):
        plt.plot([-math.pi + i * (2 * math.pi) / pmt_position_class.x_sectors,
                 -math.pi + i * (2 * math.pi) / pmt_position_class.x_sectors],
                 [0, math.pi],
                 'k--')
        if i % 5 is 0:
            plt.text(-math.pi + i * (2 * math.pi) / pmt_position_class.x_sectors,
                     -0.1,
                     i)
            plt.plot([-math.pi + i * (2 * math.pi) / pmt_position_class.x_sectors,
                      -math.pi + i * (2 * math.pi) / pmt_position_class.x_sectors],
                     [0, math.pi],
                     'k-',
                     linewidth=1.5)

    for i in range(pmt_position_class.y_sectors + 1):
        plt.plot([-math.pi, math.pi],
                 [i * math.pi / pmt_position_class.y_sectors, i * math.pi / pmt_position_class.y_sectors],
                 'k--')

        if i % 5 is 0:
            plt.text(-math.pi - 0.15,
                     i * math.pi / pmt_position_class.y_sectors,
                     i)
            plt.plot([-math.pi, math.pi],
                     [i * math.pi / pmt_position_class.y_sectors, i * math.pi / pmt_position_class.y_sectors],
                     'k-',
                     linewidth=1.5)


def draw_pattern_results(sector_pattern_array, template_radius, sub_plot):
    phi = []
    theta = []

    for i in range(len(sector_pattern_array)):
        phi.append(sector_pattern_array[i].phi)
        theta.append(sector_pattern_array[i].theta)
        scatter_circle = plt.Circle((phi[i], theta[i]), 0.1, color='white', fill=False)
        sub_plot.add_artist(scatter_circle)


def draw_reco_muons(pattern_class_array):
    phi = []
    theta = []

    for i in range(len(pattern_class_array)):
        phi.append(pattern_class_array[i].phi)
        theta.append(pattern_class_array[i].theta)
    plt.scatter(phi, theta, facecolors='none', edgecolors='white', marker='s', s=220)


def print_sector_pmts(pmt_position_class, out_path):
    for sector in range(len(pmt_position_class.sector_list_id)):

        id_np = np.asarray(pmt_position_class.sector_list_id[sector])
        phi_np = np.asarray(pmt_position_class.sector_list_phi[sector])
        theta_np = np.asarray(pmt_position_class.sector_list_theta[sector])

        fig = plt.figure(num=None, figsize=(20, 10))
        ax1 = fig.add_subplot(111)
        scatterHits = ax1.scatter(phi_np, theta_np)

        plt.savefig(out_path + "sectors/" + str(sector) + ".png")