import statusAlert, recoPreparation, color_schemes, contour_analyze, reco_from_contour, PointVecDist
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.ndimage.filters import gaussian_filter
import os


def snippet_drawer(pmt_position_class, snippet_class, muon_points, out_path, number_contour_level,
                   reco_points=None):
    statusAlert.processStatus("Iterate snippets and draw")

    '''photon data pre-processed and ordered into time snippets'''
    snippet_array = np.asarray(snippet_class.time_snippets)

    '''iterate snippets'''
    if reco_points is None:
        for snippet in range(len(snippet_array)):
            statusAlert.processStatus("processing snippet: " + str(snippet))

            '''Draw the detector picture for the certain time snippet'''
            draw_snippet_picture(pmt_position_class, snippet_class.time_snippets[snippet], muon_points, snippet,
                                 out_path, mode="absolute")
            draw_snippet_contour_plot(pmt_position_class, snippet_class.time_snippets[snippet],
                                      muon_points, snippet, out_path, number_contour_level, mode="absolute")
            if len(snippet_array) is 1:
                draw_snippet_contour_plot(pmt_position_class, snippet_class.time_snippets[snippet],
                                          muon_points, snippet, out_path, number_contour_level, mode="total")

    else:
        for snippet in range(len(snippet_array)):
            statusAlert.processStatus("processing snippet: " + str(snippet))

            '''Draw the detector picture for the certain time snippet'''
            draw_snippet_picture(pmt_position_class, snippet_class.time_snippets[snippet], muon_points, snippet,
                                 out_path, reco_points, "absolute")
            draw_snippet_contour_plot(pmt_position_class, snippet_class.time_snippets[snippet], muon_points, snippet,
                                      out_path, number_contour_level, reco_points, mode="absolute")
            if len(snippet_array) is 1:
                draw_snippet_contour_plot(pmt_position_class, snippet_class.time_snippets[snippet],
                                          muon_points, snippet, out_path, number_contour_level, reco_points, mode="total")


def snippet_drawer_difference(pmt_position_class, snippet_class, muon_points, out_path,
                              number_contour_level, reco_points=None):
    statusAlert.processStatus("Iterate snippets and draw")

    '''photon data pre-processed and ordered into time snippets'''
    snippet_array = np.asarray(snippet_class.time_snippets)

    '''iterate snippets'''
    if reco_points is None:
        for snippet in range(len(snippet_array)):
            if snippet > 0:
                snippet_diff = np.asarray(snippet_class.time_snippets[snippet]) \
                               - np.asarray(snippet_class.time_snippets[snippet-1])
                statusAlert.processStatus("processing snippet: " + str(snippet))

                '''Draw the detector picture for the certain time snippet'''
                draw_snippet_picture(pmt_position_class, snippet_diff, muon_points,
                                     snippet, out_path, "differential")
                draw_snippet_contour_plot(pmt_position_class, snippet_diff, muon_points,
                                          snippet, out_path, number_contour_level, "differential")
    else:
        for snippet in range(len(snippet_array)):
            if snippet > 0:
                snippet_diff = np.asarray(snippet_class.time_snippets[snippet]) \
                               - np.asarray(snippet_class.time_snippets[snippet - 1])
                statusAlert.processStatus("processing snippet: " + str(snippet))

                '''Draw the detector picture for the certain time snippet'''
                draw_snippet_picture(pmt_position_class, snippet_diff, muon_points, snippet,
                                     out_path, reco_points, "differential")
                draw_snippet_contour_plot(pmt_position_class, snippet_diff, muon_points, snippet,
                                          out_path, number_contour_level, reco_points, "differential")


def draw_snippet_picture(pmt_position_class, snippet_array, muon_points, snippet, out_path, reco_points=None, mode=None):
    statusAlert.processStatus("Creating graphics")
    fig = plt.figure(num=None, figsize=(20, 10))

    '''Analysis design'''
    ax1 = fig.add_subplot(111, axisbg='gainsboro', projection='aitoff')
    # scatterHits = ax1.scatter(pmt_position_class.phi_position, pmt_position_class.theta_position, marker='o', c='k')
    scatterHits = ax1.scatter(pmt_position_class.phi_position, pmt_position_class.theta_position2,
                              c=snippet_array, cmap=color_schemes.analysis_design(),
                              edgecolor='k', label='hit', )
    cb = fig.colorbar(scatterHits)

    '''Nice design'''
    # ax1 = fig.add_subplot(111)
    # scatter_hits_new = ax1.scatter(pmt_position_class.phi_hammer_aitov, pmt_position_class.theta_hammer_aitov,
    #                                c=snippet_array, cmap=color_schemes.analysis_design(),
    #                                edgecolor='k', label='hit')
    # cb = fig.colorbar(scatter_hits_new)

    # draw_sector_lines(pmt_position_class)
    draw_muon_points(muon_points, ax1)
    if reco_points is None:
        a=1
    else:
        draw_reconstructed_points(reco_points, ax1)
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

    # elif mode is "total":
    #     cb.set_label("Number of photons")
    #
    #     plt.savefig(out_path + str(file_number) + "_" + str(snippet) + ".png", bbox_inches='tight')
    #
    #     plt.close()

    else:
        print("Nothing to print")


def contour_data_reader(pmt_position_class, snippet_array, number_contour_level):
    statusAlert.processStatus("Collect contour data")
    ax = plt.subplot(111, projection='aitoff')

    phi_i = np.linspace(-math.pi, math.pi, 1200)
    theta_i = np.linspace(-math.pi / 2, math.pi / 2, 600)

    phi_position = pmt_position_class.phi_position
    theta_position2 = pmt_position_class.theta_position2

    zi = plt.mlab.griddata(phi_position, theta_position2,
                           snippet_array, phi_i, theta_i, interp='linear')
    zi = gaussian_filter(zi, 5)

    cont_plot_axes = ax.contour(phi_i, theta_i, zi, number_contour_level)

    contour_data = contour_analyze.get_contour_data(cont_plot_axes)

    return contour_data


def draw_snippet_contour_plot(pmt_position_class, snippet_array, muon_points, snippet, out_path, number_contour_level,
                              reco_points=None, mode=None):
    statusAlert.processStatus("Creating graphics")
    fig = plt.figure(num=None, figsize=(20, 10))

    '''Analysis design'''

    phi_i = np.linspace(-math.pi, math.pi, 1200)
    theta_i = np.linspace(-math.pi/2, math.pi/2, 600)

    phi_position_draw = pmt_position_class.phi_position
    theta_position2_draw = pmt_position_class.theta_position2
    # phi_position_draw = pmt_position_class.phi_shifted
    # theta_position2_draw = pmt_position_class.theta_shifted

    # threshold = max([max(snippet_array)/2, 200])
    # print(len(snippet_array))
    # snippet_array = [max(threshold, min(value, max(snippet_array))) for value in snippet_array]
    # print(len(snippet_array))
    zi = plt.mlab.griddata(phi_position_draw, theta_position2_draw,
                           snippet_array, phi_i, theta_i, interp='linear')

    zi = gaussian_filter(zi, 5)

    ax = plt.subplot(111, axisbg='gainsboro', projection='aitoff')
    ax.grid(True)
    cont_plot_axes = ax.contour(phi_i, theta_i, zi, number_contour_level, cmap=color_schemes.analysis_design())

    contour_data = contour_analyze.get_contour_data(cont_plot_axes)

    plt.ylabel("theta (deg)")
    plt.xlabel("phi (deg)")

    try:
        plt.colorbar(cont_plot_axes)
    except:
        print("No colorbar possible! ")

    draw_muon_points(muon_points, ax)
    if reco_points is None:
        pass
    else:
        draw_reconstructed_points(reco_points, ax)

    if mode is "absolute":
        try:
            os.chdir(out_path + "absolute hits contour/")
        except:
            os.makedirs(out_path + "absolute hits contour/")
        plt.savefig(out_path + "absolute hits contour/" + str(snippet) + ".png", bbox_inches='tight')

        plt.close()

    elif mode is "differential":
        try:
            os.chdir(out_path + "differential hits contour/")
        except:
            os.makedirs(out_path + "differential hits contour/")
        plt.savefig(out_path + "differential hits contour/" + str(snippet) + ".png", bbox_inches='tight')

        plt.close()

    # elif mode is "total":
    #     plt.savefig("/media/gpu/Data1/Analysis/Output/muReconstruction/run2/lala/" + str(file) + "_" + str(snippet) + "_Contour.png", bbox_inches='tight')
    #
    #     plt.close()

    else:
        print("Nothing to print")

    plt.close()


def draw_muon_points(muon_points, subplott):
    phi_muon_first = []
    theta_muon_first = []
    phi_muon_sec = []
    theta_muon_sec = []

    for muon_event in muon_points:
        for point in muon_event:
            if point.event is 0:
                phi_muon_first.append(point.phi)
                theta_muon_first.append(point.theta2)
            if point.event is 1:
                phi_muon_sec.append(point.phi)
                theta_muon_sec.append(point.theta2)

    scatter_muon_entry = subplott.scatter(phi_muon_first, theta_muon_first,
                                     facecolors='none', edgecolors='white', marker='v', s=120)
    scatter_muon_exit = subplott.scatter(phi_muon_sec, theta_muon_sec,
                                   facecolors='none', edgecolors='white', marker='^', s=120)


def draw_reconstructed_points(reco_points, subplott):
    reco_phi = []
    reco_theta = []

    for point in reco_points:
        reco_phi.append(point.x_coordinate_rad)
        reco_theta.append(point.y_coordinate_rad)
    scatter_reco_point = subplott.scatter(reco_phi, reco_theta,
                                          edgecolors='black', facecolors='none', marker='o', s=120)


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