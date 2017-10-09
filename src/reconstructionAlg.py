import statusAlert, recoPreparation, color_schemes
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


class PatternPosition:
    def __init__(self):  # this method creates the class object.
        self.hits = 0
        self.phi = 0
        self.theta = 0


def pattern_detector(pmt_position_class, snippet_class, muon_points, out_path, scanning_cut_radius, sector_threshold,
                     cut_threshold):
    statusAlert.processStatus("search patterns in snippet")

    '''photon data pre-processed and ordered into time snippets'''
    snippet_array = np.asarray(snippet_class.time_snippets)

    '''This array saves all found patterns'''
    snippet_pattern_array = []

    muon_finder = muon_entry_exit_reconstructor(pmt_position_class, snippet_class, out_path)

    '''iterate snippets'''
    for snippet in range(50):
        statusAlert.processStatus("processing snippet: " + str(snippet))

        '''sum of all photon hits in this time snippet'''
        sector_pattern_array = pattern_scan(pmt_position_class, snippet_array, snippet, scanning_cut_radius,
                                            sector_threshold, cut_threshold)

        snippet_pattern_array.append(sector_pattern_array)
        '''Draw the detector picture for the certain time snippet'''
        draw_snippet_picture(pmt_position_class, snippet_class, muon_points, snippet,
                             sector_pattern_array, muon_finder[snippet], scanning_cut_radius, out_path)


def pattern_scan(pmt_position_class, snippet_array, snippet, template_radius, sector_threshold, cut_threshold):
    '''this array saves all patterns within one snippet'''
    sector_pattern_array = []
    '''iterate sectors'''
    for sector in range(len(pmt_position_class.sector_list_id)):

        '''get data into numpy arrays'''
        id_np = np.asarray(pmt_position_class.sector_list_id[sector])
        phi_np = np.asarray(pmt_position_class.sector_list_phi[sector])
        theta_np = np.asarray(pmt_position_class.sector_list_theta[sector])

        total_hit_sum = sum(snippet_array[snippet])

        '''initialize class, which savespattern position. This class is given to the sector_pattern_array'''
        pattern_position = PatternPosition()

        '''get sum of photon hits in each sector'''
        sector_hit_sum = 0
        for i in range(len(id_np)):
            sector_hit_sum += snippet_array[snippet][id_np[i]]

        if float(sector_hit_sum) > sector_threshold * total_hit_sum:
            '''scan the template over the sector'''
            cut_phi_position = min(phi_np)
            while cut_phi_position < max(phi_np) + template_radius:
                cut_theta_position = min(theta_np)
                while cut_theta_position < max(theta_np) + template_radius:
                    '''sum of photons within the template'''
                    hit_sum = 0
                    '''scan over all pmts in the sector and check if they are inside the template cut'''
                    # faster to get pmts directly instead of iterating all in sector
                    for i in range(len(id_np)):
                        if (pow(template_radius, 2) > (pow(phi_np[i] - cut_phi_position, 2)
                                                           + pow(theta_np[i] - cut_theta_position, 2))):
                            hit_sum += snippet_array[snippet][id_np[i]]
                    '''check if photon sum is over threshold and if its the "biggest" pattern in the sector.
                     Then overwrite previous data'''

                    if hit_sum > cut_threshold * sector_hit_sum and hit_sum > pattern_position.hits:
                        pattern_position.hits = hit_sum
                        pattern_position.phi = cut_phi_position
                        pattern_position.theta = cut_theta_position

                    cut_theta_position += template_radius
                cut_phi_position += template_radius

        sector_pattern_array.append(pattern_position)
    return sector_pattern_array


def muon_entry_exit_reconstructor(pmt_position_class, snippet_class, out_path):
    statusAlert.processStatus("Searching for entry and exit points")

    '''photon data pre-processed and ordered into time snippets'''
    snippet_array = np.asarray(snippet_class.time_snippets)

    '''sector count data is saved here for each snippet'''
    snippet_sector_array = [[] for _ in range(len(snippet_array))]
    for snippet in range(len(snippet_array)):
        sector_array = [[] for _ in range(len(pmt_position_class.sector_list_id))]
        for sector in range(len(pmt_position_class.sector_list_id)):
            id_np = np.asarray(pmt_position_class.sector_list_id[sector])
            sector_hit_sum = 0
            for i in range(len(id_np)):
                sector_hit_sum += snippet_array[snippet][id_np[i]]
            sector_array[sector] = (sector_hit_sum)
        snippet_sector_array[snippet] = (sector_array)

    '''iterate over snippet, sector and the pmts of the sector to get its hit_sum'''
    snippet_sector_pattern_array = []
    for snippet in range(len(snippet_array)):

        sector_pattern_array = []
        for sector in range(len(pmt_position_class.sector_list_id)):
            pattern_position = PatternPosition()

            id_np = np.asarray(pmt_position_class.sector_list_id[sector])
            phi_np = np.asarray(pmt_position_class.sector_list_phi[sector])
            theta_np = np.asarray(pmt_position_class.sector_list_theta[sector])

            dist_sum = 0.0
            phi_sum = 0.0
            theta_sum = 0.0
            counter = 0.0
            sector_hit_sum = 0

            for i in range(len(id_np)):
                sector_hit_sum += snippet_array[snippet][id_np[i]]
                if snippet_array[snippet][id_np[i]] > 10:
                    dist_sum += math.sqrt(pow(phi_np[i]-phi_np[0], 2) + pow(theta_np[i]-theta_np[0], 2))
                    phi_sum += phi_np[i]
                    theta_sum += theta_np[i]
                    counter += 1.0
            if snippet > 0:
                if 4 > counter > 0 and dist_sum/counter < 0.1 and \
                                snippet_sector_array[snippet][sector] > 50 * snippet_sector_array[snippet-1][sector]:
                    pattern_position.phi = phi_sum / counter
                    pattern_position.theta = theta_sum / counter

            sector_pattern_array.append(pattern_position)

        snippet_sector_pattern_array.append(sector_pattern_array)

    return snippet_sector_pattern_array


    # for i in range(len(snippet_sector_array)):
    #     if i > 0:
    #         for j in range(len(snippet_sector_array[i])):
    #             if snippet_sector_array[i][j] > 150 * snippet_sector_array[i-1][j]:
    #                 print("Snippet: " + str(i))
    #                 print ("Sector: " + str(j))
    #                 print("_____________________")
    #         print(":::::::::::::::::::::::")


def draw_snippet_picture(pmt_position_class, snippet_class, muon_points, snippet, sector_pattern_array,
                         muon_finder_psnippet, template_radius, out_path):
    c = mcolors.ColorConverter().to_rgb

    color_analysis = color_schemes.make_colormap([c('crimson'), c('yellow'), 0.05,
                                        c('yellow'), c('dodgerblue'), 0,65,
                                        c('dodgerblue'), c('blue'), 0.85,
                                        c('blue'), c('darkblue'), 0.90, c('darkblue')])

    color_corporate = color_schemes.make_colormap([color_schemes.corporate_design('karminrot'),
                                                   color_schemes.corporate_design('white'), 0.5,
                                                   color_schemes.corporate_design('white')])

    statusAlert.processStatus("Creating graphics")
    fig = plt.figure(num=None, figsize=(20, 10))

    '''Analysis design'''
    ax1 = fig.add_subplot(111, axisbg='gainsboro')
    scatterHits = ax1.scatter(pmt_position_class.phi_position, pmt_position_class.theta_position,
                               c=snippet_class.time_snippets[snippet], cmap=color_analysis, label='hit')
    cb = fig.colorbar(scatterHits)

    '''Nice design'''
    # ax1 = fig.add_subplot(111)
    # scatter_hits_new = ax1.scatter(pmt_position_class.phi_hammer_aitov, pmt_position_class.theta_hammer_aitov,
    #                                c=snippet_class.time_snippets[snippet], cmap=color_analysis, label='hit')
    # cb = fig.colorbar(scatter_hits_new)

    draw_sector_lines(pmt_position_class)
    draw_muon_points(muon_points)
    # draw_pattern_results(sector_pattern_array, template_radius, ax1)
    draw_reco_muons(muon_finder_psnippet)

    plt.ylabel("theta (deg)")
    plt.xlabel("phi (deg)")

    cb.set_label("Number of Photons")

    # plt.savefig(out_path + str(snippet) + ".pdf")
    plt.savefig(out_path + str(snippet) + ".png")

    plt.close()


def draw_muon_points(muonPoints):
    phi_muon_entry = []
    theta_muon_entry = []
    phi_muon_exit = []
    theta_muon_exit = []

    for i in range(len(muonPoints)):
        for j in range(len(muonPoints[i])):
            if j == 0:
                phi_muon_entry.append(recoPreparation.calcPMTPolarPhi(muonPoints[i][j]))
                theta_muon_entry.append(recoPreparation.calcPMTPolarTheta(muonPoints[i][j]))
            else:
                phi_muon_exit.append(recoPreparation.calcPMTPolarPhi(muonPoints[i][j]))
                theta_muon_exit.append(recoPreparation.calcPMTPolarTheta(muonPoints[i][j]))
    scatter_muon_entry = plt.scatter(phi_muon_entry, theta_muon_entry,
                                     facecolors='none', edgecolors='white', marker='v', s=120)
    scatter_muon_exit = plt.scatter(phi_muon_exit, theta_muon_exit,
                                   facecolors='none', edgecolors='white', marker='^', s=120)


def draw_sector_lines(pmt_position_class):
    for i in range(pmt_position_class.x_sectors + 1):
        plt.plot([-math.pi + i * (2 * math.pi) / pmt_position_class.x_sectors,
                  -math.pi + i * (2 * math.pi) / pmt_position_class.x_sectors], [0, math.pi], 'k-')

    for i in range(pmt_position_class.y_sectors + 1):
        plt.plot([-math.pi, math.pi],
                 [i * math.pi / pmt_position_class.y_sectors, i * math.pi / pmt_position_class.y_sectors], 'k-')


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