import statusAlert, recoPreparation, color_schemes, contour_analyze, reco_from_contour, PointVecDist
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.ndimage.filters import gaussian_filter
import os


def reco_resulter(difference_array, outdir):
    result_array = []
    for event in difference_array:
        for point_diff in event:
            result_array.append(point_diff)
    n, bins, patches = plt.hist(result_array, bins=20, range=(0, 30))

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    plt.xlabel(r'$\Delta\phi$ ($^\circ$)')
    plt.ylabel('Detected points')
    plt.grid(True)

    plt.savefig(outdir + "reco_accuracy.pdf", bbox_inches='tight')
    plt.savefig(outdir + "reco_accuracy.png", bbox_inches='tight')

    plt.close()


def reco_comparer(truth_points, reco_points):
    min_differences = []
    for truth_point in truth_points:
        differences = []
        for reco_point in reco_points:
            differences.append(math.sqrt(pow(truth_point[0]-reco_point.x_coordinate_deg, 2)
                                         + pow(truth_point[1]-reco_point.y_coordinate_deg, 2)))
            differences.append(math.sqrt(pow(truth_point[0] - reco_point.x_coordinate_deg - 360, 2)
                                         + pow(truth_point[1] - reco_point.y_coordinate_deg, 2)))
            differences.append(math.sqrt(pow(truth_point[0] - reco_point.x_coordinate_deg + 360, 2)
                                         + pow(truth_point[1] - reco_point.y_coordinate_deg, 2)))
        try:
            min_differences.append(min(differences))
        except ValueError:
            pass
    return min_differences


def found_intersects(intersec_array, outdir):
    n, bins, patches = plt.hist(intersec_array, bins=10, range=(0, 10))

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    plt.xlabel("Found intersections")
    plt.ylabel("In event")
    plt.grid(True)

    eff = n[4]/n.sum()
    plt.figtext(0.8, 0.7, "Eff: %.2f" % eff)

    plt.savefig(outdir + "intersec_found.pdf", bbox_inches='tight')
    plt.savefig(outdir + "intersec_found.png", bbox_inches='tight')

    plt.close()


def found_tracks(track_array, outdir):
    n, bins, patches = plt.hist(track_array, bins=5, range=(0, 5))

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    plt.xlabel("Found track solutions")
    plt.ylabel("No.")
    plt.grid(True)

    plt.savefig(outdir + "tracks_found.pdf", bbox_inches='tight')
    plt.savefig(outdir + "tracks_found.png", bbox_inches='tight')

    plt.close()


def track_dist_hist(diff_array, outdir):
    n, bins, patches = plt.hist(diff_array, bins=100, range=(-1000, 1000))

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    plt.xlabel("DistM MC - Reco")
    plt.ylabel("No.")
    plt.grid(True)

    sum = 0
    for dist, entries in enumerate(n):
        sum += dist*entries

    rms = sum/len(n)

    plt.figtext(0.8, 0.7, "RMS: " + str(rms))
    plt.savefig(outdir + "track_reco.pdf", bbox_inches='tight')
    plt.savefig(outdir + "track_reco.png", bbox_inches='tight')

    plt.close()


def reco_accuracy(reco_data, mc_data, accuracy_array):
    for rec_muon_track in reco_data:
        diff = np.inf

        correction = calc_geometric_correction(rec_muon_track.entry_point, rec_muon_track.exit_point)
        print "Correction: "
        print correction
        correction = 0
        for mc_muon_track in mc_data:
            distance = abs(mc_muon_track.distance_track_to_center) - abs(rec_muon_track.distance_track_to_center - 2*correction)
            if abs(distance) < abs(diff):
                diff = distance
                # print "Rec:"
                # print rec_muon_track.distance_track_to_center
                # print "Truth: "
                # print mc_muon_track.distance_track_to_center
            else:
                pass
        accuracy_array.append(diff)
        print "True diff: "
        print diff
        # accuracy_array.append(mc_muon_track.distance_track_to_center)


def calc_geometric_correction(entry_point, exit_point):
    d3entry = PointVecDist.D3Vector()
    d3entry.x = entry_point.real_x
    d3entry.y = entry_point.real_y
    d3entry.z = entry_point.real_z

    d3exit = PointVecDist.D3Vector()
    d3exit.x = exit_point.real_x
    d3exit.y = exit_point.real_y
    d3exit.z = exit_point.real_z

    k = math.sqrt(pow(17700, 2) / (pow(d3entry.x, 2) + pow(d3entry.y, 2) + pow(d3entry.z, 2)))
    new_point = PointVecDist.D3Vector()
    new_point.x = d3entry.x * k
    new_point.y = d3entry.y * k
    new_point.z = d3entry.z * k

    k = math.sqrt(pow(17700, 2) / (pow(d3exit.x, 2) + pow(d3exit.y, 2) + pow(d3exit.z, 2)))
    new_point = PointVecDist.D3Vector()
    new_point.x = d3exit.x * k
    new_point.y = d3exit.y * k
    new_point.z = d3exit.z * k

    correction = PointVecDist.calc_distance_line_point(d3entry, d3exit, new_point)
    return correction







