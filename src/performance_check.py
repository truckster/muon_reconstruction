import statusAlert, recoPreparation, color_schemes, contour_analyze, reco_from_contour
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
    # if len(truth_points) is not len(reco_points):

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