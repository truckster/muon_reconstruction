import statusAlert, recoPreparation, color_schemes, contour_analyze, reco_from_contour, PointVecDist
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.colors as mcolors
from scipy.ndimage.filters import gaussian_filter
import os


def do_performance_checks(performance_class, outdir):
    point_accuracy(performance_class, outdir)
    point_accuracy_d_dep(performance_class, outdir)
    point_track_accuracy(performance_class, outdir)
    found_intersects(performance_class.found_point_list, outdir)
    found_tracks(performance_class.tracks_per_event_list, outdir)
    track_dist_hist(performance_class, outdir)
    track_angle_diff_hist(performance_class, outdir)
    merge_performance(performance_class, outdir)
    reco_accuracy_d_dep_print(performance_class, outdir)
    shower_dependence(performance_class, outdir)
    time_performance(performance_class, outdir)
    parallelity_performance(performance_class, outdir)


def point_accuracy(performance_class, outdir):
    result_array = []
    for event in performance_class.point_reco_accuracy_list:
        for point_diff in event:
            result_array.append(point_diff)
    n, bins, patches = plt.hist(result_array, bins=20, range=(0, 30))
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel(r'$\Delta\phi$ ($^\circ$)')
    plt.ylabel('Detected points')
    plt.grid(True)
    # plt.savefig(outdir + "reco_accuracy.pdf", bbox_inches='tight')
    plt.savefig(outdir + "point_accuracy.png", bbox_inches='tight')
    plt.close()


def point_accuracy_d_dep(performance_class, outdir):
    point_accuracy_list = []
    truth_track_dist_list = []
    for event_index, event in enumerate(performance_class.point_reco_accuracy_list):
        for point_diff in event:
            truth_track_dist_list.append(performance_class.mc_track_distance_list[event_index])
            point_accuracy_list.append(point_diff)

    plt.scatter(truth_track_dist_list, point_accuracy_list, marker="x")
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel("$D_{MC}$ (mm)")
    plt.ylabel("$\Delta\phi$ ($^\circ$)")
    plt.grid(True)
    # plt.savefig(outdir + "correction_d_dep.pdf", bbox_inches='tight')
    plt.savefig(outdir + "point_accuracy_over_mc_distance_scatter.png", bbox_inches='tight')
    plt.close()

    plt.hist2d(truth_track_dist_list, point_accuracy_list, bins=[20, 100], cmap=color_schemes.analysis_design(),
               range=[[0, 17700], [0, 20]])
    plt.xlabel("$D_{MC}$ (mm)")
    plt.ylabel("$\Delta\phi$ ($^\circ$)")
    plt.savefig(outdir + "point_accuracy_over_mc_distance.png", bbox_inches='tight')
    plt.close()


def point_track_accuracy(performance_class, outdir):
    point_accuracy_list = []
    track_accuracy_list = []
    for event_index, event in enumerate(performance_class.point_reco_accuracy_list):
        for point_diff in event:
            track_accuracy_list.append(performance_class.track_reco_accuracy_list_D[event_index])
            point_accuracy_list.append(point_diff)

    plt.scatter(track_accuracy_list, point_accuracy_list, marker="x")
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel("$\DeltaD$ (mm)")
    plt.ylabel(r'$\Delta\phi$ ($^\circ$)')
    plt.grid(True)
    # plt.savefig(outdir + "correction_d_dep.pdf", bbox_inches='tight')
    plt.savefig(outdir + "point_accuracy_over_rec_accuracy.png", bbox_inches='tight')
    plt.close()


def found_intersects(intersec_array, outdir):
    n, bins, patches = plt.hist(intersec_array, bins=10, range=(0, 10))
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel("Found intersections")
    plt.ylabel("Events")
    plt.grid(True)
    plt.figtext(0.7, 0.75, "Entries: %.0f" % len(intersec_array))
    eff = n[4]/n.sum()
    # plt.figtext(0.7, 0.7, "Eff: %.2f" % eff)
    # plt.savefig(outdir + "intersec_found.pdf", bbox_inches='tight')
    plt.savefig(outdir + "intersec_found.png", bbox_inches='tight')
    plt.close()


def found_tracks(track_array, outdir):
    n, bins, patches = plt.hist(track_array, bins=5, range=(0, 5))
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel("Found track solutions")
    plt.ylabel("No.")
    plt.grid(True)
    # plt.savefig(outdir + "tracks_found.pdf", bbox_inches='tight')
    plt.savefig(outdir + "tracks_found.png", bbox_inches='tight')
    plt.close()


def track_dist_hist(performance_class, outdir):
    n, bins, patches = plt.hist(performance_class.track_reco_accuracy_list_D, bins=100, range=(-1000, 1000))
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel(r"\textit{\Delta D} (mm)")
    plt.ylabel("# tracks")
    plt.grid(True)

    def gaussian(x, a, mean, sigma):
        return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))
    popt, pcov = curve_fit(gaussian, bins[:-1], n, [200, 0, 10])
    i = np.linspace(-1000, 1000, 2001)
    plt.plot(i, gaussian(i, *popt))

    plt.figtext(0.7, 0.75, "Entries: %.0f" % len(performance_class.track_reco_accuracy_list_D))
    plt.figtext(0.7, 0.7, "Center: %.2f" % popt[1])
    plt.figtext(0.7, 0.65, "FWHM: %.2f" % (2.3548*popt[2]))
    # plt.savefig(outdir + "track_reco.pdf", bbox_inches='tight')
    plt.savefig(outdir + "track_reco.png", bbox_inches='tight')

    plt.close()


def track_angle_diff_hist(performance_class, outdir):
    n, bins, patches = plt.hist(performance_class.track_reco_accuracy_list_phi, bins=100, range=(0, 2.50))
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel(r'$\Delta\phi$ (rad)')
    plt.ylabel("count")
    plt.grid(True)

    # plt.savefig(outdir + "track_reco.pdf", bbox_inches='tight')
    plt.savefig(outdir + "track_reco_angle.png", bbox_inches='tight')

    plt.close()


def points_to_track(performance_class, outdir):
    n, bins, patches = plt.hist(performance_class.found_points_for_good_reco)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel("No. of Points")
    plt.ylabel("Good reconstruction events")
    plt.grid(True)
    plt.figtext(0.7, 0.75, "Entries: %.0f" % len(performance_class.found_points_for_good_reco))
    plt.figtext(0.7, 0.7, "Efficiency: %.2f" % float((len(performance_class.found_points_for_good_reco))
                / float(len(performance_class.found_point_list))))
    plt.savefig(outdir + "pointsForTracks.png", bbox_inches='tight')
    # plt.savefig(outdir + "pointsForTracks.pdf", bbox_inches='tight')
    plt.close()


def reco_accuracy_d_dep_print(performance_class, outdir):
    plt.scatter(performance_class.mc_track_distance_list, performance_class.track_reco_accuracy_list_D, marker="x")
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel("$D_{MC}$ (mm)")
    plt.ylabel("$\Delta D_{rec-truth}$ (mm)")
    plt.grid(True)
    # plt.savefig(outdir + "scatter_reco_over_truth.pdf", bbox_inches='tight')
    plt.savefig(outdir + "d_reco_accurady_d_dep.png", bbox_inches='tight')
    plt.close()

    n, x_edges, y_edges, image = plt.hist2d(performance_class.mc_track_distance_list,
                                            performance_class.track_reco_accuracy_list_D,
                                            bins=[20, 100], cmap=color_schemes.analysis_design(),
                                            range=[[0, 14000], [-1000, 1000]])
    x_bin_center_list = (x_edges[:-1] + x_edges[1:]) / 2
    y_bin_center_list = (y_edges[:-1] + y_edges[1:]) / 2
    mean_values_list = []
    var_list = []
    for x_bin in n:
        y_bin_entry_list = []
        for y_bin_index, y_bin in enumerate(x_bin):
            for hit in range(int(y_bin)):
                y_bin_entry_list.append(y_bin_index)
        if not np.isnan(np.mean(y_bin_entry_list)):
            mean_values_list.append(y_bin_center_list[int(np.mean(y_bin_entry_list))])
            var_list.append(np.std(y_bin_entry_list)*((y_edges[0]-y_edges[-1])/len(y_edges)))
        else:
            mean_values_list.append(0.0)
            var_list.append(0.0)
    plt.errorbar(x_bin_center_list, mean_values_list, yerr=var_list, color="k", marker='x', linestyle="")
    # plt.figtext(0.17, 0.8, "Max dev.: %.2f $\pm$ %.2f" % (max(mean_values_list),
    #                                                      var_list[mean_values_list.index(max(mean_values_list))]))
    plt.xlabel("$D_{MC}$ (mm)")
    plt.ylabel("$\Delta D_{rec-truth}$ (mm)")
    plt.savefig(outdir + "d_reco_accurady_d_dep_2D.png", bbox_inches='tight')
    plt.close()

    plt.scatter(performance_class.mc_track_distance_list, performance_class.track_reco_accuracy_list_phi, marker="x")
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel("$D_{MC}$ (mm)")
    plt.ylabel("$\Delta \phi_{rec-truth}$ (deg)")
    plt.grid(True)
    # plt.savefig(outdir + "scatter_reco_over_truth.pdf", bbox_inches='tight')
    plt.savefig(outdir + "phi_reco_accurady_d_dep.png", bbox_inches='tight')
    plt.close()

    n, x_edges, y_edges, image = plt.hist2d(performance_class.mc_track_distance_list,
                                            performance_class.track_reco_accuracy_list_phi, bins=[20, 100],
                                            cmap=color_schemes.analysis_design(), range=[[0, 14000], [0, 3]])

    x_bin_center_list = (x_edges[:-1] + x_edges[1:]) / 2
    y_bin_center_list = (y_edges[:-1] + y_edges[1:]) / 2
    mean_values_list = []
    var_list = []
    for x_bin in n:
        y_bin_entry_list = []
        for y_bin_index, y_bin in enumerate(x_bin):
            for hit in range(int(y_bin)):
                y_bin_entry_list.append(y_bin_index)
        if not np.isnan(np.mean(y_bin_entry_list)):
            mean_values_list.append(y_bin_center_list[int(np.mean(y_bin_entry_list))])
            var_list.append(np.std(y_bin_entry_list) * ((y_edges[0] - y_edges[-1]) / len(y_edges)))
        else:
            mean_values_list.append(0.0)
            var_list.append(0.0)
    plt.errorbar(x_bin_center_list, mean_values_list, yerr=var_list, color="k", marker='x', linestyle="")
    # plt.figtext(0.17, 0.8, "Max dev.: %.2f $\pm$ %.2f" % (max(mean_values_list),
    #                                                      var_list[mean_values_list.index(max(mean_values_list))]))
    plt.xlabel("$D_{MC}$ (mm)")
    plt.ylabel("$\Delta \phi_{rec-truth}$ (deg)")
    plt.savefig(outdir + "phi_reco_accurady_d_dep_2D.png", bbox_inches='tight')
    plt.close()


def merge_performance(performance_class, outdir):
    fig = plt.figure(num=None, figsize=(20, 10))
    ax1 = fig.add_subplot(111)
    hist = ax1.hist2d(performance_class.unmerged, performance_class.merged, cmap=color_schemes.analysis_design(),
               range=[[0, 9], [0, 9]], bins=[9,9])
    cb = fig.colorbar(hist[3])
    plt.xlabel("Intersection points unmerged")
    plt.ylabel("Intersection points merged")
    plt.savefig(outdir + "merge_distri.png", bbox_inches='tight')
    # plt.savefig(outdir + "merge_distri.pdf", bbox_inches='tight')
    plt.close()


def shower_dependence(performance_class, outdir):
    n, bins, patches = plt.hist(performance_class.number_of_points_showering, bins=10, range=(0, 10))
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel("No. of reconstructed intersection points")
    plt.ylabel("count")
    plt.grid(True)
    plt.savefig(outdir + "shower_points.png", bbox_inches='tight')
    # plt.savefig(outdir + "shower_points.pdf", bbox_inches='tight')
    plt.close()

    n, bins, patches = plt.hist(performance_class.number_of_points_not_showering, bins=10, range=(0, 10))
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel("No. of reconstructed intersection points")
    plt.ylabel("count")
    plt.grid(True)
    plt.savefig(outdir + "no_shower_points.png", bbox_inches='tight')
    # plt.savefig(outdir + "shower_points.pdf", bbox_inches='tight')
    plt.close()


def time_performance(performance_class, outdir):
    n, bins, patches = plt.hist(performance_class.time_difference_list, bins=100)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel("$\Delta T_{rec-truth}$ (ns)")
    plt.ylabel("count")
    plt.grid(True)
    plt.savefig(outdir + "time_performance.png", bbox_inches='tight')
    # plt.savefig(outdir + "shower_points.pdf", bbox_inches='tight')
    plt.close()


def parallelity_performance(performance_class, outdir):
    n, bins, patches = plt.hist(performance_class.parallelity_list, bins=20)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel("Parallelity factor")
    plt.ylabel("count")
    plt.grid(True)
    plt.savefig(outdir + "parallelity_performance.png", bbox_inches='tight')
    # plt.savefig(outdir + "shower_points.pdf", bbox_inches='tight')
    plt.close()


def calc_geometric_correction(entry_point, exit_point):
    d3entry = entry_point.D3_vector
    d3exit = exit_point.D3_vector
    radius = 16750
    # radius = 17700

    k = math.sqrt(pow(radius, 2) / (pow(d3entry.x, 2) + pow(d3entry.y, 2) + pow(d3entry.z, 2)))
    new_entry_point = PointVecDist.D3Vector()
    new_entry_point.x = d3entry.x * k
    new_entry_point.y = d3entry.y * k
    new_entry_point.z = d3entry.z * k

    k2 = math.sqrt(pow(radius, 2) / (pow(d3exit.x, 2) + pow(d3exit.y, 2) + pow(d3exit.z, 2)))
    new_exit_point = PointVecDist.D3Vector()
    new_exit_point.x = d3exit.x * k2
    new_exit_point.y = d3exit.y * k2
    new_exit_point.z = d3exit.z * k2

    corrected_dist = PointVecDist.calc_track_dist_to_center(new_entry_point, new_exit_point)
    return corrected_dist


def calc_angle_between_lines(truth, reco):
    truth_direct = PointVecDist.D3Vector()
    truth_direct.x = truth.entry_point.D3_vector.x - truth.exit_point.D3_vector.x
    truth_direct.y = truth.entry_point.D3_vector.y - truth.exit_point.D3_vector.y
    truth_direct.z = truth.entry_point.D3_vector.z - truth.exit_point.D3_vector.z

    reco_direct = PointVecDist.D3Vector()
    reco_direct.x = reco.entry_point.x - reco.exit_point.x
    reco_direct.y = reco.entry_point.y - reco.exit_point.y
    reco_direct.z = reco.entry_point.z - reco.exit_point.z

    theta = PointVecDist.calc_angle_between_d3vecs(truth_direct, reco_direct)
    if theta > math.pi/2:
        return abs(theta-math.pi)
    else:
        return theta


def reco_accuracy(reco_data, mc_data, performance_class):
    for rec_muon_track in reco_data:
        D_diff = np.inf
        angle_diff = 0
        corrected_dist = calc_geometric_correction(rec_muon_track.entry_point, rec_muon_track.exit_point)
        # correction = 0
        for mc_muon_track in mc_data:
            distance = abs(corrected_dist) - abs(mc_muon_track.distance_track_to_center)
            if abs(distance) < abs(D_diff):
                D_diff = distance
                angle_diff = calc_angle_between_lines(rec_muon_track, mc_muon_track)
            else:
                pass
        performance_class.track_reco_accuracy_list_D.append(D_diff)
        performance_class.track_reco_accuracy_list_phi.append(angle_diff/math.pi*180.0)


def reco_accuracy_d_dep(reco_data, mc_data, performance_class):
    for rec_muon_track in reco_data:
        diff = np.inf

        correction = calc_geometric_correction(rec_muon_track.entry_point, rec_muon_track.exit_point)
        # correction = 0
        for mc_muon_track in mc_data:
            distance = abs(rec_muon_track.distance_track_to_center) - abs(mc_muon_track.distance_track_to_center)
            if abs(distance) < abs(diff):
                diff = distance
            else:
                pass
        performance_class.track_reco_accuracy_list_uncorr.append(diff)
        performance_class.mc_track_distance_list.append(mc_muon_track.distance_track_to_center)
        performance_class.reconstructed_track_distance_list.append(rec_muon_track.distance_track_to_center)
        performance_class.distance_correction.append(correction)


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

