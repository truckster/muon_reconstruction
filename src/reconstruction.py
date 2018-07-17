import recoPreparation, reconstructionAlg, statusAlert, TreeReadFunc, gauss_fit_reco, contour_analyze
import total_event_reconstruction, intersec_time_finder, point_allocate, diff_event_analysis, performance_check

import cPickle as pickle
from os import chdir, remove, path, getcwd
from glob import glob
import gc
import numpy as np

'''General script to use sub-scripts for muon reconstruction.'''
statusAlert.processStatus("Process started")

# input_path = "/home/gpu/Simulation/processed_sim/lala/"
# input_path = "/home/gpu/Simulation/processed_sim/Vortrag/"
input_path = "/home/gpu/Simulation/processed_sim/lala2/"

output_path = "/home/gpu/Analysis/muReconstruction/Output/dev/"
# output_path = "/home/gpu/Analysis/muReconstruction/Output/Vortrag/"
# output_path = "/home/gpu/Analysis/muReconstruction/Output/LPMT/"

# input_path = "/home/gpu/Simulation/presentation/y/"
# output_path = "/home/gpu/Analysis/muReconstruction/Output/presentation/y/"

'''Overwrite fit results file'''
if path.isfile(output_path + "results.txt"):
    remove(output_path + "results.txt")

reco_accuracy = []
found_point_array = []
found_track_array = []
det_points = 0

chdir(input_path)
for folder in glob("*/"):
    new_output_path = recoPreparation.create_output_path(output_path, folder, "/totalEventHist/",  input_path)
    new_output_path_fit = recoPreparation.create_output_path(output_path, folder, "/fits/",  input_path)
    total_path = recoPreparation.create_output_path(output_path, folder, "/total_event/", input_path)

    file_number = folder.split("-")[0]

    chdir(input_path + folder)
    statusAlert.processStatus("Reading file: " + str(folder))

    '''calculate PMT positions for this file'''
    PmtPositions = pickle.load(open("PMT_positions.pkl", 'rb'))

    '''collect entry and exit points of all muons in event'''
    intersec_radius = 19600
    time_resolution = 1*10**-9
    muon_points = pickle.load(open("muon_truth.pkl", 'rb'))

    '''collect information of all photons within certain time snippet and save the separately'''
    photons_in_time_window = pickle.load(open("framed_photons.pkl", 'rb'))
    photons_of_entire_event = pickle.load(open("total_event_photons.pkl", 'rb'))
    # photon_data = pickle.load(open("photon_data_array.pkl", 'rb'))
    frame_time_cut = 10
    number_contour_level = 8

    contour_array_total = pickle.load(open("Contours/contour_array_total.pkl"))
    contour_array_diff = pickle.load(open("Contours/contour_array_diff.pkl"))
    contour_array_total_dPhi = pickle.load(open("Contours/contour_array_total_phi_rotate.pkl"))
    contour_array_diff_dPhi = pickle.load(open("Contours/contour_array_diff_phi_rotate.pkl"))
    contour_array_total_dTheta = pickle.load(open("Contours/contour_array_total_theta_rotate.pkl"))
    contour_array_diff_dTheta = pickle.load(open("Contours/contour_array_diff_theta_rotate.pkl"))
    contour_array_total_shifted = pickle.load(open("Contours/contour_array_total_rotate.pkl"))
    contour_array_diff_shifted = pickle.load(open("Contours/contour_array_diff_rotate.pkl"))

    contour_entire_event = pickle.load(open("Contours/event_photons.pkl"))
    contour_entire_event_dPhi = pickle.load(open("Contours/event_photons_dPhi.pkl"))
    contour_entire_event_dTheta = pickle.load(open("Contours/event_photons_dTheta.pkl"))
    contour_entire_event_shifted = pickle.load(open("Contours/event_photons_shifted.pkl"))

    MC_positions = recoPreparation.MC_truth_writer(muon_points, output_path, folder, frame_time_cut)

    # total_event = [contour_entire_event, contour_entire_event_dPhi, contour_entire_event_dTheta]
    total_event = [contour_entire_event, contour_entire_event_shifted]
    # total_event = [contour_entire_event]
    framed_event = [contour_array_total, contour_array_total_dPhi, contour_array_total_dTheta]
    # diff_event = [contour_array_diff, contour_array_diff_dPhi, contour_array_diff_dTheta]
    diff_event = [contour_array_diff, contour_array_diff_shifted]

    '''Reconstruction by looking at entire event'''
    found_points = total_event_reconstruction.entry_exit_detector(total_event)
    '''Reconstruction by looking at differences between two frames'''
    # found_points_2 = diff_event_analysis.entry_exit_detector(diff_event)
    """Cross-check results with diffs"""
    # found_points.extend(found_points_2)
    found_points = diff_event_analysis.intersec_crosscheck(diff_event, found_points)
    found_points = reconstructionAlg.coordinate_calculation(found_points)
    found_points = reconstructionAlg.orientation_resolver(found_points)
    diff_event_analysis.point_merger(found_points, 150)

    print("Found points: " + str(len(found_points)))

    if len(found_points) > 3:
        det_points += len(found_points)

    '''Detection of entry and exit time of muons'''
    intersec_time_finder.find_times(contour_array_total, found_points)

    '''Allocate respective points'''
    # point_allocate.allocate_points(contour_array_diff, found_points, found_frames)
    track_class = point_allocate.allocate_tracks_to_points(found_points)
    total_event_reconstruction.reco_result_writer(output_path, found_points)
    found_point_array.append(len(found_points))
    try:
        found_track_array.append(len(track_class[0]))
    except:
        pass

    '''Draw all kinds of images'''
    reconstructionAlg.snippet_drawer(PmtPositions, photons_of_entire_event, muon_points, total_path,
                                     number_contour_level, found_points)
    # reconstructionAlg.snippet_drawer(PmtPositions, photons_in_time_window,
    #                                  muon_points, new_output_path, number_contour_level, found_points)
    # reconstructionAlg.snippet_drawer_difference(PmtPositions, photons_in_time_window,
    #                                             muon_points, new_output_path, number_contour_level, found_points)

    # reconstructionAlg.print_sector_pmts(PmtPositions, output_path)

    if len(found_points) is 4:
        reco_accuracy.append(performance_check.reco_comparer(MC_positions, found_points))



print(found_point_array)
print(found_track_array)
performance_check.reco_resulter(reco_accuracy, output_path)
performance_check.found_intersects(found_point_array, output_path)
print("Long events: " + str(len(reco_accuracy)))
print("Found points: " + str(det_points))
statusAlert.processStatus("Process finished")

