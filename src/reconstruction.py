import recoPreparation, reconstructionAlg, statusAlert, TreeReadFunc, gauss_fit_reco, contour_analyze
import total_event_reconstruction, intersec_time_finder, point_allocate, diff_event_analysis

import cPickle as pickle
from os import chdir, remove, path, getcwd
from glob import glob
import gc
import numpy as np

'''General script to use sub-scripts for muon reconstruction.'''
statusAlert.processStatus("Process started")

# input_path = "/home/gpu/Simulation/mult/new/"
input_path = "/home/gpu/Simulation/processed_sim/"
# input_path = "/home/gpu/Simulation/mult/test/"
# input_path = "/home/gpu/Simulation/test/"
# input_path = "/home/gpu/Simulation/test_short/"
# input_path = "/home/gpu/Simulation/single/"

output_path = "/home/gpu/Analysis/muReconstruction/Output/"
# output_path = "/home/gpu/Analysis/muReconstruction/Output/LPMT/"

# input_path = "/home/gpu/Simulation/presentation/y/"
# output_path = "/home/gpu/Analysis/muReconstruction/Output/presentation/y/"

'''Overwrite fit results file'''
if path.isfile(output_path + "results.txt"):
    remove(output_path + "results.txt")

reco_accuracy = []

chdir(input_path)
for folder in glob("*/"):
    chdir(input_path + folder)
    # TODO start new process for each file. This might solve the memory problem.
    statusAlert.processStatus("Reading file: " + str(folder))

    '''Control, which events are useful for the analysis'''

    new_output_path = recoPreparation.create_output_path(output_path, folder, "/totalEventHist/",  input_path)
    new_output_path_fit = recoPreparation.create_output_path(output_path, folder, "/fits/",  input_path)
    '''calculate PMT positions for this file'''
    PmtPositions = pickle.load(open("PMT_positions.pkl", 'rb'))

    '''collect entry and exit points of all muons in event'''
    intersec_radius = 17600
    time_resolution = 1*10**-9
    muon_points = pickle.load(open("muon_truth.pkl", 'rb'))

    '''collect information of all photons within certain time snippet and save the separately'''
    photons_in_time_window = pickle.load(open("framed_photons.pkl", 'rb'))
    photons_of_entire_event = pickle.load(open("total_event_photons.pkl", 'rb'))

    frame_time_cut = 10
    number_contour_level = 10
    contour_array_total, contour_array_diff = contour_analyze.collect_contour_data(photons_in_time_window, PmtPositions,
                                                                                   number_contour_level)
    MC_positions = recoPreparation.MC_truth_writer(muon_points, output_path, folder, frame_time_cut)

    '''Reconstruction by looking at entire event'''
    total_path = recoPreparation.create_output_path(output_path, folder, "/total_event/", input_path)
    found_points = total_event_reconstruction.entry_exit_detector(PmtPositions, photons_of_entire_event,
                                                                  number_contour_level)
    reconstructionAlg.snippet_drawer(PmtPositions, photons_of_entire_event, muon_points, total_path,
                                     number_contour_level, found_points)
    reco_positions = total_event_reconstruction.reco_result_writer(output_path, found_points)

    """Cross-check results with diffs"""
    # diff_event_analysis.intersec_crosscheck(contour_array_diff, found_points)

    '''Detection of entry and exit time of muons'''
    found_frames = intersec_time_finder.find_times(contour_array_total, found_points)
    intersec_time_finder.reco_result_writer(output_path, found_frames)

    '''Allocate respective points'''
    # point_allocate.allocate_points(contour_array_diff, found_points, found_frames)

    # '''Draw all kinds of images'''
    # reconstructionAlg.snippet_drawer(PmtPositions, photons_in_time_window,
    #                                  muon_points, new_output_path, number_contour_level, found_points)
    #
    # number_contour_level = 5
    # reconstructionAlg.snippet_drawer_difference(PmtPositions, photons_in_time_window,
    #                                             muon_points, new_output_path, number_contour_level, found_points)
    # reconstructionAlg.print_sector_pmts(PmtPositions, output_path)

    reco_accuracy.append(reconstructionAlg.reco_comparer(MC_positions, reco_positions))

    gc.collect()

reconstructionAlg.reco_resulter(reco_accuracy, output_path)
statusAlert.processStatus("Process finished")

