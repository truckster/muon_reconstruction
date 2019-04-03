import recoPreparation, reconstructionAlg, statusAlert
import total_event_reconstruction, intersec_time_finder, point_allocate, diff_event_analysis, performance_check
import cPickle as pickle
from os import chdir, remove, path
from glob import glob
import math
import image_creator
import data
import muon_analysis
import backtracking_reco
import fht_analysis
import fht_fit


'''General script to use sub-scripts for muon reconstruction.'''
statusAlert.processStatus("Process started")

input_path = "/media/gpu/Data1/Simulation/processed_sim/run2/"
# input_path = "/media/gpu/Data1/Simulation/processed_sim/runs/run21/"
# input_path = "/media/gpu/Data1/Simulation/processed_sim/dev/"
# input_path = "/media/gpu/Data1/Simulation/processed_sim/runs/test3ns/"

# output_path = "/media/gpu/Data1/Analysis/Output/muReconstruction/total/"
output_path = "/media/gpu/Data1/Analysis/Output/muReconstruction/dev/"
# output_path = "/media/gpu/Data1/Analysis/Output/muReconstruction/dev2/"
# output_path = "/media/gpu/Data1/Analysis/Output/muReconstructiThe on/test/"

'''Overwrite fit results file'''
if path.isfile(output_path + "results.txt"):
    remove(output_path + "results.txt")

statusAlert.verbosity = 0

reco_performance = data.RecoPerformanceCheckTotal()

chdir(input_path)
folder_count = 0
folder_total = 0
for folder in glob("*/"):
    folder_total += 1
for folder in glob("*/"):
    statusAlert.processStatus("%.2f %% completed" % (float(folder_count)/float(folder_total)*100))
    folder_count += 1
    graph_output_path = output_path + "/graphs/"
    new_output_path = recoPreparation.create_output_path(graph_output_path, folder, "/totalEventHist/",  input_path)
    new_output_path_fit = recoPreparation.create_output_path(graph_output_path, folder, "/fits/",  input_path)
    total_path = recoPreparation.create_output_path(graph_output_path, folder, "/total_event/", input_path)
    fht_path = recoPreparation.create_output_path(graph_output_path, folder, "/fht/", input_path)
    ov_output_path = graph_output_path + "ov/"

    run_number = folder.split("-")[0]
    file_number = folder.split("-")[1]
    name = str(run_number) + "--" + file_number

    chdir(input_path + folder)
    statusAlert.processStatus("Reading file: " + str(folder))

    '''calculate PMT positions for this file'''
    PmtPositions = pickle.load(open("PMT_positions.pkl", 'rb'))

    '''collect entry and exit points of all muons in event'''
    intersec_radius = 19500
    # intersec_radius = 16800
    time_resolution = 1*10**-9
    merge_radius = 5
    frame_time_cut = 5
    number_contour_level = 8

    # muon_points = pickle.load(open("muon_points.pkl", 'rb'))
    muon_truth = pickle.load(open("muon_truth.pkl"))[0]
    mu_intersec_radius = 17700
    mu_time_resolution = 1*10**-9

    # muon_points = muon_analysis.calc_muon_detector_intersec_points(muon_truth.muon_data, mu_intersec_radius, mu_time_resolution)
    muon_points = muon_analysis.calc_muon_detector_intersec_points2(muon_truth.muon_data, mu_intersec_radius, mu_time_resolution)
    MC_positions, mc_muon_track = recoPreparation.MC_truth_writer(muon_points, output_path, folder, frame_time_cut)
    shower_status = muon_analysis.muon_is_showering_truth(muon_truth)

    '''collect information of all photons within certain time snippet and save the separately'''
    photons_in_time_window = pickle.load(open("framed_photons.pkl", 'rb'))
    photons_of_entire_event = pickle.load(open("total_event_photons.pkl", 'rb'))
    # photon_data = pickle.load(open("photon_data_array.pkl", 'rb'))

    # '''Load raw photon data'''
    # statusAlert.processStatus("Reading raw photon data...")
    # raw_photons = pickle.load(open("raw_photons.pkl", 'rb'))
    # statusAlert.processStatus("Raw photon data read!")

    '''FHT analysis'''
    fht_frames = fht_analysis.fht_reader(photons_of_entire_event.fht_array[0], frame_time_cut)
    # fht_analysis.fht_drawer(PmtPositions, fht_frames, muon_points, fht_path)
    # fht_analysis.draw_fht_picture_total(PmtPositions, photons_of_entire_event.fht_array, muon_points, fht_path,
    #                                     reco_points=None, mode="absolute")

    contour_array_total = pickle.load(open("Contours/contour_array_total.pkl"))
    contour_array_diff = pickle.load(open("Contours/contour_array_diff.pkl"))
    contour_array_total_shifted = pickle.load(open("Contours/contour_array_total_rotate.pkl"))
    contour_array_diff_shifted = pickle.load(open("Contours/contour_array_diff_rotate.pkl"))
    contour_entire_event = pickle.load(open("Contours/event_photons.pkl"))
    contour_entire_event_shifted = pickle.load(open("Contours/event_photons_shifted.pkl"))

    '''data sets'''
    total_event = [contour_entire_event, contour_entire_event_shifted]
    diff_event = [contour_array_diff, contour_array_diff_shifted]

    '''Reconstruction by looking at entire event'''
    found_points = total_event_reconstruction.entry_exit_detector(total_event)
    '''Reconstruction by looking at differences between two frames'''
    # found_points = diff_event_analysis.entry_exit_detector(diff_event)
    """Cross-check results with diffs"""
    # found_points = diff_event_analysis.intersec_crosscheck(diff_event, found_points)
    """Process reconstructed points"""
    processed_events = reconstructionAlg.process_events(found_points, merge_radius, intersec_radius,
                                                        reco_performance, PmtPositions)

    '''Detection of entry and exit time of muons'''
    # intersec_time_finder.find_times(contour_array_total, processed_events)
    intersec_time_finder.find_times_fht(fht_frames, processed_events)

    track_array = []
    if len(processed_events) > 3:
    # if len(processed_events) == 4:
        track_array = point_allocate.allocate_tracks_to_points(processed_events, reco_performance)
        reco_performance.point_reco_accuracy_list.append(performance_check.reco_comparer(MC_positions, processed_events))
        performance_check.reco_accuracy(track_array, mc_muon_track, reco_performance)
        performance_check.reco_accuracy_d_dep(track_array, mc_muon_track, reco_performance)

        if len(track_array) is 2:
            reco_performance.found_points_for_good_reco.append(len(processed_events))

        try:
            reco_performance.tracks_per_event_list.append(len(track_array))
        except:
            pass

        # fht_fit.do(fht_frames, track_array)

    intersec_time_finder.time_finder_performance(muon_points, track_array, frame_time_cut, reco_performance)

    reco_performance.found_point_list.append(len(processed_events))
    if shower_status:
        reco_performance.number_of_points_showering.append(len(processed_events))
    else:
        reco_performance.number_of_points_not_showering.append(len(processed_events))

    '''Draw all kinds of images'''
    # image_creator.snippet_drawer(PmtPositions, photons_of_entire_event, muon_points, total_path,
    #                                  number_contour_level, processed_events)
    # image_creator.snippet_drawer(PmtPositions, photons_in_time_window, muon_points, new_output_path,
    #                                  number_contour_level, processed_events)
    # image_creator.snippet_drawer_difference(PmtPositions, photons_in_time_window,
    #                                             muon_points, new_output_path, number_contour_level, processed_events)

    statusAlert.processStatus("--------------------- event finished ---------------------")
performance_check.do_performance_checks(reco_performance, output_path)


statusAlert.project_finished()

