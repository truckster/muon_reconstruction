import recoPreparation, reconstructionAlg, statusAlert, TreeReadFunc, gauss_fit_reco, contour_analyze
import total_event_reconstruction, intersec_time_finder

from os import chdir, remove, path
from glob import glob
import gc

'''General script to use sub-scripts for muon reconstruction.'''
statusAlert.processStatus("Process started")

input_path = "/home/gpu/Simulation/mult/new/"
# input_path = "/home/gpu/Simulation/mult/test/"
# input_path = "/home/gpu/Simulation/test/"
# input_path = "/home/gpu/Simulation/single/"
output_path = "/home/gpu/Analysis/muReconstruction/Output/"
# output_path = "/home/gpu/Analysis/muReconstruction/Output/LPMT/"

# input_path = "/home/gpu/Simulation/presentation/y/"
# output_path = "/home/gpu/Analysis/muReconstruction/Output/presentation/y/"

'''Overwrite fit results file'''
if path.isfile(output_path + "results.txt"):
    remove(output_path + "results.txt")

TreeReadFunc.check_file(input_path, "mu", "TrackLengthInScint", "MuMult")
chdir(input_path)
for file in glob("*.root"):
    # TODO start new process for each file. This might solve the memory problem.
    statusAlert.processStatus("Reading file: " + str(file))

    '''Control, which events are useful for the analysis'''
    # TreeReadFunc.interestChecker(inputpath, "mu", "TrackLengthInScint", "MuMult", "MuMult")
    # TreeReadFunc.interestChecker(inputpath, "mu", "MuMult"
    # TreeReadFunc.muon_file_reader(file)

    new_output_path = recoPreparation.create_output_path(output_path, file, "/totalEventHist/",  input_path)
    new_output_path_fit = recoPreparation.create_output_path(output_path, file, "/fits/",  input_path)

    '''calculate PMT positions for this file'''
    x_sectors = 20
    y_sectors = 10
    PmtPositions = recoPreparation.calc_pmt_positions(input_path, x_sectors, y_sectors)

    '''collect entry and exit points of all muons in event'''
    intersec_radius = 17600
    time_resolution = 1*10**-9
    muon_points = recoPreparation.calc_muon_detector_intersec_points(file, intersec_radius, time_resolution)

    recoPreparation.MC_truth_writer(muon_points, output_path, file)

    '''collect information of all photons within certain time snippet and save the separately'''
    snippet_time_cut = 10
    max_snips = 30
    photons_in_time_window, photons_of_entire_event = recoPreparation.hitPMTinTimeSnippetHist2(file,
                                                                                               snippet_time_cut,
                                                                                               max_snips)

    '''Reconstruction by looking at entire event'''
    total_path = recoPreparation.create_output_path(output_path, file, "/total_event/", input_path)
    # # reconstructionAlg.snippet_drawer(PmtPositions, photons_of_entire_event, muon_points, total_path)
    found_points = total_event_reconstruction.entry_exit_detector(PmtPositions, photons_of_entire_event)
    # total_event_reconstruction.reco_result_writer(output_path, found_points)

    '''Detection of entry and exit time of muons'''
    intersec_time_finder.find_times(PmtPositions, photons_in_time_window, found_points)

    '''Draw all kinds of images'''
    # reconstructionAlg.snippet_drawer(PmtPositions, photons_in_time_window, muon_points, new_output_path)
    # reconstructionAlg.snippet_drawer_difference(PmtPositions, photons_in_time_window, muon_points, new_output_path, max_snips)
    # reconstructionAlg.print_sector_pmts(PmtPositions, output_path)

    gc.collect()

statusAlert.processStatus("Process finished")
