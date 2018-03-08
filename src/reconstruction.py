import recoPreparation,reconstructionAlg, statusAlert, TreeReadFunc, gauss_fit_reco
from os import chdir, remove, path
from glob import glob


'''General script to use sub-scripts for muon reconstruction.'''
statusAlert.processStatus("Process started")

input_path = "/home/gpu/Simulation/mult/test/"
# input_path = "/home/gpu/Simulation/lele/2mult_xxl/"
# input_path = "/home/gpu/Simulation/single/"
output_path = "/home/gpu/Analysis/muReconstruction/Output/"

# input_path = "/home/gpu/Simulation/presentation/y/"
# output_path = "/home/gpu/Analysis/muReconstruction/Output/presentation/y/"

'''Overwrite fit results file'''
if path.isfile(output_path + "results.txt"):
    remove(output_path + "results.txt")

TreeReadFunc.check_file(input_path, "mu", "TrackLengthInScint", "MuMult")
chdir(input_path)
for file in glob("*.root"):
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

    '''collect information of all photons within certain time snippet and save the separately'''
    snippet_time_cut = 5
    photons_in_time_window = recoPreparation.hitPMTinTimeSnippetHist2(file, snippet_time_cut)

    '''draw pictures of all snippets'''
    # PMTAnalysis.drawPMTHitMapHisto(
    # photons_in_time_window.time_snippets, muon_points, PmtPositions, new_output_path, 50
    # )

    '''create output file'''
    result_file = open(output_path + "results.txt", 'a')
    '''write header and real muon points'''
    result_file.write("File: " + str(file)+'\n' + "MC truth"+'\n')
    for event in muon_points:
        result_file.write("Event: " + str(event.event) + '\n')
        if event.enters is True:
            result_file.write("Entry point: " + '\n')
        if event.leaves is True:
            result_file.write("Exit point: " + '\n')
        result_file.write("Phi: " + str(event.phi) + '\n')
        result_file.write("Theta: " + str(event.theta) + '\n')
        result_file.write("Z: " + str(event.z) + '\n')
        result_file.write("Time: " + str(event.intersec_time) + '\n')
        result_file.write("------------------------------------------" + '\n'+ '\n')
        # result_file.write("Calculated snippet: " + str(event.intersec_time//(snippet_time_cut*10**-9)-12) + '\n')
        # print((event.intersec_time/time_resolution)//snippet_time_cut)
    result_file.write("----- Reconstructed Values ------" + '\n')

    '''Take data from 'snippets' for reconstruction: find all patches within one time snippet'''
    reconstructionAlg.pattern_detector(PmtPositions, photons_in_time_window, muon_points, new_output_path)
    reconstructionAlg.pattern_detector_difference(PmtPositions, photons_in_time_window, muon_points, new_output_path)
    # gauss_fit_reco.fit_function_caller(PmtPositions, photons_in_time_window, muon_points, new_output_path_fit, result_file)
    # reconstructionAlg.print_sector_pmts(PmtPositions, output_path)

    result_file.close()

statusAlert.processStatus("Process finished")
