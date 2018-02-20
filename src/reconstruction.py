import recoPreparation, PMTAnalysis, reconstructionAlg, statusAlert, TreeReadFunc, fitting_reco, gauss_fit_reco
from os import chdir
from glob import glob


'''General script to use sub-scripts for muon reconstruction.'''
statusAlert.processStatus("Process started")

input_path = "/home/gpu/Simulation/mult/new/"
# input_path = "/home/gpu/Simulation/single/"
output_path = "/home/gpu/Analysis/muReconstruction/Output/"

# input_path = "/home/gpu/Simulation/presentation/y/"
# output_path = "/home/gpu/Analysis/muReconstruction/Output/presentation/y/"

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
    muon_points = recoPreparation.muonEntryAndExitPoints(file)

    '''collect information of all photons within certain time snippet and save the separately'''
    snippet_time_cut = 5
    photons_in_time_window = recoPreparation.hitPMTinTimeSnippetHist2(file, snippet_time_cut)

    '''draw pictures of all snippets'''
    # PMTAnalysis.drawPMTHitMapHisto(
    # photons_in_time_window.time_snippets, muon_points, PmtPositions, new_output_path, 50
    # )

    '''create output file'''
    result_file = open(output_path + "results.txt", 'w')
    '''write header and real muon points'''
    result_file.write("File: " + str(file)+'\n')
    muon_sphere_points = gauss_fit_reco.calc_muon_points_sphere(muon_points)
    for coordinate in muon_sphere_points:
        result_file.write(str(coordinate)+'\n')

    '''Take data from 'snippets' for reconstruction: find all patches within one time snippet'''
    # reconstructionAlg.pattern_detector(PmtPositions, photons_in_time_window, muon_points, new_output_path)
    # fitting_reco.reco_by_fittig_gauss(PmtPositions, photons_in_time_window, muon_points, output_path)
    gauss_fit_reco.fit_function_caller(PmtPositions, photons_in_time_window, muon_points, new_output_path_fit, result_file)
    # reconstructionAlg.print_sector_pmts(PmtPositions, output_path)

    result_file.close()

statusAlert.processStatus("Process finished")
