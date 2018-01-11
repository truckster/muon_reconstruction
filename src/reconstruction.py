import recoPreparation, PMTAnalysis, reconstructionAlg, statusAlert, TreeReadFunc, fitting_reco
from os import chdir
from glob import glob


'''General script to use sub-scripts for muon reconstruction.'''
statusAlert.processStatus("Process started")

input_path = "/home/gpu/Simulation/mult/test/"
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

    new_output_path = recoPreparation.create_output_path(output_path, file, "/totalEventHist/",  input_path)

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

    '''Take data from 'snippets' for reconstruction: find all patches within one time snippet'''
    cut_radius = 0.15
    sector_threshold = 0.0
    cut_threshold = 0.4
    # reconstructionAlg.pattern_detector(PmtPositions, photons_in_time_window, muon_points, new_output_path, cut_radius,
    #                                    sector_threshold, cut_threshold)
    fitting_reco.reco_by_fittig_gauss(PmtPositions, photons_in_time_window, muon_points, output_path)
    # reconstructionAlg.print_sector_pmts(PmtPositions, output_path)

statusAlert.processStatus("Process finished")