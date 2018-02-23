import recoPreparation, PMTAnalysis, reconstructionAlg, statusAlert, TreeReadFunc, fitting_reco, gauss_fit_reco
from os import chdir, remove
from glob import glob

statusAlert.processStatus("Process started")

input_path = "/home/gpu/Simulation/mult/new/"

chdir(input_path)
for file in glob("*.root"):
    statusAlert.processStatus("Reading file: " + str(file))
    '''collect entry and exit points of all muons in event'''
    muon_points = recoPreparation.calc_muon_detector_intersec_points(file, 17600, (1*10**-9))

