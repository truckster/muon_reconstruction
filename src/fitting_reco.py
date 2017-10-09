import statusAlert, recoPreparation, color_schemes
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mco

def reco_by_fittig_gauss(pmt_position_class, snippet_class, muon_points, out_path):

    statusAlert.processStatus("search patterns in snippet")

    '''photon data pre-processed and ordered into time snippets'''
    snippet_array = np.asarray(snippet_class.time_snippets)

    '''iterate snippets'''
    for snippet in range(50):
        statusAlert.processStatus("processing snippet: " + str(snippet))
        vertical_gauss_fit(pmt_position_class)

def vertical_gauss_fit(pmt_position_class):
    '''get data into numpy arrays'''
    id_np = np.asarray(pmt_position_class.id)
    phi_np = np.asarray(pmt_position_class.phi_position)
    theta_np = np.asarray(pmt_position_class.theta_position)

    print id_np
