import TreeReadFunc, statusAlert, PointVecDist
import time, math, random, os
import ROOT as root
import numpy as np
import matplotlib.pyplot as plt
'''Scripts to prepare simulated data for muon reconstruction'''


def getKey(item):
    """important for time sorter"""
    return item[0]


class PhotonPMTHits:
    def __init__(self):  # this method creates the class object.
        self.time_snippets = []

    def add_pmt_array(self, value):
        self.time_snippets.append(value)


def hitPMTinTimeSnippetHist2(source, wtime, snippet_limit):
    """reads time and ID of each event in the given chain. wtime defines the time window for each package of condensed
        events Then it is ordered by hit time and written into arrays:time, X, Y, Z. Furthermore the number of hits per
        PMT and time windows is counted.
        """
    statusAlert.processStatus("Start snippet read for time: " + str(wtime))
    pmt_time_resolution = 3.0
    data_read = TreeReadFunc.readPhotonRecoData2(source, pmt_time_resolution)
    statusAlert.processStatus("Sort data")
    data_sorted = sorted(data_read, key=lambda x: x.hit_time, reverse=False)

    statusAlert.processStatus("Divide arrays")
    return_class_snippets = PhotonPMTHits()
    return_class_total = PhotonPMTHits()
    start_time = data_sorted[0].hit_time
    photon_index = 0
    snippet_index = 0

    pmt_array_total = [0]*17739
    while photon_index < len(data_sorted)-1 and snippet_index < snippet_limit:
        pmt_array_snippet = [0] * 17739
        hit_sum_snippet = 0
        while data_sorted[photon_index].hit_time < start_time + wtime and photon_index < len(data_sorted)-1:
            if data_sorted[photon_index].pmt_id < 17739:
                pmt_array_snippet[data_sorted[photon_index].pmt_id] += 1
                hit_sum_snippet += 1
                pmt_array_total[data_sorted[photon_index].pmt_id] += 1
            photon_index += 1
        snippet_index += 1
        if hit_sum_snippet > 0:
            return_class_snippets.add_pmt_array(pmt_array_snippet)

        start_time = data_sorted[photon_index].hit_time
    return_class_total.add_pmt_array(pmt_array_total)
    statusAlert.processStatus("Done")

    return return_class_snippets, return_class_total


class MuonIntersecPoint:
    def __init__(self):
        self.event = 0

        self.enters = False
        self.leaves = False

        self.x = 0
        self.y = 0
        self.z = 0

        self.phi = 0
        self.theta = 0
        self.theta2 = 0

        self.phi_hammer_aitoff = 0
        self.theta_hammer_aitoff = 0

        self.intersec_time = 0


def distance_to_center(x, y, z):
    dist = math.sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))
    return dist


def is_muon_in_detector(dist_to_center, radius):
    if dist_to_center < radius:
        muon_in_det = True
    else:
        muon_in_det = False

    return muon_in_det


def calc_muon_detector_intersec_points(source, radius, time_steps):
    """Calculates and gives back the entry and exit points of the muons belonging to the event."""
    statusAlert.processStatus("Calculate entry and exit point of muons in this event from monte carlo truth")
    muon_data = TreeReadFunc.read_muon_data(source)
    returnarray = []

    muon_mass = 105.6583745
    c = 299792548000

    for index, muon_event in enumerate(muon_data):
        muon_in = MuonIntersecPoint()
        muon_out = MuonIntersecPoint()

        muon_in.event = index
        muon_out.event = index

        time_current = 0

        total_momentum = math.sqrt(pow(muon_event.x_momentum_init, 2)
                                   + pow(muon_event.y_momentum_init, 2)
                                   + pow(muon_event.z_momentum_init, 2))

        x_velocity = muon_event.x_momentum_init / math.sqrt(pow(muon_mass, 2) + total_momentum**2) * c
        y_velocity = muon_event.y_momentum_init / math.sqrt(pow(muon_mass, 2) + total_momentum**2) * c
        z_velocity = muon_event.z_momentum_init / math.sqrt(pow(muon_mass, 2) + total_momentum**2) * c

        x_position_current = muon_event.x_position_init
        y_position_current = muon_event.y_position_init
        z_position_current = muon_event.z_position_init

        distance_to_det_center = distance_to_center(x_position_current, y_position_current, z_position_current)

        while z_position_current > -23000:
            muon_in_detector = is_muon_in_detector(distance_to_det_center, radius)
            x_position_current = x_position_current + time_steps * x_velocity
            y_position_current = y_position_current + time_steps * y_velocity
            z_position_current = z_position_current + time_steps * z_velocity

            time_current += time_steps
            distance_to_det_center = distance_to_center(x_position_current, y_position_current, z_position_current)
            muon_in_detector_second = is_muon_in_detector(distance_to_det_center, radius)

            if muon_in_detector is False and muon_in_detector_second is True:
                muon_in.enters = True
                muon_in.x = x_position_current
                muon_in.y = y_position_current
                muon_in.z = z_position_current

                muon_in.phi = math.atan2(y_position_current, x_position_current)
                muon_in.theta = math.acos(-z_position_current/radius)
                muon_in.theta2 = math.acos(-z_position_current/radius) - (math.pi/2)

                muon_in.phi_hammer_aitoff = (math.sqrt(8) * math.cos(muon_out.theta) * math.sin(muon_out.phi / 2)) / \
                                             (math.sqrt(1 + math.cos(muon_out.theta) * math.cos(muon_out.phi / 2)))

                muon_in.theta_hammer_aitoff = (math.sqrt(2) * math.sin(muon_out.theta)) / \
                                               (math.sqrt(1 + math.cos(muon_out.theta) * math.cos(muon_out.phi / 2)))

                muon_in.intersec_time = time_current

                returnarray.append(muon_in)

            if muon_in_detector is True and muon_in_detector_second is False:
                muon_out.leaves = True
                muon_out.x = x_position_current
                muon_out.y = y_position_current
                muon_out.z = z_position_current

                muon_out.phi = math.atan2(y_position_current, x_position_current)
                muon_out.theta = math.acos(-z_position_current/radius)
                muon_out.theta2 = math.acos(-z_position_current/radius) - (math.pi/2)

                muon_out.phi_hammer_aitoff = (math.sqrt(8) * math.cos(muon_out.theta) * math.sin(muon_out.phi / 2)) / \
                          (math.sqrt(1 + math.cos(muon_out.theta) * math.cos(muon_out.phi / 2)))

                muon_out.theta_hammer_aitoff = (math.sqrt(2) * math.sin(muon_out.theta)) / \
                            (math.sqrt(1 + math.cos(muon_out.theta) * math.cos(muon_out.phi / 2)))

                muon_out.intersec_time = time_current

                returnarray.append(muon_out)
    return returnarray


def create_output_path(outpath, file, extension,  inpath):
    try:
        os.chdir(outpath)
    except:
        os.makedirs(outpath + "/")
    os.chdir(outpath)
    output_path2 = outpath + str(file) + extension
    if os.path.isdir(output_path2) is False:
        statusAlert.processStatus("Create output path")
        os.makedirs(output_path2)
        os.chdir(output_path2)
        # os.makedirs(output_path2 + "totalEventHist/")
        # os.chdir(output_path2 + "totalEventHist/")
        statusAlert.processStatus("Done")
    os.chdir(inpath)

    return output_path2


class PMTPositions:
    def __init__(self):  # this method creates the class object.
        self.id = []
        self.x_position = []
        self.y_position = []
        self.z_position = []
        self.phi_position = []
        self.theta_position = []
        self.theta_position2 = []
        self.phi_hammer_aitov = []
        self.theta_hammer_aitov = []
        self.sector_list_id = []
        self.sector_list_phi = []
        self.sector_list_theta = []
        self.x_sectors = 0
        self.y_sectors = 0
        self.is_x_sector = []
        self.is_y_sector = []

    def add_id(self, value):
        self.id.append(value)

    def add_x_position(self, value):
        self.x_position.append(value)

    def add_y_position(self, value):
        self.y_position.append(value)

    def add_z_position(self, value):
        self.z_position.append(value)

    def add_phi_position(self, value):
        self.phi_position.append(value)

    """0deg to 180deg"""
    def add_theta_position(self, value):
        self.theta_position.append(value)

    """-90deg to 90deg"""
    def add_theta_position2(self, value):
        self.theta_position2.append(value)

    def add_phi_hammer_aitov(self, value):
        self.phi_hammer_aitov.append((value))

    def add_theta_hammer_aitov(self, value):
        self.theta_hammer_aitov.append((value))

    def add_x_sector(self, value):
        self.is_x_sector.append(value)

    def add_y_sector(self, value):
        self.is_y_sector.append(value)


def calc_pmt_positions(inpath, x_sector_num, y_sector_num):
    """Calculate xyz coordinates of PMTs from ID"""
    statusAlert.processStatus("Get pmt positions in x, y and z and calculate phi and theta positions")
    pmt_file = root.TFile.Open(os.listdir(inpath)[0])
    pmt_tree = pmt_file.Get("pmtpos")

    sector_num = x_sector_num * y_sector_num

    sector_array_id = [[] for _ in range(sector_num)]
    sector_array_phi = [[] for _ in range(sector_num)]
    sector_array_theta = [[] for _ in range(sector_num)]

    return_class = PMTPositions()
    for event in pmt_tree:
        return_class.add_id(event.__getattr__("pmtID"))
        return_class.add_x_position(event.__getattr__("x"))
        return_class.add_y_position(event.__getattr__("y"))
        return_class.add_z_position(event.__getattr__("z"))

        phi_position = math.atan2(event.__getattr__("y"), event.__getattr__("x"))
        return_class.add_phi_position(phi_position)

        radius = PointVecDist.VectorLength(event.__getattr__("x"), event.__getattr__("y"), event.__getattr__("z"))
        theta_position = math.acos(- event.__getattr__("z")/radius)
        return_class.add_theta_position(theta_position)
        return_class.add_theta_position2(theta_position - (math.pi/2))

        x_sector = int((phi_position + math.pi) / 2.01 / math.pi * x_sector_num)
        y_sector = int(theta_position/math.pi * y_sector_num)
        sector = y_sector*x_sector_num + x_sector

        sector_array_id[sector].append(event.__getattr__("pmtID"))
        sector_array_phi[sector].append(phi_position)
        sector_array_theta[sector].append(theta_position)

        theta_position -= math.pi / 2
        phi_h_a = (math.sqrt(8) * math.cos(theta_position) * math.sin(phi_position / 2)) / \
                  (math.sqrt(1 + math.cos(theta_position) * math.cos(phi_position / 2)))
        return_class.add_phi_hammer_aitov(phi_h_a)

        theta_h_a = (math.sqrt(2) * math.sin(theta_position)) / \
                    (math.sqrt(1 + math.cos(theta_position) * math.cos(phi_position / 2)))
        return_class.add_theta_hammer_aitov(theta_h_a)

        return_class.add_x_sector(x_sector)
        return_class.add_y_sector(y_sector)

    return_class.sector_list_id = sector_array_id
    return_class.sector_list_phi = sector_array_phi
    return_class.sector_list_theta = sector_array_theta
    return_class.x_sectors = x_sector_num
    return_class.y_sectors = y_sector_num

    return return_class


def MC_truth_writer(muon_points, output_path, file, time_cut):
    '''create output file'''
    result_file = ""
    '''write header and real muon points'''
    result_file += ("::::::::::::::::::::::::::::::::"+'\n')
    result_file += ("File: " + str(file)+'\n')
    result_file += ("::::::::::::::::::::::::::::::::"+'\n')
    result_file += ("---------MC truth---------"+'\n')
    times = []
    for event in muon_points:
        result_file += ("Event: " + str(event.event) + '\n')
        if event.enters is True:
            result_file += ("Entry point: " + '\n')
        if event.leaves is True:
            result_file += ("Exit point: " + '\n')
        result_file += ("Phi: " + str(event.phi) + '\n')
        result_file += ("Theta: " + str(event.theta2) + '\n')
        result_file += ("Z: " + str(event.z) + '\n')
        result_file += ("Time: " + str(event.intersec_time) + '\n')
        result_file += ("------------------------------------------" + '\n')
        times.append(event.intersec_time)

    first_intersec_time = min(times)
    frames = []
    result_file += "---------Frames of intersection (MC)---------"+'\n'
    for i in range(len(times)):
        # frames.append((times[i]-first_intersec_time)/time_cut)
        result_file += ("Frame: " + str((times[i]-first_intersec_time)/time_cut*pow(10, 9)) + '\n')
    r_file = open(output_path + "results.txt", 'a')
    r_file.write(result_file)
    r_file.close()
