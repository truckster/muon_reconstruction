import TreeReadFunc, statusAlert, PMTAnalysis, PointVecDist
import time, math, random, os
import ROOT as root
import numpy as np
import matplotlib.pyplot as plt
'''Scripts to prepare simulated data for muon reconstruction'''


def getKey(item):
    """important for time sorter"""
    return item[0]


class photon_snippet_output_c:
    def __init__(self):  # this method creates the class object.
        self.time_snippets = []

    def add_snippet_array(self, value):
        self.time_snippets.append(value)


def hitPMTinTimeSnippetHist2(source, wtime):
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
    returnClass = photon_snippet_output_c()
    startTime = data_sorted[0].hit_time
    index = 0

    while index < len(data_sorted)-1:
        PMTArray = [0] * 17739
        while data_sorted[index].hit_time < startTime + wtime:
            if data_sorted[index].pmt_id < 17739: # and datasort[index][2] is 1:
               PMTArray[data_sorted[index].pmt_id] += 1
            index += 1
        returnClass.add_snippet_array(PMTArray)
        startTime = data_sorted[index].hit_time
    statusAlert.processStatus("Done")
    return returnClass



class PhotonPositions:
    def __init__(self):  # this method creates the class object.
        self.hits = []

    def add_hits(self, value):
        self.hits.append(value)


def calcPMTPolarPhi(vector):
    phi = math.atan2(vector[2], vector[1])
    return phi


def calcPMTPolarTheta (vector):
    radius = PointVecDist.VectorLength(vector[1], vector[2], vector[3])
    theta = math.acos(-vector[3]/radius)
    # theta = math.atan(vector[3]/(math.sqrt(pow(vector[1], 2) + pow(vector[2], 2))))
    return theta


def muonEntryAndExitPoints(source):
    """Calculates and gives back the entry and exit points of the muons belonging to the event."""
    statusAlert.processStatus("Calculate entry and exit point of muons in this event from monte carlo truth")
    data = TreeReadFunc.readMuonRecoData(source)
    returnarray = []
    for i in range(len(data)):

        start = data[i][0]
        stop = data[i][1]
        direction = []
        for i in range(len(start)):
            length = math.sqrt(pow(stop[0]-start[0], 2) + pow(stop[1]-start[1], 2) + pow(stop[2]-start[2], 2))
            direction.append((stop[i]-start[i])/length)

        center = [0,0,0]
        oMinc = PointVecDist.VectorDifferenceVec(start, center)
        a = (PointVecDist.VectorTimesVec(direction, oMinc))
        b = pow(PointVecDist.VectorLengthVec(oMinc),2)-pow(19600,2)

        d1 = -a + math.sqrt(pow(a, 2) - b)
        d2 = -a - math.sqrt(pow(a, 2) - b)
        hitarray = []
        hit1 = []
        hit2 = []
        hit1.append("in")
        hit2.append("out")
        hit1 += PointVecDist.VectorSumVec(PointVecDist.VectorTimesValue(direction, d2), start)
        hit2 += PointVecDist.VectorSumVec(PointVecDist.VectorTimesValue(direction, d1), start)
        hitarray.append(hit1)
        hitarray.append(hit2)

        returnarray.append(hitarray)
    statusAlert.processStatus("Done")
    return returnarray


def create_output_path(outpath, file, extension,  inpath):
    try:
        os.chdir(outpath)
    except:
        os.makedirs(outpath + "/")
    os.chdir(outpath)
    print outpath
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
        self.phi_hammer_aitov = []
        self.theta_hammer_aitov = []
        self.sector_list_id = []
        self.sector_list_phi = []
        self.sector_list_theta = []
        self.x_sectors = 0
        self.y_sectors = 0

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

    def add_theta_position(self, value):
        self.theta_position.append(value)

    def add_phi_hammer_aitov(self, value):
        self.phi_hammer_aitov.append((value))

    def add_theta_hammer_aitov(self, value):
        self.theta_hammer_aitov.append((value))


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

    return_class.sector_list_id = sector_array_id
    return_class.sector_list_phi = sector_array_phi
    return_class.sector_list_theta = sector_array_theta
    return_class.x_sectors = x_sector_num
    return_class.y_sectors = y_sector_num

    return return_class

