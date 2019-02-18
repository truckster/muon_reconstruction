'''Read data from files and put them to root.chains. Here all data is read out from files. Later the respectively data
    is put into arrays. Maybe the chain-step can be skipped. Zombie check finds and deletes corrupted files from simulation'''
TreeReadFunc.isZombie(inputpath)
evtChainMu = TreeReadFunc.readMultRootFiles(inputpath, "mu", "mu", "TrackLengthInScint", "MuMult")
evtChainEvt = TreeReadFunc.readMultRootFiles(inputpath, "evt", "mu", "TrackLengthInScint", "MuMult")
evtChainGen = TreeReadFunc.readMultRootFiles(inputpath, "geninfo", "mu", "TrackLengthInScint", "MuMult")


def hitTimePMTs(source, limit):
    """reads time and ID of each event in the given chain. Limit defines the number of events to be read. Then it is ordered by hit time and written into arrays:time, X, Y, Z"""
    statusAlert.processStatus("Start muon printer")
    statusAlert.processStatus("Start")

    data = TreeReadFunc.readPhotonRecoData(source)
    datasort = sorted(data, key=getKey)

    statusAlert.processStatus("data read")
    PMTPos = PMTAnalysis.XYZfromID()
    array = []
    for i in range(limit):
        if datasort[i][1] < 17739:
            subarray = []
            subarray.append(datasort[i][0])
            subarray.append(PMTAnalysis.getXYZfromID(PMTPos, datasort[i][1])[0])
            subarray.append(PMTAnalysis.getXYZfromID(PMTPos, datasort[i][1])[1])
            subarray.append(PMTAnalysis.getXYZfromID(PMTPos, datasort[i][1])[2])
            array.append(subarray)


    statusAlert.processStatus("Finished muon printer")
    return array


def hitPMTinTimeSnippet(source, wtime):
    """reads time and ID of each event in the given chain. wtime defines the time window for each package of condensed
    events Then it is ordered by hit time and written into arrays:time, X, Y, Z
    """
    statusAlert.processStatus("Start snippet read for time: " + str(wtime))
    data = TreeReadFunc.readPhotonRecoData(source)
    datasort = sorted(data, key=getKey)

    statusAlert.processStatus("data read")
    PMTPos = PMTAnalysis.XYZfromID()

    statusAlert.processStatus("Divide arrays")
    returnArray = []
    startTime = datasort[0][0]
    index = 0

    while index < len(datasort)-1:
        array = []
        while datasort[index][0] < startTime + wtime:
            subarray = []
            if datasort[index][1] < 17739:
                # insert time smearing here
                subarray.append(datasort[index][0])
                subarray.append(PMTAnalysis.getXYZfromID(PMTPos, datasort[index][1])[0])
                subarray.append(PMTAnalysis.getXYZfromID(PMTPos, datasort[index][1])[1])
                subarray.append(PMTAnalysis.getXYZfromID(PMTPos, datasort[index][1])[2])
                array.append(subarray)
            index += 1
        startTime = datasort[index][0]
        returnArray.append(array)

    return returnArray



def XYZfromID():
    """Calculate xyz coordinates of PMTs from ID"""
    PMTFile = root.TFile.Open("/local/scratch0/Analysis/muReconstruction/Input/PMTPos.root")
    PMTTree = PMTFile.Get("pmtpos")
    coord = []
    for event in PMTTree:
        vec = []
        # vec.append(event.__getattr__("pmtID"))
        vec.append(event.__getattr__("x"))
        vec.append(event.__getattr__("y"))
        vec.append(event.__getattr__("z"))
        coord.append(vec)
    return coord


def getXYZfromID(array, id):
    return array[id]


def getPMTidxyz():
    PMTFile = root.TFile.Open("/local/scratch0/Analysis/muReconstruction/Input/PMTPos.root")
    PMTTree = PMTFile.Get("pmtpos")
    coord = []
    for event in PMTTree:
        vec = []
        vec.append(event.__getattr__("pmtID"))
        vec.append(event.__getattr__("x"))
        vec.append(event.__getattr__("y"))
        vec.append(event.__getattr__("z"))
        coord.append(vec)
    return coord


def drawPMTMapWithHits(xyzHitArray, outPath, name = None):
    fig = plt.figure(num=None, figsize=(20, 10))
    ax1 = fig.add_subplot(111)

    xyz = PMTAnalysis.getPMTidxyz()
    phi = []
    theta = []
    for i in range(len(xyz)):
        phi.append(calcPMTPolarPhi(xyz[i]))
        theta.append(calcPMTPolarTheta(xyz[i]))

    phiHit = []
    thetaHit = []
    for i in range(len(xyzHitArray)):
        phiHit.append(calcPMTPolarPhi(xyzHitArray[i]))
        thetaHit.append(calcPMTPolarTheta(xyzHitArray[i]))

    ax1.scatter(phi, theta, c='w', marker=',', label='not hit')
    ax1.scatter(phiHit, thetaHit, c='r', label='hit')
    if name is None:
        plt.savefig(outPath+"PMTMapWithHits.png")
    else:
        plt.savefig(outPath + "totalEvent/" + str(name) + ".png")

    plt.close()


def drawPMTs(outPath):
    statusAlert.processStatus("Start printing empty PMT positions")
    xyz = PMTAnalysis.getPMTidxyz()
    phi = []
    theta = []
    for i in range(len(xyz)):
        phi.append(calcPMTPolarPhi(xyz[i]))
        theta.append(calcPMTPolarTheta(xyz[i]))

    plt.figure(num=None, figsize=(20, 10))
    plt.scatter(phi, theta)
    plt.savefig(outPath+"PMTmap.png")
    statusAlert.processStatus("Finished printing empty PMT positions")


def toXYZConverter(inputarray):
    """Calculates the position of an event depending on the PMT id it hits"""
    PMTPos = PMTAnalysis.XYZfromID()
    returnarray = []
    for i in range(len(inputarray)):
        snippetArray = []
        for j in range(len(inputarray[i])):
            subarray = []
            subarray.append(inputarray[i][j])
            subarray.append(PMTAnalysis.getXYZfromID(PMTPos, j)[0])
            subarray.append(PMTAnalysis.getXYZfromID(PMTPos, j)[1])
            subarray.append(PMTAnalysis.getXYZfromID(PMTPos, j)[2])
            snippetArray.append(subarray)
        returnarray.append(snippetArray)
    statusAlert.processStatus("Converted")
    return returnarray


def hitPMTinTimeSnippetHist(source, wtime):
    """reads time and ID of each event in the given chain. wtime defines the time window for each package of condensed
        events Then it is ordered by hit time and written into arrays:time, X, Y, Z. Furthermore the number of hits per
        PMT and time windows is counted.
        """
    statusAlert.processStatus("Start snippet read for time: " + str(wtime))
    data = TreeReadFunc.readPhotonRecoData(source, 3.0)
    datasort = sorted(data, key=getKey)

    statusAlert.processStatus("data read")

    statusAlert.processStatus("Divide arrays")
    returnClass = photon_snippet_output_c()
    startTime = datasort[0][0]
    index = 0

    while index < len(datasort)-1:
        PMTArray = [0] * 17739
        while datasort[index][0] < startTime + wtime:
            if datasort[index][1] < 17739: # and datasort[index][2] is 1:
                PMTArray[datasort[index][1]] += 1
            index += 1
        returnClass.add_snippet_array(PMTArray)
        startTime = datasort[index][0]
    return returnClass


def to_phi_theta_converter(snippet_array):
    statusAlert.processStatus("Get angular position of photon hits")
    return_array = np.array([])
    for i in range(len(snippet_array)):
        PhotonPosAndHits = PhotonPositions()
        for j in range(len(snippet_array[i])):
            PhotonPosAndHits.add_hits(snippet_array[i])
            return_array = np.append(return_array, PhotonPosAndHits)
    statusAlert.processStatus("Done")
    return return_array


def number_photons_per_pmt(snippets):
    for i in range(len(snippets)):
        print "Snippet: " + str(i)
        print "Max PMT: " + str(np.argmax(snippets[i]))
        print "Hits: " + str(max(snippets[i]))
        print "______________________________"



def readMultRootFiles(path, tree, ctree=None, cbranch=None, citerator=None):
    """read multiple root files and combine them to a chain. optional only interesting files are added"""
    if ctree is None:
        chain = TChain(tree)
        chdir(path)
        for file in glob("*.root"):
            # statusAlert.processStatus("Reading file")
            chain.Add(file)
        return chain

    else:
        chain = TChain(tree)
        chdir(path)
        countTot = 0
        countInteresting = 0
        for file in glob("*.root"):
            countTot += 1
            if isInteresting(file, ctree, cbranch, citerator):
                countInteresting += 1
                # statusAlert.processStatus("Reading file " + str(countTot))
                chain.Add(file)
                print file + " is interesting!"
            # else:
                # print file + " is boring!"
        print str(countInteresting) + " of " + str(countTot) + " files are relevant for the analysis."
        return chain


def readPhotonRecoData(source, pmt_resolution):
    statusAlert.processStatus("     Collect photon data from file")
    array = []
    count = 0
    rfile = TFile(source)
    rtree = rfile.Get("evt")

    for event in rtree:
        for i in range(event.__getattr__("nPhotons")):
            # statusAlert.processStatus("Reading muon number " + str(count) + " and photon number " + str(i))
            subarray = []
            # print "original time: " + str(event.__getattr__("hitTime")[i])
            # print "smeared time: " + str(vetoPreparation.time_smearing(event.__getattr__("hitTime")[i], 12.0))
            subarray.append(time_smearing(event.__getattr__("hitTime")[i], pmt_resolution))
            subarray.append(event.__getattr__("pmtID")[i])
            subarray.append(event.__getattr__("isCerenkov")[i])
            array.append(subarray)
    statusAlert.processStatus("     Done")
    return array


bin_centers = bins[:-1] + 0.5 * (bins[1:] - bins[:-1])

        try:
            popt, pcov = curve_fit(gauss, bin_centers, n, p0=[1.0, 1.0, 1.0])

            plt.plot(bin_centers, gauss(bin_centers, *popt), 'r--')

            print popt
        #
        # fit_line = plt.plot(unsorted_fit_data.x_value, gauss(unsorted_fit_data.x_value, *popt))
        #
        # plt.plot(100, fit_line, 'r--', linewidth=2)
        except:
            print("Sector: " + str(vertical_sector) + ": Fit did not work!")
        plt.plot()

        plt.savefig(out_path + "test/vertical/" + str(vertical_sector) + ".png")
        plt.close()



def pattern_scan(pmt_position_class, snippet_array, snippet, template_radius, sector_threshold, cut_threshold):
    '''this array saves all patterns within one snippet'''
    sector_pattern_array = []
    '''iterate sectors'''
    for sector in range(len(pmt_position_class.sector_list_id)):

        '''get data into numpy arrays'''
        id_np = np.asarray(pmt_position_class.sector_list_id[sector])
        phi_np = np.asarray(pmt_position_class.sector_list_phi[sector])
        theta_np = np.asarray(pmt_position_class.sector_list_theta[sector])

        total_hit_sum = sum(snippet_array[snippet])

        '''initialize class, which savespattern position. This class is given to the sector_pattern_array'''
        pattern_position = PatternPosition()

        '''get sum of photon hits in each sector'''
        sector_hit_sum = 0
        for i in range(len(id_np)):
            sector_hit_sum += snippet_array[snippet][id_np[i]]

        if float(sector_hit_sum) > sector_threshold * total_hit_sum:
            '''scan the template over the sector'''
            cut_phi_position = min(phi_np)
            while cut_phi_position < max(phi_np) + template_radius:
                cut_theta_position = min(theta_np)
                while cut_theta_position < max(theta_np) + template_radius:
                    '''sum of photons within the template'''
                    hit_sum = 0
                    '''scan over all pmts in the sector and check if they are inside the template cut'''
                    # faster to get pmts directly instead of iterating all in sector
                    for i in range(len(id_np)):
                        if (pow(template_radius, 2) > (pow(phi_np[i] - cut_phi_position, 2)
                                                           + pow(theta_np[i] - cut_theta_position, 2))):
                            hit_sum += snippet_array[snippet][id_np[i]]
                    '''check if photon sum is over threshold and if its the "biggest" pattern in the sector.
                     Then overwrite previous data'''

                    if hit_sum > cut_threshold * sector_hit_sum and hit_sum > pattern_position.hits:
                        pattern_position.hits = hit_sum
                        pattern_position.phi = cut_phi_position
                        pattern_position.theta = cut_theta_position

                    cut_theta_position += template_radius
                cut_phi_position += template_radius

        sector_pattern_array.append(pattern_position)
    return sector_pattern_array


def muon_entry_exit_reconstructor(pmt_position_class, snippet_class, out_path):
    statusAlert.processStatus("Searching for entry and exit points")

    '''photon data pre-processed and ordered into time snippets'''
    snippet_array = np.asarray(snippet_class.time_snippets)

    '''sector count data is saved here for each snippet'''
    snippet_sector_array = [[] for _ in range(len(snippet_array))]
    for snippet in range(len(snippet_array)):
        sector_array = [[] for _ in range(len(pmt_position_class.sector_list_id))]
        for sector in range(len(pmt_position_class.sector_list_id)):
            id_np = np.asarray(pmt_position_class.sector_list_id[sector])
            sector_hit_sum = 0
            for i in range(len(id_np)):
                sector_hit_sum += snippet_array[snippet][id_np[i]]
            sector_array[sector] = (sector_hit_sum)
        snippet_sector_array[snippet] = (sector_array)

    '''iterate over snippet, sector and the pmts of the sector to get its hit_sum'''
    snippet_sector_pattern_array = []
    for snippet in range(len(snippet_array)):

        sector_pattern_array = []
        for sector in range(len(pmt_position_class.sector_list_id)):
            pattern_position = PatternPosition()

            id_np = np.asarray(pmt_position_class.sector_list_id[sector])
            phi_np = np.asarray(pmt_position_class.sector_list_phi[sector])
            theta_np = np.asarray(pmt_position_class.sector_list_theta[sector])

            dist_sum = 0.0
            phi_sum = 0.0
            theta_sum = 0.0
            counter = 0.0
            sector_hit_sum = 0

            for i in range(len(id_np)):
                sector_hit_sum += snippet_array[snippet][id_np[i]]
                if snippet_array[snippet][id_np[i]] > 10:
                    dist_sum += math.sqrt(pow(phi_np[i]-phi_np[0], 2) + pow(theta_np[i]-theta_np[0], 2))
                    phi_sum += phi_np[i]
                    theta_sum += theta_np[i]
                    counter += 1.0
            if snippet > 0:
                if 4 > counter > 0 and dist_sum/counter < 0.1 and \
                                snippet_sector_array[snippet][sector] > 50 * snippet_sector_array[snippet-1][sector]:
                    pattern_position.phi = phi_sum / counter
                    pattern_position.theta = theta_sum / counter

            sector_pattern_array.append(pattern_position)

        snippet_sector_pattern_array.append(sector_pattern_array)

    return snippet_sector_pattern_array


    # for i in range(len(snippet_sector_array)):
    #     if i > 0:
    #         for j in range(len(snippet_sector_array[i])):
    #             if snippet_sector_array[i][j] > 150 * snippet_sector_array[i-1][j]:
    #                 print("Snippet: " + str(i))
    #                 print ("Sector: " + str(j))
    #                 print("_____________________")
    #         print(":::::::::::::::::::::::")



def gauss_fit(fit_data, snippet, out_path):
    '''get data into numpy arrays'''
    for sector in range(len(fit_data)):
        print "Fit sector: " + str(sector)
        pmts_per_sector = fit_data[sector]
        sector_data = get_sector_histo(pmts_per_sector)

        single_fit = fit_single_gauss_in_sector(sector_data)
        double_fit = fit_double_gauss_in_sector(sector_data)

        picture_drawer(sector_data, sector, single_fit, double_fit, out_path)



# def fit_single_gauss_in_sector(sector_pmts):
#     bin_centers = 0.5 * (sector_pmts.bins[1:] + sector_pmts.bins[:-1])
#     plt.xlim([-4, 4])
#
#     single_gauss_params = [c, mu, sigma] = [1, -3, 1]
#     plsq = leastsq(res_single_gauss, single_gauss_params, args=(sector_pmts.entries, bin_centers))
#
#     # print("Fit parameters: " + str(plsq[0]))
#     # print("Fit quality parameter: " + str(plsq[1]))
#
#     plt.clf()
#
#     return plsq
#
#
# def fit_double_gauss_in_sector(sector_pmts):
#
#     bin_centers = 0.5 * (sector_pmts.bins[1:] + sector_pmts.bins[:-1])
#     plt.xlim([-4, 4])
#
#     double_gauss_params = [c1, mu1, sigma1, c2, dmu, sigma2] = [1, 0, 1, 1, 2, 1] # Initial guesses for leastsq
#
#     plsq = leastsq(res_double_gauss_dist, double_gauss_params, args=(sector_pmts.entries, bin_centers))
#
#     # print("Fit parameters: " + str(plsq[0]))
#     # print("Fit quality parameter: " + str(plsq[1])
#
#     plt.clf()
#
#     return plsq



# def picture_drawer(sector_pmts, sector, single_fit, double_fit, out_path):
#     # draw 1 actual data
#     plt.hist(sector_pmts.bins[:-1], len(sector_pmts.bins)-1, weights=sector_pmts.entries)
#
#     bin_centers = 0.5 * (sector_pmts.bins[1:] + sector_pmts.bins[:-1])
#     # plt.xlim([-4, 4])
#
#     single_gauss_fit_graph = single_gaussian(bin_centers, single_fit[0])
#     double_gauss_fit_graph = double_gaussian_dist(bin_centers, double_fit[0])
#
#     # plt.plot(bin_centers, single_gauss_fit_graph, c='k')
#     plt.plot(sector_pmts.bins[:-1], double_gauss_fit_graph, c='r')
#
#     plt.ylabel("theta (deg)")
#     plt.xlabel("phi (deg)")
#
#     plt.savefig(out_path + str(sector) + ".png")
#     plt.close()


def res_double_gauss_dist(p, y, x):
    (c1, mu1, sigma1, c2, dmu, sigma2) = p
    y_fit = double_gaussian_dist(x, p)
    err = y - y_fit
    return err


def res_single_gauss(p, x, y):
    (c, mu, sigma) = p
    y_fit = single_gaussian(x, p)
    err = y - y_fit
    return err


def single_gaussian(x, params):
    (c, mu, sigma) = params
    res = (c/(math.sqrt(2.0 * math.pi * sigma ** 2.0))) * np.exp(-(x - mu) ** 2.0 / (2.0 * sigma ** 2.0))
    return res


def double_gaussian_dist(x, params):
    (c1, mu1, sigma1, c2, dmu, sigma2) = params
    res = (c1/(math.sqrt(2.0 * math.pi * sigma1 ** 2.0))) * np.exp(-(x - mu1) ** 2.0 / (2.0 * sigma1 ** 2.0))\
          + (c2/(math.sqrt(2.0 * math.pi * sigma2 ** 2.0))) * np.exp(-(x - (mu1+dmu)) ** 2.0 / (2.0 * sigma2 ** 2.0))
    return res






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




import statusAlert, src.recoPreparation, color_schemes
from os import mkdir
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.colors as mco
from scipy.optimize import curve_fit
from scipy import optimize
from scipy.stats import gamma


def reco_by_fittig_gauss(pmt_position_class, snippet_class, muon_points, out_path):

    statusAlert.processStatus("search gauss patterns in snippet")

    '''photon data pre-processed and ordered into time snippets'''
    snippet_array = np.asarray(snippet_class.time_snippets)

    '''iterate snippets'''
    fit_results = []
    # for snippet in range(len(snippet_array)):
    for snippet in range(10):
        statusAlert.processStatus("processing snippet: " + str(snippet))

        # gauss_fit_data_vertical = combine_pmts_vertical(pmt_position_class, snippet_class, snippet)
        # mkdir(out_path + "test/vertical2/" + str(snippet))
        # fit_results.append(vertical_gauss_fit(gauss_fit_data_vertical, snippet, out_path + "test/vertical/" + str(snippet) + "/"))

        gauss_fit_data_horizontal = combine_pmts_horizontal(pmt_position_class, snippet_class, snippet)
        # mkdir(out_path + "test/horizontal2/" + str(snippet))
        fit_results.append(horizontal_gauss_fit(gauss_fit_data_horizontal, snippet, out_path + "test/horizontal/" + str(snippet) + "/"))



def vertical_gauss_fit(fit_data, snippet, out_path):
    '''get data into numpy arrays'''

    vertical_fit_results_per_snippet = []

    for vertical_sector in range(len(fit_data)):
        print "Fit vertical sector: " + str(vertical_sector)
        pmts_per_sector = fit_data[vertical_sector]
        # vertical_fit_results_per_snippet.append(fit_gauss_in_sector_vertical(pmts_per_sector, snippet, vertical_sector, out_path))
        find_all_gauss_in_data(pmts_per_sector, snippet, vertical_sector)


    return vertical_fit_results_per_snippet


def horizontal_gauss_fit(fit_data, snippet, out_path):
    '''get data into numpy arrays'''
    horizontal_fit_results_per_snippet = []
    for horizontal_sector in range(len(fit_data)):
        print "Fit horizontal sector: " + str(horizontal_sector)
        pmts_per_sector = fit_data[horizontal_sector]

        horizontal_fit_results_per_snippet.append(fit_gauss_in_sector_horizontal(pmts_per_sector, snippet, horizontal_sector, out_path))
        find_all_gauss_in_data(pmts_per_sector, snippet, horizontal_sector)

    return horizontal_fit_results_per_snippet


class GaussFitResults:
    def __init__(self):
        self.height = 0
        self.location = 0
        self.width = 0


def find_all_gauss_in_data(grouped_sector_pmt_data, snippet, sector):
    histogramm_data = []
    for pmt in grouped_sector_pmt_data:
        for hits in range(pmt.hits):
            histogramm_data.append(pmt.theta)
    n, bins, patches = plt.hist(histogramm_data, bins=60, range=(0., math.pi))

    all_gauss_found = False
    peak = 0
    while not all_gauss_found:
        gauss_parameters = fitgaussian(n)
        print(gauss_parameters)

        if gauss_parameters[0] > 10 and gauss_parameters[1] > 0:
            single_gauss_fit_parameters = GaussFitResults()
            single_gauss_fit_parameters.height = gauss_parameters[0]
            single_gauss_fit_parameters.location = gauss_parameters[1]/n.size*math.pi
            single_gauss_fit_parameters.width = gauss_parameters[2]/n.size*math.pi

            n = n[int(gauss_parameters[1]+5):]

            peak += 1
            plt.savefig("/home/gpu/Analysis/muReconstruction/Output/test/horizontal2/" + str(snippet) + "/" + str(sector) + "+" + str(peak) + ".png")
            plt.close()
        else:
            all_gauss_found = True


def fit_gauss_in_sector_vertical(sector_pmts, snippet, vertical_sector, out_path):
    hist_list = []
    for pmt in sector_pmts:
        for hits in range(pmt.hits):
            hist_list.append(pmt.theta)
    n, bins, patches = plt.hist(hist_list, bins=60, range=(0., math.pi))

    p = fitgaussian(n)
    bin_centers = bins[:-1] + 0.5 * (bins[1:] - bins[:-1])
    plt.plot(bin_centers, gauss(bin_centers, p[0], p[1] / n.size * math.pi, p[2] / n.size * math.pi), 'r--')
    plt.figtext(0.65, 0.7, ("Height: " + str(p[0])
                            + "\nCenter: " + str(p[1] / n.size * math.pi)
                            + "\nWidth: " + str(p[2] / n.size * math.pi)
                            + "\n\nHeight/Width: " + str(p[0] / (p[2] / n.size * math.pi))))

    if 12000 > p[0]/(p[2] / n.size * math.pi) > 3000 and 0.03 > (p[2] / n.size * math.pi) > 0 and n.max() > 25:
        print("Muon entry or exit found in vertical sector: " + str(vertical_sector) + " at theta: " + str(p[1] / n.size * math.pi))

    plt.savefig(out_path + str(vertical_sector) + ".png")
    plt.close()


def fit_gauss_in_sector_horizontal(sector_pmts, snippet, horizontal_sector, out_path):
    hist_list = []
    for pmt in sector_pmts:
        for hits in range(pmt.hits):
            hist_list.append(pmt.phi)
    n, bins, patches = plt.hist(hist_list, bins=60, range=(0., math.pi))


    p = fitgaussian(n)
    bin_centers = bins[:-1] + 0.5 * (bins[1:] - bins[:-1])
    plt.xlim([-4, 4])
    plt.plot(bin_centers, gauss(bin_centers, p[0], p[1] / n.size * math.pi, p[2] / n.size * math.pi), 'r--')
    plt.figtext(0.65, 0.7, ("Height:  %f \nCenter:  %f \nWidth:  %f\n\nHeight/Width: %f"
                            % (p[0], p[1] / n.size * math.pi, p[2] / n.size * math.pi,
                               p[0] / (p[2] / n.size * math.pi))))

    if 12000 > p[0]/(p[2] / n.size * math.pi) > 3000 and 0.03 > (p[2] / n.size * math.pi) > 0 and n.max() > 25:
        # print("n.max = " + str(n.max()))
        # print p[0] / (p[2] / n.size * math.pi)
        print("Muon entry or exit found in horizontal sector: " + str(horizontal_sector) + " at phi: " + str(p[1] / n.size * math.pi))

    plt.savefig(out_path + str(horizontal_sector) + ".png")
    plt.close()


# Function to be fitted
def gauss(x, height, center, width):
    return height*np.exp(-(((center-x)/width)**2)/2)


def gaussian(height, center, width):
    width = float(width)
    return lambda x: height*np.exp(-(((center-x)/width)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    X = np.arange(data.size)
    x = np.sum(X * data) / np.sum(data)
    width = np.sqrt(np.abs(np.sum((X - x) ** 2 * data) / np.sum(data)))
    height = data.max()

    # print(x/data.size * math.pi)
    # print(width/data.size * math.pi)
    # print(height)

    return height, x, width



def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) - data)
    p, success = optimize.leastsq(errorfunction, params)
    return p


class pmt_data_for_fit_class:
    def __init__(self):  # this method creates the class object.
        self.x_sector = 0
        self.y_sector = 0
        self.hits = 0
        self.phi = 0
        self.theta = 0
        self.hist_array = []


def combine_pmts_vertical(pmt_position_class, snippet_class, snippet):
    return_array = [[] for _ in range(pmt_position_class.x_sectors)]
    for pmt_id in range(len(pmt_position_class.id)):
        pmt_data = pmt_data_for_fit_class()
        pmt_data.x_sector = pmt_position_class.is_x_sector[pmt_id]
        pmt_data.y_sector = pmt_position_class.is_y_sector[pmt_id]
        pmt_data.hits = snippet_class.time_snippets[snippet][pmt_id]
        pmt_data.phi = pmt_position_class.phi_position[pmt_id]
        pmt_data.theta = pmt_position_class.theta_position[pmt_id]

        return_array[pmt_position_class.is_x_sector[pmt_id]].append(pmt_data)
    return return_array


def combine_pmts_horizontal(pmt_position_class, snippet_class, snippet):
    return_array = [[] for _ in range(pmt_position_class.y_sectors)]
    for pmt_id in range(len(pmt_position_class.id)):
        pmt_data = pmt_data_for_fit_class()
        pmt_data.x_sector = pmt_position_class.is_x_sector[pmt_id]
        pmt_data.y_sector = pmt_position_class.is_y_sector[pmt_id]
        pmt_data.hits = snippet_class.time_snippets[snippet][pmt_id]
        pmt_data.phi = pmt_position_class.phi_position[pmt_id]
        pmt_data.theta = pmt_position_class.theta_position[pmt_id]

        return_array[pmt_position_class.is_y_sector[pmt_id]].append(pmt_data)
    return return_array


def hit_pmt_in_total_event(source):
    statusAlert.processStatus("Start read of entire event: ")
    pmt_time_resolution = 3.0
    data_read = TreeReadFunc.readPhotonRecoData2(source, pmt_time_resolution)

    return_class = photon_snippet_output_c()

    pmt_array = [0] * 17739
    for index, photon_hit in enumerate(data_read):
        if photon_hit.pmt_id < 17739:
            pmt_array[photon_hit.pmt_id] += 1
    return_class.add_snippet_array(pmt_array)
    statusAlert.processStatus("Done")
    return return_class


