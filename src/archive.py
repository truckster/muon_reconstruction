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
