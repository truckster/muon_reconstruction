import matplotlib.pyplot as plt
import statusAlert, recoPreparation


def drawPMTHitMapHisto(hitArray, muonPoints, pmt_position_class, outPath, limit):
    for snippet in range(limit):
        statusAlert.processStatus("Processing histogram number " +str(snippet) + " and print detector picture")

        bary_center_coordinates = reconstructionAlg.find_photon_barycenter(pmt_position_class.phi_position, pmt_position_class.theta_position, hitArray[snippet])

        phiMuon = []
        thetaMuon = []
        for i in range(len(muonPoints)):
                for j in range(len(muonPoints[i])):
                    phiMuon.append(recoPreparation.calcPMTPolarPhi(muonPoints[i][j]))
                    thetaMuon.append(recoPreparation.calcPMTPolarTheta(muonPoints[i][j]))

        fig = plt.figure(num=None, figsize=(20, 10))
        ax1 = fig.add_subplot(111)
        scatterHits = ax1.scatter(pmt_position_class.phi_position, pmt_position_class.theta_position, c=hitArray[snippet], cmap='summer', label='hit')
        scatterMuon = ax1.scatter(phiMuon, thetaMuon, facecolors='none', edgecolors='w',  marker='o', s=120)
        scatter_photon_barycenter = ax1.scatter(bary_center_coordinates.phi_center_c,
                                                bary_center_coordinates.theta_center_c, facecolors='none', edgecolors='w',
                                                marker='s', s=120)
        plt.ylabel("theta (deg)")
        plt.xlabel("phi (deg)")
        cb=fig.colorbar(scatterHits)
        cb.set_label("Number of Photons")

        plt.savefig(outPath + str(snippet) + ".png")
        statusAlert.processStatus("Done")
        plt.close()
