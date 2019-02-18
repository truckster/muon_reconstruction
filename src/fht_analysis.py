import statusAlert, color_schemes
import reconstructionAlg
import matplotlib.pyplot as plt
import image_creator

def fht_reader(fht_array):
    print(sorted(fht_array[0]))
    print(len(fht_array[0]))


def fht_drawer(pmt_position_class, fht_array_frames, muon_points, out_path, reco_points=None):
    statusAlert.processStatus("Draw FHT picture")
    '''Draw the detector picture for the certain time snippet'''
    for frame_index, frame in enumerate(fht_array_frames):
        draw_fht_picture(pmt_position_class, frame, frame_index, muon_points, out_path, reco_points, "absolute")


def draw_fht_picture(pmt_position_class, fht_array, frame_index, muon_points, out_path, reco_points=None, mode=None):
    statusAlert.processStatus("Creating graphics")
    fig = plt.figure(num=None, figsize=(20, 10))

    '''Analysis design'''
    ax1 = fig.add_subplot(111, axisbg='gainsboro', projection='aitoff')
    # scatterHits = ax1.scatter(pmt_position_class.phi_position, pmt_position_class.theta_position, marker='o', c='k')
    scatterHits = ax1.scatter(pmt_position_class.phi_position, pmt_position_class.theta_position2,
                              c=fht_array, cmap=color_schemes.analysis_design(), edgecolor='k', label='hit', )
    cb = fig.colorbar(scatterHits)

    # draw_sector_lines(pmt_position_class)
    image_creator.draw_muon_points(muon_points, ax1)
    if reco_points is None:
        a=1
    else:
        image_creator.draw_reconstructed_points(reco_points, ax1)
    # draw_pattern_results(sector_pattern_array, template_radius, ax1)
    # draw_reco_muons(muon_finder_psnippet)

    plt.ylabel("theta (deg)")
    plt.xlabel("phi (deg)")

    plt.savefig(out_path + str(frame_index) + "fht.png", bbox_inches='tight')

    plt.close()
