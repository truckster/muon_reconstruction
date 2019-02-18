import heapq
import matplotlib.pyplot as plt
import PointVecDist


def do(events, photon_frames, pmt_positions):
    for track in events:
        print ":::::::::::::::::::::::::::::::::"
        step_size = 300
        entry_point = track.entry_point
        interesting_light = photon_frames.time_snippets[entry_point.frame+2]
        bright_pmts = heapq.nlargest(10, enumerate(interesting_light), key=lambda x: x[1])
        # print interesting_light.index(heapq.nlargest(10, interesting_light))
        dist_distribution_hist = []
        malax = 0
        for x_step in range(-2, 3):
            for y_step in range(-2, 3):
                for z_step in range(-3, 0):
                    dist_hist = []
                    pos_hist = []
                    test_pos = PointVecDist.D3Vector()
                    test_pos.x = entry_point.real_x + x_step*step_size
                    test_pos.y = entry_point.real_y + y_step*step_size
                    test_pos.z = entry_point.real_z + z_step*step_size
                    for pmt in bright_pmts:
                        pmt_pos = PointVecDist.D3Vector()
                        pmt_pos.x = pmt_positions.x_position[pmt[0]]
                        pmt_pos.y = pmt_positions.y_position[pmt[0]]
                        pmt_pos.z = pmt_positions.z_position[pmt[0]]

                        dist_spot_pmt = PointVecDist.d3_distance(test_pos, pmt_pos)
                        dist_hist.append(dist_spot_pmt)
                        pos_hist.append(pmt_pos)
                    hist = plt.hist(dist_hist, bins=10)
                    if max(hist[0]) >= malax:
                        print hist[0]
                        print x_step
                        print y_step
                        print z_step
                        malax = hist[0][0]
                        # plt.show()
                    # print x_step
                    # print y_step
                    # print z_step

