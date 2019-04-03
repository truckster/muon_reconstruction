import heapq
import matplotlib.pyplot as plt
import PointVecDist


def do(events, photons, pmt_positions):
    n_lab = 1.484
    c = 299.792 # mm/ns
    for track in events:
        print ":::::::::::::::::::::::::::::::::"
        step_size = 50
        entry_point = track.entry_point
        interesting_time = entry_point.frame*5
        interesting_light = [[] for i in range(17739)]
        # print interesting_light
        photon_num = 0
        for photon in photons:
            if interesting_time - 5 < photon.hit_time < interesting_time + 5 and photon.pmt_id < 17739:
                interesting_light[photon.pmt_id].append(photon.hit_time)
                photon_num += 1
        # print len(interesting_light)
        # print len(interesting_light[0])
        # print len(interesting_light[1])
        # print photon_num
        for x_step in range(-2, 3):
            for y_step in range(-2, 3):
                for z_step in range(-3, 0):
                    time_diff_hist = []
                    test_pos = PointVecDist.D3Vector()
                    test_pos.x = entry_point.real_x + x_step*step_size
                    test_pos.y = entry_point.real_y + y_step*step_size
                    test_pos.z = entry_point.real_z + z_step*step_size
                    # for pmt_id, pmt in enumerate(interesting_light):
                        # print pmt_id
                        # print len(pmt)
                    #     if len(pmt) > 100:
                    #         pmt_pos = PointVecDist.D3Vector()
                    #         pmt_pos.x = pmt_positions.x_position[pmt_id]
                    #         pmt_pos.y = pmt_positions.y_position[pmt_id]
                    #         pmt_pos.z = pmt_positions.z_position[pmt_id]
                    #         dist_spot_pmt = PointVecDist.d3_distance(test_pos, pmt_pos)
                    #         tof_theo = dist_spot_pmt / (c / n_lab)
                    #         for photon_hit in pmt:
                    #             time_diff_hist.append(photon_hit-tof_theo)
                    # hist = plt.hist(time_diff_hist)
                    # plt.show()
        interesting_light = 0

