"""Find the entry and exit times after finding points by looking at entire event. Here we use time frames"""
import statusAlert, reconstructionAlg
import matplotlib.path as mpath
import PointVecDist


def find_times(contour_array, intersec_points):
    # statusAlert.processStatus("Reconstruct intersection frames of muon")
    point_hit_array = [[] for _ in range(len(intersec_points))]
    for frame, contour_data in enumerate(contour_array):
        # statusAlert.processStatus("processing snippet: " + str(frame))
        for point_index, point in enumerate(intersec_points):
            current_level = 0
            point_coordinate = [point.x_coordinate_rad, point.y_coordinate_rad]
            for level_index, level in enumerate(contour_data):
                for patch_index, patch in enumerate(level):
                    # if mpath.Path(patch.contour_coordinates).contains_point(point_coordinate):
                    #     current_level = patch.height
                    if PointVecDist.point_distance2d(patch.center, point_coordinate) < 0.5:
                        current_level = patch.height
            point_hit_array[point_index].append(current_level)
    for frame, point_hit_list in enumerate(point_hit_array):
        # print point_hit_list
        intersec_points[frame].frame = point_hit_list.index(max(point_hit_list))


def find_times_fht(fht_list, intersec_points):
    for point in intersec_points:
        for frame, frame_list in enumerate(fht_list):
            if frame_list[point.closest_pmt]:
                point.frame = frame


def reco_result_writer(output_path, return_array):
    '''Add reconstructed intersection frames to file'''
    result_file = open(output_path + "results.txt", 'a')
    result_file.write("----- Reconstructed Values (Intersection frames)------" + '\n')
    for point_index, point in enumerate(return_array):
        result_file.write("Point: " + str(point_index) + '\n')
        result_file.write("Frame: " + str(point.frame) + '\n')

    result_file.write("End of event" + '\n' + '\n')
    result_file.close()


def time_finder_performance(mc_truth_list, reco_track_list, time_window, performance_class):
    truth_list = []
    reco_list = []

    for muon in mc_truth_list:
        for point in muon:
            truth_list.append(round(point.intersec_time))
    truth_list.sort()

    for track in reco_track_list:
        reco_list.append(track.entry_point.frame * time_window)
        reco_list.append(track.exit_point.frame * time_window)
        # reco_list.append(track.entry_point.frame)
        # reco_list.append(track.exit_point.frame)

    reco_list.sort()

    if len(truth_list) == len(reco_list):
        for i in range(len(truth_list)):
            performance_class.time_difference_list.append(reco_list[i] - truth_list[i])





