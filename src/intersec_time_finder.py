"""Find the entry and exit times after finding points by looking at entire event. Here we use time frames"""
import statusAlert, reconstructionAlg
import numpy as np
import matplotlib.path as mpath


def find_times(contour_array, intersec_points):
    statusAlert.processStatus("Reconstruct intersection frames of muon")
    point_hit_array = [[] for _ in range(len(intersec_points))]
    for frame, contour_data in enumerate(contour_array):
        # statusAlert.processStatus("processing snippet: " + str(frame))
        for point_index, point in enumerate(intersec_points):
            current_level=0
            # print("X: " + str(point.x_coordinate_deg))
            # print("Y: " + str(point.y_coordinate_deg))
            point_coordinate = [point.x_coordinate_rad, point.y_coordinate_rad]
            for level_index, level in enumerate(contour_data):
                for patch_index, patch in enumerate(level):
                    if mpath.Path(patch.contour_coordinates).contains_point(point_coordinate):
                        current_level = patch.height
            point_hit_array[point_index].append(current_level)

    for frame, point_hit_list in enumerate(point_hit_array):
        intersec_points[frame].frame = point_hit_list.index(max(point_hit_list))


def reco_result_writer(output_path, return_array):
    '''Add reconstructed intersection frames to file'''
    result_file = open(output_path + "results.txt", 'a')
    result_file.write("----- Reconstructed Values (Intersection frames)------" + '\n')
    for point_index, point in enumerate(return_array):
        result_file.write("Point: " + str(point_index) + '\n')
        result_file.write("Frame: " + str(point.frame) + '\n')

    result_file.write("End of event" + '\n' + '\n')
    result_file.close()
