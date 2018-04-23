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
            for level_index, level in enumerate(contour_data):
                for patch_index, patch in enumerate(level):
                    if mpath.Path(patch.contour_coordinates).contains_point(point):
                        current_level = patch.height
                        # print(current_level)
            point_hit_array[point_index].append(current_level)

    frame_array = []
    for frame, point_hit_list in enumerate(point_hit_array):
        frame_array.append(point_hit_list.index(max(point_hit_list)))

    return frame_array


def reco_result_writer(output_path, result_array):
    '''Add reconstructed intersection frames to file'''
    result_file = open(output_path + "results.txt", 'a')
    result_file.write("----- Reconstructed Values (Intersection frames)------" + '\n')

    for frame in result_array:
        result_file.write("Frame: " + str(frame) + '\n')

    result_file.write("End of event" + '\n' + '\n')
    result_file.close()
