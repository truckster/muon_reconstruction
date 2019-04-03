import statusAlert, recoPreparation, color_schemes, contour_analyze, reco_from_contour, PointVecDist
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.ndimage.filters import gaussian_filter
import os
import diff_event_analysis


def entry_exit_detector(pmt_position_class, snippet_class, muon_points, out_path):
    peak_heights = []
    contour_data_array = []
    for snippet_index, pmt_array in enumerate(snippet_class.time_snippets):
    # for snippet in range(compute_snippets):
        statusAlert.processStatus("processing snippet: " + str(snippet_index))
        contour_data = contour_data_reader(pmt_position_class, pmt_array)
        # reco_from_contour.concentric_level_finder(contour_data, snippet)
        # reco_from_contour.level_area_difference(contour_data, snippet)
        # reco_from_contour.container(contour_data, contour_raw)
        # reco_from_contour.gradient(contour_data)
        contour_data_array.append(contour_data)
        if snippet_index > 0:
            snippet_diff = np.asarray(snippet_class.time_snippets[snippet_index]) \
                           - np.asarray(snippet_class.time_snippets[snippet_index - 1])
            contour_data_diff = contour_data_reader(pmt_position_class, snippet_diff)
        peak_heights.append(contour_data[-1][0].height)

    reconstructed_points1 = reco_from_contour.peak_compare(peak_heights, contour_data_array)

    return reconstructed_points1


def reco_result_writer(output_path, result_array):
    '''create output file'''
    result_file = open(output_path + "results.txt", 'a')
    result_file.write("----- Reconstructed Values (Frames) ------" + '\n')

    for point in result_array:
        result_file.write("Point found in snippet: " + str(point.snippet) + '\n')
        result_file.write("Phi: " + str(point.phi) + '\n')
        result_file.write("Theta: " + str(point.theta) + '\n')
        result_file.write("________________________________________" + '\n')

    result_file.write("End of event" + '\n' + '\n')
    result_file.close()


def coordinate_calculation(reco_points):
    return_points = []
    for orientation_index, orientation in enumerate(reco_points):
        orientation_array = []
        for point_index, point in enumerate(orientation):
            contours = point.contour_data
            point.x_coordinate_rad = contours.center[0]
            point.y_coordinate_rad = contours.center[1]
            point.x_coordinate_deg = contours.center[0]/math.pi * 180.0
            point.y_coordinate_deg = contours.center[1]/math.pi * 180.0

            orientation_array.append(point)
        return_points.append(orientation_array)


def orientation_resolver(reco_points):
    return_points = []
    for orientation_index, orientation in enumerate(reco_points):
        for point_index, point in enumerate(orientation):
            if orientation_index is 2:
                if point.x_coordinate_deg > 0:
                   point.x_coordinate_deg -= 180
                   point.x_coordinate_rad -= math.pi
                else:
                    point.x_coordinate_deg += 180
                    point.x_coordinate_rad += math.pi
            elif orientation_index is 3:
                if point.y_coordinate_deg > 0:
                    point.y_coordinate_deg -= 90
                    point.y_coordinate_rad -= math.pi/2
                else:
                    point.y_coordinate_deg += 90
                    point.y_coordinate_rad += math.pi/2
            elif orientation_index is 1:
                if point.x_coordinate_deg > 0:
                    point.x_coordinate_deg -= 180
                    point.x_coordinate_rad -= math.pi
                else:
                    point.x_coordinate_deg += 180
                    point.x_coordinate_rad += math.pi
                if point.y_coordinate_deg > 0:
                    point.y_coordinate_deg -= 90
                    point.y_coordinate_rad -= math.pi/2
                else:
                    point.y_coordinate_deg += 90
                    point.y_coordinate_rad += math.pi/2
            else:
                pass
            return_points.append(point)
    return return_points


def process_events(reco_points, merge_radius, intersec_radius, performance_class, pmt_positions):
    coordinate_calculation(reco_points)
    pure_points = orientation_resolver(reco_points)
    diff_event_analysis.point_merger(pure_points, merge_radius)

    performance_class.unmerged.append(len(pure_points))
    pure_points = diff_event_analysis.pole_merger(pure_points)
    performance_class.merged.append(len(pure_points))

    calc_kartesian_coordinates(pure_points, intersec_radius)
    find_closest_pmt(pure_points, pmt_positions)

    return pure_points


def calc_kartesian_coordinates(found_points, intersec_radius):
    for point in found_points:
        vec = PointVecDist.D3Vector()
        point.real_x = intersec_radius * math.cos(point.y_coordinate_rad)*math.cos(point.x_coordinate_rad)
        point.real_y = intersec_radius * math.cos(point.y_coordinate_rad)*math.sin(point.x_coordinate_rad)
        point.real_z = intersec_radius * math.sin(point.y_coordinate_rad)

        vec.x = point.real_x
        vec.y = point.real_y
        vec.z = point.real_z

        point.D3_vector = vec


def find_closest_pmt(found_points, pmt_positions):
    for point in found_points:
        distance = np.inf
        pmt = 0
        for pmt_id in range(len(pmt_positions.id)):
            distance_tmp = math.sqrt(pow(point.real_x-pmt_positions.x_position[pmt_id], 2)
                                     +pow(point.real_y-pmt_positions.y_position[pmt_id], 2)
                                     +pow(point.real_y-pmt_positions.y_position[pmt_id], 2))
            if distance_tmp < distance:
                distance = distance_tmp
                pmt = pmt_id

        point.closest_pmt = pmt

