import matplotlib.pyplot as plt
import matplotlib.path as mpath
import math
import numpy as np


class ReconstructIntersecPoint:
    def __init__(self):  # this method creates the class object.
        self.snippet = 0
        self.phi = 0
        self.theta = 0


def peak_compare(peak_array, contour_data):
    return_array = []
    for index, value in enumerate(peak_array):
        if len(peak_array)-1 > index > 0:
            if value > peak_array[index-1]+30 and value > peak_array[index+1]+30:
                point_class = ReconstructIntersecPoint()
                point_class.snippet = index
                point_class.phi = contour_data[index][-1][0].center[0]
                point_class.theta = contour_data[index][-1][0].center[1]

                return_array.append(point_class)

    return return_array


def concentric_level_finder(contour_data, snippet):
    # print(contour_raw.collections[0].get_paths())
    for level_index, level in enumerate(contour_data):
        for patch_index, patch in enumerate(level):
            center_max_phi = patch.center[0]
            center_max_theta = patch.centers[1]
            for iso_hit in range(contour_data[level].iso_hit_patches):
                phi_current = contour_data[level].centers[iso_hit][0]
                theta_current = contour_data[level].centers[iso_hit][1]
                print(math.sqrt(pow(center_max_phi-phi_current, 2) + pow(center_max_theta-theta_current, 2)))


def level_area_difference(contour_data, snippet):
    for level in range(len(contour_data)):
        print("------------------------------------------")
        try:
            top_level_area = contour_data[-1].extents[0]
            for iso_hit in range(contour_data[level].iso_hit_patches):
                level_area = contour_data[level].extents[iso_hit]
                print(level_area - top_level_area)
        except IndexError:
            print("No data available")


def container(contour_data):
    for level_count, level_data in enumerate(contour_data[:-1]):
        contour_lines_raw_level = contour_raw.collections[level_count].get_paths()
        next_contour_level_center_points = contour_data[level_count+1].centers
        print("--------LEVEL!!!----------")
        for is_level_contour_count, iso_level_contour_array in enumerate(contour_lines_raw_level):
            nppath = np.asarray(iso_level_contour_array.vertices)
            path = mpath.Path(nppath)
            for next_level_center_count, next_level_center_point in enumerate(next_contour_level_center_points):
                isit = path.contains_points(next_contour_level_center_points)


