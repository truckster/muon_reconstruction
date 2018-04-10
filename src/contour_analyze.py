import matplotlib.path as mpath
import numpy as np


class ContourData:
    """class to store the contour data within a given snippet for one specific contour level"""
    def __init__(self):  # this method creates the class object.
        self.height = 0
        self.iso_hit_patches = 0
        self.centers = []
        self.extents = []
        self.level = 0
        self.contour_coordinates = []
        self.contour_path = 0


def get_contour_data(contour_output):
    """creates an array in which the contour data per level is stored for a given snippet"""
    return_data = []
    for level in range(len(contour_output.levels)):
        patch_data = ContourData()
        patch_data.height = contour_output.levels[level]
        patch_data.level = level
        patch_data.iso_hit_patches = len(contour_output.collections[level].get_paths())
        calc_contour_center(contour_output, level, patch_data)
        calc_contour_extent(contour_output, level, patch_data)

        try:
            patch_data.contour_path = contour_output.collections[level].get_paths()[0]
            patch_data.contour_coordinates = contour_output.collections[level].get_paths()[0].vertices
        except IndexError:
            pass

        return_data.append(patch_data)

    return return_data


def calc_contour_center(contour_output, level, return_class):
    for contour in contour_output.collections[level].get_paths():
        x_vals_sum = 0
        y_vals_sum = 0
        for coordinates in contour.vertices[:-1]:
            x_vals_sum += coordinates[0]
            y_vals_sum += coordinates[1]

        if len(contour.vertices) > 1:
            x_coord = x_vals_sum/float((len(contour.vertices)-1))
            y_coord = y_vals_sum/float((len(contour.vertices)-1))
        else:
            x_coord = x_vals_sum
            y_coord = y_vals_sum

        return_class.centers.append([x_coord, y_coord])


def calc_contour_extent(contour_output, level, return_class):
    for contour in contour_output.collections[level].get_paths():
        area = 0
        for pair in range(len(contour.vertices[:-1])):
            x0 = contour.vertices[pair][0]
            y0 = contour.vertices[pair][1]

            dx = contour.vertices[pair+1][0] - contour.vertices[pair][0]
            dy = contour.vertices[pair+1][1] - contour.vertices[pair][1]

            area += abs(0.5*(y0*dx - x0*dy))

        return_class.extents.append(area)


def is_center_in_next_higher_level_contained(contour_center_current_level, contour_path_next_level):
    """Function to check if the center of a given contour is enclosed by a contour of the next level"""
    # make suiting path object from numpy array
    path = mpath.Path(np.asarray(contour_path_next_level))
    # check if given center
    is_contained = path.contains_points(contour_center_current_level)

    return is_contained




