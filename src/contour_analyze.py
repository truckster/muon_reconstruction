import matplotlib.path as mpath
import numpy as np
import reconstructionAlg, statusAlert


class ContourDataSinglePatch:
    """class to store the contour data within a given snippet for one specific contour level"""
    def __init__(self):  # this method creates the class object.
        self.height = 0
        self.iso_hit_patches = 0
        self.center = []
        self.extent = []
        self.level = 0
        self.contour_coordinates = []
        self.contour_path = 0


def get_contour_data(contour_output):
    """creates an array in which the contour data per level is stored for a given snippet
        returns a 2D-Array per snippet [level][patch]=ContourDataSinglePatch_Class"""
    return_data_total_snippet = []
    for level in range(len(contour_output.levels)):
        """iterates over all contour levels in the snippet"""
        return_data_iso_hit_patches = []
        for patch in contour_output.collections[level].get_paths():
            """iterates over all iso hit patches of the given level in the contour"""
            patch_data = ContourDataSinglePatch()
            patch_data.height = contour_output.levels[level]
            patch_data.level = level
            patch_data.iso_hit_patches = len(contour_output.collections[level].get_paths())
            patch_data.center = calc_contour_center(patch)
            patch_data.extent = calc_contour_extent(contour_output, level)

            try:
                patch_data.contour_path = patch
                patch_data.contour_coordinates = patch.vertices
            except IndexError:
                pass

            return_data_iso_hit_patches.append(patch_data)

        return_data_total_snippet.append(return_data_iso_hit_patches)

    return return_data_total_snippet


def collect_contour_data(photons_in_time_window, PmtPositions):
    contour_data_array = []
    diff_contour_data_array = []
    for frame, pmt_array in enumerate(photons_in_time_window.time_snippets):
        statusAlert.processStatus("processing snippet: " + str(frame))
        contour_data_array.append(reconstructionAlg.contour_data_reader(PmtPositions, pmt_array))
        if frame > 0:
            snippet_diff = np.asarray(pmt_array) - np.asarray(photons_in_time_window.time_snippets[frame - 1])
            diff_contour_data_array.append(reconstructionAlg.contour_data_reader(PmtPositions, snippet_diff))

    return contour_data_array, diff_contour_data_array


def calc_contour_center(patch_path_data):
    x_vals_sum = 0
    y_vals_sum = 0
    for coordinates in patch_path_data.vertices[:-1]:
        x_vals_sum += coordinates[0]
        y_vals_sum += coordinates[1]

    if len(patch_path_data.vertices) > 1:
        x_coord = x_vals_sum/float((len(patch_path_data)-1))
        y_coord = y_vals_sum/float((len(patch_path_data)-1))
    else:
        x_coord = x_vals_sum
        y_coord = y_vals_sum

    return [x_coord, y_coord]


def calc_contour_extent(contour_output, level):
    for contour in contour_output.collections[level].get_paths():
        area = 0
        for pair in range(len(contour.vertices[:-1])):
            x0 = contour.vertices[pair][0]
            y0 = contour.vertices[pair][1]

            dx = contour.vertices[pair+1][0] - contour.vertices[pair][0]
            dy = contour.vertices[pair+1][1] - contour.vertices[pair][1]

            area += abs(0.5*(y0*dx - x0*dy))

        return area


def is_center_in_next_higher_level_contained(contour_center_current_level, contour_path_next_level):
    """Function to check if the center of a given contour is enclosed by a contour of the next level"""
    # make suiting path object from numpy array
    path = mpath.Path(np.asarray(contour_path_next_level))
    # check if given center
    is_contained = path.contains_points(contour_center_current_level)

    return is_contained




