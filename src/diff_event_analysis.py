import reconstructionAlg, statusAlert
import numpy as np
import matplotlib.path as mpath
import math


def intersec_crosscheck(diff_contours, reco_points):
    for frame, contour_data in enumerate(diff_contours):
        statusAlert.processStatus("processing snippet: " + str(frame))
        for point_index, point in enumerate(reco_points):
            for level_index, level in enumerate(contour_data):
                for patch_index, patch in enumerate(level):
                    if mpath.Path(patch.contour_coordinates).contains_point(point.contour_coordinates[0]):
                        print(patch.height)
            # point_hit_array[point_index].append(current_level)
