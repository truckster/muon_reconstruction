"""Allocate the entry and exit points after finding points by looking at entire event and reconstructing the time
 frames. Here we use difference timeframes"""
import statusAlert, reconstructionAlg, PointVecDist
import numpy as np
import matplotlib.path as mpath
import math


class DirectionReco:
    def __init__(self):
        self.frame = 0
        self.reference_point = []
        self.direction = []


def allocate_points(contour_diffs, intersec_points, intersec_times):
    for frame, contour in enumerate(contour_diffs):
        top_level = contour[-1]
        low_level = contour[0]
        print(frame+1)
        for patch_low in low_level:
            lengths = []
            if patch_low.height < 0:
                for patch_high in top_level:
                    lengths.append(math.sqrt(pow(patch_high.center[0] - patch_low.center[0], 2)
                                             + pow(patch_high.center[1] - patch_low.center[1], 2)))
                    same_area = have_same_surrounder(contour, patch_low, patch_high)
                if len(lengths) > 0 and same_area:
                    print(min(lengths))


def have_same_surrounder(contour_data, patch1, patch2):
    they_do = False
    for level in contour_data:
        for patch in level:
            if mpath.Path(patch.contour_coordinates).contains_point(patch2.contour_coordinates[0]) and \
                   mpath.Path(patch.contour_coordinates).contains_point(patch1.contour_coordinates[0]):
                    they_do = True

    return they_do
