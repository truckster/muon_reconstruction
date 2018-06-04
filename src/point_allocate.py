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
        for patch_low in low_level:
            lengths = []
            angles = []
            if patch_low.height < 0:
                for patch_high in top_level:
                    lengths.append(math.sqrt(pow(patch_high.center[0] - patch_low.center[0], 2)
                                             + pow(patch_high.center[1] - patch_low.center[1], 2)))
                    angles.append(np.arctan((patch_high.center[1] - patch_low.center[1])
                                            / patch_high.center[0] - patch_low.center[0]))
                    same_area = have_same_surrounder(contour, patch_low, patch_high)
                if len(lengths) > 0 and min(lengths) < 0.25 and same_area:
                    print("Frame: " + str(frame + 1))
                    print("Length: " + str(min(lengths)))
                    mindex = np.argmax(lengths)
                    print("Angle: " + str(angles[mindex]/math.pi * 180.0))


def have_same_surrounder(contour_data, patch1, patch2):
    they_do = False
    for level in contour_data:
        for patch in level:
            if mpath.Path(patch.contour_coordinates).contains_point(patch2.contour_coordinates[0]) and \
                   mpath.Path(patch.contour_coordinates).contains_point(patch1.contour_coordinates[0]):
                    they_do = True

    return they_do
