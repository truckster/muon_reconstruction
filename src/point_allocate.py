"""Allocate the entry and exit points after finding points by looking at entire event and reconstructing the time
 frames. Here we use difference timeframes"""
import statusAlert, reconstructionAlg
import numpy as np
import matplotlib.path as mpath


def allocate_points(contour_diffs, intersec_points, intersec_times):
    for frame, contour in enumerate(contour_diffs):
        print(contour[0][0].level)