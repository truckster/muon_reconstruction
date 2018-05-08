import reconstructionAlg, statusAlert
import numpy as np
import matplotlib.path as mpath
import math


def intersec_crosscheck(diff_contours, reco_points):
    better_reco_points = []
    for point_index, point in enumerate(reco_points):
        minimum = 0
        maximum = 0
        min_frame = 0
        max_frame = 0
        for frame, contour_data in enumerate(diff_contours):
            # statusAlert.processStatus("processing snippet: " + str(frame+1))
            levels = []
            for level_index, level in enumerate(contour_data):
                for patch_index, patch in enumerate(level):
                    if mpath.Path(patch.contour_coordinates).contains_point(point.contour_coordinates[0]):
                        levels.append(patch.height)
                        if patch.height < minimum:
                            minimum = patch.height
                            min_frame = frame
                        if patch.height > maximum:
                            maximum = patch.height
                            max_frame = frame

        if max_frame - min_frame < 3 and minimum < 100 and maximum - minimum > 400:
            better_reco_points.append(point)
    return better_reco_points


def point_merger(reco_points):
    print(len(reco_points))
    for point1 in reco_points:
        for point2 in reco_points:
            frame_diff = abs(point1.frame - point2.frame)
            distance = math.sqrt(pow(point1.coordinates[0]-point2.coordinates[0], 2)
                                         + pow(point1.coordinates[1]-point2.coordinates[1], 2))
            if frame_diff < 1 and distance < 0.9:
                reco_points.remove(point2)
    print(len(reco_points))
