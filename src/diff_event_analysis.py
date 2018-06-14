import reconstructionAlg, statusAlert
import numpy as np
import matplotlib.path as mpath
import math


def intersec_crosscheck(diff_contours_array, reco_points):
    """Checks if the previously detected points are relevant

    Therefore iterates through the framed data and checks if at found points the muon really enters or exits"""
    better_reco_points = []
    for orientation_index, orientation in enumerate(reco_points):
        points_of_orientation = []
        for point_index, point in enumerate(orientation):
            contours = point.contour_data
            minimum = 0
            maximum = 0
            min_frame = 0
            max_frame = 0
            diff_contours = diff_contours_array[orientation_index]
            for frame, contour_data in enumerate(diff_contours):
                # statusAlert.processStatus("processing snippet: " + str(frame+1))
                for level_index, level in enumerate(contour_data):
                    for patch_index, patch in enumerate(level):
                        if mpath.Path(patch.contour_coordinates).contains_point(contours.contour_coordinates[0]):
                            if patch.height < minimum:
                                minimum = patch.height
                                min_frame = frame
                            if patch.height > maximum:
                                maximum = patch.height
                                max_frame = frame

            if max_frame - min_frame < 3 and minimum < 100 and maximum - minimum > 150:
                point.frame = max_frame
                points_of_orientation.append(point)
        better_reco_points.append(points_of_orientation)
    return better_reco_points


def point_merger(reco_points):
    for point1 in reco_points:
        for point2 in reco_points:
            frame_diff = abs(point1.frame - point2.frame)
            distance = math.sqrt(pow(point1.coordinates[0]-point2.coordinates[0], 2)
                                 + pow(point1.coordinates[1]-point2.coordinates[1], 2))
            if point1 is not point2 and frame_diff < 1 and distance < 0.9:
                reco_points.remove(point2)
