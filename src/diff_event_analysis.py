import reconstructionAlg, statusAlert, total_event_reconstruction
import numpy as np
import matplotlib.path as mpath
import math


def entry_exit_detector(contour_data):
    real_top_patches = []
    statusAlert.processStatus("Searching entry and exit points from frames: ")
    for orientation_index, data in enumerate(contour_data):
        for frame, frame_data in enumerate(data):
            statusAlert.processStatus("Frame: " + str(frame))
            found_points = standalone_contour_lines(frame_data, orientation_index)
            real_top_patches.append(toplevel_check(found_points, frame_data, orientation_index))
            # top_levels = diff_intersec_finder(data, orientation_index)
            # real_top_patches.append(toplevel_check(top_levels, data[0], orientation_index))
    return real_top_patches


def standalone_contour_lines(contour_data, orientation):
    """Function to detect local maximum patches

    - takes contour data of entire event in an anrray with the three orientations
    - creates result class and adds the local max patch and the respective orientation
    - result class is written to list and passed"""
    local_max_patches = []
    data_1 = contour_data
    data_2 = contour_data
    for level_observed in data_1[5:]:
        for patch in level_observed:
            local_max_patch = True
            for level_others in data_2[5:]:
                for patch2 in level_others:
                    if patch2 is not patch:
                        if mpath.Path(patch.contour_coordinates).contains_point(patch2.contour_coordinates[0]):
                            # print(str(patch.level) + " contains " + str(patch2.level))
                            local_max_patch = False

            if local_max_patch:
                patch_class = total_event_reconstruction.RecoPointClass()
                patch_class.contour_data = patch
                patch_class.orientation_index = orientation
                local_max_patches.append(patch_class)
    return local_max_patches


def toplevel_check(top_level_patches, contour_data_total, orientation):
    """Function to iterate through previously detected local maximum patches and passes them to check_function"""
    real_patches = []
    for patch_index, patch in enumerate(top_level_patches):
        if is_real_toplevel_patch(patch.contour_data, contour_data_total):
            real_patches.append(patch)
    return real_patches


def is_real_toplevel_patch(patch, contour_data):
    """checks if patches are interesting for entry-/exit-point search"""
    patch_is_real_top = True

    """"check if patch excees given level threshold"""
    if patch.level < 3:
        patch_is_real_top = False
    """check if found patch has a neighbour patch at same level which is contained by the same contour level below.
    If so -> discard"""
    for patch_level_below in contour_data[patch.level-1]:
        if mpath.Path(patch_level_below.contour_coordinates).contains_point(patch.contour_coordinates[0]):
            for neighbour_patch in contour_data[patch.level]:
                if mpath.Path(patch_level_below.contour_coordinates).contains_point(neighbour_patch.contour_coordinates[0]):
                    if neighbour_patch is not patch:
                        patch_is_real_top = False
    """Same as above, but should find patches which are surrounded by a contour of same level
    Actually makes no sense to me. Maybe there was a case when this was necessary"""
    for patch_level_same in contour_data[patch.level]:
        if mpath.Path(patch_level_same.contour_coordinates).contains_point(patch.contour_coordinates[0]):
            for neighbour_patch in contour_data[patch.level]:
                if mpath.Path(patch_level_same.contour_coordinates).contains_point(neighbour_patch.contour_coordinates[0]):
                    if neighbour_patch is not patch:
                        patch_is_real_top = False
    return patch_is_real_top


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

            if max_frame - min_frame < 2 and minimum < 100 and maximum - minimum > 200:
                point.frame = max_frame
                points_of_orientation.append(point)
        better_reco_points.append(points_of_orientation)
    return better_reco_points


def point_merger(reco_points, merge_radius):
    for point1 in reco_points:
        for point2 in reco_points:
            frame_diff = abs(point1.frame - point2.frame)
            distance = math.sqrt(pow(point1.x_coordinate_deg-point2.x_coordinate_deg, 2)
                                 + pow(point1.y_coordinate_deg-point2.y_coordinate_deg, 2))

            if point1 is not point2 and frame_diff < 2 and distance < merge_radius:
                reco_points.remove(point2)
