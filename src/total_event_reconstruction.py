import reconstructionAlg, statusAlert
import numpy as np
import matplotlib.path as mpath
import math
import data


def entry_exit_detector(contour_data):
    real_top_patches = []
    for orientation_index, data in enumerate(contour_data):
        statusAlert.processStatus("Searching entry and exit points: ")

        top_levels = standalone_contour_lines(data[0], orientation_index)
        top_levels = toplevel_check(top_levels, data[0], orientation_index)
        real_top_patches.append(top_levels)

    return real_top_patches


def standalone_contour_lines(contour_data_total, orientation):
    """Function to detect local maximum patches

    - takes contour data of entire event in an anrray with the three orientations
    - creates result class and adds the local max patch and the respective orientation
    - result class is written to list and passed"""
    data_1 = contour_data_total
    data_2 = contour_data_total
    local_max_patches = []
    for level_observed in data_1:
        for patch in level_observed:
            local_max_patch = True
            for level_others in data_2:
                for patch2 in level_others:
                    if patch2 is not patch:
                        if mpath.Path(patch.contour_coordinates).contains_point(patch2.contour_coordinates[0]):
                            # print(str(patch.level) + " contains " + str(patch2.level))
                            local_max_patch = False

            if local_max_patch:
                patch_class = data.RecoPointClass()
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

    """"check if patch exceeds given level threshold"""
    # if patch.level < 2:
    #     patch_is_real_top = False

    # """check if found patch has a neighbour patch at same level which is contained by the same contour level below.
    # If so -> discard"""
    # for patch_level_below in contour_data[patch.level-2]:
    #     if mpath.Path(patch_level_below.contour_coordinates).contains_point(patch.contour_coordinates[0]):
    #         for neighbour_patch in contour_data[patch.level]:
    #             if mpath.Path(patch_level_below.contour_coordinates).contains_point(neighbour_patch.contour_coordinates[0]):
    #                 if neighbour_patch is not patch:
    #                     patch_is_real_top = False

    # """Same as above, but should find patches which are surrounded by a contour of same level
    # Actually makes no sense to me. Maybe there was a case when this was necessary"""
    # for patch_level_same in contour_data[patch.level]:
    #     if mpath.Path(patch_level_same.contour_coordinates).contains_point(patch.contour_coordinates[0]):
    #         for neighbour_patch in contour_data[patch.level]:
    #             if mpath.Path(patch_level_same.contour_coordinates).contains_point(neighbour_patch.contour_coordinates[0]):
    #                 if neighbour_patch is not patch:
    #                     patch_is_real_top = False

    """Check for non-closed contour paths"""
    if not compare_points(patch.contour_coordinates[0], patch.contour_coordinates[-1]):
        patch_is_real_top = False

    # if patch.level < 2:
    #     patch_is_real_top = False

    return patch_is_real_top


def reco_result_writer(output_path, result_array):
    '''Add reconstructed intersection points to file'''
    result_file = open(output_path + "results.txt", 'a')
    result_file.write("----- Reconstructed Values (Total)------" + '\n')

    result_file.write("Found patches: " + str(len(result_array)) + '\n' + '\n')

    for patch in result_array:
        result_file.write("Frame: " + str(patch.frame) + '\n')
        result_file.write("Phi: " + str(patch.x_coordinate_deg) + '\n')
        result_file.write("Theta: " + str(patch.y_coordinate_deg) + '\n')
        result_file.write("Track: " + str(patch.track)+ '\n')
        result_file.write("--------------------------------" + '\n')

    result_file.write("End of event" + '\n' + '\n')
    result_file.close()


def compare_points(a, b):
    if a[0] == b[0] and a[1] == b[1]:
        return True
    else:
        return False
