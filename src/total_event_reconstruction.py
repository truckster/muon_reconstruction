import reconstructionAlg, statusAlert
import numpy as np
import matplotlib.path as mpath
import math


def entry_exit_detector(contour_data, number_contour_level):
    statusAlert.processStatus("Total event: ")
    top_levels = standalone_contour_lines(contour_data)
    real_top_patches = toplevel_check(top_levels, contour_data)

    # contour_data_DPhi = reconstructionAlg.contour_data_reader(pmt_position_class, pmt_array, number_contour_level,
    #                                                           math.pi, 0)
    # contour_data_DTheta = reconstructionAlg.contour_data_reader(pmt_position_class, pmt_array, number_contour_level, 0,
    #                                                             math.pi / 2)
    # top_levels_DPhi = standalone_contour_lines(contour_data_DPhi)
    # top_levels_DTheta = standalone_contour_lines(contour_data_DTheta)
    # real_top_patches_DPhi = toplevel_check(top_levels_DPhi, contour_data_DPhi)
    # real_top_patches_DTheta = toplevel_check(top_levels_DTheta, contour_data_DTheta)

    print(len(real_top_patches))
    # print(len(real_top_patches_DPhi))
    # print(len(real_top_patches_DTheta))

    return real_top_patches


def standalone_contour_lines(contour_data_total):
    data_1 = contour_data_total
    data_2 = contour_data_total
    print(data_2)
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
                local_max_patches.append(patch)
    return local_max_patches


def toplevel_check(top_level_patches, contour_data_total):
    real_patches = []
    for patch_index, patch in enumerate(top_level_patches):
        if is_real_toplevel_patch(patch, contour_data_total):
            real_patches.append(patch)
    return real_patches


def reco_result_writer(output_path, result_array, return_array):
    '''Add reconstructed intersection points to file'''
    result_file = open(output_path + "results.txt", 'a')
    result_file.write("----- Reconstructed Values (Total)------" + '\n')

    result_file.write("Found patches: " + str(len(result_array)) + '\n' + '\n')
    positions = []

    for patch in result_array:
        reco_point = reconstructionAlg.RecoPointClass()
        result_file.write("Phi: " + str(patch.center[0]/math.pi * 180.0) + '\n')
        result_file.write("Theta: " + str(patch.center[1]/math.pi * 180.0) + '\n')
        result_file.write("--------------------------------" + '\n')
        # positions.append([patch.center[0]/math.pi * 180.0, patch.center[1]/math.pi * 180.0])
        reco_point.coordinates = [patch.center[0] / math.pi * 180.0, patch.center[1] / math.pi * 180.0]
        return_array.append(reco_point)

    result_file.write("End of event" + '\n' + '\n')
    result_file.close()

    return return_array


def is_real_toplevel_patch(patch, contour_data):
    patch_is_real_top = True
    # patch_degree = [coord /math.pi * 180.0 for coord in patch.center]
    # print(patch.level)
    # print(patch_degree)
    if patch.level < 2:
        patch_is_real_top = False
    for patch_level_below in contour_data[patch.level-1]:
        if mpath.Path(patch_level_below.contour_coordinates).contains_point(patch.contour_coordinates[0]):
            # print(patch_level_below.level)
            # patch_degree = [coord / math.pi * 180.0 for coord in patch_level_below.center]
            # print(patch_degree)
            for neighbour_patch in contour_data[patch.level]:
                if mpath.Path(patch_level_below.contour_coordinates).contains_point(neighbour_patch.contour_coordinates[0]) and patch.level < 2:
                    if neighbour_patch is not patch:
                        patch_is_real_top = False

    for patch_level_same in contour_data[patch.level]:
        if mpath.Path(patch_level_same.contour_coordinates).contains_point(patch.contour_coordinates[0]):
            for neighbour_patch in contour_data[patch.level]:
                if mpath.Path(patch_level_same.contour_coordinates).contains_point(neighbour_patch.contour_coordinates[0]):
                    if neighbour_patch is not patch:
                        patch_is_real_top = False

    #TODO: Use level check here. Non-close paths of high level are okay.
    # if not compare_points(patch.contour_coordinates[0], patch.contour_coordinates[-1]):
    #     patch_is_real_top = False

    # if patch.level < 2:
    #     patch_is_real_top = False

    return patch_is_real_top


def compare_points(a, b):
    if a[0] == b[0] and a[1] == b[1]:
        return True
    else:
        return False
