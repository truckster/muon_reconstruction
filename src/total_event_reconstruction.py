import reconstructionAlg, statusAlert
import numpy as np
import matplotlib.path as mpath


def entry_exit_detector(pmt_position_class, snippet_class):
    pmt_array = snippet_class.time_snippets[0]
    statusAlert.processStatus("Total event: ")
    contour_data = reconstructionAlg.contour_data_reader(pmt_position_class, pmt_array)
    top_levels = standalone_contour_lines(contour_data)
    real_top_patches = toplevel_check(top_levels, contour_data)

    return real_top_patches


def standalone_contour_lines(contour_data_total):
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
                local_max_patches.append(patch)
    return local_max_patches


def toplevel_check(top_level_patches, contour_data_total):
    real_patches = []
    for patch_index, patch in enumerate(top_level_patches):
        if is_real_toplevel_patch(patch, contour_data_total):
            real_patches.append(patch)
    return real_patches


def reco_result_writer(output_path, result_array):
    '''Add reconstructed intersection points to file'''
    result_file = open(output_path + "results.txt", 'a')
    result_file.write("----- Reconstructed Values (Total)------" + '\n')

    for patch in result_array:
        result_file.write("Phi: " + str(patch.center[0]) + '\n')
        result_file.write("Theta: " + str(patch.center[1]) + '\n')
        result_file.write("--------------------------------" + '\n')

    result_file.write("End of event" + '\n' + '\n')
    result_file.close()


def is_real_toplevel_patch(patch, contour_data):
    patch_is_real_top = True
    for patch_level_below in contour_data[patch.level-1]:
        if mpath.Path(patch_level_below.contour_coordinates).contains_point(patch.contour_coordinates[0]):
            for neighbour_patch in contour_data[patch.level]:
                if mpath.Path(patch_level_below.contour_coordinates).contains_point(neighbour_patch.contour_coordinates[0]):
                    if neighbour_patch is not patch:
                        patch_is_real_top = False

    return patch_is_real_top


