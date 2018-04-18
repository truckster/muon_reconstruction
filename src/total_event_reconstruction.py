import reconstructionAlg, statusAlert
import numpy as np
import matplotlib.path as mpath


def entry_exit_detector(pmt_position_class, snippet_class):
    pmt_array = snippet_class.time_snippets[0]
    statusAlert.processStatus("Total event: ")
    contour_data = reconstructionAlg.contour_data_reader(pmt_position_class, pmt_array)
    found_points = standalone_contour_lines(contour_data)
    return found_points


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

            if local_max_patch and patch.level > data_1[-6][0].level:
                local_max_patches.append(patch.center)
    return local_max_patches


def reco_result_writer(output_path, result_array):
    '''Add reconstructed intersection points to file'''
    result_file = open(output_path + "results.txt", 'a')
    result_file.write("----- Reconstructed Values (Total)------" + '\n')

    for point in result_array:
        result_file.write("Phi: " + str(point[0]) + '\n')
        result_file.write("Theta: " + str(point[1]) + '\n')
        result_file.write("--------------------------------" + '\n')

    result_file.write("End of event" + '\n' + '\n')
    result_file.close()
