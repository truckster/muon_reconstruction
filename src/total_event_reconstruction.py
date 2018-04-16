import reconstructionAlg, statusAlert
import numpy as np
import matplotlib.path as mpath


def entry_exit_detector(pmt_position_class, snippet_class):
    for snippet_index, pmt_array in enumerate(snippet_class.time_snippets):
        statusAlert.processStatus("Total event: ")
        contour_data = reconstructionAlg.contour_data_reader(pmt_position_class, pmt_array)
        standalone_contour_lines(contour_data)


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
                            print(str(patch.level) + " contains " + str(patch2.level))
                            local_max_patch = False

            if local_max_patch:
                local_max_patches.append(patch.level)
                local_max_patches.append(patch.center)
    print(local_max_patches)
