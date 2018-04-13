import reconstructionAlg, statusAlert
import numpy as np
import matplotlib.path as mpath


def entry_exit_detector(pmt_position_class, snippet_class):
    for snippet_index, pmt_array in enumerate(snippet_class.time_snippets):
        statusAlert.processStatus("Total event: ")
        contour_data = reconstructionAlg.contour_data_reader(pmt_position_class, pmt_array)
        standalone_contour_lines(contour_data)


def standalone_contour_lines(contour_data_total):
    local_max_patches = []
    for level_observed in contour_data_total:
        for patch in level_observed:
            local_max_patch = True
            for level_others in contour_data_total:
                for patch2 in level_others:
                    if patch2.level > patch.level:
                        if mpath.Path(np.asarray(patch.contour_coordinates)).\
                                contains_path(mpath.Path(np.asarray(patch2.contour_coordinates))):
                            print(str(patch.level) + " contains " + str(patch2.level))
                            local_max_patch = False

            if local_max_patch:
                local_max_patches.append(patch.level)
                local_max_patches.append(patch.center)
    print(local_max_patches)
