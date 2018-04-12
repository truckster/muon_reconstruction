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
    for level in contour_data_total:
        for patch in level:
            local_max_patch = True
            for patch2 in level:
                if (mpath.Path(np.asarray(patch.contour_coordinates)).contains_points([patch2.contour_coordinates[0]])
                and patch is not patch2):
                    local_max_patch = False
            if local_max_patch:
                local_max_patches.append(patch.center)

    print(local_max_patches)


# nppath = np.asarray(iso_level_contour_array.vertices)
# path = mpath.Path(nppath)
# for next_level_center_count, next_level_center_point in enumerate(next_contour_level_center_points):
#     isit = path.contains_points(next_contour_level_center_points)