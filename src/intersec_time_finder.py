"""Find the entry and exit times after finding points by looking at entire event. Here we use time frames"""
import statusAlert, reconstructionAlg
import numpy as np
import matplotlib.path as mpath


def find_times(pmt_position_class, snippet_class, intersec_points):
    for frame, pmt_array in enumerate(snippet_class.time_snippets):
        statusAlert.processStatus("processing snippet: " + str(frame))
        contour_data = reconstructionAlg.contour_data_reader(pmt_position_class, pmt_array)
        if frame > 0:
            snippet_diff = np.asarray(snippet_class.time_snippets[frame]) \
                           - np.asarray(snippet_class.time_snippets[frame - 1])
            contour_data_diff = reconstructionAlg.contour_data_reader(pmt_position_class, snippet_diff)
        for point in intersec_points:
            current_level = 0
            for level in contour_data:
                for patch in level:
                    if mpath.Path(patch.contour_coordinates).contains_point(point):
                        current_level = patch.height
                        # print(current_level)
                        print(patch.height)



