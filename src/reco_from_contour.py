import matplotlib.pyplot as plt


class ReconstructIntersecPoint:
    def __init__(self):  # this method creates the class object.
        self.snippet = 0
        self.phi = 0
        self.theta = 0


def peak_compare(peak_array, contour_data):
    return_array = []
    for count, value in enumerate(peak_array):
        if count > 0:
            if value > peak_array[count-1]+30 and value > peak_array[count+1]+30:
                point_class = ReconstructIntersecPoint()
                point_class.snippet = count
                point_class.phi = contour_data[count][-1].centers[0][0]
                point_class.theta = contour_data[count][-1].centers[0][1]

                return_array.append(point_class)

    return return_array


def concentric_level_finder(contour_data, snippet):
    # print(contour_raw.collections[0].get_paths())
    for level in range(len(contour_data)):
        print("------------------------------------------")
        for iso_hit in range(contour_data[level].iso_hit_patches):
            print(contour_data[level].centers[iso_hit])

