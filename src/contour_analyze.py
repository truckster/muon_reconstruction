class ContourData:
    def __init__(self):  # this method creates the class object.
        self.height = 0
        self.iso_hit_patches = 0
        self.centers = []
        self.extents = []


def get_contour_data(contour_output):
    levels = contour_output.levels
    return_data = []
    for level in range(len(levels)):
        patch_data = ContourData()
        patch_data.height = levels[level]
        patch_data.iso_hit_patches = len(contour_output.collections[level].get_paths())
        calc_contour_center(contour_output, level, patch_data)
        calc_contour_extent(contour_output, level, patch_data)
        return_data.append(patch_data)

    return return_data


def calc_contour_center(contour_output, level, return_class):
    for contour in contour_output.collections[level].get_paths():
        x_vals_sum = 0
        y_vals_sum = 0
        for coordinates in contour.vertices:
            x_vals_sum += coordinates[0]
            y_vals_sum += coordinates[1]

        x_coord = x_vals_sum/(len(contour.vertices)-1)
        y_coord = y_vals_sum/(len(contour.vertices)-1)

        return_class.centers.append([x_coord, y_coord])


def calc_contour_extent(contour_output, level, return_class):
    for contour in contour_output.collections[level].get_paths():
        area = 0
        for pair in range(len(contour.vertices[:-1])):
            x0 = contour.vertices[pair][0]
            y0 = contour.vertices[pair][1]

            dx = contour.vertices[pair+1][0] - contour.vertices[pair][0]
            dy = contour.vertices[pair+1][1] - contour.vertices[pair][1]

            area += abs(0.5*(y0*dx - x0*dy))

        return_class.extents.append(area)
