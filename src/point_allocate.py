"""Allocate the entry and exit points after finding points by looking at entire event and reconstructing the time
 frames. Here we use difference timeframes"""
import statusAlert, reconstructionAlg, PointVecDist
import numpy as np
import matplotlib.path as mpath
import math
import PointVecDist


class DirectionReco:
    def __init__(self):
        self.frame = 0
        self.reference_point = []
        self.direction = []


def allocate_points(contour_diffs, intersec_points):
    for frame, contour in enumerate(contour_diffs):
        top_level = contour[-1]
        low_level = contour[0]
        for patch_low in low_level:
            lengths = []
            angles = []
            if patch_low.height < 0:
                for patch_high in top_level:
                    lengths.append(math.sqrt(pow(patch_high.center[0] - patch_low.center[0], 2)
                                             + pow(patch_high.center[1] - patch_low.center[1], 2)))
                    angles.append(np.arctan((patch_high.center[1] - patch_low.center[1])
                                            / patch_high.center[0] - patch_low.center[0]))
                    same_area = have_same_surrounder(contour, patch_low, patch_high)
                if len(lengths) > 0 and min(lengths) < 0.25 and same_area:
                    print("Frame: " + str(frame + 1))
                    print("Length: " + str(min(lengths)))
                    mindex = np.argmax(lengths)
                    print("Angle: " + str(angles[mindex]/math.pi * 180.0))


def have_same_surrounder(contour_data, patch1, patch2):
    they_do = False
    for level in contour_data:
        for patch in level:
            if mpath.Path(patch.contour_coordinates).contains_point(patch2.contour_coordinates[0]) and \
                   mpath.Path(patch.contour_coordinates).contains_point(patch1.contour_coordinates[0]):
                    they_do = True

    return they_do


def allocate_tracks_to_points(reco_points):
    pairing_matrix = create_direction_matrix(reco_points)
    pairing_vector = allocator(pairing_matrix)
    result_class_writer(reco_points, pairing_vector)
    return pairing_vector


def result_class_writer(reco_points, pairing_vector):
    if len(pairing_vector) > 1:
        print("Several solutions found!!!")
    else:
        for result in pairing_vector:
            for pair_index, pair in enumerate(result):
                reco_points[pair.index_point_1].track = pair_index
                reco_points[pair.index_point_2].track = pair_index


def create_direction_matrix(reco_points):
    pairing_matrix = [[0] * len(reco_points) for _ in range(len(reco_points))]
    r = 17600
    for point1_index, point1 in enumerate(reco_points):
        x1 = r * math.cos(point1.y_coordinate_rad) * math.cos(point1.x_coordinate_rad)
        y1 = r * math.cos(point1.y_coordinate_rad) * math.sin(point1.x_coordinate_rad)
        z1 = r * math.sin(point1.y_coordinate_rad)

        for point2_index, point2 in enumerate(reco_points):
            x2 = r * math.cos(point2.y_coordinate_rad) * math.cos(point2.x_coordinate_rad)
            y2 = r * math.cos(point2.y_coordinate_rad) * math.sin(point2.x_coordinate_rad)
            z2 = r * math.sin(point2.y_coordinate_rad)

            distance = PointVecDist.VectorLength(x2-x1, y2-y1, z2-z1)

            if distance > 0:
                direction = [(x2-x1)/distance, (y2-y1)/distance, (z2-z1)/distance]
                pairing_matrix[point1_index][point2_index] = direction
    return pairing_matrix


class DirectionClass:
    def __init__(self):
        self.direction = []
        self.index_point_1 = 0
        self.index_point_2 = 0


def allocator(pairing_matrix):
    return_array = []
    collection = []
    for row_index, row in enumerate(pairing_matrix):
        for col_index, col in enumerate(pairing_matrix):
            if col_index > row_index:
                direction_start = DirectionClass()
                direction_start.direction = row[col_index]
                direction_start.index_point_1 = row_index
                direction_start.index_point_2 = col_index

                for alloc_round in range((len(pairing_matrix) / 2) - 1):
                    for row_index_comp, row_comp in enumerate(pairing_matrix):
                        for col_index_comp, col_comp in enumerate(pairing_matrix):
                            direction_array = [direction_start]
                            if col_index_comp > row_index_comp > row_index:
                                direction_stop = DirectionClass()
                                direction_stop.direction = row_comp[col_index_comp]
                                direction_stop.index_point_1 = row_index_comp
                                direction_stop.index_point_2 = col_index_comp
                                direction_array.append(direction_stop)

                                collection.append(direction_array)
    for coll in collection:
        if parallelism_check(coll, 0.1) and len(coll) is len(pairing_matrix)/2:
            return_array.append(coll)
    return return_array


# def allocator(pairing_matrix):
#     # pair_vector = len(pairing_matrix)*[None]
#     pair_vector = [[] for _ in range(len(pairing_matrix) / 2)]
#     first_alloc = True
#     for point1_index, point1 in enumerate(pairing_matrix):
#         for vector1_index, vector1 in enumerate(point1):
#             if not is_middle_diagonal_element(point1_index, vector1_index):
#                 for point2_index, point2 in enumerate(pairing_matrix):
#                     if point2_index > point1_index:
#                         for vector2_index, vector2 in enumerate(point2):
#                             if not is_middle_diagonal_element(point2_index, vector2_index):
#                                 if not are_mirror_matrix_element(point1_index, vector1_index, point2_index, vector2_index)\
#                                         and not are_same_matrix_element(point1_index, vector1_index, point2_index, vector2_index):
#                                     if parallelism_check(vector1, vector2, 0.1):
#                                         print(vector1)
#                                         print(vector2)
#                                         print("1: " + str(point1_index) + ":" + str(vector1_index))
#                                         print("2: " + str(point2_index) + ":" + str(vector2_index))
#                                         if first_alloc:
#                                             pair_vector[point1_index] = vector1_index
#                                             pair_vector[point2_index] = vector2_index
#                                             pair_vector[vector1_index] = point1_index
#                                             pair_vector[vector2_index] = point2_index
#                                         else:
#                                             is_alloc_still_correct(pair_vector)
#                                         first_alloc = False
#     print(pair_vector)
#     return pair_vector


def parallelism_check(vector_array, accuracy):
    vecs_are_parallel = 1
    for vec1_index, vector1 in enumerate(vector_array):
        for vec2_index, vector2 in enumerate(vector_array):
            if vec2_index > vec1_index:
                if len(vector1.direction) != len(vector2.direction):
                    print("Vectors have different lengths!")
                    return 0
                else:
                    for coordinate in range(len(vector1.direction)):
                        if abs(abs(vector1.direction[coordinate]) - abs(vector2.direction[coordinate])) < accuracy:
                            vecs_are_parallel *= 1
                        else:
                            vecs_are_parallel *= 0
    if vecs_are_parallel is 1:
        return True
    else:
        return False


def is_alloc_still_correct(pair_vector):
    for index, element in enumerate(pair_vector):
        if index is pair_vector[element]:
            pass
        else:
            print("Problem in point allocation")


def are_mirror_matrix_element(column1, row1, column2, row2):
    if column1 is row2 and column2 is row1:
        return True
    else:
        return False


def are_same_matrix_element(column1, row1, column2, row2):
    if column1 is column2 and row1 is row2:
        return True
    else:
        return False


def is_middle_diagonal_element(column, row):
    if column is row:
        return True
    else:
        return False


def are_same_vector(vector1, vector2):
    if len(vector1) != len(vector2):
        print("Vectors have different lengths!")
        return 0
    else:
        vecs_are_same = 1
        for coordinate in range(len(vector1)):
            if abs(vector1[coordinate]) == abs(vector2[coordinate]):
                vecs_are_same *= 1
            else:
                vecs_are_same *= 0
    if vecs_are_same is 1:
        return True
    else:
        return False
