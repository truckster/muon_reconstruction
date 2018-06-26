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


def allocate_points(contour_diffs, intersec_points, intersec_times):
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
    print(pairing_matrix)
    return pairing_matrix


def allocator(pairing_matrix):
    for point1_index, point1 in enumerate(pairing_matrix[:-2]):
        for vector1_index, vector1 in enumerate(point1[point1_index+1:]):
            print(":::::")
            print(vector1)
            print(":::::")
            for point2_index, point2 in enumerate(pairing_matrix[point1_index+1:-1]):
                for vector2_index, vector2 in enumerate(point2[point1_index + point2_index + 2:]):
                    print(vector2)
                    print("-------------")


def parallelism_check(vector1, vector2, accuracy):
    if len(vector1) != len(vector2):
        print("Vectors have different lengths!")
        return 0
    else:
        vecs_are_parallel = 1
        for coordinate in range(len(vector1)):
            if abs(abs(vector1[coordinate]) - abs(vector2[coordinate])) < accuracy:
                vecs_are_parallel *= 1
            else:
                vecs_are_parallel *= 0
    if vecs_are_parallel is 1:
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
