"""Allocate the entry and exit points after finding points by looking at entire event and reconstructing the time
 frames. Here we use difference timeframes"""
import statusAlert, reconstructionAlg, PointVecDist
import numpy as np
import matplotlib.path as mpath
import math
import PointVecDist
import data


def allocate_tracks_to_points(reco_points):
    pairing_matrix = create_direction_matrix(reco_points)
    pairing_vector = allocator(pairing_matrix)
    track_array = result_class_writer(reco_points, pairing_vector)
    return track_array


def result_class_writer(reco_points, pairing_vector):
    return_array = []
    # print(reco_points)
    for pair_index, pair in enumerate(pairing_vector):
        TrackC = data.TrackClass()
        reco_points[pair.index_point_1].track = pair_index
        reco_points[pair.index_point_2].track = pair_index
        TrackC.index = pair_index
        if reco_points[pair.index_point_1].frame < reco_points[pair.index_point_2].frame:
            TrackC.entry_point = reco_points[pair.index_point_1]
            TrackC.exit_point = reco_points[pair.index_point_2]
        else:
            TrackC.entry_point = reco_points[pair.index_point_2]
            TrackC.exit_point = reco_points[pair.index_point_1]

        TrackC.distance_track_to_center = PointVecDist.calc_track_dist_to_center(TrackC.entry_point.D3_vector,
                                                                                 TrackC.exit_point.D3_vector)

        return_array.append(TrackC)
    return return_array


def create_direction_matrix(reco_points):
    pairing_matrix = [[0] * len(reco_points) for _ in range(len(reco_points))]
    for point1_index, point1 in enumerate(reco_points):
        for point2_index, point2 in enumerate(reco_points):
            distance = PointVecDist.VectorLength(point2.real_x-point1.real_x,
                                                 point2.real_y-point1.real_y,
                                                 point2.real_z-point1.real_z)

            if distance > 0:
                direction = [(point2.real_x-point1.real_x)/distance,
                             (point2.real_y-point1.real_y)/distance,
                             (point2.real_z-point1.real_z)/distance]
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

    parallelism_factor = np.inf
    parallel_pair = 0
    for coll in collection:
        if get_parallelity(coll) < parallelism_factor:
            parallelism_factor = get_parallelity(coll)
            parallel_pair = coll

    return parallel_pair


def parallelism_check(vector_array, accuracy):
    vecs_are_parallel = 1
    for vec1_index, vector1 in enumerate(vector_array):
        for vec2_index, vector2 in enumerate(vector_array):
            if vec2_index > vec1_index:
                if len(vector1.direction) != len(vector2.direction):
                    print("Vectors have different lengths!")
                    return 0
                else:
                    for v1_coord, v2_coord in zip(vector1.direction,vector2.direction):
                        if abs(abs(v1_coord) - abs(v2_coord)) < accuracy:
                            vecs_are_parallel *= 1
                        else:
                            vecs_are_parallel *= 0
    if vecs_are_parallel is 1:
        return True
    else:
        return False


def get_parallelity(vector_array):
    parallelity = np.inf
    for vec1_index, vector1 in enumerate(vector_array):
        for vec2_index, vector2 in enumerate(vector_array):
            if vec2_index > vec1_index:
                parallelity = 0
                for v1_coord, v2_coord in zip(vector1.direction, vector2.direction):
                    parallelity += abs(abs(v1_coord) - abs(v2_coord))**2
                parallelity = math.sqrt(parallelity)
    return parallelity

