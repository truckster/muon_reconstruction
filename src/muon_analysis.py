import TreeReadFunc, statusAlert, PointVecDist
import time, math, random, os
import ROOT as root
import numpy as np
import matplotlib.pyplot as plt


class MuonIntersecPoint:
    def __init__(self):
        self.event = 0

        self.enters = False
        self.leaves = False

        self.x = 0
        self.y = 0
        self.z = 0

        self.phi = 0
        self.theta = 0
        self.theta2 = 0

        self.phi_hammer_aitoff = 0
        self.theta_hammer_aitoff = 0

        self.intersec_time = 0


def distance_to_center(x, y, z):
    dist = math.sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))
    return dist


def is_muon_in_detector(dist_to_center, radius):
    if dist_to_center < radius:
        muon_in_det = True
    else:
        muon_in_det = False

    return muon_in_det


def calc_muon_detector_intersec_points(data, radius, time_steps):
    """Calculates and gives back the entry and exit points of the muons belonging to the event."""
    statusAlert.processStatus("Calculate entry and exit point of muons in this event from monte carlo truth")
    returnarray = []
    muon_mass = 105.6583745
    c = 299792548000

    for index, muon_event in enumerate(data):
        muon_event_array = []
        muon_in = MuonIntersecPoint()
        muon_out = MuonIntersecPoint()
        """Track-Sphere intersection points"""
        point_array = PointVecDist.calc_line_sphere_intersect_points(muon_event, radius)
        if len(point_array) == 2:
            muon_in.event = index
            muon_out.event = index

            time_current = 0

            total_momentum = math.sqrt(pow(muon_event.x_momentum_init, 2)
                                       + pow(muon_event.y_momentum_init, 2)
                                       + pow(muon_event.z_momentum_init, 2))

            x_velocity = muon_event.x_momentum_init / math.sqrt(pow(muon_mass, 2) + total_momentum**2) * c
            y_velocity = muon_event.y_momentum_init / math.sqrt(pow(muon_mass, 2) + total_momentum**2) * c
            z_velocity = muon_event.z_momentum_init / math.sqrt(pow(muon_mass, 2) + total_momentum**2) * c

            x_position_current = muon_event.x_position_init
            y_position_current = muon_event.y_position_init
            z_position_current = muon_event.z_position_init

            distance_to_det_center = distance_to_center(x_position_current, y_position_current, z_position_current)
            while z_position_current > -23000:
                muon_in_detector = is_muon_in_detector(distance_to_det_center, radius)
                x_position_current = x_position_current + time_steps * x_velocity
                y_position_current = y_position_current + time_steps * y_velocity
                z_position_current = z_position_current + time_steps * z_velocity

                time_current += time_steps
                distance_to_det_center = distance_to_center(x_position_current, y_position_current, z_position_current)
                muon_in_detector_second = is_muon_in_detector(distance_to_det_center, radius)

                if muon_in_detector is False and muon_in_detector_second is True:
                    muon_in.enters = True
                    if point_array[0].z > point_array[1].z:
                        muon_in.x = point_array[0].x
                        muon_in.y = point_array[0].y
                        muon_in.z = point_array[0].z
                    else:
                        muon_in.x = point_array[1].x
                        muon_in.y = point_array[1].y
                        muon_in.z = point_array[1].z

                    muon_in.phi = math.atan2(muon_in.y, muon_in.x)
                    muon_in.theta = math.acos(-muon_in.z/radius)
                    muon_in.theta2 = math.acos(-muon_in.z/radius) - (math.pi/2)

                    muon_in.phi_hammer_aitoff = (math.sqrt(8) * math.cos(muon_in.theta) * math.sin(muon_in.phi / 2)) / \
                                                (math.sqrt(1 + math.cos(muon_in.theta) * math.cos(muon_in.phi / 2)))

                    muon_in.theta_hammer_aitoff = (math.sqrt(2) * math.sin(muon_in.theta)) / \
                                                  (math.sqrt(1 + math.cos(muon_in.theta) * math.cos(muon_in.phi / 2)))

                    muon_in.intersec_time = time_current

                    muon_event_array.append(muon_in)

                if muon_in_detector is True and muon_in_detector_second is False and is_muon_stopping(muon_event) is False:
                    muon_out.leaves = True
                    muon_out.x = x_position_current
                    muon_out.y = y_position_current
                    muon_out.z = z_position_current

                    if point_array[1].z > point_array[0].z:
                        muon_out.x = point_array[0].x
                        muon_out.y = point_array[0].y
                        muon_out.z = point_array[0].z
                    else:
                        muon_out.x = point_array[1].x
                        muon_out.y = point_array[1].y
                        muon_out.z = point_array[1].z

                    muon_out.phi = math.atan2(muon_out.y, muon_out.x)
                    muon_out.theta = math.acos(-muon_out.z/radius)
                    muon_out.theta2 = math.asin(muon_out.z/radius)

                    muon_out.phi_hammer_aitoff = (math.sqrt(8) * math.cos(muon_out.theta) * math.sin(muon_out.phi / 2)) / \
                              (math.sqrt(1 + math.cos(muon_out.theta) * math.cos(muon_out.phi / 2)))

                    muon_out.theta_hammer_aitoff = (math.sqrt(2) * math.sin(muon_out.theta)) / \
                                (math.sqrt(1 + math.cos(muon_out.theta) * math.cos(muon_out.phi / 2)))

                    muon_out.intersec_time = time_current

                    muon_event_array.append(muon_out)
            returnarray.append(muon_event_array)
        else:
            print("Bad track (not 2 intersection points found)")
    return returnarray


def is_muon_stopping(muon_event):
    x_position_current = muon_event.x_position_final
    y_position_current = muon_event.y_position_final
    z_position_current = muon_event.z_position_final

    if PointVecDist.VectorLength(x_position_current, y_position_current, z_position_current) > 17600:
        return False
    else:
        return True


def muon_is_showering_truth(muon_truth):
    if muon_truth[0].iso_num > 0:
        print "Shower power"
        return True
    else:
        print "Good boy"
        return False