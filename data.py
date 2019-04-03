class MuonTrackTruth:
    def __init__(self):
        self.track_index = 0
        self.entry_point = 0
        self.exit_point = 0
        self.entry_time = 0
        self.exit_time = 0
        self.distance_track_to_center = 0
        self.is_showering = False


class RecoPointClass:
    def __init__(self):
        self.frame = 0
        self.x_coordinate_deg = 0
        self.y_coordinate_deg = 0
        self.x_coordinate_rad = 0
        self.y_coordinate_rad = 0
        self.real_x = 0
        self.real_y = 0
        self.real_z = 0
        self.D3_vector = 0
        self.contour_data = 0
        self.orientation_index = 0
        self.track = None
        self.closest_pmt = 0


class TrackClass:
    def __init__(self):
        self.index = 0
        self.entry_point = 0
        self.exit_point = 0
        self.distance_track_to_center = 0


class RecoPerformanceCheckTotal:
    def __init__(self):
        self.point_reco_accuracy_list = []
        self.track_reco_accuracy_list_D = []
        self.track_reco_accuracy_list_phi = []
        self.tracks_per_event_list = []
        self.found_point_list = []
        self.found_points_for_good_reco = []
        self.unmerged = []
        self.merged = []
        self.track_reco_accuracy_list_uncorr = []
        self.mc_track_distance_list = []
        self.distance_correction = []
        self.reconstructed_track_distance_list = []
        self.number_of_points_showering = []
        self.number_of_points_not_showering = []
        self.time_difference_list = []
        self.parallelity_list = []


class MuonEventData:
    def __init__(self):
        # save MuonData class here
        self.muon_data = []

        # isotope information
        self.iso_num = 0
        self.iso_name = []
        self.iso_proc = []
        self.iso_pos_x = []
        self.iso_pos_y = []
        self.iso_pos_z = []
        self.iso_time = []


class MuonData:
    def __init__(self):
        self.x_position_init = 0
        self.y_position_init = 0
        self.z_position_init = 0

        self.x_momentum_init = 0
        self.y_momentum_init = 0
        self.z_momentum_init = 0

        self.x_position_final = 0
        self.y_position_final = 0
        self.z_position_final = 0

        self.scint_track_length = 0

        self.init_kinetic_energy = 0
        self.e_loss_rock = 0
        self.e_loss_Veto = 0
        self.e_loss_cd_water = 0
        self.e_loss_acrylic = 0
        self.e_loss_steel = 0
        self.e_loss_scint = 0

