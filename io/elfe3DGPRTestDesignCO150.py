import numpy as np

class elfe3DGPRTestDesign:
    """
    Class to generate the test design for the 3D GPR simulations.
    The class returns all parameters to run the elfe3DTestWriting class.
    """

    def __init__(
        self, num_points_per_range=11, least_samples_per_wavelength=10,
        num_layers=2, layer_thicknesses=[2.0, 8.0],
        eps_r=[1.0, 3.0, 23.0], sigma=[1e-8, 5e-3, 2.5e-2], mu_r=[1.0, 1.0, 1.0], sigma_m=[0.0, 0.0, 0.0],
        source_type=6, antenna_position=[0.0, 0.0, 0.0], current_direction=1, current=1.0, source_moment=1, num_segments=4, ricker_central_f=500e6,
        num_receivers_inline=3, num_receivers_endfire=3, num_receivers_oblique=3,
        solver_type=2, max_ref_steps=0, max_unknowns=1e7, beta_ref=0.85, accuracy_tol=3e-5, vtk=0, error_est_method=4, ref_strategy=1, output_fields_vtk=1,
        num_air_wavelengths=10.0, num_PML_layers=2, PML_layer_thickness=0.5, PML_type="log", PML_theory=0, PML_decay_type=0,
        x_e=[], y_e=[], z_e=[], experiment_name="test", m=4, box_x=4.8 , s_f=50, box_present=True, bh_f=1.0,
):
        self.c = 3e8
        self.least_samples_per_wavelength = least_samples_per_wavelength

        # Assembling dictionaries for structured access
        self.frequency_range = {
            'num_points_per_range': num_points_per_range,
        }

        self.model_layers = {
            'num_layers': int(num_layers),
            'layer_thicknesses': [float(thickness) for thickness in layer_thicknesses]
        }

        self.regions = {
            'num_regions': num_layers+1,
            'eps_r': eps_r,
            'sigma': sigma,
            'mu_r': mu_r,
            'sigma_m': sigma_m,
            'rho': None  # Assuming rho is not provided in the parameters
        }
        self.regions['rho'] = [1/(self.regions['sigma'][i]) for i in range(self.regions['num_regions'])]

        self.source_antenna = {
            'source_type': source_type,
            'antenna_position': antenna_position,
            'current_direction': current_direction,
            'current': current,
            'source_moment': source_moment,
            'num_segments': num_segments,
            'ricker_central_f': ricker_central_f,
            'm': m,
            'box_present': box_present,
        }

        self.receiver_antenna = {
            'num_receivers_inline': num_receivers_inline,
            'num_receivers_endfire': num_receivers_endfire,
            'num_receivers_oblique': num_receivers_oblique
        }

        self.solver = {
            'solver': solver_type,
            'maxRefSteps': max_ref_steps,
            'maxUnknowns': max_unknowns,
            'betaRef': beta_ref,
            'accuracyTol': accuracy_tol,
            'vtk': vtk,
            'errorEst_method': error_est_method,
            'refStrategy': ref_strategy,
            'output_fields_vtk': output_fields_vtk,
        }

        # Domain parameters
        self.domain = {'num_air_wavelengths': num_air_wavelengths}
        self.domain['x_e'] = x_e
        self.domain['y_e'] = y_e
        self.domain['z_e'] = z_e

        # PML parameters
        self.PML = {
            'num_layers': num_PML_layers,
            'layer_thickness': PML_layer_thickness,
            'PML_theory': PML_theory,
            'PML_decay_type': PML_decay_type,
            'eps_r': eps_r,
            'sigma': sigma,
            'mu_r': mu_r,
            'sigma_m': sigma_m,
            'PML_type': PML_type
        }

        # Lambda functions to calculate wavelength and frequency
        self.wavelength_from_f = lambda f, eps_r: float(self.c/(f*np.sqrt(eps_r))) # in meters
        self.frequency_from_lambda = lambda wl, eps_r: float(self.c/(wl*np.sqrt(eps_r))) # in Hz
        
        # Frequency list for source antenna
        if self.frequency_range['num_points_per_range'] > 1:
            self.source_antenna['f_list'] = [i*self.source_antenna['ricker_central_f']*(6/self.frequency_range['num_points_per_range']) for i in range(1,self.frequency_range['num_points_per_range']+1,1)]
        else:
            self.source_antenna['f_list'] = [self.source_antenna['ricker_central_f']]

        self.frequency_range['f_low'] = self.source_antenna['f_list'][0]
        self.frequency_range['f_high'] = self.source_antenna['f_list'][-1]

        # Calculate largest and smallest wavelengths for regions
        self.regions['largest_wavelengths'] = [self.wavelength_from_f(self.frequency_range['f_low'],e_r) for e_r in self.regions['eps_r']]
        # self.regions['smallest_wavelengths'] = [self.wavelength_from_f(self.frequency_range['f_high'],e_r) for e_r in self.regions['eps_r']]
        # self.regions['central_wavelengths'] = [self.wavelength_from_f(self.source_antenna['ricker_central_f'],e_r) for e_r in self.regions['eps_r']]
        
        # Calculate max element edges and volumes
        self.regions['max_element_edges'] = (np.array(self.regions['largest_wavelengths']) / self.least_samples_per_wavelength).tolist()
        print(f"odepths: {[d/4.0 for d in self.regions['max_element_edges']]}")
        self.regions['max_element_volume'] = [
            edge**3 / (6 * np.sqrt(2)) for edge in self.regions['max_element_edges']
        ]
        self.regions['max_element_volume'] = [
            f"{vol:.4e}" for vol in self.regions['max_element_volume']
        ]

        # Update domain boundaries
        self.domain.update({
            'x_min': self.domain['x_e'][0],
            'x_max': self.domain['x_e'][1],
            'y_min': self.domain['y_e'][0],
            'y_max': self.domain['y_e'][1],
            'z_min': self.domain['z_e'][0],
            'z_max': self.domain['z_e'][1]
        })

        # Update model boundaries
        self.model_layers.update({
            'x_min': self.domain['x_e'][0],
            'x_max': self.domain['x_e'][1],
            'y_min': self.domain['y_e'][0],
            'y_max': self.domain['y_e'][1],
            'z': [0.0] + [np.sum(self.model_layers['layer_thicknesses'][:i], dtype=np.float64) * -1 for i in range(1,self.model_layers['num_layers'],1)]
        })
        
        # Source antenna extents
        self.source_antenna['length'] = self.regions['largest_wavelengths'][0] / s_f if num_layers == 0 else self.regions['largest_wavelengths'][1] / s_f
        print(f"Source antenna length: {self.source_antenna['length']} m")
        self.source_antenna['x_extents'] = [self.source_antenna['antenna_position'][0]-self.source_antenna['length']/2, 
                                            self.source_antenna['antenna_position'][0]+self.source_antenna['length']/2]
        self.source_antenna['y_extents'] = [self.source_antenna['antenna_position'][1]] * 2
        self.source_antenna['z_extents'] = [self.source_antenna['antenna_position'][2]] * 2

        self.source_antenna['x_disc'] = [x for x in np.linspace(self.source_antenna['x_extents'][0], self.source_antenna['x_extents'][1], self.source_antenna['num_segments']+1)]
        self.source_antenna['y_disc'] = [y for y in np.linspace(self.source_antenna['y_extents'][0], self.source_antenna['y_extents'][1], self.source_antenna['num_segments']+1)]
        self.source_antenna['z_disc'] = [z for z in np.linspace(self.source_antenna['z_extents'][0], self.source_antenna['z_extents'][1], self.source_antenna['num_segments']+1)]

        # Source antenna box - lambda/4 around the antenna and receivers
        self.source_antenna['box_x'] = [[self.source_antenna['antenna_position'][0] - self.regions['largest_wavelengths'][0] / (4*bh_f),
                                         box_x if num_layers == 0 else box_x[0]]]
                                        #  self.source_antenna['antenna_position'][0] + self.regions['largest_wavelengths'][0] / 4]]
        self.source_antenna['box_y'] = [[self.source_antenna['antenna_position'][1] - self.regions['largest_wavelengths'][0] / (4*bh_f),
                                         self.source_antenna['antenna_position'][1] + self.regions['largest_wavelengths'][0] / (4*bh_f)]]
        self.source_antenna['box_z'] = [[self.source_antenna['antenna_position'][2] + self.regions['largest_wavelengths'][0] / (4*bh_f), 0.0]]

        if num_layers > 0:
            self.source_antenna['box_x'].append(
                [self.source_antenna['antenna_position'][0] - self.regions['largest_wavelengths'][1] / (4*bh_f),
                                        box_x[1]]
                                        #  self.source_antenna['antenna_position'][0] + self.regions['largest_wavelengths'][1] / 4]
            )
            self.source_antenna['box_y'].append(
                [self.source_antenna['antenna_position'][1] - self.regions['largest_wavelengths'][1] / (4*bh_f),
                                         self.source_antenna['antenna_position'][1] + self.regions['largest_wavelengths'][1] / (4*bh_f)]
            )
            self.source_antenna['box_z'].append(
                [0.0, self.source_antenna['antenna_position'][2] - self.regions['largest_wavelengths'][1] / (4*bh_f)]
            )
        else:
            self.source_antenna['box_x'].append(
                [self.source_antenna['antenna_position'][0] - self.regions['largest_wavelengths'][0] / (4*bh_f),
                                            box_x]
                                        #  self.source_antenna['antenna_position'][0] + self.regions['largest_wavelengths'][0] / 4]
            )
            self.source_antenna['box_y'].append(
                [self.source_antenna['antenna_position'][1] - self.regions['largest_wavelengths'][0] / (4*bh_f),
                                         self.source_antenna['antenna_position'][1] + self.regions['largest_wavelengths'][0] / (4*bh_f)]
            )
            self.source_antenna['box_z'].append(
                [0.0, self.source_antenna['antenna_position'][2] - self.regions['largest_wavelengths'][0] / (4*bh_f)]
            )
        print("Source antenna box extents:", self.source_antenna['box_x'], self.source_antenna['box_y'], self.source_antenna['box_z'])

        # Receiver antenna depth, length, and height
        self.receiver_antenna['depth'] = -1 * self.source_antenna['length'] / 4.0
        print(f"Receiver antenna depth: {self.receiver_antenna['depth']} m")
        self.receiver_antenna['length'] = self.source_antenna['length']
        # Assuming the receiver antenna length is for an equilateral triangle, the height of the triangle is given by:
        self.receiver_antenna['height'] = (self.receiver_antenna['length'] * np.sqrt(3)) / 2

        # Inline receiver antenna positions
        PML_Buff = 0.0
        if PML_type == "log":
            PML_Buff = self.PML['layer_thickness']
        else:
            PML_Buff = self.PML['layer_thickness'] * self.PML['num_layers']

        self.receiver_antenna['inline_x'] = [round(0.1 + i * (1.0 - 0.1) / (self.receiver_antenna['num_receivers_inline'] - 1), 5) for i in range(self.receiver_antenna['num_receivers_inline'])]
        self.receiver_antenna['inline_y'] = [0.0] * self.receiver_antenna['num_receivers_inline']
        self.receiver_antenna['inline_z'] = [self.receiver_antenna['depth']] * self.receiver_antenna['num_receivers_inline']
        
        # Inline receiver antenna tetrahedra positions
        self.receiver_antenna['inline_x_tet'] = [x - self.receiver_antenna['length']/2 for x in self.receiver_antenna['inline_x']]
        self.receiver_antenna['inline_x_0_tet'] = self.receiver_antenna['inline_x']
        self.receiver_antenna['inline_x_1_tet'] = [x + self.receiver_antenna['length']/2 for x in self.receiver_antenna['inline_x']]

        self.receiver_antenna['inline_y_tet'] = [-self.receiver_antenna['height']/2] * self.receiver_antenna['num_receivers_inline']
        self.receiver_antenna['inline_y_0_tet'] = [self.receiver_antenna['height']/2] * self.receiver_antenna['num_receivers_inline']
        self.receiver_antenna['inline_y_1_tet'] = [-self.receiver_antenna['height']/2] * self.receiver_antenna['num_receivers_inline']

        self.receiver_antenna['inline_z_tet'] = [0.0] * self.receiver_antenna['num_receivers_inline']

        # Endfire receiver antenna positions
        self.receiver_antenna['endfire_x'] = [0.0] * self.receiver_antenna['num_receivers_endfire']
        self.receiver_antenna['endfire_y'] = [round(0.1 + i * (1.0 - 0.1) / (self.receiver_antenna['num_receivers_endfire'] - 1), 5) for i in range(self.receiver_antenna['num_receivers_endfire'])]
        self.receiver_antenna['endfire_z'] = [self.receiver_antenna['depth']] * self.receiver_antenna['num_receivers_endfire']

        # Endfire receiver antenna tetrahedra positions
        self.receiver_antenna['endfire_x_tet'] = [self.receiver_antenna['length']/2] * self.receiver_antenna['num_receivers_endfire']
        self.receiver_antenna['endfire_x_0_tet'] = [-self.receiver_antenna['length']/2] * self.receiver_antenna['num_receivers_endfire']
        self.receiver_antenna['endfire_x_1_tet'] = [0.0] * self.receiver_antenna['num_receivers_endfire']

        self.receiver_antenna['endfire_y_tet'] = [y - self.receiver_antenna['height'] / 2 for y in self.receiver_antenna['endfire_y']]
        self.receiver_antenna['endfire_y_0_tet'] = [y - self.receiver_antenna['height'] / 2 for y in self.receiver_antenna['endfire_y']]
        self.receiver_antenna['endfire_y_1_tet'] = [y + self.receiver_antenna['height'] / 2 for y in self.receiver_antenna['endfire_y']]

        self.receiver_antenna['endfire_z_tet'] = [0.0] * self.receiver_antenna['num_receivers_endfire']

        # Oblique receiver antenna positions
        self.receiver_antenna['oblique_x'] = [round(x * np.cos(np.pi / 4), 5) for x in self.receiver_antenna['inline_x']] if self.receiver_antenna['num_receivers_oblique'] > 0 else []
        self.receiver_antenna['oblique_y'] = [round(x * np.sin(np.pi / 4), 5) for x in self.receiver_antenna['inline_x']] if self.receiver_antenna['num_receivers_oblique'] > 0 else []
        self.receiver_antenna['oblique_z'] = [self.receiver_antenna['depth']] * self.receiver_antenna['num_receivers_oblique']

        # Oblique receiver antenna tetrahedra positions
        self.receiver_antenna['oblique_x_tet'] = [x - self.receiver_antenna['height'] / 2 for x in self.receiver_antenna['oblique_x']]
        self.receiver_antenna['oblique_x_0_tet'] = [x - self.receiver_antenna['height'] / 2 for x in self.receiver_antenna['oblique_x']]
        self.receiver_antenna['oblique_x_1_tet'] = [x + self.receiver_antenna['height'] / 2 for x in self.receiver_antenna['oblique_x']]

        self.receiver_antenna['oblique_y_tet'] = [y + self.receiver_antenna['length'] / 2 for y in self.receiver_antenna['oblique_y']]
        self.receiver_antenna['oblique_y_0_tet'] = [y - self.receiver_antenna['length'] / 2 for y in self.receiver_antenna['oblique_y']]
        self.receiver_antenna['oblique_y_1_tet'] = self.receiver_antenna['oblique_y']

        self.receiver_antenna['oblique_z_tet'] = [0.0] * self.receiver_antenna['num_receivers_oblique']

        # Combine receiver antenna positions
        self.receiver_antenna['x'] = self.receiver_antenna['inline_x'] + self.receiver_antenna['endfire_x'] + self.receiver_antenna['oblique_x']
        self.receiver_antenna['y'] = self.receiver_antenna['inline_y'] + self.receiver_antenna['endfire_y'] + self.receiver_antenna['oblique_y']
        self.receiver_antenna['z'] = self.receiver_antenna['inline_z'] + self.receiver_antenna['endfire_z'] + self.receiver_antenna['oblique_z']

        # Combine receiver antenna tetrahedra positions
        self.receiver_antenna['x_tet'] = self.receiver_antenna['inline_x_tet'] + self.receiver_antenna['endfire_x_tet'] + self.receiver_antenna['oblique_x_tet']
        self.receiver_antenna['x_0_tet'] = self.receiver_antenna['inline_x_0_tet'] + self.receiver_antenna['endfire_x_0_tet'] + self.receiver_antenna['oblique_x_0_tet']
        self.receiver_antenna['x_1_tet'] = self.receiver_antenna['inline_x_1_tet'] + self.receiver_antenna['endfire_x_1_tet'] + self.receiver_antenna['oblique_x_1_tet']
        
        self.receiver_antenna['y_tet'] = self.receiver_antenna['inline_y_tet'] + self.receiver_antenna['endfire_y_tet'] + self.receiver_antenna['oblique_y_tet']
        self.receiver_antenna['y_0_tet'] = self.receiver_antenna['inline_y_0_tet'] + self.receiver_antenna['endfire_y_0_tet'] + self.receiver_antenna['oblique_y_0_tet']
        self.receiver_antenna['y_1_tet'] = self.receiver_antenna['inline_y_1_tet'] + self.receiver_antenna['endfire_y_1_tet'] + self.receiver_antenna['oblique_y_1_tet']
        
        self.receiver_antenna['z_tet'] = self.receiver_antenna['inline_z_tet'] + self.receiver_antenna['endfire_z_tet'] + self.receiver_antenna['oblique_z_tet']

        # IO file parameters
        n_layers = int(self.model_layers['num_layers'])
        f_low_mhz = int(np.ceil(self.frequency_range['f_low'] / 1e6))
        f_high_mhz = int(np.ceil(self.frequency_range['f_high'] / 1e6))

        self.io = {
            'input_folder': f"in_{experiment_name}",
            'output_folder': f"out_{experiment_name}",
            'input_base_file': "elfe3D_input.txt",
            'input_model_file': f"in/GPR_model_{experiment_name}.",
            'input_model_file_base': f"GPR_model_{experiment_name}.poly",
            'output_E_file': f"out_{experiment_name}/electric_fields",
            'output_H_file': f"out_{experiment_name}/magnetic_fields",
            'vtk_file': f"GPR_model_{experiment_name}.1.vtk",
            'postprocess_base_folder': f"{experiment_name}",
            'source_file': "source.txt",
            'region_params_file': "regionparameters.txt",
            'input_mesh_sizing_file': f"GPR_model_{n_layers}_earth_{f_low_mhz}_{f_high_mhz}_MHz_{self.PML['num_layers']}_PML.mtr",
            'output_electric_fields': "electric_fields_receiver_line.txt",
            'output_magnetic_fields': "magnetic_fields_receiver_line.txt",
        }