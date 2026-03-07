import numpy as np
import math
from itertools import groupby
import os
import bisect

class elfe3DGPRTestWritingPML:
    """
    Class to generate test input files for elfe3D GPR simulations.
    The base values for the tests were generated using the elfe3DGPRTestDesign class.
    """

    def __init__(self, frequency_range, model_layers, regions, source_antenna, receiver_antenna,
                 solver, domain, io, PML):
        """
        Initialize the elfe3DGPRTestWritingPML_inc class with the provided parameters from the elfe3DGPRTestDesign class.
        """

        self.input_base_file = io['input_base_file']
        self.input_model_file = io['input_model_file']
        self.input_model_file_base = io['input_model_file_base']
        self.output_E_file = io['output_E_file']
        self.output_H_file = io['output_H_file']
        self.source_file = io['source_file']
        self.region_params_file = io['region_params_file']
        self.input_mesh_sizing_file = io['input_mesh_sizing_file']
        self.base_folder = io['output_folder']
        
        self.f_list = source_antenna['f_list']

        self.frequency = frequency_range['f_low'] # Frequency in Hz
        self.NUM_PML = PML['num_layers']
        self.pml_thickness = PML['layer_thickness']
        self.PML_theory = PML['PML_theory']
        self.PML_decay_type = PML['PML_decay_type']
        self.num_region = 0

        # Maximum Domain Size
        self.x_min = domain['x_min']
        self.x_max = domain['x_max']
        self.y_min = domain['y_min']
        self.y_max = domain['y_max']
        self.z_min = domain['z_min']
        self.z_max = domain['z_max']

        # Earth-Layer Interfaces
        self.num_layers = model_layers['num_layers']
        self.layer_interfaces = {
            'x_min': model_layers['x_min'],
            'x_max': model_layers['x_max'],
            'y_min': model_layers['y_min'],
            'y_max': model_layers['y_max'],
            'z': model_layers['z'],
        }

        # Layer Electromagnetic Parameters
        self.region_epsilon_r = regions['eps_r']
        self.region_sigma = regions['sigma']
        self.region_mu_r = regions['mu_r']
        self.region_sigma_m = regions['sigma_m']
        self.region_rho = regions['rho']
        self.region_volumes = regions['max_element_volume']
        self.region_element_edges = regions['max_element_edges']

        # Source Parameters
        self.source_type = source_antenna['source_type']
        self.current_direction = source_antenna['current_direction']
        self.source_moment = source_antenna['source_moment']
        self.m = source_antenna['m']
        self.box_present = source_antenna['box_present']

        # Source Number and Locations
        self.num_source_segments = source_antenna['num_segments']
        self.length_antenna = source_antenna['length']
        self.source_x = source_antenna['x_extents']
        self.source_y = source_antenna['y_extents']
        self.source_z = source_antenna['z_extents']

        self.source_x_disc = source_antenna['x_disc']
        self.source_y_disc = source_antenna['y_disc']
        self.source_z_disc = source_antenna['z_disc']


        if self.box_present:
            self.source_x_box = source_antenna['box_x']
            self.source_y_box = source_antenna['box_y']
            self.source_z_box = source_antenna['box_z']

        # Receiver Locations
        self.receivers_x = receiver_antenna['x']
        self.receivers_y = receiver_antenna['y']
        self.receivers_z = receiver_antenna['z']

        self.receivers_x_tet = receiver_antenna['x_tet']
        self.receivers_x_0_tet = receiver_antenna['x_0_tet']
        self.receivers_x_1_tet = receiver_antenna['x_1_tet']

        self.receivers_y_tet = receiver_antenna['y_tet']
        self.receivers_y_0_tet = receiver_antenna['y_0_tet']
        self.receivers_y_1_tet = receiver_antenna['y_1_tet']

        self.receivers_z_tet = receiver_antenna['z_tet']
        
        # Solver Parameters
        self.solver = solver['solver']
        self.maxRefSteps = solver['maxRefSteps']
        self.maxUnknowns = int(solver['maxUnknowns'])
        self.betaRef = solver['betaRef']
        self.accuracyTol = solver['accuracyTol']
        self.vtk = solver['vtk']
        self.errorEst_method = solver['errorEst_method']
        self.refStrategy = solver['refStrategy']
        self.output_fields_vtk = solver['output_fields_vtk']

    ####################################
    # Writing the elfe3D_input.txt file
    ####################################
    def write_base_input_file(self):
        """
        Function to write the elfe3D_input file with the specified parameters.
        """

        self.input_base_file = os.path.join(self.base_folder, self.input_base_file)
        with open(self.input_base_file, 'w+') as f:
            f.write("solver                  {}\n".format(self.solver))
            f.write("model_size\n")

            PML_BUFF = self.pml_thickness * self.NUM_PML
            f.write("{}\t{}\t{}\n".format(round(self.x_min-PML_BUFF, 4), round(self.y_min-PML_BUFF, 4), round(self.z_min-PML_BUFF, 4)))
            f.write("{}\t{}\t{}\n".format(round(self.x_max+PML_BUFF, 4), round(self.y_max+PML_BUFF, 4), round(self.z_max+PML_BUFF, 4)))

            f.write("num_freq                {}\n".format(len(self.f_list)))
            for i in range(len(self.f_list)):
                f.write("{}\n".format(self.f_list[i]))
            f.write("num_rec                 {}\n".format(len(self.receivers_x)))
            for i in range(len(self.receivers_x)):
                f.write("{} {} {}\n".format(round(self.receivers_x[i], 5), round(self.receivers_y[i], 5), self.receivers_z[i]))
            f.write("output_E_file           {}\n".format(self.output_E_file))
            f.write("output_H_file           {}\n".format(self.output_H_file))
            f.write("source_type             {}\n".format(self.source_type))
            f.write("{} {} {}\n".format(round(self.source_x[0], 5), round(self.source_y[0], 5), round(self.source_z[0], 5)))
            f.write("{} {} {}\n".format(round(self.source_x[1], 5), round(self.source_y[1], 5), round(self.source_z[1], 5)))
            f.write("current_direction       {}\n".format(self.current_direction))
            f.write("source_moment           {}\n".format(self.source_moment))
            f.write("PEC_present             {}\n".format(0))
            f.write("num_PEC                 {}\n".format(0))
            f.write("model_file_name         {}\n".format(self.input_model_file))
            f.write("maxRefSteps             {}\n".format(self.maxRefSteps))
            f.write("maxUnknowns             {}\n".format(self.maxUnknowns))
            f.write("betaRef                 {}\n".format(self.betaRef))
            f.write("accuracyTol             {:.6f}\n".format(self.accuracyTol))
            f.write("vtkRef                  {}\n".format(self.vtk))
            f.write("errorEst_method         {}\n".format(self.errorEst_method))
            f.write("refStrategy             {}\n".format(self.refStrategy))
            f.write("output_fields_vtk       {}\n".format(self.output_fields_vtk))
            f.write("PML_present             {}\n".format(int(self.NUM_PML > 0)))
            f.write("PML_thickness           {}\n".format(self.pml_thickness))
            # f.write("PML_theory              {}\n".format(self.PML_theory))
            f.write("PML_decay_type          {}\n".format(self.PML_decay_type))

    ############################################
    # Evaluating the .poly file aspect - Nodes
    ############################################
    def create_rectangular_prism_nodes(self, x_min, x_max, y_min, y_max, z, marker, num_node):
        """
        Helper function for the base domain nodes generation.
        Function to create nodes for a rectangular prism.
        """

        nodes = []
        for i in [0, 1, 3, 2]:  # Counterclockwise order
            x = x_max if i & 1 else x_min
            y = y_max if i & 2 else y_min
            num_node += 1
            nodes.append([num_node, x, y, z, marker])
        return nodes, num_node

    def evaluate_all_nodes(self):
        """
        Function that defines all the nodes that will be written in the .poly file.
        Calls also the function that evaluates nodes for the PML layers.
        """

        self.num_node = 0

        # Corner Nodes for the whole domain - rectangular prism
        node_list_domain_prism, self.num_node = self.create_rectangular_prism_nodes(self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, 1, self.num_node)
        top_nodes, self.num_node = self.create_rectangular_prism_nodes(self.x_min, self.x_max, self.y_min, self.y_max, self.z_max, 1, self.num_node)
        node_list_domain_prism.extend(top_nodes)

        # Nodes for the interfaces
        node_list_of_list_of_interfaces = []
        self.marker_interfaces = [0, 1, 7, 8, 9, 10]

        for layer in range(self.num_layers):
            interface_nodes, self.num_node = self.create_rectangular_prism_nodes(
                self.layer_interfaces['x_min'], self.layer_interfaces['x_max'],
                self.layer_interfaces['y_min'], self.layer_interfaces['y_max'],
                self.layer_interfaces['z'][layer], self.marker_interfaces[layer], self.num_node
            )
            node_list_of_list_of_interfaces.append(interface_nodes)

        # Nodes for the receiver
        node_list_receiver = []
        marker_receivers = 99

        for i in range(len(self.receivers_x)):
            self.num_node += 1
            node_list_receiver.append([self.num_node, self.receivers_x_0_tet[i], self.receivers_y_0_tet[i], self.receivers_z_tet[i], marker_receivers])
            self.num_node += 1
            node_list_receiver.append([self.num_node, self.receivers_x_1_tet[i], self.receivers_y_1_tet[i], self.receivers_z_tet[i], marker_receivers])
            self.num_node += 1
            node_list_receiver.append([self.num_node, self.receivers_x_tet[i], self.receivers_y_tet[i], self.receivers_z_tet[i], marker_receivers])

        # Nodes for the source - discretized
        node_list_source = []
        marker_source = 0
        for i in range(self.num_source_segments+1):
            self.num_node += 1
            node_list_source.append([self.num_node, self.source_x_disc[i], self.source_y_disc[i], self.source_z_disc[i], marker_source])

        # Nodes for the source - box
        if self.box_present:
            self.marker_source_box = 98
            if self.source_z_box[0][0] == -1*self.source_z_box[1][1]: # hence, only one medium
                box_z = [self.source_z_box[0][0], 0.0, self.source_z_box[1][1]]
                box_x = [self.source_x_box[0][0], self.source_x_box[0][1]]
                box_y = [self.source_y_box[0][0], self.source_y_box[0][1]]
                node_list_source_box = []
                for i in range(3):
                    for j in range(2):
                        for k in range(2):
                            self.num_node += 1
                            node_list_source_box.append([self.num_node, box_x[k], box_y[j], box_z[i], self.marker_source_box])

            else:  # hence, two or more mediums
                node_list_source_box = []

                for m in range(2):
                    box_z = [self.source_z_box[m][0], self.source_z_box[m][1]]
                    box_x = [self.source_x_box[m][0], self.source_x_box[m][1]]
                    box_y = [self.source_y_box[m][0], self.source_y_box[m][1]]
                    for i in range(2):
                        for j in range(2):
                            for k in range(2):
                                self.num_node += 1
                                node_list_source_box.append([self.num_node,  box_x[k],  box_y[j], box_z[i], self.marker_source_box])

            # out of all source box nodes, find ones that have z==0 and save in new var
            self.node_list_source_box_0 = [node for node in node_list_source_box if node[3] == 0]


        self.node_list_domain_prism = node_list_domain_prism
        self.node_list_of_list_of_interfaces = node_list_of_list_of_interfaces
        self.node_list_receiver = node_list_receiver
        self.node_list_source = node_list_source
        if self.box_present:
            self.node_list_source_box = node_list_source_box

        # Nodes for the PML Layers
        if self.NUM_PML > 0:
            self.node_list_PML_2, self.node_list_PML_3, self.pml_type_2_layer_nodes  = self.evaluate_all_nodes_PML()

    def evaluate_all_nodes_PML(self):
        """
        Function to evaluate nodes for the PML layers.

        Each layer has three types of extensions:
        1. Face extensions: Simplest type, each layer has just one new prism per face shifted along its axis.
        2. Edge extensions: First layer has one prism per edge, next one has 3, and so on: (2 * layer_number - 1) prisms.
        3. Vertex extensions: First layer has one prism per vertex, next one has 7, and so on: (3 * layer_number**2 - 3 * layer_number + 1) prisms.
        """

        t = self.pml_thickness

        marker_pml_layers = [1000 for _ in range(self.NUM_PML)]
        self.marker_pml_layers = marker_pml_layers
        self.marker_pml_types = [1, 2, 3]

        pml_type_2_nodes = []
        pml_type_2_layer_nodes = []
        pml_type_3_nodes = []
        
        z_layers = [self.layer_interfaces['z'][j] for j in range(self.num_layers)]

        # Type 2:

        # Front-Right edge addition
        x = np.linspace(self.x_min, self.x_min - self.NUM_PML * t, self.NUM_PML+1)
        y = np.linspace(self.y_min, self.y_min - self.NUM_PML * t, self.NUM_PML+1)
        xyz = np.array(np.meshgrid(x, y, z_layers)).T.reshape(-1, 3)
        xyz = xyz[~((xyz[:, 0] == self.x_min) & (xyz[:, 1] == self.y_min))]

        for x, y, z in xyz:
            self.num_node += 1
            pml_type_2_layer_nodes.append([self.num_node, x, y, z, marker_pml_layers[0] + self.marker_pml_types[1]])
            pml_type_2_nodes.append(pml_type_2_layer_nodes[-1])

    
        # Front-Left edge addition
        x = np.linspace(self.x_max, self.x_max + self.NUM_PML * t, self.NUM_PML+1)
        y = np.linspace(self.y_min, self.y_min - self.NUM_PML * t, self.NUM_PML+1)
        xyz = np.array(np.meshgrid(x, y, z_layers)).T.reshape(-1, 3)
        xyz = xyz[~((xyz[:, 0] == self.x_max) & (xyz[:, 1] == self.y_min))]

        for x, y, z in xyz:
            self.num_node += 1
            pml_type_2_layer_nodes.append([self.num_node, x, y, z, marker_pml_layers[0] + self.marker_pml_types[1]])
            pml_type_2_nodes.append(pml_type_2_layer_nodes[-1])

        # Back-Left edge addition
        x = np.linspace(self.x_max, self.x_max + self.NUM_PML * t, self.NUM_PML+1)
        y = np.linspace(self.y_max, self.y_max + self.NUM_PML * t, self.NUM_PML+1)
        xyz = np.array(np.meshgrid(x, y, z_layers)).T.reshape(-1, 3)
        xyz = xyz[~((xyz[:, 0] == self.x_max) & (xyz[:, 1] == self.y_max))]

        for x, y, z in xyz:
            self.num_node += 1
            pml_type_2_layer_nodes.append([self.num_node, x, y, z, marker_pml_layers[0] + self.marker_pml_types[1]])
            pml_type_2_nodes.append(pml_type_2_layer_nodes[-1])

        # Back-Right edge addition
        x = np.linspace(self.x_min, self.x_min - self.NUM_PML * t, self.NUM_PML+1)
        y = np.linspace(self.y_max, self.y_max + self.NUM_PML * t, self.NUM_PML+1)
        xyz = np.array(np.meshgrid(x, y, z_layers)).T.reshape(-1, 3)
        xyz = xyz[~((xyz[:, 0] == self.x_min) & (xyz[:, 1] == self.y_max))]

        for x, y, z in xyz:
            self.num_node += 1
            pml_type_2_layer_nodes.append([self.num_node, x, y, z, marker_pml_layers[0] + self.marker_pml_types[1]])
            pml_type_2_nodes.append(pml_type_2_layer_nodes[-1])
    
        # Type 3

        # Front-Right-Top vertex addition
        x = np.linspace(self.x_max, self.x_max + self.NUM_PML * t, self.NUM_PML+1)
        y = np.linspace(self.y_min, self.y_min - self.NUM_PML * t, self.NUM_PML+1)
        z = np.linspace(self.z_max, self.z_max + self.NUM_PML * t, self.NUM_PML+1)
        xyz = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)[1:]

        for x, y, z in xyz:
            self.num_node += 1
            pml_type_3_nodes.append([self.num_node, x, y, z, marker_pml_layers[0] + self.marker_pml_types[2]])

        # Front-Right-Bottom vertex addition
        x = np.linspace(self.x_max, self.x_max + self.NUM_PML * t, self.NUM_PML+1)
        y = np.linspace(self.y_min, self.y_min - self.NUM_PML * t, self.NUM_PML+1)
        z = np.linspace(self.z_min, self.z_min - self.NUM_PML * t, self.NUM_PML+1)
        xyz = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)[1:]

        for x, y, z in xyz:
            self.num_node += 1
            pml_type_3_nodes.append([self.num_node, x, y, z, marker_pml_layers[0] + self.marker_pml_types[2]])

        # Front-Left-Top vertex addition
        x = np.linspace(self.x_min, self.x_min - self.NUM_PML * t, self.NUM_PML+1)
        y = np.linspace(self.y_min, self.y_min - self.NUM_PML * t, self.NUM_PML+1)
        z = np.linspace(self.z_max, self.z_max + self.NUM_PML * t, self.NUM_PML+1)
        xyz = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)[1:]

        for x, y, z in xyz:
            self.num_node += 1
            pml_type_3_nodes.append([self.num_node, x, y, z, marker_pml_layers[0] + self.marker_pml_types[2]])

        # Front-Left-Bottom vertex addition
        x = np.linspace(self.x_min, self.x_min - self.NUM_PML * t, self.NUM_PML+1)
        y = np.linspace(self.y_min, self.y_min - self.NUM_PML * t, self.NUM_PML+1)
        z = np.linspace(self.z_min, self.z_min - self.NUM_PML * t, self.NUM_PML+1)
        xyz = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)[1:]

        for x, y, z in xyz:
            self.num_node += 1
            pml_type_3_nodes.append([self.num_node, x, y, z, marker_pml_layers[0] + self.marker_pml_types[2]])

        # Back-Right-Top vertex addition
        x = np.linspace(self.x_max, self.x_max + self.NUM_PML * t, self.NUM_PML+1)
        y = np.linspace(self.y_max, self.y_max + self.NUM_PML * t, self.NUM_PML+1)
        z = np.linspace(self.z_max, self.z_max + self.NUM_PML * t, self.NUM_PML+1)
        xyz = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)[1:]

        for x, y, z in xyz:
            self.num_node += 1
            pml_type_3_nodes.append([self.num_node, x, y, z, marker_pml_layers[0] + self.marker_pml_types[2]])

        # Back-Right-Bottom vertex addition
        x = np.linspace(self.x_max, self.x_max + self.NUM_PML * t, self.NUM_PML+1)
        y = np.linspace(self.y_max, self.y_max + self.NUM_PML * t, self.NUM_PML+1)
        z = np.linspace(self.z_min, self.z_min - self.NUM_PML * t, self.NUM_PML+1)
        xyz = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)[1:]

        for x, y, z in xyz:
            self.num_node += 1
            pml_type_3_nodes.append([self.num_node, x, y, z, marker_pml_layers[0] + self.marker_pml_types[2]])

        # Back-Left-Top vertex addition
        x = np.linspace(self.x_min, self.x_min - self.NUM_PML * t, self.NUM_PML+1)
        y = np.linspace(self.y_max, self.y_max + self.NUM_PML * t, self.NUM_PML+1)
        z = np.linspace(self.z_max, self.z_max + self.NUM_PML * t, self.NUM_PML+1)
        xyz = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)[1:]

        for x, y, z in xyz:
            self.num_node += 1
            pml_type_3_nodes.append([self.num_node, x, y, z, marker_pml_layers[0] + self.marker_pml_types[2]])

        # Back-Left-Bottom vertex addition
        x = np.linspace(self.x_min, self.x_min - self.NUM_PML * t, self.NUM_PML+1)
        y = np.linspace(self.y_max, self.y_max + self.NUM_PML * t, self.NUM_PML+1)
        z = np.linspace(self.z_min, self.z_min - self.NUM_PML * t, self.NUM_PML+1)
        xyz = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)[1:]

        for x, y, z in xyz:
            self.num_node += 1
            pml_type_3_nodes.append([self.num_node, x, y, z, marker_pml_layers[0] + self.marker_pml_types[2]])

        return pml_type_2_nodes, pml_type_3_nodes, pml_type_2_layer_nodes

    ############################################
    # Evaluating the .poly file aspect - Facets
    ############################################  
    def angular_sort(self, nodes, center_axis_1, center_axis_2):
        """
        Helper function for facet generation of domain and interfaces, and indirectly, also for the PML.
        Sort nodes angularly based on their relative positions in the x-z plane, considering their angular positions.
        """
        # Calculate the centroid
        center_x = sum(node[center_axis_1] for node in nodes) / len(nodes)
        center_z = sum(node[center_axis_2] for node in nodes) / len(nodes)

        # Function to calculate the angle relative to the centroid
        def angle_from_center(node):
            x = node[center_axis_1]
            z = node[center_axis_2]
            return -math.atan2(z - center_z, x - center_x)  # Negate the angle to reverse order

        # Sort the nodes based on their angle from the center
        sorted_nodes = sorted(nodes, key=angle_from_center)

        # Select the node with the smallest numerical value (e.g., smallest node ID)
        start_node = min(sorted_nodes, key=lambda node: node[0])  # assuming the node ID is the first element

        # Find the index of the start node in the sorted list
        start_index = sorted_nodes.index(start_node)

        # Reorder the list so it starts from the selected node
        sorted_nodes = sorted_nodes[start_index:] + sorted_nodes[:start_index]

        # Return the sorted node IDs in the correct order
        return [node[0] for node in sorted_nodes]

    def get_face_nodes(self, node_list, axis_value, axis_index, sort_axis_1, sort_axis_2):
        """
        Helper function for facet generation of domain and interfaces, and indirectly, also for the PML.
        Filters and sorts nodes based on their coordinates and axes.
        """
        # Filter nodes based on the y-axis value
        filtered_nodes = [node for node in node_list if node[axis_index] == axis_value]
        # Sort the nodes based on their angular position in the x-z plane
        sorted_nodes = self.angular_sort(filtered_nodes, center_axis_1=sort_axis_1, center_axis_2=sort_axis_2)
        
        return sorted_nodes

    def create_face_string(self, marker, nodes, interfaces_nodes, axis_value, axis_index, sort_axis_1, sort_axis_2):
        """
        Helper function for facet generation of the PML layers, Type 1 and Type 2.
        Create a string representation of the facet information.
        """

        face_nodes_list = self.get_face_nodes(
            node_list=nodes, 
            axis_value=axis_value, 
            axis_index=axis_index,  
            sort_axis_1=sort_axis_1,  
            sort_axis_2=sort_axis_2
        )
        num_face_nodes = len(face_nodes_list)

        if axis_index == 3:
            interfaces_nodes_list = []
        elif interfaces_nodes is None:
            interfaces_nodes_list = []
        else:
            interfaces_nodes = [list(g) for _, g in groupby(sorted(interfaces_nodes, key=lambda n: n[3]), key=lambda n: n[3])]  # Group nodes by z-axis value
            interfaces_nodes_list = [self.get_face_nodes(
            node_list=interface_nodes, 
            axis_value=axis_value, 
            axis_index=axis_index,  
            sort_axis_1=sort_axis_1,  
            sort_axis_2=sort_axis_2
        ) for interface_nodes in interfaces_nodes]
            
        if axis_index == 3 or interfaces_nodes is None:
            num_polygon = 1
        else:
            num_polygon = 1 + self.num_layers
        face_information = [num_polygon, 0, marker]

        str_1 = " ".join(map(str, face_information))
        str_2 = " ".join(map(str, face_nodes_list))
        str_3 = "\n".join(f"{len(interface_nodes)} " + " ".join(map(str, interface_nodes)) for interface_nodes in interfaces_nodes_list)

        if interfaces_nodes is None:
            return f"{str_1}\n{num_face_nodes} {str_2}\n"
        else:
            return f"\n{str_1}\n{num_face_nodes} {str_2}\n{str_3}\n"

    def create_cuboid_faces_from_nodes(self, nodes, tol=1e-6):
        """
        Helper function for cuboid generation.
        Given a list of 8 nodes representing the vertices of a cuboid in arbitrary order,
        returns the 6 faces (each as a list of 4 nodes) with the vertices arranged in proper order.
        """
        
        # Group nodes by their y coordinate (using tolerance)
        groups = {}
        for node in nodes:
            y = node[2]
            found_key = None
            for key in groups:
                if abs(key - y) < tol:
                    found_key = key
                    break
            if found_key is None:
                groups[y] = []
                found_key = y
            groups[found_key].append(node)
        
        if len(groups) != 2:
            raise ValueError("Expected exactly 2 distinct y values for a cuboid (top and bottom faces).")
        
        # Identify the bottom and top groups (by y value)
        ys = sorted(groups.keys())
        bottom_nodes = groups[ys[0]]
        top_nodes = groups[ys[1]]
        
        if len(bottom_nodes) != 4 or len(top_nodes) != 4:
            raise ValueError("Expected exactly 4 nodes per face (bottom and top).")
        
        # Order nodes in a face using the x and z coordinates (assumes faces are horizontal)
        def order_face(face_nodes):
            cx = sum(n[1] for n in face_nodes) / len(face_nodes)
            cz = sum(n[3] for n in face_nodes) / len(face_nodes)
            # Compute angle relative to the centroid (counterclockwise)
            return sorted(face_nodes, key=lambda n: math.atan2(n[3] - cz, n[1] - cx))
        
        bottom_ordered = order_face(bottom_nodes)
        top_ordered = order_face(top_nodes)
        
        # For consistency, we assume that the order produced for both faces is such that the i-th
        # node in the bottom corresponds to the i-th node in the top. (This is true if the cuboid
        # is axis-aligned in y and the top and bottom faces are parallel.)
        
        # Build lateral faces:
        lateral_faces = []
        for i in range(4):
            next_i = (i + 1) % 4
            # Each lateral face connects an edge from the bottom to the corresponding edge from the top.
            face = [bottom_ordered[i], bottom_ordered[next_i],
                    top_ordered[next_i], top_ordered[i]]
            lateral_faces.append(face)
        
        # Define top and bottom faces.
        # (They are already ordered counterclockwise.)
        top_face = top_ordered[:]      # higher y value
        bottom_face = bottom_ordered[:]  # lower y value        

        
        # Return all faces: lateral faces first, then top and bottom (or any preferred order)
        faces = lateral_faces + [top_face, bottom_face]
        return faces

    def create_half_cuboid_faces_from_nodes(self, nodes, not_needed_node, tol=1e-6):
        """
        Helper function for cuboid generation.
        Given a list of 8 nodes representing the vertices of a cuboid in arbitrary order,
        returns the 6 faces (each as a list of 4 nodes) with the vertices arranged in proper order.
        """
        
        # Group nodes by their y coordinate (using tolerance)
        groups = {}
        for node in nodes:
            y = node[2]
            found_key = None
            for key in groups:
                if abs(key - y) < tol:
                    found_key = key
                    break
            if found_key is None:
                groups[y] = []
                found_key = y
            groups[found_key].append(node)
        
        if len(groups) != 2:
            raise ValueError("Expected exactly 2 distinct y values for a cuboid (top and bottom faces).")
        
        # Identify the bottom and top groups (by y value)
        ys = sorted(groups.keys())
        bottom_nodes = groups[ys[0]]
        top_nodes = groups[ys[1]]
        
        if len(bottom_nodes) != 4 or len(top_nodes) != 4:
            raise ValueError("Expected exactly 4 nodes per face (bottom and top).")
        
        # Order nodes in a face using the x and z coordinates (assumes faces are horizontal)
        def order_face(face_nodes):
            cx = sum(n[1] for n in face_nodes) / len(face_nodes)
            cz = sum(n[3] for n in face_nodes) / len(face_nodes)
            # Compute angle relative to the centroid (counterclockwise)
            return sorted(face_nodes, key=lambda n: math.atan2(n[3] - cz, n[1] - cx))
        
        bottom_ordered = order_face(bottom_nodes)
        top_ordered = order_face(top_nodes)
        
        # For consistency, we assume that the order produced for both faces is such that the i-th
        # node in the bottom corresponds to the i-th node in the top. (This is true if the cuboid
        # is axis-aligned in y and the top and bottom faces are parallel.)
        
        # Build lateral faces:
        lateral_faces = []
        for i in range(4):
            next_i = (i + 1) % 4
            # Each lateral face connects an edge from the bottom to the corresponding edge from the top.
            face = [bottom_ordered[i], bottom_ordered[next_i],
                    top_ordered[next_i], top_ordered[i]]
            lateral_faces.append(face)
        
        # Define top and bottom faces.
        # (They are already ordered counterclockwise.)
        top_face = top_ordered[:]      # higher y value
        bottom_face = bottom_ordered[:]  # lower y value        

        
        # Return all faces: lateral faces first, then top and bottom (or any preferred order)
        faces = lateral_faces + [top_face, bottom_face]

        # Remove the face that contains the not_needed_node by checking the node's coordinates
        needed_faces = []
        for face in faces:
            # Only append the face if none of its nodes match the not_needed_node (by coordinates with tolerance)
            if not any(
                (float(node[1]) == float(not_needed_node[0])) and
                (float(node[2]) == float(not_needed_node[1])) and
                (float(node[3]) == float(not_needed_node[2]))
                for node in face
            ):
                needed_faces.append(face)
        return needed_faces


    def evaluate_all_facets(self):
        """
        Function that defines all the facets that will be written in the .poly file.
        Calls the function that evaluates facets for the PML layers as well.
        """

        # Facets for the domain with the air-earth interface
        self.num_facet = 0

        # Front Face
        self.num_facet += 1
        num_polygon = 0
        num_face_nodes = 0

        front_face_nodes_list = []
        front_face_interfaces_nodes_list = []
        
        node_list_interfaces = []
        for node_list in self.node_list_of_list_of_interfaces:
                node_list_interfaces += node_list

        # FRONT FACE: Domain and interface
        front_face_nodes_list = self.get_face_nodes(
            node_list=self.node_list_domain_prism + node_list_interfaces, 
            axis_value=self.y_min, 
            axis_index=2,  # y-axis
            sort_axis_1=1,  # x-axis
            sort_axis_2=3   # z-axis
        )
        num_face_nodes += len(front_face_nodes_list)

        front_face_interfaces_nodes_list = [[] for _ in range(len(self.node_list_of_list_of_interfaces))]
        for i in range(len(self.node_list_of_list_of_interfaces)):
            front_face_interfaces_nodes_list[i] = self.get_face_nodes(
                node_list=self.node_list_of_list_of_interfaces[i], 
                axis_value=self.y_min, 
                axis_index=2,  # y-axis
                sort_axis_1=1,  # x-axis
                sort_axis_2=3   # z-axis
            )

        # Construct facet information
        num_polygon += 1+len(self.node_list_of_list_of_interfaces)
        front_face_information = [num_polygon, 0, 1]

        # Convert lists to strings for output
        str_1 = " ".join(map(str, front_face_information))
        str_2 = " ".join(map(str, front_face_nodes_list))
        str_3 = ""
        for i, interface_nodes in enumerate(front_face_interfaces_nodes_list):
            str_3 += f"{len(front_face_interfaces_nodes_list[i])} "
            str_3 += " ".join(map(str, interface_nodes)) + "\n"

        # Create facet string
        domain_front_face_string = f"#front\n{str_1}\n{num_face_nodes} {str_2}\n{str_3}\n"

        # Right Face
        self.num_facet += 1
        num_polygon = 0
        num_face_nodes = 0

        right_face_nodes_list = []
        right_face_interfaces_nodes_list = []

        right_face_nodes_list = self.get_face_nodes(
            node_list=self.node_list_domain_prism + node_list_interfaces, 
            axis_value=self.x_max, 
            axis_index=1,  # x-axis
            sort_axis_1=2,  # y-axis
            sort_axis_2=3   # z-axis
        )
        num_face_nodes += len(right_face_nodes_list)

        right_face_interfaces_nodes_list = []
        for node_list in self.node_list_of_list_of_interfaces:
            right_face_interfaces_nodes_list.append(self.get_face_nodes(
                node_list=node_list, 
                axis_value=self.x_max, 
                axis_index=1,  # x-axis
                sort_axis_1=2,  # y-axis
                sort_axis_2=3   # z-axis
            ))

        num_polygon += 1 + len(self.node_list_of_list_of_interfaces)
        right_face_information = [num_polygon, 0, 1]

        str_1 = " ".join(map(str, right_face_information))
        str_2 = " ".join(map(str, right_face_nodes_list))
        str_3 = ""
        for i, interface_nodes in enumerate(right_face_interfaces_nodes_list):
            str_3 += f"{len(right_face_interfaces_nodes_list[i])} "
            str_3 += " ".join(map(str, interface_nodes)) + "\n"

        domain_right_face_string = f"#right\n{str_1}\n{num_face_nodes} {str_2}\n{str_3}\n"

        # Back Face
        self.num_facet += 1
        num_polygon = 0
        num_face_nodes = 0

        back_face_nodes_list = []
        back_face_interfaces_nodes_list = []

        back_face_nodes_list = self.get_face_nodes(
            node_list=self.node_list_domain_prism + node_list_interfaces, 
            axis_value=self.y_max, 
            axis_index=2,  # y-axis
            sort_axis_1=1,  # x-axis
            sort_axis_2=3   # z-axis
        )
        num_face_nodes += len(back_face_nodes_list)

        back_face_interfaces_nodes_list = [[] for _ in range(len(self.node_list_of_list_of_interfaces))]
        for i in range(len(self.node_list_of_list_of_interfaces)):
            back_face_interfaces_nodes_list[i] = self.get_face_nodes(
                node_list=self.node_list_of_list_of_interfaces[i], 
                axis_value=self.y_max, 
                axis_index=2,  # y-axis
                sort_axis_1=1,  # x-axis
                sort_axis_2=3   # z-axis
            )

        num_polygon += 1 + len(self.node_list_of_list_of_interfaces)
        back_face_information = [num_polygon, 0, 1]

        str_1 = " ".join(map(str, back_face_information))
        str_2 = " ".join(map(str, back_face_nodes_list))
        str_3 = ""
        for i, interface_nodes in enumerate(back_face_interfaces_nodes_list):
            str_3 += f"{len(back_face_interfaces_nodes_list[i])} "
            str_3 += " ".join(map(str, interface_nodes)) + "\n"

        domain_back_face_string = f"#back\n{str_1}\n{num_face_nodes} {str_2}\n{str_3}\n"

        # Left Face
        self.num_facet += 1
        num_polygon = 0
        num_face_nodes = 0

        left_face_nodes_list = []
        left_face_interfaces_nodes_list = []

        left_face_nodes_list = self.get_face_nodes(
            node_list=self.node_list_domain_prism + node_list_interfaces, 
            axis_value=self.x_min, 
            axis_index=1,  # x-axis
            sort_axis_1=2,  # y-axis
            sort_axis_2=3   # z-axis
        )
        num_face_nodes += len(left_face_nodes_list)

        left_face_interfaces_nodes_list = [[] for _ in range(len(self.node_list_of_list_of_interfaces))]
        for i in range(len(self.node_list_of_list_of_interfaces)):
            left_face_interfaces_nodes_list[i] = self.get_face_nodes(
                node_list=self.node_list_of_list_of_interfaces[i], 
                axis_value=self.x_min, 
                axis_index=1,  # x-axis
                sort_axis_1=2,  # y-axis
                sort_axis_2=3   # z-axis
            )

        num_polygon += 1 + len(self.node_list_of_list_of_interfaces)
        left_face_information = [num_polygon, 0, 1]

        str_1 = " ".join(map(str, left_face_information))
        str_2 = " ".join(map(str, left_face_nodes_list))
        str_3 = ""
        for i, interface_nodes in enumerate(left_face_interfaces_nodes_list):
            str_3 += f"{len(left_face_interfaces_nodes_list[i])} "
            str_3 += " ".join(map(str, interface_nodes)) + "\n"

        domain_left_face_string = f"#left\n{str_1}\n{num_face_nodes} {str_2}\n{str_3}\n"

        # top
        self.num_facet += 1
        num_polygon = 1
        top_face_nodes_list = []
        
        top_face_nodes_list = self.get_face_nodes(
            node_list=self.node_list_domain_prism, 
            axis_value=self.z_min, 
            axis_index=3,  # z-axis
            sort_axis_1=1,  # x-axis
            sort_axis_2=2   # y-axis
        )

        top_face_information = [num_polygon, 0, 1]

        str_1 = " ".join(map(str, top_face_information))
        str_2 = " ".join(map(str, top_face_nodes_list))
        domain_top_face_string = f"# top and bottom\n{str_1}\n4 {str_2}\n"

        # bottom
        self.num_facet += 1
        num_polygon = 1
        bottom_face_nodes_list = []
        bottom_face_nodes_list = self.get_face_nodes(
            node_list=self.node_list_domain_prism, 
            axis_value=self.z_max, 
            axis_index=3,  # z-axis
            sort_axis_1=1,  # x-axis
            sort_axis_2=2   # y-axis
        )

        bottom_face_information = [num_polygon, 0, 1]

        str_1 = " ".join(map(str, bottom_face_information))
        str_2 = " ".join(map(str, bottom_face_nodes_list))
        domain_bottom_face_string = f"{str_1}\n4 {str_2}\n"

        self.domain_string = "# polygons - # holes - boundary marker \n" + domain_front_face_string + domain_right_face_string + domain_back_face_string + domain_left_face_string + domain_top_face_string + domain_bottom_face_string

        # source box
        if self.box_present:
            if len(self.node_list_source_box) == 12:
                source_box_1_faces = self.create_cuboid_faces_from_nodes(self.node_list_source_box[:8])
                source_box_2_faces = self.create_cuboid_faces_from_nodes(self.node_list_source_box[4:])
            elif len(self.node_list_source_box) == 16:
                source_box_1_faces = self.create_cuboid_faces_from_nodes(self.node_list_source_box[:8])
                source_box_2_faces = self.create_cuboid_faces_from_nodes(self.node_list_source_box[8:])
            else:
                source_box_1_faces = []
                source_box_2_faces = []

            # Remove faces from source_box_1_faces and source_box_2_faces where all 4 nodes have z == 0
            def filter_faces_no_z0(faces):
                faces_with_z0 = []
                faces_without_z0 = []
                for face in faces:
                    if all(abs(node[3]) < 1e-8 for node in face):
                        faces_with_z0.append(face)
                    else:
                        faces_without_z0.append(face)
                return faces_without_z0, faces_with_z0

            source_box_1_faces, source_box_1_faces_z0 = filter_faces_no_z0(source_box_1_faces)
            source_box_2_faces, source_box_2_faces_z0 = filter_faces_no_z0(source_box_2_faces)

            # Build a string for all source box faces at z=0 (from both boxes)
            source_box_z0_string = ""
            for face in source_box_1_faces_z0 + source_box_2_faces_z0:
                node_indices = [node[0] for node in face]
                source_box_z0_string += "4 " + " ".join(map(str, node_indices)) + "\n"

        # Facet for the air-earth interface with the receivers

        if any(node[3] > 0.0 for node in self.node_list_source):
            # interface
            num_polygon = 0
            self.num_facet += 1
            if self.num_layers > 0:
                interface_string = ""
                interface_face_nodes_list = []
                node_list_air_earth_interface = self.node_list_of_list_of_interfaces[0]
                
                num_polygon += 1
                for node in sorted(node_list_air_earth_interface):
                    interface_face_nodes_list.append(node[0])

                interface_n_string = " ".join(map(str, interface_face_nodes_list)) + "\n"
                str_1 = "# air earth interface\n4 "
                interface_string += str_1 + interface_n_string

            n_receivers_string = ""
            # receivers triangles
            for i in range(len(self.node_list_receiver)//3):
                num_polygon += 1
                receiver_triangle_nodes_list = [3]
                for receiver_edge in self.node_list_receiver[i*3:i*3+3]:
                    receiver_triangle_nodes_list.append(receiver_edge[0])
                receiver_triangle_nodes_string = " ".join(map(str, receiver_triangle_nodes_list)) + "\n"
                n_receivers_string += receiver_triangle_nodes_string

            num_polygon += 2 if self.box_present else 0 # for the source box
            interface_facet_information = [num_polygon, 0, 99]
            interface_base_string = " ".join(map(str, interface_facet_information)) + "\n"

            if self.box_present:
                if self.num_layers > 0:
                    self.interface_string = "# air-earth facet\n" + interface_base_string + interface_string + "# receiver triangles\n" + "".join(n_receivers_string) + "\n# source box\n" + source_box_z0_string
                    # remove this layer from the list of interfaces
                    self.node_list_of_list_of_interfaces_old = self.node_list_of_list_of_interfaces
                    self.node_list_of_list_of_interfaces = self.node_list_of_list_of_interfaces[1:]
                else:
                    self.interface_string = "# Tx-Rx facet\n" + interface_base_string + "# receiver triangles\n" + "".join(n_receivers_string) + "\n# source box\n" + source_box_z0_string
            else:
                if self.num_layers > 0:
                    self.interface_string = "# air-earth facet\n" + interface_base_string + interface_string + "# receiver triangles\n" + "".join(n_receivers_string)
                    # remove this layer from the list of interfaces
                    self.node_list_of_list_of_interfaces_old = self.node_list_of_list_of_interfaces
                    self.node_list_of_list_of_interfaces = self.node_list_of_list_of_interfaces[1:]
                else:
                    self.interface_string = "# Tx-Rx facet\n" + interface_base_string + "# receiver triangles\n" + "".join(n_receivers_string)
            # source
            self.num_facet += 1
            num_polygon = 1
            source_face_nodes_list = []
            for node in sorted(self.node_list_source):
                source_face_nodes_list.append(node[0])
            num_source_points = len(source_face_nodes_list)
            source_string = f"{num_source_points} " + " ".join(map(str, source_face_nodes_list)) + "\n"
            self.source_string = f"# source facet\n1 0 987\n{source_string}\n"

        else:
            # interface
            num_polygon = 0
            self.num_facet += 1
            if self.num_layers > 0:
                interface_string = ""
                interface_face_nodes_list = []
                node_list_air_earth_interface = self.node_list_of_list_of_interfaces[0]
                
                num_polygon += 1
                for node in sorted(node_list_air_earth_interface):
                    interface_face_nodes_list.append(node[0])

                interface_n_string = " ".join(map(str, interface_face_nodes_list)) + "\n"
                str_1 = "# air earth interface\n4 "
                interface_string += str_1 + interface_n_string

            n_receivers_string = ""
            # receivers triangles
            for i in range(len(self.node_list_receiver)//3):
                num_polygon += 1
                receiver_triangle_nodes_list = [3]
                for receiver_edge in self.node_list_receiver[i*3:i*3+3]:
                    receiver_triangle_nodes_list.append(receiver_edge[0])
                receiver_triangle_nodes_string = " ".join(map(str, receiver_triangle_nodes_list)) + "\n"
                n_receivers_string += receiver_triangle_nodes_string

            # source
            num_polygon += 1
            source_face_nodes_list = []
            for node in sorted(self.node_list_source):
                source_face_nodes_list.append(node[0])
            num_source_points = len(source_face_nodes_list)
            source_string = f"{num_source_points} " + " ".join(map(str, source_face_nodes_list)) + "\n"

            num_polygon += 2 if self.box_present else 0 # for the source box
            interface_facet_information = [num_polygon, 0, 99]
            interface_base_string = " ".join(map(str, interface_facet_information)) + "\n"

            if self.num_layers > 0:
                self.interface_string = "# air-earth facet\n" + interface_base_string + interface_string + "# receiver triangles\n" + "".join(n_receivers_string) + "# source nodes\n" + "".join(source_string)
                # remove this layer from the list of interfaces
                self.node_list_of_list_of_interfaces_old = self.node_list_of_list_of_interfaces
                self.node_list_of_list_of_interfaces = self.node_list_of_list_of_interfaces[1:]
            else:
                self.interface_string = "# Tx-Rx facet\n" + interface_base_string + "# receiver triangles\n" + "".join(n_receivers_string) + "# source nodes\n" + "".join(source_string)
            
            if self.box_present:
                self.interface_string += "# source box\n" + source_box_z0_string

        # Facets for the interfaces other than the air-earth interface
        if self.num_layers > 0: 
            self.num_facet += len(self.node_list_of_list_of_interfaces)
            marker_interfaces = [6, 7, 8, 9, 10]
            num_polygon = 0
            interfaces_n_string = ""
            i = 0
            for interface_nodes in self.node_list_of_list_of_interfaces:
                interface_n_string = ""
                num_polygon += 1
                interface_face_nodes_list = []

                for node in sorted(interface_nodes):
                    interface_face_nodes_list.append(node[0])
                interface_n_string += " ".join(map(str, interface_face_nodes_list)) + "\n"

                if num_polygon % 2 != 0:
                    str_1 = "# interface layer {}\n".format(i+1)
                else:
                    i += 1
                    str_1 = "# transition layer {}\n".format(i+1)

                interfaces_n_string += str_1 + f"1 0 {marker_interfaces[num_polygon-1]} \n4 " + interface_n_string

            self.interfaces_string = "# other interfaces and transitions\n" + interfaces_n_string
        else:
            self.interfaces_string = ""
        
        if self.box_present:
            self.source_box_string = "# source box facets\n"
            for i, face in enumerate(source_box_1_faces+source_box_2_faces):
                self.source_box_string += f"1 0 {self.marker_source_box}\n"
                self.source_box_string += "4 " + " ".join(map(str, [node[0] for node in face])) + "\n"
                self.num_facet += 1

        # PML Layers
        if self.NUM_PML > 0:
            self.evaluate_all_facets_PML()


    def evaluate_all_facets_PML(self):
        """
        Function to evaluate all three types of facets for each PML layer.
        All these facets will be rectangular prisms.
        """

        t = self.pml_thickness
        self.pml_string = "# PML Layers\n"
        
        all_z = self.all_z
        if self.num_layers > 0:
            all_face_nodes = self.node_list_domain_prism + [
                        node for node_list in self.node_list_of_list_of_interfaces_old for node in node_list
                    ] + self.node_list_PML_2 + self.node_list_PML_3
        else:
            all_face_nodes = self.node_list_domain_prism + self.node_list_PML_2 + self.node_list_PML_3
        
        xe = [self.x_min, self.x_max]
        ye = [self.y_min, self.y_max]

        # Type 1: Extend each face of the domain: 1 new face per original face
        # Assembling facets per layer
        for i in range(self.NUM_PML):
            self.pml_string += f"\n# PML Layer {i+1} - Type 1:\n"
            faces = ["Front", "Right", "Back", "Left", "Top", "Bottom"]
            
            all_face_filters = [
                ("Front", lambda node: (node[1] == self.x_min or node[1] == self.x_max) and (node[2] == self.y_min - i * t or node[2] == self.y_min - (i+1) * t) and (node[3] >= self.z_min and node[3] <= self.z_max)),
                ("Right", lambda node: (node[1] == self.x_max + i * t or node[1] == self.x_max + (i+1) * t) and (node[2] == self.y_min or node[2] == self.y_max) and (node[3] >= self.z_min and node[3] <= self.z_max)),
                ("Back", lambda node: (node[1] == self.x_min or node[1] == self.x_max) and (node[2] == self.y_max + i * t or node[2] == self.y_max + (i+1) * t) and (node[3] >= self.z_min and node[3] <= self.z_max)),
                ("Left", lambda node: (node[1] == self.x_min - (i+1) * t or node[1] == self.x_min - i * t) and (node[2] == self.y_min or node[2] == self.y_max) and (node[3] >= self.z_min and node[3] <= self.z_max)),
                ("Top", lambda node: (node[2] == self.y_min or node[2] == self.y_max) and (node[1] == self.x_min or node[1] == self.x_max) and (node[3] == self.z_max + i * t or node[3] == self.z_max + (i+1) * t)),
                ("Bottom", lambda node: (node[2] == self.y_min or node[2] == self.y_max) and (node[1] == self.x_min or node[1] == self.x_max) and (node[3] == self.z_min - i * t or node[3] == self.z_min - (i+1) * t))
            ]

            # Front face: sort by x, z for the main face, then by z for each interface
            front_face_sorting_settings = [(self.y_min - (i+1)*t, 2, 1, 3)] + [
                (z, 3, 1, 2) for z in self.layer_interfaces['z']
            ]

            right_face_sorting_settings = [(self.x_max + (i+1)*t, 1, 2, 3)] + [
                (z, 3, 1, 2) for z in self.layer_interfaces['z']
            ]

            back_face_sorting_settings = [(self.y_max + (i+1)*t, 2, 1, 3)] + [
                (z, 3, 1, 2) for z in self.layer_interfaces['z']
            ]

            left_face_sorting_settings = [(self.x_min - (i+1)*t, 1, 2, 3)] + [
                (z, 3, 1, 2) for z in self.layer_interfaces['z']
            ]

            top_face_sorting_settings = [
                (self.z_max + (i+1)*t, 3, 1, 2),  # Top face: sort by z, x
            ]

            bottom_face_sorting_settings = [
                (self.z_min - (i+1)*t, 3, 1, 2),  # Bottom face: sort by z, x
            ]

            pml_type_1_faces_nodes_list = []
            it = 0
            for _, filter_func in all_face_filters:
                it += 1
                face_nodes = [node for node in all_face_nodes if filter_func(node)]
                pml_type_1_faces_nodes_list.append([face_nodes])              

            if self.num_layers > 0:
                combined_interfaces_nodes = [node for node_list in self.node_list_of_list_of_interfaces_old for node in node_list] + self.pml_type_2_layer_nodes
            else:
                combined_interfaces_nodes = []
            pml_type_1_interfaces_nodes_list = []
            for _, filter_func in all_face_filters:
                interface_nodes = [node for node in combined_interfaces_nodes if filter_func(node)]
                pml_type_1_interfaces_nodes_list.append([interface_nodes])

            for face in range(len(faces)):
                if face == 0:
                    face_sorting_settings = front_face_sorting_settings
                    self.pml_string += f"\n# 1+l new face with interface segments for the Front face\n"
                elif face == 1:
                    face_sorting_settings = right_face_sorting_settings
                    self.pml_string += f"\n# 1+l new face with interface segments for the Right face\n"
                elif face == 2:
                    face_sorting_settings = back_face_sorting_settings
                    self.pml_string += f"\n# 1+l new face with interface segments for the Back face\n"
                elif face == 3:
                    face_sorting_settings = left_face_sorting_settings
                    self.pml_string += f"\n# 1+l new face with interface segments for the Left face\n"
                elif face == 4:
                    face_sorting_settings = top_face_sorting_settings
                    self.pml_string += f"\n# 1 new face with interface segments for the Top face\n"
                    pml_type_1_interfaces_nodes_list[face][0] = None
                elif face == 5:
                    face_sorting_settings = bottom_face_sorting_settings
                    self.pml_string += f"\n# 1 new face with interface segments for the Bottom face\n"
                    pml_type_1_interfaces_nodes_list[face][0] = None

                marker = self.marker_pml_1[i*6 + face]
                for new_face in range(len(face_sorting_settings)): 
                    self.num_facet += 1
                    self.pml_string += self.create_face_string(marker, pml_type_1_faces_nodes_list[face][0], pml_type_1_interfaces_nodes_list[face][0], *face_sorting_settings[new_face])
                
        # __________________________________________________________
        # Create 12 cuboids at each edge of the domain - Type 2
        # __________________________________________________________
        n = self.NUM_PML
        if n > 1:
            raise ValueError("The number of PML layers must be 1.")
        
        n_list = np.linspace(t, n*t, n)
        n_pairs = [(i, j) for i in n_list for j in n_list]
        n_b_list = np.linspace(0, (n-1)*t, n)
        n_b_pairs = [(i, j) for i in n_b_list for j in n_b_list]
        n_s_list = np.linspace(0, n*t, n+1)
        n_s_pairs = [(i, j) for i in n_s_list for j in n_s_list]

        all_edge_filters = [
            ("Front-Right", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[1] >= self.x_max + bi and node[1] <= self.x_max + i) and (node[2] <= self.y_min - bj and node[2] >= self.y_min - j) and (node[3] >= self.z_min and node[3] <= self.z_max))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Back-Right", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[1] >= self.x_max + bi and node[1] <= self.x_max + i) and (node[2] >= self.y_max + bj and node[2] <= self.y_max + j) and (node[3] >= self.z_min and node[3] <= self.z_max))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Back-Left", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[1] <= self.x_min - bi and node[1] >= self.x_min - i) and (node[2] >= self.y_max + bj and node[2] <= self.y_max + j) and (node[3] >= self.z_min and node[3] <= self.z_max))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Front-Left", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[1] <= self.x_min - bi and node[1] >= self.x_min - i) and (node[2] <= self.y_min - bj and node[2] >= self.y_min - j) and (node[3] >= self.z_min and node[3] <= self.z_max))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        
        ] + [("Front-Top", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[2] <= self.y_min - bj and node[2] >= self.y_min - j) and (node[1] >= self.x_min and node[1] <= self.x_max) and (node[3] >= self.z_max + bi and node[3] <= self.z_max + i))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Right-Top", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[2] >= self.y_min and node[2] <= self.y_max) and (node[1] >= self.x_max + bi and node[1] <= self.x_max + i) and (node[3] >= self.z_max + bj and node[3] <= self.z_max + j))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Back-Top", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[2] >= self.y_max + bj and node[2] <= self.y_max + j) and (node[1] >= self.x_min and node[1] <= self.x_max) and (node[3] >= self.z_max + bi and node[3] <= self.z_max + i))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Left-Top", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[2] >= self.y_min and node[2] <= self.y_max) and (node[1] <= self.x_min - bi and node[1] >= self.x_min - i) and (node[3] >= self.z_max + bj and node[3] <= self.z_max + j))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        
        ] + [("Front-Bottom", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[2] <= self.y_min - bj and node[2] >= self.y_min - j) and (node[1] >= self.x_min and node[1] <= self.x_max) and (node[3] <= self.z_min - bi and node[3] >= self.z_min - i))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Right-Bottom", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[2] >= self.y_min and node[2] <= self.y_max) and (node[1] >= self.x_max + bi and node[1] <= self.x_max + i) and (node[3] <= self.z_min - bj and node[3] >= self.z_min - j))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Back-Bottom", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[2] >= self.y_max + bj and node[2] <= self.y_max + j) and (node[1] >= self.x_min and node[1] <= self.x_max) and (node[3] <= self.z_min - bi and node[3] >= self.z_min - i))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Left-Bottom", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[2] >= self.y_min and node[2] <= self.y_max) and (node[1] <= self.x_min - bi and node[1] >= self.x_min - i) and (node[3] <= self.z_min - bj and node[3] >= self.z_min - j))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)]
        
        # Lateral cuboids

        front_right_sorting_settings = [val for pair in zip(
            [(self.x_max + i, 1, 2, 3) for i, _ in n_s_pairs], # Right faces: sort by y, z
            [(self.y_min - j, 2, 1, 3) for _, j in n_s_pairs], # Front faces: sort by x, z
            [(z, 3, 1, 2) for z in all_z], # z faces
        ) for val in pair]

        back_right_sorting_settings = [val for pair in zip(
            [(self.x_max + i, 1, 2, 3) for i, _ in n_s_pairs],  # Right face: sort by y, z
            [(self.y_max + j, 2, 1, 3) for _, j in n_s_pairs],   # Back face: sort by x, z
            [(z, 3, 1, 2) for z in all_z], # z faces
        ) for val in pair]

        back_left_sorting_settings = [val for pair in zip(
            [(self.x_min - i, 1, 2, 3) for i, _ in n_s_pairs],  # Left face: sort by y, z
            [(self.y_max + j, 2, 1, 3) for _, j in n_s_pairs],   # Back face: sort by x, z
            [(z, 3, 1, 2) for z in all_z], # z faces
            ) for val in pair]

        front_left_sorting_settings = [val for pair in zip(
            [(self.x_min - i, 1, 2, 3) for i, _ in n_s_pairs],  # Left face: sort by y, z
            [(self.y_min - j, 2, 1, 3) for _, j in n_s_pairs],   # Front face: sort by x, z
            [(z, 3, 1, 2) for z in all_z], # z faces
        ) for val in pair]

        # "Horizontal" cuboids (along the top and bottom faces)

        front_top_sorting_settings = [val for pair in zip(
            [(self.y_min - j, 2, 1, 3) for _, j in n_s_pairs],  # Front face: sort by x, z
            [(self.z_max + j, 3, 1, 2) for _, j in n_s_pairs],   # Top face: sort by z, x
            [(x, 1, 2, 3) for x in xe],  # x faces
        ) for val in pair]

        right_top_sorting_settings = [val for pair in zip(
            [(self.x_max + j, 1, 2, 3) for _, j in n_s_pairs],  # Right face: sort by y, z
            [(self.z_max + j, 3, 1, 2) for _, j in n_s_pairs],   # Top face: sort by z, x
            [(y, 2, 1, 3) for y in ye],  # y faces
            ) for val in pair]

        back_top_sorting_settings = [val for pair in zip(
            [(self.y_max + j, 2, 1, 3) for _, j in n_s_pairs],  # Back face: sort by x, z
            [(self.z_max + j, 3, 1, 2) for _, j in n_s_pairs],   # Top face: sort by z, x
            [(x, 1, 2, 3) for x in xe],  # x faces
        ) for val in pair]

        left_top_sorting_settings = [val for pair in zip(
            [(self.x_min - j, 1, 2, 3) for _, j in n_s_pairs],  # Left face: sort by y, z
            [(self.z_max + j, 3, 1, 2) for _, j in n_s_pairs],   # Top face: sort by z, x
            [(y, 2, 1, 3) for y in ye],  # y faces
            ) for val in pair]

        front_bottom_sorting_settings = [val for pair in zip(
            [(self.y_min - j, 2, 1, 3) for _, j in n_s_pairs],  # Front face: sort by x, z
            [(self.z_min - j, 3, 1, 2) for _, j in n_s_pairs],   # Bottom face: sort by z, x
            [(x, 1, 2, 3) for x in xe],  # x faces
            ) for val in pair]

        right_bottom_sorting_settings = [val for pair in zip(
            [(self.x_max + j, 1, 2, 3) for _, j in n_s_pairs],  # Right face: sort by y, z
            [(self.z_min - j, 3, 1, 2) for _, j in n_s_pairs],   # Bottom face: sort by z, x
            [(y, 2, 1, 3) for y in ye],  # y faces
            ) for val in pair]

        back_bottom_sorting_settings = [val for pair in zip(
            [(self.y_max + j, 2, 1, 3) for _, j in n_s_pairs],  # Back face: sort by x, z
            [(self.z_min - j, 3, 1, 2) for _, j in n_s_pairs],   # Bottom face: sort by z, x
            [(x, 1, 2, 3) for x in xe],  # x faces
            ) for val in pair]

        left_bottom_sorting_settings = [val for pair in zip(
            [(self.x_min - j, 1, 2, 3) for _, j in n_s_pairs],  # Left face: sort by y, z
            [(self.z_min - j, 3, 1, 2) for _, j in n_s_pairs],   # Bottom face: sort by z, x
            [(y, 2, 1, 3) for y in ye],  # y faces
            ) for val in pair]

        # PML Type 2 - Create 12 cuboids at each edge of the domain
        self.pml_string += f"\n# PML Layer - Type 2:\n"

        pml_type_2_faces_nodes_list = []
        for _, filter_func in all_edge_filters:
            edge_nodes = [node for node in all_face_nodes if filter_func(node)]
            pml_type_2_faces_nodes_list.append([edge_nodes])

        if self.num_layers > 0:
            combined_interfaces_nodes = [node for node_list in self.node_list_of_list_of_interfaces_old for node in node_list] + self.pml_type_2_layer_nodes
        else:
            combined_interfaces_nodes = []

        pml_type_2_interfaces_nodes_list = []
        for _, filter_func in all_edge_filters:
            interface_nodes = [node for node in combined_interfaces_nodes if filter_func(node)]
            pml_type_2_interfaces_nodes_list.append([interface_nodes])

        int_idx = 0
        int_markers = [99-n-self.num_layers for n in range(self.num_layers*4)]
        
        for edge in range(len(all_edge_filters)):
            # FRONT-RIGHT, BACK-RIGHT, BACK-LEFT, FRONT-LEFT, FRONT-TOP, RIGHT-TOP, BACK-TOP, LEFT-TOP, FRONT-BOTTOM, RIGHT-BOTTOM, BACK-BOTTOM, LEFT-BOTTOM
            if edge == 0:
                edge_sorting_settings = front_right_sorting_settings
                self.pml_string += f"\n# 6+l new faces with interface segments for the Front-Right edge\n"
            elif edge == len(n_pairs):
                edge_sorting_settings = back_right_sorting_settings
                self.pml_string += f"\n# 6+l new faces with interface segments for the Back-Right edge\n"
            elif edge == 2*len(n_pairs):
                edge_sorting_settings = back_left_sorting_settings
                self.pml_string += f"\n# 6+l new faces with interface segments for the Back-Left edge\n"
            elif edge == 3*len(n_pairs):
                edge_sorting_settings = front_left_sorting_settings
                self.pml_string += f"\n# 6+l new faces with interface segments for the Front-Left edge\n"
            elif edge == 4*len(n_pairs):
                edge_sorting_settings = front_top_sorting_settings
                self.pml_string += f"\n# 6 new faces with interface segments for the Front-Top edge\n"
            elif edge == 5*len(n_pairs):
                edge_sorting_settings = right_top_sorting_settings
                self.pml_string += f"\n# 6 new faces with interface segments for the Right-Top edge\n"
            elif edge == 6*len(n_pairs):
                edge_sorting_settings = back_top_sorting_settings
                self.pml_string += f"\n# 6 new faces with interface segments for the Back-Top edge\n"
            elif edge == 7*len(n_pairs):
                edge_sorting_settings = left_top_sorting_settings
                self.pml_string += f"\n# 6 new faces with interface segments for the Left-Top edge\n"
            elif edge == 8*len(n_pairs):
                edge_sorting_settings = front_bottom_sorting_settings
                self.pml_string += f"\n# 6 new faces with interface segments for the Front-Bottom edge\n"
            elif edge == 9*len(n_pairs):
                edge_sorting_settings = right_bottom_sorting_settings
                self.pml_string += f"\n# 6 new faces with interface segments for the Right-Bottom edge\n"
            elif edge == 10*len(n_pairs):
                edge_sorting_settings = back_bottom_sorting_settings
                self.pml_string += f"\n# 6 new faces with interface segments for the Back-Bottom edge\n"
            elif edge == 11*len(n_pairs):
                edge_sorting_settings = left_bottom_sorting_settings
                self.pml_string += f"\n# 6 new faces with interface segments for the Left-Bottom edge\n"

            if edge >= 4*len(n_pairs):
                pml_type_2_interfaces_nodes_list[edge][0] = None
                            
            # Use int_markers for faces corresponding to z value filters for edges less than 4*n_pairs
            if edge < 4 * len(n_pairs):
                # edge_sorting_settings contains (axis_value, axis_index, sort_axis_1, sort_axis_2) for each face
                # The z-faces are every 3rd in edge_sorting_settings, starting at index 2
                marker_list = []
                for new_face in range(len(edge_sorting_settings)):
                    axis_value, axis_index, _, _ = edge_sorting_settings[new_face]
                    # If this is a z-face (axis_index == 3), assign int_markers by z value
                    if axis_index == 3 and axis_value != self.z_min and axis_value != self.z_max:
                        marker = int_markers[int_idx]
                        int_idx += 1
                    else:
                        marker = self.marker_pml_2[edge]
                    marker_list.append(marker)
                all_face_strings = [
                    self.create_face_string(marker_list[i], pml_type_2_faces_nodes_list[edge][0], pml_type_2_interfaces_nodes_list[edge][0], *edge_sorting_settings[i])
                    for i in range(len(edge_sorting_settings))
                ]
                unique_face_strings = list(set(all_face_strings))
                for face_str in unique_face_strings:
                    self.num_facet += 1
                    self.pml_string += face_str
                    all_face_strings = [
                        self.create_face_string(marker, pml_type_2_faces_nodes_list[edge][0], pml_type_2_interfaces_nodes_list[edge][0], axis_value, axis_index, sort_axis_1, sort_axis_2)
                        for axis_value, axis_index, sort_axis_1, sort_axis_2 in edge_sorting_settings
                    ]
            else:
                marker = self.marker_pml_2[edge] 
                all_component_strings = []
                for new_face in range(len(edge_sorting_settings)):
                    face_string = self.create_face_string(
                        marker, 
                        pml_type_2_faces_nodes_list[edge][0], 
                        pml_type_2_interfaces_nodes_list[edge][0], 
                        *edge_sorting_settings[new_face]
                    )
                    all_component_strings.append(face_string)
                unique_component_strings = list(set(all_component_strings))
                for new_face in range(len(unique_component_strings)):
                    self.num_facet += 1
                    self.pml_string += unique_component_strings[new_face]

        # Type 3: Extend each vertex of the domain: n^2 new cubes per original vertex
        # Assembling facets for all layers in one loop

        n_triplets = [(i, j, k) for i in n_list for j in n_list for k in n_list]
        n_b_triplets = [(i, j, k) for i in n_b_list for j in n_b_list for k in n_b_list]
        
        self.pml_string += f"\n# PML Layer - Type 3:\n"

        all_vertex_filters = [
            ("Front-Right-Top", (self.x_max, self.y_min, self.z_max), lambda node, i=i, j=j, k=k, bi=bi, bj=bj, bk=bk: ((node[1] >= self.x_max + bi and node[1] <= self.x_max + i) and (node[2] <= self.y_min - bj and node[2] >= self.y_min - j) and (node[3] >= self.z_max + bk and node[3] <= self.z_max + k))) for (i, j, k), (bi, bj, bk) in zip(n_triplets, n_b_triplets)
        ] + [("Back-Right-Top", (self.x_max, self.y_max, self.z_max), lambda node, i=i, j=j, k=k, bi=bi, bj=bj, bk=bk: ((node[1] >= self.x_max + bi and node[1] <= self.x_max + i) and (node[2] >= self.y_max + bj and node[2] <= self.y_max + j) and (node[3] >= self.z_max + bk and node[3] <= self.z_max + k))) for (i, j, k), (bi, bj, bk) in zip(n_triplets, n_b_triplets)
        ] + [("Back-Left-Top", (self.x_min, self.y_max, self.z_max), lambda node, i=i, j=j, k=k, bi=bi, bj=bj, bk=bk: ((node[1] <= self.x_min - bi and node[1] >= self.x_min - i) and (node[2] >= self.y_max + bj and node[2] <= self.y_max + j) and (node[3] >= self.z_max + bk and node[3] <= self.z_max + k))) for (i, j, k), (bi, bj, bk) in zip(n_triplets, n_b_triplets)
        ] + [("Front-Left-Top", (self.x_min, self.y_min, self.z_max), lambda node, i=i, j=j, k=k, bi=bi, bj=bj, bk=bk: ((node[1] <= self.x_min - bi and node[1] >= self.x_min - i) and (node[2] <= self.y_min - bj and node[2] >= self.y_min - j) and (node[3] >= self.z_max + bk and node[3] <= self.z_max + k))) for (i, j, k), (bi, bj, bk) in zip(n_triplets, n_b_triplets)
        ] + [("Front-Right-Bottom", (self.x_max, self.y_min, self.z_min), lambda node, i=i, j=j, k=k, bi=bi, bj=bj, bk=bk: ((node[1] >= self.x_max + bi and node[1] <= self.x_max + i) and (node[2] <= self.y_min - bj and node[2] >= self.y_min - j) and (node[3] <= self.z_min - bk and node[3] >= self.z_min - k))) for (i, j, k), (bi, bj, bk) in zip(n_triplets, n_b_triplets)
        ] + [("Back-Right-Bottom", (self.x_max, self.y_max, self.z_min), lambda node, i=i, j=j, k=k, bi=bi, bj=bj, bk=bk: ((node[1] >= self.x_max + bi and node[1] <= self.x_max + i) and (node[2] >= self.y_max + bj and node[2] <= self.y_max + j) and (node[3] <= self.z_min - bk and node[3] >= self.z_min - k))) for (i, j, k), (bi, bj, bk) in zip(n_triplets, n_b_triplets)
        ] + [("Back-Left-Bottom", (self.x_min, self.y_max, self.z_min), lambda node, i=i, j=j, k=k, bi=bi, bj=bj, bk=bk: ((node[1] <= self.x_min - bi and node[1] >= self.x_min - i) and (node[2] >= self.y_max + bj and node[2] <= self.y_max + j) and (node[3] <= self.z_min - bk and node[3] >= self.z_min - k))) for (i, j, k), (bi, bj, bk) in zip(n_triplets, n_b_triplets)
        ] + [("Front-Left-Bottom", (self.x_min, self.y_min, self.z_min), lambda node, i=i, j=j, k=k, bi=bi, bj=bj, bk=bk: ((node[1] <= self.x_min - bi and node[1] >= self.x_min - i) and (node[2] <= self.y_min - bj and node[2] >= self.y_min - j) and (node[3] <= self.z_min - bk and node[3] >= self.z_min - k))) for (i, j, k), (bi, bj, bk) in zip(n_triplets, n_b_triplets)]

        pml_type_3_faces_nodes_list = []
        for _, _, filter_func in all_vertex_filters:
            vertex_nodes = [node for node in all_face_nodes if filter_func(node)]
            # Only keep unique vertex nodes (by coordinates)
            unique_vertex_nodes = []
            seen = set()
            for node in vertex_nodes:
                key = (round(node[1], 8), round(node[2], 8), round(node[3], 8))
                if key not in seen:
                    seen.add(key)
                    unique_vertex_nodes.append(node)
            vertex_nodes = unique_vertex_nodes
            pml_type_3_faces_nodes_list.append([vertex_nodes])

        all_cuboid_faces = []
        not_needed_nodes = [n for _, n, _ in all_vertex_filters]
        for i in range(len(pml_type_3_faces_nodes_list)):
            faces = self.create_half_cuboid_faces_from_nodes(pml_type_3_faces_nodes_list[i][0], not_needed_nodes[i])
            all_cuboid_faces.append(faces)

        # Remove duplicates from all_cuboid_faces
        unique_faces = []
        unique_markers = []
        unique_sorted_faces = []
        it = 0
        for faces in all_cuboid_faces:
            for face in faces:
                # Sort the face by ascending order of its entries (node numbers)
                sorted_face = sorted(face, key=lambda node: node[0])
                if sorted_face not in unique_sorted_faces:
                    unique_faces.append(face)
                    unique_sorted_faces.append(sorted_face)
                    unique_markers.append(self.marker_pml_3[it])
                
            it += 1
        

        for i, face in enumerate(unique_faces):
            self.pml_string += f"1 0 {unique_markers[i]}\n"
            self.pml_string += "4 " + " ".join(map(str, [node[0] for node in face])) + "\n"
            self.num_facet += 1


    #######################################################################################
    # Evaluating the .poly file aspect - Regions, as well as the regionparameters.txt file
    #######################################################################################
    def evaluate_all_regions(self):
        """
        Function to evaluate all regions in the model.
        Calls the function to evaluate all regions in the PML domain as well.
        """

        # All different material regions
        self.num_regions = 0
        regions = []
        for i in range(self.num_layers + 1):
            region_height = 0
            if i == 0:
                region_height = self.layer_interfaces['z'][0] + abs(self.layer_interfaces['z'][0] - self.z_max) / 2
                r_label = "# Air"
            elif i == self.num_layers:
                region_height = self.layer_interfaces['z'][-1] - abs(self.layer_interfaces['z'][-1] - self.z_min) / 2
                r_label = "# Last Earth Layer"
            else:
                region_height = self.layer_interfaces['z'][i - 1] - abs(self.layer_interfaces['z'][i - 1] - self.layer_interfaces['z'][i]) / 2
                r_label = f"# Earth Layer {i}"
            self.num_regions += 1
            earth_region = [
                self.num_regions, 0, 0, round(region_height, 5), self.num_regions, self.region_volumes[i], r_label,
                self.region_rho[i], self.region_mu_r[i], self.region_epsilon_r[i]
            ]
            regions.append(earth_region)
        
        # Source box region
        if self.box_present:
            if self.num_layers == 0:
                region_height = 0.0 + float(self.region_volumes[0]) * 5
                self.num_regions += 1
                vol = float(self.region_volumes[0])/self.m**3
                source_region = [
                    self.num_regions, 0, 0, round(region_height, 5), self.num_regions, f"{vol:.4e}", "# Source Box 1",
                    self.region_rho[0], self.region_mu_r[0], self.region_epsilon_r[0]
                ]
                regions.append(source_region)
                region_height = 0.0 - float(self.region_volumes[0]) * 5
                self.num_regions += 1
                source_region = [
                    self.num_regions, 0, 0, round(region_height, 5), self.num_regions, f"{vol:.4e}", "# Source Box 2",
                    self.region_rho[0], self.region_mu_r[0], self.region_epsilon_r[0]
                ]
                regions.append(source_region)
            else: # Source box region with two different materials
                region_height = 0.0 + float(self.region_volumes[0]) * 5
                self.num_regions += 1
                vol = float(self.region_volumes[0])/self.m**3
                source_region = [
                    self.num_regions, 0, 0, round(region_height, 5), self.num_regions, f"{vol:.4e}", "# Source Box 1",
                    self.region_rho[0], self.region_mu_r[0], self.region_epsilon_r[0]
                ]
                regions.append(source_region)
                region_height = 0.0 - float(self.region_volumes[1]) * 5
                self.num_regions += 1
                vol = float(self.region_volumes[1])/self.m**3
                source_region = [
                    self.num_regions, 0, 0, round(region_height, 5), self.num_regions, f"{vol:.4e}", "# Source Box 2",
                    self.region_rho[1], self.region_mu_r[1], self.region_epsilon_r[1]
                ]
                regions.append(source_region)

        self.regions = regions

        # PML regions
        if self.NUM_PML == 0:
            return
        else:
            self.PML_regions = self.evaluate_all_regions_PML()
            self.regions += self.PML_regions
    

    def evaluate_all_regions_PML(self):
        """
        Function to evaluate all PML regions.
        For each PML layer now, there will be 1 region per prism.
        Type 1: 6 prisms per layer
        Type 2: 12 prisms per layer
        Type 3: 8 prisms per layer

        Instead of region_height, this produces region_coordinates
        that correspond to the midpoint of each prism—and then
        assigns the correct [volume, rho, μr, εr] purely by that z.
        """
        if self.num_layers > 0:
            all_z = sorted([self.z_min] + self.layer_interfaces['z'] + [self.z_max])
        else:
            all_z = [self.z_min, self.z_max]
        self.all_z = all_z

        PML_regions = []
        t = self.pml_thickness
        n = self.NUM_PML

        # precompute offsets
        n_list     = np.linspace(t/2, n*t/2, n)
        n_pairs    = [(i, j) for i in n_list for j in n_list]
        n_triplets = [(i, j, k) for i in n_list for j in n_list for k in n_list]

        # all the z‐interfaces and mid‐slab heights
        midheights = [(self.all_z[j] + self.all_z[j+1]) / 2
                    for j in range(len(self.all_z)-1)]

        # prepare a lookup table of [vol, rho, mu_r, eps_r] per layer index
        materials = [
            [self.region_volumes[i],
            self.region_rho[i],
            self.region_mu_r[i],
            self.region_epsilon_r[i] ]
            for i in range(self.num_layers+1)
        ]

        # ─── Type 1: face extensions ───────────────────────────────────────────────
        self.marker_pml_1 = []
        for i in range(self.NUM_PML):
            # For each region (layer), create a face at the correct midheight for Front/Right/Back/Left
            faces = [
                ("Front", [0.0, self.y_min - (i + 1) * t / 2, 0.0]),
                ("Right", [self.x_max + (i + 1) * t / 2, 0.0, 0.0]),
                ("Back", [0.0, self.y_max + (i + 1) * t / 2, 0.0]),
                ("Left", [self.x_min - (i + 1) * t / 2, 0.0, 0.0]),
                ("Top", [0.0, 0.0, self.z_max + (i + 1) * t / 2]),
                ("Bottom", [0.0, 0.0, self.z_min - (i + 1) * t / 2])
            ]

            
            for face_name, midpoint in faces:
                r_label = f"# PML Layer {i+1} - Type 1: {face_name}"
                if face_name in ("Top", "Bottom"):
                    self.num_regions += 1
                    self.marker_pml_1.append(self.num_regions)
                    PML_regions.append([
                        self.num_regions,
                        round(midpoint[0],5),
                        round(midpoint[1],5),
                        round(midpoint[2],5),
                        self.marker_pml_1[-1],
                        None,  # placeholder for volume
                        r_label,
                        None, None, None
                    ])
                else:
                    for l in range(self.num_layers+1):
                        self.num_regions += 1
                        self.marker_pml_1.append(self.num_regions)
                        PML_regions.append([
                            self.num_regions,
                            round(midpoint[0],5),
                            round(midpoint[1],5),
                            round(midheights[l],5),
                            self.marker_pml_1[-1],
                            None,
                            r_label,
                            None, None, None
                        ])

        # ─── Type 2: edge extensions ───────────────────────────────────────────────
        self.marker_pml_2 = []
        edges = (
            [("Front-Right",  [self.x_max + i, self.y_min - j, 0.0]) for i,j in n_pairs] +
            [("Back-Right",   [self.x_max + i, self.y_max + j, 0.0]) for i,j in n_pairs] +
            [("Back-Left",    [self.x_min - i, self.y_max + j, 0.0]) for i,j in n_pairs] +
            [("Front-Left",   [self.x_min - i, self.y_min - j, 0.0]) for i,j in n_pairs] +
            [("Front-Top",    [0.0, self.y_min - j, self.z_max + i]) for i,j in n_pairs] +
            [("Right-Top",    [self.x_max + i, 0.0, self.z_max + j]) for i,j in n_pairs] +
            [("Back-Top",     [0.0, self.y_max + j, self.z_max + i]) for i,j in n_pairs] +
            [("Left-Top",     [self.x_min - i, 0.0, self.z_max + j]) for i,j in n_pairs] +
            [("Front-Bottom", [0.0, self.y_min - j, self.z_min - i]) for i,j in n_pairs] +
            [("Right-Bottom", [self.x_max + j, 0.0, self.z_min - i]) for i,j in n_pairs] +
            [("Back-Bottom",  [0.0, self.y_max + j, self.z_min - i]) for i,j in n_pairs] +
            [("Left-Bottom",  [self.x_min - j, 0.0, self.z_min - i]) for i,j in n_pairs]
        )
        for edge_name, midpoint in edges:
            r_label = f"# PML Layer - Type 2: {edge_name}"
            if edge_name in ("Front-Right","Front-Left","Back-Left","Back-Right"):
                for l in range(self.num_layers+1):
                    self.num_regions += 1
                    self.marker_pml_2.append(self.num_regions)
                    PML_regions.append([
                        self.num_regions,
                        round(midpoint[0],5),
                        round(midpoint[1],5),
                        round(midheights[l],5),
                        self.marker_pml_2[-1],
                        None,
                        r_label,
                        None, None, None
                    ])
            else:
                self.num_regions += 1
                self.marker_pml_2.append(self.num_regions)
                PML_regions.append([
                    self.num_regions,
                    round(midpoint[0],5),
                    round(midpoint[1],5),
                    round(midpoint[2],5),
                    self.marker_pml_2[-1],
                    None,
                    r_label,
                    None, None, None
                ])

        # ─── Type 3: corner extensions ─────────────────────────────────────────────
        self.marker_pml_3 = []
        vertices = (
            [("Front-Right-Top",    [self.x_max + i, self.y_min - j, self.z_max + k]) for i,j,k in n_triplets] +
            [("Back-Right-Top",     [self.x_max + i, self.y_max + j, self.z_max + k]) for i,j,k in n_triplets] +
            [("Back-Left-Top",      [self.x_min - i, self.y_max + j, self.z_max + k]) for i,j,k in n_triplets] +
            [("Front-Left-Top",     [self.x_min - i, self.y_min - j, self.z_max + k]) for i,j,k in n_triplets] +
            [("Front-Right-Bottom", [self.x_max + i, self.y_min - j, self.z_min - k]) for i,j,k in n_triplets] +
            [("Back-Right-Bottom",  [self.x_max + i, self.y_max + j, self.z_min - k]) for i,j,k in n_triplets] +
            [("Back-Left-Bottom",   [self.x_min - i, self.y_max + j, self.z_min - k]) for i,j,k in n_triplets] +
            [("Front-Left-Bottom",  [self.x_min - i, self.y_min - j, self.z_min - k]) for i,j,k in n_triplets]
        )
        for vertex_name, midpoint in vertices:
            r_label = f"# PML Layer - Type 3: {vertex_name}"
            self.num_regions += 1
            self.marker_pml_3.append(self.num_regions)
            PML_regions.append([
                self.num_regions,
                round(midpoint[0],5),
                round(midpoint[1],5),
                round(midpoint[2],5),
                self.marker_pml_3[-1],
                None,
                r_label,
                None, None, None
            ])
        print(self.layer_interfaces['z'])
        # ─── Single‐pass material assignment by z-midpoint ────────────────────────
        interfaces = self.layer_interfaces['z']
        num_layers = self.num_layers
        mat_indices = []
        for region in PML_regions:
            mz = region[3]
            # Case 1: above topmost interface (air)
            if mz >= interfaces[0]:
                mat_idx = 0
            # Case 2: below bottommost interface (deepest layer)
            elif mz < interfaces[-1]:
                mat_idx = num_layers
            # Case 3: between interfaces
            else:
                # find first interval where interfaces[i-1] > mz >= interfaces[i]
                mat_idx = None
                for i in range(1, len(interfaces)):
                    if interfaces[i-1] > mz >= interfaces[i]:
                        # layers between interfaces start at index 1
                        mat_idx = i
                        break
                # safety fallback
                if mat_idx is None:
                    raise ValueError(f"Cannot assign material index for z={mz}")
            mat_indices.append(mat_idx)
            # Assign material properties
            region[5] = self.region_volumes[mat_idx]
            region[7] = self.region_rho[mat_idx]
            region[8] = self.region_mu_r[mat_idx]
            region[9] = self.region_epsilon_r[mat_idx]

        return PML_regions

    
    ####################################
    # Writing the .poly file - Assembly
    ####################################
    def write_input_model_file(self):
        """
        Function to write the .poly file using all evaluated nodes, facets and regions."
        """
        
        self.input_model_file_base = os.path.join(self.base_folder, self.input_model_file_base)
        with open(self.input_model_file_base, 'w+') as f:            
            f.write("# Part 1 - Node List\n")
            f.write("{} 3 0 1\n".format(self.num_node))

            f.write("# Points for cube (domain)\n")
            for node in self.node_list_domain_prism:
                f.write("{} {} {} {} {}\n".format(node[0], node[1], node[2], node[3], node[4]))
            f.write("\n")

            node_list_interface = []
            f.write("# Layers\n")
            for layer in range(self.num_layers):
                node_list_interface = self.node_list_of_list_of_interfaces_old[layer]
                if layer == 0:
                    f.write("# interface air-earth\n")
                else:
                    f.write("# interface layer {}\n".format(layer-1))

                for node in node_list_interface:
                    f.write("{} {} {} {} {}\n".format(node[0], node[1], node[2], node[3], node[4]))
                
                f.write("\n")
            f.write("\n")

            f.write("# receivers\n")
            for i, node in enumerate(self.node_list_receiver):
                f.write("{} {} {} {} {}\n".format(node[0], round(node[1], 5), round(node[2], 5), round(node[3], 5), node[4]))
                if (i + 1) % 3 == 0:
                    f.write("\n")
            f.write("\n")

            f.write("# source nodes\n")
            for node in self.node_list_source:
                f.write("{} {} {} {} {}\n".format(node[0], round(node[1],5), round(node[2],5), round(node[3],5), node[4]))
            f.write("\n")

            if self.box_present:
                f.write('# source box nodes\n')
                for node in self.node_list_source_box:
                    f.write("{} {} {} {} {}\n".format(node[0], round(node[1], 5), round(node[2], 5), round(node[3], 5), node[4]))
                f.write("\n")

            if self.NUM_PML > 0:
                f.write("# PML nodes for n PML Layers\n")

                f.write("# Type 2 Interface Intersections: Count Per Lateral Edge: num_interfaces*(n+1)^2 - num_interfaces\n")
                for node in self.node_list_PML_2:
                    f.write("{} {} {} {} {}\n".format(node[0], node[1], node[2], node[3], node[4]))
                f.write("\n")
                f.write("# Type 3: Corner Cubes: Count Per Corner: (n+1)^3 - 1\n")
                for node in self.node_list_PML_3:
                    f.write("{} {} {} {} {}\n".format(node[0], node[1], node[2], node[3], node[4]))
                f.write("\n")
            
            f.write("\n# Part 2 - Facet List\n")
            f.write("{} 1\n".format(self.num_facet))
            f.write(self.domain_string)
            f.write("\n")
            f.write(self.interface_string)
            f.write("\n")
            f.write(self.interfaces_string)
            f.write("\n")
            if any(node[3] > 0.0 for node in self.node_list_source):
                f.write(self.source_string)
            f.write("\n")
            if self.box_present:
                f.write(self.source_box_string)
                f.write("\n")
            if self.NUM_PML > 0:
                f.write(self.pml_string)
                f.write("\n")

            f.write("\n# Part 3 - Hole List\n")
            f.write("0\n")

            f.write("\n# Part 4 - Region List\n")
            f.write("{}\n".format(self.num_regions))
            for region in self.regions:
                f.write("{} {} {} {} {} {} {} \n".format(region[0], region[1], region[2], region[3], region[4], region[5], region[6]))

    #####################################
    # Writing the source coordinates file
    ######################################
    def write_source_file(self):

        self.source_file = os.path.join(self.base_folder, self.source_file)
        with open(self.source_file, 'w+') as f:
            f.write("{}\n".format(2))
            for i in range(2):
                f.write("{} {} {}\n".format(round(self.source_x[i],6), self.source_y[i], self.source_z[i]))

    #########################################
    # Writing the region parameters file
    #########################################
    def write_region_params_file(self):
       
        self.region_params_file = os.path.join(self.base_folder, self.region_params_file)
        with open(self.region_params_file, 'w+') as f:
            f.write("# eleattr\n")
            f.write("{}\n".format(self.num_regions))
            f.write("# eleattr rho mu_r epsilon_r\n")
            for region in self.regions:
                f.write("{} {} {} {}\n".format(region[0], region[7], region[8], region[9]))
    
    #####################################################
    # Writing the mesh sizing .mtr file (unused for now)
    #####################################################
    def write_mesh_sizing_file(self):
        
        self.input_mesh_sizing_file = os.path.join(self.base_folder, self.input_mesh_sizing_file)
        with open(self.input_mesh_sizing_file, 'w+') as f:
            node_list_domain_prism, node_list_of_list_of_interfaces, node_list_receiver, node_list_source, num_node = self.evaluate_all_nodes_with_mesh_size()
            
            f.write("{} 1\n".format(num_node))

            for node in node_list_domain_prism:
                f.write("{} {:.6f}\n".format(node[0], round(node[-1], 7)))

            node_list_interface = []
            for layer in range(self.num_layers):
                node_list_interface = node_list_of_list_of_interfaces[layer]
                for node in node_list_interface:
                    f.write("{} {:.6f}\n".format(node[0], round(node[-1], 7)))
            
            for node in node_list_receiver:
                f.write("{} {:.6f}\n".format(node[0], round(node[-1], 7)))

            for node in node_list_source:
                f.write("{} {:.6f}\n".format(node[0], round(node[-1], 7)))

            for node in self.pml_nodes:
                f.write("{} {:.6f}\n".format(node[0], round(node[-1], 7)))
    
    ###############################################
    # Combined function to evaluate all input data
    ###############################################
    def evaluate_all_input_data(self):
        self.evaluate_all_nodes()
        self.evaluate_all_regions()
        self.evaluate_all_facets()

    ################################################
    # Combined function to write all input files
    ################################################
    def write_all_input_files(self):
        
        # Create the directory if it doesn't exist
        os.makedirs(self.base_folder, exist_ok=True)

        # Write all input files
        self.write_base_input_file()
        self.write_input_model_file()
        self.write_source_file()
        self.write_region_params_file()
        # self.write_mesh_sizing_file()
