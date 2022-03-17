classdef parameters
   % Default input class for multem
   % Contains short descriptions of all parameters and a function to
   % execute the simulation
   
   properties
      system_conf = multem_input.system_conf;                % System Configuration
      
      %%%%%%%%%%%%%% Electron-Specimen interaction model %%%%%%%%%%%%%%%%%
      
      interaction_model(1,1) uint64 {mustBeLessThanOrEqual(interaction_model,3),mustBePositive} = 1;           % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
      potential_type(1,1) uint64 {mustBeLessThanOrEqual(potential_type,6),mustBePositive} = 6;                 % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6
      operation_mode(1,1) uint64 {mustBeLessThanOrEqual(operation_mode,2),mustBePositive} = 1;                 % eOM_Normal = 1, eOM_Advanced = 2
      memory_size(1,1) uint64 {mustBeNonnegative} = 0;                                                         % memory size to be used(Mb)
      reverse_multislice(1,1) uint64 {mustBeLessThanOrEqual(reverse_multislice,1),mustBeNonnegative} = 0;      % 1: true, 0:false
      
      %%%%%%%%%%%%%%% Electron-Phonon interaction model %%%%%%%%%%%%%%%%%%
      
      pn_model(1,1) uint64 {mustBeLessThanOrEqual(pn_model,3),mustBePositive} = 1;                             % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
      pn_coh_contrib(1,1) uint64 {mustBeLessThanOrEqual(pn_coh_contrib,1),mustBeNonnegative} = 0;              % 1: true, 0:false
      pn_single_conf(1,1) uint64 {mustBeLessThanOrEqual(pn_single_conf,1),mustBeNonnegative} = 0;              % 1: true, 0:false (extract single configuration)
      pn_nconf(1,1) uint64 {mustBePositive} = 1;                                                               % true: specific phonon configuration, false: number of frozen phonon configurations
      pn_dim(1,1) uint64 {mustBePositive} = 110;                                                               % phonon dimensions (xyz)
      pn_seed(1,1) uint64 {mustBePositive} =  300183;                                                          % Random seed(frozen phonon)
      %%%%%%%%%%%%%%%%%%%%%%% Specimen information %%%%%%%%%%%%%%%%%%%%%%%
      
      spec_atoms = [];                                                                                         % Specimen atoms
      spec_dz(1,1) double {mustBeNonnegative} = 0.25;                                                          % slice thick (Angstrom)

      spec_lx(1,1) double {mustBeNonnegative} = 10;                                                            % simulation box length in x direction (Angstrom)
      spec_ly(1,1) double {mustBeNonnegative} = 10;                                                            % simulation box length in y direction (Angstrom)
      spec_lz(1,1) double {mustBeNonnegative} = 10;                                                            % simulation box length gpuDEin z direction (Angstrom)

      spec_cryst_na(1,1) uint64 {mustBeNonnegative} = 1;                                                       % number of unit cell along a
      spec_cryst_nb(1,1) uint64 {mustBeNonnegative} = 1;                                                       % number of unit cell along b
      spec_cryst_nc(1,1) uint64 {mustBeNonnegative} = 1;                                                       % number of unit cell along c
      spec_cryst_a(1,1) double {mustBeNonnegative} = 0;                                                        % length along a (Angstrom)
      spec_cryst_b(1,1) double {mustBeNonnegative} = 0;                                                        % length along b (Angstrom)
      spec_cryst_c(1,1) double {mustBeNonnegative} = 0;                                                        % length along c (Angstrom)
      spec_cryst_x0(1,1) double = 0;                                                                           % reference position along x direction (Angstrom)
      spec_cryst_y0(1,1) double = 0;                                                                           % reference position along y direction (Angstrom)

      spec_amorp = struct('z_0', 0, 'z_e', 0, 'dz', 2.0);                                                      % Include slice of amorphous layer
      %%%%%%%%%%%%%%%%%%%%%%%%% Specimen Rotation %%%%%%%%%%%%%%%%%%%%%%%%
      
      spec_rot_theta(1,1) double {mustBeNonnegative} = 0;                                                      % angle (Degree)
      spec_rot_u0(3,1) double {mustBeNonnegative} = [0 0 1];                                                   % unitary vector			
      spec_rot_center_type(1,1) uint64 {mustBeLessThanOrEqual(spec_rot_center_type,2),mustBePositive} = 1;     % 1: geometric center, 2: User define		
      spec_rot_center_p(3,1) double {mustBeNonnegative} = [0 0 0];                                             % rotation point
      %%%%%%%%%%%%%%%%%%%%%% Specimen thickness %%%%%%%%%%%%%%%%%%%%%%%%%%
      
      thick_type(1,1) uint64 {mustBeLessThanOrEqual(thick_type,3),mustBePositive} = 1;                         % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
      thick = 0;                                                                                               % Array of thickness (Angstrom)
      %%%%%%%%%%%%%%%%%%%%%%% Potential slicing %%%%%%%%%%%%%%%%%%%%%%%%%%
      
      potential_slicing(1,1) uint64 {mustBeLessThanOrEqual(potential_slicing,4),mustBePositive} = 1;          % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
      islice uint64                                                                                           % retrieve the projected potential at given slices
      %%%%%%%%%%%%%%%%%%%%%% x-y sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      nx(1,1) double {mustBePositive} = 256;                                                                  % number of pixels in x direction
      ny(1,1) double {mustBePositive} = 256;                                                                  % number of pixels in y direction
      bwl(1,1) uint64 {mustBeLessThanOrEqual(bwl,1),mustBeNonnegative} = 0;                                   % Band-width limit, 1: true, 0:false

      %%%%%%%%%%%%%%%%%%%%%%%% Simulation type %%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Simulation type
      % eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, 
      % eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, 
      % eTEMST_EWFS=51, eTEMST_EWRS=52, eTEMST_EELS=61, eTEMST_EFTEM=62, 
      % eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82, 
      % eTEMST_TFFS=91, eTEMST_TFRS=92
       simulation_type(1,1) uint64 = 52;                                                                       % Simulation type
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
      
      iw_type(1,1) uint64 {mustBeLessThanOrEqual(iw_type,4),mustBePositive} = 4;                               % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
      iw_psi double = 0;                                                                                       % User defined incident wave. It will be only used if iw_type=3
      iw_x double = 0;                                                                                         % x position 
      iw_y double = 0;                                                                                         % y position
      %%%%%%%%%%%%%%%%%%%% Microscope parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
      
      E_0(1,1) double = 300;                                                                                   % Acceleration Voltage (keV)
      theta(1,1) double = 0;                                                                                   % Polar angle (Angstrom)
      phi(1,1) double = 0;                                                                                     % Azimuthal angle (Angstrom)
      %%%%%%%%%%%%%%%%%%%%%% Illumination model %%%%%%%%%%%%%%%%%%%%%%%%%%
      
      illumination_model(1,1) uint64 {mustBeLessThanOrEqual(illumination_model,4),mustBePositive} = 3; 			% 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient, 4: Numerical integration
      temporal_spatial_incoh(1,1) uint64 {mustBeLessThanOrEqual(temporal_spatial_incoh,3),mustBePositive} = 1; 	% 1: Temporal and Spatial, 2: Temporal, 3: Spatial
      %%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
      
      cond_lens_m(1,1) double = 0.00;                       % Vortex momentum

      cond_lens_c_10(1,1) double = 14.0312; 				% [C1]      Defocus (Angstrom)
      cond_lens_c_12(1,1) double = 0.00; 					% [A1]      2-fold astigmatism (Angstrom)
      cond_lens_phi_12(1,1) double = 0.00; 					% [phi_A1]	Azimuthal angle of 2-fold astigmatism (Angstrom)

      cond_lens_c_21(1,1) double = 0.00; 					% [B2]      Axial coma (Angstrom)
      cond_lens_phi_21(1,1) double = 0.00; 					% [phi_B2]	Azimuthal angle of axial coma (Angstrom)
      cond_lens_c_23(1,1) double = 0.00;					% [A2]      3-fold astigmatism (Angstrom)
      cond_lens_phi_23(1,1) double = 0.00; 					% [phi_A2]	Azimuthal angle of 3-fold astigmatism (Angstrom)

      cond_lens_c_30(1,1) double = 1e-3; 					% [C3] 		3rd order spherical aberration (mm)
      cond_lens_c_32(1,1) double = 0.00; 					% [S3]      Axial star aberration (Angstrom)
      cond_lens_phi_32(1,1) double = 0.00; 					% [phi_S3]	Azimuthal angle of axial star aberration (Angstrom)
      cond_lens_c_34(1,1) double = 0.00; 					% [A3]      4-fold astigmatism (Angstrom)
      cond_lens_phi_34(1,1) double = 0.00; 					% [phi_A3]	Azimuthal angle of 4-fold astigmatism (Angstrom)

      cond_lens_c_41(1,1) double = 0.00; 					% [B4]      4th order axial coma (Angstrom)
      cond_lens_phi_41(1,1) double = 0.00; 					% [phi_B4]	Azimuthal angle of 4th order axial coma (Angstrom)
      cond_lens_c_43(1,1) double = 0.00; 					% [D4]      3-lobe aberration (Angstrom)
      cond_lens_phi_43(1,1) double = 0.00; 					% [phi_D4]	Azimuthal angle of 3-lobe aberration (Angstrom)
      cond_lens_c_45(1,1) double = 0.00; 					% [A4]      5-fold astigmatism (Angstrom)
      cond_lens_phi_45(1,1) double = 0.00; 					% [phi_A4]	Azimuthal angle of 5-fold astigmatism (Angstrom)

      cond_lens_c_50(1,1) double = 0.00; 					% [C5]      5th order spherical aberration (mm)
      cond_lens_c_52(1,1) double = 0.00; 					% [S5]      5th order axial star aberration (Angstrom)
      cond_lens_phi_52(1,1) double = 0.00; 					% [phi_S5]	Azimuthal angle of 5th order axial star aberration (Angstrom)
      cond_lens_c_54(1,1) double = 0.00; 					% [R5]      5th order rosette aberration (Angstrom)
      cond_lens_phi_54(1,1) double = 0.00; 					% [phi_R5]	Azimuthal angle of 5th order rosette aberration (Angstrom)
      cond_lens_c_56(1,1) double = 0.00; 					% [A5]      6-fold astigmatism (Angstrom)
      cond_lens_phi_56(1,1) double = 0.00; 					% [phi_A5]	Azimuthal angle of 6-fold astigmatism (Angstrom)

      cond_lens_inner_aper_ang(1,1) double = 0.00; 			% Inner aperture (mrad) 
      cond_lens_outer_aper_ang(1,1) double = 21.0; 			% Outer aperture (mrad)
      %%%%%%%%%% source spread function %%%%%%%%%%%%
      
      cond_lens_si_a(1,1) double = 1.0                      % Height proportion of a normalized Gaussian [0, 1]
      cond_lens_si_sigma(1,1) double = 0.0072; 				% standard deviation: For parallel ilumination(Angstrom^-1); otherwise (Angstrom)
      cond_lens_si_beta(1,1) double = 0.0;                 	% Standard deviation of the source spread function for the Exponential component: For parallel ilumination(�^-1); otherwise (�)
      cond_lens_si_rad_npts(1,1) uint64 = 4;                % # of integration points. It will be only used if illumination_model=4
      cond_lens_si_azm_npts(1,1) uint64 = 4;                % # of radial integration points. It will be only used if illumination_model=4
      %%%%%%%%% defocus spread function %%%%%%%%%%%%
      
      cond_lens_ti_a(1,1) double = 1.0;                     % Height proportion of a normalized Gaussian [0, 1]
      cond_lens_ti_sigma(1,1) double = 32.0; 				% standard deviation (Angstrom)
      cond_lens_ti_beta(1,1) double = 0.0;                 	% Standard deviation of the defocus spread for the Exponential component
      cond_lens_ti_npts(1,1) uint64 = 10;  			        % # of integration points. It will be only used if illumination_model=4
      %%%%%%%%% zero defocus reference %%%%%%%%%%%%
      
      cond_lens_zero_defocus_type(1,1) uint64 {mustBeLessThanOrEqual(cond_lens_zero_defocus_type,4),mustBePositive} = 1; 		% eZDT_First = 1, eZDT_User_Define = 4
      cond_lens_zero_defocus_plane(1,1) double = 0.00;  	% It will be only used if cond_lens_zero_defocus_type = eZDT_User_Define
      %%%%%%%%%%%%%%%%%%%%%%%% Objective lens %%%%%%%%%%%%%%%%%%%%%%%%
      
      obj_lens_m(1,1) double = 0;                           % Vortex momentum

      obj_lens_c_10(1,1) double = 14.0312;  				% [C1]      Defocus (Angstrom)
      obj_lens_c_12(1,1) double = 0.00;  					% [A1]      2-fold astigmatism (Angstrom)
      obj_lens_phi_12(1,1) double = 0.00;  					% [phi_A1]	Azimuthal angle of 2-fold astigmatism (Angstrom)

      obj_lens_c_21(1,1) double = 0.00;  					% [B2]      Axial coma (Angstrom)
      obj_lens_phi_21(1,1) double = 0.00;  					% [phi_B2]	Azimuthal angle of axial coma (Angstrom)
      obj_lens_c_23(1,1) double = 0.00;  					% [A2]      3-fold astigmatism (Angstrom)
      obj_lens_phi_23(1,1) double = 0.00;  					% [phi_A2]	Azimuthal angle of 3-fold astigmatism (Angstrom)

      obj_lens_c_30(1,1) double = 1e-03;  					% [C3] 		3rd order spherical aberration (mm)
      obj_lens_c_32(1,1) double = 0.00;  					% [S3]      Axial star aberration (Angstrom)
      obj_lens_phi_32(1,1) double = 0.00;  					% [phi_S3]	Azimuthal angle of axial star aberration (Angstrom)
      obj_lens_c_34(1,1) double = 0.00;  					% [A3]      4-fold astigmatism (Angstrom)
      obj_lens_phi_34(1,1) double = 0.00;  					% [phi_A3]	Azimuthal angle of 4-fold astigmatism (Angstrom)

      obj_lens_c_41(1,1) double = 0.00;  					% [B4]      4th order axial coma (Angstrom)
      obj_lens_phi_41(1,1) double = 0.00;  					% [phi_B4]	Azimuthal angle of 4th order axial coma (Angstrom)
      obj_lens_c_43(1,1) double = 0.00;  					% [D4]      3-lobe aberration (Angstrom)
      obj_lens_phi_43(1,1) double = 0.00;  					% [phi_D4]	Azimuthal angle of 3-lobe aberration (Angstrom)
      obj_lens_c_45(1,1) double = 0.00;  					% [A4]      5-fold astigmatism (Angstrom)
      obj_lens_phi_45(1,1) double = 0.00;  					% [phi_A4]	Azimuthal angle of 5-fold astigmatism (Angstrom)

      obj_lens_c_50(1,1) double = 0.00;  					% [C5]      5th order spherical aberration (mm)
      obj_lens_c_52(1,1) double = 0.00;  					% [S5]      5th order axial star aberration (Angstrom)
      obj_lens_phi_52(1,1) double = 0.00;  					% [phi_S5]	Azimuthal angle of 5th order axial star aberration (Angstrom)
      obj_lens_c_54(1,1) double = 0.00;  					% [R5]      5th order rosette aberration (Angstrom)
      obj_lens_phi_54(1,1) double = 0.00;  					% [phi_R5]	Azimuthal angle of 5th order rosette aberration (Angstrom)
      obj_lens_c_56(1,1) double = 0.00;  					% [A5]      6-fold astigmatism (Angstrom)
      obj_lens_phi_56(1,1) double = 0.00;  					% [phi_A5]	Azimuthal angle of 6-fold astigmatism (Angstrom)

      obj_lens_inner_aper_ang(1,1) double = 0.00;  			% Inner aperture (mrad) 
      obj_lens_outer_aper_ang(1,1) double = 24.0;  			% Outer aperture (mrad)         
     %%%%%%%%%% source spread function %%%%%%%%%%%%
      
      obj_lens_si_a(1,1) double = 1.0                       % Height proportion of a normalized Gaussian [0, 1]
      obj_lens_si_sigma(1,1) double = 0.0072; 			    % standard deviation: For parallel ilumination(Angstrom^-1); otherwise (Angstrom)
      obj_lens_si_beta(1,1) double = 0.0;                 	% Standard deviation of the source spread function for the Exponential component: For parallel ilumination(�^-1); otherwise (�)
      obj_lens_si_rad_npts(1,1) uint64 = 4;                 % # of integration points. It will be only used if illumination_model=4
      obj_lens_si_azm_npts(1,1) uint64 = 4;                 % # of radial integration points. It will be only used if illumination_model=4
      %%%%%%%%% defocus spread function %%%%%%%%%%%%
      
      obj_lens_ti_a(1,1) double = 1.0;                      % Height proportion of a normalized Gaussian [0, 1]
      obj_lens_ti_sigma(1,1) double = 32.0; 				% standard deviation (Angstrom)
      obj_lens_ti_beta(1,1) double = 0.0;                 	% Standard deviation of the defocus spread for the Exponential component
      obj_lens_ti_npts(1,1) uint64 = 10;  			        % # of integration points. It will be only used if illumination_model=4

      %%%%%%%%% zero defocus reference %%%%%%%%%%%%
      
      obj_lens_zero_defocus_type(1,1) uint64 {mustBeLessThanOrEqual(obj_lens_zero_defocus_type,4),mustBePositive} = 1; 		% eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
      obj_lens_zero_defocus_plane(1,1) double = 0.00;                                                                       % It will be only used if obj_lens_zero_defocus_type = eZDT_User_Define
      % %%%%%%%%%%%%%%%%%%%%%%%%%% STEM Detector %%%%%%%%%%%%%%%%%%%%%%%%%%
      
      detector = multem_input.detector;                                                                                     % detector object
      %%%%%%%%%%%% Scanning area for ISTEM/STEM/EELS %%%%%%%%%%%%%%%%%

      scanning_type(1,1) uint64 {mustBeLessThanOrEqual(scanning_type,2),mustBePositive} = 2;                                % eST_Line = 1, eST_Area = 2
      scanning_periodic(1,1) uint64 {mustBeLessThanOrEqual(scanning_periodic,1),mustBeNonnegative} = 0;                     % 1: true, 0:false (periodic boundary conditions)
      scanning_square_pxs(1,1) uint64 {mustBeLessThanOrEqual(scanning_square_pxs,1),mustBeNonnegative} = 1; 				% 0: false, 1: true
      scanning_ns(1,1) uint64 {mustBePositive} = 10;        % number of sampling points
      scanning_x0(1,1) double = 0.00;                       % x-starting point (Angstrom)
      scanning_y0(1,1) double = 0.00;                       % y-starting point (Angstrom)
      scanning_xe(1,1) double = 4.078;                      % x-final point (Angstrom)
      scanning_ye(1,1) double = 4.078;                      % y-final point (Angstrom)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      ped_nrot(1,1) double = 360;                           % Number of orientations
      ped_theta(1,1) double = 3.0;                          % Precession angle (degrees)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HCI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      hci_nrot(1,1) double = 360;                           % number of orientations
      hci_theta(1,1) double = 3.0;                          % Precession angle (degrees)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      eels_Z(1,1) uint64 = 79;                              % atomic type
      eels_E_loss(1,1) double = 80.0;                       % Energy loss (eV)
      eels_collection_angle(1,1) double = 100.0;   			% Collection half angle (mrad)
      eels_m_selection(1,1) double = 3;   					% selection rule
      eels_channelling_type(1,1) uint64 {mustBeLessThanOrEqual(eels_channelling_type,3),mustBePositive} = 1; 			% eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EFTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      eftem_Z(1,1) uint64 = 79;                             % atomic type
      eftem_E_loss(1,1) double = 80.0;                      % Energy loss (eV)
      eftem_collection_angle(1,1) double = 100.0;  			% Collection half angle (mrad)
      eftem_m_selection(1,1) double = 3;                    % selection rule
      eftem_channelling_type(1,1) uint64 {mustBeLessThanOrEqual(eftem_channelling_type,3),mustBePositive} = 1; 			% eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 
      
      %%%%%%%%%%%%%%%%%%%%%%% OUTPUT REGION %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% This option is not used for eTEMST_STEM and eTEMST_EELS %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      output_area_ix_0(1,1) double = 0.0; 					% x-starting in pixel
      output_area_iy_0(1,1) double = 0.0; 					% y-starting in pixel
      output_area_ix_e(1,1) double = 0.0; 					% x-final in pixel
      output_area_iy_e(1,1) double = 0.0; 					% y-final in pixel
    end
    methods
        function out_mt = ilc_incident_wave(obj)
            prms = obj.toStruct;
            clear ilc_incident_wave;
            out_mt = ilc_incident_wave(prms.system_conf, prms);
        end
        function out_mt = ilc_multem(obj)
            prms = obj.toStruct;
            clear ilc_multem;
            out_mt = ilc_multem(prms.system_conf, prms);
        end
        function out_mt = ilc_projected_potential(obj)
            prms = obj.toStruct;
            clear ilc_projected_potential;
            out_mt = ilc_projected_potential(prms.system_conf, prms);
        end
        function out_mt = ilc_propagate(obj)
            prms = obj.toStruct;
            clear ilc_propagate;
            out_mt = ilc_propagate(prms.system_conf, prms);
        end
        function out_mt = ilc_microscope_aberrations(obj)
            prms = obj.toStruct;
            clear ilc_microscope_aberrations;
            out_mt = ilc_microscope_aberrations(prms.system_conf, prms);
        end
        function out_mt = ilc_apply_ctf(obj)
            prms = obj.toStruct;
            clear ilc_propagate;
            out_mt = ilc_apply_ctf(prms.system_conf, prms);
        end
        function out_mt = ilc_transmission_function(obj)
            prms = obj.toStruct;
            clear ilc_transmission_function;
            out_mt = ilc_transmission_function(prms.system_conf, prms);
        end
        function out_mt = ilc_wave_function(obj)
            prms = obj.toStruct;
            clear ilc_wave_function;
            out_mt = ilc_wave_function(prms.system_conf, prms);
        end
        function s = toStruct(class)
            public_props = properties(class);
            s = struct();
            for fi = 1:numel(public_props)
                if isobject(class.(public_props{fi}))
                   sub_props = properties(class.(public_props{fi}));
                   for fi_s = 1:numel(sub_props)
                       val = class.(public_props{fi}).(sub_props{fi_s});
                       if isnumeric(val)
                           val = double(val);
                       end
                        s.(public_props{fi}).(sub_props{fi_s}) = val;   
                   end
                else
                    val = class.(public_props{fi});
                    if isnumeric(val)
                        val = double(val);
                    end
                    s.(public_props{fi}) = val; 
                end            
            end                  
        end 
    end
end
