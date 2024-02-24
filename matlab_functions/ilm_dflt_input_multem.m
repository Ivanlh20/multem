function [input_multem] = ilm_dflt_input_multem()
    %%%%%%%%%%%%%% Electron-Specimen interaction model %%%%%%%%%%%%%%%%%
    input_multem.interaction_model = 1;                         % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
    input_multem.potential_type = 6;                            % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.operation_mode = 1;                            % eOM_Normal = 1, eOM_Advanced = 2
    input_multem.memory_size = 0;                               % memory size to be used(Mb)
    input_multem.reverse_multislice = 0;                        % 1: true, 0:false

    %%%%%%%%%%%%%%% Electron-Phonon interaction model %%%%%%%%%%%%%%%%%%
    input_multem.pn_model = 1;                                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
    input_multem.pn_coh_contrib = 0;                            % 1: true, 0:false
    input_multem.pn_single_conf = 0;                            % 1: true, 0:false (extract single configuration)
    input_multem.pn_nconf = 1;                                  % true: specific phonon configuration, false: number of frozen phonon configurations
    input_multem.pn_dim = 110;                                  % phonon dimensions (xyz)
    input_multem.pn_seed = 300183;                              % Random seed(frozen phonon)

    %%%%%%%%%%%%%%%%%%%%%%% Specimen information %%%%%%%%%%%%%%%%%%%%%%%
    input_multem.spec_atoms = [];                               % simulation box length in x direction (Å)
    input_multem.spec_dz = 0.25;                                % slice thick (Å)

    input_multem.spec_lx = 10;                                  % simulation box length in x direction (Å)
    input_multem.spec_ly = 10;                                  % simulation box length in y direction (Å)
    input_multem.spec_lz = 10;                                  % simulation box length gpuDEin z direction (Å)

    input_multem.spec_cryst_na = 1;                             % number of unit cell along a
    input_multem.spec_cryst_nb = 1;                             % number of unit cell along b
    input_multem.spec_cryst_nc = 1;                             % number of unit cell along c
    input_multem.spec_cryst_a = 0;                              % length along a (Å)
    input_multem.spec_cryst_b = 0;                              % length along b (Å)
    input_multem.spec_cryst_c = 0;                              % length along c (Å)
    input_multem.spec_cryst_x0 = 0;                             % reference position along x direction (Å)
    input_multem.spec_cryst_y0 = 0;                             % reference position along y direction (Å)

    input_multem.spec_amorp(1).z_0 = 0;                         % Starting z position of the amorphous layer (Å)
    input_multem.spec_amorp(1).z_e = 0;                         % Ending z position of the amorphous layer (Å)
    input_multem.spec_amorp(1).dz = 2.0;                        % slice thick of the amorphous layer (Å)

    %%%%%%%%%%%%%%%%%%%%%%%%% Specimen Rotation %%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.spec_rot_theta = 0;                            % angle (º)
    input_multem.spec_rot_u0 = [0 0 1];                          % unitary vector			
    input_multem.spec_rot_center_type = 1;                       % 1: geometric center, 2: User define		
    input_multem.spec_rot_center_p = [0 0 0];                    % rotation point

    %%%%%%%%%%%%%%%%%%%%%% Specimen thickness %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.thick_type = 1;                                % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
    input_multem.thick = 0;                                     % Array of thickness (Å)

    %%%%%%%%%%%%%%%%%%%%%%% Potential slicing %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.potential_slicing = 1;                         % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4

    %%%%%%%%%%%%%%%%%%%%%% x-y sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.nx = 256;                                      % number of pixels in x direction
    input_multem.ny = 256;                                      % number of pixels in y direction
    input_multem.bwl = 0;                                       % Band-width limit, 1: true, 0:false

    %%%%%%%%%%%%%%%%%%%%%%%% Simulation type %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, 
    % eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, 
    % eTEMST_EWFS=51, eTEMST_EWRS=52, eTEMST_EELS=61, eTEMST_EFTEM=62, 
    % eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82, 
    % eTEMST_TFFS=91, eTEMST_TFRS=92
    input_multem.simulation_type = 52;    

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.iw_type = 4;                                   % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
    input_multem.iw_psi = 0;                                    % User define incident wave. It will be only used if User_Define=3
    input_multem.iw_x = 0.0;                                    % x position 
    input_multem.iw_y = 0.0;                                    % y position

    %%%%%%%%%%%%%%%%%%%% Microscope parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.E_0 = 300;                                     % Acceleration Voltage (keV)
    input_multem.theta = 0.0;                                   % Polar angle (°)
    input_multem.phi = 0.0;                                     % Azimuthal angle (°)

    %%%%%%%%%%%%%%%%%%%%%% Illumination model %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.illumination_model = 2;                        % 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient, 4: Numerical integration
    input_multem.temporal_spatial_incoh = 1;                    % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % On the optimun probe in aberration corrected ADF-STEM
    % Ultramicroscopy 111(2014) 1523-1530
    % C_{nm} Krivanek --- {A, B, C, D, R}_{n} Haider notation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.cond_lens_m = 0;                               % Vortex momentum

    input_multem.cond_lens_c_10 = 14.0312;                      % [C1]      Defocus (Å)
    input_multem.cond_lens_c_12 = 0.00;                         % [A1]      2-fold astigmatism (Å)
    input_multem.cond_lens_phi_12 = 0.00;                       % [phi_A1]	Azimuthal angle of 2-fold astigmatism (°)

    input_multem.cond_lens_c_21 = 0.00;                         % [B2]      Axial coma (Å)
    input_multem.cond_lens_phi_21 = 0.00;                       % [phi_B2]	Azimuthal angle of axial coma (°)
    input_multem.cond_lens_c_23 = 0.00;                         % [A2]      3-fold astigmatism (Å)
    input_multem.cond_lens_phi_23 = 0.00;                       % [phi_A2]	Azimuthal angle of 3-fold astigmatism (°)

    input_multem.cond_lens_c_30 = 1e-03;                        % [C3] 		3rd order spherical aberration (mm)
    input_multem.cond_lens_c_32 = 0.00;                         % [S3]      Axial star aberration (Å)
    input_multem.cond_lens_phi_32 = 0.00;                       % [phi_S3]	Azimuthal angle of axial star aberration (°)
    input_multem.cond_lens_c_34 = 0.00;                         % [A3]      4-fold astigmatism (Å)
    input_multem.cond_lens_phi_34 = 0.0;                        % [phi_A3]	Azimuthal angle of 4-fold astigmatism (°)

    input_multem.cond_lens_c_41 = 0.00;                         % [B4]      4th order axial coma (Å)
    input_multem.cond_lens_phi_41 = 0.00;                       % [phi_B4]	Azimuthal angle of 4th order axial coma (°)
    input_multem.cond_lens_c_43 = 0.00;                         % [D4]      3-lobe aberration (Å)
    input_multem.cond_lens_phi_43 = 0.00;                       % [phi_D4]	Azimuthal angle of 3-lobe aberration (°)
    input_multem.cond_lens_c_45 = 0.00;                         % [A4]      5-fold astigmatism (Å)
    input_multem.cond_lens_phi_45 = 0.00;                       % [phi_A4]	Azimuthal angle of 5-fold astigmatism (°)

    input_multem.cond_lens_c_50 = 0.00;                         % [C5]      5th order spherical aberration (mm)
    input_multem.cond_lens_c_52 = 0.00;                         % [S5]      5th order axial star aberration (Å)
    input_multem.cond_lens_phi_52 = 0.00;                       % [phi_S5]	Azimuthal angle of 5th order axial star aberration (°)
    input_multem.cond_lens_c_54 = 0.00;                         % [R5]      5th order rosette aberration (Å)
    input_multem.cond_lens_phi_54 = 0.00;                       % [phi_R5]	Azimuthal angle of 5th order rosette aberration (°)
    input_multem.cond_lens_c_56 = 0.00;                         % [A5]      6-fold astigmatism (Å)
    input_multem.cond_lens_phi_56 = 0.00;                       % [phi_A5]	Azimuthal angle of 6-fold astigmatism (°)

    input_multem.cond_lens_inner_aper_ang = 0.0;                % Inner aperture (mrad) 
    input_multem.cond_lens_outer_aper_ang = 21.0;   			% Outer aperture (mrad)

    %%%%%%%%% defocus spread function %%%%%%%%%%%%
    dsf_sigma = ilc_iehwgd_2_sigma(32);                         % from defocus spread to standard deviation
    input_multem.cond_lens_ti_a = 1.0;
    input_multem.cond_lens_ti_sigma = dsf_sigma;                % standard deviation
    input_multem.cond_lens_ti_beta = 1;
    input_multem.cond_lens_ti_npts = 5;                         % # of integration points. It will be only used if illumination_model=4

    %%%%%%%%%% source spread function %%%%%%%%%%%%
    ssf_sigma = ilc_hwhm_2_sigma(0.45);                         % half width at half maximum to standard deviation
    input_multem.cond_lens_si_a = 1.0;
    input_multem.cond_lens_si_sigma = ssf_sigma;                % standard deviation: For parallel ilumination(Å^-1); otherwise
    input_multem.cond_lens_si_beta = 1.0;
    input_multem.cond_lens_si_rad_npts = 8;                     % # of integration points. It will be only used if illumination_model=4
    input_multem.cond_lens_si_azm_npts = 12;                    % # of integration points. It will be only used if illumination_model=4

    %%%%%%%%% zero defocus reference %%%%%%%%%%%%
    input_multem.cond_lens_zero_defocus_type = 1;   			% eZDT_First = 1, eZDT_User_Define = 4
    input_multem.cond_lens_zero_defocus_plane = 0;  			% It will be only used if cond_lens_zero_defocus_type = eZDT_User_Define

    %%%%%%%%%%%%%%%%%%% condenser lens variable %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% it will be active in the future %%%%%%%%%%%%%%%%%
    % 1: Vortex momentum, 2: Defocus (Å), 2: Third order spherical aberration (mm)
    % 3: Third order spherical aberration (mm),  4: Fifth order spherical aberration (mm)
    % 5: Twofold astigmatism (Å), 2: Defocus (Å), 6: Azimuthal angle of the twofold astigmatism (°)
    % 7: Threefold astigmatism (Å),  8: Azimuthal angle of the threefold astigmatism (°)
    % 9: Inner aperture (mrad), 2: Defocus (Å), 10: Outer aperture (mrad)
    % input_multem.cdl_var_type = 0;                  			% 0:off 1: m, 2: f, 3 Cs3, 4:Cs5, 5:mfa2, 6:afa2, 7:mfa3, 8:afa3, 9:inner_aper_ang , 10:outer_aper_ang
    % input_multem.cdl_var = [-2 -1 0 1 2];           			% variable array

    %%%%%%%%%%%%%%%%%%%%%%%% Objective lens %%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.obj_lens_m = 0;                                % Vortex momentum

    input_multem.obj_lens_c_10 = 14.0312;                       % [C1]      Defocus (Å)
    input_multem.obj_lens_c_12 = 0.00;                          % [A1]      2-fold astigmatism (Å)
    input_multem.obj_lens_phi_12 = 0.00;                        % [phi_A1]	Azimuthal angle of 2-fold astigmatism (°)

    input_multem.obj_lens_c_21 = 0.00;                          % [B2]      Axial coma (Å)
    input_multem.obj_lens_phi_21 = 0.00;                        % [phi_B2]	Azimuthal angle of axial coma (°)
    input_multem.obj_lens_c_23 = 0.00;                          % [A2]      3-fold astigmatism (Å)
    input_multem.obj_lens_phi_23 = 0.00;                        % [phi_A2]	Azimuthal angle of 3-fold astigmatism (°)

    input_multem.obj_lens_c_30 = 1e-03;                         % [C3] 		3rd order spherical aberration (mm)
    input_multem.obj_lens_c_32 = 0.00;                          % [S3]      Axial star aberration (Å)
    input_multem.obj_lens_phi_32 = 0.00;                        % [phi_S3]	Azimuthal angle of axial star aberration (°)
    input_multem.obj_lens_c_34 = 0.00;                          % [A3]      4-fold astigmatism (Å)
    input_multem.obj_lens_phi_34 = 0.0;                         % [phi_A3]	Azimuthal angle of 4-fold astigmatism (°)

    input_multem.obj_lens_c_41 = 0.00;                          % [B4]      4th order axial coma (Å)
    input_multem.obj_lens_phi_41 = 0.00;                        % [phi_B4]	Azimuthal angle of 4th order axial coma (°)
    input_multem.obj_lens_c_43 = 0.00;                          % [D4]      3-lobe aberration (Å)
    input_multem.obj_lens_phi_43 = 0.00;                        % [phi_D4]	Azimuthal angle of 3-lobe aberration (°)
    input_multem.obj_lens_c_45 = 0.00;                          % [A4]      5-fold astigmatism (Å)
    input_multem.obj_lens_phi_45 = 0.00;                        % [phi_A4]	Azimuthal angle of 5-fold astigmatism (°)

    input_multem.obj_lens_c_50 = 0.00;                          % [C5]      5th order spherical aberration (mm)
    input_multem.obj_lens_c_52 = 0.00;                          % [S5]      5th order axial star aberration (Å)
    input_multem.obj_lens_phi_52 = 0.00;                        % [phi_S5]	Azimuthal angle of 5th order axial star aberration (°)
    input_multem.obj_lens_c_54 = 0.00;                          % [R5]      5th order rosette aberration (Å)
    input_multem.obj_lens_phi_54 = 0.00;                        % [phi_R5]	Azimuthal angle of 5th order rosette aberration (°)
    input_multem.obj_lens_c_56 = 0.00;                          % [A5]      6-fold astigmatism (Å)
    input_multem.obj_lens_phi_56 = 0.00;                        % [phi_A5]	Azimuthal angle of 6-fold astigmatism (°)

    input_multem.obj_lens_inner_aper_ang = 0.0;     			% Inner aperture (mrad) 
    input_multem.obj_lens_outer_aper_ang = 21.4;    			% Outer aperture (mrad)

    %%%%%%%%% defocus spread function %%%%%%%%%%%%
    dsf_sigma = ilc_iehwgd_2_sigma(32);                         % from defocus spread to standard deviation
    input_multem.obj_lens_ti_a = 1.0;
    input_multem.obj_lens_ti_sigma = dsf_sigma;                % standard deviation
    input_multem.obj_lens_ti_beta = 1;
    input_multem.obj_lens_ti_npts = 5;                         % # of integration points. It will be only used if illumination_model=4

    %%%%%%%%%% source spread function %%%%%%%%%%%%
    ssf_sigma = ilc_hwhm_2_sigma(0.45);                         % half width at half maximum to standard deviation
    input_multem.obj_lens_si_a = 1.0;
    input_multem.obj_lens_si_sigma = ssf_sigma;                % standard deviation: For parallel ilumination(Å^-1); otherwise
    input_multem.obj_lens_si_beta = 1.0;
    input_multem.obj_lens_si_rad_npts = 8;                     % # of integration points. It will be only used if illumination_model=4
    input_multem.obj_lens_si_azm_npts = 12;                    % # of integration points. It will be only used if illumination_model=4
    
    %%%%%%%%% zero defocus reference %%%%%%%%%%%%
    input_multem.obj_lens_zero_defocus_type = 3;    			% eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
    input_multem.obj_lens_zero_defocus_plane = 0;   			% It will be only used if obj_lens_zero_defocus_type = eZDT_User_Define

    %%%%%%%%%%%%%%%%%%%%%%%%%% STEM Detector %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.detector.type = 1;  % eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3

    input_multem.detector.cir(1).inner_ang = 60;    			% Inner angle(mrad) 
    input_multem.detector.cir(1).outer_ang = 180;   			% Outer angle(mrad)

    input_multem.detector.radial(1).x = 0;          			% radial detector angle(mrad)
    input_multem.detector.radial(1).fx = 0;         			% radial sensitivity value

    input_multem.detector.matrix(1).R = 0;          			% 2D detector angle(mrad)
    input_multem.detector.matrix(1).fR = 0;         			% 2D sensitivity value

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Scanning area for ISTEM/STEM/EELS %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.scanning_type = 1;                 			% eST_Line = 1, eST_Area = 2
    input_multem.scanning_periodic = 1;             			% 1: true, 0:false (periodic boundary conditions)
    input_multem.scanning_square_pxs = 1;             			% 0: false, 1: true
    input_multem.scanning_ns = 10;                  			% number of sampling points
    input_multem.scanning_x0 = 0.0;                 			% x-starting point (Å)
    input_multem.scanning_y0 = 0.0;                 			% y-starting point (Å)
    input_multem.scanning_xe = 4.078;               			% x-final point (Å)
    input_multem.scanning_ye = 4.078;               			% y-final point (Å)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.ped_nrot = 360;                    			% Number of orientations
    input_multem.ped_theta = 3.0;                   			% Precession angle (degrees)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HCI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.hci_nrot = 360;                    			% number of orientations
    input_multem.hci_theta = 3.0;                   			% Precession angle (degrees)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.eels_Z = 79;                       			% atomic type
    input_multem.eels_E_loss = 80;                  			% Energy loss (eV)
    input_multem.eels_collection_angle = 100;       			% Collection half angle (mrad)
    input_multem.eels_m_selection = 3;              			% selection rule
    input_multem.eels_channelling_type = 1;         			% eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EFTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.eftem_Z = 79;                      			% atomic type
    input_multem.eftem_E_loss = 80;                 			% Energy loss (eV)
    input_multem.eftem_collection_angle = 100;      			% Collection half angle (mrad)
    input_multem.eftem_m_selection = 3;             			% selection rule
    input_multem.eftem_channelling_type = 1;        			% eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 

    %%%%%%%%%%%%%%%%%%%%%%% OUTPUT REGION %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% This option is not use for eTEMST_STEM and eTEMST_EELS %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.output_area_ix_0 = 1;                             % x-starting in pixel
    input_multem.output_area_iy_0 = 1;                             % y-starting in pixel
    input_multem.output_area_ix_e = 1;                             % x-final in pixel
    input_multem.output_area_iy_e = 1;                             % y-final in pixel
end
