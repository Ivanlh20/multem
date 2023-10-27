function [input_multem] = ilm_dflt_input_multem()
    %%%%%%%%%%%%%% Electron-Specimen interaction model %%%%%%%%%%%%%%%%%
    input_multem.elec_spec_interact_mod = 1; % eesim_multislice = 1, eesim_phase_object = 2, eesim_weak_phase_object = 3
    input_multem.atomic_pot_parm_typ = 6; % eappt_doyle_0_4 = 1, eappt_peng_0_4 = 2, eappt_peng_0_12 = 3, eappt_kirkland_0_12 = 4, eappt_weickenmeier_0_12 = 5, eappt_lobato_0_12 = 6

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.operation_mode = 1; % eOM_Normal = 1, eOM_Advanced = 2
    input_multem.memory_size = 0; % memory size to be used(Mb)
    input_multem.reverse_multislice = 0; % 1: true, 0:false

    %%%%%%%%%%%%%%%%%%%% atomic vibrations model %%%%%%%%%%%%%%%%%%%%%%
    input_multem.atomic_vib_mod = 1; % eavm_still_atom = 1, eavm_absorptive_pot = 2, eavm_frozen_phonon = 3, eavm_user_def = 4
    input_multem.atomic_vib_coh_contrib = 0; % 1: true, 0: false
    input_multem.atomic_vib_sgl_conf = 0; % 1: true, 0:false (extract single configuration)
    input_multem.atomic_vib_nconf = 1; % true: specific phonon configuration, false: number of frozen phonon configurations
    input_multem.atomic_vib_dim = [true, true, false]; % phonon dimensions [xyz]
    input_multem.atomic_vib_seed = 300183; % Random seed(frozen phonon)

    %%%%%%%%%%%%%%%%%%%%%% specimen information %%%%%%%%%%%%%%%%%%%%%%%
    input_multem.spec_atoms = []; % simulation box length in x direction (Å)

    input_multem.spec_bs_x = 10; % simulation box size in x direction (Å)
    input_multem.spec_bs_y = 10; % simulation box size in y direction (Å)
    input_multem.spec_bs_z = 10; % simulation box size in z direction (Å)

    input_multem.spec_xtl_na = 1; % number of unit cell along a
    input_multem.spec_xtl_nb = 1; % number of unit cell along b
    input_multem.spec_xtl_nc = 1; % number of unit cell along c
    input_multem.spec_xtl_a = 0; % length along a (Å)
    input_multem.spec_xtl_b = 0; % length along b (Å)
    input_multem.spec_xtl_c = 0; % length along c (Å)
    input_multem.spec_xtl_x0 = 0; % reference position along x direction (Å)
    input_multem.spec_xtl_y0 = 0; % reference position along y direction (Å)
 
    %%%%%%%%%%%%%%%%%%%%%%%%% specimen rotation %%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.spec_rot_theta = 0; % angle (º)
    input_multem.spec_rot_u_0 = [0 0 1]; % unitary vector			
    input_multem.spec_rot_ctr_typ = 1; % 1: geometric center, 2: User define		
    input_multem.spec_rot_ctr_p = [0 0 0]; % rotation point

    %%%%%%%%%%%%%%%%%%%%%%%%% specimen slicing %%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.spec_slic(1).typ = 1; % esst_plns_proj = 1, esst_dz_proj = 2, esst_plns_sub = 3, esst_dz_sub = 4, esst_user_def = 5, esst_auto = 6
    input_multem.spec_slic(1).sli_thick = 2.0; % slice thickness
    input_multem.spec_slic(1).sel_typ = 1; % 1: by tag, 2: by z position 
    input_multem.spec_slic(1).sel_tag = 0; % tag
    input_multem.spec_slic(1).sel_Z = 0; % 0: all atomic numbers
    input_multem.spec_slic(1).sel_z_lim = [0, 0]; % [z_0, z_e]
    input_multem.spec_slic(1).z_plns = []; % vector
    %%%%%%%%%%%%%%%%%%%%%%%% specimen thickness %%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.thick_typ = 1; % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
    input_multem.thick = 0; % Array of thickness (Å)
    
    %%%%%%%%%%%%%%%%%%%%%% x-y sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.nx = 512; % number of pixels in x direction
    input_multem.ny = 512; % number of pixels in y direction
    input_multem.bwl = 0; % Band-width limit, 1: true, 0:false

    %%%%%%%%%%%%%%%% electron microscopy simulation type %%%%%%%%%%%%%%%
    % eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, 
    % eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, 
    % eTEMST_EWFS=51, eTEMST_EWRS=52, eTEMST_STEM_EELS=61, eTEMST_ISTEM_EELS=62, 
    % eTEMST_EFTEMFS=71, eTEMST_EFTEMRS=72, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72,
    % eTEMST_PPFS=81, eTEMST_PPRS=82, eTEMST_TFFS=91, eTEMST_TFRS=92
    input_multem.em_sim_typ = 52; % electron microscopy simulation type

    %%%%%%%%%%%%%%%%%%%% microscope parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.E_0 = 300; % Acceleration Voltage (keV)
    input_multem.theta = 0.0; % Polar angle (º)
    input_multem.phi = 0.0; % Azimuthal angle (º)

    %%%%%%%%%%%%%%%%%%%%%% Illumination model %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.illum_mod = 2; % 1: coherente mode, 2: Partial coherente mode, 3: transmission cross coefficient, 4: Numerical integration
    input_multem.illum_inc = 1; % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % On the optimun probe in aberration corrected ADF-STEM
    % Ultramicroscopy 111(2014) 1523-1530
    % C_{nm} Krivanek --- {A, B, C, D, R}_{n} Haider notation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.cond_lens_m = 0;           % Vortex momentum

    input_multem.cond_lens_c_10 = 14.0312;  % [C1]      Defocus (Å)
    input_multem.cond_lens_c_12 = 0.00;     % [A1]      2-fold astigmatism (Å)
    input_multem.cond_lens_phi_12 = 0.00;   % [phi_A1]	Azimuthal angle of 2-fold astigmatism (º)

    input_multem.cond_lens_c_21 = 0.00;     % [B2]      Axial coma (Å)
    input_multem.cond_lens_phi_21 = 0.00;   % [phi_B2]	Azimuthal angle of axial coma (º)
    input_multem.cond_lens_c_23 = 0.00;     % [A2]      3-fold astigmatism (Å)
    input_multem.cond_lens_phi_23 = 0.00;   % [phi_A2]	Azimuthal angle of 3-fold astigmatism (º)

    input_multem.cond_lens_c_30 = 1e-03;    % [C3] 		3rd order spherical aberration (mm)
    input_multem.cond_lens_c_32 = 0.00;     % [S3]      Axial star aberration (Å)
    input_multem.cond_lens_phi_32 = 0.00;   % [phi_S3]	Azimuthal angle of axial star aberration (º)
    input_multem.cond_lens_c_34 = 0.00;     % [A3]      4-fold astigmatism (Å)
    input_multem.cond_lens_phi_34 = 0.0;    % [phi_A3]	Azimuthal angle of 4-fold astigmatism (º)

    input_multem.cond_lens_c_41 = 0.00;     % [B4]      4th order axial coma (Å)
    input_multem.cond_lens_phi_41 = 0.00;   % [phi_B4]	Azimuthal angle of 4th order axial coma (º)
    input_multem.cond_lens_c_43 = 0.00;     % [D4]      3-lobe aberration (Å)
    input_multem.cond_lens_phi_43 = 0.00;   % [phi_D4]	Azimuthal angle of 3-lobe aberration (º)
    input_multem.cond_lens_c_45 = 0.00;     % [A4]      5-fold astigmatism (Å)
    input_multem.cond_lens_phi_45 = 0.00;   % [phi_A4]	Azimuthal angle of 5-fold astigmatism (º)

    input_multem.cond_lens_c_50 = 0.00;     % [C5]      5th order spherical aberration (mm)
    input_multem.cond_lens_c_52 = 0.00;     % [S5]      5th order axial star aberration (Å)
    input_multem.cond_lens_phi_52 = 0.00;   % [phi_S5]	Azimuthal angle of 5th order axial star aberration (º)
    input_multem.cond_lens_c_54 = 0.00;     % [R5]      5th order rosette aberration (Å)
    input_multem.cond_lens_phi_54 = 0.00;   % [phi_R5]	Azimuthal angle of 5th order rosette aberration (º)
    input_multem.cond_lens_c_56 = 0.00;     % [A5]      6-fold astigmatism (Å)
    input_multem.cond_lens_phi_56 = 0.00;   % [phi_A5]	Azimuthal angle of 6-fold astigmatism (º)

    input_multem.cond_lens_inner_aper_ang = 0.0;    % Inner aperture (mrad) 
    input_multem.cond_lens_outer_aper_ang = 21.0;   % Outer aperture (mrad)

    %%%%%%%%% defocus spread function %%%%%%%%%%%%
    input_multem.cond_lens_tp_inc_a = 1.0;      % Height proportion of a normalized Gaussian [0, 1]
    input_multem.cond_lens_tp_inc_sigma = 0.5;  % Standard deviation of the source defocus spread for the Gaussian component (Å)
    input_multem.cond_lens_tp_inc_beta = 0.1;   % Standard deviation of the source defocus spread for the Exponential component (Å)
    input_multem.cond_lens_tp_inc_npts = 10;    % Number of integration points. It will be only used if illum_mod=4
    
    %%%%%%%%%% source spread function %%%%%%%%%%%%
    input_multem.cond_lens_spt_inc_a = 1.0;        % Height proportion of a normalized Gaussian [0, 1]
    input_multem.cond_lens_spt_inc_sigma = 0.5;    % Standard deviation of the source spread function for the Gaussian component: For parallel ilumination(Å^-1);otherwise (Å)
    input_multem.cond_lens_spt_inc_beta = 0.1;     % Standard deviation of the source spread function for the Exponential component: For parallel ilumination(Å^-1);otherwise (Å)
    input_multem.cond_lens_spt_inc_rad_npts = 8;   % Number of radial integration points. It will be only used if illum_mod=4
    input_multem.cond_lens_spt_inc_azm_npts = 8;   % Number of radial integration points. It will be only used if illum_mod=4
    
    %%%%%%%%% zero defocus reference %%%%%%%%%%%%
    input_multem.cond_lens_zero_def_typ = 1;    % eZDT_First = 1, eZDT_User_Define = 2
    input_multem.cond_lens_zero_def_plane = 0;  % It will be only used if cond_lens_zero_def_typ = eZDT_User_Define

    %%%%%%%%%%%%%%%%%%%%%%%% Objective lens %%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.obj_lens_m = 0;            % Vortex momentum

    input_multem.obj_lens_c_10 = 14.0312;   % [C1]      Defocus (Å)
    input_multem.obj_lens_c_12 = 0.00;      % [A1]      2-fold astigmatism (Å)
    input_multem.obj_lens_phi_12 = 0.00;    % [phi_A1]	Azimuthal angle of 2-fold astigmatism (º)

    input_multem.obj_lens_c_21 = 0.00;      % [B2]      Axial coma (Å)
    input_multem.obj_lens_phi_21 = 0.00;    % [phi_B2]	Azimuthal angle of axial coma (º)
    input_multem.obj_lens_c_23 = 0.00;      % [A2]      3-fold astigmatism (Å)
    input_multem.obj_lens_phi_23 = 0.00;    % [phi_A2]	Azimuthal angle of 3-fold astigmatism (º)

    input_multem.obj_lens_c_30 = 1e-03;     % [C3] 		3rd order spherical aberration (mm)
    input_multem.obj_lens_c_32 = 0.00;      % [S3]      Axial star aberration (Å)
    input_multem.obj_lens_phi_32 = 0.00;    % [phi_S3]	Azimuthal angle of axial star aberration (º)
    input_multem.obj_lens_c_34 = 0.00;      % [A3]      4-fold astigmatism (Å)
    input_multem.obj_lens_phi_34 = 0.0;     % [phi_A3]	Azimuthal angle of 4-fold astigmatism (º)

    input_multem.obj_lens_c_41 = 0.00;      % [B4]      4th order axial coma (Å)
    input_multem.obj_lens_phi_41 = 0.00;    % [phi_B4]	Azimuthal angle of 4th order axial coma (º)
    input_multem.obj_lens_c_43 = 0.00;      % [D4]      3-lobe aberration (Å)
    input_multem.obj_lens_phi_43 = 0.00;    % [phi_D4]	Azimuthal angle of 3-lobe aberration (º)
    input_multem.obj_lens_c_45 = 0.00;      % [A4]      5-fold astigmatism (Å)
    input_multem.obj_lens_phi_45 = 0.00;    % [phi_A4]	Azimuthal angle of 5-fold astigmatism (º)

    input_multem.obj_lens_c_50 = 0.00;      % [C5]      5th order spherical aberration (mm)
    input_multem.obj_lens_c_52 = 0.00;      % [S5]      5th order axial star aberration (Å)
    input_multem.obj_lens_phi_52 = 0.00;    % [phi_S5]	Azimuthal angle of 5th order axial star aberration (º)
    input_multem.obj_lens_c_54 = 0.00;      % [R5]      5th order rosette aberration (Å)
    input_multem.obj_lens_phi_54 = 0.00;    % [phi_R5]	Azimuthal angle of 5th order rosette aberration (º)
    input_multem.obj_lens_c_56 = 0.00;      % [A5]      6-fold astigmatism (Å)
    input_multem.obj_lens_phi_56 = 0.00;    % [phi_A5]	Azimuthal angle of 6-fold astigmatism (º)

    input_multem.obj_lens_inner_aper_ang = 0.0;     % Inner aperture (mrad) 
    input_multem.obj_lens_outer_aper_ang = 24.0;	% Outer aperture (mrad)

    %%%%%%%%% defocus spread function %%%%%%%%%%%%
    input_multem.obj_lens_tp_inc_a = 1.0;       % Height proportion of a normalized Gaussian [0, 1]
    input_multem.obj_lens_tp_inc_sigma = 0.5;   % Standard deviation of the source defocus spread for the Gaussian component (Å)
    input_multem.obj_lens_tp_inc_beta = 0.1;    % Standard deviation of the source defocus spread for the Exponential component (Å)
    input_multem.obj_lens_tp_inc_npts = 10;     % Number of integration points. It will be only used if illum_mod=4
    
    %%%%%%%%% zero defocus reference %%%%%%%%%%%%
    input_multem.obj_lens_zero_def_typ = 3;     % eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
    input_multem.obj_lens_zero_def_plane = 0;	% It will be only used if obj_lens_zero_def_typ = eZDT_User_Define
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Scan tag for ISTEM/STEM/EELS %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.scan_pat_typ = 1;	% espt_line = 1, espt_area = 2, espt_user_def = 3
    input_multem.scan_pat_pbc = 1;	% periodic boundary conditions: 1: true, 0:false
    input_multem.scan_pat_spxs = 1; % square pixel size: 1: true, 0:false
    input_multem.scan_pat_nsp = 10; % number of sampling points
    input_multem.scan_pat_r_0 = [0.0; 0.0]; % starting point (Å)
    input_multem.scan_pat_r_e = [4.078; 4.078]; % final point (Å)
    input_multem.scan_pat_r = []; % user define scanning pattern. It will be only used if User_Define=3
    
    % if input_multem.scan_pat_typ = eST_User_Define, then the beam positions
    % must be define in input_multem.beam_pos
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.incdt_wav_typ = 4; % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
    input_multem.incdt_wav_psi = 0; % User define incident wave. It will be only used if User_Define=3

    %%%%%%%%%%%%%%%%%%%%%%%%%%% beam positions %%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.beam_pos = [0.0; 0.0]; % x-y positions
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% STEM Detector %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.detector.typ = 1; % eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3

    input_multem.detector.cir(1).inner_ang = 60; % Inner angle(mrad) 
    input_multem.detector.cir(1).outer_ang = 180; % Outer angle(mrad)

    input_multem.detector.radial(1).x = 0; % radial detector angle(mrad)
    input_multem.detector.radial(1).fx = 0; % radial sensitivity value

    input_multem.detector.matrix(1).R = 0; % 2D detector angle(mrad)
    input_multem.detector.matrix(1).fR = 0; % 2D sensitivity value
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.ped_nrot = 360; % Number of orientations
    input_multem.ped_theta = 3.0; % Precession angle (degrees)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HCI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.hci_nrot = 360; % number of orientations
    input_multem.hci_theta = 3.0; % Precession angle (degrees)

    %%%%%%%%%%%%%%%%%%%%%%%%%% STEM-EELS %%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.eels_Z = 79; % atomic type
    input_multem.eels_E_loss = 80; % energy loss (eV)
    input_multem.eels_collection_angle = 100; % Collection half angle (mrad)
    input_multem.eels_m_selection = 3; % selection rule
    input_multem.eels_channelling_type = 1; % eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% EFTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.eftem_Z = 79; % atomic type
    input_multem.eftem_E_loss = 80; % energy loss (eV)
    input_multem.eftem_collection_angle = 100; % Collection half angle (mrad)
    input_multem.eftem_m_selection = 3; % selection rule
    input_multem.eftem_channelling_type = 1; % eCT_Single_Channelling = 1, eCT_Mixed_Channelling = 2, eCT_Double_Channelling = 3 

    %%%%%%%%%%%%%%%%%%%%%%%%% output tag %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% This option is not used for eTEMST_STEM and eTEMST_STEM_EELS %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multem.output_area_ip_0 = [1; 1]; % Starting position in pixels
    input_multem.output_area_ip_e = [1; 1]; % End position in pixels
    
    %%%%%%%%%%%%%%%%%%% condenser lens variable %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% it will be active in the future %%%%%%%%%%%%%%%%%
    % 1: Vortex momentum, 2: Defocus (Å), 2: Third order spherical aberration (mm)
    % 3: Third order spherical aberration (mm), 4: Fifth order spherical aberration (mm)
    % 5: Twofold astigmatism (Å), 2: Defocus (Å), 6: Azimuthal angle of the twofold astigmatism (º)
    % 7: Threefold astigmatism (Å), 8: Azimuthal angle of the threefold astigmatism (º)
    % 9: Inner aperture (mrad), 2: Defocus (Å), 10: Outer aperture (mrad)
    % input_multem.cdl_var_type = 0; % 0:off 1: m, 2: f, 3 Cs3, 4:Cs5, 5:mfa2, 6:afa2, 7:mfa3, 8:afa3, 9:inner_aper_ang , 10:outer_aper_ang
    % input_multem.cdl_var = [-2 -1 0 1 2]; % variable array
end