function[Probe]=setProbe(x0, y0, input_multislice)
Probe.gpu = input_multislice.gpu;            % Gpu card
Probe.E_0 = input_multislice.E_0;
Probe.theta = 0.0; Probe.phi = 0; % Till ilumination (degrees)
Probe.lx = input_multislice.lx; Probe.ly = input_multislice.ly;
Probe.nx =input_multislice.nx; Probe.ny = input_multislice.ny;
Probe.x0 = x0; Probe.y0 = y0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Probe.m = input_multislice.lens_m;         %mm
Probe.f = input_multislice.lens_f;         %Angs
Probe.Cs3 = input_multislice.lens_Cs3;     %mm
Probe.Cs5 = input_multislice.lens_Cs5;     %mm
Probe.mfa2 = input_multislice.lens_mfa2; Probe.afa2 = input_multislice.lens_afa2;     %(Angs, degrees)
Probe.mfa3 = input_multislice.lens_mfa3; Probe.afa3 = input_multislice.lens_afa3;     %(Angs, degrees)
Probe.aobjl = 0.0; Probe.aobju = 21.0659;               %(mrad, mrad)
Probe.sf = input_multislice.lens_sf; Probe.nsf = input_multislice.lens_nsf;           % (Angs, number of steps)ne
Probe.beta = input_multislice.lens_beta; Probe.nbeta = input_multislice.lens_nbeta;	%(mrad, half number of steps)