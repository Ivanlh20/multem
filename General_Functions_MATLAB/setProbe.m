function[Probe]=setProbe(x0, y0, TEM)
Probe.gpu = TEM.gpu;            % Gpu card
Probe.E0 = TEM.E0;
Probe.theta = 0.0; Probe.phi = 0; % Till ilumination (degrees)
Probe.lx = TEM.lx; Probe.ly = TEM.ly;
Probe.nx =TEM.nx; Probe.ny = TEM.ny;
Probe.x0 = x0; Probe.y0 = y0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Probe.m = TEM.MC.m;         %mm
Probe.f = TEM.MC.f;         %Angs
Probe.Cs3 = TEM.MC.Cs3;     %mm
Probe.Cs5 = TEM.MC.Cs5;     %mm
Probe.mfa2 = TEM.MC.mfa2; Probe.afa2 = TEM.MC.afa2;     %(Angs, degrees)
Probe.mfa3 = TEM.MC.mfa3; Probe.afa3 = TEM.MC.afa3;     %(Angs, degrees)
Probe.aobjl = 0.0; Probe.aobju = 21.0659;               %(mrad, mrad)
Probe.sf = TEM.MC.sf; Probe.nsf = TEM.MC.nsf;           % (Angs, number of steps)ne
Probe.beta = TEM.MC.beta; Probe.nbeta = TEM.MC.nbeta;	%(mrad, half number of steps)