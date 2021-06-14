clear; clc;
load('EDX_pot.mat');

path = 'edx_data';
[fn, n_fn] = ilm_dir(path, '*.dat');
n_Z = 103;
n_coef = 29;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
keySet = 50:50:400;
n_E0 = length(keySet);
valueSet = 1:n_E0;
E0_map = containers.Map(keySet,valueSet);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f, 5d, 6p, 7s, 5f, 6d, 7p, (8s, 5g, 6f, 7d, 8p, and 9s)
% keySet = {'1s','2s','2p','3s','3p','4s','3d','4p','5s','4d','5p'};
keySet = {'1s','2s','2p','3s','3p','4s','3d','4p','4d','5p'};
n_orb = length(keySet);
valueSet = 1:n_orb;
orb_map = containers.Map(keySet,valueSet);

% s
s = [0.0, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.625, 0.75, 0.875, 1.0,...
    1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0];

% g
g_0 = 2*s;
n_g = 2^15;
n_r = 128;

% set initial paramters
clear data_table;
for iZ=1:n_Z
    for iE_0=1:n_E0
        for iorb=1:n_orb
            data_table(iZ).orb(iorb).E0(iE_0).g = zeros(n_coef, 1);
            data_table(iZ).orb(iorb).E0(iE_0).fg = zeros(n_coef, 1);
            data_table(iZ).orb(iorb).E0(iE_0).r = zeros(n_r, 1);
            data_table(iZ).orb(iorb).E0(iE_0).vp = zeros(n_r, 1);
            data_table(iZ).orb(iorb).E0(iE_0).r_max = 0;
        end
    end
end

% read data
for it = 1:n_fn
    fid = fopen(fn{it}, 'r');
    while (~feof(fid))
        out_text = fgetl(fid);
        Z_exist = strfind(out_text, ' Z =');
        if Z_exist>0
            % Z
            out_text = strtrim(split(out_text, '='));
            out_text = textscan(out_text{2}, '%f%s');
            Z = out_text{1};
            
            % orbital
            orb = char(out_text{2});
            idx_orb = orb_map(orb);

            clear vr_ss_tmp
            for iE_0=1:n_E0
                % E_0
                out_text = fgetl(fid);
                out_text = strtrim(split(out_text, '='));
                out_text = textscan(out_text{2}, '%f%s');
                E_0 = out_text{1};
                idx_E0 = E0_map(E_0);

                % fedxg
                data = textscan(fid, '%f', 'delimiter','\t', 'HeaderLines', 5);
                fg = cell2mat(data).';
                n_coef = length(fg);
                g = g_0(1:n_coef);

                data_table(Z).orb(idx_orb).E0(idx_E0).g = g;
                data_table(Z).orb(idx_orb).E0(idx_E0).fg = fg;

                disp([num2str(Z), ' - ', num2str(E_0), ' - ', orb, ' - ', num2str(n_coef)])
                disp([Z, idx_E0, idx_orb])
            end
        end
    end
    fclose(fid);
end

% saving data
save('data_table_edx.mat', 'data_table');