clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
valueSet = 50:50:400;
n_E0 = length(valueSet);
keySet = 1:n_E0;
E0_map = containers.Map(keySet,valueSet);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f, 5d, 6p, 7s, 5f, 6d, 7p, (8s, 5g, 6f, 7d, 8p, and 9s)
% keySet = {'1s','2s','2p','3s','3p','4s','3d','4p','5s','4d','5p'};
valueSet = {'1s','2s','2p','3s','3p','4s','3d','4p','4d','5p'};
n_orb = length(valueSet);
keySet = 1:n_orb;
orb_map = containers.Map(keySet,valueSet);

% set initial paramters
load('data_table_edx.mat');

n_Z = length(data_table);
n_r = 128;

% read data
str_sft_fcn = '\t\t\t';

str_sft_Z_switch = [str_sft_fcn, '\t'];
str_sft_Z_case = [str_sft_Z_switch, '\t'];
str_sft_Z_return = [str_sft_Z_case, '\t'];

str_sft_Orb_switch = [str_sft_Z_case, '\t'];
str_sft_Orb_case = [str_sft_Orb_switch, '\t'];
str_sft_Orb_return = [str_sft_Orb_case, '\t'];

str_sft_E_switch = [str_sft_Orb_case, '\t'];
str_sft_E_case = [str_sft_E_switch, '\t'];
str_sft_E_return = [str_sft_E_case, '\t'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str_tot = sprintf([str_sft_fcn, 'Edx_vr_cpu<double> load_edx_vr_parm(const dt_int32& Z, const dt_int32& orb_idx, const dt_int32& E_idx)\n']);
str_tot = sprintf([str_tot, str_sft_fcn, '{\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str_Z = sprintf([str_sft_Z_switch, 'switch (Z)\n']); 
str_Z = sprintf([str_Z, str_sft_Z_switch, '{\n']);
str_tot = sprintf([str_tot, str_Z]);

 for iZ=1:n_Z
    % check out if edx potential for a given atomic number Z is available
    Z_bb = false;
    for iorb=1:n_orb
        if data_table(iZ).orb(iorb).E0(1).fg(1)>1e-7
            Z_bb = true;
            break;
        end
    end
    
    if ~ Z_bb
        continue;
    end

    str_Z = sprintf([str_sft_Z_case, 'case ', num2str(iZ), ':\n']);
    str_tot = [str_tot, str_Z];
    
    str_Orb = sprintf([str_sft_Orb_switch 'switch (orb_idx)\n']);  
    str_Orb = sprintf([str_Orb, str_sft_Orb_switch, '{\n']);
    str_tot = [str_tot, str_Orb];
    
    for iorb=1:n_orb
        orb = orb_map(iorb);
        if data_table(iZ).orb(iorb).E0(1).fg(1)<1e-7
            continue;
        end

        str_E = sprintf([str_sft_Orb_case, 'case ', num2str(iorb-1), ':\n']);
        str_E = sprintf([str_E, str_sft_E_switch, 'switch (E_idx)\n']);
        str_E = sprintf([str_E, str_sft_E_switch, '{\n']);
        str_tot = [str_tot, str_E];
        
        for iE=1:n_E0
            E_0 = E0_map(iE);
            
            % read data
            r_max = data_table(iZ).orb(iorb).E0(iE).r_max;
            r = data_table(iZ).orb(iorb).E0(iE).r;
            vp = data_table(iZ).orb(iorb).E0(iE).vp;
            
            disp([num2str(iZ), ' - ', num2str(E_0), ' - ', orb])

            str_vp_parm = fcn_str_vp_parm(r_max, r, vp);

            str_E = sprintf([str_sft_E_case, 'case ', num2str(iE-1), ':\n']);
            str_E = sprintf([str_E, str_sft_E_return, 'return ', str_vp_parm, ';\n']);
            str_tot = [str_tot, str_E];
        end
        str_vp_parm = fcn_str_vp_parm(0.0, zeros(1, n_r), zeros(1, n_r));
        str_E = sprintf([str_sft_E_case, 'default:\n']);
        str_E = sprintf([str_E, str_sft_E_return, 'return ', str_vp_parm, ';\n']);
        str_E = sprintf([str_E, str_sft_E_switch, '}\n']);
        str_tot = [str_tot, str_E];
    end
    str_vp_parm = fcn_str_vp_parm(0.0, zeros(1, n_r), zeros(1, n_r));
    str_Orb = sprintf([str_sft_Orb_case, 'default:\n']);
    str_Orb = sprintf([str_Orb, str_sft_Orb_return, 'return ', str_vp_parm, ';\n']);
    str_Orb = sprintf([str_Orb, str_sft_Orb_switch, '}\n']);
    str_tot = [str_tot, str_Orb];
 end
 
str_vp_parm = fcn_str_vp_parm(0.0, zeros(1, n_r), zeros(1, n_r));
str_Z = sprintf([str_sft_Z_case, 'default:\n']);
str_Z = sprintf([str_Z, str_sft_Z_return, 'return ', str_vp_parm, ';\n']);
str_Z = sprintf([str_Z, str_sft_Z_switch, '}\n']);
str_tot = [str_tot, str_Z];

str_tot = sprintf([str_tot, str_sft_fcn, '}\n']);

fid = fopen('edx_vr_parm.txt','wt');
fprintf(fid, '%s', str_tot);
fclose(fid);

% saving data
% save('data_table_edx.mat', 'data_table');

function [str_vp_parm]= fcn_str_vp_parm(r_max, r, vp)
    n_r = length(r);

    % r_max
    str_r_max = num2str(r_max, '%.7e');

    % r
    str_r = '{';
    for ir = 1:n_r
        sp = ilm_ifelse(ir<n_r, ', ', '');
        str_r = [str_r, num2str(r(ir), '%.7e'), sp];
    end
    str_r = [str_r, '}'];

    % vp
    str_vp = '{';
    for ir = 1:n_r
        sp = ilm_ifelse(ir<n_r, ', ', '');
        str_vp = [str_vp, num2str(vp(ir), '%.7e'), sp];
    end
    str_vp = [str_vp, '}'];

    str_vp_parm = ['{', str_r_max, ', ', str_r, ', ', str_vp '}'];
end