clear; clc;
load('fitting_coef.mat');

fid = fopen('feg_lobato.txt', 'w');
c_l = [c_l, zeros(size(c_l, 1), 1)];
c_nl = [c_nl, zeros(size(c_nl, 1), 1)];
[nr, nc] = size(c_l);
for Z = 1:nr
    str = '';
    for ic = 1:nc
        str_cl = ['data_table[', num2str(Z-1), '].feg[0].cl[', num2str(ic-1), '] = ', num2str(c_l(Z, ic), '%.8e'), ';'];
        str_cnl = ['data_table[', num2str(Z-1), '].feg[0].cnl[', num2str(ic-1), '] = ', num2str(c_nl(Z, ic), '%.8e'), ';'];
        str = [str, '\t', str_cl, '\t', str_cnl];
    end
    disp(str);
    fprintf(fid, ['\t\t\t', str, '\n']);
end
fclose(fid);