clear all; clc;
Im = importdata('stem_0.mat');
Im = Im(100:1123, 100:1123);
figure(1);
imagesc(Im);
colormap gray;
axis image;

np_min = 35;
delta = 0.04175;

tic;
y = il_projected_intensity(Im, np_min, delta);
toc;

figure(2); clf;
subplot(2, 1, 1);
plot(y, '-r');

d_delta = 1;
tic;
[angle, psd, peaks] = il_psd(Im, np_min, d_delta);
toc;

subplot(2, 1, 2);
plot(angle, psd, '-r');
t_min = min(psd); 
t_max = max(psd);
for i=1:length(peaks)
    hold on;
    plot([peaks(i) peaks(i)], [t_min t_max], '-b');
end;

nmed = max(1, round(3/d_delta));
tt = psd - il_filter_median(psd, nmed);