function[image_av]=ilm_mean_image(data)

n_images = length(data);
[ny, nx] = size(data(1).image);
image_av = zeros(ny, nx);
for i = 1:n_images
    image_av = image_av + data(i).image;
end
image_av = image_av/n_images;