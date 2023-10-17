function [im_o]= ilm_lcwt(im_i, sigma)
    n = 2*ceil(sigma)+1;
    kn = ones(n,n);
    kn = kn/sum(kn(:));
    
    im_o = (im_i - mean(im_i, 'all'))/std(im_i, 0, 'all');
    
    im_mean_0 = conv2(im_o, kn,'same'); 
    im_mean = imgaussfilt(im_mean_0, sigma, 'FilterSize', 2*ceil(2.5*sigma)+1, 'Padding', 'symmetric');

    im_std_0 = (im_o - im_mean).^2;
    im_std_0 = conv2(im_std_0, kn, 'same'); 
    im_std_0 = max(0.01, im_std_0);
    im_std = imgaussfilt(im_std_0, sigma, 'FilterSize', 2*ceil(2.5*sigma)+1, 'Padding', 'symmetric');
    
    if 0
        figure(1);clf;
        subplot(2, 2, 1);
        imagesc(im_mean_0);
        axis image off;
        colorbar;
        subplot(2, 2, 2);
        imagesc(im_mean);
        axis image off;
        colorbar;
        subplot(2, 2, 3);
        imagesc(im_std_0);
        axis image off;
        colorbar;
        subplot(2, 2, 4);
        imagesc(im_std);
        axis image off;
        colorbar;
    end
    
    im_o = (im_o - im_mean)./im_std;
    im_o = (im_o - mean(im_o, 'all'))/std(im_o, 0, 'all');
end