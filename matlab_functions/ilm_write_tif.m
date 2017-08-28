function [] = ilm_write_tif(stack, n_stack, stack_name)
    t = Tiff(stack_name,'w');

    tagstruct.ImageLength = size(stack, 1);
    tagstruct.ImageWidth = size(stack, 2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 8;
    tagstruct.RowsPerStrip = 256;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    t.setTag(tagstruct);

    Nframes = min(n_stack, size(stack, 3));
    imgdata = uint8(255*stack(:, :, 1));
    t.write(imgdata);
    for ik = 2:Nframes
        writeDirectory(t);
        imgdata = uint8(255*stack(:, :, ik));
        t.setTag(tagstruct);
        t.write(imgdata);
    end
    t.close();
end