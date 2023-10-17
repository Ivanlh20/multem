function[rect] = ilm_bd_2_rect(lx, ly, bd_px)
    rect = [bd_px(1), bd_px(3); lx-bd_px(2), ly-bd_px(4)];
end