function [Obj, Probe, Psi] = ilm_ePIE(Obj_0, Probe_0, CBED, n_iter, alpha, beta)
    CBED = single(CBED);
    M = sqrt(ifftshift(CBED));

    ifft2_sc = sqrt(size(CBED, 1)*size(CBED, 1));
    cbed_sc = sqrt(sum(CBED, 'all'));

    probe_0_sc = cbed_sc/sqrt(sum(abs(Probe_0).^2, 'all'));
    Probe_0 = Probe_0*probe_0_sc;

    obj_0_sc = cbed_sc/sqrt(sum(abs(Obj_0.*Probe_0).^2, 'all'));
    Obj_0 = Obj_0*obj_0_sc;

    for ik=1:n_iter
        Psi_0 = Obj_0.*Probe_0;
        Psi = fft2(ifftshift(Psi_0));
        Psi = M.*exp(1j*angle(Psi));
        Psi = fftshift(ifft2(Psi*ifft2_sc));

        Obj_0_2 = abs(Obj_0).^2;
        Probe = Probe_0 + beta*conj(Obj_0).*(Psi-Psi_0)./max(Obj_0_2, [], 'all');
        Probe = Probe*cbed_sc/sqrt(sum(abs(Probe).^2, 'all'));

        Probe_0_2 = abs(Probe_0).^2;
        Obj = Obj_0 + alpha*conj(Probe_0).*(Psi-Psi_0)./max(Probe_0_2, [], 'all');
        Obj = Obj*cbed_sc/sqrt(sum(abs(Obj.*Probe).^2, 'all'));   

%         disp([sum(abs(psi_0(:)).^2), sum(abs(psi(:)).^2), sum(CBED(:))]/256^2)
%         disp([sum(abs(Probe_0(:)).^2), sum(abs(Probe(:)).^2)]/256^2)
%         disp([sum(abs(Obj_0.*Probe_0).^2, 'all'), sum(abs(Obj.*Probe).^2, 'all')]/256^2)

        if 0
            figure(2); clf;
            subplot(2, 2, 1)
            imagesc(abs(Psi_0));
            axis image;
            title('Module psi_0');
            colorbar;
            subplot(2, 2, 3)
            imagesc(angle(Psi_0));
            axis image;
            colorbar;
            title('Phase psi_0')

            subplot(2, 2, 2)
            imagesc(abs(Psi));
            axis image;
            title('Module psi');
            colorbar;
            subplot(2, 2, 4)
            imagesc(angle(Psi));
            axis image;
            colorbar;
            title('Phase psi');
        end

        if 0                        
            figure(1); clf;
            subplot(2, 4, 1)
            imagesc(abs(Obj_0));
            axis image;
            title('Module Obj');
            colorbar;
            subplot(2, 4, 5)
            imagesc(angle(Obj_0));
            axis image;
            colorbar;
            title('Phase Obj')

            subplot(2, 4, 2)
            imagesc(abs(Probe_0));
            axis image;
            title('Module Probe');
            colorbar;
            subplot(2, 4, 6)
            imagesc(angle(Probe_0));
            axis image;
            colorbar;
            title('Phase Probe');
            
            subplot(2, 4, 3)
            imagesc(abs(Obj));
            axis image;
            title('Module Obj');
            colorbar;
            subplot(2, 4, 7)
            imagesc(angle(Obj));
            axis image;
            colorbar;
            title('Phase Obj')

            subplot(2, 4, 4)
            imagesc(abs(Probe));
            axis image;
            title('Module Probe');
            colorbar;
            subplot(2, 4, 8)
            imagesc(angle(Probe));
            axis image;
            colorbar;
            title('Phase Probe');                   
        end

        Obj_0 = Obj;          
        Probe_0 = Probe;
    end


    Probe = Probe/probe_0_sc;
    Obj = Obj/obj_0_sc;
    Psi = Psi/(obj_0_sc*probe_0_sc);
end