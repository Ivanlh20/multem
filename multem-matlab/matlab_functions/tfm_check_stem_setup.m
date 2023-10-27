function tfm_check_stem_setup(input_multem, n_detector)

    fcn_check_probe(input_multem);
    fcn_check_slicing(input_multem);
    if nargin > 1
        fcn_check_detector_coverage(input_multem, n_detector);
    end

    function fcn_check_probe(input_multem)
        input_multem.system_conf.precision = 1;                           % eP_Float = 1, eP_double = 2
        input_multem.system_conf.device = 1;                              % eD_CPU = 1, eD_GPU = 2
        input_multem.system_conf.cpu_nthread = 1; 
        input_multem.system_conf.gpu_device = 0;
        
        input_multem.iw_type = 2;                       % 2 = convergent beam
        input_multem.iw_x = 0.5*input_multem.spec_lx;
        input_multem.iw_y = 0.5*input_multem.spec_ly;

        output_incident_wave = input_multem.ilc_incident_wave; 

        psi_0 = output_incident_wave.psi_0;
        figure(1); clf
        subplot(1,2,1)
        imagesc(0:output_incident_wave.dx:(input_multem.spec_lx),...
                0:output_incident_wave.dy:(input_multem.spec_ly),...
                abs(psi_0).^2);
        title('Probe profile (intensity)');
        axis image;
        set(gca,'YDir','normal');
        colormap jet;
        
        subplot(1,2,2)
        imagesc(0:output_incident_wave.dx:(input_multem.spec_lx),...
                0:output_incident_wave.dy:(input_multem.spec_ly),...
                angle(psi_0));
        title('Probe profile (Phase)');
        axis image;
        set(gca,'YDir','normal');
        colormap jet;
    end

    function fcn_check_slicing(input_multem)
        
        if ~isempty(input_multem.spec_atoms)
            % Slicing
            [~, Slice] = ilc_spec_slicing(input_multem.toStruct);   
            [nslice, ~] = size(Slice);

            figure(2); clf;
            subplot(1,2,1);
            plot(input_multem.spec_atoms(:, 2), input_multem.spec_atoms(:, 4), 'ok');   
            hold on;
            for i = 1:nslice
                plot([min(input_multem.spec_atoms(:, 2)) max(input_multem.spec_atoms(:, 2))], [Slice(i, 1) Slice(i, 1)], '-b', [min(input_multem.spec_atoms(:, 2)) max(input_multem.spec_atoms(:, 2))], [Slice(i, 2) Slice(i, 2)], '-b');    
            end
            hold off;
            title('Slice positions');
            ylabel('z');
            xlabel('x');
            axis equal;
            axis([min(input_multem.spec_atoms(:, 2))-3 max(input_multem.spec_atoms(:, 2))+3 Slice(1, 1)-3 Slice(end, 1)+3]);

            legend({'Atom positions','Slice boundaries'},'Location','southoutside')
            set(gca,'YDir','reverse')

            subplot(1,2,2);
            plot3(input_multem.spec_atoms(:,2), input_multem.spec_atoms(:,3), input_multem.spec_atoms(:,4),'or');
            hold on;
            if input_multem.simulation_type <= 12 
                plot3([input_multem.scanning_x0; input_multem.scanning_xe; input_multem.scanning_xe; input_multem.scanning_x0;input_multem.scanning_x0],...
                      [input_multem.scanning_y0; input_multem.scanning_y0; input_multem.scanning_ye; input_multem.scanning_ye;input_multem.scanning_y0],...
                      [0; 0; 0; 0; 0],'-k');
                legend({'Atom positions','STEM Scan Field'},'Location','southoutside')
            else
                z0 = get_defocus_ref(input_multem);
                [xp, yp, zp] = cylinder([0 2], 64);
                surf(input_multem.iw_x+xp, input_multem.iw_y+yp, z0+(10*zp),'FaceAlpha',0.5,'EdgeColor','none','FaceColor','g')
                surf(input_multem.iw_x+xp, input_multem.iw_y+yp, z0+(-10*zp),'FaceAlpha',0.5, 'EdgeColor','none','FaceColor','g')
                plot3(input_multem.iw_x, input_multem.iw_y, z0,'xk');
                legend({'Atom positions','Probe'},'Location','southoutside')
            end

            hold off;
            axis equal;
            title('Atom positions and Scan Field')
            zlabel('z');
            xlabel('x');
            ylabel('y');
            set(gca,'ZDir','reverse')
        else
            fprintf(2,['No atoms are specified. Check your input structure! \n']);
        end

    end

    function fcn_check_detector_coverage(input_multem, detector)

        E_0 = input_multem.E_0;
        nx = input_multem.nx;
        ny = input_multem.ny;
        lx = input_multem.spec_lx;
        ly = input_multem.spec_ly;
        ri_detector = input_multem.detector.cir(detector).inner_ang;
        ro_detector = input_multem.detector.cir(detector).outer_ang;
    
        gmax(1)=nx/(2*lx); gmax(2) = ny/(2*ly);
        gmax_abs=min(gmax(:));
        ri_detector = ilm_mrad_2_rAng(E_0,ri_detector);
        ro_detector = ilm_mrad_2_rAng(E_0,ro_detector);
        t = linspace(0, 2*pi);
        xi = ri_detector*cos(t);
        yi = ri_detector*sin(t);
        xo = ro_detector*cos(t);
        yo = ro_detector*sin(t);
        bwx = 2/3*(gmax_abs*cos(t));
        bwy = 2/3*(gmax_abs*sin(t));

        gspace=polyshape([-gmax(1) gmax(1) gmax(1) -gmax(1)],[-gmax(2) -gmax(2) gmax(2) gmax(2)]);
        d=polyshape({xi,xo},{yi,yo},'Simplify',false);
        figure(3); clf
        plot(d); hold on;
        plot(gspace);
        plot(bwx,bwy,'--r')
        hold off;
        legend({'Detector','Simulation Box','Bandwith limit'},'Location','southoutside')
        title('Detector coverage')
        axis equal

        str = ['Real space sampling = ' num2str(round(lx/nx,3)) ' by ' num2str(round(ly/ny,3)) sprintf(' %c ', 197) '\n' 'g_max = ' num2str(gmax_abs,5) '\n'];
        fprintf(str)
        if gmax_abs<ro_detector
            fprintf(2,'Reciprocal space does not cover the detector! \n')
        end
        if tfm_pn_fact(nx,1)~=nx && tfm_pn_fact(nx,2)~=nx && tfm_pn_fact(nx,3)~=nx
           fprintf(2,['Using ' num2str(nx) ' pixels in x is computationally unefficient. Consider using ' num2str(tfm_pn_fact(nx,1)) ' or ' num2str(tfm_pn_fact(nx,3)) ' instead.'  '\n'])
        end
        if tfm_pn_fact(ny,1)~=ny && tfm_pn_fact(ny,2)~=ny && tfm_pn_fact(ny,3)~=ny
           fprintf(2,['Using ' num2str(ny) ' pixels in y is computationally unefficient. Consider using ' num2str(tfm_pn_fact(ny,1)) ' or ' num2str(tfm_pn_fact(ny,3)) ' instead.'  '\n'])
        end
    end
end

function z0 = get_defocus_ref(input_multem)
    at_min = min(input_multem.spec_atoms(:,4));
    at_max = max(input_multem.spec_atoms(:,4));
    switch input_multem.cond_lens_zero_defocus_type 
        case 1
            z0 = at_min;
        case 2
            z0 = (at_max -at_min)/2+at_min;
        case 3
            z0 = at_max;
        case 4
            z0 = input_multem.cond_lens_zero_defocus_plane;
    end
end