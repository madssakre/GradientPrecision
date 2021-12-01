clear
close all

Run_this = 'One_of_each'; % sit back and wait. It will cycle through all

% % or uncomment this
% Run_this = 'Just_one'; % Choose Gradient_Type and ROI_Type below

% % and some of these
% Gradient_Type = 'standard_0';
% Gradient_Type = 'standard_1';
% Gradient_Type = 'exact';
% Gradient_Type = 'midpoint';
%
% ROI_Type = 'Gr';
% ROI_Type = 'GrBW';
% ROI_Type = 'BW';


load('Time.mat')

load('Opimization.mat')


switch Run_this
    
    case 'Just_one'
        
        load(sprintf('Space_%s.mat',ROI_Type))
        
        Optimization.Par.Grad = Gradient_Type;
        
        fprintf(2,'Running GRADPREC with a %s gradient on a %s type of ROI\n\n',Optimization.Par.Grad,ROI_Type)
        
        if contains(Optimization.Par.Grad,'exa')
            
            Optimization.par_Ncores = 2;
            % or e.g. Optimization.par_Ncores = feature('numcores');
            % or 23 as in the paper.
            if Optimization.par_Ncores > 1
                parpool;
            end
            
            fprintf(1,'The exact gradient is run with %i CPU cores. You might want to adjust this...\n\n',Optimization.par_Ncores)
            
        end
        
        
        Optimization = GRADPREC(Space,Time,Optimization);
        
        Simulation = SIMULATION([],Space,Time,[],'u',Optimization.uo,'v',Optimization.vo,'s',[],'Output','MapsSingle');
        
        figure
        colormap jet
        subplot(2,3,1)
        imagesc(Space.Mtarget_map(:,:,:,:,1),[-1,1])
        colorbar
        axis equal off
        title({'';'M_x/M_0'})
        subplot(2,3,2)
        imagesc(Space.Mtarget_map(:,:,:,:,2),[-1,1])
        title({sprintf('Target of a %s type',ROI_Type);'M_y/M_0'})
        colorbar
        axis equal off
        subplot(2,3,3)
        imagesc(Space.Mtarget_map(:,:,:,:,3),[-1,1])
        title({'';'M_z/M_0'})
        colorbar
        axis equal off
        subplot(2,3,4)
        imagesc(Simulation.M_t{1}(:,:,:,:,1),[-1,1])
        title('M_x(T)/M_0')
        colorbar
        axis equal off
        subplot(2,3,5)
        imagesc(Simulation.M_t{1}(:,:,:,:,2),[-1,1])
        title({sprintf('GRADPREC after %i iterations with a %s gradient\n\n',Optimization.ksafe-1,Optimization.Par.Grad);'M_y(T)/M_0'})
        colorbar
        axis equal off
        subplot(2,3,6)
        imagesc(Simulation.M_t{1}(:,:,:,:,3),[-1,1])
        title('M_z(T)/M_0')
        colorbar
        axis equal off
        drawnow
        
        
        
    case 'One_of_each'
        
        Gradient_Type_ = {'standard_0','standard_1','exact','midpoint','standard_0_loop','standard_1_loop'};
        
        
        ROI_Type_ = {'Gr','GrBW','BW'};
        Optimization_array = cell(3,4);
        Simulation_array = cell(3,4);
        for r = 1:3
            
            load(sprintf('Space_%s.mat',ROI_Type_{r}))
            
            for g = 1:4
                
                Optimization.Par.Grad = Gradient_Type_{g};
                
                
                
                fprintf(2,'\n\nRunning GRADPREC with a %s gradient on a %s type of ROI\n\n',Optimization.Par.Grad,ROI_Type_{r})
                
                if contains(Optimization.Par.Grad,'exa')
                    
                    Optimization.par_Ncores = 2;
                    % or e.g. Optimization.par_Ncores = feature('numcores');
                    % or 23 as in the paper.
                    if Optimization.par_Ncores > 1
                        parpool;
                    end
                    
                    fprintf(1,'The exact gradient is run with %i CPU cores. You might want to adjust this...\n\n',Optimization.par_Ncores)
                    
                end
                
                Optimization_array{r,g} = GRADPREC(Space,Time,Optimization);
                
                Simulation_array{r,g} = SIMULATION([],Space,Time,[],'u',Optimization_array{r,g}.uo,'v',Optimization_array{r,g}.vo,'s',[],'Output','MapsSingle');
                
                
                figure
                colormap jet
                subplot(2,3,1)
                imagesc(Space.Mtarget_map(:,:,:,:,1),[-1,1])
                colorbar
                axis equal off
                title({'';'M_x/M_0'})
                subplot(2,3,2)
                imagesc(Space.Mtarget_map(:,:,:,:,2),[-1,1])
                title({sprintf('Target of a %s type',ROI_Type_{r});'M_y/M_0'})
                colorbar
                axis equal off
                subplot(2,3,3)
                imagesc(Space.Mtarget_map(:,:,:,:,3),[-1,1])
                title({'';'M_z/M_0'})
                colorbar
                axis equal off
                subplot(2,3,4)
                imagesc(Simulation_array{r,g}.M_t{1}(:,:,:,:,1),[-1,1])
                title('M_x(T)/M_0')
                colorbar
                axis equal off
                subplot(2,3,5)
                imagesc(Simulation_array{r,g}.M_t{1}(:,:,:,:,2),[-1,1])
                title({sprintf('GRADPREC after %i iterations with a %s gradient\n\n',Optimization_array{r,g}.ksafe-1,Optimization.Par.Grad);'M_y(T)/M_0'})
                colorbar
                axis equal off
                subplot(2,3,6)
                imagesc(Simulation_array{r,g}.M_t{1}(:,:,:,:,3),[-1,1])
                title('M_z(T)/M_0')
                colorbar
                axis equal off
                drawnow
                
                
                
                
                
            end
            
        end
        
        %%
        
        figure
        for n = 1:3
            subplot(1,3,n)
            
            bar([1,2,3,4],[mean(Optimization_array{n,1}.Durations),mean(Optimization_array{n,2}.Durations),mean(Optimization_array{n,3}.Durations),mean(Optimization_array{n,4}.Durations)])
        end
        
end
