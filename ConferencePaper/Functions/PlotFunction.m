function fig = PlotFunction(PlotDataCell,T_ini, r, tit, leg, obsv_state, plot_lims, Ts)
    % This function plots results from a simulation step
    % To be given in the PlotDataCell format, as specified in the script
    % above

    % This script assumes single input dimension, and up to 4 states for
    % plots

    fig = figure();
    colorlist = {'blue', 'r', "#77AC30", 'g', 'k'};
    linesize = 1.5;
    [n_states, L_sim] = size(PlotDataCell{1,1});

    
    if ~exist("leg")
        leg = cell(size(PlotDataCell,2),1);
    end
    
    if ~exist("obsv_state")
        obsv_state = 1;
    end
     if ~exist("Ts")
        Ts = 1;
     end
    t = (1:L_sim)*Ts;

    if n_states <= 3
        plot_dims = [n_states+1, 1];
        u_plot_range = n_states+1;
    else
        plot_dims = [3,2];
        u_plot_range = [5,6];
    end
    %determine maxima:
    fullmat = cell2mat(PlotDataCell);

    % the "plot_lims" variable has 3 options: 
    % 1. "max" sets the plots to the max value
    % 2. plot_lims given as a matrix [x_min, x_max; u_min, u_max]
    % 3. No info given, then all is set to [-2, 2]
    if exist("plot_lims", "var")
        if strcmpi(plot_lims, "max")
            x_max = 1.2*max(fullmat(1:n_states,:),[],2); x_min = 1.2*min(fullmat(1:n_states,:),[],2);
            u_max = 1.2*max(fullmat(end,:),[],2); u_min = 1.2*min(fullmat(end,:),[],2);
        
            % Scaling up the axis limits
            x_max = x_max + (0.1.*(x_max-x_min).*ones(n_states,1));
            x_min = x_min - (0.1.*(x_max-x_min).*ones(n_states,1));
            u_max = u_max + 0.1*(u_max-u_min);
            u_min = u_min - 0.1*(u_max-u_min);
        else
            x_max = plot_lims(1,2)*ones(n_states,1); x_min = plot_lims(1,1)*ones(n_states,1);
            u_max = plot_lims(2,2); u_min = plot_lims(2,1);
        end
    else
        x_max = 2*ones(n_states,1); x_min = -x_max; u_max = 2; u_min = -2;
    end
    % Loop over each trajectory in the PlotDataCell set
    for i_controller_type = 1:size(PlotDataCell,2)
        x_traj = PlotDataCell{1,i_controller_type};
        u_traj = PlotDataCell{2,i_controller_type};
        col = colorlist{i_controller_type};

        for i_subplot = 1:n_states
            % % plotting x_i
            subplot(plot_dims(1),plot_dims(2),i_subplot)
            stairs(t, x_traj(i_subplot,:), 'Color', col,"LineWidth", linesize, "DisplayName",leg{i_controller_type}); 

            hold on; grid on
            xline(T_ini, 'black:','LineWidth', 1, "HandleVisibility", 'off');
            ylabel(strcat('$x_', num2str(i_subplot),'$'))
            ylim([x_min(i_subplot), x_max(i_subplot)]);
            xlim([0, L_sim*Ts])
        end

        subplot(plot_dims(1),plot_dims(2),u_plot_range)
        stairs(t, u_traj(1,:), 'Color', col,"LineWidth",linesize, "DisplayName",leg{i_controller_type}); 

        hold on; grid on
        xline(T_ini, 'black:','LineWidth', 1, 'HandleVisibility', 'off');
        ylabel('$u$')
        ylim([u_min, u_max]);
        xlim([0, L_sim*Ts])
        legend()
    end

    if exist("r") && ~all(r(t) == 0)
        subplot(plot_dims(1),plot_dims(2),obsv_state)
        stairs(t, r(t),'black--', LineWidth=0.75, DisplayName='r')

    end

    if exist("tit")
        sgtitle(tit);
    end

    
end