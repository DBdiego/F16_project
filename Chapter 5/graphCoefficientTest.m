load('FindF16Dynaimcs_workspace.mat')


%% original transfer function response to negative step (NOT ZOOMED IN)
%%
if 0
    an_neg_de_step = figure('Position', [300 300 1000 300]);

    step(H_an_de, 400)
    grid on;

    % Exporting figure as pdf
    set(an_neg_de_step,'Units','Inches');
    pos = get(an_neg_de_step,'Position');
    set(an_neg_de_step,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(an_neg_de_step,'./plots/an_neg_de_step.pdf','-dpdf','-r0')
end

%% original transfer function response to negative step (ZOOMED IN)
%%
if 0
    an_neg_de_step_zoom = figure('Position', [300 300 1000 300]);

    step(H_an_de, 5)
    grid on;


    % Exporting figure as pdf
    set(an_neg_de_step_zoom,'Units','Inches');
    pos = get(an_neg_de_step_zoom,'Position');
    set(an_neg_de_step_zoom,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(an_neg_de_step_zoom,'./plots/an_neg_de_step_zoom.pdf','-dpdf','-r0')

end

%% original transfer function poles & zeros
%%
if 0
    pole_zero_de_an = figure;%('Position', [300 300 1000 400]);

    pzmap(minreal(H_an_de));
    grid on;


    % Exporting figure as pdf
    set(pole_zero_de_an,'Units','Inches');
    pos = get(pole_zero_de_an,'Position');
    set(pole_zero_de_an,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(pole_zero_de_an,'./plots/pole_zero_de_an.pdf','-dpdf','-r0')

end


%% poles & zeros of all a_x positions
%%
if 0
    
    x_a_labels = [0, 5, 5.9, 6, 7, 15];
    
    zeros_array = [];
    index = 1;
    pole_zero_de_an_ax = figure;
    for x_a_label = x_a_labels
        load(strcat('./FindF16Dynaimcs_workspace_', num2str(x_a_label), '.mat'))

        disp(strcat('x_a position: ', num2str(x_a_label), 'ft.'))
        H = zpk(minreal(H_an_de))
        disp('');


        zeros_array =[zeros_array, zero(minreal(H_an_de))];
        pzmap(minreal(H_an_de));
        hold on


        index = index + 1;
    end

    legend('0 ft.', '5 ft.', '5.9 ft.', '6 ft.', '7 ft.','15 ft.', 'Location', 'southeast');
    grid ('on', 'GridAlpha', 0.4);

    hold off;

    % Exporting figure as pdf
    set(pole_zero_de_an_ax,'Units','Inches');
    pos = get(pole_zero_de_an_ax,'Position');
    set(pole_zero_de_an_ax,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(pole_zero_de_an_ax,'./plots/pole_zero_de_an_ax.pdf','-dpdf','-r0')

    % zoom version
    xlim([-320 80])
    print(pole_zero_de_an_ax,'./plots/pole_zero_de_an_ax_zoom.pdf','-dpdf','-r0')

end


%% Step responses of Tranfser Functions of a_x positions
%%
colors = ['r','b','g','k','y','m'] ;
if 1
    
    simulation_time = 0.2;
    index = 1;
    an_neg_de_step_xa = figure;
    for x_a_label = x_a_labels
        load(strcat('./FindF16Dynaimcs_workspace_', num2str(x_a_label), '.mat'))

        
        step(H_an_de, simulation_time);
        hold on


        index = index + 1;
    end

    legend('0 ft.', '5 ft.', '5.9 ft.', '6 ft.', '7 ft.','15 ft.', 'Location', 'southeast');
    grid on;

    % Exporting figure as pdf
    set(an_neg_de_step_xa,'Units','Inches');
    pos = get(an_neg_de_step_xa,'Position');
    set(an_neg_de_step_xa,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(an_neg_de_step_xa,'./plots/an_neg_de_step_xa.pdf','-dpdf','-r0')

end
