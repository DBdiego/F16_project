load('FindF16Dynaimcs_workspace.mat')

x_a_labels = [0, 5, 5.9, 6, 7, 15];
colors = ['b','r','g','k','y','m'] ;


%% 
%%
zeros_array = [];
index = 1;
figure(1);
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
grid on;

hold off;


%% Step on Tranfser Functios
%%

simulation_time = 0.2;
index = 1;
figure(2);
for x_a_label = x_a_labels
    load(strcat('./FindF16Dynaimcs_workspace_', num2str(x_a_label), '.mat'))
    

    step(H_an_de, simulation_time);
    hold on
    
    
    index = index + 1;
end

legend('0 ft.', '5 ft.', '5.9 ft.', '6 ft.', '7 ft.','15 ft.', 'Location', 'southeast');
grid on;

