

A_OL_long = [A_longitude_lo(2,2:5) ; A_longitude_lo(3,2:5) ; A_longitude_lo(4,2:5) ; A_longitude_lo(5,2:5)];
B_OL_long = [A_longitude_lo(2,  7) ; A_longitude_lo(3,  7) ; A_longitude_lo(4,  7) ; A_longitude_lo(5,  7)];
C_OL_long = [C_longitude_lo(2,2:5) ; C_longitude_lo(3,2:5) ; C_longitude_lo(4,2:5) ; C_longitude_lo(5,2:5)]; 
D_OL_long = zeros(4,1);

A_OL_lat = [A_lateral_lo(1,1) A_lateral_lo(1,4:6) ; A_lateral_lo(4,1) A_lateral_lo(4,4:6) ; A_lateral_lo(5,1) A_lateral_lo(5,4:6) ; A_lateral_lo(6,1) A_lateral_lo(6,4:6)];
B_OL_lat = [A_lateral_lo(1,8) A_lateral_lo(1,  9) ; A_lateral_lo(4,8) A_lateral_lo(4,  9) ; A_lateral_lo(5,8) A_lateral_lo(5,  9) ; A_lateral_lo(6,8) A_lateral_lo(6,  9)];
C_OL_lat = [C_lateral_lo(1,1) C_lateral_lo(1,4:6) ; C_lateral_lo(4,1) C_lateral_lo(4,4:6) ; C_lateral_lo(5,1) C_lateral_lo(5,4:6) ; C_lateral_lo(6,1) C_lateral_lo(6,4:6)];
D_OL_lat = zeros(4,2);

ss_long = ss(A_OL_long, B_OL_long, C_OL_long, D_OL_long);
ss_lat  = ss(A_OL_lat , B_OL_lat , C_OL_lat , D_OL_lat );

H_long = tf(ss_long);
H_lat  = tf(ss_lat );

a = 20.2;
s = tf('s');
H_act = a/(s+a);

H_long_actuator = H_long * H_act;
H_lat_actuator  = H_lat  * H_act;


eigenmotion_plots = 0;
four_vs_two_plots = 0;
steady_state_error_plot = 0;

if eigenmotion_plots
    %% Comparison with original model
    figure %elevator impulse responses
    [y1, te] = impulse(H_long_actuator(1), 5);
    y1_sys = impulse(all_tfs(5,2), 5);
    plot(te, y1)
    hold on
    plot(te, y1_sys, '+')
    hold on

    [y2, te] = impulse(H_long_actuator(2), 5);
    y2_sys = impulse(all_tfs(7,2), 5);
    plot(te, y2)
    hold on
    plot(te, y2_sys, '+')
    hold on

    [y3, te] = impulse(H_long_actuator(3), 5);
    y3_sys = impulse(all_tfs(8,2), 5);
    plot(te, y3)
    hold on
    plot(te, y3_sys, '+')
    hold on

    [y4, te] = impulse(H_long_actuator(4), 5);
    y4_sys = impulse(all_tfs(11,2), 5);
    plot(te, y4)
    hold on
    plot(te, y4_sys, '+')
    hold on
    title('Responses to an elevator impluse')
    legend('\theta', '\theta original', 'v', 'v original', '\alpha', '\alpha original','q', 'q original')

    clear te y1 y1_sys y2 y2_sys y3 y3_sys y4 y4_sys

    figure %aileron impulse responses
    [y1, te] = impulse(H_lat_actuator(1,1), 5);
    y1_sys = impulse(all_tfs(4,3), 5);
    plot(te, y1)
    hold on
    plot(te, y1_sys, '+')
    hold on

    [y2, te] = impulse(H_lat_actuator(2,1), 5);
    y2_sys = impulse(all_tfs(9,3), 5);
    plot(te, y2)
    hold on
    plot(te, y2_sys, '+')
    hold on

    [y3, te] = impulse(H_lat_actuator(3,1), 5);
    y3_sys = impulse(all_tfs(10,3), 5);
    plot(te, y3)
    hold on
    plot(te, y3_sys, '+')
    hold on

    [y4, te] = impulse(H_lat_actuator(4,1), 5);
    y4_sys = impulse(all_tfs(12,3), 5);
    plot(te, y4)
    hold on
    plot(te, y4_sys, '+')
    hold on
    title('Responses to an aileron impluse')
    legend('\phi', '\phi original', '\beta', '\beta original', 'p', 'p original','r', 'r original')


    clear te y1 y1_sys y2 y2_sys y3 y3_sys y4 y4_sys

    figure %rudder impulse responses
    [y1, te] = impulse(H_lat_actuator(1,2), 5);
    y1_sys = impulse(all_tfs(4,4), 5);
    plot(te, y1)
    hold on
    plot(te, y1_sys, '+')
    hold on

    [y2, te] = impulse(H_lat_actuator(2,2), 5);
    y2_sys = impulse(all_tfs(9,4), 5);
    plot(te, y2)
    hold on
    plot(te, y2_sys, '+')
    hold on

    [y3, te] = impulse(H_lat_actuator(3,2), 5);
    y3_sys = impulse(all_tfs(10,4), 5);
    plot(te, y3)
    hold on
    plot(te, y3_sys, '+')
    hold on

    [y4, te] = impulse(H_lat_actuator(4,2), 5);
    y4_sys = impulse(all_tfs(12,4), 5);
    plot(te, y4)
    hold on
    plot(te, y4_sys, '+')
    hold on
    title('Responses to a rudder impluse')
    legend('\phi', '\phi original', '\beta', '\beta original', 'p', 'p original','r', 'r original')

    %% Pole zero maps
    figure 
    pzmap(minreal(H_lat))
    title('pole zero map H_{lat}')

    figure
    pzmap(minreal(H_long))
    title('pole zero map H_{long}')

    pole(minreal(H_long))
    pole(minreal(H_lat))

    %% Step responses (not working very well)
    %% Short Period
    %%
    opt = stepDataOptions('StepAmplitude',-0.005);
    figure('name', 'short period')
    step(H_long(4),10,opt) %short period
    title('')
    ylabel('q [deg/s]')

    %% Phugoid
    %%
    opt = stepDataOptions('StepAmplitude',-0.005);
    figure('name', 'phugoid')
    subplot(2,1,1)
    step(H_long(2),300,opt) %phugoid
    title('')
    ylabel('V [ft/s]')
    subplot(2,1,2)
    step(H_long(3),300,opt) %phugoid
    title('')
    ylabel('\alpha [deg]')


    %% Dutch Roll
    %%
    opt = stepDataOptions('StepAmplitude',+0.025);
    figure('name', 'dutch roll')
    subplot(2,1,1)
    impulse(H_lat(3,2),10)%,opt) %dutch roll
    title('')
    ylabel('p [deg/s]')
    subplot(2,1,2)
    impulse(H_lat(4,2),10)%,opt) %dutch roll
    title('')
    ylabel('r [deg/s]')


    %% Spiral
    %%
    figure('name', 'Spiral')
    impulse(H_lat(1,1),100)%,opt) %Spiral
    title('')
    ylabel('\phi [deg]')


    %% Aperiodic Roll
    %%
    figure('name', 'Aperiodic Roll')
    %step(H_lat(3,1),10,opt)
    %impulse(H_lat(3,1),10,opt)
    T=10;
    dt=0.01;
    d_rud=1;
    m_rud=rad2deg(0.025);
    t=0:dt:T;
    N = length(t);
    N_rud = d_rud/dt;
    sig_rud=[m_rud*ones(1,N_rud), zeros(1,N-N_rud)];
    sig_null = zeros(1,N);
    sig_dr=[sig_rud'];
    %[A1 B1 C1 D1] = tf2ss(cell2mat(H_lat(3,1).Numerator), cell2mat(H_lat(3,1).Denominator));
    %y_dr = lsim([A1 B1 C1 D1],sig_dr,t);
    y_dr = lsim(H_lat(3,1),sig_dr,t);
    plot(t,y_dr)
    title('')
    ylabel('p [deg/s]')

    
end





%% Short Period Reduced model
A_SP = [A_OL_long(3,3:4) ; A_OL_long(4,3:4)];
B_SP = [B_OL_long(3:4)];
C_SP = [C_OL_long(3,3:4) ; C_OL_long(4,3:4)];
D_SP = zeros(2,1);

ss_SP_reduced = ss(A_SP, B_SP, C_SP, D_SP);
H_SP_reduced = tf(ss_SP_reduced);

H_SP_init = H_SP_reduced(2);


if four_vs_two_plots
    figure
    step(H_SP_reduced(2),30)
    hold on
    step(H_long(4),35)
    legend('Reduced 2 states model' , '4 states model')
    title('Short Period Reduced model')

    figure
    step(H_SP_reduced(2),10)
    hold on
    step(H_long(4),10)
    legend('Reduced 2 states model' , '4 states model')
    title('Short Period Reduced model')

end 


%% Gibson 
omega_sp = 0.03*velocity*0.3048;
damp_rat = 0.5;
T_theta2 = 1/(0.75*omega_sp);

[numerator, denominator] = tfdata(H_SP_init);
T_theta2_O = numerator{1}(2)/numerator{1}(3);

%desired poles
p1 = -damp_rat*omega_sp + j*omega_sp*sqrt(1-damp_rat^2);
p2 = -damp_rat*omega_sp - j*omega_sp*sqrt(1-damp_rat^2);


% Placing poles
K = place(A_SP,B_SP,[p1 p2]);

% induced angles
d_alpha = atan(4.572/(velocity*0.3048));
d_delta_e = K(1)*d_alpha; %Ok between -25 and +25 deg

% New short period transfer function
H_SP_new = tf(ss(A_SP-B_SP*K, B_SP, C_SP , D_SP));

% Transfer function of q
H_SP_new_q = H_SP_new(2);

% Finding Filter Gain (USE SISO TOOOOOOOOL)
%sisotool(H_SP_new_q)
Kf = -5.3885;  %-0.665452; %-0.663;%1;
F  =  Kf * (T_theta2   * s + 1)/...
           (T_theta2_O * s + 1);  %- zero(H_SP_new_q));


H_gibs = minreal(H_SP_new_q*F);

% Step response of original and filtered 
if steady_state_error_plot
    figure('name', 'Steady state Error')
    step(H_gibs, 2)
    hold on
    step(H_SP_new_q, 2)
    title('')
    ylabel('q [deg/s]')
end



%% Gibson criterion
CAP = (gD*omega_sp^2*T_theta2)/velocity

step_info = stepinfo(H_gibs);
qm_qs = step_info.Peak/1  %q_max / q_steady


%response of the sys to a 10s ramp
figure
t_ramp=0:0.01:10;
alpha=1;
ramp=alpha*t_ramp;
input_ramp = [ramp 10*ones(1,length(t_ramp)-1)];
t_sig_ramp = 0:0.01:20;
[y_ramp,t_ramp]=lsim(H_gibs,input_ramp,t_sig_ramp);
plot(t_ramp,y_ramp)
hold on
plot(t_ramp,input_ramp)
title('')

DB = findpeaks(y_ramp,'NPeaks',1) - 10; %Drop Back value
DB_qs = DB / 10 %%% = T_theta2 - (2*damp_rat/omega_sp); ?
DB_qs_des = T_theta2 - (2*damp_rat/omega_sp)









