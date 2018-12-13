%% Creating Matrices for OL analysis
A_OL_long = [A_longitude_lo(2,2:5) ; A_longitude_lo(3,2:5) ; A_longitude_lo(4,2:5) ; A_longitude_lo(5,2:5)];
B_OL_long = [A_longitude_lo(2,7) ; A_longitude_lo(3,7) ; A_longitude_lo(4,7) ; A_longitude_lo(5,7)];
C_OL_long = [C_longitude_lo(2,2:5) ; C_longitude_lo(3,2:5) ; C_longitude_lo(4,2:5) ; C_longitude_lo(5,2:5)]; 
D_OL_long = zeros(4,1);

A_OL_lat = [A_lateral_lo(1,1) A_lateral_lo(1,4:6) ; A_lateral_lo(4,1) A_lateral_lo(4,4:6) ; A_lateral_lo(5,1) A_lateral_lo(5,4:6) ; A_lateral_lo(6,1) A_lateral_lo(6,4:6)];
B_OL_lat = [A_lateral_lo(1,8) A_lateral_lo(1,9); A_lateral_lo(4,8) A_lateral_lo(4,9); A_lateral_lo(5,8) A_lateral_lo(5,9); A_lateral_lo(6,8)  A_lateral_lo(6,9)];
C_OL_lat = [C_lateral_lo(1,1) C_lateral_lo(1,4:6) ; C_lateral_lo(4,1) C_lateral_lo(4,4:6) ; C_lateral_lo(5,1) C_lateral_lo(5,4:6) ; C_lateral_lo(6,1) C_lateral_lo(6,4:6)];
D_OL_lat = zeros(4,2);

ss_long = ss(A_OL_long, B_OL_long, C_OL_long, D_OL_long);
ss_lat = ss(A_OL_lat, B_OL_lat, C_OL_lat, D_OL_lat);

H_long = tf(ss_long);
H_lat = tf(ss_lat);

a = 20.2;
s = tf('s');
H_act = a/(s+a); %actuator dynamics

H_long_actuator = H_long * H_act;
H_lat_actuator = H_lat * H_act;

%% Comparison with original model
figure %elevator impulse responses
T_sym = 5;
dt_sym=0.01;
d_sym=1;
m_sym=rad2deg(0.025);
t_sym=0:dt_sym:T_sym;
N_s = length(t_sym);
N_sym = d_sym/dt_sym;
sig_sym=[m_sym*ones(1,N_sym), zeros(1,N_s-N_sym)];

y1_sym = lsim(H_long_actuator(1),sig_sym',t_sym);
plot(t_sym,y1_sym)
hold on 
y1_sym_sys = lsim(all_tfs(5,2),sig_sym',t_sym);
plot(t_sym,y1_sym_sys,'+')
hold on

y2_sym = lsim(H_long_actuator(2),sig_sym',t_sym);
plot(t_sym,y2_sym)
hold on 
y2_sym_sys = lsim(all_tfs(7,2),sig_sym',t_sym);
plot(t_sym,y2_sym_sys,'+')
hold on

y3_sym = lsim(H_long_actuator(3),sig_sym',t_sym);
plot(t_sym,y3_sym)
hold on 
y3_sym_sys = lsim(all_tfs(8,2),sig_sym',t_sym);
plot(t_sym,y3_sym_sys,'+')
hold on

y4_sym = lsim(H_long_actuator(4),sig_sym',t_sym);
plot(t_sym,y4_sym)
hold on 
y4_sym_sys = lsim(all_tfs(11,2),sig_sym',t_sym);
plot(t_sym,y4_sym_sys,'+')
hold on

title('Responses to a .5s elevator impluse of amplitude +0.025rad')
legend('\theta', '\theta original', 'v', 'v original', '\alpha', '\alpha original','q', 'q original')

clear y1_sym y1_sym_sys y2_sym y2_sym_sys y3_sym y3_sym_sys y4_sym y4_sym_sys

figure %aileron impulse responses
y1_sym = lsim(H_lat_actuator(1,1),sig_sym',t_sym);
plot(t_sym,y1_sym)
hold on 
y1_sym_sys = lsim(all_tfs(4,3),sig_sym',t_sym);
plot(t_sym,y1_sym_sys,'+')
hold on

y2_sym = lsim(H_lat_actuator(2,1),sig_sym',t_sym);
plot(t_sym,y2_sym)
hold on 
y2_sym_sys = lsim(all_tfs(9,3),sig_sym',t_sym);
plot(t_sym,y2_sym_sys,'+')
hold on

y3_sym = lsim(H_lat_actuator(3,1),sig_sym',t_sym);
plot(t_sym,y3_sym)
hold on 
y3_sym_sys = lsim(all_tfs(10,3),sig_sym',t_sym);
plot(t_sym,y3_sym_sys,'+')
hold on

y4_sym = lsim(H_lat_actuator(4,1),sig_sym',t_sym);
plot(t_sym,y4_sym)
hold on 
y4_sym_sys = lsim(all_tfs(12,3),sig_sym',t_sym);
plot(t_sym,y4_sym_sys,'+')
hold on

title('Responses to a .5s aileron impluse of amplitude +0.025rad')
legend('\phi', '\phi original', '\beta', '\beta original', 'p', 'p original','r', 'r original')


clear y1_sym y1_sym_sys y2_sym y2_sym_sys y3_sym y3_sym_sys y4_sym y4_sym_sys

figure %rudder impulse responses
y1_sym = lsim(H_lat_actuator(1,2),sig_sym',t_sym);
plot(t_sym,y1_sym)
hold on 
y1_sym_sys = lsim(all_tfs(4,4),sig_sym',t_sym);
plot(t_sym,y1_sym_sys,'+')
hold on

y2_sym = lsim(H_lat_actuator(2,2),sig_sym',t_sym);
plot(t_sym,y2_sym)
hold on 
y2_sym_sys = lsim(all_tfs(9,4),sig_sym',t_sym);
plot(t_sym,y2_sym_sys,'+')
hold on

y3_sym = lsim(H_lat_actuator(3,2),sig_sym',t_sym);
plot(t_sym,y3_sym)
hold on 
y3_sym_sys = lsim(all_tfs(10,4),sig_sym',t_sym);
plot(t_sym,y3_sym_sys,'+')
hold on

y4_sym = lsim(H_lat_actuator(4,2),sig_sym',t_sym);
plot(t_sym,y4_sym)
hold on 
y4_sym_sys = lsim(all_tfs(12,4),sig_sym',t_sym);
plot(t_sym,y4_sym_sys,'+')
hold on

title('Responses to a .5s rudder impluse of amplitude +0.025rad')
legend('\phi', '\phi original', '\beta', '\beta original', 'p', 'p original','r', 'r original')

%% Pole zero maps
figure 
pzmap(minreal(H_lat))
title('pole zero map for H_{lat}')

figure
pzmap(minreal(H_long))
title('pole zero map for H_{long}')

pole(minreal(H_long))
pole(minreal(H_lat))

%% Step responses
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

%%Pulse shape input parametters
T=10;
dt=0.01;
d_dev=1;
m_dev=rad2deg(0.025);
t=0:dt:T;
N = length(t);
N_dev = d_dev/dt;
sig_dev=[m_dev*ones(1,N_dev), zeros(1,N-N_dev)];
sig_fin=[sig_dev'];


%% Dutch Roll
%%
figure('name', 'dutch roll')
subplot(2,1,1)
y_drp = lsim(H_lat(3,2),sig_fin,t);
plot(t,y_drp)
title('')
ylabel('p [deg/s]')
subplot(2,1,2)
y_drr = lsim(H_lat(4,2),sig_fin,t);
plot(t,y_drr)
title('')
ylabel('r [deg/s]')


%% Spiral
%%
figure('name', 'Spiral')
y_sp = lsim(H_lat(1,1),sig_fin,t);
plot(t,y_sp)
title('')
ylabel('\phi [deg]')


%% Aperiodic Roll
%%
figure('name', 'Aperiodic Roll')
y_ar = lsim(H_lat(3,1),sig_fin,t);
plot(t,y_ar)
title('')
ylabel('p [deg/s]')


%% Short Period Reduced model
A_SP = [A_OL_long(3,3:4) ; A_OL_long(4,3:4)];
B_SP = [B_OL_long(3:4)];
C_SP = [C_OL_long(3,3:4) ; C_OL_long(4,3:4)];
D_SP = zeros(2,1);

ss_SP_reduced = ss(A_SP, B_SP, C_SP, D_SP);
H_SP_reduced = tf(ss_SP_reduced);

H_SP_init = H_SP_reduced(2);


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

%% Gibson 
omega_sp = 0.03*velocity*0.3048;
damp_rat = 0.5;
T_theta2 = 1/(0.75*omega_sp);

%desired poles
p1 = -damp_rat*omega_sp + j*omega_sp*sqrt(1-damp_rat^2);
p2 = -damp_rat*omega_sp - j*omega_sp*sqrt(1-damp_rat^2);

K = place(A_SP,B_SP,[p1 p2]);

d_alpha = atan(4.572/(velocity*0.3048));
d_delta_e = K(1)*d_alpha; %Ok between -25 and +25 deg

H_SP_new = tf(ss(A_SP-B_SP*K, B_SP, C_SP , D_SP)); 

H_SP_new_q = H_SP_new(2); %Transfer function with new poles q/delta_e


Kf=1;
F = Kf*(1+T_theta2*s)/(s - zero(H_SP_new_q)); %Pre-filter to compensqte for q/delta_e zero and place the desired zero

H_gibs = minreal(H_SP_new_q*F); %final sys only gain remains to be fixed

Kf = -2.73884; % gain that permits to have no steady state error

H_gibs_new = Kf*H_gibs; % 

%% Gibson criterion
CAP = (gD*omega_sp^2*T_theta2)/velocity;

step_info = stepinfo(H_gibs_new);
qm_qs = step_info.Peak/1 ; %q_max / q_steady

%response of the sys to a 10s ramp
t_ramp=0:0.01:10;
alpha=1;
ramp=alpha*t;
input_ramp = [ramp 10*ones(1,length(t_ramp)-1)];
t_sig_ramp = 0:0.01:20;
[y_ramp,t_ramp]=lsim(H_gibs_new,input_ramp,t_sig_ramp);
plot(t_ramp,y_ramp)
hold on
plot(t_ramp,input_ramp)
title('')

DB = findpeaks(y_ramp,'NPeaks',1) - 10; %Drop Back value
DB_qs = DB / 10; %%% = T_theta2 - (2*damp_rat/omega_sp); ?
