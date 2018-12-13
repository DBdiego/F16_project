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
H_act = a/(s+a);

H_long_actuator = H_long * H_act;
H_lat_actuator = H_lat * H_act;

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

figure
pzmap(minreal(H_long))

pole(minreal(H_long))
pole(minreal(H_lat))

%% Step responses (not working very well)
opt = stepDataOptions('StepAmplitude',-0.005);
figure
step(H_long(4),10,opt) %short period
figure
subplot(2,1,1)
step(H_long(2),10,opt) %phugoid
subplot(2,1,2)
step(H_long(3),10,opt) %phugoid

opt = stepDataOptions('StepAmplitude',+0.025);
figure
subplot(2,1,1)
step(H_lat(3,2),10,opt) %dutch roll
subplot(2,1,2)
step(H_lat(4,2),10,opt) %dutch roll

figure
step(H_lat(1,1),10,opt) %Spiral

figure
step(H_lat(3,1),10,opt) %Aperiodic Roll

%% Short Period Reduced model
A_SP = [A_OL_long(3,3:4) ; A_OL_long(4,3:4)];
B_SP = [B_OL_long(3:4)];
C_SP = [C_OL_long(3,3:4) ; C_OL_long(4,3:4)];
D_SP = zeros(2,1);

ss_SP_reduced = ss(A_SP, B_SP, C_SP, D_SP);
H_SP_reduced = tf(ss_SP_reduced);

figure
step(H_SP_reduced(2),30)
hold on
step(H_long(4),30)
legend('Reduced 2 states model' , '4 states model')