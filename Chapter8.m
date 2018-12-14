%% Creating the 5 states reduced model
A_R = [A_longitude_lo(1,1:5) ; A_longitude_lo(2,1:5) ; A_longitude_lo(3,1:5) ; A_longitude_lo(4,1:5) ; A_longitude_lo(5,1:5)];
B_R = [A_longitude_lo(1,  7) ; A_longitude_lo(2,  7) ; A_longitude_lo(3,  7) ; A_longitude_lo(4,  7) ; A_longitude_lo(5,  7)];
C_R = [C_longitude_lo(1,1:5) ; C_longitude_lo(2,1:5) ; C_longitude_lo(3,1:5) ; C_longitude_lo(4,1:5) ; C_longitude_lo(5,1:5)]; 
D_R = zeros(5,1);

ss_R = ss(A_R, B_R, C_R, D_R);

H_R = tf(ss_R);

a = 20.2;
s = tf('s');
H_act = a/(s+a); %actuator dynamics

H_R_act = H_R * H_act;

