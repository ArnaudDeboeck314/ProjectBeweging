function [t_2,dt_2,ddt_2] = Voorwaartse_kinemtica(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,t,Ts,M2,t_init,t1)

F1_2_x = zeros(size(t));
F1_2_y = zeros(size(t));
F2_12_x = zeros(size(t));
F2_12_y = zeros(size(t));
F2_3_x = zeros(size(t));
F2_3_y = zeros(size(t));
F3_4_x = zeros(size(t));
F3_4_y = zeros(size(t));
F3_5_x = zeros(size(t));
F3_5_y = zeros(size(t));
F4_6_x = zeros(size(t));
F4_6_y = zeros(size(t));
F6_7_x = zeros(size(t));
F6_7_y = zeros(size(t));
F6_8_x = zeros(size(t));
F6_8_y = zeros(size(t));
F8_9_x = zeros(size(t));
F8_9_y = zeros(size(t));
F8_10_x = zeros(size(t));
F8_10_y = zeros(size(t));
F12_10_x = zeros(size(t));
F12_10_y = zeros(size(t));
F12_11_x = zeros(size(t));
F12_11_y = zeros(size(t));
F1_5_x = zeros(size(t));
F1_5_y = zeros(size(t));
F1_7_x = zeros(size(t));
F1_7_y = zeros(size(t));
F1_9_x = zeros(size(t));
F1_9_y = zeros(size(t));
F1_11_x = zeros(size(t));
F1_11_y = zeros(size(t));

ddt_2 = zeros(size(t));
ddt_3 = zeros(size(t));
ddt_4 = zeros(size(t));
ddt_5 = zeros(size(t));
ddt_6 = zeros(size(t));
ddt_7 = zeros(size(t));
ddt_8 = zeros(size(t));
ddt_9 = zeros(size(t));
ddt_10 = zeros(size(t));
ddt_11 = zeros(size(t));
ddt_12 = zeros(size(t));

dt_2 = zeros(size(t));
dt_3 = zeros(size(t));
dt_4 = zeros(size(t));
dt_5 = zeros(size(t));
dt_6 = zeros(size(t));
dt_7 = zeros(size(t));
dt_8 = zeros(size(t));
dt_9 = zeros(size(t));
dt_10 = zeros(size(t));
dt_11 = zeros(size(t));
dt_12 = zeros(size(t));

t_2 = zeros(size(t));
t_3 = zeros(size(t));
t_4 = zeros(size(t));
t_5 = zeros(size(t));
t_6 = zeros(size(t));
t_7 = zeros(size(t));
t_8 = zeros(size(t));
t_9 = zeros(size(t));
t_10 = zeros(size(t));
t_11 = zeros(size(t));
t_12 = zeros(size(t));

t_2(1) = 1.4;

dt_2(1) = 0;

optim_options = optimset('Display','off');

for k=1:size(t)
    
    [x, ~, exitflag] = fsolve('Loop_closure_ghassaei',t_init',optim_options,t_2(k),r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,t1);
    
    if (exitflag ~= 1)
        disp 'The fsole exit flag was not 1, probably no convergence!'
    end
    
    t_3(k) = x(1);
    t_4(k) = x(2);
    t_5(k) = x(3);
    t_6(k) = x(4);
    t_7(k) = x(5);
    t_8(k) = x(6);
    t_9(k) = x(7);
    t_10(k) = x(8);
    t_11(k) = x(9);
    t_12(k) = x(10);
    
    A = [0, 0, 0, 0, 0, 0, 0, 0, r11*sin(t_11(k)), -r12*sin(t_12(k));
         0, 0, 0, 0, 0, 0, 0, 0, -r11*cos(t_11(k)), r12*cos(t_12(k));
         0, 0, 0, 0, -r7*sin(t_7(k)), r8*sin(t_8(k)), -r9*sin(t_9(k)), 0, 0, 0;
         0, 0, 0, 0, r7*cos(t_7(k)), -r8*cos(t_8(k)), r9*cos(t_9(k)), 0, 0, 0;
         0, -r4*sin(t_4(k)), -r5*sin(t_5(k)), r6*sin(t_6(k)), r7*sin(t_7(k)), 0, 0, 0, 0, 0;
         0, r4*cos(t_4(k)), r5*cos(t_5(k)), -r6*cos(t_6(k)), -r7*cos(t_7(k)), 0, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 0, -r9*sin(t_9(k)), -r10*sin(t_10(k)), r11*sin(t_11(k)), 0;
         0, 0, 0, 0, 0, 0, r9*cos(t_9(k)), r10*cos(t_10(k)), -r11*cos(t_11(k)), 0;
         r3*sin(t_3(k)), 0, -r5*sin(t_5(k)), 0, 0, 0, 0, 0, 0, 0;
         -r3*cos(t_3(k)), 0, r5*cos(t_5(k)), 0, 0, 0, 0, 0, 0, 0];
        
    B = [r2*sin(t_2(k))*dt_2(k);
         -r2*cos(t_2(k))*dt_2(k);
         0;
         0;
         0;
         0;
         0;
         0;
         r2*sin(t_2(k))*dt_2(k);
         -r2*cos(t_2(k))*dt_2(k)];
    
    x = A\B;
    
    dt_3(k) = x(1);
    dt_4(k) = x(2);
    dt_5(k) = x(3);
    dt_6(k) = x(4);
    dt_7(k) = x(5);
    dt_8(k) = x(6);
    dt_9(k) = x(7);
    dt_10(k) = x(8);
    dt_11(k) = x(9);
    dt_12(k) = x(10);
    
    
    A = [1 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 m2*r2*0.5*sin(t_2(k)) 0 0 0 0 0 0 0 0 0 0;
         0 1 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -m2*r2*0.5*cos(t_2(k)) 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 1 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 m3*r3*0.5*sin(t_3(k)) 0 -m5*r5*sin(t_5(k)) 0 0 0 0 0 0 0;
         0 0 0 0 0 1 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -m3*r3*0.5*cos(t_3(k)) 0 m5*r5*cos(t_5(k)) 0 0 0 0 0 0 0;
         0 0 0 0 0 0 1 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -m4*r4*0.5*sin(t_4(k)) -m5*r5*sin(t_5(k)) 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 1 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 m4*r4*0.5*cos(t_4(k)) m5*r5*cos(t_5(k)) 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 -m5*r5*0.5*sin(t_5(k)) 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 m5*r5*0.5*cos(t_5(k)) 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 1 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -m6*r6*0.5*sin(t_6(k)) -m7*r7*sin(t_7(k)) 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 1 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 m6*r6*0.5*cos(t_6(k)) m7*r7*cos(t_7(k)) 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 -m7*r7*0.5*sin(t_7(k)) 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 m7*r7*0.5*cos(t_7(k)) 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -m7*r7*sin(t_7(k)) m8*r8*0.5*sin(t_8(k)) 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 m7*r7*cos(t_7(k)) -m8*r8*0.5*cos(t_8(k)) 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 m9*r9*sin(t_9(k)) m10*r10*0.5*sin(t_10(k)) 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 -m9*r9*cos(t_9(k)) -m10*r10*0.5*cos(t_10(k)) 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 m9*r9*0.5*sin(t_9(k)) 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -m9*r9*0.5*cos(t_9(k)) 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 m11*r11*0.5*sin(t_11(k)) 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 m11*r11*0.5*cos(t_11(k)) 0;
         0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 -1 0 0 0 0 0 0 0 0 0 m2*r2*sin(t_2(k)) 0 0 0 0 0 0 0 0 0 m12*r12*0.5*sin(t_12(k));
         0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 -1 0 0 0 0 0 0 0 0 -m2*r2*cos(t_2(k)) 0 0 0 0 0 0 0 0 0 -m12*r12*0.5*cos(t_12(k));
         0 0 r2*sin(t_2(k)) -r2*cos(t_2(k)) r2*sin(t_2(k)) -r2*cos(t_2(k)) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -I2 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0.5*r3*-sin(t_3(k)) 0.5*r3*cos(t_3(k)) 0.5*r3*sin(t_3(k)+pi) 0.5*r3*-cos(t_3(k)+pi) 0.5*r3*sin(t_3(k)+pi) 0.5*r3*-cos(t_3(k)+pi) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -I3 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0.5*r4*-sin(t_4(k)) 0.5*r4*cos(t_4(k)) 0 0 0.5*r4*sin(t_4(k)+pi) 0.5*r4*-cos(t_4(k)+pi) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -I4 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0.5*r5*-sin(t_5(k)+pi) 0.5*r5*cos(t_5(k)+pi) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5*r5*-sin(t_5(k)) 0.5*r5*cos(t_5(k)) 0 0 0 0 0 0 0 0 0 -I5 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0.5*r6*-sin(t_6(k)+pi) 0.5*r6*cos(t_6(k)+pi) 0.5*r6*sin(t_6(k)) 0.5*r6*-cos(t_6(k)) 0.5*r6*sin(t_6(k)) 0.5*r6*-cos(t_6(k)) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -I6 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0.5*r7*-sin(t_7(k)+pi) 0.5*r7*cos(t_7(k)+pi) 0 0 0 0 0 0 0 0 0 0 0 0 0.5*r7*-sin(t_7(k)) 0.5*r7*cos(t_7(k)) 0 0 0 0 0 0 0 0 0 -I7 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5*r8*-sin(t_8(k)+pi) 0.5*r8*cos(t_8(k)+pi) 0.5*r8*sin(t_8(k)) 0.5*r8*-cos(t_8(k)) 0.5*r8*sin(t_8(k)) 0.5*r8*-cos(t_8(k)) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -I8 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5*r9*-sin(t_9(k)) 0.5*r9*cos(t_9(k)) 0 0 0 0 0 0 0 0 0 0 0.5*r9*-sin(t_9(k)+pi) 0.5*r9*cos(t_9(k)+pi) 0 0 0 0 0 0 0 0 0 -I9 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5*r10*-sin(t_10(k)+pi) 0.5*r10*cos(t_10(k)+pi) 0.5*r10*-sin(t_10(k)) 0.5*r10*cos(t_10(k)) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -I10 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5*r11*-sin(t_11(k)) 0.5*r11*cos(t_11(k)) 0 0 0 0 0 0 0.5*r11*-sin(t_11(k)+pi) 0.5*r11*cos(t_11(k)+pi) 0 0 0 0 0 0 0 0 0 -I11 0;
         0 0 0.5*r12*-sin(t_12(k)+pi) 0.5*r12*cos(t_12(k)+pi) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5*r12*sin(t_12(k)) 0.5*r12*-cos(t_12(k)) 0.5*r12*sin(t_12(k)) 0.5*r12*-cos(t_12(k)) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -I12;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -r2*sin(t_2(k)) 0 0 0 0 0 0 0 0 r11*sin(t_11(k)) -r12*sin(t_12(k));
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 r2*cos(t_2(k)) 0 0 0 0 0 0 0 0 -r11*cos(t_11(k)) r12*cos(t_12(k));
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -r7*sin(t_7(k)) r8*sin(t_8(k)) -r9*sin(t_9(k)) 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 r7*cos(t_7(k)) -r8*cos(t_8(k)) r9*cos(t_9(k)) 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -r4*sin(t_4(k)) -r5*sin(t_5(k)) r6*sin(t_6(k)) r7*sin(t_7(k)) 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 r4*cos(t_4(k)) r5*cos(t_5(k)) -r6*cos(t_6(k)) -r7*cos(t_7(k)) 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -r9*sin(t_9(k)) -r10*sin(t_10(k)) r11*sin(t_11(k)) 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 r9*cos(t_9(k)) r10*cos(t_10(k)) -r11*cos(t_11(k)) 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -r2*sin(t_2(k)) r3*sin(t_3(k)) 0 -r5*sin(t_5(k)) 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 r2*cos(t_2(k)) -r3*cos(t_3(k)) 0 r5*cos(t_5(k)) 0 0 0 0 0 0 0;];
    
    
    B = [-dt_2(k)*dt_2(k)*r2*0.5*cos(t_2(k))*m2;
         -dt_2(k)*dt_2(k)*r2*0.5*sin(t_2(k))*m2;
         -dt_3(k)*dt_3(k)*r3*0.5*cos(t_3(k))*m3 + dt_5(k)*dt_5(k)*r5*cos(t_5(k))*m5;
         -dt_3(k)*dt_3(k)*r3*0.5*sin(t_3(k))*m3 + dt_5(k)*dt_5(k)*r5*sin(t_5(k))*m5;
         dt_4(k)*dt_4(k)*r4*0.5*cos(t_4(k))*m4 + dt_5(k)*dt_5(k)*r5*cos(t_5(k))*m5;
         dt_4(k)*dt_4(k)*r4*0.5*sin(t_4(k))*m4 + dt_5(k)*dt_5(k)*r5*sin(t_5(k))*m5;
         dt_5(k)*dt_5(k)*r5*0.5*cos(t_5(k))*m5;
         dt_5(k)*dt_5(k)*r5*0.5*sin(t_5(k))*m5;
         dt_6(k)*dt_6(k)*r6*0.5*cos(t_6(k))*m6 + dt_7(k)*dt_7(k)*r7*cos(t_7(k))*m7;
         dt_6(k)*dt_6(k)*r6*0.5*sin(t_6(k))*m6 + dt_7(k)*dt_7(k)*r7*sin(t_7(k))*m7;
         dt_7(k)*dt_7(k)*r7*0.5*cos(t_7(k))*m7;
         dt_7(k)*dt_7(k)*r7*0.5*sin(t_7(k))*m7;
         -dt_8(k)*dt_8(k)*r8*0.5*cos(t_8(k))*m8 + dt_7(k)*dt_7(k)*r7*cos(t_7(k))*m7;
         -dt_8(k)*dt_8(k)*r8*0.5*cos(t_8(k))*m8 + dt_7(k)*dt_7(k)*r7*sin(t_7(k))*m7;
         -dt_9(k)*dt_9(k)*r9*0.5*cos(t_9(k))*m9;
         -dt_9(k)*dt_9(k)*r9*0.5*sin(t_9(k))*m9;
         -dt_10(k)*dt_10(k)*r10*0.5*cos(t_10(k))*m10 + -dt_9(k)*dt_9(k)*r9*cos(t_9(k))*m9;
         -dt_10(k)*dt_10(k)*r10*0.5*sin(t_10(k))*m10 + -dt_9(k)*dt_9(k)*r9*sin(t_9(k))*m9;
         -dt_11(k)*dt_11(k)*r11*0.5*cos(t_11(k))*m11;
         -dt_11(k)*dt_11(k)*r11*0.5*sin(t_11(k))*m11;
         -dt_12(k)*dt_12(k)*r12*0.5*cos(t_12(k))*m12 + -dt_2(k)*dt_2(k)*r2*0.5*cos(t_2(k))*m2;
         -dt_12(k)*dt_12(k)*r12*0.5*sin(t_12(k))*m12 + -dt_2(k)*dt_2(k)*r2*0.5*sin(t_2(k))*m2;
         -M2(k);
         0;
         0;
         0;
         0;
         0;
         0;
         0;
         0;
         0;
         0;
         -r11*cos(t_11(k))*dt_11(k)^2+r12*cos(t_12(k))*dt_12(k)^2+r2*cos(t_2(k))*dt_2(k)^2;
         -r11*sin(t_11(k))*dt_11(k)^2+r12*sin(t_12(k))*dt_12(k)^2+r2*sin(t_2(k))*dt_2(k)^2;
         r7*cos(t_7(k))*dt_7(k)^2-r8*cos(t_8(k))*dt_8(k)^2+r9*cos(t_9(k))*dt_9(k)^2;
         r7*sin(t_7(k))*dt_7(k)^2-r8*sin(t_8(k))*dt_8(k)^2+r9*sin(t_9(k))*dt_9(k)^2;
         r4*cos(t_4(k))*dt_4(k)^2+r5*cos(t_5(k))*dt_5(k)^2-r6*cos(t_6(k))*dt_6(k)^2-r7*cos(t_7(k))*dt_7(k)^2;
         r4*sin(t_4(k))*dt_4(k)^2+r5*sin(t_5(k))*dt_5(k)^2-r6*sin(t_6(k))*dt_6(k)^2-r7*sin(t_7(k))*dt_7(k)^2;
         r9*cos(t_9(k))*dt_9(k)^2+r10*cos(t_10(k))*dt_10(k)^2-r11*cos(t_11(k))*dt_11(k)^2;
         r9*sin(t_9(k))*dt_9(k)^2+r10*sin(t_10(k))*dt_10(k)^2-r11*sin(t_11(k))*dt_11(k)^2;
         -r3*cos(t_3(k))*dt_3(k)^2+r5*cos(t_5(k))*dt_5(k)^2+r2*cos(t_2(k))*dt_2(k)^2;
         -r3*sin(t_3(k))*dt_3(k)^2+r5*sin(t_5(k))*dt_5(k)^2+r2*sin(t_2(k))*dt_2(k)^2;];
    
    x = A\B;
    
    F1_2_x(k) = x(1);
    F1_2_y(k) = x(2);
    F2_12_x(k) = x(3);
    F2_12_y(k) = x(4);
    F2_3_x(k) = x(5);
    F2_3_y(k) = x(6);
    F3_4_x(k) = x(7);
    F3_4_y(k) = x(8);
    F3_5_x(k) = x(9);
    F3_5_y(k) = x(10);
    F4_6_x(k) = x(11);
    F4_6_y(k) = x(12);
    F6_7_x(k) = x(13);
    F6_7_y(k) = x(14);
    F6_8_x(k) = x(15);
    F6_8_y(k) = x(16);
    F8_9_x(k) = x(17);
    F8_9_y(k) = x(18);
    F8_10_x(k) = x(19);
    F8_10_y(k) = x(20);
    F12_10_x(k) = x(21);
    F12_10_y(k) = x(22);
    F12_11_x(k) = x(23);
    F12_11_y(k) = x(24);
    F1_5_x(k) = x(25);
    F1_5_y(k) = x(26);
    F1_7_x(k) = x(27);
    F1_7_y(k) = x(28);
    F1_9_x(k) = x(29);
    F1_9_y(k) = x(30);
    F1_11_x(k) = x(31);
    F1_11_y(k) = x(32);
    ddt_2(k) = x(33);
    ddt_3(k) = x(34);
    ddt_4(k) = x(35);
    ddt_5(k) = x(36);
    ddt_6(k) = x(37);
    ddt_7(k) = x(38);
    ddt_8(k) = x(39);
    ddt_9(k) = x(40);
    ddt_10(k) = x(41);
    ddt_11(k) = x(42);
    ddt_12(k) = x(43);
    
    t_init = [t_3(k) + Ts*dt_3(k), t_4(k) + Ts*dt_4(k), t_5(k) + Ts*dt_5(k), t_6(k) + Ts*dt_6(k), t_7(k) + Ts*dt_7(k),
              t_8(k) + Ts*dt_8(k), t_9(k) + Ts*dt_9(k), t_10(k) + Ts*dt_10(k), t_11(k) + Ts*dt_11(k), t_12(k) + Ts*dt_12(k)]; 
    
    t_2(k+1) = t_2(k) + Ts*dt_2(k);
    dt_2(k+1) = dt_2(k) + Ts*ddt_2(k);
end
end