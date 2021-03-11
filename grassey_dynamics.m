function [F1_2_x,F1_2_y,F2_12_x,F2_12_y,F2_3_x,F2_3_y,F3_4_x,F3_4_y,F3_5_x,F3_5_y,F4_6_x,F4_6_y,F6_7_x,F6_7_y,F6_8_x,F6_8_y,F8_9_x, ...
          F8_9_y,F8_10_x,F8_10_y,F12_10_x,F12_10_y,F12_11_x,F12_11_y,F1_5_x,F1_5_y,F1_7_x,F1_7_y,F1_9_x,F1_9_y,F1_11_x,F1_11_y,M_2] = ...
          grassey_dynamics(t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,ddt2,ddt3,ddt4,ddt5,ddt6,ddt7,ddt8,ddt9,ddt10,ddt11,ddt12,dt2,dt3,dt4,dt5,dt6,dt7,dt8,dt9,dt10,dt11,dt12, ...
                           r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,t)

F1_2_x = zeros(size(phi2));
F1_2_y = zeros(size(phi2));
F2_12_x = zeros(size(phi2));
F2_12_y = zeros(size(phi2));
F2_3_x = zeros(size(phi2));
F2_3_y = zeros(size(phi2));
F3_4_x = zeros(size(phi2));
F3_4_y = zeros(size(phi2));
F3_5_x = zeros(size(phi2));
F3_5_y = zeros(size(phi2));
F4_6_x = zeros(size(phi2));
F4_6_y = zeros(size(phi2));
F6_7_x = zeros(size(phi2));
F6_7_y = zeros(size(phi2));
F6_8_x = zeros(size(phi2));
F6_8_y = zeros(size(phi2));
F8_9_x = zeros(size(phi2));
F8_9_y = zeros(size(phi2));
F8_10_x = zeros(size(phi2));
F8_10_y = zeros(size(phi2));
F12_10_x = zeros(size(phi2));
F12_10_y = zeros(size(phi2));
F12_11_x = zeros(size(phi2));
F12_11_y = zeros(size(phi2));
F1_5_x = zeros(size(phi2));
F1_5_y = zeros(size(phi2));
F1_7_x = zeros(size(phi2));
F1_7_y = zeros(size(phi2));
F1_9_x = zeros(size(phi2));
F1_9_y = zeros(size(phi2));
F1_11_x = zeros(size(phi2));
F1_11_y = zeros(size(phi2));
M_2 = zeros(size(phi2));

a2_x = -dt2*dt2*r2*0.5*cos(t2) + -ddt2*r2*0.5*sin(t2);
a2_y = -dt2*dt2*r2*0.5*sin(t2) + ddt2*r2*0.5*cos(t2);

a3_x = a2_x + -dt3*dt3*r3*0.5*cos(t3) + -ddt3*r3*0.5*sin(t3);
a3_y = a2_y + -dt3*dt3*r3*0.5*sin(t3) + ddt3*r3*0.5*cos(t3);

a5_x = -dt5*dt5*r5*0.5*cos(t5) + -ddt5*r5*0.5*sin(t5);
a5_y = -dt5*dt5*r5*0.5*sin(t5) + ddt5*r5*0.5*cos(t5);

a4_x = a5_x + -dt4*dt4*r4*0.5*cos(t4) + -ddt4*r4*0.5*sin(t4);
a4_y = a5_y + -dt4*dt4*r4*0.5*sin(t4) + ddt4*r4*0.5*cos(t4);

a7_x = -dt7*dt7*r7*0.5*cos(t7) + -ddt7*r7*0.5*sin(t7);
a7_y = -dt7*dt7*r7*0.5*sin(t7) + ddt7*r7*0.5*cos(t7);

a6_x = a7_x + -dt6*dt6*r6*0.5*cos(t6) + -ddt6*r6*0.5*sin(t6);
a6_y = a7_y + -dt6*dt6*r2*0.5*sin(t6) + ddt6*r6*0.5*cos(t6);

a8_x = a7_x + -dt8*dt8*r3*0.5*cos(t8) + -ddt8*r8*0.5*sin(t8);
a8_y = a7_y + -dt8*dt8*r8*0.5*sin(t8) + ddt8*r8*0.5*cos(t8);

a9_x = -dt9*dt9*r9*0.5*cos(t9) + -ddt9*r9*0.5*sin(t9);
a9_y = -dt9*dt9*r9*0.5*sin(t9) + ddt9*r9*0.5*cos(t9);

a10_x = a9_x + -dt10*dt10*r10*0.5*cos(t10) + -ddt10*r10*0.5*sin(t10);
a10_y = a9_y + -dt10*dt10*r10*0.5*sin(t10) + ddt10*r10*0.5*cos(t10);

a11_x = -dt11*dt11*r11*0.5*cos(t11) + -ddt11*r11*0.5*sin(t11);
a11_y = -dt11*dt11*r11*0.5*sin(t11) + ddt11*r11*0.5*cos(t11);

a12_x = a2_x + -dt12*dt12*r12*0.5*cos(t12) + -ddt12*r12*0.5*sin(t12);
a12_y = a2_y + -dt12*dt12*r12*0.5*sin(t12) + ddt12*r12*0.5*cos(t12);
            
t_size = size(t,1);

for k=1:t_size
    
    A = [1 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 1 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 1 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 1 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 1 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 1 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 1 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 1 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0;
         0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0;
         0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 -1 0 0 0 0 0 0 0 0 0;
         0 0 r2*sin(t2(k)) r2*-cos(t2(k)) r2*sin(t2(k)) r2*-cos(t2(k)) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
         0 0 0 0 -sin(t3(k)) cos(t3(k)) sin(t3(k)+pi) -cos(t3(k)+pi) sin(t3(k)+pi) -cos(t3(k)+pi) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 -sin(t4(k)) cos(t4(k)) 0 0 sin(t4(k)+pi) -cos(t4(k)+pi) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 -sin(t5(k)+pi) cos(t5(k)+pi) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -sin(t5(k)) cos(t5(k)) 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 -sin(t6(k)+pi) cos(t6(k)+pi) sin(t6(k)) -cos(t6(k)) sin(t6(k)) -cos(t6(k)) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 -sin(t7(k)+pi) cos(t7(k)+pi) 0 0 0 0 0 0 0 0 0 0 0 0 -sin(t7(k)) cos(t7(k)) 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 -sin(t8(k)+pi) cos(t8(k)+pi) sin(t8(k)) -cos(t8(k)) sin(t8(k)) -cos(t8(k)) 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -sin(t9(k)) cos(t9(k)) 0 0 0 0 0 0 0 0 0 0 -sin(t9(k)+pi) cos(t9(k)+pi) 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -sin(t10(k)+pi) cos(t10(k)+pi) -sin(t10(k)) cos(t10(k)) 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -sin(t11(k)) cos(t11(k)) 0 0 0 0 0 0 -sin(t11(k)+pi) cos(t11(k)+pi) 0;
         0 0 -sin(t12(k)+pi) cos(t12(k)+pi) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 sin(t12(k)) -cos(t12(k)) sin(t12(k)) -cos(t12(k)) 0 0 0 0 0 0 0 0 0];
    
    B = [m2*a2_x(k);
         m2*a2_y(k);
         m3*a3_x(k);
         m3*a3_y(k);
         m4*a4_x(k);
         m4*a4_y(k);
         m5*a5_x(k);
         m5*a5_y(k);
         m6*a6_x(k);
         m6*a6_y(k);
         m7*a7_x(k);
         m7*a7_y(k);
         m8*a8_x(k);
         m8*a8_y(k);
         m9*a9_x(k);
         m9*a9_y(k);
         m10*a10_x(k);
         m10*a10_y(k);
         m11*a11_x(k);
         m11*a11_y(k);
         m12*a12_x(k);
         m12*a12_y(k);
         I2*ddt2(k);
         I3*ddt3(k)*2/r3;
         I4*ddt4(k)*2/r4;
         I5*ddt5(k)*2/r5;
         I6*ddt6(k)*2/r6;
         I7*ddt7(k)*2/r7;
         I8*ddt8(k)*2/r8;
         I9*ddt9(k)*2/r9;
         I10*ddt10(k)*2/r10;
         I11*ddt11(k)*2/r11;
         I12*ddt12(k)*2/r12];
    
    x = A/B;
    
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
    M_2(k) = x(33);
    
end


end