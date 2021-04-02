function [F1_2_x,F1_2_y,F2_12_x,F2_12_y,F2_3_x,F2_3_y,F3_4_x,F3_4_y,F3_5_x,F3_5_y,F4_6_x,F4_6_y,F6_7_x,F6_7_y,F6_8_x,F6_8_y,F8_9_x, ...
          F8_9_y,F8_10_x,F8_10_y,F12_10_x,F12_10_y,F12_11_x,F12_11_y,F1_5_x,F1_5_y,F1_7_x,F1_7_y,F1_9_x,F1_9_y,F1_11_x,F1_11_y,M_2,M_2c,L_x] = ...
          ghassaei_dynamics(t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,ddt2,ddt3,ddt4,ddt5,ddt6,ddt7,ddt8,ddt9,ddt10,ddt11,ddt12,dt2,dt3,dt4,dt5,dt6,dt7,dt8,dt9,dt10,dt11,dt12, ...
                           r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,t)

F1_2_x = zeros(size(t2));
F1_2_y = zeros(size(t2));
F2_12_x = zeros(size(t2));
F2_12_y = zeros(size(t2));
F2_3_x = zeros(size(t2));
F2_3_y = zeros(size(t2));
F3_4_x = zeros(size(t2));
F3_4_y = zeros(size(t2));
F3_5_x = zeros(size(t2));
F3_5_y = zeros(size(t2));
F4_6_x = zeros(size(t2));
F4_6_y = zeros(size(t2));
F6_7_x = zeros(size(t2));
F6_7_y = zeros(size(t2));
F6_8_x = zeros(size(t2));
F6_8_y = zeros(size(t2));
F8_9_x = zeros(size(t2));
F8_9_y = zeros(size(t2));
F8_10_x = zeros(size(t2));
F8_10_y = zeros(size(t2));
F12_10_x = zeros(size(t2));
F12_10_y = zeros(size(t2));
F12_11_x = zeros(size(t2));
F12_11_y = zeros(size(t2));
F1_5_x = zeros(size(t2));
F1_5_y = zeros(size(t2));
F1_7_x = zeros(size(t2));
F1_7_y = zeros(size(t2));
F1_9_x = zeros(size(t2));
F1_9_y = zeros(size(t2));
F1_11_x = zeros(size(t2));
F1_11_y = zeros(size(t2));
M_2 = zeros(size(t2));

M_2c = zeros(size(t2));

L_x = zeros(size(t2));

v2_x = zeros(size(t2));
v2_y = zeros(size(t2));
v3_x = zeros(size(t2));
v3_y = zeros(size(t2));
v4_x = zeros(size(t2));
v4_y = zeros(size(t2));
v5_x = zeros(size(t2));
v5_y = zeros(size(t2));
v6_x = zeros(size(t2));
v6_y = zeros(size(t2));
v7_x = zeros(size(t2));
v7_y = zeros(size(t2));
v8_x = zeros(size(t2));
v8_y = zeros(size(t2));
v9_x = zeros(size(t2));
v9_y = zeros(size(t2));
v10_x = zeros(size(t2));
v10_y = zeros(size(t2));
v11_x = zeros(size(t2));
v11_y = zeros(size(t2));
v12_x = zeros(size(t2));
v12_y = zeros(size(t2));

a2_x = zeros(size(t2));
a2_y = zeros(size(t2));
a3_x = zeros(size(t2));
a3_y = zeros(size(t2));
a5_x = zeros(size(t2));
a5_y = zeros(size(t2));
a4_x = zeros(size(t2));
a4_y = zeros(size(t2));
a7_x = zeros(size(t2));
a7_y = zeros(size(t2));
a6_x = zeros(size(t2));
a6_y = zeros(size(t2));
a8_x = zeros(size(t2));
a8_y = zeros(size(t2));
a9_x = zeros(size(t2));
a9_y = zeros(size(t2));
a10_x = zeros(size(t2));
a10_y = zeros(size(t2));
a11_x = zeros(size(t2));
a11_y = zeros(size(t2));
a12_x = zeros(size(t2));
a12_y = zeros(size(t2));
            
t_size = size(t,1);

for k=1:t_size

v2_x(k) = -dt2(k)*r2*0.5*sin(t2(k));
v2_y(k) = dt2(k)*r2*0.5*cos(t2(k));
v5_x(k) = dt5(k)*r5*0.5*sin(t5(k));
v5_y(k) = -dt5(k)*r5*0.5*cos(t5(k));
v7_x(k) = dt7(k)*r7*0.5*sin(t7(k));
v7_y(k) = -dt7(k)*r7*0.5*cos(t7(k));
v9_x(k) = -dt9(k)*r9*0.5*sin(t9(k));
v9_y(k) = dt9(k)*r9*0.5*cos(t9(k));
v11_x(k) = -dt11(k)*r11*0.5*sin(t11(k));
v11_y(k) = dt11(k)*r11*0.5*cos(t11(k));
v3_x(k) = v5_x(k)*2 - dt3(k)*r3*0.5*sin(t3(k));
v3_y(k) = v5_y(k)*2 + dt3(k)*r3*0.5*cos(t3(k));
v4_x(k) = v5_x(k)*2 + dt4(k)*r4*0.5*sin(t4(k));
v4_y(k) = v5_y(k)*2 - dt4(k)*r4*0.5*cos(t4(k));
v6_x(k) = v7_x(k)*2 + dt6(k)*r6*0.5*sin(t6(k));
v6_y(k) = v7_y(k)*2 - dt6(k)*r6*0.5*cos(t6(k));
v8_x(k) = v7_x(k)*2 - dt8(k)*r8*0.5*sin(t8(k));
v8_y(k) = v7_y(k)*2 + dt8(k)*r8*0.5*cos(t8(k));
v10_x(k) = v9_x(k)*2 - dt10(k)*r10*0.5*sin(t10(k));
v10_y(k) = v9_y(k)*2 + dt10(k)*r10*0.5*cos(t10(k));
v12_x(k) = v2_x(k)*2 - dt12(k)*r12*0.5*sin(t12(k));
v12_y(k) = v2_y(k)*2 + dt12(k)*r12*0.5*cos(t12(k));
    
a2_x(k) = -dt2(k)*dt2(k)*r2*0.5*cos(t2(k)) + -ddt2(k)*r2*0.5*sin(t2(k));
a2_y(k) = -dt2(k)*dt2(k)*r2*0.5*sin(t2(k)) + ddt2(k)*r2*0.5*cos(t2(k));
a5_x(k) = dt5(k)*dt5(k)*r5*0.5*cos(t5(k)) + ddt5(k)*r5*0.5*sin(t5(k));
a5_y(k) = dt5(k)*dt5(k)*r5*0.5*sin(t5(k)) + -ddt5(k)*r5*0.5*cos(t5(k));
a7_x(k) = dt7(k)*dt7(k)*r7*0.5*cos(t7(k)) + ddt7(k)*r7*0.5*sin(t7(k));
a7_y(k) = dt7(k)*dt7(k)*r7*0.5*sin(t7(k)) + -ddt7(k)*r7*0.5*cos(t7(k));
a9_x(k) = -dt9(k)*dt9(k)*r9*0.5*cos(t9(k)) + -ddt9(k)*r9*0.5*sin(t9(k));
a9_y(k) = -dt9(k)*dt9(k)*r9*0.5*sin(t9(k)) + ddt9(k)*r9*0.5*cos(t9(k));
a11_x(k) = -dt11(k)*dt11(k)*r11*0.5*cos(t11(k)) + -ddt11(k)*r11*0.5*sin(t11(k));
a11_y(k) = -dt11(k)*dt11(k)*r11*0.5*sin(t11(k)) + ddt11(k)*r11*0.5*cos(t11(k));
a3_x(k) = a5_x(k)*2 + -dt3(k)*dt3(k)*r3*0.5*cos(t3(k)) + -ddt3(k)*r3*0.5*sin(t3(k));
a3_y(k) = a5_y(k)*2 + -dt3(k)*dt3(k)*r3*0.5*sin(t3(k)) + ddt3(k)*r3*0.5*cos(t3(k));
a4_x(k) = a5_x(k)*2 + dt4(k)*dt4(k)*r4*0.5*cos(t4(k)) + ddt4(k)*r4*0.5*sin(t4(k));
a4_y(k) = a5_y(k)*2 + dt4(k)*dt4(k)*r4*0.5*sin(t4(k)) + -ddt4(k)*r4*0.5*cos(t4(k));  
a6_x(k) = a7_x(k)*2 + dt6(k)*dt6(k)*r6*0.5*cos(t6(k)) + ddt6(k)*r6*0.5*sin(t6(k));
a6_y(k) = a7_y(k)*2 + dt6(k)*dt6(k)*r6*0.5*sin(t6(k)) + -ddt6(k)*r6*0.5*cos(t6(k));
a8_x(k) = a7_x(k)*2 + -dt8(k)*dt8(k)*r8*0.5*cos(t8(k)) + -ddt8(k)*r8*0.5*sin(t8(k));
a8_y(k) = a7_y(k)*2 + -dt8(k)*dt8(k)*r8*0.5*sin(t8(k)) + ddt8(k)*r8*0.5*cos(t8(k));
a10_x(k) = a9_x(k)*2 + -dt10(k)*dt10(k)*r10*0.5*cos(t10(k)) + -ddt10(k)*r10*0.5*sin(t10(k));
a10_y(k) = a9_y(k)*2 + -dt10(k)*dt10(k)*r10*0.5*sin(t10(k)) + ddt10(k)*r10*0.5*cos(t10(k));
a12_x(k) = a2_x(k)*2 + -dt12(k)*dt12(k)*r12*0.5*cos(t12(k)) + -ddt12(k)*r12*0.5*sin(t12(k));
a12_y(k) = a2_y(k)*2 + -dt12(k)*dt12(k)*r12*0.5*sin(t12(k)) + ddt12(k)*r12*0.5*cos(t12(k));
    
    
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
         0 0 r2*sin(t2(k)) -r2*cos(t2(k)) r2*sin(t2(k)) -r2*cos(t2(k)) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
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
    M_2(k) = x(33);
    
    %controle dynamica

     P2 = m2*(a2_x(k)*v2_x(k) + a2_y(k)*v2_y(k)) + I2*ddt2(k)*dt2(k);
     P3 = m3*(a3_x(k)*v3_x(k) + a3_y(k)*v3_y(k)) + I3*ddt3(k)*dt3(k);
     P4 = m4*(a4_x(k)*v4_x(k) + a4_y(k)*v4_y(k)) + I4*ddt4(k)*dt4(k);
     P5 = m5*(a5_x(k)*v5_x(k) + a5_y(k)*v5_y(k)) + I5*ddt5(k)*dt5(k);
     P6 = m6*(a6_x(k)*v6_x(k) + a6_y(k)*v6_y(k)) + I6*ddt6(k)*dt6(k);
     P7 = m7*(a7_x(k)*v7_x(k) + a7_y(k)*v7_y(k)) + I7*ddt7(k)*dt7(k);
     P8 = m8*(a8_x(k)*v8_x(k) + a8_y(k)*v8_y(k)) + I8*ddt8(k)*dt8(k);
     P9 = m9*(a9_x(k)*v9_x(k) + a9_y(k)*v9_y(k)) + I9*ddt9(k)*dt9(k);
     P10 = m10*(a10_x(k)*v10_x(k) + a10_y(k)*v10_y(k)) + I10*ddt10(k)*dt10(k);
     P11 = m11*(a11_x(k)*v11_x(k) + a11_y(k)*v11_y(k)) + I11*ddt11(k)*dt11(k);
     P12 = m12*(a12_x(k)*v12_x(k) + a12_y(k)*v12_y(k)) + I12*ddt12(k)*dt12(k);
    
    
     M_2c(k) =(P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + P10 + P11 + P12)/dt2(k);
    
     L_x(k) = -(m2*a2_x(k) + m3*a3_x(k) + m4*a4_x(k) + m5*a5_x(k) + m6*a7_x(k) + m7*a7_x(k) + m8*a8_x(k) + m9*a9_x(k) + m10*a10_x(k) + m11*a11_x(k) + m12*a12_x(k)); 
     
end

end