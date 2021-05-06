%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Linkage system: Linkage System of Ghassaei 
%
% Seppe Vilain, Arnoud Deboeck
%
% Kinematic analysis check of position, velocity and acceleration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function  kinematics_check(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,...
            t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,dt2,dt3,dt4,dt5,dt6,dt7,dt8,dt9,dt10,...
            dt11,dt12,ddt2,ddt3,ddt4,ddt5,ddt6,ddt7,ddt8,ddt9,ddt10,ddt11,ddt12, t, plot_kinCheck)
        
%%%%%%%%%%%%
% Position %
%%%%%%%%%%%%

% CALCULATION OF POSITION FROM DIFFERENT POINTS IN COMPLEX NOTATION

% initialize iteration variables
% 4 loops   AHCBA
%           ADEFA
%           AFGA
%           AGHA

%loop1
A1 = zeros(size(t));
A2 = zeros(size(t));
A3 = zeros(size(t));
A4 = zeros(size(t));

B1 = zeros(size(t));
B2 = zeros(size(t));
B3= zeros(size(t));

C1 = zeros(size(t));
C2 = zeros(size(t));
C3 = zeros(size(t));

D1 = zeros(size(t));

H1 = zeros(size(t));
H2 = zeros(size(t));


t_size = size(t,1);
for k=1:t_size
    
    %loop1
    B1(k) = A1(k) + r1 * exp(j*t1);
    C1(k) = B1(k) + r2 * exp(j*t2(k));
    H1(k) = C1(k) + r12 * exp(j*t12(k));
    A2(k) = H1(k) - r11 * exp(j*t11(k));
   
    H2(k) = A1(k) + r11 * exp(j*t11(k));
    C2(k) = H2(k) - r12 * exp(j*t12(k));
    B2(k) = C2(k) - r2 * exp(j*t2(k));
    A3(k) = B2(k) - r1 * exp(j*t1);
    
    %loop 2
    B3(k) = A1(k) + r1 * exp(j*t1);
    C3(k) = B3(k) + r2 * exp(j*t2(k));
    D1(k) = C3(k) - r3 * exp(j*t3(k));
    A4(k) = D1(k) + r5 * exp(j*t5(k));

end

% ##############
% #### plot ####
% ##############

if plot_kinCheck
    
    figure
    % POINT A
    subplot(3,4,1)
    plot(t,real(A1)-real(A2))
    xlabel('time [s]')
    ylabel('real(A1)-real(A2)')
    axis tight
    
    subplot(3,4,2)
    plot(t,imag(A1)-imag(A2))
    xlabel('time [s]')
    ylabel('imag(A1)-imag(A2)')
    axis tight
    
    subplot(3,4,3)
    plot(t,real(A1)-real(A3))
    xlabel('time [s]')
    ylabel('real(A1)-real(A3)')
    axis tight
    
    subplot(3,4,4)
    plot(t,imag(A1)-imag(A3))
    xlabel('time [s]')
    ylabel('imag(A1)-imag(A3)')
    axis tight
    
    subplot(3,4,5)
    plot(t,real(A1)-real(A4))
    xlabel('time [s]')
    ylabel('real(A1)-real(A4)')
    axis tight
    
    subplot(3,4,6)
    plot(t,imag(A1)-imag(A4))
    xlabel('time [s]')
    ylabel('imag(A1)-imag(A4)')
    axis tight
    
    %point B
    subplot(3,4,7)
    plot(t,real(B1)-real(B2))
    xlabel('time [s]')
    ylabel('real(B1)-real(B2)')
    axis tight
    
    subplot(3,4,8)
    plot(t,imag(B1)-imag(B2))
    xlabel('time [s]')
    ylabel('imag(B1)-imag(B2)')
    axis tight
    
    %point C
    subplot(3,4,9)
    plot(t,real(C1)-real(C2))
    xlabel('time [s]')
    ylabel('real(C1)-real(C2)')
    axis tight
    
    subplot(3,4,10)
    plot(t,imag(C1)-imag(C2))
    xlabel('time [s]')
    ylabel('imag(C1)-imag(C2)')
    axis tight
    
    %point H
    subplot(3,4,11)
    plot(t,real(H1)-real(H2))
    xlabel('time [s]')
    ylabel('real(H1)-real(H2)')
    axis tight
    
    subplot(3,4,12)
    plot(t,imag(H1)-imag(H2))
    xlabel('time [s]')
    ylabel('imag(H1)-imag(H2)')
    axis tight
   
    
    
end   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  VELOCITY and ACCELERATION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%corner velocity 
omega2 = [zeros(size(t2)) zeros(size(t2)) dt2];
omega3 = [zeros(size(t2)) zeros(size(t2)) dt3];
omega4 = [zeros(size(t2)) zeros(size(t2)) dt4];
omega5 = [zeros(size(t2)) zeros(size(t2)) dt5];
omega6 = [zeros(size(t2)) zeros(size(t2)) dt6];
omega7 = [zeros(size(t2)) zeros(size(t2)) dt7];
omega11 = [zeros(size(t2)) zeros(size(t2)) dt11];
omega12 = [zeros(size(t2)) zeros(size(t2)) dt12];

%corner acceleration
alpha2 = [zeros(size(t2)) zeros(size(t2)) ddt2];
alpha3 = [zeros(size(t2)) zeros(size(t2)) ddt3];
alpha4 = [zeros(size(t2)) zeros(size(t2)) ddt4];
alpha5 = [zeros(size(t2)) zeros(size(t2)) ddt5];
alpha6 = [zeros(size(t2)) zeros(size(t2)) ddt6];
alpha7 = [zeros(size(t2)) zeros(size(t2)) ddt7];
alpha11 = [zeros(size(t2)) zeros(size(t2)) ddt11];
alpha12 = [zeros(size(t2)) zeros(size(t2)) ddt12];

%length
BC = [r2*cos(t2) r2*sin(t2) zeros(size(t2))];
CD = [r3*cos(t3) r3*sin(t3) zeros(size(t2))];
DE = [r4*cos(t4) r4*sin(t4) zeros(size(t2))];
AD = [r5*cos(t5) r5*sin(t5) zeros(size(t2))];
EF = [r6*cos(t6) r6*sin(t6) zeros(size(t2))];
AF = [r7*cos(t7) r7*sin(t7) zeros(size(t2))];
AH = [r11*cos(t11) r11*sin(t11) zeros(size(t2))];
HC = [r12*cos(t12) r12*sin(t12) zeros(size(t2))];

%Calculate velocity of a point starting from 2 different points

vel_A1 = 0 + cross(omega2,BC) - cross(omega3,CD) + cross(omega5,AD) ;
vel_A2 = 0 + cross(omega2,BC) + cross(omega12,HC) - cross(omega11,AH);
vel_A3 = - cross(omega5,AD) - cross(omega4,DE) + cross(omega6,EF) + cross(omega7, AF);

vel_B1 = cross(omega2,BC) - cross(omega3,CD) + cross(omega5,AD) +0 ;
vel_B2 = cross(omega2,BC) + cross(omega12,HC) - cross(omega11,AH) + 0;

vel_C1 = cross(omega11,AH) - cross(omega12,HC);
vel_C2 = 0 + cross(omega2,BC);
vel_C3 = -cross(omega5,AD) + cross(omega3,CD);



% plot errors 
if plot_kinCheck
    figure
    subplot(3,2,1)
    plot(t,vel_A1 - vel_A2)
    xlabel('t [s]')
    ylabel('vel_{A1} - vel_{A2} [m/s] ')
    axis tight

    subplot(3,2,2)
    plot(t,vel_A1 - vel_A3)
    xlabel('t [s]')
    ylabel('vel_{A1} - vel_{A3} [m/s] ')
    axis tight

    subplot(3,2,3)
    plot(t,vel_B1 - vel_B2)
    xlabel('t [s]')
    ylabel('vel_{B1} - vel_{B2} [m/s] ')
    axis tight

    subplot(3,2,4)
    plot(t,vel_C1 - vel_C2)
    xlabel('t [s]')
    ylabel('vel_{C1} - vel_{C2} [m/s] ')
    axis tight

    subplot(3,2,5)
    plot(t,vel_C1 - vel_C3)
    xlabel('t [s]')
    ylabel('vel_{C1} - vel_{C3} [m/s] ')
    axis tight
end

%acceleration 

acc_H1 = cross(omega11,cross(omega11,AH)) + cross(alpha11,AH); %
acc_H2 = cross(omega2,cross(omega2,BC)) + cross(alpha2,BC) + cross(omega12,cross(omega12,HC)) + cross(alpha12,HC);%

acc_C1 = acc_H1 - cross(omega12,cross(omega12,HC)) - cross(alpha12,HC);%
acc_C2 = cross(omega2,cross(omega2,BC)) + cross(alpha2,BC);%
acc_C3 = -cross(omega5,cross(omega5,AD)) - cross(alpha5,AD) + cross(omega3,cross(omega3,CD)) + cross(alpha3,CD);

acc_A1 = acc_C2 + cross(omega12,cross(omega12,HC)) + cross(alpha12,HC) - acc_H1;
acc_A2 = acc_C2 - cross(omega3,cross(omega3,CD)) - cross(alpha3,CD) + cross(omega5,cross(omega5,AD)) + cross(alpha5,AD);
acc_A3 = -cross(omega5,cross(omega5,AD)) - cross(alpha5,AD) - cross(omega4,cross(omega4,DE)) - cross(alpha4,DE) + cross(omega6,cross(omega6,EF)) + cross(alpha6,EF) + cross(omega7,cross(omega7,AF)) + cross(alpha7,AF);

acc_B1 = acc_C1 + cross(omega2,cross(omega2,BC)) + cross(alpha2,BC);%
acc_B2 = acc_C2 - cross(omega12,cross(omega12,HC)) - cross(alpha12,HC) + acc_H1; %

%plot errors 
if plot_kinCheck
    figure
    subplot(3,2,1)
    plot(t,acc_A1 - acc_A2)
    xlabel('t [s]')
    ylabel('acc_{A1} - acc_{A2} [m/s^2] ')
    axis tight

    subplot(3,2,2)
    plot(t,acc_A1 - acc_A3)
    xlabel('t [s]')
    ylabel('acc_{A1} - acc_{A3} [m/s^2] ')
    axis tight

    subplot(3,2,3)
    plot(t,acc_B1 - acc_B2)
    xlabel('t [s]')
    ylabel('acc_{B1} - acc_{B2} [m/s^2] ')
    axis tight

    subplot(3,2,4)
    plot(t,acc_C1 - acc_C2)
    xlabel('t [s]')
    ylabel('acc_{C1} - acc_{C2} [m/s^2] ')
    axis tight

    subplot(3,2,5)
    plot(t,acc_C1 - acc_C3)
    xlabel('t [s]')
    ylabel('acc_{C1} - acc_{C3} [m/s^2] ')
    axis tight
    
    subplot(3,2,6)
    plot(t,acc_H1 - acc_H2)
    xlabel('t [s]')
    ylabel('acc_{H1} - acc_{H2} [m/s^2] ')
    axis tight

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYTICAL SOLUTION VS NUMERICAL DERIVATIVE OF VELOCITY AND ACCELERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize iteration variables
difft3 = zeros(size(t));
difft4 = zeros(size(t));
difft5 = zeros(size(t));
difft6 = zeros(size(t));
difft11 = zeros(size(t));
difft12 = zeros(size(t));
ddifft3 = zeros(size(t));
ddifft4 = zeros(size(t));
ddifft5 = zeros(size(t));
ddifft6 = zeros(size(t));
ddifft11 = zeros(size(t));
ddifft12 = zeros(size(t));
errdt3 = zeros(size(t));
errdt4 = zeros(size(t));
errdt5 = zeros(size(t));
errdt6 = zeros(size(t));
errdt11 = zeros(size(t));
errdt12 = zeros(size(t));
errddt3 = zeros(size(t));
errddt4 = zeros(size(t));
errddt5 = zeros(size(t));
errddt6 = zeros(size(t));
errddt11 = zeros(size(t));
errddt12 = zeros(size(t));

% finite difference calculation of first and second derivative difft_i and ddifft_i
% f'(x) = (f(x+1)-f(x-1))/(2*Ts)
% f''(x) = (f(x+1)-2*f(x)+f(x-1))/(Ts^2)
Ts = t(2) - t(1);      % timestep
for k=2:t_size-1       % (2:n-1) instead of (1:n) because f(-1) and f(n+1) do not exist
    difft3(k) = (t3(k+1) - t3(k-1)) / (2*Ts); 
    difft4(k) = (t4(k+1) - t4(k-1)) / (2*Ts); 
    difft5(k) = (t5(k+1) - t5(k-1)) / (2*Ts);
    difft6(k) = (t6(k+1) - t6(k-1)) / (2*Ts);
    difft11(k) = (t11(k+1) - t11(k-1)) / (2*Ts);
    difft12(k) = (t12(k+1) - t12(k-1)) / (2*Ts);
    ddifft3(k) = (t3(k+1) - 2*t3(k) + t3(k-1)) / (Ts^2); 
    ddifft4(k) = (t4(k+1) - 2*t4(k) + t4(k-1)) / (Ts^2); 
    ddifft5(k) = (t5(k+1) - 2*t5(k) + t5(k-1)) / (Ts^2);
    ddifft6(k) = (t6(k+1) - 2*t6(k) + t6(k-1)) / (Ts^2);
    ddifft11(k) = (t11(k+1) - 2*t11(k) + t11(k-1)) / (Ts^2);
    ddifft12(k) = (t12(k+1) - 2*t12(k) + t12(k-1)) / (Ts^2);
end

% extrapolation for difft_i(1), difft_i(n), ddifft_i(1) and ddifft_i(n)
% y = y_1 + ((x-x_1)/(x_2-x_1))*(y_2-y_1) met (x-x_1)/(x_2-x_1)=-1
% y = 2*y_1 - y_2
difft3(1) = 2*difft3(2) - difft3(3);
difft4(1) = 2*difft4(2) - difft4(3); 
difft5(1) = 2*difft5(2) - difft5(3);
difft6(1) = 2*difft6(2) - difft6(3);
difft11(1) = 2*difft11(2) - difft11(3);
difft12(1) = 2*difft12(2) - difft12(3);
difft3(t_size) = 2*difft3(t_size-1) - difft3(t_size-2);
difft4(t_size) = 2*difft4(t_size-1) - difft4(t_size-2);
difft5(t_size) = 2*difft5(t_size-1) - difft5(t_size-2);
difft6(t_size) = 2*difft6(t_size-1) - difft6(t_size-2);
difft11(t_size) = 2*difft11(t_size-1) - difft11(t_size-2); 
difft12(t_size) = 2*difft12(t_size-1) - difft12(t_size-2);
ddifft3(1) = 2*ddifft3(2) - ddifft3(3);
ddifft4(1) = 2*ddifft4(2) - ddifft4(3); 
ddifft5(1) = 2*ddifft5(2) - ddifft5(3);
ddifft6(1) = 2*ddifft6(2) - ddifft6(3);
ddifft11(1) = 2*ddifft11(2) - ddifft11(3);
ddifft12(1) = 2*ddifft12(2) - ddifft12(3);
ddifft3(t_size) = 2*ddifft3(t_size-1) - ddifft3(t_size-2);
ddifft4(t_size) = 2*ddifft4(t_size-1) - ddifft4(t_size-2);
ddifft5(t_size) = 2*ddifft5(t_size-1) - ddifft5(t_size-2);
ddifft6(t_size) = 2*ddifft6(t_size-1) - ddifft6(t_size-2);
ddifft11(t_size) = 2*ddifft11(t_size-1) - ddifft11(t_size-2); 
ddifft12(t_size) = 2*ddifft12(t_size-1) - ddifft12(t_size-2);

% calculation of errdt_i and errddt_i
for k=1:t_size
    errdt3(k) = difft3(k) - dt3(k);
    errdt4(k) = difft4(k) - dt4(k); 
    errdt5(k) = difft5(k) - dt5(k);
    errdt6(k) = difft6(k) - dt6(k);
    errdt11(k) = difft11(k) - dt11(k);
    errdt12(k) = difft12(k) - dt12(k);
    errddt3(k) = ddifft3(k) - ddt3(k);
    errddt4(k) = ddifft4(k) - ddt4(k); 
    errddt5(k) = ddifft5(k) - ddt5(k);
    errddt6(k) = ddifft6(k) - ddt6(k);
    errddt11(k) = ddifft11(k) - ddt11(k);
    errddt12(k) = ddifft12(k) - ddt12(k);
end

% **********************
% *** plot figures ***
% **********************
if plot_kinCheck
    figure
    subplot(231)
    plot(t,errdt3)
    xlabel('t [s]')
    ylabel('errdt3 [rad/s]') 
    subplot(232)
    plot(t,errdt4)
    xlabel('t [s]')
    ylabel('errdt4 [rad/s]')
    subplot(233)
    plot(t,errdt5)
    xlabel('t [s]')
    ylabel('errdt5 [rad/s]') 
    subplot(234)
    plot(t,errdt6)
    xlabel('t [s]')
    ylabel('errdt6 [rad/s]')
    subplot(235)
    plot(t,errdt11)
    xlabel('t [s]')
    ylabel('errdt11 [rad/s]')
    subplot(236)
    plot(t,errdt12)
    xlabel('t [s]')
    ylabel('errdt12 [rad/s]')
    
    figure
    subplot(231)
    plot(t,errddt3)
    xlabel('t [s]')
    ylabel('errddt3 [rad/s^2]') 
    subplot(232)
    plot(t,errddt4)
    xlabel('t [s]')
    ylabel('errddt4 [rad/s^2]')
    subplot(233)
    plot(t,errddt5)
    xlabel('t [s]')
    ylabel('errddt5 [rad/s^2]') 
    subplot(234)
    plot(t,errddt6)
    xlabel('t [s]')
    ylabel('errddt6 [rad/s^2]')
    subplot(235)
    plot(t,errddt11)
    xlabel('t [s]')
    ylabel('errddt11 [rad/s^2]')
    subplot(236)
    plot(t,errddt12)
    xlabel('t [s]')
    ylabel('errddt12 [rad/s^2]')
end
        
end

