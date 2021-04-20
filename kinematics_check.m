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
    subplot(4,6,1)
    plot(t,real(A1)-real(A2))
    xlabel('time [s]')
    ylabel('real(A1)-real(A2)')
    axis tight
    
    subplot(4,6,2)
    plot(t,imag(A1)-imag(A2))
    xlabel('time [s]')
    ylabel('imag(A1)-imag(A2)')
    axis tight
    
    subplot(4,6,3)
    plot(t,real(A1)-real(A3))
    xlabel('time [s]')
    ylabel('real(A1)-real(A3)')
    axis tight
    
    subplot(4,6,4)
    plot(t,imag(A1)-imag(A3))
    xlabel('time [s]')
    ylabel('imag(A1)-imag(A3)')
    axis tight
    
    subplot(4,6,5)
    plot(t,real(A1)-real(A4))
    xlabel('time [s]')
    ylabel('real(A1)-real(A4)')
    axis tight
    
    subplot(4,6,6)
    plot(t,imag(A1)-imag(A4))
    xlabel('time [s]')
    ylabel('imag(A1)-imag(A4)')
    axis tight
    
    %point B
    subplot(4,6,7)
    plot(t,real(B1)-real(B2))
    xlabel('time [s]')
    ylabel('real(B1)-real(B2)')
    axis tight
    
    subplot(4,6,8)
    plot(t,imag(B1)-imag(B2))
    xlabel('time [s]')
    ylabel('imag(B1)-imag(B2)')
    axis tight
    
    %point C
    subplot(2,6,9)
    plot(t,real(C1)-real(C2))
    xlabel('time [s]')
    ylabel('real(C1)-real(C2)')
    axis tight
    
    subplot(2,6,10)
    plot(t,imag(C1)-imag(C2))
    xlabel('time [s]')
    ylabel('imag(C1)-imag(C2)')
    axis tight
    
    %point H
    subplot(2,6,11)
    plot(t,real(H1)-real(H2))
    xlabel('time [s]')
    ylabel('real(H1)-real(H2)')
    axis tight
    
    subplot(2,6,12)
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

acc_H1 = cross(omega11,cross(omega11,AH)) + cross(alpha11,AH);
acc_H2 = cross(omega2,cross(omega2,BC)) + cross(alpha2,BC) - cross(omega12,cross(omega12,HC)) - cross(alpha12,HC);

acc_C1 = acc_H1 + cross(omega12,cross(omega12,HC)) + cross(alpha12,HC);
acc_C2 = cross(omega2,cross(omega2,BC)) + cross(alpha2,BC);
acc_C3 = cross(omega4,cross(omega4,AD)) + cross(alpha4,AD) - cross(omega3,cross(omega3,CD)) - cross(alpha3,CD);

acc_A1 = acc_C2 - cross(omega12,cross(omega12,HC)) - cross(alpha12,HC) - acc_H1;
acc_A2 = acc_C2 + cross(omega3,cross(omega3,CD)) + cross(alpha3,CD) - cross(omega4,cross(omega4,AD)) - cross(alpha4,AD);
acc_A3 = cross(omega4,cross(omega4,AD)) + cross(alpha4,AD)+ cross(omega5,cross(omega5,DE)) + cross(alpha5,DE) + cross(omega6,cross(omega6,EF)) + cross(alpha6,EF) - cross(omega11,cross(omega11,AF)) - cross(alpha11,AF);

acc_B1 = acc_C1 - cross(omega2,cross(omega2,BC)) - cross(alpha2,BC);
acc_B2 = acc_C2 - cross(omega12,cross(omega12,HC)) - cross(alpha12,HC) - acc_H1; 

%plot errors 
if plot_kinCheck
    figure
    subplot(5,2,1)
    plot(t,acc_A1 - acc_A2)
    xlabel('t [s]')
    ylabel('acc_{A1} - acc_{A2} [m/s^2] ')
    axis tight

    subplot(5,2,2)
    plot(t,acc_A1 - acc_A3)
    xlabel('t [s]')
    ylabel('acc_{A1} - acc_{A3} [m/s^2] ')
    axis tight

    subplot(5,2,3)
    plot(t,acc_B1 - acc_B2)
    xlabel('t [s]')
    ylabel('acc_{B1} - acc_{B2} [m/s^2] ')
    axis tight

    subplot(5,2,4)
    plot(t,acc_C1 - acc_C2)
    xlabel('t [s]')
    ylabel('acc_{C1} - acc_{C2} [m/s^2] ')
    axis tight

    subplot(5,2,5)
    plot(t,acc_C1 - acc_C3)
    xlabel('t [s]')
    ylabel('acc_{C1} - acc_{C3} [m/s^2] ')
    axis tight
end
        
end

