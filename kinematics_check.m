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
D2 = zeros(size(t));

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
    
    
end    
        
end

