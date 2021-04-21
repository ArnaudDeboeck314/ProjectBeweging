% #####################################################
%
% Linkage system: Linkage system of Ghassaei
% 
% Seppe Vilain, Arnoud Deboeck
% 
% #####################################################
function [t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,dt3,dt4,dt5,dt6,dt7,dt8,dt9,dt10,...
            dt11,dt12,ddt3,ddt4,ddt5,ddt6,ddt7,ddt8,ddt9,ddt10,ddt11,ddt12] ...
            = kinematics_ghassaei(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,t1,t2,dt2,ddt2,t_init,t,plot_kin)

t3 = zeros(size(t));
t4 = zeros(size(t));
t5 = zeros(size(t));
t6 = zeros(size(t));
t7 = zeros(size(t));
t8 = zeros(size(t));
t9 = zeros(size(t));
t10 = zeros(size(t));
t11 = zeros(size(t));
t12 = zeros(size(t));

dt3 = zeros(size(t));
dt4 = zeros(size(t));
dt5 = zeros(size(t));
dt6 = zeros(size(t));
dt7 = zeros(size(t));
dt8 = zeros(size(t));
dt9 = zeros(size(t));
dt10 = zeros(size(t));
dt11 = zeros(size(t));
dt12 = zeros(size(t));

ddt3 = zeros(size(t));
ddt4 = zeros(size(t));
ddt5 = zeros(size(t));
ddt6 = zeros(size(t));
ddt7 = zeros(size(t));
ddt8 = zeros(size(t));
ddt9 = zeros(size(t));
ddt10 = zeros(size(t));
ddt11 = zeros(size(t));
ddt12 = zeros(size(t));

% fsolve options (help fsolve, help optimset)
optim_options = optimset('Display','off');

% *** loop over positions ***
Ts = t(2) - t(1);    % timestep 
t_size = size(t,1);  % number of simulation steps

for k=1:t_size
    
     % *** position analysis ***
    
    % fsolve solves the non-linear set of equations
    % loop closure equations: see loop_closure_eqs.m
    % argument loop_closure_eqs: file containing closure equations
    % argument [..]': initial values of unknown angles phi3 and phi4
    % argument optim options: parameters for fsolve
    % argument phi2(k): input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
    % argument a1 ... phi1: constants
    % return value x: solution for the unknown angles phi3 and phi4
    % return exitflag: indicates convergence of algorithm
    [x, ~, exitflag] = fsolve('Loop_closure_ghassaei',t_init',optim_options,t2(k),r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,t1);
    
    if (exitflag ~= 1)
        disp 'The fsole exit flag was not 1, probably no convergence!'
    end
    
    % save the results of the fsolve above
    t3(k) = x(1);
    t4(k) = x(2);
    t5(k) = x(3);
    t6(k) = x(4);
    t7(k) = x(5);
    t8(k) = x(6);
    t9(k) = x(7);
    t10(k) = x(8);
    t11(k) = x(9);
    t12(k) = x(10);
    
    % ###Velocity Analysis###
    
    A = [0, 0, 0, 0, 0, 0, 0, 0, r11*sin(t11(k)), -r12*sin(t12(k));
         0, 0, 0, 0, 0, 0, 0, 0, -r11*cos(t11(k)), r12*cos(t12(k));
         0, 0, 0, 0, -r7*sin(t7(k)), r8*sin(t8(k)), -r9*sin(t9(k)), 0, 0, 0;
         0, 0, 0, 0, r7*cos(t7(k)), -r8*cos(t8(k)), r9*cos(t9(k)), 0, 0, 0;
         0, -r4*sin(t4(k)), -r5*sin(t5(k)), r6*sin(t6(k)), r7*sin(t7(k)), 0, 0, 0, 0, 0;
         0, r4*cos(t4(k)), r5*cos(t5(k)), -r6*cos(t6(k)), -r7*cos(t7(k)), 0, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 0, -r9*sin(t9(k)), -r10*sin(t10(k)), r11*sin(t11(k)), 0;
         0, 0, 0, 0, 0, 0, r9*cos(t9(k)), r10*cos(t10(k)), -r11*cos(t11(k)), 0;
         r3*sin(t3(k)), 0, -r5*sin(t5(k)), 0, 0, 0, 0, 0, 0, 0;
         -r3*cos(t3(k)), 0, r5*cos(t5(k)), 0, 0, 0, 0, 0, 0, 0];
        
    B = [r2*sin(t2(k))*dt2(k);
         -r2*cos(t2(k))*dt2(k);
         0;
         0;
         0;
         0;
         0;
         0;
         r2*sin(t2(k))*dt2(k);
         -r2*cos(t2(k))*dt2(k)];
    
    x = A\B;
    
    %results
    dt3(k) = x(1);
    dt4(k) = x(2);
    dt5(k) = x(3);
    dt6(k) = x(4);
    dt7(k) = x(5);
    dt8(k) = x(6);
    dt9(k) = x(7);
    dt10(k) = x(8);
    dt11(k) = x(9);
    dt12(k) = x(10);
    
    % ###Acceleration Analysis###
    
    A = [0, 0, 0, 0, 0, 0, 0, 0, r11*sin(t11(k)), -r12*sin(t12(k));
         0, 0, 0, 0, 0, 0, 0, 0, -r11*cos(t11(k)), r12*cos(t12(k));
         0, 0, 0, 0, -r7*sin(t7(k)), r8*sin(t8(k)), -r9*sin(t9(k)), 0, 0, 0;
         0, 0, 0, 0, r7*cos(t7(k)), -r8*cos(t8(k)), r9*cos(t9(k)), 0, 0, 0;
         0, -r4*sin(t4(k)), -r5*sin(t5(k)), r6*sin(t6(k)), r7*sin(t7(k)), 0, 0, 0, 0, 0;
         0, r4*cos(t4(k)), r5*cos(t5(k)), -r6*cos(t6(k)), -r7*cos(t7(k)), 0, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 0, -r9*sin(t9(k)), -r10*sin(t10(k)), r11*sin(t11(k)), 0;
         0, 0, 0, 0, 0, 0, r9*cos(t9(k)), r10*cos(t10(k)), -r11*cos(t11(k)), 0;
         r3*sin(t3(k)), 0, -r5*sin(t5(k)), 0, 0, 0, 0, 0, 0, 0;
         -r3*cos(t3(k)), 0, r5*cos(t5(k)), 0, 0, 0, 0, 0, 0, 0;];

    
    B = [-r11*cos(t11(k))*dt11(k)^2+r12*cos(t12(k))*dt12(k)^2+r2*cos(t2(k))*dt2(k)^2+r2*sin(t2(k))*ddt2(k);
         -r11*sin(t11(k))*dt11(k)^2+r12*sin(t12(k))*dt12(k)^2+r2*sin(t2(k))*dt2(k)^2-r2*cos(t2(k))*ddt2(k);
         r7*cos(t7(k))*dt7(k)^2-r8*cos(t8(k))*dt8(k)^2+r9*cos(t9(k))*dt9(k)^2;
         r7*sin(t7(k))*dt7(k)^2-r8*sin(t8(k))*dt8(k)^2+r9*sin(t9(k))*dt9(k)^2;
         r4*cos(t4(k))*dt4(k)^2+r5*cos(t5(k))*dt5(k)^2-r6*cos(t6(k))*dt6(k)^2-r7*cos(t7(k))*dt7(k)^2;
         r4*sin(t4(k))*dt4(k)^2+r5*sin(t5(k))*dt5(k)^2-r6*sin(t6(k))*dt6(k)^2-r7*sin(t7(k))*dt7(k)^2;
         r9*cos(t9(k))*dt9(k)^2+r10*cos(t10(k))*dt10(k)^2-r11*cos(t11(k))*dt11(k)^2;
         r9*sin(t9(k))*dt9(k)^2+r10*sin(t10(k))*dt10(k)^2-r11*sin(t11(k))*dt11(k)^2;
         -r3*cos(t3(k))*dt3(k)^2+r5*cos(t5(k))*dt5(k)^2+r2*cos(t2(k))*dt2(k)^2+r2*sin(t2(k))*ddt2(k);
         -r3*sin(t3(k))*dt3(k)^2+r5*sin(t5(k))*dt5(k)^2+r2*sin(t2(k))*dt2(k)^2-r2*cos(t2(k))*ddt2(k);
         ];
    


     
    x = A\B;
    
    % results
    ddt3(k) = x(1);
    ddt4(k) = x(2);
    ddt5(k) = x(3);
    ddt6(k) = x(4);
    ddt7(k) = x(5);
    ddt8(k) = x(6);
    ddt9(k) = x(7);
    ddt10(k) = x(8);
    ddt11(k) = x(9);
    ddt12(k) = x(10);
    
    % change the initial values by the new calculated initial values for
    % the next iteration step.
    t_init = [t3(k) + Ts*dt3(k), t4(k) + Ts*dt4(k), t5(k) + Ts*dt5(k), t6(k) + Ts*dt6(k), t7(k) + Ts*dt7(k),
              t8(k) + Ts*dt8(k), t9(k) + Ts*dt9(k), t10(k) + Ts*dt10(k), t11(k) + Ts*dt11(k), t12(k) + Ts*dt12(k)]; 
    
end

% *** create movie ***

% point P = fixed
A = 0;
% point S = fixed
B = r1*exp(j*t1);
% define which positions we want as frames in our movie
frames = 100;    % number of frames in movie
delta = floor(t_size/frames); % time between frames
index_vec = [1:delta:t_size]';

% Create a window large enough for the whole mechanisme in all positions, to prevent scrolling.
% This is done by plotting a diagonal from (x_left, y_bottom) to (x_right, y_top), setting the
% axes equal and saving the axes into "movie_axes", so that "movie_axes" can be used for further
% plots.
x_left = -r9;
y_bottom = -( r7 + r6);
x_right = r1*cos(t1);
y_top = r11;

figure(10)
hold on
plot([x_left, x_right], [y_bottom, y_top]);
axis equal;
movie_axes = axis;   %save current axes into movie_axes

% draw and save movie frame
for m=1:length(index_vec)
    index = index_vec(m);
    C = B + r2*exp(j*t2(index));
    H = C + r12*exp(j*(t12(index)));
    F = A + r7 * exp(j*(t7(index)+pi));
    E = F + r6 * exp(j*(t6(index)+pi));
    D = E + r4 * exp(j*(t4(index)));
    G = A + r9 * exp(j*(t9(index)));
    
    loop1 = [A, B, C, H, A, F, E, D, A, B, C, D, A, G, F, G, H];
    
    figure(10)
    clf
    hold on
    plot(real(loop1),imag(loop1),'-o')
    
    axis(movie_axes);     % set axes as in movie_axes
    Movie(m) = getframe;  % save frame to a variable Film
end

% save movie
save ghassaei_movie Movie
close(10)

% #### PLOTS ####

if plot_kin % only plot when asked (option  for time purpose)
    index = 1; 
    
    A = 0;
    B = r1*exp(j*t1);
    C = B + r2*exp(j*t2(index));
    H = C + r12*exp(j*(t12(index)));
    F = A + r7 * exp(j*(t7(index)+pi));
    E = F + r6 * exp(j*(t6(index)+pi));
    D = E + r4 * exp(j*(t4(index)));
    G = A + r9 * exp(j*(t9(index)));
    
    
    figure
    assembly=[A, B, C, H, A, F, E, D, A, B, C, D, A, G, F, G, H];
    plot(real(assembly),imag(assembly),'ro-')
    xlabel('[m]')
    ylabel('[m]')
    title('assembly')
    axis equal
    
    %###### corner plot #####
    figure
    subplot(11,1,1)
    plot(t,t2)
    ylabel('\phi_2 [rad]')
    
    subplot(11,1,2)
    plot(t,t3)
    ylabel('\phi_3 [rad]')
    
    subplot(11,1,3)
    plot(t,t4)
    ylabel('\phi_4 [rad]')
    
    subplot(11,1,4)
    plot(t,t5)
    ylabel('\phi_5 [rad]')
    
    subplot(11,1,5)
    plot(t,t6)
    ylabel('\phi_6 [rad]')
    
    subplot(11,1,6)
    plot(t,t7)
    ylabel('\phi_7 [rad]')
    
    subplot(11,1,7)
    plot(t,t8)
    ylabel('\phi_8 [rad]')
    
    subplot(11,1,8)
    plot(t,t9)
    ylabel('\phi_9 [rad]')
    
    subplot(11,1,9)
    plot(t,t10)
    ylabel('\phi_{10} [rad]')
    
    subplot(11,1,10)
    plot(t,t11)
    ylabel('\phi_{11} [rad]')
    
    subplot(11,1,11)
    plot(t,t12)
    ylabel('\phi_{12} [rad]')
    
    %###### corner velocity plot #####
    figure 
     subplot(11,1,1)
    plot(t,dt2)
    ylabel('d\phi_2 [rad/s]')
    
    subplot(11,1,2)
    plot(t,dt3)
    ylabel('d\phi_3 [rad/s]')
    
    subplot(11,1,3)
    plot(t,dt4)
    ylabel('d\phi_4 [rad/s]')
    
    subplot(11,1,4)
    plot(t,dt5)
    ylabel('d\phi_5 [rad/s]')
    
    subplot(11,1,5)
    plot(t,dt6)
    ylabel('d\phi_6 [rad/s]')
    
    subplot(11,1,6)
    plot(t,dt7)
    ylabel('d\phi_7 [rad/s]')
    
    figure
    
    subplot(11,1,7)
    plot(t,dt8)
    ylabel('d\phi_8 [rad/s]')
    
    subplot(11,1,8)
    plot(t,dt9)
    ylabel('d\phi_9 [rad/s]')
    
    subplot(11,1,9)
    plot(t,dt10)
    ylabel('d\phi_{10} [rad/s]')
    
    subplot(11,1,10)
    plot(t,dt11)
    ylabel('d\phi_{11} [rad/s]')
    
    subplot(11,1,11)
    plot(t,dt12)
    ylabel('d\phi_{12} [rad/s]')
    
    %###### corner acceleration plot #####
    figure 
    subplot(11,1,1)
    plot(t,ddt2)
    ylabel('dd\phi_2 [rad/s^2]')
    
    subplot(11,1,2)
    plot(t,ddt3)
    ylabel('dd\phi_3 [rad/s^2]')
    
    subplot(11,1,3)
    plot(t,ddt4)
    ylabel('dd\phi_4 [rad/s^2]')
    
    subplot(11,1,4)
    plot(t,ddt5)
    ylabel('dd\phi_5 [rad/s^2]')
    
    subplot(11,1,5)
    plot(t,ddt6)
    ylabel('dd\phi_6 [rad/s^2]')
    
    subplot(11,1,6)
    plot(t,ddt7)
    ylabel('dd\phi_7 [rad/s^2]')
    
    figure
    subplot(11,1,7)
    plot(t,ddt8)
    ylabel('dd\phi_8 [rad/s^2]')
    
    subplot(11,1,8)
    plot(t,ddt9)
    ylabel('dd\phi_9 [rad/s^2]')
    
    subplot(11,1,9)
    plot(t,dt10)
    ylabel('dd\phi_{10} [rad/s^2]')
    
    subplot(11,1,10)
    plot(t,ddt11)
    ylabel('dd\phi_{11} [rad/s^2]')
    
    subplot(11,1,11)
    plot(t,ddt12)
    ylabel('dd\phi_{12} [rad/s^2]')
end
