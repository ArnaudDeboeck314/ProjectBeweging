load('exc_geen_veer_17,25.mat');

P1 = zeros(1,36000);

for i = 1:36000
    P1(i) = normalforce_tot(i)*sin(pressure_angle(i))*(60 + S(i))*w;
end

load('exc_9mm_veer_17,25.mat');

P2 = zeros(1,36000);

for i = 1:36000
    P2(i) = normalforce_tot(i) * cos(pressure_angle(i)) * (tan(pressure_angle(i)) * (sqrt((60^2) - (exc^2)) + S(i)) + exc) * w;
end

P = mean(P1) + theta*0;

M_gem = P/w;

M_L = P1/w;

M = M_gem-M_L;

En = 0;

delta_theta = theta(2) - theta(1);

for i = 5132:14221
    En = En + M(i)*delta_theta;
end

I = abs(En)/(0.1*w^2)
