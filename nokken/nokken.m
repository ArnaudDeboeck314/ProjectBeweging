% Seppe Vilain, Arnaud 

clear
close all

%figuren plotten 
singlerise=0;
multirise=1;
contactkracht=0;

%data inladen uit de gegevens geproduceerd met het matcam.m file
excOutput=load('exc_9mm_veer_17,25.mat');
output=load('exc_geen_veer_17,25.mat');


zeta=0.071;  % dempingscoefficient uit gegevens

%##     nokken profiel  ##
%##     vermogen        ##
%##     Vliegwiel       ##



%##     single rise     ##

% afleiding van formules staan in het verslag !

% bepalen van de kritische heffing
wdegree=excOutput.w*180/pi; %omzetten van de hoeksnelheid naar [°/s]

b1=45;  %beta eerste cycloïde
b2=60;  %beta tweede cycloïde
b3=105; %beta derde cycloïde

t11=b1/wdegree; %eerste cycloïde 
t12=b2/wdegree; %tweede cycloïde
t13=b3/wdegree; %derde cycloïde

h1=10e-3;   %S1 in meter (S1 = 10 mm)
h2=20e-3;   %S2 in meter (S2 = 10 mm)
h3=-30e-3;   %S3 in meter (S3 = 10 mm)

stijging1=h1/b1;
stijging2=h2/b2;    %grootste
stijging3=h3/b3;

t1=t12;         %tijd van de grootste heffing
beta=b2/180*pi; %beta van de grootste stijging in radialen 

% bereken van equivalente volgerstijfheid k_f
lambda=0.75/zeta;
tn=t1/lambda;
omegan=(2*pi)/tn;
kf=excOutput.mass*omegan^2;

% numerieke simulatie van de single rise (formules uit de cursus)
teller = (2*pi*lambda)^2;
noemer = [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];
sys = tf(teller, noemer); 


tau_single = 0:1/6000:5-1/6000; %array met de waarden tussen 1 en 5-1/6000 en met step size 1/6000 [1x30000]

% theta = tau - (1/(2pi)) * sin(2*pi*tau) met theta(tau>=0) = 1
theta= rectangularPulse(-0.1,1,tau_single).*((excOutput.w*tau_single*t1/beta) - (1/(2*pi))*sin(2*pi*(excOutput.w*tau_single*t1/beta))) + heaviside(tau_single-1);

theta0=0;
theta_dot0=0;

[A,B,C,D]=tf2ss(teller,noemer);

X0=[1/C(2)*theta_dot0; 1/C(2)*theta0];

gamma_single=lsim(A,B,C,D,theta,tau_single,X0);

difference=theta(:)-gamma_single(:);

if singlerise
    figure('Name','Single rise')
    subplot(3,1,1)
    plot(tau_single,theta)
    title('Input excitatie')
    xlabel('\tau [-]')
    ylabel('\theta [-]')
    axis([0 2 0 1])
    subplot(3,1,2)
    plot(tau_single,gamma_single)
    title('Output responsie')
    xlabel('\tau [-]')
    ylabel('\gamma [-]')
    axis([0 2 0 1])
    subplot(3,1,3)
    plot(tau_single,difference)
    title('Verschil tussen excitatie en responsie')
    xlabel('\tau [-]')
    ylabel('\theta-\gamma [-]')
    axis([0 2 -inf inf])
end

% benaderende analyse
N = 3;
Q = (2*pi)^2;
Amp = Q/((2*pi*lambda)^N);
[yup,ylow]=envelope(difference,1,'peak');
yexp = Amp*exp(-zeta*(2*pi*lambda)*(tau_single-1));

if singlerise
    figure('Name','Benaderende singe rise')
    subplot(2,1,1)
    plot(tau_single,difference,tau_single,yup,tau_single,ylow,tau_single,Amp)
    title('Single Rise Envelope')
    xlabel('\tau [-]')
    ylabel('\theta-\gamma [-]')
    xlim([1 3])
    ylim([-0.0002 0.0002])
    legend('excitatie-responsie','bovenenvelop','onderenvelop')
    subplot(2,1,2)
    plot(tau_single,difference,tau_single,Amp,tau_single, yexp)
    title('Exponentieel omhullende')
    xlabel('\tau [-]')
    ylabel('\theta-\gamma [-]')
    xlim([1 3])
    ylim([-0.0002 0.0002])
    legend('excitatie-responsie','exponentieel omhullende')
end

%##     Multi Rise      ##

T=0.25;
lambdatilde=T/tn;

teller=(2*pi*lambdatilde)^2;
noemer=[1, 2*zeta*2*pi*lambdatilde, (2*pi*lambdatilde)^2];
sys=tf(teller, noemer);

tau_multi = 0:1/36000:25-1/36000;
theta=repmat(output.S/30,1,25);

gamma_multi=lsim(sys,theta,tau_multi);

difference=theta(:)-gamma_multi(:);

verschil=zeros(24*36000,1);
for i=1:24*36000
    verschil(i)=difference(i+36000)-difference(i);
end

if multirise
    figure('Name','Multi rise')
    subplot(311)
    plot(tau_multi,theta)
    title('Input excitatie')
    xlabel('\tau')
    ylabel('\theta')
    xlim([24 25])
    subplot(312)
    plot(tau_multi,gamma_multi)
    title('Output responsie')
    xlabel('\tau')
    ylabel('\gamma')
    xlim([24 25])
    subplot(313)
    plot(tau_multi,difference);
    title('Verschil tussen excitatie en responsie')
    xlabel('\tau')
    ylabel('\gamma-\theta')
    xlim([24 25])
    
    figure('Name','Verschil tussen perioden')
    plot(tau_multi(1:length(verschil)),verschil)
    title('Verschil tussen perioden')
    xlabel('\tau')
    ylabel('[\gamma(\tau+1)-\theta(\tau+1)]-[\gamma(\tau)-\theta(\tau)]')

    figure('Name','Verschil tussen Single - en Multi Rise')
    plot(tau_single(1:6000),1/3+2/3*gamma_single(1:6000)-gamma_multi(873001:879000))
    title('Verschil tussen single en multi rise (90°-150°)')
    xlabel('\tau [-]')
    ylabel('\gamma_{single} - \gamma_{multi} [-]')
end


%##     Contactkracht    ##

F=zeros(1,36000);
N=zeros(1,36000);
for i=1:36000
    Nextra(i)=-3.299*10^7*(gamma_multi(i+24*36000)-theta(i+24*36000))*0.03/cos(excOutput.pressure_angle(i));
    N(i)=excOutput.normalforce_tot(i)+Nextra(i);
end

if contactkracht
    figure('Name','Contactkracht')
    plot(Nextra)
    hold on
    plot(N)
    hold off
    title('Influence of N_{extra} on N_{total}')
    xlim([0 3.6e4])
    ylabel('Normal Force [N]')
    legend('F_{extra}','N_{total}')
end





