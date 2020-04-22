
%% AIRCRAFT FLIGHT CONDITION 'LANDING'.
V     = 125.675;
m     = 5535;
twmuc = 89.7;
KY2   = 1.4996;
c     = 2.057;
S     = 30;
lh    = 5.5;
g     = 9.80665;

%% TURBULENCE PARAMETERS
disp(' ');
sigma = 1;
Lg    = 150;

sigmaug_V = sigma/V;
sigmaag   = sigma/V;

%% AIRCRAFT SYMMETRIC AERODYNAMIC DERIVATIVES : 
CX0 = 0.0000;     CZ0  =-0.2292;     Cm0  =  0.0000;
CXu =-0.0032;     CZu  =-0.4592;     Cmu  =  0.0236;
CXa = 0.1692;     CZa  =-5.7874;     Cma  = -0.7486;
CXq = -0.0450;    CZq  =-4.5499;     Cmq  = -7.4647;
CXd = 0.0000;     CZd  =-0.5798;     Cmd  = -1.4440;
CXfa= 0.0000;     CZfa =-0.3980;     Cmfa = -4.2255;
                  CZfug= 0.0000;     Cmfug= -Cm0*lh/c;
                  CZfag= CZfa-CZq;   Cmfag=  Cmfa-Cmq;

%% CALCULATION OF AIRCRAFT SYMMETRIC STABILITY DERIVATIVES
xu   = (V/c)*(CXu/twmuc);
xa   = (V/c)*(CXa/twmuc);
xt   = (V/c)*(CZ0/twmuc);
xq   = 0;
xd   = (V/c)*(CXd/twmuc);
xug  = xu;
xfug = 0;
xag  = xa;
xfag = 0;

zu   = (V/c)*( CZu/(twmuc-CZfa));
za   = (V/c)*( CZa/(twmuc-CZfa));
zt   = (V/c)*(-CX0/(twmuc-CZfa));
zq   = (V/c)*((CZq+twmuc)/(twmuc-CZfa));
zd   = (V/c)*( CZd/(twmuc-CZfa));
zug  = zu;
zfug = (V/c)*( CZfug/(twmuc-CZfa));
zag  = za;
zfag = (V/c)*( CZfag/(twmuc-CZfa));

mu   = (V/c)*(( Cmu+CZu*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
ma   = (V/c)*(( Cma+CZa*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mt   = (V/c)*((-CX0*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mq   = (V/c)*(Cmq+Cmfa*(twmuc+CZq)/(twmuc-CZfa))/(twmuc*KY2);
md   = (V/c)*((Cmd+CZd*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mug  = mu;
mfug = (V/c)*(Cmfug+CZfug*Cmfa/(twmuc-CZfa))/(twmuc*KY2);
mag  = ma;
mfag = (V/c)*(Cmfag+CZfag*Cmfa/(twmuc-CZfa))/(twmuc*KY2);

%% Define time basis
dt = 0.01; T = 200;
t  = [0:dt:T];
N = length(t);
fs = 1/dt;

%% Uncontrolled system
% STATE- AND INPUT MATRICES
Au=[xu xa xt 0    xug                  xag       0;
   zu za zt zq   zug-zfug*V/Lg*(c/V)  zag       zfag*(c/V);
   0  0  0  V/c  0                    0         0;
   mu ma mt mq   mug-mfug*V/Lg*(c/V)  mag       mfag*(c/V);
   0  0  0  0   -V/Lg                 0         0;
   0  0  0  0    0                    0         1;
   0  0  0  0    0                   -(V/Lg)^2 -2*V/Lg];

Bu=...
 [xd 0                                 0;
  zd zfug*(c/V)*sigmaug_V*sqrt(2*V/Lg) zfag*(c/V)*sigmaag*sqrt(3*V/Lg);
  0  0                                 0;
  md mfug*(c/V)*sigmaug_V*sqrt(2*V/Lg) mfag*(c/V)*sigmaag*sqrt(3*V/Lg);
  0  sigmaug_V*sqrt(2*V/Lg)            0;
  0  0                                 sigmaag*sqrt(3*V/Lg);
  0  0                                 (1-2*sqrt(3))*sigmaag*sqrt((V/Lg)^3)];

Cuv = zeros(1,7);
Duv = zeros(1,3);
sysu = ss(Au,Bu,Cuv,Duv);

disp('Eigenvalues uncontrolled aircraft model')
disp(eig(Au))
figure(1)
pzplot(sysu)
a = findobj(gca,'type','line');
for i = 1:length(a)
    set(a(i),'markersize',12) %change marker size
    set(a(i), 'linewidth',2)  %change linewidth
end

%% Controlled system
% GAIN FACTORS FOR AUTOPILOT CHAPTER 7 AND DEFINITION OF MATRIX At
% (approximately "Dead Beat" damping)
Kt = -0.21; Kq = 0;      % gains on "theta" and "q"
K  = [0 0 Kt Kq 0 0 0];  % feedback matrix
Ac = Au-Bu(:,1)*K;         % new A matrix = (A - BK) because of feedback
Bc = Bu;
Cc = zeros(1,7);
Dc = zeros(1,3);
sysc = ss(Ac,Bc,Cc,Dc);

disp('Eigenvalues of uncontrolled aircraft model')
disp(eig(Ac))
figure(2)
pzplot(sysc)
a = findobj(gca,'type','line');
for i = 1:length(a)
    set(a(i),'markersize',12) %change marker size
    set(a(i), 'linewidth',2)  %change linewidth
end

%% Reduced system
% Ar = [CZa+(CZfa-2*twmuc)*(c/V) CZq+2*twmuc;
%       Cma+Cmfa*(c/V) Cmq-2*twmuc*KY2*(c/V)];
% Br = [CXa CZfag;
%       Cma Cmfag];

Ar = Ac([2 4 5 6 7],[2 4 5 6 7]);
Br = Bc([2 4 5 6 7],:);
Cr = zeros(1,5);
Dr = zeros(1,3);
sysr = ss(Ar,Br,Cr,Dr);

disp('Eigenvalues for reduced aircraft model')
disp(eig(Ar))
figure(3)
pzplot(sysr)
a = findobj(gca,'type','line');
for i = 1:length(a)
    set(a(i),'markersize',12) %change marker size
    set(a(i), 'linewidth',2)  %change linewidth
end

%% Time domain analyses Controlled system

% White noise input
delta_e = zeros(1,N);
w_1 = randn(1,N)/sqrt(dt);
w_3 = randn(1,N)/sqrt(dt);

u = [delta_e' w_1' w_3'];

% Defining C and D matrix with n_z   
Cc = eye(4,7);      Cc(5,:) = (V/g)*(Ac(3,:)-Ac(2,:));
Dc = zeros(5,3);    Dc(5,:) = (V/g)*(Bc(3,:)-Bc(2,:));

% Simulating
yc = lsim(ss(Ac,Bc,Cc,Dc),u,t);

% Plotting the results
figure(4)
plot(t,yc(:,1))
xlabel('time [sec]')
ylabel('Airspeed deviation [-]')

figure(5)
plot(t,yc(:,2))
xlabel('time [sec]')
ylabel('Angle of attack [rad]')

figure(6)
plot(t,yc(:,3))
xlabel('time [sec]')
ylabel('Pitch angle [rad]')

figure(7)
plot(t,yc(:,4))
xlabel('time [sec]')
ylabel('Pitch rate [rad/s]')

figure(8)
plot(t,yc(:,5))
xlabel('time [sec]')
ylabel('load factor [m/s^2]')

%% Time domain analyses Reduced system

% Simulating
Cr = eye(2,5); Dr = zeros(2,3);
yr = lsim(ss(Ar,Br,Cr,Dr),u,t);

%Plotting results
figure(9)
plot(t,yr(:,1))
xlabel('time [sec]')
ylabel('Angle of attack [rad]')

figure(10)
plot(t,yr(:,2))
xlabel('time [sec]')
ylabel('Pitch rate [rad/s]')

%% Spectral analysis Controlled system

%frequency axis
Nomega = 2001; w = logspace(-2,3,Nomega);

% Calculation of analytic PSD
mag = bode(Ac,Bc,Cc(1,:),Dc(1,:),2,w) + bode(Ac,Bc,Cc(1,:),Dc(1,:),3,w); Suuc = sum(mag.*mag,2);
mag = bode(Ac,Bc,Cc(2,:),Dc(2,:),2,w) + bode(Ac,Bc,Cc(2,:),Dc(2,:),3,w); Saac = sum(mag.*mag,2);
mag = bode(Ac,Bc,Cc(3,:),Dc(3,:),2,w) + bode(Ac,Bc,Cc(3,:),Dc(3,:),3,w); Sttc = sum(mag.*mag,2);
mag = bode(Ac,Bc,Cc(4,:),Dc(4,:),2,w) + bode(Ac,Bc,Cc(4,:),Dc(4,:),3,w); Sqqc = sum(mag.*mag,2);
mag = bode(Ac,Bc,Cc(5,:),Dc(5,:),2,w) + bode(Ac,Bc,Cc(5,:),Dc(5,:),3,w); Snnc = sum(mag.*mag,2);

% Calculation of experimental PSD
u_v = dt*fft(yc(:,1));       Pu_v = real((1/T)*u_v.*conj(u_v));
alpha = dt*fft(yc(:,2));     Palpha = real((1/T)*alpha.*conj(alpha));
theta = dt*fft(yc(:,3));     Ptheta = real((1/T)*theta.*conj(theta));
q = dt*fft(yc(:,4));         Pq = real((1/T)*q.*conj(q));
n_z = dt*fft(yc(:,5));       Pn_z = real((1/T)*n_z.*conj(n_z));

% Calculation of smoothed PSD
Psu_v = Pu_v;
Psalpha = Palpha;
Pstheta = Ptheta;
Psq = Pq;
Psn_z = Pn_z;

for k=2:length(Psalpha)-1
    Psu_v(k) = 0.25*Pu_v(k-1) + 0.5*Pu_v(k) + 0.25*Pu_v(k+1);
    Psalpha(k) = 0.25*Palpha(k-1) + 0.5*Palpha(k) + 0.25*Palpha(k+1);
    Pstheta(k) = 0.25*Ptheta(k-1) + 0.5*Ptheta(k) + 0.25*Ptheta(k+1);
    Psq(k) = 0.25*Pq(k-1) + 0.5*Pq(k) + 0.25*Pq(k+1);
    Psn_z(k) = 0.25*Pn_z(k-1) + 0.5*Pn_z(k) + 0.25*Pn_z(k+1);
end

% Plotting the result
omega = 2*pi*fs*(0:(N/2)-1)/N;

figure(11)
loglog(omega,Pu_v(1:round(N/2)-1),omega,Psu_v(1:round(N/2)-1),w,Suuc,'linewidth',1.5);
xlim(10.^[-1,2])
xlabel('omega [rad/s]')
ylabel('Suu [rad^2]')
legend('Experimental PSD','Smoothed experimental PSD','Analytic PSD','Location','southwest')

figure(12)
loglog(omega,Palpha(1:round(N/2)-1),omega,Psalpha(1:round(N/2)-1),w,Saac,'linewidth',1.5);
xlim(10.^[-1,2])
xlabel('omega [rad/s]')
ylabel('Saa [rad^2]')
legend('Experimental PSD','Smoothed experimental PSD','Analytic PSD','Location','southwest')

figure(13)
loglog(omega,Ptheta(1:round(N/2)-1),omega,Pstheta(1:round(N/2)-1),w,Sttc,'linewidth',1.5);
xlim(10.^[-1,2])
xlabel('omega [rad/s]')
ylabel('Stt [rad^2]')
legend('Experimental PSD','Smoothed experimental PSD','Analytic PSD','Location','southwest')

figure(14)
loglog(omega,Pq(1:round(N/2)-1),omega,Psq(1:round(N/2)-1),w,Sqqc,'linewidth',1.5);
xlim(10.^[-1,2])
xlabel('omega [rad/s]')
ylabel('Sqq [rad^2]')
legend('Experimental PSD','Smoothed experimental PSD','Analytic PSD','Location','southwest')

figure(15)
loglog(omega,Pn_z(1:round(N/2)-1),omega,Psn_z(1:round(N/2)-1),w,Snnc,'linewidth',1.5);
xlim(10.^[-1,2])
xlabel('omega [rad/s]')
ylabel('Snn [rad^2]')
legend('Experimental PSD','Smoothed experimental PSD','Analytic PSD','Location','southwest')

%% Spectral Analysis Reduced system

% Calculation of analytic PSD
mag = bode(Ar,Br,Cr(1,:),Dr(1,:),2,w) + bode(Ar,Br,Cr(1,:),Dr(1,:),3,w); Saar = sum(mag.*mag,2);
mag = bode(Ar,Br,Cr(2,:),Dr(2,:),2,w) + bode(Ar,Br,Cr(2,:),Dr(2,:),3,w); Sqqr = sum(mag.*mag,2);

% Calculation of experimental PSD
alpha = dt*fft(yr(:,1));    Palpha = real((1/T)*alpha.*conj(alpha));
q = dt*fft(yr(:,2));        Pq = real((1/T)*q.*conj(q));

% Calculation of smoothed PSD
Psalpha = Palpha;
Psq = Pq;

for k=2:length(Psalpha)-1
    Psalpha(k) = 0.25*Palpha(k-1) + 0.5*Palpha(k) + 0.25*Palpha(k+1);
    Psq(k) = 0.25*Pq(k-1) + 0.5*Pq(k) + 0.25*Pq(k+1);
end

% Plotting the result
figure(16)
loglog(omega,Palpha(1:round(N/2)-1),omega,Psalpha(1:round(N/2)-1),w,Saar,'linewidth',1.5);
xlim(10.^[-1,2])
xlabel('omega [rad/s]')
ylabel('Saa [rad^2]')
legend('Experimental PSD','Smoothed experimental PSD','Analytic PSD','Location','southwest')

figure(17)
loglog(omega,Pq(1:round(N/2)-1),omega,Psq(1:round(N/2)-1),w,Sqqr,'linewidth',1.5);
xlim(10.^[-1,2])
xlabel('omega [rad/s]')
ylabel('Sqq [rad^2]')
legend('Experimental PSD','Smoothed experimental PSD','Analytic PSD','Location','southwest')

%% Variances Controlled system

con_vartable = zeros(5,3);

% Calculation of analytic variance
dw = diff(w);
dw(length(dw)+1) = 0;

con_vartable(1,1) = sum(Suuc'.*dw)/pi;
con_vartable(2,1) = sum(Saac'.*dw)/pi;
con_vartable(3,1) = sum(Sttc'.*dw)/pi;
con_vartable(4,1) = sum(Sqqc'.*dw)/pi;
con_vartable(5,1) = sum(Snnc'.*dw)/pi;

% Calculation of variance with impulse response method
x0 = Bc(:,2) + Bc(:,3);
yci = lsim(ss(Ac,Bc,Cc,Dc),u,t,x0);

con_vartable(1,2) = cov(yci(:,1));
con_vartable(2,2) = cov(yci(:,2));
con_vartable(3,2) = cov(yci(:,3));
con_vartable(4,2) = cov(yci(:,4));
con_vartable(5,2) = cov(yci(:,5));

% Calculation of variance with var.m
con_vartable(1,3) = var(yc(:,1));
con_vartable(2,3) = var(yc(:,2));
con_vartable(3,3) = var(yc(:,3)); 
con_vartable(4,3) = var(yc(:,4));
con_vartable(5,3) = var(yc(:,5));

disp('Variances of controlled aircraft model')
disp(con_vartable)

%% Variances Reduced system

red_vartable = zeros(2,3);

% Calculation of analytic variance
red_vartable(1,1) = sum(Saar'.*dw)/pi;
red_vartable(2,1) = sum(Sqqr'.*dw)/pi;

% Calculation of variance with impulse response method
x0 = Br(:,2) + Br(:,3);
yri = lsim(ss(Ar,Br,Cr,Dr),u,t,x0);

red_vartable(1,2) = cov(yri(:,1));
red_vartable(2,2) = cov(yri(:,2));

% Calculation of variance with var.m
red_vartable(1,3) = var(yr(:,1));
red_vartable(2,3) = var(yr(:,2));

disp('Variances of reduced aircraft model')
disp(red_vartable)