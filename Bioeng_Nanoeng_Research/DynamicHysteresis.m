function dhys
global sigma epsilong f

mu0 = 4*pi*1e-7;
Kb = 1.38e-23;
T = 300;

fi = [170 210 250 320 380];          % can change perimeter

for ifx = 1:5

f = fi(ifx)*1.e3;        %Hz

Hi = [10 15 20 25 30 40 50 60];          % can change perimeter

for ixx = 1:8
    
    Keff = 13500;   %J/m3
    r = 15e-9;      %m
    V = 4/3*pi*r^3; %m3

    Ms = 440000;    %A/m

    H0 = Hi(ixx)*1.e3;    %A/m

    sigma = Keff*V/(Kb*T);
    epsilong = H0*mu0*Ms*V/(Kb*T);
    Hk = 2*Keff/(Ms*mu0);

    nT = 30;
    options = odeset('RelTol',1e-5,'AbsTol',1e-6,'MaxStep',0.1/f);
    [X,P] = ode23t(@dpdt,[0 nT/f],0.5,options);

    M = Ms*(2*P-1);
    H = H0*cos(X*2*pi*f);

    %Search for index of the last nx cycle
    nx = 10;
    nx2 = length(M);
    t1 = (nT-nx)/f;
    i = 1;
    while X(i) < t1;
        i = i+1;
    end
    nx1 = i;

    %Integrate for the area
    A = 0;
    for i = nx1:nx2-1
        A = A + (H(i)-H(i+1))*(M(i+1)+M(i))/2;
    end

    [fi(ifx) Hi(ixx)]
    SAR = mu0*f*A/nx/5170000

    figure
    subplot(3,1,1);
    plot(X, H)

    subplot(3,1,2);
    plot(X, M);

    subplot(3,1,3);
    plot(X, P);

    lx = length(X);
    lxh = double(int16(lx/2));

    figure
    plot(H(lxh:lx)*mu0,M(lxh:lx))
    xlim([-0.05 0.05])

    px = [num2str(fi(ifx)) 'kHz\'];
    fx = num2str(Hi(ixx));
    save([px 'X' fx],'X');
    save([px 'H' fx],'H');
    save([px 'M' fx],'M');
    save([px 'P' fx],'P');

end

end

end

function dy=dpdt(t,y)
global sigma epsilong f

epsilongT = epsilong*cos(t*2*pi*f);

E1 = - epsilongT;
E2 = epsilongT;

atheta = -epsilongT/(2*sigma);
if atheta > 1
    atheta = 1;
else if atheta < -1
        atheta = -1;
    end
end

theta3 = acos(atheta);
E3 = sigma*(sin(theta3))^2-epsilongT*cos(theta3);

v1 = 1.e12*exp(E1-E3);
v2 = 1.e12*exp(E2-E3);

dy = (1-y)*v2-y*v1;

end
          
