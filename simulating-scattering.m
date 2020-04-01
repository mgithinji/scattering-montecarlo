%% Morgan Githinji (ENEE313H) - Monte Carlo Simulation

% declaring varibles
SumE = [];
SumV = [];

% defining constants
T = 300; % units: K
q = 1.6e-19;% units: C
kB_eV = 8.617e-5; % units: eV/K
kB = 1.381e-23; % units: J/K
h_bar = 1.05457e-34; % units: Js
m0 = 9.1e-31; % units: kg
m_l = 0.98*m0; % units: kg
m_t = 0.19*m0; % units: kg
alpha = 0.5/q; % units: 1/J
rho = 2330; % units: kg/m^3
D_ac = 9.0*q; % units: J 
v_s = 9e3; % units: m/s
Gamma = 2e14; % units: 1/s
dt = 1e-15; % units: s

% calculating m_d and effective mass
m = 1/((1/3)*((1/m_l)+(2/m_t)));
m_d = (m_l*(m_t^2))^(1/3);

% parameters
Zf = 4;
Zg = 1;

DtK_g1_eV = 0.5e8; % units: eV/cm
DtK_g1 = (0.5e8)*q*100; % units: J/m
DtK_g2_eV = 0.8e8; % units: eV/cm
DtK_g2 = (0.8e8)*q*100; % units: J/m
DtK_g3_eV = 11e8; % units: eV/cm
DtK_g3 = (11e8)*q*100; % units: J/m
DtK_f1_eV = 0.3e8; % units: eV/cm
DtK_f1 = (0.3e8)*q*100; % units: J/m
DtK_f2_eV = 2.0e8; % units: eV/cm
DtK_f2 = (2.0e8)*q*100; % units: J/m
DtK_f3_eV = 2.0e8; % units: eV/cm
DtK_f3 = (2.0e8)*q*100; % units: J/m

T_op_g1 = 140; % units: K
T_op_g2 = 215; % units: K
T_op_g3 = 720; % units: K
T_op_f1 = 220; % units: K
T_op_f2 = 550; % units: K
T_op_f3 = 685; % units: K

w_op_g1 = kB*T_op_g1/h_bar; % units: 1/s
w_op_g2 = kB*T_op_g2/h_bar; % units: 1/s
w_op_g3 = kB*T_op_g3/h_bar; % units: 1/s
w_op_f1 = kB*T_op_f1/h_bar; % units: 1/s
w_op_f2 = kB*T_op_f2/h_bar; % units: 1/s
w_op_f3 = kB*T_op_f3/h_bar; % units: 1/s

N_op_g1 = 1/(exp(h_bar*w_op_g1/(kB*T))-1);
N_op_g2 = 1/(exp(h_bar*w_op_g2/(kB*T))-1);
N_op_g3 = 1/(exp(h_bar*w_op_g3/(kB*T))-1);
N_op_f1 = 1/(exp(h_bar*w_op_f1/(kB*T))-1);
N_op_f2 = 1/(exp(h_bar*w_op_f2/(kB*T))-1);
N_op_f3 = 1/(exp(h_bar*w_op_f3/(kB*T))-1);

%% simulation loop
F_mag = [10^5 10^6 10^7 10^9 10^13 10^18 10^20];
a = 1;
F = F_mag(a); % units: V/m
while F < F_mag(7)
    F = F_mag(a); % units: V/m
    
    % initializing
    E = (3/2)*kB*T;
    k_mag = sqrt(2*m*E*(1+alpha*E))/h_bar;
    kx = k_mag/sqrt(3);
    ky = k_mag/sqrt(3);
    kz = k_mag/sqrt(3);
    
    % initializing other stuff
    Fx = F/sqrt(3); % units: V/m
    Fy = F/sqrt(3); % units: V/m
    Fz = F/sqrt(3); % units: V/m
    
    % allocating arrays for Energy and Velocity
    Energy = zeros(1,6e5);
    Velocity = zeros(1,6e5);
    Velocity_x = zeros(1,6e5);
    Velocity_y = zeros(1,6e5);
    Velocity_z = zeros(1,6e5);
    
    i = 1; % updates our arrays to be put in histogram
    time = 1; % keeping track of time
    num_flights = 100000; % this is my "max time"
    
    while time < num_flights
        if i~=1
            Energy (i) = (sqrt(1+4*alpha*(((h_bar^2)*k_mag^2)/(2*m)))-1)/(2*alpha);
            E = Energy (i);
            Velocity_x(i) = (h_bar*kx/m)/(1+(2*alpha*E));
            Velocity_y(i) = (h_bar*ky/m)/(1+(2*alpha*E));
            Velocity_z(i) = (h_bar*kz/m)/(1+(2*alpha*E));
            Velocity(i) = h_bar*dot([Velocity_x(i),Velocity_y(i),Velocity_z(i)],[1/sqrt(3),1/sqrt(3),1/sqrt(3)])/(m*(1+2*alpha*E));
            i = i + 1;
        end
        
        % getting random flight time
        r1 = rand(); % uniformly random number from 0 to 1
        tau = -log(r1)/Gamma; % random flight time
        
        while tau > dt
            tau = tau - dt;
            kx = kx - (q*Fx/h_bar)*dt;
            ky = ky - (q*Fy/h_bar)*dt;
            kz = kz - (q*Fz/h_bar)*dt;
            k_mag = sqrt(kx^2 + ky^2 + kz^2);
            Energy (i) = (sqrt(1+4*alpha*(((h_bar^2)*k_mag^2)/(2*m)))-1)/(2*alpha);
            E = Energy (i);
            Velocity_x(i) = (h_bar*kx/m)/(1+(2*alpha*E));
            Velocity_y(i) = (h_bar*ky/m)/(1+(2*alpha*E));
            Velocity_z(i) = (h_bar*kz/m)/(1+(2*alpha*E));
            Velocity(i) = h_bar*dot([Velocity_x(i),Velocity_y(i),Velocity_z(i)],[-1/sqrt(3),-1/sqrt(3),-1/sqrt(3)])/(m*(1+2*alpha*E));
            i = i + 1;
        end
        
        % calculating scattering rates
        % acoustic scattering
        Sac = ((sqrt(2)*(m_d^(3/2))*(kB*T*D_ac^2))/(pi*(h_bar^4)*(v_s^2)*rho))*sqrt(E)*(1+2*alpha*E)*sqrt(1+alpha*E);
        % optical scattering
        E_p = E+h_bar*w_op_g1; % g1a
        Sop_g1a = (((DtK_g1^2)*(m_d^(3/2))*Zg)/(sqrt(2)*pi*rho*(h_bar^3)*w_op_g1))*N_op_g1*sqrt(E_p+alpha*(E_p)^2)*(1+2*alpha*E_p);
        E_p = E-h_bar*w_op_g1; % g1e
        Sop_g1e = (((DtK_g1^2)*(m_d^(3/2))*Zg)/(sqrt(2)*pi*rho*(h_bar^3)*w_op_g1))*(N_op_g1 + 1)*sqrt(E_p+alpha*(E_p)^2)*(1+2*alpha*E_p);
        Sop_g1e = isreal(Sop_g1e)*Sop_g1e; % zero out the imaginary stuff
        E_p = E+h_bar*w_op_g2; % g2a
        Sop_g2a = (((DtK_g2^2)*(m_d^(3/2))*Zg)/(sqrt(2)*pi*rho*(h_bar^3)*w_op_g2))*N_op_g2*sqrt(E_p+alpha*(E_p)^2)*(1+2*alpha*E_p);
        E_p = E-h_bar*w_op_g2; % g2e
        Sop_g2e = (((DtK_g2^2)*(m_d^(3/2))*Zg)/(sqrt(2)*pi*rho*(h_bar^3)*w_op_g2))*(N_op_g2 + 1)*sqrt(E_p+alpha*(E_p)^2)*(1+2*alpha*E_p);
        Sop_g2e = isreal(Sop_g2e)*Sop_g2e; % zero out the imaginary stuff
        E_p = E+h_bar*w_op_g3; % g3a
        Sop_g3a = (((DtK_g3^2)*(m_d^(3/2))*Zg)/(sqrt(2)*pi*rho*(h_bar^3)*w_op_g3))*N_op_g3*sqrt(E_p+alpha*(E_p)^2)*(1+2*alpha*E_p);
        E_p = E-h_bar*w_op_g3; % g3e
        Sop_g3e = (((DtK_g3^2)*(m_d^(3/2))*Zg)/(sqrt(2)*pi*rho*(h_bar^3)*w_op_g3))*(N_op_g3 + 1)*sqrt(E_p+alpha*(E_p)^2)*(1+2*alpha*E_p);
        Sop_g3e = isreal(Sop_g3e)*Sop_g3e; % zero out the imaginary stuff
        E_p = E+h_bar*w_op_f1; % f1a
        Sop_f1a = (((DtK_f1^2)*(m_d^(3/2))*Zf)/(sqrt(2)*pi*rho*(h_bar^3)*w_op_f1))*N_op_f1*sqrt(E_p+alpha*(E_p)^2)*(1+2*alpha*E_p);
        E_p = E-h_bar*w_op_f1; % f1e
        Sop_f1e = (((DtK_f1^2)*(m_d^(3/2))*Zf)/(sqrt(2)*pi*rho*(h_bar^3)*w_op_f1))*(N_op_f1 + 1)*sqrt(E_p+alpha*(E_p)^2)*(1+2*alpha*E_p);
        Sop_f1e = isreal(Sop_f1e)*Sop_f1e; % zero out the imaginary stuff
        E_p = E+h_bar*w_op_f2; % f2a
        Sop_f2a = (((DtK_f2^2)*(m_d^(3/2))*Zf)/(sqrt(2)*pi*rho*(h_bar^3)*w_op_f2))*N_op_f2*sqrt(E_p+alpha*(E_p)^2)*(1+2*alpha*E_p);
        E_p = E-h_bar*w_op_f2; % f2e
        Sop_f2e = (((DtK_f2^2)*(m_d^(3/2))*Zf)/(sqrt(2)*pi*rho*(h_bar^3)*w_op_f2))*(N_op_f2 + 1)*sqrt(E_p+alpha*(E_p)^2)*(1+2*alpha*E_p);
        Sop_f2e = isreal(Sop_f2e)*Sop_f2e; % zero out the imaginary stuff
        E_p = E+h_bar*w_op_f3; % f3a
        Sop_f3a = (((DtK_f3^2)*(m_d^(3/2))*Zf)/(sqrt(2)*pi*rho*(h_bar^3)*w_op_f3))*N_op_f3*sqrt(E_p+alpha*(E_p)^2)*(1+2*alpha*E_p);
        E_p = E-h_bar*w_op_f3; % f3e
        Sop_f3e = (((DtK_f3^2)*(m_d^(3/2))*Zf)/(sqrt(2)*pi*rho*(h_bar^3)*w_op_f3))*(N_op_f3 + 1)*sqrt(E_p+alpha*(E_p)^2)*(1+2*alpha*E_p);
        Sop_f3e = isreal(Sop_f3e)*Sop_f3e; % zero out the imaginary stuff
        
        % get values for S1...S13
        S1 = Sac;
        S2 = Sop_g1a;
        S3 = Sop_g1e;
        S4 = Sop_g2a;
        S5 = Sop_g2e;
        S6 = Sop_g3a;
        S7 = Sop_g3e;
        S8 = Sop_f1a;
        S9 = Sop_f1e;
        S10 = Sop_f2a;
        S11 = Sop_f2e;
        S12 = Sop_f3a;
        S13 = Sop_f3e;
        
        % get values for L1...L13
        L1 = S1;
        L2 = S1 + S2;
        L3 = S1 + S2 + S3;
        L4 = S1 + S2 + S3 + S4;
        L5 = S1 + S2 + S3 + S4 + S5;
        L6 = S1 + S2 + S3 + S4 + S5 + S6;
        L7 = S1 + S2 + S3 + S4 + S5 + S6 + S7;
        L8 = S1 + S2 + S3 + S4 + S5 + S6 + S7 + S8;
        L9 = S1 + S2 + S3 + S4 + S5 + S6 + S7 + S8 + S9;
        L10 = S1 + S2 + S3 + S4 + S5 + S6 + S7 + S8 + S9 + S10;
        L11 = S1 + S2 + S3 + S4 + S5 + S6 + S7 + S8 + S9 + S10 + S11;
        L12 = S1 + S2 + S3 + S4 + S5 + S6 + S7 + S8 + S9 + S10 + S11 + S12;
        L13 = S1 + S2 + S3 + S4 + S5 + S6 + S7 + S8 + S9 + S10 + S11 + S12 + S13;
        
        flag = 0;
        % randomly choose scattering mechanism and update E' and k'
        r2 = rand(); % uniformly random number from 0 to 1
        if (r2>=0)&&(r2<=L1/Gamma) % acoustic
            E_p = E;
            k_p = sqrt(2*m*E_p*(1+alpha*E_p))/h_bar;
        elseif (r2>=L1/Gamma)&&(r2<=L2/Gamma)
            E_p = E+h_bar*w_op_g1; % g1a
            k_p = sqrt(2*m*E_p*(1+alpha*E_p))/h_bar;
        elseif (r2>=L2/Gamma)&&(r2<=L3/Gamma)
            E_p = E-h_bar*w_op_g1; % g1e
            k_p = sqrt(2*m*E_p*(1+alpha*E_p))/h_bar;
        elseif (r2>=L3/Gamma)&&(r2<=L4/Gamma)
            E_p = E+h_bar*w_op_g2; % g2a
            k_p = sqrt(2*m*E_p*(1+alpha*E_p))/h_bar;
        elseif (r2>=L4/Gamma)&&(r2<=L5/Gamma)
            E_p = E-h_bar*w_op_g2; % g2e
            k_p = sqrt(2*m*E_p*(1+alpha*E_p))/h_bar;
        elseif (r2>=L5/Gamma)&&(r2<=L6/Gamma)
            E_p = E+h_bar*w_op_g3; % g3a
            k_p = sqrt(2*m*E_p*(1+alpha*E_p))/h_bar;
        elseif (r2>=L6/Gamma)&&(r2<=L7/Gamma)
            E_p = E-h_bar*w_op_g3; % g3e
            k_p = sqrt(2*m*E_p*(1+alpha*E_p))/h_bar;
        elseif (r2>=L7/Gamma)&&(r2<=L8/Gamma)
            E_p = E+h_bar*w_op_f1; % f1a
            k_p = sqrt(2*m*E_p*(1+alpha*E_p))/h_bar;
        elseif (r2>=L8/Gamma)&&(r2<=L9/Gamma)
            E_p = E-h_bar*w_op_f1; % f1e
            k_p = sqrt(2*m*E_p*(1+alpha*E_p))/h_bar;
        elseif (r2>=L9/Gamma)&&(r2<=L10/Gamma)
            E_p = E+h_bar*w_op_f2; % f2a
            k_p = sqrt(2*m*E_p*(1+alpha*E_p))/h_bar;
        elseif (r2>=L10/Gamma)&&(r2<=L11/Gamma)
            E_p = E-h_bar*w_op_f2; % f2e
            k_p = sqrt(2*m*E_p*(1+alpha*E_p))/h_bar;
        elseif (r2>=L11/Gamma)&&(r2<=L12/Gamma)
            E_p = E+h_bar*w_op_f3; % f3a
            k_p = sqrt(2*m*E_p*(1+alpha*E_p))/h_bar;
        elseif (r2>=L12/Gamma)&&(r2<=L13/Gamma)
            E_p = E-h_bar*w_op_f3; % f3e
            k_p = sqrt(2*m*E_p*(1+alpha*E_p))/h_bar;
        else
            E_p = E;
            k_p = k_mag;
            flag = 1;
        end
        
        if flag == 0
            % calculating the values of kx' ky' and kz'
            r3 = rand(); % uniformly random number from 0 to 1
            r4 = rand(); % uniformly random number from 0 to 1
            phi_p = 2*pi*r3;
            theta_p = acos(1-2*r4);
            kx_p = k_p*sin(theta_p)*cos(phi_p);
            ky_p = k_p*sin(theta_p)*sin(phi_p);
            kz_p = k_p*cos(theta_p);
            % update state E = E' and k = k'
            E = E_p;
            kx = kx_p;
            ky = ky_p;
            kz = kz_p;
            k_mag = sqrt(kx^2 + ky^2 + kz^2);
        end
        time = time + 1;
    end
    % chopping unused array spaces
    Energy = Energy(1:i-1);
    Velocity = Velocity(1:i-1);
    
    % plotting histograms
    figure
    histogram(Energy/q)
    str = sprintf('Frequency of Energy for F = %d V/m',F);
    title(str)
    xlabel('Energy (eV)')
    figure
    histogram(Velocity)
    str = sprintf('Frequency of Velocity for F = %d V/m',F);
    title(str)
    xlabel('Velocity (m/s)')
    
    incrementing variable that changes the field magnitude
    a = a + 1;
    % appending array for avg energy and velocity
    SumE = [SumE mean(Energy)];
    SumV = [SumV mean(Velocity)];
end

%% plotting averages
figure
loglog(F_mag,SumE/q); 
title('Average Energy as Field Varies')
xlabel('Field (V/m)')
ylabel('Average Energy (eV)')
figure
loglog(F_mag,SumV*100);
title('Average Velocity as Field Varies')
xlabel('Field (V/cm)')
ylabel('Average Velocity (eV)')
%set(gca,'YScale','log')