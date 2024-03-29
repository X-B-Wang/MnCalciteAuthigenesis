%%
close all; clear all;
format short;

%% units
% weight: g
% concentration: 1 mol/m3 = 1e-3 mol/L = mmol/L = 1e-3 mmol/cm3 = 1 umol/cm3
% time: thousand years, kyr
% length: m

%% Parameter settings
% sediment properties
phi = 0.60; % porosity, 60%

% chemical properties
Raero = 1e-9 * 1e-3 * 60*60*24*365*1000; % dissolved oxygen penetration
RMMR = 5e-11 * 1e-3 * 60*60*24*365*1000; % MMR rate
Rox = 1e-4 * 1e-3 * 60*60*24*365*1000; % oxidation kinetics of Mn2+, 1e-5 ~ 1e-4 1/M/s
Rcal = 2e-6 * 1e-3 * 60*60*24*365*1000; % precipitation rate of calcite
Ksp = 3.36e-9 * 1e6; % calcite solubility product, mol2/L2

KMn = 5; % partition coefficient of Mn2+

DMn = 6.57e-10 * 60*60*24*365*1000; % diffusivity of Mn2+, m2/s
DO2 = 2.13e-9 * 60*60*24*365*1000;
DCO3 = 8.67e-10 * 60*60*24*365*1000;
DHCO3 = 1.11e-9 * 60*60*24*365*1000;
DH2CO3 = 1.79e-9 * 60*60*24*365*1000;
DCa = 7.53e-10 * 60*60*24*365*1000;

DMn = DMn / (1 - 2 * log(phi)); % diffusivity corrected by porosity
DO2 = DO2 / (1 - 2 * log(phi));
DDIC = mean(DCO3 + DHCO3 + DH2CO3);
DDIC = DDIC / (1 - 2 * log(phi));
DCa = DCa / (1 - 2 * log(phi));

% boundary conditions
O0 = 250e-3; % dissolved O2 conc, saturated=250 uM
DIC0 = 20e-3; % DIC conc, 0.6-138 uM, mean=37 uM, median=19 uM, set as 20 uM
Ca0 = 0.2; % Ca2+ conc, mean=0.38 mM, median=0.14 mM, set as 0.2 mM
Mn0 = 2 * 55*1e-6*1e2/3*800; % Mn-solid conc 'in porewater', Mn in sediment is 2 mM/kg~110 ppm
% Mn0 = 0.0007 * rho/(55+16)*(1-phi)/phi*1e3; % 0.07 wt.% MnO in UCC
CH0 = 800; % 3 wt.% TOC, (g/g)/(g/mol)*(g/cm3)=mol/cm3
SR = 0.03 * 2.5; % sedimentation rate, 5 cm/kyr
MMR_depth = 5e-3; % oxygen tolerance, where MMR starts


%% Correction zone
% modifying the parameters to study site
Mn0 = Mn0 * 0.6; % HCl-leaching, reactive Mn
Ca0 = Ca0 * 6;
DIC0 = DIC0 * 1;
CH0 = CH0 * 0.1;
O0 = O0 *0.08;


%% numerical variables
maxD = 10; dz = 0.01;
zz = 0:dz:maxD; zz = zz';
N = length(zz)-1; nii = 100;
eps = 1e-9; % 1 nM


%% Geochemical simulation
s = SR;
% dissolved oxygen: simple
O2 = O0 * exp((s - sqrt(s^2 + 4 * Raero * CH0 * DO2)) / (2 * DO2) * zz);
RR = RMMR + 0 * zz;
for i = 1:N+1
    if O2(i) > MMR_depth * 2
        RR(i) = 0;
    elseif O2(i) > MMR_depth
        RR(i) = interp1([1 2], [RMMR 0], O2(i)/MMR_depth);
    end
end


% MMR
Mns = Mn0 + 0 * zz;
CH = CH0 + 0 * zz;
for ii = 1:10
    CHI = 0 * zz;
    ch = CH .* RR / s;
    for i = 1:N+1
        CHI(i) = ISimpson(ch, i, dz);
    end
    Mns_prev = Mns;
    Mns = Mn0 * exp(-CHI);

    MnsI = 0 * zz;
    mns_o = Mns .* RR / 2 / s + Raero * O2 / s;
    for i = 1:N+1
        MnsI(i) = ISimpson(mns_o, i, dz);
    end
    CH_prev = CH;
    CH = CH0 * exp(-MnsI);
    if mean(abs(CH_prev - CH)) < eps && mean(abs(Mns_prev - Mns)) < eps
        break;
    end
    if ii == nii
        'MMR!'
    end
end

% DIC
DIC = DIC0 + 0 * zz;
Ca = Ca0 + 0 * zz;
Mn = 0 + 0 * zz;
for jj = 1:nii
    SSA = RR / 2 .* Mns .* CH + Raero * O2 .* CH;
    calciteFinal = 0;
    for ii = 1:nii
        % Thomas method
        G = -2 * DDIC/dz^2 * eye(N);
        for i=1:N-1
            G(i,i+1) = DDIC/dz^2 - s/2/dz;
            G(i+1,i) = DDIC/dz^2 + s/2/dz;
            G(i,i) = G(i,i) - Rcal*Ca(i+1) - KMn*Rcal*Mn(i+1);
        end
        G(N,N) = -DDIC/dz^2 - s/2/dz - Rcal*Ca(N+1) - KMn*Rcal*Mn(N+1); % free flux at the bottom
        
        RHS = -SSA(2:end);
        RHS(1) = RHS(1) - DIC0 * (DDIC/dz^2 + s/2/dz);
        f=RHS;
        n=length(f);
        v=zeros(n,1);
        y=v;
        w=G(1,1);
        y(1)=f(1)/w;
        for j=2:n
            v(j-1)=G(j-1,j)/w;
            w=G(j,j)-G(j,j-1)*v(j-1);
            y(j)=(f(j)-G(j,j-1)*y(j-1))/w;
        end
        for j=n-1:-1:1
            y(j)=y(j)-v(j)*y(j+1);
        end
        DIC_prev = DIC;
        DIC = [DIC0; y];
        
        % Resolve [Ca2+] by the analytical expression
        Ca_prev = Ca;
        vv = (s - sqrt(s^2 + 4*Rcal*DCa*DIC)) / (2*DCa);
        VI = 0 * zz;
        for i = 1:N+1
            VI(i) = ISimpson(vv, i, dz);
        end
        for i=1:N+1
            Ca(i) = Ca0 * exp(VI(i));
        end
        if mean(abs(Ca_prev - Ca)) < eps && mean(abs(DIC_prev - DIC)) < eps
            calciteFinal = 1;
            break;
        end
        if ii == nii
            'Calcite!'
        end
    end

    SSA = RR .* Mns .* CH;
    Mn_prev = Mn;
    % Thomas method
    G = -2 * DMn/dz^2 * eye(N);
    for i=1:N-1
        G(i,i+1) = DMn/dz^2 - s/2/dz;
        G(i+1,i) = DMn/dz^2 + s/2/dz;
        G(i,i) = G(i,i) - KMn*Rcal*DIC(i+1) - Rox*O2(i+1);
    end
    G(N,N) = -DMn/dz^2 - s/2/dz - KMn*Rcal*DIC(i+1) - Rox*O2(N+1); % free flux at the bottom
    
    RHS = -SSA(2:end);
    f=RHS;
    n=length(f);
    v=zeros(n,1);
    y=v;
    w=G(1,1);
    y(1)=f(1)/w;
    for j=2:n
        v(j-1)=G(j-1,j)/w;
        w=G(j,j)-G(j,j-1)*v(j-1);
        y(j)=(f(j)-G(j,j-1)*y(j-1))/w;
    end
    for j=n-1:-1:1
        y(j)=y(j)-v(j)*y(j+1);
    end
    Mn = [0; y];

    if mean(abs(Mn_prev - Mn)) < eps && calciteFinal
        break;
    end

    if jj == nii
        'Fail!'
    end
end

calcite = Rcal * Ca .* DIC;
Mn_cal = KMn * Rcal * Mn .* DIC;


pN = 7;

mk = 1;
subplot(1, pN, mk);
plot(O2 * 1e3, zz, 'lineWidth', 1.5); hold on;
set(gca, 'yDir', 'reverse', 'lineWidth', 1, 'fontSize', 10);
xlabel('[O_2] (μM)'); ylabel('depth (cm)');

mk = mk + 1;
subplot(1, pN, mk);
plot(Mns, zz, 'lineWidth', 1.5); hold on;
set(gca, 'yDir', 'reverse', 'lineWidth', 1, 'fontSize', 10);
xlabel('[Mn(SR)] (mM)');

mk = mk + 1;
subplot(1, pN, mk);
plot(CH, zz, 'lineWidth', 1.5); hold on;
set(gca, 'yDir', 'reverse', 'lineWidth', 1, 'fontSize', 10);
xlabel('[CH_2O] (mM)');

mk = mk + 1;
subplot(1, pN, mk);
plot(DIC * 1e3, zz, 'lineWidth', 1.5); hold on;
set(gca, 'yDir', 'reverse', 'lineWidth', 1, 'fontSize', 10);
xlabel('[DIC] (μM)');

mk = mk + 1;
subplot(1, pN, mk);
plot(Ca, zz, 'lineWidth', 1.5); hold on;
set(gca, 'yDir', 'reverse', 'lineWidth', 1, 'fontSize', 10);
xlabel('[Ca^{2+}] (mM)');

mk = mk + 1;
subplot(1, pN, mk);
plot(log10(Ca .* DIC / Ksp), zz, 'lineWidth', 1.5); hold on;
% plot([log10(0.6), log10(0.6)], [0, maxD]); hold on;
set(gca, 'yDir', 'reverse', 'lineWidth', 1, 'fontSize', 10);
xlabel('log(\Omega)');

mk = mk + 1;
subplot(1, pN, mk);
plot(Mn * 1e3, zz, 'lineWidth', 1.5); hold on;
set(gca, 'yDir', 'reverse', 'lineWidth', 1, 'fontSize', 10);
xlabel('[Mn^{2+}] (μM)');


Mncarb_ppm = ISimpson(Mn_cal, N+1, dz)/ISimpson(calcite, N+1, dz) *55/100 * 1000000


function ret = ISimpson(f, n, h) % Simpson
    ret = (f(1) + f(n) + 2*sum(f(2:2:n-1)) + 4*sum(f(3:2:n-1))) * h / 3;
end