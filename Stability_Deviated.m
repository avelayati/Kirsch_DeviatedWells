%%
%%% Arian Velayati, PhD
%%%% This script is used to calculate stresses at the wall of deviated
%%%% wellbores and identifying stress magnitudes and locations for shear
%%%% failure (breakouts)and tensile fracs

clc; clear; close;

%% Inputs

Azi = 70; % Wellbore azimuth
Inc = 90; % Wellbore inclination
Sp = [70 0 0;0 67 0;0 0 45];% Principal in-situ stress tensor
Shmin_azi = 70; % Azimuth of min horizontal stress
SHmax_azi = 160;
Pp = 32; % MPa (Pore Pressure)
Pw = 32; %MPa; Mud pressure

% Inputs for kirsch
v = 0.3; % Poisson ratio
u = 0.8; % Friction coefficient

%% Transformation matrix : links the geographical coordinate system and the wellbore coordinate system


Rgw = [-cosd(Azi)*cosd(Inc) -sind(Azi)*cosd(Inc) sind(Inc)
        sind(Azi) -cosd(Azi) 0;
        cosd(Azi)*sind(Inc) sind(Azi)*sind(Inc) cosd(Inc)]; % Transformation matrix

    
% Director angles
I = input('Select the faulting regime: NF = 1; SS = 2; RF = 3:  ');

if I == 1
% Normal Faulting
a = Shmin_azi;
B = 90;
Gamma = 0;
    elseif I == 2
a = SHmax_azi;
B = 0;
Gamma = 90;
    else
a = SHmax_azi;
B = 0;
Gamma = 0; 
end

% RPG corresponding change (as a result of principal in-situ stresses Sp)
% of coordinate matrix to the geographical coordinate system 
RPG = [cosd(a)*cosd(B), sind(a)*cosd(B), -sind(B);
       cosd(a)*sind(B)*sind(Gamma)-sind(a)*cosd(Gamma), sind(a)*sind(B)*sind(Gamma)+cosd(a)*cosd(Gamma), cosd(B)*sind(Gamma);
       cosd(a)*sind(B)*cosd(Gamma) + sind(a)*sind(Gamma), sind(a)*sind(B)*cosd(Gamma)-cosd(a)*sind(Gamma), cosd(B)*cosd(Gamma)];

SG = RPG' * Sp * RPG;  % corresponding geographical location stress tensor

%% Method 1: Finding wellbore coordinate system stresses from geographical stress tensot

SW = Rgw*SG*Rgw';

%% Method 2 : Finding wellbore coordinate system stresses from the principal stress tensor Sp

SW2 = Rgw*RPG'*Sp*RPG*Rgw';

%% Kirsch eqn for isotropic rock with far-field shear stresses (wellbore wall)

% Effective stresses
sig11 = SW(1,1)-Pp; sig22 = SW(2,2)-Pp;sig33 = SW(3,3)-Pp; sig12 = SW(1,2); sig13 = SW(1,3); sig23 =  SW(2,3); 

% s12 s13 and s23 are needed in order to account for principal stresses not coinciding with the wellbore orientation

% Stresses at the wellbore wall 

teta = 0;
Delp = Pw - Pp;
sigrr  = Delp;
sigtt  = sig11 + sig22 -2*(sig11-sig22)*cosd(2*teta)-4*sig12*sind(2*teta)-Delp;
tau_tz = 2*(sig23*cosd(2*teta)-sig13*sind(teta));
sigzz = sig33 - 2*v*(sig11-sig22)*cosd(2*teta) - 4*v*sig12*sind(2*teta);
S_k0 = [sigrr sigtt tau_tz sigzz teta];

teta = 1:360;
S_kout = S_k0;

for i = 1:length(teta)
    
Delp = Pw - Pp;
sigrr  = Delp;
sigtt  = sig11 + sig22 -2*(sig11-sig22)*cosd(2*teta(i))-4*sig12*sind(2*teta(i))-Delp;
tau_tz = 2*(sig23*cosd(2*teta(i))-sig13*sind(teta(i)));
sigzz = sig33 - 2*v*(sig11-sig22)*cosd(2*teta(i)) - 4*v*sig12*sind(2*teta(i));

S_k = [sigrr sigtt tau_tz sigzz teta(i)];

S_kout = [S_kout;S_k];

end

%% Required UCS to avert shear failure

q = (sqrt(u^2 + 1) + u)^2 ; % Anisotropy factor in Mohr Cuolomb failure criterion

UCS = zeros(361,1);

for i = 1:361
    UCS(i) = max(S_kout(i,1:4))-q*min(S_kout(i,1:4));
end

T = find(UCS == max(UCS)); disp(' angle (teta) with the highest required ucs is : (Potential breakouts here) '); disp(T(1)-1)
disp(' The required UCS for this angle is : '); disp(max(UCS))