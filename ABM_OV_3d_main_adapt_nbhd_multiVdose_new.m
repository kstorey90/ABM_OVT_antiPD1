%ABM_OV_main.m
%Simulate ABM modeling GBM treatment with OVT and antiPD-1
% changed tau to 1 (from 0.5)
% include a second viral dose at a specified time
% stores virus in state matrix

%allows to start center dosing, and switch to adaptive


%% Create and initialize the states
%clear variables; close all; clc;
CELLSIZE       = 0.0014;%0.0018;    %[cm]
days           = 80;%40;     %should be integer
totTime        = 24*days;     %[hours]
totRep         = 1;%10;
n              = 100;    %grid size (change later)
L              = n*CELLSIZE;            %domain length [cm]
%total          = zeros(n,n,n,4,totTime/24+1,totRep);

rng('shuffle');

%%%%% Change if you want to save/not save data throughout!
saveOn = 0;  % 1 to save state data throughout sim, 0 to save only final allStates

%%%Set OV on(1) or off(0):
OV=1;   %I don't think I'm using this anymore (just change dose size to 0)

%%antiPD1 on(1) or off(0):
antiPD1 = 1;%1;
APdose = 3; %anti-PD-1 dose in mg/m^2

%Set PD-1/PD-L1 checkpoint on/off:
PD1L1 = 1;


a_AT = 0.06;%0.07;%0.1;%0.0016;      %Rate of tumor cell-mediated proliferation of antiT adap imm cells (baseline: 0.025)

k_TA = 0.02;%1/24;%0.1;%1/24;        %killing rate of tumor cells by antiT adap (now using k=0.03)

g_inn1 = 5e-5;%1e-5;%5e-6;     %par used in P_inn - viral-med activ of innate (using baseline: 5e-5)

%CHECK ON THIS:
beta = 2.5e-8;%1e-7;%2.5e-8;              %viral infection rate

delta_A = 0.0019;%0.008;%0.0019;%0.01;%0.05;               % decay rate of anti-PD-1 (had been using 0.008)

numCent=1;%6;%1;%5;%3;  %the number of first center doses to give, before switching to dosingType
dosingType='adap';%'adap';%'center';%'168';
adapType ='max'; %'max2';%'max';%'max_split3';

doseTimes = [0 10 13 30 33 60];%[0 10 13 30 33];%[0 10 13 30 50];%[0 10 13 55];%[0 10 13 55 85];%[0 10 13 55];%[0 4];%[0 10];      %times of all doses
numDoses  = size(doseTimes,2);      %total number of OV doses
Vdose_split = [1 1 1 1 1 1];% 1];%[1 1 1 1 1];%[1 3 3];%[1 6 8];%[1 1 1];%[1 4 2];%1   %number of parts to split each viral dose into (number of entries = numDoses)
%Vdose_split = [6 9];
           %radius (# of sites) of viral dose)
   
 % for center dosing:          
doseLocs =  zeros(numDoses,3);          

%Don't need this anymore:
doseLocs1 = [0 0 0];        %shifts from the center for the dose locations
doseLocs2 = [0 0 0];%[10 0 0; -10 0 0; 0 10 0; 0 -10 0; 0 0 10; 0 0 -10];%[10 0 0;];%%   
%[10 0 0; -10 0 0; 0 10 0; 0 -10 0];
doseLocs3 = [0 0 0];%[7 7 7; -7 -7 7; 7 -7 7; -7 7 7; 7 7 7; -7 -7 -7; 7 -7 -7; -7 7 -7];%[-10 0 0];%%[0 0 10; 0 0 -10];
            %[5 5 5; -5 -5 5; 5 -5 5; -5 5 5; 5 5 5; -5 -5 -5; 5 -5 -5; -5 5 -5];% % If splitting dose1:
% doseLocs1 = [30 0 0; -30 0 0; 0 30 0; 0 -30 0; 0 0 30; 0 0 -30]; 
% doseLocs2 = [0 0 0; 20 20 20; -20 -20 20; 20 -20 20; -20 20 20; 20 20 20; -20 -20 -20; 20 -20 -20; -20 20 -20];
doseLocs4 = [0 0 0];%[10 0 0; -10 0 0; 0 10 0; 0 -10 0; 0 0 10; 0 0 -10];

dens_nbhd_sz = 5;   %size of the neighborhood to determine viral dosing (baseline=5)

doseSize = 10^6;%10^5;%10^7;        %viral dose (in tumor sphere)
r_v = 2; 

max_time=totTime/24+1;
%vol=zeros((max_time-1),10,totRep);



g_inn2 = 5e-5;%5e-3;%5*10^(-2);             %par used in P_inn - innate pos feedback activ

gamma_AV = 2;%10;       %par used in prob that new adap imm cell will be antiviral


g_div = 1e-5;%1e-6;%1*10^(-4);         %par used in probability of reducing cell div counter


tau = 1; %0.5;%1;      %time step for ABM (in hours)

m_Z = 2;%1.1;          %times more likely for innate imm cell to move toward an inf cell (should be >1)
m_AT = 2;%1.5;         %times more likely for antitumor adap imm cell to move toward a tumor cell (should be >1)
m_AV = 2;%1.5;         %times more likely for antiviral adap imm cell to move toward an inf cell (should be >1)

r_movZ = 0.05;%0.1;       %rate of movement of innate immune cell to an empty neighboring site
r_movAT = 0.01;%0.05;       %rate of movement of an adapT immune cell to an empty neighboring site
r_movAV = 0.01;%0.05;       %rate of movement of an adapT immune cell to an empty neighboring site


%K_YQ = 1.5649e-22;    % inhib of T cells by PD-1/PD-L1
K_YQ = 1e-15;%1.296e-9;

%mu_PA = delta_A/(1.022e-10);    % blocking rate of PD-1
mu_PA = delta_A/(2*2.36e-5);%delta_A/(5*2.36e-5);%delta_A/(9*2.36e-5);


%dir_name='figures/antiPD1_inf5e-8';
if PD1L1
    if antiPD1
        %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_ginn' num2str(g_inn1) '_gdiv' num2str(g_div)]; %'_oneVdose']; %'_burst2d'];
        %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_ginn' num2str(g_inn1) '_gdiv' num2str(g_div) '_tau' num2str(tau) '_antiPD1_t4']; %'_oneVdose']; %'_burst2d'];
        %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_aAT' num2str(a_AT)]; %'_oneVdose']; %'_burst2d'];
        if numDoses==1
            dir_name = ['figures_3d/antiPD1_nbhd_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_ginn1_' num2str(g_inn1)];  %_parts' num2str(Vdose_split)]; %'_oneVdose']; %'_burst2d'];
        else
            %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_dose_' dosingType '_tau'  num2str(tau) '_adapt'];
            
            %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_dose_' dosingType '_tau'  num2str(tau) '_aAT'  num2str(a_AT)];
            %adaptive:
            %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end))  '_adapt_density_aAT_' num2str(a_AT)];
            %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end))  '_adapt_density_kTA_' num2str(k_TA)];
            %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end))  '_adapt_density_ginn1_' num2str(g_inn1)];
            if strcmp(dosingType,'center')
                %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_dose_' dosingType '_aAT_'  num2str(a_AT)];
                % when varying the third dose:
                %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_Third_t' num2str(doseTimes(3)) '_dose_' dosingType '_aAT_'  num2str(a_AT)];
                %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_Third_t' num2str(doseTimes(3)) '_dose_' dosingType '_ginn1_'  num2str(g_inn1)];
                %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_Third_t' num2str(doseTimes(3)) '_dose_' dosingType '_beta_'  num2str(beta)];
                dir_name = ['figures_3d/upd_nbhd_antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_Third_t' num2str(doseTimes(3)) '_dose_' dosingType '_aAT_'  num2str(a_AT) '_nbhdSz_' num2str(dens_nbhd_sz)];
                %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_Third_t' num2str(doseTimes(3)) '_dose_' dosingType '_aAT_'  num2str(a_AT) '_kTA_' num2str(k_TA)];

            else
                if totRep > 1
                    dir_name = ['figures_3d/upd_nbhd8_antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_Third_t' num2str(doseTimes(3)) '_dose_cent_' num2str(numCent) '_' dosingType '_' num2str(adapType) '_aAT_'  num2str(a_AT) '_reps' num2str(totRep)];
                else
                %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_dose_' dosingType '_' num2str(adapType) '_aAT_'  num2str(a_AT)];
                % when varying the third dose:
                %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_Third_t' num2str(doseTimes(3)) '_dose_' dosingType '_' num2str(adapType) '_aAT_'  num2str(a_AT)];
                %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_Third_t' num2str(doseTimes(3)) '_dose_' dosingType '_' num2str(adapType) '_ginn1_'  num2str(g_inn1)];
                %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_Third_t' num2str(doseTimes(3)) '_dose_' dosingType '_' num2str(adapType) '_beta_'  num2str(beta)];
                dir_name = ['figures_3d/upd_nbhd_antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_Third_t' num2str(doseTimes(3)) '_dose_cent_' num2str(numCent) '_' dosingType '_' num2str(adapType) '_aAT_'  num2str(a_AT) '_mu_' num2str(mu_PA)]; %'_kTA_' num2str(k_TA)];
                %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_Third_t' num2str(doseTimes(3)) '_dose_' dosingType '_' num2str(adapType) '_aAT_'  num2str(a_AT) '_kTA_' num2str(k_TA)];
                end
            end
                %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_dose_' dosingType '_tau'  num2str(tau) '_adapt_ginn1_'  num2str(g_inn1)];
            %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_dose_' dosingType '_tau'  num2str(tau) '_adapt_aAT_'  num2str(a_AT)];
       
            %%%if varying aAT:
            %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_small_tumor_aAT'  num2str(a_AT)];
            %%%if varying ginn1:
            %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_dose_center_ginn1_' num2str(g_inn1) ];
            %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_imm_mig_ginn2_' num2str(g_inn2) ]; %'_parts' num2str(Vdose_split)]; 
            %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_1Vdose_parts6_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_parts' num2str(Vdose_split)]; 
            %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_dosesplit142_imm_mig_ginn1_' num2str(g_inn1)];
            %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_dose_center_kTA'  num2str(k_TA)];
            %dir_name = ['figures_3d/antiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_dose_center_aAT'  num2str(a_AT)];
        end
            %dir_name = ['figures_3d/testing'];
    
    else
        %dir_name = ['figures_3d/noantiPD1_inf' num2str(beta) '_t' num2str(days) '_ginn' num2str(g_inn1)]; %'_oneVdose']; %'_burst2d'];
        %dir_name = ['figures_3d/noantiPD1_inf' num2str(beta) '_t' num2str(days) '_gdiv' num2str(g_div) '_tau' num2str(tau)]; %'_oneVdose']; %'_burst2d'];
        %dir_name = ['figures_3d/noantiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_imm_mig']; %'_parts' num2str(Vdose_split)]; 
        %dir_name = ['figures_3d/noantiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_small_tumor_kTA'  num2str(k_TA)];
        if strcmp(dosingType,'center')
            dir_name = ['figures_3d/upd_noantiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_Third_t' num2str(doseTimes(3)) '_dose_' dosingType '_aAT_'  num2str(a_AT) ];%'_beta_'  num2str(beta)];
        else 
            dir_name = ['figures_3d/upd_noantiPD1_inf' num2str(beta) '_t' num2str(days) '_' num2str(numDoses) 'Vdose_t' num2str(doseTimes(end)) '_Third_t' num2str(doseTimes(3)) '_dose' dosingType '_' num2str(adapType) '_aAT_'  num2str(a_AT)]; %'_beta_'  num2str(beta)]; 
        end
    end
else
    dir_name = ['figures_3d/nocheckpt_inf' num2str(beta) '_t' num2str(days) '_ginn' num2str(g_inn1)]; %'_oneVdose']; %'_burst2d'];
end
mkdir(dir_name)
%filename=['/nocheckpt' OV];




%%Cell cycle times
r     = 6;%10;%5;%3;%28;%2;      %radius of initial susc. cells
Sm = 36.1;%18.3;        %mean tumor cell cycle time (baseline: 36.1--converted from growth rate)
Ssd = 2;%3;                %sd tumor cell cycle time (estimate)

%%Viral pars
%diffusion
%d         = 1;                     %viral diffusion rate
D         = 3.6e-4;                %viral diffusion coeff [cm^2/s]
delx      = CELLSIZE/L;            %RDE space step or 1/n
delt      = delx^2/(6*D);          %RDE time step
tol       = 1e-1;%1.0e-5;                %FD method convergence tolerance
courant   = (D*delt)/delx^2;       %Courant number (must be <= 0.5)
cInf      = 0;                %background viral concentration [mol/cm^3]
bc        = 0;%1;                     %viral boundary conditions
%other
omega = 0.025;%0.001;%0.025;              %viral clearance rate


%%other pars
b_T = 50;               %inf cell burst size (# of viral particles)
alpha_I = 1;            %viral part. uptaken by tumor cells during infec
delta_T = (1/18)*tau;   %death rate of infected cells 



%innate

k_VZ = 0.005;           % killing rate of virus by innate imm cells (/cell/hr)
k_I = 0.02;             % killing rate of infected cells by innate imm cells (/cell/hr)
a_Z = 2.4e-6;           %rate of infected cell-med prolif of innate imm cells (/cell/hr)
delta_Z = 0.008;        %death rate for innate imm cells (death time given by log(2)/delta_z)
sigma_Z = 5;            %stdev for time until death of innate imm cells
kappa_Z = 10;            %number of possible inf cell kills before death



%adaptive
%g_ad = 10^(-3);         %par used in P_ad - innate imm-med recruitment of adap
a_TZ = 0.05;        %rate of recruitment of adap by innate imm




delta_AT = 3.75e-4;        %death rate for antitumor adap imm cells (death time given by log(2)/delta_z)
delta_AV = 5.54e-3;        %death rate for antiviral adap imm cells (death time given by log(2)/delta_z)
sigma_AT = 10;            %stdev for time until death of antiT adap imm cells
sigma_AV = 10;            %stdev for time until death of antiV adap imm cells
kappa_AT = 10;%15;      %number of possible tumor cell kills before antitumor T cell death
kappa_AV = 10;%15;      %number of possible inf cell kills before antiviral T cell death

k_IA = 1/24;        %killing rate of inf cells by antiV adap
k_VA = 10^(-5);     %killing rate of virus by antiV adap

a_AV = 0.025;       %Rate of inf cell-mediated proliferation of antiV adap imm cells


%%PD-1/PD-L1
%rho_P = 8.9992e-15;  %molar conc of PD-1 per T cell
rho_P = 1.259e-11;
%rho_L = 1.793e-14;   %molar conc of PD-L1 per T cell
rho_L = 2.51e-11;   
eps_T = 10;         %exp of PD-L1 on tumor cells vs T cells
eps_Z = 10;         %exp of PD-L1 on innate imm cells vs T cells


%%anti-PD-1


%APconc = pi^(1/3)*(3/4*(0.139*APdose+0.064))^(2/3); %dose of A (in micromol/cm^2)
APconc = 0.139*APdose+0.064; %dose of A (in micromol/cm^3)

numTypes = 6; %added antiPD1

%%States:
% 'S'  = susceptible tumor cells (1 in 1st entry)
% 'I'  = infected tumor cells    (2 in 1st entry)
% 'II' = innate immune cells     (1 in 2nd entry)
% 'AT' = adaptive antitumor immune cells  (1 in 3rd entry)
% 'AV' = adaptive antiviral immune cells  (2 in 3rd entry)
% 'E'  = empty (0)

meanDZ = log(2)/delta_Z;    %mean death time for innate imm cells
meanDAT = log(2)/delta_AT;    %mean death time for antiT imm cells
meanDAV = log(2)/delta_AV;    %mean death time for antiV imm cells

% store parameter values (for reference later) 
% add more here if I change values
paramfile = [dir_name '/params.mat'];
save(paramfile,'beta','g_inn1','g_inn2','gamma_AV','g_div','a_AT','a_TZ', 'doseTimes', 'm_Z', 'm_AT', 'm_AV', 'r_movZ', 'r_movAT', 'r_movAV', 'doseSize', 'doseLocs1','doseLocs2','doseLocs3','r', 'delta_A', 'k_TA','tau', 'K_YQ', 'mu_PA')    

finalPops = zeros(totRep,numTypes+1);   %save the final pops for all sims
allSimSusc = zeros(totTime/tau+1,totRep); %save the susc pop for all t, for all sims

tic %timer

for rep=1:totRep
    
    totTypes = zeros(1,numTypes);   %Stores a vector of total number of each cell type
                                    %1=Susc, 2=Inf, 3=Innate, 4=AdapT, 
                                    %5=AdapV, 6=antiPD1
                                    
    
    
    %virus    = zeros(n,n,n);               %initial virus concentration 
    
    
    %% Initialize grid with cells
    
    h     = n/2;    %center of the circle
    state = zeros(n,n,n,4);  %store cells at each site %In third dim:
                                 %first entry: tumor cells, 
                                 %second entry: innate imm cells,
                                 %third entry: adap immune cells
                                 %fourth entry: virus
    for x=1:n
        for y=1:n
            for z=1:n
                if (x-h)^2+(y-h)^2+(z-h)^2<=r^2
                    state(x,y,z,1)=1; %add susc cell
                    totTypes(1) = totTypes(1)+1;
                end
            end
        end
    end
    
    % Randomly initialize cell cycle times (0,Sm)
    divisionCounter = zeros(n,n,n);
    divisionCounter(state(:,:,:,1)==1) = abs(floor(Sm.*rand(size(divisionCounter(state(:,:,:,1)==1)))));
    
    
    % Innate immune cell death counter
    innDeathCounter = zeros(n,n,n,2);     % store time ctr until natural death in (:,:,1) and
                                        % inf cell kill ctr until death in (:,:,2)
                                        
    % Adaptive antitumor immune cell death counter
    adTDeathCounter = zeros(n,n,n,2);     % store time ctr until natural death in (:,:,1) and
                                        % inf cell kill ctr until death in (:,:,2)
     % Adaptive antiviral immune cell death counter
    adVDeathCounter = zeros(n,n,n,2);
    

    AllTypes = zeros(totTime/tau, 2+numTypes);    % store totTypes+virus for all time steps (time in col 1)                            
    curvirus=sum(sum(sum(state(:,:,:,4))));
    AllTypes(1,:) = [0 totTypes curvirus];            
    
    %% Store states at each time point
    t=0;
%     allStates = zeros(n,n,n,4,totTime/24+1);
%     allStates(:,:,:,1:3,t+1) = state;
%     allStates(:,:,:,4,t+1) = virus;
    t=t+1;
    
    Vdose=0; %keep track of number of viral doses given
    
    %I CHANGED THIS TO START ANTI-PD-1 EARLIER:
    lastAPdose = -10*24;%-7*24;
    
    %% Update grid
    for iTime=tau:tau:totTime
        iTime
        
        if Vdose>0
            %infected cell lysis
            tmp = ones(n,n,n);
            tmp(state(:,:,:,1)==2) = rand(size(tmp(state(:,:,:,1)==2)));
            lysis_mat = zeros(n,n,n);
            lysis_mat(tmp<delta_T) = rand(size(lysis_mat(tmp<delta_T)));
            %Iterate over cells being removed 
            while nnz(lysis_mat)~=0
                [i,j,k] = ind2sub([n n n],find(lysis_mat==max(max(max(lysis_mat))))); %indices of cell with highest rand
                lysis_mat(i,j,k)=0;
                state(i,j,k,1)=0;
                totTypes(2) = totTypes(2)-1;
                state(i,j,k,4) = state(i,j,k,4) + b_T;
            end
            
            
            %cell infection by virus
            %p_inf = totTypes(1)*beta*tau*virus;
            p_inf = totTypes(1)*beta*tau*state(:,:,:,4);
            infec_mat=ones(n,n,n);
            infec_mat(state(:,:,:,1)==1) = rand(size(infec_mat(state(:,:,:,1)==1)));
            for i=1:n
                for j=1:n
                    for k=1:n
                        if (infec_mat(i,j,k)<p_inf(i,j,k) && p_inf(i,j,k)<1)
                            state(i,j,k,1)=2;
                            totTypes(1) = totTypes(1)-1;
                            totTypes(2) = totTypes(2)+1;
                            %virus(i,j,k) = max(virus(i,j,k) - alpha_I,0);  %virus uptaken during infec
                            state(i,j,k,4) = max(state(i,j,k,4) - alpha_I,0);
                            %iTime
                        end
                    end
                end 
            end

            
            %Include viral diffusion step here
            %totalT=0;
            %virus = virus - omega*tau*virus;   %natural clearance
            state(:,:,:,4) = state(:,:,:,4) - omega*tau*state(:,:,:,4);
            tmp = zeros(n,n,n);
            for tt=1:tau/delt
                %totalT=totalT+delt;
                tmp(1,:,:) = bc;
                tmp(:,1,:) = bc;
                tmp(:,:,1) = bc;
                tmp(n,:,:) = bc;
                tmp(:,n,:) = bc;
                tmp(:,:,n) = bc;
                for i = 2:n-1
                    for j = 2:n-1
                        for k = 2:n-1
                            %tmp(i,j,k) = virus(i,j,k)+courant*(virus(i-1,j,k)+virus(i+1,j,k)+virus(i,j-1,k)+virus(i,j+1,k)+virus(i,j,k-1)+virus(i,j,k+1)-6*virus(i,j,k));%-delt*consumptionRate(i,j);
                            tmp(i,j,k) = state(i,j,k,4)+courant*(state(i-1,j,k,4)+state(i+1,j,k,4)+state(i,j-1,k,4)+state(i,j+1,k,4)+state(i,j,k-1,4)+state(i,j,k+1,4)-6*state(i,j,k,4));
                        end  
                    end
                end
               
                state(:,:,:,4)=tmp;
            end
        end
        
        %natural depletion of antiPD1
        totTypes(6) = totTypes(6)*(1-delta_A*tau);
        
        %% Viral dose 
        if (Vdose<numDoses && iTime/24>=doseTimes(Vdose+1)) 
            numSitesDose = (2*r_v-1)^3+6;       %number of sites dose is dist'd over
                                                %r_v=radius of initial susc. cells
            if Vdose< numCent
                doseSplit=Vdose_split(1);
                h     = n/2*ones(Vdose_split(1),3) + doseLocs1;
               %h     = [n/2 n/2 n/2];    %center of the circle
               %doseSplit=1;
            else
                
                curDose=Vdose+1;
                if strcmp(dosingType,'adap')
                    
                    
                    
                    statemat1=zeros(n,n,n);     %new matrix to contain 1 for all tumor cells
                    statemat1(state(:,:,:,1)~=0)=1;
                    
                    nbhd_size = ones(dens_nbhd_sz,dens_nbhd_sz,dens_nbhd_sz);
                    nbhd_dens = convn(statemat1,nbhd_size,'same');
                    maxDens = max(nbhd_dens(:));
                    
                   
                    [xMaxes, yMaxes, zMaxes] = ind2sub(size(nbhd_dens),find(nbhd_dens == maxDens));   %find row,col index for max dens

                    if strcmp(adapType,'max')
                        doseSplit=size(xMaxes,1);     %number of sites with maxdens
                    elseif strcmp(adapType,'max2')
                        nbhd_dens(nbhd_dens==maxDens)=0;
                        maxDens = max(nbhd_dens(:));
                        [newx, newy, newz] = ind2sub(size(nbhd_dens),find(nbhd_dens == maxDens));
                        xMaxes = [xMaxes; newx];
                        yMaxes = [yMaxes; newy];
                        zMaxes = [zMaxes; newz];
                        doseSplit=size(xMaxes,1);     %number of sites with maxdens
                    else %fix top 3
                        doseSplit=Vdose_split(curDose);
                        while length(xMaxes) < doseSplit %get the next highest if we don't have enough locs
                            nbhd_dens(nbhd_dens==maxDens)=0;
                            maxDens = max(nbhd_dens(:));
                            [newx, newy, newz] = ind2sub(size(nbhd_dens),find(nbhd_dens == maxDens));
                            xMaxes = [xMaxes; newx];
                            yMaxes = [yMaxes; newy];
                            zMaxes = [zMaxes; newz];
                        end

                    end 


                        %doseSplit=size(rowsOfMaxes,1);     %number of sites with maxdens
                        h     = zeros(doseSplit,3);         %h matrix to store sites for doses
                        for dosenum=1:doseSplit
                           h(dosenum,:) = [xMaxes(dosenum), yMaxes(dosenum), zMaxes(dosenum)];
                        end

                        doseLocsFile = [dir_name '/doseLocs' num2str(curDose) '_sim' num2str(rep) '.mat'];
                        save(doseLocsFile, 'h')
                else     
                    doseSplit=Vdose_split(curDose);
                    
                    curDoseLocs = doseLocs(curDose,:);
                    h     = (n/2)*ones(Vdose_split(curDose),3) + curDoseLocs;
                end
                
%                 if ~isempty(rowT)
%                     maxDim=abs(rowT(floor(size(rowT,1)/2))-n/2);
%                     doseLocs2 = [maxDim 0 0; -maxDim 0 0; 0 maxDim 0; 0 -maxDim 0; 0 0 maxDim; 0 0 -maxDim];
%                 else
%                     doseLocs2 = [10 0 0; -10 0 0; 0 10 0; 0 -10 0; 0 0 10; 0 0 -10];
%                 end
%                     
%                 h     = n/2*ones(Vdose_split(2),3) + doseLocs2; 
                    
            end

               
            pfuPerSite = doseSize/(doseSplit*numSitesDose);      %initial pfu per lattice site 
            
            
            for x=1:n
                for y=1:n
                    for z=1:n
                        if any((x-h(:,1)).^2+(y-h(:,2)).^2+(z-h(:,3)).^2<=r_v^2)
                            %virus(x,y,z)=virus(x,y,z)+pfuPerSite; %add virus to site
                            state(x,y,z,4)=state(x,y,z,4)+pfuPerSite;
                        end
                    end
                end
            end

            Vdose = Vdose + 1;
        end
        
        %%antiPD1 dose
        if antiPD1
            if iTime >= lastAPdose + 14*24
                totTypes(6) = totTypes(6) + APconc;
                lastAPdose = iTime;
            end
        end
        
        
        %Update cell cycle time
        %totalFree = checkTumorNeighbors(state,n);
        
        % Reduce division counter
        p_div = exp(-g_div*(totTypes(1)+totTypes(2)));
        tmp = ones(n,n,n);
        tmp(state(:,:,:,1)==1) = rand(size(tmp(state(:,:,:,1)==1)));
        divCtrFlag=zeros(n,n,n);
        divCtrFlag(tmp<p_div) = rand(size(divCtrFlag(tmp<p_div)));
        while nnz(divCtrFlag)~=0
            [i,j,k] = ind2sub([n n n],find(divCtrFlag==max(max(max(divCtrFlag)))));
            divCtrFlag(i,j,k)=0;
            divisionCounter(i,j,k) = divisionCounter(i,j,k)-tau;
        end
        
        
        
        
        
        %%% Innate immune events
        if totTypes(3) > 0
            
            %Reduce death counter
            innDeathCounter(state(:,:,:,2)==1) = innDeathCounter(state(:,:,:,2)==1) - tau;
%             for i=1:n
%                 for j=1:n
%                     for k=1:n 
%                         if state(i,j,k,2)==1
%                             innDeathCounter(i,j,k,1) = innDeathCounter(i,j,k,1) - tau;
%                         end
%                     end
%                 end
%             end
            
            % innate immune killing of virus 
            for i=1:n
                for j=1:n
                    for k=1:n 
                        if (state(i,j,k,2)==1) 
                            state(i,j,k,4)=state(i,j,k,4) - floor(k_VZ*tau*state(i,j,k,4));
                            %red. clock?
                        end
                    end
                end 
            end
            
            % innate imm killing of inf cells
            p_IZ = k_I*tau;    %probability of innate imm cell killing inf cell
            tmp = ones(n,n,n);
            tmp(state(:,:,:,2)==1 & state(:,:,:,1)==2) = rand(size(tmp(state(:,:,:,2)==1 & state(:,:,:,1)==2)));    
            virInnFlag=zeros(n,n,n);
            virInnFlag(tmp<p_IZ) = rand(size(virInnFlag(tmp<p_IZ)));
            while nnz(virInnFlag)~=0
                [i,j,k] = ind2sub([n n n],find(virInnFlag==max(max(max(virInnFlag)))));
                virInnFlag(i,j,k)=0;
                state(i,j,k,1)=0;
                totTypes(2) = totTypes(2)-1;
                innDeathCounter(i,j,k,2) = innDeathCounter(i,j,k,2) - 1;
            end
            
            
            %innate imm death
            for i=1:n
                for j=1:n
                    for k=1:n 
                        if (state(i,j,k,2)==1 && (innDeathCounter(i,j,k,1)<=0 || innDeathCounter(i,j,k,2)<=0))
                            state(i,j,k,2)=0;
                            totTypes(3) = totTypes(3) - 1;
                            innDeathCounter(i,j,k,:)=[0; 0];
                        end
                    end
                end
            end
            
            %innate imm movement
            nbrTotals = checkAllNeighbors(state,n);    %check all neighbor states
            innNbrs = nbrTotals(:,:,:,3); 
            tmp = ones(n,n,n);
            tmp(state(:,:,:,2)==1 & innNbrs<26) = rand(size(tmp(state(:,:,:,2)==1 & innNbrs<26)));
            innMovFlag = zeros(n,n,n);
            innMovFlag(tmp<r_movZ*tau) = rand(size(innMovFlag(tmp<r_movZ*tau)));
            while nnz(innMovFlag)~=0
                [i,j,k] = ind2sub([n n n],find(innMovFlag==max(max(max(innMovFlag)))));
                innMovFlag(i,j,k)=0;
                [state,innDeathCounter] = immuneMigrate(state, i, j, k, 1, m_Z, innDeathCounter);    
            end
%             for i=1:n
%                 for j=1:n
%                     for k=1:n 
%                         if innMovFlag(i,j,k)< r_movZ*tau
%                             %innate immune cell moves
%                             [state,innDeathCounter] = immuneMigrate(state, i, j, k, 1, m_Z, innDeathCounter);   
%                         end
%                     end
%                 end
%             end
        
            
            %innate imm prolif
            p_innProl = a_Z*tau*totTypes(2);    %Prob of innate imm prolif at site occupied by both innate imm, inf tumor
            tmp = ones(n,n,n);
            tmp(state(:,:,:,1)==2 & state(:,:,:,2)==1) = rand(size(tmp(state(:,:,:,1)==2 & state(:,:,:,2)==1)));
            innProlifFlag=zeros(n,n,n);
            innProlifFlag(tmp<p_innProl) = rand(size(innProlifFlag(tmp<p_innProl)));
            while nnz(innProlifFlag)~=0
                [i,j,k] = ind2sub([n n n],find(innProlifFlag==max(max(max(innProlifFlag)))));
                innProlifFlag(i,j,k)=0;
                [state,totTypes,innDeathCounter]=immuneDivide(state, i, j, k, 1, totTypes, innDeathCounter, meanDZ, sigma_Z, kappa_Z);   
            end
            
            
        end
        
        %%%PD-1/PD-L1
        if PD1L1
            if totTypes(6)>0
                PD1 = rho_P*(totTypes(4)+totTypes(5));
                APblock = mu_PA*tau*PD1*totTypes(6);
                PD1 = max(PD1-APblock,0); %blocking by antiPD1  
                totTypes(6) = max(totTypes(6)-APblock,0);   %depletion of antiPD1
            else
                PD1 = rho_P*(totTypes(4)+totTypes(5));
            end
            PDL1 = rho_L*(totTypes(4)+totTypes(5)+eps_T*(totTypes(1)+totTypes(2))+eps_Z*totTypes(3));
            F_PL = 1/(1+PD1*PDL1/K_YQ);
        else
            F_PL = 1;
        end
        
        %%% Antitumor adaptive immune events
        if totTypes(4) > 0
            
            
            %Reduce death counter
            adTDeathCounter(state(:,:,:,3)==1) = adTDeathCounter(state(:,:,:,3)==1) - tau;
            
            
            % adapT imm killing of tumor cells
            p_TA = k_TA*tau;    %probability of antiT adap imm cell killing tumor cell
            tmp = ones(n,n,n);
            tmp(state(:,:,:,3)==1 & state(:,:,:,1)~=0) = rand(size(tmp(state(:,:,:,3)==1 & state(:,:,:,1)~=0)));
            tumAdapFlag=zeros(n,n,n);
            tumAdapFlag(tmp<p_TA) = rand(size(tumAdapFlag(tmp<p_TA)));
            while nnz(tumAdapFlag)~=0
                [i,j,k] = ind2sub([n n n],find(tumAdapFlag==max(max(max(tumAdapFlag)))));
                tumAdapFlag(i,j,k)=0;
                typeij = state(i,j,k,1);
                state(i,j,k,1)=0;
                totTypes(typeij) = totTypes(typeij)-1;
                adTDeathCounter(i,j,k,2) = adTDeathCounter(i,j,k,2) - 1;
            end
            
            
            
            % antiT adap imm death
            aTdeathFlag=zeros(n,n,n);
            aTdeathFlag(state(:,:,:,3)==1 & (adTDeathCounter(:,:,:,1)<=0 | adTDeathCounter(:,:,:,2)<=0)) = rand(size(aTdeathFlag(state(:,:,:,3)==1 & (adTDeathCounter(:,:,:,1)<=0 | adTDeathCounter(:,:,:,2)<=0))));
            while nnz(aTdeathFlag)~=0
                [i,j,k] = ind2sub([n n n],find(aTdeathFlag==max(max(max(aTdeathFlag))))); %indices of cell with highest rand
                aTdeathFlag(i,j,k)=0;
                state(i,j,k,3)=0;
                totTypes(4) = totTypes(4) - 1;
                adTDeathCounter(i,j,k,:)=[0; 0];
            end

            
            %adap movement
            nbrTotals = checkAllNeighbors(state,n);    %check all neighbor states
            adNbrs = nbrTotals(:,:,:,4)+nbrTotals(:,:,:,5); 
            tmp = ones(n,n,n);
            tmp(state(:,:,:,3)==1 & innNbrs<26) = rand(size(tmp(state(:,:,:,3)==1 & innNbrs<26)));
            adTMovFlag = zeros(n,n,n);
            adTMovFlag(tmp<r_movAT*tau) = rand(size(adTMovFlag(tmp<r_movAT*tau)));
            while nnz(adTMovFlag)~=0
                [i,j,k] = ind2sub([n n n],find(adTMovFlag==max(max(max(adTMovFlag)))));
                adTMovFlag(i,j,k)=0;
                [state,adTDeathCounter] = immuneMigrate(state, i, j, k, 2, m_AT, adTDeathCounter);    
            end

            
            %adap proliferation
            p_adTProl = a_AT*tau*F_PL;    %Prob of antiT adap imm prolif at site occupied by both adap imm, tumor
            tmp = ones(n,n,n);
            tmp(state(:,:,:,1)~=0 & state(:,:,:,3)==1) = rand(size(tmp(state(:,:,:,1)~=0 & state(:,:,:,3)==1)));
            adTProlifFlag=zeros(n,n,n);
            adTProlifFlag(tmp<p_adTProl) = rand(size(adTProlifFlag(tmp<p_adTProl)));
            while nnz(adTProlifFlag)~=0
                [i,j,k] = ind2sub([n n n],find(adTProlifFlag==max(max(max(adTProlifFlag)))));
                adTProlifFlag(i,j,k)=0;
                [state,totTypes,adTDeathCounter]=immuneDivide(state, i, j, k, 2, totTypes, adTDeathCounter, meanDAT, sigma_AT, kappa_AT); 
            end
            
          
        end
        
        
        
        %%% Antiviral adaptive immune events
        if totTypes(5) > 0
            
            
            %Reduce death counter
            adVDeathCounter(state(:,:,:,3)==2) = adVDeathCounter(state(:,:,:,3)==2) - tau;
            
            % ADD adap immune killing of virus 
            for i=1:n
                for j=1:n
                    for k=1:n 
                        if (state(i,j,k,3)==2) 
                            state(i,j,k,4)=state(i,j,k,4) - floor(k_VA*tau*state(i,j,k,4));
                            %red. clock?
                        end
                    end
                end 
            end
 
            
            %%adapV killing of inf cells
            p_VA = k_VA*tau;    %probability of antiV adap imm cell killing inf cell
            tmp = ones(n,n,n);
            tmp(state(:,:,:,3)==2 & state(:,:,:,1)==2) = rand(size(tmp(state(:,:,:,3)==2 & state(:,:,:,1)==2))); 
            infAdapFlag=zeros(n,n,n);
            infAdapFlag(tmp<p_VA) = rand(size(infAdapFlag(tmp<p_VA)));
            while nnz(infAdapFlag)~=0
                [i,j,k] = ind2sub([n n n],find(infAdapFlag==max(max(max(infAdapFlag)))));
                infAdapFlag(i,j,k)=0;
                state(i,j,k,1)=0;
                totTypes(2) = totTypes(2)-1;
                adVDeathCounter(i,j,k,2) = adVDeathCounter(i,j,k,2) - 1;
            end
            
            
            
            % antiV adap imm death
           
            aVdeathFlag=zeros(n,n,n);
            aVdeathFlag(state(:,:,:,3)==2 & (adVDeathCounter(:,:,:,1)<=0 | adVDeathCounter(:,:,:,2)<=0)) = rand(size(aVdeathFlag(state(:,:,:,3)==2 & (adVDeathCounter(:,:,:,1)<=0 | adVDeathCounter(:,:,:,2)<=0))));
            while nnz(aVdeathFlag)~=0
                [i,j,k] = ind2sub([n n n],find(aVdeathFlag==max(max(max(aVdeathFlag))))); %indices of cell with highest rand
                aVdeathFlag(i,j,k)=0;
                state(i,j,k,3)=0;
                totTypes(5) = totTypes(5) - 1;
                adVDeathCounter(i,j,k,:)=[0; 0];
            end

            
            %adap movement
            nbrTotals = checkAllNeighbors(state,n);    %check all neighbor states
            adNbrs = nbrTotals(:,:,:,4)+nbrTotals(:,:,:,5); 
            tmp = ones(n,n,n);
            tmp(state(:,:,:,3)==2 & innNbrs<26) = rand(size(tmp(state(:,:,:,3)==2 & innNbrs<26)));
            adVMovFlag = zeros(n,n,n);
            adVMovFlag(tmp<r_movAV*tau) = rand(size(adVMovFlag(tmp<r_movAV*tau)));
            while nnz(adVMovFlag)~=0
                [i,j,k] = ind2sub([n n n],find(adVMovFlag==max(max(max(adVMovFlag)))));
                adVMovFlag(i,j,k)=0;
                [state,adVDeathCounter] = immuneMigrate(state, i, j, k, 3, m_AV, adVDeathCounter);    
            end
            
           
            %adap proliferation
            p_adVProl = a_AV*tau*F_PL;    %Prob of antiV adap imm prolif at site occupied by both adap imm, inf tumor
            tmp = ones(n,n,n);
            tmp(state(:,:,:,1)==2 & state(:,:,:,3)==2) = rand(size(tmp(state(:,:,:,1)==2 & state(:,:,:,3)==2)));
            adVProlifFlag=zeros(n,n,n);
            adVProlifFlag(tmp<p_adVProl) = rand(size(adVProlifFlag(tmp<p_adVProl)));
            while nnz(adVProlifFlag)~=0
                [i,j,k] = ind2sub([n n n],find(adVProlifFlag==max(max(max(adVProlifFlag)))));
                adVProlifFlag(i,j,k)=0;
                [state,totTypes,adVDeathCounter]=immuneDivide(state, i, j, k, 3, totTypes, adVDeathCounter, meanDAV, sigma_AV, kappa_AV); 
            end
            
           
            
            
        end
        
        %UNCOMMENT BELOW
        
        %Add innate immune cells
        nbrTotals = checkAllNeighbors(state,n);    %check all neighbor states
        innNbrs = nbrTotals(:,:,:,3);                 %number of innate imm neighbors
        %p_inn = 1-exp(-g_inn1*virus-g_inn2*innNbrs);    %probability of new (activ)innate imm cell
        p_inn = 1-exp(-g_inn1*state(:,:,:,4)-g_inn2*innNbrs);    %probability of new (activ)innate imm cell
        inn_mat=ones(n,n,n);
        inn_mat(state(:,:,:,2)==0) = rand(size(inn_mat(state(:,:,:,2)==0)));
            for i=1:n
                for j=1:n
                    for k=1:n 
                        if (inn_mat(i,j,k)<p_inn(i,j,k))
                            state(i,j,k,2)=1;
                            totTypes(3) = totTypes(3)+1;
                            innDeathCounter(i,j,k,1) = initializeCycleDur(meanDZ,sigma_Z);
                            innDeathCounter(i,j,k,2) = kappa_Z;
                        end
                    end
                end 
            end
            
         
            
        %innate cells recruit adaptive imm cells
        %p_ad = 1-exp(-g_ad*totTypes(3));       %probability of new adap imm cell (scalar)
        if totTypes(1)+totTypes(2)>0
            adap_mat=ones(n,n,n);
            adap_mat(state(:,:,:,3)==0 & state(:,:,:,2)==1) = rand(size(adap_mat(state(:,:,:,3)==0 & state(:,:,:,2)==1)));
            p_ad=tau*a_TZ*F_PL;  %probability of recruitment
            for i=1:n
                for j=1:n
                    for k=1:n 
                        if (adap_mat(i,j,k)<p_ad)  %determine whether a recruitment occurs
                            %Prob that it is antitumor:
                            newATProp = exp(-gamma_AV*totTypes(2)/(totTypes(1)+totTypes(2)));
                            if rand < newATProp  
                                %new antitumor T cell
                                state(i,j,k,3)=1;
                                totTypes(4) = totTypes(4)+1;
                                adTDeathCounter(i,j,k,1) = initializeCycleDur(meanDAT,sigma_AT);
                                adTDeathCounter(i,j,k,2) = kappa_AT;
                            else
                                %new antiviral T cell
                                state(i,j,k,3)=2;
                                totTypes(5) = totTypes(5)+1;
                                adVDeathCounter(i,j,k,1) = initializeCycleDur(meanDAT,sigma_AT);
                                adVDeathCounter(i,j,k,2) = kappa_AV;
                            end     
                        end
                    end
                end
            end
        end
        
        
        
        
        %Cells divide
        % If division counter is zero then attempt division
        divisionFlag = zeros(n,n,n);
        divisionFlag(state(:,:,:,1)==1 & divisionCounter<=0) = rand(size(divisionFlag(state(:,:,:,1)==1 & divisionCounter<=0)));
       
        predivCounter=divisionCounter;
        % Iterate over dividing cells starting from the largest rand
        while nnz(divisionFlag)~=0
            [i,j,k] = ind2sub([n n n],find(divisionFlag==max(max(max(divisionFlag))))); %indices of cell with highest rand
            divisionFlag(i,j,k)=0;
            divisionCounter(i,j,k)=initializeCycleDur(Sm,Ssd);
            [state,divisionCounter,divisionFlag] = divide(state,divisionCounter,divisionFlag,i,j,k,n,Sm,Ssd);
            totTypes(1) = totTypes(1)+1;
        end
        
        
        
        %% Store results
        if mod(iTime,24)==0
%             allStates(:,:,:,1:3,t+1) = state;
%             allStates(:,:,:,4,t+1) = virus;
            if saveOn
                filename = [dir_name '/state_t' num2str(t)];
                save(filename,'state')
                filename = [dir_name '/AllTypes_t' num2str(t)];
                save(filename,'AllTypes')
                %counterfile = [dir_name '/counters_t' num2str(t)];
                %save(counterfile,'divisionCounter','innDeathCounter','adTDeathCounter','adVDeathCounter')
            
            end
            t=t+1;
            
        end
        
        %store totTypes
        curvirus=sum(sum(sum(state(:,:,:,4))));
        AllTypes(iTime*(1/tau)+1,:) = [iTime totTypes curvirus];
        

    end
     toc
     
%      total(:,:,:,:,:,rep) = allStates;
%      filename = [dir_name '/allStates'];
%      save(filename,'allStates')
     
     filename = [dir_name '/AllTypes_sim' num2str(rep)];
     save(filename,'AllTypes')
     
     allSimSusc(:,rep) = AllTypes(:,2);
     
     %state_end=zeros(n,n,n,4);
     %state_end(:,:,:,1:3)=state;
     %state_end(:,:,:,4)=virus;
     state_end=state;
     
     finalPops(rep,:) = [totTypes curvirus];


    if saveOn
        filename = [dir_name '/state_end_sim' num2str(rep)];
        save(filename,'state_end')
      
        %% 2-D cross sections (at z=n/2)
        state2d=state_end(:,:,n/2,:);
        % virus 
        virus=state2d(:,:,4);
        figure('name','Virus')
        imh_virus = imagesc(virus);
        cbh_virus = colorbar;
        filename = [dir_name '/cross2d_xy_virus'];
        savefig(filename)
        filename = [filename '.jpg'];
        saveas(gcf,filename)

        % tumor cells
        tumor=state2d(:,:,1,1);
        plotState=zeros(n,n);
        plotState(tumor==1)=1; %Susceptible
        plotState(tumor==2)=2; %Infected
        plotState(tumor==0)=3; %empty
        map = [0,0,1; 1,0,0 ; 1,1,1];
        figure('name','Tumor Cells')
        imh_state = imagesc(plotState);
        colormap(map)
        cbh_state = colorbar('YTick',[1.45 2.28 3],'YTickLabel',{'Susc','Inf','Empty'});
        xlabel('x')
        ylabel('y')
        filename = [dir_name '/cross2d_xy_tumor_cells'];
        savefig(filename)
        filename = [filename '.jpg'];
        saveas(gcf,filename)

        % innate imm cells
        innate=state2d(:,:,1,2);
        plotState=zeros(n,n);
        plotState(innate==1)=1; %Innate
        plotState(innate==0)=2; %empty
        map = [1,0,1 ; 1,1,1];
        figure('name','Innate immune cells')
        imh_state = imagesc(plotState);
        colormap(map)
        cbh_state = colorbar('YTick',[1.45 3],'YTickLabel',{'Innate','Empty'});
        xlabel('x')
        ylabel('y')
        filename = [dir_name '/cross2d_xy_innate'];
        savefig(filename)
        filename = [filename '.jpg'];
        saveas(gcf,filename)

        % adaptive imm cells
        adaptive=state2d(:,:,1,3);
        plotState=zeros(n,n);
        plotState(adaptive==1)=1; %adapT
        plotState(adaptive==2)=2; %adapV
        plotState(adaptive==0)=3; %empty
        map = [0,1,0; 0.66,0,0 ; 1,1,1];
        figure('name','Adaptive immune cells')
        imh_state = imagesc(plotState);
        colormap(map)
        cbh_state = colorbar('YTick',[1.45 2.28 3],'YTickLabel',{'antitumor adap','antiviral adap','Empty'});
        filename = [dir_name '/cross2d_xy_adaptive'];
        savefig(filename)
        filename = [filename '.jpg'];
        saveas(gcf,filename)




        %% 3D plots
        tumormat=state(:,:,:,1);
        suscind=find(tumormat==1);
        infind=find(tumormat==2);
        [S1, S2, S3]=ind2sub([n n n],suscind);
        [I1, I2, I3]=ind2sub([n n n],infind);
        figure('name','Tumor cells')
        scatter3(S1,S2,S3,'MarkerFaceColor',[0 0 1])
        hold on
        scatter3(I1,I2,I3,'MarkerFaceColor',[1 0 0])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        legend('Susc.','Inf.')
        %hold off
        filename = [dir_name '/plot3d_tumor'];
        savefig(filename)
        filename = [filename '.jpg'];
        saveas(gcf,filename)

        %innate imm
        innatemat=state(:,:,:,2);
        innind=find(innatemat==1);
        [In1, In2, In3]=ind2sub([n n n],innind);
        figure('name','Innate immune cells')
        scatter3(In1,In2,In3,'MarkerFaceColor',[1 0 1])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        filename = [dir_name '/plot3d_innate'];
        savefig(filename)
        filename = [filename '.jpg'];
        saveas(gcf,filename)

        %adap imm
        adaptivemat=state(:,:,:,3);
        ATind=find(adaptivemat==1);
        AVind=find(adaptivemat==2);
        [AT1, AT2, AT3]=ind2sub([n n n],ATind);
        [AV1, AV2, AV3]=ind2sub([n n n],AVind);
        figure('name','Adaptive immune cells')
        scatter3(AT1,AT2,AT3,'MarkerFaceColor',[0 1 0])
        hold on
        scatter3(AV1,AV2,AV3,'MarkerFaceColor',[0.66 0 0])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        legend('Antitumor','Antiviral')    
        filename = [dir_name '/plot3d_adaptive'];
        savefig(filename)
        filename = [filename '.jpg'];
        saveas(gcf,filename)
    % 

     %% Plots of total pop'n vs time
        figure('name','Time: susc cells')
        plot(AllTypes(:,1),AllTypes(:,2),'linewidth',3)
        xlabel('hours')
        ylabel('Susceptible tumor cells')
        set(gca,'FontSize',16);
        filename = [dir_name '/time_susc'];
        savefig(filename)
        filename = [filename '.jpg'];
        saveas(gcf,filename)
    %     
        figure('name','Time: infected cells')
        plot(AllTypes(:,1),AllTypes(:,3),'linewidth',3,'Color',[1 0 0])
        xlabel('hours')
        ylabel('Infected tumor cells')
        set(gca,'FontSize',16);
        filename = [dir_name '/time_infected'];
        savefig(filename)
        filename = [filename '.jpg'];
        saveas(gcf,filename)
    %     
        figure('name','Time: innate immune cells')
        plot(AllTypes(:,1),AllTypes(:,4),'linewidth',3,'Color',[1 0 1])
        xlabel('hours')
        ylabel('Innate immune cells')
        set(gca,'FontSize',16);
        filename = [dir_name '/time_innate'];
        savefig(filename)
        filename = [filename '.jpg'];
        saveas(gcf,filename)
    %     
        figure('name','Time: antitumor adaptive immune')
        plot(AllTypes(:,1),AllTypes(:,5),'linewidth',3,'Color',[0 1 0])
        xlabel('hours')
        ylabel('Antitumor adaptive immune cells')
        set(gca,'FontSize',16);
        filename = [dir_name '/time_antitumor_adap'];
        savefig(filename)
        filename = [filename '.jpg'];
        saveas(gcf,filename)

        figure('name','Time: antiviral adaptive immune')
        plot(AllTypes(:,1),AllTypes(:,6),'linewidth',3,'Color',[0.66 0 0])
        xlabel('hours')
        ylabel('Antiviral adaptive immune cells')
        set(gca,'FontSize',16);
        filename = [dir_name '/time_antiviral_adap'];
        savefig(filename)
        filename = [filename '.jpg'];
        saveas(gcf,filename)

        figure('name','Time: virus')
        plot(AllTypes(:,1),AllTypes(:,8),'linewidth',3,'Color',[0 .75 .75])
        xlabel('hours')
        ylabel('Virus (pfu)')
        set(gca,'FontSize',16);
        filename = [dir_name '/time_virus'];
        savefig(filename)
        filename = [filename '.jpg'];
        saveas(gcf,filename)

        if antiPD1
            figure('name','Time: anti-PD-1')
            plot(AllTypes(:,1),AllTypes(:,7),'linewidth',3)
            xlabel('hours')
            ylabel('Anti-PD-1 (\mu mol/cm^3)')
            set(gca,'FontSize',16);
            filename = [dir_name '/time_antiPD1'];
            savefig(filename)
            filename = [filename '.jpg'];
            saveas(gcf,filename)
        end
    end
end

if totRep>1
    filename = [dir_name '/finalPops'];
    save(filename,'finalPops')
    filename = [dir_name '/allSimSusc'];
    save(filename,'allSimSusc')
end
