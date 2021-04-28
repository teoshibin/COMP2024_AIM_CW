%% SETTING JADE+DE - Template
clear;

%-------- setExperiment --------------------
bbobpath = getBBOBRootPath('bbob-local');
setExperiment.datapath = fullfile(bbobpath, 'results', 'DEAE5C8'); % path with the experiment data
setExperiment.name = 'DEAE CCOVFAC study';     % name the experiment e.g.: 'DEAE with new feature'
setExperiment.comment = 'DE pop5DIM F~U(0.5,1) CR=0.5 AE CCOVFAC 8';
setExperiment.maxfunevals = '5e3'; %  e.g.: '10000 * dim'

setExperiment.dim = [2,3,5,10,20,40]; % e.g.: [2,3,5,10,20,40]
setExperiment.benchmark_fun = 'benchmarks(''FunctionIndices'')'; %e.g.: 'benchmarks(''FunctionIndices'')' or 'benchmarksnoisy(''FunctionIndices'')'
setExperiment.useFunctions = 'all';% e.g. 'all' or for noiseless:[1-24] or for noisy[101-130];
setExperiment.instance = [1:15];     % #instances = #trials e.g.: [1:15] for BBOB 2012: [1:5 21:30]

%-------- setDE --------------------------
setDE.popsize = '5*DIM'; %e.g.: '6*DIM' e.g.: '4 + floor(3 * log(DIM))'
setDE.maxpopsize = inf;%e.g.: 80

setDE.minvarcondition =  1e-10;%e.g.: 1e-10;       % minimum variance which don't restart the algorithm
setDE.stuckcond_noImp = inf;%e.g.:50;          % maximum interations without improvement after which algorithm is restarted
setDE.stuckcond_lowVar = 30;%e.gl: 30; 

setDE.mutation_op = 'best';   % type of mutation operator: ['rand' or 'best' or 'average' or 'JA']
setDE.crossover_op = 'bin';   % type of mutation operator: ['bin' or 'exp' of 'JA' or 'non']
setDE.CR = 0.5;    %e.g.: 0.5

% switch using statistics graphs
setDE.statsDE = false; % show  in a graph at the end of every trial
setJA.statsJA = false; % show mu_CR, mu_F, p_AE in a graph at the end of every trial

%-------- setJA ---------------------------
setJA.adaptJA_muF = false;  % use adaptation of mu_F (mutation param.)
setJA.adaptJA_muCR = false; % use adaptation of mu_CR (crossover param.)
setJA.adaptJA_pAE = false;  % use adaptation of p_AE 

setJA.JA_init_mu_F = 0.5;        % init value for mu_F
setJA.JA_init_mu_CR = 0.5;       % init value for mu_CR
setJA.JA_init_p_AE = 0.0;        % init value for p_AE

setJA.JA_pArch = 0.1;   % how much best individuals to store in archive (percentage, recommanded values [0.05 - 0.20]; 0 = switch off)
setJA.JA_pMut = 0.1;    % from 100*JA_pMut% of best ind. are randomly picked ind. in JA-mutation 
                        % (percentage, recommanded same as p_pArch but its min. value is fixed to 1 individual)
setJA.JA_c_CR = 0.1; % muCR learning rate (recommanded values [0.2 - 0.05]) i.e. life span from 5 to 20 generations
setJA.JA_c_F = 0.1;  % muF learning rate (recommanded values [0.2 - 0.05]) i.e. life span from 5 to 20 generations
setJA.JA_c_AE = 0.1; % pAE learning rate (for start same value as c_CR and c_F)

setJA.JA_p_AE_min = 0.05;   % pAE min value
setJA.JA_p_AE_max = 0.95;   % pAE max value

%-------- setAE ---------------------------
setDE.useAE = true;    % To use DE: set false, to use DE+AE: set true

setAE.mu = 'ceil(popsize/2)'; %e.g.: 'ceil(popsize/2)'
setAE.CCOVFAC = 8;
setAE.c_p = '1/sqrt(AE.DIM)';    % learning rate for the evolution path e.g.: [] e.g.: '1/sqrt(AE.DIM)'
setAE.c_1 = '0.1';   % learning rate for rank-ONE-update (evolution path) e.g.: [] e.g.: '0.1' e.g.: '(AE.alp_c*0.2) / ((DIM+1.3)^2+AE.mu_w)'
setAE.c_mu = '0.1';   % learning rate for rank-MU-update (covariance matrix estimation) e.g.: [] e.g.: '0.1' e.g.: 'alp_mu = 0.2; (AE.alp_c * 0.2 *(AE.mu_w-2+(1/AE.mu_w))) / ((DIM+2)^2+alp_mu*AE.mu_w)'

%% Copy the m-files to the results
copypath = fullfile(setExperiment.datapath,'matlab','');
mkdir(copypath);
copyfile('*.m', copypath);
% newdiary(fullfile(copypath,'diary'));


% try
   experimentJADEAE_fn(setExperiment, setDE, setAE, setJA);
% catch exception
%    fprintf('The experiment: "%s" has FAILED:\n\n',setExperiment.name);
%    disp(exception.message);
%    
%    
%    filename_1 = fullfile(setExperiment.datapath, '00_FAILURE.txt');
%    file_1 = fopen(filename_1,'w');
%    fprintf(file_1,'The experiment:  "%s" has FAILED:\n\n',setExperiment.name);
%    fprintf(file_1,'%s\n',exception.message);
%    fclose(file_1);
%    
%    filename_2 = fullfile(setExperiment.datapath, 'log.txt');
%    file_2 = fopen(filename_2,'a');
%    fprintf(file_2,'The experiment:  "%s" has FAILED:\n\n',setExperiment.name);
%    fprintf(file_2,'%s\n',exception.message);
%    fclose(file_2);
% 
% end

diary off;
