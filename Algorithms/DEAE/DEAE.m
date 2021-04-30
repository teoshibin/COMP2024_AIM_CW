% Differential Evolution with Adaptive Encoding Algorithm
%
% DE with adaptive encoding is known to cater for rotational invariance 
% especially when it comes to crossover of non-separable functions. 
% This differs such that the offspring creation takes places 
% in a different space provided by a linear transformation. 
%
% This code is a refactored version of the source code available online.
% Please find the online source code citation as shown below.
%%%--------------------------------------------------------------------------------------------------%%%
%    Title: DEAE Source Code
%    Author: Pošík and Klemš
%    Date: 2012
%    Code version: 1
%    Availability: https://coco.gforge.inria.fr/doku.php?id=bbob-2012-algorithms
%%%--------------------------------------------------------------------------------------------------%%%
% (Version 1)[Source Code].https://coco.gforge.inria.fr/doku.php?id=bbob-2012-algorithms

function [exitcode] = DEAE(FUN, DIM, ftarget, maxfunevals)
    % JADEAE_fn(FUN, {dim, setDE, setAE, setJA}, ftarget, maxfunevals)
    
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
    %------------------------------------------
        
    
    % init return values
    exitcode = 'should not return this';
    
    maxfunevals = min(1e8 * DIM, maxfunevals); % set maximum number of fitness evaluations
    
    popsize = min(eval(setDE.popsize), setDE.maxpopsize);   % set parental population size e.g. 4 + floor(3 * log(DIM))
    offspringsize = popsize; % set offspring size
    
    offspring = zeros(DIM,popsize); %initialize auxilary data structure
    
    % Setup restarts
    stuckcount = 0;
    minvarcondition = setDE.minvarcondition;          % minimum variance which don't restart the algorithm
    stuckcond_noImp = setDE.stuckcond_noImp;          % maximum interations without improvement after which algorithm is restarted
    stuckcond_lowVar = setDE.stuckcond_lowVar;        % max. iter without improvement + low diversity condition
    
    crossover_op = setDE.crossover_op;   % type of mutation operator: ['bin' or 'exp']
    mutation_op = setDE.mutation_op;   % type of mutation operator: ['rand' or 'best' or 'average']
    
    %F = 0.8;             % set parameter
    CR = setDE.CR;      % set cross-over probability
    
    useAE = setDE.useAE;       % USE ADAPTIVE ENCODING
    
    statsDE = setDE.statsDE;   % show variance in population, best solution and time per generations in a graph at the end of every trial
    statsJA = setJA.statsJA;   % show mu_CR, mu_F, p_AE in a graph at the end of every trial
    
    JA_pArch = setJA.JA_pArch;          % how much best individuals to store (percentage, 0 = switch off)
    JA_pMut = setJA.JA_pMut;% from 100*JA_pMut% of best ind. are randomly picked ind. in JA-mutation
    % (percentage, recommanded same as p_pArch but its min. value is fixed to 1 individual)
    JA_c_CR = setJA.JA_c_CR;    % mu_CR learning rate
    JA_c_F = setJA.JA_c_F;      % mu_F learning rate
    JA_c_AE = setJA.JA_c_AE;    % p_AE learning rate
    JA_p_AE_min = setJA.JA_p_AE_min;   % pAE min value
    JA_p_AE_max = setJA.JA_p_AE_max;   % pAE max value
    
    pop_de = 10 * rand(DIM, popsize) - 5;  % initialize population randomly within [-5, 5] in every dimension
    popfit = feval(FUN, pop_de);
    [popfit, index] = sort(popfit);        % sort population
    pop_de = pop_de(:,index);              % sort population
    fbest = popfit(1);
    fbestold = fbest;
    
    maxiter = ceil( (maxfunevals-popsize)/popsize );
    
    % init archive - if nBestToArchive = 0 - then archive is switched off
    nBestToArchive = min(max(0, ceil(popsize * JA_pArch)), popsize);
    archive = [];
    auxArchive = [];
    newmembers = [];
    newmembersfit = [];
    
    % init JA-mutation
    adaptJA_muF = setJA.adaptJA_muF;    % use adaptation of mu_F
    adaptJA_muCR = setJA.adaptJA_muCR;  % use adaptation of mu_CR
    adaptJA_pAE = setJA.adaptJA_pAE;    % use adaptation of p_AE
    
    mu_F = setJA.JA_init_mu_F;    % init value for mu_F
    mu_CR = setJA.JA_init_mu_CR;	% init value for mu_CR
    p_AE = setJA.JA_init_p_AE;    % init value for p_AE
    
    nBestToMutate = min(max(1, ceil(popsize * JA_pMut)), popsize);
    randvectCR = ones(1,popsize)*mu_CR;
    randvectF = ones(1,popsize)*mu_F;
    
    if useAE
        muae = eval(setAE.mu);
%         AE = AdaptiveEncoding(DIM, popsize, eval(setAE.mu), pop_de);
        AE = AEupdate([], pop_de(:,1:muae), setAE.CCOVFAC);
%         if ~isempty(setAE.c_1)
%             AE.c_1 = eval(setAE.c_1);
%         end
%         if ~isempty(setAE.c_mu)
%             AE.c_mu = eval(setAE.c_mu);
%         end
%         if ~isempty(setAE.c_p)
%             AE.c_p = eval(setAE.c_p);
%         end
    end
    
    if (strcmp(mutation_op, 'average'))
        % set recombination weight
        mut_weigh = zeros(popsize,1);
        mut_mu = ceil(popsize/2);
        
        ln_sum = sum(log(1:mut_mu));
        for i = 1:mut_mu
            mut_weigh(i) = (log(mut_mu+1) - log(i)) / (mut_mu*log(mut_mu+1) - ln_sum);
        end
        
    end
    
    if statsDE
        stat_var = zeros(1,maxiter);
        stat_best = zeros(1,maxiter);
        stat_time = zeros(1,maxiter);
        tic
    end
    
    if (statsJA == true)
        stat_muCR = zeros(1, maxiter+1);
        stat_muF = zeros(1, maxiter+1);
        stat_pAE = zeros(1, maxiter+1);
        
        stat_muCR(1) = mu_CR;
        stat_muF(1) = mu_F;
        stat_pAE(1) = p_AE;
    end
    
    for iter = 1:maxiter
        
        % ENCODING 1st part - except of case when we encode only some individuals
        if (useAE && ~adaptJA_pAE)
            isEncoded = ones(popsize,1);
%             pop_en = AE.encode(pop_de);
            pop_en = AE.invB * pop_de;
        else
            isEncoded = zeros(popsize, 1);
            pop_en = pop_de;
        end
        
        % MUTATION - rotationaly independent
        F = rand/2 + 0.5; % randomize F for every generation within range [0.5, 1]
        switch mutation_op
            case 'rand'
                randvectMU = ceil(rand(3,popsize)*popsize);   % may improve when applaying different numbers in every column
                for i = 1:popsize
                    %F = rand/2 + 0.5; % randomize F for every vector within range [0.5, 1]
                    offspring(:,i) = pop_en(:,randvectMU(1,i)) + F *( pop_en(:,randvectMU(2,i)) - pop_en(:,randvectMU(3,i)) );
                end
            case 'best'
                randvectMU = ceil(rand(2,popsize)*popsize);
                randvectMU(2,randvectMU(1,:) ==randvectMU(2,:)) = mod(randvectMU(2,randvectMU(1,:) ==randvectMU(2,:))+1,popsize);
                randvectMU(2,randvectMU(2,:)==0) = popsize;
                for i = 1:popsize
                    %F = rand/2 + 0.5; % randomize F for every vector within range [0.5, 1]
                    offspring(:,i) = pop_en(:,1) + F *( pop_en(:,randvectMU(1,i)) - pop_en(:,randvectMU(2,i)) );
                end
            case 'average'
                avg_ind = pop_en*mut_weigh;
                
                randvectMU = ceil(rand(2,popsize)*popsize);
                randvectMU(2,randvectMU(1,:) ==randvectMU(2,:)) = mod(randvectMU(2,randvectMU(1,:) ==randvectMU(2,:))+1,popsize);
                randvectMU(2,randvectMU(2,:)==0) = popsize;
                for i = 1:popsize
                    %F = rand/2 + 0.5; % randomize F for every vector within range [0.5, 1]
                    offspring(:,i) = avg_ind + F *( pop_en(:,randvectMU(1,i)) - pop_en(:,randvectMU(2,i)) );
                end
            case 'JA'
                % random F
                generateSize = popsize*2;
                countInRange = 0;
                randvectF = zeros(1,popsize);
                for i = 1:50 % not probable to take so many loops
                    randc = 0.1 * tan(pi*(rand(1,generateSize)-0.5)) + mu_F;
                    randc = randc(randc > 0);
                    if (countInRange + length(randc) < popsize)
                        randvectF(countInRange + 1 : countInRange + length(randc)) = randc;
                        countInRange = countInRange + length(randc);
                    else
                        randvectF((countInRange+1) : popsize) = randc(1:popsize-countInRange );
                        break;
                    end
                end
                randvectF(randvectF > 1) = 1;
                F = repmat(randvectF, DIM, 1);
                
                % random R1 and R2
                sequence = 1:popsize;
                
                randvectR1 = ceil(rand(1,popsize) * popsize);
                randvectR1(randvectR1 == sequence) = mod(randvectR1(randvectR1 == sequence) + 1, popsize);
                randvectR1(randvectR1 == 0) = popsize;
                
                archiveCount = size(archive,2);
                sizePandA = popsize + archiveCount;
                randvectR2 = ceil(rand(1,popsize) * sizePandA);
                randvectR2(randvectR2 == sequence) = mod(randvectR2(randvectR2 == sequence) + 1, sizePandA);
                randvectR2(randvectR2 == 0) = sizePandA;
                randvectR2(randvectR2 == randvectR1) = mod(randvectR2(randvectR2 == randvectR1) + 1, sizePandA);
                randvectR2(randvectR2 == 0) = sizePandA;
                
                % random from nBestToMutate individuals from population
                randvectB = randi([1, nBestToMutate], 1, popsize);
                
                % append archive
                if(archiveCount > 0)
                    aux = [pop_en, archive(:,1:archiveCount)];
                else
                    aux = pop_en;
                end
                
                % create donors
                offspring = pop_en + F .* ( pop_en(:,randvectB) - pop_en ) + F .* (pop_en(:,randvectR1) - aux(:,randvectR2));
        end
        
        % ENCODING 2nd part - if we encode only some individuals
        if (useAE && adaptJA_pAE && 0 < p_AE)
            isEncoded = rand(popsize, 1) <= p_AE;
%             pop_en(:,isEncoded) = AE.encode(pop_en(:,isEncoded));
            pop_en(:,isEncoded) = AE.invB * pop_en(:,isEncoded);
%             offspring(:,isEncoded) = AE.encode(offspring(:,isEncoded));
            offspring(:,isEncoded) = AE.invB * offspring(:,isEncoded);
        end
        
        % CROSS-OVER - rotationaly dependent
        switch crossover_op
            case 'bin'
                randvect = rand(DIM,popsize);
                offspring(randvect > CR) = pop_en(randvect > CR);    % new element from offspring is accepted when rand(0,1) < CR
            case 'exp'
                for i = 1:popsize
                    beginpos = randi([1,DIM], 1);
                    maxL = randi([1,DIM], 1);
                    pbbL = rand(maxL, 1);
                    L = find(pbbL>CR);
                    if isempty(L)
                        L = maxL;
                    else
                        L = L(1);
                    end
                    uvector = pop_en(:,i);
                    uvector(beginpos:min(beginpos+L-1,DIM)) = offspring(beginpos:min(beginpos+L-1,DIM),i);
                    uvector(1:mod(beginpos+L-1,DIM)) = offspring(1:mod(beginpos+L-1,DIM),i);
                    offspring(:,i) = uvector;
                end
            case 'JA'
                randvectCR = mu_CR + 0.1*randn(1,popsize);
                randvectCR(randvectCR < 0) = 0;
                randvectCR(randvectCR > 1) = 1;
                randmatCR2 = repmat(randvectCR, DIM, 1);
                randmatCR = rand(DIM, popsize);
                randmatCR = randmatCR < randmatCR2;
                
                randvectJ  = randi(DIM,1,popsize);
                randmatJ = zeros(DIM, popsize);
                randmatJ(sub2ind(size(randmatJ),randvectJ,1:popsize)) = 1;
                
                offspring(~(randmatCR | randmatJ)) = pop_en(~(randmatCR | randmatJ));
            case 'non'
                %don't do crossover
                
        end
        
        % OFFSPRING DECODING AND EVALUATION
        if (useAE && ~adaptJA_pAE)
%         if (useAE && p_AE == 1) % decode all
%             offspring = AE.decode(offspring);
            offspring = AE.B * offspring;
        elseif(useAE) % decode only encoded
%             offspring(:,isEncoded) = AE.decode(offspring(:,isEncoded));
            offspring(:,isEncoded) = AE.B * offspring(:,isEncoded);
        end
        offspringfit = feval(FUN,offspring);
        
        % SELECTION
        % find better individuals
        ind = find(popfit>offspringfit);
        
        % save into archive
        %    if (nBestToArchive > 0)
        %       for i = ind
        %          if (archiveCount < nBestToArchive)
        %             archiveCount = archiveCount + 1;
        %             randPosArchive = archiveCount;
        %          else
        %             randPosArchive = randi(archiveCount);
        %          end
        %
        %          archive(:,randPosArchive) = pop_de(:,i);
        %       end
        %    end
        if (nBestToArchive > 0 && ~isempty(ind))
            archiveCount = size(archive,2);
            auxArchive(:,1:(archiveCount + length(ind))) = [archive, pop_de(:,ind)];
            
            p = randperm(archiveCount + length(ind));
            p = p(1:min(nBestToArchive, archiveCount + length(ind)));
            
            archive = auxArchive(:,p);
        end
        
        % overwrite by better individuals from offspring
        pop_de(:,ind) = offspring(:,ind);
        popfit(ind) = offspringfit(ind);
        newmembers = [newmembers offspring(:,ind)];
        newmembersfit = [newmembersfit offspringfit(ind)];
        
        % ADAPT ADAPTABLES
        if (adaptJA_muCR)
            meanA = mean(randvectCR(ind));
            if(~isnan(meanA))
                mu_CR = (1-JA_c_CR)*mu_CR + JA_c_CR*meanA;
            end
        end
        if (adaptJA_muF)
            meanL = sum(sum(randvectF(ind).^2))./sum(sum(randvectF(ind)));
            if(~isnan(meanL))
                mu_F = (1-JA_c_F)*mu_F + JA_c_F*meanL;
            end
        end
        if (adaptJA_pAE)
            meanAE = mean(isEncoded(ind));
            if(~isnan(meanAE))
                p_AE = (1-JA_c_AE)*p_AE + JA_c_AE*meanAE;
                p_AE = max(p_AE, JA_p_AE_min);
                p_AE = min(p_AE, JA_p_AE_max);
            end
        end
        
        % ORDER POPULATION BELONG FITNESS - because of adaptive encoding
        %                                   and average mutation
        [popfit, index] = sort(popfit);
        pop_de = pop_de(:,index);
        fbest = popfit(1);
        
        [newmembersfit, index] = sort(newmembersfit);
        newmembers = newmembers(:,index);
        
        % STOPPING CRITERIA
        % interim stopping (algorithm is stucked in local optima)
        if fbest < fbestold
            % Improvement detected
            stuckcount = 0;
            fbestold = fbest;
        else
            stuckcount = stuckcount + 1;
            if stuckcount >= stuckcond_noImp
                exitcode = 'no improvement';
%                 break;
            end
        end
        %interim stopping (too low diversity)
        if (stuckcount > stuckcond_lowVar) && (sum(var(pop_de,1,2))/DIM < minvarcondition)
%         if (max(std(pop_de,1,2)) < minvarcondition)
            exitcode = 'low variance';
            break;
        end
        % ftarget is reached
        if feval(FUN, 'fbest') < ftarget  % task achieved
            exitcode = 'solution found';
            break;
        end
        
        % DEBBUG AND STATS
        if statsDE
            stat_var(iter) = sum(var(pop_de,0,2));
            stat_best(iter) = feval(FUN, 'fbest') - ftarget;
            stat_time(iter) = toc;
        end
        
        if statsJA
            stat_muCR(iter+1) = mu_CR;
            stat_muF(iter+1)  = mu_F;
            stat_pAE(iter+1)  = p_AE;
        end
        
        % ADAPT AE
        if useAE
%             AE = AE.update(pop_de);
            if length(newmembersfit) >= muae,
                AE = AEupdate(AE,newmembers(:,1:muae));
                newmembers = [];
                newmembersfit = [];
            end
        end
    end
    
    %x = feval(FUN, 'fbest');
    
    if statsDE
        figure(1);
        plot(stat_var(1:iter-1));
        title('Variance in population');
        figure(2);
        plot(stat_best(1:iter-1));
        titstr = strcat ('Best found solution: ', num2str(feval(FUN, 'fbest') - ftarget));
        title(titstr);
        figure(3);
        plot(stat_time(1:iter-1));
        title('Computation time')
        
        %input('Press space to continue: ');
        %pause(1);
    end
    
    if statsJA
        figure(4);
        plot(stat_muCR(1:iter),'r');
        hold on
        title('Development of adaptation: mu_{CR} mu_{F} p_{AE}');
        plot(stat_muF(1:iter),'g');
        plot(stat_pAE(1:iter),'b');
        legend('mu_{CR}', 'mu_{F}', 'p_{AE}');
        hold off
        %input('Press space to continue: ');
        %pause(0);
    end
