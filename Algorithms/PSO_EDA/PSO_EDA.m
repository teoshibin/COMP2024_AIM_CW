% Hybrid of Estimation Distribution Algorithm and Partical Swarm
% Optimization
%  
% EDA-PSO
%
% This code is a refactored version of the source code available online.
% Please find the online source code citation as shown below.
%%%--------------------------------------------------------------------------------------------------%%%
%    Title: Black-Box Optimization Benchmarking for Noiseless
%           Function Testbed using an EDA and PSO Hybrid (PDF)
%    Author: El-Abd and Kamel
%    Date: 2009
%    Code version: 1
%    Availability: https://coco.gforge.inria.fr/doku.php?id=bbob-2009-algorithms
%%%--------------------------------------------------------------------------------------------------%%%

function PSO_EDA(FUN, DIM, ftarget, maxfunevals)

    % Set algorithm parameters
    popsize = 40;
    c1 = 2;
    c2 = 2;
    xbound = 5;
    vbound = 5;
    p = 0.5;
    num_PSO = 0;
    tot_PSO = 0;
    num_EDA = 0;
    tot_EDA = 0;
    C_AVS = 1;
    C_AVS_Max = 10;
    C_AVS_Min = 0.1;
    E_Dec = 0.9;
    E_Inc = 1 / E_Dec;
    
    % Allocate memory and initialize
    xmin = -xbound * ones(1,DIM);
    xmax = xbound * ones(1,DIM);
    vmin = -vbound * ones(1,DIM);
    vmax = vbound * ones(1,DIM);
    x = 2 * xbound * rand(popsize,DIM) - xbound;
    v = 2 * vbound * rand(popsize,DIM) - vbound;
    pbest = x;
    
    % update pbest and gbest
    cost_x = feval(FUN, x');
    cost_p = cost_x;
    [cost,index] = min(cost_p);
    gbest = pbest(index,:);
    maxfunevals = min(1e5 * DIM, maxfunevals);
    maxiterations = ceil(maxfunevals/popsize);
    
    for iter = 2 : maxiterations
        % PSO part
        % Update inertia weight
        w = 0.9 - 0.8*(iter-2)/(maxiterations-2);

        prev_num_EDA = num_EDA;

        % Update velocity
        candidate_v = w*v + c1*rand(popsize,DIM).*(pbest-x) + c2*rand(popsize,DIM).* (repmat(gbest,popsize,1)-x);

        % Clamp veloctiy
        s = candidate_v < repmat(vmin,popsize,1);
        candidate_v = ~s.*candidate_v + s.*repmat(vmin,popsize,1);
        b = candidate_v > repmat(vmax,popsize,1);
        candidate_v = ~b.*candidate_v + b.*repmat(vmax,popsize,1);

        % Update position
        candidate_x = x + candidate_v;

        % Clamp position - Absorbing boundaries
        % Set candidate x to the boundary
        s = candidate_x < repmat(xmin,popsize,1);
        candidate_x = ~s.*candidate_x + s.*repmat(xmin,popsize,1);
        b = candidate_x > repmat(xmax,popsize,1);
        candidate_x = ~b.*candidate_x + b.*repmat(xmax,popsize,1);

        % Clamp position - Absorbing boundaries
        % Set candidate v to zero
        b = s | b;
        candidate_v = ~b.*candidate_v + b.*zeros(popsize,DIM);

        % EDA part
        % Calculate Gaussian model where Mus and Sigma
        % are based on the best half of the swarm
        [cost_p_sorted,indices] = sort(cost_p);
        Mue = mean(pbest(indices(1:popsize/2),:));
        Sigma = std(pbest(indices(1:popsize/2),:));
        % Calculate correlation
        D = pbest(indices(1:popsize/2),:) - repmat(Mue,popsize/2,1);
        Distances = sum(abs(D'));
        Fitness = cost_p_sorted(1:popsize/2);
        Rho = corr(Distances', Fitness', 'type', 'spearman');
        if(Rho>-0.55)
            Sigma = sqrt(C_AVS) * Sigma;
        end

        % Generate candidate solution
        candidate_EDA_x = repmat(Mue,popsize,1) + repmat(Sigma,popsize,1).*randn(popsize,DIM);

        % Depending on which component to choose
        % select candidates for consideration and
        % update the equivelant counters
        r = rand(popsize,1)<p;
        R = repmat(r,1,DIM);
        candidates = ~R.*candidate_EDA_x + R.*candidate_x;
        candidates_fitness = feval(FUN, candidates');
        tot_PSO = tot_PSO + sum(r);
        tot_EDA = tot_EDA + popsize - sum(r);

        % Update x if candidates are better
        c = candidates_fitness<cost_x;
        C = repmat(c',1,DIM);
        x = ~C.*x + C.*candidates;
        v = ~(R&C).*v + (R&C).*candidate_v;
        cost_x = ~c.*cost_x + c.*candidates_fitness;

        % Update success counters
        num_PSO = num_PSO + sum(r&c');
        num_EDA = num_EDA + sum((~r)&c');

        % Update pbest if necessary
        c = cost_x<cost_p;
        C = repmat(c',1,DIM);
        pbest = ~C.*pbest + C.*x;
        cost_p = ~c.*cost_p + c.*cost_x;
        % Update gbest if necessary
        [cost,index] = min(cost_p);
        gbest = pbest(index,:);

        % Update C_AVS based on the
        % EDA component performance
        if(num_EDA>prev_num_EDA)
            C_AVS = E_Inc * C_AVS;
        else
            C_AVS = E_Dec * C_AVS;
        end
        if((C_AVS<C_AVS_Min)||(C_AVS>C_AVS_Max))
            C_AVS = C_AVS_Max;
        end

        % Update probability using
        % percentage of improvements
        PSO_Imp_Perc = num_PSO / tot_PSO;
        EDA_Imp_Perc = num_EDA / tot_EDA;
        p = PSO_Imp_Perc/(PSO_Imp_Perc+EDA_Imp_Perc);

        % Exit if target is reached
        if feval(FUN, 'fbest') < ftarget
            break;
        end
    end 
end
