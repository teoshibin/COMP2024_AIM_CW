% Partical Swarm Optimization with population-based incremental learning [1]
% 
% PSO_Bounds (PSO with PIBL)
%
% This code is a refactored version of the source code available online.
% Please find the online source code citation as shown below.
%%%--------------------------------------------------------------------------------------------------%%%
%    Title: Black-Box Optimization Benchmarking for Noiseless
%           Function Testbed using PSO_Bounds (PDF)
%    Author: El-Abd and Kamel
%    Date: 2009
%    Code version: 1
%    Availability: https://coco.gforge.inria.fr/doku.php?id=bbob-2009-algorithms
%%%--------------------------------------------------------------------------------------------------%%%

function PSO_Bounds(FUN, DIM, ftarget, maxfunevals)

    % Set algorithm parameters
    popsize = 40;
    c1 = 2;
    c2 = 2;
    xbound = 5;
    vbound = 5;
    pmin = 0.2;
    pmax = 0.8;
    alpha = 0.05;
    T = 0.0001;
    
    % Allocate memory and initialize
    xmin = -xbound * ones(1,DIM);
    xmax = xbound * ones(1,DIM);
    vmin = -vbound * ones(1,DIM);
    vmax = vbound * ones(1,DIM);
    x = 2 * xbound * rand(popsize,DIM) - xbound;
    v = 2 * vbound * rand(popsize,DIM) - vbound;
    pbest = x;
    p = 0.5 * ones(1,DIM);
    
    % update pbest and gbest
    cost_p = feval(FUN, pbest');
    [cost,index] = min(cost_p);
    gbest = pbest(index,:);
    maxfunevals = min(1e5 * DIM, maxfunevals);
    maxiterations = ceil(maxfunevals/popsize);
    
    for iter = 2 : maxiterations

        % Update inertia weight
        w = 0.9 - 0.8*(iter-2)/(maxiterations-2);

        % Update velocity
        v = w*v + c1*rand(popsize,DIM).*(pbest-x) + c2*rand(popsize,DIM).*(repmat(gbest,popsize,1)-x);

        % Clamp veloctiy
        s = v < repmat(vmin,popsize,1);
        v = ~s.*v + s.*repmat(vmin,popsize,1);
        b = v > repmat(vmax,popsize,1);
        v = ~b.*v + b.*repmat(vmax,popsize,1);

        % Update position
        x = x + v;

        % Clamp position - Absorbing boundaries
        % Set x to the boundary
        s = x < repmat(xmin,popsize,1);
        x = ~s.*x + s.*repmat(xmin,popsize,1);
        b = x > repmat(xmax,popsize,1);
        x = ~b.*x + b.*repmat(xmax,popsize,1);

        % Clamp position - Absorbing boundaries
        % Set v to zero
        b = s | b;
        v = ~b.*v + b.*zeros(popsize,DIM);
        % Update pbest and gbest if necessary
        cost_x = feval(FUN, x');
        s = cost_x<cost_p;
        cost_p = ~s.*cost_p + s.*cost_x;
        s = repmat(s',1,DIM);
        pbest = ~s.*pbest + s.*x;
        [cost,index] = min(cost_p);
        gbest = pbest(index,:);

        % Update dimension probability
        probability = sum(pbest>repmat(((xmin+xmax)/2),popsize,1));
        p = (1-alpha)*p + alpha*(probability/popsize);

        % Update bounds if necessary
        % Shift xmax
        pmn = p<pmin;
        xmax = ~pmn.*xmax + pmn.*(xmax - (xmax-xmin)/2);

        % Shift xmin
        pmx = p>pmax;
        xmin = ~pmx.*xmin + pmx.*(xmin + (xmax-xmin)/2);

        % In either case, set p to 0.5
        pm = pmn | pmx;
        p = ~pm.*p + pm*0.5;

        % Re-initialize if necessary
        t = (xmax-xmin)<(2*T*xbound);
        xmin = ~t.*xmin + t*-xbound;
        xmax = ~t.*xmax + t*xbound;
        vmax = ~t.*((xmax-xmin)/2) + t.*vbound;
        vmin = -vmax;
        t = repmat(t,popsize,1);
        v = (1-t).*v + t.*(2 * vbound * rand(popsize,DIM) - vbound);

        % Exit if target is reached
        if feval(FUN, 'fbest') < ftarget
            break;
        end
    end
end


