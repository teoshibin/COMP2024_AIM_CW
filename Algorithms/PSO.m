 function PSO(FUN, DIM, ftarget, maxfunevals)
% MY_OPTIMIZER(FUN, DIM, ftarget, maxfunevals)
% samples new points uniformly randomly in [-5,5]^DIM
% and evaluates them on FUN until ftarget of maxfunevals
% is reached, or until 1e8 * DIM fevals are conducted. 
% Relies on FUN to keep track of the best point. 

% FUN - benchmark function
% DIM - dimansions
% ftarget - close enough to solution (not for the algorithm but for termination)
% maxfunevals - maximum fitness function evaluation

    % Set algorithm parameters
    popsize = 40;   % population size
    c1 = 1.4944;    % local exploration
    c2 = 1.4944;    % global exploitation
    w = 0.792;      % local exploitation ~ global exploration
    xbound = 5;     % for [-5, 5]^D
    vbound = 5;     % for velocity boundary which is [-5, 5]^D as well
    
    % Initialize
    xmin = -xbound * ones(1,DIM);   % min search space
    xmax = xbound * ones(1,DIM);    % max search space
    vmin = -vbound * ones(1,DIM);   % minimum velocity
    vmax = vbound * ones(1,DIM);    % maximum velocity
    x = 2 * xbound * rand(popsize,DIM) - xbound; % swarm positions
    v = 2 * vbound * rand(popsize,DIM) - vbound; % swarm velocity
    pbest = x; % initial swarm personal best = initial positions
    
    % init pbest and gbest
    cost_p = feval(FUN, x');    % calculate population fitness
    [~,index] = min(cost_p);        % get minimum population fitness
    gbest = pbest(index,:);         % set gbest as best population fitness
    maxfunevals = min(1e5 * DIM, maxfunevals); % maxfunevals cannot be larger than 1*10^5 * DIM
    
    % each individual will complete one eval each iter
    % so evalperiter * iter = totaleval where evalperiter = popsize
    maxiterations = ceil(maxfunevals/popsize); 
    
    for iter = 2 : maxiterations
        
        % Update inertia weight dynamically
        % w = 0.9 - 0.8*(iter-2)/(maxiterations-2);

        % Calculate new velocity
        % v = w*v + c1*r*pbest_distance + c2*r*gbest_distance
        v = w*v + c1*rand(popsize,DIM).*(pbest-x) + c2*rand(popsize,DIM).*(repmat(gbest,popsize,1)-x);

        % Restrict veloctiy
        lt = v < repmat(vmin,popsize,1); % find logical matrix for v that is smaller than min v
        v = ~lt.*v + lt.*repmat(vmin,popsize,1); % replace v with min v that is smaller than min v
        gt = v > repmat(vmax,popsize,1); % find logical matrix for v that is larger than max v
        v = ~gt.*v + gt.*repmat(vmax,popsize,1); % replace v with max v that is larger than max v

        % Update position
        x = x + v;

        % Restrict position
        lt = x < repmat(xmin,popsize,1);
        x = ~lt.*x + lt.*repmat(xmin,popsize,1);
        gt = x > repmat(xmax,popsize,1);
        x = ~gt.*x + gt.*repmat(xmax,popsize,1);

        % Reset restructed particals' velocity
        glt = lt | gt; % get logical matrix for v that is hitting the boundary
        v = ~glt.*v + glt.*zeros(popsize,DIM); % reset v of particals that hit the boundary
        
        % Update pbest and gbest if necessary
        cost_x = feval(FUN, x');
        lt = cost_x < cost_p;
        cost_p = ~lt.*cost_p + lt.*cost_x; % replace lower cost with orignal personal best cost
        lt = repmat(lt',1,DIM); % repmat to the size of pbest where there will be DIM coordinates for each position
        pbest = ~lt.*pbest + lt.*x; % replace new pbest that is lower in cost
        [~,index] = min(cost_p);
        gbest = pbest(index,:); % update gbest

        % Exit if target is reached
        if feval(FUN, 'fbest') < ftarget
            break;
        end
    end
end 


  
