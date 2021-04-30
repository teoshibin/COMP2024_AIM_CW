% Binary Genenatic Algorithm
% with tournement selection, two-point crossover & mutation with even
% probibility
%
% GA
%
% This code is coded by us with online referencing.
% Please find the online referenced source code citation as shown below.
%%%--------------------------------------------------------------------------------------------------%%%
%    Title: Binary and Real-Coded Genetic Algorithms
%    Date: 2015
%    Code version: 1
%    Availability: https://www.mathworks.com/matlabcentral/fileexchange/52856-binary-and-real-coded-genetic-algorithms
%%%--------------------------------------------------------------------------------------------------%%%

function GA(FUN, DIM, ftarget, maxfunevals)
    % GA
    % representation - a population with popsize of individual solutions,
    %                  an individual with DIM of chromosomes, a chromosome 
    %                  with 32 bit logical values that represents a float number
    % selection      - tournement
    % crossover      - two point
    % mutation       - even probability with bit length
    % termination    - fbest < ftarget
    
    % Initialise Parameters for Binary Encoding
    % 1 sign bit, 3 integer bit, 28 decimal bits
    BIT_SIZE = 32;
    INT_BIT = 3;
    DEC_BIT = BIT_SIZE - INT_BIT - 1;

    %% Problem Definition
    MAX_VALUE = 5;
    MIN_VALUE= -5;

    %% GA Parameters

    maxfunevals = min(1e8 * DIM, maxfunevals); 
    counteval = 0;

    nPop=floor(sqrt(maxfunevals)*4); 	% Population Size

    pc=1.0;                 % Crossover Percentage
    nc=2*round(pc*nPop/2);  % Number of Offsprings (also Parents)
    
    mu=2/(BIT_SIZE * DIM);  % Mutation Rate -- P_m in report
    
    TournamentSize = 3;
   

    %% Initialization
    empty_individual.Position=zeros(DIM,1);
    empty_individual.Binary=zeros(DIM, BIT_SIZE);
    empty_individual.Cost=[];

    pop=repmat(empty_individual,nPop,1);
 
    % Evaluating the fitness function of the initial population
    for i=1:nPop
        % Initialize Population
        pop(i).Position = MIN_VALUE + (MAX_VALUE-MIN_VALUE) * rand(DIM, 1);
        
        % Initialize Chromosomes
        pop(i).Binary = ga_encode(pop(i).Position, BIT_SIZE, INT_BIT, DEC_BIT);

        % Evaluation
        pop(i).Cost = feval(FUN, pop(i).Position);
    end
    counteval = counteval + nPop;

    % Sort Population
    Costs=[pop.Cost];
    [~, SortOrder]=sort(Costs);
    pop = pop(SortOrder);

    %% Main Loop
    while counteval < maxfunevals
     
        popc=repmat(empty_individual,nc/2,2);
        for k=1:nc/2
            % Selection
            i1=tournament_selection(pop,TournamentSize);
            i2=tournament_selection(pop,TournamentSize);
            

            % Select Parents
            p1=pop(i1);
            p2=pop(i2);

            % Perform Crossover
            [popc(k,1).Position, popc(k,2).Position] = crossover(p1.Position,p2.Position);
            
            % Perform Mutation
            popc(k,1).Binary = ga_encode(popc(k,1).Position, BIT_SIZE, INT_BIT, DEC_BIT);
            popc(k,1).Position = ga_decode(mutate(popc(k,1).Binary, mu), INT_BIT, DEC_BIT);
            
            popc(k,2).Binary = ga_encode(popc(k,2).Position, BIT_SIZE, INT_BIT, DEC_BIT);
            popc(k,2).Position = ga_decode(mutate(popc(k,2).Binary, mu), INT_BIT, DEC_BIT);
            
            % Evaluate Offsprings
            popc(k,1).Cost=feval(FUN,popc(k,1).Position);
            popc(k,2).Cost=feval(FUN,popc(k,2).Position);
        end
        popc=popc(:);
        counteval = counteval + nc;
        
        % Create Merged Population
        pop=[pop 
             popc]; %#ok

        % Sort Population
        Costs = [pop.Cost];
        [~, SortOrder] = sort(Costs);
        pop = pop(SortOrder);
        
        % Truncation
        pop=pop(1:nPop);
        Costs=[];

        % Show Iteration Information
        % fprintf("Iteration: %d --- Best: %f\n", it, fgeneric('fbest'));

        % Exit if target is reached
        if feval(FUN, 'fbest') < ftarget
            break;
        end
    end
end 



function x = ga_decode(binary_x, int_bit, dec_bit)
    x(:,1) = binary_x(:,2:end) * pow2([int_bit-1:-1:0 -(1:dec_bit)].');
    x(binary_x(:,1)==1) = x(binary_x(:,1)==1) * -1;
end

function binary_x = ga_encode(x, bit_size, int_bit, dec_bit)
    binary_x = zeros(numel(x), bit_size);
    binary_x(:,1) = x < 0;
    binary_x(:,2:end) = fix(abs(rem(x.*pow2(-(int_bit-1):dec_bit),2)));
end

function i = tournament_selection(pop,m)

    nPop=numel(pop);

    S=randsample(nPop,m);
    
    spop=pop(S);
    
    scosts=[spop.Cost];
    
    [~, j]=min(scosts);
    
    i=S(j);
end

function [y1, y2] = crossover(x1, x2)
    nVar = numel(x1);
    
    x1 = x1';
    x2 = x2';
    
    cc=randsample(nVar,2);
    c1=min(cc);
    c2=max(cc);
    
    y1=[x1(1:c1) x2(c1+1:c2) x1(c2+1:end)];
    y2=[x2(1:c1) x1(c1+1:c2) x2(c2+1:end)];
    
    y1 = y1';
    y2 = y2';
end

function y = mutate(x,mu)
    [n, m] = size(x);
    temp = rand(n,m);
    y = abs(x - (temp < mu));
end
  
