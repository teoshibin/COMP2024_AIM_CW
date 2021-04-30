
addpath(genpath(pwd));
datasets_name = ["Bin_GA" "PSO" "PSO_Bounds" "PSO_EDA" "DE" "JADE" "DEAE"];

empty_opt.algName = "";
empty_opt.comments = "";
opt = repmat(empty_opt, 1, length(datasets_name));

opt(1).algName = 'Bin_GA';
opt(2).algName = 'PSO';
opt(3).algName = 'PSO_Bounds';
opt(4).algName = 'PSO_EDA';
opt(5).algName = 'DE';
opt(6).algName = 'JADE';
opt(7).algName = 'DEAE CCOVFAC study';

opt(1).comments = 'Binary Genetic Algorithm using Tournament Selection and Double Crossoover Mutation'; % Bin_GA
opt(2).comments = 'Partical Swarm Optimization'; % PSO
opt(3).comments = 'Population-Based Incremental Learning PSO'; % PSO_Bounds
opt(4).comments = 'PSO with Estimation of Distribution Evolutionary Algorithm'; % PSO_EDA
opt(5).comments = 'Differential Evolution'; % DE
opt(6).comments = 'Adaptive Differential Evolution with current-to-pbest mutation strategy'; % JADE 
opt(7).comments = 'DE pop5DIM F~U(0.5,1) CR=0.5 AE CCOVFAC 8'; % DEAE

maxfunevals = '1e3 * dim'; % INCREMENT maxfunevals successively to larger value(s)
minfunevals = 'dim + 2';  % PUT MINIMAL SENSIBLE NUMBER OF EVALUATIONS for a restart
maxrestarts = 1e4;        % SET to zero for an entirely deterministic algorithm

for ialgo = 1:length(datasets_name)
    
    fprintf("\n\n" + opt(ialgo).algName + "\n");
    
    datapath = convertStringsToChars('../Datasets/' + datasets_name(ialgo)); % path with the experiment data
    
    t0 = clock;
    rng('default');

    for dim = [2,3,5,10,20,40] % small dimensions first, for CPU reasons
        for ifun = benchmarks('FunctionIndices')  % or benchmarksnoisy(...)
            for iinstance = [1:15]  % first 15 function instances
                fgeneric('initialize', ifun, iinstance, datapath, opt(ialgo)); 

                % independent restarts until maxfunevals or ftarget is reached
                for restarts = 0:maxrestarts
                    
                    switch(ialgo)
                        case 1
                            % Running Bin_GA Algorithm
                            GA('fgeneric', dim, fgeneric('ftarget'), ...
                                eval(maxfunevals) - fgeneric('evaluations'));
                        case 2
                            % Running PSO Algorithm
                            PSO('fgeneric', dim, fgeneric('ftarget'), ...
                                eval(maxfunevals) - fgeneric('evaluations'));
                        case 3
                            % Running PSO_Bounds Algorithm
                            PSO_Bounds('fgeneric', dim, fgeneric('ftarget'), ...
                                eval(maxfunevals) - fgeneric('evaluations'));
                        case 4
                            % Running PSO_EDA Algorithm
                            PSO_EDA('fgeneric', dim, fgeneric('ftarget'), ...
                                eval(maxfunevals) - fgeneric('evaluations'));
                        case 5
                            % Running DE Algorithm
                            DE('fgeneric', dim, fgeneric('ftarget'), ...             
                                eval(maxfunevals) - fgeneric('evaluations'));
                        case 6
                            % Running JADE Algorithm
                            JADE('fgeneric', dim, fgeneric('ftarget'), ...             
                                eval(maxfunevals) - fgeneric('evaluations'));
                        case 7
                            % Running DEAE Algorithm
                            DEAE('fgeneric', dim, fgeneric('ftarget'), ...          
                                eval(maxfunevals) - fgeneric('evaluations'));
                    end

                    if fgeneric('fbest') < fgeneric('ftarget') || ...
                       fgeneric('evaluations') + eval(minfunevals) > eval(maxfunevals)
                       break;
                    end  
                end

                fprintf(['  f%d in %d-D, instance %d: FEs=%d with %d restarts,' ...
                        ' fbest-ftarget=%.4e, elapsed time [h]: %.2f\n'], ...
                       ifun, dim, iinstance, ...
                       fgeneric('evaluations'), ...
                       restarts, ...
                       fgeneric('fbest') - fgeneric('ftarget'), ...
                       etime(clock, t0)/60/60);

                fgeneric('finalize');
            end
            disp(['      date and time: ' num2str(clock, ' %.0f')]);
        end
        fprintf('---- dimension %d-D done ----\n', dim);
    end
end
