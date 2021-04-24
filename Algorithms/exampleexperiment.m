% runs an entire experiment for benchmarking MY_OPTIMIZER
% on the noise-free testbed. fgeneric.m and benchmarks.m
% must be in the path of Matlab/Octave
% CAPITALIZATION indicates code adaptations to be made

addpath(pwd);  % should point to fgeneric.m etc.
datapath = '../DEAE';  % different folder for each experiment
opt.algName = 'DEAE';
opt.comments = 'Differential Evolution with Adaptive Encoding Algorithm';
maxfunevals = '5e4'; % 10*dim is a short test-experiment taking a few minutes 
                          % INCREMENT maxfunevals successively to larger value(s)
minfunevals = 'dim + 2';  % PUT MINIMAL SENSIBLE NUMBER OF EVALUATIONS for a restart
maxrestarts = 1e4;        % SET to zero for an entirely deterministic algorithm

more off;  % in octave pagination is on by default

t0 = clock;
% rand('state', sum(100 * t0)); % THIS IS DEPRECIATED MATLAB CODE
rng(sum(100 * t0));

for dim = [2,3,5,10,20,40]  % small dimensions first, for CPU reasons
    for ifun = benchmarks('FunctionIndices')  % or benchmarksnoisy(...)
        for iinstance = [1:15]  % first 15 function instances
            fgeneric('initialize', ifun, iinstance, datapath, opt); 

            % independent restarts until maxfunevals or ftarget is reached
            for restarts = 0:maxrestarts
                % Running PSO Algorithm
                % PSO('fgeneric', dim, fgeneric('ftarget'), ...
                %             eval(maxfunevals) - fgeneric('evaluations'));
                
                % Running DE Algorithm
                % DE('fgeneric', dim, fgeneric('ftarget'), ...             
                %        eval(maxfunevals) - fgeneric('evaluations'));
                    
                % Running JADE Algorithm
                % JADE('fgeneric', dim, fgeneric('ftarget'), ...             
                %       eval(maxfunevals) - fgeneric('evaluations'));
                
                % Running DEAE Algorithm
                 DEAE('fgeneric', dim, fgeneric('ftarget'), ...             
                        eval(maxfunevals) - fgeneric('evaluations'));
                    
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

