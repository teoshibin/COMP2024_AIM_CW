function[] = experimentJADEAE_fn(setExperiment, setDE, setAE, setJA)
% runs an entire experiment for benchmarking MY_OPTIMIZER
% on the noise-free testbed. fgeneric.m and benchmarks.m
% must be in the path of Matlab/Octave
% CAPITALIZATION indicates code adaptations to be made

%clear all

%addpath('PUT_PATH_TO_BBOB/matlab');  % should point to fgeneric.m etc.
datapath = setExperiment.datapath;                        % DIFFERENT FOLDER FOR EACH EXPERIMENT
opt.algName = setExperiment.name;                     % ALGORIHTM NAME
opt.comments = setExperiment.comment;
maxfunevals = setExperiment.maxfunevals;  % 10*dim is a short test-experiment taking a few minutes 
                                          % INCREMENT maxfunevals successively to larger value(s)
minfunevals = 'dim + 2';  % PUT MINIMAL SENSIBLE NUMBER OF EVALUATIONS for a restart
maxrestarts = 1e4;        % SET to zero for an entirely deterministic algorithm

if strcmp(setExperiment.useFunctions,'all')
   funList = eval(setExperiment.benchmark_fun);
else 
   funList = setExperiment.useFunctions;
end

% init log file
mkdir(datapath);
filename_1 = fullfile(datapath, 'log.txt');
file_1 = fopen(filename_1,'w+');
fclose(file_1);

t0 = clock;
rng('default');

for dim = setExperiment.dim  % small dimensions first, for CPU reasons
   for ifun = funList;  % or benchmarksnoisy(...)
      for iinstance = setExperiment.instance  % first 15 function instances
         fgeneric('initialize', ifun, iinstance, datapath, opt); 
      
      % independent restarts until maxfunevals or ftarget is reached
         for restarts = 0:maxrestarts
            if restarts > 0  % write additional restarted info
               fgeneric('restart', exitcode)
            end
            exitcode = JADEAEHansen_fn('fgeneric', {dim, setDE, setAE, setJA}, fgeneric('ftarget'), ...              % HERE CHANGE ALGORITHM
                        eval(maxfunevals) - fgeneric('evaluations'));
            if fgeneric('fbest') < fgeneric('ftarget') || ...
               fgeneric('evaluations') + eval(minfunevals) > eval(maxfunevals)
               break;
            end  
         end

         disp(sprintf(['  f%d in %d-D, instance %d: FEs=%d with %d restarts,' ...
                    ' fbest-ftarget=%.4e, elapsed time [h]: %.2f'], ...
                   ifun, dim, iinstance, ...
                   fgeneric('evaluations'), ...
                   restarts, ...
                   fgeneric('fbest') - fgeneric('ftarget'), ...
                   etime(clock, t0)/60/60));

         
         file_1 = fopen(filename_1,'a');
         fprintf(file_1,['f%d in %d-D, instance %d: FEs=%d with %d restarts,' ...
                    ' fbest-ftarget=%.4e, elapsed time [h]: %.2f\n'], ...
                   ifun, dim, iinstance, ...
                   fgeneric('evaluations'), ...
                   restarts, ...
                   fgeneric('fbest') - fgeneric('ftarget'), ...
                   etime(clock, t0)/60/60);
         fclose(file_1);
         fgeneric('finalize');
      end
      disp(['      date and time: ' num2str(clock, ' %.0f')]);
   end
	disp(sprintf('---- dimension %d-D done ----', dim));
end
