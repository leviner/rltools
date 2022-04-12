function Pool =  Start_parallel(varargin)
% -----------------------------------------------
%  Start_parallel([Number_of_workers]) 
%
%   Version 3
%
%
%   If Number_of_workers is 0 or missing
%
%   Will default to max_workers below which may not be maximum
%
% ------------------------------------------------


Nproc = varargin{1};
if Nproc == 0
   
    fprintf('   No workers specified, no running in parallel...\n');
    delete(gcp('nocreate'));   % stop existing parpool, if running
    Pool = 0;
    return
    
end


max_workers = 56;

MyCluster = parcluster(parallel.defaultClusterProfile);
maxdesired_workers = MyCluster.NumWorkers;
 
if (exist('Pool','var') && Pool.NumWorkers == Nproc)
           fprintf('Parallel Pool set to %d workers',Pool.NumWorkers);
           return
end



Pool =  gcp('nocreate');
 
% No command line arguments        
if (nargin == 0)
    
    % if parallel not already started, start with default
    if (isempty(Pool) )
        try
          Pool = parpool(maxdesired_workers);
        catch 
            fprintf('Reset maximum workers to default...\n');
            Pool = parpool;
        end
    end
    
  
else                % Command line arguments exist
    
    Nproc = varargin{1};
    
    % if already set, check if Npoc is same as Pool
    if (exist('Pool.NumWorkers','var') && (Pool.NumWorkers == Nproc) )
           fprintf('Parallel Pool is set to %d workers\n',Pool.NumWorkers);
           return
    end

    % Check if within range, reset max number of workers
    if (Nproc > max_workers || Nproc < 1)
        Nproc = max_workers;
        fprintf('   Reset to maximum number of workers...\n');
    end
  
    % arguments included  
    if (isempty(Pool))                  % first time
        try
          Pool = parpool(Nproc);
        catch        
           fprintf('\t Resetting maximum workers to default...\n');
            Pool = parpool;
        end
        
    else                                % already started
        
        if Pool.NumWorkers ~= Nproc
            delete (Pool)
            try
                 Pool = parpool(Nproc);
            catch
                fprintf('\t Resetting maximum workers to default...\n');
                Pool = parpool;
            end
        end
        
    end
      
end
