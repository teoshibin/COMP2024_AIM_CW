function ae = AEupdate(ae, pop, ccovfac, iterdiag)
% ae = AEupdate(ae, pop, ccovfac=1, iterdiag=0)
% 
% AEupdate (Adaptive Encoding) implements a CMA-like update for a
% linear encoding that is, in principle, applicable to any
% continuous domain iterative (search) algorithm.
% 
% INPUT: 
%   AE: an "adaptive encoding object". If AE is empty, AEupdate 
%      returns a new adaptive encoding "object". ae.B contains
%      the linear "(de-)coding", i.e. the fitness functions  
%      should be called like fitness(ae.B * x). ae.invB is the 
%      respective inverse. 
%   POP: DIM by POPSIZE array of most recent best
%      solutions. DIM is the search space dimensionality, 
%      the number of solutions POPSIZE might well be one. 
%      Only for non-uniform weights (see code below) an ordering
%      of the solutions regarding their fitness is required. 
%   CCOVFAC: factor for learning rate, default is one (about
%      ten times slower than default for CMA-ES). CCOVFAC depends
%      on the algorithm to which AEupdate is applied to and needs 
%      to be identified in advance, see reference. 
%   ITERDIAG: number of initial iterations, where the transformation 
%      is diagonal and therefore no rotation is applied. 
%      iterdiag=-1 invokes 200*DIM/POPSIZE. 
% OUTPUT: new or updated "adaptive encoding object". 
% 
% Example: 
%    x = rand(10,100);        % initial N by pop-size array x
%    opt = myoptim([], x);    % initialize optimizer, to be implemented
%    ae = AEupdate([], x);    % initialize encoding
%    for iter = 1:1e3
%       x = myoptim(opt, x, ae.B); % optimizers iteration step on x
%                                  % using fitness(ae.B * x)
%       x = ae.B * x;         % decode x
%       ae = AEupdate(ae, x); % adapt encoding 
%       x = ae.invB * x;      % encoding has changed, this could be 
%                             % the first row in the loop 
%    end
%    return ae.B * x(:,ibest) % return best solution, "decoded"
%
% REFERENCE: Hansen, N. (2008). Adaptive Encoding: How to Render
%    Search Coordinate System Invariant. In Rudolph et al. (eds.)
%    Parallel Problem Solving from Nature, PPSN X,
%    Proceedings, Springer. http://hal.inria.fr/inria-00287351/en/

version = '0.13'; 
last_change = '09/17/2009'; 

% TODO: * initial B optionally with random matrix? DONE (V0.13)
%       * find a better mapping of B between iterations? 
%         not sorting by EV but match scalar products 
%         column-wise? 

if nargin < 2 || isempty(pop)
  error('need two arguments, first can be empty');
end
N = size(pop, 1);

if nargin < 3 || isempty(ccovfac)
  ccovfac = 1;
end

if nargin < 4 || isempty(iterdiag)
  iterdiag = 0; 
end

if iterdiag < 0 
  iterdiag = 200*N/size(pop, 2);
end

  % initialize "object" 
  if isempty(ae) 
    % parameter setting
    ae.N = N; 
    ae.mu = size(pop, 2);
    ae.weights = ones(ae.mu, 1) / ae.mu; 
    ae.mucov = ae.mu; % for computing c1 and cmu
%     if 11 < 3  % non-uniform weights, assumes a correct ordering of
               % input arguments
      ae.weights = log(ae.mu+1)-log(1:ae.mu)'; 
      ae.weights = ae.weights/sum(ae.weights); 
      ae.mucov = 1/sum(ae.weights.^2); 
%     end
    ae.alpha_p = 1;                    % 
    ae.c1 = ccovfac*0.2/(ae.mucov+(N+1.3).^2);
    ae.cmu = ccovfac*0.2*(ae.mucov-2+1/ae.mucov)/(0.2*ae.mucov+(N+2).^2);
    ae.cc = 1/sqrt(N); 
    ae.diagonaliterations = iterdiag; 

    % initialization
    ae.countiter = 0; 
    ae.pc = zeros(N,1); 
    ae.xmean = pop * ae.weights;
    ae.C = eye(N); 
    ae.diagD = ones(N,1); 
    ae.B = eye(N);
    % generate if desired an arbitrary orthogonal initial transformation
    if 11 < 2  
      for i = 1:N
        v = randn(N,1);
        for j = 1:i-1
          v = v - (v'*ae.B(:,j)) * ae.B(:,j);
        end
        ae.B(:,i) = v / norm(v);  % might fail with a very small probability
      end
    end
    ae.Bo = ae.B;     % assuming the initial B is orthogonal
    ae.invB = ae.B';  % assuming the initial B is orthogonal
    return
  end % initialize object
  
  % check consistency
  if N ~= ae.N
    error(['dimensionality changed from ' num2str(ae.N) ' to ' ...
           num2str(N)]);
  end

  % begin 
  ae.countiter = ae.countiter + 1; 
  ae.xold = ae.xmean; 
  ae.xmean = pop * ae.weights;

  arnorm = sqrt(sum((ae.invB*(pop - repmat(ae.xold,1,ae.mu))).^2,1));

  alpha0 = sqrt(N) / sqrt(sum((ae.invB*(ae.xmean-ae.xold)).^2));
  if any(~isreal(alpha0)),
      error('imaginary');
  end

  alphai = sqrt(N) ./ arnorm; 
  alphai = sqrt(N) * min(1./median(arnorm), 2./arnorm); 
  if ~isfinite(alpha0)
    alpha0 = 1;
  end
  alphai(~isfinite(alphai)) = 1;

  % adapt the encoding
  z = alpha0 * (ae.xmean-ae.xold); 

  % for cumulation of z, no use if z is normalized anyway 
  if 11 < 3 && norm(ae.invB*z) < 0.1*norm(ae.invB*ae.pc)
    % in case of fast step size decrease
    % so far did rarely show up! 
    disp(['z is small' num2str([ norm(ae.invB*z) norm(ae.invB*ae.pc)])]);
    z = z / norm(ae.invB*z) * 0.1*norm(ae.invB*ae.pc);
  end
  ae.pc = (1-ae.cc)*ae.pc + sqrt(ae.cc*(2-ae.cc)) * z;
  if any(~isreal(ae.pc)),
      error('imaginary');
  end

  % S = z * z';
  S = ae.pc * ae.pc'; % 10-D: factor three on fcigar 

  zmu =  repmat(alphai, N,1) .* (pop - repmat(ae.xold, 1, ae.mu));
  Cmu = zmu * diag(ae.weights) * zmu';  
  % trace normalization seemly puts too much emphasis in large EVs, the 
  % small directions die out
  % Mmu = trace(ae.C) / trace(Mmu) * Mmu; 

     % norm_fac = N*norm(invB*ae.pc)^-2; disp(norm_fac)
     %        ae.C = (1-ae.ccov) * ae.C + ae.ccov/ae.mucov * ae.alpha_p * S ...
     %            + ae.ccov * (1-1/ae.mucov) * Cmu;

  ae.C = (1-(ae.c1+ae.cmu)) * ae.C + ae.c1 * ae.alpha_p * S ...
       + ae.cmu * Cmu;
  if ae.countiter <= ae.diagonaliterations
    ae.C = diag(diag(ae.C)); 
    ae.Bo = eye(N); 
    EV = ae.C; % like second output arg of eig() 
  else
    ae.C = (triu(ae.C)+triu(ae.C,1)');
    [ae.Bo, EV] = eig(ae.C);
            % limit condition of C to 1e14 + 1
            if max(diag(EV)) > 1e14*min(diag(EV))
                tmp = max(diag(EV))/1e14 - min(diag(EV));
                ae.C = ae.C + tmp*eye(N); EV = EV + tmp*eye(N);
            end
    if any(~isreal(ae.Bo(:))) || any(EV(:)<0),
      error('imaginary');
    end
  end
  [EV idx] = sort(diag(EV)); % vector of eigenvalues
  ae.diagD = sqrt(EV); 
  if any(~isreal(ae.diagD(:))),
      error('imaginary');
  end
  ae.Bo = ae.Bo(:,idx); 

        if 11 < 3 % just a check: fails to map Bo to corresponding C from above
          idx = randperm(N);
          EV = EV(idx);
          ae.Bo = ae.Bo(:,idx); 
        end
        % B = B*sqrt(EV)*B'; % worse on fcigar as below, because no
                             % additional scaling can be exploited

  ae.B = ae.Bo * diag(ae.diagD);
  ae.invB = diag(1./ae.diagD) * ae.Bo'; 
  if any(~isreal(ae.invB(:))),
      error('imaginary');
  end
  if 11 < 3  % orthogonal version versus full B 
    ae.B = ae.Bo;  
    ae.invB = ae.Bo';
  end

% ---------------------------------------------------------------
% Adaptive encoding. To be used under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).  
% Author (copyright): Nikolaus Hansen, 2008.  
% e-mail: nikolaus.hansen AT inria.fr 
% URL:http://www.lri.fr/~hansen
% ---------------------------------------------------------------

