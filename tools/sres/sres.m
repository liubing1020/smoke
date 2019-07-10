function [xb,Statistics,Gm] = SRES_CONTINUOUS(fcn,mm,lu,lambda,G,mu,pf,varphi),
% SRES Evolution Strategy using Stochastic Ranking
% usage:
%        [xb,Stats,Gm] = sres(fcn,mm,lu,lambda,G,mu,pf,varphi) ;
% where
%        fcn       : name of function to be optimized (string)
%        mm        : 'max' or 'min' (for maximization or minimization)
%        lu        : parameteric constraints (lower and upper bounds)
%        lambda    : population size (number of offspring) (100 to 200)
%        G         : maximum number of generations
%        mu        : parent number (mu/lambda usually 1/7)
%        pf        : pressure on fitness in [0 0.5] try around 0.45
%        varphi    : expected rate of convergence (usually 1)
%
%        xb        : best feasible individual found
%        Stats     : [min(f(x)) mean(f(x)) number_feasible(x)]
%        Gm        : the generation number when "xb" was found

% Copyright (C) 1998-1999 Thomas Philip Runarsson (e-mail: tpr@verk.hi.is)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.

% randomize seed
  rand('seed',sum(100*clock)) ;
  if strcmp(lower(mm),'max'), mm = -1 ; else, mm = 1 ; end

% Initialize Population
  n = size(lu,2) ;
  x = ones(lambda,1)*lu(1,:)+rand(lambda,n).*(ones(lambda,1)*(lu(2,:)-lu(1,:))) ;

% Selection index vector
  sI = (1:mu)'*ones(1,ceil(lambda/mu)) ; sI = sI(1:lambda) ;

% Initial parameter settings
  eta = ones(lambda,1)*(lu(2,:)-lu(1,:))/sqrt(n) ;
  tau  = varphi/(sqrt(2*sqrt(n))) ;
  tau_ = varphi/(sqrt(2*n)) ;
  ub = ones(lambda,1)*lu(2,:) ;
  lb = ones(lambda,1)*lu(1,:) ;
  eta_u = eta(1,:) ;
  BestMin = inf ;
  nretry = 10 ;
  xb = [] ;

% Start Generation loop ...
  for g=1:G,
	print = ['generation:',num2str(g)];
	disp(print); 
  % fitness evaluation
    %[f,phi] = feval(fcn,x) ; f = mm*f ;
    [f,phi] = fcn(x) ; f = mm*f ;
    Feasible = find((sum((phi>0),2)<=0)) ;

  % Performance / statistics
    if ~isempty(Feasible),
      [Min(g),MinInd] = min(f(Feasible)) ;
      MinInd = Feasible(MinInd) ;
      Mean(g) = mean(f(Feasible)) ;
    else,
      Min(g) = NaN ; Mean(g) = NaN ;
    end
    NrFeas(g) = length(Feasible) ;

  % Keep best individual found
    if (Min(g)<BestMin) & ~isempty(Feasible)
      xb = x(MinInd,:) ;
      BestMin = Min(g) ;
      Gm = g ;
    end

  % Compute penalty function "quadratic loss function" (or any other)
    phi(find(phi<=0)) = 0 ;
    phi = sum(phi.^2,2) ;

  % Selection using stochastic ranking (see srsort.c)
    I = srsort(f,phi,pf) ;
    x = x(I(sI),:) ; eta = eta(I(sI),:) ;

  % Update eta (traditional technique with global intermediate recombination)
    eta = arithx(eta) ;
    eta = eta.*exp(tau_*randn(lambda,1)*ones(1,n)+tau*randn(lambda,n)) ;
    
  % Upper bound on eta (used?)
    for i=1:n,   
      I = find(eta(:,i)>eta_u(i)) ; 
      eta(I,i) = eta_u(i)*ones(size(I)) ;
    end

  % Mutation
    x_ = x ; % make a copy of the individuals for repeat ...
    x = x + eta.*randn(lambda,n) ;

  % If variables are out of bounds retry "nretry" times 
    I = find((x>ub) | (x<lb)) ;
    retry = 1 ;
    while ~isempty(I)
      x(I) = x_(I) + eta(I).*randn(length(I),1) ;
      I = find((x>ub) | (x<lb)) ;
      if (retry>nretry), break ; end
      retry = retry + 1 ;
    end
    % ignore failures
    if ~isempty(I),
      x(I) = x_(I) ;
    end
  end

% Check Output
  if isempty(xb),
    [dummy,MinInd] = min(phi) ;
    xb = x(MinInd,:) ;
    Gm = g ;
    disp('warning: solution is infeasible') ;
  end
  if nargout > 1,
    Statistics = [mm*[Min' Mean'] NrFeas'] ;
  end
