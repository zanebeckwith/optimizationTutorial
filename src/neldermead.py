#!/usr/bin/env python

"""
This file implements a function to determine
the minimum of the input function 'func' using 
the Nelder-Mead (downhill simplex) method.

Originally written for the UNC-CH CAP-REU
program, summer of 2013.

Author: Zane Beckwith
Adapted from code provided by Dr. Fabian Heitsch.

input :
	initial_simplex	: list, with n+1 members, each member
		is a list of n values, where n is the dimension of func
	func						: the objective function. Must accept n-vector
	max_iterations	: maximum number of iterations used
	tol  						: tolerance (applied on f)
	iverb						: 0 = no diagnostics
			 						  1 = print out progress
			 						  2 = show map with each simplex
output :
	xmin 						: location of the minimum 
										(center of gravity of final simplex).
"""

def minimize(initial_simplex, max_iterations, tol, iverb):
  n    = size(initial_simplex)
  fcur = zeros(1,n(2))   % n+1 function values (current)
  fsor = zeros(1,n(2))   % sorted function values
  scur = initial_simplex         % initialize with first simplex
  ssor = zeros(n(1),n(2)) % sorted simplex positions.
  xmin = zeros(n(1),1)
  xref = zeros(n(1),1) % reflected point
  xexp = zeros(n(1),1) % expanded point
  xcon = zeros(n(1),1) % contracted point
  x0   = zeros(n(1),1)
  x1   = zeros(n(1),1)
  fref = 0.0
  fexp = 0.0
  fcon = 0.0
  dx   = 0.0  % keeps track of largest distance between simplex points
  dxs  = 0.0
  istol= 0
  done = 0    % keeps track of whether we're done within one iteration
  it   = 0
  jmax = 0

% these are the coefficients for 
  alpha= 1.0 % reflection
  gamma= 2.0 % expansion
  rho  = 0.5 % contraction
  sigma= 0.5 % reduction

% draw the starting simplex
  if (iverb == 2)
    color('white')
    plot([initial_simplex(1,1) initial_simplex(1,2)],[initial_simplex(2,1) initial_simplex(2,2)],'w')
    plot([initial_simplex(1,2) initial_simplex(1,3)],[initial_simplex(2,2) initial_simplex(2,3)],'w')
    plot([initial_simplex(1,3) initial_simplex(1,1)],[initial_simplex(2,3) initial_simplex(2,1)],'w')
  end
  pause(1.0)

  while ((it < max_iterations) && (istol == 0))
    done = 0 % start again
%   step 1: calculate f(x) for all vertices, and sort
    for j=1:n(2)
      fcur(j) = fxy(scur(:,j))
    end
    [fsor indx] = sort(fcur)  % fsor = fcur(indx) 
    ssor        = scur(:,indx)% re-index the positions 
%   step2: calculate the center of gravity x0 for all points except the worst x(n+1)
    sumf = sum(fsor(1:n(2)-1))
    for i=1:n(1)
      x0(i) = sum(fsor(1:n(2)-1).*ssor(i,1:n(2)-1))/sumf
    end 
    if (iverb > 0)
      fprintf(1,'minimize: it=%4d: centroid   : f = %13.5e\n',it+1,fxy(x0))
    end
%   step3: reflection
    xref = x0+alpha.*(x0-ssor(:,n(2)))
    fref = fxy(xref)
    if ((fsor(1) <= fref) && (fref < fsor(n(1))))
      if (iverb > 0)
        fprintf(1,'minimize: it=%4d: reflection : f = %13.5e\n',it+1,fref)
      end
      ssor(:,n(2)) = xref(:)
      done         = 1
    end
%   step4: expansion
    if (done == 0)
      if (fref < fsor(1))
        xexp = x0+gamma.*(x0-ssor(:,n(2)))
        fexp = fxy(xexp)
        if (iverb > 0)                                                                                                              
          fprintf(1,'minimize: it=%4d: expansion  : f = %13.5e\n',it+1,fexp)
        end                                                                                                                         
        if (fexp < fref)
          ssor(:,n(2)) = xexp(:)
          fsor(n(2))   = fexp
        else
          ssor(:,n(2)) = xref(:)
          fsor(n(2))   = fref
        end
        done      = 1
      end
    end
%   step5: contraction
    if (done == 0) 
      if (fref < fsor(n(2)-1))
        fprintf(1,'minimize: invalid branch in step 5\n')
        stop
      end
      xcon = ssor(:,n(2))+rho.*(x0-ssor(:,n(2)))
      fcon = fxy(xcon)
      if (fcon < fsor(n(2)))
        if (iverb > 0)                                                                                                                
          fprintf(1,'minimize: it=%4d: contraction: f = %13.5e\n',it+1,fcon)                                                                        
        end
        ssor(:,n(2)) = xcon(:)
        fsor(n(2))   = fcon
        done         = 1
      end
    end
%   step6: reduction
    if (done == 0) 
      if (iverb > 0)
        fprintf(1,'minimize: it=%4d: reduction\n',it+1)
      end
      for j=2:n(2)
        ssor(:,j) = ssor(:,1)+sigma.*(ssor(:,j)-ssor(:,1))
        fsor(j)   = fxy(ssor(:,j))
      end
    end
%  still need to compare current f etc (all steps).
    scur(:,:) = ssor(:,:) % need to write back the modified simplex as new input
    it = it+1
%  calculate largest distance between points: this is not recommended for real applications...
    dx = -1e10
    for i=1:n(1)
      for   j=i+1:n(2)
        dxs = sqrt(sum((ssor(:,i)-ssor(:,j)).^2))
        if (dxs > dx)
          dx   = dxs
          jmax = j
        end
      end
    end
    if (dx < tol) 
      istol = 1
      xmin  = ssor(:,jmax)
    end
    fprintf(1,'minimize: it=%4d: |dx| = %13.5e\n',it,dx)
    
%   for visualization: plot most recent simplex
    if (iverb == 2)
      color('white')
      plot([ssor(1,1) ssor(1,2)],[ssor(2,1) ssor(2,2)],'w') 
      plot([ssor(1,2) ssor(1,3)],[ssor(2,2) ssor(2,3)],'w')
      plot([ssor(1,3) ssor(1,1)],[ssor(2,3) ssor(2,1)],'w')
    end
    pause(1.0)
  end
  if (it == max_iterations)
    fprintf(1,'minimize: reached maximum number of iterations')
  end
end
