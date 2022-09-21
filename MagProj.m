%                     MagProj.m
%             written on July 24, 2001 by 
%                   Russell Luke
%   Inst. Fuer Numerische und Angewandte Mathematik
%                Universitaet Gottingen
%
% DESCRIPTION:  Projection operator onto a magnitude constraint
%
% INPUT:        constr = a nonnegative array that is the magnitude
%                        constraint
%               u = the function to be projected onto constr 
%                   (can be complex)
%
% OUTPUT:       unew = the projection
%
% USAGE: unew = MagProj(constr,u)
%
% Nonstandard Matlab function calls:  

function unew = MagProj(constr,u)

EPSILON = 3e-15;

% Fancy way of computing the projections.  NOT EXACT:
% it's really the gradient of a perturbed squared distance
% function to the constraint set as described in Luke-Burke-Lyon,
% SIREV 2002.

modsq_u = conj(u).*u;
denom  = (modsq_u+EPSILON).^(.5);
denom2 = modsq_u+EPSILON;
r_eps  = (modsq_u./denom) - constr;
dr_eps = (denom2 +EPSILON)./(denom2.*denom);
unew   = (1- dr_eps.*r_eps).*u;







