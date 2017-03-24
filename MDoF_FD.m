% === Mass-Spring array, Finite Difference Approach

% This code solves the equations of motion for a lumped mass-spring array
% using the finite difference method. 
% It can serve as a (very!) crude first example of a physics-based
% sound synthesis algorithm. 
% The masses and the stiffness coefficient are chosen randomly, as well as
% the initial conditions (u0 and u1). 
% The sampling rate chosen should be high enough to avoid instability. A
% stability check is added to prevent the code from running with unsuitable
% parameters.
%
% If you want to experiment some aleatoric music, try running this code
% several times, like this:
%   for II=1:10
%   MDoF_FD, pause(rand())
%   end
%
% More information on the FD implementation can be found in:
% S. Bilbao, Numerical Sound Synthesis, Wiley, Ch. 3.4
%
% INPUT PARAMETERS:
%   DT: duration of the simulation (in seconds)
%   SR: sampling rate (in Hz)
%   N:  number of degrees of freedom (masses)
%   pt: output location (from 1 to N)
%   M:  vector of the mass values (default random values)
%   C:  vector of the stiffness coefficients (default random values)
%
% AUTHOR: Alberto Torin, University of Edinburgh
% Under Creative Commons By Attribution license.
% ==========================================================

% == Initial Parameters
DT = 2;
SR = 96000;

N = 60;
pt = 4;

M = 0.5+rand(N,1);
C = (2*pi*440)^2 * rand(N+1,1);


% == CODE STARTS HERE 
% convert vectors into diagonal matrices
M = sparse(diag(M));
C = sparse(diag(C));

% create matrix A and A transpose (see lecture slides)
A = spdiags([-ones(N,1), ones(N,1)], -1:0, N+1, N);
AT = A';

% create stiffness matrix
K = AT*C*A;

% divide by masses matrix
B = M\K;

% Check stability of the scheme
if sqrt(eigs(B, 1, 'lm'))/2 > SR
    fprintf('Stability condition violated\n')
    return
end

% Check output location
if pt > N
    fprintf('Output location should not be larger than N\n')
    return
end

% find time step and number of samples
t = 1/SR;
NF = floor(DT*SR);

% == Initial conditions and preallocate output
u0 = rand(N,1).^3;
u1 = rand(N,1).^3;

out = zeros(NF,1);

% == Main Loop: solve FD scheme
for nn = 1:NF
   u = -t^2*B*u1 + 2*u1 - u0;
   
   out(nn) = u(pt);
   
   u0 = u1; u1 = u;
    
end

% == Plot output sound
plot((0:NF-1)*t, out/max(abs(out)))
xlabel('time (s)')
title('Output sound')

% == fade out to avoid click
ramp = linspace(1,0,1000)';
fade = ones(NF,1);
fade(end-999:end) = ramp;

% == play the output
soundsc(out.*fade, SR)