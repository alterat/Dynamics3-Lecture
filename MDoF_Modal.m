% === Mass-Spring array, Modal Approach

% This code solves the equations of motion for a lumped mass-spring array
% using the modal method. 
% It can serve as a (very!) crude first example of a physics-based
% sound synthesis algorithm. 
% The vector of amplitudes is chosen randomly, as well as the frequency of
% the first mode.
%
% If you want to experiment some aleatoric music, try running this code
% several times, like this:
%   for II=1:10
%   MDoF_Modal; pause(rand())
%   end
%
% INPUT PARAMETERS:
%   DT: duration of the simulation (in seconds)
%   SR: sampling rate (in Hz)
%   N:  number of degrees of freedom (masses)
%   pt: output location (from 1 to N)
%   M:  vector of the mass values (default random values)
%   C:  vector of the stiffness coefficients (default random values)
%   Ampli: vector of amplitudes for the modes (length N)
%   phi: vector of phases (length N)
%
% AUTHOR: Alberto Torin, University of Edinburgh
% Under Creative Commons By Attribution license.
% ==========================================================

% == Initial parameters
DT = 2;
SR = 96000;

N = 70;
pt = 23;

M = 0.5 + rand(N,1);
C = (2*pi*440)^2 * rand(N+1,1);

Ampli = 10*rand(N,1).^2;
phi = zeros(N,1);


% == CODE STARTS HERE 
% find number of samples and create time vector
NF = floor(DT*SR);
time = (0:NF-1)/SR;

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

% find eigenvectors and eigenvalues
% ATTENTION: this operation requires creating a full matrix. This can
% require a lot of memory when N is large!
[V,E] = eig(full(B));
evals = sqrt(diag(E));

% create a matrix of the argument for the sine function
T = kron(evals, time) + kron(phi, ones(numel(time),1)');

% find the relative amplitudes of the various modes at the specified point
Amps = V(pt, :)'.*Ampli;

% create the individual sound signals for the various modes
snd = sparse(diag(Amps)) * sin(T); 

% create the final signal by summing the individual ones
y = sum(snd,1)';

% == Plot the initial portion of the final signal
plot(y(1:2000));

% == fade out to avoid click
ramp = linspace(1,0,1000)';
fade = ones(NF,1);
fade(end-999:end) = ramp;

% == play the sound
soundsc(fade.*y, SR)
