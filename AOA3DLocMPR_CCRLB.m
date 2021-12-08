function [ CRLB ] = AOA3DLocMPR_CCRLB( srcLoc, senPos, varargin )
% [ CRLB ] = AOA3DLocMPR_CCRLB( srcLoc,senPos,Qa,Qs )
%
% Evaluate the CRLB for AOA localization in MPR using the Constrained CRLB
% formulation
%
% Input:
%   srcLoc: (N x 1), source location (N = localization dimension)
%   senPos: (N x M), sensor location (M = number of sensors)
%   varargin:
%       Qa: (2M x 2M), AOA covariance matrix
%       Qs: (MN x MN), sensor position covariance matrix
%
% Output:
%   CRLB: CRLB matrix
%
% Reference:
% Y. Sun, K. C. Ho, and Q. Wan, "Eigenspace solution for AOA localization
% in modified polar representation," IEEE Trans. Signal Process.,
% vol. 68, pp. 2256-2271, 2020.
%
% Yimao Sun, K. C. Ho   06-28-2020
%
%       Copyright (C) 2020
%       Computational Intelligence Signal Processing Laboratory
%       University of Missouri
%       Columbia, MO 65211, USA
%       hod@missouri.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,M] = size(senPos);

thetaM = atan2(srcLoc(2)-senPos(2,:),srcLoc(1)-senPos(1,:))';
d = sqrt(sum((srcLoc(1:2)-senPos(1:2,:)).^2,1))';
r = sqrt(sum((srcLoc-senPos).^2,1))';
r0 = norm(srcLoc);

DrAng = zeros(2*M,N+1);
for i = 1:M
    DrAng(i,:) = [-r0/d(i)^2*(srcLoc(2)-senPos(2,i)), ...
        r0/d(i)^2*(srcLoc(1)-senPos(1,i)), ...
        0, ...
        r0/d(i)^2*(senPos(1,i)*(srcLoc(2))-srcLoc(1)*(senPos(2,i)))
        ];
    DrAng(i+M,:) = [-r0/r(i)^2*(srcLoc(3)-senPos(3,i))*cos(thetaM(i)), ...
        -r0/r(i)^2*(srcLoc(3)-senPos(3,i))*sin(thetaM(i)), ...
        r0*d(i)/r(i)^2, ...
        -r0/r(i)^2*(senPos(3,i))*d(i) + ...
        r0/r(i)^2*(srcLoc(3)-senPos(3,i)) * ((senPos(1,i))*cos(thetaM(i))+(senPos(2,i))*sin(thetaM(i)))
        ];
end

L = length(varargin);
if L == 1
    Q = varargin{1};
elseif L == 2
    phiM = atan2(srcLoc(3)-senPos(3,:),sqrt(sum((srcLoc(1:2)-senPos(1:2,:)).^2,1)))';
    Ba = diag([d;r])/r0;
    C1 = zeros(M,M*N); C2 = zeros(M,M*N);
    for m = 1:M
        alpha1 = [sin(thetaM(m));-cos(thetaM(m));0];
        C1(m,(1:N)+(m-1)*N) = -alpha1'/r0;
        alpha2 = [cos(thetaM(m))*sin(phiM(m));sin(thetaM(m))*sin(phiM(m));-cos(phiM(m))];
        C2(m,(1:N)+(m-1)*N) = -alpha2'/r0;
    end
    Ca = [C1;C2];
    Q = varargin{1} + Ba\Ca*varargin{2}*Ca'/Ba;
else
    error('Too many input.');
end

FIM = DrAng'/Q*DrAng;

psi = [srcLoc;1]/r0;
F = 2*psi'*diag([ones(1,N),0]);
U = null(F);

theta = atan2(srcLoc(2),srcLoc(1));
phi = atan2(srcLoc(3),norm(srcLoc(1:2),2));
D = [-sin(theta)/cos(phi), cos(theta)/cos(phi),    0,        0;
    -cos(theta)*sin(phi), -sin(theta)*sin(phi),   cos(phi),  0;
    0,                    0,                     0,          1];

CRLB = D*U/(U'*FIM*U)*U'*D';