function [ Cov ] = AOA3DLocMPR_CovBR( srcLoc, senPos, Qa, Qs )
% [ Cov ] = AOA3DLocMPR_CovBR( srcLoc, senPos, Qa, Qs )
%
% Evaluate the theoretical covariance matrix of the MPR source location 
% estimate by the BR method using AOA measurements
%
% Input:
%   srcLoc: (3 x 1), source location in Cartesian (localization dimension = 3)
%   senPos: (3 x M), sensor location (M = number of sensors)
%   Qa:     (2M x 2M), AOA covariance matrix
%   Qs:     (3M x 3M), sensor position covariance matrix
%
% Output:
%   Cov: (3 x 3), covariance matrix of source location estimate in MPR by BR
%
% Reference:
% Y. Sun, K. C. Ho, and Q. Wan, "Eigenspace solution for AOA localization
% in modified polar representation," IEEE Trans. Signal Process.,
% vol. 68, pp. 2256-2271, 2020.
%
% Yimao Sun, K. C. Ho   03-28-2021
%
%       Copyright (C) 2020
%       Computational Intelligence Signal Processing Laboratory
%       University of Missouri
%       Columbia, MO 65211, USA
%       hod@missouri.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,M] = size(senPos);

go = 1/norm(srcLoc);
u0 = srcLoc/norm(srcLoc);
theta = atan2(u0(2),u0(1));
phi = atan2(u0(3),sqrt(u0(1)^2+u0(2)^2));

thetaM = atan2(srcLoc(2)-senPos(2,:),srcLoc(1)-senPos(1,:))';
phiM = atan2(srcLoc(3)-senPos(3,:),sqrt(sum((srcLoc(1:2)-senPos(1:2,:)).^2,1)))';
psi = [u0;go];
S = diag([ones(N,1);0]);

GA1 = [sin(thetaM),-cos(thetaM),zeros(M,1)];
Ga1 = -diag(GA1*(senPos));
GA2 = [cos(thetaM).*sin(phiM),sin(thetaM).*sin(phiM),-cos(phiM)];
Ga2 = -diag(GA2*(senPos));
G1 = [GA1,Ga1;GA2,Ga2];

b1 = sqrt(sum((u0(1:2)-go*(senPos(1:2,:))).^2,1))';
b2 = sqrt(sum((u0-go*(senPos)).^2,1))';
Ba = diag([b1;b2]);

for m = 1:M
    thetaTmp = atan2(u0(2)-go*senPos(2,m), u0(1)-go*senPos(1,m));
    phiTmp = atan2(u0(3)-go*senPos(3,m), norm(u0(1:2)-go*senPos(1:2,m),2));
    alpha1 = [sin(thetaTmp);-cos(thetaTmp);0];
    C1(m,(1:N)+(m-1)*N) = -alpha1'*go;
    alpha2 = [cos(thetaTmp)*sin(phiTmp);sin(thetaTmp)*sin(phiTmp);-cos(phiTmp)];
    C2(m,(1:N)+(m-1)*N) = -alpha2'*go;
end
Ca = [C1;C2];

W1 = inv(Ba*Qa*Ba' + Ca*Qs*Ca');
Q = blkdiag(Qa,Qs);

H1 = [cos(thetaM),sin(thetaM),zeros(M,1)];
h1 = -diag(H1*senPos);
a1 = -reshape(GA1',[],1);
A1 = [H1,h1;zeros(M,N),zeros(M,1);zeros(M*N,N),a1];

H21 = [-sin(thetaM).*sin(phiM), cos(thetaM).*sin(phiM), zeros(M,1)];
H22 = [cos(thetaM).*cos(phiM), sin(thetaM).*cos(phiM), sin(phiM)];
h21 = -diag(H21*senPos);
h22 = -diag(H22*senPos);
a2 = -reshape(GA2',[],1);

A2 = [H21,          h21;
      H22,          h22;
      zeros(M*N,N), a2];

A = [A1;A2];

W11 = W1(1:M,1:M);
W12 = W1(1:M,M+1:2*M);
W21 = W1(M+1:2*M,1:M);
W22 = W1(M+1:2*M,M+1:2*M);
k1 = N-1;                                
k2 = 2*N-1;                              
QW = kron(ones(k1,k1),Q).* ...
    [kron(ones(k2,k2),W11), kron(ones(k2,k2),W12);
     kron(ones(k2,k2),W21),  kron(ones(k2,k2),W22)];
Omega = A'*QW*A;

J = G1'*W1*G1;

kappa = Omega*psi;
U = null(kappa');
P = eye(N+1) - psi*psi'*S;
C = P*U/(U'*J*U)*U'*P';


D1 = [-sin(theta)/cos(phi),  cos(theta)/cos(phi),    0,          0;
      -cos(theta)*sin(phi), -sin(theta)*sin(phi),    cos(phi),   0;
      0,                     0,                      0,          1];
  
Cov = D1*C*D1';

end

