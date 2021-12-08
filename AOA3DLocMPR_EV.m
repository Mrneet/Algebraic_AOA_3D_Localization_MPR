function [ mprSol ] = AOA3DLocMPR_EV( theta, phi, senPos, Qa, Qs )
% [ mprSol ] = AOA3DLocMPR_EV( senPos, Qa, Qs, theta, phi )
%
% Estimate the source location in MPR by the Eigenvector (EV) method
%
% Input:
%   theta:  (M x 1), noisy azimuth angle measurements
%   phi:    (M x 1), noisy elevation angle measurements
%   senPos: (3 x M), positions of sensors, each column is a sensor
%           position(3D) and first column is for the reference sensor;
%   Qa:     (2M x 2M), AOA covariance matrix
%   Qs:     (3M x 3M), sensor position covariance matrix
%
% Output:
%   mprSol: (3 x 1), EV solution in MPR
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

iterNum = 2;    % number of iterations to improve the weighting matrix

[N,M] = size(senPos);

A1 = [sin(theta),-cos(theta),zeros(M,1)];
a1 = -diag(A1*(senPos));
A2 = [cos(theta).*sin(phi),sin(theta).*sin(phi),-cos(phi)];
a2 = -diag(A2*(senPos));
A = [A1;A2];
a = [a1;a2];

C1 = zeros(M,M*N); C2 = zeros(M,M*N);
W = inv(Qa);

for it = 1:iterNum
    O = W-W*a*a'*W/(a'*W*a);
    [V,D] = eig(A'*O*A);
    [~,IX] = sort(sum(D));
    
    u_sol = zeros(N,2);g_sol = zeros(1,2);
    J = zeros(1,2);

    u_sol(:,1) = V(:,IX(1));
    g_sol(1) = -a'*W*A*u_sol(:,1)/(a'*W*a);
    thetaM = zeros(M,1); phiM = zeros(M,1);
    
    for j = 1:M
        thetaM(j,1) = atan2(u_sol(2,1)-g_sol(1)*senPos(2,j), u_sol(1,1)-g_sol(1)*senPos(1,j));
        phiM(j,1) = atan2(u_sol(3,1)-g_sol(1)*senPos(3,j), norm(u_sol(1:2,1)-g_sol(1)*senPos(1:2,j),2));
    end
    J(1) = sum((theta-thetaM).^2+(phi-phiM).^2);

    u_sol(:,2) = -V(:,IX(1));
    g_sol(2) = -a'*W*A*u_sol(:,2)/(a'*W*a);
    
    for j = 1:M
        thetaM(j,1) = atan2(u_sol(2,2)-g_sol(2)*senPos(2,j), u_sol(1,2)-g_sol(2)*senPos(1,j));
        phiM(j,1) = atan2(u_sol(3,2)-g_sol(2)*senPos(3,j), norm(u_sol(1:2,2)-g_sol(2)*senPos(1:2,j),2));
    end
    J(2) = sum((theta-thetaM).^2+(phi-phiM).^2);

    [~,ind] = min(J);
    u_bar = u_sol(:,ind);
    g_est = -a'*W*A*u_bar/(a'*W*a);
    
    % update W
    b1 = sqrt(sum((u_bar(1:2)-g_est*(senPos(1:2,:))).^2,1))'/norm(u_bar,2);
    b2 = sqrt(sum((u_bar-g_est*(senPos)).^2,1))'/norm(u_bar,2);
    B1 = diag((b1));
    B2 = diag((b2));
    Ba = blkdiag(B1,B2);
    
    for m = 1:M
        thetaTmp = atan2(u_bar(2)-g_est*senPos(2,m), u_bar(1)-g_est*senPos(1,m));
        phiTmp = atan2(u_bar(3)-g_est*senPos(3,m), norm(u_bar(1:2)-g_est*senPos(1:2,m),2));
        alpha1 = [sin(thetaTmp);-cos(thetaTmp);0];
        C1(m,(1:N)+(m-1)*N) = -alpha1'*g_est;
        alpha2 = [cos(thetaTmp)*sin(phiTmp);sin(thetaTmp)*sin(phiTmp);-cos(phiTmp)];
        C2(m,(1:N)+(m-1)*N) = -alpha2'*g_est;
    end
    Ca = [C1;C2];
    
    W = inv(Ba*Qa*Ba' + Ca*Qs*Ca');
end

theta = atan2(u_bar(2),u_bar(1));
phi = atan2(u_bar(3),sqrt(u_bar(1)^2+u_bar(2)^2));
mprSol = [theta;phi;g_est];

