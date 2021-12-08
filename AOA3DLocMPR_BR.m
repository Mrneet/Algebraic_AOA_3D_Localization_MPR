function [ mprSol ] = AOA3DLocMPR_BR( theta, phi, senPos, Qa, Qs )
% [ mprSol ] = AOA3DLocMPR_BR( theta, phi, senPos, Qa, Qs )
%
% Estimate the source location in MPR by the Bias Reduction (BR) method
%
% Input:
%   theta:  (M x 1), noisy azimuth angle measurements
%   phi:    (M x 1), noisy elevation angle measurements
%   senPos: (3 x M), noisy positions of sensors, each column is a sensor
%           position(3D) and first column is for the reference sensor;
%   Qa:     (2M x 2M), AOA covariance matrix
%   Qs:     (3M x 3M), sensor position covariance matrix
%
% Output:
%   mprSol:     (3 x 1), BR solution in MPR
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


iterNum = 2;

[N,M] = size(senPos);

GA1 = [sin(theta),-cos(theta),zeros(M,1)];
Ga1 = -diag(GA1*(senPos));
GA2 = [cos(theta).*sin(phi),sin(theta).*sin(phi),-cos(phi)];
Ga2 = -diag(GA2*(senPos));
G1 = [GA1,Ga1;GA2,Ga2];

H1 = [cos(theta),sin(theta),zeros(M,1)];
h1 = -diag(H1*senPos);
a1 = -reshape(GA1',[],1);
A1 = [H1,h1;zeros(M,N),zeros(M,1);zeros(M*N,N),a1];

H21 = [-sin(theta).*sin(phi), cos(theta).*sin(phi), zeros(M,1)];
H22 = [cos(theta).*cos(phi), sin(theta).*cos(phi), sin(phi)];
h21 = -diag(H21*senPos);
h22 = -diag(H22*senPos);
a2 = -reshape(GA2',[],1);

A2 = [H21,          h21;
      H22,          h22;
      zeros(M*N,N), a2];

A = [A1;A2];



C1 = zeros(M,M*N); C2 = zeros(M,M*N);
W1 = inv(Qa);
Q = blkdiag(Qa,Qs);

for it = 1:iterNum
    W11 = W1(1:M,1:M);
    W12 = W1(1:M,M+1:2*M);
    W21 = W1(M+1:2*M,1:M);
    W22 = W1(M+1:2*M,M+1:2*M);
    k1 = N-1;                                
    k2 = 2*N-1;                              
    QW = kron(ones(k1,k1),Q).* ...
        [kron(ones(k2,k2),W11), kron(ones(k2,k2),W12);
         kron(ones(k2,k2),W21),  kron(ones(k2,k2),W22)];
    R = G1'*W1*G1;
    Omega = A'*QW*A;
    
    [V,D2] = eig(R,Omega);
    [~,IX] = sort(sum(D2));
    u_sol = [V(1:N,IX(1)), -V(1:N,IX(1))]/norm(V(1:N,IX(1)),2);
    g_sol = [V(N+1,IX(1)), -V(N+1,IX(1))]/norm(V(1:N,IX(1)),2);

    J = zeros(1,2);
    for i = 1:2
        thetaM = zeros(M,1); phiM = zeros(M,1);
        for j = 1:M
            thetaM(j,1) = atan2(u_sol(2,i)-g_sol(i)*senPos(2,j), u_sol(1,i)-g_sol(i)*senPos(1,j));
            phiM(j,1) = atan2(u_sol(3,i)-g_sol(i)*senPos(3,j), norm(u_sol(1:2,i)-g_sol(i)*senPos(1:2,j),2));
        end
        J(i) = sum((theta-thetaM).^2+(phi-phiM).^2);
    end
    [~,ind] = min(J);
    u_bar = u_sol(:,ind);
    g_est = g_sol(ind);

    % update W
    b1 = sqrt(sum((u_bar(1:2)-g_est*(senPos(1:2,:))).^2,1))'/norm(u_bar,2);
    b2 = sqrt(sum((u_bar-g_est*(senPos)).^2,1))'/norm(u_bar,2);
    Ba = diag([b1;b2]);
    
    for m = 1:M
        thetaTmp = atan2(u_bar(2)-g_est*senPos(2,m), u_bar(1)-g_est*senPos(1,m));
        phiTmp = atan2(u_bar(3)-g_est*senPos(3,m), norm(u_bar(1:2)-g_est*senPos(1:2,m),2));
        alpha1 = [sin(thetaTmp);-cos(thetaTmp);0];
        C1(m,(1:N)+(m-1)*N) = -alpha1'*g_est;
        alpha2 = [cos(thetaTmp)*sin(phiTmp);sin(thetaTmp)*sin(phiTmp);-cos(phiTmp)];
        C2(m,(1:N)+(m-1)*N) = -alpha2'*g_est;
    end
    Ca = [C1;C2];
    
    W1 = inv(Ba*Qa*Ba' + Ca*Qs*Ca');
end

theta = atan2(u_bar(2),u_bar(1));
phi = atan2(u_bar(3),sqrt(u_bar(1)^2+u_bar(2)^2));
mprSol = [theta;phi;g_est];

