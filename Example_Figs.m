% This program gives an example of how to call the routines for AOA
% localization in MPR using the EV and BR methods that can reproduce the
% simulation figures in the paper
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

clear all; 
%clc;
warning off

% =========================================
%  im = 1: performance vs noise power (Figs. 3-6)
%  im = 2: perforamnce vs sensor position errors (Figs. 7-10)
%  im = 3: performance vs source range (Figs. 11-14)
% =========================================
im = 1;

% -- settings --
senPosTrue = [              % sensor positions
    0      0      0
    -13.1   22.3  -32.4
    34.6   44.4   15.5
    30.5  -29.4   18.4
    9.9   29.2   -8.4
    -7.4   48.0   -3.5
    4.5   -0.9   21.7
    -47.2   -7.3   29.1
    ]';
thetaSrc = 22.13*pi/180;    % source azimuth
phiSrc = 15.41*pi/180;      % source elevation

[N,M] = size(senPosTrue);

models = ['nse';'err';'rag'];

model = models(im,:);
switch model
    case 'nse'
        % ******* vs. noise power config, Fig. 3-6 *******
        nsePwr = -70:10:20;         % 10log(rad^2)
        srcRange = 200;             % m
        errLvl = -20;
    case 'err'
        % ******* vs. sensor position error level config, Fig. 7-10 *******
        nsePwr = -50;               % 10log(rad^2)
        srcRange = 200;             % m
        errLvl = -60:10:40;
    case 'rag'
        % ******* vs. range config, Fig. 11-14 *******
        nsePwr = -20;               % 10log(rad^2)
        srcRange = [20,50:50:800];  % m
        errLvl = -20;
end

T = 1000;     % number of ensemble runs

% AOA measurements noise
rng('Default');
nseTmp = zeros(2*M,T);
for m = 1:T
    nseTmp(:,m) = randn(2*M,1);
end
nseTmp = nseTmp - mean(nseTmp,2);

% sensor position errors
errTmp = zeros(M*N,1);
for m = 1:T
    errTmp(:,m) = randn(M*N,1);
end
errTmp = errTmp - mean(errTmp,2);

Ra=eye(2*M);
aErr = [1,80,80,50,5,3,20,10];
Rs = kron(eye(N),diag(aErr(1:M)));

rng(1);
for ii = 1:2
    a1 = rand(3*M,1);
    a2 = rand(M,1);
    Rt = diag(roundn(a2/mean(a2),-2));
    a3 = rand(M,1);
    Rp = diag(roundn(a3/mean(a3),-2));
    Ra = blkdiag(Rt,Rp);
end

totalTime = zeros(1,2);

% -- perform processing --
% over source range
for ir = 1:length(srcRange)
    gSrc = 1/srcRange(ir);
    u0 = [cos(thetaSrc)*cos(phiSrc); sin(thetaSrc)*cos(phiSrc); sin(phiSrc)];
    
    srcLoc = srcRange(ir)*u0;
    
    % true AOA measurements
    thetaTrue = atan2(srcLoc(2)-senPosTrue(2,:),srcLoc(1)-senPosTrue(1,:))';
    phiTrue = atan2(srcLoc(3)-senPosTrue(3,:),sqrt(sum((srcLoc(1:2)-senPosTrue(1:2,:)).^2,1)))';
    
    % over noise power
    for in = 1:length(nsePwr)
        Qa = 10^(nsePwr(in)/10)*Ra;
        
        % over sensor position error level
        for is = 1:length(errLvl)
            disp(['source range: ',num2str(srcRange(ir)),'m, 20log10(\sigma_a(rad)): ', num2str(nsePwr(in)),', 20log10(\sigma_s(m)): ',num2str(errLvl(is)),' ...']);
            
            Qs = 10^(errLvl(is)/10)*Rs;
            
            % CRLB
            CRB_MPR = AOA3DLocMPR_CCRLB( srcLoc,senPosTrue,Qa,Qs );
            crlb_t(ir,in,is) = CRB_MPR(1,1);
            crlb_p(ir,in,is) = CRB_MPR(2,2);
            crlb_g(ir,in,is) = CRB_MPR(3,3);
            
            % Theoretical covariance matrices
            Cov_EV  = AOA3DLocMPR_CovEV( srcLoc, senPosTrue, Qa, Qs );
            Cov_BR = AOA3DLocMPR_CovBR( srcLoc, senPosTrue, Qa, Qs );
            covEV_t(ir,in,is) = Cov_EV(1,1);
            covEV_p(ir,in,is) = Cov_EV(2,2);
            covEV_g(ir,in,is) = Cov_EV(3,3);
            covBR_t(ir,in,is) = Cov_BR(1,1);
            covBR_p(ir,in,is) = Cov_BR(2,2);
            covBR_g(ir,in,is) = Cov_BR(3,3);
            
            [pos11,pos12,pos13,pos14,pos15,pos16] = deal(zeros(N,T));
            
            % over Monte-Carlo runs
            for m = 1:T
                angNse = sqrtm(Qa)*nseTmp(:,m);
                theta = thetaTrue + angNse(1:M);
                phi = phiTrue + angNse((1:M)+M);
                
                senErr = reshape(chol(Qs)'*errTmp(:,m),N,M);
                senPos = senPosTrue + senErr;
                
                % EV solution
                tic;
                mprSol1 = AOA3DLocMPR_EV( theta, phi, senPos, Qa, Qs );
                th(1,m) = mprSol1(1);
                ph(1,m) = mprSol1(2);
                g(1,m) = mprSol1(3);
                totalTime(1) = totalTime(1) + toc;
                
                % BR solution
                tic;
                mprSol2 = AOA3DLocMPR_BR( theta, phi, senPos, Qa, Qs );
                th(2,m) = mprSol2(1);
                ph(2,m) = mprSol2(2);
                g(2,m) = mprSol2(3);
                totalTime(2) = totalTime(2) + toc;
            end
            
            % performance statistics
            algNum = size(g,1);
            for ia = 1:algNum
                mse_t(ir,in,is,ia) = mean((thetaSrc-th(ia,:)).^2);
                mse_p(ir,in,is,ia) = mean((phiSrc-ph(ia,:)).^2);
                mse_g(ir,in,is,ia) = mean((gSrc-g(ia,:)).^2);
                bias_t(ir,in,is,ia) = abs(mean(th(ia,:)) - thetaSrc);
                bias_p(ir,in,is,ia) = abs(mean(ph(ia,:)) - phiSrc);
                bias_g(ir,in,is,ia) = (mean(g(ia,:)) - gSrc);
            end
        end
    end
end

% -- plot results --
symbs = ['o','^'];
name = {'EV-MPR','BR-MPR'};

switch model
    case 'nse'
        xlabtext = '10log(\sigma_a^2(rad^2))';
        xdata = nsePwr;
        yl_mse = [-70, 20;
            -100,-10];
        yl_bias = [-120,10;
            -160,0];
    case 'err'
        xlabtext = '10log(\sigma_s^2(m^2))';
        xdata = errLvl;
        yl_mse = [-65, 20;
            -95,-30];
        yl_bias = [-130,10;
            -170,-40];
    case 'rag'
        xlabtext = 'Range(m)';
        xdata = srcRange;
        yl_mse = [-35, 5;
            -65,-25];
        yl_bias = [-100,0;
            -120,0];
end

% MSE
figure(1); clf;
for ia = 1:algNum
    plot(xdata,10*log10(reshape(mse_t(:,:,:,ia)+mse_p(:,:,:,ia),[],1)),symbs(ia),'linewidth',1.5,'DisplayName',name{ia});hold on;grid on;
end
plot(xdata,10*log10(reshape(crlb_t(:,:,:)+crlb_p(:,:,:),[],1)),'-','linewidth',1.5,'DisplayName','CRLB');
% plot(xdata,10*log10(reshape(covEV_t(:,:,:)+covEV_p(:,:,:),[],1)),'-','linewidth',1.5,'DisplayName','covEV');
% plot(xdata,10*log10(reshape(covBR_t(:,:,:)+covBR_p(:,:,:),[],1)),'-','linewidth',1.5,'DisplayName','covBR');
xlabel(xlabtext,'FontSize',13);ylabel('10log(MSE(\theta,\phi)(rad^2))','FontSize',13);
leg1 = legend('show','Location','Northwest');
set(leg1,'FontSize',11);
ylim(yl_mse(1,:));

figure(2); clf;
for ia = 1:algNum
    plot(xdata,10*log10(reshape(mse_g(:,:,:,ia),[],1)),symbs(ia),'linewidth',1.5,'DisplayName',name{ia});hold on;grid on;
end
plot(xdata,10*log10(reshape(crlb_g(:,:,:),[],1)),'-','linewidth',1.5,'DisplayName','CRLB');
% plot(xdata,10*log10(reshape(covEV_g(:,:,:),[],1)),'-','linewidth',1.5,'DisplayName','covEV');
% plot(xdata,10*log10(reshape(covBR_g(:,:,:),[],1)),'-','linewidth',1.5,'DisplayName','covBR');
xlabel(xlabtext,'FontSize',13);ylabel('10log(MSE(g)(1/m^2))','FontSize',13);
leg1 = legend('show','Location','Northwest');
set(leg1,'FontSize',11);
ylim(yl_mse(2,:));

% bias
figure(3); clf;
for ia = 1:algNum
    plot(xdata,20*log10(reshape(sqrt(bias_t(:,:,:,ia).^2+bias_p(:,:,:,ia).^2),[],1)),symbs(ia),'linewidth',1.5,'DisplayName',name{ia});hold on;grid on;
end
xlabel(xlabtext,'FontSize',13);ylabel('20log(Bias(\theta,\phi)(rad))','FontSize',13);
leg1 = legend('show','Location','Northwest');
set(leg1,'FontSize',11);
ylim(yl_bias(1,:));

figure(4); clf;
for ia = 1:algNum
    plot(xdata,20*log10(reshape(bias_g(:,:,:,ia),[],1)),symbs(ia),'linewidth',1.5,'DisplayName',name{ia});hold on;grid on;
end
xlabel(xlabtext,'FontSize',13);ylabel('20log(Bias(g)(1/m))','FontSize',13);
leg1 = legend('show','Location','Northwest');
set(leg1,'FontSize',11);
ylim(yl_bias(2,:));

clear mse_t mse_p mse_g bias_t bias_p bias_g crlb_t crlb_p crlb_g;


