function [Phi]=Phase_Generate(P,N)
% This function aims at generating random pilot phase matrix
% Phi=exp(1i*2*pi*rand(P,N));% randomly generation
% 
% X_temp=randn(P,P);
% [U,~,~]=svd(X_temp);
% Phi=U(1:P,1:N);

% Phase matrix for ALS based algorithm
TEMP=dftmtx(N);
Phi=1/sqrt(N)*TEMP(1:P,1:N);

% % Phase matrix for AMP based algorithm (seems impractical)
% Phi=zeros(P,N);
% for p=1:P
%     non_index=randi([1,N],round(0.8*N),1);
%     Phi(p,non_index)=exp(1i*2*pi*rand);
% end