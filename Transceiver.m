function [trans_X,trans_X_inv]=Transceiver(M,T) 
% This function aims at generating signals to be transmitted at BS
% X_temp=randn(T,T);
% [U,~,~]=svd(X_temp);
% trans_X=U(1:M,:);
% P_sum=T; %The total energy

if M<T
   TEMP=dftmtx(T);
   trans_X=1/sqrt(T)*TEMP(1:M,1:T);
   trans_X_inv=trans_X'*inv(trans_X*trans_X');
else
   TEMP=dftmtx(M);
   trans_X=1/sqrt(M)*TEMP(1:M,1:T);
   trans_X_inv=inv(trans_X'*trans_X)*trans_X';
end
end