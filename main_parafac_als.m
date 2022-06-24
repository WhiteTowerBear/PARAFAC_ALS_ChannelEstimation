%========================================================a
%  Author: Wei Li;
%  Date:   11/17/ 2019,  SUTD;
%  Version: V2.0 
%  Note: This main function is ALS-based PARAFAC channel estimation method

%The references to theorems and equations refer to the following paper:
%
% L. Wei, C. Huang, G. C. Alexandropoulos, C. Yuen, Z. Zhang and M. Debbah, 
% "Channel Estimation for RIS-Empowered Multi-User MISO Wireless 
% Communications," in IEEE Transactions on Communications, vol. 69, 
% no. 6, pp. 4144-4157, June 2021.

%License: If you in any way use this code for research that results in 
% publications, please cite our original article listed above.
%%========================================================
clc
clear all 
fid=fopen('result.txt','a+');
%%========================================================
%  Set the parameters
M=16; % The number of antennas at BS
K=16; % The number of users
N=16; % The total RIS elements
% T=N+1;
T=16; % The numebr of time slots
N1=16; % The number of elements placed on sub-RIS
sub_block=N/N1; % The number of sub-RIS
P=16;  % The number of pilot phase matrices
var_channel=1;    % The variance of reflecting channel 

FRAME=100;
iter=15;
diff_temp=[];
snr=1;


for SNR= 0:5:30  
    SNR
    mc_sum_Hs=0;
    mc_sum_Hr=0;
    mc_sum_H=0;
    mc_sum_Hr=0;
    mmse_sum_Hs=0;
    mmse_sum_Hr=0;
    CRB_total=0;
    CRB_aa_total=0;
    CRB_bb_total=0;
   
    for frame=1:1:FRAME
    % Generate transmitted signals X
    [X,X_inv]=Transceiver(M,T); 
    var_noise=10^(-0.1*SNR);
    %%========================================================
    %  Pass channel
    Hr= binornd(1,1,K,N).*(sqrt(var_channel/2)*(randn(K,N)+1i*randn(K,N))); 
    Hs=binornd(1,1,N,M).*(sqrt(var_channel/2)*(randn(N,M)+1i*randn(N,M)));
  
    Hs(:,1)=1;
    % Generate phase matrix Phi
    [Phi]=Phase_Generate(P,N1);
    % Received signals
    for tt=1:sub_block
        [Hs_est1((tt-1)*N1+1:tt*N1,:),Hr_est1(:,(tt-1)*N1+1:tt*N1)]=parafac(K,T,P,M,N1,X,X_inv,Hr(:,(tt-1)*N1+1:tt*N1),Hs((tt-1)*N1+1:tt*N1,:),iter,var_noise,var_channel,Phi);
    end
%%========================================================%
%%
% Simulation result
      mc_sum_Hs=mc_sum_Hs+norm(Hs-Hs_est1,'fro')^2/(norm(Hs,'fro')^2);
      mc_sum_Hr=mc_sum_Hr+norm(Hr-Hr_est1,'fro')^2/(norm(Hr,'fro')^2);
      mc_sum_H=mc_sum_H+norm(Hr*Hs-Hr_est1*Hs_est1,'fro')^2/(norm(Hr*Hs,'fro')^2);
    end
    diff_Hs(1,snr)=mc_sum_Hs/FRAME
    diff_Hr(1,snr)=mc_sum_Hr/FRAME
    diff_H(1,snr)=mc_sum_H/FRAME 
    snr=snr+1;
end
fprintf(fid,'\n \n M=%d  K=%d  N=%d  T=%d  P=%d  %d Monte Carlo simulations \n',M,K,N,T,P,FRAME);
fprintf(fid,' %s\r\n ',datestr(now,31));
fprintf(fid,'\n NMSE of H is \n');
fprintf(fid, ' %d  ', diff_H);
fclose(fid);
figure
semilogy(0:5:30,diff_Hr,'-bo','linewidth',2)
hold on
semilogy(0:5:30,diff_Hs,'-ys','linewidth',2)
hold on
semilogy(0:5:30,diff_H,'-ms','linewidth',2)


legend({'Hr (ALS)','Hs (ALS)','Cascaded channel (ALS)'});
title(['M=',num2str(M),',K=',num2str(K),',T=',num2str(T), ',N=',num2str(N),',P=',num2str(P)]);
xlabel('SNR (dB)');
ylabel('NMSE');
grid on
saveas(gcf,'CRB1.fig')
