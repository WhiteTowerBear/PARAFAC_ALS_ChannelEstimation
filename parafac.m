function [Hs_est1,Hr_est1]=parafac(K,T,P,M,N1,X,X_inv,Hr,Hs,iter,var_noise,var_channel,Phi)
% This function estimates two channels based on PARAFAC
%The references to theorems and equations refer to the following paper:
%
% L. Wei, C. Huang, G. C. Alexandropoulos, C. Yuen, Z. Zhang and M. Debbah, 
% "Channel Estimation for RIS-Empowered Multi-User MISO Wireless 
% Communications," in IEEE Transactions on Communications, vol. 69, 
% no. 6, pp. 4144-4157, June 2021.

%License: If you in any way use this code for research that results in 
% publications, please cite our original article listed above.

noise=zeros(K,T,P);
    rec_y=zeros(K,T,P);
    for p=1:P
        noise(:,:,p)=sqrt(var_noise/2)*(randn(K,T)+1i*randn(K,T));
        rec_y(:,:,p)=Hr*diag(Phi(p,:))*Hs*X+noise(:,:,p);
        rec_y_TEMP(:,:,p)=rec_y(:,:,p)*X_inv;
    end
%     %%========================================================
%     % Receiver
%     % Mode-1, mode-2, mode-3 unfolded forms of PARAFAC    
    for m=1:M
        for p=1:P
            for k=1:K
                Z_KP_M((p-1)*K+k,m)=rec_y_TEMP(k,m,p); % MODE2
                Z_PM_K((m-1)*P+p,k)=rec_y_TEMP(k,m,p); % MODE1
            end
        end
    end
  %%  
    % Initialization
    Hs_est=zeros(N1,M,iter+1);
    Hr_est=zeros(K,N1,iter+1);
    Hr_est(:,:,1)=sqrt(var_channel/2)*(randn(K,N1)+1i*randn(K,N1));
    
    A1=zeros(P*M,N1,iter);
    A1_inv=zeros(N1,P*M,iter);
    A2=zeros(K*P,N1,iter);
    A2_inv=zeros(N1,K*P,iter);
    
    %%========================================================
    %% PARAFAC 
    % Iteration between channels
    com_Hs=ones(iter+1,1);
    for i=2:iter        
        A2(:,:,i-1)=kr(Phi,Hr_est(:,:,i-1));
        A2_inv(:,:,i-1)=inv(A2(:,:,i-1)'*A2(:,:,i-1))*A2(:,:,i-1)';
        Hs_est(:,:,i)=A2_inv(:,:,i-1)*Z_KP_M;
        for n=1:N1 % for scaling ambiguity removal later:
            Hs_est(n,:,i) = Hs_est(n,:,i) / Hs_est(n,1,i);
        end
        A1(:,:,i)=kr(Hs_est(:,:,i).',Phi);
        A1_inv(:,:,i)=inv(A1(:,:,i)'*A1(:,:,i))*A1(:,:,i)';
        Hr_est(:,:,i)=(A1_inv(:,:,i)*Z_PM_K).';
        fit=0;
        for p=1:P
            fit=fit+norm(Z_KP_M-kr(Phi,Hr_est(:,:,i))*Hs_est(:,:,i),'fro')^2;
        end
        com_Hs(i,1)=fit;
        delta(i,1)=(com_Hs(i-1,1)-com_Hs(i,1))/com_Hs(i,1);
        if abs(delta(i,1))<1e-5
            break
        end
    end   
    Hs_est1=Hs_est(:,:,i);
    Hr_est1=Hr_est(:,:,i);
end