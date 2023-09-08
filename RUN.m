% Load the database with timeVector (Nx1), powerVector (Nx1) and Aux (N*m)
% Aux refers to weather features values at each sample time
% timeVector is a datetime vetor
% powerVector is a numeric vector
% Aux is a table

%Set prediction time (refered in publication as t0) in variable tP
%tp = ...;
P.Pr = [tP tP+1];
P.Tr = [timeVector(1) P.Pr(1)];

%Select usable weather features indices
%P.Aux = [...];
P.maxP = max(powerVector);

%Set method parameters
%P.alpha = ...;
%P.M = ...;
%P.L = ...;

%Get dataset configuration 
[ind_Tr,ind_Ev,Tt,Yt,Te,Ye,Yapp] = getConf(P,timeVector,powerVector);

%Build mass functions (BPA stands for basic probability assignment)
[BPAs,Yf,E] = getBPAs(P,ind_Tr,ind_Ev,Yt,Ye,Yapp,Aux);

%Yager Combination method
Yp = combine(P,BPAs);

%Evaluate prediction
MAPE = mean(abs(Yp-Ye)./Ye);
MAE = mean(abs(Yp-Ye));
MSE = mean((Yp-Ye).^2);
RMS = sqrt(MSE);
NRMSE = RMS/mean(Ye);

%Plot results
subplot(2,1,1)
plot(Te,Ye,Te,Yf)
legend([{'N/A'} arrayfun(@(i) sprintf('MAPE = %6.2f',nanmean(abs(Yf(:,i)...
    -Ye)./Ye)),1:size(Yf,2),'UniformOutput',false)])
ylabel('Power (in kW)')
title('Before combination')
subplot(2,1,2)
plot(Te,Ye,Te,Yp)
legend({'N/A',sprintf('MAPE = %6.2f',MAPE)})
title('After combination')
ylabel('Power (in kW)')