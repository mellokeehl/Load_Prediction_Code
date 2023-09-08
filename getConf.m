function [ind_Tr,ind_Ev,Tt,Yt,Te,Ye,Yapp] = getConf(P,timeVector,powerVector)
ind_Tr = find(timeVector >= P.Tr(1) & timeVector < P.Tr(2));
ind_Ev = find(timeVector >= P.Pr(1) & timeVector <= P.Pr(2));
Yt = powerVector(ind_Tr);
Tt = timeVector(ind_Tr);
Ye = powerVector(ind_Ev);
Te = timeVector(ind_Ev);
Yapp = nan(length(powerVector),1);
Yapp(ind_Tr) = Yt;