function [BPAs,Yf,E] = getBPAs(P,ind_Tr,ind_Ev,Yt,Ye,Yapp,Aux)
iF = 1;
%% LAG

E = zeros(1,2);
Yf = zeros(length(Ye),length(P.Aux)+1);

tail = P.L;
y = arrayfun(@(i) Yt(end-24*(2*tail-i+1)+1:end-24*(tail+1-i)),1:tail+1,'UniformOutput',false);
y = horzcat(y{:});
E(1) = mean(abs(mean(y(:,1:tail),2)-y(:,tail+1))); %1 Lag 24h x 7

for i=1:length(ind_Ev)
    y = arrayfun(@(lag) Yapp(ind_Ev(i)-lag*24),1:tail);
    Yf(i,iF) = mean(y(~isnan(y))); %1 Lag 24h x 14
end
iF = iF+1;

%% AUX

for j=1:length(P.Aux)

    x = Aux{ind_Tr,P.Aux(j)};
    y = Yt(~isnan(x));
    x = x(~isnan(x));
    w = exp(-P.alpha*((x-Aux{ind_Tr(end+1-24*7:end),P.Aux(j)}')/std(x)).^2);
    E(iF) = mean(abs(y'*w./sum(w)-y(end+1-24*7:end)'));
    for i=1:length(ind_Ev)
        xe = Aux{ind_Ev(i),P.Aux(j)};
        w = exp(-P.alpha*((x-xe)/std(x)).^2);
        Yf(i,iF) = y'*w/sum(w);    
    end

    iF = iF+1;

end

%% BUILD BPAs
%E = (E+1).^5/sum((E+1).^5)*mean(E);
BPAs = cell(length(ind_Ev),length(P.Aux)+1);
for i=1:length(ind_Ev)
    x = linspace(0,3,P.M+1)';
    x = x(2:end);
    %x = (0.5:1:3)';
    for j=1:length(P.Aux)+1
        BPAs{i,j}.sets = [Yf(i,j)-x*abs(E(j)) Yf(i,j)+x*abs(E(j))];
        BPAs{i,j}.masses =  (2*x*abs(E(j)))/sum(2*x*E(j));
    end
end
