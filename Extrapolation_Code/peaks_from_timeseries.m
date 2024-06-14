function [M, T]=peaks_from_timeseries(x, t)
Mend=ceil(length(x)/2);
M=zeros(1,Mend);
T=zeros(1,Mend);
ind=0;
if ~exist('t', 'var')
    t=1:length(x);
end
ii=[false, diff(x)==0];
x(ii)=[];
t(ii)=[];
for i=2:length(x)-1
    if x(i-1)<=x(i) && x(i)>x(i+1)
        ind=ind+1;
        M(ind)=x(i);
        T(ind)=t(i);
    end
end;
M=M(1:ind);
T=T(1:ind);
