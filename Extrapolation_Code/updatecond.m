function [condition] = updatecond(cond) 

cond = +cond; % Convert logical values to numeric

dc = diff(cond);    % Finding intermitten data 
[~,poz,val] = find(dc);
if ~isempty(poz);
    leng = diff([0 poz length(cond)]);
    N = length(leng);
    if val(1) == -1 % first block contains ones
        T = max(leng(1:2:N)); % ones go odd
    else
        T = max(leng(2:2:N)); % first block contains zeros, ones go even
    end
    n = find(leng == T); % position of the biggest
    cond = 0*cond;
    cond( sum(leng(1:n-1))+1 : sum(leng(1:n-1))+T ) = 1;
    
%     plot(cond,'o');
%     hold on; plot(cond1,'.r')
end

condition = logical(cond);

end
