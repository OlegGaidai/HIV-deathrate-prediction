function [ q0, b0, a0, c0, r ] = guess(bar_lev, tail_mark, process, ACER)

b0 = mean2(process);
if b0 > tail_mark
    b0 = tail_mark;
end
c0 = 2;

[r slope intersept] = regression((bar_lev-b0).^c0,log(ACER)); % y = a*x + b linear regression
a0 = -slope;
q0 = exp(intersept);

end
