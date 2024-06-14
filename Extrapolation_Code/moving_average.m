function [CI_smooth] = moving_average(CI, lag)

% ---------------------------------------------------------------------- %
%                                                                        %
%   Function that smooths the provided confidence intervals by           %
%   moving average                                                       %
%                                                                        %
%   Inputs:                                                              %
%     CI: Confidence intervals                                           %
%     lag: Number of points to be taken into aeraing                     %
%                                                                        %
%   Outputs:                                                             %
%    CI_smooth:   Smoothed confidence intervals with indicated           %
%                 polynomial                                             %
%                                                                        %

CI_smooth = nan(1,length(CI));

for kk=((lag+1)/2):(length(CI)-(lag-1)/2)
    start = kk-((lag-1)/2);
    stop = kk+((lag-1)/2);
    CI_smooth(kk) = sum(CI(start:stop))/lag;
end
