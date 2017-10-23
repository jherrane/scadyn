function [ C_sca ] = get_Csca( S11, nang )
%GET_CSCA Summary of this function goes here
%   Detailed explanation goes here
n = size(S11,1);
times = mod(n,nang) + 1;
C_sca = zeros(times,1);
for i = 1:times
    for j = 1:nang
        C_sca(i) = C_sca(i) + S11(j+(times-1)*i);
    end
    C_sca(i) = C_sca(i)/nang;
end

end

