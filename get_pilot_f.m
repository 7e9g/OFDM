function [output,H] = get_pilot_f(input,pilot_interval)
%将数据，与块状导频分开
[N,NL] = size(input);
H = zeros(N,ceil(NL/(1+pilot_interval)));
output = zeros(N,NL-ceil(NL/(1+pilot_interval)));
jj=1;
kk=1;
for ii = 1:NL
    if mod(ii,pilot_interval+1)==1
        H(:,jj) = input(:,ii);
        jj=jj+1;
    else
        output(:,kk) = input(:,ii);
        kk=kk+1;
    end
end

end