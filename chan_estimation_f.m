function output = chan_estimation_f(input,H,pilot_seq,pilot_interval)
%信道估计
[N,NL] = size(input);
ML = size(H,2);
H_out = zeros(N,ML);
output = zeros(N,NL);
for ii=1:ML
    H_out(:,ii)=H(:,ii)./pilot_seq;
end
for jj=1:ML

    h = inv(diag(H_out(:,jj),0));
    if jj==ML
        Y = input(:,((jj-1)*pilot_interval+1):end);
        output(:,((jj-1)*pilot_interval+1):end)=h*Y;
    else
        Y = input(:,((jj-1)*pilot_interval+1):jj*pilot_interval);
        output(:,((jj-1)*pilot_interval+1):jj*pilot_interval)=h*Y;
    end
end
    
end