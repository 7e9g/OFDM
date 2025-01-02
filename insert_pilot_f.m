%函数的功能为添加导频信号
function output = ...
    insert_pilot_f(input,pilot_seq,pilot_inter,is_pilot_k)
    
if is_pilot_k==1
    [N,NL] = size(input);
    NL_p = NL+ceil(NL/pilot_inter);
    output = zeros(N,NL_p);
    j = 1;
    for i=1:NL_p
        if mod(i,pilot_inter+1)==1
            output(:,i)=pilot_seq;
        else
            output(:,i)=input(:,j);
            j = j+1;
        end
    end
else
end

end