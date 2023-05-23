function [LM,bin_LM_len,I] = LocationMap(I,a,b)
[A,B] = size(I);
flow_map = zeros(1,floor((A-2)/a)*floor((B-2)/b));

n = 1;
for i = 1:floor((A-2)/a)
    for j = 1:floor((B-2)/b)
        
        X = I(a*(i-1)+1:a*i,b*(j-1)+1:b*j);
        X = X(:);
        
        ind = find(X==255,1,'first');
        ind2 = find(X==0,1,'first');
        
        if isempty(ind) && isempty(ind2)
            flow_map(n) = 0;
        else
            flow_map(n) = 1;
        end
         
        n = n + 1;
        
        
    end
end

    LM = flow_map;
    
    %after compression
    cPos_x = cell(1,1);%ËãÊõ±àÂëÑ¹Ëõ
    cPos_x{1} = flow_map;
    if isempty(flow_map)
        bin_LM_len = 0;
        
    else
        bin_LM_len = sum(flow_map>0);
    end
    
end