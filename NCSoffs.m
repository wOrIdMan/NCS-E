function offs = NCSoffs(pop,sigma,lb,ub,gen,maxGen)
[m,d] = size(pop);
offs = zeros(m,d);
for i = 1:m
    parenti = i;
    randv = (randn(1,d).*sigma(parenti,:))';
    v  = pop(parenti,:) + randv';

%     v(v<lb) = 2*lb(v<lb)-v(v<lb);
%     v(v>ub) = 2*ub(v>ub)-v(v>ub);
%     v(v<lb) = lb(v<lb);
%     v(v>ub) = ub(v>ub);
    
    w = find(v < lb);
    if ~isempty(w)
        l=rand;
        if l < 0.5
            v(w) = 2 * lb(w) -  v(w);
            w1 = find( v(w) > ub(w));
            if ~isempty(w1)
                v(1, w(w1)) = ub(w(w1));
            end
        else
            if gen<maxGen*0.6
                v(w) =  ub(w);
            else
                v(w) =  lb(w);
            end
        end
    end
    
    y = find(v > ub);
    if ~isempty(y)
        l=rand;
        if l<0.5
            
            v(y) =  2 * ub(y) - v(y);
            y1 = find(v(y) < lb(y));
            if ~isempty(y1)
                v(y(y1)) = ub(y(y1));
            end

        else
           if gen < maxGen*0.4
               v(y) =  lb(y);   
           else
               v(y) =  ub(y);  
           end
        end
    end



    offs(i,:) = v;
end
