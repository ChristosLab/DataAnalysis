function out = nn_segment(q)
a=diff(q);
b=find([a inf]>1);
c=diff([0 b]); 
d=cumsum(c); 
out = [[d - c + 1]', c'];
end