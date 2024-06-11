q = parallel.pool.DataQueue;
afterEach(q, @calccc);
parfor i = 1:10
    send(q, i)
end

function out = calccc(in)
out = mod(in, 2);
disp(out)
end