function f = cost(x)
    p = [1,0;
        0,1;
        -1,-1
        1,1];
    r = [1,1,sqrt(3),1.4];
    f = zeros(5,1);
    for i = 1:4
        f(i+1) = norm(x-p(i,:)') - r(i);
    end
end