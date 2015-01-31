function [NC] = FourierConstants(C, e)
    K = [];
    m = [];
    i = 1;
    while abs(sum(K .^ 2) - 1) > e
        k = C(i);
        if abs(k) > e
            m(end + 1) = i;
            K(end + 1) = k;
        end
        i = i + 1;
    end
    NC = zeros(2, length(m));
    NC(1, :) = m;
    NC(2, :) = K;