function N = Normalize(x, f)
    I = trapz(x, f .^ 2);
    N = f / sqrt(I);
end