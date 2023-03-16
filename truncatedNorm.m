function sumsquares = truncatedNorm(value, qc)

N = 2^qc-1;
sumsquares = 0;

for n = 0:1:N

    sumsquares = sumsquares+((value^n)/sqrt(factorial(n)))^2;

end

end
    