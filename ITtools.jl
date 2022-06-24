# binary entropy function
function H_b(p)
    @assert(p >= 0 && p <= 1)
    if (p == 0 || p == 1)
        return 0
    else
        return -p*log2(p) - (1-p)*log2(1-p)
    end
end