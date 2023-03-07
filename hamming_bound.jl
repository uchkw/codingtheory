# Hamming bound

function volume(n::Int64, r::Int64)
    s = 1
    for i=1:r
        s += binomial(big(n),i)
    end
    return s
end

#%%
# Hamming bound: upper bound of minimum distance
function get_dmin_ub(n::Int64, k::Int64)::Int64
    dmin = 1
    for d=2:n-k
        if log2(volume(n, floor(Int64, (d-1)/2))) > n-k
            println("dmin_ub= $dmin")
            break
        else
            dmin = d
        end
    end
    return dmin    
end

#%%

# get_dmin_ub(511,493)
# get_dmin_ub(512,493)
# get_dmin_ub(144,128)
