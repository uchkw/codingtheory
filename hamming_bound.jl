# Hamming bound

function volume(n::Int64, r::Int64)
    s = 1
    for i=1:r
        s += binomial(n,i)
    end
    return s
end

#%%
function get_achievable_dmin(n::Int64, k::Int64)::Int64
    dmin = 1
    for d=2:n-k
        if log2(volume(n, floor(Int64, (d-1)/2))) > n-k
            println("dmin= $dmin")
            break
        else
            dmin = d
        end
    end
    return dmin    
end

#%%

get_achievable_dmin(511,493)
get_achievable_dmin(512,493)
get_achievable_dmin(144,128)
