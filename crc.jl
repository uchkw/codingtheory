# misdetection probability
function misdetburst(b::Integer, r::Integer, q)
    # b: burst length
    # r: degree of CRC
    # q: field size of the CRC polynomial
    if b <= r
        return 0
    elseif b == r+1
        return 1.0 / ((Float64(q)-1)*Float64(q)^(r-1))
    else
        return 1.0 / (Float64(q)^r)
    end
end
