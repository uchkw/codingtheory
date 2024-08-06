# Personal Libraries for GaloisFields.jl
using GaloisFields
using Polynomials
using Formatting
isdefined(Main, :F2) || const F2 = @GaloisField 2

exsize(ff::DataType) = Int(floor(log2(length(ff))))
fieldsize(a::GaloisFields.AbstractExtensionField)::Int = length(typeof(a))
# assuming that the field is a power of 2
function exsize(a::GaloisFields.AbstractExtensionField)::Int
    fs = fieldsize(a)
    exsz = 1
    while true
        fs >>= 1
        if isone(fs)
            return exsz
        end
        exsz += 1
    end
end

# local で宣言されたGaloisFieldに対してこの関数で原始元を取得できない 2024-04-12
# eval はグローバル変数にしか有効じゃない．
primitiveroot(FF::DataType) = eval(FF.parameters[3])

# dec to binary converter
function de2bi(d::Int; width::Int = 0)::Array{Int, 1}
    bw = 1
    d > 0 && (bw = Int(floor(log2(d))) + 1)
    bw > width && (width = bw)
    b = zeros(Int, width)
    i = 1
    while d > 0
        b[i] = d & 0x1
        d >>= 1
        i += 1
    end
    return b
end

de2bi(d::UInt; width::Int = 0) = de2bi(Int(d); width)
de2bi(d::UInt8; width::Int = 0) = de2bi(Int(d); width)
de2bi(d::UInt16; width::Int = 0) = de2bi(Int(d); width)

de2f2(d::Int; width::Int = 0)::Array{F2, 1} = map(F2, de2bi(d; width=width))
de2f2poly(d::Int; width::Int = 0)::Polynomial{F2, :x} = Polynomial(map(F2, de2bi(d; width=width)))

function bi2de(b::Array{F,1}) where F <: GaloisFields.AbstractGaloisField
    sum([b[i].n * 2^(i - 1) for i = eachindex(b)])
end
# convert binary vector to decimal. you can also specify bit width
bi2de(b::Array{Int8,1}) = sum([b[i] * 2^(i - 1) for i = 1:length(b)])
bi2de(b::Array{Int64,1}) = sum([b[i] * 2^(i - 1) for i = 1:length(b)])
bi2de(α::F) where F <: GaloisFields.AbstractGaloisField = bi2de(bvec(α))

# get the GF(2) vector corresponding to the element 
function bvec(α::Fe)::Vector{F2} where Fe <: GaloisFields.AbstractExtensionField
    return [ iszero((α.n >> (i-1)) & 1) ? zero(F2) : one(F2) for i in 1:exsize(α) ]
end

function bvec(x::T, bw::Int=0)::Array{F2, 1} where T <: Integer
    if bw == 0
        bw = ceil(Int, log2(x))
    end
    bx = de2bi(x, width=bw)
    return map(F2, bx)
end

function bvec(p::Polynomial{F2,:x}, bw::Int=0)::Array{F2, 1}
    if iszero(bw) || iszero(bw-length(p))
        return p.coeffs
    else
        return vcat(p.coeffs,zeros(F2, bw-length(p)))
    end
end

# get binary weights of the F2 vector
wt(v::Array{F2,1})::Int = sum([v[i].n for i = eachindex(v)])

function F2p(b::Array{Fb, 1}, α::Fe) where Fb <: GaloisFields.AbstractGaloisField where Fe <: GaloisFields.AbstractExtensionField
    @assert length(b) == log2(length(typeof(α)))
    basis = [ α^i for i in 0:length(b)-1 ]
    return sum(basis .* b)
end

function F2p(Fe::DataType, b::Array{Fb, 1}) where Fb <: GaloisFields.AbstractGaloisField
    @assert length(b) == exsize(Fe)
    α = primitiveroot(Fe)
    basis = [ α^i for i in 0:length(b)-1 ]
    return sum(basis .* b)
end


function hex(v::Array{F2,1}, wordsize=32)
    nwd = ceil(Int, length(v)/wordsize)
    if length(v) < wordsize
        wordsize = length(v)
    end
    d = ceil(Int,wordsize/4)
    fs = FormatSpec("0$d"*"x")
    vv = zeros(nwd)
    for w in 1:nwd
        for i in 1:wordsize
            idx = wordsize*(w-1)+i
            idx > length(v) && break
            vv[w] += v[idx].n * 2^(i - 1)
        end
    end

    ret = ""
    for w in nwd:-1:1
        ret *= fmt(fs,vv[w])
    end
    return ret
end

function fromF2mattoIntmat(M::Array{F, 2}) where F <: GaloisFields.AbstractGaloisField
    Mi = zeros(Int, size(M))
    for j in 1:size(M)[2]
        for i in 1:size(M)[1]
            Mi[i,j] = Int(M[i,j] == F(1))
        end
    end
    return Mi
end

#%%
# Assuming GF size of 2^m
function trace(a::F)::F2 where F <: GaloisFields.AbstractGaloisField
    r = a
    for i in 1:exsize(F)-1
        r += a^(2^i)
    end
    return ifelse(iszero(r), F2(0), F2(1))
end

function getTrOneElem(a::F)::F where F <: GaloisFields.AbstractGaloisField
    for i in 0:fieldsize(a)-2
        if isone(trace(a^i))
            return a^i
        end
    end
end

function cuberoot(a::F, α::F)::F where F <: GaloisFields.AbstractGaloisField
    if iszero(a)
        return zero(F)
    else
        N = length(F) - 1
        if N % 3 != 0
            i = log(a,α)
            for j in 0:N-1
                if i == (3*j % N)
                    return α^j
                end
            end
        else
            return a^(N ÷ 3)
        end
    end        
end

function Base. *(a::Vector{F}, b::Vector{F2})::F where F <: GaloisFields.AbstractGaloisField
    @assert length(a) == length(b)
    result = F(0) # FF型の初期値を設定
    @inbounds for i in 1:length(a)
        result += a[i] * b[i]
    end
    return result
end
function Base. *(a::Vector{F2}, b::Vector{F})::F where F <: GaloisFields.AbstractGaloisField
    return Base.*(b, a)
end

function Base.inv(M::Array{F, 2}) where F <: GaloisFields.AbstractGaloisField
    N = zeros(F, size(M))
    matsize = size(M)[1]
    for i in 1:matsize
        N[i,i] = F(1)
    end
    MN = hcat(M, N)
    for j in 1:matsize
        if MN[j,j] == F(0)
            # pivot search
            for i in j+1:matsize
                if MN[i,j] != F(0)
                    MN[j,:], MN[i,:] = MN[i,:], MN[j,:]
                    break
                end
            end
        end
        if MN[j,j] != F(1) # pivotを1に
            MN[j,:] *= inv(MN[j,j])
        end
        if MN[j,j] == F(0) # error
            println("No pivot at $j column")
            return nothing
        end
        # make all elements under the pivot zero on the jth column
        for i in j+1:matsize
            if MN[i,j] != F(0)
                MN[i,:] += MN[i,j] * MN[j,:]
            end
        end
        # make all elements under the pivot zero on the jth column
        for i in j-1:-1:1
            if MN[i,j] != F(0)
                MN[i,:] += MN[i,j] * MN[j,:]
            end
        end
    end
    return MN[:,matsize+1:end]
end

getconwaypolynomialofdegree(d::Int) = Polynomial(GaloisFields.conwaypolynomial(2,d))

function Base.log(a::F, α::F) where F <: GaloisFields.AbstractGaloisField
    if iszero(a)
        return Inf
    end
    for i in 0:length(F)-2
        if a == α^i
            return i
        end
    end
end

function Base.log(a::F; base::F=α) where F <: GaloisFields.AbstractGaloisField
    if iszero(a)
        return Inf
    end
    for i in 0:length(F)-2
        if a == base^i
            return i
        end
    end
end


#エラーが出る．2024-04-12
function Base.log(a::F) where F <: GaloisFields.AbstractGaloisField
    if iszero(a)
        return Inf
    end
    p = primitiveroot(F) # 体から原始元を取得できなみたい．なぜだ
    for i in 0:length(F)-2
        if a == p^i
            return i
        end
    end
end

function getconjugates(β::F)::Array{F,1} where F <: GaloisFields.AbstractExtensionField
    ret = [β]
    i = 1
    x = β^2
    while x != β
        push!(ret, x)
        i += 1
        x = β^(2^i)
    end
    ret
end

function getroots(px::Polynomial{F}, α::F)::Array{F,1} where F <: GaloisFields.AbstractExtensionField
    N = fieldsize(α)-2
    ret = Vector{F}()
    for i in 0:N
        if iszero(px(α^i))
            push!(ret, α^i)
        end
    end
    return ret
end

# from Polynomials v.0.6.1
# return the polynomial with roots r
function poly(r::AbstractVector{T}, var::Polynomials.SymbolLike=:x) where {T}
    n = length(r)
    c = zeros(T, n+1)
    c[1] = one(T)
    for j = 1:n
        for i = j:-1:1
            c[i+1] = c[i+1]-r[j]*c[i]
        end
    end
    return Polynomial(reverse(c), var)
end

function getminimumpolynomial(β::F)::Polynomial{F} where F <: GaloisFields.AbstractExtensionField
    poly(getconjugates(β))
end

# orderを調べる
function getorder(genpoly::Polynomial{F}) where F <: GaloisFields.AbstractGaloisField
    F2 = @GaloisField 2
    A = makecompanionmatrix(Array{F2}(genpoly.coeffs[1:end-1]))
    s = zeros(F2, size(A)[1])
    s[1] = 1
    o = 1
    n = A*s
    while n != s
        # @show n,s;
        o += 1
        n = A*n
    end
    return o
end

order(x) = getorder(x)

# orderを調べる
function getorder(x::F) where F <: GaloisFields.AbstractGaloisField
    for i in 2:length(typeof(x))
        if isequal(x^i, F(1))
            return i
        end
    end
    return -1
end

function onehotvector(len::Int, pos::Int, v::F)::Array{F,1} where F <: GaloisFields.AbstractGaloisField
    ret = zeros(F, len)
    ret[pos] = v
    ret
end

# MSBを落とす
function makecompanionmatrix(g::Array{F,1})::Array{F,2} where F <: GaloisFields.AbstractGaloisField
    A = zeros(F, length(g), length(g))
    A[:,end] = copy(g)
    for i in 2:length(g)
        A[i,i-1] = F(1)
    end
    A
end

# primitiveroot関数を使っているため，ローカルスコープで定義されたGaloisFieldに使用できない
function makecompanionmatrix(a::F)::Array{F2,2} where F <: GaloisFields.AbstractGaloisField
    d = length(bvec(a))
    A = zeros(F2, d, d)
    v = a
    p = primitiveroot(typeof(a))
    for i in 1:d
        A[:,i] = bvec(v)
        v *= p
    end
    A
end

# a はコンパニオン行列の最初の列の要素, αは基底となる原始元
function makecompanionmatrix(a::F, α::F)::Array{F2,2} where F <: GaloisFields.AbstractGaloisField
    d = length(bvec(a))
    A = zeros(F2, d, d)
    v = a
    for i in 1:d
        A[:,i] = bvec(v)
        v *= α
    end
    A
end

function isprimitive(p::Polynomial{F2})::Bool
    if getorder(p) == 2^Polynomials.degree(p) - 1
        return true
    end
    return false
end

# generate the list of primitive polynomials of argument degree
function generateprimitivepolynomials(degree::Integer)::Array{Polynomial{F2}, 1}
    polylist = []
    for i in 2^degree+1:2:2^(degree+1) - 1
        p = de2f2poly(i)
        if isprimitive(p)
            push!(polylist, p)
        end
    end
    return polylist
end

# return number of non-zero coefficients
function numterms(f::Polynomial{F}) where F <: GaloisFields.AbstractGaloisField
    return sum(f.coeffs .!= 0)
end

######## この辺以下は別のファイルにしたほうがいいかも．Coding theoryにより過ぎてる
function getbvec(A::Array{F,2}, i::Int)::Array{F,1} where F <: GaloisFields.AbstractGaloisField
    b = zeros(F, size(A)[1])
    b[1] = 1
    for j in 1:i
        b = A * b
    end
    return b
end

# b1からp個A回転させた列ベクトル並べた行列を作る．
function makepreprocmatrix(A::Array{F,2}, b1::Array{F,1}, p::Int) where F <: GaloisFields.AbstractGaloisField
    @assert size(A)[1] == length(b1)
    Bp = zeros(F, size(A)[1], p)
    Bp[:,1] = b1
    for j in 2:p
        Bp[:,j] = A * Bp[:,j-1]
    end
    Bp
end

# パリティ検査行列のバイナリ表現を生成する．
function makeHmat(genpoly::Polynomial{F}, Hlen::Int = 0)::Array{F2, 2} where F <: GaloisFields.AbstractGaloisField
    g = Array{F2}(genpoly.coeffs[1:end-1])
    A = makecompanionmatrix(g)
    b = getbvec(A, 0)
    if Hlen == 0 
        Hlen = getorder(genpoly)
    end
    H = makepreprocmatrix(A, b, Hlen)
    return H
end

# m=8だとオフセットを1(z=0)にできない
# 2bit訂正時に使用するY行列を生成する
function generateYvec(α::Fe, z::Int=0)::Array{Fe,1} where Fe <: GaloisFields.AbstractExtensionField
    if iszero(trace(α^z))
        for i in 0:length(α)-2
            if isone(trace(α^i))
                z = i
                break
            end
        end
    end
    dim = exsize(α)
    Y = zeros(Fe, dim)
    for i in 1:dim
        im1 = i - 1
        basis = α^im1 + trace(α^im1) * α^z
        for j in 0:length(typeof(α))-2
            yi = α^j
            y = yi + yi^2
            if isequal(y, basis)
                Y[i] = yi
                break
            end
        end
    end
    return Y
end

"""
Return polynomial coefficients with increasing degree order

"""
function logcoeffs(f::F, b::Fb) where F <: GaloisFields.AbstractGaloisField where Fb <: GaloisFields.AbstractGaloisField
    return map(x -> log(x, b), f.coeffs)
end
function logcoeffs(f::F) where F <: GaloisFields.AbstractGaloisField
    return map(x -> log(x), f.coeffs)
end


"""
Return polynomial coefficients with increasing degree order

Examples
≡≡≡≡≡≡≡≡≡≡

julia> logcoeffs(Polynomial([α^10, 1]), α)
2-element Array{Int64,1}:
 10
  0
"""
function logcoeffs(f::Polynomial{F}, b::F) where F <: GaloisFields.AbstractGaloisField
    return map(x -> log(x, b), f.coeffs)
end
function logcoeffs(f::Polynomial{F}) where F <: GaloisFields.AbstractGaloisField
    return map(x -> log(x), f.coeffs)
end

# べき表示にするのにもっといい方法ないものか．．．
function gfpretty(a::F, α::F) where F <: GaloisFields.AbstractGaloisField
    p = log(a, α)
    return "α^"*string(p)
end
function gfpretty(a::F) where F <: GaloisFields.AbstractGaloisField
    p = log(a)
    return "α^"*string(p)
end
function gfpretty(v::Vector{F}) where F <: GaloisFields.AbstractGaloisField
    map(gfpretty, v)
end
# show polynomial representation of binary vector
function gfpretty(b::Array{Int8,1})
    for i = 1:length(b)
        if b[i] != 0
            i == 1 ? print("1") : print("x")
            if i > 2
                print("^", i - 1)
            end
            if i < length(b)
                print(" + ")
            end
        end
    end
    print("\n")
end


function extract_degrees(polynomial::AbstractString)
    degrees = Vector{Int64}()

    # 項ごとに分割して処理する
    terms = split(polynomial, "+")
    for term in terms
        # 次数の部分を取り出す
        match_result = match(r"x\^(\d+)", term)
        if match_result !== nothing
            degree = parse(Int64, match_result.captures[1])
            push!(degrees, degree)
        else
            if occursin("x",term)
                push!(degrees, 1)
            elseif occursin("1",term)
                push!(degrees, 0)
            end
        end
    end

    # 次数が高い順にソートして返す
    return sort(degrees, rev=true)
end

function string2coefvec(polynomial::AbstractString)
    degrees = extract_degrees(polynomial)
    degreearray = zeros(Int,maximum(degrees)+1)
    for i in degrees
        degreearray[i+1] = 1
    end
    return degreearray
end

function string2F2poly(polynomial::AbstractString)::Polynomial{F2, :x}
    f2array = map(F2,string2coefvec(polynomial))
    return Polynomial(f2array)
end

# return inverse matrix
# if the matrix is singular, return zero matrix
function Base.inv(_A::Matrix{FE})::Array{FE} where FE <: GaloisFields.AbstractExtensionField
    A = copy(_A)
    Asize = size(A, 1)
    @assert size(A) == (Asize, Asize)
    B = zeros(FE, Asize, Asize)
    for i in 1:Asize
        B[i, i] = FE(1)
    end
    
    # Forward elimination
    for i in 1:Asize
        # Find pivot
        pivot_row = i
        # @show pivot_row i A[pivot_row, i]
        if iszero(A[pivot_row, i])
            ofs = findfirst(!iszero, A[pivot_row+1:end, i])
            if isnothing(ofs)
                return zeros(FE)  # Singular matrix, no inverse
            end
            pivot_row += ofs
        end
        A[i, :], A[pivot_row, :] = A[pivot_row, :], A[i, :]
        B[i, :], B[pivot_row, :] = B[pivot_row, :], B[i, :]
        
        # Normalize pivot row
        B[i, :] *= inv(A[i, i])
        A[i, :] *= inv(A[i, i])
        
        # Eliminate below
        for j in i+1:Asize
            B[j, :] += B[i, :] * A[j, i]
            A[j, :] += A[i, :] * A[j, i]
        end
    end
    
    # Backward substitution
    for i in Asize:-1:2
        for j in i-1:-1:1
            B[j, :] += B[i, :] * A[j, i]
            A[j, :] += A[i, :] * A[j, i]
        end
    end
    
    return B
end

function isprimitive(a::FF)::Bool where FF <: GaloisFields.AbstractExtensionField
    for i in 1:length(FF)-2
        if isone(a^i)
            return false
        end
    end
    return true
end

function primitiveelements(a::FF)::Array{FF, 1} where FF <: GaloisFields.AbstractExtensionField
    pelements = [a]
    for i in 2:length(FF)-2
        if isprimitive(a^i)
            push!(pelements, a^i)
        end
    end
    return pelements
end
