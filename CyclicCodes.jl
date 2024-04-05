isdefined(Main, :F2) || "GaloisFieldUtils.jl"
# Matlabチックな関数
hex2dec(h::String) = parse(Int, h, base = 16)

dec2hex(d::Int) = string(d, base = 16)

# 生成行列を巡回シフトして生成行列を生成する
function makeCircGmat(gx::Polynomial{F2,:x}, n::Int)::Array{F2, 2}
    k = n - Polynomials.degree(gx)
    G = zeros(F2, (k,n))
    g = bvec(gx, n)
    for i in 1:k
        G[i,:] = circshift(g,i-1)
    end
    return G
end

# [q, r] = a/b
# 引数の多項式はde2biで入力されることを想定
function gfdeconv(
    a::Array{Int8,1},
    b::Array{Int8,1},
)::Tuple{Array{Int8,1},Array{Int8,1}}
    deg_b = length(b)
    deg_a = length(a)
    if deg_a < deg_b
        return ([Int8(0)], a)
    end
    q = zeros(Int8, deg_a - deg_b + 1)
    q_i = length(q)
    for d = deg_a:-1:deg_b
        if a[d] != 0
            q[q_i] = 1
            a[d-deg_b+1:d-1] = a[d-deg_b+1:d-1] .⊻ b[1:end-1]
        else
            q[q_i] = 0
        end
        q_i -= 1
    end
    return (q, a[1:deg_b-1])
end

struct LFSR
    poly::Int
    reg_size::Int
    reg_mask::Int
    reg_fdbk::Int
end

function generate_lfsr(poly::Int)
    reg_size = Int(floor(log2(poly)))
    msb_pos = 1 << reg_size
    fdbk = poly - msb_pos
    if fdbk == 0
        println("No feedback polynomial! It must be wrong.")
        exit()
    end
    lfsr = LFSR(poly, reg_size, msb_pos - 1, fdbk)
end

# LFSRを１クロック進めて最上位ビットと新たなレジスタ状態を返す．
function clock(reg, lfsr::LFSR)
    reg <<= 1
    obit = reg >> lfsr.reg_size
    obit == 1 && (reg = (reg & lfsr.reg_mask) ⊻ lfsr.reg_fdbk)
    return obit, reg
end

# LFSRの周期を返す．
function get_period(lfsr::LFSR)
    reg = 1
    obit, reg = clock(reg, lfsr)
    i = 1
    while reg != 1
        obit, reg = clock(reg, lfsr)
        i += 1
    end
    return i
end

# LFSRの系列を返す．
function get_sequence(lfsr::LFSR)
    reg = 1
    seq = zeros(Int8, 1)
    seq[1], reg = clock(reg, lfsr)
    i = 1
    while reg != 1
        obit, reg = clock(reg, lfsr)
        push!(seq, obit)
        i += 1
    end
    return seq
end

# # 上で定義したシフトレジスタを使って検査行列を作成
function make_parity_check_matrix(lfsr::LFSR)
    period = get_period(lfsr)
    H = zeros(Int8, lfsr.reg_size, period)
    reg = 1
    H[end, end] = 1
    obit, reg = clock(reg, lfsr)
    for j = (period-1):-1:1
        bvec = de2bi(reg)
        for i = 1:length(bvec)
            H[end+1-i, j] = bvec[i]
        end
        obit, reg = clock(reg, lfsr)
    end
    return H
end

using LinearAlgebra
# 末尾mbitが単位行列であること前提
function make_generator_matrix(H::Array{Int8,2})
    m, n = size(H)
    k = n - m
    return hcat(Matrix{Int8}(I, k, k), transpose(H[:, 1:k]))
end

using SymPy
# MacWilliam変換を計算
# A = 重み分布, k = Aの重み分布を持つ符号の情報長, n = Aの重み分布を持つ符号の符号長
mw_transform(A::Sym, x::Sym, k::Int, n::Int) =
    (A((1 - x) / (1 + x)) * (1 + x)^n) / 2^k

function get_dualcodes_wd(hex_poly::String, n::Int)
    px, rx = gfdeconv(de2bi(hex2dec(hex_poly)), de2bi(3)) #(x+1)で割る
    (rx[1] == 1) && (px = de2bi(hex2dec(hex_poly)))
    px_lfsr = generate_lfsr(bi2de(px))
    seq = get_sequence(px_lfsr) # pxの周期を確認．CCITTなら32767
    n_max = length(seq)
    if rx[1] == 0
        W = zeros(BigInt, n)
        N = sum(seq[1:n])
        W[N] = 1
        for i = 2:n_max
            w_end = n + i - 1
            (w_end > n_max) && (w_end %= n_max)
            x = seq[w_end] # binary digit entering the window
            if seq[i-1] ⊻ x == 1
                N += x == 1 ? 1 : -1
            end
            W[N] += 1
        end
        @assert W[n] == 0
        B = zeros(BigInt, n)
        # All 1はWに含まれないことを仮定
        B[n] = 1 # (x+1)を持つ場合はall-1がある．
        for i = 1:n-1
            B[i] = W[i] + W[n-i]
        end
        return B
    else
        B = zeros(BigInt, n)
        N = sum(seq[1:n])
        B[N] = 1
        for i = 2:n_max
            w_end = n + i - 1
            (w_end > n_max) && (w_end %= n_max)
            x = seq[w_end] # binary digit entering the window
            if seq[i-1] ⊻ x == 1
                N += x == 1 ? 1 : -1
            end
            B[N] += 1
        end
        return B
    end
end

function prob_ud_using_dual(p::Float64, n::Int, r::Int, B::Array{BigInt,1})
    pud = BigFloat(0.0)
    for i = n:-1:1
        (B[i] > 0) && (pud += B[i] * BigFloat(1 - 2p)^i)
    end
    pud += BigFloat(1.0) # all-0 codewordの寄与で初期化
    pud *= BigFloat(0.5)^r
    pud -= BigFloat(1 - p)^n
    return pud
end

# errnum 次の全ての係数が 1 の有限体の係数を持つ多項式を返す．
function get_biterr_poly(errnum::Int, FE::DataType=FF)::Polynomial{FE}
    return Polynomial([FF(1) for i in 1:errnum], :x)
end

# errnum 次の全ての係数が 1 の有限体の係数を持つ多項式を返す．
function get_rnderr_poly(errnum::Int, FE::DataType=FF)::Polynomial{FE}
    return Polynomial(rand(FF, errnum), :x)
end

# Obtain the syndrome from the error or received polynomial
function calc_syndrome(ex::Polynomial{FE}, num::Int, x::FE=α)::Array{FE, 1} where FE <: GaloisFields.AbstractExtensionField
    return ex.([x^i for i in 0:num-1])
end

# Determine the error location polynomial by the PGZ algorithm
function get_elp_by_pgz(S::Array{FE, 1})::Polynomial{FE} where FE <: GaloisFields.AbstractExtensionField
    t = length(S) ÷ 2
    v = t
    Smat = [S[i+j] for i in 1:v, j in 0:v-1]
    invS = inv(Smat)
    while iszero(invS)
        v -= 1
        if iszero(v)
            return FE(0)
        end
        Smat = [S[i+j] for i in 1:v, j in 0:v-1]
        invS = inv(Smat)
    end
    v = size(Smat,2)
    Svec = [S[v+i] for i in 1:v]
    σ = invS * Svec
    append!(σ, FE(1)) # add maximum degree coefficient "1"
    return Polynomial(σ, :x)
end
