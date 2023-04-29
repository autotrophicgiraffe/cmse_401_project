function _norm!(A, x)
    m, n = size(A)
    rows = rowvals(A)
    vals = nonzeros(A)
    Threads.@threads for i = 1:n
        @simd for j in nzrange(A, i)
          vals[j] *= x[rows[j]]*x[i]
       end
    end
    nothing
end

using SparseArrays

n = 10^6;
A = sprand(n,n,0.0001);
x = randn(Float64, n);

_norm!(A, x)
@time _norm!(A,x)
  # 0.688902 seconds (5 allocations: 224 bytes) # without multithreading
  # 0.189398 seconds (5 allocations: 224 bytes) # with multithreading
