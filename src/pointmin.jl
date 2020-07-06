function pointmin(_I::AbstractArray{T}) where T
    if isempty(_I)
        return false
    end
    I = _I
    Fx=zeros(T,size(I))
    Fy=zeros(T,size(I))
    Fz=zeros(T,size(I))

    J=zeros(T, tuple([size(I,i)+2 for i = 1:ndims(I)]...) )
    @views J .= maximum(I[:])

    @views J[2:end-1,2:end-1] = I
    Ne = [-1 -1; -1  0; -1  1; 0 -1; 0  1; 1 -1;  1  0; 1  1]
    for i = 1:maximum(size(Ne))
        @views In=J[2+Ne[i,1]:end-1+Ne[i,1],2+Ne[i,2]:end-1+Ne[i,2]]
        check = In .< I
        I[check] = In[check]
        @views D = Ne[i,:]
        D = D ./ sqrt(sum(D.^2))
        @views Fx[check] .= D[1]
        @views Fy[check] .= D[2]
    end
    return (Fx,Fy)
end