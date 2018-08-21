function pointmin(_I :: Array{Float64,2})
    if isempty(_I)
        return false
    end
    I = copy(_I)
    Fx=zeros(size(I))
    Fy=zeros(size(I))
    Fz=zeros(size(I))

    J=zeros( tuple([size(I,i)+2 for i = 1:ndims(I)]...) )
    # J[:,:]=maximum(I[:])
    J .= maximum(I[:])

    J[2:end-1,2:end-1]=I
    Ne = [-1 -1; -1  0; -1  1; 0 -1; 0  1; 1 -1;  1  0; 1  1]
    for i=1:maximum(size(Ne))
        In=J[2+Ne[i,1]:end-1+Ne[i,1],2+Ne[i,2]:end-1+Ne[i,2]]
        check = In.<I #TODO Check this is what is intended (In<I)
        I[check]= copy(In[check])
        D=Ne[i,:]
        D=D./sqrt(sum(D.^2))
        Fx[check] .= D[1]
        Fy[check] .= D[2]
    end

    return (Fx,Fy)
end

function pointmin(_I :: Array{Float64,3})
    I = copy(_I)
    Fx=zeros(size(I))
    Fy=zeros(size(I))
    Fz=zeros(size(I))

    J=zeros( tuple([size(I,i)+2 for i = 1:ndims(I)]...) )
    J[:,:]=maximum(I[:])

    J[2:end-1,2:end-1,2:end-1]=I
    Ne=[-1 -1 -1; -1 -1  0; -1 -1  1; -1  0 -1; -1  0  0; -1  0  1; -1  1 -1; -1  1  0; -1  1  1;        
    0 -1 -1;  0 -1  0;  0 -1  1;  0  0 -1;            0  0  1;  0  1 -1;  0  1  0;  0  1  1; 
        1 -1 -1;  1 -1  0;  1 -1  1;  1  0 -1;  1  0  0;  1  0  1;  1  1 -1;  1  1  0;  1  1  1];
    
    for i=1:maximum(size(Ne))
        In = J[2+Ne[i,1]:end-1+Ne[i,1],2+Ne[i,2]:end-1+Ne[i,2],2+Ne[i,3]:end-1+Ne[i,3]]
        check = In.<I
        I[check]= In[check]
        D=Ne[i,:]
        D=D./sqrt(sum(D.^2))
        Fx[check]= D[1]
        Fy[check]= D[2]
        Fz[check]= D[3]
    end

    return (Fx,Fy,Fz)
end

