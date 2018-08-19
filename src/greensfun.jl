# HierarchicalMatrix of GreensFuns on TensorSpace of HierarchicalSpaces

converttoPiecewiseSpace(H::Space) = H
converttoPiecewiseSpace(H::HierarchicalSpace) = PiecewiseSpace(H)

function partition(ss::ApproxFun.TensorSpace{Tuple{HS1,HS2}}) where {HS1<:HierarchicalSpace,HS2<:HierarchicalSpace}
    ss11,ss12 = partition(factor(ss.space,1))
    ss21,ss22 = partition(factor(ss.space,2))
    (ss11⊗ss21,ss12⊗ss22),(converttoPiecewiseSpace(ss12)⊗converttoPiecewiseSpace(ss21),converttoPiecewiseSpace(ss11)⊗converttoPiecewiseSpace(ss22))
end

function partition(ss::CauchyWeight{O,Tuple{HS1,HS2}}) where {O,HS1<:HierarchicalSpace,HS2<:HierarchicalSpace}
    ss11,ss12 = partition(factor(ss.space,1))
    ss21,ss22 = partition(factor(ss.space,2))
    (CauchyWeight(ss11⊗ss21,O),CauchyWeight(ss12⊗ss22,O)),(CauchyWeight(converttoPiecewiseSpace(ss12)⊗converttoPiecewiseSpace(ss21),O),CauchyWeight(converttoPiecewiseSpace(ss11)⊗converttoPiecewiseSpace(ss22),O))
end

function GreensFun(f::DFunction,ss::AbstractProductSpace{Tuple{HS1,HS2}};method::Symbol=:lowrank,kwds...) where {HS1<:HierarchicalSpace,HS2<:HierarchicalSpace}
    (ss11,ss22),(ss21,ss12) = partition(ss)
    meth1 = method == :Cholesky || method == :lowrank ? :standard : method
    G11 = GreensFun(f,ss11;method=method,kwds...)
    G22 = GreensFun(f,ss22;method=method,kwds...)
    G21 = GreensFun(LowRankFun(f,ss21;method=meth1,kwds...))
    G12 = method == :Cholesky ? transpose(G21) : GreensFun(LowRankFun(f,ss12;method=meth1,kwds...))
    return HierarchicalMatrix((G11,G22),(G21,G12))
end

blocksize(H::HierarchicalMatrix{F,G}) where {F<:GreensFun,G<:GreensFun} = map(ncomponents,domain(H).domains)

function domain(H::HierarchicalMatrix{F,G}) where {F<:GreensFun,G<:GreensFun}
    H11,H22 = diagonaldata(H)
    H21,H12 = offdiagonaldata(H)
    m1,n2 = factor(domain(H12),1),factor(domain(H12),2)
    m2,n1 = factor(domain(H21),1),factor(domain(H21),2)
    @assert (m1,n1) == (factor(domain(H11),1),factor(domain(H11),2))
    @assert (m2,n2) == (factor(domain(H22),1),factor(domain(H22),2))
    (m1∪m2)*(n1∪n2)
end

for TYP in (:(ApproxFun.DefiniteLineIntegralWrapper),:DefiniteLineIntegral)
    @eval function Base.getindex(⨍::$TYP,H::HierarchicalMatrix{G,GreensFun{L,T}}) where {G<:GreensFun,L<:LowRankFun,T}
        H11,H22 = diagonaldata(H)
        wsp = components(domainspace(⨍))
        if ncomponents(factor(domain(H11),2)) ≥ 2
            ⨍1 = DefiniteLineIntegral(PiecewiseSpace(wsp[1:ncomponents(factor(domain(H11),2))]))
        else
            ⨍1 = DefiniteLineIntegral(wsp[1])
        end
        if ncomponents(factor(domain(H22),2)) ≥ 2
            ⨍2 = DefiniteLineIntegral(PiecewiseSpace(wsp[end-ncomponents(factor(domain(H22),2))+1:end]))
        else
            ⨍2 = DefiniteLineIntegral(wsp[end])
        end
        HierarchicalOperator((qr(⨍1[H11]),qr(⨍2[H22])),
                             map(LowRankIntegralOperator,offdiagonaldata(H)))
    end
end

qr(H::HierarchicalOperator) = H # trivial no-op for now.
