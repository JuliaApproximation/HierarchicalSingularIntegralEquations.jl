
## HierarchicalDomain

function Base.show(io::IO, H::HierarchicalDomain{S,T}) where {S,T}
    print(io,"$(nlevels(H))-level HierarchicalDomain{$S,$T}:\n")
    show(io,UnionDomain(collectdata(H)))
end

## HierarchicalSpace

function Base.show(io::IO, H::HierarchicalSpace{S,T}) where {S,T}
    print(io,"$(nlevels(H))-level HierarchicalSpace{$S,$T}:\n")
    show(io,PiecewiseSpace(collectdata(H)))
end

## HierarchicalMatrix{F<:GreensFun,G<:GreensFun}
# Base.writemime because HierarchicalMatrix <: AbstractArray

## HierarchicalOperator{U<:Operator,V<:AbstractLowRankOperator}

Base.alignment(io::IO, x::Infinity) = (1,0)
function Base.show(io::IO, ::MIME"text/plain", H::HierarchicalMatrix{F,GreensFun{L,T}}) where {F<:GreensFun,L<:LowRankFun,T}
    print(io,"$(nlevels(H))-level HierarchicalMatrix of GreensFun's with blockwise ranks:\n")
    Base.print_matrix(io,blockrank(H),"["," ","]")
end
function Base.show(io::IO, ::MIME"text/plain", H::HierarchicalOperator{U,V}) where {U<:Operator,V<:AbstractLowRankOperator}
    print(io,"$(nlevels(H))-level HierarchicalOperator with blockwise ranks:\n")
    Base.print_matrix(io,blockrank(H),"["," ","]")
end
