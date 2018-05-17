__precompile__()

module QuickShiftClustering

isinstalled(a) = isa(Pkg.installed(a), VersionNumber)
isinstalled("PyPlot") && begin
    using PyPlot
end

using NearestNeighbors, ProgressMeter, Distances, FunctionalData

export quickshift, quickshiftlabels, quickshiftplot

mutable struct QuickShift
    rootind
    links
    sigma
end

function dist(a::Array{Float32,2}, i::Int, b::Array{Float32,2}, j::Int)
    sum = 0f0
    for d = 1:size(a,1)
        sum += (a[d,i]-b[d,j])^2
    end
    sum
end

function link(i, inds, G, data)
    mindist = typemax(eltype(data))
    minind = 0
    for n = inds
        if G[n]>G[i]
            d = dist(data,i,data,n)
            if d < mindist
                mindist = d
                minind = n::Int
            end
        end
    end
    minind, mindist
end

function gauss(data, n, ind, factor1, factor2)
    s = 0f0
    for m = ind
        s += exp(-dist(data,n,data,m)*factor1)
    end
    factor2 * s
end

quickshift(data, a...) = quickshift(convert(Array{Float32,2},data), a...)
quickshift(data::Array{Float32,2}, sigma) = quickshift(data, convert(Float32,sigma))
function quickshift(data::Array{Float32,2}, sigma::Float32=convert(Float32,median(pairwise(Euclidean(), randsample(data,1000)))/100))
    tree = KDTree(data)
    # @show sigma
    N = len(data)
    factor1 = 1f0 / (2*sigma^2) ::Float32
    factor2 = 1/(2*pi*sigma^2*N)
    G = zeros(Float32, 1, N)
    nninds = [Array(Int,0) for n in 1:N]
    @showprogress 1 "Computing kernel distances ... " for n = 1:N
        knnind = @p knn tree vec(at(data,n)) 100 true | fst
        ind = length(knnind) > 10 ? knnind : 1:N
        nninds[n] = ind
        G[n] = gauss(data, n, ind, factor1, factor2)
    end
    # println("median lenghts:",(@p map nninds length | flatten | median))
    links = [Any[] for i in 1:N]
    rootind = -1
    inflength = typemax(eltype(G))
    Nrange = 1:N
    minind = 0
    mindist = inflength
    @showprogress 1 "Linking ... " for i in 1:N
        for inds = [nninds[i]; Nrange]
            minind, mindist = link(i, inds, G, data)
            if minind != 0
                break
            end
        end
        if mindist == inflength
            if rootind < 0
                rootind = i
            else
                push!(links[rootind], (sqrt(mindist), i))
            end
        else
            push!(links[minind], (sqrt(mindist), i))
        end
    end
    QuickShift(rootind, links, sigma)
end

function quickshiftlabels(a::QuickShift, maxlength = 10*a.sigma)
    labels = zeros(Int32,length(a.links))
    cut_internal(a.rootind, a.links, labels, maxlength, 1, 2) 
    labels
end

function cut_internal(ind, links, labels, maxlength, label, maxlabel) 
    labels[ind] = label
    for x in links[ind]
        if x[1] > maxlength
            maxlabel += 1
        end
        maxlabel = max(label, cut_internal(x[2], links, labels, maxlength, x[1] > maxlength ? maxlabel : label, maxlabel)) 
    end
    maxlabel
end

function quickshiftplot(a::QuickShift, data, labels)
    if !isdefined(:PyPlot)
        error("quickshiftplot needs PyPlot installed and loaded using 'using PyPlot'")
    end
    if size(data,1) != 2
        error("quickshiftplot only works on 2D data, i.e. size(data,1)==2")
    end

    for i = 1:len(a.links)
        from = data[:,i]
        for x = a.links[i]
            to = data[:,x[2]]
            p = hcat(from,to)'
            plot(fst(p),snd(p),"b-")
        end

    end
    scatter(data[1,:],data[2,:], c = Array{Int, 1}(labels), edgecolor = "none")
end

end
