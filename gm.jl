struct DynamicalSystem
    initial::Array{Float64,1}
    map::Function
end
Base.iterate(ds::DynamicalSystem, pt=ds.initial ) = (pt, ds.map(pt) )

function push_points(initial::Array{Float64,1} ,
                    map::Function ,
                    drop::Int64  ,
                    num::Int64 ,
                    plot_x::Array{Float16,1} , plot_y::Array{Float16,1})
    for pt in Iterators.take(
                Iterators.drop(DynamicalSystem(initial,map ),drop ) , num)
        push!(plot_x,pt[1])
        push!(plot_y,pt[2])
    end
end

@enum GMType f1g1 f1g2 f2g1 f2g2

struct GumowskiMira
    gmtype::GMType
    α::Float64
    σ::Float64
    μ::Float64

    # g1,g2 : 有界非線形性 (bounded nonlinerity)
    g1::Function
    g2::Function

    f1::Function
    f2::Function

    function GumowskiMira(gmtype::GMType ;α::Float64=0. ,σ::Float64=0. ,μ::Float64)
        g1=x->μ*x + (2.0(1.0-μ)x^2 )/(1.0 + x^2)
        g2=x->μ*x + (1.0-μ)x^2*exp((1.0 - x^2)/4.0)
        f1=g->function(pt::Array{Float64,1})
            xx =  pt[2] + g(pt[1])
            yy = -pt[1] + g(xx)
            [xx,yy]
        end
        f2=g->function(pt::Array{Float64,1})
            xx =  pt[2] + α*pt[2]*( 1.0 - σ*pt[2]^2 ) + g(pt[1])
            yy = -pt[1] + g(xx)
            [xx,yy]
        end
        new(gmtype,α,σ,μ,g1,g2,f1,f2)
    end

end

gmmap(gm::GumowskiMira)=
                    if      (gm.gmtype==f1g1)
                        gm.f1(gm.g1)
                    elseif (gm.gmtype==f1g2)
                        gm.f1(gm.g2)
                    elseif (gm.gmtype==f2g1)
                        gm.f2(gm.g1)
                    elseif (gm.gmtype==f2g2)
                        gm.f2(gm.g2)
                    else
                        error("gmtypeが不正")
                    end


#gm = GumowskiMira(f2g1,α=0.008 , σ=0.05, μ=-0.496) # 「神話の鳥 (mythic bird)」3枚羽の翼
gm = GumowskiMira(f2g1,α=0.009 , σ=0.05, μ=-0.801) # 「神話の鳥 (mythic bird)」5枚羽の翼
#gm = GumowskiMira(f2g2,α=0.0083, σ=0.1 , μ=-0.38 ) #
#gm = GumowskiMira(f2g2,α=0.01  , σ=0.1 , μ=0.8 ) #
#gm = GumowskiMira(f1g1,μ=0.39 ) #
#gm = GumowskiMira(f1g1,μ=0.365 ) #
#gm = GumowskiMira(f1g1,μ=0.34 ) #

initial = [1.1,1.1]
plot_x=Float16[]
plot_y=Float16[]
push_points(initial,gmmap(gm),20000,200000,plot_x,plot_y)

using Plots
pyplot()
scatter(plot_x,plot_y,markersize=1,markercolor=:black)
