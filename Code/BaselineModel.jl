
# what structure here? Maybe make so we can simulate across sites??

mutable struct Parameters
	np::NamedTuple{(:ψ, :δ, :pi_d, :ρ, :pi_k, :ζc, :ζf, :βw, :vω, :scale, :αl, :αθ, :k0, :ρa),NTuple{14,Int64}}
	lb::NamedTuple{(:ψ, :δ, :pi_d, :ρ, :pi_k, :ζc, :ζf, :βw, :vω, :scale, :αl, :αθ, :k0, :ρa),Tuple{Array{Int64,1},Array{Int64,1},Int64,Int64,Int64,Int64,Int64,Array{Int64,1},Array{Int64,1},Int64,Array{Float64,1},Array{Float64,1},Int64,Float64}}
	ub::NamedTuple{(:ψ, :δ, :pi_d, :ρ, :pi_k, :ζc, :ζf, :βw, :vω, :scale, :αl, :αθ, :k0, :ρa),Tuple{Array{Int64,1},Array{Int64,1},Int64,Int64,Int64,Int64,Int64,Array{Int64,1},Array{Int64,1},Int64,Array{Float64,1},Array{Float64,1},Int64,Float64}}

    αc::Float64 #<- coefficient on consumption
    αθ::Float64 #<- coefficient on skills
    αH::Float64 #<- disutility of work
    αP::Float64 #<- disutility of program participation (add something here?)
    δI::Array{Float64,2} #<- Cobb Douglas share on investment (skill by age)
    δθ::Array{Float64,2} #<- Cobb Douglas share on own skills (skill by age)
    p̃::Array{Float64,1} #<- relative price of investment good (home vs exterior care)
    ϵ::Array{Float64,1} #<- elasticity of investments across periods

    wage::Array{Float64,1} #<- wage at each site
    β::Float64 #<- discounting (won't matter in lots of sites)
end

# How do we want to do this!?



# function to take parameters and simulate a bunch of things
