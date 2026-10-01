# ------------------------------------------------------------------------
#                        Particles     
# ------------------------------------------------------------------------
mutable struct powder_struct
    id::Int64
    position::SVector{2, Float64}
    velocity::SVector{2, Float64}
    acceleration::SVector{2, Float64}
    radius::Float64
    mass::Float64

    rigidbody::Int64
    softbody::Int64

    active::Int64
    collision::Int64
    gravity::Int64

    material::String
end

mutable struct gas_struct
    id::Int64
    position::SVector{2, Float64}
    velocity::SVector{2, Float64}
    acceleration::SVector{2, Float64}
    radius::Float64 # half the grid size
    mass::Float64

    rigidbody::Int64
    softbody::Int64

    active::Int64
    collision::Int64
    gravity::Int64
    sph::Int64

    time_active::Int64
    lifetime::Int64

    material::String
end

mutable struct liquid_struct
    id::Int64
    position::SVector{2, Float64}
    velocity::SVector{2, Float64}
    acceleration::SVector{2, Float64}
    radius::Float64
    mass::Float64

    rigidbody::Int64
    softbody::Int64

    density::Float64
    pressure::Float64
    target_density::Float64
    stiff_coef::Float64       
    viscosity_coef::Float64 

    active::Int64
    collision::Int64
    gravity::Int64
    sph::Int64

    color_id::Int64
    material::String
end

mutable struct solid_struct
    id::Int64
    position::SVector{2, Float64}
    velocity::SVector{2, Float64}
    acceleration::SVector{2, Float64}
    radius::Float64
    mass::Float64

    rigidbody::Int64
    softbody::Int64

    active::Int64
    collision::Int64
    gravity::Int64

    material::String
end

# ------------------------------------------------------------------------
#                        Rigidbodies     
# ------------------------------------------------------------------------
mutable struct rigidbody_struct
    id::Int
    particle_indices::Vector{Int}
    cm::SVector{2, Float64}
    V::SVector{2, Float64}
    ω::SVector{3, Float64}
    M::Float64
    bonds::Vector{Tuple{Int,Int}}   
    break_threshold::Float64
end

# ------------------------------------------------------------------------
#                        Softbodies     
# ------------------------------------------------------------------------
mutable struct softbody_struct
    particle_indices::Vector{Int}
    constraints::Vector{Tuple{Int,Int,Float64}}   # (local_i, local_j, rest_length)
    stiffness::Float64
    pinned::Vector{Bool}
end