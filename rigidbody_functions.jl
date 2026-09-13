# -----------------------------------------------------------------------------
#                           RigidBody physics
# ----------------------------------------------------------------------------- 
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

function rigidbody_physics(particles, rigidbodies)

    for rb in rigidbodies

        # translation
        #F_gravity = @SVector zeros(2)
        F_gravity = calculate_gravity(rb.cm, rb.M, 0, nothing)

        rb.V += (F_gravity / rb.M) * dt      
        translation = rb.V * dt
        new_cm = rb.cm + translation         

        # rotation
        angle = rb.ω[3] * dt
        cos_a = cos(angle)
        sin_a = sin(angle)

        for idx in rb.particle_indices
            p = particles[idx]

            r = p.position - rb.cm

            r_rot = SVector(cos_a*r[1] - sin_a*r[2], sin_a*r[1] + cos_a*r[2])

            p.position = new_cm + r_rot

            r_new = p.position - new_cm
            p.velocity = rb.V + SVector(-rb.ω[3]*r_new[2], rb.ω[3]*r_new[1])
        end

        rb.cm = new_cm
    end
end

function calculate_inertia(particles, rb)
    inertia = 0.0
    for i in rb.particle_indices
        p = particles[i]
        r = p.position - rb.cm
        inertia += p.mass * (r[1]^2 + r[2]^2)
    end
    return max(inertia, 1.0)
end

function calculate_inertia_tensor(particles, rb)

    I_tensor = zeros(2,2)
    I3 = Matrix{Float64}(I, 2, 2)  # identity

    for i in rb.particle_indices
        
        p = particles[i]
        r = p.position .- rb.cm   

        r2 = dot(r, r)
        rrT = r * transpose(r)
        I_tensor .+= p.mass .* (r2 .* I3 .- rrT)
    end

    return I_tensor
end

function calculate_center_of_mass(particles)
    total_mass = 0.0
    cm_x = 0.0
    cm_y = 0.0
    
    for p in particles
        total_mass += p.mass
        cm_x += p.mass * p.position[1]
        cm_y += p.mass * p.position[2]
    end
    
    return SVector(cm_x / total_mass, cm_y / total_mass), total_mass
end

function create_cube!(particles, rigidbodies, id, offset, v_init, ω_init, m, n)
    particle_radius = grid_size/2
    particle_diam = 2 * particle_radius

    positions = SVector{2,Float64}[]
    for row in 0:(m-1)
        for col in 0:(n-1)
            push!(positions, SVector(col * particle_diam, row * particle_diam))
        end
    end

    indices = Int[]

    for pos in positions
        p = solid_struct(
            length(particles)+1,
            offset .+ pos,
            @SVector(zeros(2)),
            @SVector(zeros(2)),
            particle_radius,
            3.0,
            id,
            0,
            1,              # active
            1,              # collision
            1,              # gravity
            "solid"
        )
        push!(particles, p)
        push!(indices, length(particles))
    end

    # Calculate center of mass
    cube_particles = [particles[i] for i in indices]
    cm, total_mass = calculate_center_of_mass(cube_particles)

    # Set initial velocities
    for i in indices
        r = particles[i].position - cm
        particles[i].velocity = v_init + SVector(-ω_init[1]*r[2], ω_init[1]*r[1])
    end

    bonds = build_grid_bonds(positions, particle_diam, indices)

    rb = rigidbody_struct(
        id,
        indices, #global indices of particles in the rigidbody
        cm,
        SVector(v_init[1], v_init[2]),
        SVector(0.0, 0.0, ω_init[1]),
        total_mass,
        bonds,
        99
    )
    push!(rigidbodies, rb)
end


function build_grid_bonds(local_positions, spacing, indices)
    bonds = Tuple{Int,Int}[]
    n = length(local_positions)
    for a in 1:n
        for b in a+1:n
            d = norm(local_positions[a] - local_positions[b])
            if d < spacing * 1.1   # adjacent (no diagonals)
                push!(bonds, (indices[a], indices[b]))   # global ids
            end
        end
    end
    return bonds
end

# check the two groups of bounds
# subtract and recalculate the info of the first
# create a new body with the rest
function split_rigidbody!(particles, rigidbodies, rb, broken_bond)

    a_global, b_global = broken_bond

    # ensure the broken bond endpoints are part of this rigidbody
    particle_set = Set(rb.particle_indices)
    if !(a_global in particle_set) || !(b_global in particle_set)
        return
    end

    # build adjacency among current rb particles and filter invalid bonds
    adj = Dict{Int, Vector{Int}}()
    for pid in rb.particle_indices
        adj[pid] = Int[]
    end
    valid_bonds = Tuple{Int,Int}[]
    for bond in rb.bonds
        x, y = bond
        if x in particle_set && y in particle_set
            push!(valid_bonds, bond)
            push!(adj[x], y)
            push!(adj[y], x)
        end
    end
    rb.bonds = valid_bonds

    # flood-fill (BFS) from a_global to find its connected component
    visited = Set{Int}()
    queue = [a_global]
    push!(visited, a_global)
    while !isempty(queue)
        cur = popfirst!(queue)
        for nb in adj[cur]
            if !(nb in visited)
                push!(visited, nb)
                push!(queue, nb)
            end
        end
    end

    # if b_global is still connected, nothing to split
    if b_global in visited
        return
    end

    # partition into global-id groups
    group_a = collect(visited)
    group_b = [pid for pid in rb.particle_indices if !(pid in visited)]

    original_bonds = rb.bonds

    # ---- group_a: either stays as rb, or becomes a free particle if alone ----
    if length(group_a) == 1
        lone_index = group_a[1]
        particles[lone_index].rigidbody = 0
    else
        rb.particle_indices = group_a
        new_bonds_a = Tuple{Int,Int}[]
        for bond in original_bonds
            x, y = bond
            if (x in visited) && (y in visited)
                push!(new_bonds_a, (x, y))
            end
        end
        rb.bonds = new_bonds_a

        piece_particles_a = [particles[i] for i in rb.particle_indices]
        cm_a, mass_a = calculate_center_of_mass(piece_particles_a)
        rb.cm = cm_a
        rb.M = mass_a

        for i in rb.particle_indices
            particles[i].rigidbody = rb.id
        end
    end

    # ---- group_b: either becomes a free particle, or a brand new rigidbody ----
    if length(group_b) == 1
        lone_index = group_b[1]
        particles[lone_index].rigidbody = 0
    else
        new_particle_indices_b = group_b

        new_bonds_b = Tuple{Int,Int}[]
        for bond in original_bonds
            x, y = bond
            if !(x in visited) && !(y in visited)
                push!(new_bonds_b, (x, y))
            end
        end

        piece_particles_b = [particles[i] for i in new_particle_indices_b]
        cm_b, mass_b = calculate_center_of_mass(piece_particles_b)

        new_id = length(rigidbodies) + 1
        for i in new_particle_indices_b
            particles[i].rigidbody = new_id
        end

        new_rigidbody = rigidbody_struct(new_id, new_particle_indices_b, cm_b, rb.V, rb.ω, mass_b, new_bonds_b, rb.break_threshold)
        push!(rigidbodies, new_rigidbody)
    end
end