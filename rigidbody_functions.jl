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
    return inertia 
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

    bonds = build_grid_bonds(positions, particle_diam)

    rb = rigidbody_struct(
        id,
        indices,
        cm,
        SVector(v_init[1], v_init[2]),
        SVector(0.0, 0.0, ω_init[1]),
        total_mass,
        bonds,
        50
    )
    push!(rigidbodies, rb)
end


function build_grid_bonds(local_positions, spacing)
    bonds = Tuple{Int,Int}[]
    n = length(local_positions)
    for a in 1:n
        for b in a+1:n
            d = norm(local_positions[a] - local_positions[b])
            if d < spacing * 1.1   # adjacent cell (up/down/left/right), not diagonal
                push!(bonds, (a, b))
            end
        end
    end
    return bonds
end

    # check the two groups of bounds
    # subtract and recalculate the info of the first
    # create a new body with the rest

function split_rigidbody!(particles, rigidbodies, rb, broken_bond)

    a, b = broken_bond
    n = length(rb.particle_indices)

    if a > n || b > n
        return
    end

    valid_bonds = Tuple{Int,Int}[]
    for bond in rb.bonds
        if bond[1] <= n && bond[2] <= n
            push!(valid_bonds, bond)
        end
    end
    rb.bonds = valid_bonds

    in_group_a = falses(n)
    in_group_a[a] = true

    changed = true
    while changed
        changed = false
        for bond in rb.bonds
            x = bond[1]
            y = bond[2]
            if in_group_a[x] && !in_group_a[y]
                in_group_a[y] = true
                changed = true
            elseif in_group_a[y] && !in_group_a[x]
                in_group_a[x] = true
                changed = true
            end
        end
    end

    if in_group_a[b]
        return
    end

    group_a = Int[]
    group_b = Int[]
    for k in 1:n
        if in_group_a[k]
            push!(group_a, k)
        else
            push!(group_b, k)
        end
    end

    original_indices = rb.particle_indices
    original_bonds = rb.bonds

    # ---- group_a: either stays as rb, or becomes a free particle if alone ----

    if length(group_a) == 1
        lone_index = original_indices[group_a[1]]
        particles[lone_index].rigidbody = 0
    else
        new_particle_indices_a = Int[]
        for k in group_a
            push!(new_particle_indices_a, original_indices[k])
        end
        rb.particle_indices = new_particle_indices_a

        new_bonds_a = Tuple{Int,Int}[]
        for bond in original_bonds
            x = bond[1]
            y = bond[2]
            if in_group_a[x] && in_group_a[y]
                push!(new_bonds_a, (x, y))
            end
        end
        rb.bonds = new_bonds_a

        piece_particles_a = []
        for i in rb.particle_indices
            push!(piece_particles_a, particles[i])
        end
        cm_a, mass_a = calculate_center_of_mass(piece_particles_a)
        rb.cm = cm_a
        rb.M = mass_a

        for i in rb.particle_indices
            particles[i].rigidbody = rb.id
        end
    end

    # ---- group_b: either becomes a free particle, or a brand new rigidbody ----

    if length(group_b) == 1
        lone_index = original_indices[group_b[1]]
        particles[lone_index].rigidbody = 0
    else
        new_particle_indices_b = Int[]
        for k in group_b
            push!(new_particle_indices_b, original_indices[k])
        end

        remap = Dict{Int,Int}()
        for new_local in 1:length(group_b)
            old_local = group_b[new_local]
            remap[old_local] = new_local
        end

        new_bonds_b = Tuple{Int,Int}[]
        for bond in original_bonds
            x = bond[1]
            y = bond[2]
            if !in_group_a[x] && !in_group_a[y]
                new_x = remap[x]
                new_y = remap[y]
                push!(new_bonds_b, (new_x, new_y))
            end
        end

        piece_particles_b = []
        for i in new_particle_indices_b
            push!(piece_particles_b, particles[i])
        end
        cm_b, mass_b = calculate_center_of_mass(piece_particles_b)

        new_id = length(rigidbodies) + 1
        for i in new_particle_indices_b
            particles[i].rigidbody = new_id
        end

        new_rigidbody = rigidbody_struct(new_id, new_particle_indices_b, cm_b, rb.V, rb.ω, mass_b, new_bonds_b, rb.break_threshold)
        push!(rigidbodies, new_rigidbody)
    end
end