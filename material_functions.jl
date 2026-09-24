function material_transformations(p1, p2, particles, gas, id_grid, cell_of_particle)
   
    if p1.material == "powder" && p2.material == "liquid"

        new_gas = gas_struct(length(particles)+1, p1.position, SVector(rand(),0.0), @SVector(zeros(2)),
                        grid_size/2, 0.1, 
                        0, 0, 
                        1, 1, 1, 1, 
                        0, 300, "gas")
        transform_particle!(particles, gas, id_grid, cell_of_particle, p1, p2, new_gas)
        return
    elseif p1.material == "liquid" && p2.material == "powder"

        new_gas = gas_struct(length(particles)+1, p2.position, SVector(rand(),0.0), @SVector(zeros(2)),
                        grid_size/2, 0.1, 
                        0, 0, 
                        1, 1, 1, 1, 
                        0, 300, "gas")
        transform_particle!(particles, gas, id_grid, cell_of_particle, p1, p2, new_gas)
        return
    end
end

function transform_particle!(particles, target_array, id_grid, cell_of_particle, p, p2, new_particle)

    if p.active == 0 || p2.active == 0
        return   # already transformed earlier this same scan
    end

    px = Int(floor(p.position[1] / grid_size)) + 1
    py = Int(floor(p.position[2] / grid_size)) + 1

    # remove p from its grid cell 
    erase_particle!(p, id_grid, cell_of_particle)

    # remove p2 from its grid cell 
    erase_particle!(p2, id_grid, cell_of_particle)

    # spawn the new gas particle
    push!(target_array, new_particle)
    push!(particles, new_particle)
    push!(cell_of_particle, (px, py))

    if !haskey(id_grid, (px, py))
        id_grid[(px, py)] = Int[]
    end
    push!(id_grid[(px, py)], new_particle.id)
end