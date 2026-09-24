# -----------------------------------------------------------------------------
#                           Create Scene
# ----------------------------------------------------------------------------- 
function create_scene1()

    particles = Union{liquid_struct, solid_struct, gas_struct, powder_struct}[]
    liquid = liquid_struct[]
    liquid2 = liquid_struct[]
    solid = solid_struct[]
    gas = gas_struct[]
    powder = powder_struct[]
    rigidbodies = rigidbody_struct[]
    softbodies = softbody_struct[]

    # boundary of world
    for i in 1:pixel_size_x
        for j in 1:pixel_size_y

            if (i >= 1 && i <=5) || (i>=pixel_size_x-4 && i<=pixel_size_x) || (j>=1 && j<=5) || (j>=pixel_size_y-4 && j<=pixel_size_y)

                pos_x = (i - 1) * grid_size + grid_size/2
                pos_y = (j - 1) * grid_size + grid_size/2

                p = solid_struct(
                    length(particles) + 1,
                    [pos_x, pos_y],           
                    [0.0, 0.0],     # velocity
                    [0.0, 0.0],     # acceleration
                    grid_size/2,    # radius
                    10000.0,      # mass

                    0,              # rigidbody
                    0,

                    0,              # active
                    1,              # collision
                    0,              # gravity

                    "solid"         # material
                )
                push!(solid, p)
                push!(particles, p)

            end

            u_left = 300
            u_right = 350
            u_bottom = 50
            u_top = 150
            u_thickness = 3
            
            if (i >= u_left && i <= u_left + u_thickness && j >= u_bottom && j <= u_top) ||      # Left wall
            (i >= u_right - u_thickness && i <= u_right && j >= u_bottom && j <= u_top) ||       # Right wall
            (i >= u_left && i <= u_right && j >= u_bottom && j <= u_bottom + u_thickness)      # Bottom wall
            #(i >= u_left && i <= u_right && j >= u_top && j <= u_top + u_thickness)              # top wall

                
                pos_x = (i - 1) * grid_size + grid_size/2
                pos_y = (j - 1) * grid_size + grid_size/2

                p = solid_struct(
                    length(particles) + 1,
                    [pos_x, pos_y],           
                    [0.0, 0.0],     # velocity
                    [0.0, 0.0],     # acceleration
                    grid_size/2,    # radius
                    1000000.0,      # mass

                    0,              # rigidbody
                    0,

                    0,              # active
                    1,              # collision
                    0,              # gravity

                    "solid"         # material
                )
                push!(solid, p)
                push!(particles, p)
            end
        end
    end

    # some liquid
    for i in 1:10

        x = 325
        y = 55 + 10*rand()

        p = liquid_struct(length(particles)+1, SVector(x,y), @SVector(zeros(2)), @SVector(zeros(2)),
                       grid_size/2, 0.1, 0, 0,
                       0.4, 0.0, 0.4, 0.1, 0.1,   # density, pressure, target_density, stiff_coef, viscosity_coef
                       1, 1, 1, 1, 
                       1, "liquid")
        push!(liquid, p)
        push!(particles, p)
    end
    for i in 1:0

        x = 325
        y = 55 + 10*rand()

        p = liquid_struct(length(particles)+1, SVector(x,y), @SVector(zeros(2)), @SVector(zeros(2)),
                       grid_size/2, 0.1, 0, 0,
                       0.2, 0.0, 0.2, 0.1, 0.1,   # density, pressure, target_density, stiff_coef, viscosity_coef
                       1, 1, 1, 1, 
                       1, "liquid")
        push!(liquid2, p)
        push!(particles, p)
    end

    # some gas
    for i in 1:0

        x = 325
        y = 55 + 80*rand()

        p = gas_struct(length(particles)+1, SVector(x,y), SVector(rand(),0.0), @SVector(zeros(2)),
                        grid_size/2, 0.1, 
                        0, 0, 
                        1, 1, 1, 1, 
                        0, 300, "gas")
        push!(gas, p)
        push!(particles, p)
    end

    # some powder
    for i in 1:10

        x = 325
        y = 65 + 80*rand()

        p = powder_struct(
            length(particles) + 1,
            [x,y],           
            [0.0, 0.0],     # velocity
            [0.0, 0.0],     # acceleration
            grid_size/2,    # radius
            10.0,            # mass

            0,              # rigidbody
            0,

            1,              # active 
            1,              # collision
            1,              # gravity

            "powder"        
        )
        push!(powder, p)
        push!(particles, p)
    end

    create_cube!(particles, rigidbodies, length(rigidbodies)+1, [100.0, 8.0], [0.0, 0.0], [0.0],15, 3)
    create_cube!(particles, rigidbodies, length(rigidbodies)+1, [100-6, 25.0], [0.0, 0.0], [0.0],2, 27)
    create_cube!(particles, rigidbodies, length(rigidbodies)+1, [100+12, 8.0], [0.0, 0.0], [0.0],15, 3)

    create_sphere!(particles, rigidbodies, length(rigidbodies)+1, [100.0, 200.0], [0.0, 0.0], [0.0], 10)
    #create_sphere!(particles, rigidbodies, length(rigidbodies)+1, [328.0, 200.0], [0.0, 0.0], [0.0], 10)
    
    create_rope!(particles, softbodies, length(softbodies)+1, [250.0, 180.0], 15, 0.1, grid_size)
    create_rope2!(particles, softbodies, length(softbodies)+1, [200.0, 180.0], 15, 0.1, grid_size)

    return particles, liquid, liquid2, gas, powder, solid, rigidbodies, softbodies
end

function create_scene2()

    particles = Union{liquid_struct, solid_struct, gas_struct, powder_struct}[]
    liquid = liquid_struct[]
    liquid2 = liquid_struct[]
    solid = solid_struct[]
    gas = gas_struct[]
    powder = powder_struct[]
    rigidbodies = rigidbody_struct[]
    softbodies = softbody_struct[]

    # boundary of world
    for i in 1:pixel_size_x
        for j in 1:pixel_size_y

            if (i >= 1 && i <=5) || (i>=pixel_size_x-4 && i<=pixel_size_x) || (j>=1 && j<=5) || (j>=pixel_size_y-4 && j<=pixel_size_y)

                pos_x = (i - 1) * grid_size + grid_size/2
                pos_y = (j - 1) * grid_size + grid_size/2

                p = solid_struct(length(particles) + 1, [pos_x, pos_y], [0.0, 0.0], [0.0, 0.0], 
                        grid_size/2, 10000.0, 0, 0,
                        0, 1, 0,   
                        "solid" )
                push!(solid, p)
                push!(particles, p)

            end

            u_left = 230
            u_right = 270
            u_bottom = 50
            u_top = 250
            u_thickness = 3
            
            if (i >= u_left && i <= u_left + u_thickness && j >= u_bottom && j <= u_top) ||      # Left wall
            (i >= u_right - u_thickness && i <= u_right && j >= u_bottom && j <= u_top) ||       # Right wall
            (i >= u_left && i <= u_right && j >= u_bottom && j <= u_bottom + u_thickness)      # Bottom wall
            #(i >= u_left && i <= u_right && j >= u_top && j <= u_top + u_thickness)              # top wall

                
                pos_x = (i - 1) * grid_size + grid_size/2
                pos_y = (j - 1) * grid_size + grid_size/2

                p = solid_struct(
                    length(particles) + 1, [pos_x, pos_y], [0.0, 0.0],  [0.0, 0.0],     
                            grid_size/2, 1000000.0, 0, 0,
                            0, 1, 0,   
                            "solid")
                push!(solid, p)
                push!(particles, p)
            end
        end
    end

    # some liquid
    for i in 1:500

        x = 230 + rand()*40
        y = 50 + rand()*200

        p = liquid_struct(length(particles)+1, SVector(x,y), @SVector(zeros(2)), @SVector(zeros(2)),
                       grid_size/2, 0.1, 0, 0,
                       0.4, 0.0, 0.4, 0.1, 0.1,   # density, pressure, target_density, stiff_coef, viscosity_coef
                       1, 1, 1, 1, 
                       1, "liquid")
        push!(liquid, p)
        push!(particles, p)
    end
    

    return particles, liquid, liquid2, gas, powder, solid, rigidbodies, softbodies
end