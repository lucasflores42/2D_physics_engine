#import Pkg
#Pkg.add(["StaticArrays", "Plots", "LinearAlgebra", "GLMakie"])
using Plots, LinearAlgebra, StaticArrays #, GLMakie

# -----------------------------------------------------------------------------
#                           Parameters
# ----------------------------------------------------------------------------- 
# 270 height and 480 width, total 129,600 pixels.
const grid_size = 1.0
const pixel_size_x = 480
const pixel_size_y = 270 
const box_size_x = pixel_size_x * grid_size
const box_size_y = pixel_size_y * grid_size

const tmax = 1000.0
const dt = 0.01

include("sph_functions.jl")
include("rigidbody_functions.jl")
include("softbody_functions.jl")
include("particle_functions.jl")
include("collision_functions.jl")
include("other_functions.jl")
include("material_functions.jl")
include("scenes_functions.jl")
include("electromag_functions.jl")
include("quantum_functions.jl")
include("thermo_functions.jl")
include("relativity_functions.jl")


# -----------------------------------------------------------------------------
#                           Simulation step
# ----------------------------------------------------------------------------- 
function simulation_step(particles, liquid, liquid2, gas, powder, solid, rigidbodies, softbodies, id_grid, cell_of_particle)

    particle_physics(particles, liquid, liquid2, gas, powder, solid, id_grid, cell_of_particle)
    rigidbody_physics(particles, rigidbodies)
    softbody_physics(particles, softbodies)

    collision_physics!(particles, rigidbodies, powder, liquid, gas, id_grid, cell_of_particle)
    #clamp_particles(particles)

    update_grids!(particles, id_grid, cell_of_particle)
end

# -----------------------------------------------------------------------------
#                           Visualization
# ----------------------------------------------------------------------------- 
function visualization(particles, id_grid, step)

    material_grid = build_material_grid(particles, id_grid)

    colors = cgrad([:white, :brown, :blue, :green, :gray, :orange], 6, categorical=true)

    plt = heatmap(material_grid', color=colors, clims=(0,5),
                  xlim=(0, box_size_x), ylim=(0, box_size_y),
                  title="Time $(round(step, digits=2))s",
                  xlabel="X", ylabel="Y",
                  size=(1920, 1080), aspect_ratio=:equal, legend=false)

    return plt
end

# -----------------------------------------------------------------------------
#                           Main Simulation
# ----------------------------------------------------------------------------- 
function main()
    t = 0.0
    step = 0

    particles, liquid, liquid2, gas, powder, solid, rigidbodies, softbodies = create_scene2()
    id_grid, cell_of_particle = init_grids(particles)

    while t < tmax
        step_time = @elapsed simulation_step(
            particles, liquid, liquid2, gas, powder, solid,
            rigidbodies, softbodies, id_grid, cell_of_particle
        )

        step += 1

        render_time = 0.0
        if step % 10 == 0
            render_time = @elapsed begin
                plt = visualization(particles, id_grid, t)
                display(plt)
            end
        end

        t += dt

        println("t = $(round(t, digits=2))s | step: $(round(step_time*1000, digits=2))ms | render: $(round(render_time*1000, digits=2))ms")
    end
end

main()
