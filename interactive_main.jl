using Plots, LinearAlgebra, StaticArrays, GLMakie

include("2D_physics_engine.jl")

Makie.inline!(false)

particles, liquid, gas, powder, solid, rigidbodies, softbodies = create_scene()
id_grid, cell_of_particle = init_grids(particles)

fig = Figure(size = (900, 600))

layout_botoes = fig[1, 1] = GridLayout(tellwidth = true, tellheight = false)
ax = Axis(fig[1, 2], title = "Ferramenta Atual: NADA", aspect = DataAspect())
deregister_interaction!(ax, :rectanglezoom)

GLMakie.xlims!(ax, 0, box_size_x)
GLMakie.ylims!(ax, 0, box_size_y)

material_grid_obs = Observable(build_material_grid(particles, id_grid))
colors = cgrad([:white, :brown, :blue, :gray, :orange], 5, categorical=true)
GLMakie.heatmap!(ax, material_grid_obs, colorrange = (0,4), colormap = colors)

btn1 = Button(layout_botoes[1, 1], label = "Liquid", buttoncolor = :blue, labelcolor = :white)
btn2 = Button(layout_botoes[2, 1], label = "Solid", buttoncolor = :brown, labelcolor = :white)
btn3 = Button(layout_botoes[3, 1], label = "Gas", buttoncolor = :gray, labelcolor = :white)
btn4 = Button(layout_botoes[4, 1], label = "Powder", buttoncolor = :orange, labelcolor = :white)

tipo_atual = Observable(0)

on(btn1.clicks) do _; tipo_atual[] = 1; ax.title = "Desenhando: Liquid"; end
on(btn2.clicks) do _; tipo_atual[] = 2; ax.title = "Desenhando: Solid"; end
on(btn3.clicks) do _; tipo_atual[] = 3; ax.title = "Desenhando: Gas"; end
on(btn4.clicks) do _; tipo_atual[] = 4; ax.title = "Desenhando: Powder"; end

material_names = Dict(1 => "liquid", 2 => "solid", 3 => "gas", 4 => "powder")

mouse_pressionado = Observable(false)

function adicionar_particula_se_ativo()
    if tipo_atual[] != 0
        posicao = mouseposition(ax.scene)

        if 0 <= posicao[1] <= box_size_x && 0 <= posicao[2] <= box_size_y
            px = Int(floor(posicao[1]/grid_size)) + 1
            py = Int(floor(posicao[2]/grid_size)) + 1
            if 1 <= px <= pixel_size_x && 1 <= py <= pixel_size_y
                spawn_particle!(particles, liquid, gas, powder, solid, id_grid, cell_of_particle,
                                 px, py, material_names[tipo_atual[]])
            end
        end
    end
end

on(ax.scene.events.mousebutton) do event
    if event.button == Mouse.left
        if event.action == Mouse.press
            mouse_pressionado[] = true
            adicionar_particula_se_ativo()
        elseif event.action == Mouse.release
            mouse_pressionado[] = false
        end
    end
end

on(ax.scene.events.mouseposition) do _
    if mouse_pressionado[]
        adicionar_particula_se_ativo()
    end
end

screen = display(fig)

t = 0.0
step = 0
@async begin
    while t < tmax && screen.window_open[]
        try
            simulation_step(particles, liquid, gas, powder, solid, rigidbodies, softbodies, id_grid, cell_of_particle)

            global step += 1
            if step % 3 == 0
                material_grid_obs[] = build_material_grid(particles, id_grid)
            end

            global t += dt
        catch e
            println("SIMULATION FAILED: ", e)
            break
        end
        sleep(0.001)
        yield()
    end
    println("loop ended, t=$t, step=$step")
end

if !isinteractive()
    wait(screen)
end