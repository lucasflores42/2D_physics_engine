
# função contração do espaço
function space_contraction(particles, p)

    gamma = 1 / (1 - v^2/c^2)^(1/2)
    Lo = p.radius   # change p.radius to a vector?
    L = Lo / gamma

    p.radius = L
end





