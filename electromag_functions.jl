
# definir corrente ou estado das particulas como eletrificado ou n

# corrente induzida por campo magnetico
function faraday_law(particles, p)

    flux = integral B. dA
    epsilon = - flux
    current = epsilon / resistence 
end

# campo magnetico induzido por corrente
function ampere_law(particles, p)

    integral B. dl = mu_o * current
end


