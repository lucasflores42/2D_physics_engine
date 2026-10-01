# ------------------------------------------------------------------------
#                        World     
# ------------------------------------------------------------------------
# 270 height and 480 width, total 129,600 pixels.
const grid_size = 1.0
const pixel_size_x = 480
const pixel_size_y = 270 
const box_size_x = pixel_size_x * grid_size
const box_size_y = pixel_size_y * grid_size

const tmax = 1000.0
const dt = 0.01

# ------------------------------------------------------------------------
#                        Collisions     
# ------------------------------------------------------------------------
const restitution_x = 0.5
const restitution_y = 0.3
const restitution_angular = 0.3
const collision_min_distance = grid_size #* sqrt(2)
const max_velocity = 50.0
const max_angular_velocity = 20.0
const friction_coef = 0.3

# ------------------------------------------------------------------------
#                        Particles     
# ------------------------------------------------------------------------
const gravity_coef = 0.1
const colision_restitution_coefficient = 0.5
const collision_min_distance = grid_size #* sqrt(2)

# ------------------------------------------------------------------------
#                        Relativity     
# ------------------------------------------------------------------------
const c = 3 * 10^8 

# ------------------------------------------------------------------------
#                        Fluid dynamics     
# ------------------------------------------------------------------------
const smoothing_length = 0.2
const surface_tension = 0.15 
const sph_cell_range = Int(ceil(3 * smoothing_length / grid_size))