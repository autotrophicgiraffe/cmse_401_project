using NLsolve
using LinearAlgebra, SparseArrays, IterativeSolvers, Roots
using AlgebraicMultigrid
using DelimitedFiles
# using Plots
using GLMakie
using SpecialFunctions
using Printf

elapsed_time = @elapsed begin
# file provide the geometry information (level-set function)
# source function and dirichelet boundary conditions
include("disk.jl")

# source function for the cut cell solver
include("cut_cell_solver_explicit.jl")

# use the line by line strategy whenever possible, if true
line_by_line_opt = true
# check the distance between the two boundary points used
# in the RBF, if true
check_dist_opt = true#true #false #true
# string to distinguish different test
test_name = "_disk_heat"
save_dir  = "data/"

# computational domain
x_left = -1.1;
y_left = -1.1;
Lx = 2.2;
time_final = 10.2*0.5*0.1999562886971;
graphical_output = false
if false
    #for n_ref in (100)
    for n_ref in (100,200,400,800,1600)
        use_explicit = true 
        CFL = 0.7
        θ = 0.0
        maxerr, h, iter_num,u,err,x,y = cut_cell_solver(
            n_ref+1,n_ref+1,CFL,θ,x_left,Lx,y_left,time_final,
            line_by_line_opt,check_dist_opt,use_explicit,graphical_output)
        @printf " %8.3e & %8.3e & & " h maxerr 
    
        use_explicit = false 
        CFL = 2.0
        θ = 0.5
        maxerr, h, iter_num,u,err,x,y = cut_cell_solver(
            n_ref+1,n_ref+1,CFL,θ,x_left,Lx,y_left,time_final,
            line_by_line_opt,check_dist_opt,use_explicit,graphical_output)
        @printf "%8.3e & & %.2f " maxerr iter_num
        
        CFL = 2.0
        θ = 0.25
        maxerr, h, iter_num,u,err,x,y = cut_cell_solver(
            n_ref+1,n_ref+1,CFL,θ,x_left,Lx,y_left,time_final,
            line_by_line_opt,check_dist_opt,use_explicit,graphical_output)
        @printf "%8.3e & & %.2f " maxerr iter_num
        
        CFL = 0.85
        θ = 1.0/12.0
        maxerr, h, iter_num,u,err,x,y = cut_cell_solver(
            n_ref+1,n_ref+1,CFL,θ,x_left,Lx,y_left,time_final,
            line_by_line_opt,check_dist_opt,use_explicit,graphical_output)
        @printf "%8.3e & & %.2f \n" maxerr iter_num
    end
end

graphical_output = false
n_ref = 400
use_explicit = true 
CFL = 0.7
θ = 0.0
maxerr, h, iter_num,u,err,x,y = cut_cell_solver(
    n_ref+1,n_ref+1,CFL,θ,x_left,Lx,y_left,time_final,
    line_by_line_opt,check_dist_opt,use_explicit,graphical_output)

using CairoMakie

f = Figure(fontsize = 18,
           resolution = (600,500))
ax = Axis(f[1, 1],
     title = L"u(x,y,T\,)",
     xlabel = "x",
     ylabel = "y")
co = contourf!(x, y, u, levels = LinRange(-0.27,0.27,11), colormap=:jet)
contour!(co,x, y, u, levels = LinRange(-0.27,0.27,10), color = :black)
Colorbar(f[1, 2], co, width = 20)
colsize!(f.layout, 1, Aspect(1, 1.0))
limits!(ax, -1, 1, -1, 1) 
save("u_wave_exp.pdf", f)

f = Figure(fontsize = 18,
           resolution = (600,500))
ax = Axis(f[1, 1],
     title = L"e(x,y,T\,)",
     xlabel = "x",
     ylabel = "y")
co = contourf!(x, y, err, levels = 11, colormap=:jet)
contour!(co,x, y, err, levels = 10, color = :black)
Colorbar(f[1, 2], co, width = 20)
colsize!(f.layout, 1, Aspect(1, 1.0))
limits!(ax, -1, 1, -1, 1) 
save("err_wave_exp.pdf", f)


use_explicit = false 
CFL = 2.0
θ = 0.5
maxerr, h, iter_num,u,err,x,y = cut_cell_solver(
    n_ref+1,n_ref+1,CFL,θ,x_left,Lx,y_left,time_final,
        line_by_line_opt,check_dist_opt,use_explicit,graphical_output)

f = Figure(fontsize = 18,
           resolution = (600,500))
ax = Axis(f[1, 1],
     title = L"u(x,y,T\,)",
     xlabel = "x",
     ylabel = "y")
co = contourf!(x, y, u, levels = LinRange(-0.27,0.27,11), colormap=:jet)
contour!(co,x, y, u, levels = LinRange(-0.27,0.27,10), color = :black)
Colorbar(f[1, 2], co, width = 20)
colsize!(f.layout, 1, Aspect(1, 1.0))
limits!(ax, -1, 1, -1, 1) 
save("u_wave_1_2.pdf", f)

f = Figure(fontsize = 18,
           resolution = (600,500))
ax = Axis(f[1, 1],
     title = L"e(x,y,T\,)",
     xlabel = "x",
     ylabel = "y")
co = contourf!(x, y, err, levels = 11, colormap=:jet)
contour!(co,x, y, err, levels = 10, color = :black)
Colorbar(f[1, 2], co, width = 20)
colsize!(f.layout, 1, Aspect(1, 1.0))
limits!(ax, -1, 1, -1, 1) 
save("err_wave_1_2.pdf", f)

CFL = 2.0
θ = 0.25
maxerr, h, iter_num,u,err,x,y = cut_cell_solver(
    n_ref+1,n_ref+1,CFL,θ,x_left,Lx,y_left,time_final,
    line_by_line_opt,check_dist_opt,use_explicit,graphical_output)

f = Figure(fontsize = 18,
           resolution = (600,500))
ax = Axis(f[1, 1],
     title = L"u(x,y,T\,)",
     xlabel = "x",
     ylabel = "y")
co = contourf!(x, y, u, levels = LinRange(-0.27,0.27,11), colormap=:jet)
contour!(co,x, y, u, levels = LinRange(-0.27,0.27,10), color = :black)
Colorbar(f[1, 2], co, width = 20)
colsize!(f.layout, 1, Aspect(1, 1.0))
limits!(ax, -1, 1, -1, 1) 
save("u_wave_1_4.pdf", f)

f = Figure(fontsize = 18,
           resolution = (600,500))
ax = Axis(f[1, 1],
     title = L"e(x,y,T\,)",
     xlabel = "x",
     ylabel = "y")
co = contourf!(x, y, err, levels = 11, colormap=:jet)
contour!(co,x, y, err, levels = 10, color = :black)
Colorbar(f[1, 2], co, width = 20)
colsize!(f.layout, 1, Aspect(1, 1.0))
limits!(ax, -1, 1, -1, 1) 
save("err_wave_1_4.pdf", f)


CFL = 0.85
θ = 1.0/12.0
maxerr, h, iter_num,u,err,x,y = cut_cell_solver(
    n_ref+1,n_ref+1,CFL,θ,x_left,Lx,y_left,time_final,
    line_by_line_opt,check_dist_opt,use_explicit,graphical_output)

f = Figure(fontsize = 18,
           resolution = (600,500))
ax = Axis(f[1, 1],
     title = L"u(x,y,T\,)",
     xlabel = "x",
     ylabel = "y")
co = contourf!(x, y, u, levels = LinRange(-0.27,0.27,11), colormap=:jet)
contour!(co,x, y, u, levels = LinRange(-0.27,0.27,10), color = :black)
Colorbar(f[1, 2], co, width = 20)
colsize!(f.layout, 1, Aspect(1, 1.0))
limits!(ax, -1, 1, -1, 1) 
save("u_wave_1_12.pdf", f)

f = Figure(fontsize = 18,
           resolution = (600,500))
ax = Axis(f[1, 1],
     title = L"e(x,y,T\,)",
     xlabel = "x",
     ylabel = "y")
co = contourf!(x, y, err, levels = 11, colormap=:jet)
contour!(co,x, y, err, levels = 10, color = :black)
Colorbar(f[1, 2], co, width = 20)
colsize!(f.layout, 1, Aspect(1, 1.0))
limits!(ax, -1, 1, -1, 1) 
save("err_wave_1_12.pdf", f)
end
println("Elapsed time: $elapsed_time seconds")
