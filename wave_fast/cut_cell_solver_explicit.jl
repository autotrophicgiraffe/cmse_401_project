using ForwardDiff

##########################
#
#NUMBER OF THREADS
num_threads = Threads.nthreads()
#
##########################

############################################################################
# If not provided, use auto-differentiation to compute the gradient of the
# level set function
############################################################################
# greadient of the level set function
level_gradient = x -> ForwardDiff.gradient(level, x);
############################################################################
# define radial basis function
############################################################################
function radial_basis(r)
    return r^3
end

function polynomial_tail(z)
    return [1;z[1];z[2]]
end
############################################################################
# correction to the diagonal element based on the RBF
# InterpMatrix : determine the coefficient for the RBF+polynomial
# z_gp             : ghost point
# z_interir    : 2 by 1 vector saving the interiro point
# z_tau1           : boundary point closest to the ghost point z_gp
# z_tau2           : boundary point closet to the interior point z_interior
# h                        : mesh size
# N_poly_basis : number of basis for the polynomial tail
############################################################################
function correct_bc_rbf!(InterpMatrix,z_gp,z_interior,z_tau1,z_tau2,
                         h,N_poly_basis)
    N_tau = 2
    update_RBF_interp_matrix!(InterpMatrix,z_interior,[z_tau1;z_tau2],N_tau)
    C_local = inv(InterpMatrix)
    # Compute the coefficients for the boundary point in the RBF interpolation,
    # which is [rbf_interoir(z_gp);rbf_bc[z_gp];polynomial tails]
    gp_rbf_dof = [radial_basis(h);
                  radial_basis(norm(z_tau1-z_gp));
                  radial_basis(norm(z_tau2-z_gp));
                  polynomial_tail(z_gp)]
    # compute the correction for the discrete Laplace operator
    coeff_gp = dot(C_local[:,1],gp_rbf_dof)
    if(coeff_gp>1.0)
        println("Diagonal dominance may be broken!")
        return
    end
    # compute the correction for the righthand side
    bc_values = [dirichlet_bc(z_tau1);dirichlet_bc(z_tau2)]
    bc_correction =  dot(gp_rbf_dof,C_local[:,2:3]*bc_values)
    return coeff_gp,bc_correction
end
##############################################################################
# find closeset point on the curve
##############################################################################
function find_closest_point(x0,y0)
    function f1!(F, z)
        grad1,grad2 = level_gradient(z);
        F[1] = grad2*(z[1]-x0)-grad1*(z[2]-y0)
        F[2] = level(z)
    end
    z_initial = [x0;y0]
    res = nlsolve(f1!,z_initial,ftol=1e-12)
    # solution to the equation, obtain the coordinate
    # of the point on the bounday
    return res.zero
end
##############################################################################
# find z_tau2 such that z0z_tau2 and z0z_tau1 is 45 degrees,
# and z_tau2 is on the curve with z0 = (x0,y0),
# with initial guess z_initial
##############################################################################
function find_45_degree_point(x0,y0,z_tau1,z_initial)
    function h1!(F, z)
        v11 = z[1]-x0
        v12 = z[2]-y0
        v21 = z_tau1[1]-x0
        v22 = z_tau1[2]-y0
        r1 = sqrt( v11^2+v12^2 )
        r2 = sqrt( v21^2+v22^2 )
        F[1] = v11*v21+v12*v22-r1*r2*sqrt(2)/2
        F[2] = level(z)
    end
    res = nlsolve(h1!,z_initial,ftol=1e-12)
    return res.zero
end


##############################################################################
# function to update the linear syste which determines the RBF interpolation
# InterpMatrix: save the matrix
# p_ij: interior point ij
# p_tau: boundary points
# N_tau: number of boundary points
##############################################################################
function update_RBF_interp_matrix!(InterpMatrix,p_ij,p_tau_array,N_tau,N_poly_basis=3)
    # RBF part
        Threads.@threads for j = 1:N_tau
            rbf_1j = radial_basis(norm(p_ij-p_tau_array[2*j-1:2*j]))
            InterpMatrix[1,j+1] = rbf_1j
            InterpMatrix[j+1,1] = rbf_1j
        end
        Threads.@threads for i = 1:N_tau-1
            for j = i+1:N_tau
                rbf_ij = radial_basis(norm(p_tau_array[2*j-1:2*j]-p_tau_array[2*i-1:2*i]))
                InterpMatrix[i+1,j+1] = rbf_ij
                InterpMatrix[j+1,i+1] = rbf_ij
            end
        end
    # polynomial part, linear polynomial
    ind = N_tau+2:N_tau+N_poly_basis+1
    values = polynomial_tail(p_ij)
    InterpMatrix[1,ind] = values[:]
    InterpMatrix[ind,1] = values[:]
        for i= 1:N_tau
	    values = polynomial_tail(p_tau_array[2*i-1:2*i])
            InterpMatrix[i+1,ind] = values[:]
            InterpMatrix[ind,i+1] = values[:]
        end
end

############################################################################
# solver function
# nx                       : number of grid points in the x direction
# ny                       : number of grid points in the y direction
# time_final       : Final time
# line_by_line_opt : true, use line by line strategy whenever possible
#                                        flase, use rbf everywhere
# check_dist_opt   : true, check the distance for the two boundary points used
#                                        in the rbf
#                                        false, not check
############################################################################
function cut_cell_solver(nx,ny,CFL,θ,x_start,Lx,y_start,time_final,
                         line_by_line_opt,check_dist_opt=true,
                         use_explicit=true,graphical_output=false)
    ############################################################
    # Step 1: Assemble the mask matrix to detect the geometry
    ############################################################
    # Assume a bounding box with lower left corner [0,0] and upper right corner [Lx,Ly]
    # Use a fixed h and make Ly = h*ny
    h =  Lx / (nx-1)
    Ly = Lx

    epsilon = 0.025
    # mask matrix saving the geometry info,
    # interior of the computational domain 1 (different from interior points and
    # interior ghost point)
    # otherwise 0
    mask = zeros(Int64,nx,ny)
    x = zeros(nx)
    y = zeros(ny)
    # Set up grids
    Threads.@threads for i = 1:nx
        x[i] = x_start+h*(i-1)
    end
    Threads.@threads for i = 1:ny
        y[i] = y_start+h*(i-1)
    end
    # Evaluate the levelset function to get mask
    lval = zeros(Float64,nx,ny)
    Threads.@threads for j = 1:ny
        for i = 1:nx
            l = level([x[i];y[j]])
            lval[i,j] = l
            if l < 0
                mask[i,j] = 1
            end
        end
    end
    # Get interior boundary points
    # maskin: 1 interior points, -1 interior ghost points
    # maskgp: 1 interior boundary points
    maskin = copy(mask)
    maskgp = zeros(Int64,nx,ny)
    Threads.@threads for j = 2:ny-1
        for i = 2:nx-1
            if mask[i,j]==1
                ms = mask[i+1,j] + mask[i-1,j] + mask[i,j+1] + mask[i,j-1]
                if ms < 4
                    maskin[i,j] = -1
                    maskgp[i,j] = 1
                end
            end
        end
    end
    # We now have three arrays:
    # mask = 1 where levelset is < 0
    # maskin = -1 on ghost points and 1 inside
    # maskgp = 1 on ghost points
    # Count number of degrees of freedom, i.e. interior points
    N = sum(mask-maskgp)
    # Count number of ghostpoints
    Ngp = sum(maskgp)
    # We create a N x 2 list with 2D index coordinates for each interior point.
    # We also create a 2D array that maps cartesian index to linear.
    cart2lin = zeros(Int64,nx,ny)
    lin2cart = zeros(Int64,2,N)
    # index for the interior points
    # cart2lin: >0, index for the interior point
    #                   =0, not interior points
    ilin = 1
    for j = 1:ny    
	for i = 1:nx
            if maskin[i,j] == 1
                lin2cart[1,ilin] = i
                lin2cart[2,ilin] = j
                cart2lin[i,j] = ilin
		ilin += 1
            end
        end
    end
    ############################################################
    # Step 2: Assemble the discrete Laplacian operator
    ############################################################
    # "Assemble" the Laplacian A
    # obtain the sparsity pattern for the interior part
    # no boundary point is assigned here
    # index for the array of rows, cols and values
    innz = 1
    rows = zeros(Int64,5*N)
    cols = zeros(Int64,5*N)
    vals = zeros(Float64,5*N)
    # sweep over all the interior points
    #Threads.@threads for k in 1:num_threads
    for ilin = 1:N
        i = lin2cart[1,ilin]
        j = lin2cart[2,ilin]
        # cart2lin == 0
        if cart2lin[i-1,j] != 0
            rows[innz] = ilin
            cols[innz] = cart2lin[i-1,j]
            vals[innz] = coefficientbeta(x[i]-0.5*h, y[j])
            innz +=1
        end
        if cart2lin[i+1,j] != 0
            rows[innz] = ilin
            cols[innz] = cart2lin[i+1,j]
            vals[innz] = coefficientbeta(x[i]+0.5*h, y[j])
            innz +=1
        end
        if cart2lin[i,j-1] != 0
            rows[innz] = ilin
            cols[innz] = cart2lin[i,j-1]
            vals[innz] = coefficientbeta(x[i], y[j]-0.5*h)
            innz +=1
        end
        if cart2lin[i,j+1] != 0
            rows[innz] = ilin
            cols[innz] = cart2lin[i,j+1]
            vals[innz] = coefficientbeta(x[i], y[j]+0.5*h)
            innz +=1
        end
        rows[innz] = ilin
        cols[innz] = ilin
        vals[innz] = -coefficientbeta(x[i]-0.5*h, y[j])-coefficientbeta(x[i]+0.5*h, y[j])-
            coefficientbeta(x[i], y[j]-0.5*h)-coefficientbeta(x[i], y[j]+0.5*h)
        innz +=1
    end
    #end
    innz -=1
    A = sparse(rows[1:innz],cols[1:innz],vals[1:innz])
    Nbddof = 0
    # How many DOFs has at least one gp as neighbour?
    # Loop to count
    Threads.@threads for ilin = 1:N
        i = lin2cart[1,ilin]
        j = lin2cart[2,ilin]
        if maskin[i-1,j]==-1 || maskin[i+1,j]==-1 || maskin[i,j-1]==-1 || maskin[i,j+1]==-1
            Nbddof += 1
	end
    end
    Nbddof = Nbddof[]   
    print("   ", Nbddof)

    # This array holds the distances from a boundary DOF gridpoint to
    # the boundary on the
    # left, right, below, above
    # the distance is given as a positive number >= 1
    # zero if there is no boundary
    bc_stuff = zeros(Float64,3,Nbddof) # distance from the boundary
    bd_idx       = zeros(Int64,Nbddof)     # index of the boundary dof (1D format)
    InterpMatrix = zeros(Float64,6,6)
    N_poly_basis = 3
    N_tau = 2
    ind = N_tau+2:N_tau+N_poly_basis+1

    # Define the right hand side
    bc_correction_array = zeros(N)
    # source for the manufactured solution with 0 boundary conditions
    j = 1
    for i = 1:N
        i1 = lin2cart[1,i]
        j1 = lin2cart[2,i]
        # interior point having a ghost neighbour
        if maskin[i1-1,j1]==-1 || maskin[i1+1,j1]==-1 || maskin[i1,j1-1]==-1 || maskin[i1,j1+1]==-1
            bd_idx[j] = i
            if maskin[i1-1,j1] == -1
                # line by line can not be used or not using that strategy at all
                if( ( (i1-2>0)&&(mask[i1-2,j1]==1) )||(!line_by_line_opt) )
                    ##########################################
                    # Ministep 1: locate the boundary point,
                    # used to interpolate the ghost point
                    ##########################################
                    x0 = x[i1-1]
                    y0 = y[j1]
                    # find point closet to x0,y0
                    z_tau1 = find_closest_point(x0,y0)
                    # find point closet to the interior point
                    z_tau2 = find_closest_point(x0+h,y0)
                    # if z_tau1 and z_tau2 are too close, find the boundary point z_tau2
                    # such that the angle between the direction z_tau2 z_gp and
                    # the direction z_tau1z_gp is 45 degrees
                    if( (norm(z_tau1-z_tau2)<epsilon*h)&&(check_dist_opt) )
                        z_initial = [x0-h;y0-h]
                        z_tau2 = find_45_degree_point(x0,y0,z_tau1,z_initial)
                    end
                    ####################################################
                    # Ministep 2: correct the diagonal term based on the RBF
                    ####################################################
                    coeff_gp,bc_correction = correct_bc_rbf!(InterpMatrix,[x0;y0],[x0+h;y0],z_tau1,z_tau2,h,N_poly_basis)
                    A[i,i] = A[i,i] + coeff_gp*coefficientbeta(x0+0.5*h,y0)
                    #b[i] = b[i] - bc_correction*coefficientbeta(x0+0.5*h,y0)
                    bc_correction_array[i] = bc_correction*coefficientbeta(x0+0.5*h,y0)
                else
                    x0 = x[i1-1]
                    y0 = y[j1]
                    f1(z) = level([z;y0])
                    # fzero, f1 target function, x0 initial guess
                    r = fzero(f1,x0,order=1)
                    gamma = (x0+h-r)/h
                    A[i,i] = A[i,i] + (gamma-1.0)/gamma*coefficientbeta(x0+0.5*h,y0)
                    #b[i] = b[i] - dirichlet_bc([r;y0])/gamma*coefficientbeta(x0+0.5*h,y0)
                    bc_correction_array[i] = dirichlet_bc([r;y0])/gamma*coefficientbeta(x0+0.5*h,y0)
                end
            end
            if maskin[i1+1,j1] == -1
                # line by line can not be used or not using that strategy at all
                if( ( (i1+2<nx+1)&&(mask[i1+2,j1]==1) )||(!line_by_line_opt) )
                    ##########################################
                    # Ministep 1: locate the boundary point,
                    # used to interpolate the ghost point
                    ##########################################
                    x0 = x[i1+1]
                    y0 = y[j1]
                    # find the closet point to the boundary ghost point
                    z_tau1 = find_closest_point(x0,y0)
                    # find point closet to the interior point
                    z_tau2 = find_closest_point(x0-h,y0)
                    # if the distance of z_tau1 and z_tau2 are two close, find a new
                    # z_tau2
                    if( (norm(z_tau1-z_tau2)<epsilon*h)&&(check_dist_opt) )
                        z_initial = [x0+h;y0+h]
                        z_tau2 = find_45_degree_point(x0,y0,z_tau1,z_initial)
                    end
                    ####################################################
                    # Ministep 2: correct the diagonal term based on the RBF
                    ####################################################
                    coeff_gp,bc_correction = correct_bc_rbf!(InterpMatrix,[x0;y0],[x0-h;y0],z_tau1,z_tau2,h,N_poly_basis)
                    A[i,i] = A[i,i] + coeff_gp*coefficientbeta(x0-0.5*h,y0)
                    #b[i] = b[i] - bc_correction*coefficientbeta(x0-0.5*h,y0)
                    bc_correction_array[i] = bc_correction*coefficientbeta(x0-0.5*h,y0)
                else
                    x0 = x[i1+1]
                    y0 = y[j1]
                    f2(z) = level([z;y0])
                    r = fzero(f2,x0,order=1)
                    gamma = -(x0-h-r)/h
                    A[i,i] = A[i,i] + (gamma-1.0)/gamma*coefficientbeta(x0-0.5*h,y0)
                    #b[i] = b[i] - dirichlet_bc([r;y0])/gamma*coefficientbeta(x0-0.5*h,y0)
                    bc_correction_array[i] = dirichlet_bc([r;y0])/gamma*coefficientbeta(x0-0.5*h,y0)
                end
            end
            if maskin[i1,j1-1] == -1
                # line by line can not be used or not using that strategy at all
                if( ( (j1-2>0)&&(mask[i1,j1-2]==1) )||(!line_by_line_opt) )
                    ##########################################
                    # Ministep 1: locate the boundary point,
                    # used to interpolate the ghost point
                    ##########################################
                    x0 = x[i1]
                    y0 = y[j1-1]
                    # find the closet point to the boundary ghost point
                    z_tau1 = find_closest_point(x0,y0)
                    # find the closet point to the boundary ghost point
                    z_tau2 = find_closest_point(x0,y0+h)
                    # if the distance of z_tau1 and z_tau2 are two close, find a new
                    # z_tau2
                    if( (norm(z_tau1-z_tau2)<epsilon*h)&&(check_dist_opt) )
                        z_initial = [x0-h;y0-h]
                        z_tau2 = find_45_degree_point(x0,y0,z_tau1,z_initial)
                    end
                    ####################################################
                    # Ministep 2: correct the diagonal term based on the RBF
                    ####################################################
                    coeff_gp,bc_correction = correct_bc_rbf!(InterpMatrix,[x0;y0],[x0;y0+h],z_tau1,z_tau2,h,N_poly_basis)
                    A[i,i] = A[i,i] + coeff_gp*coefficientbeta(x0,y0+0.5*h)
                    #b[i] = b[i] - bc_correction*coefficientbeta(x0,y0+0.5*h)
                    bc_correction_array[i] = bc_correction*coefficientbeta(x0,y0+0.5*h)
                else
                    x0 = x[i1]
                    y0 = y[j1-1]
                    f3(z) = level([x0;z])
                    r = fzero(f3,y0,order=1)
                    gamma = (y0+h-r)/h
                    A[i,i] = A[i,i] + (gamma-1.0)/gamma*coefficientbeta(x0,y0+0.5*h)
                    #b[i] = b[i] - dirichlet_bc([x0;r])/gamma*coefficientbeta(x0,y0+0.5*h)
                    bc_correction_array[i] = dirichlet_bc([x0;r])/gamma*coefficientbeta(x0,y0+0.5*h)
                end
            end
            if maskin[i1,j1+1] == -1
                # line by line can not be used or not using that strategy at all
                if( ( (j1+2<ny+1)&&(mask[i1,j1+2]==1) )||(!line_by_line_opt) )
                    ##########################################
                    # Ministep 1: locate the boundary point,
                    # used to interpolate the ghost point
                    ##########################################
                    # define the nonlienar equation to find
                    # the point on the physical boundary
                    x0 = x[i1]
                    y0 = y[j1+1]
                    # find the closet point to the boundary ghost point
                    z_tau1 = find_closest_point(x0,y0)
                    # find the closet point to the boundary ghost point
                    z_tau2 = find_closest_point(x0,y0-h)
                    # if the distance of z_tau1 and z_tau2 are two close, find a new
                    # z_tau2
                    if( (norm(z_tau1-z_tau2)<epsilon*h)&&(check_dist_opt) )
                        z_initial = [x0+h;y0+h]
                        z_tau2 .= find_45_degree_point(x0,y0,z_tau1,z_initial)
                    end
                    ##########################################################
                    # Ministep 2: correct the diagonal term based on the RBF
                    ##########################################################
                    coeff_gp,bc_correction = correct_bc_rbf!(InterpMatrix,[x0;y0],[x0;y0-h],z_tau1,z_tau2,h,N_poly_basis)
                    A[i,i] = A[i,i] + coeff_gp*coefficientbeta(x0,y0-0.5*h)
                    #b[i] = b[i] - bc_correction*coefficientbeta(x0,y0-0.5*h)
                    bc_correction_array[i] = bc_correction*coefficientbeta(x0,y0-0.5*h)
                else
                    x0 = x[i1]
                    y0 = y[j1+1]
                    f4(z) = level([x0;z])
                    r = fzero(f4,y0,order=1)
                    gamma = -(y0-h-r)/h
                    A[i,i] = A[i,i] + (gamma-1.0)/gamma*coefficientbeta(x0,y0-0.5*h)
                    #b[i] = b[i] - dirichlet_bc([x0;r])/gamma*coefficientbeta(x0,y0-0.5*h)
                    bc_correction_array[i] = dirichlet_bc([x0;r])/gamma*coefficientbeta(x0,y0-0.5*h)
                end
            end
            j += 1
        end
    end
    ############################################################
    # Step 3: Set up the time integrator
    ############################################################
    ##################################
    # Step 3a: set up time step size
    ##################################
    dt = CFL*h
    N_time_step = Int(ceil(time_final/dt))
    if(N_time_step>0)
        dt = time_final/N_time_step
    end
    ##################################
    # Step 3b: initial conditions
    ##################################
    cm = zeros(N)
    c = zeros(N)
    cp = zeros(N)
    err_lin = zeros(N)
    #for i = 1:N
    Threads.@threads for i = 1:N
        i1 = lin2cart[1,i]
        j1 = lin2cart[2,i]
        c[i] = initial_condition([x[i1];y[j1]])
        cm[i] = exact_solution([x[i1];y[j1]],-dt)
        #cm[i] = c[i]
    end
    ##################################
    # Step 3c: set up the operators
    ##################################

    #θ = 1.0/2.0
    Laplace =1.0/h^2*A
    Id = sparse(I,N,N)
    PoissonOperator = Id-θ*dt^2/h^2*A
    CurrentStepOperator = 2.0*Id+(1.0-2.0*θ)*dt^2/h^2*A
    ml = ruge_stuben(PoissonOperator,strength = Classical(1.0))
    # set up preconditioner
    p = aspreconditioner(ml);
    # Correct initial data
    #cp = Laplace*c
    #cm += dt^2/2.0*cp
    #if use_explicit
    #else
        # Fourth order accurate for u_t = 0 and zero Dirichlet
    #    cp = Laplace*Laplace*c
    #    cm += dt^4/24.0*cp
    #end
    
    ##################################
    # time marching
    ##################################
    u = zeros(nx,ny)
    U = Observable(u)
    err = zeros(nx,ny)
    ERR = Observable(err)

    function explicit_step!(it)
        tt = (it-1.0)*dt
        
        Threads.@threads for i = 1:N
            i1 = lin2cart[1,i]
            j1 = lin2cart[2,i]
            cp[i] = 0*dt^2*source_function([x[i1];y[j1]],tt)
        end
        
        # Should contain BC...?
        cp += dt^2*Laplace*c
        cp += 2.0*c
        cp -= cm
        # Swap
        cm = copy(c)
        c = copy(cp)

        fill!(u,NaN)
        fill!(err,NaN)
       
        Threads.@threads for i = 1:N
            i1 = lin2cart[1,i]
            j1 = lin2cart[2,i]
            #if x[i1] > 0 || y[j1] > 0
            u[i1,j1] = c[i]
            err[i1,j1] = c[i]-exact_solution([x[i1];y[j1]],it*dt)
            #end
        end
        err = err'
        ERR[] = err
        u = u'
        U[] = u
    end

    function implicit_step!(it)
        ##################################
        # time marching
        ##################################
        b = zeros(N)
        tt = (it-1)*dt
        b = CurrentStepOperator*c
        b -= PoissonOperator*cm
        log = gmres!(cp,PoissonOperator, b, Pl = p, log=true, reltol=1e-12)
        iter_num = log[2].iters
        cm = copy(c)
        c = copy(cp)
        fill!(u,NaN)
        fill!(err,NaN)
        Threads.@threads for i = 1:N
	    i1 = lin2cart[1,i]
            j1 = lin2cart[2,i]
            #if x[i1] > 0 || y[j1] > 0
            u[i1,j1] = c[i]
            err[i1,j1] = c[i]-exact_solution([x[i1];y[j1]],it*dt)
            #end
        end
        err = err'
        ERR[] = err
        u = u'
        U[] = u
        return iter_num
    end
    
    iter_num = 0        
    # Setup figure
    fig = Figure(resolution = (1400, 500))
    ax, sf = heatmap(fig[1,1], x,y,U,colorrange=(-0.5, 0.5),colormap=:jet)
    Colorbar(fig[1,2],sf)
    ax,sf=heatmap(fig[1,3], x,y,ERR,colorrange=(-0.001,0.001),colormap=:jet)
    Colorbar(fig[1,4],sf)
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    colsize!(fig.layout, 3, Aspect(1, 1.0))
    if graphical_output
        record(fig, "movie.mp4", profile = "high", compression = 1) do io
            for it = 1:N_time_step
                if use_explicit
                    explicit_step!(it)
                else
                    implicit_step!(it)
                end
                recordframe!(io)
            end
        end
    else
        for it = 1:N_time_step
            if use_explicit
                explicit_step!(it)
            else
                iter_num += implicit_step!(it)
            end
        end
    end
    iter_num /=N_time_step
    
    #############################################################
    # present the results
    #############################################################
    u_diff = zeros(nx,ny)
    #for i = 1:N
    Threads.@threads for i = 1:N
	i1 = lin2cart[1,i]
        j1 = lin2cart[2,i]
        u[i1,j1] = c[i]
        u_diff[i1,j1] = c[i]-exact_solution([x[i1];y[j1]],time_final)
    end
    err = maximum(abs.(u_diff))
    fill!(u_diff,NaN)
    Threads.@threads for i = 1:N
    #for i = 1:N
        i1 = lin2cart[1,i]
        j1 = lin2cart[2,i]
        u_diff[i1,j1] = c[i]-exact_solution([x[i1];y[j1]],time_final)
    end
    return err, h, iter_num, u, u_diff, x, y

end
