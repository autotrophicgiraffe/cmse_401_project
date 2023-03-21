# cmse_401_project
Spring Semester 2023 Project Repository for CMSE 401 Parallel Programming

To run the scripts, the system you are using needs to have Julia installed and added to the PATH.

#########################################

#		INSTALL JULIA		#

#########################################

On Linux/Ubuntu:
-   Download latest version on Julia Website (x86 for me)
-   Unzip with tar. tar -xvf julia-1.8.5-linux-x86_64.tar.gz
-   Add to path by editing following file: nano ~/.bashrc
-   Add to bottom of file: export PATH="$PATH:/location/julia-version/bin" 
-   Run file again: source ~/.bashrc

This enables you to utilize Julia on your system.

To run a julia script:
-   julia name_of_script.jl

For HPCC
-   module load julia/1.5.2 (or most recent version)

#########################################

#		ADD PACKAGES		#

#########################################

For most of the files in the repo, the following commands will need to be run to install packages
-   julia
-   import Pkg 
-   Pkg.add("IterativeSolvers")
-   Pkg.add("Plots")
-   Pkg.add("Roots")
-   Pkg.add("AlgebraicMultigrid")
-   Pkg.add("Arpack")
-   Pkg.add("NearestNeighbors")
-   Pkg.add("PyPlot")

If there are still packages that are not installed, install using same format as above.

Test if you can run the code by navigating to the code folder within the repo and running:
-   julia amg_eb.jl
May be warnings about assignments... go ahead and ignore these.
Should run and produce several .txt files. If no h.txt file exists then something has gone wrong.

#########################################

#		RESOURCES		#

#########################################

Installing Julia on Ubuntu:

https://www.digitalocean.com/community/tutorials/how-to-install-julia-programming-language-on-ubuntu-22-04


