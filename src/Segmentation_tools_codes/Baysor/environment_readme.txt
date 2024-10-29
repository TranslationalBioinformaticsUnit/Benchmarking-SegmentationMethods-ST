# INSTALLING BAYSOR in cmd

# Load Julia (install Julia if not present in your system)
module load julia/1.8.4 # type: ignore

# Install Baysor using Julia
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git")); Pkg.build()' # type: ignore