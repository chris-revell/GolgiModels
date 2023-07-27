FROM julia:1.9.2-bullseye

WORKDIR /app/GolgiModels

EXPOSE 9384

# Copy files
COPY Project.toml ./Project.toml
COPY src ./src
COPY supplementary/GolgiCompartmentModel_reduced.png ./supplementary/
COPY supplementary/GolgiCompartmentModel.png ./supplementary/

# Run on container build: Install deps
RUN julia -e 'using Pkg; Pkg.develop(path="./")'
RUN julia --project=. -i -e "using Pkg; Pkg.instantiate(); using GolgiModels; golgiApp()"

# Run on container startup
CMD julia --project=. -i -e "using GolgiModels; golgiApp()"


# FROM julia:1.9.2-bullseye

# WORKDIR /app/GolgiModels

# EXPOSE 9384

# # Copy files
# # COPY dummyProject.toml ./Project.toml
# COPY Project.toml ./Project.toml
# COPY src ./src
# COPY supplementary/GolgiCompartmentModel_reduced.png ./supplementary/
# COPY supplementary/GolgiCompartmentModel.png ./supplementary/

# # Run on container build: Install deps
# RUN julia -e 'using Pkg; Pkg.develop(path="./")'
# # RUN julia --project=. -e 'using Pkg;Pkg.add("Catalyst");Pkg.add("DifferentialEquations");Pkg.add("DrWatson");Pkg.add("FileIO");Pkg.add("Symbolics");'
# # RUN julia --project=. -e 'using Pkg;Pkg.add("Format");Pkg.add("FromFile");Pkg.add("GeometryBasics");Pkg.add("Images");Pkg.add("JSServe");Pkg.add("UnPack");'
# # RUN julia --project=. -e 'using Pkg;Pkg.add("LinearAlgebra");Pkg.add("Makie");Pkg.add("OrdinaryDiffEq");Pkg.add("PrecompileTools");Pkg.add("SparseArrays");Pkg.add("WGLMakie")'

# # Pkg.instantiate(); using GolgiModels"

# # # Run on container startup
# CMD julia --project=. -i -e "using Pkg; Pkg.instantiate(); using GolgiModels; golgiApp()"
