FROM julia:1.9.2-bullseye

WORKDIR /app/GolgiModels

EXPOSE 9384

# Copy files
COPY Project.toml ./
COPY scripts ./scripts
COPY src ./src
COPY supplementary ./supplementary

# Run on container build: Install deps
RUN julia -e 'using Pkg; Pkg.develop(path="./")'
RUN julia --project=. -e "import Pkg; Pkg.instantiate(); using GolgiModels"

# Run on container startup
CMD julia --project=. -i -e "using GolgiModels; golgiApp()"
