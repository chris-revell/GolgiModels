FROM julia:1.9.2-bullseye

WORKDIR /app

EXPOSE 9384

# Copy files
COPY Project.toml .

# Install deps
RUN julia --project=. -e "import Pkg; Pkg.instantiate(); Pkg.precompile()"

# Run apps
COPY scripts ./scripts
COPY src ./src
COPY supplementary ./supplementary

# TODO: run the script and create a sysimage to speed up load time of the container

CMD julia --project=. -i -e 'include("scripts/interactiveToolWithDwellTimes.jl")'
