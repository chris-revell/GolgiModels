FROM julia:latest
COPY . /app
WORKDIR /app
CMD julia scripts/interactiveCatalystBothNonLinear.jl 
