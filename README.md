# iceFEM.jl

Contains the Julia package developed during Lift-Off Fellowship at UoA. This package implements the coupled problem of wave-induced motion of ice-shelves and icebergs. 

## How to get it

To download the pacakge, clone the repository, change directory and start Julia using 

```shell
git clone git@github.com:Balaje/iceFEM.jl.git
cd iceFEM.jl
julia --project=.
```

## Run a simple example

```julia
using iceFEM

# Define the dimensional parameters of the ice and fluid
ρᵢ = 922.5; # Density of Ice
E = 2e9; ν = 0.33;
L = 40000; h = 200;
ρₒ = 1025; # Density of Ocean
xg = 0; k0 = 0; # Position of Grounding line and spring constant of bedrock (0 if not available)
H = 500; # Ocean Depth
g = 9.8; # Acceleration due to gravity
ice = Ice(ρᵢ, E, ν, L, h);
fluid = Fluid(ρₒ, k0, g, H, xg);
```

Now we solve the frequency domain problem

```julia
# Solve the frequency domain problem
ω = 2π/500;
Aₚ = g/(1im*ω);
sol1 = solve(ice, fluid, ω , FreeClamped(), FiniteDepth(4)); # Solve the finite depth problem using 4 modes
R1 = sol1.aₘ[1]/Aₚ;
sol2 = solve(ice, fluid, ω , FreeClamped(), ShallowWater()); # Solve the shallow water problem
R2 = sol2.a₀[1]/Aₚ;
Aₚ = g/(1im*ω);

# Show the reflection coeffcients
@show R1 abs(R1);
@show R2 abs(R2);
```

We can obtain the displacement profiles and plot them
``` julia
using Plots
x = 0:0.01:sol1.ndp.geo[1]; 
U1 = u₁(x, sol1);
U2 = u₁(x, sol2);
plt = plot(x, abs.(U₁));
plot!(plt, x, abs.(U₂))
```
Similarly we can obtain the slope/shear force/bending moment
``` julia
# Slope ...
DxU1 = ∂ₓu₁(x, sol1);
DxU2 = ∂ₓu₁(x, sol2);
# Bending Moment ...
Dx2U1 = ∂ₓ²u₁(x, sol1);
Dx2U2 = ∂ₓ²u₁(x, sol2);
# Shear Force ...
Dx3U1 = ∂ₓ³u₁(x, sol1);
Dx3U2 = ∂ₓ³u₁(x, sol2);
```

More examples will be added slowly
