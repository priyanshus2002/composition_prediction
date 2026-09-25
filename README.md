# Optimal Working Fluid Composition for Internal Cascade Heat Pumps

This repository contains the optimization code used in for finding optimal working fluid compositions for internal and external cascade heat pump configurations. The optimizer minimizes total exergy destruction across all heat exchangers.


## Structure

```
.
├── load_model.jl           # Single entry point — include this in any script
├── trial.jl                # Example usage script
└── files/
    ├── stream.jl           # Properties and create_stream structs; thermodynamic variable definitions
    ├── utils.jl            # Constraint-building utilities (separator, throttle, exergy, etc.)
    ├── mixture_examples.jl # Library of 44 named multicomponent mixtures
    ├── configurations.jl   # Heat pump configuration structs and build! methods
    └── helper.jl           # optimize_model! and Properties(::MixtureProperties) constructor
```

## Usage

```julia
using Revise
include("load_model.jl")

m = Model()
m[:prop] = Properties(mixture_property("Mixture_18"))

config = InternalCascade_1so1si(
    duty        = 1.0,
    source_temp = [320.0, 300.0],   # [T_high, T_low] of source, K
    sink_temp   = [400.0, 420.0],   # [T_low, T_high] of sink, K
    min_TP      = 0.5               # minimum throttle pressure, bar
)

obj = build!(m, config)
optimize_model!(m, obj)

to_print(m[:streams][1])
```

## Supported Configurations

All configurations are subtypes of `HeatPumpConfig` and follow the same `build!(m, config)` interface.

| Type | Topology | Sources | Sinks | Streams |
|------|----------|---------|-------|---------|
| `InternalCascade_1so1si` | Internal | 1 | 1 | 9 |
| `InternalCascade_1so2si` | Internal | 1 | 2 | 10 |
| `InternalCascade_2so2si` | Internal | 2 | 2 | 9 |
| `InternalCascade_3so3si` | Internal | 3 | 3 | 14 |
| `ExternalCascade_1so2si` | External | 1 | 2 | 6 |
| `ExternalCascade_1so1si` | External (N-stage) | 1 | 1 | 4N |

The naming convention `Xso` / `Ysi` denotes the number of sources and sinks. **Internal** cascade configurations create two-phase streams and phase separators within the working fluid loop; **External** cascade configurations have no two-phase streams. `ExternalCascade_1so1si` accepts an integer `N` controlling the number of cascade stages, with each stage contributing 4 streams. After solving, all stream variables are accessible via `m[:streams]`.

### N-stage external cascade example

```julia
using Revise
include("load_model.jl")

m = Model()
m[:prop] = Properties(mixture_property("Mixture_18"))

config = ExternalCascade_1so1si(
    duty        = 1.0,
    source_temp = [320.0, 300.0],   # [T_high, T_low] of source, K
    sink_temp   = [400.0, 420.0],   # [T_low, T_high] of sink, K
    min_TP      = 0.5,              # minimum throttle pressure, bar
    N           = 2                 # number of cascade stages
)

obj = build!(m, config)
optimize_model!(m, obj)

# 4N streams total; loop i occupies streams 4(i-1)+1 through 4i
for i in 1:n_streams(config)
    println("=== Stream $i ===")
    to_print(m[:streams][i])
    println()
end
```

## Mixture Library

`mixture_examples.jl` defines 44 hydrocarbon mixtures (Mixture_1 through Mixture_44) spanning binary through six-component systems (C3–C8 range). Each entry stores:

- `α` — relative volatilities
- `λ`, `λ_P` — averaged latent heat and latent heat as a function of pressure
- `antoine` — Antoine equation coefficients `[A, B, C]`
- `pRef` — reference pressure (bar)
- `n` — number of components

Retrieve a mixture with `mixture_property("Mixture_18")`, then construct model properties with `Properties(mix)`.

## T-Q Plots

After solving, call `tq_plot(m, config)` to generate a temperature–heat duty diagram for every heat exchanger in the configuration. Each heat exchanger gets its own subplot; hot and cold side curves are drawn in red and blue respectively.

```julia
fig = tq_plot(m, config)
save("tq.png", fig)
```

`tq_data(m, config)` returns the underlying `Vector{HXData}` if you want to build custom plots. Each `HXData` entry has fields `label`, `T_hot`, `T_cold` (both as `[T_in, T_out]`), and `Q`.

Requires [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl).

## Solver Settings

`optimize_model!` defaults to a 600 s time limit and 1% MIP gap. These can be overridden:

```julia
optimize_model!(m, obj; time_limit=1200, mip_gap=0.005)
```

## References

```bibtex
@article{singh2026exergy,
  title={An Exergy-Based Approach for Finding the Composition of Working Fluids in Heat Pump Configurations},
  author={Singh, Priyanshu and Agrawal, Rakesh},
  journal={},
  pages={},
  year={in prep},
  publisher={}
}

@article{singh2026internal,
  title={Internal cascade heat pumps for high temperature lift applications},
  author={Singh, Priyanshu and Nogaja, Akash Sanjay and Agrawal, Rakesh},
  journal={Applied Thermal Engineering},
  pages={130567},
  year={2026},
  publisher={Elsevier}
}

@article{singh2026onecompressor,
  title={A One-Compressor Cascade for Heat Pumping Between Multiple Heat Sources and Heat Sinks},
  author={Singh, Priyanshu and Nogaja, Akash Sanjay and Agrawal, Rakesh},
  journal={Submitted for publication},
  year={2026b}
}

@article{nogaja2022identifying,
  title={Identifying heat-integrated energy-efficient multicomponent distillation configurations},
  author={Nogaja, Akash Sanjay and Mathew, Tony Joseph and Tawarmalani, Mohit and Agrawal, Rakesh},
  journal={Industrial \& Engineering Chemistry Research},
  volume={61},
  number={37},
  pages={13984--13995},
  year={2022},
  publisher={ACS Publications}
}

@article{mathew2023relaxing,
  title={Relaxing the constant molar overflow assumption in distillation optimization},
  author={Mathew, Tony Joseph and Tawarmalani, Mohit and Agrawal, Rakesh},
  journal={AIChE Journal},
  volume={69},
  number={9},
  pages={e18125},
  year={2023},
  publisher={Wiley Online Library}
}

@article{mathew2021simple,
  title={A simple criterion for feasibility of heat integration between distillation streams based on relative volatilities},
  author={Mathew, Tony Joseph and Tumbalam Gooty, Radhakrishna and Tawarmalani, Mohit and Agrawal, Rakesh},
  journal={Industrial \& Engineering Chemistry Research},
  volume={60},
  number={28},
  pages={10286--10302},
  year={2021},
  publisher={ACS Publications}
}




```
