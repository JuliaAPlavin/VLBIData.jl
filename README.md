# VLBIData.jl 📡

Comprehensive data types, functions, and calculations for visibilities from astronomical interferometers. While primarily focused on radio VLBI, its functionality can also be useful for other interferometric contexts.

## ✨ Features

- **🏗️ Core data structures**: `Antenna`, `Baseline`, `UV` coordinates, and visibility specifications (`VisSpec`, `VisAmpSpec`)
- **📊 Visibility tables**: Simple, flexible data tables requiring only `spec` and `value` fields, with optional fields including `datetime`, `freq_spec`, `stokes`, `count`
- **🔄 Closure quantities**: `ClosurePhaseSpec` and `ClosureAmpSpec` with automatic computation from visibility tables
- **⚖️ Visibility error rescaling**: Multiple methods to correct thermal noise estimates when original weights are improperly scaled
- **📈 Data averaging**: Flexible time-based averaging with proper uncertainty propagation
- **🌊 Polarization handling**: Full support for Stokes parameters and conversion to coherency matrices

## 🔗 Related Packages

- **[VLBIFiles.jl](https://github.com/JuliaAPlavin/VLBIFiles.jl)**: File I/O for VLBI data formats (images, uvfits, models, etc.)
- **[InterferometricModels.jl](https://github.com/JuliaAPlavin/InterferometricModels.jl)**: Source models for interferometric observations  
- **[VLBIPlots.jl](https://github.com/JuliaAPlavin/VLBIPlots.jl)**: Plotting and visualization for VLBI data and models

## 🚧 Development

The data model and interface are designed to be flexible and will evolve in future updates (with appropriate version bumps). Contributions to improve or completely revamp the definitions are very welcome!

> [!NOTE]
> Earlier versions of `VLBIData` (v0.1, v0.2, v0.3) primarily handled file I/O. Since v0.4, `VLBIData` has been dedicated to data model definitions and visibility calculations. For file I/O functionality, please use [VLBIFiles.jl](https://github.com/JuliaAPlavin/VLBIFiles.jl) that relies upon the definitions provided here.
