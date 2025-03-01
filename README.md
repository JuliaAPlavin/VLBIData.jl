# VLBIData.jl

Common data types, functions, and calculations for visibilities from astronomical interferometers. While it is primarily focused on radio VLBI, its functionality may also be useful in other contexts.

The current data model and specifics of the interface are flexible and will likely evolve in future updates (with the version number updated accordingly). Contributions to improve or even completely revamp the definitions are welcome!

> [!NOTE]
> Earlier versionds of `VLBIData` (v0.1, v0.2, v0.3) primarily handled file I/O. Since v0.4, `VLBIData` has been dedicated to data model definitions and visibility calculations. For file I/O functionality, please use [VLBIFiles.jl](https://github.com/JuliaAPlavin/VLBIFiles.jl) that relies upon the definitions provided here.
