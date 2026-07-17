"""
    SDM

Loader for ASDM (Science Data Model) datasets, e.g. EVLA/ALMA raw observations.

Use:

```julia
using AURIS
ds   = SDM.open_sdm("/path/to/sdm")
main = ds.mains[1]                              # one Main row == one BDF file
bdf  = SDM.open_bdf(ds, main)
amps = SDM.load_amplitudes(bdf; stokes=:I)      # (n_time, n_chan, n_baseline)
close(bdf)
sp   = SDM.load_syspower(ds)                     # binary SysPower table
```

Split across three files:

  * `tables.jl`    – SDM XML tables → `Components` structs, and `open_sdm`.
  * `bdf.jl`       – Binary Data Format visibility reader (`open_bdf`, …).
  * `syspower.jl`  – the MIME-wrapped binary `SysPower` table (`load_syspower`).
"""
module SDM

using Mmap
using ...Components
using ...Components: Main   # disambiguate the SDM `Main` table from `Core.Main`

include("tables.jl")
include("bdf.jl")
include("syspower.jl")

export SDMDataset, open_sdm, bdf_path, asdm_time_to_unix
export BDFFile, open_bdf, close_bdf, load_amplitudes, integration_times
export load_syspower

end 
