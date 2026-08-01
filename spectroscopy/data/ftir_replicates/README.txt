ATR-FTIR replicates, recorded on a diamond ATR by James Sturgis, October 2025.

  Glucose.1-3.dpt      D-glucose, three replicate acquisitions
  CelluloseX.0-2.dpt   microcrystalline cellulose, three replicates
  H2O.0-2.dpt          water, three replicates -- the reference to subtract

Enough to demonstrate the whole ATR-FTIR workflow: load a folder, group by
sample, average the replicates, subtract the water contribution, baseline
correct, normalise and pick peaks.

Tab separated, no header, one x/y pair per line, CRLF line endings -- exactly
as Bruker OPUS writes them.
