# rotornet
Matlab scripts to model RotorNet in the 2017 Sigcomm paper

Artifacts from "RotorNet: A Scalable, Low-complexity, Optical
Datacenter Network," Sigcomm 2017

09/20/2017

Requirements: Matlab 2016b or newer
	(may work with older Matlab versions or with GNU Octave, but
    not tested)

Description:

- Matlab scipts and datasets are arranged into directories, one for
  each graph in the paper.
- Each directory contains a plotting script "PlotFigX.m", where X
  ranges from 7-12.
- At the top of each script "PlotFigX.m", setting the variable
  "rerunsims" = 0 will produce the graph from the existing .mat files
  without rerunning the simulation function.
- Setting "rerunsims" = 1 will rerun the simulation function, generate
  one or more new .mat files, and produce the graph based on those new
  .mat files.
