# tnp_tools

A repository for performing tag-and-probe analyses to derive simulation correction ("scale") factors. This repository is particularly designed for fits where a large amount of manual intervention is needed, which is common when there is a large background component demanding fitting functions with a large number of parameters.

Note this README is a work in proress.

## Introduction 

In high-energy physics, a common task is to derive corrections for simulated samples to correct differences in selections caused by imperfections in the detector simulation. This can be done by comparing the efficiency of the selection in simulation and data. However, one may not be able to isolate a 100% pure sample in data due to contamination from background. The tag-and-probe method is commonly used to derive a correction, even in the presence of background contamination in data by fitting a variable with a well-known shape such as the Z boson resonance and performing a simultaneous fit to signal and background events in data and simulation.

## Setting the environment

Before running this code, one needs ROOT, python3, and CorrectionLib (TODO: determine minimal versions). For those working on the UCSB servers, there is a file `set_env.sh` that can be used to set up the environment.

## Running the code



