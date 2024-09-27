# tnp_tools

A repository for performing tag-and-probe analyses to derive simulation correction ("scale") factors. This repository is particularly designed for fits where a large amount of manual intervention is needed for fits.

Note this README is a work in proress.

## Introduction 

In high-energy physics, a common task is to derive corrections for simulated samples to correct differences in selections caused by imperfections in the detector simulation. This can be done by comparing the efficiency of the selection in simulation and data. However, one may not be able to isolate a 100% pure sample in data due to contamination from background. The tag-and-probe method is commonly used to derive a correction, even in the presence of background contamination in data by performing a signal+background fit to a distribution in which the the signal shape is well-known (such as the Z boson mass) to extract the amounts of signal and background events in data, from which efficiencies can be measured.

## Setting the environment

Before running this code, one needs ROOT, python3, and CorrectionLib (TODO: determine minimal versions). For those working on the UCSB servers, there is a file `set_env.sh` that can be used to set up the environment.

## Running the code - basic use

The main analysis is performed with the TnpAnalyzer class. One creates an instance of this class, then calls the `set_input_files`, `set_fitting_variable`, `set_measurement_variable`, `set_preselection`, `add_nd_binning`, and `add_model` methods in order to configure the analysis, then calls `run_interactive`. See the example in [scripts/example.py](scripts/example.py) for a full working example of a T&P analysis script.

When running the analysis script, one will be presented with interactive prompts. Use the command `(h)elp` to access the list of available commands. The typical workflow is to run `(p)roduce` from the main interactive session to generate the histograms to be used in the tag-and-probe analysis, then run `(f)it` for each histogram to extract the amounts of signal and bakcground in each bin. Finally, one runs the `(o)utput` command to generate the summary plots and JSON file with the measured efficiencies. 

Note that running the `(f)it` command from the main interactive session brings up an interactive fitting session. The interactive fitting session lets the user visualize the fit shape as the fit parameters are manually adjusted before running the standard Minuit fit. Again, use `(h)elp` to see the available commands in the interactive session.

## Running the code - typical use

While the basic use instructions can be used to perform generic tag-and-probe analyses, the `rms_sf_analyzer` is preconfigured for typical use cases.

Section in progress.

## The structure of the code

Section in progress.

