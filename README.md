# tnp_tools

A repository for performing tag-and-probe analyses to derive simulation correction ("scale") factors. This repository is particularly designed for fits where a large amount of manual intervention is needed for fits.

Note this README is a work in proress.

## Introduction 

In high-energy physics, a common task is to derive corrections for simulated samples to correct differences in selections caused by imperfections in the detector simulation. This can be done by comparing the efficiency of the selection in simulation and data. However, one may not be able to isolate a 100% pure sample in data due to contamination from background. The tag-and-probe method is commonly used to derive a correction, even in the presence of background contamination in data by performing a signal+background fit to a distribution in which the the signal shape is well-known (such as the Z boson mass) to extract the amounts of signal and background events in data, from which efficiencies can be measured.

## Setting the environment

Before running this code, one needs ROOT, pdflatex, python3, and CorrectionLib (TODO: determine minimal versions). For those working on the UCSB servers, there is a file `set_env.sh` that can be used to set up the environment.

## Running the code - typical use

### Driver script

For typical use cases, the `rms_sf_analyzer` class is used to set up the analysis. An example can be found in [scripts/elid_wpl.py](scripts/elid_wpl.py). After initializing the RmsSFAnalyzer with a unique name, one sets the year, input files, fitting variable, measurement variable, preselection cuts, and binning. Then one calls the `run_interactive` method.

When this script is run, one will be sent to an interactive analysis command line tool. You can use `(h)elp` to see the available options and `(q)uit` to exit. The typical top-level commands for an analysis are

```
p       #produces histograms to fit
f nom   #starts interactive fitter on data using default model
f mc    #starts interactive fitter on MC using alternative signal-only model
f alts  #starts interactive fitter on data using alternative signal model
f altb  #starts interactive fitter on data using alternative background model
f altsb #starts interactive fitter on data using alt signal+background model
o       #produces output once all fits have been performed
```

`(p)roduce` must be run before any `(f)it` commands and all fits must be done
before running `(o)utput`. Furthermore, the `mc` fit must be done before the
`alts` and `altsb` fits because the MC fit is used to constraint fit parameters in the alternative signal model. The actual MC efficiencies for the nominal and alternative MC samples are generated without a fit since there is no background to consider.

The `(f)it` command defaultly starts with the passing leg of the first bin, but if you need to do a fit for a particular bin and leg, you can use `(f)it <sample> <bin> <pf>` where `<bin>` is the bin number and `<pf>` is `(p)ass` or `(f)ail`.

The interactive fitting tool allows one to visualize the fitting procedure (using ex. X11, VNC, or similar). The interactive fitting tool similarly provides a set of commands than can be seen with `(h)elp`. When started, the tool will automatically perform a fit. If this default fit is satisfactory, one can confirm it using either `(q)uit` or `(n)ext`. `(q)uit` exits the interactive fitter and returns to the top level menu while `(n)ext` goes to the next fit, progressing the bin number or wrapping around and going to the first failing bin after the last passing bin. To exit without saving the fit, use `(q)uit(!)`.

If the initial fit is not satisfactory, the tool allows fine user control to help make the fit converge. `(l)ist` lists the parameters of the fitting model and their current values. One can manually set the parameters using `(s)et <parameter> <value>`. After setting the parameters to reasonable values by hand, one can call `(f)it` to redo the fit. 

## Running the code - general use

Information in this section can be used to write more general tag-and-probe analyses.

The main analysis is performed with the TnpAnalyzer class. One creates an instance of this class, then calls the `set_input_files`, `set_fitting_variable`, `set_measurement_variable`, `set_preselection`, `add_nd_binning`, and `add_model` methods in order to configure the analysis, then calls `run_interactive`. See the example in [scripts/example.py](scripts/example.py) for a full working example of a T&P analysis script.

When running the analysis script, one will be presented with interactive prompts. Use the command `(h)elp` to access the list of available commands. The typical workflow is to run `(p)roduce` from the main interactive session to generate the histograms to be used in the tag-and-probe analysis, then run `(f)it` for each histogram to extract the amounts of signal and bakcground in each bin. Finally, one runs the `(o)utput` command to generate the summary plots and JSON file with the measured efficiencies. 

Note that running the `(f)it` command from the main interactive session brings up an interactive fitting session. The interactive fitting session lets the user visualize the fit shape as the fit parameters are manually adjusted before running the standard Minuit fit. Again, use `(h)elp` to see the available commands in the interactive session.

## The structure of the code

Section in progress.

