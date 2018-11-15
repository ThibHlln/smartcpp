[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI Version](https://badge.fury.io/py/smartcpp.svg)](https://pypi.python.org/pypi/smartcpp)


# SMARTcpp - a C++ accelerator extension of the rainfall-runoff SMART for Python

SMARTcpp is an open-source C++ accelerator extension for the hydrological catchment model SMART in Python. It is licensed under GNU GPL-3.0 (see [licence file](https://github.com/ThibHlln/smartcpp/blob/master/LICENCE.md) provided). SMART (Soil Moisture Accounting and Routing for Transport) is a top-down rainfall-runoff model composed of a soil moisture accounting component and a linear routing component. It requires rainfall and potential evapotranspiration time series, it features a set of ten parameters, and it yields a discharge time series. This C++ extension is giving access to the calculation of the states, processes, and outputs of the model for one simulation time step or for all time steps. SMARTcpp is intended to be used in combination with a wrapping script in Python.

## How to Install

SMARTcpp is available on PyPI for Python 2.7 and Python 3.6 (both for macOS and Windows 64 bits), so you can simply use pip:

    python -m pip install smartcpp

You can also use a link to the GitHub repository directly:

	python -m pip install git+https://github.com/ThibHlln/smartcpp.git

Alternatively, you can download the source code (*i.e.* the GitHub repository) and, from the downloaded directory itself, run the command:

    python setup.py install

## Model Specifications

### Model Inputs

* rainfall time series [mm/time step]
* potential evapotranspiration time series [mm/time step]

### Model Parameters

* T: rainfall aerial correction coefficient [-]
* C: evaporation decay parameter [-]
* H: quick runoff coefficient [-]
* D: drain flow parameter - fraction of saturation excess diverted to drain flow [-]
* S: soil outflow coefficient [-]
* Z: effective soil depth [mm]
* SK: surface routing parameter [hours]
* FK: inter flow routing parameter [hours]
* GK: groundwater routing parameter [hours]
* RK: river channel routing parameter [hours]

### Model Outputs

* discharge time series at catchment outlet [m<sup>3</sup>/s]
* groundwater contribution to catchment runoff [-]

### References

Mockler, E., O’Loughlin, F., and Bruen, M.: Understanding hydrological flow paths in conceptual catchment models using uncertainty and sensitivity analysis, *Computers & Geosciences*, 90, 66–77,[doi:10.1016/j.cageo.2015.08.015](https://dx.doi.org/10.1016/j.cageo.2015.08.015), 2016

## Version History

* 0.2.0 [15 Nov 2018]: Version capable of running all steps
    * Adds method allsteps that brings the loop into C++
* 0.1.2 [18 Jul 2018]: Version with proper PyPI display
	* Fixes display issue of README.md on PyPI
* 0.1.1 [18 Jul 2018]: Version with Python 3.x compatibility
	* Adds compatibility with Python 3.x extensions
* 0.1.0 [09 Jul 2018]: First version of SMARTcpp

## Acknowledgment

This tool was developed with the financial support of Ireland's Environmental Protection Agency (Grant Number 2014-W-LS-5).
