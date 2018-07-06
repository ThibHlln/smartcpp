# SMARTcpp - a C++ extension of the rainfall-runoff SMART for Python

SMARTcpp is an open-source C++ extension for the hydrological catchment model SMART in Python. It is licensed under GNU GPL-3.0 (see [licence file](LICENCE.md) provided). SMART is a top-down rainfall-runoff model composed of a soil moisture accounting component and a routing component. It requires rainfall and potential evapotranspiration time series, it features a set of ten parameters, and it yields a discharge time series. This C++ extension is giving access to the calculation of the states, processes, and outputs of the model for one simulation time-step. It is intended to be used in combination with a warpping script in Python, where the loop through the simulation time series is defined.

Mockler, E., O’Loughlin, F., and Bruen, M.: Understanding hydrological flow paths in conceptual catchment models using uncertainty and sensitivity analysis, *Computers & Geosciences*, 90, 66–77,[doi:10.1016/j.cageo.2015.08.015](https://dx.doi.org/10.1016/j.cageo.2015.08.015), 2016

## Model Specifications

### Model Inputs

* aerial rainfall time series [mm/time step]
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
* actual evapotranpiration time series [mm/time step]

## Version History

* 0.1.0 [06 Jul 2018]: First version of SMARTcpp

## Acknowledgment

This tool was developed with the financial support of Ireland's Environmental Protection Agency (Grant Number 2014-W-LS-5).
