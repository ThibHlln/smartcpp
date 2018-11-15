// -*- coding: utf-8 -*-

// This file is part of SMARTcpp - An open-source C++ extension of the SMART model in Python
// Copyright (C) 2018  Thibault Hallouin (1), Eva Mockler (1,2), Michael Bruen (1)
//
// (1) Dooge Centre for Water Resources Research, University College Dublin, Ireland
// (2) Environmental Protection Agency, Ireland
//
// SMARTcpp is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SMARTcpp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SMARTcpp. If not, see <http://www.gnu.org/licenses/>.

#include <Python.h>
#include <numpy/arrayobject.h>

#include <iostream>
#include <cmath>

class Catchment{
    public:
        // class members
        double out_aeva;
        double out_q_h2o_ove;
        double out_q_h2o_dra;
        double out_q_h2o_int;
        double out_q_h2o_sgw;
        double out_q_h2o_dgw;
        double s_v_h2o_ove;
        double s_v_h2o_dra;
        double s_v_h2o_int;
        double s_v_h2o_sgw;
        double s_v_h2o_dgw;
        double s_v_h2o_ly1;
        double s_v_h2o_ly2;
        double s_v_h2o_ly3;
        double s_v_h2o_ly4;
        double s_v_h2o_ly5;
        double s_v_h2o_ly6;
        double pr_eff_rain_to_ove;
        double pr_eff_rain_to_dra;
        double pr_eff_rain_to_int;
        double pr_eff_rain_to_sgw;
        double pr_eff_rain_to_dgw;
        // class constructor
        Catchment(){
            out_aeva = 0.0;
            out_q_h2o_ove = 0.0;
            out_q_h2o_dra = 0.0;
            out_q_h2o_int = 0.0;
            out_q_h2o_sgw = 0.0;
            out_q_h2o_dgw = 0.0;
            s_v_h2o_ove = 0.0;
            s_v_h2o_dra = 0.0;
            s_v_h2o_int = 0.0;
            s_v_h2o_sgw = 0.0;
            s_v_h2o_dgw = 0.0;
            s_v_h2o_ly1 = 0.0;
            s_v_h2o_ly2 = 0.0;
            s_v_h2o_ly3 = 0.0;
            s_v_h2o_ly4 = 0.0;
            s_v_h2o_ly5 = 0.0;
            s_v_h2o_ly6 = 0.0;
            pr_eff_rain_to_ove = 0.0;
            pr_eff_rain_to_dra = 0.0;
            pr_eff_rain_to_int = 0.0;
            pr_eff_rain_to_sgw = 0.0;
            pr_eff_rain_to_dgw = 0.0;
        }
};

class River{
    public:
        // class members
        double out_q_riv;
        double s_v_riv;
        // class constructor
        River(){
            out_q_riv = 0.0;
            s_v_riv = 0.0;
        }
};

static Catchment onestep_catchment(
        double area_m2, double time_delta_sec,
        double c_in_rain, double c_in_peva,
        double c_p_t, double c_p_c, double c_p_h, double c_p_d, double c_p_s, double c_p_z, double c_p_sk, double c_p_fk, double c_p_gk,
        double c_s_v_h2o_ove, double c_s_v_h2o_dra, double c_s_v_h2o_int, double c_s_v_h2o_sgw, double c_s_v_h2o_dgw,
        double c_s_v_h2o_ly1, double c_s_v_h2o_ly2, double c_s_v_h2o_ly3, double c_s_v_h2o_ly4, double c_s_v_h2o_ly5, double c_s_v_h2o_ly6){

    /*
    Catchment Constants
    _ area_m2                   catchment area [m2]
    _ time_gap_min              time gap between two simulation time steps [minutes]

    Catchment model * c_ *
    _ Hydrology
    ___ Inputs * in_ *
    _____ c_in_rain             precipitation as rain [mm/time step]
    _____ c_in_peva             potential evapotranspiration [mm/time step]
    ___ Parameters * p_ *
    _____ c_p_t                 T: rainfall aerial correction coefficient
    _____ c_p_c                 C: evaporation decay parameter
    _____ c_p_h                 H: quick runoff coefficient
    _____ c_p_d                 D: drain flow parameter - fraction of saturation excess diverted to drain flow
    _____ c_p_s                 S: soil outflow coefficient
    _____ c_p_z                 Z: effective soil depth [mm]
    _____ c_p_sk                SK: surface routing parameter [hours]
    _____ c_p_fk                FK: inter flow routing parameter [hours]
    _____ c_p_gk                GK: groundwater routing parameter [hours]
    ___ States * s_ *
    _____ c_s_v_h2o_ove         volume of water in overland store [m3]
    _____ c_s_v_h2o_dra         volume of water in drain store [m3]
    _____ c_s_v_h2o_int         volume of water in inter store [m3]
    _____ c_s_v_h2o_sgw         volume of water in shallow groundwater store [m3]
    _____ c_s_v_h2o_dgw         volume of water in deep groundwater store [m3]
    _____ c_s_v_h2o_ly1         volume of water in first soil layer store [m3]
    _____ c_s_v_h2o_ly2         volume of water in second soil layer store [m3]
    _____ c_s_v_h2o_ly3         volume of water in third soil layer store [m3]
    _____ c_s_v_h2o_ly4         volume of water in fourth soil layer store [m3]
    _____ c_s_v_h2o_ly5         volume of water in fifth soil layer store [m3]
    _____ c_s_v_h2o_ly6         volume of water in sixth soil layer store [m3]
    ___ Processes * pr_ *
    _____ c_pr_eff_rain_to_ove  effective rainfall converted into overland flow runoff [mm]
    _____ c_pr_eff_rain_to_dra  effective rainfall converted into to drain flow runoff [mm]
    _____ c_pr_eff_rain_to_int  effective rainfall converted into to interflow runoff [mm]
    _____ c_pr_eff_rain_to_sgw  effective rainfall converted into to shallow groundwater flow runoff [mm]
    _____ c_pr_eff_rain_to_dgw  effective rainfall converted into to deep groundwater flow runoff [mm]
    ___ Outputs * out_ *
    _____ c_out_aeva            actual evapotranspiration [m3/s]
    _____ c_out_q_h2o_ove       overland flow [m3/s]
    _____ c_out_q_h2o_dra       drain flow [m3/s]
    _____ c_out_q_h2o_int       inter flow [m3/s]
    _____ c_out_q_h2o_sgw       shallow groundwater flow [m3/s]
    _____ c_out_q_h2o_dgw       deep groundwater flow [m3/s]
    */

    Catchment c;

    // 1. Hydrology
    // 1.0. Define internal constants
    double nb_soil_layers = 6.0;  // number of layers in soil column [-]

    // 1.1. Convert non-SI units
    c_p_sk *= 3600.0;  // convert hours in seconds
    c_p_fk *= 3600.0;  // convert hours in seconds
    c_p_gk *= 3600.0;  // convert hours in seconds

    // 1.2. Hydrological calculations

    // /!\ all calculations in mm equivalent until further notice

    // calculate capacity Z and level LVL of each layer (assumed equal) from effective soil depth
    double list_z_lyr[7];
    list_z_lyr[0] = 0.0;  // artificial null value added to keep script clear later
    list_z_lyr[1] = c_p_z / nb_soil_layers;  // Soil Layer 1
    list_z_lyr[2] = c_p_z / nb_soil_layers;  // Soil Layer 2
    list_z_lyr[3] = c_p_z / nb_soil_layers;  // Soil Layer 3
    list_z_lyr[4] = c_p_z / nb_soil_layers;  // Soil Layer 4
    list_z_lyr[5] = c_p_z / nb_soil_layers;  // Soil Layer 5
    list_z_lyr[6] = c_p_z / nb_soil_layers;  // Soil Layer 6

    double list_lvl_lyr[7];
    list_lvl_lyr[0] = 0.0;  // artificial null value added to keep script clear later
    list_lvl_lyr[1] = c_s_v_h2o_ly1 / area_m2 * 1000;  // Soil Layer 1
    list_lvl_lyr[2] = c_s_v_h2o_ly2 / area_m2 * 1000;  // Soil Layer 2
    list_lvl_lyr[3] = c_s_v_h2o_ly3 / area_m2 * 1000;  // Soil Layer 3
    list_lvl_lyr[4] = c_s_v_h2o_ly4 / area_m2 * 1000;  // Soil Layer 4
    list_lvl_lyr[5] = c_s_v_h2o_ly5 / area_m2 * 1000;  // Soil Layer 5
    list_lvl_lyr[6] = c_s_v_h2o_ly6 / area_m2 * 1000;  // Soil Layer 6

    // calculate cumulative level of water in all soil layers at beginning of time step (i.e. soil moisture)
    double lvl_total_start = list_lvl_lyr[1] + list_lvl_lyr[2] + list_lvl_lyr[3] + list_lvl_lyr[4] + list_lvl_lyr[5] + list_lvl_lyr[6];

    // apply parameter T to rainfall data (aerial rainfall correction)
    double rain = c_in_rain * c_p_t;
    // calculate excess rainfall
    double excess_rain = rain - c_in_peva;
    // initialise actual evapotranspiration variable
    double aeva = 0.0;
    // initialise effective rainfall to runoff pathways
    double c_pr_eff_rain_to_ove, c_pr_eff_rain_to_dra, c_pr_eff_rain_to_int, c_pr_eff_rain_to_sgw, c_pr_eff_rain_to_dgw;

    if (excess_rain >= 0.0) {  // excess rainfall available for runoff and infiltration
        // actual evapotranspiration = potential evapotranspiration
        aeva += c_in_peva;

        // calculate surface runoff using quick runoff parameter H and relative soil moisture content
        double h_prime = c_p_h * (lvl_total_start / c_p_z);
        c_pr_eff_rain_to_ove = h_prime * excess_rain;  // excess rainfall contribution to quick surface runoff store
        excess_rain -= c_pr_eff_rain_to_ove;  // remainder that infiltrates

        // calculate percolation through soil layers (from top layer [1] to bottom layer [6])
        for (int i = 1; i <= 6; i++) {
            double space_in_lyr = list_z_lyr[i] - list_lvl_lyr[i];
            if (excess_rain <= space_in_lyr) {
                list_lvl_lyr[i] += excess_rain;
                excess_rain = 0.0;
            } else {
                list_lvl_lyr[i] = list_z_lyr[i];
                excess_rain -= space_in_lyr;
            }
        }

        // calculate saturation excess from remaining excess rainfall after filling layers (if not 0)
        c_pr_eff_rain_to_dra = c_p_d * excess_rain;  // sat. excess contribution (if not 0) to quick interflow runoff store
        c_pr_eff_rain_to_int = (1.0 - c_p_d) * excess_rain;  // sat. excess contribution (if not 0) to slow interflow runoff store

        // calculate leak from soil layers (i.e. piston flow becoming active during rainfall events)
        double s_prime = c_p_s * (lvl_total_start / c_p_z);
        // leak to interflow
        for (int i = 1; i <= 6; i++) {  // soil moisture outflow reducing exponentially downwards
            double leak_interflow = list_lvl_lyr[i] * pow(s_prime, i);
            if (leak_interflow < list_lvl_lyr[i]) {
                c_pr_eff_rain_to_int += leak_interflow;  // soil moisture outflow contribution to slow interflow runoff store
                list_lvl_lyr[i] -= leak_interflow;
            }
        }
        // leak to shallow groundwater flow
        c_pr_eff_rain_to_sgw = 0.0;
        for (int i = 1; i <= 6; i++) {  // soil moisture outflow reducing linearly downwards
            double leak_shallow_flow = list_lvl_lyr[i] * (s_prime / i);
            if (leak_shallow_flow < list_lvl_lyr[i]) {
                c_pr_eff_rain_to_sgw += leak_shallow_flow;  // soil moisture outflow contribution to slow shallow GW runoff store
                list_lvl_lyr[i] -= leak_shallow_flow;
            }
        }
        // leak to deep groundwater flow
        c_pr_eff_rain_to_dgw = 0.0;
        for (int i = 1; i <= 6; i++) {  // soil moisture outflow reducing exponentially upwards
            double leak_deep_flow = list_lvl_lyr[i] * pow(s_prime, 7 - i);
            if (leak_deep_flow < list_lvl_lyr[i]) {
                c_pr_eff_rain_to_dgw += leak_deep_flow;  // soil moisture outflow contribution to slow deep GW runoff store
                list_lvl_lyr[i] -= leak_deep_flow;
            }
        }
    } else {  // no excess rainfall (i.e. potential evapotranspiration not satisfied by available rainfall)
        c_pr_eff_rain_to_ove = 0.0;  // no effective rainfall contribution to quick overland flow runoff store
        c_pr_eff_rain_to_dra = 0.0;  // no effective rainfall contribution to quick drain flow runoff store
        c_pr_eff_rain_to_int = 0.0;  // no effective rainfall contribution to quick + leak interflow runoff store
        c_pr_eff_rain_to_sgw = 0.0;  // no effective rainfall contribution to shallow groundwater flow runoff store
        c_pr_eff_rain_to_dgw = 0.0;  // no effective rainfall contribution to deep groundwater flow runoff store

        double deficit_rain = excess_rain * (-1.0);  // excess is negative => excess is actually a deficit
        aeva += rain;
        for (int i = 1; i <= 6; i++) {  // attempt to satisfy PE from soil layers (from top layer [1] to bottom layer [6]
            if (list_lvl_lyr[i] >= deficit_rain) {  // i.e. all moisture required available in this soil layer
                list_lvl_lyr[i] -= deficit_rain;  // soil layer is reduced by the moisture required
                aeva += deficit_rain;  // this moisture contributes to the actual evapotranspiration
                deficit_rain = 0.0; // the full moisture still required has been met
            } else { // i.e. not all moisture required available in this soil layer
                aeva += list_lvl_lyr[i];  // takes what is available in this layer for evapotranspiration
                // effectively reduce the evapotranspiration demand for the next layer using parameter C
                // i.e. the more you move down through the soil layers, the less AET can meet PET (exponentially)
                deficit_rain = c_p_c * (deficit_rain - list_lvl_lyr[i]);
                list_lvl_lyr[i] = 0.0; // soil layer is now empty
            }
        }
    }

    // /!\ all calculations in S.I. units now (i.e. mm converted into cubic metres)

    // calculate actual evapotranspiration as a flux
    c.out_aeva = aeva / 1e3 * area_m2 / time_delta_sec;  // [m3/s]

    // route overland flow (quick surface runoff)
    c.out_q_h2o_ove = c_s_v_h2o_ove / c_p_sk;  // [m3/s]
    c.s_v_h2o_ove = c_s_v_h2o_ove + (c_pr_eff_rain_to_ove / 1e3 * area_m2) - (c.out_q_h2o_ove * time_delta_sec);  // [m3] - [m3]
    if (c.s_v_h2o_ove < 0.0)
        c.s_v_h2o_ove = 0.0;
    // route drain flow (quick interflow runoff)
    c.out_q_h2o_dra = c_s_v_h2o_dra / c_p_sk;  // [m3/s]
    c.s_v_h2o_dra = c_s_v_h2o_dra + (c_pr_eff_rain_to_dra / 1e3 * area_m2) - (c.out_q_h2o_dra * time_delta_sec);  // [m3] - [m3]
    if (c.s_v_h2o_dra < 0.0)
        c.s_v_h2o_dra = 0.0;
    // route interflow (slow interflow runoff)
    c.out_q_h2o_int = c_s_v_h2o_int / c_p_fk;  // [m3/s]
    c.s_v_h2o_int = c_s_v_h2o_int + (c_pr_eff_rain_to_int / 1e3 * area_m2) - (c.out_q_h2o_int * time_delta_sec);  // [m3] - [m3]
    if (c.s_v_h2o_int < 0.0)
        c.s_v_h2o_int = 0.0;
    // route shallow groundwater flow (slow shallow GW runoff)
    c.out_q_h2o_sgw = c_s_v_h2o_sgw / c_p_gk;  // [m3/s]
    c.s_v_h2o_sgw = c_s_v_h2o_sgw + (c_pr_eff_rain_to_sgw / 1e3 * area_m2) - (c.out_q_h2o_sgw * time_delta_sec);  // [m3] - [m3]
    if (c.s_v_h2o_sgw < 0.0)
        c.s_v_h2o_sgw = 0.0;
    // route deep groundwater flow (slow deep GW runoff)
    c.out_q_h2o_dgw = c_s_v_h2o_dgw / c_p_gk;  // [m3/s]
    c.s_v_h2o_dgw = c_s_v_h2o_dgw + (c_pr_eff_rain_to_dgw / 1e3 * area_m2) - (c.out_q_h2o_dgw * time_delta_sec);  // [m3] - [m3]
    if (c.s_v_h2o_dgw < 0.0)
        c.s_v_h2o_dgw = 0.0;

    // update soil layer levels
    c.s_v_h2o_ly1 = list_lvl_lyr[1] / 1000 * area_m2;
    c.s_v_h2o_ly2 = list_lvl_lyr[2] / 1000 * area_m2;
    c.s_v_h2o_ly3 = list_lvl_lyr[3] / 1000 * area_m2;
    c.s_v_h2o_ly4 = list_lvl_lyr[4] / 1000 * area_m2;
    c.s_v_h2o_ly5 = list_lvl_lyr[5] / 1000 * area_m2;
    c.s_v_h2o_ly6 = list_lvl_lyr[6] / 1000 * area_m2;

    // store internal process variables
    c.pr_eff_rain_to_ove = c_pr_eff_rain_to_ove;
    c.pr_eff_rain_to_dra = c_pr_eff_rain_to_dra;
    c.pr_eff_rain_to_int = c_pr_eff_rain_to_int;
    c.pr_eff_rain_to_sgw = c_pr_eff_rain_to_sgw;
    c.pr_eff_rain_to_dgw = c_pr_eff_rain_to_dgw;

    return c;
}

static River onestep_river(
        double time_delta_sec,
        double r_in_q_riv,
        double r_p_rk,
        double r_s_v_riv){

    /*
    River model * r_ *
    _ Hydrology
    ___ Inputs * in_ *
    _____ r_in_q_riv      flow at inlet [m3/s]
    ___ Parameters * p_ *
    _____ r_p_rk          linear factor k for water where Storage = k.Flow [hours]
    ___ States * s_ *
    _____ r_s_v_riv       volume of water in store [m3]
    ___ Outputs * out_ *
    _____ r_out_q_riv     flow at outlet [m3/s]
    */

    River r;

    // 1. Hydrology
    // 1.0. Define internal constants
    r_p_rk *= 3600.0;  // convert hours in seconds

    // 1.1. Hydrological calculations

    // calculate outflow, at current time step
    r.out_q_riv = r_s_v_riv / r_p_rk;
    // calculate storage in temporary variable, for next time step
    double r_s_v_h2o_old = r_s_v_riv;
    double r_s_v_h2o_temp = r_s_v_h2o_old + (r_in_q_riv - r.out_q_riv) * time_delta_sec;
    // check if storage has gone negative
    if (r_s_v_h2o_temp < 0.0) {  // temporary cannot be used
        // constrain outflow: allow maximum outflow at 95% of what was in store
        r.out_q_riv = 0.95 * (r_in_q_riv + r_s_v_h2o_old / time_delta_sec);
        // calculate final storage with constrained outflow
        r.s_v_riv = r_s_v_riv + (r_in_q_riv - r.out_q_riv) * time_delta_sec;
    } else {
        r.s_v_riv = r_s_v_h2o_temp;  // temporary storage becomes final storage
    }

    return r;
}


static PyObject *smartcpp_allsteps(PyObject *self, PyObject *args) {
    /* interface variables */
    double area_m2, time_delta_sec;  // constants
    int simu_length; // constants
    PyArrayObject *nd_rain, *nd_peva;  // inputs
    PyArrayObject *nd_parameters;  // model parameters
    PyArrayObject *nd_initial;  // initial conditions
    PyArrayObject *nd_discharge, *nd_history;  // outputs
    double groundwater_component;  // outputs
    int report_type, report_gap;  // reporting information

    /* internal variables */
    npy_intp dim_q[1], dim_h[1];  // arrays for the lengths of the dimensions of the discharge and history arrays
    int marker;  // a marker to know when a reporting step has been reached
    int counter;  // a counter to know which index of the discharge time series to save at
    double riv, gw_ro, all_ro;  // tallies for calculation of discharge and groundwater contribution

    /* Parse the arguments coming from the Python function call */
    if (!PyArg_ParseTuple(args, "ddiO!O!O!O!ii",
                          &area_m2, &time_delta_sec, &simu_length,
                          &PyArray_Type, &nd_rain, &PyArray_Type, &nd_peva,
                          &PyArray_Type, &nd_parameters, &PyArray_Type, &nd_initial,
                          &report_type, &report_gap)) {
        return NULL;
    }

    /* Create C arrays with the same dimensions as numpy arrays */
    int dim = simu_length / report_gap;
    double *arr_discharge = new double[dim];
    double arr_history[19];
    double arr_parameters[10];

    /* Get initial conditions and model parameter values in C arrays */
    for (int i = 0; i < 19; i++)
    {
        arr_history[i] = *(double *)PyArray_GETPTR1(nd_initial, i);
    }
    for (int i = 0; i < 10; i++)
    {
        arr_parameters[i] = *(double *)PyArray_GETPTR1(nd_parameters, i);
    }

    /* Loop through the data set and call the one step functions */
    Catchment c;
    River r;
    marker = 1;  // initialise the marker
    counter = 0;  // initialise the counter
    riv = 0.0;  // initialise tally for river flow
    gw_ro = 0.0;  // initialise tally for groundwater runoff
    all_ro = 0.0;  // initialise tally for total runoff

    for (int i = 1; i <= simu_length; i++)
    {
        /* Run the catchment model for one time step */
        c = onestep_catchment(
                area_m2, time_delta_sec,
                *(double *)PyArray_GETPTR1(nd_rain, i - 1),
                *(double *)PyArray_GETPTR1(nd_peva, i - 1),
                arr_parameters[0], arr_parameters[1], arr_parameters[2], arr_parameters[3], arr_parameters[4],
                arr_parameters[5], arr_parameters[6], arr_parameters[7], arr_parameters[8],
                arr_history[7], arr_history[8], arr_history[9], arr_history[10],
                arr_history[11], arr_history[12], arr_history[13], arr_history[14],
                arr_history[15], arr_history[16], arr_history[17]);

        /* Run the river model for one time step */
        r = onestep_river(
                time_delta_sec,
                c.out_q_h2o_ove + c.out_q_h2o_dra + c.out_q_h2o_int + c.out_q_h2o_sgw + c.out_q_h2o_dgw,
                arr_parameters[9],
                arr_history[18]);

        /* Store the current time step state values for next time step */
        arr_history[7] = c.s_v_h2o_ove;
        arr_history[8] = c.s_v_h2o_dra;
        arr_history[9] = c.s_v_h2o_int;
        arr_history[10] = c.s_v_h2o_sgw;
        arr_history[11] = c.s_v_h2o_dgw;
        arr_history[12] = c.s_v_h2o_ly1;
        arr_history[13] = c.s_v_h2o_ly2;
        arr_history[14] = c.s_v_h2o_ly3;
        arr_history[15] = c.s_v_h2o_ly4;
        arr_history[16] = c.s_v_h2o_ly5;
        arr_history[17] = c.s_v_h2o_ly6;
        arr_history[18] = r.s_v_riv;

        /* Save the necessary output values in C arrays */
        if (report_type == 1)  // summary of the simulation steps included in one reporting step
        {
            if (marker == report_gap)  // reached a reporting step
            {
                riv += r.out_q_riv;
                gw_ro += c.out_q_h2o_sgw + c.out_q_h2o_dgw;
                all_ro += c.out_q_h2o_ove + c.out_q_h2o_dra + c.out_q_h2o_int + c.out_q_h2o_sgw + c.out_q_h2o_dgw;
                arr_discharge[counter] = riv / report_gap;
                riv = 0;
                marker = 0;
                counter++;
            }
            else  // between two reporting steps
            {
                riv += r.out_q_riv;
                gw_ro += (c.out_q_h2o_sgw + c.out_q_h2o_dgw);
                all_ro += (c.out_q_h2o_ove + c.out_q_h2o_dra + c.out_q_h2o_int + c.out_q_h2o_sgw + c.out_q_h2o_dgw);
            }
        }
        else  // raw extraction of the value at the corresponding to the reporting step
        {
            if (marker == report_gap)  // reached a reporting step
            {
                arr_discharge[counter] = r.out_q_riv;
                gw_ro += (c.out_q_h2o_sgw + c.out_q_h2o_dgw);
                all_ro += (c.out_q_h2o_ove + c.out_q_h2o_dra + c.out_q_h2o_int + c.out_q_h2o_sgw + c.out_q_h2o_dgw);
                marker = 0;
                counter++;
            }
            // else // do nothing, no reporting required
        }

        marker++;
    }

    groundwater_component = gw_ro / all_ro;

    /* Convert C arrays for discharge and storage to a numpy arrays */
    dim_q[0] = simu_length / report_gap;
    nd_discharge = (PyArrayObject *) PyArray_SimpleNewFromData(1, dim_q, NPY_DOUBLE, arr_discharge);
    dim_h[0] = 19;
    nd_history = (PyArrayObject *) PyArray_SimpleNewFromData(1, dim_h, NPY_DOUBLE, arr_history);

    /* Return a tuple containing the discharge time series, the groundwater component, and the outputs+states history */
    PyObject *tuple_result = PyTuple_New(3);
    PyTuple_SetItem(tuple_result, 0, PyArray_Return(nd_discharge));
    PyTuple_SetItem(tuple_result, 1, Py_BuildValue("d", groundwater_component));
    PyTuple_SetItem(tuple_result, 2, PyArray_Return(nd_history));
    return tuple_result;
}


static PyObject *smartcpp_onestep(PyObject *self, PyObject *args) {
    double area_m2, time_delta_sec; // constants
    double c_in_rain, c_in_peva;  // inputs
    double c_p_t, c_p_c, c_p_h, c_p_d, c_p_s, c_p_z, c_p_sk, c_p_fk, c_p_gk, r_p_rk;  // parameters
    double c_s_v_h2o_ove, c_s_v_h2o_dra, c_s_v_h2o_int, c_s_v_h2o_sgw, c_s_v_h2o_dgw, c_s_v_h2o_ly1, c_s_v_h2o_ly2, c_s_v_h2o_ly3, c_s_v_h2o_ly4, c_s_v_h2o_ly5, c_s_v_h2o_ly6, r_s_v_riv;  // states

    if (!PyArg_ParseTuple(args, "dddddddddddddddddddddddddd",
                          &area_m2, &time_delta_sec, &c_in_rain, &c_in_peva,
                          &c_p_t, &c_p_c, &c_p_h, &c_p_d, &c_p_s, &c_p_z, &c_p_sk, &c_p_fk, &c_p_gk, &r_p_rk,
                          &c_s_v_h2o_ove, &c_s_v_h2o_dra, &c_s_v_h2o_int, &c_s_v_h2o_sgw, &c_s_v_h2o_dgw, &c_s_v_h2o_ly1, &c_s_v_h2o_ly2, &c_s_v_h2o_ly3, &c_s_v_h2o_ly4, &c_s_v_h2o_ly5, &c_s_v_h2o_ly6, &r_s_v_riv)) {
        return NULL;
    }

    /* Calculations for the catchment runoff */
    Catchment c;

    c = onestep_catchment(
            area_m2, time_delta_sec, c_in_rain, c_in_peva,
            c_p_t, c_p_c, c_p_h, c_p_d, c_p_s, c_p_z, c_p_sk, c_p_fk, c_p_gk,
            c_s_v_h2o_ove, c_s_v_h2o_dra, c_s_v_h2o_int, c_s_v_h2o_sgw, c_s_v_h2o_dgw, c_s_v_h2o_ly1, c_s_v_h2o_ly2, c_s_v_h2o_ly3, c_s_v_h2o_ly4, c_s_v_h2o_ly5, c_s_v_h2o_ly6);

    /* Calculations for the river routing */
    River r;

    r = onestep_river(
            time_delta_sec,
            c.out_q_h2o_ove + c.out_q_h2o_dra + c.out_q_h2o_int + c.out_q_h2o_sgw + c.out_q_h2o_dgw,
            r_p_rk,
            r_s_v_riv);

    return Py_BuildValue("ddddddddddddddddddd",
                         c.out_aeva, c.out_q_h2o_ove, c.out_q_h2o_dra, c.out_q_h2o_int, c.out_q_h2o_sgw, c.out_q_h2o_dgw, r.out_q_riv,
                         c.s_v_h2o_ove, c.s_v_h2o_dra, c.s_v_h2o_int, c.s_v_h2o_sgw, c.s_v_h2o_dgw,
                         c.s_v_h2o_ly1, c.s_v_h2o_ly2, c.s_v_h2o_ly3, c.s_v_h2o_ly4, c.s_v_h2o_ly5, c.s_v_h2o_ly6,
                         r.s_v_riv);
}


static PyObject *smartcpp_onestep_c(PyObject *self, PyObject *args) {
    double area_m2, time_delta_sec; // constants
    double c_in_rain, c_in_peva;  // inputs
    double c_p_t, c_p_c, c_p_h, c_p_d, c_p_s, c_p_z, c_p_sk, c_p_fk, c_p_gk;  // parameters
    double c_s_v_h2o_ove, c_s_v_h2o_dra, c_s_v_h2o_int, c_s_v_h2o_sgw, c_s_v_h2o_dgw, c_s_v_h2o_ly1, c_s_v_h2o_ly2, c_s_v_h2o_ly3, c_s_v_h2o_ly4, c_s_v_h2o_ly5, c_s_v_h2o_ly6;  // states

    if (!PyArg_ParseTuple(args, "dddddddddddddddddddddddd",
                          &area_m2, &time_delta_sec, &c_in_rain, &c_in_peva,
                          &c_p_t, &c_p_c, &c_p_h, &c_p_d, &c_p_s, &c_p_z, &c_p_sk, &c_p_fk, &c_p_gk,
                          &c_s_v_h2o_ove, &c_s_v_h2o_dra, &c_s_v_h2o_int, &c_s_v_h2o_sgw, &c_s_v_h2o_dgw, &c_s_v_h2o_ly1, &c_s_v_h2o_ly2, &c_s_v_h2o_ly3, &c_s_v_h2o_ly4, &c_s_v_h2o_ly5, &c_s_v_h2o_ly6)) {
        return NULL;
    }

    /* Calculations for the catchment runoff */
    Catchment c;

    c = onestep_catchment(
            area_m2, time_delta_sec, c_in_rain, c_in_peva,
            c_p_t, c_p_c, c_p_h, c_p_d, c_p_s, c_p_z, c_p_sk, c_p_fk, c_p_gk,
            c_s_v_h2o_ove, c_s_v_h2o_dra, c_s_v_h2o_int, c_s_v_h2o_sgw, c_s_v_h2o_dgw, c_s_v_h2o_ly1, c_s_v_h2o_ly2, c_s_v_h2o_ly3, c_s_v_h2o_ly4, c_s_v_h2o_ly5, c_s_v_h2o_ly6);

    return Py_BuildValue("dddddddddddddddddddddd",
                         c.out_aeva, c.out_q_h2o_ove, c.out_q_h2o_dra, c.out_q_h2o_int, c.out_q_h2o_sgw, c.out_q_h2o_dgw,
                         c.s_v_h2o_ove, c.s_v_h2o_dra, c.s_v_h2o_int, c.s_v_h2o_sgw, c.s_v_h2o_dgw,
                         c.s_v_h2o_ly1, c.s_v_h2o_ly2, c.s_v_h2o_ly3, c.s_v_h2o_ly4, c.s_v_h2o_ly5, c.s_v_h2o_ly6,
                         c.pr_eff_rain_to_ove, c.pr_eff_rain_to_dra, c.pr_eff_rain_to_int,
                         c.pr_eff_rain_to_sgw, c.pr_eff_rain_to_dgw);
}


static PyObject *smartcpp_onestep_r(PyObject *self, PyObject *args) {
    double time_delta_sec; // constants
    double r_in_q_riv;  // inputs
    double r_p_rk;  // parameters
    double r_s_v_riv;  // states

    if (!PyArg_ParseTuple(args, "dddd",
                          &time_delta_sec, &r_in_q_riv,
                          &r_p_rk,
                          &r_s_v_riv)) {
        return NULL;
    }

    /* Calculations for the river routing */
    River r;

    r = onestep_river(
            time_delta_sec,
            r_in_q_riv,
            r_p_rk,
            r_s_v_riv);

    return Py_BuildValue("dd",
                         r.out_q_riv,
                         r.s_v_riv);
}

static char smartcpp_docstring[] =
        "This module provides access to the Rainfall-Runoff Model SMART.\n";

static char allsteps_docstring[] =
        "Calculates SMART discharge and groundwater component for all time steps.\n";

static char onestep_docstring[] =
        "Calculates SMART variables for one time step (Catchment Runoff + River Routing).\n";

static char onestep_c_docstring[] =
        "Calculates SMART variables for one time step (Catchment Runoff only).\n";

static char onestep_r_docstring[] =
        "Calculates SMART variables for one time step (River Routing only).\n";

static PyMethodDef smartcpp_methods[] = {
        { "allsteps", smartcpp_allsteps, METH_VARARGS, allsteps_docstring },
        { "onestep", smartcpp_onestep, METH_VARARGS, onestep_docstring },
        { "onestep_c", smartcpp_onestep_c, METH_VARARGS, onestep_c_docstring },
        { "onestep_r", smartcpp_onestep_r, METH_VARARGS, onestep_r_docstring },
        { NULL, NULL, 0, NULL }
};

// For Python 2.x
#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initsmartcpp(void) {
    import_array();
    Py_InitModule3( "smartcpp", smartcpp_methods, smartcpp_docstring );
}
#endif

// For Python 3.x
#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef smartcpp =
        {
                PyModuleDef_HEAD_INIT,
                "smartcpp", /* name of module */
                smartcpp_docstring,          /* module documentation, may be NULL */
                -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
                smartcpp_methods
        };

PyMODINIT_FUNC PyInit_smartcpp(void)
{
    import_array();
    return PyModule_Create(&smartcpp);
}
#endif
