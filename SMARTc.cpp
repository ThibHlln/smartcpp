#include <Python.h>

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
    // initialise runoff pathways
    double overland_flow, drain_flow, inter_flow, shallow_flow, deep_flow;

    if (excess_rain >= 0.0) {  // excess rainfall available for runoff and infiltration
        // actual evapotranspiration = potential evapotranspiration
        aeva += c_in_peva;

        // calculate surface runoff using quick runoff parameter H and relative soil moisture content
        double h_prime = c_p_h * (lvl_total_start / c_p_z);
        overland_flow = h_prime * excess_rain;  // excess rainfall contribution to quick surface runoff store
        excess_rain -= overland_flow;  // remainder that infiltrates

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
        drain_flow = c_p_d * excess_rain;  // sat. excess contribution (if not 0) to quick interflow runoff store
        inter_flow = (1.0 - c_p_d) * excess_rain;  // sat. excess contribution (if not 0) to slow interflow runoff store

        // calculate leak from soil layers (i.e. piston flow becoming active during rainfall events)
        double s_prime = c_p_s * (lvl_total_start / c_p_z);
        // leak to interflow
        for (int i = 1; i <= 6; i++) {  // soil moisture outflow reducing exponentially downwards
            double leak_interflow = list_lvl_lyr[i] * pow(s_prime, i);
            if (leak_interflow < list_lvl_lyr[i]) {
                inter_flow += leak_interflow;  // soil moisture outflow contribution to slow interflow runoff store
                list_lvl_lyr[i] -= leak_interflow;
            }
        }
        // leak to shallow groundwater flow
        shallow_flow = 0.0;
        for (int i = 1; i <= 6; i++) {  // soil moisture outflow reducing linearly downwards
            double leak_shallow_flow = list_lvl_lyr[i] * (s_prime / i);
            if (leak_shallow_flow < list_lvl_lyr[i]) {
                shallow_flow += leak_shallow_flow;  // soil moisture outflow contribution to slow shallow GW runoff store
                list_lvl_lyr[i] -= leak_shallow_flow;
            }
        }
        // leak to deep groundwater flow
        deep_flow = 0.0;
        for (int i = 1; i <= 6; i++) {  // soil moisture outflow reducing exponentially upwards
            double leak_deep_flow = list_lvl_lyr[i] * pow(s_prime, 7 - i);
            if (leak_deep_flow < list_lvl_lyr[i]) {
                deep_flow += leak_deep_flow;  // soil moisture outflow contribution to slow deep GW runoff store
                list_lvl_lyr[i] -= leak_deep_flow;
            }
        }
    } else {  // no excess rainfall (i.e. potential evapotranspiration not satisfied by available rainfall)
        overland_flow = 0.0;  // no soil moisture contribution to quick overland flow runoff store
        drain_flow = 0.0;  // no soil moisture contribution to quick drain flow runoff store
        inter_flow = 0.0;  // no soil moisture contribution to quick + leak interflow runoff store
        shallow_flow = 0.0;  // no soil moisture contribution to shallow groundwater flow runoff store
        deep_flow = 0.0;  // no soil moisture contribution to deep groundwater flow runoff store

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
    c.s_v_h2o_ove = c_s_v_h2o_ove + (overland_flow / 1e3 * area_m2) - (c.out_q_h2o_ove * time_delta_sec);  // [m3] - [m3]
    if (c.s_v_h2o_ove < 0.0)
        c.s_v_h2o_ove = 0.0;
    // route drain flow (quick interflow runoff)
    c.out_q_h2o_dra = c_s_v_h2o_dra / c_p_sk;  // [m3/s]
    c.s_v_h2o_dra = c_s_v_h2o_dra + (drain_flow / 1e3 * area_m2) - (c.out_q_h2o_dra * time_delta_sec);  // [m3] - [m3]
    if (c.s_v_h2o_dra < 0.0)
        c.s_v_h2o_dra = 0.0;
    // route interflow (slow interflow runoff)
    c.out_q_h2o_int = c_s_v_h2o_int / c_p_fk;  // [m3/s]
    c.s_v_h2o_int = c_s_v_h2o_int + (inter_flow / 1e3 * area_m2) - (c.out_q_h2o_int * time_delta_sec);  // [m3] - [m3]
    if (c.s_v_h2o_int < 0.0)
        c.s_v_h2o_int = 0.0;
    // route shallow groundwater flow (slow shallow GW runoff)
    c.out_q_h2o_sgw = c_s_v_h2o_sgw / c_p_gk;  // [m3/s]
    c.s_v_h2o_sgw = c_s_v_h2o_sgw + (shallow_flow / 1e3 * area_m2) - (c.out_q_h2o_sgw * time_delta_sec);  // [m3] - [m3]
    if (c.s_v_h2o_sgw < 0.0)
        c.s_v_h2o_sgw = 0.0;
    // route deep groundwater flow (slow deep GW runoff)
    c.out_q_h2o_dgw = c_s_v_h2o_dgw / c_p_gk;  // [m3/s]
    c.s_v_h2o_dgw = c_s_v_h2o_dgw + (deep_flow / 1e3 * area_m2) - (c.out_q_h2o_dgw * time_delta_sec);  // [m3] - [m3]
    if (c.s_v_h2o_dgw < 0.0)
        c.s_v_h2o_dgw = 0.0;

    // update soil layer levels
    c.s_v_h2o_ly1 = list_lvl_lyr[1] / 1000 * area_m2;
    c.s_v_h2o_ly2 = list_lvl_lyr[2] / 1000 * area_m2;
    c.s_v_h2o_ly3 = list_lvl_lyr[3] / 1000 * area_m2;
    c.s_v_h2o_ly4 = list_lvl_lyr[4] / 1000 * area_m2;
    c.s_v_h2o_ly5 = list_lvl_lyr[5] / 1000 * area_m2;
    c.s_v_h2o_ly6 = list_lvl_lyr[6] / 1000 * area_m2;

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


static PyObject *SMARTc_onestep( PyObject *self, PyObject *args ) {
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

static char SMARTc_docstring[] =
        "This module provides access to the Rainfall-Runoff Model SMART.\n";

static char onestep_docstring[] =
        "Calculates SMART variables for one time step.\n";

static PyMethodDef SMARTc_methods[] = {
        { "onestep", SMARTc_onestep, METH_VARARGS, onestep_docstring },
        { NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC initSMARTc(void) {
    Py_InitModule3( "SMARTc", SMARTc_methods, SMARTc_docstring );
}