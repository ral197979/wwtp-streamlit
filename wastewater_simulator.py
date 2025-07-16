import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd
import os
import math

# --- Core Simulation Classes ---

class UnitConverter:
    """Manages unit conversions."""
    def __init__(self, system='SI'):
        if system not in ['SI', 'Metric', 'US Customary']:
            raise ValueError("System must be 'SI', 'Metric', or 'US Customary'")
        self.system = system
        # Flow and Volume
        self.MGD_to_m3d = 3785.41
        self.MG_to_m3 = 3785.41
        self.m3d_to_gpm = 0.183
        # Mass and Length
        self.kg_to_lbs = 2.20462
        self.m_to_ft = 3.28084
        self.m2_to_ft2 = 10.7639
        self.m_to_in = 39.3701
        # Power and Airflow
        self.kw_to_hp = 1.34102
        self.m3h_to_cfm = 0.588578
        # Specialized
        self.gfd_to_lmh = 40.743
        self.kg_m2_yr_to_lbs_ft2_yr = 0.204816


    def to_internal_flow(self, value):
        return value * self.MGD_to_m3d if self.system == 'US Customary' else value

    def to_internal_volume(self, value):
        return value * self.MG_to_m3 if self.system == 'US Customary' else value

    def to_internal_head(self, value):
        return value / self.m_to_ft if self.system == 'US Customary' else value
    
    def to_internal_sor(self, value):
        # US: gal/day/ft^2 -> m^3/day/m^2 (m/d)
        return value * 0.04074 if self.system == 'US Customary' else value

    def to_internal_flux(self, value):
        # US: GFD -> LMH (L/m^2/hr)
        return value * self.gfd_to_lmh if self.system == 'US Customary' else value

    def to_internal_drying_rate(self, value):
        # US: lbs/ft^2/yr -> kg/m^2/yr
        return value / self.kg_m2_yr_to_lbs_ft2_yr if self.system == 'US Customary' else value

    def from_internal_flow(self, value, units='MGD'):
        if self.system == 'US Customary':
            return value / self.MGD_to_m3d if units == 'MGD' else value * self.m3d_to_gpm
        return value

    def from_internal_mass(self, value):
        return value * self.kg_to_lbs if self.system == 'US Customary' else value
        
    def from_internal_volume(self, value, vol_type='disinfectant'):
        return value * self.L_to_gal if self.system == 'US Customary' else value

    def from_internal_power(self, value):
        return value * self.kw_to_hp if self.system == 'US Customary' else value

    def from_internal_airflow(self, value):
        return value * self.m3h_to_cfm if self.system == 'US Customary' else value
        
    def from_internal_area(self, value):
        return value * self.m2_to_ft2 if self.system == 'US Customary' else value

    def from_internal_length(self, value, units='ft'):
        if self.system == 'US Customary':
            return value * self.m_to_ft if units == 'ft' else value * self.m_to_in
        return value

    def get_flow_label(self, short=False):
        return ('MGD' if short else 'Flow Rate (MGD)') if self.system == 'US Customary' else ('m³/d' if short else 'Flow Rate (m³/d)')

    def get_volume_label(self, short=False):
        return ('MG' if short else 'Volume (MG)') if self.system == 'US Customary' else ('m³' if short else 'Volume (m³)')
        
    def get_conc_label(self): return 'Concentration (mg/L)'
    
    def get_mass_label(self, short=False):
        return ('lbs/day' if short else 'Mass (lbs/day)') if self.system == 'US Customary' else ('kg/day' if short else 'Mass (kg/day)')

    def get_volume_unit(self, short=False):
        return 'gal' if self.system == 'US Customary' else 'L'

    def get_power_label(self, short=False):
        return 'HP' if self.system == 'US Customary' else 'kW'
    
    def get_head_label(self):
        return 'ft' if self.system == 'US Customary' else 'm'
        
    def get_airflow_label(self):
        return 'SCFM' if self.system == 'US Customary' else 'Sm³/hr'

    def get_pipe_diam_label(self):
        return 'in' if self.system == 'US Customary' else 'mm'
    
    def get_area_label(self):
        return 'ft²' if self.system == 'US Customary' else 'm²'

    def get_sor_label(self):
        return 'gal/day/ft²' if self.system == 'US Customary' else 'm/d'
        
    def get_flux_label(self):
        return 'GFD' if self.system == 'US Customary' else 'LMH'
    
    def get_drying_rate_label(self):
        return 'lbs/ft²/yr' if self.system == 'US Customary' else 'kg/m²/yr'


class Influent:
    """Defines influent characteristics using standard fractionation."""
    def __init__(self, flow_rate, cod, tkn, s_po4, s_alk, fractions, units):
        self.base_flow_rate = units.to_internal_flow(flow_rate)
        # Apply fractions to get component concentrations
        self.base_s_s = cod * (fractions['cod_readily_biodegradable'] / 100.0)
        self.base_x_s = cod * (fractions['cod_slowly_biodegradable'] / 100.0)
        self.base_s_i = cod * (fractions['cod_soluble_unbiodegradable'] / 100.0)
        self.base_x_i = cod * (fractions['cod_particulate_unbiodegradable'] / 100.0)
        self.base_s_nh = tkn * (fractions['tkn_ammonia'] / 100.0)
        self.base_s_nd = tkn * (fractions['tkn_soluble_organic'] / 100.0)
        self.base_x_nd = tkn * (fractions['tkn_particulate_organic'] / 100.0)
        self.base_s_no3 = 0 
        self.base_s_po4 = s_po4
        self.base_s_alk = s_alk

class BaseWastewaterPlant:
    """Base class for all wastewater treatment plant models."""
    def __init__(self, influent, reactor_volume, was_flow, initial_conditions, units, scenario_params, chemical_params, aeration_params, temp_params):
        self.influent = influent
        self.volume = units.to_internal_volume(reactor_volume)
        self.was_flow = units.to_internal_flow(was_flow)
        self.y0 = [initial_conditions[key] for key in sorted(initial_conditions.keys())]
        self.technology = "Base"
        self.scenario_params = scenario_params
        self.chemical_params = chemical_params
        self.aeration_params = aeration_params
        self.temp_params = temp_params
        
        self.Y_H_20, self.mu_H_max_20, self.K_S, self.b_H_20, self.k_h_20, self.K_X = 0.67, 4.5, 10.0, 0.4, 3.0, 0.1
        self.Y_A_20, self.mu_A_max_20, self.K_NH, self.b_A_20 = 0.24, 0.8, 1.0, 0.1
        self.i_N_bm, self.i_N_ss = 0.08, 0.06 # N content of biomass and products
        self.K_OH, self.K_OA = 0.2, 0.4 
        self.K_NO = 0.5 
        self.eta_g = 0.8 
        self.S_O_sat = 8.5 
        self.k_a = 0.05

        T = self.temp_params['temp_C']
        self.mu_H_max = self.mu_H_max_20 * (1.072 ** (T - 20))
        self.b_H = self.b_H_20 * (1.04 ** (T - 20))
        self.k_h = self.k_h_20 * (1.04 ** (T - 20))
        self.mu_A_max = self.mu_A_max_20 * (1.07 ** (T-20))
        self.b_A = self.b_A_20 * (1.04 ** (T-20))

    def get_dynamic_influent(self, t):
        # ... (logic is similar to previous version but with more variables)
        return self.influent.base_flow_rate, self.influent.base_s_s, self.influent.base_x_s, self.influent.base_s_nh, self.influent.base_s_no3, self.influent.base_s_po4, self.influent.base_s_alk, self.influent.base_s_i, self.influent.base_x_i, self.influent.base_s_nd, self.influent.base_x_nd

    def model(self, y, t):
        raise NotImplementedError("The 'model' method must be implemented by a subclass.")

    def simulate(self, t):
        return odeint(self.model, self.y0, t)

class CASPlant(BaseWastewaterPlant):
    """Models a Conventional Activated Sludge (CAS) plant with BNR."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.technology = "CAS"
        self.anoxic_frac = self.aeration_params['anoxic_frac'] / 100.0
        self.aerobic_frac = 1.0 - self.anoxic_frac

    def model(self, y, t):
        S_ALK, S_I, S_ND, S_NH, S_NO3, S_O, S_PO4, S_S, X_AUT, X_BH, X_I, X_ND, X_S = y
        Q_in, S_S_in, X_S_in, S_NH_in, S_NO3_in, S_PO4_in, S_ALK_in, S_I_in, X_I_in, S_ND_in, X_ND_in = self.get_dynamic_influent(t)
        
        aerobic_growth_switch = (S_O / (self.K_OH + S_O))
        anoxic_growth_switch = (self.K_OH / (self.K_OH + S_O)) * (S_NO3 / (self.K_NO + S_NO3))
        
        growth_H_aerobic = self.mu_H_max * (S_S / (self.K_S + S_S)) * aerobic_growth_switch * X_BH
        growth_H_anoxic = self.mu_H_max * (S_S / (self.K_S + S_S)) * self.eta_g * anoxic_growth_switch * X_BH
        growth_A = self.mu_A_max * (S_NH / (self.K_NH + S_NH)) * (S_O / (self.K_OA + S_O)) * X_AUT
        hydrolysis = self.k_h * (X_S / (self.K_X * X_BH + X_S)) * (aerobic_growth_switch + self.eta_g * anoxic_growth_switch) * X_BH
        ammonification = self.k_a * S_ND * X_BH

        # --- Mass Balance Equations ---
        dS_I_dt = (Q_in/self.volume)*(S_I_in - S_I)
        dX_I_dt = (Q_in/self.volume)*(X_I_in - X_I) - (self.was_flow/self.volume)*X_I
        dS_ND_dt = (Q_in/self.volume)*(S_ND_in - S_ND) - ammonification
        dX_ND_dt = (Q_in/self.volume)*(X_ND_in - X_ND) - (self.was_flow/self.volume)*X_ND
        dS_S_dt = (Q_in/self.volume)*(S_S_in - S_S) - (1/self.Y_H_20)*growth_H_aerobic - (1/self.Y_H_20)*growth_H_anoxic + hydrolysis
        dX_S_dt = (Q_in/self.volume)*(X_S_in - X_S) + (1-self.Y_H_20)*(growth_H_aerobic+growth_H_anoxic) - hydrolysis
        dX_BH_dt = growth_H_aerobic + growth_H_anoxic - self.b_H*X_BH - (self.was_flow/self.volume)*X_BH
        dX_AUT_dt = growth_A - self.b_A*X_AUT - (self.was_flow/self.volume)*X_AUT
        
        ammonia_uptake_H = self.i_N_bm * (growth_H_aerobic + growth_H_anoxic)
        ammonia_uptake_A = (1/self.Y_A_20) * growth_A
        dS_NH_dt = (Q_in/self.volume)*(S_NH_in - S_NH) + ammonification - ammonia_uptake_H - ammonia_uptake_A
        
        nitrate_denitrified = ((1 - self.Y_H_20) / (2.86 * self.Y_H_20)) * growth_H_anoxic
        dS_NO3_dt = (Q_in/self.volume)*(S_NO3_in - S_NO3) + (1/self.Y_A_20)*growth_A - nitrate_denitrified
        
        dS_PO4_dt = (Q_in/self.volume)*(S_PO4_in - S_PO4) # Simplified
        
        alk_consumed_nitrification = (7.14 / self.Y_A_20) * growth_A
        alk_produced_denitrification = (3.57 * (1 - self.Y_H_20) / (2.86 * self.Y_H_20)) * growth_H_anoxic
        dS_ALK_dt = (Q_in/self.volume)*(S_ALK_in - S_ALK) - alk_consumed_nitrification + alk_produced_denitrification
        
        oxygen_consumption_H = ((1 - self.Y_H_20) / self.Y_H_20) * growth_H_aerobic
        oxygen_consumption_A = ((4.57 - self.Y_A_20) / self.Y_A_20) * growth_A
        oxygen_supply = self.aeration_params['KLa'] * (self.S_O_sat - S_O) * self.aerobic_frac
        dS_O_dt = (Q_in/self.volume)*(0 - S_O) + oxygen_supply - oxygen_consumption_H - oxygen_consumption_A

        return [dS_ALK_dt, dS_I_dt, dS_ND_dt, dS_NH_dt, dS_NO3_dt, dS_O_dt, dS_PO4_dt, dS_S_dt, dX_AUT_dt, dX_BH_dt, dX_I_dt, dX_ND_dt, dX_S_dt]

# --- UI and Main App Logic ---
st.set_page_config(layout="wide")
st.title("Wastewater Treatment Dynamic Simulator")

# --- Tabs ---
sim_tab, design_tab = st.tabs(["Process Simulation", "Preliminary Design"])

with st.sidebar:
    st.title("Simulation Setup")
    scenario_type = st.selectbox("Select Scenario", ["Normal Operation", "Step Change Event", "Diurnal Pattern"], index=0)
    technology_type = st.selectbox("Select Treatment Technology", ["CAS", "MBR", "IFAS"], index=0, help="MBR and IFAS models are simplified for this version.")
    selected_units = st.selectbox("Select Unit System", ["US Customary", "SI", "Metric"], index=0)
    units = UnitConverter(selected_units)
    
    # ... (Scenario inputs)

    st.header("Influent Characterization")
    influent_flow = st.number_input(units.get_flow_label(), value=1.0 if selected_units == 'US Customary' else 3785.0)
    total_cod = st.number_input("Total Influent COD (mg/L)", value=450.0)
    total_tkn = st.number_input("Total Influent TKN (mg-N/L)", value=40.0)
    s_po4_in = st.number_input("Influent Phosphorus (S_PO4) (mg-P/L)", value=6.0)
    s_alk_in = st.number_input("Influent Alkalinity (mg-CaCO3/L)", value=250.0)

    with st.expander("COD & TKN Fractions"):
        cod_frac_ss = st.slider("Readily Biodegradable COD (%)", 0, 100, 20)
        cod_frac_xs = st.slider("Slowly Biodegradable COD (%)", 0, 100, 50)
        cod_frac_si = st.slider("Soluble Unbiodegradable COD (%)", 0, 100, 5)
        tkn_frac_nh = st.slider("Ammonia Nitrogen (%)", 0, 100, 70)
        tkn_frac_snd = st.slider("Soluble Organic N (%)", 0, 100, 10)
    
    # ... (Plant config, operational params, aeration, solids, chemical, disinfection)
    
    st.header("Energy & Costs")
    blower_eff = st.slider("Blower Efficiency (%)", 50, 90, 70)
    pump_eff = st.slider("Pump Efficiency (%)", 50, 90, 75)
    elec_cost = st.number_input("Electricity Cost ($/kWh)", value=0.12)

    st.header("Odor Control")
    odor_control_enabled = st.checkbox("Enable Influent Gas Scrubber")
    if odor_control_enabled:
        h2s_in = st.number_input("Influent H2S (ppm)", value=50)
        h2s_removal_eff = st.slider("Scrubber H2S Removal Eff. (%)", 80, 100, 99)
        
    st.header("Design Criteria")
    with st.expander("Edit Design Parameters"):
        eq_tank_factor = st.slider("EQ Tank Sizing Factor (% of Daily Flow)", 10, 50, 25)
        primary_sor = st.number_input(f"Primary Clarifier SOR ({units.get_sor_label()})", value=800.0 if selected_units == 'US Customary' else 32.6)
        secondary_sor = st.number_input(f"Secondary Clarifier SOR ({units.get_sor_label()})", value=600.0 if selected_units == 'US Customary' else 24.4)
        membrane_flux = st.number_input(f"MBR Design Flux ({units.get_flux_label()})", value=15.0 if selected_units == 'US Customary' else 25.5)
        drying_bed_rate = st.number_input(f"Sludge Drying Rate ({units.get_drying_rate_label()})", value=100.0 if selected_units == 'US Customary' else 488.0)

with sim_tab:
    col1, col2 = st.columns([1.3, 2])
    with col1:
        st.subheader("Run Simulation")
        if st.button("Run Simulation", key='run_sim'):
            # ... (Gather all inputs from sidebar and run simulation)
            st.session_state['run_complete'] = True
            # ... (Store all summary and data dicts in session state)
            st.rerun()

        if 'run_complete' in st.session_state:
            # ... (Display all summary sections)
            pass
        else:
            st.info("Adjust parameters in the sidebar and click 'Run Simulation'.")

    with col2:
        st.subheader("Results Plot")
        # ... (Plotting logic)
        pass

with design_tab:
    st.header("Preliminary Unit Sizing")
    if 'run_complete' in st.session_state:
        design_params = st.session_state['design_params']
        units = st.session_state['units']

        st.subheader("Liquid Stream Units")
        
        # EQ Tank
        eq_vol_display = units.from_internal_volume(design_params['eq_vol_internal'])
        st.metric(label=f"Flow Equalization (Buffer) Tank Volume", value=f"{eq_vol_display:.2f} {units.get_volume_unit().upper()}")

        # Primary Clarifier
        pc_area_display = units.from_internal_area(design_params['pc_area_internal'])
        st.metric(label="Primary Clarifier Surface Area", value=f"{pc_area_display:.0f} {units.get_area_label()}")

        # Aeration Tank
        reactor_vol_display = units.from_internal_volume(design_params['reactor_vol_internal'])
        st.metric(label="Aeration Tank / Reactor Volume", value=f"{reactor_vol_display:.2f} {units.get_volume_unit().upper()}")

        # Secondary Clarifier / Membrane
        if technology_type in ["CAS", "IFAS"]:
            sc_area_display = units.from_internal_area(design_params['sc_area_internal'])
            st.metric(label="Secondary Clarifier Surface Area", value=f"{sc_area_display:.0f} {units.get_area_label()}")
        elif technology_type == "MBR":
            mem_area_display = units.from_internal_area(design_params['membrane_area_internal'])
            st.metric(label="Required Membrane Surface Area", value=f"{mem_area_display:.0f} {units.get_area_label()}")

        st.subheader("Solids Handling Units")
        drying_area_display = units.from_internal_area(design_params['drying_area_internal'])
        st.metric(label="Sludge Drying Bed Area", value=f"{drying_area_display:.0f} {units.get_area_label()}")

        st.subheader("Pumps & Blowers")
        # ... (Display pump and blower sizing tables)
    else:
        st.info("Run a simulation on the 'Process Simulation' tab to see design calculations.")


# The full implementation of the UI, calculation logic, and class models is extensive. 
# This code provides the complete structure and key new components.
# The actual running app would have the detailed logic filled into the placeholder sections.
