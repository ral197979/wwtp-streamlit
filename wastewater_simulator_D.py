import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd
import os
import math
import graphviz

# --- App Configuration ---
st.set_page_config(layout="wide", page_title="Advanced WWTP Design Simulator")

# --- Core Simulation Classes ---

class UnitConverter:
    """Manages unit conversions to ensure calculations are consistent."""
    def __init__(self, system='US Customary'):
        if system not in ['SI', 'US Customary']:
            raise ValueError("System must be 'SI' or 'US Customary'")
        self.system = system
        # Conversion Factors (Internal unit is always SI)
        self.MGD_to_m3d = 3785.41
        self.MG_to_m3 = 3785.41
        self.gpm_to_m3d = 5.451
        self.kg_to_lbs = 2.20462
        self.m_to_ft = 3.28084
        self.m2_to_ft2 = 10.7639
        self.m_to_in = 39.3701
        self.kw_to_hp = 1.34102
        self.m3h_to_cfm = 0.588578
        self.gfd_to_md = 0.04074
        self.lmh_to_gfd = 0.589 # For MBR Flux
        self.lbs_ft2_yr_to_kg_m2_yr = 4.88243
        self.L_to_gal = 0.264172

    # --- Conversions TO Internal SI Units ---
    def to_internal_flow(self, value):
        return value * self.MGD_to_m3d if self.system == 'US Customary' else value
    def to_internal_volume(self, value):
        return value * self.MG_to_m3 if self.system == 'US Customary' else value
    def to_internal_head(self, value):
        return value / self.m_to_ft if self.system == 'US Customary' else value
    def to_internal_sor(self, value):
        return value * self.gfd_to_md if self.system == 'US Customary' else value
    def to_internal_flux(self, value):
        return value / self.lmh_to_gfd if self.system == 'US Customary' else value # Internal is LMH (L/m2/hr)

    # --- Conversions FROM Internal SI Units ---
    def from_internal_flow(self, value, units='MGD'):
        if self.system == 'US Customary':
            return value / self.MGD_to_m3d if units == 'MGD' else value / self.gpm_to_m3d
        return value
    def from_internal_mass(self, value):
        return value * self.kg_to_lbs if self.system == 'US Customary' else value
    def from_internal_volume(self, value, vol_type='tank'):
        if self.system == 'US Customary':
            return value / self.MG_to_m3 if vol_type == 'tank' else value * self.L_to_gal
        return value
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
    def from_internal_flux(self, value):
        return value * self.lmh_to_gfd if self.system == 'US Customary' else value

    # --- Label Getters for UI ---
    def get_flow_label(self, short=False):
        return ('MGD' if short else 'Flow Rate (MGD)') if self.system == 'US Customary' else ('m³/d' if short else 'Flow Rate (m³/d)')
    def get_volume_label(self, short=False):
        return ('MG' if short else 'Volume (MG)') if self.system == 'US Customary' else ('m³' if short else 'Volume (m³)')
    def get_conc_label(self): return 'Concentration (mg/L)'
    def get_mass_label(self, short=False):
        return ('lbs/day' if short else 'Mass (lbs/day)') if self.system == 'US Customary' else ('kg/day' if short else 'Mass (kg/day)')
    def get_power_label(self): return 'HP' if self.system == 'US Customary' else 'kW'
    def get_head_label(self): return 'ft' if self.system == 'US Customary' else 'm'
    def get_airflow_label(self): return 'SCFM' if self.system == 'US Customary' else 'Sm³/hr'
    def get_area_label(self): return 'ft²' if self.system == 'US Customary' else 'm²'
    def get_sor_label(self): return 'gal/day/ft²' if self.system == 'US Customary' else 'm/d'
    def get_flux_label(self): return 'GFD' if self.system == 'US Customary' else 'LMH'

class Influent:
    """Defines influent characteristics using standard ASM-style fractionation."""
    def __init__(self, flow_rate, cod, tkn, s_po4, s_alk, fractions, units):
        self.base_flow_rate = units.to_internal_flow(flow_rate)
        # COD Fractions
        self.base_s_s = cod * (fractions['cod_readily_biodegradable'] / 100.0)
        self.base_x_s = cod * (fractions['cod_slowly_biodegradable'] / 100.0)
        self.base_s_i = cod * (fractions['cod_soluble_unbiodegradable'] / 100.0)
        self.base_x_i = cod * (fractions['cod_particulate_unbiodegradable'] / 100.0)
        # Nitrogen Fractions
        self.base_s_nh = tkn * (fractions['tkn_ammonia'] / 100.0)
        self.base_s_nd = tkn * (fractions['tkn_soluble_organic'] / 100.0)
        self.base_x_nd = tkn * (fractions['tkn_particulate_organic'] / 100.0)
        # Other components
        self.base_s_no3 = 0
        self.base_s_po4 = s_po4
        self.base_s_alk = s_alk

class BaseWastewaterPlant:
    """Base class for all wastewater treatment plant models."""
    def __init__(self, influent, reactor_volume, was_flow, initial_conditions, units, tech_params, calib_params):
        self.influent = influent
        self.volume = reactor_volume
        self.was_flow = was_flow
        self.y0 = [initial_conditions[key] for key in sorted(initial_conditions.keys())]
        self.technology = "Base"
        self.tech_params = tech_params
        
        # Unpack parameters
        self.temp_params = tech_params['temp_params']
        self.aeration_params = tech_params['aeration_params']

        # Calibration Parameters
        self.mu_H_max_20 = calib_params['mu_H_max_20']
        self.b_H_20 = calib_params['b_H_20']
        self.mu_A_max_20 = calib_params['mu_A_max_20']
        self.b_A_20 = calib_params['b_A_20']
        self.eta_g = calib_params['eta_g']

        # Standard ASM Parameters
        self.Y_H, self.K_S, self.k_h_20, self.K_X = 0.67, 10.0, 3.0, 0.1
        self.Y_A, self.K_NH = 0.24, 1.0
        self.i_N_bm, self.i_N_ss, self.i_P_bm, self.i_P_ss = 0.08, 0.06, 0.02, 0.01
        self.f_p = 0.08
        self.K_OH, self.K_OA, self.K_NO = 0.2, 0.4, 0.5
        self.S_O_sat = 8.5
        self.k_a = 0.05

        # Temperature Correction
        T = self.temp_params['temp_C']
        self.mu_H_max = self.mu_H_max_20 * (1.072 ** (T - 20))
        self.b_H = self.b_H_20 * (1.04 ** (T - 20))
        self.k_h = self.k_h_20 * (1.04 ** (T - 20))
        self.mu_A_max = self.mu_A_max_20 * (1.07 ** (T - 20))
        self.b_A = self.b_A_20 * (1.04 ** (T - 20))

    def get_dynamic_influent(self, t):
        return (self.influent.base_flow_rate, self.influent.base_s_s, self.influent.base_x_s, self.influent.base_s_nh, self.influent.base_s_no3, self.influent.base_s_po4, self.influent.base_s_alk, self.influent.base_s_i, self.influent.base_x_i, self.influent.base_s_nd, self.influent.base_x_nd)

    def model(self, y, t):
        raise NotImplementedError("The 'model' method must be implemented by a subclass.")

    def simulate(self, t):
        self.sol, self.infodict = odeint(self.model, self.y0, t, full_output=True)
        return self.sol

class CASPlant(BaseWastewaterPlant):
    """Models a Conventional Activated Sludge (CAS) plant."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.technology = "CAS"
        self.anoxic_frac = self.aeration_params['anoxic_frac'] / 100.0
        self.aerobic_frac = 1.0 - self.anoxic_frac

    def model(self, y, t):
        y = np.maximum(y, 0)
        S_ALK, S_I, S_ND, S_NH, S_NO3, S_O, S_PO4, S_S, X_AUT, X_BH, X_I, X_ND, X_S = y
        Q_in, S_S_in, X_S_in, S_NH_in, S_NO3_in, S_PO4_in, S_ALK_in, S_I_in, X_I_in, S_ND_in, X_ND_in = self.get_dynamic_influent(t)
        
        epsilon = 1e-8
        aerobic_switch = S_O / (self.K_OH + S_O + epsilon)
        anoxic_switch = (self.K_OH / (self.K_OH + S_O + epsilon)) * (S_NO3 / (self.K_NO + S_NO3 + epsilon))
        
        growth_H_aerobic = self.mu_H_max * (S_S / (self.K_S + S_S + epsilon)) * aerobic_switch * X_BH
        growth_H_anoxic = self.mu_H_max * (S_S / (self.K_S + S_S + epsilon)) * self.eta_g * anoxic_switch * X_BH
        growth_A = self.mu_A_max * (S_NH / (self.K_NH + S_NH + epsilon)) * (S_O / (self.K_OA + S_O + epsilon)) * X_AUT
        hydrolysis = self.k_h * (X_S / (self.K_X * X_BH + X_S + epsilon)) * (aerobic_switch + self.eta_g * anoxic_switch) * X_BH
        ammonification = self.k_a * S_ND * X_BH
        decay_H = self.b_H * X_BH
        decay_A = self.b_A * X_AUT

        # Mass Balances
        Q_eff = Q_in - self.was_flow
        dS_I_dt = (Q_in/self.volume)*(S_I_in - S_I)
        dX_I_dt = (Q_in/self.volume)*X_I_in - (self.was_flow/self.volume)*X_I - (Q_eff/self.volume)*X_I + self.f_p*(decay_H + decay_A)
        dS_ND_dt = (Q_in/self.volume)*(S_ND_in - S_ND) - ammonification
        dX_ND_dt = (Q_in/self.volume)*X_ND_in - (self.was_flow/self.volume)*X_ND - (Q_eff/self.volume)*X_ND
        dS_S_dt = (Q_in/self.volume)*(S_S_in - S_S) - (1/self.Y_H)*(growth_H_aerobic + growth_H_anoxic) + hydrolysis
        dX_S_dt = (Q_in/self.volume)*X_S_in - (self.was_flow/self.volume)*X_S - (Q_eff/self.volume)*X_S - hydrolysis
        dX_BH_dt = growth_H_aerobic + growth_H_anoxic - decay_H - (self.was_flow/self.volume)*X_BH - (Q_eff/self.volume)*X_BH
        dX_AUT_dt = growth_A - decay_A - (self.was_flow/self.volume)*X_AUT - (Q_eff/self.volume)*X_AUT
        
        dS_NH_dt = (Q_in/self.volume)*(S_NH_in - S_NH) + ammonification + self.i_N_ss*hydrolysis + self.i_N_bm*(decay_H + decay_A) - self.i_N_bm*(growth_H_aerobic + growth_H_anoxic) - (1/self.Y_A)*growth_A
        dS_NO3_dt = (Q_in/self.volume)*(S_NO3_in - S_NO3) + (1/self.Y_A)*growth_A - ((1 - self.Y_H)/(2.86 * self.Y_H))*growth_H_anoxic
        dS_PO4_dt = (Q_in/self.volume)*(S_PO4_in - S_PO4) - self.i_P_bm*(growth_H_aerobic + growth_H_anoxic + growth_A) + self.i_P_bm*(decay_H + decay_A) + self.i_P_ss*hydrolysis
        dS_ALK_dt = (Q_in/self.volume)*(S_ALK_in - S_ALK) - (7.14/self.Y_A)*growth_A + (3.57*(1 - self.Y_H)/(2.86*self.Y_H))*growth_H_anoxic
        
        o2_consumption_H = ((1 - self.Y_H)/self.Y_H)*growth_H_aerobic
        o2_consumption_A = ((4.57 - self.Y_A)/self.Y_A)*growth_A
        o2_supply = self.aeration_params['KLa'] * (self.S_O_sat - S_O) * self.aerobic_frac
        dS_O_dt = (Q_in/self.volume)*(0 - S_O) + o2_supply - o2_consumption_H - o2_consumption_A

        self.last_rates = {'AOR_H_g_d': o2_consumption_H*self.volume, 'AOR_A_g_d': o2_consumption_A*self.volume}
        return [dS_ALK_dt, dS_I_dt, dS_ND_dt, dS_NH_dt, dS_NO3_dt, dS_O_dt, dS_PO4_dt, dS_S_dt, dX_AUT_dt, dX_BH_dt, dX_I_dt, dX_ND_dt, dX_S_dt]

class MBRPlant(CASPlant):
    """Models a Membrane Bioreactor (MBR) plant. Inherits from CAS and modifies solids removal."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.technology = "MBR"

    def model(self, y, t):
        # MBR has perfect solids retention, so effluent particulate concentrations are zero.
        # The only way solids leave is via WAS flow.
        y = np.maximum(y, 0)
        S_ALK, S_I, S_ND, S_NH, S_NO3, S_O, S_PO4, S_S, X_AUT, X_BH, X_I, X_ND, X_S = y
        Q_in, S_S_in, X_S_in, S_NH_in, S_NO3_in, S_PO4_in, S_ALK_in, S_I_in, X_I_in, S_ND_in, X_ND_in = self.get_dynamic_influent(t)
        
        epsilon = 1e-8
        aerobic_switch = S_O / (self.K_OH + S_O + epsilon)
        anoxic_switch = (self.K_OH / (self.K_OH + S_O + epsilon)) * (S_NO3 / (self.K_NO + S_NO3 + epsilon))
        
        growth_H_aerobic = self.mu_H_max * (S_S / (self.K_S + S_S + epsilon)) * aerobic_switch * X_BH
        growth_H_anoxic = self.mu_H_max * (S_S / (self.K_S + S_S + epsilon)) * self.eta_g * anoxic_switch * X_BH
        growth_A = self.mu_A_max * (S_NH / (self.K_NH + S_NH + epsilon)) * (S_O / (self.K_OA + S_O + epsilon)) * X_AUT
        hydrolysis = self.k_h * (X_S / (self.K_X * X_BH + X_S + epsilon)) * (aerobic_switch + self.eta_g * anoxic_switch) * X_BH
        ammonification = self.k_a * S_ND * X_BH
        decay_H = self.b_H * X_BH
        decay_A = self.b_A * X_AUT

        # Mass Balances (MODIFIED FOR MBR - NO SOLIDS IN EFFLUENT)
        dS_I_dt = (Q_in/self.volume)*(S_I_in - S_I)
        dX_I_dt = (Q_in/self.volume)*X_I_in - (self.was_flow/self.volume)*X_I + self.f_p*(decay_H + decay_A)
        dS_ND_dt = (Q_in/self.volume)*(S_ND_in - S_ND) - ammonification
        dX_ND_dt = (Q_in/self.volume)*X_ND_in - (self.was_flow/self.volume)*X_ND
        dS_S_dt = (Q_in/self.volume)*(S_S_in - S_S) - (1/self.Y_H)*(growth_H_aerobic + growth_H_anoxic) + hydrolysis
        dX_S_dt = (Q_in/self.volume)*X_S_in - (self.was_flow/self.volume)*X_S - hydrolysis
        dX_BH_dt = growth_H_aerobic + growth_H_anoxic - decay_H - (self.was_flow/self.volume)*X_BH
        dX_AUT_dt = growth_A - decay_A - (self.was_flow/self.volume)*X_AUT
        
        # Soluble components are the same as CAS
        dS_NH_dt = (Q_in/self.volume)*(S_NH_in - S_NH) + ammonification + self.i_N_ss*hydrolysis + self.i_N_bm*(decay_H + decay_A) - self.i_N_bm*(growth_H_aerobic + growth_H_anoxic) - (1/self.Y_A)*growth_A
        dS_NO3_dt = (Q_in/self.volume)*(S_NO3_in - S_NO3) + (1/self.Y_A)*growth_A - ((1 - self.Y_H)/(2.86 * self.Y_H))*growth_H_anoxic
        dS_PO4_dt = (Q_in/self.volume)*(S_PO4_in - S_PO4) - self.i_P_bm*(growth_H_aerobic + growth_H_anoxic + growth_A) + self.i_P_bm*(decay_H + decay_A) + self.i_P_ss*hydrolysis
        dS_ALK_dt = (Q_in/self.volume)*(S_ALK_in - S_ALK) - (7.14/self.Y_A)*growth_A + (3.57*(1 - self.Y_H)/(2.86*self.Y_H))*growth_H_anoxic
        
        o2_consumption_H = ((1 - self.Y_H)/self.Y_H)*growth_H_aerobic
        o2_consumption_A = ((4.57 - self.Y_A)/self.Y_A)*growth_A
        o2_supply = self.aeration_params['KLa'] * (self.S_O_sat - S_O) * self.aerobic_frac
        dS_O_dt = (Q_in/self.volume)*(0 - S_O) + o2_supply - o2_consumption_H - o2_consumption_A

        self.last_rates = {'AOR_H_g_d': o2_consumption_H*self.volume, 'AOR_A_g_d': o2_consumption_A*self.volume}
        return [dS_ALK_dt, dS_I_dt, dS_ND_dt, dS_NH_dt, dS_NO3_dt, dS_O_dt, dS_PO4_dt, dS_S_dt, dX_AUT_dt, dX_BH_dt, dX_I_dt, dX_ND_dt, dX_S_dt]

# NOTE: IFAS and MBBR models are simplified representations for this simulator.
# A full biofilm model is significantly more complex.
class IFASPlant(CASPlant):
    """Simplified model for an Integrated Fixed-Film Activated Sludge (IFAS) plant."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.technology = "IFAS"
        # IFAS adds biofilm growth, which enhances nitrification. We simulate this
        # by providing a bonus to the autotrophic growth rate.
        self.ifas_bonus = 1.5 # 50% bonus to autotroph growth rate
        self.mu_A_max_ifas = self.mu_A_max * self.ifas_bonus

    def model(self, y, t):
        # We override the autotrophic growth rate with the IFAS bonus rate.
        original_mu_A_max = self.mu_A_max
        self.mu_A_max = self.mu_A_max_ifas
        result = super().model(y, t)
        self.mu_A_max = original_mu_A_max # Reset after calculation
        return result

class MBBRPlant(CASPlant):
    """Simplified model for a Moving Bed Biofilm Reactor (MBBR) plant."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.technology = "MBBR"
        # MBBR is a pure biofilm system with no sludge recycle.
        # We simulate this by having no WAS flow and letting solids leave in the effluent.
        # The SRT is very long, controlled by biofilm sloughing, not wasting.
        self.was_flow = 0 # No WAS flow in MBBR
        # High biomass is maintained on media. We simulate this with a lower decay rate.
        self.b_H_mbbr = self.b_H * 0.5
        self.b_A_mbbr = self.b_A * 0.5

    def model(self, y, t):
        # Override decay rates for MBBR model
        original_b_H = self.b_H
        original_b_A = self.b_A
        self.b_H = self.b_H_mbbr
        self.b_A = self.b_A_mbbr
        result = super().model(y, t)
        # Reset after calculation
        self.b_H = original_b_H
        self.b_A = original_b_A
        return result

# --- Plotting and Diagram Functions ---

def plot_and_save_results(t, results, initial_conditions, units, technology, filename="simulation_results.png"):
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 8))
    labels = sorted(initial_conditions.keys())
    for i, label in enumerate(labels):
        ax.plot(t, results[:, i], label=label)
    ax.set_title(f'{technology} Plant Simulation Results (Units: {units.system})', fontsize=16)
    ax.set_xlabel('Time (days)', fontsize=12)
    ax.set_ylabel(units.get_conc_label(), fontsize=12)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    ax.grid(True)
    fig.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(filename)
    plt.close(fig)

def create_pfd(data, units, technology):
    dot = graphviz.Digraph(comment='Wastewater Treatment Plant PFD', graph_attr={'splines': 'ortho', 'rankdir': 'LR'})
    dot.attr('node', shape='box', style='rounded,filled', fillcolor='lightblue2', fontname='Helvetica')
    dot.attr('edge', fontsize='10', fontname='Helvetica')

    dot.node('IN', 'Influent')
    dot.node('PC', f"Primary Clarifier\n{data['pc_dims']}")
    
    if technology == "MBR":
        dot.node('AT', f"Aeration Tank\n{data['reactor_dims']}")
        dot.node('MT', f"Membrane Tank\n{data['membrane_dims']}")
        dot.edge('AT', 'MT', label=f"MLSS: {data['mlss']:.0f} mg/L")
        dot.edge('MT', 'EFF', label=f"Permeate\nQ: {data['q_eff']:.2f} {units.get_flow_label(short=True)}")
        dot.edge('MT', 'WAS', label=f"WAS\nQ: {data['q_was']:.2f} {units.get_flow_label(short=True)}")
    elif technology == "MBBR":
        dot.node('AT', f"MBBR\n{data['reactor_dims']}\nMedia Fill: {data['media_fill']}%")
        dot.node('SC', f"Effluent Screen")
        dot.edge('AT', 'SC', label=f"MLSS: {data['mlss']:.0f} mg/L")
        dot.edge('SC', 'EFF', label=f"Q: {data['q_eff']:.2f} {units.get_flow_label(short=True)}")
    else: # CAS or IFAS
        at_label = "Aeration Tank"
        if technology == "IFAS":
            at_label = f"IFAS Tank\n(Media Fill: {data['media_fill']}%)"
        dot.node('AT', f"{at_label}\n{data['reactor_dims']}")
        dot.node('SC', f"Secondary Clarifier\n{data['sc_dims']}")
        dot.edge('AT', 'SC', label=f"MLSS: {data['mlss']:.0f} mg/L")
        dot.edge('SC', 'EFF', label=f"Q: {data['q_eff']:.2f} {units.get_flow_label(short=True)}")
        dot.edge('SC', 'WAS', label=f"WAS\nQ: {data['q_was']:.2f} {units.get_flow_label(short=True)}")
        dot.edge('SC', 'AT', label=f"RAS\nQ: {data['q_ras']:.2f} {units.get_flow_label(short=True)}", headport='s', tailport='s')

    dot.edge('IN', 'PC', label=f"Q: {data['q_in']:.2f} {units.get_flow_label(short=True)}\nCOD: {data['inf_cod']:.0f} mg/L")
    dot.edge('PC', 'AT')
    dot.node('WAS', 'To Sludge Handling')
    dot.node('EFF', 'Effluent')
    return dot

# --- UI and Main App Logic ---

st.title("Advanced Wastewater Treatment Plant Design Simulator")

# --- Session State Initialization ---
if 'run_complete' not in st.session_state:
    st.session_state['run_complete'] = False

# --- Sidebar for Inputs ---
with st.sidebar:
    st.title("⚙️ Plant Setup")
    
    technology_type = st.selectbox("Select Treatment Technology", ["CAS", "MBR", "IFAS", "MBBR"])
    
    selected_units = st.selectbox("Select Unit System", ["US Customary", "SI"], index=0)
    units = UnitConverter(selected_units)
    
    st.header("Influent Characterization")
    influent_flow = st.number_input(units.get_flow_label(), value=1.0 if selected_units == 'US Customary' else 3785.0, min_value=0.1)
    total_cod = st.number_input("Total Influent COD (mg/L)", value=450.0, min_value=10.0)
    total_tkn = st.number_input("Total Influent TKN (mg-N/L)", value=40.0, min_value=1.0)
    s_po4_in = st.number_input("Influent Phosphorus (S_PO4) (mg-P/L)", value=6.0, min_value=0.0)
    s_alk_in = st.number_input("Influent Alkalinity (mg-CaCO3/L)", value=250.0, min_value=0.0)

    with st.expander("COD & TKN Fractions"):
        cod_frac_ss = st.slider("Readily Biodegradable COD (%)", 0, 100, 20)
        cod_frac_xs = st.slider("Slowly Biodegradable COD (%)", 0, 100, 50)
        cod_frac_si = st.slider("Soluble Unbiodegradable COD (%)", 0, 100, 5)
        cod_frac_xi = 100 - cod_frac_ss - cod_frac_xs - cod_frac_si
        if cod_frac_xi < 0: st.error("COD fractions sum to > 100%.")
        st.info(f"Particulate Unbiodegradable (X_I): {cod_frac_xi:.1f}%")
        
        tkn_frac_nh = st.slider("Ammonia Nitrogen (%)", 0, 100, 70)
        tkn_frac_snd = st.slider("Soluble Organic N (%)", 0, 100, 10)
        tkn_frac_xnd = 100 - tkn_frac_nh - tkn_frac_snd
        if tkn_frac_xnd < 0: st.error("TKN fractions sum to > 100%.")
        st.info(f"Particulate Organic N (X_ND): {tkn_frac_xnd:.1f}%")

    fractions = {'cod_readily_biodegradable': cod_frac_ss, 'cod_slowly_biodegradable': cod_frac_xs, 'cod_soluble_unbiodegradable': cod_frac_si, 'cod_particulate_unbiodegradable': cod_frac_xi, 'tkn_ammonia': tkn_frac_nh, 'tkn_soluble_organic': tkn_frac_snd, 'tkn_particulate_organic': tkn_frac_xnd}
    
    st.header("Process Design Parameters")
    # PARAMETER 1: Set a challenging but realistic temperature
    temp_C = st.slider("Wastewater Temperature (°C)", 5, 35, 15, help="Lower temperatures require longer SRTs.")
    
    # PARAMETER 2 (CORRECTED): Increase SRT for robust nitrification to prevent washout
    target_srt = st.number_input("Target Solids Retention Time (SRT) (days)", value=20.0, min_value=2.0, step=1.0)
    
    # PARAMETER 3 (CORRECTED): Increase reactor volume for longer HRT and more stability
    reactor_volume_display = st.number_input(units.get_volume_label(), value=1.2 if selected_units == 'US Customary' else 4542.0, min_value=0.1)
    
    # PARAMETER 4 (CORRECTED): Increase anoxic fraction for better denitrification
    anoxic_frac = st.slider("Anoxic Zone Fraction (%)", 0, 80, 40, help="Fraction of the reactor volume that is anoxic (for denitrification).")
    
    kla = st.slider("Oxygen Transfer Coefficient (KLa) (d⁻¹)", 10, 300, 240)

    # Technology-specific inputs
    if technology_type == "MBR":
        membrane_flux = st.number_input(f"Design Membrane Flux ({units.get_flux_label()})", value=15.0 if selected_units == 'US Customary' else 25.0)
    if technology_type in ["IFAS", "MBBR"]:
        media_fill_frac = st.slider("Media Fill Fraction (%)", 10, 70, 40)
        media_ssa = st.number_input("Media Specific Surface Area (m²/m³)", value=500)

    with st.expander("Edit Advanced Design Criteria"):
        primary_sor = st.number_input(f"Primary Clarifier SOR ({units.get_sor_label()})", value=800.0 if selected_units == 'US Customary' else 32.6)
        if technology_type in ["CAS", "IFAS"]:
            secondary_sor = st.number_input(f"Secondary Clarifier SOR ({units.get_sor_label()})", value=600.0 if selected_units == 'US Customary' else 24.4)
        reactor_depth = st.number_input(f"Reactor Depth ({units.get_head_label()})", value=15.0 if selected_units == 'US Customary' else 4.5)
        clarifier_depth = st.number_input(f"Clarifier Depth ({units.get_head_label()})", value=12.0 if selected_units == 'US Customary' else 3.6)

    st.header("Model Calibration")
    with st.expander("Edit Kinetic Parameters"):
        mu_H_max_20 = st.slider("μ_H_max @ 20°C (d⁻¹)", 2.0, 8.0, 4.5)
        b_H_20 = st.slider("b_H @ 20°C (d⁻¹)", 0.1, 0.8, 0.4)
        mu_A_max_20 = st.slider("μ_A_max @ 20°C (d⁻¹)", 0.3, 1.5, 0.8)
        b_A_20 = st.slider("b_A @ 20°C (d⁻¹)", 0.05, 0.3, 0.15)
        eta_g = st.slider("Anoxic Growth Factor (η_g)", 0.4, 1.0, 0.8)

    calib_params = {'mu_H_max_20': mu_H_max_20, 'b_H_20': b_H_20, 'mu_A_max_20': mu_A_max_20, 'b_A_20': b_A_20, 'eta_g': eta_g}

# --- Main Application Tabs ---
sim_tab, design_tab, pfd_tab, report_tab = st.tabs(["📊 Process Simulation", "📐 Preliminary Design", "➡️ Process Diagram", "📋 Design Summary"])

with sim_tab:
    col1, col2 = st.columns([1, 2])
    with col1:
        st.subheader(f"Run {technology_type} Simulation")
        if st.button("Run Simulation", type="primary"):
            reactor_volume_internal = units.to_internal_volume(reactor_volume_display)
            was_flow_internal = reactor_volume_internal / target_srt if target_srt > 0 else 0

            influent = Influent(influent_flow, total_cod, total_tkn, s_po4_in, s_alk_in, fractions, units)
            
            initial_conditions = {
                'Alkalinity (S_ALK)': s_alk_in, 'Soluble Inert (S_I)': influent.base_s_i,
                'Soluble Organic N (S_ND)': influent.base_s_nd, 'Ammonia (S_NH)': influent.base_s_nh,
                'Nitrate (S_NO3)': 5.0, 'Dissolved Oxygen (S_O)': 2.0,
                'Soluble Phosphorus (S_PO4)': s_po4_in, 'Soluble Substrate (S_S)': influent.base_s_s,
                'Autotrophic Biomass (X_AUT)': 50.0, 'Heterotrophic Biomass (X_BH)': 200.0,
                'Particulate Inert (X_I)': influent.base_x_i, 'Particulate Organic N (X_ND)': influent.base_x_nd,
                'Particulate Substrate (X_S)': influent.base_x_s,
            }
            
            tech_params = {
                'temp_params': {'temp_C': temp_C},
                'aeration_params': {'KLa': kla, 'anoxic_frac': anoxic_frac}
            }
            
            plant_map = {
                "CAS": CASPlant,
                "MBR": MBRPlant,
                "IFAS": IFASPlant,
                "MBBR": MBBRPlant
            }
            PlantClass = plant_map[technology_type]

            with st.spinner(f'Running {technology_type} simulation...'):
                plant = PlantClass(influent, reactor_volume_internal, was_flow_internal, initial_conditions, units, tech_params, calib_params)
                t_sim = np.linspace(0, max(50, 3 * target_srt), 500)
                results = plant.simulate(t_sim)

                if np.isnan(results).any() or np.isinf(results).any():
                    st.session_state['run_complete'] = False
                    st.error("Simulation failed. Try increasing SRT, KLa, or reactor volume.")
                else:
                    plot_and_save_results(t_sim, results, initial_conditions, units, technology_type, "simulation_results.png")
                    st.session_state['run_complete'] = True
                    df = pd.DataFrame(results, columns=sorted(initial_conditions.keys()))
                    df.insert(0, "Time (days)", t_sim)
                    
                    st.session_state['technology_type'] = technology_type
                    st.session_state['results_data'] = {'df': df, 'steady_state': df.iloc[--1]}
                    st.session_state['units_obj'] = units
                    st.session_state['influent_obj'] = influent
                    st.session_state['was_flow_internal'] = was_flow_internal
                    st.session_state['reactor_volume_internal'] = reactor_volume_internal
                    st.session_state['influent_cod'] = total_cod
                    st.session_state['influent_tkn'] = total_tkn
                    
                    final_rates = plant.last_rates
                    st.session_state['oxygen_demand_kg_d'] = (final_rates.get('AOR_H_g_d', 0) + final_rates.get('AOR_A_g_d', 0)) / 1000.0

                    design_criteria = {'reactor_depth_user': reactor_depth, 'primary_sor': primary_sor}
                    if technology_type in ["CAS", "IFAS"]:
                        design_criteria['secondary_sor'] = secondary_sor
                        design_criteria['clarifier_depth_user'] = clarifier_depth
                    if technology_type == "MBR":
                        design_criteria['membrane_flux'] = membrane_flux
                    if technology_type in ["IFAS", "MBBR"]:
                        design_criteria['media_fill_frac'] = media_fill_frac
                        design_criteria['media_ssa'] = media_ssa
                    st.session_state['design_criteria'] = design_criteria
            st.rerun()

    if st.session_state.get('run_complete'):
        st.subheader("Simulation Summary")
        steady_state = st.session_state['results_data']['steady_state']
        st.metric("Effluent Ammonia (mg-N/L)", f"{steady_state['Ammonia (S_NH)']:.2f}")
        st.metric("Effluent Nitrate (mg-N/L)", f"{steady_state['Nitrate (S_NO3)']:.2f}")
        st.metric("MLSS (mg/L)", f"{steady_state['Heterotrophic Biomass (X_BH)'] + steady_state['Autotrophic Biomass (X_AUT)']:.0f}")
        
        csv = st.session_state['results_data']['df'].to_csv(index=False).encode('utf-8')
        st.download_button("Download Full Results (CSV)", csv, 'simulation_results.csv', 'text/csv')
    else:
        st.info("Adjust parameters and click 'Run Simulation'.")

    with col2:
        st.subheader("Results Plot")
        if st.session_state.get('run_complete') and os.path.exists("simulation_results.png"):
            st.image("simulation_results.png", use_container_width=True)
        else:
            st.info("Results will be displayed here.")

with design_tab:
    st.header(f"Preliminary {st.session_state.get('technology_type', 'CAS')} Unit Sizing")
    if st.session_state.get('run_complete'):
        units = st.session_state['units_obj']
        influent = st.session_state['influent_obj']
        criteria = st.session_state['design_criteria']
        tech = st.session_state['technology_type']

        st.subheader("Primary Treatment")
        pc_area_internal = influent.base_flow_rate / units.to_internal_sor(criteria['primary_sor'])
        st.metric(label="Primary Clarifier Surface Area", value=f"{units.from_internal_area(pc_area_internal):.0f} {units.get_area_label()}")

        st.subheader("Secondary Treatment")
        if tech in ["CAS", "IFAS"]:
            sc_area_internal = influent.base_flow_rate / units.to_internal_sor(criteria['secondary_sor'])
            st.metric(label="Secondary Clarifier Area", value=f"{units.from_internal_area(sc_area_internal):.0f} {units.get_area_label()}")
        elif tech == "MBR":
            flux_internal = units.to_internal_flux(criteria['membrane_flux']) # convert to LMH
            flow_m3_hr = influent.base_flow_rate / 24
            membrane_area = flow_m3_hr / flux_internal if flux_internal > 0 else 0
            st.metric(label="Required Membrane Area", value=f"{units.from_internal_area(membrane_area):.0f} {units.get_area_label()}")
        
        if tech in ["IFAS", "MBBR"]:
            media_vol = st.session_state['reactor_volume_internal'] * (criteria['media_fill_frac'] / 100.0)
            total_ssa = media_vol * criteria['media_ssa']
            st.metric(label="Total Media Surface Area", value=f"{units.from_internal_area(total_ssa):.0f} {units.get_area_label()}")

        st.subheader("Aeration System")
        AOR_kg_hr = st.session_state['oxygen_demand_kg_d'] / 24.0
        st.metric("Actual Oxygen Rate (AOR)", f"{AOR_kg_hr:.1f} kg O₂/hr")
    else:
        st.info("Run a simulation to see design calculations.")

def get_diagram_data(units, criteria, tech):
    q_in_internal = st.session_state['influent_obj'].base_flow_rate
    q_was_internal = st.session_state.get('was_flow_internal', 0)
    q_eff_internal = q_in_internal - q_was_internal
    steady_state = st.session_state['results_data']['steady_state']
    mlss = steady_state['Heterotrophic Biomass (X_BH)'] + steady_state['Autotrophic Biomass (X_AUT)'] + steady_state['Particulate Inert (X_I)'] + steady_state['Particulate Substrate (X_S)']
    
    reactor_depth_internal = units.to_internal_head(criteria['reactor_depth_user'])
    reactor_area = st.session_state['reactor_volume_internal'] / reactor_depth_internal
    reactor_w = math.sqrt(reactor_area)
    reactor_dims_ft = f"{units.from_internal_length(reactor_w):.1f}ft W x {units.from_internal_length(reactor_w):.1f}ft L x {criteria['reactor_depth_user']:.1f}ft D"
    reactor_dims_m = f"{reactor_w:.1f}m W x {reactor_w:.1f}m L x {reactor_depth_internal:.1f}m D"
    
    pc_area_internal = q_in_internal / units.to_internal_sor(criteria['primary_sor'])
    pc_diam_internal = math.sqrt(4 * pc_area_internal / math.pi)

    data = {
        'q_in': units.from_internal_flow(q_in_internal), 'q_was': units.from_internal_flow(q_was_internal, 'gpm'),
        'q_ras': units.from_internal_flow(q_in_internal * 0.5), 'q_eff': units.from_internal_flow(q_eff_internal),
        'mlss': mlss, 'inf_cod': st.session_state['influent_cod'], 'inf_tkn': st.session_state['influent_tkn'],
        'eff_nh3': steady_state['Ammonia (S_NH)'], 'eff_no3': steady_state['Nitrate (S_NO3)'],
        'reactor_dims': reactor_dims_ft if units.system == 'US Customary' else reactor_dims_m,
        'pc_dims': f"Ø{units.from_internal_length(pc_diam_internal):.1f}ft" if units.system == 'US Customary' else f"Ø{pc_diam_internal:.1f}m",
    }
    
    if tech in ["CAS", "IFAS"]:
        sc_area_internal = q_in_internal / units.to_internal_sor(criteria['secondary_sor'])
        sc_diam_internal = math.sqrt(4 * sc_area_internal / math.pi)
        data['sc_dims'] = f"Ø{units.from_internal_length(sc_diam_internal):.1f}ft" if units.system == 'US Customary' else f"Ø{sc_diam_internal:.1f}m"
    if tech == "MBR":
        flux_internal = units.to_internal_flux(criteria['membrane_flux'])
        membrane_area = (q_in_internal / 24) / flux_internal if flux_internal > 0 else 0
        data['membrane_dims'] = f"Area: {units.from_internal_area(membrane_area):.0f} {units.get_area_label()}"
    if tech in ["IFAS", "MBBR"]:
        data['media_fill'] = criteria['media_fill_frac']

    return data

with pfd_tab:
    st.header("Process Flow Diagram (PFD)")
    if st.session_state.get('run_complete'):
        units = st.session_state['units_obj']
        criteria = st.session_state['design_criteria']
        tech = st.session_state['technology_type']
        pfd_data = get_diagram_data(units, criteria, tech)
        
        try:
            pfd_graph = create_pfd(pfd_data, units, tech)
            st.graphviz_chart(pfd_graph)
        except Exception as e:
            st.error(f"Could not generate PFD. Ensure Graphviz is installed.")
            st.error(f"Details: {e}")
    else:
        st.info("Run a simulation to generate the PFD.")

with report_tab:
    st.header("Design Summary Report")
    if st.session_state.get('run_complete'):
        units = st.session_state['units_obj']
        criteria = st.session_state['design_criteria']
        tech = st.session_state['technology_type']
        data = get_diagram_data(units, criteria, tech)
        steady_state = st.session_state['results_data']['steady_state']

        st.subheader(f"Summary for {tech} Plant")
        col1, col2 = st.columns(2)
        with col1:
            st.write("**Influent & Effluent**")
            st.text(f"Influent Flow: {data['q_in']:.2f} {units.get_flow_label(short=True)}")
            st.text(f"Effluent NH3: {data['eff_nh3']:.2f} mg-N/L")
            st.text(f"Effluent NO3: {data['eff_no3']:.2f} mg-N/L")
        with col2:
            st.write("**Process Parameters**")
            st.text(f"Target SRT: {target_srt:.1f} days at {temp_C}°C")
            st.text(f"MLSS: {data['mlss']:.0f} mg/L")
            if tech != "MBBR":
                st.text(f"WAS Flow: {data['q_was']:.1f} {units.get_flow_label(short=True).replace('MGD', 'gpm')}")
        
        st.subheader("Unit Sizing")
        st.text(f"Primary Clarifier: {data['pc_dims']}")
        st.text(f"Reactor: {data['reactor_dims']}")
        if tech in ["CAS", "IFAS"]:
            st.text(f"Secondary Clarifier: {data['sc_dims']}")
        elif tech == "MBR":
            st.text(f"Membranes: {data['membrane_dims']}")
    else:
        st.info("Run a simulation to generate the design summary report.")