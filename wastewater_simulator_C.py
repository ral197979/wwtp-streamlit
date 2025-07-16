import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd
import os
import math
import graphviz

# --- App Configuration ---
st.set_page_config(layout="wide")

# --- Core Simulation Classes ---

class UnitConverter:
    """Manages unit conversions."""
    def __init__(self, system='SI'):
        if system not in ['SI', 'Metric', 'US Customary']:
            raise ValueError("System must be 'SI', 'Metric', or 'US Customary'")
        self.system = system
        # Conversion Factors
        self.MGD_to_m3d = 3785.41
        self.MG_to_m3 = 3785.41
        self.m3d_to_gpm = 0.183
        self.kg_to_lbs = 2.20462
        self.m_to_ft = 3.28084
        self.m2_to_ft2 = 10.7639
        self.m_to_in = 39.3701
        self.kw_to_hp = 1.34102
        self.m3h_to_cfm = 0.588578
        self.gfd_to_lmh = 40.743
        self.kg_m2_yr_to_lbs_ft2_yr = 0.204816
        self.L_to_gal = 0.264172

    def to_internal_flow(self, value):
        return value * self.MGD_to_m3d if self.system == 'US Customary' else value
    def to_internal_volume(self, value):
        return value * self.MG_to_m3 if self.system == 'US Customary' else value
    def to_internal_head(self, value):
        return value / self.m_to_ft if self.system == 'US Customary' else value
    def to_internal_sor(self, value):
        return value / 24.54 # US: gal/day/ft² to m/d. SI: m³/m²/d is already m/d
    def to_internal_flux(self, value):
        return value * self.gfd_to_lmh if self.system == 'US Customary' else value
    def to_internal_drying_rate(self, value):
        return value / self.kg_m2_yr_to_lbs_ft2_yr if self.system == 'US Customary' else value
    def from_internal_flow(self, value, units='MGD'):
        if self.system == 'US Customary':
            return value / self.MGD_to_m3d if units == 'MGD' else value * self.m3d_to_gpm
        return value
    def from_internal_mass(self, value):
        return value * self.kg_to_lbs if self.system == 'US Customary' else value
    def from_internal_volume(self, value, vol_type='disinfectant'):
        if self.system == 'US Customary':
            return value * self.L_to_gal if vol_type == 'disinfectant' else value / self.MG_to_m3
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
    def __init__(self, influent, reactor_volume, was_flow, initial_conditions, units, scenario_params, chemical_params, aeration_params, temp_params, calib_params):
        self.influent = influent
        self.volume = reactor_volume 
        self.was_flow = was_flow 
        self.y0 = [initial_conditions[key] for key in sorted(initial_conditions.keys())]
        self.technology = "Base"
        self.scenario_params = scenario_params
        self.chemical_params = chemical_params
        self.aeration_params = aeration_params
        self.temp_params = temp_params
        
        self.mu_H_max_20 = calib_params['mu_H_max_20']
        self.b_H_20 = calib_params['b_H_20']
        self.mu_A_max_20 = calib_params['mu_A_max_20']
        self.b_A_20 = calib_params['b_A_20']
        self.eta_g = calib_params['eta_g']

        self.Y_H_20, self.K_S, self.k_h_20, self.K_X = 0.67, 10.0, 3.0, 0.1
        self.Y_A_20, self.K_NH = 0.24, 1.0
        self.i_N_bm, self.i_N_ss = 0.08, 0.06 
        self.i_P_bm, self.i_P_ss = 0.02, 0.01 
        self.f_p = 0.2 
        self.K_OH, self.K_OA = 0.2, 0.4 
        self.K_NO = 0.5 
        self.S_O_sat = 8.5 
        self.k_a = 0.05

        T = self.temp_params['temp_C']
        self.mu_H_max = self.mu_H_max_20 * (1.072 ** (T - 20))
        self.b_H = self.b_H_20 * (1.04 ** (T - 20))
        self.k_h = self.k_h_20 * (1.04 ** (T - 20))
        self.mu_A_max = self.mu_A_max_20 * (1.07 ** (T-20))
        self.b_A = self.b_A_20 * (1.04 ** (T-20))

    def get_dynamic_influent(self, t):
        flow, s_s, x_s, s_nh, s_no3, s_po4, s_alk, s_i, x_i, s_nd, x_nd = (
            self.influent.base_flow_rate, self.influent.base_s_s, 
            self.influent.base_x_s, self.influent.base_s_nh, 
            self.influent.base_s_no3, self.influent.base_s_po4, 
            self.influent.base_s_alk, self.influent.base_s_i, 
            self.influent.base_x_i, self.influent.base_s_nd, 
            self.influent.base_x_nd
        )
        return flow, s_s, x_s, s_nh, s_no3, s_po4, s_alk, s_i, x_i, s_nd, x_nd

    def model(self, y, t):
        raise NotImplementedError("The 'model' method must be implemented by a subclass.")

    def simulate(self, t):
        self.sol, self.infodict = odeint(self.model, self.y0, t, full_output=True)
        return self.sol

class CASPlant(BaseWastewaterPlant):
    """Models a Conventional Activated Sludge (CAS) plant with BNR."""
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
        
        aerobic_growth_switch = (S_O / (self.K_OH + S_O + epsilon))
        anoxic_growth_switch = (self.K_OH / (self.K_OH + S_O + epsilon)) * (S_NO3 / (self.K_NO + S_NO3 + epsilon))
        
        growth_H_aerobic = self.mu_H_max * (S_S / (self.K_S + S_S + epsilon)) * aerobic_growth_switch * X_BH
        growth_H_anoxic = self.mu_H_max * (S_S / (self.K_S + S_S + epsilon)) * self.eta_g * anoxic_growth_switch * X_BH
        growth_A = self.mu_A_max * (S_NH / (self.K_NH + S_NH + epsilon)) * (S_O / (self.K_OA + S_O + epsilon)) * X_AUT
        hydrolysis = self.k_h * (X_S / (self.K_X * X_BH + X_S + epsilon)) * (aerobic_growth_switch + self.eta_g * anoxic_growth_switch) * X_BH
        ammonification = self.k_a * S_ND * X_BH
        decay_H = self.b_H * X_BH
        decay_A = self.b_A * X_AUT

        dS_I_dt = (Q_in/self.volume)*(S_I_in - S_I)
        dX_I_dt = (Q_in/self.volume)*(X_I_in - X_I) + self.f_p * (decay_H + decay_A) - (self.was_flow/self.volume)*X_I
        dS_ND_dt = (Q_in/self.volume)*(S_ND_in - S_ND) - ammonification
        dX_ND_dt = (Q_in/self.volume)*(X_ND_in - X_ND) - (self.was_flow/self.volume)*X_ND
        dS_S_dt = (Q_in/self.volume)*(S_S_in - S_S) - (1/self.Y_H_20)*growth_H_aerobic - (1/self.Y_H_20)*growth_H_anoxic + hydrolysis
        dX_S_dt = (Q_in/self.volume)*(X_S_in - X_S) - hydrolysis - (self.was_flow/self.volume)*X_S
        dX_BH_dt = growth_H_aerobic + growth_H_anoxic - decay_H - (self.was_flow/self.volume)*X_BH
        dX_AUT_dt = growth_A - decay_A - (self.was_flow/self.volume)*X_AUT
        
        ammonia_from_hydrolysis = self.i_N_ss * hydrolysis
        ammonia_from_decay = self.i_N_bm * (decay_H + decay_A)
        ammonia_uptake_H = self.i_N_bm * (growth_H_aerobic + growth_H_anoxic)
        ammonia_uptake_A = (1/self.Y_A_20) * growth_A
        dS_NH_dt = (Q_in/self.volume)*(S_NH_in - S_NH) + ammonification + ammonia_from_hydrolysis + ammonia_from_decay - ammonia_uptake_H - ammonia_uptake_A
        
        nitrate_denitrified = ((1 - self.Y_H_20) / (2.86 * self.Y_H_20)) * growth_H_anoxic
        dS_NO3_dt = (Q_in/self.volume)*(S_NO3_in - S_NO3) + (1/self.Y_A_20)*growth_A - nitrate_denitrified
        
        p_uptake = self.i_P_bm * (growth_H_aerobic + growth_H_anoxic + growth_A)
        p_release = (self.i_P_bm * (decay_H + decay_A)) + (self.i_P_ss * hydrolysis)
        dS_PO4_dt = (Q_in/self.volume)*(S_PO4_in - S_PO4) - p_uptake + p_release
        
        alk_consumed_nitrification = (7.14 / self.Y_A_20) * growth_A
        alk_produced_denitrification = (3.57 * (1 - self.Y_H_20) / (2.86 * self.Y_H_20)) * growth_H_anoxic
        dS_ALK_dt = (Q_in/self.volume)*(S_ALK_in - S_ALK) - alk_consumed_nitrification + alk_produced_denitrification
        
        oxygen_consumption_H = ((1 - self.Y_H_20) / self.Y_H_20) * growth_H_aerobic
        oxygen_consumption_A = ((4.57 - self.Y_A_20) / self.Y_A_20) * growth_A
        oxygen_supply = self.aeration_params['KLa'] * (self.S_O_sat - S_O) * self.aerobic_frac
        dS_O_dt = (Q_in/self.volume)*(0 - S_O) + oxygen_supply - oxygen_consumption_H - oxygen_consumption_A

        self.last_rates = {
            'AOR_H_g_d': oxygen_consumption_H * self.volume,
            'AOR_A_g_d': oxygen_consumption_A * self.volume
        }

        return [dS_ALK_dt, dS_I_dt, dS_ND_dt, dS_NH_dt, dS_NO3_dt, dS_O_dt, dS_PO4_dt, dS_S_dt, dX_AUT_dt, dX_BH_dt, dX_I_dt, dX_ND_dt, dX_S_dt]

def plot_and_save_results(t, results, initial_conditions, units, technology, filename="simulation_results.png"):
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 8))
    labels = sorted(initial_conditions.keys())
    for i, label in enumerate(labels):
        ax.plot(t, results[:, i], label=label)
    ax.set_title(f'{technology} Plant Simulation (Units: {units.system})')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel(units.get_conc_label())
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.grid(True)
    fig.tight_layout()
    plt.savefig(filename)
    plt.close(fig)

def create_pfd(data, units):
    dot = graphviz.Digraph(comment='Wastewater Treatment Plant PFD', graph_attr={'splines': 'ortho', 'rankdir': 'LR'})
    dot.attr('node', shape='box', style='rounded,filled', fillcolor='lightblue2', fontname='Helvetica')
    dot.attr('edge', fontsize='10', fontname='Helvetica')

    dot.node('IN', 'Influent')
    # ADDITION: EQ Tank and Pumps
    dot.node('EQ', f"Equalization Tank\n{data['eq_dims']}")
    dot.node('EQ_Pumps', 'EQ Pumps')
    dot.node('PC', f"Primary Clarifier\n{data['pc_dims']}")
    dot.node('AT', f"Aeration Tank\n{data['reactor_dims']}")
    dot.node('SC', f"Secondary Clarifier\n{data['sc_dims']}")
    dot.node('EFF', 'Effluent')
    dot.node('WAS', 'To Sludge Handling')
    
    # ADDITION: Update edges for new EQ system
    dot.edge('IN', 'EQ', label=f"Q: {data['q_in']:.2f} {units.get_flow_label(short=True)}\nCOD: {data['inf_cod']:.0f} mg/L\nTKN: {data['inf_tkn']:.0f} mg/L")
    dot.edge('EQ', 'EQ_Pumps')
    dot.edge('EQ_Pumps', 'PC')
    dot.edge('PC', 'AT')
    dot.edge('AT', 'SC', label=f"Q (ML): {data['q_ml']:.2f} {units.get_flow_label(short=True)}\nMLSS: {data['mlss']:.0f} mg/L")
    dot.edge('SC', 'EFF', label=f"Q: {data['q_eff']:.2f} {units.get_flow_label(short=True)}\nNH₃: {data['eff_nh3']:.2f} mg/L\nNO₃: {data['eff_no3']:.2f} mg/L")
    dot.edge('SC', 'WAS', label=f"Q: {data['q_was']:.2f} {units.get_flow_label(short=True)}")
    dot.edge('SC', 'AT', label=f"RAS\nQ: {data['q_ras']:.2f} {units.get_flow_label(short=True)}", headport='s', tailport='s')

    return dot

def create_pid(data, units):
    """
    Creates a P&ID with styling to mimic AutoCAD Plant 3D.
    """
    pid = graphviz.Digraph(
        'P&ID',
        graph_attr={
            'splines': 'ortho',
            'rankdir': 'LR',
            'nodesep': '1.0',
            'fontname': 'Arial'
        }
    )
    pid.attr('node', style='filled', fillcolor='white', color='black', fontname='Arial', fontsize='10', penwidth='1.5')
    pid.attr('edge', fontname='Arial', fontsize='9', color='black', penwidth='2.0')

    # --- Main Equipment ---
    # ADDITION: EQ Tank
    pid.node('EQ', f"EQUALIZATION TANK\nTK-100\n{data['eq_dims']}", shape='box', height='2.5', width='3.5')
    pid.node('PC', f"PRIMARY CLARIFIER\nCL-101\n{data['pc_dims']}", shape='box', height='2.5', width='3.5')
    pid.node('AT', f"AERATION TANK\nTK-101\n{data['reactor_dims']}", shape='box', height='2.5', width='3.5')
    pid.node('SC', f"SECONDARY CLARIFIER\nCL-102\n{data['sc_dims']}", shape='box', height='2.5', width='2.5')
    pid.node('EFF', 'TO EFFLUENT\nDISCHARGE', shape='rarrow', style='filled', fillcolor='gray90', color='black')

    # --- Rotating Equipment ---
    # ADDITION: EQ Pumps
    with pid.subgraph(name='cluster_p0') as s:
        s.attr(style='invis')
        s.node('P0', 'P-100', shape='circle', fixedsize='true', width='0.8')
        s.node('P0_label', 'EQ\nPUMPS', shape='plain', fontsize='9')
    with pid.subgraph(name='cluster_p1') as s:
        s.attr(style='invis')
        s.node('P1', 'P-101', shape='circle', fixedsize='true', width='0.8')
        s.node('P1_label', 'INFLUENT\nPUMP', shape='plain', fontsize='9')
    with pid.subgraph(name='cluster_b1') as s:
        s.attr(style='invis')
        s.node('B1', 'B-101', shape='circle', fixedsize='true', width='0.8')
        s.node('B1_label', 'AERATION\nBLOWER', shape='plain', fontsize='9')
    with pid.subgraph(name='cluster_p2') as s:
        s.attr(style='invis')
        s.node('P2', 'P-102', shape='circle', fixedsize='true', width='0.8')
        s.node('P2_label', 'RAS\nPUMP', shape='plain', fontsize='9')
    with pid.subgraph(name='cluster_p3') as s:
        s.attr(style='invis')
        s.node('P3', 'P-103', shape='circle', fixedsize='true', width='0.8')
        s.node('P3_label', 'WAS\nPUMP', shape='plain', fontsize='9')

    # --- Instrumentation & Controls ---
    pid.node('LT100', 'LT\n100', shape='circle', fixedsize='true', width='0.6') # EQ Level Transmitter
    pid.node('AT1', 'AT\n101', shape='circle', fixedsize='true', width='0.6') # Analyzer Transmitter
    pid.node('LT1', 'LT\n101', shape='circle', fixedsize='true', width='0.6') # Level Transmitter
    pid.node('V1', 'FCV-101', shape='invtriangle', style='filled', fillcolor='black', fixedsize='true', width='0.3', label='')

    # --- Piping & Signal Lines ---
    pid.edge('P1', 'EQ', label=f" Q: {data['q_in']:.2f} {units.get_flow_label(short=True)}", tailport='e', headport='w')
    pid.edge('EQ', 'P0', tailport='s', headport='n')
    pid.edge('P0', 'PC', tailport='e', headport='w')
    pid.edge('PC', 'AT', tailport='e', headport='w') # Assuming influent pump P1 now feeds EQ tank
    pid.edge('AT', 'SC', tailport='e', headport='w')
    pid.edge('B1', 'V1', tailport='e', headport='n')
    pid.edge('V1', 'AT', label=f" AIR: {data['airflow']:.0f} {units.get_airflow_label()}", arrowhead='normal', tailport='s', headport='s')
    pid.edge('SC', 'P2', tailport='s', headport='n')
    pid.edge('P2', 'AT', label=f" RAS: {data['q_ras']:.2f} {units.get_flow_label(short=True)}", tailport='s', headport='s')
    pid.edge('SC', 'P3', tailport='s', headport='n')
    pid.edge('SC', 'EFF', tailport='e', headport='w')

    # Instrument Lines
    pid.edge('EQ', 'LT100', penwidth='1.0', arrowhead='none', headport='w') # EQ Level
    pid.edge('AT', 'AT1', penwidth='1.0', arrowhead='none', headport='w')
    pid.edge('AT', 'LT1', penwidth='1.0', arrowhead='none', headport='w')

    # Control Signals
    pid.edge('LT100', 'P0', style='dashed', penwidth='1.0', arrowhead='open', label=' LSL/LSH', color='gray50')
    pid.edge('AT1', 'B1', style='dashed', penwidth='1.0', arrowhead='open', label=' DO Signal (AIC)', color='gray50')

    return pid

# --- UI and Main App Logic ---

st.title("Wastewater Treatment Dynamic Simulator")

# --- Tabs ---
sim_tab, design_tab, pfd_tab, pid_tab, external_tab = st.tabs(["Process Simulation", "Preliminary Design", "Process Diagram", "P&ID", "External Simulator"])

with st.sidebar:
    st.title("Simulation Setup")
    scenario_type = st.selectbox("Select Scenario", ["Normal Operation"], index=0, help="Dynamic scenarios are simplified in this version.")
    technology_type = st.selectbox("Select Treatment Technology", ["CAS"], index=0, help="Other models are for demonstration only.")
    selected_units = st.selectbox("Select Unit System", ["US Customary", "SI"], index=0)
    units = UnitConverter(selected_units)
    
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
        cod_frac_xi = 100 - cod_frac_ss - cod_frac_xs - cod_frac_si
        st.info(f"Particulate Unbiodegradable (X_I): {cod_frac_xi:.1f}%")
        tkn_frac_nh = st.slider("Ammonia Nitrogen (%)", 0, 100, 70)
        tkn_frac_snd = st.slider("Soluble Organic N (%)", 0, 100, 10)
        tkn_frac_xnd = 100 - tkn_frac_nh - tkn_frac_snd
        st.info(f"Particulate Organic N (X_ND): {tkn_frac_xnd:.1f}%")

    fractions = {
        'cod_readily_biodegradable': cod_frac_ss, 'cod_slowly_biodegradable': cod_frac_xs,
        'cod_soluble_unbiodegradable': cod_frac_si, 'cod_particulate_unbiodegradable': cod_frac_xi,
        'tkn_ammonia': tkn_frac_nh, 'tkn_soluble_organic': tkn_frac_snd,
        'tkn_particulate_organic': tkn_frac_xnd
    }
    
    st.header("Process Design")
    # MODIFICATION: Increased default SRT to a very robust value for nitrification
    target_srt = st.number_input("Target Solids Retention Time (SRT) (days)", value=20.0, min_value=1.0, step=0.5)
    reactor_volume_display = st.number_input(units.get_volume_label(), value=1.0 if selected_units == 'US Customary' else 3785.0)
    eq_tank_factor = st.slider("EQ Tank Sizing Factor (% of Daily Flow)", 10, 50, 25)
    temp_C = st.slider("Wastewater Temperature (°C)", 5, 35, 20)
    anoxic_frac = st.slider("Anoxic Zone Fraction (%)", 0, 80, 30)
    # MODIFICATION: Increased default KLa to a very robust value for nitrification
    kla = st.slider("Oxygen Transfer Coefficient (KLa) (d⁻¹)", 10, 300, 250)

    with st.expander("Edit Design Criteria"):
        pump_eff = st.slider("Pump Efficiency (%)", 50, 90, 75)
        eq_depth = st.number_input(f"EQ Tank Depth ({units.get_head_label()})", value=15.0 if selected_units == 'US Customary' else 4.5)
        eq_pump_tdh = st.number_input(f"EQ Pump Total Dynamic Head (TDH) ({units.get_head_label()})", value=30.0 if selected_units == 'US Customary' else 9.0)
        primary_sor = st.number_input(f"Primary Clarifier SOR ({units.get_sor_label()})", value=800.0 if selected_units == 'US Customary' else 32.6)
        secondary_sor = st.number_input(f"Secondary Clarifier SOR ({units.get_sor_label()})", value=600.0 if selected_units == 'US Customary' else 24.4)
        reactor_depth = st.number_input(f"Reactor Depth ({units.get_head_label()})", value=15.0 if selected_units == 'US Customary' else 4.5)
        clarifier_depth = st.number_input(f"Clarifier Depth ({units.get_head_label()})", value=12.0 if selected_units == 'US Customary' else 3.6)

    st.header("Model Calibration")
    with st.expander("Edit Kinetic Parameters"):
        mu_H_max_20 = st.slider("Heterotroph Max. Growth Rate (d⁻¹)", 2.0, 8.0, 4.5, help="μ_H_max @ 20°C")
        b_H_20 = st.slider("Heterotroph Decay Rate (d⁻¹)", 0.1, 0.8, 0.4, help="b_H @ 20°C")
        mu_A_max_20 = st.slider("Autotroph Max. Growth Rate (d⁻¹)", 0.3, 1.5, 0.8, help="μ_A_max @ 20°C")
        b_A_20 = st.slider("Autotroph Decay Rate (d⁻¹)", 0.05, 0.3, 0.1, help="b_A @ 20°C")
        eta_g = st.slider("Anoxic Growth Factor (η_g)", 0.4, 1.0, 0.8, help="Correction factor for heterotroph growth in anoxic vs. aerobic conditions.")

    calib_params = {'mu_H_max_20': mu_H_max_20, 'b_H_20': b_H_20, 'mu_A_max_20': mu_A_max_20, 'b_A_20': b_A_20, 'eta_g': eta_g}

with sim_tab:
    col1, col2 = st.columns([1.3, 2])
    with col1:
        st.subheader("Run Simulation")
        if st.button("Run Simulation", key='run_sim', type="primary"):
            reactor_volume_internal = units.to_internal_volume(reactor_volume_display)
            was_flow_internal = reactor_volume_internal / target_srt

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

            plant_args = (influent, reactor_volume_internal, was_flow_internal, initial_conditions, units, 
                          {'type': scenario_type}, {}, 
                          {'KLa': kla, 'anoxic_frac': anoxic_frac}, {'temp_C': temp_C}, calib_params)
            
            with st.spinner('Running simulation...'):
                plant = CASPlant(*plant_args)
                t_sim = np.linspace(0, 50, 500)
                results = plant.simulate(t_sim)

                if np.isnan(results).any() or np.isinf(results).any():
                    st.session_state['run_complete'] = False
                    st.error("Simulation failed due to numerical instability. Please adjust input parameters (e.g., SRT, KLa, or reactor volume) and try again.")
                else:
                    plot_and_save_results(t_sim, results, initial_conditions, units, technology_type)
                    st.session_state['run_complete'] = True
                    df = pd.DataFrame(results, columns=sorted(initial_conditions.keys()))
                    df.insert(0, "Time (days)", t_sim)
                    
                    st.session_state['results_data'] = {'df': df, 'steady_state': df.iloc[-1]}
                    st.session_state['units_obj'] = units
                    st.session_state['influent_obj'] = influent
                    st.session_state['was_flow_internal'] = was_flow_internal
                    st.session_state['reactor_volume_internal'] = reactor_volume_internal
                    st.session_state['influent_cod'] = total_cod
                    st.session_state['influent_tkn'] = total_tkn
                    
                    final_rates = plant.last_rates
                    oxygen_demand_g_d = final_rates['AOR_H_g_d'] + final_rates['AOR_A_g_d']
                    st.session_state['oxygen_demand_kg_d'] = oxygen_demand_g_d / 1000.0

                    st.session_state['design_criteria'] = {
                        'primary_sor': primary_sor, 'secondary_sor': secondary_sor,
                        'reactor_depth_user': reactor_depth, 'clarifier_depth_user': clarifier_depth,
                        'eq_tank_factor': eq_tank_factor, 'pump_eff': pump_eff,
                        'eq_depth_user': eq_depth, 'eq_pump_tdh_user': eq_pump_tdh
                    }
            st.rerun()

        if st.session_state.get('run_complete'):
            st.subheader("Simulation Summary")
            steady_state = st.session_state['results_data']['steady_state']
            st.metric("Effluent Ammonia (mg-N/L)", f"{steady_state['Ammonia (S_NH)']:.2f}")
            st.metric("Effluent Nitrate (mg-N/L)", f"{steady_state['Nitrate (S_NO3)']:.2f}")
            st.metric("Effluent Phosphorus (mg-P/L)", f"{steady_state['Soluble Phosphorus (S_PO4)']:.2f}")
            
            df_to_download = st.session_state['results_data']['df']
            csv = df_to_download.to_csv(index=False).encode('utf-8')
            st.download_button(label="Download Full Results (CSV)", data=csv, file_name='simulation_results.csv', mime='text/csv')
        else:
            st.info("Adjust parameters in the sidebar and click 'Run Simulation'.")

    with col2:
        st.subheader("Results Plot")
        if st.session_state.get('run_complete'):
            if os.path.exists("simulation_results.png"):
                st.image("simulation_results.png", use_container_width=True)
            else:
                st.error("Plot file could not be generated.")
        else:
            st.info("Results will be displayed here after running the simulation.")

with design_tab:
    st.header("Preliminary Unit Sizing")
    if st.session_state.get('run_complete'):
        units = st.session_state['units_obj']
        influent = st.session_state['influent_obj']
        criteria = st.session_state['design_criteria']
        oxygen_demand_kg_d = st.session_state['oxygen_demand_kg_d']
        
        # --- ADDITION: EQ System Sizing ---
        st.subheader("Equalization System")
        eq_vol_internal = influent.base_flow_rate * (criteria['eq_tank_factor'] / 100.0)
        eq_depth_internal = units.to_internal_head(criteria['eq_depth_user'])
        eq_area_internal = eq_vol_internal / eq_depth_internal
        eq_w_internal = math.sqrt(eq_area_internal)
        
        # Pump Power Calculation
        flow_m3s = influent.base_flow_rate / (24 * 3600) # m³/s
        head_m = units.to_internal_head(criteria['eq_pump_tdh_user'])
        pump_eff_frac = criteria['pump_eff'] / 100.0
        rho_g = 1000 * 9.81 # kg/m³ * m/s²
        
        power_watts = (rho_g * flow_m3s * head_m) / pump_eff_frac
        power_kw = power_watts / 1000.0
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric(label="EQ Tank Volume", value=f"{units.from_internal_volume(eq_vol_internal, 'tank'):.2f} {units.get_volume_label(short=True)}")
        with col2:
            st.metric(label="EQ Tank Area", value=f"{units.from_internal_area(eq_area_internal):.0f} {units.get_area_label()}")
        with col3:
            st.metric(label="EQ Pump Power (per pump)", value=f"{units.from_internal_power(power_kw):.1f} {units.get_power_label()}")


        st.subheader("Liquid Stream Units")
        pc_area_internal = influent.base_flow_rate / units.to_internal_sor(criteria['primary_sor'])
        sc_area_internal = influent.base_flow_rate / units.to_internal_sor(criteria['secondary_sor'])
        
        col1, col2 = st.columns(2)
        with col1:
            st.metric(label="Primary Clarifier Surface Area", value=f"{units.from_internal_area(pc_area_internal):.0f} {units.get_area_label()}")
        with col2:
            st.metric(label="Secondary Clarifier Surface Area", value=f"{units.from_internal_area(sc_area_internal):.0f} {units.get_area_label()}")

        st.subheader("Aeration System")
        alpha, beta, DO_setpoint, blower_eff = 0.65, 0.95, 2.0, 70
        C_sat_field = 9.0 
        
        AOR_kg_hr = oxygen_demand_kg_d / 24.0
        sotr_correction = alpha * (((beta * C_sat_field) - DO_setpoint) / 9.09)
        SOTR_kg_hr = AOR_kg_hr / sotr_correction
        airflow_m3_hr = SOTR_kg_hr / 0.15
        airflow_display = units.from_internal_airflow(airflow_m3_hr)
        blower_power_kw = SOTR_kg_hr * 0.5 / (blower_eff / 100.0)
        blower_power_display = units.from_internal_power(blower_power_kw)
        
        col1, col2, col3 = st.columns(3)
        with col1:
             st.metric("Actual Oxygen Rate (AOR)", f"{AOR_kg_hr:.1f} kg O₂/hr")
        with col2:
            st.metric("Required Airflow", f"{airflow_display:.0f} {units.get_airflow_label()}")
        with col3:
            st.metric("Estimated Blower Power", f"{blower_power_display:.1f} {units.get_power_label()}")

    else:
        st.info("Run a simulation on the 'Process Simulation' tab to see design calculations.")

def get_diagram_data(units, criteria):
    """Helper function to create the data dictionary for diagrams."""
    # EQ Tank Dims
    eq_vol_internal = st.session_state['influent_obj'].base_flow_rate * (criteria['eq_tank_factor'] / 100.0)
    eq_depth_internal = units.to_internal_head(criteria['eq_depth_user'])
    eq_area_internal = eq_vol_internal / eq_depth_internal
    eq_w_internal = math.sqrt(eq_area_internal)
    eq_dims_m = f"{eq_w_internal:.1f}m W x {eq_w_internal:.1f}m L x {eq_depth_internal:.1f}m D"
    eq_dims_ft = f"{units.from_internal_length(eq_w_internal):.1f}ft W x {units.from_internal_length(eq_w_internal):.1f}ft L x {criteria['eq_depth_user']:.1f}ft D"

    # Reactor Dims
    reactor_depth_internal = units.to_internal_head(criteria['reactor_depth_user'])
    reactor_area = st.session_state['reactor_volume_internal'] / reactor_depth_internal
    reactor_w = math.sqrt(reactor_area)
    reactor_dims_m = f"{reactor_w:.1f}m W x {reactor_w:.1f}m L x {reactor_depth_internal:.1f}m D"
    reactor_dims_ft = f"{units.from_internal_length(reactor_w):.1f}ft W x {units.from_internal_length(reactor_w):.1f}ft L x {criteria['reactor_depth_user']:.1f}ft D"

    sc_area_internal = st.session_state['influent_obj'].base_flow_rate / units.to_internal_sor(criteria['secondary_sor'])
    sc_diam_internal = math.sqrt(4 * sc_area_internal / math.pi)
    sc_dims_m = f"Ø{sc_diam_internal:.1f}m"
    sc_dims_ft = f"Ø{units.from_internal_length(sc_diam_internal):.1f}ft"
    
    q_in_internal = st.session_state['influent_obj'].base_flow_rate
    q_was_internal = st.session_state['was_flow_internal']
    q_eff_internal = q_in_internal - q_was_internal
    
    steady_state = st.session_state['results_data']['steady_state']
    mlss = steady_state['Heterotrophic Biomass (X_BH)'] + steady_state['Autotrophic Biomass (X_AUT)'] + steady_state['Particulate Inert (X_I)'] + steady_state['Particulate Substrate (X_S)']

    q_ras_internal = q_in_internal * 0.5

    data = {
        'q_in': units.from_internal_flow(q_in_internal),
        'q_was': units.from_internal_flow(q_was_internal),
        'q_ras': units.from_internal_flow(q_ras_internal),
        'q_eff': units.from_internal_flow(q_eff_internal),
        'q_ml': units.from_internal_flow(q_in_internal + q_ras_internal),
        'mlss': mlss,
        'inf_cod': st.session_state['influent_cod'], 
        'inf_tkn': st.session_state['influent_tkn'],
        'eff_nh3': steady_state['Ammonia (S_NH)'],
        'eff_no3': steady_state['Nitrate (S_NO3)'],
        'eq_dims': eq_dims_ft if units.system == 'US Customary' else eq_dims_m,
        'reactor_dims': reactor_dims_ft if units.system == 'US Customary' else reactor_dims_m,
        'sc_dims': sc_dims_ft if units.system == 'US Customary' else sc_dims_m,
        'pc_dims': "(Simplified)",
    }
    return data

with pfd_tab:
    st.header("Process Flow Diagram")
    if st.session_state.get('run_complete'):
        units = st.session_state['units_obj']
        criteria = st.session_state['design_criteria']
        pfd_data = get_diagram_data(units, criteria)
        
        try:
            pfd_graph = create_pfd(pfd_data, units)
            st.graphviz_chart(pfd_graph)
        except Exception as e:
            st.error(f"Could not generate PFD. This usually means the Graphviz system software is not installed.")
            st.info("For local deployment, install Graphviz (e.g., `brew install graphviz` on macOS, `sudo apt-get install graphviz` on Debian/Ubuntu). This is pre-installed on Streamlit Community Cloud.")
            st.error(f"Details: {e}")
    else:
        st.info("Run a simulation on the 'Process Simulation' tab to generate the PFD.")

with pid_tab:
    st.header("Piping & Instrumentation Diagram (P&ID)")
    if st.session_state.get('run_complete'):
        units = st.session_state['units_obj']
        criteria = st.session_state['design_criteria']
        pid_data = get_diagram_data(units, criteria)
        
        AOR_kg_hr = st.session_state['oxygen_demand_kg_d'] / 24.0
        sotr_correction = 0.65 * (((0.95 * 9.0) - 2.0) / 9.09)
        SOTR_kg_hr = AOR_kg_hr / sotr_correction
        airflow_m3_hr = SOTR_kg_hr / 0.15
        pid_data['airflow'] = units.from_internal_airflow(airflow_m3_hr)
        
        try:
            pid_graph = create_pid(pid_data, units)
            st.graphviz_chart(pid_graph)
        except Exception as e:
            st.error(f"Could not generate P&ID. Ensure Graphviz is installed on your system.")
            st.error(f"Details: {e}")
    else:
        st.info("Run a simulation on the 'Process Simulation' tab to generate the P&ID.")

with external_tab:
    st.header("External Simulator")
    st.markdown("This tab embeds an external Streamlit application for comparison and reference.")
    st.components.v1.iframe("https://asm-simulator.streamlit.app/", height=800, scrolling=True)