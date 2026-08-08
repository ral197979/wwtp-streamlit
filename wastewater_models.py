# Part 1 of 2
import streamlit as st
import pandas as pd
import numpy as np
import math
import graphviz
from fpdf import FPDF
from PIL import Image

st.set_page_config(layout="wide", page_title="WWTP Design Tool")

# --- Unit Conversion Factors ---
unit_systems = {
    "SI (m³/d, kg)": {"flow": 1, "mass": 1, "volume": 1},
    "US Customary (MGD, lb)": {"flow": 0.00026417, "mass": 2.20462, "volume": 0.26417},
    "Metric (L/d, g)": {"flow": 1000, "mass": 1000, "volume": 1000}
}

# --- Sidebar Inputs ---
st.sidebar.title("Design Parameters")

unit_choice = st.sidebar.selectbox("Unit System", list(unit_systems.keys()))
unit = unit_systems[unit_choice]

st.sidebar.markdown("### 👥 Population & Flow")
population = st.sidebar.number_input("Population Served", 100, 1_000_000, 10000)
gpcd = st.sidebar.slider("Flow per Capita (gpcd)", 50, 300, 150)

st.sidebar.markdown("### 🧪 Influent Parameters")
inf_bod = st.sidebar.number_input("BOD₅ (mg/L)", 100, 500, 250)
inf_tss = st.sidebar.number_input("TSS (mg/L)", 100, 500, 250)
inf_tkn = st.sidebar.number_input("TKN (mg/L)", 20, 80, 40)
inf_tp = st.sidebar.number_input("TP (mg/L)", 2, 15, 6)

st.sidebar.markdown("### 🧬 Effluent Targets")
eff_bod = st.sidebar.number_input("Effluent BOD₅", 5, 60, 30)
eff_tss = st.sidebar.number_input("Effluent TSS", 5, 60, 30)
eff_tn = st.sidebar.number_input("Effluent TN", 5, 20, 10)
eff_tp = st.sidebar.number_input("Effluent TP", 0.1, 4.0, 1.0)

nutrient_removal = st.sidebar.checkbox("Include Nutrient Removal (N & P)", value=True)
chem_choice = st.sidebar.selectbox("P Removal Chemical", ["Alum", "Ferric"])

srt = st.sidebar.slider("SRT (days)", 3, 30, 10)
mlss = st.sidebar.slider("MLSS (mg/L)", 1500, 5000, 3000)
mvss = st.sidebar.slider("MLVSS (mg/L)", 1000, 4000, 2500)
ras_ratio = st.sidebar.slider("RAS Ratio", 0.2, 1.0, 0.5)

# --- Design Tabs ---
tabs = st.tabs(["Activated Sludge", "IFAS", "MBBR", "GIS Layout"])

# === Activated Sludge Tab ===
with tabs[0]:
    st.header("Activated Sludge Process Design")

    Q_m3d = population * gpcd * 3.78541 / 1000
    Q_peak = Q_m3d * 2.5
    Q_peak_h = Q_peak / 24 * 4

    BOD_kg = Q_m3d * inf_bod / 1000
    TSS_kg = Q_m3d * inf_tss / 1000
    TKN_kg = Q_m3d * inf_tkn / 1000
    TP_kg = Q_m3d * inf_tp / 1000

    HRT = 6 if not nutrient_removal else 10
    V_reactor = Q_m3d * HRT / 24
    MLVSS_kg = V_reactor * mvss / 1000
    FM = BOD_kg / MLVSS_kg

    O2_BOD = 1.42 * BOD_kg
    O2_N = 4.57 * TKN_kg if nutrient_removal else 0
    O2_total = O2_BOD + O2_N
    blower_kWh_d = O2_total / 1.8
    blower_kW = blower_kWh_d / 24

    EQ_vol = Q_peak_h * 3
    pump_kW = Q_peak_h * 6 * 9.81 / (3600 * 0.65)

    clar_area = Q_peak / 65
    clar_dia = math.sqrt((clar_area * 4) / math.pi)

    WAS_kg = BOD_kg * 0.5
    sludge_m3d = WAS_kg / 50

    ras_flow = Q_m3d * ras_ratio
    was_flow = sludge_m3d

    P_removed = max(0, TP_kg - (Q_m3d * eff_tp / 1000))
    chem_factor = 5.5 if chem_choice == "Alum" else 4.5
    chem_kg = P_removed * chem_factor

    capex = 500 * V_reactor + 800 * clar_area + 300 * EQ_vol
    opex = blower_kWh_d * 0.12 * 365 + sludge_m3d * 100 * 365 + chem_kg * 0.5 * 365 + 50000

    sludge_age = MLVSS_kg / WAS_kg if WAS_kg > 0 else 0
    solids_loading = (Q_m3d * mlss) / clar_area

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Flow & Loading")
        st.write(f"**Avg Flow:** {Q_m3d:.1f} m³/d")
        st.write(f"**BOD Load:** {BOD_kg:.1f} kg/d")
        st.write(f"**Reactor Vol:** {V_reactor:.1f} m³")
        st.write(f"**F/M Ratio:** {FM:.3f}")
        st.write(f"**O₂ Demand:** {O2_total:.1f} kg/d")
        st.write(f"**Sludge Age (SRT):** {sludge_age:.1f} days")
    with col2:
        st.subheader("Equipment")
        st.write(f"**Clarifier Area:** {clar_area:.1f} m²")
        st.write(f"**Solids Loading:** {solids_loading:.1f} kg/m²·d")
        st.write(f"**Pump Power:** {pump_kW:.2f} kW")
        st.write(f"**Blower Power:** {blower_kW:.1f} kW")
        st.write(f"**Sludge Volume (WAS):** {sludge_m3d:.1f} m³/d")
        st.write(f"**WAS Flow Rate:** {was_flow:.2f} m³/d")
        st.write(f"**RAS Flow Rate:** {ras_flow:.2f} m³/d")
        if nutrient_removal:
            st.write(f"**P Removed:** {P_removed:.2f} kg/d")
            st.write(f"**{chem_choice} Dose:** {chem_kg:.2f} kg/d")
        st.write(f"**CAPEX:** ${capex:,.0f}")
        st.write(f"**OPEX:** ${opex:,.0f}/yr")

# (The rest remains unchanged for IFAS, MBBR, PDF export)
