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

st.sidebar.markdown("### 🧼 Effluent Targets")
eff_bod = st.sidebar.number_input("Effluent BOD₅", 5, 60, 30)
eff_tss = st.sidebar.number_input("Effluent TSS", 5, 60, 30)
eff_tn = st.sidebar.number_input("Effluent TN", 5, 20, 10)
eff_tp = st.sidebar.number_input("Effluent TP", 0.1, 4.0, 1.0)

nutrient_removal = st.sidebar.checkbox("Include Nutrient Removal (N & P)", value=True)
chem_choice = st.sidebar.selectbox("P Removal Chemical", ["Alum", "Ferric"])

srt = st.sidebar.slider("SRT (days)", 3, 30, 10)
mlss = st.sidebar.slider("MLSS (mg/L)", 1500, 5000, 3000)

# --- Design Tabs ---
tabs = st.tabs(["Activated Sludge", "IFAS", "MBBR", "GIS Layout"])

# === Activated Sludge Tab ===
with tabs[0]:
    st.header("Activated Sludge Process Design")

    # Flow
    Q_m3d = population * gpcd * 3.78541 / 1000
    Q_peak = Q_m3d * 2.5
    Q_peak_h = Q_peak / 24 * 4

    # Loads
    BOD_kg = Q_m3d * inf_bod / 1000
    TSS_kg = Q_m3d * inf_tss / 1000
    TKN_kg = Q_m3d * inf_tkn / 1000
    TP_kg = Q_m3d * inf_tp / 1000

    # Reactor
    HRT = 6 if not nutrient_removal else 10
    V_reactor = Q_m3d * HRT / 24
    MLVSS_kg = V_reactor * mlss / 1000
    FM = BOD_kg / MLVSS_kg

    # Oxygen
    O2_BOD = 1.42 * BOD_kg
    O2_N = 4.57 * TKN_kg if nutrient_removal else 0
    O2_total = O2_BOD + O2_N
    blower_kWh_d = O2_total / 1.8
    blower_kW = blower_kWh_d / 24

    # Equalization
    EQ_vol = Q_peak_h * 3
    pump_kW = Q_peak_h * 6 * 9.81 / (3600 * 0.65)

    # Clarifier
    clar_area = Q_peak / 65
    clar_dia = math.sqrt((clar_area * 4) / math.pi)

    # Sludge
    WAS_kg = BOD_kg * 0.5
    sludge_m3d = WAS_kg / 50

    # P Removal
    P_removed = max(0, TP_kg - (Q_m3d * eff_tp / 1000))
    chem_factor = 5.5 if chem_choice == "Alum" else 4.5
    chem_kg = P_removed * chem_factor

    # Cost
    capex = 500 * V_reactor + 800 * clar_area + 300 * EQ_vol
    opex = blower_kWh_d * 0.12 * 365 + sludge_m3d * 100 * 365 + chem_kg * 0.5 * 365 + 50000

    # Outputs
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Flow & Loading")
        st.write(f"**Avg Flow:** {Q_m3d:.1f} m³/d")
        st.write(f"**BOD Load:** {BOD_kg:.1f} kg/d")
        st.write(f"**Reactor Vol:** {V_reactor:.1f} m³")
        st.write(f"**F/M Ratio:** {FM:.3f}")
        st.write(f"**O₂ Demand:** {O2_total:.1f} kg/d")
    with col2:
        st.subheader("Equipment")
        st.write(f"**Clarifier Area:** {clar_area:.1f} m²")
        st.write(f"**Pump Power:** {pump_kW:.2f} kW")
        st.write(f"**Blower Power:** {blower_kW:.1f} kW")
        st.write(f"**Sludge Volume:** {sludge_m3d:.1f} m³/d")
        if nutrient_removal:
            st.write(f"**P Removed:** {P_removed:.2f} kg/d")
            st.write(f"**{chem_choice} Dose:** {chem_kg:.2f} kg/d")
        st.write(f"**CAPEX:** ${capex:,.0f}")
        st.write(f"**OPEX:** ${opex:,.0f}/yr")
# === IFAS Tab ===
with tabs[1]:
    st.header("IFAS Process Design")
    st.info("This tab will reuse influent & effluent values, and simulate a hybrid fixed + suspended growth system.")

    biofilm_surface_area = st.slider("Biofilm Surface Area (m²/m³)", 400, 1000, 600)
    fill_fraction = st.slider("Media Fill Fraction (%)", 20, 70, 40)

    media_vol = V_reactor * (fill_fraction / 100)
    total_surface_area = media_vol * biofilm_surface_area

    st.write(f"**Media Volume:** {media_vol:.1f} m³")
    st.write(f"**Biofilm Surface Area:** {total_surface_area:.1f} m²")

# === MBBR Tab ===
with tabs[2]:
    st.header("MBBR Process Design")
    st.info("This tab simulates pure attached-growth using MBBR media.")

    reactor_vol_mbbr = st.slider("Reactor Volume (m³)", 500, 5000, 2000)
    biofilm_surface_area_mbbr = st.slider("Biofilm Area (m²/m³)", 300, 1000, 500)
    fill_fraction_mbbr = st.slider("Fill Fraction (%)", 30, 70, 50)

    media_vol_mbbr = reactor_vol_mbbr * fill_fraction_mbbr / 100
    surface_area_mbbr = media_vol_mbbr * biofilm_surface_area_mbbr

    st.write(f"**Media Volume:** {media_vol_mbbr:.1f} m³")
    st.write(f"**Biofilm Surface Area:** {surface_area_mbbr:.1f} m²")

# === GIS Layout Tab ===
with tabs[3]:
    st.header("🗺️ GIS Plant Layout")
    st.markdown("Upload a PNG or JPG site layout of your WWTP. Then visualize the layout with unit labels.")

    uploaded_file = st.file_uploader("Upload Site Layout Image (PNG or JPG)", type=["png", "jpg", "jpeg"])
    if uploaded_file:
        image = Image.open(uploaded_file)
        st.image(image, caption="Plant Layout", use_column_width=True)
        st.markdown("**Suggested Unit Positions:**")
        st.markdown("- ✅ EQ Tank → Left\n- ✅ Aeration Basin → Center\n- ✅ Clarifier → Right\n- ✅ Chemical Dosing → Above")
    else:
        st.warning("Please upload a layout image to view the GIS overlay.")

# === PDF Export ===
def generate_pdf(summary: dict):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", "B", 16)
    pdf.cell(0, 10, "WWTP Design Report", ln=True)
    pdf.set_font("Arial", "", 12)
    for k, v in summary.items():
        pdf.cell(0, 8, f"{k}: {v}", ln=True)
    return pdf.output(dest="S").encode("latin-1")

if tabs[0]:  # Activated Sludge tab
    summary_dict = {
        "Flow (m³/d)": f"{Q_m3d:.1f}",
        "Reactor Volume": f"{V_reactor:.1f} m³",
        "Blower Power": f"{blower_kW:.1f} kW",
        "Sludge Volume": f"{sludge_m3d:.1f} m³/d",
        "Pump Power": f"{pump_kW:.2f} kW",
        "Clarifier Area": f"{clar_area:.1f} m²",
        "P Removed": f"{P_removed:.2f} kg/d" if nutrient_removal else "N/A",
        "Chemical Dose": f"{chem_kg:.2f} kg/d" if nutrient_removal else "N/A",
        "CAPEX": f"${capex:,.0f}",
        "OPEX": f"${opex:,.0f}/yr"
    }
    pdf_bytes = generate_pdf(summary_dict)
    st.download_button("📄 Download PDF Report", data=pdf_bytes, file_name="WWTP_Design_Report.pdf", mime="application/pdf")

# === Footer ===
st.markdown("---")
st.caption("Developed using Streamlit · Based on EPA & Metcalf & Eddy Design Standards · © 2025")

