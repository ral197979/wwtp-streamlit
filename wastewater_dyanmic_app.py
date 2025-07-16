import streamlit as st
import math

# --- Page Configuration ---
st.set_page_config(
    page_title="Wastewater Process Design Simulator",
    page_icon="💧",
    layout="wide"
)

# --- Conversion Factors ---
MGD_to_M3D = 3785.41
M3_to_GAL = 264.172
M2_to_FT2 = 10.7639
M_to_FT = 3.28084
KG_to_LBS = 2.20462
GFD_to_MMD = 1 / 24.5425 # Gal/ft²/day to m³/m²/day

# --- Calculation Logic ---
def calculate_design(inputs, unit_system):
    """
    Performs all wastewater design calculations based on user inputs.
    All internal calculations are done in Metric units.
    
    Args:
        inputs (dict): A dictionary containing all the input values from the sidebar.
        unit_system (str): The selected unit system ("Metric" or "US Conventional").
        
    Returns:
        dict: A dictionary containing all the calculated results in Metric units.
    """
    
    # --- Convert inputs to Metric if US Conventional is selected ---
    if unit_system == "US Conventional (MGD, gal, lbs)":
        Q = inputs['qAvg'] * MGD_to_M3D
        maxSOR = inputs['maxSOR'] * GFD_to_MMD
    else: # Metric
        Q = inputs['qAvg']
        maxSOR = inputs['maxSOR']

    # Unpack other inputs (they are unit-independent or already handled)
    S0 = inputs['s0']
    NH3N_in = inputs['nh3n']
    peakFactor = inputs['peakFactor']
    Y = inputs['y']
    kd = inputs['kd']
    SRT = inputs['srt']
    MLSS = inputs['mlss']
    mlvssRatio = inputs['mlvssRatio'] / 100
    Se = inputs['seBOD']
    fillFraction = inputs['fillFraction'] / 100
    carrierSSA = inputs['carrierSSA']
    salr = inputs['salr']

    results = {}

    # 1. Aeration Tank & Sludge (Calculations in Metric)
    results['bodLoad_kg'] = Q * S0 / 1000 if Q > 0 and S0 > 0 else 0
    results['MLVSS'] = MLSS * mlvssRatio
    
    denominator_tank_vol = results['MLVSS'] * (1 + kd * SRT)
    if denominator_tank_vol > 0:
        results['tankVolume_m3'] = (Q * Y * (S0 - Se) * SRT) / denominator_tank_vol
    else:
        results['tankVolume_m3'] = 0

    results['hrt_hr'] = (results['tankVolume_m3'] / Q) * 24 if Q > 0 else 0
    
    biomassInTank_kg = results['tankVolume_m3'] * results['MLVSS'] / 1000
    results['was_kg'] = biomassInTank_kg / SRT if SRT > 0 else 0
    
    denominator_fm = results['MLVSS'] * results['tankVolume_m3'] / 1000
    results['fmRatio'] = results['bodLoad_kg'] / denominator_fm if denominator_fm > 0 else 0

    # 2. Secondary Clarifier (Calculations in Metric)
    results['qPeak_m3d'] = Q * peakFactor
    results['clarifierArea_m2'] = results['qPeak_m3d'] / maxSOR if maxSOR > 0 else 0
    results['peakSOR_mmd'] = results['qPeak_m3d'] / results['clarifierArea_m2'] if results['clarifierArea_m2'] > 0 else 0
    results['clarifierDiameter_m'] = math.sqrt(4 * results['clarifierArea_m2'] / math.pi) if results['clarifierArea_m2'] > 0 else 0

    # 3. MBBR Design (Calculations in Metric)
    results['nh3nLoad_kg'] = Q * NH3N_in / 1000
    results['requiredCarrierArea_m2'] = (results['nh3nLoad_kg'] * 1000) / salr if salr > 0 else 0
    results['carrierVolume_m3'] = results['requiredCarrierArea_m2'] / carrierSSA if carrierSSA > 0 else 0
    results['mbbrTankVolume_m3'] = results['carrierVolume_m3'] / fillFraction if fillFraction > 0 else 0
    
    return results

# --- UI Layout ---

# Header
st.title("💧 Wastewater Process Design Simulator")
st.markdown("This application provides a preliminary design for an activated sludge process with MBBR integration for nitrification. Adjust the parameters in the sidebar to see the design calculations update in real-time.")

# Sidebar for Inputs
st.sidebar.header("Inputs & Assumptions")

unit_system = st.sidebar.radio(
    "Select Unit System",
    ("Metric (m³, L, kg)", "US Conventional (MGD, gal, lbs)")
)

inputs = {}
if unit_system == "Metric (m³, L, kg)":
    inputs['qAvg'] = st.sidebar.number_input("Average Daily Flow (Q, m³/day)", min_value=0.0, value=5000.0, step=100.0)
    inputs['maxSOR'] = st.sidebar.number_input("Max Clarifier SOR (m³/m²·d)", min_value=1.0, value=48.0, step=1.0)
else: # US Conventional
    inputs['qAvg'] = st.sidebar.number_input("Average Daily Flow (Q, MGD)", min_value=0.0, value=1.3, step=0.1)
    inputs['maxSOR'] = st.sidebar.number_input("Max Clarifier SOR (gal/ft²·d)", min_value=1.0, value=1200.0, step=10.0)

# Unit-independent inputs
inputs['s0'] = st.sidebar.number_input("Influent BOD₅ (S₀, mg/L)", min_value=0.0, value=250.0, step=10.0)
inputs['nh3n'] = st.sidebar.number_input("Influent NH₃-N (mg/L)", min_value=0.0, value=40.0, step=5.0)
inputs['peakFactor'] = st.sidebar.number_input("Peaking Factor (for Qpeak)", min_value=1.0, value=2.5, step=0.1)
inputs['y'] = st.sidebar.number_input("Biomass Yield (Y, mg VSS/mg BOD₅)", min_value=0.0, value=0.5, step=0.05, format="%.2f")
inputs['kd'] = st.sidebar.number_input("Endogenous Decay Coeff. (kₔ, 1/day)", min_value=0.0, value=0.06, step=0.01, format="%.3f")
inputs['srt'] = st.sidebar.number_input("Target SRT (days)", min_value=1.0, value=12.0, step=1.0)
inputs['mlss'] = st.sidebar.number_input("Target MLSS (mg/L)", min_value=0.0, value=3500.0, step=100.0)
inputs['mlvssRatio'] = st.sidebar.number_input("MLVSS/MLSS Ratio (%)", min_value=0.0, value=80.0, step=1.0)
inputs['seBOD'] = st.sidebar.number_input("Effluent BOD₅ (Se, mg/L)", min_value=0.0, value=10.0, step=1.0)
inputs['fillFraction'] = st.sidebar.number_input("MBBR Carrier Fill Fraction (%)", min_value=0.0, max_value=100.0, value=50.0, step=1.0)
inputs['carrierSSA'] = st.sidebar.number_input("Carrier Specific Surface Area (m²/m³)", min_value=0.0, value=500.0, step=10.0)
inputs['salr'] = st.sidebar.number_input("Surface Area Loading Rate (SALR, g NH₃-N/m²·d)", min_value=0.0, value=1.2, step=0.1)

# --- Main Panel for Outputs ---
st.header("Calculated Design Parameters")

# Perform calculations
results = calculate_design(inputs, unit_system)

# Create columns for the result cards
col1, col2, col3 = st.columns(3)

# Display results based on selected unit system
if unit_system == "Metric (m³, L, kg)":
    with col1:
        st.subheader("📊 Aeration Tank & Sludge")
        st.metric(label="BOD₅ Load (kg/day)", value=f"{results['bodLoad_kg']:.2f}")
        st.metric(label="MLVSS Target (mg/L)", value=f"{results['MLVSS']:.0f}")
        st.metric(label="Required Tank Volume (m³)", value=f"{results['tankVolume_m3']:.0f}")
        st.metric(label="HRT (hours)", value=f"{results['hrt_hr']:.2f}")
        st.metric(label="Sludge Wasted (WAS, kg/day)", value=f"{results['was_kg']:.2f}")
        st.metric(label="F:M Ratio", value=f"{results['fmRatio']:.2f}")

    with col2:
        st.subheader("💧 Secondary Clarifier")
        st.metric(label="Peak Flow (m³/day)", value=f"{results['qPeak_m3d']:.0f}")
        st.metric(label="Required Surface Area (m²)", value=f"{results['clarifierArea_m2']:.1f}")
        
        peak_sor_val = results['peakSOR_mmd']
        max_sor_val = inputs['maxSOR']
        st.metric(label="Peak SOR (m³/m²·d)", value=f"{peak_sor_val:.2f}")
        if peak_sor_val > max_sor_val:
            st.error(f"Warning: Peak SOR ({peak_sor_val:.2f}) exceeds max ({max_sor_val:.2f}).")
        else:
            st.success(f"Peak SOR is within the design limit.")

        st.metric(label="Calculated Diameter (m)", value=f"{results['clarifierDiameter_m']:.2f}")

    with col3:
        st.subheader("🦠 MBBR for Nitrification")
        st.metric(label="Ammonia Load (kg/day)", value=f"{results['nh3nLoad_kg']:.2f}")
        st.metric(label="Required Carrier Surface Area (m²)", value=f"{results['requiredCarrierArea_m2']:.0f}")
        st.metric(label="Required Carrier Volume (m³)", value=f"{results['carrierVolume_m3']:.1f}")
        st.metric(label="Required MBBR Tank Volume (m³)", value=f"{results['mbbrTankVolume_m3']:.1f}")

else: # US Conventional
    with col1:
        st.subheader("📊 Aeration Tank & Sludge")
        st.metric(label="BOD₅ Load (lbs/day)", value=f"{results['bodLoad_kg'] * KG_to_LBS:,.2f}")
        st.metric(label="MLVSS Target (mg/L)", value=f"{results['MLVSS']:.0f}")
        st.metric(label="Required Tank Volume (gallons)", value=f"{results['tankVolume_m3'] * M3_to_GAL:,.0f}")
        st.metric(label="HRT (hours)", value=f"{results['hrt_hr']:.2f}")
        st.metric(label="Sludge Wasted (WAS, lbs/day)", value=f"{results['was_kg'] * KG_to_LBS:,.2f}")
        st.metric(label="F:M Ratio", value=f"{results['fmRatio']:.2f}")

    with col2:
        st.subheader("💧 Secondary Clarifier")
        st.metric(label="Peak Flow (MGD)", value=f"{results['qPeak_m3d'] / MGD_to_M3D:.2f}")
        st.metric(label="Required Surface Area (ft²)", value=f"{results['clarifierArea_m2'] * M2_to_FT2:,.1f}")
        
        peak_sor_val = results['peakSOR_mmd'] / GFD_to_MMD
        max_sor_val = inputs['maxSOR']
        st.metric(label="Peak SOR (gal/ft²·d)", value=f"{peak_sor_val:,.2f}")
        if peak_sor_val > max_sor_val:
            st.error(f"Warning: Peak SOR ({peak_sor_val:,.2f}) exceeds max ({max_sor_val:,.2f}).")
        else:
            st.success(f"Peak SOR is within the design limit.")

        st.metric(label="Calculated Diameter (ft)", value=f"{results['clarifierDiameter_m'] * M_to_FT:.2f}")

    with col3:
        st.subheader("🦠 MBBR for Nitrification")
        st.metric(label="Ammonia Load (lbs/day)", value=f"{results['nh3nLoad_kg'] * KG_to_LBS:,.2f}")
        st.metric(label="Required Carrier Surface Area (ft²)", value=f"{results['requiredCarrierArea_m2'] * M2_to_FT2:,.0f}")
        st.metric(label="Required Carrier Volume (ft³)", value=f"{results['carrierVolume_m3'] * M3_to_GAL / 7.481:,.1f}") # 7.481 gal/ft³
        st.metric(label="Required MBBR Tank Volume (ft³)", value=f"{results['mbbrTankVolume_m3'] * M3_to_GAL / 7.481:,.1f}")


# --- Instructions and Notes ---
with st.expander("About this Calculator and Formulas Used"):
    st.markdown("""
        ### How to Use
        1.  **Select Unit System:** Choose between Metric and US Conventional units in the sidebar.
        2.  **Adjust Inputs:** Use the controls in the sidebar to enter your specific project data and design assumptions.
        3.  **Review Outputs:** The calculated design parameters will update automatically in the main panel in the chosen units.
        4.  **Check Warnings:** Pay attention to any success, warning, or error messages, as they indicate if key design parameters are within typical limits.

        ### Key Formulas (Calculated in Metric)
        - **Aeration Tank Volume (V):** `V = (Q * Y * (S₀ - Se) * SRT) / (MLVSS * (1 + kₔ * SRT))`
        - **Sludge Wasted (Px):** `Px = (V * MLVSS) / SRT`
        - **F:M Ratio:** `F:M = (Q * S₀) / (V * MLVSS)`
        - **Clarifier Area (A):** `A = Q_peak / SOR_max`
        - **MBBR Required Surface Area (A_carrier):** `A_carrier = (NH₃-N Load * 1000) / SALR`
        ---
        *Disclaimer: This is a tool for preliminary design and educational purposes. Final engineering designs should be conducted and verified by a qualified professional.*
    """)
