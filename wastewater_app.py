import streamlit as st
import math

# --- Page Configuration ---
st.set_page_config(
    page_title="Wastewater Process Design Simulator",
    page_icon="💧",
    layout="wide"
)

# --- Calculation Logic ---
def calculate_design(inputs):
    """
    Performs all wastewater design calculations based on user inputs.
    
    Args:
        inputs (dict): A dictionary containing all the input values from the sidebar.
        
    Returns:
        dict: A dictionary containing all the calculated results.
    """
    # Unpack inputs for easier use
    Q = inputs['qAvg']
    S0 = inputs['s0']
    NH3N_in = inputs['nh3n']
    peakFactor = inputs['peakFactor']
    Y = inputs['y']
    kd = inputs['kd']
    SRT = inputs['srt']
    MLSS = inputs['mlss']
    mlvssRatio = inputs['mlvssRatio'] / 100
    Se = inputs['seBOD']
    maxSOR = inputs['maxSOR']
    fillFraction = inputs['fillFraction'] / 100
    carrierSSA = inputs['carrierSSA']
    salr = inputs['salr']

    results = {}

    # 1. Aeration Tank & Sludge
    results['bodLoad'] = Q * S0 / 1000 if Q > 0 and S0 > 0 else 0
    results['MLVSS'] = MLSS * mlvssRatio
    
    # Check for division by zero
    denominator_tank_vol = results['MLVSS'] * (1 + kd * SRT)
    if denominator_tank_vol > 0:
        results['tankVolume'] = (Q * Y * (S0 - Se) * SRT) / denominator_tank_vol
    else:
        results['tankVolume'] = 0

    results['hrt'] = (results['tankVolume'] / Q) * 24 if Q > 0 else 0
    
    biomassInTank = results['tankVolume'] * results['MLVSS'] / 1000
    results['was'] = biomassInTank / SRT if SRT > 0 else 0
    
    denominator_fm = results['MLVSS'] * results['tankVolume'] / 1000
    results['fmRatio'] = results['bodLoad'] / denominator_fm if denominator_fm > 0 else 0

    # 2. Secondary Clarifier
    results['qPeak'] = Q * peakFactor
    results['clarifierArea'] = results['qPeak'] / maxSOR if maxSOR > 0 else 0
    results['peakSOR'] = results['qPeak'] / results['clarifierArea'] if results['clarifierArea'] > 0 else 0
    results['clarifierDiameter'] = math.sqrt(4 * results['clarifierArea'] / math.pi) if results['clarifierArea'] > 0 else 0

    # 3. MBBR Design
    results['nh3nLoad'] = Q * NH3N_in / 1000
    results['requiredCarrierArea'] = (results['nh3nLoad'] * 1000) / salr if salr > 0 else 0
    results['carrierVolume'] = results['requiredCarrierArea'] / carrierSSA if carrierSSA > 0 else 0
    results['mbbrTankVolume'] = results['carrierVolume'] / fillFraction if fillFraction > 0 else 0
    
    return results

# --- UI Layout ---

# Header
st.title("💧 Wastewater Process Design Simulator")
st.markdown("This application provides a preliminary design for an activated sludge process with MBBR integration for nitrification. Adjust the parameters in the sidebar to see the design calculations update in real-time.")

# Sidebar for Inputs
st.sidebar.header("Inputs & Assumptions")

inputs = {
    'qAvg': st.sidebar.number_input("Average Daily Flow (Q, m³/day)", min_value=0.0, value=5000.0, step=100.0),
    's0': st.sidebar.number_input("Influent BOD₅ (S₀, mg/L)", min_value=0.0, value=250.0, step=10.0),
    'nh3n': st.sidebar.number_input("Influent NH₃-N (mg/L)", min_value=0.0, value=40.0, step=5.0),
    'peakFactor': st.sidebar.number_input("Peaking Factor (for Qpeak)", min_value=1.0, value=2.5, step=0.1),
    'y': st.sidebar.number_input("Biomass Yield (Y, mg VSS/mg BOD₅)", min_value=0.0, value=0.5, step=0.05, format="%.2f"),
    'kd': st.sidebar.number_input("Endogenous Decay Coeff. (kₔ, 1/day)", min_value=0.0, value=0.06, step=0.01, format="%.3f"),
    'srt': st.sidebar.number_input("Target Solids Retention Time (SRT, days)", min_value=1.0, value=12.0, step=1.0),
    'mlss': st.sidebar.number_input("Target MLSS (mg/L)", min_value=0.0, value=3500.0, step=100.0),
    'mlvssRatio': st.sidebar.number_input("MLVSS/MLSS Ratio (%)", min_value=0.0, value=80.0, step=1.0),
    'seBOD': st.sidebar.number_input("Effluent BOD₅ (Se, mg/L)", min_value=0.0, value=10.0, step=1.0),
    'maxSOR': st.sidebar.number_input("Max Secondary Clarifier SOR (m³/m²·d)", min_value=1.0, value=48.0, step=1.0),
    'fillFraction': st.sidebar.number_input("MBBR Carrier Fill Fraction (%)", min_value=0.0, max_value=100.0, value=50.0, step=1.0),
    'carrierSSA': st.sidebar.number_input("Carrier Specific Surface Area (m²/m³)", min_value=0.0, value=500.0, step=10.0),
    'salr': st.sidebar.number_input("Surface Area Loading Rate (SALR, g NH₃-N/m²·d)", min_value=0.0, value=1.2, step=0.1)
}

# --- Main Panel for Outputs ---
st.header("Calculated Design Parameters")

# Perform calculations
results = calculate_design(inputs)

# Create columns for the result cards
col1, col2, col3 = st.columns(3)

with col1:
    st.subheader("📊 Aeration Tank & Sludge")
    st.metric(label="BOD₅ Load (kg/day)", value=f"{results['bodLoad']:.2f}")
    st.metric(label="MLVSS Target (mg/L)", value=f"{results['MLVSS']:.0f}")
    st.metric(label="Required Tank Volume (m³)", value=f"{results['tankVolume']:.0f}")
    st.metric(label="Hydraulic Retention Time (HRT, hours)", value=f"{results['hrt']:.2f}")
    st.metric(label="Sludge Wasted (WAS, kg/day)", value=f"{results['was']:.2f}")
    st.metric(label="Food to Microorganism (F:M) Ratio", value=f"{results['fmRatio']:.2f}")

with col2:
    st.subheader("💧 Secondary Clarifier")
    st.metric(label="Peak Flow (Qpeak, m³/day)", value=f"{results['qPeak']:.0f}")
    st.metric(label="Required Surface Area (m²)", value=f"{results['clarifierArea']:.1f}")
    
    # Conditional Formatting for SOR
    peak_sor_val = results['peakSOR']
    max_sor_val = inputs['maxSOR']
    st.metric(label="Peak Surface Overflow Rate (SOR, m³/m²·d)", value=f"{peak_sor_val:.2f}")
    if peak_sor_val > max_sor_val:
        st.error(f"Warning: Peak SOR ({peak_sor_val:.2f}) exceeds maximum ({max_sor_val:.2f}). Consider increasing clarifier area.")
    elif peak_sor_val > max_sor_val * 0.9:
        st.warning(f"Note: Peak SOR ({peak_sor_val:.2f}) is close to the maximum limit ({max_sor_val:.2f}).")
    else:
        st.success(f"Peak SOR ({peak_sor_val:.2f}) is within the design limit ({max_sor_val:.2f}).")

    st.metric(label="Calculated Diameter (m)", value=f"{results['clarifierDiameter']:.2f}")

with col3:
    st.subheader("🦠 MBBR for Nitrification")
    st.metric(label="Ammonia Load (kg/day)", value=f"{results['nh3nLoad']:.2f}")
    st.metric(label="Required Carrier Surface Area (m²)", value=f"{results['requiredCarrierArea']:.0f}")
    st.metric(label="Required Carrier Volume (m³)", value=f"{results['carrierVolume']:.1f}")
    st.metric(label="Required MBBR Tank Volume (m³)", value=f"{results['mbbrTankVolume']:.1f}")

# --- Instructions and Notes ---
with st.expander("About this Calculator and Formulas Used"):
    st.markdown("""
        ### How to Use
        1.  **Adjust Inputs:** Use the controls in the sidebar on the left to enter your specific project data and design assumptions.
        2.  **Review Outputs:** The calculated design parameters will update automatically in the main panel.
        3.  **Check Warnings:** Pay attention to any success, warning, or error messages, as they indicate if key design parameters are within typical limits.

        ### Key Formulas
        - **Aeration Tank Volume (V):** `V = (Q * Y * (S₀ - Se) * SRT) / (MLVSS * (1 + kₔ * SRT))`
        - **Sludge Wasted (Px):** `Px = (V * MLVSS) / SRT`
        - **F:M Ratio:** `F:M = (Q * S₀) / (V * MLVSS)`
        - **Clarifier Area (A):** `A = Q_peak / SOR_max`
        - **MBBR Required Surface Area (A_carrier):** `A_carrier = (NH₃-N Load) / SALR`
        ---
        *Disclaimer: This is a tool for preliminary design and educational purposes. Final engineering designs should be conducted and verified by a qualified professional.*
    """)

