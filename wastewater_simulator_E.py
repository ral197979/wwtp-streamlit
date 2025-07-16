import streamlit as st
import math
import csv
import io

# --- Unit Conversion Constants ---
# These constants are used to convert between Metric, US, and Imperial units.
M3_TO_MGD = 0.264172  # Cubic meters to Million US Gallons
M3_TO_IGD = 0.219969  # Cubic meters to Million Imperial Gallons
KG_TO_LBS = 2.20462   # Kilograms to Pounds
M2_TO_FT2 = 10.7639   # Square meters to Square feet
M_TO_FT = 3.28084     # Meters to Feet
L_TO_GAL_US = 0.264172 # Liters to US Gallons
L_TO_GAL_IMP = 0.219969# Liters to Imperial Gallons
M3_M2_D_TO_GPD_FT2 = 24.5424 # m³/m²/day to US gpd/ft²
M3_M2_D_TO_IGPD_FT2 = 20.4355 # m³/m²/day to Imperial gpd/ft²
KG_M3_D_TO_LBS_1000FT3_D = 62.428 # kg/m³/day to lbs/1000ft³/day

class WWTP_Design_Calculator:
    """
    A class to perform preliminary design calculations for a Municipal 
    Wastewater Treatment Plant (WWTP).

    This enhanced version includes multiple treatment technologies (CAS, MBR, MBBR, IFAS)
    and separates design criteria from the calculation logic for easier modification.
    """
    def __init__(self):
        """Initializes the calculator with default parameters and design criteria."""
        # --- User-configurable input parameters ---
        self.params = {
            'current_population': 50000,
            'growth_rate': 2.0,  # in percent
            'design_life': 30,   # in years
            'water_consumption': 350,  # Base unit: L/c/d
            'wastewater_fraction': 0.85, # Fraction of water consumption that becomes wastewater
            'bod_per_capita': 50, # g/c/d
            'tss_per_capita': 70, # g/c/d
            'tkn_per_capita': 12, # g/c/d
            'use_calculated_pf': True, # Use Harmon's peak factor formula
            'peak_factor_manual': 2.5,
        }

        # --- Technology-specific design criteria ---
        self.design_criteria = {
            'primary': {
                'bod_removal': 0.35, # Fractional BOD removal
                'hlr_design': 40,  # Hydraulic Loading Rate in m³/m²/d
                'hrt_design': 2.0, # Hydraulic Retention Time in hours
            },
            'CAS': { # Conventional Activated Sludge
                'theta_c': 10, 'mlss_conc': 3000, 'y': 0.5, 'kd': 0.06,
                's_effluent': 10, 'clarifier_hlr': 24,
            },
            'MBR': { # Membrane Bioreactor
                'theta_c': 15, 'mlss_conc': 10000, 'y': 0.5, 'kd': 0.06,
                's_effluent': 5, 'membrane_flux': 25,
            },
            'MBBR': { # Moving Bed Biofilm Reactor
                'bod_salr': 5.0, 'media_ssa': 500, 'media_fill': 0.55, 'clarifier_hlr': 32,
            },
            'IFAS': { # Integrated Fixed-Film Activated Sludge
                'theta_c': 8, 'mlss_conc': 2000, 'y': 0.5, 'kd': 0.06,
                's_effluent': 10, 'load_split_fixed_film': 0.6, 'bod_salr': 6.0,
                'media_ssa': 500, 'media_fill': 0.30, 'clarifier_hlr': 28,
            }
        }
        
        self.units = 'Metric'
        self.treatment_type = 'CAS'
        self.calculations = {}
        self.sizing = {}
        self.unit_config = {}
        self.set_units(self.units) # Initialize with default metric units

    def set_units(self, unit_system='Metric'):
        """Sets the unit system and updates the configuration for labels."""
        self.units = unit_system
        if unit_system == 'US':
            self.unit_config = {
                'flow': 'MGD', 'volume_large': 'MG', 'volume_small': 'ft³', 'area': 'ft²', 
                'depth': 'ft', 'load': 'lbs/day', 'consumption': 'gal/c/d',
                'hlr': 'gpd/ft²', 'olr': 'lbs/1000ft³/day', 'flux': 'gfd'
            }
        elif unit_system == 'Imperial':
            self.unit_config = {
                'flow': 'MIGD', 'volume_large': 'MIG', 'volume_small': 'ft³', 'area': 'ft²',
                'depth': 'ft', 'load': 'lbs/day', 'consumption': 'gal/c/d',
                'hlr': 'gpd/ft²', 'olr': 'lbs/1000ft³/day', 'flux': 'gpd'
            }
        else: # Metric
            self.unit_config = {
                'flow': 'm³/day', 'volume_large': 'm³', 'volume_small': 'm³', 'area': 'm²',
                'depth': 'm', 'load': 'kg/day', 'consumption': 'L/c/d',
                'hlr': 'm³/m²/day', 'olr': 'kg/m³/day', 'flux': 'L/m²/hr'
            }

    def run_calculations(self):
        """Executes all design calculations based on current parameters."""
        self._calculate_design_basis()
        self._calculate_sizing()

    def _calculate_design_basis(self):
        """Calculates future population, flows, and pollutant loads."""
        p0 = self.params['current_population']
        r = self.params['growth_rate'] / 100
        n = self.params['design_life']
        wc_L = self.params['water_consumption']
        projected_population = round(p0 * (1 + r)**n)
        # Corrected: Convert L/day to m³/day by dividing by 1000
        avg_daily_flow_m3 = (projected_population * wc_L * self.params['wastewater_fraction']) / 1000

        if self.params['use_calculated_pf']:
            p_thousands = projected_population / 1000
            peak_factor = (18 + math.sqrt(p_thousands)) / (4 + math.sqrt(p_thousands)) if p_thousands > 0 else 2.5
        else:
            peak_factor = self.params['peak_factor_manual']
        
        peak_hourly_flow_m3_d = avg_daily_flow_m3 * peak_factor

        self.calculations = {
            'projected_population': projected_population,
            'avg_daily_flow_m3': avg_daily_flow_m3,
            'peak_hourly_flow_m3_d': peak_hourly_flow_m3_d,
            'peak_factor': peak_factor,
            'bod_load_kg': (projected_population * self.params['bod_per_capita']) / 1000,
            'tss_load_kg': (projected_population * self.params['tss_per_capita']) / 1000,
            'tkn_load_kg': (projected_population * self.params['tkn_per_capita']) / 1000,
        }
        self.calculations['influent_bod'] = (self.calculations['bod_load_kg'] * 1000) / avg_daily_flow_m3 if avg_daily_flow_m3 > 0 else 0
        self.calculations['influent_tss'] = (self.calculations['tss_load_kg'] * 1000) / avg_daily_flow_m3 if avg_daily_flow_m3 > 0 else 0
        self.calculations['influent_tkn'] = (self.calculations['tkn_load_kg'] * 1000) / avg_daily_flow_m3 if avg_daily_flow_m3 > 0 else 0

    def _calculate_sizing(self):
        """Calculates the sizing for each treatment unit by dispatching to the correct method."""
        self._size_primary_treatment()
        bod_removal_primary = self.design_criteria['primary']['bod_removal']
        self.sizing['secondary_bod_load_kg'] = self.calculations['bod_load_kg'] * (1 - bod_removal_primary)
        sizing_method = getattr(self, f'_size_secondary_{self.treatment_type.lower()}', None)
        if sizing_method:
            sizing_method()

    def _size_primary_treatment(self):
        crit = self.design_criteria['primary']
        q = self.calculations['avg_daily_flow_m3']
        clarifier_area = q / crit['hlr_design'] if crit['hlr_design'] > 0 else 0
        clarifier_volume = (q / 24) * crit['hrt_design']
        self.sizing['primary'] = {
            'clarifier_area_m2': clarifier_area,
            'clarifier_volume_m3': clarifier_volume,
            'clarifier_depth_m': clarifier_volume / clarifier_area if clarifier_area > 0 else 0,
            'hydraulic_loading_rate_m3_m2_d': crit['hlr_design'],
        }

    def _size_secondary_cas(self):
        crit = self.design_criteria['CAS']
        q = self.calculations['avg_daily_flow_m3']
        s0 = self.calculations['influent_bod'] * (1 - self.design_criteria['primary']['bod_removal'])
        
        numerator = q * crit['y'] * (s0 - crit['s_effluent']) * crit['theta_c']
        denominator = crit['mlss_conc'] * (1 + crit['kd'] * crit['theta_c'])
        aeration_volume = numerator / denominator if denominator > 0 else 0
        
        px_vss = (crit['y'] * q * (s0 - crit['s_effluent'])) / (1 + crit['kd'] * crit['theta_c'])
        
        # Corrected: Calculate O2 demand based on BOD removed, not influent load
        bod_removed_kg = q * (s0 - crit['s_effluent']) / 1000 # Convert g/day to kg/day
        o2_demand = bod_removed_kg - (1.42 * px_vss) + (self.calculations['tkn_load_kg'] * 4.57)

        self.sizing['secondary'] = {
            'type': 'Conventional Activated Sludge (CAS)', 'aeration_volume_m3': aeration_volume,
            'aeration_hrt_hr': (aeration_volume / q) * 24 if q > 0 else 0,
            'fm_ratio': self.sizing['secondary_bod_load_kg'] / (aeration_volume * (crit['mlss_conc'] / 1000)) if aeration_volume > 0 else 0,
            'sludge_production_kg_d': px_vss, 'oxygen_demand_kg_d': o2_demand,
            'secondary_clarifier_area_m2': q / crit['clarifier_hlr'],
            'organic_loading_rate_kg_m3_d': self.sizing['secondary_bod_load_kg'] / aeration_volume if aeration_volume > 0 else 0,
            'secondary_clarifier_hlr_m3_m2_d': crit['clarifier_hlr'],
        }

    def _size_secondary_mbr(self):
        crit = self.design_criteria['MBR']
        q = self.calculations['avg_daily_flow_m3']
        s0 = self.calculations['influent_bod'] * (1 - self.design_criteria['primary']['bod_removal'])
        numerator = q * crit['y'] * (s0 - crit['s_effluent']) * crit['theta_c']
        denominator = crit['mlss_conc'] * (1 + crit['kd'] * crit['theta_c'])
        aeration_volume = numerator / denominator if denominator > 0 else 0
        flow_lph = (q * 1000) / 24
        membrane_area = flow_lph / crit['membrane_flux'] if crit['membrane_flux'] > 0 else 0
        self.sizing['secondary'] = {
            'type': 'Membrane Bioreactor (MBR)', 'aeration_volume_m3': aeration_volume,
            'aeration_hrt_hr': (aeration_volume / q) * 24 if q > 0 else 0,
            'membrane_area_m2': membrane_area, 'membrane_flux_lmh': crit['membrane_flux'],
            'organic_loading_rate_kg_m3_d': self.sizing['secondary_bod_load_kg'] / aeration_volume if aeration_volume > 0 else 0,
        }

    def _size_secondary_mbbr(self):
        crit = self.design_criteria['MBBR']
        bod_load_g_d = self.sizing['secondary_bod_load_kg'] * 1000
        media_area_req = bod_load_g_d / crit['bod_salr'] if crit['bod_salr'] > 0 else 0
        reactor_volume = media_area_req / (crit['media_ssa'] * crit['media_fill']) if (crit['media_ssa'] * crit['media_fill']) > 0 else 0
        q = self.calculations['avg_daily_flow_m3']
        self.sizing['secondary'] = {
            'type': 'Moving Bed Biofilm Reactor (MBBR)', 'reactor_volume_m3': reactor_volume,
            'reactor_hrt_hr': (reactor_volume / q) * 24 if q > 0 else 0,
            'media_surface_area_req_m2': media_area_req,
            'organic_loading_rate_kg_m3_d': self.sizing['secondary_bod_load_kg'] / reactor_volume if reactor_volume > 0 else 0,
            'secondary_clarifier_area_m2': q / crit['clarifier_hlr'],
            'secondary_clarifier_hlr_m3_m2_d': crit['clarifier_hlr'],
        }

    def _size_secondary_ifas(self):
        crit = self.design_criteria['IFAS']
        q = self.calculations['avg_daily_flow_m3']
        bod_load_fixed_film_kg = self.sizing['secondary_bod_load_kg'] * crit['load_split_fixed_film']
        media_area_req = (bod_load_fixed_film_kg * 1000) / crit['bod_salr'] if crit['bod_salr'] > 0 else 0
        s0_total = self.calculations['influent_bod'] * (1 - self.design_criteria['primary']['bod_removal'])
        s0_suspended = s0_total * (1 - crit['load_split_fixed_film'])
        numerator = q * crit['y'] * (s0_suspended - crit['s_effluent']) * crit['theta_c']
        denominator = crit['mlss_conc'] * (1 + crit['kd'] * crit['theta_c'])
        volume_for_suspended = numerator / denominator if denominator > 0 else 0
        volume_for_fixed = media_area_req / (crit['media_ssa'] * crit['media_fill']) if (crit['media_ssa'] * crit['media_fill']) > 0 else 0
        total_volume = max(volume_for_suspended, volume_for_fixed)
        self.sizing['secondary'] = {
            'type': 'Integrated Fixed-Film Activated Sludge (IFAS)', 'reactor_volume_m3': total_volume,
            'reactor_hrt_hr': (total_volume / q) * 24 if q > 0 else 0,
            'media_surface_area_req_m2': media_area_req,
            'organic_loading_rate_kg_m3_d': self.sizing['secondary_bod_load_kg'] / total_volume if total_volume > 0 else 0,
            'secondary_clarifier_area_m2': q / crit['clarifier_hlr'],
            'secondary_clarifier_hlr_m3_m2_d': crit['clarifier_hlr'],
        }

    def get_display_value(self, metric_value, unit_type):
        if metric_value is None or math.isinf(metric_value) or math.isnan(metric_value):
            return "N/A"
        val = metric_value
        if self.units == 'US':
            if unit_type == 'flow': val *= M3_TO_MGD
            elif unit_type == 'load': val *= KG_TO_LBS
            elif unit_type == 'area': val *= M2_TO_FT2
            elif unit_type == 'depth': val *= M_TO_FT
            elif unit_type == 'volume_small': val *= (M_TO_FT**3)
            elif unit_type == 'hlr': val *= M3_M2_D_TO_GPD_FT2
            elif unit_type == 'olr': val *= KG_M3_D_TO_LBS_1000FT3_D
            elif unit_type == 'flux': val *= (24.5424 / 1000 * 24) # LMH to GFD
        elif self.units == 'Imperial':
            if unit_type == 'flow': val *= M3_TO_IGD
            elif unit_type == 'load': val *= KG_TO_LBS
            elif unit_type == 'area': val *= M2_TO_FT2
            elif unit_type == 'depth': val *= M_TO_FT
            elif unit_type == 'volume_small': val *= (M_TO_FT**3)
            elif unit_type == 'hlr': val *= M3_M2_D_TO_IGPD_FT2
            elif unit_type == 'olr': val *= KG_M3_D_TO_LBS_1000FT3_D
            elif unit_type == 'flux': val *= (20.4355 / 1000 * 24) # LMH to IGPD
        return f"{val:,.2f}"

    def get_csv_data(self):
        """Generates the full results as a CSV formatted string."""
        output = io.StringIO()
        writer = csv.writer(output)
        headers = ["Category", "Parameter", "Value", "Units"]
        writer.writerow(headers)
        
        data = []
        data.append(["Input Parameters", "Unit System", self.units, ""])
        data.append(["Input Parameters", "Treatment Technology", self.treatment_type, ""])
        for key, value in self.params.items():
            data.append(["Input Parameters", key.replace('_', ' ').title(), value, ""])
        calc = self.calculations
        data.append(["Design Basis", "Projected Population", f"{calc.get('projected_population', 0):.0f}", "people"])
        data.append(["Design Basis", "Average Daily Flow", self.get_display_value(calc.get('avg_daily_flow_m3'), 'flow'), self.unit_config.get('flow')])
        data.append(["Design Basis", "Peak Factor", f"{calc.get('peak_factor', 0):.2f}", ""])
        data.append(["Design Basis", "BOD Load", self.get_display_value(calc.get('bod_load_kg'), 'load'), self.unit_config.get('load')])
        data.append(["Design Basis", "TSS Load", self.get_display_value(calc.get('tss_load_kg'), 'load'), self.unit_config.get('load')])
        data.append(["Design Basis", "TKN Load", self.get_display_value(calc.get('tkn_load_kg'), 'load'), self.unit_config.get('load')])
        
        prim = self.sizing.get('primary', {})
        sec = self.sizing.get('secondary', {})
        data.append(["Primary Sizing", "Clarifier Area", self.get_display_value(prim.get('clarifier_area_m2'), 'area'), self.unit_config.get('area')])
        data.append(["Primary Sizing", "HLR", self.get_display_value(prim.get('hydraulic_loading_rate_m3_m2_d'), 'hlr'), self.unit_config.get('hlr')])
        data.append(["Secondary Sizing", "Type", sec.get('type', 'N/A'), ""])
        vol_key = 'reactor_volume_m3' if 'reactor_volume_m3' in sec else 'aeration_volume_m3'
        if vol_key in sec:
            data.append(["Secondary Sizing", "Reactor Volume", self.get_display_value(sec.get(vol_key), 'volume_small'), self.unit_config.get('volume_small')])
        
        writer.writerows(data)
        return output.getvalue()

# --- Streamlit App UI ---
st.set_page_config(page_title="WWTP Design Calculator", layout="wide")
st.title("Wastewater Treatment Plant Preliminary Design Calculator")

# Instantiate calculator and store in session state
if 'calculator' not in st.session_state:
    st.session_state.calculator = WWTP_Design_Calculator()

calc = st.session_state.calculator

# --- Sidebar for Inputs ---
st.sidebar.header("Design Settings")

# General Settings
calc.set_units(st.sidebar.selectbox("Unit System", ["Metric", "US", "Imperial"]))
calc.treatment_type = st.sidebar.selectbox("Secondary Treatment Technology", ["CAS", "MBR", "MBBR", "IFAS"])

# Population and Flow Parameters
st.sidebar.subheader("Influent Characteristics")
calc.params['current_population'] = st.sidebar.number_input("Current Population", value=calc.params['current_population'], step=1000)
calc.params['design_life'] = st.sidebar.slider("Design Life (years)", 5, 50, calc.params['design_life'])
calc.params['growth_rate'] = st.sidebar.slider("Annual Growth Rate (%)", 0.0, 5.0, calc.params['growth_rate'], 0.1)
wc_unit = calc.unit_config['consumption']
calc.params['water_consumption'] = st.sidebar.number_input(f"Avg. Water Consumption ({wc_unit})", value=calc.params['water_consumption'])
calc.params['wastewater_fraction'] = st.sidebar.slider("Wastewater Fraction", 0.0, 1.0, calc.params['wastewater_fraction'], 0.05)
calc.params['bod_per_capita'] = st.sidebar.number_input("BOD per Capita (g/c/d)", value=calc.params['bod_per_capita'])
calc.params['tss_per_capita'] = st.sidebar.number_input("TSS per Capita (g/c/d)", value=calc.params['tss_per_capita'])
calc.params['tkn_per_capita'] = st.sidebar.number_input("TKN per Capita (g/c/d)", value=calc.params['tkn_per_capita'])

# Advanced Criteria Expander
with st.sidebar.expander("Edit Advanced Design Criteria"):
    for tech, criteria in calc.design_criteria.items():
        st.subheader(f"{tech.upper()} Criteria")
        for key, value in criteria.items():
            new_val = st.number_input(f"{tech} - {key}", value=float(value), key=f"{tech}_{key}", format="%.4f")
            calc.design_criteria[tech][key] = new_val

# --- Main Page for Results ---
if st.sidebar.button("Run Calculation", type="primary"):
    calc.run_calculations()
    st.session_state.results_generated = True

if 'results_generated' in st.session_state and st.session_state.results_generated:
    st.header("Design Results")
    
    # --- Design Basis ---
    st.subheader("1. Design Basis")
    c = calc.calculations
    col1, col2, col3 = st.columns(3)
    col1.metric("Projected Population", f"{c.get('projected_population', 0):,.0f} people")
    col2.metric(f"Average Daily Flow ({calc.unit_config['flow']})", calc.get_display_value(c.get('avg_daily_flow_m3'), 'flow'))
    col3.metric(f"Peak Hourly Flow ({calc.unit_config['flow']})", f"{calc.get_display_value(c.get('peak_hourly_flow_m3_d'), 'flow')} (PF: {c.get('peak_factor',0):.2f})")
    
    col1, col2, col3 = st.columns(3)
    col1.metric(f"BOD Load ({calc.unit_config['load']})", calc.get_display_value(c.get('bod_load_kg'), 'load'))
    col2.metric(f"TSS Load ({calc.unit_config['load']})", calc.get_display_value(c.get('tss_load_kg'), 'load'))
    col3.metric(f"TKN Load ({calc.unit_config['load']})", calc.get_display_value(c.get('tkn_load_kg'), 'load'))

    # --- Sizing ---
    st.subheader("2. Equipment Sizing & Loading")
    prim = calc.sizing.get('primary', {})
    sec = calc.sizing.get('secondary', {})

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("#### Primary Treatment")
        st.text(f"  Clarifier Area:   {calc.get_display_value(prim.get('clarifier_area_m2'), 'area')} {calc.unit_config['area']}")
        st.text(f"  Clarifier Volume: {calc.get_display_value(prim.get('clarifier_volume_m3'), 'volume_small')} {calc.unit_config['volume_small']}")
        st.text(f"  HLR:              {calc.get_display_value(prim.get('hydraulic_loading_rate_m3_m2_d'), 'hlr')} {calc.unit_config['hlr']}")

    with col2:
        st.markdown(f"#### Secondary Treatment ({sec.get('type', 'N/A')})")
        vol_key = 'reactor_volume_m3' if 'reactor_volume_m3' in sec else 'aeration_volume_m3'
        if vol_key in sec:
            st.text(f"  Reactor Volume:   {calc.get_display_value(sec.get(vol_key), 'volume_small')} {calc.unit_config['volume_small']}")
        if 'oxygen_demand_kg_d' in sec:
             st.text(f"  Oxygen Demand:    {calc.get_display_value(sec.get('oxygen_demand_kg_d'), 'load')} {calc.unit_config['load']}")
        if 'membrane_area_m2' in sec:
            st.text(f"  Membrane Area:    {calc.get_display_value(sec.get('membrane_area_m2'), 'area')} {calc.unit_config['area']}")
        if 'media_surface_area_req_m2' in sec:
            st.text(f"  Media Area Req'd: {calc.get_display_value(sec.get('media_surface_area_req_m2'), 'area')} {calc.unit_config['area']}")
        if 'secondary_clarifier_area_m2' in sec:
            st.text(f"  Sec. Clarifier Area: {calc.get_display_value(sec.get('secondary_clarifier_area_m2'), 'area')} {calc.unit_config['area']}")

    st.download_button(
       label="Download Full Results as CSV",
       data=calc.get_csv_data(),
       file_name='WWTP_Design_Results.csv',
       mime='text/csv',
    )
else:
    st.info("Adjust the settings in the sidebar and click 'Run Calculation' to see the results.")
