"""
labmate-mcp Streamlit Web App
Chemistry tools accessible via browser
"""

import streamlit as st
import requests
from io import BytesIO
import base64
import re
import os

# Page config
st.set_page_config(
    page_title="Labmate",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .stApp { max-width: 1200px; margin: 0 auto; }
    .result-box { 
        background-color: #f0f2f6; 
        padding: 1rem; 
        border-radius: 0.5rem; 
        margin: 1rem 0;
    }
</style>
""", unsafe_allow_html=True)

# Sidebar navigation
st.sidebar.title("🧪 Labmate")
st.sidebar.markdown("*Your AI lab companion*")

tool_category = st.sidebar.radio(
    "Category",
    ["📷 Image → SMILES", "🔬 Compound Lookup", "📚 Literature Search", 
     "⚗️ Named Reactions", "🧮 Calculations", "🔮 Predictions"]
)

# API key management
with st.sidebar.expander("🔑 API Keys"):
    rowan_key = st.text_input("Rowan API Key", type="password", 
                              value=os.environ.get("ROWAN_API_KEY", ""))
    rxn_key = st.text_input("IBM RXN API Key", type="password",
                            value=os.environ.get("RXN_API_KEY", ""))
    st.caption("Keys are not stored. Get free keys at rowan.ai and rxn.res.ibm.com")

# Session state
if 'smiles' not in st.session_state:
    st.session_state.smiles = ""


# ============================================================
# IMAGE TO SMILES (Featured!)
# ============================================================
if tool_category == "📷 Image → SMILES":
    st.title("📷 Image → SMILES → Predictions")
    st.markdown("Upload a chemical structure image (photo, screenshot, hand-drawn) and convert to SMILES")
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("1. Upload Image")
        uploaded_file = st.file_uploader(
            "Choose an image", 
            type=['png', 'jpg', 'jpeg', 'webp'],
            help="Works with printed structures, screenshots, and hand-drawn molecules"
        )
        
        if uploaded_file:
            st.image(uploaded_file, caption="Uploaded structure", use_container_width=True)
            
            if st.button("🔍 Extract SMILES", type="primary"):
                with st.spinner("Analyzing structure with DECIMER..."):
                    try:
                        # Using Hugging Face DECIMER model
                        API_URL = "https://api-inference.huggingface.co/models/yuvraj17/DECIMER-V2"
                        
                        response = requests.post(
                            API_URL,
                            data=uploaded_file.getvalue(),
                            timeout=60
                        )
                        
                        if response.status_code == 200:
                            result = response.json()
                            if isinstance(result, list) and len(result) > 0:
                                smiles = result[0].get("generated_text", "")
                                st.session_state.smiles = smiles
                                st.success(f"✅ Extracted SMILES!")
                            else:
                                st.error("Could not extract structure")
                        elif response.status_code == 503:
                            st.warning("⏳ Model is loading... try again in 20 seconds")
                        else:
                            st.error(f"API error: {response.status_code}")
                            
                    except Exception as e:
                        st.error(f"Error: {str(e)}")
                        st.info("Try manual SMILES input below")
    
    with col2:
        st.subheader("2. Results & Predictions")
        
        # Manual input fallback
        manual_smiles = st.text_input(
            "SMILES (extracted or manual):",
            value=st.session_state.smiles,
            placeholder="CC(=O)OC1=CC=CC=C1C(=O)O"
        )
        
        if manual_smiles:
            st.session_state.smiles = manual_smiles
            st.code(manual_smiles)
            
            # Show structure from PubChem
            try:
                img_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{requests.utils.quote(manual_smiles)}/PNG"
                st.image(img_url, caption="Structure (PubChem)", width=200)
            except:
                pass
            
            st.markdown("**🚀 Use this SMILES:**")
            
            c1, c2 = st.columns(2)
            with c1:
                if st.button("🔬 Lookup compound"):
                    st.session_state['goto'] = 'lookup'
                    st.rerun()
                if st.button("🔮 Predict pKa"):
                    st.session_state['goto'] = 'predict'
            with c2:
                if st.button("⚗️ Retrosynthesis"):
                    st.session_state['goto'] = 'retro'
                if st.button("📋 Copy SMILES"):
                    st.code(manual_smiles)


# ============================================================
# COMPOUND LOOKUP
# ============================================================
elif tool_category == "🔬 Compound Lookup":
    st.title("🔬 Compound Lookup")
    st.markdown("Search by name, CAS number, SMILES, or formula")
    
    query = st.text_input(
        "Enter compound:", 
        value=st.session_state.smiles if st.session_state.get('goto') == 'lookup' else "",
        placeholder="e.g., aspirin, 50-78-2, CC(=O)OC1=CC=CC=C1C(=O)O"
    )
    
    if st.button("Search", type="primary") and query:
        with st.spinner("Searching PubChem..."):
            try:
                # Try by name first, then SMILES
                if query.startswith("C") and any(c in query for c in "()=#"):
                    # Looks like SMILES
                    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{requests.utils.quote(query)}/property/Title,MolecularFormula,MolecularWeight,CanonicalSMILES,InChIKey,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount/JSON"
                else:
                    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{query}/property/Title,MolecularFormula,MolecularWeight,CanonicalSMILES,InChIKey,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount/JSON"
                
                response = requests.get(url, timeout=10)
                
                if response.status_code == 200:
                    props = response.json()['PropertyTable']['Properties'][0]
                    cid = props.get('CID')
                    
                    col1, col2 = st.columns([1, 2])
                    
                    with col1:
                        if cid:
                            img_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG"
                            st.image(img_url, caption="2D Structure", width=250)
                    
                    with col2:
                        st.subheader(props.get('Title', query))
                        st.markdown(f"**Formula:** {props.get('MolecularFormula', 'N/A')}")
                        st.markdown(f"**MW:** {props.get('MolecularWeight', 'N/A')} g/mol")
                        st.markdown(f"**SMILES:** `{props.get('CanonicalSMILES', 'N/A')}`")
                        st.markdown(f"**InChIKey:** `{props.get('InChIKey', 'N/A')}`")
                        
                        st.divider()
                        st.markdown("**Lipinski Parameters**")
                        
                        metrics = st.columns(4)
                        with metrics[0]:
                            xlogp = props.get('XLogP', 'N/A')
                            st.metric("XLogP", xlogp, "✓" if isinstance(xlogp, (int, float)) and xlogp < 5 else "")
                        with metrics[1]:
                            st.metric("MW", f"{props.get('MolecularWeight', 0):.0f}")
                        with metrics[2]:
                            st.metric("HBD", props.get('HBondDonorCount', 'N/A'))
                        with metrics[3]:
                            st.metric("HBA", props.get('HBondAcceptorCount', 'N/A'))
                else:
                    st.warning(f"Compound '{query}' not found")
                    
            except Exception as e:
                st.error(f"Error: {str(e)}")


# ============================================================
# LITERATURE SEARCH
# ============================================================
elif tool_category == "📚 Literature Search":
    st.title("📚 Literature Search")
    st.markdown("Search papers via Semantic Scholar")
    
    query = st.text_input("Search query:", placeholder="e.g., asymmetric catalysis Maruoka")
    
    col1, col2, col3 = st.columns(3)
    with col1:
        year_from = st.number_input("Year from:", min_value=1900, max_value=2026, value=2020)
    with col2:
        year_to = st.number_input("Year to:", min_value=1900, max_value=2026, value=2026)
    with col3:
        limit = st.number_input("Results:", min_value=5, max_value=50, value=10)
    
    if st.button("Search", type="primary") and query:
        with st.spinner("Searching Semantic Scholar..."):
            try:
                url = "https://api.semanticscholar.org/graph/v1/paper/search"
                params = {
                    "query": query,
                    "limit": limit,
                    "fields": "title,authors,year,citationCount,openAccessPdf,tldr,url,externalIds",
                    "year": f"{year_from}-{year_to}"
                }
                
                response = requests.get(url, params=params, timeout=15)
                
                if response.status_code == 200:
                    papers = response.json().get('data', [])
                    
                    if papers:
                        st.success(f"Found {len(papers)} papers")
                        
                        for paper in papers:
                            authors = ", ".join([a['name'] for a in paper.get('authors', [])[:3]])
                            if len(paper.get('authors', [])) > 3:
                                authors += " et al."
                            
                            with st.expander(f"📄 {paper.get('title', 'No title')} ({paper.get('year', 'N/A')})"):
                                st.markdown(f"**Authors:** {authors}")
                                st.markdown(f"**Citations:** {paper.get('citationCount', 0)}")
                                
                                if paper.get('tldr'):
                                    st.info(f"**TL;DR:** {paper['tldr'].get('text', '')}")
                                
                                doi = paper.get('externalIds', {}).get('DOI')
                                if doi:
                                    st.markdown(f"**DOI:** [{doi}](https://doi.org/{doi})")
                                
                                if paper.get('openAccessPdf'):
                                    st.markdown(f"[📥 Open Access PDF]({paper['openAccessPdf']['url']})")
                    else:
                        st.warning("No papers found")
                else:
                    st.error(f"API error: {response.status_code}")
                    
            except Exception as e:
                st.error(f"Error: {str(e)}")


# ============================================================
# NAMED REACTIONS
# ============================================================
elif tool_category == "⚗️ Named Reactions":
    st.title("⚗️ Named Reactions Database")
    st.markdown("Browse 200+ named organic reactions")
    
    REACTIONS = {
        "Suzuki Coupling": {
            "category": "Cross-Coupling",
            "substrate": "Aryl/vinyl halide + boronic acid",
            "catalyst": "Pd(PPh₃)₄, Pd(dppf)Cl₂",
            "base": "K₂CO₃, Cs₂CO₃, K₃PO₄",
            "solvent": "THF, dioxane, DMF/H₂O",
            "temp": "60-100°C",
            "tips": "Degassing is critical. Use fresh boronic acid."
        },
        "Buchwald-Hartwig": {
            "category": "Cross-Coupling",
            "substrate": "Aryl halide + amine",
            "catalyst": "Pd₂(dba)₃ + ligand (XPhos, BrettPhos)",
            "base": "NaOtBu, Cs₂CO₃",
            "solvent": "Toluene, dioxane",
            "temp": "80-110°C",
            "tips": "Ligand choice crucial. BrettPhos for 1° amines."
        },
        "Grignard Reaction": {
            "category": "Addition",
            "substrate": "Carbonyl + RMgX",
            "catalyst": "None (stoichiometric)",
            "base": "N/A",
            "solvent": "Et₂O, THF (anhydrous!)",
            "temp": "0°C to reflux",
            "tips": "Absolutely dry conditions. Activate Mg with I₂."
        },
        "Wittig Reaction": {
            "category": "Olefination",
            "substrate": "Aldehyde/ketone + phosphonium ylide",
            "catalyst": "None",
            "base": "nBuLi, NaHMDS",
            "solvent": "THF, Et₂O",
            "temp": "-78°C to RT",
            "tips": "Stabilized ylides → E. Non-stabilized → Z."
        },
        "Diels-Alder": {
            "category": "Cycloaddition",
            "substrate": "Diene + dienophile",
            "catalyst": "Lewis acids optional",
            "base": "N/A",
            "solvent": "Toluene, CH₂Cl₂, neat",
            "temp": "RT to 200°C",
            "tips": "Endo rule for stereochemistry."
        },
        "Swern Oxidation": {
            "category": "Oxidation",
            "substrate": "Alcohol → aldehyde/ketone",
            "catalyst": "None",
            "base": "Et₃N",
            "solvent": "CH₂Cl₂",
            "temp": "-78°C (critical!)",
            "tips": "Keep cold until Et₃N added. No overoxidation."
        },
        "Mitsunobu": {
            "category": "Substitution",
            "substrate": "Alcohol + nucleophile",
            "catalyst": "None",
            "base": "DEAD/DIAD + PPh₃",
            "solvent": "THF, CH₂Cl₂",
            "temp": "0°C to RT",
            "tips": "pKa of nucleophile < 11. Clean inversion."
        },
        "Heck Reaction": {
            "category": "Cross-Coupling",
            "substrate": "Aryl halide + alkene",
            "catalyst": "Pd(OAc)₂, Pd(PPh₃)₄",
            "base": "Et₃N, K₂CO₃",
            "solvent": "DMF, MeCN",
            "temp": "80-140°C",
            "tips": "β-hydride elimination gives alkene."
        },
    }
    
    search = st.text_input("Search:", placeholder="e.g., Suzuki, coupling, oxidation")
    
    filtered = {k: v for k, v in REACTIONS.items() 
                if not search or search.lower() in k.lower() or search.lower() in v['category'].lower()}
    
    st.caption(f"Showing {len(filtered)} reactions")
    
    for name, rxn in filtered.items():
        with st.expander(f"**{name}** ({rxn['category']})"):
            col1, col2 = st.columns(2)
            with col1:
                st.markdown(f"**Substrate:** {rxn['substrate']}")
                st.markdown(f"**Catalyst:** {rxn['catalyst']}")
                st.markdown(f"**Base:** {rxn['base']}")
            with col2:
                st.markdown(f"**Solvent:** {rxn['solvent']}")
                st.markdown(f"**Temperature:** {rxn['temp']}")
            st.info(f"💡 {rxn['tips']}")


# ============================================================
# CALCULATIONS
# ============================================================
elif tool_category == "🧮 Calculations":
    st.title("🧮 Chemistry Calculations")
    
    calc_type = st.selectbox("Calculator:", [
        "Molarity",
        "Dilution (C₁V₁ = C₂V₂)",
        "Molecular Weight",
        "Percent Yield"
    ])
    
    if calc_type == "Molarity":
        col1, col2, col3 = st.columns(3)
        with col1:
            mass = st.number_input("Mass (g):", min_value=0.0, format="%.4f")
        with col2:
            mw = st.number_input("MW (g/mol):", min_value=0.0, format="%.2f")
        with col3:
            volume = st.number_input("Volume (L):", min_value=0.0, format="%.4f")
        
        if st.button("Calculate", type="primary") and mw > 0 and volume > 0:
            moles = mass / mw
            molarity = moles / volume
            st.success(f"**Molarity: {molarity:.4f} M** ({moles:.4f} mol)")
    
    elif calc_type == "Dilution (C₁V₁ = C₂V₂)":
        col1, col2 = st.columns(2)
        with col1:
            c1 = st.number_input("C₁:", min_value=0.0, format="%.4f")
            v1 = st.number_input("V₁:", min_value=0.0, format="%.4f")
        with col2:
            c2 = st.number_input("C₂:", min_value=0.0, format="%.4f")
            v2 = st.number_input("V₂:", min_value=0.0, format="%.4f")
        
        solve = st.radio("Solve for:", ["V₂", "V₁", "C₂"], horizontal=True)
        
        if st.button("Calculate", type="primary"):
            if solve == "V₂" and c2 > 0:
                st.success(f"**V₂ = {(c1*v1)/c2:.4f}**")
            elif solve == "V₁" and c1 > 0:
                st.success(f"**V₁ = {(c2*v2)/c1:.4f}**")
            elif solve == "C₂" and v2 > 0:
                st.success(f"**C₂ = {(c1*v1)/v2:.4f}**")
    
    elif calc_type == "Molecular Weight":
        WEIGHTS = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.065, 
                   'P': 30.974, 'F': 18.998, 'Cl': 35.453, 'Br': 79.904, 'I': 126.90}
        
        formula = st.text_input("Formula:", placeholder="C6H12O6")
        
        if st.button("Calculate", type="primary") and formula:
            matches = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
            mw = sum(WEIGHTS.get(el, 0) * (int(n) if n else 1) for el, n in matches if el)
            st.success(f"**MW: {mw:.3f} g/mol**")
    
    elif calc_type == "Percent Yield":
        col1, col2 = st.columns(2)
        with col1:
            actual = st.number_input("Actual (g):", min_value=0.0, format="%.4f")
        with col2:
            theoretical = st.number_input("Theoretical (g):", min_value=0.0, format="%.4f")
        
        if st.button("Calculate", type="primary") and theoretical > 0:
            pct = (actual / theoretical) * 100
            st.success(f"**Yield: {pct:.1f}%**")


# ============================================================
# PREDICTIONS
# ============================================================
elif tool_category == "🔮 Predictions":
    st.title("🔮 Property Predictions")
    
    smiles = st.text_input(
        "SMILES:", 
        value=st.session_state.smiles,
        placeholder="c1ccccc1O"
    )
    
    pred_type = st.selectbox("Prediction:", [
        "pKa", "NMR shifts", "Solubility", "Retrosynthesis"
    ])
    
    if st.button("Predict", type="primary") and smiles:
        if pred_type == "Retrosynthesis" and not rxn_key:
            st.warning("⚠️ Requires IBM RXN API key")
            st.info("Get free key: [rxn.res.ibm.com](https://rxn.res.ibm.com)")
        elif pred_type in ["pKa", "NMR shifts", "Solubility"] and not rowan_key:
            st.warning("⚠️ Requires Rowan API key")
            st.info("Get free key: [rowan.ai](https://rowan.ai)")
        else:
            st.info(f"🚧 {pred_type} prediction coming soon! Use CLI: `labmate-mcp`")


# Footer
st.sidebar.divider()
st.sidebar.markdown("""
**labmate-mcp** v7.3.1  
[GitHub](https://github.com/JonasRackl/labmate-mcp) · [PyPI](https://pypi.org/project/labmate-mcp/)
""")
