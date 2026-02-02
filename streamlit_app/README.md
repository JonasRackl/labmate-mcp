# Labmate Web App 🧪

A Streamlit web interface for labmate-mcp chemistry tools.

## Features

- 📷 **Image → SMILES** — Upload a structure image and extract SMILES using DECIMER AI
- 🔬 **Compound Lookup** — Search PubChem by name, CAS, SMILES, or formula
- 📚 **Literature Search** — Search Semantic Scholar for papers
- ⚗️ **Named Reactions** — Browse named organic reactions with conditions
- 🧮 **Calculations** — Molarity, dilution, MW, percent yield
- 🔮 **Predictions** — pKa, NMR, solubility (requires API keys)

## Run Locally
```bash
pip install -r requirements.txt
streamlit run app.py
```

## Deploy to Streamlit Cloud (Free)

1. Push to GitHub (done!)
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Connect repo → select `streamlit_app/app.py`
4. Deploy!

---

Part of [labmate-mcp](https://github.com/JonasRackl/labmate-mcp) · MIT License
