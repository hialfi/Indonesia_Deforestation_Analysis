# Indonesia_Deforestation_Analysis

## Project Overview
This repository contains R scripts and datasets used to analyze drivers of deforestation in Indonesia (2008–2024) using PCA and Bayesian modeling.

## Data Description

### DATA/Indonesia_Deforestation_MergedData.csv
- **TreeCoverLoss**: Annual tree cover loss (hectares) from Global Forest Watch (UMD)
- **PopulationTotal**: Total population (World Bank WDI)
- **ForestArea**: Forest area in hectares (WDI)
- **ProtectedArea**: Protected forest area (WDI)
- **UrbanGrowth**: Urban population growth (%)
- **GDPGrowth**: GDP growth rate (%)
- **AgriValue**: Agricultural value added (% of GDP)
- **EnergyUsePerCapita**: Energy use per person
- **FDIInflow**: Foreign direct investment inflow (USD)
- **ODAUSD**: Official development assistance (USD)
- **Flood, Fire, Drought**: Number of verified natural disaster events per year
- **Budget**: Environmental budget from Statistics Indonesia (APBN)

> All variables are merged into a panel dataset covering 2008–2024.

### DATA/APBLingkungan20082024.csv
- Central government Environmental Function Budget (annual) in IDR
- Used in models incorporating fiscal policy as a predictor
