# python preprocess_chemprop.py --filepath ./CO2_Soleimani/augmentation/StratifiedKFold/input_t*.csv --column_names Augmented_SMILES,T_K,P_Mpa --target exp_CO2_sol_g_g --augment True

# python preprocess_chemprop.py --filepath ./CO2_Soleimani/manual_frag/StratifiedKFold/input_t*.csv --column_names Polymer_manual_aug_str,T_K,P_Mpa --target exp_CO2_sol_g_g --augment True

python preprocess_chemprop.py --filepath ./CO2_Soleimani/manual_frag_str/StratifiedKFold/input_t*.csv --column_names Polymer_manual_str,T_K,P_Mpa --target exp_CO2_sol_g_g --augment False

# python preprocess_chemprop.py --filepath ./CO2_Soleimani/manual_frag/StratifiedKFold/input_t*.csv --column_names Polymer_SMILES,T_K,P_Mpa --target exp_CO2_sol_g_g --augment False

# python preprocess_chemprop.py --filepath ./OPV_Min/augmentation/KFold/input_t*.csv --column_names Augmented_SMILES --target calc_PCE_percent --augment True

# python preprocess_chemprop.py --filepath ./OPV_Min/manual_frag/KFold/input_t*.csv --column_names DA_manual_aug_str --target calc_PCE_percent --augment True

python preprocess_chemprop.py --filepath ./OPV_Min/manual_frag_str/KFold/input_t*.csv --column_names DA_manual_str --target calc_PCE_percent --augment False

# python preprocess_chemprop.py --filepath ./OPV_Min/SMILES/KFold/input_t*.csv --column_names DA_SMILES --target calc_PCE_percent --augment False

# python preprocess_chemprop.py --filepath ./PV_Wang/SMILES/StratifiedKFold/input_t*.csv --column_names PS_SMILES,Contact_angle,Thickness_um,Solvent,Solvent_solubility_parameter_Mpa_sqrt,xw_wt_percent,Temp_C,Permeate_pressure_mbar --target J_Total_flux_kg_m_2h_1 --augment False

# python preprocess_chemprop.py --filepath ./PV_Wang/manual_frag/StratifiedKFold/input_t*.csv --column_names PS_manual_aug_str,Contact_angle,Thickness_um,Solvent,Solvent_solubility_parameter_Mpa_sqrt,xw_wt_percent,Temp_C,Permeate_pressure_mbar --target J_Total_flux_kg_m_2h_1 --augment True

python preprocess_chemprop.py --filepath ./PV_Wang/manual_frag_str/StratifiedKFold/input_t*.csv --column_names PS_manual_str,Contact_angle,Thickness_um,Solvent,Solvent_solubility_parameter_Mpa_sqrt,xw_wt_percent,Temp_C,Permeate_pressure_mbar --target J_Total_flux_kg_m_2h_1 --augment False

# python preprocess_chemprop.py --filepath ./PV_Wang/augmentation/StratifiedKFold/input_t*.csv --column_names Augmented_SMILES,Contact_angle,Thickness_um,Solvent,Solvent_solubility_parameter_Mpa_sqrt,xw_wt_percent,Temp_C,Permeate_pressure_mbar --target J_Total_flux_kg_m_2h_1 --augment True

# python preprocess_chemprop.py --filepath ./Swelling_Xu/augmentation/StratifiedKFold/input_t*.csv --column_names Augmented_SMILES --target SD --augment True

# python preprocess_chemprop.py --filepath ./Swelling_Xu/manual_frag/StratifiedKFold/input_t*.csv --column_names PS_manual_aug_str --target SD --augment True

# python preprocess_chemprop.py --filepath ./Swelling_Xu/SMILES/StratifiedKFold/input_t*.csv --column_names PS_SMILES --target SD --augment False

python preprocess_chemprop.py --filepath ./Swelling_Xu/manual_frag_str/StratifiedKFold/input_t*.csv --column_names PS_manual_str --target SD --augment False