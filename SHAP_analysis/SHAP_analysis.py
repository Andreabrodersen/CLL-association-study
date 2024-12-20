#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 14:25:26 2024

@author: AmoBro (Andrea Moen Brodersen)
"""



import os
import pandas as pd
import matplotlib.pyplot as plt
import shap
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import seaborn as sns
# Load the dataset (replace the directory with your own)
rf_df_flow_single = pd.read_csv('/Rf_df_single_nona_177.csv', index_col=0)


# Set the column ranges
start_col = 41
end_col = len(rf_df_flow_single.columns)

# Create a folder to save the plots(add your own directory)
os.makedirs("/Test_plots_shap", exist_ok=True)
os.makedirs("/Train_plots_shap", exist_ok=True)
# Initialize an empty DataFrame to store variable importance data
df_ImpData_single = pd.DataFrame()
all_shap_values = pd.DataFrame()
all_shap_values_train = pd.DataFrame()


for col_index in range(start_col, end_col):
    # Get the column name for the current col_index
    col_name = rf_df_flow_single.columns[col_index]
    
    # Create random forest model
    X = rf_df_flow_single.iloc[:,0:41] # Features
    y = rf_df_flow_single[col_name]  # Target
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    rf_fit = RandomForestRegressor(n_estimators=500)
    rf_fit.fit(X_train, y_train)
    
    # Calculate SHAP values
    
    
    explainer = shap.TreeExplainer(rf_fit)
    shap_values = explainer.shap_values(X_test)
    shap_values_train = explainer.shap_values(X_train)
    
    # Save SHAP values to a DataFrame
    shap_df = pd.DataFrame(shap_values, columns=X_test.columns)
    shap_df["Drug"] = col_name
    
    #Shap values train
    shap_df_train = pd.DataFrame(shap_values_train, columns=X_train.columns)
    shap_df_train["Drug"] = col_name
    
    # Append SHAP values to the main DataFrame
    shap_df = pd.DataFrame(shap_values, columns=X_test.columns)
    shap_df["Drug"] = col_name
    
    # Append SHAP values to the main DataFrame_train
    shap_df_train = pd.DataFrame(shap_values_train, columns=X_train.columns)
    shap_df_train["Drug"] = col_name
    
    # Append SHAP values to the main DataFrame
    all_shap_values = pd.concat([all_shap_values, shap_df])
    
    # Append SHAP values to the main DataFrame_train
    all_shap_values_train = pd.concat([all_shap_values_train, shap_df_train])

    # Save SHAP values to a CSV file (add your own directory)
    csv_file_path = "/all_shap_values.csv"
    all_shap_values.to_csv(csv_file_path, index=False)
    
    csv_file_path = "/all_shap_values_train.csv"
    all_shap_values_train.to_csv(csv_file_path, index=False)
    
    # Create and save SHAP summary plot (add your own directory path after //)
    plt.figure(figsize=( 6, 12))
    shap.summary_plot(shap_values, X_test, max_display=42, show=False)
    plt.savefig(f"//Test_plots_shap/shap_summary_plot_{col_name}.pdf")
    plt.close()
    
    # Create and save SHAP summary plot
    plt.figure(figsize=( 6, 12))
    shap.summary_plot(shap_values_train, X_train, max_display=42, show=False)
    plt.savefig(f"//Train_plots_shap/shap_summary_plot_train_{col_name}.pdf")
    plt.close()
    
    # Combine SHAP values with variable importance data
    imp_data = pd.DataFrame(rf_fit.feature_importances_, index=X.columns, columns=["%IncMSE"])
    imp_data["Var.Names"] = imp_data.index
    imp_data["Drug"] = col_name
    imp_data = imp_data.sort_values(by="%IncMSE")
    imp_data["Var.Names"] = pd.Categorical(imp_data["Var.Names"], categories=imp_data["Var.Names"])
    
    # Create and save the plot #aqdd your own directory) 
    plt.figure(figsize=(8, 6))
    sns.barplot(x="%IncMSE", y="Var.Names", data=imp_data, palette="Blues")
    plt.title(f"Plot for {col_name}")
    plt.savefig(f"/Test_plots_shap20240802/plot_column_{col_name}.png")
    plt.close()
    
    # Append the variable importance data to the main DataFrame
    df_ImpData_single = pd.concat([df_ImpData_single, imp_data])
    
    # Visualize the direction of associations using SHAP values
    plt.figure(figsize=(8, 6))
    sns.barplot(x=imp_data["%IncMSE"], y=imp_data["Var.Names"], palette="Blues")
    plt.title(f"SHAP Values for {col_name}")
    plt.savefig(f"/Test_plots_shap20240802/shap_values_plot_{col_name}.png")
    plt.close()

# Group and summarize the variable importance data
df_ImpData_sum_single = df_ImpData_single.groupby("Var.Names").agg(
    Mean=("%IncMSE", "mean"),
    SD=("%IncMSE", "std"),
    Max=("%IncMSE", "max"),
    Which_max=("Drug", lambda x: x[df_ImpData_single["%IncMSE"].idxmax()])
).reset_index()
df_ImpData_sum_single["tval"] = df_ImpData_sum_single["Mean"] / df_ImpData_sum_single["SD"]










#for combos------
#read file (add your own directory path)

rf_df_flow_combo = pd.read_csv('/Rf_df_combo_nona_177.csv', index_col=0)


# Set the column ranges
start_col = 41
end_col = len(rf_df_flow_combo.columns)

# Create a folder to save the plots(add your preffered directory path) 
os.makedirs("/plots_combo_shap", exist_ok=True)
os.makedirs("/plots_combo_shap_train", exist_ok=True)


# Initialize an empty DataFrame to store variable importance data
df_ImpData_single = pd.DataFrame()
all_shap_values_combo = pd.DataFrame()
all_shap_values_combo_train = pd.DataFrame()

for col_index in range(start_col, end_col):
    # Get the column name for the current col_index
    col_name = rf_df_flow_combo.columns[col_index]
    
    # Create random forest model
    X = rf_df_flow_combo.iloc[:,0:41] # Features
    y = rf_df_flow_combo[col_name]  # Target
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    rf_fit = RandomForestRegressor(n_estimators=500)
    rf_fit.fit(X_train, y_train)
    
    # Calculate SHAP values
    explainer = shap.TreeExplainer(rf_fit)
    shap_values = explainer.shap_values(X_test)
    
    #train
    shap_values_train = explainer.shap_values(X_train)
    
    
    # Save SHAP values to a DataFrame
    shap_df = pd.DataFrame(shap_values, columns=X_test.columns)
    shap_df["Drug"] = col_name
    
    # Append SHAP values to the main DataFrame
    shap_df = pd.DataFrame(shap_values, columns=X_test.columns)
    shap_df["Drug"] = col_name
    
    # Append SHAP values to the main DataFrame
    all_shap_values = pd.concat([all_shap_values, shap_df])
    
    
    # Save SHAP values to a DataFrame_TRAIN
    shap_df_train = pd.DataFrame(shap_values, columns=X_train.columns)
    shap_df_train["Drug"] = col_name
    
    # Append SHAP values to the main DataFrame_TRAIN
    shap_df_train = pd.DataFrame(shap_values, columns=X_test.columns)
    shap_df_train["Drug"] = col_name
    
    # Append SHAP values to the main DataFrame_TRAIN
    all_shap_values_combo = pd.concat([all_shap_values_combo, shap_df])
    
    # Append SHAP values to the main DataFrame_TRAIN
    all_shap_values_combo_train = pd.concat([all_shap_values_combo_train, shap_df_train])

    # Save SHAP values to a CSV file (add your own directory)
    csv_file_path = "/all_shap_values_combo.csv"
    all_shap_values_combo.to_csv(csv_file_path, index=False)
    
    
    # Save SHAP values to a CSV file (add your own directory) 
    csv_file_path = "/all_shap_values_combo_train.csv"
    all_shap_values_combo_train.to_csv(csv_file_path, index=False)
    
    # Create and save SHAP summary plot(add your own directory)
    plt.figure(figsize=( 6, 12))
    shap.summary_plot(shap_values, X_test, max_display=42, show=False)
    plt.savefig(f"/plots_combo_shap/shap_summary_plot_{col_name}.pdf")
    plt.close()
    
    # Create and save SHAP summary plot (add your own directory)
    plt.figure(figsize=( 6, 12))
    shap.summary_plot(shap_values_train, X_train, max_display=42, show=False)
    plt.savefig(f"/plots_combo_shap_train/shap_summary_plot_{col_name}.pdf")
    plt.close()
    
    # Combine SHAP values with variable importance data
    imp_data = pd.DataFrame(rf_fit.feature_importances_, index=X.columns, columns=["%IncMSE"])
    imp_data["Var.Names"] = imp_data.index
    imp_data["Drug"] = col_name
    imp_data = imp_data.sort_values(by="%IncMSE")
    imp_data["Var.Names"] = pd.Categorical(imp_data["Var.Names"], categories=imp_data["Var.Names"])
    
    # Create and save the plot (add your own directory)
    plt.figure(figsize=(8, 6))
    sns.barplot(x="%IncMSE", y="Var.Names", data=imp_data, palette="Blues")
    plt.title(f"Plot for {col_name}")
    plt.savefig(f"/plots_combo_shap/plot_column_{col_name}.png")
    plt.close()
    
    # Append the variable importance data to the main DataFrame
    df_ImpData_single = pd.concat([df_ImpData_single, imp_data])
    
    # Visualize the direction of associations using SHAP values
    plt.figure(figsize=(8, 6))
    sns.barplot(x=imp_data["%IncMSE"], y=imp_data["Var.Names"], palette="Blues")
    plt.title(f"SHAP Values for {col_name}")
    plt.savefig(f"plots_combo_shap/shap_values_plot_{col_name}.png")
    plt.close()

# Group and summarize the variable importance data
df_ImpData_sum_single = df_ImpData_single.groupby("Var.Names").agg(
    Mean=("%IncMSE", "mean"),
    SD=("%IncMSE", "std"),
    Max=("%IncMSE", "max"),
    Which_max=("Drug", lambda x: x[df_ImpData_single["%IncMSE"].idxmax()])
).reset_index()
df_ImpData_sum_single["tval"] = df_ImpData_sum_single["Mean"] / df_ImpData_sum_single["SD"]
















