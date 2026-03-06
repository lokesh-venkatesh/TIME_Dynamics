if __name__ == '__main__':
    os.makedirs('data', exist_ok=True)
    os.makedirs('data/more_data', exist_ok=True)
    os.makedirs("data/figure_3i-m", exist_ok=True)
    p1p2_simulation_results = run_sweep({'p1': np.linspace(0.2, 0.4, 25),
                                         'p2': np.linspace(0.05, 0.08, 25)}, output_prefix='sweep_p1p2')

    # Example usage:
    
    # 1. Run with default parameters
    # df_default = run_default()
    # df_default.to_csv('simulation_results.csv', index=False)
    # print(f"Default run completed, saved to results_default.csv")
    
    # 2. Modify specific parameters BY NAME
    # df_mod = run_modified(['p1', 'p2'], [0.3, 0.1])
    # df_mod = run_modified(['p1'], [0.25], t_final=1000)
    # df_mod.to_csv('results_for_figure_2i.csv', index=False)
    
    # 3. Parameter sweep by name (single parameter)
    # results = run_sweep({'p1': np.linspace(0.2, 0.4, 100)}, output_prefix='sweep_p1')
    # print(results)
    
    # 4. Parameter sweep (multiple parameters)
    # results = run_sweep({'betaM2': np.linspace(1e-16, 1e-14, 3), 
    #                      'betaTc': np.linspace(1e-9, 1e-7, 3)}, 
    #                     output_prefix='sweep_multi')

    # Extract p1 and p2 values and create grid
    p1_values = sorted(set([key[0][1] for key in p1p2_simulation_results.keys()]))
    p2_values = sorted(set([key[1][1] for key in p1p2_simulation_results.keys()]))

    print(f"p1 values: {len(p1_values)}, p2 values: {len(p2_values)}")

    # Create 2D arrays for each species at steady state
    S_map = np.zeros((len(p1_values), len(p2_values)))
    SR_map = np.zeros((len(p1_values), len(p2_values)))
    C_map = np.zeros((len(p1_values), len(p2_values)))
    CR_map = np.zeros((len(p1_values), len(p2_values)))

    # Fill the maps with final steady-state values
    for key, df in p1p2_simulation_results.items():
        p1_idx = p1_values.index(key[0][1])
        p2_idx = p2_values.index(key[1][1])
        
        S_map[p1_idx, p2_idx] = df['S'].iloc[-1]
        SR_map[p1_idx, p2_idx] = df['SR'].iloc[-1]
        C_map[p1_idx, p2_idx] = df['C'].iloc[-1]
        CR_map[p1_idx, p2_idx] = df['CR'].iloc[-1]

    """FILLED IN CODE IS BELOW"""

    # Convert to pandas DataFrames with proper indexing for neat labels
    S_df = pd.DataFrame(S_map, index=np.round(p1_values, 4), columns=np.round(p2_values, 4))
    SR_df = pd.DataFrame(SR_map, index=np.round(p1_values, 4), columns=np.round(p2_values, 4))
    C_df = pd.DataFrame(C_map, index=np.round(p1_values, 4), columns=np.round(p2_values, 4))
    CR_df = pd.DataFrame(CR_map, index=np.round(p1_values, 4), columns=np.round(p2_values, 4))

    dfs_data = [
        (S_df, 'S', '3a'),
        (SR_df, 'SR', '3b'),
        (C_df, 'C', '3c'),
        (CR_df, 'CR', '3d')
    ]

    for df_data, species_name, figure_label in dfs_data:
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Get min and max for this specific plot
        vmin = df_data.values.min()
        vmax = df_data.values.max()
        
        # Create heatmap with individual color scaling
        sns.heatmap(df_data, ax=ax, cmap='viridis', vmin=vmin, vmax=vmax, 
                    cbar=True, cbar_kws={'label': 'Concentration', 
                                        'ticks': np.linspace(vmin, vmax, 5)})
        
        # Format colorbar tick labels to show min/max with appropriate precision
        cbar = ax.collections[0].colorbar
        cbar.ax.set_yticklabels([f'{val:.2e}' for val in np.linspace(vmin, vmax, 5)])
        
        ax.set_xlabel('p2', fontsize=12, fontweight='bold')
        ax.set_ylabel('p1', fontsize=12, fontweight='bold')
        ax.set_title(species_name, fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(f'data/figure_3a-d/figure_{figure_label}.png', dpi=300, bbox_inches='tight')
        plt.show()

    if os.path.exists('__pycache__'):
        shutil.rmtree('__pycache__')