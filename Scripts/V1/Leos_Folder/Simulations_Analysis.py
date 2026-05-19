#!/usr/bin/python3
import gromacs
import subprocess
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import pyedr
import glob
import itertools
################################################################################
class Simulations_Analysis(object):
    def __init__(self):
        pass
################################################################################
    def get_dipole_model(self): 
    
        gromacs.environment.flags['capture_output'] = True
        input_str = '0\n'
        tmp = gromacs.tools.Dipoles(f=f'nvt_vacuum2',
                                s=f'nvt_vacuum2',
                                input=input_str)
        chk, dipole_output, stderr = tmp.run()

        f = open(f'Dipole.txt', 'w')
        f.write(dipole_output)
        f.close()

        dipole_lines = dipole_output.splitlines()  
        dipole_line = dipole_lines[8]  
        dipole_model = float(dipole_line.split()[2])
        return dipole_model
######################################################
    def get_dipole_model_liquid(self): 
    
        gromacs.environment.flags['capture_output'] = True
        input_str = '0\n'
        tmp = gromacs.tools.Dipoles(f=f'Pure_QMMM_md3.trr',
                                s=f'Pure_QMMM_md3.tpr',
                                input=input_str)
        chk, dipole_output, stderr = tmp.run()

        f = open(f'Dipole2.txt', 'w')
        f.write(dipole_output)
        f.close()

        dipole_lines = dipole_output.splitlines()  
        dipole_line = dipole_lines[8]  
        Model_Dipole = float(dipole_line.split()[2])  
        epsilon_line = dipole_lines[-1]  
        Epsilon = float(epsilon_line.split()[2])  
        
        return Model_Dipole, Epsilon
################################################################################
    def Sim_Density(self):
        gromacs.environment.flags['capture_output'] = True
        input_str = 'Density 0\n'
        tmp = gromacs.tools.Energy(f=f'Pure_QMMM_md3.edr',
                               s=f'Pure_QMMM_md3.tpr',
                               input=input_str)
        chk, density_output, stderr = tmp.run()

        f = open(f'Density.txt', 'w')
        f.write(density_output)
        f.close()

        density_lines = density_output.splitlines()
        density_line = density_lines[9]  
        density = float(density_line.split()[1])
        
        return density  
################################################################################
    def Sim_Enth(self,Gas_Potential, Liq_Potential, T, N):

        R = 8.314e-3 
        return Gas_Potential - (Liq_Potential / N) + R * T   
                
################################################################################
    def Pot_Gas(self):

        gas_edr = 'vacuum_potential.edr'
        gas_data = pyedr.edr_to_dict(gas_edr)
        equil_time_ps = 100
        gas_mask = gas_data['Time'] >= equil_time_ps
        
        return np.mean(gas_data['Potential'][gas_mask])     
################################################################################
    def Pot_Liq(self):
        liq_edr = f'Pure_QMMM_md3.edr'
        liq_data = pyedr.edr_to_dict(liq_edr)
        equil_time_ps = 100
        liq_mask = liq_data['Time'] >= equil_time_ps        
        return np.mean(liq_data['Potential'][liq_mask])
############################################################################################################
    def process_dipole(self, scee_df, mu_Vacuum, PCM1, PCM2, settings):
        delta_mu = scee_df['muL_SCEE'] * (PCM1 / PCM2) - mu_Vacuum
        mu_liquid = delta_mu + settings.molecule.gas_dipole
        
        df = scee_df[['config', 'dipole_l', 'dipole_m', 'dipole_h', 'muL_SCEE']].copy()
        df['delta_mu'] = delta_mu
        df['mu_liquid'] = mu_liquid
    
        df.to_csv('Dipole_Calculations.csv', index=False)
        return df    
############################################################################################################
    def read_data(self):
        csvfile_list = glob.glob(f'Simulations/replica_*/*/*/Dipole_Calculations.csv')
        df = pd.DataFrame()
        for i, csvfile in enumerate(csvfile_list, start=1):
            tmp = pd.read_csv(csvfile)
            df = pd.concat([df, tmp])
        return df
######################################################
    def _iqr_masks(self,x, k):
        x = np.asarray(x, dtype=float)
        q1, q3 = np.percentile(x, [25, 75])
        iqr = q3 - q1
        lo = q1 - k * iqr
        hi = q3 + k * iqr
        out = (x < lo) | (x > hi)
        keep = ~out
        return keep, out, (q1, q3, iqr, lo, hi)
######################################################
    def _normal_pdf(self,x, mu, sigma):
        return (1.0 / (sigma * np.sqrt(2.0 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
######################################################
    def plot_state_point(self, df_pool, T, p, mu_Vacuum, k=2.0, bin_width=0.25, x_ceil=10.0):
        x_all = df_pool['mu_liquid'].astype(float).to_numpy()
        finite_mask = np.isfinite(x_all)
 
        if finite_mask.sum() < 5:
            print(f"Warning: <5 finite values at T={T}K, p={p}Bar — skipping")
            return finite_mask 
 
        x_finite = x_all[finite_mask]
        keep_iqr, out_iqr, (q1, q3, iqr, lo, hi) = self._iqr_masks(x_finite, k=k)
        out_phys = x_finite < mu_Vacuum
        out_combined = out_iqr | out_phys
 
        stats_p1 = self._plot_distribution_pass(
            x_finite, out_combined,
            T=T, p=p, pass_label='raw',
            bin_width=bin_width, x_ceil=x_ceil,
        )
        stats_p1.update({
            'T': T, 'p': p, 'pass': 'raw',
            'N_outliers_phys': int(out_phys.sum()),
            'N_outliers_stat': int(out_iqr.sum()),
            'mu_Vacuum_threshold': mu_Vacuum,
        })
 
        keep_p1_inner = ~out_combined
        x_p1_kept = x_finite[keep_p1_inner]
 
        if x_p1_kept.size < 5:
            print(f"Warning: <5 values after pass 1 at T={T}K, p={p}Bar — "
                  f"skipping pass 2")

            keep_full = np.zeros(len(x_all), dtype=bool)
            keep_full[finite_mask] = keep_p1_inner
            self._append_stats_rows([stats_p1])
            return keep_full
 
        keep_p2_inner, out_p2_inner, _ = self._iqr_masks(x_p1_kept, k=k)
 
        stats_p2 = self._plot_distribution_pass(
            x_p1_kept, out_p2_inner,
            T=T, p=p, pass_label='filtered',
            bin_width=bin_width, x_ceil=x_ceil,
        )
        stats_p2.update({
            'T': T, 'p': p, 'pass': 'filtered',
            'N_outliers_phys': 0,
            'N_outliers_stat': int(out_p2_inner.sum()),
            'mu_Vacuum_threshold': mu_Vacuum,
        })
 
        self._append_stats_rows([stats_p1, stats_p2])

        keep_full = np.zeros(len(x_all), dtype=bool)

        finite_indices = np.where(finite_mask)[0]
        p1_kept_indices = finite_indices[keep_p1_inner]

        p2_kept_indices = p1_kept_indices[keep_p2_inner]
        keep_full[p2_kept_indices] = True
 
        return keep_full
######################################################
    def _plot_distribution_pass(self, x_all, outlier_mask, T, p, pass_label, bin_width=0.25, x_ceil=10.0):
        x_keep = x_all[~outlier_mask]
        x_out  = x_all[outlier_mask]
 
        if x_keep.size < 2:
            print(f"Warning: <2 kept values at T={T}K, p={p}Bar pass={pass_label} "
                  f"— cannot fit Gaussian")
            return {
                'N_total':    int(x_all.size),
                'N_kept':     int(x_keep.size),
                'N_outliers': int(x_out.size),
                'mean':       float(x_keep.mean()) if x_keep.size else np.nan,
                'std':        np.nan,
                'median':     float(np.median(x_keep)) if x_keep.size else np.nan,
            }
 
        mu    = x_keep.mean()
        sigma = x_keep.std(ddof=1) 

        q1, q3 = np.percentile(x_keep, [25, 75])
 
        xmin = 0.0
        x_within_ceil = x_all[x_all <= x_ceil]
        if x_within_ceil.size > 0:
            xmax = np.ceil(x_within_ceil.max() / bin_width) * bin_width
        else:
            xmax = x_ceil
        xmax = min(xmax, x_ceil)
 
        bin_edges = np.arange(xmin, xmax + bin_width, bin_width)
        bw = bin_edges[1] - bin_edges[0]
 

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.hist(x_keep, bins=bin_edges, alpha=0.9, color='purple', label="Data")
        if x_out.size:
            x_out_visible = x_out[(x_out >= xmin) & (x_out <= xmax)]
            if x_out_visible.size:
                ax.hist(x_out_visible, bins=bin_edges, alpha=0.9,
                        color='orange', label="Outliers")
 

        x_line = np.linspace(xmin - bin_width, xmax + bin_width, 400)
        y_line = self._normal_pdf(x_line, mu, sigma) * x_keep.size * bw
        ax.plot(x_line, y_line, color='#008080', linewidth=2,
                label="Gaussian fit")
 
        ax.set_xlim(xmin, xmax)
        raw_spacing = (xmax - xmin) / 8
        tick_spacing = max(round(raw_spacing * 2) / 2, 0.5)
        ax.set_xticks(np.arange(xmin, xmax + tick_spacing, tick_spacing))
 
        ax.set_xlabel(r'$\mu_{liquid}$ / D', fontsize=16)
        ax.set_ylabel('Frequency', fontsize=16)
        ax.set_title(f'T = {T} K, p = {p} Bar', fontsize=14)
        ax.legend(frameon=False, fontsize=12)
        ax.tick_params(axis='both', which='major', labelsize=12)
 
        plt.tight_layout()
        fname = f'dipole_distribution_{pass_label}_T{T}K_p{p}Bar'
        plt.savefig(f'{fname}.png', bbox_inches="tight", dpi=300)
        plt.savefig(f'{fname}.pdf', bbox_inches="tight", dpi=300)
        plt.close(fig)
 
        return {
            'N_total':    int(x_all.size),
            'N_kept':     int(x_keep.size),
            'N_outliers': int(x_out.size),
            'mean':       round(float(mu), 6),
            'std':        round(float(sigma), 6),
            'median':     round(float(np.median(x_keep)), 6),
            'q1':         round(float(q1), 6),
            'q3':         round(float(q3), 6),
        }
######################################################
    def _append_stats_rows(self, rows):
    
        path = 'dipole_stats.csv'
        df_new = pd.DataFrame(rows)
        file_exists = os.path.exists(path)
        df_new.to_csv(
            path,
            mode='a' if file_exists else 'w',
            header=not file_exists,
            index=False,
        )
######################################################
    def write_filtered_per_leaf(self, dipole_results, keep_masks):
        cols_to_blank = ['muL_SCEE', 'delta_mu', 'mu_liquid',
                         'dipole_l', 'dipole_m', 'dipole_h']    
        
        for (T, p), keep_mask in keep_masks.items():
            df_pool = dipole_results[
                (dipole_results['T'] == T) & (dipole_results['p'] == p)
            ].reset_index(drop=True)
 
            if len(df_pool) != len(keep_mask):
                raise ValueError(
                    f"keep_mask length {len(keep_mask)} doesn't match "
                    f"pool length {len(df_pool)} at T={T}K, p={p}Bar"
                )
 

            for replica in df_pool['replica'].unique():
                replica_mask = (df_pool['replica'] == replica).to_numpy()
                replica_reject = ~keep_mask & replica_mask
 
                df_leaf = df_pool[replica_mask].copy().reset_index(drop=True)
                leaf_rejected = replica_reject[replica_mask]
 
                cols_present = [c for c in cols_to_blank
                               if c in df_leaf.columns]
                df_leaf.loc[leaf_rejected, cols_present] = np.nan
 
                leaf_dir = f'Simulations/replica_{replica}/{T}K/{p}Bar'
                out_path = os.path.join(leaf_dir,
                                        'Dipole_Calculations_filtered.csv')
                df_leaf.to_csv(out_path, index=False)
######################################################
    def _se_replica(self, x, ddof=1):
    
        x = np.asarray(x, dtype=float)
        x = x[np.isfinite(x)]
        
        if x.size < 2:
            return np.nan
        
        return x.std(ddof=ddof) / np.sqrt(x.size)
######################################################
    def build_collection(self, dipole_results,leaf_summary,keep_masks,settings, mu_Vacuum, PCM1, PCM2, dipole_model,gas_potential):

        n_exp_sq = settings.molecule.refractive_index ** 2
        qr1 = settings.advanced.charge_scaling.qr1
        qr2 = settings.advanced.charge_scaling.qr2
        qr3 = settings.advanced.charge_scaling.qr3
        
        per_replica_rows = []
        system_rows = []
        
        for T, p in itertools.product(settings.state.temperatures, 
                                       settings.state.pressures):
            df_pool = dipole_results[
                (dipole_results['T'] == T) & (dipole_results['p'] == p)
            ].reset_index(drop=True)
            
            keep_mask = keep_masks.get((T, p))
            if keep_mask is None:
                print(f"Warning: no keep_mask for T={T}K p={p}Bar, skipping")
                continue
            
            df_filtered = df_pool[keep_mask].copy()
            n_kept_total = len(df_filtered)
            n_configs_total = len(df_pool)
            
            leaf_subset = leaf_summary[
                (leaf_summary['T'] == T) & (leaf_summary['p'] == p)
            ]
            
            per_replica_accumulator = {
                'muL_SCEE_mean': [],
                'delta_mu_mean': [],
                'mu_liquid_mean': [],
                'SE_muL_SCEE': [],
                'SE_mu_liquid': [],
                'SE_dielectric': [],
                'dielectric_correction': [],
                'Liquid_Dipole_Moment_Model': [],
                'Model_Epsilon': [],
                'Simulation_Density': [],
                'Enthalpy_of_Vaporisation': [],
            }
            
            for replica in df_filtered['replica'].unique():
                replica_group = df_filtered[df_filtered['replica'] == replica]
                n_kept_replica = len(replica_group)
                
                leaf_row = leaf_subset[leaf_subset['replica'] == replica]
                if len(leaf_row) != 1:
                    print(f"Warning: replica {replica} at T={T}K p={p}Bar has "
                          f"{len(leaf_row)} leaf_summary rows, expected 1; skipping")
                    continue
                leaf_row = leaf_row.iloc[0]
                
                mean_muL_SCEE_r = replica_group['muL_SCEE'].mean()
                mean_delta_mu_r = replica_group['delta_mu'].mean()
                mean_mu_liquid_r = replica_group['mu_liquid'].mean()
                
                SE_muL_SCEE_r = self._se_replica(replica_group['muL_SCEE'].to_numpy())
                SE_mu_liquid_r = self._se_replica(replica_group['mu_liquid'].to_numpy())
                
                liq_dipole_model_r = leaf_row['Liquid Dipole Moment Model']
                model_eps_r = leaf_row['Model Epsilon']
                sim_density_r = leaf_row['Simulation Density']
                liq_potential_r = leaf_row['Liquid Potential']
                enth_of_vap_r = leaf_row['Enthalpy of Vaporisation']
                n_configs_total_r = len(df_pool[df_pool['replica'] == replica])
                
                if mean_mu_liquid_r > liq_dipole_model_r:
                    diel_corr_r = (
                        n_exp_sq
                        + (mean_mu_liquid_r / liq_dipole_model_r) 
                        * (model_eps_r + 1)
                    )
                    
                    SE_diel_r = ((model_eps_r + 1) / liq_dipole_model_r) * SE_mu_liquid_r
                else:
                    diel_corr_r = np.nan
                    SE_diel_r = np.nan
                
                num1_r = qr1 * liq_dipole_model_r
                num2_r = qr2 * liq_dipole_model_r
                num3_r = qr3 * liq_dipole_model_r
                
                per_replica_rows.append({
                    'System': settings.molecule.system_title,
                    'replica': replica,
                    'T': T,
                    'p': p,
                    'N_configs_total': n_configs_total_r,
                    'N_configs_kept': n_kept_replica,
                    'mean_muL_SCEE': mean_muL_SCEE_r,
                    'SE_muL_SCEE': SE_muL_SCEE_r,
                    'mean_delta_mu': mean_delta_mu_r,
                    'mean_mu_liquid': mean_mu_liquid_r,
                    'SE_mu_liquid': SE_mu_liquid_r,
                    'Simulation Density': sim_density_r,
                    'Liquid Dipole Moment Model': liq_dipole_model_r,
                    'Model Epsilon': model_eps_r,
                    'Scaled Charge 1': num1_r,
                    'Scaled Charge 2': num2_r,
                    'Scaled Charge 3': num3_r,
                    'Liquid Potential': liq_potential_r,
                    'Enthalpy of Vaporisation': enth_of_vap_r,
                    'Dielectric Constant Correction': diel_corr_r,
                    'SE Dielectric Constant Correction': SE_diel_r,
                })
                
                per_replica_accumulator['muL_SCEE_mean'].append(mean_muL_SCEE_r)
                per_replica_accumulator['delta_mu_mean'].append(mean_delta_mu_r)
                per_replica_accumulator['mu_liquid_mean'].append(mean_mu_liquid_r)
                per_replica_accumulator['SE_muL_SCEE'].append(SE_muL_SCEE_r)
                per_replica_accumulator['SE_mu_liquid'].append(SE_mu_liquid_r)
                per_replica_accumulator['SE_dielectric'].append(SE_diel_r)
                per_replica_accumulator['dielectric_correction'].append(diel_corr_r)
                per_replica_accumulator['Liquid_Dipole_Moment_Model'].append(liq_dipole_model_r)
                per_replica_accumulator['Model_Epsilon'].append(model_eps_r)
                per_replica_accumulator['Simulation_Density'].append(sim_density_r)
                per_replica_accumulator['Enthalpy_of_Vaporisation'].append(enth_of_vap_r)
            
            n_replicas_actual = len(per_replica_accumulator['muL_SCEE_mean'])
            if n_replicas_actual == 0:
                continue
            
            system_mean_muL_SCEE = np.mean(per_replica_accumulator['muL_SCEE_mean'])
            system_mean_delta_mu = np.mean(per_replica_accumulator['delta_mu_mean'])
            system_mean_mu_liquid = np.mean(per_replica_accumulator['mu_liquid_mean'])
            system_mean_liq_dipole = np.mean(per_replica_accumulator['Liquid_Dipole_Moment_Model'])
            system_mean_model_eps = np.mean(per_replica_accumulator['Model_Epsilon'])
            system_mean_density = np.mean(per_replica_accumulator['Simulation_Density'])
            system_mean_enth = np.mean(per_replica_accumulator['Enthalpy_of_Vaporisation'])
            
            system_num1 = qr1 * system_mean_liq_dipole
            system_num2 = qr2 * system_mean_liq_dipole
            system_num3 = qr3 * system_mean_liq_dipole
            
            def rms_combine(se_list):
                arr = np.asarray(se_list, dtype=float)
                arr = arr[np.isfinite(arr)]
                if arr.size == 0:
                    return np.nan
                return np.sqrt(np.sum(arr**2)) / len(arr)
            
            system_SE_muL_SCEE = rms_combine(per_replica_accumulator['SE_muL_SCEE'])
            system_SE_mu_liquid = rms_combine(per_replica_accumulator['SE_mu_liquid'])
            system_SE_dielectric = rms_combine(per_replica_accumulator['SE_dielectric'])
            
            diel_corrs = np.asarray(per_replica_accumulator['dielectric_correction'])
            n_diel_valid = np.sum(np.isfinite(diel_corrs))
            system_diel_corr = np.nanmean(diel_corrs) if n_diel_valid > 0 else np.nan
            
            system_rows.append({
                'System': settings.molecule.system_title,
                'T': T,
                'p': p,
                'N_replicas': n_replicas_actual,
                'Number of Configurations': n_configs_total,
                'Number of Configurations Filtered': n_kept_total,
                
                'Experimental Density': settings.molecule.density,
                'Experimental Gas Dipole Moment': settings.molecule.gas_dipole,
                'Experimental Refractive Index': settings.molecule.refractive_index,
                'Experimental Epsilon': settings.molecule.dielectric_constant,

                'Scaling value 1': qr1,
                'Scaling value 2': qr2,
                'Scaling value 3': qr3,

                'Scaled Charge 1': system_num1,
                'Scaled Charge 2': system_num2,
                'Scaled Charge 3': system_num3,

                'Vacuum Dipole Moment Model': dipole_model,
                'mu_Vacuum': mu_Vacuum,
                'PCM1': PCM1,
                'PCM2': PCM2,
                'cal_diconst': settings.molecule.calculated_dielectric,

                'Simulation Density': system_mean_density,
                'Liquid Dipole Moment Model': system_mean_liq_dipole,
                'Model Epsilon': system_mean_model_eps,
                'Enthalpy of Vaporisation': system_mean_enth,

                'mean_muL_SCEE': system_mean_muL_SCEE,
                'SE_muL_SCEE': system_SE_muL_SCEE,
                'mean_delta_mu': system_mean_delta_mu,
                'mean_mu_liquid': system_mean_mu_liquid,
                'SE_mu_liquid': system_SE_mu_liquid,

                'Dielectric Constant Correction': system_diel_corr,
                'SE Dielectric Constant Correction': system_SE_dielectric,
            })
        
        per_replica_df = pd.DataFrame(per_replica_rows)
        system_df = pd.DataFrame(system_rows)
        return per_replica_df, system_df
######################################################



