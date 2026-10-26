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
    # Single definition so the file Run_SCEE deletes is the one we append to
    STATS_PATH = './Results/Dipole_Stats.csv'

    def __init__(self):
        pass
################################################################################
    @staticmethod
    def md3_prefix():
        """Name stem of the production run in the current leaf folder.

        New runs write Pure_QMMM_md3.*; older runs used <Molecule>_QMMM_md3.*
        (e.g. Phenol_QMMM_md3). Use whichever exists here.
        """
        if os.path.exists('Pure_QMMM_md3.tpr'):
            return 'Pure_QMMM_md3'
        found = sorted(f[:-len('.tpr')] for f in glob.glob('*_QMMM_md3.tpr')
                       if not f.startswith('Mixture_'))
        if not found:
            raise FileNotFoundError(
                f"No *_QMMM_md3.tpr production run found in {os.getcwd()}")
        if len(found) > 1:
            print(f'WARNING: several production runs found {found}; using {found[0]}')
        return found[0]
################################################################################
    @staticmethod
    def leaf_path(replica, T, p):
        """Leaf folder for a replica/state point.

        New runs use e.g. 330.15K/1Bar; older runs used 330.15K/1.0Bar.
        Prefer a folder that actually holds a production run, so an empty
        folder left behind by a crashed run can't hide the real data.
        Falls back to the new-style name if neither exists yet.
        """
        new = f'Simulations/replica_{replica}/{T:g}K/{p:g}Bar'
        old = f'Simulations/replica_{replica}/{float(T)}K/{float(p)}Bar'
        for cand in (new, old):
            if glob.glob(os.path.join(cand, '*_QMMM_md3.tpr')):
                return cand
        for cand in (new, old):
            if os.path.isdir(cand):
                return cand
        return new
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

        # Search for "Average  =   2.3301  Std. Dev. =   0.1303  Error =   0.0003"
        dipole_model = None
        for line in dipole_output.splitlines():
            stripped = line.strip()
            if stripped.startswith('Average') and stripped.split()[1] == '=':
                tokens = stripped.split()
                return float(tokens[2])
        if dipole_model is None:
            raise RuntimeError("Could not find 'Average' dipole line in gmx dipoles output")
        return dipole_model
######################################################
    def get_dipole_model_liquid(self): 
    
        gromacs.environment.flags['capture_output'] = True
        input_str = '0\n'
        run = self.md3_prefix()
        tmp = gromacs.tools.Dipoles(f=f'{run}.trr',
                                s=f'{run}.tpr',
                                input=input_str)
        chk, dipole_output, stderr = tmp.run()

        with open('Dipole2.txt', 'w') as f:
            f.write(dipole_output)
        
        Model_Dipole = None
        Epsilon = None
    
        for line in dipole_output.splitlines():
            stripped = line.strip()
            if stripped.startswith('Average') and stripped.split()[1] == '=' and Model_Dipole is None:
                Model_Dipole = float(stripped.split()[2])
            elif stripped.startswith('Epsilon'):
                Epsilon = float(stripped.split('=')[1].strip())
    
        if Model_Dipole is None:
            raise RuntimeError("Could not find 'Average  =' dipole line in gmx dipoles output")
        if Epsilon is None:
            raise RuntimeError("Could not find 'Epsilon =' line in gmx dipoles output")
            
        return Model_Dipole, Epsilon
################################################################################
    def Sim_Density(self):
        gromacs.environment.flags['capture_output'] = True
        input_str = 'Density 0\n'
        run = self.md3_prefix()
        tmp = gromacs.tools.Energy(f=f'{run}.edr',
                               s=f'{run}.tpr',
                               input=input_str)
        chk, density_output, stderr = tmp.run()

        with open('Density.txt', 'w') as f:
            f.write(density_output)
    
        # Search for the Density line (rather than fixed line number)
        for line in density_output.splitlines():
            if line.strip().startswith('Density'):
                tokens = line.split()
                # Format: "Density   997.219   0.14   8.65933   -1.03743   (kg/m^3)"
                return float(tokens[1])
    
        raise ValueError(
            "Could not find 'Density' line in gmx energy output. "
            "Check Density.txt for the actual output."
        )
        
        return density  
################################################################################
    def Sim_Enth(self,Gas_Potential, Liq_Potential, T, N):

        R = 8.314e-3 
        return Gas_Potential - (Liq_Potential / N) + R * T   
                
################################################################################
    def Pot_Gas(self):

        gas_edr = 'nvt_vacuum2.edr'
        gas_data = pyedr.edr_to_dict(gas_edr)
        equil_time_ps = 100
        gas_mask = gas_data['Time'] >= equil_time_ps
        
        return np.mean(gas_data['Potential'][gas_mask])     
################################################################################
    def Pot_Liq(self):
        liq_edr = f'{self.md3_prefix()}.edr'
        liq_data = pyedr.edr_to_dict(liq_edr)
        equil_time_ps = 100
        liq_mask = liq_data['Time'] >= equil_time_ps        
        return np.mean(liq_data['Potential'][liq_mask])
############################################################################################################
    def process_dipole(self, scee_df, mu_Vacuum, PCM1, PCM2, settings):
        SCEE=scee_df['muL_SCEE']
        delta_mu = SCEE * (PCM1 / PCM2) - mu_Vacuum
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
    def _build_pass_row(self, mu_liq, mu_scee, mu_delta, **base_fields):
        """Build a Dipole_Stats row with metrics grouped across the three quantities.
        
        Returns a dict laid out as: base fields, then mean/std/median/q1/q3 each
        with three columns (muL_SCEE, delta_mu, mu_liquid).
        """
        def _safe_stat(x, fn, *args, **kwargs):
            if x.size == 0:
                return np.nan
            if fn is np.std and x.size < 2:
                return np.nan
            return float(fn(x, *args, **kwargs))
        
        return {
            **base_fields,
            
            'mean_muL_SCEE':    _safe_stat(mu_scee,  np.mean),
            'mean_delta_mu':    _safe_stat(mu_delta, np.mean),
            'mean_mu_liquid':   _safe_stat(mu_liq,   np.mean),
            
            'std_muL_SCEE':     _safe_stat(mu_scee,  np.std, ddof=1),
            'std_delta_mu':     _safe_stat(mu_delta, np.std, ddof=1),
            'std_mu_liquid':    _safe_stat(mu_liq,   np.std, ddof=1),
            
            'median_muL_SCEE':  _safe_stat(mu_scee,  np.median),
            'median_delta_mu':  _safe_stat(mu_delta, np.median),
            'median_mu_liquid': _safe_stat(mu_liq,   np.median),
            
            'q1_muL_SCEE':      _safe_stat(mu_scee,  np.percentile, 25),
            'q1_delta_mu':      _safe_stat(mu_delta, np.percentile, 25),
            'q1_mu_liquid':     _safe_stat(mu_liq,   np.percentile, 25),
            
            'q3_muL_SCEE':      _safe_stat(mu_scee,  np.percentile, 75),
            'q3_delta_mu':      _safe_stat(mu_delta, np.percentile, 75),
            'q3_mu_liquid':     _safe_stat(mu_liq,   np.percentile, 75),
        }
######################################################
    def plot_state_point(self, df_pool, T, p, mu_Vacuum, mu_gas_exp,
                         k, max_iqr_passes, min_pass_outliers, min_kept_configs,
                         bin_width=0.25, x_ceil=10.0, system_title=None):
        # Pull all three dipole quantities from the pool
        mu_liquid_all = df_pool['mu_liquid'].astype(float).to_numpy()
        muL_SCEE_all  = df_pool['muL_SCEE'].astype(float).to_numpy()
        delta_mu_all  = df_pool['delta_mu'].astype(float).to_numpy()
        
        # Finite across all three (belt-and-braces against any NaN mismatch)
        finite_mask = (np.isfinite(mu_liquid_all) & 
                       np.isfinite(muL_SCEE_all) & 
                       np.isfinite(delta_mu_all))
        
        if finite_mask.sum() < 5:
            print(f"Warning: <5 finite values at T={T}K, p={p}Bar — skipping")
            return {
                'keep_mask': finite_mask,
                'below_floor': True,
                'n_passes_run': 0,
            }
        
        # Working arrays restricted to finite configs
        mu_liq_finite   = mu_liquid_all[finite_mask]
        mu_scee_finite  = muL_SCEE_all[finite_mask]
        mu_delta_finite = delta_mu_all[finite_mask]
        
        # Identify pass-1 outliers on mu_liquid (the filter-driving quantity).
        # IQR only: the experimental gas-phase dipole is no longer used as a
        # lower cut-off, because a one-sided cut biases the mean upwards.
        keep_iqr1, out_iqr1, _ = self._iqr_masks(mu_liq_finite, k=k)
        out_pass1 = out_iqr1
        # Diagnostic only (nothing is removed on this basis)
        n_below_gas = int((mu_liq_finite < mu_gas_exp).sum())
        
        below_floor_flag = False
        stats_rows = []
        
        # ----- Row 1: raw baseline (all finite, before any pass applied) -----
        stats_rows.append(self._build_pass_row(
            mu_liq=mu_liq_finite,
            mu_scee=mu_scee_finite,
            mu_delta=mu_delta_finite,
            N_total=int(mu_liq_finite.size),
            N_outliers=int(out_pass1.sum()),
            T=T, p=p, **{'pass': 'raw'},
            N_outliers_stat=int(out_iqr1.sum()),
            N_below_gas_exp=n_below_gas,
            mu_gas_exp=mu_gas_exp,
            mu_Vacuum_QM=mu_Vacuum,
            below_floor=False,
        ))
        
        # ----- Apply pass 1 -----
        surviving_idx = np.where(~out_pass1)[0]
        mu_liq_current   = mu_liq_finite[surviving_idx]
        mu_scee_current  = mu_scee_finite[surviving_idx]
        mu_delta_current = mu_delta_finite[surviving_idx]
        
        if mu_liq_current.size < min_kept_configs:
            print(f"WARNING: T={T}K p={p}Bar: after pass 1, {mu_liq_current.size} configs "
                  f"remain — below floor of {min_kept_configs}.")
            below_floor_flag = True
        
        # ----- Iterative IQR passes -----
        pass_num = 1
        while pass_num < max_iqr_passes and mu_liq_current.size >= 5:
            keep_next, out_next, _ = self._iqr_masks(mu_liq_current, k=k)
            n_removed = int(out_next.sum())
            
            # Record this pass's "going in" stats with the outliers IT identifies
            stats_rows.append(self._build_pass_row(
                mu_liq=mu_liq_current,
                mu_scee=mu_scee_current,
                mu_delta=mu_delta_current,
                N_total=int(mu_liq_current.size),
                N_outliers=n_removed,
                T=T, p=p, **{'pass': f'pass{pass_num}'},
                N_outliers_stat=n_removed,
                N_below_gas_exp=int((mu_liq_current < mu_gas_exp).sum()),
                mu_gas_exp=mu_gas_exp,
                mu_Vacuum_QM=mu_Vacuum,
                below_floor=below_floor_flag,
            ))
            
            # Early termination if pass found too few outliers
            if n_removed < min_pass_outliers:
                break
            
            # Floor check before applying
            if mu_liq_current.size - n_removed < min_kept_configs:
                print(f"WARNING: T={T}K p={p}Bar: pass {pass_num+1} would remove "
                      f"{n_removed} configs, breaching floor of {min_kept_configs}. "
                      f"Not applying.")
                below_floor_flag = True
                break
            
            # Apply this pass
            surviving_idx = surviving_idx[keep_next]
            mu_liq_current   = mu_liq_finite[surviving_idx]
            mu_scee_current  = mu_scee_finite[surviving_idx]
            mu_delta_current = mu_delta_finite[surviving_idx]
            pass_num += 1
        
        # ----- Final 'filtered' row: after all applied passes -----
        stats_rows.append(self._build_pass_row(
            mu_liq=mu_liq_current,
            mu_scee=mu_scee_current,
            mu_delta=mu_delta_current,
            N_total=int(mu_liq_current.size),
            N_outliers=0,
            T=T, p=p, **{'pass': 'filtered'},
            N_outliers_stat=0,
            N_below_gas_exp=int((mu_liq_current < mu_gas_exp).sum()),
            mu_gas_exp=mu_gas_exp,
            mu_Vacuum_QM=mu_Vacuum,
            below_floor=below_floor_flag,
        ))
        
        self._append_stats_rows(stats_rows)
        
        # ----- Single distribution plot: kept (purple) vs all removed (orange) -----
        all_removed_mask = np.ones(mu_liq_finite.size, dtype=bool)
        all_removed_mask[surviving_idx] = False
        mu_liq_outliers = mu_liq_finite[all_removed_mask]
        self._plot_distribution(mu_liq_current, mu_liq_outliers, T=T, p=p,
                                bin_width=bin_width, x_ceil=x_ceil,
                                system_title=system_title, mu_gas=mu_gas_exp)
        
        # ----- Build keep_full mask over the original df_pool indices -----
        finite_indices = np.where(finite_mask)[0]
        kept_global_indices = finite_indices[surviving_idx]
        keep_full = np.zeros(len(mu_liquid_all), dtype=bool)
        keep_full[kept_global_indices] = True
        
        return {
            'keep_mask': keep_full,
            'below_floor': below_floor_flag,
            'n_passes_run': pass_num,
        }

######################################################
    def _append_stats_rows(self, rows):
    
        path = self.STATS_PATH
        df_new = pd.DataFrame(rows)
        file_exists = os.path.exists(path)
        df_new.to_csv(
            path,
            mode='a' if file_exists else 'w',
            header=not file_exists,
            index=False,
        )
######################################################
    def _plot_distribution(self, x_kept, x_outliers, T, p, bin_width=0.25, x_ceil=10.0,
                           system_title=None, mu_gas=None):
        if x_kept.size < 2:
            print(f"Warning: <2 kept values at T={T}K, p={p}Bar — cannot plot")
            return
    
        mu    = x_kept.mean()
        sigma = x_kept.std(ddof=1)
    
        xmin = 0.0
        x_all_visible = np.concatenate([x_kept, x_outliers])
        x_within_ceil = x_all_visible[x_all_visible <= x_ceil]
        if x_within_ceil.size > 0:
            xmax = np.ceil(x_within_ceil.max() / bin_width) * bin_width
        else:
            xmax = x_ceil
        xmax = min(xmax, x_ceil)
    
        bin_edges = np.arange(xmin, xmax + bin_width, bin_width)
        bw = bin_edges[1] - bin_edges[0]
    
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.hist(x_kept, bins=bin_edges, alpha=0.9, color='purple', label="Data")
        if x_outliers.size:
            x_out_visible = x_outliers[(x_outliers >= xmin) & (x_outliers <= xmax)]
            if x_out_visible.size:
                ax.hist(x_out_visible, bins=bin_edges, alpha=0.9,
                        color='orange', label="Outliers")
    
        x_line = np.linspace(xmin - bin_width, xmax + bin_width, 400)
        y_line = self._normal_pdf(x_line, mu, sigma) * x_kept.size * bw
        ax.plot(x_line, y_line, color='#008080', linewidth=2, label="Gaussian fit")
        if mu_gas is not None:
            ax.axvline(mu_gas, color='black', linestyle='--', linewidth=1.5,
                       label=rf'$\mu_{{gas}}$ = {mu_gas:.2f} D')
    
        ax.set_xlim(xmin, xmax)
        raw_spacing = (xmax - xmin) / 8
        tick_spacing = max(round(raw_spacing * 2) / 2, 0.5)
        ax.set_xticks(np.arange(xmin, xmax + tick_spacing, tick_spacing))
    
        ax.set_xlabel(r'$\mu_{liquid}$ / D', fontsize=16)
        ax.set_ylabel('Frequency', fontsize=16)
        state = f'T = {T:g} K, p = {p:g} bar'
        ax.set_title(f'{system_title}: {state}' if system_title else state, fontsize=14)
        ax.legend(frameon=False, fontsize=12)
        ax.tick_params(axis='both', which='major', labelsize=12)
    
        plt.tight_layout()
        fname = f'dipole_distribution_T{T}K_p{p}Bar'
        plt.savefig(f'./Results/{fname}.png', bbox_inches="tight", dpi=300)
        plt.savefig(f'./Results/{fname}.pdf', bbox_inches="tight", dpi=300)
        plt.close(fig)
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
 
                leaf_dir = self.leaf_path(replica, T, p)
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
    def build_collection(self, dipole_results, leaf_summary,
                         keep_masks, filter_results, settings,
                         mu_Vacuum, dipole_model, PCM1, PCM2, gas_stats=None):

        n_exp_sq = settings.molecule.refractive_index ** 2
        mu_gas_exp = settings.molecule.gas_dipole
        frac = settings.advanced.optimal_dipole.polarisation_fraction
        gas_stats = gas_stats or {}

        def _gas(kind, key):
            return gas_stats.get(kind, {}).get(key, np.nan)
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
            filter_info = filter_results[(T, p)]
            below_floor = filter_info['below_floor']
            n_passes_run = filter_info['n_passes_run']
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
                'mu_optimal': [],
                'SE_mu_optimal': [],
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
                
                # Dielectric constant correction (Jorge & Lue 2019, Eq. 6):
                #   eps_corr = n^2 + (mu_L / mu_model)^2 * (eps_model - 1)
                # No special case for mu_L == mu_model: the formula is continuous,
                # and at k = 1 it gives eps_model + (n^2 - 1), not eps_model.
                k_r = mean_mu_liquid_r / liq_dipole_model_r
                diel_corr_r = n_exp_sq + k_r**2 * (model_eps_r - 1)
                # Propagated from SE(mu_L) only; eps_model and mu_model are
                # treated as exact (the replica spread is a fuller estimate).
                SE_diel_r = (2 * mean_mu_liquid_r * (model_eps_r - 1)
                             / liq_dipole_model_r**2) * SE_mu_liquid_r

                # Optimal dipole moment:
                #   mu_opt = mu_G_exp + frac * (mu_L - mu_G_SCEE)
                # with mu_L = mu_liquid and mu_G_SCEE = calculated gas-phase dipole
                # (configuration-averaged Vacuum value, mu_Vacuum).
                mu_opt_r = mu_gas_exp + frac * (mean_mu_liquid_r - mu_Vacuum)
                # Configuration-sampling uncertainty only (mu_Vacuum treated as fixed)
                SE_mu_opt_r = frac * SE_mu_liquid_r
                
                num1_r = qr1 * liq_dipole_model_r
                num2_r = qr2 * liq_dipole_model_r
                num3_r = qr3 * liq_dipole_model_r
                
                per_replica_rows.append({
                    'System': settings.molecule.system_title,
                    'replica': replica,
                    'T': T,
                    'p': p,
                    'below_floor': below_floor,
                    'N_configs_total': n_configs_total_r,
                    'N_configs_kept': n_kept_replica,
                    'mean_muL_SCEE': mean_muL_SCEE_r,
                    'SE_muL_SCEE': SE_muL_SCEE_r,
                    'mean_delta_mu': mean_delta_mu_r,
                    'mean_mu_liquid': mean_mu_liquid_r,
                    'SE_mu_liquid': SE_mu_liquid_r,
                    'mu_optimal': mu_opt_r,
                    'SE_mu_optimal': SE_mu_opt_r,
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
                per_replica_accumulator['mu_optimal'].append(mu_opt_r)
                per_replica_accumulator['SE_mu_optimal'].append(SE_mu_opt_r)
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
            system_mu_opt = np.mean(per_replica_accumulator['mu_optimal'])
            system_SE_mu_opt = rms_combine(per_replica_accumulator['SE_mu_optimal'])

            def replica_spread(vals):
                # SD of the replica means / sqrt(N): the between-replica error
                arr = np.asarray(vals, dtype=float)
                arr = arr[np.isfinite(arr)]
                return arr.std(ddof=1) / np.sqrt(arr.size) if arr.size > 1 else np.nan
            
            diel_corrs = np.asarray(per_replica_accumulator['dielectric_correction'])
            n_diel_valid = np.sum(np.isfinite(diel_corrs))
            system_diel_corr = np.nanmean(diel_corrs) if n_diel_valid > 0 else np.nan
            
            system_rows.append({
                'System': settings.molecule.system_title,
                'T': T,
                'p': p,
                'N_replicas': n_replicas_actual,
                'below_floor': below_floor,
                'N_passes_run': n_passes_run,
                'iqr_k': settings.advanced.filtering.iqr_k,
                'max_iqr_passes': settings.advanced.filtering.max_iqr_passes,
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
                'mu_Vacuum_QM': mu_Vacuum,
                'SE_mu_Vacuum_QM': _gas('Vacuum', 'SE'),
                'N_Vacuum_configs': _gas('Vacuum', 'N'),
                'PCM1': PCM1,
                'SE_PCM1': _gas('PCM1', 'SE'),
                'PCM2': PCM2,
                'SE_PCM2': _gas('PCM2', 'SE'),
                'Polarisation fraction': frac,
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
                'Replica-spread SE Dielectric Constant Correction':
                    replica_spread(per_replica_accumulator['dielectric_correction']),

                'mu_optimal': system_mu_opt,
                'SE_mu_optimal': system_SE_mu_opt,
                'Replica-spread SE mu_optimal':
                    replica_spread(per_replica_accumulator['mu_optimal']),
            })
        
        per_replica_df = pd.DataFrame(per_replica_rows)
        system_df = pd.DataFrame(system_rows)
        return per_replica_df, system_df
######################################################



