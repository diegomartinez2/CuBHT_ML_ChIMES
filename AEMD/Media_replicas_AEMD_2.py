#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Media_replicas.py
Modificado para promediar Th, Tc y DeltaT de múltiples réplicas de LAMMPS.
Salida: Step, Time_ps, Th, Tc, DeltaT, Std_DeltaT
"""

import os
import pandas as pd
import numpy as np
import glob
import re

def process_results():
    # 1. Buscar todas las carpetas de ejecución
    all_runs = glob.glob("run_S*")

    if not all_runs:
        print("No se encontraron carpetas que coincidan con 'run_S*'.")
        return

    # 2. Agrupar carpetas por casos (ignorando la réplica _R*)
    cases = {}
    for run_dir in all_runs:
        match = re.match(r"(run_S.+_DT.+_TS.+)_R\d+", run_dir)
        if match:
            case_name = match.group(1)
            if case_name not in cases:
                cases[case_name] = []
            cases[case_name].append(run_dir)

    print(f"Detectados {len(cases)} casos únicos para promediar.")

    # 3. Procesar cada caso
    for case_id, replica_dirs in cases.items():
        data_dT = []
        data_Th = []
        data_Tc = []
        time_mapping = {}

        print(f"Procesando {case_id} ({len(replica_dirs)} réplicas)...")

        for rep_dir in replica_dirs:
            folder_basename = os.path.basename(rep_dir)
            file_basename = folder_basename.replace("run_", "delta_T_") + ".dat"
            file_path = os.path.join(rep_dir, file_basename)

            if os.path.exists(file_path):
                try:
                    # El fix de LAMMPS: 1.Step 2.Count 3.Time_ps 4.Th 5.Tc 6.DeltaT
                    df = pd.read_csv(file_path, sep='\s+', skiprows=2,
                                     names=['Step', 'Count', 'Time_ps', 'Th', 'Tc', 'DeltaT'])

                    if not df.empty:
                        # Guardamos las series usando Step como índice
                        data_dT.append(df.set_index('Step')['DeltaT'])
                        data_Th.append(df.set_index('Step')['Th'])
                        data_Tc.append(df.set_index('Step')['Tc'])

                        # Guardamos la relación Step -> Time_ps de la primera réplica válida
                        if not time_mapping:
                            time_mapping = df.set_index('Step')['Time_ps'].to_dict()

                except Exception as e:
                    print(f"  [!] Error leyendo {file_path}: {e}")
            else:
                print(f"  [?] Archivo no encontrado: {file_path}")

        # 4. Calcular estadísticas
        if data_dT:
            # Combinamos y promediamos
            df_dT = pd.concat(data_dT, axis=1)
            df_Th = pd.concat(data_Th, axis=1)
            df_Tc = pd.concat(data_Tc, axis=1)

            mean_dT = df_dT.mean(axis=1)
            std_dT = df_dT.std(axis=1)
            mean_Th = df_Th.mean(axis=1)
            mean_Tc = df_Tc.mean(axis=1)

            final_df = pd.DataFrame({
                'Step': mean_dT.index,
                'Mean_Th': mean_Th.values,
                'Mean_Tc': mean_Tc.values,
                'Mean_DeltaT': mean_dT.values,
                'Std_DeltaT': std_dT.values,
                'N_Replicas': df_dT.count(axis=1).values
            })

            # Añadir la columna de tiempo en ps mapeando los Steps
            final_df['Time_ps'] = final_df['Step'].map(time_mapping)

            # Reordenar columnas según el formato solicitado y compatible con fit_aemd.py
            # Salida: Step, Time_ps, Th, Tc, DeltaT, Std_DeltaT
            final_df = final_df[['Step', 'Time_ps', 'Mean_Th', 'Mean_Tc', 'Mean_DeltaT', 'Std_DeltaT']]

            # 5. Guardar el resultado
            output_name = f"averaged_{case_id.replace('run_', '')}.dat"

            # Guardamos con tabulaciones.
            # Nota: fit_aemd.py leerá Th, Tc y dT promediados correctamente.
            final_df.to_csv(output_name, sep='\t', index=False)
            print(f"  [OK] Generado: {output_name}")
        else:
            print(f"  [!] No se pudieron recolectar datos para {case_id}")

if __name__ == "__main__":
    process_results()
