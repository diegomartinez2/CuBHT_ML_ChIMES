#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Media_replicas.py
Modificado para el nuevo formato de salida de LAMMPS:
Fix: Step, N_count, Time_ps, Th, Tc, DeltaT
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
        all_data = []
        time_mapping = {} # Para guardar la relación Step -> Time(ps)
        print(f"Procesando {case_id} ({len(replica_dirs)} réplicas)...")

        for rep_dir in replica_dirs:
            # Reconstruir el nombre del archivo generado por LAMMPS
            folder_basename = os.path.basename(rep_dir)
            # Asumimos que el archivo se llama delta_T_AEMD.dat dentro de cada carpeta
            # o sigue el patrón delta_T_S...R1.dat según tu script de Bash
            file_basename = folder_basename.replace("run_", "delta_T_") + ".dat"
            file_path = os.path.join(rep_dir, file_basename)

            if os.path.exists(file_path):
                try:
                    # El fix de LAMMPS: 1.Step 2.Count 3.Time_ps 4.Th 5.Tc 6.DeltaT
                    # Saltamos 2 líneas de cabecera (# ...)
                    df = pd.read_csv(file_path, sep='\s+', skiprows=2,
                                     names=['Step', 'Count', 'Time_ps', 'Th', 'Tc', 'DeltaT'])

                    if not df.empty:
                        # Guardamos DeltaT usando Step como índice
                        all_data.append(df.set_index('Step')['DeltaT'])

                        # Guardamos la relación Step -> Time_ps de la primera réplica válida
                        if not time_mapping:
                            time_mapping = df.set_index('Step')['Time_ps'].to_dict()

                except Exception as e:
                    print(f"  [!] Error leyendo {file_path}: {e}")
            else:
                print(f"  [?] Archivo no encontrado: {file_path}")

        # 4. Calcular estadísticas
        if all_data:
            combined = pd.concat(all_data, axis=1)

            mean_series = combined.mean(axis=1)
            std_series = combined.std(axis=1)

            final_df = pd.DataFrame({
                'Step': mean_series.index,
                'Mean_Th': mean_Th.values,
                'Mean_Tc': mean_Tc.values,
                'Mean_DeltaT': mean_series.values,
                'Std_DeltaT': std_series.values,
                'N_Replicas': combined.count(axis=1).values
            })

            # Añadir la columna de tiempo en ps mapeando los Steps
            final_df['Time_ps'] = final_df['Step'].map(time_mapping)

            # Reordenar columnas para que sea fácil de leer
            final_df = final_df[['Step', 'Time_ps',  'Mean_Th', 'Mean_Tc', 'Mean_DeltaT', 'Std_DeltaT', 'N_Replicas']]

            # 5. Guardar el resultado
            output_name = f"averaged_{case_id.replace('run_', '')}.dat"
            # Guardamos con formato compatible para que fit_aemd.py lo lea fácil
            # fit_aemd.py espera: step time Th Tc dT
            # Aquí guardamos una versión simplificada pero útil
            final_df.to_csv(output_name, sep='\t', index=False)
            print(f"  [OK] Generado: {output_name}")
        else:
            print(f"  [!] No se pudieron recolectar datos para {case_id}")

if __name__ == "__main__":
    process_results()
