#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Media_replicas_AEMD_v2.py
Promedia resultados de múltiples réplicas de LAMMPS.
Entrada esperada: TimeStep, c_th, c_tc, v_delta_T (skiprows=2)
Salida solicitada: step, time_ps, Th, Tc, DeltaT, Std_DeltaT
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

    # 2. Agrupar carpetas por casos (ej: run_S100_DT50_TS0.5_R1, R2...)
    cases = {}
    for run_dir in all_runs:
        # Busca el patrón del nombre del caso ignorando el sufijo de la réplica _R*
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

        # Extraer el timestep del nombre del directorio (ej: TS0.5 -> 0.5)
        # Esto es útil para calcular el tiempo en ps si no viene en el archivo
        ts_match = re.search(r"_TS([\d.]+)", case_id)
        timestep_val = float(ts_match.group(1)) if ts_match else 0.001

        print(f"Procesando {case_id} ({len(replica_dirs)} réplicas)...")

        for rep_dir in replica_dirs:
            folder_basename = os.path.basename(rep_dir)
            file_basename = folder_basename.replace("run_", "delta_T_") + ".dat"
            file_path = os.path.join(rep_dir, file_basename)

            if os.path.exists(file_path):
                try:
                    # MODIFICACIÓN: Ajuste de nombres de entrada según tu solicitud
                    # Suponiendo que el archivo tiene: TimeStep, c_th, c_tc, v_delta_T
                    df = pd.read_csv(file_path, sep='\s+', skiprows=2,
                                     names=['TimeStep', 'c_th', 'c_tc', 'v_delta_T'])

                    if not df.empty:
                        # Usamos TimeStep como índice para alinear réplicas
                        data_dT.append(df.set_index('TimeStep')['v_delta_T'])
                        data_Th.append(df.set_index('TimeStep')['c_th'])
                        data_Tc.append(df.set_index('TimeStep')['c_tc'])

                except Exception as e:
                    print(f"  [!] Error leyendo {file_path}: {e}")
            else:
                print(f"  [?] Archivo no encontrado: {file_path}")

        # 4. Calcular estadísticas
        if data_dT:
            # Combinar series de réplicas
            df_dT = pd.concat(data_dT, axis=1)
            df_Th = pd.concat(data_Th, axis=1)
            df_Tc = pd.concat(data_Tc, axis=1)

            # Promedios y desviación estándar
            mean_dT = df_dT.mean(axis=1)
            std_dT = df_dT.std(axis=1).fillna(0) # fillna(0) si solo hay 1 réplica
            mean_Th = df_Th.mean(axis=1)
            mean_Tc = df_Tc.mean(axis=1)

            # Construir DataFrame final con los nombres de salida solicitados
            final_df = pd.DataFrame({
                'step': mean_dT.index,
                'Th': mean_Th.values,
                'Tc': mean_Tc.values,
                'DeltaT': mean_dT.values,
                'Std_DeltaT': std_dT.values
            })

            # Calcular time_ps: step * timestep / 1000 (convertir fs a ps si el TS está en fs)
            # Si tu TS ya está en ps, quita el / 1000
            final_df['time_ps'] = final_df['step'] * timestep_val / 1000.0

            # Reordenar columnas exactamente como pediste:
            # step, time_ps, Th, Tc, DeltaT, Std_DeltaT
            final_df = final_df[['step', 'time_ps', 'Th', 'Tc', 'DeltaT', 'Std_DeltaT']]

            # 5. Guardar el resultado
            output_name = f"averaged_{case_id.replace('run_', '')}.dat"
            final_df.to_csv(output_name, sep='\t', index=False)
            print(f"  [OK] Generado: {output_name}")
        else:
            print(f"  [!] No se pudieron recolectar datos para {case_id}")

if __name__ == "__main__":
    process_results()
