#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Media_replicas.py
#  Modified for flat directory structure (run_S*_DT*_TS*_R*)
#
"""
Este script extrae las salidas de datos obtenidas con el script de automatización Bash,
agrupa las réplicas que tienen el mismo SIZE, DT y TS, calcula la media y la
desviación estándar, y guarda un archivo promediado por cada caso.
"""

import os
import pandas as pd
import numpy as np
import glob
import re

def process_results():
    # 1. Buscar todas las carpetas de ejecución
    # El patrón busca carpetas que empiecen por "run_S"
    all_runs = glob.glob("run_S*")

    if not all_runs:
        print("No se encontraron carpetas que coincidan con 'run_S*'.")
        return

    # 2. Agrupar carpetas por casos (ignorando la réplica _R*)
    # Usamos un diccionario donde la clave es el "caso" y el valor es la lista de carpetas
    cases = {}
    for run_dir in all_runs:
        # Extraer el prefijo del caso usando regex: run_S{SIZE}_DT{DT}_TS{TS}
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
        print(f"Procesando {case_id} ({len(replica_dirs)} réplicas)...")

        for rep_dir in replica_dirs:
            # Extraer info de la carpeta para reconstruir el nombre del archivo .dat
            # Estructura: run_S444_DT5_TS0.5_R1 -> delta_T_S444_DT5_TS0.5_R1.dat
            folder_basename = os.path.basename(rep_dir)
            file_basename = folder_basename.replace("run_", "delta_T_") + ".dat"
            file_path = os.path.join(rep_dir, file_basename)

            if os.path.exists(file_path):
                try:
                    # Lectura del archivo de LAMMPS (skiprows=2 para saltar comentarios)
                    # Formato esperado: TimeStep v_delta_T
                    df = pd.read_csv(file_path, sep='\s+', skiprows=2, names=['Step', 'DeltaT'])

                    if not df.empty:
                        all_data.append(df.set_index('Step')['DeltaT'])
                except Exception as e:
                    print(f"  [!] Error leyendo {file_path}: {e}")
            else:
                print(f"  [?] Archivo no encontrado: {file_path}")

        # 4. Calcular estadísticas si hay datos
        if all_data:
            # Alinear por Step (índice) por si acaso hay discrepancias de longitud
            combined = pd.concat(all_data, axis=1)

            mean_series = combined.mean(axis=1)
            std_series = combined.std(axis=1)

            final_df = pd.DataFrame({
                'Mean_DeltaT': mean_series,
                'Std_DeltaT': std_series,
                'N_Replicas': combined.count(axis=1)
            }).reset_index()

            # 5. Guardar el resultado
            # Se guarda un archivo .dat por cada caso en el directorio raíz
            output_name = f"averaged_{case_id.replace('run_', '')}.dat"
            final_df.to_csv(output_name, sep='\t', index=False)
            print(f"  [OK] Generado: {output_name}")
        else:
            print(f"  [!] No se pudieron recolectar datos para {case_id}")

if __name__ == "__main__":
    process_results()
