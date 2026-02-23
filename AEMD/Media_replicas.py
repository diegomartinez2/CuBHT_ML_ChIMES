#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Media_replicas.py
#
#  Copyright 2026 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
# ---------------------------
# Importación de los módulos
# ---------------------------
"""
Este script extrae las salidas de datos obtenidas con run_AEMD_study3.sh
calcula la media de las repeticiones con los mismos DT y ST, y prepara
para el ajuste (Ajuste.py)
"""
import os
import pandas as pd
import numpy as np
import glob

def process_results():
    base_dir = "RESULTS"
    # Buscamos todas las combinaciones de parámetros (exceptuando la carpeta REP)
    # Estructura: RESULTS/SIZE_*/DT_*/TS_*/
    case_patterns = glob.glob(os.path.join(base_dir, "SIZE_*", "DT_*", "TS_*"))

    for case_path in case_patterns:
        all_data = []

        # Buscar todas las réplicas dentro de este caso
        replica_dirs = glob.glob(os.path.join(case_path, "REP_*"))

        if not replica_dirs:
            continue

        for rep_dir in replica_dirs:
            # Normalmente es "delta_T_AEMD.dat", pero en nuestro caso es "delta_T_S{}_DT{}_TS{}_R{}.dat"
            file_path = os.path.join(rep_dir, "delta_T_AEMD.dat")
            if os.path.exists(file_path):
                # LAMMPS ave/time tiene 2 líneas de comentario iniciales
                # El formato suele ser: TimeStep Number-of-bins v_delta_T
                # pero en este caso es de TimeStep v_delta_T
                try:
                    #df = pd.read_csv(file_path, sep='\s+', skiprows=2, names=['Step', 'Count', 'DeltaT'])
                    df = pd.read_csv(file_path, sep='\s+', skiprows=2, names=['Step', 'DeltaT'])
                    all_data.append(df.set_index('Step')['DeltaT']) # TimeStep v_delta_T
                except Exception as e:
                    print(f"Error leyendo {file_path}: {e}")

        if all_data:
            # Concatenar réplicas por el índice (Step)
            combined = pd.concat(all_data, axis=1)

            # Calcular media y desviación estándar
            mean_series = combined.mean(axis=1)
            std_series = combined.std(axis=1)

            # Crear DataFrame de salida
            final_df = pd.DataFrame({
                'Mean_DeltaT': mean_series,
                'Std_DeltaT': std_series,
                'N_Replicas': combined.count(axis=1)
            }).reset_index()

            # Guardar en la carpeta del caso (fuera de REP_X)
            #if os.path.exists(output_name):
                #print(f"Saltando: {output_name} ya existe.")
                #continue
            output_name = os.path.join(case_path, "averaged_delta_T.dat")
            final_df.to_csv(output_name, sep='\t', index=False)
            print(f"Generado: {output_name}")

if __name__ == "__main__":
    process_results()
