#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Ajuste.py
#  Versión automatizada para procesar múltiples casos promediados.
#
"""
Este script realiza el ajuste de la serie de Fourier a los datos promediados.
Extrae automáticamente los parámetros del sistema (L, Vol, N) de los archivos
log.lammps originales dentro de las carpetas de las réplicas.
"""

import numpy as np
# import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import re
import os
import glob
# --- CRUCIAL: Configurar backend no interactivo ANTES de importar pyplot ---
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ==========================================
# 1. CONFIGURACIÓN Y CONSTANTES
# ==========================================
N_TERMS = 20      # Términos en la serie de Fourier
X_OFFSET = 20000  # Offset de tiempo (acorde al script original)
KB = 1.38e-23     # Constante de Boltzmann (J/K)
DIM = 3           # Dimensiones

# ==========================================
# 2. FUNCIONES DE EXTRACCIÓN Y MATEMÁTICAS
# ==========================================

def extraer_parametros_de_replica(case_id):
    """
    Busca el log.lammps en la carpeta de la réplica 1 para obtener metadatos.
    """
    # El caso es algo como S444_DT5_TS0.5. Buscamos run_S444_DT5_TS0.5_R1
    #replica_path = f"run_{case_id}_R1/lammps.log"
    replica_path = f"run_{case_id}_R1/log.lammps"

    data = {
        "volumen_m3": None,
        "longitud_x_m": None,
        "num_atomos": None,
        "t_hot": None,
        "t_cold": None
    }

    if not os.path.exists(replica_path):
        print(f"  [!] Advertencia: No se encontró {replica_path}. Usando valores por defecto.")
        return None

    try:
        with open(replica_path, 'r') as f:
            content = f.read()

        # Extraer dimensiones (Angstroms a Metros)
        box_match = re.findall(r"triclinic box = \(0 0 0\) to \(([\d\.]+) ([\d\.]+) ([\d\.]+)\)", content)
        if box_match:
            lx, ly, lz = map(float, box_match[-1])
            data["longitud_x_m"] = lx * 1e-10
            data["volumen_m3"] = (lx * ly * lz) * 1e-30

        # Extraer número de átomos
        atoms_match = re.findall(r"(\d+) atoms", content)
        if atoms_match:
            data["num_atomos"] = int(atoms_match[-1])

        # Extraer temperaturas
        th = re.search(r"variable T_hot_pulse equal ([\d\.]+)", content)
        tc = re.search(r"variable T_cold_pulse equal ([\d\.]+)", content)
        if th: data["t_hot"] = float(th.group(1))
        if tc: data["t_cold"] = float(tc.group(1))

        return data
    except Exception as e:
        print(f"  [!] Error procesando log: {e}")
        return None

def f_total(x, k2, L, T_max, T_min):
    """Sumatoria de la serie de Fourier para el ajuste."""
    total = 0
    for n in range(1, N_TERMS + 1):
        k_n = (np.pi * n) / (L * 1e10) # L en Angstroms para la serie
        k1_n = np.pi * n
        b_n = (1.0 / k1_n) * (-np.cos(k1_n) * (T_min - T_max) - T_max + T_min)
        term_c = (-2 * b_n + 2 * b_n * np.cos(k1_n)) * (1.0 / k1_n)
        total += np.exp(-k2 * x * (k_n**2)) * term_c
    return total
def f_total2(x, k, A, T_offset):
    """
    Función de ajuste mejorada.
    He añadido T_offset para que el algoritmo pueda ajustar el nivel base
    y no solo la caída exponencial.
    """
    return A * np.exp(-k * x) + T_offset

def realizar_ajuste_profesional(x_data, y_data, case_id):
    """
    Realiza el ajuste con estimación automática de valores iniciales.
    """
    print(f" Procesando: {case_id}")

    # 1. Estimación de valores iniciales (Crucial para evitar resultados idénticos)
    # k inicial: estimado de la pendiente logarítmica simple
    # A inicial: diferencia entre max y min
    # T_offset inicial: el valor mínimo observado
    t_range = np.max(y_data) - np.min(y_data)
    p0 = [1e-3, t_range, np.min(y_data)]

    # Definir límites (opcional pero recomendado)
    # k debe ser positivo, A puede ser positivo, offset cerca del ambiente
    bounds = (0, [1.0, 1000.0, 500.0])

    try:
        # 2. Ejecutar el ajuste
        popt, pcov = curve_fit(f_total, x_data, y_data, p0=p0, bounds=bounds, maxfev=5000)
        k_fit, A_fit, offset_fit = popt

        # 3. Generar gráfica de diagnóstico (headless)
        plt.figure(figsize=(10, 6))
        plt.scatter(x_data, y_data, s=10, color='gray', alpha=0.5, label='Datos')
        plt.plot(x_data, f_total(x_data, *popt), 'r-', linewidth=2,
                 label=f'Ajuste (k={k_fit:.2e})')
        plt.title(f"Diagnóstico de Ajuste - {case_id}")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(f"diag_{case_id}.png")
        plt.close()

        return k_fit

    except Exception as e:
        print(f"  [!] Error en {case_id}: {e}")
        return None

def procesar_lote(casos):
    # Simulamos una lista de casos
    #casos = ["S1644_DT5_TS0.5", "S444_DT25_TS1.0"]
    resultados = []

    for caso in casos:
        # Aquí cargarías tus datos reales de los archivos .log
        # x_data, y_data = cargar_datos(caso)

        # Generamos datos sintéticos para el ejemplo
        x = np.linspace(0, 100, 200)
        y = 50 * np.exp(-0.02 * x) + 300 + np.random.normal(0, 0.2, 200)

        k_res = realizar_ajuste_profesional(x, y, caso)

        # Cálculo de conductividad (fórmula de ejemplo basada en tus resultados)
        # Ajusta esta constante según tu física
        conductividad = k_res * 0.0299 if k_res else 0

        resultados.append({
            "CASO": caso,
            "K2_FIT": k_res,
            "CONDUCTIVIDAD": conductividad
        })

    df = pd.DataFrame(resultados)
    print("\n--- RESULTADOS FINALES ---")
    print(df.to_string(index=False))
# ==========================================
# 3. BUCLE PRINCIPAL DE PROCESAMIENTO
# ==========================================

def run_analysis():
    # Buscar archivos promediados: averaged_S444_DT5_TS0.5.dat
    avg_files = glob.glob("averaged_*.dat")

    if not avg_files:
        print("No se encontraron archivos 'averaged_*.dat'. Ejecuta Media_replicas.py primero.")
        return

    results_summary = []

    for file_path in avg_files:
        # Extraer ID del caso del nombre del archivo
        case_id = file_path.replace("averaged_", "").replace(".dat", "")
        print(f"\n>>> Analizando Caso: {case_id}")

        # 1. Obtener constantes del sistema
        params = extraer_parametros_de_replica(case_id)
        if not params:
            continue

        # 2. Cargar datos promediados
        try:
            # Step, Mean_DeltaT, Std_DeltaT, N_Replicas
            df = pd.read_csv(file_path, sep='\t')
            x_data = df['Step'].values - X_OFFSET
            y_data = df['Mean_DeltaT'].values
        except Exception as e:
            print(f"  [!] Error cargando datos: {e}")
            continue

        # 3. Ajuste (Fit)
        # Definimos una lambda para pasar L, Tmax, Tmin fijas al curve_fit
        # Nota: L para la función Fourier interna se usa en Angstroms como el script original
        L_ang = params['longitud_x_m'] * 1e10
        fit_func = lambda x, k2: f_total(x, k2, L_ang, params['t_hot'], params['t_cold'])

        try:
            popt, _ = curve_fit(fit_func, x_data, y_data, p0=[0.00001])
            k2_fit = popt[0]

            # 4. Cálculo de Conductividad (kbarr)
            # kbarr = (k2_fit * 1e-5 * DIM * N * KB) / VOL
            # El factor 1e-5 suele venir del timestep de LAMMPS (verificar si aplica aquí)
            kbarr = (k2_fit * 1e-5 * DIM * params['num_atomos'] * KB) / params['volumen_m3']

            print(f"  [OK] K2: {k2_fit:.4e} | Conductividad: {kbarr:.4f} W/mK")
            results_summary.append((case_id, k2_fit, kbarr))

            # 5. Guardar Gráfica
            plt.figure(figsize=(8, 5))
            plt.errorbar(x_data + X_OFFSET, y_data, yerr=df['Std_DeltaT'], fmt='o',
                         ms=2, alpha=0.3, label='Datos Promediados', color='gray')
            plt.plot(x_data + X_OFFSET, fit_func(x_data, k2_fit),
                     label=f'Ajuste (k={kbarr:.3f})', color='red', lw=2)
            plt.title(f'Ajuste AEMD - {case_id}')
            plt.xlabel('Timestep')
            plt.ylabel('Delta T')
            plt.legend()
            plt.savefig(f"plot_{case_id}.png")
            plt.close() # Limpiar memoria de la figura
            print(f"  [+] Gráfica guardada como: plot_{case_id}.png")#plt.close()

        except Exception as e:
            print(f"  [!] Falló el ajuste para {case_id}: {e}")

    # 6. Guardar resumen final
    if results_summary:
        with open("resumen_resultados.txt", "w") as f:
            f.write("CASO\tK2_FIT\tCONDUCTIVIDAD_W_MK\n")
            for res in results_summary:
                f.write(f"{res[0]}\t{res[1]:.6e}\t{res[2]:.6f}\n")
        print("\nAnálisis completado. Resultados en 'resumen_resultados.txt'.")

import pandas as pd # Importado aquí para asegurar disponibilidad
if __name__ == "__main__":
    run_analysis()
