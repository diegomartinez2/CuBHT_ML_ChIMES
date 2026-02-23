import numpy as np
import os
import re
import pandas as pd
from scipy.optimize import curve_fit

# --- CONFIGURACION FISICA ---
N_TERMS = 10  # Numero de terminos en la serie
L_ANG = 500.0 # Longitud en Angstroms (ajustar segun tu muestra)

def f_total(x, k2, L, T_max, T_min):
    """
    Sumatoria de la serie de Fourier original para el ajuste.
    k2 es la difusividad que queremos hallar.
    """
    total = 0
    # L se convierte internamente si es necesario, aqui usamos la logica original
    for n in range(1, N_TERMS + 1):
        k_n = (np.pi * n) / (L * 1e-10)
        k1_n = np.pi * n
        # Coeficientes segun el modelo original
        b_n = (1.0 / k1_n) * (-np.cos(k1_n) * (T_min - T_max) - T_max + T_min)
        term_c = (-2 * b_n + 2 * b_n * np.cos(k1_n)) * (1.0 / k1_n)
        total += np.exp(-k2 * x * (k_n**2)) * term_c
    return total

def extraer_datos_nombre(nombre_archivo):
    """Extrae DT y TS del nombre del archivo."""
    dt = float(re.search(r'DT(\d+\.?\d*)', nombre_archivo).group(1)) if "DT" in nombre_archivo else 1.0
    ts = float(re.search(r'TS(\d+\.?\d*)', nombre_archivo).group(1)) if "TS" in nombre_archivo else 1.0
    return dt, ts

def procesar_archivo(ruta):
    try:
        # Cargar datos (Tiempo, Temperatura)
        data = np.loadtxt(ruta, skiprows=1)
        x_data = data[:, 0]
        y_data = data[:, 1]

        # Definir temperaturas base para la serie
        t_max = np.max(y_data)
        t_min = np.min(y_data)

        # Creamos una funcion lambda que encapsula las constantes (L, Tmax, Tmin)
        # para que curve_fit solo tenga que buscar k2
        fit_func = lambda x, k2: f_total(x, k2, L_ANG, t_max, t_min)

        # Ajuste: p0 es el valor inicial de k2
        popt, _ = curve_fit(fit_func, x_data, y_data, p0=[1e-18], maxfev=10000)
        k2_fit = popt[0]

        # Calculo de conductividad final
        dt, ts = extraer_datos_nombre(ruta)
        # Aplicamos tu formula de conversion (ejemplo)
        conductividad = k2_fit * (dt / ts) * 1e15 # Ajustar factor de escala

        return {
            "CASO": os.path.basename(ruta),
            "K2_FIT": k2_fit,
            "CONDUCTIVIDAD": conductividad,
            "T_MAX": t_max,
            "T_MIN": t_min
        }
    except Exception as e:
        print(f"Error en {ruta}: {e}")
        return None

def main():
    # Lista de archivos .log en la carpeta actual
    archivos = [f for f in os.listdir('.') if f.endswith('.log')]
    resultados = []

    for archivo in archivos:
        res = procesar_archivo(archivo)
        if res:
            resultados.append(res)

    # Mostrar resultados
    if resultados:
        df = pd.DataFrame(resultados)
        print("\n" + "="*60)
        print(f"{'CASO':<25} | {'K2_FIT':<12} | {'CONDUCTIVIDAD':<12}")
        print("-" * 60)
        for _, row in df.iterrows():
            print(f"{row['CASO']:<25} | {row['K2_FIT']:<12.6e} | {row['CONDUCTIVIDAD']:<12.6f}")
        print("="*60)

if __name__ == "__main__":
    main()
