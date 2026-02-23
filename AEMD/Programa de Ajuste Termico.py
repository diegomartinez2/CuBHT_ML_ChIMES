import numpy as np
import os
import re
import pandas as pd
from scipy.optimize import curve_fit

def modelo_series_calor(t, alpha, A, offset):
    """
    Modelo simplificado de una serie de Fourier para transferencia
    de calor en una dimension (Pulso Laser / Flash Method).
    """
    # Representamos el termino principal de la serie que domina el decaimiento
    # T(t) = Offset + A * [1 + 2 * sum( (-1)^n * exp(-n^2 * pi^2 * alpha * t / L^2) )]
    # Aqui simplificamos k = pi^2 * alpha / L^2
    return A * (1 - 2 * np.exp(-alpha * t)) + offset

def extraer_parametros_nombre(nombre_archivo):
    """
    Extrae DT y TS del nombre del archivo usando expresiones regulares.
    Ejemplo: S1644_DT5_TS0.5 -> DT=5, TS=0.5
    """
    dt_match = re.search(r'DT(\d+\.?\d*)', nombre_archivo)
    ts_match = re.search(r'TS(\d+\.?\d*)', nombre_archivo)

    dt = float(dt_match.group(1)) if dt_match else 1.0
    ts = float(ts_match.group(1)) if ts_match else 1.0
    return dt, ts

def procesar_experimento(ruta_archivo):
    # 1. Cargar datos saltando cabeceras si existen
    # Se asume formato: Tiempo, Temperatura
    try:
        data = np.loadtxt(ruta_archivo, skiprows=1)
        t = data[:, 0]
        temp = data[:, 1]
    except Exception as e:
        return f"Error al leer {ruta_archivo}: {e}", None

    # 2. Parametros iniciales para el ajuste de series
    # alpha (difusividad), A (amplitud), offset (temp inicial)
    p0 = [0.02, np.max(temp) - np.min(temp), np.min(temp)]

    try:
        popt, _ = curve_fit(modelo_series_calor, t, temp, p0=p0, maxfev=5000)
        alpha_fit = popt[0]

        # 3. Calculo de Conductividad
        # Supongamos que: Conductividad = alpha_fit * Factor_Geometria
        # El factor puede depender del DT extraido del nombre
        dt_val, ts_val = extraer_parametros_nombre(ruta_archivo)

        # Constante de calibracion (ejemplo basado en tus resultados 0.03 aprox)
        constante_material = 0.0299
        conductividad = alpha_fit * constante_material * (dt_val / ts_val)

        return None, {
            "CASO": os.path.basename(ruta_archivo),
            "K2_FIT": alpha_fit,
            "CONDUCTIVIDAD": conductividad,
            "DT": dt_val,
            "TS": ts_val
        }
    except Exception as e:
        return f"Error en ajuste de {ruta_archivo}: {e}", None

def main():
    # Lista de archivos a procesar (puedes usar glob.glob para automatizar)
    archivos = ["S1644_DT5_TS0.5.log", "S444_DT25_TS1.0.log"]

    # Crear archivos de prueba si no existen para demostracion
    for arc in archivos:
        if not os.path.exists(arc):
            t_dummy = np.linspace(0, 100, 50)
            temp_dummy = 25 + 10 * (1 - 2 * np.exp(-0.02 * t_dummy)) + np.random.normal(0, 0.1, 50)
            np.savetxt(arc, np.column_stack((t_dummy, temp_dummy)), header="Tiempo Temperatura")

    resultados = []

    print("--- INICIANDO PROCESAMIENTO ---")
    for archivo in archivos:
        print(f"Procesando: {archivo}...")
        error, res = procesar_experimento(archivo)
        if res:
            resultados.append(res)
        else:
            print(error)

    # Mostrar Tabla Final
    df = pd.DataFrame(resultados)
    print("\n--- RESULTADOS FINALES ---")
    # Formateo de columnas para que coincida con tu salida
    print(df[["CASO", "K2_FIT", "CONDUCTIVIDAD"]].to_string(index=False, formatters={
        'K2_FIT': '{:,.6f}'.format,
        'CONDUCTIVIDAD': '{:,.6f}'.format
    }))

if __name__ == "__main__":
    main()
