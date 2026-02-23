import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# Configuración para servidores sin pantalla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def f_total(x, k, A, T_offset):
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

# Ejemplo de ejecución masiva
def procesar_lote_ejemplo():
    # Simulamos una lista de casos
    casos = ["S1644_DT5_TS0.5", "S444_DT25_TS1.0"]
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

if __name__ == "__main__":
    procesar_lote_ejemplo()
