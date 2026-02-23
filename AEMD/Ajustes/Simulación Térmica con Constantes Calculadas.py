import numpy as np
import os
import re
from scipy.optimize import curve_fit

# --- CONFIGURACION DEL SISTEMA ---
N_TERMS = 10
L_ANG = 500.0  # Longitud en Angstroms
AREA_SECCION = 20.0 * 20.0 # Area en Angstroms cuadrados (ejemplo)
MASA_ATOMICA_UMAS = 28.085 # Ejemplo: Silicio (en daltons/u)
KB = 1.380649e-23 # Constante de Boltzmann (J/K)
NA = 6.022140e23  # Numero de Avogadro

def calcular_propiedades_material(masa_u, longitud_ang, area_ang):
    """
    Calcula la densidad y capacidad calorifica estimada.
    En sistemas atomicos a temp. ambiente, Cv approx 3 * kB por atomo (Dulong-Petit).
    """
    # 1. Volumen en m^3
    volumen_m3 = (longitud_ang * area_ang) * 1e-30

    # 2. Masa en kg (asumiendo que conocemos el numero de atomos, ej: 1000)
    # Aqui un ejemplo simplificado: densidad tipica del material si no hay conteo de atomos
    # O mejor: calcular masa total si sabes el numero de atomos en la simulacion.
    n_atomos = 1000 # Sustituir por valor real de tu simulacion
    masa_kg = (n_atomos * masa_u) / (NA * 1000)

    rho = masa_kg / volumen_m3 # kg/m^3

    # 3. Cp (Capacidad calorifica especifica J / kg*K)
    # Segun Dulong-Petit: C_molar = 3 * R = 3 * 8.314 J/mol*K
    cp = (3 * 8.314) / (masa_u / 1000) # J/kg*K

    return rho * cp # Esto es la capacidad calorifica volumetrica (J/m^3*K)

def f_total(x, k2, L, T_max, T_min):
    """Modelo de Fourier original."""
    total = 0
    L_m = L * 1e-10
    for n in range(1, N_TERMS + 1):
        k_n = (np.pi * n) / L_m
        k1_n = np.pi * n
        b_n = (1.0 / k1_n) * (-np.cos(k1_n) * (T_min - T_max) - T_max + T_min)
        term_c = (-2 * b_n + 2 * b_n * np.cos(k1_n)) * (1.0 / k1_n)
        total += np.exp(-k2 * x * (k_n**2)) * term_c
    return total

def procesar_archivo(ruta):
    try:
        data = np.loadtxt(ruta, skiprows=1)
        x_data, y_data = data[:, 0], data[:, 1]

        t_max, t_min = np.max(y_data), np.min(y_data)

        # Ajuste de difusividad (alpha / k2)
        fit_func = lambda x, k2: f_total(x, k2, L_ANG, t_max, t_min)
        popt, _ = curve_fit(fit_func, x_data, y_data, p0=[1e-7], maxfev=5000)
        alpha_fit = popt[0] # Unidades: m^2/s si x esta en segundos

        # CALCULO DE LA CONSTANTE DEL MATERIAL
        # En lugar de hardcodear, calculamos C_volumetrica
        c_vol = calcular_propiedades_material(MASA_ATOMICA_UMAS, L_ANG, AREA_SECCION)

        # Extraer parametros de simulacion DT/TS
        dt_val = float(re.search(r'DT(\d+\.?\d*)', ruta).group(1)) if "DT" in ruta else 1.0
        ts_val = float(re.search(r'TS(\d+\.?\d*)', ruta).group(1)) if "TS" in ruta else 1.0

        # Conductividad final: Kappa = alpha * C_vol * (escalado temporal)
        conductividad = alpha_fit * c_vol * (dt_val / ts_val)

        return {
            "Archivo": os.path.basename(ruta),
            "Difusividad_m2s": alpha_fit,
            "Conductividad_WmK": conductividad
        }
    except Exception as e:
        return f"Error: {e}"

# Lista de archivos para demostracion (suponiendo archivos .log en carpeta)
print("Procesando archivos y calculando constantes materiales...")
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
