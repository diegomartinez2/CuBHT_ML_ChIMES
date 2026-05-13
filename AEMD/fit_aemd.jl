#!/usr/bin/env julia

"""
fit_aemd.jl
Traducción del script de ajuste AEMD para extraer difusividad (α) y conductividad térmica (κ).
"""

using ArgParse
using DelimitedFiles
using Statistics
using LinearAlgebra
using Plots

# --- Funciones de Lectura ---

function read_geom(fname="geom.info")
    !isfile(fname) && error("Archivo de geometría no encontrado: $fname")

    lines = readlines(fname)
    usable_lines = filter(ln -> !isempty(strip(ln)) && !startswith(strip(ln), "#"), lines)

    isempty(usable_lines) && error("No se encontraron datos utilizables en $fname")

    parts = split(usable_lines[end])
    length(parts) < 3 && error("Se esperaban al menos 3 columnas en $fname")

    N = floor(Int, parse(Float64, parts[1]))
    V_A3 = parse(Float64, parts[2])
    Lz_A = parse(Float64, parts[3])

    if N <= 0 || V_A3 <= 0.0 || Lz_A <= 0.0
        error("Valores físicos no válidos: N=$N, V_A3=$V_A3, Lz_A=$Lz_A")
    end

    return N, V_A3, Lz_A
end

function read_media_aemd_deltaT(fname="media_aemd_deltaT.dat")
    !isfile(fname) && error("Archivo ΔT no encontrado: $fname")

    # readdlm lee matrices ignorando comentarios por defecto si se especifica
    data = readdlm(fname, Float64, comments=true, comment_char='#')

    isempty(data) && error("No se pudieron leer datos de $fname")

    ncol = size(data, 2)
    # Julia usa indexación basada en 1
    if ncol == 6
        return data[:,1], data[:,2], nothing, nothing, data[:,5], data[:,6]
    elseif ncol == 5
        return data[:,1], nothing, data[:,2], data[:,3], data[:,4], data[:,5]
    elseif ncol == 4
        return data[:,1], data[:,2], nothing, nothing, data[:,3], data[:,4]
    elseif ncol == 3
        return data[:,1], nothing, nothing, nothing, data[:,2], data[:,3]
    else
        error("Se esperaban 3, 4, 5 o 6 columnas en $fname, se obtuvieron $ncol")
    end
end

# --- Modelo Físico y Ajuste ---

function model_deltaT(t_ps, alpha_A2_per_ps, dT0, Lz_A; nterms=20)
    nterms < 1 && error("nterms debe ser >= 1")

    # Broadcast de operaciones para manejar arrays de tiempo
    result = zeros(length(t_ps))
    odd = [2.0 * i + 1.0 for i in 0:(nterms-1)]
    coef = 1.0 ./ (odd.^2)
    lam2 = (odd .* pi ./ Lz_A).^2

    # Cálculo eficiente de la serie
    for i in eachindex(t_ps)
        expo = exp.(-lam2 .* alpha_A2_per_ps .* t_ps[i])
        result[i] = sum(coef .* expo)
    end

    return (8.0 * dT0 / (pi^2)) .* result
end

function fit_alpha(t_ps, dT, Lz_A; nterms=20)
    length(t_ps) < 10 && error("No hay suficientes puntos para el ajuste.")

    dT_abs0 = abs(dT[1])
    (iszero(dT_abs0) || !isfinite(dT_abs0)) && error("ΔT(0) es cero o no finito.")

    # Estimación inicial de tau (cruce 1/e)
    target = dT_abs0 / exp(1)
    idx = findfirst(x -> abs(x) < target, dT)
    tau_ps = isnothing(idx) ? t_ps[end] : t_ps[idx]
    tau_ps = max(tau_ps, 1e-6)

    alpha_guess = (Lz_A^2) / (pi^2 * tau_ps)

    # Búsqueda en rejilla (Grid Search)
    alphas_coarse = alpha_guess .* 10 .^ range(-2.0, 2.0, length=121)
    best = (mse=Inf, alpha=0.0, dT0=0.0)

    for a in alphas_coarse
        base = model_deltaT(t_ps, a, 1.0, Lz_A, nterms=nterms)
        denom = dot(base, base)
        denom <= 1e-30 && continue

        dT0 = dot(base, dT) / denom
        pred = dT0 .* base
        mse = mean((pred .- dT).^2)

        if mse < best.mse
            best = (mse=mse, alpha=a, dT0=dT0)
        end
    end

    # Refinamiento local
    a0 = best.alpha
    for _ in 1:6
        alphas_fine = a0 .* 10 .^ range(-0.25, 0.25, length=41)
        best_local = (mse=Inf, alpha=0.0, dT0=0.0)

        for a in alphas_fine
            base = model_deltaT(t_ps, a, 1.0, Lz_A, nterms=nterms)
            denom = dot(base, base)
            denom <= 1e-30 && continue

            dT0 = dot(base, dT) / denom
            pred = dT0 .* base
            mse = mean((pred .- dT).^2)

            if mse < best_local.mse
                best_local = (mse=mse, alpha=a, dT0=dT0)
            end
        end
        best_local.mse == Inf && break
        best = best_local
        a0 = best.alpha
    end

    return best.alpha, best.dT0, best.mse
end

# --- Función Principal ---

function main()
    s = ArgParseSettings(description = "Ajuste AEMD en Julia.")
    @add_arg_table! s begin
        "--geom"
            default = "geom.info"
            help = "Archivo de geometría LAMMPS."
        "--media"
            default = "media_aemd_deltaT.dat"
            help = "Archivo de datos ΔT."
        "--dt_fs"
            arg_type = Float64
            default = 1.0
            help = "Timestep en fs."
        "--nterms"
            arg_type = Int
            default = 20
            help = "Términos de la serie."
        "--tmin_ps"
            arg_type = Float64
            default = 0.0
        "--tmax_ps"
            arg_type = Float64
            default = nothing
        "--plotfile"
            default = "aemd_fit.png"
    end

    args = parse_args(s)

    # Lectura de datos
    N, V_A3, Lz_A = read_geom(args["geom"])
    step, time_ps, Th, Tc, dT, Std_DeltaT = read_media_aemd_deltaT(args["media"])

    # Eje de tiempo
    t_ps = isnothing(time_ps) ? step .* (args["dt_fs"] * 1e-3) : time_ps

    # Filtrado (Máscara)
    mask = isfinite.(dT) .& isfinite.(t_ps) .& (t_ps .>= args["tmin_ps"])
    if !isnothing(args["tmax_ps"])
        mask .&= (t_ps .<= args["tmax_ps"])
    end

    t_ps_fit = t_ps[mask]
    dT_fit = dT[mask]

    # Ordenar por tiempo
    p = sortperm(t_ps_fit)
    t_ps_fit = t_ps_fit[p]
    dT_fit = dT_fit[p]

    # Ajuste
    alpha_A2ps, dT0_fit, mse = fit_alpha(t_ps_fit, dT_fit, Lz_A, nterms=args["nterms"])

    # Conversiones SI
    alpha_SI = alpha_A2ps * 1.0e-8
    kB = 1.380649e-23
    V_m3 = V_A3 * 1.0e-30
    n_dens = N / V_m3
    Cvol = 3.0 * n_dens * kB
    kappa_WmK = alpha_SI * Cvol

    # Gráfico con Plots.jl
    t_smooth = range(minimum(t_ps_fit), maximum(t_ps_fit), length=600)
    dT_smooth = model_deltaT(t_smooth, alpha_A2ps, dT0_fit, Lz_A, nterms=args["nterms"])

    p_plot = scatter(t_ps, dT, label="ΔT data (all)", markersize=2, alpha=0.3, color=:blue)
    scatter!(p_plot, t_ps_fit, dT_fit, label="ΔT data (fit window)", markersize=3, alpha=0.8, color=:red)
    plot!(p_plot, t_smooth, dT_smooth, label="Fit (series model)", linewidth=2, color=:black)

    xlabel!(p_plot, "time (ps)")
    ylabel!(p_plot, "ΔT (K)")
    title!(p_plot, "AEMD: ΔT(t) data and fit")
    savefig(p_plot, args["plotfile"])

    # Resultados
    println("=== AEMD fit results (Julia) ===")
    @printf("N atoms             : %d\n", N)
    @printf("Volume (Å^3)        : %.6e\n", V_A3)
    @printf("Lz (Å)              : %.6f\n", Lz_A)
    @printf("alpha (m^2/s)       : %.6e\n", alpha_SI)
    @printf("kappa (W/m/K)       : %.6f\n", kappa_WmK)
    @printf("Saved plot          : %s\n", args["plotfile"])
end

# Necesario para el formateo de strings estilo C
using Printf

main()
