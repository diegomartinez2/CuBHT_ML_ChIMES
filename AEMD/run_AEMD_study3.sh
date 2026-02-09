#!/bin/bash

###############################################################################
# SCRIPT DE AUTOMATIZACIÓN PARA ESTUDIO SISTEMÁTICO AEMD (LAMMPS)
###############################################################################

# 1. LISTAS DE PARÁMETROS
SIZES=("4 4 4" "8 4 4" "16 2 2" "16 4 4")
DELTA_TS=(5 15 25)
TIMESTEPS=(0.5 1.0 1.5)
REPLICAS=3

INPUT_BASE="Thermo_AEMD_CuBHT.in"
LAMMPS_EXEC="mpirun -np 70 /home/dmartinezgutierrez/lmp_mpi_chimes_diego_intel_env_new"

# 2. BUCLES DE EJECUCIÓN
for SIZE in "${SIZES[@]}"; do
    SIZE_NAME=$(echo $SIZE | tr -d ' ')

    for DT_VAL in "${DELTA_TS[@]}"; do
        T_HOT=$(awk -v dt="$DT_VAL" 'BEGIN { print 300.0 + dt }')
        T_COLD=$(awk -v dt="$DT_VAL" 'BEGIN { print 300.0 - dt }')

        for STEP in "${TIMESTEPS[@]}"; do
            for ((R=1; R<=REPLICAS; R++)); do

                WORK_DIR="run_S${SIZE_NAME}_DT${DT_VAL}_TS${STEP}_R${R}"
                OUTPUT_DAT="delta_T_S${SIZE_NAME}_DT${DT_VAL}_TS${STEP}_R${R}.dat"

                mkdir -p "$WORK_DIR"

                # Semillas aleatorias
                SEED1=$((RANDOM + R * 111))
                SEED2=$((RANDOM + R * 222))
                SEED3=$((RANDOM + R * 333))

                # Cálculo de NUEVO_VALOR (Tdamp)
                # Forzamos punto decimal y aseguramos que no sea 0
                NUEVO_VALOR=$(LC_NUMERIC=C awk -v s="$STEP" 'BEGIN { printf "%d", s * 100 }')

                echo "==> Configurando: $WORK_DIR | Tdamp: $NUEVO_VALOR | T_hot: $T_HOT"

                # --- MODIFICACIÓN DEL INPUT (SED REFORZADO) ---
                # 1. Usamos [[:space:]]* para manejar espacios variables
                # 2. Usamos [0-9.]* para capturar valores numéricos previos
                sed -e "s/^replicate.*/replicate $SIZE/" \
                    -e "s/^variable dt equal.*/variable dt equal $STEP/" \
                    -e "s/^\(fix NVT_eq all nvt temp \${T_eq} \${T_eq}\).*/\1 $NUEVO_VALOR/" \
                    -e "s/^variable T_hot_pulse equal.*/variable T_hot_pulse equal $T_HOT/" \
                    -e "s/^variable T_cold_pulse equal.*/variable T_cold_pulse equal $T_COLD/" \
                    -e "s/\(velocity all create \${T_eq}\) [0-9]*/\1 $SEED1/" \
                    -e "s/\(velocity HOT create \${T_hot_pulse}\) [0-9]*/\1 $SEED2/" \
                    -e "s/\(velocity COLD create \${T_cold_pulse}\) [0-9]*/\1 $SEED3/" \
                    -e "s|fix output all ave/time \(.*\) file .*|fix output all ave/time \1 file $OUTPUT_DAT|" \
                    "$INPUT_BASE" > "$WORK_DIR/input.in"

                # Copiar dependencias
                cp params.txt AA-optimized.lmp "$WORK_DIR/" 2>/dev/null

                # Ejecución
                cd "$WORK_DIR"
                if [ -s "input.in" ]; then
                    $LAMMPS_EXEC -in input.in > lammps.log 2>&1
                else
                    echo "ERROR: Falló la creación de input.in"
                fi
                cd ..

                echo "==> Finalizado: $WORK_DIR"
            done
        done
    done
done

echo "Estudio sistemático completado."
