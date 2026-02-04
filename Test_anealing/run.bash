# ===== 5 MINUTOS TOTAL =====
echo "FASE 1: Minimización..."
lmp_mpi_chimes -in min_delicado.in > min.log

echo "FASE 2: Annealing..."
lmp_mpi_chimes -in anneal_80K.in > anneal.log

echo "FASE 3: V_eq..."
lmp_mpi_chimes -in npt_final.in > npt.log

echo "RESULTADO:"
cat v_equilibrio.txt
tail -5 npt.log
