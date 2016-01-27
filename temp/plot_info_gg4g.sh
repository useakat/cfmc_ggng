#!/bin/bash

plot_unwgt=("pt_g1" "pt_g2" "pt_g3" "pt_g4" "y_g1" "y_g2" "y_g3" "y_g4" \
#"phi_g1" "phi_g2" "phi_g3" "phi_g4" 
"dr_g1g2" "dr_g1g3" "dr_g1g4" "dr_g2g3" "dr_g2g4" "dr_g3g4" "m_g1g2" "m_g1g3" "m_g1g4" "m_g2g3" "m_g2g4" "m_g3g4") 
#"kt_g1g2" "kt_g1g3" "kt_g1g4" "kt_g2g3" "kt_g2g4" "kt_g3g4")

xlabel=("pT (g1)" "pT (g2)" "pT (g3)" "pT (g4)" "y (g1)" "y (g2)" "y (g3)" "y (g4)" \
#"phi (g1)" "phi (g2)" "phi (g3)" "phi (g4)" 
"dR (g1g2)" "dR (g1g3)" "dR (g1g4)" "dR (g2g3)" "dR (g2g4)" "dR (g3g4)" "m (g1g2)" "m (g1g3)" "m (g1g4)" "m (g2g3)" \
"m (g2g4)" "m (g3g4)")
#"kT (g1g2)" "kT (g1g3)" "kT (g1g4)" "kT (g2g3)" "kT (g2g4)" "kT (g3g4)")

binsize=("2." "2." "2." "2." "0.2" "0.2" "0.2" "0.2" \
#"0.13" "0.13" "0.13" "0.13" \
"0.1" "0.1" "0.1" "0.1" "0.1" "0.1" \
"2." "2." "2." "2." "2." "2.")
# "2." "2." "2." "2." "2." "2.")

xshift=(0 0 0 0  -5 -5 -5 -5 \
#0 0 0 0  \
0 0 0 0 0 0  0 0 0 0 0 0)
#  0 0 0 0 0 0)