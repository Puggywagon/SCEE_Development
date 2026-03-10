#!/bin/bash


# This script builds the figures for the manuscript.
#
# Dependencies:
# - '../../TIP4P_Results/5ns/Analysis/data.csv'
# - '../../Replica_Dipoles/Dipoles2.csv'
# - '../../TIP4P_Results/pcm_dipoles.csv'
#
# - './Dipoles2_lines.csv'
# - './NIST_vaporpressure.csv'
# - './TIP4P_vaporpressure.csv'
# - './NIST_isotherms.csv'


echo "plot_states.py"
#python ./plot_states.py --show False

echo "plot_density.py"
#python ./plot_density.py --show False

echo "plot_dipolestats.py"
#python ./plot_dipolestats.py --show False
python ./plot_dipoledist.py --show False
#python ./plot_dipoledistliq.py --show False
#python ./dipdist_nHB.py --show False
#python ./dipdist_donor.py --show False
#python ./dipdist_acceptor.py --show False
#python ./mu_nHB.py --show False
#python ./plot_SCEE_PCM.py --show False

echo "plot_rdf.py --T 298"
#python ./plot_rdf.py --T 298 --show False

echo "plot_rdf.py --T 500"
#python ./plot_rdf.py --T 500 --show False

echo "plot_rdf.py --T 700"
#python ./plot_rdf.py --T 700 --show False

echo "plot_rdf.py --T 600"
#python ./plot_rdf.py --T 600 --show False

echo "plot_rdf.py --T 800"
#python ./plot_rdf.py --T 800 --show False

echo "plot_rdf.py --T 400"
#python ./plot_rdf.py --T 400 --show False

echo "plot_rdf.py --T 900"
#python ./plot_rdf.py --T 900 --show False

echo "plot_rdf.py --T 1000"
#python ./plot_rdf.py --T 1000 --show False

echo "plot_hydrogenbond.py"
#python ./plot_hydrogenbond.py --show False

echo "plot_murho.py"
#python ./plot_murho.py --show False

echo "plot_mun.py"
#python ./plot_mun.py --show False

