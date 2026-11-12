#!/bin/sh
# usage: ulp_report.sh <build-dir>   (run from the repo root; CI's "ULP report" list)
B=$1; T=tests/data; rc=0
run() { exe=$1; shift; echo "### $exe $*"; "$B/tests/$exe" "$@" || { echo "### rc=$? $exe"; rc=1; }; }
run test_exp_dd $T/exp_dd_reference.txt
run test_log_dd $T/log_dd_reference.txt
run test_dd_special $T/dd_special_reference.txt $T/log1pmx_ulo_reference.txt $T/expm1_dd_reference.txt
run test_erf_ulp $T/erf_reference.txt
run test_erfc_ulp $T/erfc_reference.txt
run test_lgamma_ulp $T/lgamma_reference.txt
run test_erfinv_ulp $T/erfinv_reference.txt $T/erfcinv_reference.txt
run test_gamma_ulp $T/gamma_p_reference.txt $T/gamma_q_reference.txt
run test_beta_ulp $T/beta_p_reference.txt $T/beta_q_reference.txt
run test_digamma_ulp $T/digamma_reference.txt
run test_trigamma_ulp $T/trigamma_reference.txt
run test_gammainv_ulp $T/gammainv_p_reference.txt $T/gammainv_q_reference.txt
run test_betainv_ulp $T/betainv_p_reference.txt $T/betainv_q_reference.txt $T/betainv_backward_reference.txt
run test_bessel_ulp $T/i0_reference.txt $T/i1_reference.txt $T/i0e_reference.txt $T/i1e_reference.txt
run test_lbeta_ulp $T/lbeta_reference.txt
run test_trig_ulp $T/cos_reference.txt $T/sin_reference.txt
run test_exp_ulp $T/exp_reference.txt $T/expm1_reference.txt
run test_log_ulp $T/log_reference.txt $T/log1p_reference.txt $T/log1pmx_reference.txt
exit $rc
