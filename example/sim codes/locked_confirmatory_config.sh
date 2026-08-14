#!/usr/bin/env bash
# Fixed SILK simulation configuration shared by the one-replication validation
# and the 50-replication confirmatory run. Hyperparameters were fixed using a
# separate development design that is not used for confirmatory claims. Neither
# cluster entry point performs outcome-based tuning.

export SILK_CONFIG_ID="silk-confirmatory-2026-08-13-v2"
export SILK_HYPERPARAMETER_SELECTION="separate-development-design"

export SILK_N_TRAIN_GRID="200,400,1000"
export SILK_TP_N_TRAIN="400"
export SILK_N_TEST="1000"
export SILK_HORIZONS="1,2,3,4"

export SILK_BIOMARKER_KERNELS="gaussian,laplace,matern32"
export SILK_BIOMARKER_BW="median"
export SILK_BIOMARKER_BW_SCALE="1"
export SILK_SHIFT_GRID_STEP="0.2"
export SILK_PROFILE_SEARCH_RADIUS="2.0"
export SILK_ALT_LOCAL_RADIUS="2.0"
export SILK_N_STARTS="3"
export SILK_N_FOLDS="5"
export SILK_ENABLE_JMBAYES2="true"
export SILK_JMBAYES_N_CHAINS="1"
export SILK_JMBAYES_N_ITER="1000"
export SILK_JMBAYES_N_BURNIN="500"
export SILK_JMBAYES_PRED_N_SAMPLES="200"
export SILK_RSF_NUM_TREES="300"
export SILK_RSF_MIN_NODE_SIZE="15"
export SILK_DEEPSURV_EPOCHS="80"
export SILK_DEEPSURV_BATCH_SIZE="128"
export SILK_DEEPSURV_VALIDATION_FRACTION="0.20"
export SILK_DEEPSURV_PATIENCE="10"
export SILK_TIMEERROR_N_IMPUTE="5"

export REPORT_N_TRAIN="400"
export REPORT_N_BOOT="10000"
export REPORT_BOOT_SEED="20260530"
