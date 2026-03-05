#!/bin/bash

RUN=$1
ISCOL4=$2
EXP=2026-03_RARiS
RAW_DIR=~/data/${EXP}/physics
ROOT_DIR=~/data/${EXP}_root


python /opt/eudaq/user/epic/scripts/LeCroyDump.py ${RAW_DIR}/run${RUN}.raw ${ROOT_DIR}/run${RUN}.root --output-type root 
root make_event_waveforms.C\(\"${ROOT_DIR}/run${RUN}.root\",\"waveform_${RUN}.root\",\"events\",-1,-10.,-50.,20\)
root make_hit_position.C\(\"waveform_${RUN}.root\",\"pos_${RUN}.root\",8.,14.5,0.5,600,13.,19.,0.,300.,false,-6.,6.,$ISCOL2\)
#root 'draw_hit_position.C("pos_${RUN}.root","thr_lin")'
