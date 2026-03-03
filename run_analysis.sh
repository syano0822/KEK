#!/bin/bash

RUNNUMBER=$1
root 'make_event_waveforms.C("${RUNNUMBER}.root","waveform_${RUNNUMBER}.root","events",-1,-10.,-50.,20)'
root 'make_hit_position.C("waveform_$RUNNUMBER.root","pos_${RUNNUMBER}.root",8.,14.5,0.5,900,13.,19.,0.,300.,true)'
root 'draw_hit_position.C("$RUNNUMBER","thr_lin")'
