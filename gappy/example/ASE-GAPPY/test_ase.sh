#!/bin/bash
echo ''
echo '**TESTING FIX CELL OPTIMIZATION'
echo ''
python run_opt.py
echo ''
echo '**TESTING NVT MD'
echo ''
python run_md.py
echo ''
echo '***FINISHED***'
echo ''

