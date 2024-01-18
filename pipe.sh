#!/bin/sh

chmod +x *.py *.sh
pip3 install -r requirements.txt

echo "[0/4] syncing"
./sync.sh

echo "[1/4] preparing data"
./utils.py

echo "[2/4] image plotting"
./plot.py 10-image
./plot.py 10-image i,r,g

echo "[3/4] source detection"
./source.py 10-
./source.py

echo "[4/4] isophote fitting"
./isophote.py 10-
