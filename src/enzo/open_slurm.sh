#!/bin/bash

srun --job-name=compile --pty --nodes=1  --ntasks-per-node=8 -p gpu --gpus=1 bash
