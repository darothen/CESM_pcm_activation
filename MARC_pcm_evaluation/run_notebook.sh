#!/usr/bin/env bash

jupyter-nbconvert \
   --ExecutePreprocessor.timeout=600 \
   --to notebook \
   --execute $1
