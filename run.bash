#!/bin/bash
gfortran -fcheck=all -fbounds-check mod_kinds.F mod_scalars.F vandera.F  && ./a.out
