#!/bin/bash

awk 'OFS="\t" {print $1":"$2"-"$3, $1, $2, $3, "."}' $1
