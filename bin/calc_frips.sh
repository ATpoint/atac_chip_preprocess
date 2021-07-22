#!/bin/bash

signal=$(grep -w 'Assigned' $1 | cut -f2)
noise=$(grep -w 'Unassigned_NoFeatures' $1 | cut -f2)
bc <<< "scale=6;${signal}/(${signal}+${noise})"

