#!/usr/bin/env bash

mkdir backup_inlists
cp inlist* backup_inlists

if [[ "$(uname)" == "Darwin" ]]; then
  # Macs require a different syntax for sed
  sed -i '' 's/log_center_density_limit/log_center_density_upper_limit/' inlist*
  sed -i '' 's/log_center_temp_limit/log_center_temp_upper_limit/' inlist*
  sed -i '' 's/center_entropy_limit/center_entropy_upper_limit/' inlist*
  sed -i '' 's/max_entropy_limit/max_entropy_upper_limit/' inlist*
  sed -E -i '' '/inlist[1-5]/s/([1-5])([_a-z]*) =/\2\(\1\) =/' inlist*
else
  sed 's/log_center_density_limit/log_center_density_upper_limit/' -i inlist*
  sed 's/log_center_temp_limit/log_center_temp_upper_limit/' -i inlist*
  sed 's/center_entropy_limit/center_entropy_upper_limit/' -i inlist*
  sed 's/max_entropy_limit/max_entropy_upper_limit/' -i inlist*
  sed -E '/inlist[1-5]/s/([1-5])([_a-z]*) =/\2\(\1\) =/' -i inlist*
fi
