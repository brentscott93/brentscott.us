#!/bin/bash 
cd ~/rstats/brentscott.us
git pull
r -e "blogdown::build_dir('~/rstats/brentscott.us/content/run', force = TRUE)"
git add -A
git commit -m "cron auto-update run log"
git push