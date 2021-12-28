#!/bin/bash 
cd ~/rstats/brentscott.us
git pull
r -e "blogdown::build_dir('~/rstats/brentscott.us/content/run', ignore = /^2021.Rmd/)"
git add -A
git commit -m "cron auto-update run log"
git push