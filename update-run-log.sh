#!/bin/bash 
cd ~/rstats/brentscott.us
git pull
r -e "blogdown::build_site(local = FALSE, run_hugo = FALSE, build_rmd = 'content/run/2021/2021.Rmd')"
git add -A
git commit -m "cron auto-update run log"
git push