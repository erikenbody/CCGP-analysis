#!/bin/bash

src="/scratch2/erik/CCGP-reruns/projects"
dest="eenbody@courtyard.gi.ucsc.edu:/public/groups/corbettlab/eenbody/ccgp/CCGP-module"

for dir in $(find ${src} -type d -name "CCGP")
do
  project=$(dirname $(dirname $(dirname $dir))) # to get the $PROJECT directory
  project_name=$(basename $project) # to get the name of $PROJECT directory
  rsync -avz ${dir} ${dest}/${project_name}/
done
