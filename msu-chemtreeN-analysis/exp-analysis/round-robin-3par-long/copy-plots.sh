#!/bin/zsh
find . -name '*.pdf' -print  | xargs -t tar -cf - | (cd ~/Dropbox/Second_folder/round-robin-3par-long; tar xf - )
