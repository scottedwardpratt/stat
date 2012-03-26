#!/bin/zsh
find . -name '*.pdf' -print  | xargs -t tar -cf - | (cd ~/Dropbox/chemTreeN-project/round-robin-2; tar xf - )
