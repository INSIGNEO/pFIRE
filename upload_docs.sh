#!/bin/bash

(
basedir="$HOME/repos/pfire/"
docdir="doc"
ghpdir="gh-pages"


cd $basedir
(
cd $docdir
zip -r _static/pfire_tutorial_files.zip tutorial_files
)
export PATH="${PATH}:${basedir}/bin"
sphinx-build -an $docdir $ghpdir
cd $ghpdir
git add --all
git commit -m "$(date)"
git push upstream gh-pages
)
