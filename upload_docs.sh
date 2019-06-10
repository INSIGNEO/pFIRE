#!/bin/bash

(
basedir="$PWD"
docdir="${basedir}/doc/"
ghpdir="${basedir}/gh-pages"


cd $docdir
rm _static/pfire_tutorial_files.zip
zip --exclude=tutorial_files/faces_1/sad2happy.xdmf.h5 \
    -r \
    _static/pfire_tutorial_files.zip \
    tutorial_files
export PATH="${PATH}:${basedir}/bin:${basedir}/reglab"
sphinx-build -an . $ghpdir
cd $ghpdir
git add --all
git commit -m "$(date)"
git push upstream gh-pages
)
