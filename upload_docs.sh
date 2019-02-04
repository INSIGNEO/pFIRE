#!/bin/bash

(
cd $HOME/repos/pfire
export PATH="${PATH}:${PWD}/bin"
sphinx-build doc/ gh-pages/
cd gh-pages/
git add --all
git commit -m "$(date)"
git push upstream gh-pages
)
