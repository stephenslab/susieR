#!/bin/bash
#
# This script will be executed every time you run "git commit". It
# will commit changes made to package DESCRIPTION by the pre-commit hook
#
# To use this script, copy it to the .git/hooks directory of your
# local repository to filename `post-commit`, and make it executable.
#
ROOT_DIR=`git rev-parse --show-toplevel`
# Only commit DESCRIPTION file when it is not staged (due to changes by pre-commit hook)
if [[ -z `git diff HEAD` ]] || [[ ! -f $ROOT_DIR/DESCRIPTION ]] || [[ -z `git diff $ROOT_DIR/DESCRIPTION` ]]; then
    exit 0
else
    git add $ROOT_DIR/DESCRIPTION
    git commit --amend -C HEAD --no-verify
    echo "Amend current commit to incorporate version bump"
fi
