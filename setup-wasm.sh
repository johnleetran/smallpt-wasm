#!/bin/bash

pushd

function goBackDir {
	popd
}

trap goBackDir EXIT
mkdir -p build

set -eux
# cp index.html build/index.html
# cp pre.js  build/pre.js

cd build
emcmake cmake ..
emcmake make

### web
#em++ -g0 libsmallptwasm.a -o smallptwasm.html --bind -lidbfs.js -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB #--pre-js pre.js

### nodejs
em++ -g0 libsmallptwasm.a -o smallptwasm.js --bind -lnodefs.js -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB #--pre-js pre.js

### example from action parsing project
#em++ -g0 libfilesystem.a -o filesystem.js --bind -lnodefs.js -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB --pre-js pre.js