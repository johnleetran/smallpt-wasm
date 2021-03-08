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
cp app.js  build/app.js

cd build
emcmake cmake ..
emcmake make

### web
#em++ -g0 libsmallptwasm.a -o smallptwasm.html --bind -lidbfs.js -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB #--pre-js pre.js

### nodejs
# em++ -g0 libsmallptwasm.a -o smallptwasm.js --bind -lnodefs.js -s "EXTRA_EXPORTED_RUNTIME_METHODS=['UTF8ToString']" -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -s DISABLE_EXCEPTION_CATCHING=0 -s FORCE_FILESYSTEM=1 #--pre-js pre.js

### nodejs standalone, meaning manual loading
em++ libsmallptwasm.a -o smallptwasm.js -s WASM=1 --bind -lnodefs.js -s "EXTRA_EXPORTED_RUNTIME_METHODS=['UTF8ToString']" -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -s DISABLE_EXCEPTION_CATCHING=0 -s FORCE_FILESYSTEM=1 -s WASM=1 -s MODULARIZE=1 #--pre-js pre.js



### example from action parsing project
#em++ -g0 libfilesystem.a -o filesystem.js --bind -lnodefs.js -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB --pre-js pre.js