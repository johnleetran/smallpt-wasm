# smallptWASM

## build

```
chmod +x ./setup-wasm.sh
./setup-wasm.sh
```

## to run
```
cd build
node app.js 4 "$(cat ../sample-meshes/cornell.json)"
```

## TODO
* bind filesystem in node.js and write output to file system [DONE]
* node.js to read the sample-meshes config file [BS - the file needs to be compiled into the wasm.]
