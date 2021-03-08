const fs = require('fs');
let smallptwasm = require('./build/smallptwasm')
async function main(spp, sceneFile) {
    let a = await smallptwasm();
    let sceneData = fs.readFileSync(sceneFile);
    await a.renderJS(spp, sceneData);
}

let spp = process.argv[2] || 4
let sceneFile = process.argv[3] || "./sample-meshes/cornell.json"
main(spp, sceneFile);
