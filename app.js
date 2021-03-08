const fs = require('fs');
let smallptwasm = require('./smallptwasm')
async function main() {
    let a = await smallptwasm();
    let sceneData = fs.readFileSync("../sample-meshes/cornell.json");
    await a.renderJS(4, sceneData);
    console.log("here the bean yo")
}

main();
