// ---------- Utilities ----------
const ZERO_C = math.complex(0,0);
function c(re, im=0) { return math.complex(re, im); }
const cre = math.re, cim = math.im;

// ---------- Gate matrices (2x2) ----------
const SQRT1_2 = 1/Math.sqrt(2);
const GATES = {
  X: [[c(0,0), c(1,0)], [c(1,0), c(0,0)]],
  Y: [[c(0,0), c(0,-1)], [c(0,1), c(0,0)]],
  Z: [[c(1,0), c(0,0)], [c(0,0), c(-1,0)]],
  H: [[c(SQRT1_2,0), c(SQRT1_2,0)], [c(SQRT1_2,0), c(-SQRT1_2,0)]],
  S: [[c(1,0), c(0,0)], [c(0,0), c(0,1)]],                      // diag(1, i)
  Sdg: [[c(1,0), c(0,0)], [c(0,0), c(0,-1)]],                   // diag(1, -i)
  T: [[c(1,0), c(0,0)], [c(0,0), math.exp(math.multiply(c(0,1), Math.PI/4))]], // diag(1, e^{iπ/4})
  Tdg: [[c(1,0), c(0,0)], [c(0,0), math.exp(math.multiply(c(0,-1), Math.PI/4))]],
};

// Parameterized single-qubit gates
function Rx(theta){
  const th = theta/2;
  return [
    [c(Math.cos(th),0), c(0,-Math.sin(th))],
    [c(0,-Math.sin(th)), c(Math.cos(th),0)]
  ];
}
function Ry(theta){
  const th = theta/2;
  return [
    [c(Math.cos(th),0), c(-Math.sin(th),0)],
    [c(Math.sin(th),0), c(Math.cos(th),0)]
  ];
}
function Rz(theta){
  const th = theta/2;
  return [
    [math.exp(math.multiply(c(0,-1), th)), c(0,0)],
    [c(0,0), math.exp(math.multiply(c(0,1), th))]
  ];
}
function Phase(phi){
  return [[c(1,0), c(0,0)], [c(0,0), math.exp(math.multiply(c(0,1), phi))]];
}

// ---------- DOM elements ----------
const btnSet = document.getElementById('btnSet');
const afterSet = document.getElementById('afterSet');
const numQInput = document.getElementById('numQ');
const basisSelect = document.getElementById('basisState');
const basisRepresentation = document.getElementById('basisRepresentation');

const showBasisDensity = document.getElementById('showBasisDensity');
const basisDensityFormat = document.getElementById('basisDensityFormat');
const basisDensityContainer = document.getElementById('basisDensityContainer');
const btnToggleResultFormat = document.getElementById('btnToggleResultFormat');

const gateType = document.getElementById('gateType');
const singleTargetDiv = document.getElementById('singleTargetDiv');
const targetQ = document.getElementById('targetQ');

const angleDiv = document.getElementById('angleDiv');
const angleDeg = document.getElementById('angleDeg');
const angleLabel = document.getElementById('angleLabel');

const cnotDiv = document.getElementById('cnotDiv');
const controlQ = document.getElementById('controlQ');
const targetQ2 = document.getElementById('targetQ2');

const swapDiv = document.getElementById('swapDiv');
const swapA = document.getElementById('swapA');
const swapB = document.getElementById('swapB');

const ccnotDiv = document.getElementById('ccnotDiv');
const cc_c1 = document.getElementById('cc_c1');
const cc_c2 = document.getElementById('cc_c2');
const cc_t = document.getElementById('cc_t');

const btnAddGate = document.getElementById('btnAddGate');
const btnUndo = document.getElementById('btnUndo');
const btnClearGates = document.getElementById('btnClearGates');
const gatesListDiv = document.getElementById('gatesList');
const btnRun = document.getElementById('btnRun');

const resultsDiv = document.getElementById('results');
const blochSpheresDiv = document.getElementById('blochSpheres');

// ---------- App state ----------
let nQ = 2;
let stateVec = []; // array of math.complex
let gateSequence = []; // each gate object: {type, params, angle?}
let resultDensityAsTable = false; // toggles final/reduced density format

// ---------- Setup handlers ----------
btnSet.addEventListener('click', onSet);
gateType.addEventListener('change', onGateTypeChange);
btnAddGate.addEventListener('click', onAddGate);
btnUndo.addEventListener('click', onUndo);
btnClearGates.addEventListener('click', onClearGates);
btnRun.addEventListener('click', onRun);
basisSelect.addEventListener('change', updateBasisRepresentation);

showBasisDensity.addEventListener('change', () => {
  if (showBasisDensity.checked) renderAllBasisDensity(nQ);
  else basisDensityContainer.innerHTML = '';
});
basisDensityFormat.addEventListener('change', () => {
  if (showBasisDensity.checked) renderAllBasisDensity(nQ);
});
btnToggleResultFormat.addEventListener('click', () => {
  resultDensityAsTable = !resultDensityAsTable;
  // re-render results if already present
  if (stateVec && stateVec.length) {
    const rho = outerProduct(stateVec);
    const reducedList = [];
    for (let q=0; q<nQ; q++){
      reducedList.push(partialTrace(rho, nQ, q));
    }
    displayResults(stateVec, rho, reducedList);
  }
});

// initialize UI
onGateTypeChange();

// ---------- Functions ----------
function onSet(){
  nQ = parseInt(numQInput.value);
  if (!(nQ >=1 && nQ <=5)) { alert("Choose n between 1 and 5"); return; }
  populateBasis(nQ);
  populateQubitSelectors(nQ);
  initState(nQ);
  afterSet.classList.remove('hidden');
  gateSequence = [];
  renderGateList();
  resultsDiv.innerHTML = "<div class='small'>Initial state set. Add gates and click Run.</div>";
  blochSpheresDiv.innerHTML = "";
  updateBasisRepresentation();

  // If user had "show basis density" checked, refresh it
  if (showBasisDensity.checked) renderAllBasisDensity(nQ);
}

function populateBasis(n){
  basisSelect.innerHTML = "";
  for (let i=0;i< (1<<n); i++){
    const bits = i.toString(2).padStart(n, '0');
    const opt = document.createElement('option');
    opt.value = bits;
    // Compose tensor-product style: |b0⟩ ⊗ |b1⟩ ⊗ ... (q0 leftmost)
    const tensorParts = bits.split('').map(b => `|${b}⟩`).join(' ⊗ ');
    opt.text = `|${bits}⟩    (${tensorParts})`;
    opt.title = `${tensorParts}  — ordering q0 (left) to q${n-1} (right)`;
    basisSelect.appendChild(opt);
  }
  // default to |00...0⟩
  basisSelect.value = '0'.repeat(n);
}

function updateBasisRepresentation(){
  // show a clearer tensor-product representation for currently selected basis state
  const bits = basisSelect.value || '0'.repeat(nQ);
  const parts = bits.split('').map((b, idx) => {
    return `q${idx}: |${b}⟩`;
  }).join('  ⊗  ');
  basisRepresentation.textContent = `Tensor product (leftmost = q0): ${parts}`;
}

function populateQubitSelectors(n){
  const sels = [targetQ, controlQ, targetQ2, swapA, swapB, cc_c1, cc_c2, cc_t];
  sels.forEach(s => s.innerHTML = '');
  for (let i=0;i<n;i++){
    const createOpt = (label) => { const o=document.createElement('option'); o.value=i; o.text='q'+i; return o; };
    sels.forEach(s => s.appendChild(createOpt()));
  }
}

function initState(n){
  const dim = 1<<n;
  stateVec = Array(dim).fill(0).map(()=>c(0,0));
  const initIndex = parseInt(basisSelect.value || "0", 2);
  stateVec[initIndex] = c(1,0);
}

function onGateTypeChange(){
  const type = gateType.value;
  // hide all
  singleTargetDiv.classList.add('hidden');
  cnotDiv.classList.add('hidden');
  swapDiv.classList.add('hidden');
  ccnotDiv.classList.add('hidden');
  angleDiv.classList.add('hidden');

  // show relevant
  if (['X','Y','Z','H','S','Sdg','T','Tdg','Rx','Ry','Rz','Phase'].includes(type)){
    singleTargetDiv.classList.remove('hidden');
  }
  if (['Rx','Ry','Rz','Phase'].includes(type)){
    angleDiv.classList.remove('hidden');
    angleLabel.textContent = (type==='Phase') ? 'φ (degrees):' : 'θ (degrees):';
  }
  if (type === 'CNOT' || type === 'CZ'){
    cnotDiv.classList.remove('hidden');
  }
  if (type === 'SWAP'){
    swapDiv.classList.remove('hidden');
  }
  if (type === 'CCNOT'){
    ccnotDiv.classList.remove('hidden');
  }
}

function onAddGate(){
  const type = gateType.value;
  let gate = { type, params: [] };

  if (['X','Y','Z','H','S','Sdg','T','Tdg','Rx','Ry','Rz','Phase'].includes(type)){
    const t = parseInt(targetQ.value);
    gate.params = [t];
    if (['Rx','Ry','Rz','Phase'].includes(type)){
      gate.angle = (parseFloat(angleDeg.value) || 0) * Math.PI/180; // store radians
    }
  } else if (type === 'CNOT' || type === 'CZ'){
    const c = parseInt(controlQ.value), t = parseInt(targetQ2.value);
    if (c === t) { alert("Control and target must be different"); return; }
    gate.params = [c, t];
  } else if (type === 'SWAP'){
    const a = parseInt(swapA.value), b = parseInt(swapB.value);
    if (a === b) { alert("Choose two different qubits"); return; }
    gate.params = [a, b];
  } else if (type === 'CCNOT'){
    const c1 = parseInt(cc_c1.value), c2 = parseInt(cc_c2.value), t = parseInt(cc_t.value);
    const set = new Set([c1,c2,t]);
    if (set.size < 3) { alert("Controls and target must be all different"); return; }
    if (nQ < 3) { alert("CCNOT needs at least 3 qubits"); return; }
    gate.params = [c1, c2, t];
  }

  gateSequence.push(gate);
  renderGateList();
}

function onUndo(){
  gateSequence.pop();
  renderGateList();
}

function onClearGates(){
  gateSequence = [];
  renderGateList();
}

function renderGateList(){
  gatesListDiv.innerHTML = "";
  if (gateSequence.length === 0){
    gatesListDiv.innerHTML = "<div class='small'>No gates added yet.</div>";
    return;
  }
  gateSequence.forEach((g,i)=>{
    const d = document.createElement('div');
    d.className = "gate-item";

    const left = document.createElement('div');
    left.className = 'gate-left';
    let desc = `${i+1}. ${g.type}`;
    if (g.params?.length){
      desc += ` (${g.params.join(',')})`;
    }
    if (g.angle !== undefined){
      const deg = (g.angle*180/Math.PI).toFixed(2);
      desc += `, ${deg}°`;
    }
    left.textContent = desc;

    const right = document.createElement('div');
    right.className = 'gate-right';

    const up = document.createElement('button');
    up.textContent = '↑';
    up.onclick = ()=>{ if(i>0){ [gateSequence[i-1],gateSequence[i]]=[gateSequence[i],gateSequence[i-1]]; renderGateList(); } };

    const down = document.createElement('button');
    down.textContent = '↓';
    down.onclick = ()=>{ if(i<gateSequence.length-1){ [gateSequence[i+1],gateSequence[i]]=[gateSequence[i],gateSequence[i+1]]; renderGateList(); } };

    const rm = document.createElement('button');
    rm.textContent = 'Remove';
    rm.className = 'rm';
    rm.onclick = ()=>{ gateSequence.splice(i,1); renderGateList(); };

    right.appendChild(up);
    right.appendChild(down);
    right.appendChild(rm);

    d.appendChild(left);
    d.appendChild(right);
    gatesListDiv.appendChild(d);
  });
}

// ---------- Quantum ops ----------
function applySingleQubitGate(psi, n, target, U){
  const dim = psi.length;
  const out = Array(dim).fill(0).map(()=>c(0,0));
  for (let i=0;i<dim;i++){
    const bin = i.toString(2).padStart(n, '0');
    // note: target is index in bitstring where leftmost is 0
    const bit = parseInt(bin[target]);
    for (let j=0;j<2;j++){
      const newBin = bin.substring(0,target) + j.toString() + bin.substring(target+1);
      const idx = parseInt(newBin, 2);
      const coeff = U[j][bit];
      out[idx] = math.add(out[idx], math.multiply(coeff, psi[i]));
    }
  }
  return out;
}

function applyCNOT(psi, n, control, target){
  const dim = psi.length;
  const out = Array(dim).fill(0).map(()=>c(0,0));
  for (let i=0;i<dim;i++){
    const bin = i.toString(2).padStart(n, '0');
    if (bin[control] === '1'){
      const flippedBit = bin[target] === '1' ? '0' : '1';
      const newBin = bin.substring(0,target) + flippedBit + bin.substring(target+1);
      const idx = parseInt(newBin, 2);
      out[idx] = math.add(out[idx], psi[i]);
    } else {
      out[i] = math.add(out[i], psi[i]);
    }
  }
  return out;
}

function applyCZ(psi, n, control, target){
  const dim = psi.length;
  const out = Array(dim).fill(0).map(()=>c(0,0));
  for (let i=0;i<dim;i++){
    const bin = i.toString(2).padStart(n, '0');
    const phase = (bin[control]==='1' && bin[target]==='1') ? c(-1,0) : c(1,0);
    out[i] = math.add(out[i], math.multiply(phase, psi[i]));
  }
  return out;
}

function applySWAP(psi, n, a, b){
  if (a===b) return psi.slice();
  const dim = psi.length;
  const out = Array(dim).fill(0).map(()=>c(0,0));
  for (let i=0;i<dim;i++){
    const bin = i.toString(2).padStart(n, '0');
    if (bin[a] === bin[b]){
      out[i] = math.add(out[i], psi[i]);
    } else {
      const swapped = bin.substring(0, Math.min(a,b))
        + (a<b ? bin[b] : bin[a])
        + bin.substring(Math.min(a,b)+1, Math.max(a,b))
        + (a<b ? bin[a] : bin[b])
        + bin.substring(Math.max(a,b)+1);
      const idx = parseInt(swapped, 2);
      out[idx] = math.add(out[idx], psi[i]);
    }
  }
  return out;
}

function applyCCNOT(psi, n, c1, c2, t){
  const dim = psi.length;
  const out = Array(dim).fill(0).map(()=>c(0,0));
  for (let i=0;i<dim;i++){
    const bin = i.toString(2).padStart(n, '0');
    if (bin[c1]==='1' && bin[c2]==='1'){
      const flippedBit = bin[t] === '1' ? '0' : '1';
      const newBin = bin.substring(0,t) + flippedBit + bin.substring(t+1);
      const idx = parseInt(newBin, 2);
      out[idx] = math.add(out[idx], psi[i]);
    } else {
      out[i] = math.add(out[i], psi[i]);
    }
  }
  return out;
}

function outerProduct(psi){
  const dim = psi.length;
  const rho = Array(dim).fill(0).map(()=>Array(dim).fill(0).map(()=>c(0,0)));
  for (let i=0;i<dim;i++){
    for (let j=0;j<dim;j++){
      rho[i][j] = math.multiply(psi[i], math.conj(psi[j]));
    }
  }
  return rho;
}

function partialTrace(rho, n, target){
  const dim = rho.length;
  let red = [[c(0,0), c(0,0)], [c(0,0), c(0,0)]];
  for (let i=0;i<dim;i++){
    for (let j=0;j<dim;j++){
      const ib = i.toString(2).padStart(n,'0');
      const jb = j.toString(2).padStart(n,'0');
      let equalOther = true;
      for (let k=0;k<n;k++){ if (k!==target && ib[k]!==jb[k]) { equalOther=false; break; } }
      if (equalOther){
        const a = parseInt(ib[target]); const b = parseInt(jb[target]);
        red[a][b] = math.add(red[a][b], rho[i][j]);
      }
    }
  }
  return red;
}

function densityToBloch(red){
  const rho01 = red[0][1];
  const x = 2 * cre(rho01);
  const y = -2 * cim(rho01);
  const z = cre(red[0][0]) - cre(red[1][1]);
  return {x:x, y:y, z:z};
}

// ---------- Run simulation ----------
function onRun(){
  initState(nQ);
  let psi = stateVec.slice();

  for (const g of gateSequence){
    if (g.type in GATES){
      psi = applySingleQubitGate(psi, nQ, g.params[0], GATES[g.type]);
    } else if (g.type === 'Rx'){
      psi = applySingleQubitGate(psi, nQ, g.params[0], Rx(g.angle));
    } else if (g.type === 'Ry'){
      psi = applySingleQubitGate(psi, nQ, g.params[0], Ry(g.angle));
    } else if (g.type === 'Rz'){
      psi = applySingleQubitGate(psi, nQ, g.params[0], Rz(g.angle));
    } else if (g.type === 'Phase'){
      psi = applySingleQubitGate(psi, nQ, g.params[0], Phase(g.angle));
    } else if (g.type === 'CNOT'){
      psi = applyCNOT(psi, nQ, g.params[0], g.params[1]);
    } else if (g.type === 'CZ'){
      psi = applyCZ(psi, nQ, g.params[0], g.params[1]);
    } else if (g.type === 'SWAP'){
      psi = applySWAP(psi, nQ, g.params[0], g.params[1]);
    } else if (g.type === 'CCNOT'){
      psi = applyCCNOT(psi, nQ, g.params[0], g.params[1], g.params[2]);
    }
  }
  stateVec = psi;

  const rho = outerProduct(stateVec);
  const reducedList = [];
  for (let q=0; q<nQ; q++){
    const red = partialTrace(rho, nQ, q);
    reducedList.push(red);
  }

  displayResults(stateVec, rho, reducedList);
  drawAllBloch(reducedList);

  // if user wants basis density matrices visible, refresh them (they are independent of gate sequence)
  if (showBasisDensity.checked) renderAllBasisDensity(nQ);
}

// ---------- Display & plotting ----------
function displayResults(psi, rho, reducedList){
  resultsDiv.innerHTML = '';
  const dim = psi.length;

  // ----- Final state amplitudes -----
  let s = "<div class='result-block'><h3>Final state amplitudes (nonzero)</h3>";
  s += "<div style='overflow-x:auto;'><pre>";
  for (let i=0;i<dim;i++){
    const mag = Math.hypot(cre(psi[i]), cim(psi[i]));
    if (mag > 1e-9){
      const amp = `${cre(psi[i]).toFixed(6)}${cim(psi[i])>=0?'+':'-'}${Math.abs(cim(psi[i])).toFixed(6)}j`;
      s += `|${i.toString(2).padStart(nQ,'0')}> : ${amp}\n`;
    }
  }
  s += "</pre></div></div>";

  // ----- Full density matrix -----
  s += "<div class='result-block'><h3>Full density matrix ρ</h3>";
  if (resultDensityAsTable){
    s += "<div style='overflow-x:auto;'>" + formatMatrixAsTable(rho) + "</div>";
  } else {
    s += "<div style='overflow-x:auto;'><pre>" + formatComplexMatrix(rho) + "</pre></div>";
  }
  s += "</div>";

  // ----- Reduced density matrices -----
  for (let q=0;q<reducedList.length;q++){
    s += "<div class='result-block'>";
    s += `<h3>Reduced ρ (qubit ${q})</h3>`;
    if (resultDensityAsTable){
      s += "<div style='overflow-x:auto;'>" + formatMatrixAsTable(reducedList[q]) + "</div>";
    } else {
      s += "<div style='overflow-x:auto;'><pre>" + formatComplexMatrix(reducedList[q]) + "</pre></div>";
    }
    const bloch = densityToBloch(reducedList[q]);
    s += `<div class='small'>Bloch vector: (${bloch.x.toFixed(6)}, ${bloch.y.toFixed(6)}, ${bloch.z.toFixed(6)})</div>`;
    s += "</div>";
  }

  resultsDiv.innerHTML = s;
}


function formatComplexMatrix(mat){
  return mat.map(row => row.map(c=>`${cre(c).toFixed(6)}${cim(c)>=0?'+':'-'}${Math.abs(cim(c)).toFixed(6)}j`).join('\t')).join('\n');
}

function formatMatrixAsTable(mat){
  let html = "<table class='density-table'>";
  for (let row of mat){
    html += "<tr>";
    for (let cell of row){
      html += `<td>${cre(cell).toFixed(6)}${cim(cell)>=0?'+':'-'}${Math.abs(cim(cell)).toFixed(6)}j</td>`;
    }
    html += "</tr>";
  }
  html += "</table>";
  return html;
}

// ---------- Basis density rendering (new) ----------
function renderAllBasisDensity(n){
  // show density matrices for all basis states |b> where b runs 0..2^n-1
  // Format controlled by basisDensityFormat.value ('aligned' or 'table')
  const format = basisDensityFormat.value || 'aligned';
  const total = 1<<n;
  let html = '';
  for (let i=0;i<total;i++){
    const bits = i.toString(2).padStart(n,'0');
    html += `<div class="basis-density-block">`;
    html += `<div class="basis-density-title">|${bits}⟩  —  ( ${bits.split('').map(b => `|${b}⟩`).join(' ⊗ ')} )</div>`;

    // build full density matrix for this basis state: it's a matrix with 1 at (i,i) and 0 elsewhere
    const dim = 1<<n;
    const rho = Array(dim).fill(0).map(()=>Array(dim).fill(c(0,0)));
    rho[i][i] = c(1,0);

    if (format === 'table') {
      html += formatMatrixAsTable(rho);
    } else {
      // aligned pre format — reusing formatComplexMatrix but need to limit width for readability
      html += `<div class="aligned-pre">${formatComplexMatrix(rho)}</div>`;
    }

    html += `</div>`;
  }
  basisDensityContainer.innerHTML = html;
}

// ---------- Bloch plotting ----------
function drawAllBloch(reducedList){
  blochSpheresDiv.innerHTML = '';
  for (let q=0; q<reducedList.length; q++){
    const block = document.createElement('div');
    block.id = 'bloch_'+q; block.style.width = '350px'; block.style.height = '350px';
    blochSpheresDiv.appendChild(block);
    const bloch = densityToBloch(reducedList[q]);
    plotBloch(block.id, bloch, q);
  }
}

function plotBloch(containerId, bloch, q){
  const U = 30, V = 30;
  let xs = [], ys = [], zs = [];
  for (let i=0;i<=U;i++){
    const rowx=[], rowy=[], rowz=[];
    const theta=Math.PI*i/U;
    for (let j=0;j<=V;j++){
      const phi=2*Math.PI*j/V;
      rowx.push(Math.sin(theta)*Math.cos(phi));
      rowy.push(Math.sin(theta)*Math.sin(phi));
      rowz.push(Math.cos(theta));
    }
    xs.push(rowx); ys.push(rowy); zs.push(rowz);
  }
  const sphere = { type: 'surface', x: xs, y: ys, z: zs, opacity: 0.25, showscale:false, colorscale:'Blues' };
  const axes = [
    { type:'scatter3d', mode:'lines', x:[-1,1], y:[0,0], z:[0,0], line:{width:2,color:'gray'} },
    { type:'scatter3d', mode:'lines', x:[0,0], y:[-1,1], z:[0,0], line:{width:2,color:'gray'} },
    { type:'scatter3d', mode:'lines', x:[0,0], y:[0,0], z:[-1,1], line:{width:2,color:'gray'} }
  ];
  const vx = bloch.x, vy = bloch.y, vz = bloch.z;
  const vec = { type:'scatter3d', mode:'lines+markers', x:[0,vx], y:[0,vy], z:[0,vz],
                line:{width:6,color:'red'}, marker:{size:4,color:'red'} };
  const data = [sphere].concat(axes).concat([vec]);
  const layout = { title: `Qubit ${q}`, margin:{l:0,r:0,b:0,t:30},
                  scene: { aspectmode:'cube', xaxis:{range:[-1,1]}, yaxis:{range:[-1,1]}, zaxis:{range:[-1,1]} } };
  Plotly.newPlot(containerId, data, layout, {displayModeBar:false});
}

function drawAllBloch(reducedList){
  blochSpheresDiv.innerHTML = '';
  for (let q=0; q<reducedList.length; q++){
    // Container for this qubit
    const qubitContainer = document.createElement('div');
    qubitContainer.style.display = 'inline-block';
    qubitContainer.style.verticalAlign = 'top';
    qubitContainer.style.margin = '10px';

    // Sphere container
    const block = document.createElement('div');
    block.id = 'bloch_'+q; 
    block.style.width = '350px'; 
    block.style.height = '350px';
    qubitContainer.appendChild(block);

    const bloch = densityToBloch(reducedList[q]);

    // Defer plotting
    setTimeout(() => {
      plotBloch(block.id, bloch, q);
    }, 0);

    // Properties panel
    const propDiv = document.createElement('div');
    propDiv.className = 'bloch-properties';
    const btn = document.createElement('button');
    btn.textContent = 'Properties';
    btn.onclick = () => showBlochProperties(bloch, propDiv, q);
    propDiv.appendChild(btn);

    qubitContainer.appendChild(propDiv);

    // Add the whole container to main div
    blochSpheresDiv.appendChild(qubitContainer);
  }
}


function showBlochProperties(bloch, container, q){
  // remove old checkboxes except button
  const btn = container.querySelector('button');
  container.innerHTML = '';
  container.appendChild(btn);

  // Determine properties
  const len = Math.sqrt(bloch.x**2 + bloch.y**2 + bloch.z**2);
  const props = [
    {name:'Pure state', value: len > 0.9999},
    {name:'Mixed state', value: len <= 0.9999},
    {name:'More |0⟩ than |1⟩', value: bloch.z > 0.01},
    {name:'More |1⟩ than |0⟩', value: bloch.z < -0.01},
    {name:'Superposition along X-axis', value: Math.abs(bloch.x) > 0.01},
    {name:'Superposition along Y-axis', value: Math.abs(bloch.y) > 0.01},
    {name:'Along +Z axis', value: bloch.z > 0.99},
    {name:'Along -Z axis', value: bloch.z < -0.99}
  ];

  // Append checkboxes
  props.forEach(p=>{
    const label = document.createElement('label');
    const checkbox = document.createElement('input');
    checkbox.type = 'checkbox';
    checkbox.checked = p.value;
    checkbox.disabled = true;
    label.appendChild(checkbox);
    label.appendChild(document.createTextNode(p.name));
    container.appendChild(label);
  });
}
