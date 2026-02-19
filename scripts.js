// scripts.js — nw62 Needleman-Wunsch BLOSUM62 Aligner (FULL VERSION)

const aminoAcids = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'];

const blosum62 = [
    [4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0],
    [-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3],
    [-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3],
    [0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2],
    [-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2],
    [0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3],
    [-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1],
    [-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2],
    [-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2],
    [1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2],
    [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1],
    [0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4]
];

class NeedlemanSearch {
    constructor({ gapPenalty = -4 } = {}) { this.gapPenalty = gapPenalty; }
    blosumScore(a, b) {
        const i = aminoAcids.indexOf(a.toUpperCase());
        const j = aminoAcids.indexOf(b.toUpperCase());
        if (i === -1 || j === -1) return this.gapPenalty;
        return blosum62[i][j];
    }
    classify(score) {
        if (score >= 2) return 'Conservative';
        if (score >= 0) return 'Neutral';
        return 'Non-conservative';
    }
    align(seq1, seq2) {
        const m = seq1.length, n = seq2.length;
        const M = Array.from({ length: m + 1 }, () => Array(n + 1).fill(0));
        for (let i = 0; i <= m; i++) M[i][0] = i * this.gapPenalty;
        for (let j = 0; j <= n; j++) M[0][j] = j * this.gapPenalty;
        for (let i = 1; i <= m; i++) for (let j = 1; j <= n; j++) {
            const diag = M[i-1][j-1] + this.blosumScore(seq1[i-1], seq2[j-1]);
            const up = M[i-1][j] + this.gapPenalty;
            const left = M[i][j-1] + this.gapPenalty;
            M[i][j] = Math.max(diag, up, left);
        }
        let align1 = '', align2 = '';
        let i = m, j = n;
        while (i > 0 && j > 0) {
            const score = M[i][j];
            const sDiag = M[i-1][j-1], sUp = M[i-1][j], sLeft = M[i][j-1];
            if (score === sDiag + this.blosumScore(seq1[i-1], seq2[j-1])) {
                align1 = seq1[i-1] + align1; align2 = seq2[j-1] + align2; i--; j--;
            } else if (score === sUp + this.gapPenalty) {
                align1 = seq1[i-1] + align1; align2 = '-' + align2; i--;
            } else {
                align1 = '-' + align1; align2 = seq2[j-1] + align2; j--;
            }
        }
        while (i > 0) { align1 = seq1[i-1] + align1; align2 = '-' + align2; i--; }
        while (j > 0) { align1 = '-' + align1; align2 = seq2[j-1] + align2; j--; }
        return { score: M[m][n], alignment1: align1, alignment2: align2 };
    }
}

// DOM
const seq1Input   = document.getElementById('seq1');
const seq2Input   = document.getElementById('seq2');
const gapInput    = document.getElementById('gapPenalty');
const alignBtn    = document.getElementById('alignBtn');
const clearBtn    = document.getElementById('clearBtn');
const resetBtn    = document.getElementById('resetConfig');
const saveBtn     = document.getElementById('saveBtn');
const demoBtn     = document.getElementById('demoBtn');
const resultsDiv  = document.getElementById('results');
const vizDiv      = document.getElementById('alignmentViz');
const scoreEl     = document.getElementById('score');
const idEl        = document.getElementById('identity');
const simEl       = document.getElementById('similarity');
const gapsEl      = document.getElementById('gaps');
const warningEl   = document.getElementById('warning');
const variantTable = document.getElementById('variantTable');
const historyDiv  = document.getElementById('history');
const countEl     = document.getElementById('count');
const exampleContainer = document.getElementById('exampleBtns');

let currentResult = null;
let history = JSON.parse(localStorage.getItem('nw62_history')) || [];

// Quick examples
const examples = [
    { label: "Short test 1", s1: "AIHV", s2: "-V-I" },
    { label: "Short test 2", s1: "GVYY", s2: "D---" }
];
examples.forEach(ex => {
    const btn = document.createElement('button');
    btn.className = "btn btn-muted btn-sm";
    btn.textContent = ex.label;
    btn.onclick = () => { seq1Input.value = ex.s1; seq2Input.value = ex.s2; alignAndShow(); };
    exampleContainer.appendChild(btn);
});

// Demo sequences
const SPIKE1 = `MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT`;

const SPIKE2 = `MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHVISGTNGTKRFDNPVLPFNDGVYFASIEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLDHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPIIVREPEDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFDEVFNATRFASVYAWNRKRISNCVADYSVLYNLAPFFTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKVSGNYNYLYRLFRKSNLKPFERDISTEIYQAGNKPCNGVAGFNCYFPLRSYSFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLKGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQGVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEYVNNSYECDIPIGAGICASYQTQTKSHRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLKRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKYFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFKGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNHNAQALNTLVKQLSSKFGAISSVLNDIFSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT`;

// Demo button
demoBtn.addEventListener('click', () => {
    seq1Input.value = SPIKE1;
    seq2Input.value = SPIKE2;
    alignAndShow();
});

// Helper functions
function buildMatchLine(a1, a2) {
    let identity = 0, similarity = 0;
    for (let i = 0; i < a1.length; i++) {
        if (a1[i] === a2[i]) { identity++; similarity++; }
        else if (a1[i] !== '-' && a2[i] !== '-') {
            const sc = new NeedlemanSearch().blosumScore(a1[i], a2[i]);
            if (sc > 0) similarity++;
        }
    }
    return {
        identity: a1.length ? ((identity / a1.length) * 100).toFixed(1) : '0',
        similarity: a1.length ? ((similarity / a1.length) * 100).toFixed(1) : '0'
    };
}

function renderWrappedAlignment(align1, align2) {
    const BLOCK = 60;
    let html = '<div>';
    let pos1 = 0, pos2 = 0;
    for (let start = 0; start < align1.length; start += BLOCK) {
        const end = Math.min(start + BLOCK, align1.length);
        const b1 = align1.substring(start, end);
        const b2 = align2.substring(start, end);
        let matchB = '';
        let htmlB1 = '', htmlB2 = '';
        for (let k = 0; k < b1.length; k++) {
            const c1 = b1[k], c2 = b2[k];
            if (c1 === c2 && c1 !== '-') matchB += '|';
            else if (c1 !== '-' && c2 !== '-') matchB += (new NeedlemanSearch().blosumScore(c1, c2) > 0 ? ':' : '.');
            else matchB += ' ';
            if (c1 !== c2 && c1 !== '-' && c2 !== '-') {
                const sc = new NeedlemanSearch().blosumScore(c1, c2);
                const cls = sc >= 2 ? 'sub-conservative' : sc >= 0 ? 'sub-neutral' : 'sub-noncons';
                htmlB1 += `<span class="${cls}">${c1}</span>`;
                htmlB2 += `<span class="${cls}">${c2}</span>`;
            } else if (c1 === '-' || c2 === '-') {
                htmlB1 += (c1 === '-') ? `<span class="sub-gap">-</span>` : c1;
                htmlB2 += (c2 === '-') ? `<span class="sub-gap">-</span>` : c2;
            } else {
                htmlB1 += c1; htmlB2 += c2;
            }
        }
        const start1 = pos1 + 1;
        const end1 = pos1 + (b1.match(/[^-]/g) || []).length;
        const start2 = pos2 + 1;
        const end2 = pos2 + (b2.match(/[^-]/g) || []).length;
        html += `
            <div class="aln-row aln-row-gap">
                <span class="aln-pos aln-pos-left">${start1}</span>
                <span class="aln-seq">${htmlB1}</span>
                <span class="aln-pos aln-pos-right">${end1}</span>
            </div>
            <div class="aln-row aln-row-match"><span class="aln-pos"></span><span class="aln-seq">${matchB}</span><span class="aln-pos"></span></div>
            <div class="aln-row aln-row-block">
                <span class="aln-pos aln-pos-left">${start2}</span>
                <span class="aln-seq">${htmlB2}</span>
                <span class="aln-pos aln-pos-right">${end2}</span>
            </div>
        `;
        pos1 = end1; pos2 = end2;
    }
    html += '</div>';
    return html;
}

function getVariantBlocks(aln1, aln2) {
    const blocks = [];
    let i = 0;
    while (i < aln1.length) {
        if (aln1[i] === aln2[i]) { i++; continue; }
        const start = i;
        while (i < aln1.length && aln1[i] !== aln2[i]) i++;
        const block1 = aln1.substring(start, i);
        const block2 = aln2.substring(start, i);
        const len = i - start;
        const posStart = start + 1;
        const posEnd = i;
        let gaps1 = 0, gaps2 = 0, hasSub = 0;
        for (let k = 0; k < block1.length; k++) {
            if (block1[k] === '-') gaps1++;
            if (block2[k] === '-') gaps2++;
            if (block1[k] !== '-' && block2[k] !== '-' && block1[k] !== block2[k]) hasSub++;
        }
        let type;
        if (hasSub === 0) {
            if (gaps1 > 0 && gaps2 === 0) type = 'Deletion';
            else if (gaps2 > 0 && gaps1 === 0) type = 'Insertion';
            else type = 'Indel';
        } else if (gaps1 === 0 && gaps2 === 0) {
            type = 'Substitution Block';
        } else {
            type = 'Complex (Indel + Sub)';
        }
        blocks.push({ posStart, posEnd, block1, block2, len, type, startIdx: start, hasSub });
    }
    return blocks;
}

function alignAndShow() {
    let s1 = seq1Input.value.trim().toUpperCase().replace(/[^ARNDCQEGHILKMFPSTWYV-]/g, '');
    let s2 = seq2Input.value.trim().toUpperCase().replace(/[^ARNDCQEGHILKMFPSTWYV-]/g, '');
    const gap = parseInt(gapInput.value) || -4;
    if (!s1 || !s2) { alert("Please enter both sequences"); return; }
    const searcher = new NeedlemanSearch({ gapPenalty: gap });
    currentResult = searcher.align(s1, s2);

    vizDiv.innerHTML = renderWrappedAlignment(currentResult.alignment1, currentResult.alignment2);

    const stats = buildMatchLine(currentResult.alignment1, currentResult.alignment2);
    scoreEl.textContent = currentResult.score;
    idEl.textContent = stats.identity;
    simEl.textContent = stats.similarity;

    let gapCount = 0;
    for (let k = 0; k < currentResult.alignment1.length; k++) if (currentResult.alignment1[k] === '-' || currentResult.alignment2[k] === '-') gapCount++;
    const gapsPct = ((gapCount / currentResult.alignment1.length) * 100).toFixed(1);
    gapsEl.textContent = `${gapCount} (${gapsPct}%)`;

    const lenRatio = Math.abs(s1.length - s2.length) / Math.max(s1.length, s2.length);
    warningEl.classList.toggle('hidden', lenRatio <= 0.6);

    const blocks = getVariantBlocks(currentResult.alignment1, currentResult.alignment2);
    let html = `<thead><tr><th class="text-center" style="width:2rem"></th><th>Positions</th><th>Seq1</th><th>Seq2</th><th>Len</th><th>Type</th></tr></thead><tbody>`;
    blocks.forEach((b, idx) => {
        html += `<tr class="variant-row" data-index="${idx}">
            <td class="text-center"><span class="chevron">›</span></td>
            <td class="mono">${b.posStart}–${b.posEnd}</td>
            <td class="mono">${b.block1}</td>
            <td class="mono">${b.block2}</td>
            <td class="mono text-center">${b.len}</td>
            <td>${b.type}</td>
        </tr>`;
    });
    if (blocks.length === 0) html += `<tr><td colspan="6" class="text-center text-dim" style="padding:2rem 0">Perfect match — no variants</td></tr>`;
    html += `</tbody>`;
    variantTable.innerHTML = html;

    document.querySelectorAll('.variant-row').forEach(row => {
        row.addEventListener('click', () => {
            const index = parseInt(row.getAttribute('data-index'));
            const block = blocks[index];
            const nextRow = row.nextElementSibling;
            if (nextRow && nextRow.classList.contains('details-row')) {
                nextRow.remove();
                row.classList.remove('expanded');
                return;
            }
            let detailHTML = `<tr class="details-row"><td colspan="6"><div class="detail-inner">`;
            if (block.hasSub > 0) {
                detailHTML += `<table class="w-full text-xs"><thead><tr class="text-dim"><th class="text-left">Pos</th><th class="text-left">From → To</th><th class="text-left">Score</th><th class="text-left">Classification</th></tr></thead><tbody>`;
                for (let k = block.startIdx; k < block.startIdx + block.len; k++) {
                    const a = currentResult.alignment1[k];
                    const b = currentResult.alignment2[k];
                    if (a !== b && a !== '-' && b !== '-') {
                        const sc = new NeedlemanSearch().blosumScore(a, b);
                        const cls = sc >= 2 ? 'sub-conservative' : sc >= 0 ? 'sub-neutral' : 'sub-noncons';
                        detailHTML += `<tr class="${cls}"><td class="mono">${k+1}</td><td class="mono">${a} → ${b}</td><td class="mono">${sc}</td><td>${new NeedlemanSearch().classify(sc)}</td></tr>`;
                    }
                }
                detailHTML += `</tbody></table>`;
            } else {
                detailHTML += `<p class="text-sm text-dim">Pure indel block — no amino-acid substitutions to classify.</p>`;
            }
            detailHTML += `</div></td></tr>`;
            row.insertAdjacentHTML('afterend', detailHTML);
            row.classList.add('expanded');
        });
    });

    resultsDiv.classList.remove('hidden');

    currentResult.seq1 = s1;
    currentResult.seq2 = s2;
    currentResult.gapPenalty = gap;
    currentResult.identity = stats.identity;
    currentResult.similarity = stats.similarity;
}

function saveComparison() {
    if (!currentResult) return;
    history.unshift({ ...currentResult });
    if (history.length > 20) history.pop();
    localStorage.setItem('nw62_history', JSON.stringify(history));
    renderHistory();
}

function renderHistory() {
    historyDiv.innerHTML = '';
    countEl.textContent = `(${history.length})`;
    if (history.length === 0) {
        historyDiv.innerHTML = '<p class="text-dim text-center" style="padding:2rem 0">No saved comparisons yet.</p>';
        return;
    }
    history.forEach((item, idx) => {
        const div = document.createElement('div');
        div.className = "history-item";
        div.innerHTML = `<div class="mono text-sm">${item.seq1.substring(0,20)}${item.seq1.length>20?'…':''} vs ${item.seq2.substring(0,20)}${item.seq2.length>20?'…':''}</div><div class="text-right"><div class="font-bold stat-green">${item.score}</div><div class="text-xs text-dim">${item.identity}% id</div></div>`;
        div.onclick = () => { seq1Input.value = item.seq1; seq2Input.value = item.seq2; gapInput.value = item.gapPenalty; alignAndShow(); };
        historyDiv.appendChild(div);
    });
}

// Event listeners
alignBtn.addEventListener('click', alignAndShow);
saveBtn.addEventListener('click', saveComparison);
clearBtn.addEventListener('click', () => { seq1Input.value = ''; seq2Input.value = ''; resultsDiv.classList.add('hidden'); });
resetBtn.addEventListener('click', () => gapInput.value = '-4');
document.addEventListener('keydown', e => { if (e.key === 'Enter' && !e.shiftKey) { e.preventDefault(); alignAndShow(); } });

// Init
renderHistory();
