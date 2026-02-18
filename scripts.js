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

// Theme
const themeToggle = document.getElementById('themeToggle');
const themeIcon   = document.getElementById('themeIcon');
function setTheme(isDark) {
    document.documentElement.classList.toggle('dark', isDark);
    themeIcon.textContent = isDark ? '🌙' : '☀️';
    localStorage.setItem('theme', isDark ? 'dark' : 'light');
}
themeToggle.addEventListener('click', () => setTheme(!document.documentElement.classList.contains('dark')));
if (localStorage.getItem('theme') === 'dark' || (!localStorage.getItem('theme') && window.matchMedia('(prefers-color-scheme: dark)').matches)) setTheme(true);

// Quick examples
const examples = [
    { label: "Short test 1", s1: "AIHV", s2: "-V-I" },
    { label: "Short test 2", s1: "GVYY", s2: "D---" }
];
examples.forEach(ex => {
    const btn = document.createElement('button');
    btn.className = "px-4 py-2 text-sm font-medium bg-gray-200 dark:bg-gray-700 hover:bg-gray-300 dark:hover:bg-gray-600 rounded-xl transition-colors";
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
    let html = '<div class="font-mono text-base leading-none tracking-[0.5px]">';
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
                const cls = sc >= 2 ? 'bg-emerald-200 dark:bg-emerald-900/60 text-emerald-800 dark:text-emerald-200' : sc >= 0 ? 'bg-amber-200 dark:bg-amber-900/60 text-amber-800 dark:text-amber-200' : 'bg-rose-200 dark:bg-rose-900/60 text-rose-800 dark:text-rose-200';
                htmlB1 += `<span class="${cls}">${c1}</span>`;
                htmlB2 += `<span class="${cls}">${c2}</span>`;
            } else if (c1 === '-' || c2 === '-') {
                const gCls = 'bg-sky-100 dark:bg-sky-900/60 text-sky-700 dark:text-sky-300';
                htmlB1 += (c1 === '-') ? `<span class="${gCls}">-</span>` : c1;
                htmlB2 += (c2 === '-') ? `<span class="${gCls}">-</span>` : c2;
            } else {
                htmlB1 += c1; htmlB2 += c2;
            }
        }
        const start1 = pos1 + 1;
        const end1 = pos1 + (b1.match(/[^-]/g) || []).length;
        const start2 = pos2 + 1;
        const end2 = pos2 + (b2.match(/[^-]/g) || []).length;
        html += `
            <div class="flex items-baseline mb-0.5">
                <span class="w-14 text-right pr-4 text-gray-400 text-xs">${start1}</span>
                <span class="flex-1">${htmlB1}</span>
                <span class="w-14 text-left pl-4 text-gray-400 text-xs">${end1}</span>
            </div>
            <div class="flex items-baseline mb-1 text-amber-500"><span class="w-14"></span><span class="flex-1">${matchB}</span><span class="w-14"></span></div>
            <div class="flex items-baseline mb-7">
                <span class="w-14 text-right pr-4 text-gray-400 text-xs">${start2}</span>
                <span class="flex-1">${htmlB2}</span>
                <span class="w-14 text-left pl-4 text-gray-400 text-xs">${end2}</span>
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
        let type, cls;
        if (hasSub === 0) {
            if (gaps1 > 0 && gaps2 === 0) { type = 'Deletion'; cls = 'bg-sky-100 dark:bg-sky-900/50 text-sky-800'; }
            else if (gaps2 > 0 && gaps1 === 0) { type = 'Insertion'; cls = 'bg-sky-100 dark:bg-sky-900/50 text-sky-800'; }
            else { type = 'Indel'; cls = 'bg-sky-100 dark:bg-sky-900/50 text-sky-800'; }
        } else if (gaps1 === 0 && gaps2 === 0) {
            type = 'Substitution Block'; cls = 'bg-amber-100 dark:bg-amber-900/50 text-amber-800';
        } else {
            type = 'Complex (Indel + Sub)'; cls = 'bg-purple-100 dark:bg-purple-900/50 text-purple-800';
        }
        blocks.push({ posStart, posEnd, block1, block2, len, type, cls, startIdx: start, hasSub });
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
    let html = `<thead><tr class="bg-gray-100 dark:bg-gray-800"><th class="w-8"></th><th>Positions</th><th>Seq1</th><th>Seq2</th><th>Len</th><th>Type</th></tr></thead><tbody>`;
    blocks.forEach((b, idx) => {
        html += `<tr class="variant-row cursor-pointer hover:bg-gray-50 dark:hover:bg-gray-800 border-b border-gray-100 dark:border-gray-700" data-index="${idx}">
            <td class="text-center"><span class="chevron inline-block w-4">›</span></td>
            <td class="font-mono">${b.posStart}–${b.posEnd}</td>
            <td class="font-mono">${b.block1}</td>
            <td class="font-mono">${b.block2}</td>
            <td class="font-mono text-center">${b.len}</td>
            <td>${b.type}</td>
        </tr>`;
    });
    if (blocks.length === 0) html += `<tr><td colspan="6" class="text-center py-8 text-gray-500">Perfect match — no variants</td></tr>`;
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
            let detailHTML = `<tr class="details-row bg-gray-50 dark:bg-gray-900"><td colspan="6" class="p-0"><div class="p-4">`;
            if (block.hasSub > 0) {
                detailHTML += `<table class="w-full text-xs"><thead><tr class="text-gray-500"><th class="text-left">Pos</th><th class="text-left">From → To</th><th class="text-left">Score</th><th class="text-left">Classification</th></tr></thead><tbody>`;
                for (let k = block.startIdx; k < block.startIdx + block.len; k++) {
                    const a = currentResult.alignment1[k];
                    const b = currentResult.alignment2[k];
                    if (a !== b && a !== '-' && b !== '-') {
                        const sc = new NeedlemanSearch().blosumScore(a, b);
                        const cls = sc >= 2 ? 'bg-emerald-100 dark:bg-emerald-900/40 text-emerald-700 dark:text-emerald-300' : sc >= 0 ? 'bg-amber-100 dark:bg-amber-900/40 text-amber-700 dark:text-amber-300' : 'bg-rose-100 dark:bg-rose-900/40 text-rose-700 dark:text-rose-300';
                        detailHTML += `<tr class="${cls}"><td class="font-mono">${k+1}</td><td class="font-mono">${a} → ${b}</td><td class="font-mono">${sc}</td><td>${new NeedlemanSearch().classify(sc)}</td></tr>`;
                    }
                }
                detailHTML += `</tbody></table>`;
            } else {
                detailHTML += `<p class="text-sm text-gray-600 dark:text-gray-400">Pure indel block — no amino-acid substitutions to classify.</p>`;
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
        historyDiv.innerHTML = '<p class="text-gray-500 text-center py-8">No saved comparisons yet.</p>';
        return;
    }
    history.forEach((item, idx) => {
        const div = document.createElement('div');
        div.className = "flex justify-between items-center bg-gray-50 dark:bg-gray-800 p-4 rounded-xl cursor-pointer hover:bg-gray-100 dark:hover:bg-gray-700";
        div.innerHTML = `<div class="font-mono text-sm">${item.seq1.substring(0,20)}${item.seq1.length>20?'…':''} vs ${item.seq2.substring(0,20)}${item.seq2.length>20?'…':''}</div><div class="text-right"><div class="font-bold text-green-600">${item.score}</div><div class="text-xs text-gray-500">${item.identity}% id</div></div>`;
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