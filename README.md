# Needleman-Wunsch BLOSUM62 Aligner

A 100% client-side global pairwise amino-acid alignment tool using the Needleman-Wunsch algorithm with BLOSUM62 scoring. Features VectorBuilder-style alignment visualization, expandable grouped variant blocks with per-substitution classification, quick examples, and a SARS-CoV-2 spike variant demo.

## Features
- **Global pairwise alignment** using the Needleman-Wunsch algorithm with BLOSUM62 substitution matrix.
- **Configurable gap penalty** with one-click reset.
- **Visual alignment output** with color-coded substitution types (conservative, neutral, non-conservative, gap).
- **Grouped variant blocks** — click any row to expand and inspect individual substitutions.
- **Quick examples** and a built-in SARS-CoV-2 spike variant demo.
- **Save & compare** multiple alignments in a session.
- **Export-friendly** — runs entirely in the browser with zero data transmission.

## Usage
1. Open `index.html` in a modern browser.
2. Paste two amino-acid sequences (or click a quick example / demo button).
3. Adjust the gap penalty if desired and click **Align Sequences**.
4. Inspect the alignment visualization, stats, and grouped variant blocks.

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

### Attribution

If you use this tool (in research, publications, pipelines, derivatives, or any other work), please credit me:

**Matt Moser**
GitHub: [https://github.com/mosermd-bio/nw-b62](https://github.com/mosermd-bio/nw-b62)

Example citation (for papers/posters/talks):
> Matt Moser. (2026). Needleman-Wunsch BLOSUM62 Aligner. GitHub repository. https://github.com/mosermd-bio/nw-b62

Thanks for your support — happy to discuss collaborations or improvements!
-Matt
