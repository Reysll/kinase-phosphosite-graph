# Follow-up prompt for a new session

Paste the block below into a fresh Claude Code session in this repo. It is
written to be self-contained, so a session with no memory of the earlier work
can pick up safely.

Keep this file updated as the work moves. Everything it asserts is recorded in
`thesis/STATE.md`, which is the authoritative recovery point.

---

## Paste this

I am continuing a master's thesis and journal manuscript on this repo, a
heterogeneous graph network for kinase-substrate association prediction in
liver cancer. Earlier sessions did the analysis and validation work. Start by
reading, in this order:

1. `thesis/STATE.md` — the authoritative recovery point. 16 findings with the
   output file and exact command for each, 8 retracted claims, the advisor's
   2026-10-15 direction, and 7 open decisions.
2. `CLAUDE.md` — project instructions, environment, commands and current scale.
3. `thesis/references.md` — citation verification status.
4. `thesis/ADVISOR_BRIEF.md` — note its Options A-D are superseded by the
   advisor direction section in STATE.md.
5. `thesis/chapters/` — Chapter II and III drafts. Chapters I, IV, V, VI and
   VII are not written yet and are on hold.

### Standing constraints, which override any default behaviour

- **Do not modify any project code, data or results.** That means `src/`,
  `src_prediction/`, `data/`, `outputs/`, `outputs_prediction/` and `scripts/`
  are read-only. Write everything new under `thesis/`. The only exception on
  record was a single legend tweak to `build_ks_volcano.py` before this rule
  took effect.
- **Run no experiments or analyses** unless I ask in that message. The analysis
  phase is closed.
- **Accuracy.** Every number, metric, dataset size, parameter and result must
  come from the code, config or output files. Cite the file path in your notes
  to me, not in the thesis. Never invent results, sample sizes or performance
  values. If a value is missing or a run is incomplete, insert a visible
  placeholder like `[TODO: value from X run]`. If something I describe is not
  supported by the repo, tell me instead of writing around it.
- **Citations.** APA 7, prefer 2020-2025. Only cite sources you can verify.
  Give a DOI or a freely accessible link, PubMed Central when the DOI is
  paywalled. Flag anything you could not verify instead of guessing.
- **Style.** Plain, direct academic prose. No filler, no flattery, no
  AI-sounding transitions. No em dashes. No Oxford comma.
- **Framing.** Do not use the word "reversal" or present the kinase-level
  result as overturning the trial-level result. They are complementary
  analyses answering different questions.
- **Blocked topics.** Do not write about the context-conditioned kinase
  feature, the 15.2% vs 60% feasibility gate, the relational-only baseline,
  the drop-K diagnostic or the lbfgs flag. None are in the repo and I have not
  confirmed a source.
- **Interpreter.** `.venv312/Scripts/python`, never bare `python`.
- Thesis formatting target is UTRGV: Letter, Times New Roman 12pt black, bold
  only on headings 1-4, double-spaced body, 2in top margin on section-header
  pages and 1in elsewhere, page numbers centred 0.5in from the bottom, roman
  iii at the Abstract, arabic 1 at Chapter I, each chapter starting a new page.
  Use `[TODO:]` placeholders for committee names and dates.

### Where things stand

The validation round found that the project's headline results do not hold.
A degree-preserving null reproduces the Liver FC advantage in about 96% of
trials and 100% of the necroptosis cluster shift, seed noise exceeds the
context effect, kinase identity alone explains 99.97% of rank variance, and a
substrate-count popularity baseline beats every trained model. Eight prior
claims are retracted. All of this is documented with commands in STATE.md.

The advisor has since redirected the work: journal paper instead of a
conference, drop prediction performance as the headline metric in favour of
within-category comparative kinase ranking, and produce statistically
significant liver-cancer-associated substrates for a collaborator to validate
experimentally. Within-category ranking is already implemented and does not
support the hypothesis. No per-site differential test exists in the codebase
yet, so the collaborator deliverable is genuinely new work.

### Open work, in rough priority order

1. Reply to the advisor's 2026-10-15 email. Not drafted.
2. A graph-construction and integration walkthrough for the meeting he asked
   for. Not drafted. The two things he will most want to see are that 90.6% of
   the generic graph is pathway co-membership cliques, and that correlations
   are pairwise-complete with a mode of 6 overlapping samples.
3. Scope the per-site differential analysis for the collaborator's validation
   targets. This is the deliverable most robust to everything that failed
   validation.
4. Chapters I, IV, V, VI and VII, once the framing settles.
5. Checkpoint 4, the `.docx` build with python-docx plus PDF export and a
   page-by-page formatting check. Open decision 7 in STATE.md flags that
   Chapter II contains LaTeX equations and python-docx has no LaTeX support,
   so the conversion approach needs deciding first.
6. Unresolved provenance: `data/raw/high_confidence_score.csv` has no recorded
   source and must not be cited as STRING, and the liver file has 18 tumor
   plus 18 control samples against the 159-patient cohort Gao et al. report.

### Git

Work is on branch `docs/accuracy-fixes`, pushed to `origin`. It is not merged
to `main`, and whether to merge is open decision 4. Large regenerable outputs
under `thesis/analysis/out/` are gitignored and reproducible from the commands
in STATE.md.

Tell me what you have read and what you think the first move should be before
you change anything.
