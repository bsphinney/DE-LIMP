# Signal to Noise — Who stays with the chaperone?

- **Show:** Signal to Noise
- **Title:** Who stays with the chaperone? Hsp104 and Ssa1 pulldowns in heat-shocked yeast (synthetic example)
- **Hosts:** MAYA (cell biologist, voice Kore) · LEO (statistician, voice Charon)
- **Sources:** synthetic_report.md (a synthetic test fixture; no real study)
- **Written by:** a test author, in the exact layout of a real Signal to Noise script (2026-09-25). AI-generated; the written report is the record.

## Pronunciation

| Written | Spoken |
|---|---|
| Hsp104 | H S P 104 |
| Hsp26 | H S P 26 |
| Hsp42 | H S P 42 |
| Ssa1 | S S A 1 |
| Sis1 | S I S 1 |
| Ydj1 | Y D J 1 |
| KRT1 | K R T 1 |
| KRT10 | K R T 10 |
| DPC-Quant | D P C quant |
| PASEF | pasef |
| FDR | F D R |
| 0421 | zero four two one |

## Claims beyond the report

Nothing in the transcript states a fact that is not in the sources. Everything below is the script author's framing, so check it if strict fidelity matters:

- Metaphors and running jokes: "the bouncer", "who stays for the after-party", and "judging a singer from 4 notes".
- Background descriptions of the baits from general knowledge, not from the report: Hsp104 is "a disaggregase" and Ssa1 is "an Hsp70".
- Teaching lines that are general explanation rather than results from this report:
  - proteins are cut into peptides first because whole proteins are hard to measure;
  - a 1% FDR means about 1 in 100 precursor identifications is expected to be wrong;
  - proteins are reported as groups because some peptides are shared between related proteins;
  - empirical Bayes borrows information across proteins, which is why 4 replicates are workable.
- The hosts are synthetic voices.

<!-- TRANSCRIPT START -->

**MAYA:** Leo, picture a yeast cell at 42 degrees, proteins unfolding everywhere, and one machine whose whole job is to pull the clumps apart.

**LEO:** And I'm picturing a spreadsheet. 24 samples, 3,412 proteins, 5 comparisons, and a column called PropObs.

**MAYA:** Welcome to Signal to Noise. Quick note first: this is an AI-generated discussion, our voices are synthetic, and the written report is the record.

**LEO:** Today it's submission 0421, a synthetic example, so nothing here is a real result.

---

**LEO:** Set the table. What did they do?

**MAYA:** 2 baits, Hsp104 and Ssa1, plus IgG controls, from yeast early and late in a heat shock. 4 replicates per group, 6 groups.

**LEO:** The IgG is the denominator for everything. Remember that.

---

**MAYA:** Before the results, walk me through how the samples were measured. I have never touched a mass spectrometer.

**LEO:** Proteins are cut into peptides, the peptides are separated on a 44 minute chromatography gradient, and the mass spec weighs them and their fragments. On a timsTOF this was dia-PASEF: the instrument fragments everything in wide windows instead of picking peptides one at a time.

**MAYA:** And the software?

**LEO:** DIA-NN 2.6.1, a library-free search with match-between-runs, at a 1% FDR. That gave 41,236 precursors, which roll up into 3,412 protein groups. Groups, because some peptides are shared between related proteins.

**MAYA:** And 1% FDR means?

**LEO:** About 1 in 100 precursor identifications is expected to be wrong. That is the price of admission.

---

**MAYA:** First question any biologist asks. Did the pulldowns work?

**LEO:** Yes. Hsp104 is the top hit in its own pulldown, log2FC 9.84, adjusted p of 3.1e-12, hundreds of times more than in IgG. Ssa1 tops its pulldown too.

**MAYA:** The is-the-phone-plugged-in check.

---

**MAYA:** Early heat shock, Hsp104 against IgG: 412 significant proteins, 367 of them up.

**LEO:** Late: 158. Hold that drop, we will come back to it.

**MAYA:** And Sis1 and Ydj1 come down with Ssa1, which is lovely.

---

**LEO:** Here is the caveat. The Late IgG runs are thin: 2,104 and 2,388 proteins, against 2,961 for Early IgG.

**MAYA:** So a thinner control makes everything look enriched.

**LEO:** And Hsp104 itself drops -1.1 in Late, its 96 partners by a median of -0.8. They move together, so it is mostly bait recovery.

---

**MAYA:** Did anything late survive?

**LEO:** Hsp26 and Hsp42 rise in the Late Hsp104 pulldowns, log2FC 2.35 and 1.92, and Hsp42 was seen in fewer than half of the runs. Leads, not conclusions.

---

**MAYA:** Detective story?

**LEO:** Keratins, KRT1 and KRT10, in 3 of 24 runs. Removed as contaminants before quantification.

---

**MAYA:** Nerd moment. Go.

**LEO:** Between 18 and 41 percent of each sample's values are inferred, median 27 percent. The detection-probability model fills in what was never observed. And empirical Bayes borrows variance information across thousands of proteins, which is why 4 replicates are workable. Then Benjamini-Hochberg turns the p-values into adjusted p, because we tested 3,412 proteins at once.

---

**MAYA:** So what should the lab do with this?

**LEO:** Open README.html, then Analysis_Report.html, then the DE_dpc tables. Tier the hits by the Detected columns: detected in 4 of 4 Hsp104 runs goes on the slide, detected in 1 of 4 goes on the follow-up list.

**MAYA:** And validate Hsp26 and Hsp42 first. That's Signal to Noise. I'm Maya.

**LEO:** I'm Leo. Check your controls.

<!-- TRANSCRIPT END -->
