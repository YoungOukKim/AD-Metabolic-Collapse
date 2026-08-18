# Changelog

This file records what changed in the deposited material during the revision
round, and why. Commit messages carry the mechanics; this file carries the
reasoning, because several entries are corrections to earlier deposited files
and a reader needs to know that the earlier version said something different.

Entries are newest first. Version numbers refer to the revision package supplied
to the journal, not to git tags.

---

AD-Metabolic-Collapse -- repository v24
=======================================

1. output/revision2/Supporting_Data_Values_v8.xlsx     [replaces v7]

   Two sheets still held material the manuscript had corrected.

   Fig4_data carried the superseded extraction. Recomputing Figure 4 from that
   sheet returned R2 = 0.303 and p = 5.9e-08, the values that were printed in
   the figure before it was redrawn, rather than the 0.290 and 1.3e-07 the
   manuscript reports. The Data availability statement says the workbook holds
   the source data for every figure, so a reader checking Figure 4 against its
   own source would have reproduced the superseded numbers. The sheet is
   rebuilt from data/sample/donor_level_summary.csv and
   output/revision2/S1_donor_stage_full.csv and now returns 0.2896 and 1.3e-07.

   SuppFig2_data panel A was headed "Event timeline" and listed "MCT4 decline
   onset (10%) = 0.3", "PTGDS/EAAT2/ATP1A2 compensatory peak" and "PTGDS
   collapse". Every one of those is a claim the revision withdrew: the words
   onset, event timeline, compensatory and collapse were removed from the paper
   under reviewer 1's fifth point, and the manuscript states in as many words
   that the earlier report of MCT4 onset at Bin 0.3 is corrected to Bin 0.5,
   with Bin 0.3 belonging to TFRC. A reader who read that sentence and then
   opened the source data would have found the corrected claim still standing
   there. The panel is rewritten as positions along the axis, and the note
   records what the earlier version said. The donor block lower in the same
   sheet was also the superseded extraction and is replaced.

2. Table 1 of the manuscript, the ABC row
   The gap between the MCT4 and astrocytic V-ATPase rank correlations read
   +0.14, computed from the superseded values. On the corrected extraction it
   is 0.318 - 0.192 = 0.126, so +0.13. The Braak gap is unchanged at +0.08.

---

AD-Metabolic-Collapse -- repository v23
=======================================

1. figures/README.md                                   [new]
   README.md                                           [figures note added]

   The figures/ directory holds four files where the paper has ten display
   items, and nothing said why. The answer is that the directory is the output
   of the three Python redraw scripts, while the R scripts write to
   output/figures/, which is empty in a fresh clone because they have not been
   run. A reader counting files would reasonably conclude that six figures were
   missing.

   figures/README.md now states which figure comes from which script, and which
   can be regenerated from the deposited data alone: Figures 2, 4 and 6 from the
   Python scripts and Figures 1, 3 and Supplemental Figures 1 and 2 from the R
   scripts, all from deposited CSVs; Figures 5 and 7 need ADNI proteomics, which
   requires ADNI access; and Supplemental Figure 3 is a schematic whose printed
   numbers are itemised in the Supporting Data Values workbook.

   Fig4_Redraw.py was also missing from the file listing in README.md and is
   added.

---

AD-Metabolic-Collapse -- repository v22
=======================================

1. figures/Figure4.{png,tif}                          [replaced]
   R/figures/Fig4_Redraw.py                           [new]

   Figure 4 was the one display item that had not been regenerated from the
   corrected extraction. It was still drawn from the version in which log1p had
   been applied twice, so panel B printed R2 = 0.303 and p = 5.9e-08 inside the
   plot while its own legend, the body text and Table 1 gave 0.290 and 1.3e-07.
   A reader comparing the number in the panel with the number in the sentence
   below it would have found them different.

   Fig4_Redraw.py rebuilds all three panels from data/sample/donor_level_summary.csv
   and output/revision2/S1_donor_stage_full.csv, and asserts every published
   value before drawing. The values now printed:

     ANLS              Braak -0.201  CERAD -0.239  ABC -0.307
     MCT4              Braak -0.339  CERAD -0.380  ABC -0.318
     V-ATPase (astro)  Braak -0.260  CERAD -0.203  ABC -0.192
     panel B  R2 0.2896, p 1.3e-07, beta -0.0737
     panel C  Mann-Whitney p 0.0020, dementia n = 42

   Nothing in the reading changes: every sign and ordering holds and MCT4 is
   still the strongest and most consistent of the three. One significance band
   moves, MCT4 against Braak from *** to ** as p crosses 0.001.

2. R/replication/82_rosmap_braak_matched_donors.R      [paths]
   The script hard-coded D:/work and told the reader to send output files back
   to us. Both were left over from the run that produced the result and are not
   meaningful to anyone else. The path now comes from the AD_WORK environment
   variable and defaults to the working directory, and the note about sending
   files back is replaced by a statement of why the per-donor file is not
   redistributable.

---

AD-Metabolic-Collapse -- repository v21
=======================================

The first push of the revision round. Everything below v21 in this file was
prepared between the initial deposit and this push, so a reader comparing the
two commits sees a large diff at once; this entry explains the parts of it that
are not simply additions.

1. data/sample/  -- four files change value, not shape

   astro_bin_means.csv, neuron_bin_means.csv, donor_level_summary.csv and
   astro_subtype_trajectories.csv all differ from the versions deposited with
   the original submission. Row counts, bin keys and donor counts are unchanged;
   the numbers themselves moved.

   The cause is an extraction error found during the revision. An earlier
   version of the extraction applied log1p to a matrix that SEA-AD already
   stores normalised, so every value passed through that transform twice. Of the
   352 values shared between the two versions of astro_bin_means.csv, 348 are
   smaller in the old file, which is the signature of the second compression.

   The manuscript reports the corrected values. Table 5, for example, gives
   astrocytic SLC16A3 as 0.0633 early and 0.0360 late; the corrected files give
   0.0633 and 0.0360, and the previously deposited files gave 0.0468 and 0.0269.
   The percentage changes are almost unaffected, which is why the error survived
   the first round: it moves correlations and regression slopes, not ratios.

   The old files are left in the history rather than rewritten. Anyone checking
   a number from the submitted version against this repository should use the
   commit that this entry belongs to or later.

2. output/ is tracked from this commit onward

   The initial .gitignore excluded output/ entirely. The manuscript's Code
   availability statement now points readers to output/revision2/ for the
   pre-registered verdict rules written before the data were opened
   (S12_RULES.txt) and the verdicts as recorded at run time (S12_VERDICTS.txt),
   so that directory has to be present. Files derived from restricted datasets
   at individual or per-donor resolution stay excluded; the table in README.md
   lists them.

3. output/tables/rosmap_donor_mct4.csv is not committed

   It holds per-donor astrocytic MCT4 for 430 ROSMAP donors keyed by ROSMAP
   identifier. ROSMAP is released under a RADC data use agreement and per-donor
   derived values are not redistributable. The script that builds it is
   committed, and so is the three-row result it feeds
   (output/tables/82_rosmap_braak_matched.csv).

4. Two scripts renamed

     Fig5_CSF_Validation.R  ->  Fig5_CSF_Orthogonal_Test.R
     SuppFig2_Temporal.R    ->  SuppFig2_Axis_Position.R

   Both old names asserted something the revision withdrew: the CSF analysis is
   an orthogonal test rather than validation of the mechanism, and the axis is
   cross-sectional so no within-individual ordering is claimed. Headers inside
   Fig4_Clinical.R, Fig6_CrossCohort_Mediation.R and its Python counterpart were
   changed for the same reason.

5. CHANGES_v15.txt through CHANGES_v19.txt are consolidated into this file.

---

AD-Metabolic-Collapse -- repository v19
=======================================

Change from v18. One defect, in a sheet added at v18.

1. output/revision2/Supporting_Data_Values_v7.xlsx     [replaces v6]

   The SuppFig3_data sheet added at v18 was wrong in two ways.

   It labelled -0.8% as the neuronal V-ATPase change. -0.8% is the astrocytic
   value; the neuronal value is -5.4%. Stats_Summary carries both correctly and
   the figure itself prints both, so the error was confined to that one sheet.

   The sheet also ended with a note stating that every number printed on
   Supplemental Figure 3 was listed above. It was not. Eight printed values were
   missing: HK2 -35%, LDHA -21%, CP -45%, FTH1 -19%, the 93.5% of astrocyte
   subtypes declining in the same direction, the TFRC network degree of 16, the
   neuronal V-ATPase -5.4% and the 7x ratio quoted beside it.

   The sheet is rebuilt from the figure itself, box by box, and now lists every
   printed value with the sheet it comes from. The three statements on the figure
   that are not numbers -- lactate delivery inferred rather than measured,
   neuronal ATP unmeasured, and the V1A-Tau correlation not used as evidence --
   are listed as well, so a reader can see that nothing on the figure is
   unaccounted for.

   No other sheet changed and no value elsewhere in the workbook moved.

---

AD-Metabolic-Collapse -- repository v18
=======================================

Changes from v17. These follow a full audit of the revision package against the
reviewers' letter, item by item. No analysis changed and no reported value moved.

1. output/revision2/Supporting_Data_Values_v6.xlsx     [replaces v5]
   Six sheets added. The Data availability statement says the workbook holds the
   source data for every figure and table; before this it did not, because the
   analyses added during the revision had no sheet:

     S4_acid_base_panel     eighteen-gene acid-base and MCT-chaperone panel
     S5_mito_composites     forty-seven mitochondrial genes, astrocytes and neurons
     S8_staging_axes        the four neuropathological staging axes
     Table4_ROSMAP          the cross-cohort table
     S82_rosmap_matched     the ROSMAP Braak refit on the matched 380 donors
     SuppFig3_data          Supplemental Figure 3 is a schematic and has no
                            measurements of its own; this sheet lists every number
                            printed on it and the sheet each is taken from.

2. Three papers the reviewers cited by PMID were absent from the reference list
   and are now cited and discussed in the manuscript. They are recorded here
   because two of them bear directly on the acid-base analysis deposited in
   output/revision2/S4_acidbase_axis.csv:

     Stridh et al., J. Physiol. 590, 2333-2351 (2012)
       carbonic anhydrase II binds MCT1 and augments astrocytic lactate flux
       non-catalytically, acting as a proton-collecting antenna
     Forero-Quintero et al., J. Biol. Chem. 294, 593-607 (2019)
       membrane-anchored carbonic anhydrase IV reaches the transporters through
       the chaperones CD147 and GP70
     Zimmer et al., J. Neurochem. 169, e70294 (2025)
       the biphasic metabolic trajectory in which presymptomatic glial-driven
       hypermetabolism gives way to hypometabolism as tau spreads

   The manuscript reference list is now 82 entries, in the journal's style, and
   every entry has been resolved against PubMed or the published record.

3. One sentence in the manuscript stated a null result as a positive finding and
   used a temporal verb the reviewers had asked us to remove. It read that the
   protein-level Energy/Demand ratio did not differ across groups (p = 0.285)
   "confirming that quantitative protein abundance is preserved ... energetic
   vulnerability precedes structural loss". It now reports the null as a null.
   Nothing in this repository depends on that sentence; it is recorded so that
   the deposited record and the manuscript do not diverge.

---

AD-Metabolic-Collapse -- repository v17
=======================================

Changes from v16. These follow an external audit of the revision package. Every
item below is a defect the audit found in the deposited material, or one found
while checking the audit's findings. No analysis changed and no reported value
in the manuscript moved except one rounding correction, item 4.

1. output/revision2/Table3_verification.csv          [corrected]
   The "Correction 2" entry said that transferrin receptor "no longer exceeds
   its null" and that "verdict M-1 stands on ferritin light chain alone". That
   is the opposite of what the manuscript, the response, S12_VERDICTS.txt,
   S12_M_mediation_competition.csv and Supplemental Table 1c all record. The
   entry predates the reconciliation of the verdict criterion and is corrected
   in place, with a sentence saying that an earlier version said the opposite.

2. README.md                                          [corrected]
   Carried the same stale statement in the empirical-p paragraph. Corrected the
   same way. Found while checking item 1; the audit had not flagged it.

3. figures/Figure6.{png,tif}                          [replaced]
   R/figures/Fig6_Redraw.py                           [new]
   The deposited Figure 6 was an earlier version of the figure: panel B was
   MCT4-as-mediator with the reverse direction below, and panel C was a single
   detection-matched band (2,090 genes, MCT4 percentile 0.33). The printed
   figure has neither: panel B is the six-mediator competition across three
   neuronal outcomes, and panel C is three matched-null distributions (FTH1
   1,598 genes / 0.13; CP 2,378 / 0.21; MCT4 2,027 / 0.44). A reader running the
   deposited code would have obtained a figure the paper does not contain.

   Fig6_Redraw.py rebuilds the printed figure from the deposited inputs and
   asserts every published value before drawing: the three band sizes and
   percentiles, the three early-to-late changes, the mediated fractions and the
   95th percentiles for MCT4, TFRC and FTL at neuronal V-ATPase, and the 18 of 32
   sign agreements in panel D. R/figures/Fig6_CrossCohort_Mediation.py is marked
   SUPERSEDED and retained as the record of the earlier version.

   One defect in the printed figure is fixed at the same time. The panel B
   subtitle read "an asterisk marks an empirical p < 0.05", but transferrin
   receptor at neuronal V-ATPase carries a mark at p = 0.053. The verdict rule is
   not a p threshold; it is whether the observed value lies above the 95th
   percentile of that gene's own matched null. The subtitle now says that, which
   is what the marks have always meant. The manuscript legend is corrected to
   match.

4. Table 5 of the manuscript, CTSB late-bin mean
   Read 0.1693; recomputation from data/sample/astro_bin_means.csv gives
   0.169249, so 0.1692. Corrected in the manuscript. All twelve rows of Table 5
   were then recomputed from the deposited bin means and agree.

5. output/tables/couplings.csv                        [columns added]
   Written by P2_run_all.R, whose pcorZ uses n - 2 degrees of freedom, while the
   manuscript uses n - k - 2 throughout. The n - 2 columns are unchanged and the
   n - k - 2 values are added alongside as A_p_nk2, C1_p_nk2, C2_p_nk2, so the
   two conventions can be compared rather than one being quietly overwritten.

6. output/tables/partial_df_audit.csv                 [column renamed]
   The column p_in_manuscript held the value as submitted, which is no longer
   what the manuscript says. Renamed p_in_submitted_version. No value changed.

Not changed, and why:

   The ROSMAP rows of Table 4 are fitted on 380 donors (pseudotime) and 430
   (Braak), so the manuscript sentence "in the same donors ... changing the axis,
   and nothing else" was not literally true. The manuscript text and the Table 4
   note now give the donor count on each row and claim only that the discrepancy
   is confined to the pseudotime axis. Recomputing Braak on the matched 380
   donors requires the ROSMAP per-donor data, which is not redistributable and is
   not in this repository; the script that does it is
   82_rosmap_braak_matched_donors.R, run separately against the raw data.

7. output/revision2/Supplemental_Table_1_v13.xlsx     [replaces v12]
   Two new parts, both answering a reviewer request that had been answered only
   in the response letter and not in the paper:

     ST1i  acid-base and MCT-chaperone panel, eighteen genes, from
           S4_acidbase_axis.csv. Of the seventeen genes other than MCT4, none
           meets the joint condition MCT4 meets. CA4, CA14, CA9 and EMB are
           below the detection floor and are marked not measurable.
     ST1j  mitochondrial gene sets in astrocytes and excitatory neurons,
           forty-seven genes, from S5_mito_composites.csv. Every p in the table
           is uncorrected and the family-wise result is stated with it.

   The manuscript now carries both as Discussion paragraphs. They were
   previously in the response letter only, which meant a reviewer looking in the
   paper for the answer to a question they had asked would not have found it.

8. ROSMAP Braak refitted on the matched donor set             [new result]
   R/replication/82_rosmap_braak_matched_donors.R              [new]
   output/tables/82_rosmap_braak_matched.csv                   [new]
   output/tables/rosmap_donor_mct4.csv                         [new]
   tables/Table4_ROSMAP.csv, output/tables/Table4_ROSMAP.csv    [row added]
   figures/Figure6.{png,tif}                                    [redrawn]

   The manuscript claimed the two ROSMAP axes were compared "in the same donors
   and the same nuclei". They were not: pseudotime is defined for 380 of the 430
   donors carrying a Braak stage. The Braak slope was therefore refitted on
   exactly those 380 donors. It stays null:

     pseudotime, 380 donors                beta = +0.02142   p = 9.59e-08
     Braak, all 430 donors                 beta = -0.00070   p = 0.4596
     Braak, the same 380 donors            beta = -0.00080   p = 0.4315

   The claim is now literally true and Table 4 carries the matched row, which
   Figure 6A picks up automatically. The 50 donors that drop out are not a
   pathologically distinct group (Braak mean 3.48 against 3.63 in those kept).

   The script runs three gates before printing anything new: it reproduces both
   published ROSMAP values to 5e-4, checks that the pseudotime donors are a
   subset of the Braak donors, and checks that the restricted set is 380. It
   reads astrocytes.h5Seurat directly with rhdf5, reusing the reader in
   47_braak_common_axis.R, so no new dependency is introduced.

9. ROSMAP nuclei count in the manuscript                      [corrected]
   The Methods and Results paired "228,925 astrocytic nuclei" with "430 donors".
   228,925 is the size of the deposited file; the 430 donors carrying a Braak
   stage contribute 219,291 nuclei, and those are the nuclei that enter the
   donor-level analysis. Both figures are now given, with which is which. This
   is the same kind of error as the SEA-AD count corrected earlier in this
   revision, where the atlas total had been quoted in place of the analysed set.

---

AD-Metabolic-Collapse -- repository v16
=======================================

Changes from v15. No analysis changed and no reported number moved. The changes
are to the manuscript's section order and to the labels of the lettered parts of
Supplemental Table 1, and this repository is updated to match.

1. Methods moved to follow the Discussion, per this journal's Article format.
   References are numbered by first appearance, so the move renumbered 56 of the
   79 references and the reference list was reordered to match. No citation was
   added, removed or reassigned. Nothing in this repository depends on reference
   numbering, so no file here changes for that reason; it is recorded because it
   is the reason the part letters below moved.

2. The lettered parts of Supplemental Table 1 were relabelled so that they run in
   order of first appearance in the manuscript. Only four moved:

     staging axes             ST1h -> ST1e
     detection-matched null   ST1e -> ST1f
     sixteen-pathway screen   ST1f -> ST1g
     four-compartment split   ST1g -> ST1h

   Parts a to d are unchanged. Part a is the donor-level partial correlations and
   part b the Bin 0.1 sensitivity analysis; those are the two parts printed in the
   manuscript's Supplementary Material block, and their letters did not move.

3. output/revision2/Supplemental_Table_1_v12.xlsx        [replaces v11]
   Sheet names follow the new letters. No cell value changed.

4. CHANGES_v15.txt item 5 named two parts by the letters they carried at the time
   (ST1h for the staging axes, ST1e for the detection-matched null). Those letters
   have moved, so that entry now names the parts in words instead. The entry is
   corrected rather than deleted: it described real work, and a deposited record
   that silently drops a stale statement is worse than one that corrects it.

Unaffected by the relabelling, and checked:
  output/tables/SuppTable_1a_ANLS.csv        -- part a, unchanged
  output/tables/SuppTable_1b_recomputed.csv  -- part b, unchanged
  R/audit/78_supp_table_1b.R                 -- part b, unchanged

---

AD-Metabolic-Collapse -- repository v15
=======================================

Changes from v14. Nothing in the analysis layer changes; no reported number moves.

1. R/figures/Fig2D_Specificity_Control.R
   Facet titles no longer carry the internal model codes "[A]" and "[C1]".
   Those strings collided with the figure's own panel letters A and C, so a
   reader could read "[A]" inside panel D as panel A. The titles are now
   "CPS-adjusted" and "CPS + genome-wide control". The d$control column is
   unchanged, so output/tables/hk_control_by_gene.csv is unchanged.

2. R/figures/Fig2_Redraw.py                                            [new]
   The script of record for Figure 2 as printed. Reads the same two deposited
   inputs the R scripts use and asserts every published value before drawing:
   zero-order r +0.517; partial r +0.474 at p = 5.86e-06 on n - k - 2 degrees
   of freedom; the five panel B partials; specificity gaps 0.133 and 0.052.
   Fixes, in addition to item 1, four rendering defects: the clipped panel C
   x-axis label; the panel D signal label struck through by the mean line; the
   panel B tag overlapping the legend swatch, the legend having been moved
   inside the axes; and the panel tags, which are now placed in figure
   coordinates rather than axes coordinates so that the D tag begins at the
   same x as the A tag despite the panels differing in width.

   The panel tags are drawn with antialiasing off. The letterform, size and
   stroke width are identical to the ggplot2 tags in the other figures of this
   paper (95 x 91 px, median stroke 23 px at fontsize 30), but Agg lays a
   heavier grey edge than the original device: with antialiasing on the tag
   gained about 2% ink and one pixel of height, which reads as a bolder letter
   beside the ggplot2-rendered tags. With it off the total ink is within 0.6%
   of the original. Only the tags are affected; all other text keeps
   antialiasing.
   seed 42, as everywhere else in this repository.

   The R scripts are retained. They are not byte-identical to the published
   rendering -- panel D jitter comes from a different generator and font
   metrics differ between ggplot2 and matplotlib -- but every plotted value is
   the same.

3. figures/Figure2.png, figures/Figure2.tif
   Replaced with the corrected rendering. 4500 x 3099 px, 300 dpi, RGB, LZW.

4. R/revision2/91_revision2_S12_all.R
   The deposited S12_M_mediation_competition.csv carries two columns,
   p_empirical and p_empirical_uncorrected, but the script wrote only one. The
   script now computes both, so code and output agree:

     p_empirical_uncorrected = k / n
     p_empirical             = (k + 1) / (n + 1)

   The verdict column follows k / n, which is the form in which rule M-1 was
   pre-registered. The (k + 1) / (n + 1) correction was adopted at revision
   because a permutation p of exactly zero is not attainable, and it is the
   value reported in the paper; it is not used to reverse a verdict. Of the
   eighteen mediator-outcome tests it moves exactly one across 0.05 --
   transferrin receptor at neuronal V-ATPase, 0.0494 to 0.0530 -- and in the
   direction that favours this study's own claim. That is stated in the paper,
   in Supplemental Table 1c and here.

5. output/revision2/Supplemental_Table_1_v11.xlsx        [replaces v10]
   New part for the four neuropathological staging axes (ADNC, Thal, CERAD,
   Braak): astrocytic MCT4 with the neuronal identity control on the same donors
   and the same cells, the donor-count gate and the provenance of the candidate
   set. Built from output/revision2/S8_strata.csv and S8_axis_counts.csv. Notes
   updated in the mediation-competition part (footnotes on the marginal TFRC
   call and on the p convention) and in the detection-matched-null part
   (percentile convention).

   NOTE: the part letters cited here were reassigned in v16; see CHANGES_v16.txt.
   The lettered parts are named in words above so that this record stays true.
