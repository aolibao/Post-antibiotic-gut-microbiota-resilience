#!/home/qqt/anaconda3/bin/python

from pathlib import Path
import shutil
import re

ROOT = Path("/home/qqt/20061009response")
WORK = ROOT / "reviewer_response" / "manuscript_revision"


def replace_between(text: str, start: str, end: str, replacement: str) -> str:
    i = text.find(start)
    if i < 0:
        raise RuntimeError(f"Start marker not found: {start}")
    j = text.find(end, i + len(start))
    if j < 0:
        raise RuntimeError(f"End marker not found: {end}")
    return text[:i] + replacement.rstrip() + "\n\n" + text[j:]


manuscript = (WORK / "Manuscript_original_extracted.md").read_text(encoding="utf-8")
manuscript = manuscript.replace(
    "Post-antibiotic gut microbiota resilience is determined by baseline\n"
    "community structure and enhanced by tryptophan supplementation",
    "Baseline community structure is associated with post-antibiotic gut microbiota recovery, with an early tryptophan-associated ecological signal in mice",
)
manuscript = manuscript.replace(
    "**Gut microbiome at baseline drives recovery capacity following\n"
    "ciprofloxacin exposure**",
    "**Baseline gut microbiome is associated with recovery patterns following ciprofloxacin exposure**",
)
manuscript = manuscript.replace(
    "**Fig. 2. Gut microbiome at baseline drives recovery capacity following\n"
    "ciprofloxacin exposure**",
    "**Fig. 2. Baseline gut microbiome is associated with recovery patterns following ciprofloxacin exposure**",
)
manuscript = manuscript.replace(
    "recovery follows an orderly,\nnon-random process determined by baseline microbiota characteristics.",
    "recovery follows a structured, non-random process associated with baseline microbiota characteristics.",
)
manuscript = manuscript.replace(
    "post-exposure and drives gut microbial diversity loss",
    "post-exposure and coincides with gut microbial diversity loss",
)
manuscript = manuscript.replace(
    "recovery (B6-C6).",
    "recovery (B6-C6). The original cohort study reported approval under Stanford IRB protocol #25268 and written informed consent. Age was available for 45 of 60 exposed participants and sex for 58 of 60; BMI, diet, bowel preparation and concomitant-medication variables were not available in the supplied public metadata. The institutional determination for the present public-data secondary analysis will be confirmed by the authors before submission.",
    1,
)

manuscript = replace_between(
    manuscript,
    "**Abstract**",
    "**Key words**",
    """**Abstract**

Antibiotic exposure perturbs the gut microbiota, but recovery varies markedly among individuals. We analysed 60 ciprofloxacin-exposed participants with dense longitudinal sampling and five unexposed participants as a descriptive reference. Recovery was evaluated with endpoint diversity and six within-participant measures referenced to each participant's own baseline. The principal ecological patterns were consistent across these definitions: participants with stronger recovery maintained higher diversity during peak disturbance and showed more connected late-recovery ecological networks, whereas persistent *Bacteroides*-dominated configurations were associated with lower diversity. Dirichlet multinomial mixture modelling showed that the recovery phenotype was not a simple relabelling of baseline enterotype. Human shotgun profiles further identified NAD-related inferred microbial gene-content patterns, which are interpreted as associative functional potential rather than measured pathway activity or flux. In a mouse experiment assessed by 16S rRNA gene sequencing, tryptophan pretreatment showed a consistent directional signal toward earlier diversity recovery, with the largest effect during early recovery and convergence at the late time point. The mouse data do not provide metagenomic or metabolomic evidence for a tryptophan-to-NAD mechanism. Together, these findings define reproducible ecological features of post-antibiotic recovery and identify tryptophan-associated early community recovery as a candidate preclinical observation requiring mechanistic validation.""",
)

manuscript = replace_between(
    manuscript,
    "This study used a large scale longitudinal human cohort",
    "**Methods**",
    """This study used a densely sampled longitudinal human cohort to characterize inter-individual heterogeneity in gut microbiota recovery following ciprofloxacin exposure. We evaluated recovery using both the original endpoint-diversity grouping and multiple within-participant displacement metrics, and used DMM modelling to distinguish recovery phenotype from baseline community type. Human shotgun functional profiles were interpreted as inferred microbial gene-content potential. The mouse arm used 16S rRNA gene sequencing to test whether nutritional interventions were associated with post-antibiotic ecological recovery; it was not designed to provide direct metagenomic, metabolomic or causal validation of the NAD-related hypothesis.""",
)

manuscript = replace_between(
    manuscript,
    "**Taxonomic Profiling**",
    "**Functional Analysis**",
    """**Taxonomic Profiling**

Raw reads were inspected using FastQC 0.12.1 and trimmed with Trimmomatic 0.40. Host-read removal was performed with KneadData package 0.12.1 (the installed command-line banner reports v0.12.0) using Bowtie2 2.5.4 against the local GRCh38.p14 index (NCBI GCF_000001405.40); reads aligning to the human reference were discarded. Taxonomic classification was performed with Kraken2 2.1.1 against the downloaded NCBI RefSeq Kraken2 database version dated 28 December 2024, containing 14,931 bacterial reference genomes, followed by genus- and species-level abundance estimation with Bracken 2.5.0. The database was downloaded as a released reference resource and was not constructed locally. The primary taxonomic analysis rarefied bacterial species-level counts to 2,693,000 reads per sample with random seed 20261009. An explicitly labelled unrarefied branch was retained as a sequencing-depth sensitivity analysis.""",
)

manuscript = replace_between(
    manuscript,
    "**Functional Analysis**",
    "**Microbial Diversity Analysis**",
    """**Functional Analysis**

Clean paired and unpaired reads from KneadData were concatenated by sample and profiled with HUMAnN 3.9. Taxonomic prescreening used MetaPhlAn 4.0.6 with the retained `mpa_vJun23_CHOCOPhlAnSGB_202403` marker-database snapshot. HUMAnN nucleotide alignment used the ChocoPhlAn pangenome database `v201901_v31`, and translated search used UniRef90 `201901b` with the installed DIAMOND executable (CLI banner 2.0.15). The configured prescreen threshold was 0.01; nucleotide and translated query-coverage thresholds were 90%, subject-coverage thresholds were 50%, and the translated-search E-value threshold was 1.0. Stratified output was retained. Gene-family abundances were regrouped to KEGG Orthology and enzyme commission categories and normalized to copies per million. HUMAnN outputs are described throughout as inferred microbial gene-content potential; they do not measure pathway activity, NAD concentration, metabolite abundance or metabolic flux.""",
)

manuscript = replace_between(
    manuscript,
    "**Microbial Diversity Analysis**",
    "**Microbial community dynamic modeling using MDSINE2**",
    """**Microbial Diversity Analysis**

Alpha diversity was quantified using Shannon diversity and observed richness. Bray-Curtis and Jaccard dissimilarities were used to assess beta diversity, and PERMANOVA with 999 permutations was used for prespecified community-structure contrasts. Analyses were conducted with phyloseq 1.44.0.

**Recovery definitions and full-cohort sensitivity analyses**

The original endpoint-Shannon HDR/LDR grouping was retained for comparison. Recovery was additionally defined by six within-participant measures referenced to the mean of baseline visits A1-A4: endpoint Shannon deviation, endpoint Bray-Curtis distance, endpoint Jaccard distance, and time-normalized displacement areas under the curve for the same three measures. For each definition, the 15 participants with the smallest displacement were assigned to the high-recovery group and the 15 with the largest displacement to the low-recovery group. The full cohort of 60 exposed participants was also analysed without dichotomization using continuous recovery outcomes, baseline Shannon diversity as a continuous predictor, and heteroscedasticity-consistent HC3 confidence intervals. Cut-off sensitivity analyses used 12, 15 and 20 participants per extreme group.

**Dirichlet multinomial mixture analysis**

Baseline community types were estimated from participant-level A1-A4 genus count profiles using DirichletMultinomial 1.44.0 in R 4.3.1. The audit matrix retained 2,173 detected genus-level features. For stable mixture estimation, near-zero components with cohort mean relative abundance below 0.01% were excluded, leaving 143 genera, and profiles were standardized to 100,000 multinomial counts per participant. Models with two to six components were compared using the Laplace criterion. Associations between selected DMM type and the original HDR/LDR grouping were assessed by Fisher's exact test, and continuous recovery was modelled jointly with baseline Shannon diversity using HC3 inference.""",
)

manuscript = replace_between(
    manuscript,
    "**Microbial community dynamic modeling using MDSINE2**",
    "**Microbial Network Analysis**",
    """**Microbial community dynamic modelling using MDSINE2**

MDSINE2 was used as a rapid sensitivity analysis of relative-abundance dynamics under ciprofloxacin exposure. The model followed the generalized Lotka-Volterra framework with microbial interaction modules and perturbation effects. Because matched qPCR measurements were unavailable and the sensitivity runs used a shortened single-chain configuration, the resulting interaction and perturbation summaries are interpreted descriptively on the relative-abundance scale and are not presented as evidence of formal MCMC convergence or absolute-abundance dynamics.""",
)

manuscript = replace_between(
    manuscript,
    "**Microbial Network Analysis**",
    "***In vivo* mouse models**",
    """**Microbial Network Analysis**

FastSpar was used for SparCC correlation analysis of genus-level community profiles with 500 bootstrap iterations. Multiple-testing-adjusted *P* values below 0.05 defined retained edges, and networks were visualized using ggClusterNet 2.0. To evaluate mathematical coupling between taxon relative abundance and Shannon diversity, genus abundances were additionally CLR-transformed and Shannon diversity was recalculated after excluding all species assigned to the focal genus and renormalizing the remaining community. Participant-label permutations and a single global Benjamini-Hochberg correction across 124 prevalence-filtered named genera and 16 visits were used. Longitudinal mixed models separated between-participant and within-participant CLR effects and adjusted for visit and sequencing depth.""",
)

manuscript = replace_between(
    manuscript,
    "***In vivo* mouse models**",
    "**Quantitative RT‒PCR**",
    """***In vivo* mouse models**

Six-week-old C57BL/6 mice received tryptophan (2 g/kg/day), nicotinamide (160 mg/kg/day), or regular drinking water for seven days. Mice in the intervention groups then either continued supplementation or received supplementation together with ciprofloxacin (12.5 mg/kg/day) for 72 hours; the antibiotic-control group received water followed by ciprofloxacin alone. The animal experiments were approved by the Experimental Animal Ethics Committee of Jinan University (IACUC-20251216-09).

**16S rRNA gene sequencing and ecological analysis**

The mouse arm generated full-length bacterial 16S rRNA gene amplicon data only. Libraries were sequenced on an Oxford Nanopore platform, processed using QIIME 2 2020.6, and assigned taxonomy with a naive Bayes classifier trained on SILVA 132. Contaminant ASVs were removed before ecological analysis. These data support taxonomic composition, alpha-diversity and longitudinal community-recovery analyses, but they do not directly measure metagenomic KO abundance, metabolic pathways, metabolites or flux. No mouse shotgun metagenomic analysis is reported. Animal-level Shannon diversity and evenness were analysed with prespecified Welch contrasts, 95% confidence intervals, Hedges' *g*, and Benjamini-Hochberg correction within each contrast family.""",
)

manuscript = replace_between(
    manuscript,
    "**Statistics Analysis**",
    "**Result**",
    """**Statistical Analysis**

The participant or individual animal was the independent sampling unit. Continuous group contrasts used Welch tests unless a paired analysis was specified; effect estimates are reported with 95% confidence intervals and standardized effect sizes. Distance-matrix contrasts used PERMANOVA with 999 permutations. Correlation analyses used Spearman coefficients, participant-label permutations where indicated, and Benjamini-Hochberg false-discovery-rate correction over the defined testing family. Full-cohort linear models used HC3 heteroscedasticity-consistent inference. Exact unadjusted and adjusted *P* values are reported; adjusted *q* values below 0.05 were considered statistically significant.""",
)

manuscript = replace_between(
    manuscript,
    "**Tryptophan pretreatment facilitates the recovery of gut microbiota",
    "**Discussion**",
    """**Sensitivity analyses distinguish recovery phenotype from baseline enterotype**

Across six within-participant recovery definitions, the principal diversity and network trends were directionally consistent with the original endpoint-Shannon grouping: high-recovery participants maintained higher diversity during peak disturbance and showed more connected late-recovery networks. Individual taxa and MDSINE2 parameters varied among definitions and are therefore interpreted as sensitivity results rather than exact replications. In the full cohort, baseline Shannon diversity did not clearly predict continuous Bray-Curtis recovery (standardized beta=0.153, 95% CI -0.103 to 0.409, *P*=0.236). DMM modelling selected two baseline community types characterized by *Bacteroides* and *Segatella*. DMM type was not associated with the original HDR/LDR classification (Fisher's exact *P*=0.215) and did not explain continuous day-5 recovery after adjustment for baseline Shannon diversity (HC3 *P*=0.497), showing that the recovery phenotype was not a simple relabelling of baseline enterotype.

**Tryptophan pretreatment is associated with an early ecological recovery signal after ciprofloxacin exposure**

The mouse 16S data were used to evaluate ecological recovery rather than metagenomic function. The estimated Shannon difference between ABX-Trp and ABX-Trp-During was positive at days 7, 14 and 28 (+0.402, +0.372 and +0.0225), with the largest standardized effect at day 14 (Hedges' *g*=1.58). The direct timing family did not reach statistical significance after Benjamini-Hochberg correction (*q*=0.454), and the confidence intervals crossed zero. Compared with the antibiotic-control group at day 14, ABX-Trp showed a +0.769 Shannon difference and a large standardized effect (*g*=2.01; BH *q*=0.103). At day 28, the ABX-Trp versus control evenness difference was small (+0.00449; *g*=0.51; BH *q*=0.547), consistent with convergence during later recovery. We therefore interpret the mouse result as a directional, time-localized ecological signal of earlier community reorganization, not as statistically established superiority of pretreatment or evidence of a tryptophan-to-NAD mechanism (Fig. 6).""",
)

manuscript = replace_between(
    manuscript,
    "Through integrated metagenomic and functional metabolic analyses",
    "Previous studies on post-antibiotic resilience",
    """Human shotgun metagenomic profiles identified differences in NAD-related inferred microbial gene-content potential between recovery groups. These HUMAnN-derived features are carried by the taxa that compose each community and do not measure pathway activity, NAD concentration or metabolic flux. After adjustment for baseline taxonomic composition, richness and read depth, neither NAD salvage pathway II nor NAD de novo biosynthesis I was independently associated with the primary continuous recovery outcome (BH *q*=0.977 for both). We therefore interpret the NAD-related pattern as an ecological correlate embedded within community composition and as a source of hypotheses, rather than as an established mediator.

The mouse experiment provides a complementary ecological observation but not functional validation. Only 16S rRNA gene amplicon data were generated, so KO abundance and NAD pathway activity cannot be recovered from this experiment. Within the resolution of the 16S data, tryptophan pretreatment showed a consistent directional association with higher Shannon diversity during early recovery, while the groups converged at the later time point. This temporal pattern is compatible with accelerated early community reorganization. Tryptophan could act through several non-exclusive routes, including altered substrate availability, indole/AhR signalling, host-mediated effects or community restructuring; distinguishing these possibilities will require new shotgun metagenomic, targeted metabolomic or functional measurements.""",
)

manuscript = replace_between(
    manuscript,
    "Previous studies on post-antibiotic resilience have focused",
    "**Conclusion**",
    """Previous studies on post-antibiotic resilience have focused on host genetics, diet, probiotics and faecal microbiota transplantation. The present study adds a densely sampled ecological analysis of 60 antibiotic-exposed participants and shows that recovery-associated diversity and network patterns remain directionally consistent across multiple within-participant definitions. The persistent *Bacteroides*-associated configuration and the more even, connected recovery pattern may help identify ecological states for prospective validation.

The mouse 16S experiment extends this ecological framework by suggesting that intervention timing may influence early community reorganization. Because the adjusted timing contrasts were not statistically significant and no mouse metagenomic or metabolomic data were generated, tryptophan is presented as a candidate intervention for prospectively powered study rather than as an established therapeutic or NAD-directed mechanism.""",
)

manuscript = replace_between(
    manuscript,
    "**Conclusion**",
    "**Study funding**",
    """**Conclusion**

Baseline community structure was associated with heterogeneous ecological recovery after ciprofloxacin exposure, and the principal diversity and network patterns were consistent across multiple within-participant recovery definitions. DMM analysis showed that this recovery phenotype was not reducible to baseline enterotype. Human NAD-related profiles represent inferred microbial functional potential and should not be interpreted as direct pathway activity. In the mouse 16S experiment, tryptophan pretreatment showed a directional early-recovery signal with later convergence, supporting further mechanistic study without establishing a tryptophan-to-NAD causal pathway. These findings provide an ecological framework for studying post-antibiotic recovery and for designing prospectively powered intervention experiments.""",
)

manuscript = replace_between(
    manuscript,
    "**Data availability**",
    "**References**",
    """**Data availability**

The mouse 16S rRNA gene sequencing data are publicly available under BioProjects PRJNA1472735 and PRJNA1470252. The public human shotgun metagenomic cohort is available under BioProject PRJNA974858. The mouse experiment did not generate shotgun metagenomic data. The final code repository and persistent identifier will be added after deposition.""",
)

manuscript = replace_between(
    manuscript,
    "**Fig. 6. Tryptophan pretreatment facilitates the recovery of gut",
    "",  # handled below
    "",
) if False else manuscript

fig6_start = manuscript.find("**Fig. 6. Tryptophan pretreatment facilitates the recovery of gut")
if fig6_start < 0:
    raise RuntimeError("Fig. 6 legend marker not found")
manuscript = manuscript[:fig6_start] + """**Fig. 6. Animal-level 16S ecological analysis of treatment timing and recovery**

**(A)** Animal-level Shannon diversity at recovery days 7, 14 and 28; lines and points summarize group means with 95% confidence intervals.

**(B)** Raw animal-level Shannon observations at day 14.

**(C)** Direct ABX-Trp versus ABX-Trp-During Shannon contrasts at recovery days 7, 14 and 28; positive estimates favor pretreatment and bars show 95% confidence intervals.

**(D)** Day-28 evenness contrasts between each intervention group and the antibiotic-control group; bars show 95% confidence intervals. Benjamini-Hochberg-adjusted results are reported in the source analysis table. The mouse arm used 16S rRNA gene sequencing only; no metagenomic pathway or KO panels are presented."""

manuscript_path = WORK / "Manuscript_revised_true_evidence.md"
manuscript_path.write_text(manuscript, encoding="utf-8")


supplement = (WORK / "Supplement_original_extracted.md").read_text(encoding="utf-8")
supplement = supplement.replace(
    "**Post-antibiotic gut microbiota resilience is determined by baseline\n"
    "community structure and enhanced by tryptophan supplementation**",
    "**Baseline community structure is associated with post-antibiotic gut microbiota recovery, with an early tryptophan-associated ecological signal in mice**",
)
supplement = supplement.replace(
    "sFig. 4. Baseline NAD metabolic bottlenecks associated with de novo\n"
    "biosynthesis",
    "sFig. 4. Baseline NAD-related inferred microbial functional potential",
)
supplement = supplement.replace(
    "**Supplementary Fig. 4. Baseline NAD metabolic bottlenecks associated\n"
    "with de novo biosynthesis**",
    "**Supplementary Fig. 4. Baseline NAD-related inferred microbial functional potential**",
)
supplement = supplement.replace(
    "**(A) **UpSet plots show the overlap of genera across the baseline,\n"
    "exposure, and recovery phases. We included genera that had a significant\n"
    "positive correlation with the Shannon index at one or more time points\n"
    "within each phase (*P* \\< 0.05 and r \\> 0). The bars show the\n"
    "intersection sizes for the LDR (left) and HDR (right) groups.",
    "**(A)** UpSet plots show intersections of quality-controlled Bracken genus-level features with positive Shannon correlations across baseline, exposure and recovery phases (*P* < 0.05 and rho > 0). The values 862 and 1,294 are feature-correlation intersection counts in the corresponding LDR and HDR displays, respectively. All quality-controlled genus-level features, including low-abundance and provisionally named features, were retained to avoid omitting antibiotic-sensitive organisms; these values are not counts of conventional named human-gut genera.",
)
supplement = supplement.replace(
    "**(B)** Line graphs illustrate the individual recovery-phase\n"
    "trajectories of Shannon index deviations from the baseline. The baseline\n"
    "is defined as the mean Shannon index across A1–A4, and the y-axis\n"
    "represents the Shannon index at a specific time point (B6–C6) minus the\n"
    "baseline mean. The dashed line denotes the baseline level (0).",
    "**(B)** Line graphs illustrate individual recovery-phase trajectories of Shannon-index deviation from baseline. Baseline is the mean Shannon index across A1–A4, and the y-axis is the Shannon index at B6–C6 minus that baseline mean; the dashed line denotes zero deviation. The trajectory classes are Repeat-Recover, Dynamic Non-Recovery and Dynamic Sustained-Recovery. For the pathway contrast shown in Fig. 4A, Recover denotes Repeat-Recover (n=9) and Non-recover denotes Dynamic Non-Recovery (n=31); Dynamic Sustained-Recovery was not included in that binary comparison.",
)

source_mouse = ROOT / "reviewer_response" / "R1.9_mouse_reanalysis" / "figures" / "R1_9_mouse_alpha_reanalysis_word.png"
target_mouse = WORK / "source_media" / "supplement" / "media" / "sFigure5_revised_16S.png"
shutil.copy2(source_mouse, target_mouse)

s5_img = supplement.find('<img\nsrc="reviewer_response/manuscript_revision/source_media/supplement/media/image5.tiff"')
if s5_img < 0:
    raise RuntimeError("Supplementary Fig. 5 image marker not found")
supplement = supplement[:s5_img] + """<img src="reviewer_response/manuscript_revision/source_media/supplement/media/sFigure5_revised_16S.png" style="width:6.5in" />

**Supplementary Fig. 5. Animal-level 16S ecological recovery after nutritional intervention**

**(A)** Animal-level Shannon diversity trajectories during recovery.

**(B)** Raw animal-level Shannon observations at day 14.

**(C)** Direct ABX-Trp versus ABX-Trp-During contrasts across recovery days.

**(D)** Day-28 evenness contrasts between intervention groups and the antibiotic-control group. The mouse experiment generated 16S rRNA gene amplicon data only; no predicted or measured KO abundance is presented."""

# Pandoc's DOCX extraction uses raw HTML image elements. Convert all of them
# to ordinary Markdown images so both DOCX and XeLaTeX outputs embed the figures.
supplement = re.sub(
    r'<img\s+src="([^"]+)"[^>]*/>',
    r'![](\1)',
    supplement,
    flags=re.S,
)

(WORK / "Supplement_revised_true_evidence.md").write_text(supplement, encoding="utf-8")

print(manuscript_path)
print(WORK / "Supplement_revised_true_evidence.md")
