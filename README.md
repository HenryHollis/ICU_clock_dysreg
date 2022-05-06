# ICU_clock_dysreg
This is R 4.1.2 code used in the analysis of circadian clock dysregulation among people undergoing critical care.
<ul>
  <li>Differential expression calculated with "run_all_diff_expr_edgeR.R"</li>
  <li>Genes list such as those that are ubiquitously differentially expressed across tissues, etc., are created with Generate_FDR_LFC_tables.R</li>
  <li>Modulated genes between AD-ICU enriched using enrichR_script</li>
  <li>CCD and dCCD calculated with PermTissueTesting</li>
  <li> Note: in getPatientStats.R, prostate and testis only have male subjects and as such, R shifts the column values one to the left for those tissues</li>
</ul>
