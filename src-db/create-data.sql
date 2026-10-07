--Based on assembly GRCh38.p14
--https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
INSERT INTO prada.chromosome(number, sizebp) VALUES (1,248956422);
INSERT INTO prada.chromosome(number, sizebp) VALUES (2,242193529);
INSERT INTO prada.chromosome(number, sizebp) VALUES (3,198295559);
INSERT INTO prada.chromosome(number, sizebp) VALUES (4,190214555);
INSERT INTO prada.chromosome(number, sizebp) VALUES (5,181538259);
INSERT INTO prada.chromosome(number, sizebp) VALUES (6,170805979);
INSERT INTO prada.chromosome(number, sizebp) VALUES (7,159345973);
INSERT INTO prada.chromosome(number, sizebp) VALUES (8,145138636);
INSERT INTO prada.chromosome(number, sizebp) VALUES (9,138394717);
INSERT INTO prada.chromosome(number, sizebp) VALUES (10,133797422);
INSERT INTO prada.chromosome(number, sizebp) VALUES (11,135086622);
INSERT INTO prada.chromosome(number, sizebp) VALUES (12,133275309);
INSERT INTO prada.chromosome(number, sizebp) VALUES (13,114364328);
INSERT INTO prada.chromosome(number, sizebp) VALUES (14,107043718);
INSERT INTO prada.chromosome(number, sizebp) VALUES (15,101991189);
INSERT INTO prada.chromosome(number, sizebp) VALUES (16,90338345);
INSERT INTO prada.chromosome(number, sizebp) VALUES (17,83257441);
INSERT INTO prada.chromosome(number, sizebp) VALUES (18,80373285);
INSERT INTO prada.chromosome(number, sizebp) VALUES (19,58617616);
INSERT INTO prada.chromosome(number, sizebp) VALUES (20,64444167);
INSERT INTO prada.chromosome(number, sizebp) VALUES (21,46709983);
INSERT INTO prada.chromosome(number, sizebp) VALUES (22,50818468);
INSERT INTO prada.chromosome(number, sizebp) VALUES (23,156040895); --X
INSERT INTO prada.chromosome(number, sizebp) VALUES (24,57227415); --Y
--INSERT INTO prada.chromosome(number, sizebp) VALUES (25,NULL); --XY
--INSERT INTO prada.chromosome(number, sizebp) VALUES (26,NULL); --MT

WITH prev_chrom AS (SELECT CASE
WHEN c.number <23 THEN 'chr'||c.number
WHEN c.number =23 THEN 'chrX'
WHEN c.number =24 THEN 'chrY'
WHEN c.number =25 THEN 'chrXY'
WHEN c.number =26 THEN 'chrMT'
END AS name, c.number
FROM prada.chromosome c
)
UPDATE prada.chromosome SET name=p.name FROM prev_chrom p WHERE chromosome.number=p.number;

--Harmonise the pgx_gene custom data with gencode gene codes
UPDATE prada.pgx_gene pg SET "Gene"='GBA1' WHERE "Gene"='GBA';


--insert drug data
--pharmgkb cancer - cancer_drugs_pharmgkb.tsv
INSERT INTO prada.drug(name,type,weight) SELECT cdp.drug, 'cancer',0.5 FROM prada.cancer_drugs_pharmgkb cdp;
--Set initial set of drug selections
UPDATE prada.drug d SET "selected_for_analysis"=1 WHERE "name"='amitriptyline';
UPDATE prada.drug d SET "selected_for_analysis"=1 WHERE "name"='citalopram';
UPDATE prada.drug d SET "selected_for_analysis"=1 WHERE "name"='escitalopram';
UPDATE prada.drug d SET "selected_for_analysis"=1 WHERE "name"='fluvoxamine';
UPDATE prada.drug d SET "selected_for_analysis"=1 WHERE "name"='paroxetine';
UPDATE prada.drug d SET "selected_for_analysis"=1 WHERE "name"='sertraline';
UPDATE prada.drug d SET "selected_for_analysis"=1 WHERE "name"='venlafaxine';
UPDATE prada.drug d SET "selected_for_analysis"=1 WHERE "name"='vortioxetine';

--insert IDP GWAS data
INSERT INTO prada.variant(type,snp,chr,bp,bp2,mdd_p,mdd_beta,mdd_beta_se,mdd_beta_n) SELECT d.type, d.snp, d.chr, d.bp, d.bp2, d.mdd_p, d.mdd_beta, d.mdd_beta_se, d.mdd_n FROM prada.idp38_import d;

--insert MDD CNV Data
INSERT INTO prada.variant(type,snp,chr,bp,bp2,mdd_p) SELECT d.type, d.snp, d.chr, d.bp, d.bp2, d.mdd_p FROM prada."mmddcnv0index.b38.bed" d;

--insert PRADA drug recommendations
--version 2 format, excel sheet by giuseppe and danyang
--version 2: improved data import. correct import of floating point numbers of start and target dose for example.
DROP TABLE IF EXISTS t_drug_recommendations_import;
CREATE TEMP TABLE IF NOT EXISTS t_drug_recommendations_import AS
SELECT --over diplotypes, should have one row per gene
row_number() OVER (PARTITION BY hcpgx.recommendation,hcpgx.guidelineid ,hcpgx.drug_name,hcpgx.gene_name ORDER BY hcpgx.prada_ehrpriority_num DESC, hcpgx.prada_cpiclevel_num DESC, hcpgx.prada_pgkbcalevel_num DESC, hcpgx.diplotype ASC) rn,
hcpgx.recommendation,
hcpgx.guidelineid,
hcpgx.drug_name,
hcpgx.implications,
hcpgx.gene_name,
d.start_dose,
d.target_dose,
CASE --prada_titration_speed
	WHEN d.titration_speed = '0' THEN 0
	WHEN d.titration_speed = '1' THEN 2
	ELSE 1
END titration_speed,
CASE 
	WHEN strpos(d.switch_as_1st_choice,hcpgx.drug_name)!=0 THEN array_append(ARRAY[]::text[],hcpgx.drug_name)
	WHEN strpos(d.switch_as_2nd_choice,hcpgx.drug_name)!=0 THEN array_append(ARRAY[]::text[],hcpgx.drug_name)
	WHEN strpos(d.switch_as_1st_choice,'tricyclic')!=0 THEN array_append(ARRAY[]::text[],'tricyclic')
	WHEN strpos(d.switch_as_2nd_choice,'tricyclic')!=0 THEN array_append(ARRAY[]::text[],'tricyclic')
	WHEN strpos(d.switch_as_1st_choice,'tertiary amines')!=0 THEN array_append(ARRAY[]::text[],'tertiary amines') --there are some positive recommendations for these rules also, encode?
	WHEN strpos(d.switch_as_2nd_choice,'tertiary amines')!=0 THEN array_append(ARRAY[]::text[],'tertiary amines')
	ELSE ARRAY[]::text[]
END prada_switch_drug,
CASE WHEN strpos(d.switch_as_1st_choice,hcpgx.drug_name)!=0 THEN 1 ELSE 0 END prada_switch1_drug, --prada_switch1_drug
CASE WHEN strpos(d.switch_as_1st_choice,hcpgx.gene_name)!=0 THEN 1 ELSE 0 END prada_switch1_gene, --prada_switch1_gene
CASE WHEN strpos(d.switch_as_2nd_choice,hcpgx.drug_name)!=0 THEN 1 ELSE 0 END prada_switch2_drug, --prada_switch2_drug
CASE WHEN strpos(d.switch_as_2nd_choice,hcpgx.gene_name)!=0 THEN 1 ELSE 0 END prada_switch2_gene, --prada_switch2_gene
d.tdm,
hcpgx.drugrecommendation
FROM postgres.prada_cpic_pgx_variables_coded_20260925 d
LEFT OUTER JOIN prada.harmonised_combined_pgx hcpgx ON hcpgx.recommendation =d.recommendation AND d.drug_name  =hcpgx.drug_name AND d.guidelineid = hcpgx.guidelineid AND d.implications::jsonb = hcpgx.implications;
CREATE UNIQUE INDEX t_drug_recommendations_import_u ON t_drug_recommendations_import (recommendation,guidelineid ,drug_name,gene_name,rn);
--SELECT * FROM t_drug_recommendations_import;

--Import data first into t_drug_recommendations_import
--DELETE FROM prada.recommendation *;  
INSERT INTO prada.recommendation(
recommendation,
guideline,
drugid,
implications,
gene_name,
prada_start_dose,
prada_target_dose,
prada_titration_speed,
prada_switch_drug,
prada_switch1_drug,
prada_switch1_gene,
prada_switch2_drug,
prada_switch2_gene,
prada_tdm,
drugrecommendation,
version
)
SELECT 
m.recommendation,
m.guidelineid,
m.drug_name,
m.implications,
m.gene_name,
m.start_dose,
m.target_dose,
m.titration_speed,
m.prada_switch_drug,
m.prada_switch1_drug,
m.prada_switch1_gene,
m.prada_switch2_drug,
m.prada_switch2_gene,
m.tdm,
m.drugrecommendation,
2 --version
FROM t_drug_recommendations_import m
--WHERE hcpgx.recommendation=5991110
WHERE m.recommendation IS NOT NULL AND m.rn =1
--ORDER BY m.recommendation,m.guidelineid ,m.drug_name
;





--Experiments with transferring recommendations across CPIC database versions

----SELECT * FROM prada.recommendation rec WHERE rec.recommendation = 5991110;
--SELECT
--rec.recommendation,
--rec.guideline,
--rec.drugid,
--rec.implications::text,
--rec.gene_name,
--rec.prada_start_dose,
--rec.prada_target_dose,
--rec.prada_titration_speed,
--rec.prada_switch_drug[1],
--rec.prada_switch1_drug,
--rec.prada_switch1_gene,
--rec.prada_switch2_drug,
--rec.prada_switch2_gene,
--rec.prada_tdm,
--rec.drugrecommendation
--FROM prada.recommendation rec;


----version 1 format 2, export across cpic versions - not ready yet (does not work)!!!
--DROP TABLE IF EXISTS t_drug_recommendations_import;
--CREATE TEMP TABLE IF NOT EXISTS t_drug_recommendations_import AS
--SELECT --over diplotypes, should have one row per gene
--row_number() OVER (PARTITION BY hcpgx.recommendation,hcpgx.guidelineid ,hcpgx.drug_name,hcpgx.gene_name ORDER BY hcpgx.prada_ehrpriority_num DESC, hcpgx.prada_cpiclevel_num DESC, hcpgx.prada_pgkbcalevel_num DESC, hcpgx.diplotype ASC) rn,
--hcpgx.recommendation,
--d.recommendation recommendation_old,
--hcpgx.guidelineid,
--hcpgx.drug_name,
--hcpgx.implications,
--hcpgx.gene_name,
--d.prada_start_dose,
--d.prada_target_dose,
--d.prada_titration_speed,
--d.prada_switch_drug,
--d.prada_switch1_drug, --prada_switch1_drug
--d.prada_switch1_gene, --prada_switch1_gene
--d.prada_switch2_drug, --prada_switch2_drug
--d.prada_switch2_gene, --prada_switch2_gene
--d.prada_tdm,
--hcpgx.drugrecommendation
--FROM postgres.recommendation d -- we can't use the original recommendation keys
--LEFT OUTER JOIN prada.harmonised_combined_pgx hcpgx ON d.drugid =hcpgx.drug_name AND d.guideline = hcpgx.guidelineid AND d.implications = hcpgx.implications;
--CREATE UNIQUE INDEX t_drug_recommendations_import_u ON t_drug_recommendations_import (recommendation,recommendation_old,guidelineid ,drug_name,gene_name,rn);
--
--INSERT INTO prada.recommendation(
--recommendation,
--guideline,
--drugid,
--implications,
--gene_name,
--prada_start_dose,
--prada_target_dose,
--prada_titration_speed,
--prada_switch_drug,
--prada_switch1_drug,
--prada_switch1_gene,
--prada_switch2_drug,
--prada_switch2_gene,
--prada_tdm,
--drugrecommendation,
--version
--)
--SELECT 
--m.recommendation,
--m.guidelineid,
--m.drug_name,
--m.implications,
--m.gene_name,
--m.start_dose,
--m.target_dose,
--m.titration_speed,
--m.prada_switch_drug,
--m.prada_switch1_drug,
--m.prada_switch1_gene,
--m.prada_switch2_drug,
--m.prada_switch2_gene,
--m.tdm,
--m.drugrecommendation,
--1 --version
--FROM t_drug_recommendations_import m
----WHERE hcpgx.recommendation=5991110
--WHERE m.recommendation IS NOT NULL AND m.rn =1
----ORDER BY m.recommendation,m.guidelineid ,m.drug_name
;
