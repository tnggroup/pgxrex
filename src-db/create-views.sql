--Expanded version of cpic.diplotype
--DROP MATERIALIZED VIEW prada.cpic_genetics;
CREATE MATERIALIZED VIEW prada.cpic_genetics --OR REPLACE
AS SELECT
g.chr,
g.symbol genesymbol,
g.genesequenceid,
g.proteinsequenceid,
g.mrnasequenceid,
g.hgncid,
g.ncbiid,
g.ensemblid,
g.pharmgkbid,
g.frequencymethods,
g.lookupmethod,
g.notesondiplotype,
g.url,
g.functionmethods,
g.notesonallelenaming,
grd.diplotype,
grd.diplotypekey,
grd.frequency diplotype_frequency,
grl.description,
grl.lookupkey grl_lookupkey,
gr.result,
gr.activityscore,
gr.ehrpriority,
gr.consultationtext,

(
	CASE
		WHEN gr.ehrpriority='Abnormal/Priority/High Risk' THEN 5
		WHEN gr.ehrpriority='Priority/High-Risk' THEN 4
		WHEN gr.ehrpriority='Normal Risk' THEN 3
		WHEN gr.ehrpriority='Normal/Routine/ Low Risk' THEN 2
		WHEN gr.ehrpriority='Normal/Routine/Low Risk' THEN 2
		WHEN gr.ehrpriority='Routine/Low-Risk' THEN 2
		ELSE 1
	END
) AS prada_ehrpriority_num,

(
	CASE
		WHEN gr.genesymbol IS NOT NULL THEN
			json_build_object(gr.genesymbol,
		    CASE
		        WHEN g.lookupmethod = 'ACTIVITY_SCORE'::text THEN grl.totalactivityscore
		        ELSE gr.result
		    END
		    )
		ELSE NULL
	END
) AS lookupkey
        
FROM
cpic.gene g
LEFT OUTER JOIN cpic.gene_result gr ON gr.genesymbol::text = g.symbol::text --ON g.symbol = gr.genesymbol
LEFT OUTER JOIN cpic.gene_result_lookup grl ON grl.phenotypeid =gr.id
LEFT OUTER JOIN cpic.gene_result_diplotype grd ON grd.functionphenotypeid = grl.id;


CREATE INDEX cpic_genetics_i ON prada.cpic_genetics (chr,genesymbol,ensemblid,prada_ehrpriority_num, diplotype, diplotypekey, grl_lookupkey);

--REFRESH MATERIALIZED VIEW prada.cpic_genetics;

--101,956 rows with inner join
--102,084 rows with full join

--SELECT * FROM prada.cpic_genetics WHERE diplotype IS NULL;
--SELECT DISTINCT genesymbol FROM prada.cpic_genetics;

--DROP MATERIALIZED VIEW prada.harmonised_cpic_pgx;
CREATE MATERIALIZED VIEW prada.harmonised_cpic_pgx --OR REPLACE
AS 
WITH a AS (
SELECT *,row_number() OVER (PARTITION BY chr, genesymbol, diplotype ORDER BY g1.allele) rn FROM (SELECT g.*,jsonb_object_keys(g.diplotypekey::jsonb -> g.genesymbol) allele FROM prada.cpic_genetics g) g1
),
g AS (
SELECT g.*, a1.allele::text allele1, a2.allele::text allele2 FROM prada.cpic_genetics g
LEFT OUTER JOIN a a1 ON g.chr=a1.chr AND g.genesymbol=a1.genesymbol AND g.diplotypekey=a1.diplotypekey AND a1.rn=1
LEFT OUTER JOIN a a2 ON g.chr=a2.chr AND g.genesymbol=a2.genesymbol AND g.diplotypekey=a2.diplotypekey AND a2.rn=2
)
SELECT
d.name drug_name,
d.drugid cpic_drugid,
d.pharmgkbid,
d.rxnormid,
d.drugbankid,
d.atcid,
d.umlscui,
d.flowchart,
d.version,
d.guidelineid,
p.cpiclevel,
p.pgkbcalevel,
p.usedforrecommendation,
gl.name guideline_name,
gl.url guideline_url,
(
--TRANSLATION BUSINESS RULES FOR DISCREPANCIES IN GENE SYMBOLS BETWEEN CPIC AND GENCODE
	CASE 
		WHEN p.genesymbol='GBA' THEN 'GBA1'
		ELSE p.genesymbol
	END
) AS genesymbol,
(
--TRANSLATION BUSINESS RULES FOR DISCREPANCIES IN ENSEMBLID BETWEEN CPIC AND GENCODE
	CASE 
		WHEN g.ensemblid='ENSG00000112096' THEN 'ENSG00000291237' --SOD2
		ELSE g.ensemblid
	END
) AS ensemblid,
g.chr,
g.diplotype,
g.description,
g.result,
g.activityscore,
g.ehrpriority,
g.consultationtext,
g.lookupkey,
g.diplotypekey,
g.diplotype_frequency,
g.prada_ehrpriority_num,
(
	CASE
		WHEN p.cpiclevel='A' THEN 8
		WHEN p.cpiclevel='A/B' THEN 7
		WHEN p.cpiclevel='B' THEN 6
		WHEN p.cpiclevel='B/C' THEN 5
		WHEN p.cpiclevel='C' THEN 4
		WHEN p.cpiclevel='C/D' THEN 3
		WHEN p.cpiclevel='D' THEN 2
		ELSE 1
	END
) AS prada_cpiclevel_num,
(
	CASE
		WHEN p.pgkbcalevel='1A' THEN 7
		WHEN p.pgkbcalevel='1B' THEN 6
		WHEN p.pgkbcalevel='2A' THEN 5
		WHEN p.pgkbcalevel='2B' THEN 4
		WHEN p.pgkbcalevel='3' THEN 3
		WHEN p.pgkbcalevel='4' THEN 2
		ELSE 1
	END
) AS prada_pgkbcalevel_num,
(
	CASE
		WHEN p.usedforrecommendation = TRUE THEN 2
		ELSE 1
	END
) AS prada_usedforrecommendation_num,
g.allele1,
g.allele2,
al1.frequency allele1_frequency,
al2.frequency allele2_frequency

--r.implications,
--r.drugrecommendation,
--r.classification,
--r.population,
--r.comments

FROM 
g 
LEFT OUTER JOIN cpic.pair p ON p.genesymbol = g.genesymbol AND p.removed = FALSE --AND p.usedforrecommendation = TRUE
LEFT OUTER JOIN cpic.drug d ON d.drugid = p.drugid
LEFT OUTER JOIN cpic.guideline gl ON d.guidelineid = gl.id
LEFT OUTER JOIN cpic.allele al1 ON g.genesymbol=al1.genesymbol AND g.allele1=al1.name
LEFT OUTER JOIN cpic.allele al2 ON g.genesymbol=al2.genesymbol AND g.allele2=al2.name
;

CREATE INDEX harmonised_cpic_pgx_i ON prada.harmonised_cpic_pgx (drug_name,cpic_drugid,pharmgkbid,chr,genesymbol,ensemblid);
CREATE INDEX harmonised_cpic_pgx_i2 ON prada.harmonised_cpic_pgx (pharmgkbid,chr,genesymbol,ensemblid,diplotype, diplotypekey);
CREATE INDEX harmonised_cpic_pgx_i3 ON prada.harmonised_cpic_pgx (prada_ehrpriority_num,prada_cpiclevel_num,prada_pgkbcalevel_num,prada_usedforrecommendation_num);

--REFRESH MATERIALIZED VIEW prada.harmonised_cpic_pgx;

-- 2,571,684 rows with LOJ pair, drugs, guidelines, no usedforrecommendation condition

-- OLD inner join cpic_genetics
-- 1,340,988 rows with full join
-- 1,340,482 rows with LOJ
-- 1,340,482 rows with innner join

--SELECT * FROM  prada.harmonised_cpic_pgx tp WHERE tp.drug_name='azathioprine' ORDER by drug_name;
--SELECT * FROM  prada.harmonised_cpic_pgx tp ORDER by drug_name;
--SELECT * FROM  prada.harmonised_cpic_pgx tp WHERE guideline_name IS NULL ORDER by drug_name;

--SELECT * FROM  prada.harmonised_cpic_pgx p WHERE p.genesymbol = 'SOD2';
/*
SELECT DISTINCT genesymbol
	FROM prada.harmonised_cpic_pgx px
	ORDER BY genesymbol;
	
SELECT DISTINCT cpiclevel
	FROM prada.harmonised_cpic_pgx px
	ORDER BY cpiclevel;
	
SELECT DISTINCT pgkbcalevel
	FROM prada.harmonised_cpic_pgx px
	ORDER BY pgkbcalevel;

SELECT DISTINCT ehrpriority
	FROM prada.harmonised_cpic_pgx px
	ORDER BY ehrpriority;
	

*/

--DROP VIEW prada.harmonised_cpic_pgx_ancestry_eur;
CREATE OR REPLACE VIEW prada.harmonised_cpic_pgx_ancestry_eur
AS
WITH pgx AS (
SELECT DISTINCT ON (g.chr, g.genesymbol, g.diplotype) 
g.chr, g.genesymbol, g.diplotype, g.result, 
(g.diplotype_frequency ->> 'European')::numeric diplotype_frequency_pop, (g.allele1_frequency ->> 'European')::numeric allele1_frequency_pop, (g.allele2_frequency ->> 'European')::numeric allele2_frequency_pop
FROM prada.harmonised_cpic_pgx g)
SELECT pgx.*,
COALESCE(pgx.diplotype_frequency_pop, 2*pgx.allele1_frequency_pop*pgx.allele2_frequency_pop) consensus_allele_frequency
FROM pgx;

CREATE OR REPLACE VIEW prada.harmonised_cpic_pgx_ancestry_eas
AS
WITH pgx AS (
SELECT DISTINCT ON (g.chr, g.genesymbol, g.diplotype) 
g.chr, g.genesymbol, g.diplotype, g.result, 
(g.diplotype_frequency ->> 'East Asian')::numeric diplotype_frequency_pop, (g.allele1_frequency ->> 'East Asian')::numeric allele1_frequency_pop, (g.allele2_frequency ->> 'East Asian')::numeric allele2_frequency_pop
FROM prada.harmonised_cpic_pgx g)
SELECT pgx.*,
COALESCE(pgx.diplotype_frequency_pop, 2*pgx.allele1_frequency_pop*pgx.allele2_frequency_pop) consensus_allele_frequency
FROM pgx;

CREATE OR REPLACE VIEW prada.harmonised_cpic_pgx_ancestry_amlathis
AS
WITH pgx AS (
SELECT DISTINCT ON (g.chr, g.genesymbol, g.diplotype) 
g.chr, g.genesymbol, g.diplotype, g.result, 
(g.diplotype_frequency ->> 'Latino')::numeric diplotype_frequency_pop, (g.allele1_frequency ->> 'Latino')::numeric allele1_frequency_pop, (g.allele2_frequency ->> 'Latino')::numeric allele2_frequency_pop
FROM prada.harmonised_cpic_pgx g)
SELECT pgx.*,
COALESCE(pgx.diplotype_frequency_pop, 2*pgx.allele1_frequency_pop*pgx.allele2_frequency_pop) consensus_allele_frequency
FROM pgx;

CREATE OR REPLACE VIEW prada.harmonised_cpic_pgx_ancestry_cassas
AS
WITH pgx AS (
SELECT DISTINCT ON (g.chr, g.genesymbol, g.diplotype) 
g.chr, g.genesymbol, g.diplotype, g.result, 
(g.diplotype_frequency ->> 'Central/South Asian')::numeric diplotype_frequency_pop, (g.allele1_frequency ->> 'Central/South Asian')::numeric allele1_frequency_pop, (g.allele2_frequency ->> 'Central/South Asian')::numeric allele2_frequency_pop
FROM prada.harmonised_cpic_pgx g)
SELECT pgx.*,
COALESCE(pgx.diplotype_frequency_pop, 2*pgx.allele1_frequency_pop*pgx.allele2_frequency_pop) consensus_allele_frequency
FROM pgx;

CREATE OR REPLACE VIEW prada.harmonised_cpic_pgx_ancestry_afr
AS
WITH pgx AS (
SELECT DISTINCT ON (g.chr, g.genesymbol, g.diplotype) 
g.chr, g.genesymbol, g.diplotype, g.result, 
(g.diplotype_frequency ->> 'Sub-Saharan African')::numeric diplotype_frequency_pop, (g.allele1_frequency ->> 'Sub-Saharan African')::numeric allele1_frequency_pop, (g.allele2_frequency ->> 'Sub-Saharan African')::numeric allele2_frequency_pop
FROM prada.harmonised_cpic_pgx g)
SELECT pgx.*,
COALESCE(pgx.diplotype_frequency_pop, 2*pgx.allele1_frequency_pop*pgx.allele2_frequency_pop) consensus_allele_frequency
FROM pgx;

CREATE OR REPLACE VIEW prada.harmonised_cpic_pgx_ancestry_amcarafr
AS
WITH pgx AS (
SELECT DISTINCT ON (g.chr, g.genesymbol, g.diplotype) 
g.chr, g.genesymbol, g.diplotype, g.result, 
(g.diplotype_frequency ->> 'African American/Afro-Caribbean')::numeric diplotype_frequency_pop, (g.allele1_frequency ->> 'African American/Afro-Caribbean')::numeric allele1_frequency_pop, (g.allele2_frequency ->> 'African American/Afro-Caribbean')::numeric allele2_frequency_pop
FROM prada.harmonised_cpic_pgx g)
SELECT pgx.*,
COALESCE(pgx.diplotype_frequency_pop, 2*pgx.allele1_frequency_pop*pgx.allele2_frequency_pop) consensus_allele_frequency
FROM pgx;

CREATE OR REPLACE VIEW prada.harmonised_cpic_pgx_ancestry_gme
AS
WITH pgx AS (
SELECT DISTINCT ON (g.chr, g.genesymbol, g.diplotype) 
g.chr, g.genesymbol, g.diplotype, g.result, 
(g.diplotype_frequency ->> 'Near Eastern')::numeric diplotype_frequency_pop, (g.allele1_frequency ->> 'Near Eastern')::numeric allele1_frequency_pop, (g.allele2_frequency ->> 'Near Eastern')::numeric allele2_frequency_pop
FROM prada.harmonised_cpic_pgx g)
SELECT pgx.*,
COALESCE(pgx.diplotype_frequency_pop, 2*pgx.allele1_frequency_pop*pgx.allele2_frequency_pop) consensus_allele_frequency
FROM pgx;


--DROP MATERIALIZED VIEW prada.harmonised_gene;
CREATE MATERIALIZED VIEW prada.harmonised_gene --OR REPLACE
AS 
SELECT
--cpicg.genesymbol,
--cpicg.ensemblid,
--cpicg.chr cpicchr,
--pgx."Gene",
g.chr,
c.number AS chrn,
g.bp1,
g.bp2,
g.gene_name,
g.gene_id,
g.gene_id0,
--pgx."CPIC" cpiclevel2,
--pgx."PharmGKB" pgkbcalevel2,
(pgx."Gene" IS NOT NULL) AS in_pgx,
(cpicg.genesymbol IS NOT NULL) AS in_cpic,
(cpic.genesymbol IS NOT NULL) AS has_cpiclevel,
(pgx."DPWG")::boolean in_dpwg,
(pgx."Twist")::boolean in_twist,
(pgx."PharmVar")::boolean in_pharmvar,
(pgx."CMRG")::boolean in_cmrg,
(pgx."PharmCAT")::boolean in_pharmcat
FROM
prada.gencode_gene g
LEFT OUTER JOIN (
SELECT DISTINCT genesymbol,ensemblid,chr
	FROM prada.harmonised_cpic_pgx WHERE genesymbol IS NOT NULL AND ensemblid IS NOT NULL AND chr IS NOT NULL
) AS cpicg ON (
	(cpicg.ensemblid=g.gene_id OR ((cpicg.ensemblid IS NULL OR g.gene_id IS NULL) AND g.gene_name = cpicg.genesymbol))
	AND g.chr=cpicg.chr
	)
LEFT OUTER JOIN (
SELECT DISTINCT genesymbol,ensemblid,chr
	FROM prada.harmonised_cpic_pgx 
	WHERE genesymbol IS NOT NULL AND ensemblid IS NOT NULL AND chr IS NOT NULL
	AND cpiclevel IS NOT NULL
) AS cpic ON (
	(cpic.ensemblid=g.gene_id OR ((cpic.ensemblid IS NULL OR g.gene_id IS NULL) AND g.gene_name = cpic.genesymbol))
	AND g.chr=cpic.chr
	)
LEFT OUTER JOIN prada.pgx_gene pgx ON g.gene_name=pgx."Gene" --OR cpicg.genesymbol=pgx."Gene" --AND g.chr=pgx.chr -- we don't have chr data here
INNER JOIN prada.chromosome c ON g.chr = c.name OR cpicg.chr = c.name
;

CREATE UNIQUE INDEX harmonised_gene_u ON prada.harmonised_gene (chrn,gene_id); -- chr,gene_name is not unique 
CREATE UNIQUE INDEX harmonised_gene_i ON prada.harmonised_gene (chrn,bp1,bp2,gene_name,gene_id);

--REFRESH MATERIALIZED VIEW prada.harmonised_gene;

--SELECT * FROM prada.harmonised_gene WHERE in_cpic=TRUE;
--SELECT * FROM prada.harmonised_gene WHERE in_pgx=TRUE;
--SELECT * FROM prada.harmonised_gene WHERE bp1 IS NULL;
--SELECT * FROM prada.harmonised_gene WHERE in_cpic=TRUE
--SELECT * FROM prada.harmonised_gene WHERE bp1 IS NULL;
--SELECT * FROM prada.gencode_gene gg WHERE gg.gene_name = 'SOD2';
--SELECT * FROM prada.harmonised_gene WHERE gene_name= 'Metazoa_SRP';
--SELECT DISTINCT chr FROM prada.harmonised_gene;

--CPIC recommendations are linked in here
--DROP VIEW prada.harmonised_combined_pgx;
CREATE OR REPLACE VIEW prada.harmonised_combined_pgx
AS SELECT
g.chr,
g.chrn,
g.bp1,
g.bp2,
g.gene_name,
g.gene_id,
g.gene_id0,
g.in_pgx,
g.in_cpic,
g.has_cpiclevel,
g.in_dpwg,
g.in_twist,
g.in_pharmvar,
g.in_cmrg,
g.in_pharmcat,
p.drug_name,
p.pharmgkbid,
p.rxnormid,
p.drugbankid,
p.atcid,
p.umlscui,
p.flowchart,
p.version,
p.guidelineid,
d.class drug_class,
(
CASE
	WHEN d.weight IS NULL THEN 1.0
	ELSE d.weight
END
) AS drug_weight,
p.cpiclevel,
p.pgkbcalevel,
p.usedforrecommendation,
p.guideline_name,
p.guideline_url,
p.diplotype,
p.description,
p.result,
p.activityscore,
p.ehrpriority,
p.consultationtext,
p.prada_ehrpriority_num,
p.prada_cpiclevel_num,
p.prada_pgkbcalevel_num,
p.prada_usedforrecommendation_num,
p.lookupkey,
p.diplotypekey,
p.diplotype_frequency,
p.allele1,
p.allele2,
p.allele1_frequency,
p.allele2_frequency,
r.id AS recommendation,
r.implications,
r.drugrecommendation,
r.phenotypes,
r.classification,
r.population,
r.comments
FROM 
prada.harmonised_gene g
LEFT OUTER JOIN prada.harmonised_cpic_pgx p ON --we trust the ensemblid and chr mapping because of the harmonised cpic pgx view
	p.ensemblid=g.gene_id
	AND g.chr=p.chr
LEFT OUTER JOIN prada.drug d ON (p.rxnormid=d.rxnormid OR (d.rxnormid IS NULL AND p.drug_name=d.name))
LEFT OUTER JOIN cpic.recommendation r ON p.cpic_drugid=r.drugid AND p.guidelineid = r.guidelineid AND p.lookupkey::jsonb <@  r.lookupkey 
WHERE g.in_pgx=TRUE OR g.in_cpic=TRUE OR g.has_cpiclevel=TRUE OR g.in_dpwg=TRUE OR g.in_twist=TRUE OR g.in_pharmvar=TRUE OR g.in_cmrg=TRUE OR g.in_pharmcat=TRUE OR p.ensemblid IS NOT NULL
--AND d.type !='cancer'
;

--SELECT * FROM  prada.harmonised_combined_pgx tp WHERE tp.drug_name='sertraline';
--SELECT * FROM  prada.harmonised_combined_pgx tp WHERE tp.drug_name='bupropion';
--SELECT * FROM  prada.harmonised_combined_pgx tp WHERE tp.drug_class IS NOT NULL;
--SELECT * FROM  prada.harmonised_combined_pgx tp ORDER by drug_name;
/*SELECT * FROM  prada.harmonised_combined_pgx tp ORDER by prada_ehrpriority_num DESC, prada_cpiclevel_num DESC, prada_pgkbcalevel_num DESC
LIMIT 100;*/
--SELECT DISTINCT chr,bp1,bp2,gene_name,gene_id FROM  prada.harmonised_combined_pgx tp WHERE prada_cpiclevel_num > 0 OR prada_pgkbcalevel_num > 0 OR prada_pgkbcalevel_num > 0 ORDER by chr,bp1,bp2;
--SELECT DISTINCT chr,bp1,bp2,gene_name,gene_id,drug_class FROM  prada.harmonised_combined_pgx tp WHERE tp.drug_class IS NOT NULL AND (prada_cpiclevel_num > 0 OR prada_pgkbcalevel_num > 0) ORDER by gene_name,chr,bp1,bp2; --antidepressant with some pgx literture support


--DROP MATERIALIZED VIEW prada.harmonised_combined_pgx_gene;
CREATE MATERIALIZED VIEW prada.harmonised_combined_pgx_gene
AS
WITH a AS (
SELECT
p.gene_name, p.gene_id, (p.chrn::integer) AS chr,
avg(p.bp1) bp1,
avg(p.bp2) bp2,
max(p.in_pgx::int) in_pgx_num,
max(p.in_cpic::int) in_cpic_num,
max(p.has_cpiclevel::int) has_cpiclevel_num,
max(p.in_dpwg::int) in_dpwg_num,
max(p.in_twist::int) in_twist_num,
max(p.in_pharmvar::int) in_pharmvar_num,
max(p.in_cmrg::int) in_cmrg_num,
max(p.in_pharmcat::int) in_pharmcat_num,
avg(prada_ehrpriority_num) prada_ehrpriority_avg,
avg(prada_cpiclevel_num) prada_cpiclevel_avg,
avg(prada_pgkbcalevel_num) prada_pgkbcalevel_avg,
avg(prada_usedforrecommendation_num) prada_usedforrecommendation_avg,
avg(p.drug_weight) prada_drug_weight_avg,
count(*) num_diplotype
FROM prada.harmonised_combined_pgx p
GROUP BY p.gene_name, p.gene_id, p.chrn
) 
SELECT
a.*,
--prada gene prioritissation business rule score
(1 + COALESCE(in_pgx_num,0)+COALESCE(in_cpic_num,0)+COALESCE(has_cpiclevel_num,0)+COALESCE(in_dpwg_num,0)+COALESCE(in_twist_num,0)+COALESCE(in_pharmvar_num,0)+COALESCE(in_cmrg_num,0)+COALESCE(in_pharmcat_num,0)) --coverage across databases
*COALESCE(prada_ehrpriority_avg,1)*COALESCE(prada_cpiclevel_avg,1)*COALESCE(prada_pgkbcalevel_avg,1)*COALESCE(prada_usedforrecommendation_avg,1) --database information level
*COALESCE(prada_drug_weight_avg,1) --Drug weight - less priority on cancer drugs etc
prada_gene_score
FROM a;

CREATE UNIQUE INDEX harmonised_combined_pgx_gene_u ON prada.harmonised_combined_pgx_gene (chr,gene_id); -- chr,gene_name is not unique 
CREATE UNIQUE INDEX harmonised_combined_pgx_gene_i ON prada.harmonised_combined_pgx_gene (prada_gene_score,chr,bp1,bp2,gene_name,gene_id);

--REFRESH MATERIALIZED VIEW prada.harmonised_combined_pgx_gene;

--SELECT * FROM prada.harmonised_combined_pgx_gene ORDER BY prada_gene_score DESC;
--SELECT * FROM prada.harmonised_combined_pgx_gene ORDER BY gene_name;

--DROP MATERIALIZED VIEW prada.harmonised_combined_drug;
CREATE MATERIALIZED VIEW prada.harmonised_combined_drug
AS
WITH a AS
(
SELECT
p.drug_name,
p.rxnormid,
count(*) num_diplotype
FROM prada.harmonised_combined_pgx p
GROUP BY p.drug_name, p.rxnormid
) SELECT 
d.name drug_name,
d."type",
d."class",
d.weight,
d2.drugid,
d2.pharmgkbid,
d2.rxnormid,
d2.drugbankid,
d2.atcid,
d2.umlscui,
d2.flowchart,
d2.version,
d2.guidelineid,
a.num_diplotype
FROM prada.drug d
LEFT OUTER JOIN cpic.drug d2 ON (d.rxnormid=d2.rxnormid OR (d.rxnormid IS NULL AND d.name=d2.name)) --This does not join on the CPIC drugid
LEFT OUTER JOIN a ON (a.rxnormid=d.rxnormid OR (d.rxnormid IS NULL AND a.drug_name=d.name))
;
CREATE UNIQUE INDEX harmonised_combined_drug_u ON prada.harmonised_combined_drug (drug_name);
CREATE UNIQUE INDEX harmonised_combined_drug_u2 ON prada.harmonised_combined_drug (rxnormid);
CREATE INDEX harmonised_combined_drug_i ON prada.harmonised_combined_drug (type,class);

--REFRESH MATERIALIZED VIEW prada.harmonised_combined_drug;

--SELECT * FROM prada.harmonised_combined_drug ORDER BY type,class,num_diplotype DESC, drug_name;
--SELECT * FROM prada.harmonised_combined_drug WHERE type = 'antidepressant' ORDER BY type,class,num_diplotype DESC, drug_name;
--SELECT * FROM prada.drug WHERE type ='antidepressant';



