BEGIN TRANSACTION;

-- DROP TABLE prada.chromosome;
CREATE TABLE prada.chromosome
(
    number integer NOT NULL,
    name text,
    sizebp integer NOT NULL,
    CONSTRAINT chromosome_pkey PRIMARY KEY (number)
);
COMMENT ON TABLE prada.chromosome IS 'Custom project data on chromosomes.';
CREATE UNIQUE INDEX chromosome_name_u ON prada.chromosome (name);


-- prada.gencode_gene
-- imported through R-script
CREATE INDEX gencode_gene_i ON prada.gencode_gene (chr,bp1,bp2,gene_name,gene_id);

-- prada.pgx_gene
-- imported manually from the ONT gene-list (as csv, values only from excel)
CREATE INDEX pgx_gene_i ON prada.pgx_gene ("Gene"); --Let's keep this non-unique for the future


-- DROP TABLE prada.variant;
CREATE TABLE prada.variant
(
	type smallint NOT NULL, --1: GWAS, 2: CNV
	snp text NOT NULL,
    chr smallint NOT NULL,
    bp integer NOT NULL,
    bp2 integer,
    mdd_p numeric,
    mdd_beta numeric,
    mdd_beta_se numeric,
    mdd_beta_n numeric,
    CONSTRAINT variant_pkey PRIMARY KEY (type,snp)
);
COMMENT ON TABLE prada.variant IS 'Custom project data on variants. To hold for example anchor information';
CREATE INDEX variant_i ON prada.variant (type,snp,chr,bp,bp2);
CREATE INDEX variant_i2 ON prada.variant (mdd_p,mdd_beta,mdd_beta_n,mdd_beta_se);

-- DROP TABLE prada.drug;
CREATE TABLE prada.drug
(
	rxnormid varchar(20), --should match the cpic rxnormid
	name varchar(100), --in case we don't have the id
    type text DEFAULT 'generic', --this is to separate for example cancer drugs from antidepressants
    class text, --this is to distinguish classes of for example antidepressants based on metabolising pathways
    weight double precision DEFAULT 1.0
    --CONSTRAINT drug_pkey PRIMARY KEY (rxnormid,name) --we have null values in rxnormid so won't fit into pkey
);
COMMENT ON TABLE prada.drug IS 'Custom project data on drugs. To prioritise drugs to include in the analyses.';
CREATE UNIQUE INDEX drug_u ON prada.drug (rxnormid,name);
ALTER TABLE prada.drug
  ADD CONSTRAINT drug_pkey PRIMARY KEY (name);
ALTER TABLE prada.drug
	ADD COLUMN selected_for_analysis smallint;

-- DROP TABLE prada.recommendation;
CREATE TABLE prada.recommendation
(
	recommendation integer NOT NULL,
	guideline integer NOT NULL,
	drugid text NOT NULL,
	gene_name text NOT NULL, --one row per gene!
	implications jsonb NOT NULL,
	prada_start_dose numeric,
	prada_target_dose numeric,
	prada_titration_speed smallint,
	prada_switch_drug text[],
	prada_switch1_drug smallint,
	prada_switch1_gene smallint,
	prada_switch2_drug smallint,
	prada_switch2_gene smallint,
	prada_tdm smallint,
	drugrecommendation text,
	time_change TIMESTAMP WITH TIME ZONE NOT NULL DEFAULT now(),
	version integer
);
COMMENT ON TABLE prada.recommendation IS 'Custom project data on clinical recommendations. Complement to primarily the cpic.recommendation table.';
--CREATE UNIQUE INDEX recommendation_u ON prada.recommendation (recommendation,guideline,drugid,gene_name);
ALTER TABLE prada.recommendation
  ADD CONSTRAINT recommendation_pkey PRIMARY KEY (recommendation,guideline,drugid,gene_name);


COMMIT;