-- Manual SQL equivalent of Alembic migration 0023
-- Run this if alembic is unavailable or you prefer direct psql.
--
-- Usage:
--   PGPASSWORD="<password>" psql -h <host> -U bioresilient -d bioresilient -f migration_0023_manual.sql

BEGIN;

-- Add convergence_weight column (Fix 3: separate from phylop_score)
ALTER TABLE evolution_score
    ADD COLUMN IF NOT EXISTS convergence_weight FLOAT;

-- Add convergence_pval column (Fix 2: permutation null model p-value)
ALTER TABLE evolution_score
    ADD COLUMN IF NOT EXISTS convergence_pval FLOAT;

-- Back-fill convergence_weight from phylop_score for genes that have
-- convergence signal but whose phylop_score was set to the convergence weight
-- (i.e. before the UCSC PhyloP enrichment step ran on that gene).
UPDATE evolution_score
SET    convergence_weight = phylop_score
WHERE  convergence_count > 0
  AND  phylop_score IS NOT NULL
  AND  convergence_weight IS NULL;

-- Record this migration in alembic_version so alembic stays in sync
UPDATE alembic_version SET version_num = '0023';

COMMIT;

-- Verify
SELECT
    COUNT(*)                                     AS total_genes,
    COUNT(convergence_weight)                    AS with_convergence_weight,
    COUNT(convergence_pval)                      AS with_convergence_pval,
    ROUND(AVG(convergence_weight)::numeric, 3)   AS avg_convergence_weight
FROM evolution_score;
