import { useEffect, useState } from 'react'
import { Link, useParams } from 'react-router-dom'
import { motion } from 'framer-motion'
import {
  RadarChart, Radar, PolarGrid, PolarAngleAxis, ResponsiveContainer, Tooltip,
} from 'recharts'
import { ArrowLeft, ExternalLink, Shield, Pill, Dna, Activity, Info, FileText, Loader2, Syringe, AlertTriangle, Layers } from 'lucide-react'
import { api } from '@/lib/api'
import { TierBadge, ScoreBar, PageHeader, Spinner, EmptyState } from '@/components/ui'
import { formatFloat, cn } from '@/lib/utils'

const SCORE_META = {
  convergence:  { label: 'Resilience Signal',          tooltip: 'How many independent resilient lineages share this divergence' },
  selection:    { label: 'Evolutionary Pressure',      tooltip: 'Strength of positive natural selection (dN/dS)' },
  expression:   { label: 'Expression Difference',      tooltip: 'How differently this gene is expressed in resilient species' },
  disease:      { label: 'Disease Relevance',          tooltip: 'Association with diseases from OpenTargets, GWAS, gnomAD' },
  druggability: { label: 'Druggability',               tooltip: 'Binding pocket quality and known drug interactions' },
  safety:       { label: 'Safety Profile',             tooltip: 'How safe this target is based on network position and essentiality' },
  regulatory:   { label: 'Regulatory Divergence',      tooltip: 'Non-coding DNA changes near promoter (AlphaGenome Track B)' },
}

function Tooltip2({ text, children }) {
  return (
    <span className="group relative inline-flex items-center gap-1">
      {children}
      <Info className="w-3 h-3 text-text-muted/40 cursor-help" />
      <span className="absolute z-20 bottom-full left-1/2 -translate-x-1/2 mb-2 px-2.5 py-1.5 bg-elevated
                       border border-white/10 rounded-lg text-xs text-text-primary whitespace-nowrap
                       opacity-0 group-hover:opacity-100 pointer-events-none transition-opacity shadow-lg">
        {text}
      </span>
    </span>
  )
}

function ScoreRadar({ sub }) {
  const data = Object.entries(SCORE_META).map(([k, m]) => ({
    subject: m.label.split(' ')[0],
    value: Math.round((sub?.[k] ?? 0) * 100),
    fullMark: 100,
  }))
  return (
    <ResponsiveContainer width="100%" height={260}>
      <RadarChart data={data}>
        <PolarGrid stroke="rgba(255,255,255,0.08)" />
        <PolarAngleAxis dataKey="subject" tick={{ fill: '#64748b', fontSize: 11 }} />
        <Radar name="Score" dataKey="value" stroke="#00d4ff" fill="#00d4ff" fillOpacity={0.12} strokeWidth={1.5} />
        <Tooltip
          contentStyle={{ background: '#1a2235', border: '1px solid rgba(255,255,255,0.08)', borderRadius: 8, fontSize: 12 }}
          formatter={v => [`${v}`, 'Score']}
        />
      </RadarChart>
    </ResponsiveContainer>
  )
}

function Section({ title, icon: Icon, children }) {
  return (
    <div className="card space-y-4">
      <div className="flex items-center gap-2 pb-2 border-b border-white/5">
        {Icon && <Icon className="w-4 h-4 text-accent" />}
        <p className="font-semibold text-text-primary">{title}</p>
      </div>
      {children}
    </div>
  )
}

function MotifAlignment({ humanSeq, animalSeq }) {
  const len = Math.max(humanSeq?.length ?? 0, animalSeq?.length ?? 0)
  const h = (humanSeq ?? '').padEnd(len, ' ')
  const a = (animalSeq ?? '').padEnd(len, ' ')
  return (
    <div className="space-y-2">
      <div>
        <p className="text-[10px] text-text-muted mb-0.5">Human</p>
        <div className="font-mono text-xs tracking-widest flex flex-wrap gap-0.5">
          {h.split('').map((ch, i) => (
            <span
              key={`h-${i}`}
              className={ch !== a[i] && a[i] !== ' ' ? 'text-amber-400 bg-amber-500/20 rounded px-0.5' : 'text-text-primary'}
            >
              {ch}
            </span>
          ))}
        </div>
      </div>
      <div>
        <p className="text-[10px] text-accent mb-0.5">Resilient species</p>
        <div className="font-mono text-xs tracking-widest flex flex-wrap gap-0.5">
          {a.split('').map((ch, i) => (
            <span
              key={`a-${i}`}
              className={ch !== h[i] && ch !== ' ' ? 'text-amber-400 bg-amber-500/20 rounded px-0.5' : 'text-text-muted'}
            >
              {ch}
            </span>
          ))}
        </div>
      </div>
    </div>
  )
}

function KV({ label, value, mono, tooltip }) {
  const val = value ?? '—'
  const content = <span className={cn('text-sm text-text-primary', mono && 'font-mono text-accent')}>{val}</span>
  return (
    <div className="flex items-start justify-between gap-4 py-1.5 border-b border-white/3 last:border-0">
      <p className="text-xs text-text-muted shrink-0">
        {tooltip ? <Tooltip2 text={tooltip}>{label}</Tooltip2> : label}
      </p>
      {content}
    </div>
  )
}

export default function GeneDetail() {
  const { id } = useParams()
  const [gene, setGene] = useState(null)
  const [scores, setScores] = useState(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)
  const [narrativeLoading, setNarrativeLoading] = useState(false)
  const [narrativeOverride, setNarrativeOverride] = useState(null)

  useEffect(() => {
    Promise.all([
      api.getCandidate(id),
      api.getScores(id).catch(() => null),
    ]).then(([g, s]) => {
      setGene(g)
      setScores(s)
      setLoading(false)
    }).catch(e => {
      setError(e.message)
      setLoading(false)
    })
  }, [id])

  const refreshGene = () => {
    api.getCandidate(id).then(setGene).catch(() => {})
  }

  const onGenerateNarrative = () => {
    setNarrativeLoading(true)
    api.generateNarrative({ gene_id: id, force: false })
      .then((res) => {
        setNarrativeOverride(res.narrative)
        refreshGene()
      })
      .finally(() => setNarrativeLoading(false))
  }

  if (loading) return <div className="flex items-center justify-center py-40"><Spinner className="w-8 h-8" /></div>
  if (error || !gene) return <EmptyState title="Gene not found" subtitle={error ?? 'Check the gene ID'} />

  const ev = gene.evolution ?? {}
  const sub = scores?.sub_scores ?? {}
  const narrative = narrativeOverride ?? gene.narrative ?? null

  return (
    <div className="px-8 py-8 space-y-6">
      {/* Header */}
      <div>
        <Link to="/candidates" className="inline-flex items-center gap-1.5 text-xs text-text-muted hover:text-text-primary mb-4 group">
          <ArrowLeft className="w-3.5 h-3.5 group-hover:-translate-x-0.5 transition-transform" /> All candidates
        </Link>
        <div className="flex items-start justify-between">
          <div className="flex items-center gap-4">
            <div>
              <h1 className="text-3xl font-bold font-mono text-accent">{gene.gene_symbol}</h1>
              <p className="text-sm text-text-muted mt-0.5">{gene.human_protein
                ? <a href={`https://www.uniprot.org/uniprot/${gene.human_protein}`} target="_blank" rel="noreferrer"
                     className="hover:text-accent flex items-center gap-1">{gene.human_protein} <ExternalLink className="w-3 h-3" /></a>
                : 'No UniProt accession'
              }</p>
            </div>
            <TierBadge tier={gene.tier} />
          </div>
          <div className="text-right">
            <p className="label-muted">Composite Score</p>
            <p className="text-4xl font-bold text-text-primary">{(gene.composite_score * 100).toFixed(1)}</p>
          </div>
        </div>
      </div>

      {/* Research Summary */}
      <div className="card">
          <div className="flex items-center gap-2 pb-2 border-b border-white/5">
            <FileText className="w-4 h-4 text-accent" />
            <p className="font-semibold text-text-primary">Research Summary</p>
          </div>
          {narrativeLoading ? (
            <div className="flex items-center gap-2 py-4 text-text-muted">
              <Loader2 className="w-4 h-4 animate-spin" />
              <span className="text-sm">Generating summary…</span>
            </div>
          ) : narrative ? (
            <div className="prose prose-invert prose-sm max-w-none text-text-primary whitespace-pre-wrap pt-2">
              {narrative}
            </div>
          ) : (
            <button
              type="button"
              onClick={onGenerateNarrative}
              className="mt-2 px-4 py-2 rounded-lg bg-accent/15 text-accent border border-accent/30 hover:bg-accent/25 text-sm font-medium transition-colors"
            >
              Generate Summary
            </button>
          )}
        </div>

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        {/* Score radar */}
        <div className="card lg:col-span-1">
          <p className="font-semibold text-text-primary mb-2">Score Breakdown</p>
          <ScoreRadar sub={sub} />
          <div className="space-y-2 mt-2">
            {Object.entries(SCORE_META).map(([k, m]) => (
              <div key={k} className="flex items-center justify-between gap-2">
                <p className="text-xs text-text-muted shrink-0 w-28">
                  <Tooltip2 text={m.tooltip}>{m.label}</Tooltip2>
                </p>
                <ScoreBar value={sub[k] ?? 0} className="flex-1" />
              </div>
            ))}
          </div>
        </div>

        {/* Evolution details */}
        <div className="space-y-4 lg:col-span-2">
          <Section title="Evolutionary Evidence" icon={Activity}>
            <div className="grid grid-cols-2 gap-x-8">
              <KV label="dN/dS ratio" value={formatFloat(ev.dnds_ratio)} mono
                tooltip="Ratio of non-synonymous to synonymous substitution rate. >1 = positive selection." />
              <KV label="Selection p-value" value={ev.dnds_pvalue != null ? ev.dnds_pvalue.toExponential(2) : '—'} mono />
              <KV label="Convergent lineages" value={ev.convergence_count} mono
                tooltip="Number of independent lineages showing the same motif change." />
              <KV label="PhyloP score" value={formatFloat(ev.phylop_score)} mono
                tooltip="Conservation score from UCSC. Higher = more conserved." />
              <KV label="Selection model" value={ev.selection_model} />
              <KV label="Branches selected"
                value={ev.branches_under_selection?.length ? ev.branches_under_selection.join(', ') : '—'} />
              {scores?.fel_sites != null && (
                <KV label="FEL sites (pervasive)" value={scores.fel_sites} mono
                  tooltip="Sites under pervasive positive selection (FEL p<0.05). Complements MEME episodic detection." />
              )}
              {scores?.busted_pvalue != null && (
                <KV label="BUSTED p-value" value={scores.busted_pvalue.toExponential(2)} mono
                  tooltip="Gene-wide episodic selection test p-value (BUSTED). <0.05 = gene has experienced positive selection." />
              )}
              {ev.relax_k != null && (
                <KV label="RELAX k" value={formatFloat(ev.relax_k)} mono
                  tooltip="RELAX selection intensity parameter. k>1 = stronger selection in resilient branches; k<1 = relaxed selection." />
              )}
              {ev.relax_pvalue != null && (
                <KV label="RELAX p-value" value={ev.relax_pvalue.toExponential(2)} mono
                  tooltip="RELAX p-value for branch-specific selection rate shift. <0.05 = significant rate change." />
              )}
            </div>
          </Section>

          {/* Orthologs */}
          {gene.orthologs?.length > 0 && (
            <Section title="Species Orthologs" icon={Dna}>
              <div className="overflow-x-auto">
                <table className="w-full text-sm">
                  <thead>
                    <tr className="border-b border-white/5">
                      <th className="label-muted text-left py-2 pr-4">Species</th>
                      <th className="label-muted text-left py-2 pr-4">Protein ID</th>
                      <th className="label-muted text-left py-2 pr-4">Identity</th>
                      <th className="label-muted text-left py-2">Divergent motifs</th>
                    </tr>
                  </thead>
                  <tbody>
                    {gene.orthologs.map(o => (
                      <tr key={o.id} className="border-b border-white/3 hover:bg-elevated">
                        <td className="py-2 pr-4 text-text-primary font-medium capitalize">
                          {o.species_id?.replace(/_/g, ' ')}
                        </td>
                        <td className="py-2 pr-4 font-mono text-xs text-text-muted">{o.protein_id ?? '—'}</td>
                        <td className="py-2 pr-4">
                          {o.sequence_identity_pct != null ? (
                            <span className={cn('font-mono text-xs',
                              o.sequence_identity_pct < 85 ? 'text-accent' : 'text-text-muted',
                            )}>
                              {o.sequence_identity_pct.toFixed(1)}%
                            </span>
                          ) : '—'}
                        </td>
                        <td className="py-2">
                          <span className={cn('text-xs', o.motifs?.length ? 'text-warning' : 'text-text-muted/40')}>
                            {o.motifs?.length ?? 0}
                          </span>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>

              {/* Motif details */}
              {gene.orthologs.some(o => o.motifs?.length) && (
                <div className="mt-4 space-y-3">
                  <p className="label-muted">Divergent Motifs</p>
                  {gene.orthologs.flatMap(o =>
                    (o.motifs ?? []).map(m => (
                      <div key={m.id} className="bg-base rounded-lg px-4 py-3 space-y-2 border border-white/5">
                        <div className="flex items-center gap-3 flex-wrap">
                          <span className="text-xs text-text-muted">Position {m.start_pos}–{m.end_pos}</span>
                          <span className="text-xs text-text-muted capitalize">{o.species_id?.replace(/_/g, ' ')}</span>
                          {m.in_functional_domain && (
                            <span className="px-1.5 py-0.5 rounded text-[10px] font-medium bg-emerald-500/15 text-emerald-400 border border-emerald-500/25"
                              title={`Pfam/InterPro domain: ${m.domain_name || 'domain'}`}>
                              {m.domain_name ? m.domain_name.slice(0, 20) : 'Functional domain'}
                            </span>
                          )}
                          {m.convergent_aa_count >= 2 && (
                            <span className="px-1.5 py-0.5 rounded text-[10px] font-medium bg-violet-500/15 text-violet-400 border border-violet-500/25"
                              title={`True convergent substitution: same amino acid change in ${m.convergent_aa_count} independent lineages`}>
                              Conv ×{m.convergent_aa_count}
                            </span>
                          )}
                          {m.consequence_score != null && (
                            <span
                              className={cn(
                                'px-1.5 py-0.5 rounded text-[10px] font-mono font-medium border',
                                m.consequence_score >= 0.7
                                  ? 'bg-red-500/15 text-red-400 border-red-500/25'
                                  : m.consequence_score >= 0.4
                                  ? 'bg-amber-500/15 text-amber-400 border-amber-500/25'
                                  : 'bg-white/5 text-text-muted border-white/10',
                              )}
                              title="AlphaMissense pathogenicity score for this divergence position"
                            >
                              AM {m.consequence_score.toFixed(2)}
                            </span>
                          )}
                          {m.esm1v_score != null && (
                            <span
                              className={cn(
                                'px-1.5 py-0.5 rounded text-[10px] font-mono border',
                                m.esm1v_score < -1.0
                                  ? 'bg-blue-500/15 text-blue-400 border-blue-500/25'
                                  : 'bg-white/5 text-text-muted border-white/10',
                              )}
                              title="ESM-1v log-likelihood ratio — negative values indicate unusual/adaptive substitution"
                            >
                              ESM {m.esm1v_score.toFixed(2)}
                            </span>
                          )}
                          {m.motif_direction && m.motif_direction !== 'neutral' && (
                            <span
                              className={cn(
                                'px-1.5 py-0.5 rounded text-[10px] font-semibold border',
                                m.motif_direction === 'gain_of_function'
                                  ? 'bg-green-500/15 text-green-400 border-green-500/25'
                                  : m.motif_direction === 'loss_of_function'
                                  ? 'bg-orange-500/15 text-orange-400 border-orange-500/25'
                                  : 'bg-red-500/15 text-red-400 border-red-500/25',
                              )}
                              title={
                                m.motif_direction === 'gain_of_function' ? 'Predicted gain-of-function substitution' :
                                m.motif_direction === 'loss_of_function' ? 'Predicted protective loss-of-function (LoF-tolerant gene)' :
                                'Likely pathogenic — LoF-intolerant gene'
                              }
                            >
                              {m.motif_direction === 'gain_of_function' ? '↑ GoF' :
                               m.motif_direction === 'loss_of_function' ? '↓ LoF' : '⚠ Pathogenic'}
                            </span>
                          )}
                          {m.divergence_score != null && (
                            <span className="text-xs font-mono text-accent ml-auto">
                              Δ {(m.divergence_score * 100).toFixed(1)}%
                            </span>
                          )}
                        </div>
                        <div>
                          <p className="text-[10px] text-text-muted mb-1.5">Alignment (divergent positions in amber)</p>
                          <div className="bg-elevated rounded px-3 py-2">
                            <MotifAlignment humanSeq={m.human_seq} animalSeq={m.animal_seq} />
                          </div>
                        </div>
                      </div>
                    ))
                  )}
                </div>
              )}
            </Section>
          )}

          {/* Phase 2 data — shown if populated */}
          {gene.disease && (
            <Section title="Disease Association" icon={Shield}>
              <div className="grid grid-cols-2 gap-x-8">
                <KV label="Disease" value={gene.disease.disease_name} />
                <KV label="OpenTargets score" value={formatFloat(gene.disease.opentargets_score)} mono />
                <KV label="GWAS p-value" value={gene.disease.gwas_pvalue?.toExponential(2)} mono />
                <KV label="gnomAD pLI" value={formatFloat(gene.disease.gnomad_pli)} mono
                  tooltip="Loss-of-function intolerance score. >0.9 = essential gene." />
                <KV label="Mouse KO phenotype" value={gene.disease.mouse_ko_phenotype} />
                    {gene.disease.protective_variant_count != null && gene.disease.protective_variant_count > 0 && (
                  <>
                    <KV label="Protective variants" value={gene.disease.protective_variant_count} mono
                      tooltip="Rare human variants (MAF<1%) at positions matching animal divergence direction — the PCSK9 paradigm." />
                    <KV label="Best protective trait" value={gene.disease.best_protective_trait} />
                    {gene.disease.protective_variant_pvalue != null && (
                      <KV label="Protective GWAS p-value"
                        value={gene.disease.protective_variant_pvalue.toExponential(2)} mono
                        tooltip="GWAS p-value for the most significant protective phenotype association." />
                    )}
                  </>
                )}
                {gene.disease.lit_score != null && (
                  <KV label="Literature score" value={gene.disease.lit_score.toFixed(3)} mono
                    tooltip={`PubMed citations in resilience literature (${gene.disease.lit_pmid_count ?? 0} papers). Higher = better-validated gene.`} />
                )}
              </div>
            </Section>
          )}

          {gene.drug_target && (
            <Section title="Druggability" icon={Pill}>
              <div className="grid grid-cols-2 gap-x-8">
                <KV label="Binding pockets" value={gene.drug_target.pocket_count} mono />
                <KV label="Top pocket score" value={formatFloat(gene.drug_target.top_pocket_score)} mono />
                <KV label="Druggability tier" value={gene.drug_target.druggability_tier} />
                <KV label="ChEMBL target" value={gene.drug_target.chembl_target_id} mono />
                <KV label="Known drugs" value={gene.drug_target.existing_drugs?.join(', ')} />
                {gene.drug_target.p2rank_score != null && (
                  <>
                    <KV label="P2Rank score" value={gene.drug_target.p2rank_score.toFixed(3)} mono
                      tooltip="ML pocket binding probability from P2Rank. Closer to 1.0 = more druggable." />
                    <KV label="P2Rank pockets" value={gene.drug_target.p2rank_pocket_count} mono />
                  </>
                )}
              </div>
            </Section>
          )}

          {gene.gene_therapy && (
            <Section title="Gene Therapy Feasibility" icon={Syringe}>
              <div className="grid grid-cols-2 gap-x-8">
                <KV label="Gene size (bp)" value={gene.gene_therapy.gene_size_bp?.toLocaleString()} mono
                  tooltip="Coding sequence size. AAV payload limit is ~4.7 kb." />
                <KV label="AAV compatible" value={gene.gene_therapy.aav_compatible ? 'Yes' : 'No'}
                  tooltip="Whether the gene fits within AAV packaging constraints." />
                <KV label="CRISPR guide sites" value={gene.gene_therapy.crispr_sites} mono
                  tooltip="Number of valid CRISPR guide RNA sites with acceptable on-target scores." />
                <KV label="Off-target risk" value={gene.gene_therapy.offtarget_risk}
                  tooltip="Predicted CRISPR off-target risk level." />
                {gene.gene_therapy.tissue_tropism?.length > 0 && (
                  <KV label="AAV tissue tropism" value={gene.gene_therapy.tissue_tropism.join(', ')}
                    tooltip="AAV serotypes with tropism for the relevant tissue." />
                )}
              </div>
            </Section>
          )}

          {gene.safety && (
            <Section title="Safety Profile" icon={AlertTriangle}>
              <div className="grid grid-cols-2 gap-x-8">
                <KV label="Essential gene" value={gene.safety.is_essential ? 'Yes' : 'No'}
                  tooltip="pLI > 0.9 indicates loss-of-function intolerance." />
                <KV label="Hub gene risk" value={gene.safety.hub_risk ? 'Yes' : 'No'}
                  tooltip="Network degree > 50 in STRING protein-protein interaction network." />
                <KV label="Network degree" value={gene.safety.network_degree} mono
                  tooltip="Number of protein-protein interactions in STRING." />
                <KV label="Protein family size" value={gene.safety.family_size} mono
                  tooltip="Large families increase off-target risk from cross-reactivity." />
                {gene.safety.depmap_score != null && (
                  <KV label="DepMap essentiality" value={gene.safety.depmap_score.toFixed(3)} mono
                    tooltip="DepMap CRISPR chronos score. More negative = more essential across cancer cell lines. <-0.5 = broadly essential." />
                )}
                {gene.safety.gtex_tissue_count != null && (
                  <KV label="GTEx tissues expressed" value={gene.safety.gtex_tissue_count} mono
                    tooltip="Number of tissues with TPM > 1 in GTEx. >30 = ubiquitous expression, harder to target safely." />
                )}
                {gene.safety.gtex_max_tpm != null && (
                  <KV label="GTEx max TPM" value={gene.safety.gtex_max_tpm.toFixed(1)} mono
                    tooltip="Maximum median TPM across all GTEx tissues." />
                )}
                {gene.safety.phewas_hits && Object.keys(gene.safety.phewas_hits).length > 0 && (
                  <div className="col-span-2 mt-2">
                    <p className="text-xs text-text-muted mb-1">PheWAS associations</p>
                    <div className="flex flex-wrap gap-1.5">
                      {Object.entries(gene.safety.phewas_hits).slice(0, 10).map(([trait, pval]) => (
                        <span key={trait} className="px-2 py-0.5 rounded text-[10px] bg-white/5 text-text-muted border border-white/8">
                          {trait}: {Number(pval).toExponential(1)}
                        </span>
                      ))}
                    </div>
                  </div>
                )}
              </div>
            </Section>
          )}

          {gene.regulatory?.length > 0 && (
            <Section title="Regulatory Divergence (AlphaGenome)" icon={Layers}>
              <p className="text-xs text-text-muted mb-3">Non-coding sequence changes near promoter regions predicted by AlphaGenome.</p>
              <div className="overflow-x-auto">
                <table className="w-full text-sm">
                  <thead>
                    <tr className="border-b border-white/5">
                      <th className="label-muted text-left py-2 pr-4">Species</th>
                      <th className="label-muted text-left py-2 pr-4">Promoter Divergence</th>
                      <th className="label-muted text-left py-2 pr-4">Expression log2FC</th>
                      <th className="label-muted text-left py-2 pr-4">Lineages</th>
                      <th className="label-muted text-left py-2">Score</th>
                    </tr>
                  </thead>
                  <tbody>
                    {gene.regulatory.map((r, i) => (
                      <tr key={i} className="border-b border-white/3 hover:bg-elevated">
                        <td className="py-2 pr-4 text-text-primary capitalize">{r.species_id?.replace(/_/g, ' ')}</td>
                        <td className="py-2 pr-4 font-mono text-xs text-text-muted">{formatFloat(r.promoter_divergence)}</td>
                        <td className="py-2 pr-4 font-mono text-xs text-text-muted">{formatFloat(r.expression_log2fc)}</td>
                        <td className="py-2 pr-4 font-mono text-xs text-text-muted">{r.lineage_count ?? '—'}</td>
                        <td className="py-2 font-mono text-xs text-accent">{formatFloat(r.regulatory_score)}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </Section>
          )}
        </div>
      </div>
    </div>
  )
}
