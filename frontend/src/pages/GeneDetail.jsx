import { useEffect, useState } from 'react'
import { Link, useParams } from 'react-router-dom'
import { motion } from 'framer-motion'
import {
  RadarChart, Radar, PolarGrid, PolarAngleAxis, ResponsiveContainer, Tooltip,
} from 'recharts'
import { ArrowLeft, ExternalLink, Shield, Pill, Dna, Activity, Info } from 'lucide-react'
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

  if (loading) return <div className="flex items-center justify-center py-40"><Spinner className="w-8 h-8" /></div>
  if (error || !gene) return <EmptyState title="Gene not found" subtitle={error ?? 'Check the gene ID'} />

  const ev = gene.evolution ?? {}
  const sub = scores?.sub_scores ?? {}

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
                        <div className="flex items-center gap-3">
                          <span className="text-xs text-text-muted">Position {m.start_pos}–{m.end_pos}</span>
                          <span className="text-xs text-text-muted capitalize">{o.species_id?.replace(/_/g, ' ')}</span>
                          {m.divergence_score != null && (
                            <span className="text-xs font-mono text-accent ml-auto">
                              Δ {(m.divergence_score * 100).toFixed(1)}%
                            </span>
                          )}
                        </div>
                        <div className="grid grid-cols-2 gap-3">
                          <div>
                            <p className="text-[10px] text-text-muted mb-1">Human</p>
                            <p className="font-mono text-xs text-text-primary tracking-widest bg-elevated rounded px-2 py-1.5">
                              {m.human_seq}
                            </p>
                          </div>
                          <div>
                            <p className="text-[10px] text-accent mb-1">Resilient species</p>
                            <p className="font-mono text-xs text-accent tracking-widest bg-accent/5 rounded px-2 py-1.5">
                              {m.animal_seq}
                            </p>
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
              </div>
            </Section>
          )}
        </div>
      </div>
    </div>
  )
}
