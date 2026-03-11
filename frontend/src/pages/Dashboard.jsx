import { useEffect, useState } from 'react'
import { Link } from 'react-router-dom'
import { motion } from 'framer-motion'
import {
  BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, Cell,
} from 'recharts'
import { FlaskConical, Activity, ArrowRight, CheckCircle2, AlertCircle, Clock } from 'lucide-react'
import { api } from '@/lib/api'
import { StatCard, PageHeader, Spinner, TierBadge } from '@/components/ui'
import { relativeTime, formatFloat } from '@/lib/utils'

const FADE_UP = (i) => ({
  initial: { opacity: 0, y: 20 },
  animate: { opacity: 1, y: 0 },
  transition: { delay: i * 0.07, duration: 0.4, ease: 'easeOut' },
})

const TIER_COLOR = { Validated: '#10b981', Tier1: '#7c3aed', Tier2: '#f59e0b', Tier3: '#334155' }

function StatusPill({ status }) {
  const map = {
    idle:     { label: 'Not started',  icon: Clock,        cls: 'text-text-muted bg-white/5' },
    running:  { label: 'Running',      icon: Activity,     cls: 'text-accent bg-accent/10 animate-pulse' },
    complete: { label: 'Complete',     icon: CheckCircle2, cls: 'text-success bg-success/10' },
    failed:   { label: 'Failed',       icon: AlertCircle,  cls: 'text-danger bg-danger/10' },
    stopped:  { label: 'Stopped',      icon: AlertCircle,  cls: 'text-warning bg-warning/10' },
  }
  const { label, icon: Icon, cls } = map[status] ?? map.idle
  return (
    <span className={`inline-flex items-center gap-1.5 px-3 py-1 rounded-full text-xs font-semibold ${cls}`}>
      <Icon className="w-3 h-3" />{label}
    </span>
  )
}

export default function Dashboard() {
  const [candidates, setCandidates] = useState(null)
  const [pipelineStatus, setPipelineStatus] = useState(null)
  const [traits, setTraits] = useState([])
  const [traitId, setTraitId] = useState('')
  const [loading, setLoading] = useState(true)

  useEffect(() => {
    api.getTraits().then(setTraits).catch(() => [])
  }, [])

  useEffect(() => {
    setLoading(true)
    const params = { limit: 20 }
    if (traitId) params.trait_id = traitId
    Promise.all([
      api.getCandidates(params).catch(() => []),
      api.getPipelineStatus().catch(() => null),
    ]).then(([cands, status]) => {
      setCandidates(cands ?? [])
      setPipelineStatus(status)
      setLoading(false)
    })
  }, [traitId])

  const validated = candidates?.filter(c => c.tier === 'Validated').length ?? 0
  const tier1 = candidates?.filter(c => c.tier === 'Tier1').length ?? 0
  const tier2 = candidates?.filter(c => c.tier === 'Tier2').length ?? 0
  const tier3 = candidates?.filter(c => c.tier === 'Tier3').length ?? 0
  const top10 = candidates?.slice(0, 10) ?? []

  const completedSteps = pipelineStatus
    ? Object.values(pipelineStatus.steps || {}).filter(s => s.status === 'complete').length
    : 0
  const totalSteps = pipelineStatus ? Object.keys(pipelineStatus.steps || {}).length : 18

  return (
    <div className="px-8 py-8 space-y-8">
      <PageHeader
        title="Dashboard"
        subtitle="Overview of your BioResilient analysis"
      >
        {traits.length > 0 && (
          <select
            value={traitId}
            onChange={e => setTraitId(e.target.value)}
            className="px-3 py-2 bg-surface border border-white/8 rounded-lg text-sm text-text-primary focus:outline-none focus:border-accent/50"
          >
            <option value="">Default analysis</option>
            {traits.map(t => (
              <option key={t.id} value={t.id}>{t.label}</option>
            ))}
          </select>
        )}
        <Link to="/pipeline" className="btn-secondary text-sm">
          <Activity className="w-4 h-4" /> Pipeline
        </Link>
        <Link to="/candidates" className="btn-primary text-sm">
          <FlaskConical className="w-4 h-4" /> View Candidates
        </Link>
      </PageHeader>

      {loading ? (
        <div className="flex items-center justify-center py-32"><Spinner className="w-8 h-8" /></div>
      ) : (
        <>
          {/* Tier summary cards */}
          <div className="grid grid-cols-2 lg:grid-cols-5 gap-4">
            {[
              { label: 'Validated', value: validated, sub: 'Evolutionary + human genetic evidence', accent: true },
              { label: 'Tier 1 Candidates', value: tier1, sub: 'High-confidence targets', accent: true },
              { label: 'Tier 2 Candidates', value: tier2, sub: 'Promising — needs validation' },
              { label: 'Tier 3 Candidates', value: tier3, sub: 'Low priority' },
              { label: 'Pipeline Progress', value: `${completedSteps}/${totalSteps}`, sub: 'Steps complete' },
            ].map((card, i) => (
              <motion.div key={card.label} {...FADE_UP(i)}>
                <StatCard {...card} />
              </motion.div>
            ))}
          </div>

          {/* Pipeline status banner */}
          {pipelineStatus && (
            <motion.div {...FADE_UP(4)} className="card flex items-center justify-between gap-4">
              <div className="flex items-center gap-4">
                <StatusPill status={pipelineStatus.status} />
                <div>
                  <p className="text-sm font-medium text-text-primary">
                    {pipelineStatus.status === 'running'
                      ? `Running: ${pipelineStatus.steps?.[pipelineStatus.current_step]?.label ?? pipelineStatus.current_step}`
                      : pipelineStatus.status === 'complete'
                      ? 'All analysis complete'
                      : pipelineStatus.status === 'failed'
                      ? `Failed at: ${pipelineStatus.current_step}`
                      : 'No pipeline run yet'}
                  </p>
                  {pipelineStatus.started_at && (
                    <p className="text-xs text-text-muted">
                      Started {relativeTime(pipelineStatus.started_at)}
                      {pipelineStatus.elapsed_total_s && ` · ${(pipelineStatus.elapsed_total_s / 60).toFixed(1)} min total`}
                    </p>
                  )}
                </div>
              </div>
              <Link to="/pipeline" className="btn-ghost text-sm">
                View details <ArrowRight className="w-4 h-4" />
              </Link>
            </motion.div>
          )}

          {/* Top 10 chart */}
          {top10.length > 0 && (
            <motion.div {...FADE_UP(5)} className="card">
              <div className="flex items-center justify-between mb-5">
                <div>
                  <p className="font-semibold text-text-primary">Top Candidates</p>
                  <p className="text-xs text-text-muted">Ranked by composite score</p>
                </div>
                <Link to="/candidates" className="text-xs text-accent hover:underline flex items-center gap-1">
                  View all <ArrowRight className="w-3 h-3" />
                </Link>
              </div>
              <ResponsiveContainer width="100%" height={240}>
                <BarChart data={top10} layout="vertical" margin={{ left: 60, right: 20 }}>
                  <XAxis type="number" domain={[0, 1]} tickFormatter={v => `${(v * 100).toFixed(0)}`}
                    tick={{ fill: '#64748b', fontSize: 11 }} axisLine={false} tickLine={false} />
                  <YAxis type="category" dataKey="gene_symbol" tick={{ fill: '#f1f5f9', fontSize: 12, fontFamily: 'JetBrains Mono' }}
                    width={56} axisLine={false} tickLine={false} />
                  <Tooltip
                    cursor={{ fill: 'rgba(255,255,255,0.03)' }}
                    contentStyle={{ background: '#1a2235', border: '1px solid rgba(255,255,255,0.08)', borderRadius: 8, fontSize: 12 }}
                    formatter={(v, n, p) => [`${(v * 100).toFixed(1)}`, 'Score']}
                  />
                  <Bar dataKey="composite_score" radius={[0, 4, 4, 0]}>
                    {top10.map((c) => (
                      <Cell key={c.gene_id} fill={TIER_COLOR[c.tier] ?? '#334155'} />
                    ))}
                  </Bar>
                </BarChart>
              </ResponsiveContainer>
            </motion.div>
          )}

          {/* Quick candidate list */}
          {top10.length > 0 && (
            <motion.div {...FADE_UP(6)} className="card">
              <p className="font-semibold text-text-primary mb-4">Recent Top Targets</p>
              <div className="space-y-2">
                {top10.slice(0, 5).map((c, i) => (
                  <Link
                    key={c.gene_id}
                    to={`/candidates/${c.gene_id}`}
                    className="flex items-center justify-between px-4 py-3 rounded-lg hover:bg-elevated transition-colors group"
                  >
                    <div className="flex items-center gap-3">
                      <span className="text-xs text-text-muted w-5">{i + 1}</span>
                      <span className="font-mono text-sm text-accent font-medium">{c.gene_symbol}</span>
                      <TierBadge tier={c.tier} />
                    </div>
                    <div className="flex items-center gap-3">
                      <span className="text-sm font-medium text-text-primary">
                        {(c.composite_score * 100).toFixed(1)}
                      </span>
                      <ArrowRight className="w-4 h-4 text-text-muted opacity-0 group-hover:opacity-100 transition-opacity" />
                    </div>
                  </Link>
                ))}
              </div>
            </motion.div>
          )}

          {top10.length === 0 && (
            <motion.div {...FADE_UP(5)} className="card text-center py-12">
              <FlaskConical className="w-10 h-10 text-text-muted/30 mx-auto mb-3" />
              <p className="text-text-muted font-medium">No candidates yet</p>
              <p className="text-sm text-text-muted/60 mt-1">Run the pipeline to generate scored gene candidates</p>
              <Link to="/pipeline" className="btn-primary mt-4 inline-flex">
                <Activity className="w-4 h-4" /> Go to Pipeline
              </Link>
            </motion.div>
          )}
        </>
      )}
    </div>
  )
}
