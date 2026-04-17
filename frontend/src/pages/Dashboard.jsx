import { useEffect, useState } from 'react'
import { Link } from 'react-router-dom'
import { motion } from 'framer-motion'
import {
  BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, Cell,
} from 'recharts'
import {
  FlaskConical, Activity, ArrowRight, CheckCircle2,
  AlertCircle, Clock, TrendingUp, Microscope, Dna, Shield,
} from 'lucide-react'
import { api } from '@/lib/api'
import { StatCard, PageHeader, Spinner, TierBadge, ScoreBar } from '@/components/ui'
import { relativeTime, formatFloat } from '@/lib/utils'

const stagger = (i) => ({
  initial: { opacity: 0, y: 18 },
  animate: { opacity: 1, y: 0 },
  transition: { delay: i * 0.06, duration: 0.38, ease: [0.22, 1, 0.36, 1] },
})

const TIER_COLORS = {
  Validated: '#047857',
  Tier1:     '#7C3AED',
  Tier2:     '#B45309',
  Tier3:     '#94A3B8',
}

function StatusPill({ status }) {
  const map = {
    idle:     { label: 'Not started', icon: Clock,        cls: 'pill-idle' },
    running:  { label: 'Running',     icon: Activity,     cls: 'pill-running' },
    complete: { label: 'Complete',    icon: CheckCircle2, cls: 'pill-complete' },
    failed:   { label: 'Failed',      icon: AlertCircle,  cls: 'pill-failed' },
    stopped:  { label: 'Stopped',     icon: AlertCircle,  cls: 'pill-idle' },
  }
  const { label, icon: Icon, cls } = map[status] ?? map.idle
  return (
    <span className={cls}>
      <Icon className="w-3 h-3" />
      {label}
    </span>
  )
}

function CustomTooltip({ active, payload }) {
  if (!active || !payload?.length) return null
  const d = payload[0]
  return (
    <div className="bg-white rounded-xl shadow-card-lg border border-border px-4 py-3">
      <p className="text-xs text-ink-3 mb-0.5">{d.payload.gene_symbol}</p>
      <p className="text-sm font-bold text-ink">{(d.value * 100).toFixed(1)}</p>
      <TierBadge tier={d.payload.tier} />
    </div>
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
  const tier1     = candidates?.filter(c => c.tier === 'Tier1').length ?? 0
  const tier2     = candidates?.filter(c => c.tier === 'Tier2').length ?? 0
  const tier3     = candidates?.filter(c => c.tier === 'Tier3').length ?? 0
  const top10     = candidates?.slice(0, 10) ?? []

  const completedSteps = pipelineStatus
    ? Object.values(pipelineStatus.steps || {}).filter(s => s.status === 'complete').length
    : 0
  const totalSteps = pipelineStatus ? (Object.keys(pipelineStatus.steps || {}).length || 30) : 30

  return (
    <div className="px-8 py-8 space-y-7">
      <PageHeader title="Dashboard" subtitle="Your BioResilient analysis at a glance">
        {traits.length > 0 && (
          <select
            value={traitId}
            onChange={e => setTraitId(e.target.value)}
            className="select w-52"
          >
            <option value="">Default analysis</option>
            {traits.map(t => <option key={t.id} value={t.id}>{t.label}</option>)}
          </select>
        )}
        <Link to="/pipeline" className="btn-secondary">
          <Activity className="w-4 h-4" /> Pipeline
        </Link>
        <Link to="/candidates" className="btn-primary">
          <FlaskConical className="w-4 h-4" /> View Candidates
        </Link>
      </PageHeader>

      {loading ? (
        <div className="flex items-center justify-center py-40">
          <Spinner className="w-8 h-8" />
        </div>
      ) : (
        <>
          {/* ── Stats Row ── */}
          <div className="grid grid-cols-2 lg:grid-cols-4 xl:grid-cols-5 gap-4">
            {[
              {
                label: 'Validated Targets',
                value: validated,
                sub: 'Evolutionary + human genetic support',
                accent: true,
                icon: Shield,
                color: 'text-bio',
              },
              {
                label: 'Tier 1 Candidates',
                value: tier1,
                sub: 'High-confidence, ready to test',
                icon: TrendingUp,
                color: 'text-violet',
              },
              {
                label: 'Tier 2 Candidates',
                value: tier2,
                sub: 'Promising, needs validation',
                icon: Microscope,
                color: 'text-amber',
              },
              {
                label: 'Tier 3',
                value: tier3,
                sub: 'Lower priority targets',
                icon: FlaskConical,
              },
              {
                label: 'Pipeline Progress',
                value: `${completedSteps}/${totalSteps}`,
                sub: 'Steps complete',
                icon: Activity,
                color: 'text-accent',
              },
            ].map((card, i) => (
              <motion.div key={card.label} {...stagger(i)}>
                <StatCard {...card} />
              </motion.div>
            ))}
          </div>

          {/* ── Pipeline Status Banner ── */}
          {pipelineStatus && (
            <motion.div {...stagger(4)} className="card">
              <div className="flex items-center justify-between gap-4">
                <div className="flex items-center gap-4">
                  <StatusPill status={pipelineStatus.status} />
                  <div>
                    <p className="text-sm font-semibold text-ink">
                      {pipelineStatus.status === 'running'
                        ? `Running: ${pipelineStatus.current_step ?? 'initialising…'}`
                        : pipelineStatus.status === 'complete'
                        ? 'All analysis complete — results ready'
                        : pipelineStatus.status === 'failed'
                        ? `Failed at: ${pipelineStatus.current_step}`
                        : 'No pipeline run in progress'}
                    </p>
                    {pipelineStatus.started_at && (
                      <p className="text-xs text-ink-3 mt-0.5">
                        Started {relativeTime(pipelineStatus.started_at)}
                        {pipelineStatus.elapsed_total_s &&
                          ` · ${(pipelineStatus.elapsed_total_s / 60).toFixed(1)} min total`}
                      </p>
                    )}
                  </div>
                </div>
                {/* Progress bar */}
                <div className="hidden md:flex items-center gap-3 flex-1 max-w-xs">
                  <div className="flex-1 score-bar">
                    <div
                      className="score-bar-fill"
                      style={{ width: `${totalSteps ? (completedSteps / totalSteps) * 100 : 0}%` }}
                    />
                  </div>
                  <span className="text-xs font-mono text-ink-3 shrink-0">
                    {totalSteps ? Math.round((completedSteps / totalSteps) * 100) : 0}%
                  </span>
                </div>
                <Link to="/pipeline" className="btn-ghost shrink-0">
                  Details <ArrowRight className="w-4 h-4" />
                </Link>
              </div>
            </motion.div>
          )}

          {top10.length > 0 && (
            <div className="grid grid-cols-1 xl:grid-cols-5 gap-6">
              {/* ── Bar Chart ── */}
              <motion.div {...stagger(5)} className="card xl:col-span-3">
                <div className="flex items-center justify-between mb-6">
                  <div>
                    <p className="font-bold text-ink text-[15px]">Top Candidates</p>
                    <p className="text-xs text-ink-3 mt-0.5">Composite score · Phase 1 + 2</p>
                  </div>
                  <Link to="/candidates" className="btn-ghost text-xs">
                    View all <ArrowRight className="w-3.5 h-3.5" />
                  </Link>
                </div>
                <ResponsiveContainer width="100%" height={260}>
                  <BarChart data={top10} layout="vertical" margin={{ left: 4, right: 12, top: 0, bottom: 0 }}>
                    <XAxis
                      type="number" domain={[0, 1]}
                      tickFormatter={v => `${(v * 100).toFixed(0)}`}
                      tick={{ fill: '#9E9890', fontSize: 11, fontFamily: 'Plus Jakarta Sans' }}
                      axisLine={false} tickLine={false}
                    />
                    <YAxis
                      type="category" dataKey="gene_symbol" width={72}
                      tick={{ fill: '#17140F', fontSize: 12, fontFamily: 'JetBrains Mono', fontWeight: 500 }}
                      axisLine={false} tickLine={false}
                    />
                    <Tooltip content={<CustomTooltip />} cursor={{ fill: 'rgba(27,69,212,0.04)', radius: 6 }} />
                    <Bar dataKey="composite_score" radius={[0, 6, 6, 0]} maxBarSize={18}>
                      {top10.map((c) => (
                        <Cell key={c.gene_id} fill={TIER_COLORS[c.tier] ?? '#94A3B8'} />
                      ))}
                    </Bar>
                  </BarChart>
                </ResponsiveContainer>

                {/* Legend */}
                <div className="flex items-center gap-5 mt-4 pt-4 border-t border-border">
                  {Object.entries(TIER_COLORS).map(([tier, color]) => (
                    <div key={tier} className="flex items-center gap-1.5">
                      <div className="w-2.5 h-2.5 rounded-sm" style={{ background: color }} />
                      <span className="text-xs text-ink-3">{tier}</span>
                    </div>
                  ))}
                </div>
              </motion.div>

              {/* ── Quick target list ── */}
              <motion.div {...stagger(6)} className="card xl:col-span-2">
                <div className="flex items-center justify-between mb-5">
                  <p className="font-bold text-ink text-[15px]">Top 5 Targets</p>
                  <Link to="/candidates" className="btn-ghost text-xs">
                    All <ArrowRight className="w-3.5 h-3.5" />
                  </Link>
                </div>
                <div className="space-y-1">
                  {top10.slice(0, 5).map((c, i) => (
                    <Link
                      key={c.gene_id}
                      to={`/candidates/${c.gene_id}`}
                      className="flex items-center justify-between px-3 py-3 rounded-xl hover:bg-canvas transition-colors group"
                    >
                      <div className="flex items-center gap-3">
                        <span className="text-[11px] font-bold text-ink-3 w-5 text-center">{i + 1}</span>
                        <div>
                          <p className="font-mono text-sm font-semibold text-ink group-hover:text-accent transition-colors">
                            {c.gene_symbol}
                          </p>
                          <TierBadge tier={c.tier} />
                        </div>
                      </div>
                      <div className="flex items-center gap-2">
                        <div className="w-20">
                          <ScoreBar value={c.composite_score} />
                        </div>
                        <ArrowRight className="w-3.5 h-3.5 text-ink-3 opacity-0 group-hover:opacity-100 transition-opacity" />
                      </div>
                    </Link>
                  ))}
                </div>
              </motion.div>
            </div>
          )}

          {top10.length === 0 && !loading && (
            <motion.div {...stagger(5)}>
              <div className="card text-center py-16">
                <div className="w-16 h-16 rounded-2xl bg-accent-light flex items-center justify-center mx-auto mb-4">
                  <Dna className="w-7 h-7 text-accent" strokeWidth={1.5} />
                </div>
                <p className="text-[17px] font-bold text-ink">No candidates yet</p>
                <p className="text-sm text-ink-3 mt-1.5 max-w-sm mx-auto">
                  Run the analysis pipeline to generate scored gene candidates across species.
                </p>
                <Link to="/pipeline" className="btn-primary mt-5 inline-flex">
                  <Activity className="w-4 h-4" /> Go to Pipeline
                </Link>
              </div>
            </motion.div>
          )}
        </>
      )}
    </div>
  )
}
