import { useEffect, useState } from 'react'
import { Link } from 'react-router-dom'
import { motion } from 'framer-motion'
import {
  BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, Cell,
} from 'recharts'
import { Route, ExternalLink } from 'lucide-react'
import { api } from '@/lib/api'
import { PageHeader, Spinner, EmptyState } from '@/components/ui'
import { cn, formatFloat } from '@/lib/utils'

const ENRICHMENT_COLOR = (score) => {
  if (score >= 0.7) return '#7c3aed'
  if (score >= 0.4) return '#f59e0b'
  return '#475569'
}

export default function PathwayConvergencePage() {
  const [pathways, setPathways] = useState([])
  const [loading, setLoading] = useState(true)
  const [minCandidates, setMinCandidates] = useState(2)

  useEffect(() => {
    setLoading(true)
    api.getPathwayConvergence({ top: 100, min_candidates: minCandidates })
      .then(d => { setPathways(d ?? []); setLoading(false) })
      .catch(() => setLoading(false))
  }, [minCandidates])

  const chartData = pathways.slice(0, 20).map(p => ({
    name: (p.pathway_name || p.pathway_id || '').slice(0, 35),
    score: p.pathway_score ?? 0,
    candidates: p.candidate_count ?? 0,
  }))

  return (
    <div className="px-8 py-8 space-y-6">
      <PageHeader
        title="Pathway Convergence"
        subtitle="Biological pathways enriched for convergent candidate genes, scored by hypergeometric enrichment weighted by evolutionary signal"
      />

      <div className="flex items-center gap-3">
        <span className="text-xs text-ink-3">Min candidates per pathway</span>
        <input
          type="range" min={1} max={10} step={1} value={minCandidates}
          onChange={e => setMinCandidates(Number(e.target.value))}
          className="w-24 accent-accent"
        />
        <span className="font-mono text-xs text-accent w-4">{minCandidates}</span>
      </div>

      {loading ? (
        <div className="flex items-center justify-center py-32"><Spinner className="w-8 h-8" /></div>
      ) : pathways.length === 0 ? (
        <EmptyState icon={Route} title="No pathway data" subtitle="Run the pipeline through step 11d to generate pathway convergence scores." />
      ) : (
        <>
          {/* Top pathways chart */}
          <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }} className="card">
            <p className="font-semibold text-ink mb-2">Top Enriched Pathways</p>
            <p className="text-xs text-ink-3 mb-4">Ranked by combined enrichment + evolutionary weight score</p>
            <ResponsiveContainer width="100%" height={Math.max(320, chartData.length * 32)}>
              <BarChart data={chartData} layout="vertical" margin={{ left: 250, right: 20, top: 10, bottom: 10 }}>
                <XAxis type="number" domain={[0, 'auto']}
                  tick={{ fill: '#64748b', fontSize: 12 }} axisLine={false} tickLine={false} />
                <YAxis type="category" dataKey="name"
                  tick={{ fill: '#1a1f2e', fontSize: 12, fontWeight: 500 }}
                  width={240} axisLine={false} tickLine={false} />
                <Tooltip
                  cursor={{ fill: 'rgba(255,255,255,0.05)' }}
                  contentStyle={{ background: '#1a2235', border: '1px solid rgba(255,255,255,0.12)', borderRadius: 8, fontSize: 13, padding: '8px 12px' }}
                  formatter={(v, name) => [v.toFixed(3), name === 'score' ? 'Pathway Score' : 'Candidates']}
                  labelStyle={{ color: '#f1f5f9' }}
                />
                <Bar dataKey="score" radius={[0, 6, 6, 0]} isAnimationActive={true}>
                  {chartData.map((d, i) => (
                    <Cell key={i} fill={ENRICHMENT_COLOR(d.score)} />
                  ))}
                </Bar>
              </BarChart>
            </ResponsiveContainer>
          </motion.div>

          {/* Full pathway table */}
          <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }} transition={{ delay: 0.1 }} className="card overflow-x-auto p-0">
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b border-border">
                  <th className="px-4 py-3 text-left label font-medium">#</th>
                  <th className="px-4 py-3 text-left label font-medium">Pathway</th>
                  <th className="px-4 py-3 text-left label font-medium">ID</th>
                  <th className="px-4 py-3 text-left label font-medium">Candidates</th>
                  <th className="px-4 py-3 text-left label font-medium">Background</th>
                  <th className="px-4 py-3 text-left label font-medium">log10 p-value</th>
                  <th className="px-4 py-3 text-left label font-medium">Evol. Weight</th>
                  <th className="px-4 py-3 text-left label font-medium">Score</th>
                  <th className="px-4 py-3 text-left label font-medium">Genes</th>
                </tr>
              </thead>
              <tbody>
                {pathways.map((p, i) => (
                  <tr key={p.pathway_id} className={cn(
                    'border-b border-white/3 transition-colors hover:bg-surface-2',
                    i % 2 === 0 ? 'bg-surface' : 'bg-canvas/40',
                  )}>
                    <td className="px-4 py-3 text-ink-3 text-xs">{i + 1}</td>
                    <td className="px-4 py-3 text-ink font-medium">
                      {p.pathway_name || p.pathway_id}
                    </td>
                    <td className="px-4 py-3">
                      {p.pathway_id?.startsWith('R-') ? (
                        <a href={`https://reactome.org/content/detail/${p.pathway_id}`}
                          target="_blank" rel="noreferrer"
                          className="text-xs text-ink-3 hover:text-accent font-mono inline-flex items-center gap-1">
                          {p.pathway_id} <ExternalLink className="w-3 h-3" />
                        </a>
                      ) : (
                        <span className="text-xs text-ink-3 font-mono">{p.pathway_id}</span>
                      )}
                    </td>
                    <td className="px-4 py-3 font-mono text-xs text-accent">{p.candidate_count}</td>
                    <td className="px-4 py-3 font-mono text-xs text-ink-3">{p.gene_count}</td>
                    <td className="px-4 py-3 font-mono text-xs text-ink-3">{formatFloat(p.log_pvalue, 2)}</td>
                    <td className="px-4 py-3 font-mono text-xs text-ink-3">{formatFloat(p.evolutionary_weight, 2)}</td>
                    <td className="px-4 py-3 font-mono text-xs font-medium text-accent">{formatFloat(p.pathway_score, 3)}</td>
                    <td className="px-4 py-3">
                      <div className="flex flex-wrap gap-1">
                        {(p.gene_symbols || []).slice(0, 8).map(s => (
                          <span key={s} className="px-1.5 py-0.5 rounded text-[10px] font-mono bg-canvas text-ink-3 border border-border">
                            {s}
                          </span>
                        ))}
                        {(p.gene_symbols || []).length > 8 && (
                          <span className="text-[10px] text-ink-3">+{p.gene_symbols.length - 8}</span>
                        )}
                      </div>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </motion.div>
        </>
      )}
    </div>
  )
}
