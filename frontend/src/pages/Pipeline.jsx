import { useEffect, useRef, useState, useCallback } from 'react'
import { motion, AnimatePresence } from 'framer-motion'
import {
  CheckCircle2, Circle, Loader2, AlertCircle, Play, Square,
  Terminal, RefreshCw, ExternalLink, Settings2,
  X, Dna, Shield,
} from 'lucide-react'
import { api } from '@/lib/api'
import { PageHeader, Spinner } from '@/components/ui'
import { cn, relativeTime } from '@/lib/utils'

// ─── Nextflow step1 → step15 (must match nextflow/main.nf STEP_ORDER prefix) ─
/** IDs in DAG execution order through inclusive ``step15``. */
const STEP_ORDER_TO_15 = [
  'step1', 'step2', 'step3', 'step3b', 'step3c',
  'step4', 'step4b', 'step4c', 'step4d',
  'step5', 'step3d', 'step6', 'step6b', 'step6c', 'step7', 'step7b',
  'step8', 'step8b', 'step9', 'step9b',
  'step10b', 'step11', 'step11b', 'step11c', 'step11d',
  'step12', 'step12b', 'step13', 'step14', 'step14b', 'step15',
]

const STEP_LABELS = {
  step1: 'Validate environment & tools',
  step2: 'Download proteomes (cancer-resistance species set)',
  step3: 'OrthoFinder — ortholog groups',
  step3b: 'Load orthologs into database',
  step3c: 'Nucleotide conservation (MAFFT)',
  step4: 'Align proteins & divergent regions',
  step4b: 'Pfam + AlphaMissense consequences',
  step4c: 'ESM-1v variant effect',
  step4d: 'Variant direction (GoF / LoF)',
  step5: 'Species phylogeny (IQ-TREE)',
  step3d: 'Phylo conservation (PHAST / phyloP)',
  step6: 'PAML branch-site selection',
  step6b: 'HyPhy FEL + BUSTED',
  step6c: 'HyPhy RELAX',
  step7: 'Convergence & conservation scoring',
  step7b: 'Convergent amino acids (true positives)',
  step8: 'Expression evidence (Open Targets / DepMap)',
  step8b: 'Bgee tissue expression',
  step9: 'Phase 1 composite scores',
  step9b: 'Structural context (AlphaFold + pockets)',
  step10b: 'Regulatory divergence (AlphaGenome)',
  step11: 'Disease annotation (OpenTargets / GWAS)',
  step11b: 'Rare variants (gnomAD)',
  step11c: 'Literature (PubMed)',
  step11d: 'Pathway convergence',
  step12: 'Druggability (AlphaFold + fpocket)',
  step12b: 'Druggability (P2Rank)',
  step13: 'Gene therapy (AAV / CRISPR)',
  step14: 'Safety screen (PheWAS / STRING)',
  step14b: 'Safety (DepMap + GTEx)',
  step15: 'Final re-score (Phase 2 weights)',
}

/** UI list: same order as Nextflow will run (step1 … step15). */
const PIPELINE_STEPS_TO_15 = STEP_ORDER_TO_15.map((id) => ({
  id,
  label: STEP_LABELS[id] || id,
}))

/** Preset test species for cancer_resistance (config/trait_presets.json). */
const CANCER_RESISTANCE_SPECIES_IDS = [
  'naked_mole_rat',
  'blind_mole_rat',
  'bowhead_whale',
  'african_elephant',
  'greenland_shark',
]

const PHENOTYPE_FIXED = 'cancer_resistance'
const UNTIL_STEP_FIXED = 'step15'

// ─── Step components ───────────────────────────────────────────────────────

function StepIcon({ status }) {
  if (status === 'complete') return <CheckCircle2 className="w-4 h-4 text-bio" strokeWidth={2} />
  if (status === 'running')  return <Loader2 className="w-4 h-4 text-accent animate-spin" />
  if (status === 'failed')   return <AlertCircle className="w-4 h-4 text-danger" />
  return <Circle className="w-4 h-4 text-ink-3/25" strokeWidth={1.5} />
}

function StepRow({ step, info, isLast }) {
  const status = info?.status ?? 'pending'
  return (
    <div className="flex gap-3">
      {/* timeline connector */}
      <div className="flex flex-col items-center">
        <div className={cn(
          'w-7 h-7 rounded-full flex items-center justify-center shrink-0 transition-all duration-300',
          status === 'complete' ? 'bg-bio-bg'        :
          status === 'running'  ? 'bg-accent-light ring-2 ring-accent/20' :
          status === 'failed'   ? 'bg-danger-bg'     : 'bg-canvas',
        )}>
          <StepIcon status={status} />
        </div>
        {!isLast && (
          <div className={cn('w-px flex-1 mt-1', status === 'complete' ? 'bg-bio-ring' : 'bg-border')} />
        )}
      </div>

      {/* content */}
      <div className={cn('pb-4 flex-1 min-w-0', isLast && 'pb-0')}>
        <div className="flex items-center justify-between gap-2 min-h-7">
          <p className={cn(
            'text-sm font-medium transition-colors leading-tight',
            status === 'running'  ? 'text-accent' :
            status === 'complete' ? 'text-ink'    :
            status === 'failed'   ? 'text-danger' : 'text-ink-3',
          )}>
            {step.label}
          </p>
          {info?.elapsed_s != null && (
            <span className="text-[11px] font-mono text-ink-3 shrink-0">{info.elapsed_s}s</span>
          )}
        </div>
        {status === 'running' && (
          <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }} className="flex items-center gap-1.5 mt-0.5">
            <span className="w-1 h-1 rounded-full bg-accent animate-pulse" />
            <p className="text-[11px] text-accent font-medium">Processing…</p>
          </motion.div>
        )}
        {status === 'failed' && (
          <p className="text-[11px] text-danger mt-0.5">Failed — see log below</p>
        )}
      </div>
    </div>
  )
}

// ─── Run Config Drawer (cancer_resistance · step1 → step15 only) ───────────

function RunConfigDrawer({ open, onClose, species, onLaunch, isRunning }) {
  const [resumeFrom, setResumeFrom] = useState('step1')
  const [dryRun, setDryRun] = useState(false)
  const [loading, setLoading] = useState(false)

  const cancerSpecies = species.filter((sp) => CANCER_RESISTANCE_SPECIES_IDS.includes(sp.id))

  const handleLaunch = async () => {
    setLoading(true)
    try {
      await onLaunch({
        resume_from: resumeFrom,
        until_step: UNTIL_STEP_FIXED,
        dry_run: dryRun,
        phenotype: PHENOTYPE_FIXED,
        species_ids: CANCER_RESISTANCE_SPECIES_IDS,
      })
      onClose()
    } catch (e) {
      alert(e.message)
    }
    setLoading(false)
  }

  const startIdx = PIPELINE_STEPS_TO_15.findIndex((s) => s.id === resumeFrom)
  const stepsRemaining = startIdx < 0 ? PIPELINE_STEPS_TO_15.length : PIPELINE_STEPS_TO_15.length - startIdx

  return (
    <AnimatePresence>
      {open && (
        <>
          <motion.div
            initial={{ opacity: 0 }} animate={{ opacity: 1 }} exit={{ opacity: 0 }}
            className="fixed inset-0 bg-ink/20 backdrop-blur-sm z-30"
            onClick={onClose}
          />

          <motion.div
            initial={{ x: '100%' }} animate={{ x: 0 }} exit={{ x: '100%' }}
            transition={{ type: 'spring', damping: 28, stiffness: 280 }}
            className="fixed right-0 top-0 h-full w-[480px] bg-white shadow-card-lg z-40 flex flex-col overflow-hidden"
          >
            <div className="flex items-center justify-between px-6 py-5 border-b border-border">
              <div className="flex items-center gap-3">
                <div className="w-9 h-9 rounded-xl bg-accent-light flex items-center justify-center">
                  <Settings2 className="w-4.5 h-4.5 text-accent" />
                </div>
                <div>
                  <p className="font-bold text-ink text-[15px]">Configure Run</p>
                  <p className="text-xs text-ink-3">Cancer resistance · Nextflow step1 → step15</p>
                </div>
              </div>
              <button type="button" onClick={onClose} className="btn-icon">
                <X className="w-4 h-4" />
              </button>
            </div>

            <div className="flex-1 overflow-y-auto px-6 py-6 space-y-8">
              <section className="flex gap-3 p-4 rounded-xl border border-accent/20 bg-accent-light/40">
                <div className="w-9 h-9 rounded-lg bg-accent text-white flex items-center justify-center shrink-0">
                  <Shield className="w-4 h-4" />
                </div>
                <div>
                  <p className="text-sm font-semibold text-ink">Cancer resistance</p>
                  <p className="text-xs text-ink-3 leading-relaxed mt-1">
                    Species set matches <span className="font-mono">config/trait_presets.json</span> (registry
                    still adds human baseline / controls in the pipeline).
                  </p>
                </div>
              </section>

              {cancerSpecies.length > 0 && (
                <section>
                  <p className="label mb-3">Species in this run</p>
                  <ul className="space-y-2">
                    {cancerSpecies.map((sp) => (
                      <li
                        key={sp.id}
                        className="flex items-center gap-2 text-sm text-ink-2 border border-border rounded-lg px-3 py-2 bg-canvas"
                      >
                        <CheckCircle2 className="w-4 h-4 text-bio shrink-0" strokeWidth={2} />
                        <span className="font-medium text-ink">{sp.common_name ?? sp.id}</span>
                        <span className="text-xs font-mono text-ink-3 ml-auto">{sp.id}</span>
                      </li>
                    ))}
                  </ul>
                </section>
              )}

              <section>
                <p className="label mb-3">Resume from (Nextflow --from_step)</p>
                <p className="text-[11px] text-ink-3 mb-2">Run ends at <span className="font-mono">{UNTIL_STEP_FIXED}</span> inclusive.</p>
                <div className="flex flex-wrap gap-1.5 max-h-48 overflow-y-auto">
                  {PIPELINE_STEPS_TO_15.map((s) => (
                    <button
                      type="button"
                      key={s.id}
                      onClick={() => setResumeFrom(s.id)}
                      className={cn(
                        'px-2.5 py-1 rounded-lg text-xs font-mono font-medium transition-all duration-150',
                        resumeFrom === s.id
                          ? 'bg-accent text-white shadow-accent/20 shadow-sm'
                          : 'bg-canvas text-ink-3 hover:bg-canvas-2',
                      )}
                    >
                      {s.id}
                    </button>
                  ))}
                </div>
              </section>

              <section>
                <label className="flex items-center gap-3 cursor-pointer p-4 rounded-xl border border-border hover:border-border-2 hover:bg-canvas transition-all">
                  <div
                    role="switch"
                    aria-checked={dryRun}
                    tabIndex={0}
                    onClick={() => setDryRun((v) => !v)}
                    onKeyDown={(e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); setDryRun((v) => !v) } }}
                    className={cn(
                      'w-10 h-6 rounded-full flex items-center transition-all duration-200 px-0.5',
                      dryRun ? 'bg-accent' : 'bg-border-2',
                    )}
                  >
                    <div className={cn(
                      'w-5 h-5 rounded-full bg-white shadow-sm transition-transform duration-200',
                      dryRun ? 'translate-x-4' : 'translate-x-0',
                    )} />
                  </div>
                  <div>
                    <p className="text-sm font-semibold text-ink">Dry run mode</p>
                    <p className="text-xs text-ink-3">Validate configuration without executing steps</p>
                  </div>
                </label>
              </section>
            </div>

            <div className="px-6 py-5 border-t border-border bg-canvas">
              <div className="flex items-center justify-between gap-4">
                <div className="text-xs text-ink-3">
                  <span className="font-semibold text-ink">{CANCER_RESISTANCE_SPECIES_IDS.length}</span> focal species ·{' '}
                  <span className="font-semibold text-ink">{stepsRemaining}</span> stages from resume → {UNTIL_STEP_FIXED}
                </div>
                <div className="flex gap-2">
                  <button type="button" onClick={onClose} className="btn-secondary">Cancel</button>
                  <button type="button" onClick={handleLaunch} disabled={loading || isRunning} className="btn-primary">
                    {loading ? <Spinner className="w-4 h-4" /> : <Play className="w-4 h-4" />}
                    {dryRun ? 'Dry Run' : 'Launch Pipeline'}
                  </button>
                </div>
              </div>
            </div>
          </motion.div>
        </>
      )}
    </AnimatePresence>
  )
}

// ─── Main Pipeline Page ────────────────────────────────────────────────────

export default function PipelinePage() {
  const [status, setStatus]             = useState(null)
  const [logs, setLogs]                 = useState([])
  const [loading, setLoading]           = useState(true)
  const [showConfig, setShowConfig]     = useState(false)
  const [actionLoading, setActionLoading] = useState(false)
  const [species, setSpecies]           = useState([])
  const logsRef  = useRef(null)
  const esRef    = useRef(null)

  const fetchStatus = useCallback(async () => {
    try {
      const s = await api.getPipelineStatus()
      setStatus(s)
    } catch {}
    setLoading(false)
  }, [])

  useEffect(() => {
    api.getSpecies().then(setSpecies).catch(() => [])
  }, [])

  // Poll every 5s while running
  useEffect(() => {
    fetchStatus()
    const iv = setInterval(() => {
      if (status?.status === 'running') fetchStatus()
    }, 5000)
    return () => clearInterval(iv)
  }, [fetchStatus, status?.status])

  // SSE log stream
  useEffect(() => {
    esRef.current?.close()
    const es = new EventSource('/pipeline/logs?tail=300')
    esRef.current = es
    es.onmessage = (e) => {
      if (e.data.trim()) {
        setLogs(prev => [...prev.slice(-400), e.data])
        setTimeout(() => {
          if (logsRef.current) logsRef.current.scrollTop = logsRef.current.scrollHeight
        }, 50)
      }
    }
    return () => es.close()
  }, [])

  const handleLaunch = async (params) => {
    setActionLoading(true)
    try {
      await api.startPipeline(params)
      await fetchStatus()
    } catch (e) {
      throw e
    } finally {
      setActionLoading(false)
    }
  }

  const handleStop = async () => {
    if (!confirm('Stop the running pipeline?')) return
    setActionLoading(true)
    try { await api.stopPipeline() } catch {}
    await fetchStatus()
    setActionLoading(false)
  }

  const isRunning = status?.status === 'running'

  const completedSteps = status
    ? PIPELINE_STEPS_TO_15.filter((s) => status.steps?.[s.id]?.status === 'complete').length
    : 0

  return (
    <div className="px-8 py-8 space-y-7">
      <PageHeader
        title="Pipeline Control"
        subtitle="Cancer resistance · Nextflow step1 → step15 (same order as main.nf)"
      >
        {status?.seqera_run_url && (
          <a
            href={status.seqera_run_url}
            target="_blank"
            rel="noopener noreferrer"
            className="btn-secondary text-sm"
          >
            <ExternalLink className="w-4 h-4" />
            Seqera Platform
          </a>
        )}
        {isRunning ? (
          <button onClick={handleStop} disabled={actionLoading} className="btn-danger">
            {actionLoading ? <Spinner className="w-4 h-4" /> : <Square className="w-4 h-4" />}
            Stop Pipeline
          </button>
        ) : (
          <>
            <button onClick={() => setShowConfig(true)} className="btn-secondary">
              <Settings2 className="w-4 h-4" />
              Configure
            </button>
            <button
              type="button"
              onClick={() => handleLaunch({
                resume_from: 'step1',
                until_step: UNTIL_STEP_FIXED,
                dry_run: false,
                phenotype: PHENOTYPE_FIXED,
                species_ids: CANCER_RESISTANCE_SPECIES_IDS,
              })}
              disabled={actionLoading}
              className="btn-primary"
            >
              {actionLoading ? <Spinner className="w-4 h-4" /> : <Play className="w-4 h-4" />}
              Quick Run
            </button>
          </>
        )}
      </PageHeader>

      {/* Status Banner */}
      {status && (
        <motion.div
          initial={{ opacity: 0, y: 8 }} animate={{ opacity: 1, y: 0 }}
          className="card"
        >
          <div className="flex items-center justify-between gap-6">
            <div className="flex items-center gap-4">
              {isRunning ? (
                <div className="flex items-center gap-2 pill-running">
                  <span className="w-2 h-2 rounded-full bg-accent animate-pulse" />
                  Running
                </div>
              ) : status.status === 'complete' ? (
                <div className="pill-complete"><CheckCircle2 className="w-3 h-3" /> Complete</div>
              ) : status.status === 'failed' ? (
                <div className="pill-failed"><AlertCircle className="w-3 h-3" /> Failed</div>
              ) : (
                <div className="pill-idle"><Circle className="w-3 h-3" /> Idle</div>
              )}
              <div>
                <p className="text-sm font-semibold text-ink">
                  {isRunning
                    ? `Running: ${status.current_step ?? 'initialising…'}`
                    : status.status === 'complete'
                    ? `All stages through ${UNTIL_STEP_FIXED} complete (${PIPELINE_STEPS_TO_15.length} steps)`
                    : status.status === 'failed'
                    ? `Failed at ${status.current_step}`
                    : 'Ready to run'}
                </p>
                {status.started_at && (
                  <p className="text-xs text-ink-3 mt-0.5">
                    Started {relativeTime(status.started_at)}
                    {status.elapsed_total_s && ` · ${(status.elapsed_total_s / 60).toFixed(1)} min elapsed`}
                  </p>
                )}
              </div>
            </div>

            {/* Progress bar */}
            {completedSteps > 0 && (
              <div className="hidden md:flex items-center gap-3 flex-1 max-w-xs">
                <div className="flex-1 score-bar">
                  <div className="score-bar-fill" style={{ width: `${(completedSteps / PIPELINE_STEPS_TO_15.length) * 100}%` }} />
                </div>
                <span className="text-xs font-mono text-ink-3 shrink-0">
                  {completedSteps}/{PIPELINE_STEPS_TO_15.length}
                </span>
              </div>
            )}

            <button onClick={fetchStatus} className="btn-icon ml-auto shrink-0">
              <RefreshCw className="w-3.5 h-3.5" />
            </button>
          </div>
        </motion.div>
      )}

      {/* Main grid */}
      <div className="grid grid-cols-1 xl:grid-cols-5 gap-6">
        {/* ── Step Timelines (left) ── */}
        <div className="xl:col-span-2 space-y-5">
          <div className="card">
            <div className="flex items-center gap-3 mb-5">
              <div className="w-8 h-8 rounded-lg bg-accent-light flex items-center justify-center">
                <Dna className="w-4 h-4 text-accent" />
              </div>
              <div>
                <p className="font-bold text-ink text-[15px]">Active pipeline</p>
                <p className="text-xs text-ink-3">
                  Order matches <span className="font-mono">nextflow/main.nf</span> through{' '}
                  <span className="font-mono">{UNTIL_STEP_FIXED}</span>
                </p>
              </div>
              <div className="ml-auto text-xs font-mono text-ink-3">
                {completedSteps}/{PIPELINE_STEPS_TO_15.length}
              </div>
            </div>
            {loading ? (
              <Spinner className="mx-auto my-4" />
            ) : (
              PIPELINE_STEPS_TO_15.map((step, i) => (
                <StepRow
                  key={step.id}
                  step={step}
                  info={status?.steps?.[step.id]}
                  isLast={i === PIPELINE_STEPS_TO_15.length - 1}
                />
              ))
            )}
          </div>
        </div>

        {/* ── Log Viewer (right) ── */}
        <div className="xl:col-span-3 card flex flex-col" style={{ minHeight: '82vh' }}>
          <div className="flex items-center gap-3 mb-4 pb-4 border-b border-border">
            <div className="w-8 h-8 rounded-lg bg-canvas flex items-center justify-center">
              <Terminal className="w-4 h-4 text-ink-2" />
            </div>
            <div>
              <p className="font-bold text-ink text-[14px]">Live Log</p>
              <p className="text-xs text-ink-3">Streaming output from pipeline executor</p>
            </div>
            <div className="ml-auto flex items-center gap-3">
              {status?.seqera_run_url && (
                <a
                  href={status.seqera_run_url}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="flex items-center gap-1.5 text-xs text-ink-3 hover:text-ink transition-colors"
                >
                  <ExternalLink className="w-3 h-3" />
                  Seqera
                </a>
              )}
              {isRunning && (
                <div className="flex items-center gap-1.5 text-xs text-accent font-semibold">
                  <span className="w-2 h-2 rounded-full bg-accent animate-pulse" />
                  Live
                </div>
              )}
            </div>
          </div>

          <div
            ref={logsRef}
            className="flex-1 overflow-y-auto rounded-xl bg-canvas border border-border p-4 log-output text-ink-2 leading-relaxed"
            style={{ maxHeight: '75vh' }}
          >
            {logs.length === 0 ? (
              <div className="flex flex-col items-center justify-center h-full text-center py-12">
                <Terminal className="w-8 h-8 text-ink-3/30 mb-3" />
                <p className="text-sm text-ink-3 font-medium">
                  {isRunning ? 'Waiting for output…' : 'No output yet'}
                </p>
                <p className="text-xs text-ink-3/60 mt-1">Start the pipeline to see live logs here</p>
              </div>
            ) : (
              logs.map((line, i) => (
                <p
                  key={i}
                  className={cn(
                    'whitespace-pre-wrap break-all',
                    line.includes('ERROR') || line.includes('FAILED') ? 'text-danger' :
                    line.includes('WARNING') ? 'text-warning' :
                    line.includes('complete') || line.includes('✓') || line.includes('SUCCESS') ? 'text-bio' :
                    line.includes('Step ') || line.includes('Running') ? 'text-accent' :
                    'text-ink-2',
                  )}
                >
                  {line}
                </p>
              ))
            )}
          </div>
        </div>
      </div>

      {/* Config Drawer (Sprint 1 - phenotype + species) */}
      <RunConfigDrawer
        open={showConfig}
        onClose={() => setShowConfig(false)}
        species={species}
        onLaunch={handleLaunch}
        isRunning={isRunning}
      />
    </div>
  )
}
