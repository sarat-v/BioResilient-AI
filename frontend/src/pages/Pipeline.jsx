import { useEffect, useRef, useState, useCallback } from 'react'
import { motion, AnimatePresence } from 'framer-motion'
import {
  CheckCircle2, Circle, Loader2, AlertCircle, Play, Square,
  ChevronDown, Terminal, RefreshCw,
} from 'lucide-react'
import { api } from '@/lib/api'
import { PageHeader, Spinner } from '@/components/ui'
import { cn, relativeTime } from '@/lib/utils'

const ALL_STEPS = [
  { id: 'step1',  label: 'Validate environment',                    phase: 1 },
  { id: 'step2',  label: 'Download proteomes',                      phase: 1 },
  { id: 'step3',  label: 'Find orthologs across species',           phase: 1 },
  { id: 'step3b', label: 'Load orthologs into database',            phase: 1 },
  { id: 'step4',  label: 'Align sequences & find divergent regions',phase: 1 },
  { id: 'step5',  label: 'Build evolutionary tree',                 phase: 1 },
  { id: 'step6',  label: 'Detect evolutionary selection pressure',  phase: 1 },
  { id: 'step7',  label: 'Score convergence & conservation',        phase: 1 },
  { id: 'step8',  label: 'Analyse gene expression differences',     phase: 1 },
  { id: 'step9',  label: 'Compute composite candidate scores',      phase: 1 },
  { id: 'step10', label: 'API ready',                               phase: 1 },
  { id: 'step10b',label: 'Regulatory divergence (AlphaGenome)',     phase: 2 },
  { id: 'step11', label: 'Disease association lookup',              phase: 2 },
  { id: 'step12', label: 'Druggability assessment',                 phase: 2 },
  { id: 'step13', label: 'Gene therapy feasibility',                phase: 2 },
  { id: 'step14', label: 'Safety pre-screen',                       phase: 2 },
  { id: 'step15', label: 'Final re-score with all data',            phase: 2 },
  { id: 'step16', label: 'Pipeline complete',                       phase: 2 },
]

function StepIcon({ status }) {
  if (status === 'complete') return <CheckCircle2 className="w-4 h-4 text-success" />
  if (status === 'running')  return <Loader2 className="w-4 h-4 text-accent animate-spin" />
  if (status === 'failed')   return <AlertCircle className="w-4 h-4 text-danger" />
  return <Circle className="w-4 h-4 text-white/10" />
}

function StepRow({ step, info, isLast }) {
  const status = info?.status ?? 'pending'
  return (
    <div className="flex gap-3">
      {/* connector line */}
      <div className="flex flex-col items-center">
        <div className={cn(
          'w-8 h-8 rounded-full flex items-center justify-center shrink-0 transition-colors duration-300',
          status === 'complete' ? 'bg-success/10' :
          status === 'running'  ? 'bg-accent/10 ring-1 ring-accent/30' :
          status === 'failed'   ? 'bg-danger/10' : 'bg-white/5',
        )}>
          <StepIcon status={status} />
        </div>
        {!isLast && <div className={cn('w-px flex-1 mt-1', status === 'complete' ? 'bg-success/30' : 'bg-white/5')} />}
      </div>

      {/* content */}
      <div className={cn('pb-4 flex-1', isLast && 'pb-0')}>
        <div className="flex items-center justify-between min-h-8">
          <p className={cn(
            'text-sm font-medium transition-colors',
            status === 'running'  ? 'text-accent' :
            status === 'complete' ? 'text-text-primary' :
            status === 'failed'   ? 'text-danger' : 'text-text-muted',
          )}>
            {step.label}
          </p>
          {info?.elapsed_s != null && (
            <span className="text-xs font-mono text-text-muted">{info.elapsed_s}s</span>
          )}
        </div>
        {status === 'running' && (
          <motion.p initial={{ opacity: 0 }} animate={{ opacity: 1 }} className="text-xs text-accent/60 mt-0.5">
            In progress…
          </motion.p>
        )}
        {status === 'failed' && (
          <p className="text-xs text-danger/70 mt-0.5">Step failed — check logs below</p>
        )}
      </div>
    </div>
  )
}

export default function PipelinePage() {
  const [status, setStatus] = useState(null)
  const [logs, setLogs] = useState([])
  const [loading, setLoading] = useState(true)
  const [resumeFrom, setResumeFrom] = useState('step1')
  const [dryRun, setDryRun] = useState(false)
  const [showResume, setShowResume] = useState(false)
  const [actionLoading, setActionLoading] = useState(false)
  const logsRef = useRef(null)
  const esRef = useRef(null)

  const fetchStatus = useCallback(async () => {
    try {
      const s = await api.getPipelineStatus()
      setStatus(s)
    } catch {}
    setLoading(false)
  }, [])

  // Poll status every 3s while running
  useEffect(() => {
    fetchStatus()
    const interval = setInterval(() => {
      if (status?.status === 'running') fetchStatus()
    }, 3000)
    return () => clearInterval(interval)
  }, [fetchStatus, status?.status])

  // SSE log stream
  useEffect(() => {
    if (esRef.current) esRef.current.close()
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

  const handleRun = async () => {
    setActionLoading(true)
    try {
      await api.startPipeline({ resume_from: resumeFrom, dry_run: dryRun })
      await fetchStatus()
    } catch (e) {
      alert(e.message)
    }
    setActionLoading(false)
    setShowResume(false)
  }

  const handleStop = async () => {
    if (!confirm('Stop the running pipeline?')) return
    setActionLoading(true)
    try { await api.stopPipeline() } catch {}
    await fetchStatus()
    setActionLoading(false)
  }

  const isRunning = status?.status === 'running'

  const phase1Steps = ALL_STEPS.filter(s => s.phase === 1)
  const phase2Steps = ALL_STEPS.filter(s => s.phase === 2)

  return (
    <div className="px-8 py-8 space-y-8">
      <PageHeader
        title="Pipeline Control"
        subtitle="Run the analysis pipeline and monitor progress"
      >
        {isRunning ? (
          <button onClick={handleStop} disabled={actionLoading} className="btn-secondary text-sm text-danger border-danger/40">
            {actionLoading ? <Spinner className="w-4 h-4" /> : <Square className="w-4 h-4" />} Stop
          </button>
        ) : (
          <>
            <button onClick={() => setShowResume(s => !s)} className="btn-ghost text-sm">
              <ChevronDown className={cn('w-4 h-4 transition-transform', showResume && 'rotate-180')} />
              Resume from…
            </button>
            <button onClick={handleRun} disabled={actionLoading} className="btn-primary text-sm">
              {actionLoading ? <Spinner className="w-4 h-4" /> : <Play className="w-4 h-4" />}
              Run Full Pipeline
            </button>
          </>
        )}
      </PageHeader>

      {/* Resume options */}
      <AnimatePresence>
        {showResume && !isRunning && (
          <motion.div
            initial={{ opacity: 0, height: 0 }} animate={{ opacity: 1, height: 'auto' }} exit={{ opacity: 0, height: 0 }}
            className="mx-8 card overflow-hidden"
          >
            <p className="label-muted mb-3">Resume from step</p>
            <div className="flex flex-wrap gap-2 mb-4">
              {ALL_STEPS.map(s => (
                <button
                  key={s.id}
                  onClick={() => setResumeFrom(s.id)}
                  className={cn(
                    'px-3 py-1.5 rounded-lg text-xs font-medium transition-colors',
                    resumeFrom === s.id ? 'bg-accent text-base' : 'bg-white/5 text-text-muted hover:bg-white/10',
                  )}
                >
                  {s.id}
                </button>
              ))}
            </div>
            <div className="flex items-center gap-4">
              <label className="flex items-center gap-2 text-sm text-text-muted cursor-pointer">
                <input type="checkbox" checked={dryRun} onChange={e => setDryRun(e.target.checked)}
                  className="rounded border-white/20 bg-white/5" />
                Dry run (no actual execution)
              </label>
              <button onClick={handleRun} disabled={actionLoading} className="btn-primary text-sm ml-auto">
                <Play className="w-4 h-4" /> Resume from {resumeFrom}
              </button>
            </div>
          </motion.div>
        )}
      </AnimatePresence>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 px-0">
        {/* Step timeline */}
        <div className="space-y-4">
          {/* Phase 1 */}
          <div className="card">
            <div className="flex items-center justify-between mb-5">
              <div>
                <p className="font-semibold text-text-primary">Phase 1 — Core Analysis</p>
                <p className="text-xs text-text-muted">Sequence divergence & evolutionary selection</p>
              </div>
              <button onClick={fetchStatus} className="btn-ghost p-2">
                <RefreshCw className="w-3.5 h-3.5" />
              </button>
            </div>
            {loading ? <Spinner className="mx-auto" /> : (
              <div>
                {phase1Steps.map((step, i) => (
                  <StepRow
                    key={step.id}
                    step={step}
                    info={status?.steps?.[step.id]}
                    isLast={i === phase1Steps.length - 1}
                  />
                ))}
              </div>
            )}
          </div>

          {/* Phase 2 */}
          <div className="card">
            <div className="mb-5">
              <p className="font-semibold text-text-primary">Phase 2 — Deep Annotation</p>
              <p className="text-xs text-text-muted">Disease, druggability, therapy & safety</p>
            </div>
            {loading ? <Spinner className="mx-auto" /> : (
              <div>
                {phase2Steps.map((step, i) => (
                  <StepRow
                    key={step.id}
                    step={step}
                    info={status?.steps?.[step.id]}
                    isLast={i === phase2Steps.length - 1}
                  />
                ))}
              </div>
            )}
          </div>
        </div>

        {/* Live log viewer */}
        <div className="card flex flex-col" style={{ minHeight: '70vh' }}>
          <div className="flex items-center gap-2 mb-3">
            <Terminal className="w-4 h-4 text-accent" />
            <p className="font-medium text-sm text-text-primary">Live Log</p>
            {isRunning && (
              <span className="ml-auto flex items-center gap-1.5 text-xs text-accent">
                <span className="w-1.5 h-1.5 rounded-full bg-accent animate-pulse" />
                Live
              </span>
            )}
          </div>
          <div
            ref={logsRef}
            className="flex-1 overflow-y-auto rounded-lg bg-base p-3 font-mono text-xs text-text-muted leading-relaxed space-y-0.5"
            style={{ maxHeight: '65vh' }}
          >
            {logs.length === 0 ? (
              <p className="text-text-muted/40 text-center py-10">
                {isRunning ? 'Waiting for output…' : 'No log output yet. Start the pipeline to see logs here.'}
              </p>
            ) : (
              logs.map((line, i) => (
                <p
                  key={i}
                  className={cn(
                    'whitespace-pre-wrap break-all',
                    line.includes('ERROR') || line.includes('FAILED') ? 'text-danger' :
                    line.includes('WARNING') ? 'text-warning' :
                    line.includes('complete') || line.includes('✓') ? 'text-success' :
                    line.includes('Step ') ? 'text-accent' : '',
                  )}
                >
                  {line}
                </p>
              ))
            )}
          </div>
        </div>
      </div>
    </div>
  )
}
