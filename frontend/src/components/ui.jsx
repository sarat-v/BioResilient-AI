import { cn, tierBadgeClass } from '@/lib/utils'

export function TierBadge({ tier }) {
  if (!tier) return null
  return <span className={tierBadgeClass(tier)}>{tier}</span>
}

export function ScoreBar({ value, className }) {
  const pct = Math.round((value ?? 0) * 100)
  return (
    <div className={cn('flex items-center gap-2', className)}>
      <div className="flex-1 h-1.5 bg-white/5 rounded-full overflow-hidden">
        <div className="score-bar-fill h-full" style={{ width: `${pct}%` }} />
      </div>
      <span className="text-xs font-mono text-text-muted w-8 text-right">{pct}</span>
    </div>
  )
}

export function StatCard({ label, value, sub, accent }) {
  return (
    <div className="card flex flex-col gap-1">
      <p className="label-muted">{label}</p>
      <p className={cn('text-3xl font-bold', accent ? 'text-accent' : 'text-text-primary')}>
        {value ?? '—'}
      </p>
      {sub && <p className="text-xs text-text-muted">{sub}</p>}
    </div>
  )
}

export function PageHeader({ title, subtitle, children }) {
  return (
    <div className="flex items-start justify-between px-8 pt-8 pb-4">
      <div>
        <h1 className="text-2xl font-bold text-text-primary">{title}</h1>
        {subtitle && <p className="text-sm text-text-muted mt-0.5">{subtitle}</p>}
      </div>
      {children && <div className="flex items-center gap-2">{children}</div>}
    </div>
  )
}

export function Spinner({ className }) {
  return (
    <div className={cn('w-5 h-5 rounded-full border-2 border-accent/20 border-t-accent animate-spin', className)} />
  )
}

export function EmptyState({ icon: Icon, title, subtitle }) {
  return (
    <div className="flex flex-col items-center justify-center gap-3 py-20 text-center">
      {Icon && <Icon className="w-10 h-10 text-text-muted/40" />}
      <p className="text-base font-medium text-text-muted">{title}</p>
      {subtitle && <p className="text-sm text-text-muted/60">{subtitle}</p>}
    </div>
  )
}
