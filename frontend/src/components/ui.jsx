import { cn } from '@/lib/utils'

export function TierBadge({ tier }) {
  if (!tier) return null
  const map = {
    Validated: 'badge-validated',
    Tier1:     'badge-tier1',
    Tier2:     'badge-tier2',
    Tier3:     'badge-tier3',
  }
  return <span className={cn('badge', map[tier] ?? 'badge-tier3')}>{tier}</span>
}

export function ScoreBar({ value, className }) {
  const pct = Math.round((value ?? 0) * 100)
  return (
    <div className={cn('flex items-center gap-2.5', className)}>
      <div className="flex-1 score-bar">
        <div className="score-bar-fill" style={{ width: `${pct}%` }} />
      </div>
      <span className="text-xs font-mono text-ink-3 w-7 text-right tabular-nums">{pct}</span>
    </div>
  )
}

export function StatCard({ label, value, sub, accent, color }) {
  return (
    <div className="card flex flex-col gap-1.5 hover:shadow-card-md transition-shadow duration-200">
      <p className="label">{label}</p>
      <p className={cn(
        'text-3xl font-bold tracking-tight',
        color ? color : (accent ? 'text-accent' : 'text-ink'),
      )}>
        {value ?? '—'}
      </p>
      {sub && <p className="text-xs text-ink-3">{sub}</p>}
    </div>
  )
}

export function PageHeader({ title, subtitle, children, backTo }) {
  return (
    <div className="flex items-start justify-between px-8 pt-8 pb-2">
      <div>
        {backTo && (
          <a href={backTo} className="text-xs text-ink-3 hover:text-ink-2 flex items-center gap-1 mb-2 transition-colors">
            ← Back
          </a>
        )}
        <h1 className="text-[26px] font-bold text-ink tracking-tight">{title}</h1>
        {subtitle && <p className="text-sm text-ink-3 mt-0.5 font-medium">{subtitle}</p>}
      </div>
      {children && (
        <div className="flex items-center gap-2 flex-wrap justify-end">
          {children}
        </div>
      )}
    </div>
  )
}

export function Spinner({ className }) {
  return (
    <div className={cn(
      'w-5 h-5 rounded-full border-2 border-border-2 border-t-accent animate-spin',
      className,
    )} />
  )
}

export function EmptyState({ icon: Icon, title, subtitle, action }) {
  return (
    <div className="flex flex-col items-center justify-center gap-3 py-20 text-center">
      {Icon && (
        <div className="w-14 h-14 rounded-2xl bg-canvas flex items-center justify-center mb-1">
          <Icon className="w-6 h-6 text-ink-3" strokeWidth={1.5} />
        </div>
      )}
      <p className="text-base font-semibold text-ink">{title}</p>
      {subtitle && <p className="text-sm text-ink-3 max-w-sm">{subtitle}</p>}
      {action && <div className="mt-2">{action}</div>}
    </div>
  )
}

export function SectionHeader({ title, subtitle, action }) {
  return (
    <div className="flex items-center justify-between mb-5">
      <div>
        <p className="font-bold text-ink text-[15px]">{title}</p>
        {subtitle && <p className="text-xs text-ink-3 mt-0.5">{subtitle}</p>}
      </div>
      {action && <div>{action}</div>}
    </div>
  )
}

export function Pill({ children, className }) {
  return (
    <span className={cn(
      'inline-flex items-center px-2.5 py-1 rounded-lg text-xs font-semibold bg-canvas text-ink-2',
      className,
    )}>
      {children}
    </span>
  )
}

export function Divider({ className }) {
  return <div className={cn('border-t border-border', className)} />
}

export function InfoRow({ label, value, mono }) {
  return (
    <div className="flex items-start justify-between gap-4 py-2.5 border-b border-border last:border-0">
      <span className="text-sm text-ink-3 font-medium shrink-0">{label}</span>
      <span className={cn('text-sm text-ink text-right', mono && 'font-mono')}>{value ?? '—'}</span>
    </div>
  )
}
