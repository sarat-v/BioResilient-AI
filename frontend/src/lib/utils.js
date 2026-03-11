import { clsx } from 'clsx'

export function cn(...inputs) {
  return clsx(inputs)
}

export function formatScore(v) {
  if (v == null) return '—'
  return (v * 100).toFixed(0)
}

export function formatFloat(v, decimals = 3) {
  if (v == null) return '—'
  return Number(v).toFixed(decimals)
}

export function tierBadgeClass(tier) {
  if (tier === 'Validated') return 'badge-validated'
  if (tier === 'Tier1') return 'badge-tier1'
  if (tier === 'Tier2') return 'badge-tier2'
  return 'badge-tier3'
}

export function stepStatusColor(status) {
  if (status === 'complete') return 'text-success'
  if (status === 'running') return 'text-accent'
  if (status === 'failed') return 'text-danger'
  return 'text-text-muted'
}

export function relativeTime(iso) {
  if (!iso) return null
  const d = new Date(iso)
  const diff = Math.floor((Date.now() - d.getTime()) / 1000)
  if (diff < 60) return `${diff}s ago`
  if (diff < 3600) return `${Math.floor(diff / 60)}m ago`
  if (diff < 86400) return `${Math.floor(diff / 3600)}h ago`
  return d.toLocaleDateString()
}
