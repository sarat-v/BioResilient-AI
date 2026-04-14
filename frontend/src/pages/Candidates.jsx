import { useEffect, useMemo, useState } from 'react'
import { Link, useSearchParams } from 'react-router-dom'
import { motion, AnimatePresence } from 'framer-motion'
import {
  useReactTable, getCoreRowModel, getSortedRowModel, getFilteredRowModel,
  flexRender,
} from '@tanstack/react-table'
import {
  ChevronUp, ChevronDown, ChevronsUpDown, Search, FlaskConical,
  Download, Filter, X, SlidersHorizontal,
} from 'lucide-react'
import { api } from '@/lib/api'
import { TierBadge, ScoreBar, PageHeader, Spinner, EmptyState } from '@/components/ui'
import { formatFloat, cn } from '@/lib/utils'

// ─── Column definitions ────────────────────────────────────────────────────

function SortIcon({ sorted }) {
  if (sorted === 'asc')  return <ChevronUp className="w-3 h-3 text-accent" />
  if (sorted === 'desc') return <ChevronDown className="w-3 h-3 text-accent" />
  return <ChevronsUpDown className="w-3 h-3 opacity-25" />
}

const COLUMNS = [
  {
    id: 'rank', header: '#',
    cell: ({ row }) => (
      <span className="text-xs font-mono text-ink-3 tabular-nums">{row.index + 1}</span>
    ),
    size: 40, enableSorting: false,
  },
  {
    accessorKey: 'gene_symbol', header: 'Gene',
    cell: ({ row, getValue }) => (
      <Link
        to={`/candidates/${row.original.gene_id}`}
        className="font-mono text-sm font-bold text-accent hover:underline underline-offset-2"
      >
        {getValue()}
      </Link>
    ),
    size: 100,
  },
  {
    accessorKey: 'tier', header: 'Tier',
    cell: ({ getValue }) => <TierBadge tier={getValue()} />,
    size: 90,
  },
  {
    accessorKey: 'composite_score', header: 'Score',
    cell: ({ getValue }) => (
      <div className="w-28">
        <ScoreBar value={getValue()} />
      </div>
    ),
    size: 150,
  },
  {
    accessorKey: 'convergence_score', header: 'Resilience',
    cell: ({ getValue }) => (
      <span className="text-xs font-mono text-ink-2 tabular-nums">{formatFloat(getValue(), 2)}</span>
    ),
    size: 100,
  },
  {
    accessorKey: 'selection_score', header: 'Evol. Pressure',
    cell: ({ getValue }) => (
      <span className="text-xs font-mono text-ink-2 tabular-nums">{formatFloat(getValue(), 2)}</span>
    ),
    size: 110,
  },
  {
    accessorKey: 'disease_name', header: 'Top Disease',
    cell: ({ getValue }) => {
      const v = getValue()
      return v
        ? <span className="text-xs text-ink-2 leading-tight">{v}</span>
        : <span className="text-ink-3/40 text-xs">—</span>
    },
    size: 160,
    enableSorting: false,
  },
  {
    accessorKey: 'druggability_tier', header: 'Druggability',
    cell: ({ getValue }) => {
      const v = getValue()
      const color = v === 'high' ? 'text-bio bg-bio-bg ring-bio-ring' :
                    v === 'medium' ? 'text-amber bg-amber-bg ring-amber-ring' :
                    v === 'low' ? 'text-slate bg-slate-bg ring-slate-ring' : null
      return color
        ? <span className={cn('badge text-xs capitalize', color)}>{v}</span>
        : <span className="text-ink-3/40 text-xs">—</span>
    },
    size: 110,
  },
  {
    accessorKey: 'human_protein', header: 'UniProt',
    cell: ({ getValue }) => getValue()
      ? (
        <a
          href={`https://www.uniprot.org/uniprot/${getValue()}`}
          target="_blank"
          rel="noreferrer"
          className="text-xs text-ink-3 hover:text-accent font-mono underline-offset-2 hover:underline transition-colors"
        >
          {getValue()}
        </a>
      )
      : <span className="text-ink-3/40 text-xs">—</span>,
    size: 90,
  },
  {
    id: 'actions', header: '',
    cell: ({ row }) => (
      <Link
        to={`/candidates/${row.original.gene_id}`}
        className="text-xs text-ink-3 hover:text-accent transition-colors"
      >
        Details →
      </Link>
    ),
    size: 70, enableSorting: false,
  },
]

// ─── Filters ───────────────────────────────────────────────────────────────

const TIER_OPTIONS = ['All', 'Validated', 'Tier2', 'Tier3']
// Note: Tier1 removed — no candidates matched Tier1 criteria in this run

const DIRECTION_OPTIONS = [
  { value: '', label: 'All directions' },
  { value: 'loss_of_function',  label: '↓ Loss-of-function' },
  { value: 'gain_of_function',  label: '↑ Gain-of-function' },
  { value: 'likely_pathogenic', label: '⚠ Likely pathogenic' },
  { value: 'neutral',           label: 'Neutral' },
]

// ─── Page ──────────────────────────────────────────────────────────────────

export default function CandidatesPage() {
  const [searchParams]  = useSearchParams()
  const speciesIdFromUrl = searchParams.get('species') || undefined

  const [allData, setAllData]               = useState([])
  const [loading, setLoading]               = useState(true)
  const [tierFilter, setTierFilter]         = useState('All')
  const [minScore, setMinScore]             = useState(0)
  const [search, setSearch]                 = useState('')
  const [pathwayFilter, setPathwayFilter]   = useState('')
  const [pathwayOptions, setPathwayOptions] = useState({ go_terms: [], pathway_ids: [] })
  const [sorting, setSorting]               = useState([{ id: 'composite_score', desc: true }])
  const [directionFilter, setDirectionFilter] = useState('')
  const [traitId, setTraitId]               = useState('')
  const [traits, setTraits]                 = useState([])
  const [showFilters, setShowFilters]       = useState(false)

  useEffect(() => {
    api.getPathways().then(d => setPathwayOptions(d || { go_terms: [], pathway_ids: [] })).catch(() => {})
    api.getTraits().then(setTraits).catch(() => [])
  }, [])

  useEffect(() => {
    setLoading(true)
    const params = { limit: 1000 }
    if (speciesIdFromUrl) params.species_id = speciesIdFromUrl
    if (pathwayFilter) params.pathway = pathwayFilter
    if (traitId) params.trait_id = traitId
    if (directionFilter) params.variant_direction = directionFilter
    api.getCandidates(params)
      .then(d => { setAllData(d ?? []); setLoading(false) })
      .catch(() => setLoading(false))
  }, [speciesIdFromUrl, pathwayFilter, traitId, directionFilter])

  const filtered = useMemo(() => allData.filter(c => {
    if (tierFilter !== 'All' && c.tier !== tierFilter) return false
    if ((c.composite_score ?? 0) < minScore) return false
    if (search && !c.gene_symbol.toLowerCase().includes(search.toLowerCase())) return false
    return true
  }), [allData, tierFilter, minScore, search])

  const table = useReactTable({
    data: filtered,
    columns: COLUMNS,
    state: { sorting },
    onSortingChange: setSorting,
    getCoreRowModel: getCoreRowModel(),
    getSortedRowModel: getSortedRowModel(),
    getFilteredRowModel: getFilteredRowModel(),
  })

  const activeFilterCount = [
    tierFilter !== 'All', minScore > 0, directionFilter, pathwayFilter, traitId,
  ].filter(Boolean).length

  return (
    <div className="px-8 py-8 space-y-6">
      <PageHeader
        title="Gene Candidates"
        subtitle={`${filtered.length}${allData.length !== filtered.length ? ` of ${allData.length}` : ''} genes · Validated (highest confidence), Tier 2 (moderate), Tier 3 (exploratory) · click any to see full annotation`}
      >
        <a
          href={api.getCandidatesExportUrl({ tier: tierFilter !== 'All' ? tierFilter : undefined })}
          download="candidates.csv"
          className="btn-secondary text-sm"
        >
          <Download className="w-4 h-4" /> Export CSV
        </a>
      </PageHeader>

      {/* ── Filter bar ── */}
      <div className="flex flex-wrap items-center gap-3">
        {/* Search */}
        <div className="relative">
          <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-3.5 h-3.5 text-ink-3" />
          <input
            value={search}
            onChange={e => setSearch(e.target.value)}
            placeholder="Search gene…"
            className="input pl-9 w-44"
          />
        </div>

        {/* Tier selector */}
        <div className="flex items-center rounded-xl border border-border bg-white shadow-card overflow-hidden">
          {TIER_OPTIONS.map(t => (
            <button
              key={t}
              onClick={() => setTierFilter(t)}
              className={cn(
                'px-3 py-2 text-xs font-semibold transition-all duration-150',
                tierFilter === t
                  ? 'bg-accent text-white'
                  : 'text-ink-3 hover:bg-canvas hover:text-ink-2',
              )}
            >
              {t}
            </button>
          ))}
        </div>

        {/* Advanced filter toggle */}
        <button
          onClick={() => setShowFilters(v => !v)}
          className={cn('btn-secondary relative', activeFilterCount > 0 && 'border-accent text-accent')}
        >
          <SlidersHorizontal className="w-4 h-4" />
          Filters
          {activeFilterCount > 0 && (
            <span className="absolute -top-1.5 -right-1.5 w-4 h-4 rounded-full bg-accent text-white text-[10px] font-bold flex items-center justify-center">
              {activeFilterCount}
            </span>
          )}
        </button>
      </div>

      {/* ── Advanced filter panel ── */}
      <AnimatePresence>
        {showFilters && (
          <motion.div
            initial={{ opacity: 0, height: 0 }} animate={{ opacity: 1, height: 'auto' }} exit={{ opacity: 0, height: 0 }}
            transition={{ duration: 0.2 }}
            className="overflow-hidden"
          >
            <div className="card border-accent/20">
              <div className="flex items-center justify-between mb-4">
                <p className="font-semibold text-ink text-[14px] flex items-center gap-2">
                  <Filter className="w-4 h-4 text-accent" /> Advanced Filters
                </p>
                <button onClick={() => setShowFilters(false)} className="btn-icon w-7 h-7">
                  <X className="w-3.5 h-3.5" />
                </button>
              </div>

              <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-4">
                {/* Trait preset */}
                {traits.length > 0 && (
                  <div>
                    <p className="label mb-1.5">Research Phenotype</p>
                    <select value={traitId} onChange={e => setTraitId(e.target.value)} className="select">
                      <option value="">Default analysis</option>
                      {traits.map(t => <option key={t.id} value={t.id}>{t.label}</option>)}
                    </select>
                  </div>
                )}

                {/* Pathway */}
                <div>
                  <p className="label mb-1.5">Biological Process</p>
                  <select value={pathwayFilter} onChange={e => setPathwayFilter(e.target.value)} className="select">
                    <option value="">All pathways</option>
                    {pathwayOptions.go_terms?.slice(0, 50).map(go => <option key={go} value={go}>{go}</option>)}
                    {pathwayOptions.pathway_ids?.slice(0, 30).map(pid => <option key={pid} value={pid}>{pid}</option>)}
                  </select>
                </div>

                {/* Variant direction */}
                <div>
                  <p className="label mb-1.5">Variant Direction</p>
                  <select value={directionFilter} onChange={e => setDirectionFilter(e.target.value)} className="select">
                    {DIRECTION_OPTIONS.map(d => <option key={d.value} value={d.value}>{d.label}</option>)}
                  </select>
                </div>

                {/* Min score */}
                <div>
                  <p className="label mb-1.5">Min Score: <span className="text-accent">{(minScore * 100).toFixed(0)}</span></p>
                  <input
                    type="range" min={0} max={1} step={0.05} value={minScore}
                    onChange={e => setMinScore(Number(e.target.value))}
                    className="w-full accent-accent"
                  />
                </div>
              </div>

              {activeFilterCount > 0 && (
                <button
                  onClick={() => { setTierFilter('All'); setMinScore(0); setDirectionFilter(''); setPathwayFilter(''); setTraitId('') }}
                  className="mt-4 text-xs text-danger hover:underline font-medium"
                >
                  Clear all filters
                </button>
              )}
            </div>
          </motion.div>
        )}
      </AnimatePresence>

      {/* ── Table ── */}
      {loading ? (
        <div className="flex items-center justify-center py-40">
          <Spinner className="w-8 h-8" />
        </div>
      ) : filtered.length === 0 ? (
        <EmptyState
          icon={FlaskConical}
          title="No candidates match your filters"
          subtitle="Try adjusting the tier filter or minimum score"
          action={
            <button onClick={() => { setTierFilter('All'); setMinScore(0); setSearch('') }} className="btn-secondary text-sm">
              Clear filters
            </button>
          }
        />
      ) : (
        <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }} className="card p-0 overflow-hidden">
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                {table.getHeaderGroups().map(hg => (
                  <tr key={hg.id} className="border-b border-border bg-canvas">
                    {hg.headers.map(header => (
                      <th
                        key={header.id}
                        style={{ width: header.getSize() }}
                        className={cn(
                          'px-4 py-3 text-left label font-semibold whitespace-nowrap select-none',
                          header.column.getCanSort() && 'cursor-pointer hover:text-ink-2',
                        )}
                        onClick={header.column.getToggleSortingHandler()}
                      >
                        <span className="flex items-center gap-1">
                          {flexRender(header.column.columnDef.header, header.getContext())}
                          {header.column.getCanSort() && (
                            <SortIcon sorted={header.column.getIsSorted()} />
                          )}
                        </span>
                      </th>
                    ))}
                  </tr>
                ))}
              </thead>
              <tbody>
                {table.getRowModel().rows.map((row, i) => (
                  <tr
                    key={row.id}
                    className={cn(
                      'border-b border-border/60 transition-colors duration-75 hover:bg-accent-light/40',
                      i % 2 === 0 ? 'bg-white' : 'bg-canvas/40',
                    )}
                  >
                    {row.getVisibleCells().map(cell => (
                      <td key={cell.id} className="px-4 py-2.5">
                        {flexRender(cell.column.columnDef.cell, cell.getContext())}
                      </td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>

          <div className="px-4 py-3 border-t border-border bg-canvas flex items-center justify-between">
            <p className="text-xs text-ink-3">
              Showing <span className="font-semibold text-ink">{filtered.length}</span> candidates
            </p>
            <p className="text-xs text-ink-3">
              Click any gene to view full Phase 1 + 2 annotation
            </p>
          </div>
        </motion.div>
      )}
    </div>
  )
}
