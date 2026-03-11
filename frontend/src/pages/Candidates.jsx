import { useEffect, useMemo, useState } from 'react'
import { Link, useSearchParams } from 'react-router-dom'
import { motion } from 'framer-motion'
import {
  useReactTable, getCoreRowModel, getSortedRowModel, getFilteredRowModel,
  flexRender,
} from '@tanstack/react-table'
import { ChevronUp, ChevronDown, ChevronsUpDown, Search, FlaskConical, Download } from 'lucide-react'
import { api } from '@/lib/api'
import { TierBadge, ScoreBar, PageHeader, Spinner, EmptyState } from '@/components/ui'
import { formatFloat, cn } from '@/lib/utils'

const TOOLTIP_LABELS = {
  convergence_score:  'Resilience signal across species',
  selection_score:    'Evolutionary pressure (dN/dS)',
  expression_score:   'Expression difference in resilient species',
}

function SortIcon({ sorted }) {
  if (sorted === 'asc')  return <ChevronUp className="w-3 h-3" />
  if (sorted === 'desc') return <ChevronDown className="w-3 h-3" />
  return <ChevronsUpDown className="w-3 h-3 opacity-30" />
}

const COLUMNS = [
  {
    id: 'rank', header: '#',
    cell: ({ row }) => <span className="text-text-muted text-xs">{row.index + 1}</span>,
    size: 40, enableSorting: false,
  },
  {
    accessorKey: 'gene_symbol', header: 'Gene',
    cell: ({ row, getValue }) => (
      <Link to={`/candidates/${row.original.gene_id}`} className="font-mono text-sm text-accent hover:underline font-medium">
        {getValue()}
      </Link>
    ),
    size: 100,
  },
  {
    accessorKey: 'tier', header: 'Tier',
    cell: ({ getValue }) => <TierBadge tier={getValue()} />,
    size: 80,
  },
  {
    accessorKey: 'composite_score', header: 'Score',
    cell: ({ getValue }) => (
      <div className="w-28">
        <ScoreBar value={getValue()} />
      </div>
    ),
    size: 140,
  },
  {
    accessorKey: 'convergence_score', header: () => (
      <span title={TOOLTIP_LABELS.convergence_score}>Resilience Signal</span>
    ),
    cell: ({ getValue }) => <span className="font-mono text-xs text-text-muted">{formatFloat(getValue(), 2)}</span>,
    size: 120,
  },
  {
    accessorKey: 'selection_score', header: () => (
      <span title={TOOLTIP_LABELS.selection_score}>Evol. Pressure</span>
    ),
    cell: ({ getValue }) => <span className="font-mono text-xs text-text-muted">{formatFloat(getValue(), 2)}</span>,
    size: 110,
  },
  {
    accessorKey: 'expression_score', header: () => (
      <span title={TOOLTIP_LABELS.expression_score}>Expression Δ</span>
    ),
    cell: ({ getValue }) => <span className="font-mono text-xs text-text-muted">{formatFloat(getValue(), 2)}</span>,
    size: 100,
  },
  {
    accessorKey: 'human_protein', header: 'UniProt',
    cell: ({ getValue }) => getValue()
      ? <a href={`https://www.uniprot.org/uniprot/${getValue()}`} target="_blank" rel="noreferrer"
          className="text-xs text-text-muted hover:text-accent font-mono">{getValue()}</a>
      : <span className="text-text-muted/30">—</span>,
    size: 90,
  },
  {
    id: 'actions', header: '',
    cell: ({ row }) => (
      <Link to={`/candidates/${row.original.gene_id}`} className="text-xs text-text-muted hover:text-accent">
        Details →
      </Link>
    ),
    size: 70, enableSorting: false,
  },
]

const TIER_OPTIONS = ['All', 'Validated', 'Tier1', 'Tier2', 'Tier3']
const DIRECTION_OPTIONS = [
  { value: '', label: 'All directions' },
  { value: 'loss_of_function', label: '↓ Loss-of-function' },
  { value: 'gain_of_function', label: '↑ Gain-of-function' },
  { value: 'likely_pathogenic', label: '⚠ Likely pathogenic' },
  { value: 'neutral', label: 'Neutral' },
]

export default function CandidatesPage() {
  const [searchParams] = useSearchParams()
  const speciesIdFromUrl = searchParams.get('species') || undefined
  const [allData, setAllData] = useState([])
  const [loading, setLoading] = useState(true)
  const [tierFilter, setTierFilter] = useState('All')
  const [minScore, setMinScore] = useState(0)
  const [search, setSearch] = useState('')
  const [pathwayFilter, setPathwayFilter] = useState('')
  const [pathwayOptions, setPathwayOptions] = useState({ go_terms: [], pathway_ids: [] })
  const [sorting, setSorting] = useState([{ id: 'composite_score', desc: true }])
  const [directionFilter, setDirectionFilter] = useState('')
  const [traitId, setTraitId] = useState('')
  const [traits, setTraits] = useState([])

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

  const filtered = useMemo(() => {
    return allData.filter(c => {
      if (tierFilter !== 'All' && c.tier !== tierFilter) return false
      if ((c.composite_score ?? 0) < minScore) return false
      if (search && !c.gene_symbol.toLowerCase().includes(search.toLowerCase())) return false
      return true
    })
  }, [allData, tierFilter, minScore, search])

  const table = useReactTable({
    data: filtered,
    columns: COLUMNS,
    state: { sorting },
    onSortingChange: setSorting,
    getCoreRowModel: getCoreRowModel(),
    getSortedRowModel: getSortedRowModel(),
    getFilteredRowModel: getFilteredRowModel(),
  })

  return (
    <div className="px-8 py-8 space-y-6">
      <PageHeader
        title="Candidates"
        subtitle={`${filtered.length} genes${allData.length !== filtered.length ? ` of ${allData.length}` : ''} — click any gene to see full details`}
      />

      {/* Filters */}
      <div className="flex flex-wrap items-center gap-3">
        {/* Search */}
        <div className="relative">
          <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-3.5 h-3.5 text-text-muted" />
          <input
            value={search}
            onChange={e => setSearch(e.target.value)}
            placeholder="Search gene…"
            className="pl-8 pr-3 py-2 bg-surface border border-white/8 rounded-lg text-sm text-text-primary
                       placeholder:text-text-muted focus:outline-none focus:border-accent/50 w-44"
          />
        </div>

        {/* Trait preset */}
        {traits.length > 0 && (
          <select
            value={traitId}
            onChange={e => setTraitId(e.target.value)}
            className="px-3 py-2 bg-surface border border-white/8 rounded-lg text-sm text-text-primary focus:outline-none focus:border-accent/50 min-w-[150px]"
          >
            <option value="">Default analysis</option>
            {traits.map(t => (
              <option key={t.id} value={t.id}>{t.label}</option>
            ))}
          </select>
        )}

        {/* Pathway / biological process */}
        <select
          value={pathwayFilter}
          onChange={e => setPathwayFilter(e.target.value)}
          className="px-3 py-2 bg-surface border border-white/8 rounded-lg text-sm text-text-primary focus:outline-none focus:border-accent/50 min-w-[180px]"
        >
          <option value="">All pathways</option>
          {pathwayOptions.go_terms?.slice(0, 50).map(go => (
            <option key={go} value={go}>{go}</option>
          ))}
          {pathwayOptions.pathway_ids?.slice(0, 30).map(pid => (
            <option key={pid} value={pid}>{pid}</option>
          ))}
        </select>

        {/* Variant direction filter */}
        <select
          value={directionFilter}
          onChange={e => setDirectionFilter(e.target.value)}
          className="px-3 py-2 bg-surface border border-white/8 rounded-lg text-sm text-text-primary focus:outline-none focus:border-accent/50 min-w-[160px]"
        >
          {DIRECTION_OPTIONS.map(d => (
            <option key={d.value} value={d.value}>{d.label}</option>
          ))}
        </select>

        {/* Tier filter */}
        <div className="flex rounded-lg border border-white/8 overflow-hidden">
          {TIER_OPTIONS.map(t => (
            <button
              key={t}
              onClick={() => setTierFilter(t)}
              className={cn(
                'px-3 py-2 text-xs font-medium transition-colors',
                tierFilter === t ? 'bg-accent/15 text-accent' : 'text-text-muted hover:bg-white/5',
              )}
            >
              {t}
            </button>
          ))}
        </div>

        {/* Min score slider */}
        <div className="flex items-center gap-2">
          <span className="text-xs text-text-muted">Min score</span>
          <input
            type="range" min={0} max={1} step={0.05} value={minScore}
            onChange={e => setMinScore(Number(e.target.value))}
            className="w-24 accent-accent"
          />
          <span className="font-mono text-xs text-accent w-8">{(minScore * 100).toFixed(0)}</span>
        </div>

        {/* Export CSV */}
        <a
          href={api.getCandidatesExportUrl({ tier: tierFilter !== 'All' ? tierFilter : undefined })}
          download="candidates.csv"
          className="inline-flex items-center gap-1.5 px-3 py-2 rounded-lg border border-white/8 text-xs text-text-muted hover:text-accent hover:border-accent/30 transition-colors"
        >
          <Download className="w-3.5 h-3.5" /> Export CSV
        </a>
      </div>

      {loading ? (
        <div className="flex items-center justify-center py-32"><Spinner className="w-8 h-8" /></div>
      ) : filtered.length === 0 ? (
        <EmptyState icon={FlaskConical} title="No candidates match" subtitle="Try adjusting the filters" />
      ) : (
        <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }} className="card overflow-x-auto p-0">
          <table className="w-full text-sm">
            <thead>
              {table.getHeaderGroups().map(hg => (
                <tr key={hg.id} className="border-b border-white/5">
                  {hg.headers.map(header => (
                    <th
                      key={header.id}
                      style={{ width: header.getSize() }}
                      className={cn(
                        'px-4 py-3 text-left label-muted font-medium select-none',
                        header.column.getCanSort() && 'cursor-pointer hover:text-text-primary',
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
                    'border-b border-white/3 transition-colors hover:bg-elevated',
                    i % 2 === 0 ? 'bg-surface' : 'bg-white/[0.01]',
                  )}
                >
                  {row.getVisibleCells().map(cell => (
                    <td key={cell.id} className="px-4 py-3">
                      {flexRender(cell.column.columnDef.cell, cell.getContext())}
                    </td>
                  ))}
                </tr>
              ))}
            </tbody>
          </table>
        </motion.div>
      )}
    </div>
  )
}
