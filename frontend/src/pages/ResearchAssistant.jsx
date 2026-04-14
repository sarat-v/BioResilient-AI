import { useCallback, useEffect, useState } from 'react'
import { Link } from 'react-router-dom'
import { motion } from 'framer-motion'
import { Search, FileText, Loader2, Dna } from 'lucide-react'
import { api } from '@/lib/api'
import { PageHeader, Spinner } from '@/components/ui'
import { cn } from '@/lib/utils'

const DEBOUNCE_MS = 300

export default function ResearchAssistant() {
  const [query, setQuery] = useState('')
  const [hits, setHits] = useState([])
  const [searching, setSearching] = useState(false)
  const [selected, setSelected] = useState(null)
  const [narrative, setNarrative] = useState(null)
  const [narrativeLoading, setNarrativeLoading] = useState(false)
  const [traits, setTraits] = useState([])
  const [traitId, setTraitId] = useState('')
  const [pathwayOptions, setPathwayOptions] = useState({ go_terms: [], pathway_ids: [] })
  const [pathwayFilter, setPathwayFilter] = useState('')

  useEffect(() => {
    api.getTraits().then(setTraits).catch(() => [])
    api.getPathways().then(d => setPathwayOptions(d || { go_terms: [], pathway_ids: [] })).catch(() => {})
  }, [])

  const runSearch = useCallback((q) => {
    if (!q?.trim()) {
      setHits([])
      return
    }
    setSearching(true)
    api.searchGenes(q.trim())
      .then((data) => { setHits(data ?? []); setSearching(false) })
      .catch(() => setSearching(false))
  }, [])

  const onSearchChange = (e) => {
    const v = e.target.value
    setQuery(v)
    if (v.trim().length < 2) {
      setHits([])
      return
    }
    clearTimeout(onSearchChange._t)
    onSearchChange._t = setTimeout(() => runSearch(v), DEBOUNCE_MS)
  }

  const onGenerateNarrative = (geneId) => {
    setNarrativeLoading(true)
    setNarrative(null)
    api.generateNarrative({ gene_id: geneId, force: false })
      .then((res) => {
        setNarrative(res.narrative)
        setNarrativeLoading(false)
      })
      .catch(() => setNarrativeLoading(false))
  }

  return (
    <div className="px-8 py-8 space-y-6">
      <PageHeader
        title="Research Assistant"
        subtitle="Search genes by symbol, ID, or accession across all candidates. AI generates plain-English scientific narratives of each gene's role in disease resilience."
      />

      {traits.length > 0 && (
        <div className="flex items-center gap-2">
          <span className="text-sm text-ink-3">Trait:</span>
          <select
            value={traitId}
            onChange={e => setTraitId(e.target.value)}
            className="px-3 py-2 bg-surface border border-border rounded-lg text-sm text-ink focus:outline-none focus:border-accent/50"
          >
            <option value="">Default</option>
            {traits.map(t => (
              <option key={t.id} value={t.id}>{t.label}</option>
            ))}
          </select>
        </div>
      )}

      {/* Browse by pathway */}
      {(pathwayOptions.go_terms?.length > 0 || pathwayOptions.pathway_ids?.length > 0) && (
        <div className="card">
          <div className="flex items-center gap-2 mb-3">
            <Dna className="w-4 h-4 text-accent" />
            <p className="font-semibold text-ink">Browse by Biological Process</p>
          </div>
          <p className="text-xs text-ink-3 mb-3">Select a GO term or Reactome pathway to filter candidates.</p>
          <div className="flex items-center gap-3 flex-wrap">
            <select
              value={pathwayFilter}
              onChange={e => setPathwayFilter(e.target.value)}
              className="px-3 py-2 bg-base border border-border rounded-lg text-sm text-ink focus:outline-none focus:border-accent/50 min-w-[220px]"
            >
              <option value="">All pathways</option>
              {pathwayOptions.go_terms?.slice(0, 60).map(go => (
                <option key={go} value={go}>{go}</option>
              ))}
              {pathwayOptions.pathway_ids?.slice(0, 40).map(pid => (
                <option key={pid} value={pid}>{pid}</option>
              ))}
            </select>
            {pathwayFilter && (
              <a
                href={`/candidates?pathway=${encodeURIComponent(pathwayFilter)}`}
                className="px-4 py-2 rounded-lg bg-accent/15 text-accent border border-accent/30 hover:bg-accent/25 text-sm font-medium transition-colors"
              >
                View candidates →
              </a>
            )}
          </div>
        </div>
      )}

      {/* Search */}
      <div className="card max-w-2xl">
        <div className="relative">
          <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-ink-3" />
          <input
            value={query}
            onChange={onSearchChange}
            placeholder="Search by gene symbol (e.g. TP53), NCBI ID, or UniProt accession…"
            className="w-full pl-10 pr-4 py-3 bg-base border border-border rounded-lg text-ink
                       placeholder:text-ink-3 focus:outline-none focus:border-accent/50"
          />
          {searching && (
            <span className="absolute right-3 top-1/2 -translate-y-1/2">
              <Loader2 className="w-4 h-4 text-accent animate-spin" />
            </span>
          )}
        </div>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Results list */}
        <div className="card">
          <p className="font-semibold text-ink mb-3">Results</p>
          {hits.length === 0 && !searching && (
            <p className="text-sm text-ink-3 py-4">
              {query.trim().length < 2 ? 'Type at least 2 characters to search.' : 'No genes found.'}
            </p>
          )}
          <ul className="space-y-1 max-h-96 overflow-y-auto">
            {hits.map((g) => (
              <li key={g.gene_id}>
                <button
                  type="button"
                  onClick={() => {
                    setSelected(g)
                    setNarrative(null)
                  }}
                  className={cn(
                    'w-full text-left px-4 py-2.5 rounded-lg text-sm transition-colors',
                    selected?.gene_id === g.gene_id
                      ? 'bg-accent/15 text-accent border border-accent/30'
                      : 'hover:bg-canvas text-ink border border-transparent',
                  )}
                >
                  <span className="font-mono font-medium">{g.gene_symbol}</span>
                  {g.human_protein && (
                    <span className="ml-2 text-ink-3 text-xs">{g.human_protein}</span>
                  )}
                </button>
              </li>
            ))}
          </ul>
        </div>

        {/* Detail + narrative */}
        <div className="space-y-4">
          {selected && (
            <motion.div
              initial={{ opacity: 0, x: 8 }}
              animate={{ opacity: 1, x: 0 }}
              className="card space-y-4"
            >
              <div className="flex items-center justify-between">
                <div>
                  <p className="font-mono text-lg font-semibold text-accent">{selected.gene_symbol}</p>
                  <p className="text-xs text-ink-3">
                    {selected.human_gene_id && `NCBI ${selected.human_gene_id}`}
                    {selected.human_protein && ` · UniProt ${selected.human_protein}`}
                  </p>
                </div>
                <Link
                  to={`/candidates/${selected.gene_id}`}
                  className="text-xs text-accent hover:underline"
                >
                  Full detail →
                </Link>
              </div>
              <div>
                <div className="flex items-center gap-2 mb-2">
                  <FileText className="w-4 h-4 text-accent" />
                  <span className="font-semibold text-ink">Research Summary</span>
                </div>
                {narrativeLoading ? (
                  <div className="flex items-center gap-2 py-4 text-ink-3">
                    <Loader2 className="w-4 h-4 animate-spin" />
                    <span className="text-sm">Generating summary…</span>
                  </div>
                ) : narrative ? (
                  <div className="prose prose-invert prose-sm max-w-none text-ink whitespace-pre-wrap">
                    {narrative}
                  </div>
                ) : (
                  <button
                    type="button"
                    onClick={() => onGenerateNarrative(selected.gene_id)}
                    className="px-4 py-2 rounded-lg bg-accent/15 text-accent border border-accent/30
                               hover:bg-accent/25 text-sm font-medium transition-colors"
                  >
                    Generate Summary
                  </button>
                )}
              </div>
            </motion.div>
          )}
          {!selected && hits.length > 0 && (
            <p className="text-sm text-ink-3 py-8">Select a gene to see details and generate a summary.</p>
          )}
        </div>
      </div>
    </div>
  )
}
