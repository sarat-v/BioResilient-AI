import { useEffect, useState } from 'react'
import { Link } from 'react-router-dom'
import { motion } from 'framer-motion'
import { TreePine, ExternalLink } from 'lucide-react'
import { api } from '@/lib/api'
import { PageHeader, Spinner, EmptyState } from '@/components/ui'

// Key phenotype → display label + colour
const PHENOTYPE_STYLE = {
  cancer_resistance:  { label: 'Cancer resistance',  cls: 'bg-success/10 text-emerald-300 ring-success/20' },
  longevity:          { label: 'Longevity',           cls: 'bg-accent/10 text-cyan-300 ring-accent/20' },
  regeneration:       { label: 'Regeneration',        cls: 'bg-purple-400/10 text-purple-300 ring-purple-500/20' },
  hibernation:        { label: 'Hibernation',         cls: 'bg-blue-400/10 text-blue-300 ring-blue-500/20' },
  viral_resistance:   { label: 'Viral resistance',    cls: 'bg-orange-400/10 text-orange-300 ring-orange-500/20' },
  ischaemia_resistance: { label: 'Ischaemia resistance', cls: 'bg-red-400/10 text-red-300 ring-red-400/20' },
  baseline:           { label: 'Baseline (human)',    cls: 'bg-white/5 text-text-muted ring-white/10' },
}

// Cute animal emojis as stand-ins for photos
const SPECIES_EMOJI = {
  naked_mole_rat:       '🐭',
  bowhead_whale:        '🐋',
  axolotl:              '🦎',
  ground_squirrel:      '🐿️',
  little_brown_bat:     '🦇',
  damaraland_mole_rat:  '🐀',
  greenland_shark:      '🦈',
  african_elephant:     '🐘',
  mouse_lemur:          '🐒',
  spiny_mouse:          '🐁',
  rougheye_rockfish:    '🐟',
  human:                '🧬',
}

function PhenotypePill({ phenotype }) {
  const style = PHENOTYPE_STYLE[phenotype] ?? { label: phenotype, cls: 'bg-white/5 text-text-muted ring-white/10' }
  return (
    <span className={`inline-flex items-center px-2 py-0.5 rounded-full text-[11px] font-medium ring-1 ${style.cls}`}>
      {style.label}
    </span>
  )
}

function SpeciesCard({ species, index }) {
  return (
    <motion.div
      initial={{ opacity: 0, y: 20 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ delay: index * 0.06, duration: 0.4, ease: 'easeOut' }}
      className="card hover:border-white/10 hover:bg-elevated transition-all duration-200 flex flex-col gap-4"
    >
      {/* Header */}
      <div className="flex items-start gap-4">
        <div className="w-14 h-14 rounded-xl bg-gradient-to-br from-accent/10 to-accent-warm/10 flex items-center justify-center text-3xl shrink-0 border border-white/5">
          {SPECIES_EMOJI[species.id] ?? '🔬'}
        </div>
        <div className="flex-1 min-w-0">
          <p className="font-semibold text-text-primary capitalize leading-tight">
            {species.id.replace(/_/g, ' ')}
          </p>
          <p className="text-xs text-text-muted italic mt-0.5 truncate">{species.scientific_name}</p>
          <p className="text-xs text-text-muted/60 font-mono mt-0.5">TaxID {species.taxid}</p>
        </div>
      </div>

      {/* Phenotypes */}
      <div className="flex flex-wrap gap-1.5">
        {(species.phenotypes ?? []).map(p => (
          <PhenotypePill key={p} phenotype={p} />
        ))}
      </div>

      {/* Details row */}
      <div className="flex items-center justify-between pt-2 border-t border-white/5 text-xs text-text-muted">
        <span className="capitalize">{species.lineage_group ?? 'Unknown lineage'}</span>
        {species.genome_assembly && (
          <a
            href={`https://www.ncbi.nlm.nih.gov/datasets/genome/${species.genome_assembly}/`}
            target="_blank"
            rel="noreferrer"
            className="flex items-center gap-1 hover:text-accent transition-colors"
          >
            {species.genome_assembly} <ExternalLink className="w-3 h-3" />
          </a>
        )}
      </div>

      {/* Link to candidates for this species */}
      <Link
        to={`/candidates?species=${species.id}`}
        className="text-xs text-accent/60 hover:text-accent transition-colors flex items-center gap-1"
      >
        View candidate genes → 
      </Link>
    </motion.div>
  )
}

// Group species by lineage
function groupByLineage(speciesList) {
  const groups = {}
  for (const s of speciesList) {
    const g = s.lineage_group ?? 'Other'
    if (!groups[g]) groups[g] = []
    groups[g].push(s)
  }
  return groups
}

export default function SpeciesPage() {
  const [speciesList, setSpeciesList] = useState([])
  const [loading, setLoading] = useState(true)

  useEffect(() => {
    api.getSpecies()
      .then(d => { setSpeciesList(d ?? []); setLoading(false) })
      .catch(() => setLoading(false))
  }, [])

  const groups = groupByLineage(speciesList)
  const nonHuman = speciesList.filter(s => s.id !== 'human')

  return (
    <div className="px-8 py-8 space-y-8">
      <PageHeader
        title="Species Registry"
        subtitle={`${nonHuman.length} resilient species compared against the human baseline`}
      />

      {loading ? (
        <div className="flex items-center justify-center py-32"><Spinner className="w-8 h-8" /></div>
      ) : speciesList.length === 0 ? (
        <EmptyState
          icon={TreePine}
          title="No species data yet"
          subtitle="Run the pipeline to load species data"
        />
      ) : (
        Object.entries(groups).sort(([a], [b]) => a === 'Primates' ? 1 : a.localeCompare(b)).map(([lineage, members]) => (
          <div key={lineage}>
            <div className="flex items-center gap-3 mb-4">
              <p className="text-sm font-semibold text-text-muted uppercase tracking-wider">{lineage}</p>
              <div className="flex-1 h-px bg-white/5" />
              <span className="text-xs text-text-muted/50">{members.length}</span>
            </div>
            <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4 gap-4">
              {members.map((s, i) => (
                <SpeciesCard key={s.id} species={s} index={i} />
              ))}
            </div>
          </div>
        ))
      )}
    </div>
  )
}
