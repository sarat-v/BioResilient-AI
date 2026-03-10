// Centralised API helper — calls FastAPI backend via VITE_API_BASE_URL env var
const BASE = import.meta.env.VITE_API_BASE_URL || ''

async function get(path, params = {}) {
  const url = new URL(BASE + path, window.location.origin)
  Object.entries(params).forEach(([k, v]) => v != null && url.searchParams.set(k, v))
  const res = await fetch(url)
  if (!res.ok) throw new Error(`${res.status} ${res.statusText}`)
  return res.json()
}

async function post(path, body = {}) {
  const res = await fetch(BASE + path, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  })
  if (!res.ok) throw new Error(`${res.status} ${res.statusText}`)
  return res.json()
}

export const api = {
  // Candidates
  getCandidates: (params) => get('/candidates', params),
  getCandidate: (id) => get(`/candidates/${id}`),
  getCandidatesExportUrl: (params = {}) => {
    const url = new URL((import.meta.env.VITE_API_BASE_URL || '') + '/candidates/export', window.location.origin)
    Object.entries(params).forEach(([k, v]) => v != null && url.searchParams.set(k, v))
    return url.toString()
  },
  // Species
  getSpecies: () => get('/species'),
  getSpeciesById: (id) => get(`/species/${id}`),
  // Scores
  getScores: (geneId) => get(`/scores/${geneId}`),
  // Pipeline
  getPipelineStatus: () => get('/pipeline/status'),
  startPipeline: (body) => post('/pipeline/run', body),
  stopPipeline: () => post('/pipeline/stop'),
  // Research assistant
  searchGenes: (q) => get('/research/search', { q }),
  generateNarrative: (body) => post('/research/narrative', body),
  getPathways: () => get('/research/pathways'),
  getTraits: () => get('/research/traits'),
}
