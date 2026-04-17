// Centralised API helper — calls FastAPI backend via VITE_API_BASE_URL env var
const BASE = import.meta.env.VITE_API_BASE_URL || ''

function _authHeader() {
  const token = localStorage.getItem('br_token')
  return token ? { Authorization: `Bearer ${token}` } : {}
}

function _handleUnauthorized() {
  localStorage.removeItem('br_token')
  localStorage.removeItem('br_user')
  window.location.href = '/login'
}

async function get(path, params = {}) {
  const url = new URL(BASE + path, window.location.origin)
  Object.entries(params).forEach(([k, v]) => v != null && url.searchParams.set(k, v))
  const res = await fetch(url, { headers: _authHeader() })
  if (res.status === 401) { _handleUnauthorized(); throw new Error('Session expired') }
  if (!res.ok) throw new Error(`${res.status} ${res.statusText}`)
  return res.json()
}

async function post(path, body = {}) {
  const res = await fetch(BASE + path, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json', ..._authHeader() },
    body: JSON.stringify(body),
  })
  if (res.status === 401) { _handleUnauthorized(); throw new Error('Session expired') }
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
  // body: { resume_from, until_step?, dry_run, phenotype, species_ids }
  stopPipeline: () => post('/pipeline/stop'),
  // Research assistant
  searchGenes: (q) => get('/research/search', { q }),
  generateNarrative: (body) => post('/research/narrative', body),
  getPathways: () => get('/research/pathways'),
  getTraits: () => get('/research/traits'),
  getPathwayConvergence: (params) => get('/research/pathway-convergence', params),
}
