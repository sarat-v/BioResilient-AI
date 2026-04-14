import { useState } from 'react'
import { useNavigate } from 'react-router-dom'
import { motion } from 'framer-motion'
import { Dna, AlertCircle } from 'lucide-react'
import { useAuth } from '@/context/AuthContext'

const BASE = import.meta.env.VITE_API_BASE_URL || ''

export default function Login() {
  const { login } = useAuth()
  const navigate = useNavigate()

  const [email, setEmail] = useState('')
  const [password, setPassword] = useState('')
  const [error, setError] = useState('')
  const [loading, setLoading] = useState(false)

  async function handleSubmit(e) {
    e.preventDefault()
    setError('')
    setLoading(true)

    try {
      const res = await fetch(`${BASE}/auth/login`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ email, password }),
      })

      if (!res.ok) {
        const data = await res.json().catch(() => ({}))
        throw new Error(data.detail || 'Login failed')
      }

      const data = await res.json()
      login(data.access_token, data.user)
      navigate('/', { replace: true })
    } catch (err) {
      setError(err.message)
    } finally {
      setLoading(false)
    }
  }

  return (
    <div className="min-h-screen bg-canvas flex items-center justify-center px-4">
      {/* Background gradients matching the app */}
      <div className="fixed inset-0 pointer-events-none">
        <div
          className="absolute inset-0"
          style={{
            backgroundImage: `
              radial-gradient(ellipse 800px 600px at 0% 0%, rgba(27,69,212,0.035) 0%, transparent 70%),
              radial-gradient(ellipse 600px 800px at 100% 100%, rgba(124,58,237,0.025) 0%, transparent 70%)
            `,
          }}
        />
      </div>

      <motion.div
        initial={{ opacity: 0, y: 12 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.28, ease: 'easeOut' }}
        className="relative w-full max-w-[380px]"
      >
        {/* Logo */}
        <div className="flex items-center gap-3 justify-center mb-8">
          <div
            className="w-10 h-10 rounded-xl flex items-center justify-center shrink-0"
            style={{ background: 'linear-gradient(135deg, #1B45D4, #7C3AED)' }}
          >
            <Dna className="w-5 h-5 text-white" strokeWidth={2} />
          </div>
          <div>
            <p className="text-[16px] font-bold text-ink tracking-tight leading-none">BioResilient</p>
            <p className="text-[12px] text-ink-3 font-medium mt-0.5 leading-none">AI Platform</p>
          </div>
        </div>

        {/* Card */}
        <div className="card">
          <h1 className="text-[20px] font-bold text-ink mb-1">Sign in</h1>
          <p className="text-sm text-ink-3 mb-6">Access the research pipeline and candidate data</p>

          <form onSubmit={handleSubmit} className="flex flex-col gap-4">
            <div>
              <label className="label mb-1.5 block">Email</label>
              <input
                type="email"
                required
                autoFocus
                autoComplete="email"
                value={email}
                onChange={e => setEmail(e.target.value)}
                placeholder="you@lab.org"
                className="w-full px-3.5 py-2.5 rounded-xl border border-border bg-canvas text-sm text-ink
                           placeholder:text-ink-3 focus:outline-none focus:ring-2 focus:ring-accent/30 focus:border-accent
                           transition-all duration-150"
              />
            </div>

            <div>
              <label className="label mb-1.5 block">Password</label>
              <input
                type="password"
                required
                autoComplete="current-password"
                value={password}
                onChange={e => setPassword(e.target.value)}
                placeholder="••••••••"
                className="w-full px-3.5 py-2.5 rounded-xl border border-border bg-canvas text-sm text-ink
                           placeholder:text-ink-3 focus:outline-none focus:ring-2 focus:ring-accent/30 focus:border-accent
                           transition-all duration-150"
              />
            </div>

            {error && (
              <div className="flex items-center gap-2 px-3.5 py-2.5 rounded-xl bg-red-50 border border-red-200">
                <AlertCircle className="w-4 h-4 text-red-500 shrink-0" />
                <p className="text-sm text-red-700">{error}</p>
              </div>
            )}

            <button
              type="submit"
              disabled={loading}
              className="btn-primary justify-center mt-1 disabled:opacity-60 disabled:cursor-not-allowed"
            >
              {loading ? (
                <span className="flex items-center gap-2">
                  <svg className="w-4 h-4 animate-spin" viewBox="0 0 24 24" fill="none">
                    <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                    <path className="opacity-75" fill="currentColor"
                      d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                  </svg>
                  Signing in…
                </span>
              ) : 'Sign in'}
            </button>
          </form>
        </div>

        <p className="text-center text-xs text-ink-3 mt-5">
          Access is restricted. Contact your administrator to request an account.
        </p>
      </motion.div>
    </div>
  )
}
