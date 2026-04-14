import { NavLink } from 'react-router-dom'
import { motion } from 'framer-motion'
import {
  LayoutDashboard, FlaskConical, Dna, TreePine, Activity,
  FileText, Route, Zap, LogOut, User,
} from 'lucide-react'
import { cn } from '@/lib/utils'
import { useAuth } from '@/context/AuthContext'

const NAV = [
  { to: '/',           icon: LayoutDashboard, label: 'Dashboard'  },
  { to: '/pipeline',   icon: Activity,        label: 'Pipeline'   },
  { to: '/candidates', icon: FlaskConical,    label: 'Candidates' },
  { to: '/research',   icon: FileText,        label: 'Research'   },
  { to: '/pathways',   icon: Route,           label: 'Pathways'   },
  { to: '/species',    icon: TreePine,        label: 'Species'    },
]

export default function Layout({ children }) {
  const { user, logout } = useAuth()
  return (
    <div className="flex h-screen overflow-hidden bg-canvas">
      {/* Sidebar */}
      <aside className="w-[220px] shrink-0 flex flex-col bg-white border-r border-border">
        {/* Logo */}
        <div className="px-5 pt-7 pb-6">
          <div className="flex items-center gap-3">
            <div
              className="w-9 h-9 rounded-xl flex items-center justify-center shrink-0"
              style={{ background: 'linear-gradient(135deg, #1B45D4, #7C3AED)' }}
            >
              <Dna className="w-[18px] h-[18px] text-white" strokeWidth={2} />
            </div>
            <div>
              <p className="text-[14px] font-bold text-ink tracking-tight leading-none">BioResilient</p>
              <p className="text-[11px] text-ink-3 font-medium mt-0.5 leading-none">AI Platform</p>
            </div>
          </div>
        </div>

        <div className="mx-4 mb-3 border-t border-border" />

        {/* Nav */}
        <nav className="flex-1 px-3 py-1 space-y-0.5 overflow-y-auto">
          {NAV.map(({ to, icon: Icon, label }) => (
            <NavLink
              key={to}
              to={to}
              end={to === '/'}
              className={({ isActive }) =>
                cn('nav-item group', isActive && 'nav-active')
              }
            >
              {({ isActive }) => (
                <>
                  <Icon
                    className={cn('w-4 h-4 shrink-0 transition-colors', isActive ? 'text-accent' : 'text-ink-3 group-hover:text-ink-2')}
                    strokeWidth={isActive ? 2.2 : 1.8}
                  />
                  <span className="flex-1">{label}</span>
                  {isActive && (
                    <motion.div
                      layoutId="nav-dot"
                      className="w-1.5 h-1.5 rounded-full bg-accent"
                    />
                  )}
                </>
              )}
            </NavLink>
          ))}
        </nav>

        {/* Bottom info */}
        <div className="px-4 pb-6 pt-3">
          <div className="mx-1 mb-3 border-t border-border" />

          {/* Pipeline badge */}
          <div className="flex items-center gap-2.5 px-2 py-2 rounded-xl bg-canvas mb-1.5">
            <div className="w-7 h-7 rounded-lg bg-accent-light flex items-center justify-center shrink-0">
              <Zap className="w-3.5 h-3.5 text-accent" />
            </div>
            <div className="min-w-0">
              <p className="text-[11px] font-semibold text-ink truncate">Phase 1 + 2 Pipeline</p>
              <p className="text-[10px] text-ink-3">30 steps · AWS Batch</p>
            </div>
          </div>

          {/* Logged-in user + logout */}
          <button
            onClick={logout}
            className="w-full flex items-center gap-2.5 px-2 py-2 rounded-xl hover:bg-surface-2 group transition-colors duration-150 text-left"
            title="Sign out"
          >
            <div className="w-7 h-7 rounded-lg bg-surface-2 border border-border flex items-center justify-center shrink-0">
              <User className="w-3.5 h-3.5 text-ink-3" />
            </div>
            <div className="min-w-0 flex-1">
              <p className="text-[11px] font-semibold text-ink truncate">{user?.name}</p>
              <p className="text-[10px] text-ink-3 truncate">{user?.email}</p>
            </div>
            <LogOut className="w-3.5 h-3.5 text-ink-3 opacity-0 group-hover:opacity-100 transition-opacity shrink-0" />
          </button>
        </div>
      </aside>

      {/* Main content area */}
      <main className="flex-1 overflow-y-auto min-w-0">
        <div className="page-enter">
          {children}
        </div>
      </main>
    </div>
  )
}
