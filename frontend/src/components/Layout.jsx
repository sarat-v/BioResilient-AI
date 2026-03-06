import { NavLink } from 'react-router-dom'
import { LayoutDashboard, FlaskConical, Dna, TreePine, Activity } from 'lucide-react'
import { cn } from '@/lib/utils'

const NAV = [
  { to: '/',          icon: LayoutDashboard, label: 'Dashboard'  },
  { to: '/pipeline',  icon: Activity,        label: 'Pipeline'   },
  { to: '/candidates',icon: FlaskConical,    label: 'Candidates' },
  { to: '/species',   icon: TreePine,        label: 'Species'    },
]

export default function Layout({ children }) {
  return (
    <div className="flex h-screen overflow-hidden">
      {/* Sidebar */}
      <aside className="w-56 shrink-0 flex flex-col bg-surface border-r border-white/5">
        {/* Logo */}
        <div className="px-5 py-6 border-b border-white/5">
          <div className="flex items-center gap-2.5">
            <div className="w-8 h-8 rounded-lg bg-gradient-accent flex items-center justify-center">
              <Dna className="w-4 h-4 text-white" />
            </div>
            <div>
              <p className="text-sm font-semibold text-text-primary leading-tight">BioResilient</p>
              <p className="text-[10px] text-text-muted leading-tight">AI Platform</p>
            </div>
          </div>
        </div>

        {/* Nav */}
        <nav className="flex-1 px-3 py-4 space-y-1">
          {NAV.map(({ to, icon: Icon, label }) => (
            <NavLink
              key={to}
              to={to}
              end={to === '/'}
              className={({ isActive }) =>
                cn(
                  'flex items-center gap-3 px-3 py-2.5 rounded-lg text-sm font-medium transition-all duration-150',
                  isActive
                    ? 'bg-accent/10 text-accent border border-accent/20'
                    : 'text-text-muted hover:text-text-primary hover:bg-white/5',
                )
              }
            >
              <Icon className="w-4 h-4 shrink-0" />
              {label}
            </NavLink>
          ))}
        </nav>

        <div className="px-5 py-4 border-t border-white/5">
          <p className="text-[10px] text-text-muted">Phase 1 + 2 Pipeline</p>
          <p className="text-[10px] text-text-muted opacity-50">v1.0</p>
        </div>
      </aside>

      {/* Main content */}
      <main className="flex-1 overflow-y-auto">
        {children}
      </main>
    </div>
  )
}
