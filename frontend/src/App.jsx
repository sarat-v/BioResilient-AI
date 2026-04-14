import { Routes, Route, Navigate } from 'react-router-dom'
import { AnimatePresence, motion } from 'framer-motion'
import { AuthProvider, useAuth } from '@/context/AuthContext'
import Layout from '@/components/Layout'
import Login from '@/pages/Login'
import Dashboard from '@/pages/Dashboard'
import PipelinePage from '@/pages/Pipeline'
import CandidatesPage from '@/pages/Candidates'
import GeneDetail from '@/pages/GeneDetail'
import SpeciesPage from '@/pages/Species'
import ResearchAssistant from '@/pages/ResearchAssistant'
import PathwayConvergencePage from '@/pages/PathwayConvergence'

const PAGE_TRANSITION = {
  initial: { opacity: 0, x: 8 },
  animate: { opacity: 1, x: 0 },
  exit: { opacity: 0, x: -8 },
  transition: { duration: 0.2 },
}

function Animated({ children }) {
  return <motion.div {...PAGE_TRANSITION}>{children}</motion.div>
}

function AppInner() {
  const { isAuthenticated } = useAuth()

  return (
    <Routes>
      <Route
        path="/login"
        element={isAuthenticated ? <Navigate to="/" replace /> : <Login />}
      />
      <Route
        path="/*"
        element={
          isAuthenticated ? (
            <Layout>
              <AnimatePresence mode="wait">
                <Routes>
                  <Route path="/" element={<Animated><Dashboard /></Animated>} />
                  <Route path="/pipeline" element={<Animated><PipelinePage /></Animated>} />
                  <Route path="/candidates" element={<Animated><CandidatesPage /></Animated>} />
                  <Route path="/candidates/:id" element={<Animated><GeneDetail /></Animated>} />
                  <Route path="/research" element={<Animated><ResearchAssistant /></Animated>} />
                  <Route path="/pathways" element={<Animated><PathwayConvergencePage /></Animated>} />
                  <Route path="/species" element={<Animated><SpeciesPage /></Animated>} />
                </Routes>
              </AnimatePresence>
            </Layout>
          ) : (
            <Navigate to="/login" replace />
          )
        }
      />
    </Routes>
  )
}

export default function App() {
  return (
    <AuthProvider>
      <AppInner />
    </AuthProvider>
  )
}
